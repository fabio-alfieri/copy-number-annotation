# rm(list = ls())
gc(full = TRUE)

library(tidyverse)
library(pROC)
library(PRROC)
library(ggpubr)

# -----------------------------
# Settings
# -----------------------------

input_dir <- "annotation-matrix"
output_base <- "aggregateBins_all_files_lessPointsANDmean"

if (!dir.exists(output_base)) {
  dir.create(output_base)
}

cutoffs <- c(0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9)

genome_build <- "hg19"

# -----------------------------
# Load cytobands / centromeres
# -----------------------------

cytoband <- read_tsv(
  "cytoBand.txt.gz",
  col_names = c("chr", "start", "end", "band", "gieStain"),
  show_col_types = FALSE
)

centromeres <- cytoband %>%
  filter(gieStain == "acen") %>%
  group_by(chr) %>%
  summarise(
    cen_start = min(start),
    cen_end = max(end),
    .groups = "drop"
  )

# -----------------------------
# Find all annotation files
# -----------------------------

annotation_files <- list.files(
  input_dir,
  pattern = "^annotation_(ampl|del)_.+\\.tsv$",
  full.names = TRUE
)

annotation_files


# -----------------------------
# Load Jubran et al. ROC / PR curves
# -----------------------------

# we run Kubran et al. ML tool modified to extract Gain_versus_Neutral__curve_points_ROC_PR.csv and Loss_versus_Neutral__curve_points_ROC_PR.csv

uri_curve_files <- list(
  ampl = "../Uri-tool/AneuploidyML/article_results/ROC_PR_export/Gain_versus_Neutral__curve_points_ROC_PR.csv",
  del  = "../Uri-tool/AneuploidyML/article_results/ROC_PR_export/Loss_versus_Neutral__curve_points_ROC_PR.csv"
)

get_jubran_xgboost_curve <- function(cna_type, curve_type = c("ROC", "PR")) {
  
  curve_type <- match.arg(curve_type)
  curve_file <- uri_curve_files[[cna_type]]
  
  if (is.null(curve_file) || !file.exists(curve_file)) {
    warning("Jubran curve file not found for cna_type = ", cna_type)
    return(NULL)
  }
  
  curve_df <- read_csv(curve_file, show_col_types = FALSE) %>%
    filter(curve == curve_type)
  
  median_fold_df <- curve_df %>%
    distinct(comparison, model, fold, auc) %>%
    group_by(comparison, model) %>%
    mutate(
      mean_auc = mean(auc, na.rm = TRUE),
      distance_from_mean = abs(auc - mean_auc)
    ) %>%
    slice_min(distance_from_mean, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  plot_df <- curve_df %>%
    inner_join(
      median_fold_df %>% select(comparison, model, fold, auc, mean_auc),
      by = c("comparison", "model", "fold", "auc")
    ) %>%
    filter(model == "XGBoost")
  
  if (curve_type == "ROC") {
    plot_df <- plot_df %>%
      transmute(
        fpr = x,
        tpr = y,
        method = "Jubran et al., 2024",
        auc = auc
      )
  }
  
  if (curve_type == "PR") {
    plot_df <- plot_df %>%
      transmute(
        recall = x,
        precision_interp = y,
        method = "Jubran et al., 2024",
        auc = auc
      )
  }
  
  return(plot_df)
}

# -----------------------------
# Run ROC curve PENNE + Jubran
# -----------------------------

run_aggregate_roc_pr <- function(file_path) {
  
  message("Processing: ", file_path)
  
  file_name <- basename(file_path)
  
  # Parse CNA type and scale from filename
  # Example: annotation_del_Small-scale.tsv
  cna_type <- str_match(file_name, "^annotation_(ampl|del)_(.+)\\.tsv$")[, 2]
  scale_name <- str_match(file_name, "^annotation_(ampl|del)_(.+)\\.tsv$")[, 3]
  
  safe_scale <- scale_name %>%
    str_replace_all("-", "") %>%
    str_replace_all(" ", "_")
  
  out_folder <- file.path(
    output_base,
    paste0(safe_scale, "_", cna_type, "_aggregateBins")
  )
  
  if (!dir.exists(out_folder)) {
    dir.create(out_folder, recursive = TRUE)
  }
  
  # -----------------------------
  # Read data
  # -----------------------------
  
  arm.level <- read.csv(file_path, sep = "\t")
  
  required_cols <- c("chr", "start", "end", "type", "observed", "prediction")
  
  missing_cols <- setdiff(required_cols, colnames(arm.level))
  
  if (length(missing_cols) > 0) {
    warning(
      "Skipping ", file_name,
      " because these columns are missing: ",
      paste(missing_cols, collapse = ", ")
    )
    return(NULL)
  }
  
  df <- arm.level %>%
    filter(!is.na(observed), !is.na(prediction)) %>%
    mutate(
      chr = as.character(chr),
      chr = ifelse(str_detect(chr, "^chr"), chr, paste0("chr", chr)),
      midpoint = (start + end) / 2,
      bin_width = end - start + 1
    )
  
  # -----------------------------
  # Assign bins to chromosome arms
  # -----------------------------
  
  df_arm <- df %>%
    left_join(centromeres, by = "chr") %>%
    mutate(
      arm = case_when(
        midpoint < cen_start ~ "p",
        midpoint > cen_end ~ "q",
        TRUE ~ NA_character_
      ),
      chr_arm = paste0(chr, arm)
    ) %>%
    filter(!is.na(arm))
  
  if (nrow(df_arm) == 0) {
    warning("Skipping ", file_name, ": no bins could be assigned to chromosome arms.")
    return(NULL)
  }
  
  # -----------------------------
  # Aggregate predictions / observations
  # -----------------------------
  
  arm_agg <- df_arm %>%
    group_by(type, chr, arm, chr_arm) %>%
    summarise(
      start = min(start, na.rm = TRUE),
      end = max(end, na.rm = TRUE),
      n_bins = n(),
      total_bp = sum(bin_width, na.rm = TRUE),
      
      # observed_arm = weighted.mean(observed, w = bin_width, na.rm = TRUE),
      # prediction_arm = weighted.mean(prediction, w = bin_width, na.rm = TRUE),
      
      observed_arm = mean(observed, na.rm = TRUE),
      prediction_arm = mean(prediction, na.rm = TRUE),
      
      residual_arm = observed_arm - prediction_arm,
      
      .groups = "drop"
    )
  
  write.csv(
    arm_agg,
    file.path(out_folder, paste0("aggregated_arm_dataset_", cna_type, "_", safe_scale, ".csv")),
    row.names = FALSE
  )
  
  # -----------------------------
  # Observed vs predicted plot
  # -----------------------------
  
  p_scatter <- ggplot(arm_agg, aes(x = observed_arm, y = prediction_arm)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~ type, scales = "free") +
    theme_classic(base_size = 13) +
    labs(
      x = "Observed arm-level CNA frequency",
      y = "Predicted arm-level CNA frequency",
      title = paste0("Observed vs predicted - ", cna_type, " - ", scale_name)
    )
  
  ggsave(
    filename = file.path(out_folder, paste0("Observed_vs_predicted_", cna_type, "_", safe_scale, ".pdf")),
    plot = p_scatter,
    width = 8,
    height = 5
  )
  
  # -----------------------------
  # ROC / PR across cutoffs
  # -----------------------------
  
  roc_all <- list()
  pr_all <- list()
  summary_all <- list()
  
  thin_roc_points <- function(roc_df, n_points = 15) {
    
    # Always work on unique FPR/TPR points
    roc_clean <- roc_df %>%
      arrange(fpr, tpr) %>%
      group_by(fpr) %>%
      summarise(
        tpr = max(tpr, na.rm = TRUE),
        threshold = threshold[which.max(tpr)],
        specificity = specificity[which.max(tpr)],
        sensitivity = sensitivity[which.max(tpr)],
        cutoff_quantile = first(cutoff_quantile),
        roc_auc = first(roc_auc),
        cna_type = first(cna_type),
        scale = first(scale),
        file = first(file),
        .groups = "drop"
      ) %>%
      arrange(fpr)
    
    # Define approximately equally spaced FPR positions
    target_fpr <- seq(0, 1, length.out = n_points)
    
    # For each target FPR, keep the closest actual ROC point
    keep_idx <- map_int(target_fpr, function(x) {
      which.min(abs(roc_clean$fpr - x))
    }) %>%
      unique()
    
    roc_thinned <- roc_clean[keep_idx, ] %>%
      arrange(fpr, tpr)
    
    # Make sure endpoints are retained
    endpoints <- roc_clean %>%
      filter(
        row_number() == 1 | row_number() == n()
      )
    
    roc_thinned <- bind_rows(endpoints, roc_thinned) %>%
      distinct(fpr, tpr, .keep_all = TRUE) %>%
      arrange(fpr, tpr)
    
    return(roc_thinned)
  }
  
  for (observed_arbitrary_cutoff in cutoffs) {
    
    message("  cutoff: ", observed_arbitrary_cutoff)
    
    arm_agg_bin <- arm_agg %>%
      # group_by(type) %>%
      mutate(
        observed_cutoff_value = quantile(observed_arm, observed_arbitrary_cutoff, na.rm = TRUE),
        observed_bin = observed_arm >= observed_cutoff_value
      ) %>%
      ungroup()
    
    if (length(unique(arm_agg_bin$observed_bin)) < 2) {
      warning(
        "Skipping cutoff ", observed_arbitrary_cutoff,
        " for ", file_name,
        " because only one class is present."
      )
      next
    }
    
    # -----------------------------
    # ROC
    # -----------------------------
    
    roc_obj <- roc(
      response = arm_agg_bin$observed_bin,
      predictor = arm_agg_bin$prediction_arm,
      levels = c(FALSE, TRUE),
      direction = "<",
      quiet = TRUE
    )
    
    roc_auc <- as.numeric(auc(roc_obj))
    
    roc_df_full <- coords(
      roc_obj,
      x = "all",
      ret = c("threshold", "specificity", "sensitivity"),
      transpose = FALSE
    ) %>%
      as.data.frame() %>%
      mutate(
        fpr = 1 - specificity,
        tpr = sensitivity,
        cutoff_quantile = observed_arbitrary_cutoff,
        roc_auc = roc_auc,
        cna_type = cna_type,
        scale = scale_name,
        file = file_name
      ) %>%
      arrange(fpr, tpr) %>%
      distinct(fpr, tpr, .keep_all = TRUE)
    
    # Thin only for plotting
    roc_df <- thin_roc_points(
      roc_df_full
    )
    
    roc_all[[as.character(observed_arbitrary_cutoff)]] <- roc_df
    
    # -----------------------------
    # PR
    # -----------------------------
    
    pr_obj <- pr.curve(
      scores.class0 = arm_agg_bin$prediction_arm[arm_agg_bin$observed_bin == TRUE],
      scores.class1 = arm_agg_bin$prediction_arm[arm_agg_bin$observed_bin == FALSE],
      curve = TRUE
    )
    
    pr_auc <- pr_obj$auc.integral
    baseline_pr <- mean(arm_agg_bin$observed_bin)
    
    pr_df <- data.frame(
      recall = pr_obj$curve[, 1],
      precision = pr_obj$curve[, 2],
      threshold = pr_obj$curve[, 3]
    ) %>%
      arrange(recall) %>%
      group_by(recall) %>%
      summarise(
        precision = max(precision, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(recall) %>%
      mutate(
        precision_interp = rev(cummax(rev(precision))),
        cutoff_quantile = observed_arbitrary_cutoff,
        pr_auc = pr_auc,
        baseline_pr = baseline_pr,
        cna_type = cna_type,
        scale = scale_name,
        file = file_name
      )
    
    pr_all[[as.character(observed_arbitrary_cutoff)]] <- pr_df
    
    # -----------------------------
    # Summary
    # -----------------------------
    
    summary_all[[as.character(observed_arbitrary_cutoff)]] <- data.frame(
      file = file_name,
      cna_type = cna_type,
      scale = scale_name,
      cutoff_quantile = observed_arbitrary_cutoff,
      n_positive = sum(arm_agg_bin$observed_bin == TRUE),
      n_negative = sum(arm_agg_bin$observed_bin == FALSE),
      positive_fraction = mean(arm_agg_bin$observed_bin),
      roc_auc = roc_auc,
      pr_auc = pr_auc,
      baseline_pr = baseline_pr
    )
  }
  
  roc_all_df <- bind_rows(roc_all)
  pr_all_df <- bind_rows(pr_all)
  summary_df <- bind_rows(summary_all)
  
  if (nrow(roc_all_df) == 0 || nrow(pr_all_df) == 0) {
    warning("No ROC/PR results for ", file_name)
    return(NULL)
  }
  
  # -----------------------------
  # Plot all cutoffs together + Jubran et al. only for Arm/Chromosome
  # -----------------------------
  
  roc_all_df <- roc_all_df %>%
    mutate(
      cutoff_label = paste0(
        "q", cutoff_quantile,
        " | AUC=", round(roc_auc, 3)
      )
    )
  
  pr_all_df <- pr_all_df %>%
    mutate(
      cutoff_label = paste0(
        "q", cutoff_quantile,
        " | AUPRC=", round(pr_auc, 3)
      )
    )
  
  overlay_jubran <- scale_name %in% c("Arm-level", "Chromosome-level")
  
  jubran_roc_df <- if (overlay_jubran) {
    get_jubran_xgboost_curve(cna_type, curve_type = "ROC")
  } else {
    NULL
  }
  
  jubran_pr_df <- if (overlay_jubran) {
    get_jubran_xgboost_curve(cna_type, curve_type = "PR")
  } else {
    NULL
  }
  
  # -----------------------------
  # ROC plot
  # -----------------------------
  
  p_ROC_all <- ggplot(
    roc_all_df,
    aes(x = fpr, y = tpr, color = cutoff_label)
  ) +
    geom_step(linewidth = 0.9, direction = "hv", alpha = 0.85) +
    geom_abline(linetype = "dashed", color = "gray60") +
    # coord_equal() +
    theme_classic(base_size = 13) +
    labs(
      title = paste0("ROC curves - ", cna_type, " - ", scale_name),
      x = "False positive rate",
      y = "True positive rate",
      color = NULL
    ) +
    theme(
      legend.position = "right"
    )
  
  if (!is.null(jubran_roc_df) && nrow(jubran_roc_df) > 0) {
    p_ROC_all <- p_ROC_all +
      geom_line(
        data = jubran_roc_df,
        aes(x = fpr, y = tpr, color = method),
        linewidth = 2.5,
        inherit.aes = FALSE
      )
  }
  
  # -----------------------------
  # PR plot
  # -----------------------------
  
  p_PR_all <- ggplot(
    pr_all_df,
    aes(x = recall, y = precision_interp, color = cutoff_label)
  ) +
    geom_step(linewidth = 1, direction = "hv", alpha = 0.85) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_classic(base_size = 13) +
    labs(
      title = paste0("Precision-recall curves - ", cna_type, " - ", scale_name),
      x = "Recall",
      y = "Precision",
      color = NULL
    ) +
    theme(
      legend.position = "right"
    )
  
  if (!is.null(jubran_pr_df) && nrow(jubran_pr_df) > 0) {
    p_PR_all <- p_PR_all +
      geom_line(
        data = jubran_pr_df,
        aes(x = recall, y = precision_interp, color = method),
        linewidth = 2.3,
        inherit.aes = FALSE
      )
  }
  
  combined_plot <- ggarrange(
    p_ROC_all,
    p_PR_all,
    ncol = 2,
    common.legend = FALSE
  )
  
  pdf(
    file.path(out_folder, paste0("ROC_PR_all_cutoffs_", cna_type, "_", safe_scale, ".pdf")),
    width = 15,
    height = 5.5
  )
  print(combined_plot)
  dev.off()
  
  png(
    file.path(out_folder, paste0("ROC_PR_all_cutoffs_", cna_type, "_", safe_scale, ".png")),
    width = 1200, height = 600
  )
  print(combined_plot)
  dev.off()
  
   # -----------------------------
  # Save tables
  # -----------------------------
  
  write.csv(
    summary_df,
    file.path(out_folder, paste0("ROC_PR_all_cutoffs_summary_", cna_type, "_", safe_scale, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    roc_all_df,
    file.path(out_folder, paste0("ROC_all_cutoffs_points_", cna_type, "_", safe_scale, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    pr_all_df,
    file.path(out_folder, paste0("PR_all_cutoffs_points_", cna_type, "_", safe_scale, ".csv")),
    row.names = FALSE
  )
  
  return(summary_df)
}


all_summaries <- map(annotation_files, run_aggregate_roc_pr)

all_summaries_df <- bind_rows(all_summaries)
