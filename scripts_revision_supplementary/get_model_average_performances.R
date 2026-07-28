rm(list=ls())
gc(full=T)

packages <- c("stringr", "tidyr", "ggplot2", "tidyverse")

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

model.outputs <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/canonical", full.names = TRUE, pattern = "Small-scale")

patterns <- c("Arm-level",
              "Chromosome-level",
              "Mid-length",
              "Small-scale",
              "NoCluster", 
              "ampl", "del",
              "LOCO", 
              as.character(0:2000),
              "LOCTO",
              c("BRCA", "COADREAD", "ESCA", 
                "GBMLGG", "KIRC", "KIRP", 
                "LUAD", "LUSC", "OV", 
                "PAAD", "STAD"),
              "NSS", 
              "OO",
              "OS",
              "PRO", 
              "new", 
              "XL", 
              "PROXY", 
              "LOFO",
              "LAFO",
              "Chromosome", 
              "Centromere", 
              "Centromere", 
              "distance.to.telomere", 
              "distance.to.centromere", 
              "Length",
              "NO",
              "OG",
              "Type",
              "CS",
              "CSP",
              "scaled",
              "full",
              "mean",
              "TEST",
              "TTS",
              "covThr",
              "zero",
              "0.5",
              "0.8",
              "ON_PCAWG",
              "Hic",
              "repliseq",
              "1e05",
              "1e06",
              "4e06",
              "5e06")

avg_performances <- lapply(X = model.outputs, 
                           FUN = function(model_path){
                             
                             load(model_path)
                             
                             parts <- strsplit(x = model_path, split = "_")[[1]]
                             model_name <- paste(parts[parts %in% patterns], collapse = "_")
                             
                             # Extract whatever precedes "Output.regressor" in the basename as the
                             # coverage-cutoff label (e.g. "1e05", "1e06"). If nothing precedes it,
                             # this is the standard model (3Mbp cutoff).
                             base_name <- basename(model_path)
                             cutoff_prefix <- str_extract(base_name, "^.*(?=_?Output\\.regressor)")
                             cutoff_label <- if (is.na(cutoff_prefix) || cutoff_prefix == "") {
                               "3e06"
                             } else {
                               str_remove(cutoff_prefix, "_$")
                             }
                             
                             if(("LOCO" %in% parts) & ("scaled" %in% parts) & ("full" %in% parts)){
                               chr <- as.numeric(parts[9])
                               model <- Output.regressor[[1]][[chr]]
                             } else if (("LOCO" %in% parts) & ("scaled" %in% parts) & ("mean" %in% parts)){
                               chr <- as.numeric(parts[9])
                               model <- Output.regressor[[1]][[chr]]
                             } else if (("LOCO" %in% parts) & ("scaled" %in% parts) & ("TEST" %in% parts)){
                               chr <- as.numeric(parts[9])
                               model <- Output.regressor[[1]][[chr]]  
                             } else if (("LOCO" %in% parts) & ("scaled" %in% parts)){
                               chr <- as.numeric(parts[7])
                               model <- Output.regressor[[1]][[chr]]
                             } else if (("LOCO" %in% parts) & ("OV" %in% parts)){
                               chr <- as.numeric(parts[6])
                               cancer_type <- parts[9]
                               model <- Output.regressor[[1]][[chr]]
                             } else if (("LOCO" %in% parts) & ("KIRC" %in% parts)){
                               chr <- as.numeric(parts[6])
                               cancer_type <- parts[9]
                               model <- Output.regressor[[1]][[chr]]
                             } else if(("LOCO" %in% parts) & ("NSS" %in% parts)) {
                               chr <- as.numeric(parts[6])
                               model <- Output.regressor[[1]][[chr]]
                             } else if (("LOCO" %in% parts) & !("NSS" %in% parts)){
                               chr <- as.numeric(parts[5])
                               model <- Output.regressor[[1]][[chr]]
                             } else if("LOCTO" %in% parts){
                               ct <- parts[5]
                               model <- Output.regressor[[1]][[ct]]
                             } else if("TTS" %in% parts){
                               label <- paste0(parts[5], "_", parts[6])
                               model <- Output.regressor[[1]][[label]]
                             } else if (("CSP" %in% parts) & ("scaled" %in% parts)){
                               label <- paste0(parts[7], "_", parts[8])
                               model <- Output.regressor[[1]][[label]]
                             } else if ("CSP" %in% parts){
                               label <- paste0(parts[5], "_", parts[6])
                               model <- Output.regressor[[1]][[label]]
                             } else if (("Hic" %in% parts) & ("repliseq" %in% parts)){
                               # cancer type sits at the start of the *basename*, e.g.
                               # ".../HiC_repliseq_models//BRCA_Output.regressor_w_Hic_and_repliseq_..."
                               # parts[1] is unreliable here because the full path itself contains
                               # underscores (e.g. "HiC_repliseq_models"), so extract from basename instead.
                               base_parts <- strsplit(base_name, split = "_")[[1]]
                               cancer_type <- base_parts[base_parts %in% patterns][1]
                               
                               model <- Output.regressor[[1]]
                               model_name <- paste0(cancer_type, "_HiC_repliseq_", model_name)
                             } else{
                               model <- Output.regressor[[1]]
                             }
                             
                             rm(Output.regressor)
                             gc()
                             
                             r2_array <- unlist(lapply(model, function(x){x$Performance$R2}))
                             rmse_array <- unlist(lapply(model, function(x){x$Performance$RMSE}))
                             pearson_array <- unlist(lapply(model, function(x){x$Performance$Pearson}))
                             spearman_array <- unlist(lapply(model, function(x){x$Performance$Spearman}))
                             
                             # Mean and SD of y_train per fold, then averaged across folds
                             ytrain_mean_array <- unlist(lapply(model, function(x){ mean(x$y_train, na.rm = TRUE) }))
                             ytrain_sd_array   <- unlist(lapply(model, function(x){ sd(x$y_train,   na.rm = TRUE) }))
                             
                             # NEW: number of training rows per fold (from X_train), averaged across folds
                             ntrain_array <- unlist(lapply(model, function(x){ nrow(x$X_train) }))
                             
                             # NEW: number (and proportion) of zeros in y_train per fold, then averaged across folds
                             nzero_train_array <- unlist(lapply(model, function(x){ sum(x$y_train == 0, na.rm = TRUE) }))
                             prop_zero_train_array <- unlist(lapply(model, function(x){
                               mean(x$y_train == 0, na.rm = TRUE)
                             }))
                             
                             out.list <- list(avg_r2 = mean(r2_array),
                                              avg_rmse = mean(rmse_array),
                                              avg_pearson = mean(pearson_array),
                                              avg_spearman = mean(spearman_array),
                                              sd_r2 = sd(r2_array),
                                              sd_rmse = sd(rmse_array),
                                              sd_pearson = sd(pearson_array),
                                              sd_spearman = sd(spearman_array),
                                              mean_ytrain = mean(ytrain_mean_array),
                                              sd_ytrain   = mean(ytrain_sd_array),
                                              n_train     = mean(ntrain_array),   # NEW: avg training-set size across folds
                                              n_zero_train    = mean(nzero_train_array),      # NEW: avg number of zeros in y_train across folds
                                              prop_zero_train = mean(prop_zero_train_array),  # NEW: avg proportion of zeros in y_train across folds
                                              cutoff = cutoff_label
                             )
                             
                             out.list <- list(out.list)
                             
                             names(out.list) <- model_name
                             print(out.list[[1]]$avg_r2)
                             
                             print(paste0("processed ", model_name))
                             
                             return(out.list)
                             
                           })

avg_performance_table <- do.call(what = rbind, 
                                 lapply(X = avg_performances, 
                                        FUN = function(sublist){
                                          
                                          row.name <- names(sublist)
                                          vals <- unlist(sublist[[1]])
                                          row.df <- t(setNames(data.frame(vals), row.name))
                                          
                                        }))

avg_performance_table <- as.data.frame(avg_performance_table)

# Ensure cutoff stays a plain character column (rbind of mixed types can coerce it oddly)
avg_performance_table$cutoff <- as.character(avg_performance_table$cutoff)

avg_performance_table <- avg_performance_table[!(rownames(avg_performance_table) %in% to.discard),]

avg_performance_path <- "/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/avg_performances_complete_CSP_scaled.tsv"

write.table(x = avg_performance_table, 
            file = avg_performance_path,
            quote = F, sep = "\t", row.names = T, col.names = T)


### heatmap comparing cancer types

library(ggplot2)
library(dplyr)
library(stringr)

cancer_types <- c("BRCA", "COADREAD", "ESCA", "GBMLGG", "KIRC", 
                  "KIRP", "LUAD", "LUSC", "OV", "PAAD", "STAD")

plot_df <- as.data.frame(avg_performance_table) %>%
  tibble::rownames_to_column("id") %>%
  mutate(
    train = str_extract(id, paste(cancer_types, collapse = "|")),
    test  = str_extract(id, paste0("(?<=_)(", paste(cancer_types, collapse = "|"), ")(?=_Length)"))
  )

# Plot
ggplot(plot_df, aes(x = test, y = train, fill = avg_r2)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(
    label = sprintf("%.2f", avg_r2),
    color = avg_r2 < 0.5
  ), size = 2.8, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#f7fbff", "#6baed6", "#08519c"),
    values = c(0, 0.3, 1),
    limits = c(0, 1),
    name   = "R²"
  ) +
  scale_color_manual(values = c("TRUE" = "grey20", "FALSE" = "white"), guide = "none") +
  scale_x_discrete(limits = cancer_types) +
  scale_y_discrete(limits = cancer_types) +
  labs(
    x     = "Test cancer type",
    y     = "Train cancer type",
    title = "Cross-cancer prediction performance (R²)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    panel.grid        = element_blank(),
    plot.title        = element_text(face = "bold"),
    legend.key.height = unit(1.5, "cm")
  ) +
  coord_fixed()

### barplot comparing canonical, OO, OS, faceted by ampl/del
library(tidyverse)
library(patchwork)

# ── Parse avg_performance_table directly ──────────────────────────────────────
df <- avg_performance_table %>%
  as.data.frame() %>%
  rownames_to_column("model") %>%
  mutate(across(-c(model, cutoff), as.numeric)) %>%
  mutate(
    modality = case_when(
      str_starts(model, "OS_") ~ "OS",
      str_starts(model, "OO_") ~ "OO",
      TRUE                     ~ "Full"
    ),
    cnv_type = case_when(
      str_detect(model, "_ampl_covThr_zero$") ~ "Amplification",
      str_detect(model, "_del_covThr_zero$")  ~ "Deletion",
      TRUE ~ NA_character_
    ),
    class = model %>%
      str_remove("^(OS_|OO_)") %>%
      str_remove("^Length_") %>%
      str_remove("_(ampl|del)_covThr_zero$"),
    class    = factor(class, levels = c("Arm-level", "Chromosome-level",
                                        "Mid-length", "Small-scale", "NoCluster")),
    modality = factor(modality, levels = c("OS", "OO", "Full")),
    cnv_type = factor(cnv_type, levels = c("Amplification", "Deletion"))
  )

pal <- c("OS" = "#E07B54", "OO" = "#5B8DB8", "Full" = "#5BAD72")

# ── Plot: Grouped bar chart with SD + value labels, faceted by ampl/del ─────
p_bar <- ggplot(df, aes(x = class, y = avg_r2, fill = modality)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65,
           color = "white", linewidth = 0.3) +
  geom_errorbar(aes(ymin = avg_r2 - sd_r2, ymax = avg_r2 + sd_r2),
                position = position_dodge(width = 0.75),
                width = 0.25, linewidth = 0.5, color = "grey30") +
  geom_text(aes(label = sprintf("%.2f", avg_r2),
                y = avg_r2 + sd_r2),
            position = position_dodge(width = 0.75),
            vjust = -0.5, size = 2.7, color = "grey30") +
  scale_fill_manual(values = pal, name = "Modality") +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0),
                     breaks = seq(0, 1, 0.2)) +
  facet_wrap(~ cnv_type, nrow = 2) +
  labs(x = NULL, y = expression(Average~R^2)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(angle = 25, hjust = 1),
    legend.position    = "right",
    strip.text         = element_text(face = "bold", size = 12)
  )

combined_plot <- p_bar

# ── Save ────────────────────────────────────────────────────────────────────
ggsave(
  filename = "/home/ieo7429/Scrivania/results_regressor_gab/canonical_OO_OS_comparison_allclasses.pdf",
  plot     = combined_plot,
  width    = 10,
  height   = 12
)

### LOCO scaled vs unscaled vs scaled_full

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)

# ── 1. Parse row names ────────────────────────────────────────────────────────
avg_performance_table <- as.data.frame(avg_performance_table)
avg_performance_table$model <- rownames(avg_performance_table)

parsed <- avg_performance_table %>%
  mutate(
    condition = case_when(
      str_detect(model, "LOCO_scaled_full")                  ~ "Scaled Full",
      str_detect(model, "LOCO_scaled_mean")                  ~ "Scaled Mean",
      str_detect(model, "LOCO_scaled_TEST")                  ~ "Scaled Test",
      str_detect(model, "LOCO_scaled")                       ~ "Scaled",
      str_detect(model, "LOCO_\\d+_.+_[A-Z]+_ampl$")        ~ "LOCO_TTS",  # cancer type = uppercase
      TRUE                                                   ~ "Unscaled"
    ),
    chr = case_when(
      str_detect(model, "LOCO_scaled_full")                  ~ str_extract(model, "(?<=LOCO_scaled_full_)\\d+"),
      str_detect(model, "LOCO_scaled_mean")                  ~ str_extract(model, "(?<=LOCO_scaled_mean_)\\d+"),
      str_detect(model, "LOCO_scaled_TEST")                  ~ str_extract(model, "(?<=LOCO_scaled_TEST_)\\d+"),
      str_detect(model, "LOCO_scaled")                       ~ str_extract(model, "(?<=LOCO_scaled_)\\d+"),
      str_detect(model, "LOCO_\\d+_.+_[A-Z]+_ampl$")        ~ str_extract(model, "(?<=LOCO_)\\d+"),
      TRUE                                                   ~ str_extract(model, "(?<=LOCO_)\\d+")
    ),
    chr = as.integer(chr)
  )

# Keep only chromosomes present in ALL THREE conditions
shared_chrs <- parsed %>%
  group_by(chr) %>%
  filter(n_distinct(condition) == 3) %>%
  pull(chr) %>%
  unique()

plot_df <- parsed %>%
  filter(chr %in% shared_chrs) %>%
  mutate(
    condition = factor(condition, levels = c("Unscaled", "Scaled", "Scaled Full", "Scaled Mean", "Scaled Test", "LOCO_TTS")),
    chr       = factor(chr)
  )

# ── 2. Generic paired plot function (now 3 conditions) ───────────────────────
make_paired_plot <- function(df, metric, y_label, title_label) {
  ggplot(df, aes(x = condition, y = .data[[metric]], group = chr)) +
    geom_line(aes(colour = chr), alpha = 0.55, linewidth = 0.6) +
    geom_point(aes(colour = chr, fill = chr),
               size = 3, shape = 21, stroke = 0.4, alpha = 0.85) +
    geom_boxplot(aes(group = condition),
                 width = 0.25, outlier.shape = NA,
                 fill = NA, colour = "grey30",
                 linewidth = 0.7, alpha = 0.8) +
    scale_colour_viridis_d(option = "turbo", name = "Chr") +
    scale_fill_viridis_d(option = "turbo",   name = "Chr") +
    theme_bw(base_size = 11) +
    labs(title = title_label, x = NULL, y = y_label) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.key.size    = unit(0.4, "cm"),
      legend.text        = element_text(size = 8),
      plot.title         = element_text(face = "bold", size = 11)
    )
}

p1 <- make_paired_plot(plot_df, "avg_r2", "R²", "R²")


#### compute correlation between sd of target variable and performance

metrics       <- c("avg_r2", "avg_rmse", "avg_pearson", "avg_spearman")
delta_metrics <- c("delta_r2", "delta_rmse", "delta_pearson", "delta_spearman")
test_desc     <- c("mean_test_mean", "mean_test_sd", "mean_test_iqr", "pooled_test_sd")

# Isolate unscaled as the reference
unscaled_ref <- plot_df %>%
  filter(condition == "Unscaled") %>%
  dplyr::select(chr, all_of(metrics)) %>%
  rename_with(~ paste0(.x, "_unscaled"), all_of(metrics))

# All non-unscaled conditions
delta_df <- plot_df %>%
  filter(condition == "Scaled") %>%
  inner_join(unscaled_ref, by = "chr") %>%
  mutate(
    delta_r2       =  avg_r2       - avg_r2_unscaled,
    delta_rmse     = -(avg_rmse    - avg_rmse_unscaled),   # sign-flip: lower RMSE = better
    delta_pearson  =  avg_pearson  - avg_pearson_unscaled,
    delta_spearman =  avg_spearman - avg_spearman_unscaled
  ) %>%
  dplyr::select(chr, condition, starts_with("delta_"))

output.files <- list.files(path = "/mnt/ghost/fabiogokce/gokce/Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/MLTables_updated/InteractomeINSIDER",
                           full.names = TRUE, 
                           recursive = FALSE, 
                           pattern = "ampl_or_del_0.1Mbp_covThr_zero")

output.files <- output.files[1] # remove location

print(output.files)

name <- "ampl_or_del_0.1Mbp_covThr_zero"
level <- "0.1Mbp"
bincov <- "covThr_zero"
folder <- "InteractomeINSIDER"

chromosome <- read.table(file = '/mnt/ghost/fabiogokce/gokce/Codes/Codes-CNAs/MethodII/fabio/Data/chromosome.tsv', header = T)
centromeres <- read.table(file = '/mnt/ghost/fabiogokce/gokce/Codes/Codes-CNAs/MethodII/fabio/Data/centomere.tsv', header = T)

tmp <- full_join(chromosome, centromeres, by = 'chr') %>% dplyr::select(chr, Centromere_Length, Chromosome_Length, Centromere_Type)
tmp$chr <- as.character(tmp$chr)

load(output.files)
clus.group <- gsub("_ampl_or_del_0.1Mbp_covThr_zero.RData","",basename(output.files))

clusters <- names(ML.Tables[["ampl_or_del"]])[1] # select only Mid-Length and small models

ml_table_backup <- ML.Tables[["ampl_or_del"]][[clusters]]
ml_table_backup$bin <- as.character(ml_table_backup$bin)
print(clusters)

ml_table_backup$chr <- do.call(rbind, str_split(ml_table_backup$bin,'_'))[,1]
ml_table_backup <- left_join(ml_table_backup, tmp, by = 'chr')

per_chr <- ml_table_backup %>% group_by(chr) %>% summarise(sd_target = sd(ampl_score), mean_target = mean(ampl_score), coef = (mean_target / sd_target))

plot_df_final <- merge(delta_df, per_chr, by = "chr")

library(ggplot2)
library(ggrepel)

# ── Pearson + Spearman + p-values for the annotation ─────────────────────────
r_pearson  <- cor.test(plot_df_final$delta_r2, plot_df_final$coef,
                       method = "pearson")
r_spearman <- cor.test(plot_df_final$delta_r2, plot_df_final$coef,
                       method = "spearman", exact = FALSE)

anno_text <- sprintf(
  "Pearson r = %.2f  (p = %.3f)\nSpearman ρ = %.2f  (p = %.3f)",
  r_pearson$estimate,  r_pearson$p.value,
  r_spearman$estimate, r_spearman$p.value
)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(plot_df_final,
            aes(x = mean_target, y = delta_r2, label = chr)) +
  
  # zero-delta reference
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  
  # OLS ribbon + line
  geom_smooth(method = "lm", se = TRUE,
              colour = "#3A6EA5", fill = "#3A6EA5",
              alpha = 0.12, linewidth = 0.9) +
  
  # points coloured by direction of delta
  geom_point(aes(fill = delta_r2 > 0),
             shape = 21, size = 4, stroke = 0.5,
             colour = "white", alpha = 0.9) +
  scale_fill_manual(values = c("TRUE"  = "#4A9D6F",   # scaling helped
                               "FALSE" = "#C0503A"),  # scaling hurt
                    labels = c("TRUE"  = "Improved",
                               "FALSE" = "Worsened"),
                    name   = "Scaling effect") +
  
  # chromosome labels
  geom_text_repel(size = 3, colour = "grey25",
                  max.overlaps = 20, seed = 42,
                  box.padding = 0.4) +
  
  # correlation annotation
  annotate("text",
           x    = max(plot_df_final$sd_target) * 0.98,
           y    = max(plot_df_final$delta_r2)  * 0.98,
           label = anno_text,
           hjust = 1, vjust = 1, size = 3.2,
           colour = "grey20",
           family = "mono") +
  
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  
  labs(
    title    = "Δ R² (scaled − unscaled) vs. test-set mean per chromosome",
    subtitle = "Green = scaling improved R²  ·  Red = scaling worsened R²",
    x        = "Mean of y_test (held-out chromosome)",
    y        = expression(Delta~R^2~(scaled - unscaled))
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(colour = "grey50", size = 10),
    legend.position  = "bottom"
  )

# ggsave("delta_r2_vs_sd_target.pdf", width = 7, height = 6)


###

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)

# ── 1. Classify rows ──────────────────────────────────────────────────────────

avg_performance_table <- avg_performance_table %>%
  tibble::rownames_to_column("model_id") %>%
  mutate(
    model_type = case_when(
      model_id == "Length_Mid-length_ampl_covThr_zero"                ~ "Canonical Full",
      model_id == "TTS_OV_OV_Length_Mid-length_ampl_covThr_zero"      ~ "OV Full",
      model_id == "TTS_KIRC_KIRC_Length_Mid-length_ampl_covThr_zero"  ~ "KIRC Full",
      str_detect(model_id, "^LOCO_scaled_")                          ~ "LOCO Full Scaled",
      str_detect(model_id, "_KIRC_ampl_covThr_zero$")                ~ "LOCO KIRC",
      str_detect(model_id, "_OV_ampl_covThr_zero$")                  ~ "LOCO OV",
      str_detect(model_id, "^LOCO_")                                 ~ "LOCO Full Unscaled",
      TRUE                                                            ~ "Other"
    )
  )

# ── 2. Aggregate LOCO groups (mean ± sd across chromosomes) ───────────────────

loco_agg <- avg_performance_table %>%
  filter(str_detect(model_type, "^LOCO")) %>%
  group_by(model_type) %>%
  summarise(
    across(c(avg_r2, avg_rmse, avg_pearson, avg_spearman),
           list(mean = ~ mean(.x, na.rm = TRUE),
                sd   = ~ sd(.x,   na.rm = TRUE)),
           .names = "{.col}__{.fn}"),
    n = n(),
    .groups = "drop"
  )

print(loco_agg)

# ── 3. Build a unified summary table for plotting ─────────────────────────────

single_models <- avg_performance_table %>%
  filter(model_type %in% c("Canonical Full", "OV Full", "KIRC Full")) %>%
  transmute(
    model_type,
    avg_r2__mean      = avg_r2,       avg_r2__sd      = 0,
    avg_rmse__mean    = avg_rmse,     avg_rmse__sd    = 0,
    avg_pearson__mean = avg_pearson,  avg_pearson__sd = 0,
    avg_spearman__mean= avg_spearman, avg_spearman__sd= 0,
    n = 1
  )

plot_data <- bind_rows(single_models, loco_agg) %>%
  mutate(
    model_type = factor(model_type, levels = c(
      "Canonical Full", "LOCO Full Unscaled", "LOCO Full Scaled",
      "OV Full",        "LOCO OV",
      "KIRC Full",      "LOCO KIRC"
    ))
  )

# ── 4. One grouped bar chart per metric ───────────────────────────────────────

metrics <- c("avg_r2", "avg_pearson", "avg_spearman", "avg_rmse")
metric_labels <- c(
  avg_r2       = "Average R²",
  avg_pearson  = "Average Pearson r",
  avg_spearman = "Average Spearman ρ",
  avg_rmse     = "Average RMSE"
)

palette <- c(
  "Canonical Full"     = "#2166AC",
  "LOCO Full Unscaled" = "#74ADD1",
  "LOCO Full Scaled"   = "#ABD9E9",
  "OV Full"            = "#D73027",
  "LOCO OV"            = "#F4A582",
  "KIRC Full"          = "#1A9641",
  "LOCO KIRC"          = "#78C679"
)

plots <- lapply(metrics, function(m) {
  df <- plot_data %>%
    transmute(
      model_type,
      mean = get(paste0(m, "__mean")),
      sd   = get(paste0(m, "__sd"))
    )
  
  ggplot(df, aes(x = model_type, y = mean, fill = model_type)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
    geom_errorbar(
      aes(ymin = mean - sd, ymax = mean + sd),
      width = 0.25, linewidth = 0.6, colour = "grey30"
    ) +
    scale_fill_manual(values = palette, guide = "none") +
    labs(
      title   = metric_labels[m],
      x       = NULL,
      y       = metric_labels[m],
      caption = "Error bars = ±1 SD across chromosomes (LOCO models only)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 35, hjust = 1),
      plot.title   = element_text(face = "bold"),
      plot.caption = element_text(colour = "grey50", size = 8)
    )
})

# ── 5. Arrange and display ────────────────────────────────────────────────────

combined <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]]) +
  plot_annotation(
    title    = "Model Performance: Canonical vs OV vs KIRC — Full vs LOCO",
    subtitle = "LOCO bars show mean across chromosomes 1–22",
    theme    = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave("/home/ieo7429/Scrivania/loco_performance_comparison.pdf", combined, width = 16, height = 10)

combined

####

### CS_TTS analysis

avg_performance_table <- as.data.frame(avg_performance_table)

library(dplyr)
library(stringr)
library(tibble)

df <- avg_performance_table %>%
  rownames_to_column("model") %>%
  mutate(
    cancer = case_when(
      str_detect(model, "^CS_ESCA_") ~ "ESCA",
      str_detect(model, "^CS_LUSC_") ~ "LUSC",
      str_detect(model, "^CS_OV_")   ~ "OV",
      TRUE                           ~ "PanCancer"
    ),
    
    chr = str_extract(model, "\\d+(?=_Length)")
  ) %>%
  mutate(
    chr = factor(chr, levels = as.character(1:22))
  )


ggplot(
  df %>% filter(!is.na(avg_r2)),
  aes(x = chr, y = avg_r2, fill = cancer)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_errorbar(
    aes(
      ymin = avg_r2 - sd_r2,
      ymax = avg_r2 + sd_r2
    ),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  theme_bw() +
  labs(
    x = "Chromosome",
    y = expression(R^2),
    fill = "Cancer type"
  ) +
  theme(
    axis.text.x = element_text(angle = 90)
  )


#### inspect chr3 performances

df <- avg_performance_table %>%
  rownames_to_column(var = "row_name") %>%
  mutate(cancer = sub("^CS_(.+)_3_Length_Mid-length_ampl_covThr_zero$", "\\1", row_name))


ggplot(df, aes(x = reorder(cancer, avg_r2), y = avg_r2)) +
  geom_col(fill = "#4C8DC4", width = 0.7, na.rm = TRUE) +
  geom_errorbar(
    aes(ymin = avg_r2 - sd_r2, ymax = avg_r2 + sd_r2),
    width = 0.3, linewidth = 0.7, color = "grey30", na.rm = TRUE
  ) +
  geom_text(
    aes(label = ifelse(is.na(avg_r2), "NA", round(avg_r2, 3))),
    hjust = -0.15, size = 3.2, color = "grey20"
  ) +
  coord_flip(ylim = c(0, 1.08)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Average R² by Cancer Type",
    subtitle = "Mid-length amplification signal | error bars = ±1 SD",
    x        = NULL,
    y        = expression(Average~R^2)
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey50", size = 10),
    axis.text.y   = element_text(size = 11),
    panel.grid.major.x = element_line(color = "grey90")
  )


##### COMPARING small scale at different coverage thresholds

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

df <- avg_performance_table %>%
  rownames_to_column("model") %>%
  mutate(
    variant = case_when(
      grepl("_ampl_", model) ~ "ampl",
      grepl("_del_", model)  ~ "del",
      TRUE ~ NA_character_
    ),
    covThr = sub(".*_covThr_", "", model)
  ) %>%
  # order covThr nicely
  mutate(covThr = factor(covThr, levels = c("zero", "0.5", "0.8")))

# reshape means
means_long <- df %>%
  dplyr::select(model, variant, covThr, avg_r2, avg_rmse, avg_pearson, avg_spearman) %>%
  pivot_longer(cols = starts_with("avg_"), names_to = "metric", values_to = "mean") %>%
  mutate(metric = sub("avg_", "", metric))

# reshape sds
sds_long <- df %>%
  dplyr::select(model, variant, covThr, sd_r2, sd_rmse, sd_pearson, sd_spearman) %>%
  pivot_longer(cols = starts_with("sd_"), names_to = "metric", values_to = "sd") %>%
  mutate(metric = sub("sd_", "", metric))

plot_df <- left_join(means_long, sds_long, by = c("model", "variant", "covThr", "metric"))

ggplot(plot_df, aes(x = covThr, y = mean, color = covThr)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.15) +
  facet_grid(metric ~ variant, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    x = "Coverage threshold",
    y = "Value (mean ± SD)",
    title = "Performance metrics by coverage threshold, split by ampl/del"
  )

#### chr-specific barplot

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# ── 1. Extract chromosome number and order factor levels numerically ──────────

plot_df <- avg_performance_table %>%
  tibble::rownames_to_column("model_id") %>%
  mutate(
    chr = str_extract(model_id, "^CS_[0-9]+") %>% str_remove("^CS_") %>% as.integer()
  ) %>%
  mutate(
    chr = factor(chr, levels = sort(unique(chr)))
  )

# ── 2. Barplot with error bars (mean ± sd) ────────────────────────────────────

p <- ggplot(plot_df, aes(x = chr, y = avg_r2, fill = avg_r2)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = avg_r2 - sd_r2, ymax = avg_r2 + sd_r2),
    width = 0.25, linewidth = 0.5, colour = "grey30"
  ) +
  scale_fill_gradientn(
    colors = c("#f7fbff", "#6baed6", "#08519c"),
    values = c(0, 0.5, 1),
    limits = c(0, 1),
    name   = "R²"
  ) +
  labs(
    x     = "Chromosome",
    y     = expression(paste("Average ", R^2)),
    title = "Chromosome-specific model performance (R²)",
    caption = "Error bars = ±1 SD across folds"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 0, hjust = 0.5),
    plot.title   = element_text(face = "bold"),
    plot.caption = element_text(colour = "grey50", size = 8),
    legend.position = "none"
  )

p

# ── 3. Save ─────────────────────────────────────────────────────────────────

ggsave(
  filename = "/home/ieo7429/Scrivania/chr_specific_r2_barplot.pdf",
  plot     = p,
  width    = 12,
  height   = 6
)

##### NEW: avg_r2 (and sd_r2) vs. training-set size (n_train)

library(dplyr)
library(ggplot2)
library(tibble)

df_size <- avg_performance_table %>%
  as.data.frame() %>%
  rownames_to_column("model") %>%
  mutate(
    avg_r2  = as.numeric(avg_r2),
    sd_r2   = as.numeric(sd_r2),
    n_train = as.numeric(n_train)
  ) %>%
  filter(!is.na(avg_r2), !is.na(n_train)) %>%
  # order bars by training-set size, left to right
  mutate(model = factor(model, levels = model[order(n_train)]))

p_size <- ggplot(df_size, aes(x = model)) +
  # R² bars + error bars
  geom_col(
    aes(y = avg_r2, fill = "clustered"),
    width = 0.6
  ) +
  geom_errorbar(
    aes(ymin = pmax(avg_r2 - sd_r2, 0), ymax = avg_r2 + sd_r2),
    width = 0.2, colour = "grey30", linewidth = 0.4
  ) +
  # training-set size shown as a number above each bar, not as a bar
  geom_text(
    aes(y = avg_r2 + sd_r2, label = scales::comma(n_train)),
    vjust = -0.5, size = 2.6, colour = "grey20", angle = 0, hjust = 0
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("clustered" = "steelblue")
  ) +
  scale_y_continuous(
    name   = expression(Average~R^2~(mean%+-%SD)),
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "Model performance (R²) vs. training-set size",
    x     = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank(),
    plot.title      = element_text(face = "bold"),
    legend.position = "none"
  )

p_size

ggsave(
  filename = "/home/ieo7429/Scrivania/results_regressor_gab/avg_r2_vs_ntrain.pdf",
  plot = p_size, width = 14, height = 7
)

##### COMPARING Small-scale performance across coverage cutoffs (3Mbp vs 1e05 vs 1e06)

library(dplyr)
library(tibble)
library(ggplot2)

df <- avg_performance_table %>%
  rownames_to_column("model") %>%
  mutate(
    across(c(avg_r2, sd_r2, mean_ytrain, sd_ytrain),
           ~as.numeric(na_if(as.character(.), "<NA>"))),
    avg_r2 = tidyr::replace_na(avg_r2, 0),
    sd_r2 = tidyr::replace_na(sd_r2, 0),
    mean_ytrain = tidyr::replace_na(mean_ytrain, 0),
    sd_ytrain = tidyr::replace_na(sd_ytrain, 0),
    type = ifelse(grepl("ampl", model), "Amplification", "Deletion"),
    cutoff = factor(cutoff, levels = c("1e05", "1e06", "3e06", "4e06", "5e06"))
  )

p <- ggplot(df, aes(cutoff, avg_r2)) +
  geom_col(
    aes(fill = "Model performance"),
    width = 0.7
  ) +
  geom_errorbar(
    aes(ymin = pmax(avg_r2 - sd_r2, 0),
        ymax = avg_r2 + sd_r2),
    width = 0.15
  ) +
  facet_wrap(~type) +
  scale_fill_manual(
    name = NULL,
    values = c("Model performance" = "steelblue")
  ) +
  scale_y_continuous(name = expression("Average " * R^2)) +
  labs(x = "Cutoff") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top"
  )

ggsave(
  filename = "/home/ieo7429/Scrivania/results_regressor_gab/avg_performance_plot.pdf",
  plot = p,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)

##### NEW: avg_r2 vs. proportion of zeros in y_train, and vs. sd of y_train, across cutoffs, split by ampl/del
##### NAs in avg_r2 are set to 0 rather than dropped, so every point appears

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

df_zero <- avg_performance_table %>%
  as.data.frame() %>%
  rownames_to_column("model") %>%
  mutate(
    across(c(avg_r2, sd_r2, n_zero_train, prop_zero_train, sd_ytrain),
           ~ as.numeric(na_if(as.character(.), "<NA>"))),
    avg_r2 = tidyr::replace_na(avg_r2, 0),
    sd_r2  = tidyr::replace_na(sd_r2, 0),
    type   = ifelse(grepl("ampl", model), "Amplification", "Deletion"),
    cutoff = factor(cutoff, levels = c("1e05", "1e06", "3e06", "4e06", "5e06"))
  )

## Panel A: avg_r2 vs. proportion of zeros in y_train, coloured by cutoff, faceted by type
p_zero_scatter <- ggplot(df_zero, aes(x = prop_zero_train, y = avg_r2, colour = cutoff)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_errorbar(aes(ymin = pmax(avg_r2 - sd_r2, 0), ymax = avg_r2 + sd_r2),
                width = 0.01, alpha = 0.6) +
  facet_wrap(~type) +
  scale_colour_viridis_d(option = "turbo", name = "Cutoff") +
  labs(
    x     = "Proportion of zeros in y_train",
    y     = expression(Average~R^2),
    title = "Average R² vs. proportion of zeros in training target",
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey50", size = 10)
  )

p_zero_scatter

ggsave(
  filename = "/home/ieo7429/Scrivania/results_regressor_gab/avg_r2_vs_propzero_train.pdf",
  plot = p_zero_scatter, width = 9, height = 5
)

## Panel B: avg_r2 vs. sd of y_train, coloured by cutoff, faceted by type, with correlation annotation
cor_labels <- df_zero %>%
  group_by(type) %>%
  summarise(
    r = cor(sd_ytrain, avg_r2, use = "complete.obs"),
    p = cor.test(sd_ytrain, avg_r2)$p.value,
    x = max(sd_ytrain, na.rm = TRUE) * 0.7,
    y = max(avg_r2, na.rm = TRUE) * 0.95,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r = %.2f, p = %.3f", r, p))

p_sd_scatter <- ggplot(df_zero, aes(x = sd_ytrain, y = avg_r2)) +
  geom_point(aes(colour = cutoff), size = 3, alpha = 0.9) +
  geom_errorbar(aes(ymin = pmax(avg_r2 - sd_r2, 0), ymax = avg_r2 + sd_r2, colour = cutoff),
                width = 0.01, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.5, linetype = "dashed") +
  geom_text(data = cor_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 4, fontface = "italic") +
  facet_wrap(~type, scales = "free_x") +
  scale_colour_viridis_d(option = "turbo", name = "Cutoff") +
  labs(
    x     = "SD of y_train",
    y     = expression(Average~R^2),
    title = "Average R² vs. SD of training target"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey50", size = 10)
  )

p_sd_scatter

ggsave(
  filename = "/home/ieo7429/Scrivania/results_regressor_gab/avg_r2_vs_sdytrain.pdf",
  plot = p_sd_scatter, width = 9, height = 5
)
