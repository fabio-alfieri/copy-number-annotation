library(GenomicRanges)
library(IRanges)
library(data.table)
library(ggplot2)

# -------------------------
# Input files
# -------------------------

manual_file <- '../SPICE_loci_Kaufman/all_460_loci.tsv'
# ampl_file   <- 'annotation-matrix/annotation_ampl_Small-scale.tsv'
# del_file    <- 'annotation-matrix/annotation_del_Small-scale.tsv'
ampl_file   <- 'annotation-matrix/annotation_ampl_Mid-length.tsv'
del_file    <- 'annotation-matrix/annotation_del_Mid-length.tsv'

out_dir <- "../SPICE_loci_Kaufman/manual_loci_vs_selected_midlength_bins"
dir.create(out_dir, showWarnings = FALSE)

# Include both exact and "likely" positive/negative selection labels?
keep_likely <- TRUE

# -------------------------
# Helper functions
# -------------------------

plot_overlap_significance_by_cancer_type <- function(
    stats_dt,
    event_name,
    kaufman_class_name,
    out_file
) {
  
  plot_dt <- copy(stats_dt[
    event == event_name &
      kaufman_class == kaufman_class_name
  ])
  
  plot_dt[, neg_log10_padj := -log10(p_adj_BH)]
  plot_dt[is.infinite(neg_log10_padj), neg_log10_padj := max(
    neg_log10_padj[is.finite(neg_log10_padj)],
    na.rm = TRUE
  ) + 1]
  
  plot_dt[, significant := p_adj_BH < 0.05]
  
  cancer_order <- plot_dt[
    ,
    .(best = max(neg_log10_padj, na.rm = TRUE)),
    by = cancer_type
  ][order(best)]$cancer_type
  
  plot_dt[, cancer_type := factor(cancer_type, levels = cancer_order)]
  
  p <- ggplot(
    plot_dt,
    aes(
      x = cancer_type,
      y = neg_log10_padj,
      fill = selection_class
    )
  ) +
    geom_col(position = "dodge") +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    coord_flip() +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(
        event_name,
        ": enrichment in Kaufman ",
        kaufman_class_name,
        " loci"
      ),
      x = "Cancer type",
      y = expression(-log[10](BH-adjusted~p)),
      fill = "Small-scale annotation",
      caption = "Dashed line indicates BH-adjusted p = 0.05."
    )
  
  ggsave(out_file, p, width = 8, height = 5)
  return(p)
}

add_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

classify_selection <- function(annot_final, keep_likely = TRUE) {
  annot_final <- as.character(annot_final)
  
  selection_class <- rep(NA_character_, length(annot_final))
  
  if (keep_likely) {
    selection_class[grepl("Positive Selection", annot_final, ignore.case = TRUE)] <- "Positive"
    selection_class[grepl("Negative Selection", annot_final, ignore.case = TRUE)] <- "Negative"
  } else {
    selection_class[annot_final == "Positive Selection"] <- "Positive"
    selection_class[annot_final == "Negative Selection"] <- "Negative"
  }
  
  selection_class
}

make_manual_granges <- function(manual_raw) {
  manual <- as.data.frame(manual_raw)
  
  manual$chr <- add_chr(manual$chrom)
  manual$locus_n <- seq_len(nrow(manual))
  
  manual$locus_id <- paste0(
    "locus_", sprintf("%03d", manual$locus_n), "_",
    manual$type, "_",
    manual$chr, "_",
    manual$start_ci, "_",
    manual$end_ci
  )
  
  GRanges(
    seqnames = manual$chr,
    ranges = IRanges(
      start = as.integer(manual$start_ci),
      end   = as.integer(manual$end_ci)
    ),
    locus_id = manual$locus_id,
    locus_n = manual$locus_n,
    manual_type = manual$type,
    pos = manual$pos,
    nr_events_in_data = manual$nr_events_in_data,
    p_value = manual$p_value,
    genes = manual$genes
  )
}

make_selected_small_granges <- function(small_raw, event_name, keep_likely = TRUE) {
  small <- as.data.frame(small_raw)
  
  small$chr <- add_chr(small$chr)
  small$event <- event_name
  
  # In your annotation files, the column "type" is the cancer type
  small$cancer_type <- small$type
  
  small$selection_class <- classify_selection(
    small$annot_final,
    keep_likely = keep_likely
  )
  
  # Keep only positive or negative selected bins
  small <- small[!is.na(small$selection_class), ]
  
  small$bin_uid <- paste(
    small$event,
    small$cancer_type,
    small$binID,
    small$chr,
    small$start,
    small$end,
    sep = "_"
  )
  
  GRanges(
    seqnames = small$chr,
    ranges = IRanges(
      start = as.integer(small$start),
      end   = as.integer(small$end)
    ),
    bin_uid = small$bin_uid,
    binID = small$binID,
    event = small$event,
    cancer_type = small$cancer_type,
    observed = small$observed,
    prediction = small$prediction,
    residual = small$residual,
    pval.residual = small$pval.residual,
    pred_is_correct = small$pred_is_correct,
    region_is_diploid = small$region_is_diploid,
    top1 = small$top1,
    shap_top1 = small$shap_top1,
    is_candidate = small$is_candidate,
    sas_ns = small$sas_ns,
    is_strong_candidate = small$is_strong_candidate,
    annot_final = small$annot_final,
    selection_class = small$selection_class
  )
}

overlap_manual_with_small <- function(manual_gr, small_gr, event_name) {
  
  hits <- findOverlaps(
    query = manual_gr,
    subject = small_gr,
    ignore.strand = TRUE
  )
  
  if (length(hits) == 0) {
    return(data.table())
  }
  
  q <- queryHits(hits)
  s <- subjectHits(hits)
  
  overlap_gr <- pintersect(manual_gr[q], small_gr[s])
  
  res <- data.table(
    event = event_name,
    
    locus_id = mcols(manual_gr)$locus_id[q],
    locus_n = mcols(manual_gr)$locus_n[q],
    manual_type = mcols(manual_gr)$manual_type[q],
    genes = mcols(manual_gr)$genes[q],
    
    chr = as.character(seqnames(manual_gr[q])),
    manual_start = start(manual_gr[q]),
    manual_end = end(manual_gr[q]),
    manual_width_bp = width(manual_gr[q]),
    
    small_bin_uid = mcols(small_gr)$bin_uid[s],
    binID = mcols(small_gr)$binID[s],
    cancer_type = mcols(small_gr)$cancer_type[s],
    
    small_start = start(small_gr[s]),
    small_end = end(small_gr[s]),
    small_width_bp = width(small_gr[s]),
    
    overlap_start = start(overlap_gr),
    overlap_end = end(overlap_gr),
    overlap_width_bp = width(overlap_gr),
    
    frac_manual_locus_covered = width(overlap_gr) / width(manual_gr[q]),
    frac_small_bin_covered = width(overlap_gr) / width(small_gr[s]),
    
    annot_final = mcols(small_gr)$annot_final[s],
    selection_class = mcols(small_gr)$selection_class[s],
    
    observed = mcols(small_gr)$observed[s],
    prediction = mcols(small_gr)$prediction[s],
    residual = mcols(small_gr)$residual[s],
    pval.residual = mcols(small_gr)$pval.residual[s],
    top1 = mcols(small_gr)$top1[s],
    shap_top1 = mcols(small_gr)$shap_top1[s],
    is_candidate = mcols(small_gr)$is_candidate[s],
    is_strong_candidate = mcols(small_gr)$is_strong_candidate[s]
  )
  
  res[]
}

summarize_by_locus <- function(overlaps_dt, event_name) {
  
  if (nrow(overlaps_dt) == 0) {
    return(data.table())
  }
  
  summary_dt <- overlaps_dt[
    ,
    .(
      n_selected_overlap_rows = .N,
      n_unique_selected_bins = uniqueN(small_bin_uid),
      n_cancer_types = uniqueN(cancer_type),
      
      n_positive_rows = sum(selection_class == "Positive"),
      n_negative_rows = sum(selection_class == "Negative"),
      
      n_positive_cancer_types = uniqueN(cancer_type[selection_class == "Positive"]),
      n_negative_cancer_types = uniqueN(cancer_type[selection_class == "Negative"]),
      
      max_frac_manual_locus_covered = max(frac_manual_locus_covered, na.rm = TRUE),
      
      selection_classes_seen = paste(sort(unique(selection_class)), collapse = "; "),
      cancer_types_seen = paste(sort(unique(cancer_type)), collapse = "; "),
      annotations_seen = paste(sort(unique(annot_final)), collapse = "; ")
    ),
    by = .(
      event,
      locus_id,
      locus_n,
      chr,
      manual_start,
      manual_end,
      manual_type,
      genes
    )
  ]
  
  summary_dt[order(locus_n)]
}

# -------------------------
# Load data
# -------------------------

manual_raw <- fread(manual_file)
ampl_raw   <- fread(ampl_file)
del_raw    <- fread(del_file)

manual_gr <- make_manual_granges(manual_raw)

ampl_selected_gr <- make_selected_small_granges(
  ampl_raw,
  event_name = "amplification",
  keep_likely = keep_likely
)

del_selected_gr <- make_selected_small_granges(
  del_raw,
  event_name = "deletion",
  keep_likely = keep_likely
)

# -------------------------
# Overlap analyses
# -------------------------

ampl_overlaps <- overlap_manual_with_small(
  manual_gr = manual_gr,
  small_gr = ampl_selected_gr,
  event_name = "amplification"
)

del_overlaps <- overlap_manual_with_small(
  manual_gr = manual_gr,
  small_gr = del_selected_gr,
  event_name = "deletion"
)

ampl_locus_summary <- summarize_by_locus(ampl_overlaps, "amplification")
del_locus_summary  <- summarize_by_locus(del_overlaps, "deletion")

# Manual loci with no selected overlap
manual_loci_dt <- data.table(
  locus_id = mcols(manual_gr)$locus_id,
  locus_n = mcols(manual_gr)$locus_n,
  chr = as.character(seqnames(manual_gr)),
  manual_start = start(manual_gr),
  manual_end = end(manual_gr),
  manual_type = mcols(manual_gr)$manual_type,
  genes = mcols(manual_gr)$genes
)

ampl_no_overlap <- manual_loci_dt[!locus_id %in% ampl_locus_summary$locus_id]
del_no_overlap  <- manual_loci_dt[!locus_id %in% del_locus_summary$locus_id]


# ============================================================
# Per-cancer-type significance of overlap with Kaufman loci
# Fisher exact test using all small-scale bins as universe
# ============================================================

make_all_small_granges <- function(small_raw, event_name) {
  
  small <- as.data.table(small_raw)
  
  small[, chr := add_chr(chr)]
  small[, event := event_name]
  small[, cancer_type := type]
  
  small[, selection_class := classify_selection(
    annot_final,
    keep_likely = keep_likely
  )]
  
  small[, bin_uid := paste(
    event,
    cancer_type,
    binID,
    chr,
    start,
    end,
    sep = "_"
  )]
  
  GRanges(
    seqnames = small$chr,
    ranges = IRanges(
      start = as.integer(small$start),
      end   = as.integer(small$end)
    ),
    bin_uid = small$bin_uid,
    binID = small$binID,
    event = small$event,
    cancer_type = small$cancer_type,
    annot_final = small$annot_final,
    selection_class = small$selection_class
  )
}

ampl_all_gr <- make_all_small_granges(
  ampl_raw,
  event_name = "amplification"
)

del_all_gr <- make_all_small_granges(
  del_raw,
  event_name = "deletion"
)

# Fisher's test per cancer type

run_fisher_overlap_by_cancer_type <- function(
    all_bins_gr,
    manual_gr,
    event_name,
    kaufman_class = c("All", "OG", "TSG"),
    selection_classes = c("Positive", "Negative")
) {
  
  kaufman_class <- match.arg(kaufman_class)
  
  # Optionally restrict Kaufman/manual loci to OG or TSG
  manual_sub <- manual_gr
  
  if (kaufman_class != "All") {
    manual_sub <- manual_sub[mcols(manual_sub)$manual_type == kaufman_class]
  }
  
  cancer_types <- sort(unique(mcols(all_bins_gr)$cancer_type))
  
  results <- list()
  idx <- 1
  
  for (ct in cancer_types) {
    
    bins_ct <- all_bins_gr[mcols(all_bins_gr)$cancer_type == ct]
    
    # Does each bin overlap Kaufman/manual loci?
    overlaps_kaufman <- countOverlaps(
      bins_ct,
      manual_sub,
      ignore.strand = TRUE
    ) > 0
    
    for (sel in selection_classes) {
      
      # Is each bin selected positive/negative?
      is_selected <- mcols(bins_ct)$selection_class == sel
      is_selected[is.na(is_selected)] <- FALSE
      
      tab <- table(
        selected = factor(is_selected, levels = c(FALSE, TRUE)),
        overlaps_kaufman = factor(overlaps_kaufman, levels = c(FALSE, TRUE))
      )
      
      fisher_res <- fisher.test(tab, alternative = "greater")
      
      n_selected <- sum(is_selected)
      n_selected_overlapping <- sum(is_selected & overlaps_kaufman)
      
      n_background <- sum(!is_selected)
      n_background_overlapping <- sum(!is_selected & overlaps_kaufman)
      
      pct_selected_overlapping <- ifelse(
        n_selected == 0,
        NA_real_,
        100 * n_selected_overlapping / n_selected
      )
      
      pct_background_overlapping <- ifelse(
        n_background == 0,
        NA_real_,
        100 * n_background_overlapping / n_background
      )
      
      results[[idx]] <- data.table(
        event = event_name,
        cancer_type = ct,
        kaufman_class = kaufman_class,
        selection_class = sel,
        
        n_total_bins = length(bins_ct),
        n_selected_bins = n_selected,
        n_selected_bins_overlapping_kaufman = n_selected_overlapping,
        pct_selected_bins_overlapping_kaufman = pct_selected_overlapping,
        
        n_background_bins = n_background,
        n_background_bins_overlapping_kaufman = n_background_overlapping,
        pct_background_bins_overlapping_kaufman = pct_background_overlapping,
        
        odds_ratio = unname(fisher_res$estimate),
        p_value = fisher_res$p.value
      )
      
      idx <- idx + 1
    }
  }
  
  res <- rbindlist(results, fill = TRUE)
  res[, p_adj_BH := p.adjust(p_value, method = "BH")]
  res[]
}

# run the test
stats_ampl_all <- run_fisher_overlap_by_cancer_type(
  all_bins_gr = ampl_all_gr,
  manual_gr = manual_gr,
  event_name = "amplification",
  kaufman_class = "All"
)

stats_ampl_OG <- run_fisher_overlap_by_cancer_type(
  all_bins_gr = ampl_all_gr,
  manual_gr = manual_gr,
  event_name = "amplification",
  kaufman_class = "OG"
)

stats_ampl_TSG <- run_fisher_overlap_by_cancer_type(
  all_bins_gr = ampl_all_gr,
  manual_gr = manual_gr,
  event_name = "amplification",
  kaufman_class = "TSG"
)

stats_del_all <- run_fisher_overlap_by_cancer_type(
  all_bins_gr = del_all_gr,
  manual_gr = manual_gr,
  event_name = "deletion",
  kaufman_class = "All"
)

stats_del_OG <- run_fisher_overlap_by_cancer_type(
  all_bins_gr = del_all_gr,
  manual_gr = manual_gr,
  event_name = "deletion",
  kaufman_class = "OG"
)

stats_del_TSG <- run_fisher_overlap_by_cancer_type(
  all_bins_gr = del_all_gr,
  manual_gr = manual_gr,
  event_name = "deletion",
  kaufman_class = "TSG"
)

overlap_stats_by_cancer_type <- rbindlist(list(
  stats_ampl_all,
  stats_ampl_OG,
  stats_ampl_TSG,
  stats_del_all,
  stats_del_OG,
  stats_del_TSG
), fill = TRUE)

overlap_stats_by_cancer_type[, p_adj_BH_global := p.adjust(p_value, method = "BH")]


direction_matched_stats <- overlap_stats_by_cancer_type[
  (event == "amplification" & kaufman_class == "OG") |
    (event == "deletion" & kaufman_class == "TSG")
]

p_sig_ampl_OG <- plot_overlap_significance_by_cancer_type(
  stats_dt = overlap_stats_by_cancer_type,
  event_name = "amplification",
  kaufman_class_name = "OG",
  out_file = file.path(out_dir, "plot_significance_amplification_Kaufman_OG_by_cancer_type.pdf")
)

p_sig_del_TSG <- plot_overlap_significance_by_cancer_type(
  stats_dt = overlap_stats_by_cancer_type,
  event_name = "deletion",
  kaufman_class_name = "TSG",
  out_file = file.path(out_dir, "plot_significance_deletion_Kaufman_TSG_by_cancer_type.pdf")
)

# ============================================================
# Plots
# ============================================================

plot_locus_counts <- function(overlaps_dt, event_name, out_file) {
  
  plot_dt <- unique(
    overlaps_dt[, .(locus_id, manual_type, selection_class)]
  )
  
  plot_dt <- plot_dt[
    ,
    .N,
    by = .(manual_type, selection_class)
  ]
  
  p <- ggplot(
    plot_dt,
    aes(x = manual_type, y = N, fill = selection_class)
  ) +
    geom_col(position = "dodge") +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(event_name, ": manual loci overlapping selected small-scale bins"),
      x = "Manual locus class",
      y = "Number of manual loci",
      fill = "Small-scale annotation"
    )
  
  ggsave(out_file, p, width = 6, height = 4)
  return(p)
}

plot_by_cancer_type <- function(overlaps_dt, event_name, out_file) {
  
  plot_dt <- unique(
    overlaps_dt[, .(locus_id, manual_type, cancer_type, selection_class)]
  )
  
  plot_dt <- plot_dt[
    ,
    .N,
    by = .(cancer_type, selection_class)
  ]
  
  p <- ggplot(
    plot_dt,
    aes(x = reorder(cancer_type, N), y = N, fill = selection_class)
  ) +
    geom_col(position = "dodge") +
    coord_flip() +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(event_name, ": selected overlaps by cancer type"),
      x = "Cancer type",
      y = "Number of manual loci with overlap",
      fill = "Small-scale annotation"
    )
  
  ggsave(out_file, p, width = 7, height = 5)
  return(p)
}

plot_by_cancer_type_OG_TSG_with_p <- function(
    overlaps_dt,
    stats_dt,
    event_name,
    out_file
) {
  
  plot_dt <- unique(
    overlaps_dt[, .(
      locus_id,
      manual_type,
      cancer_type,
      selection_class
    )]
  )
  
  plot_dt <- plot_dt[
    ,
    .N,
    by = .(cancer_type, manual_type, selection_class)
  ]
  
  event_key <- ifelse(
    event_name == "Amplification",
    "amplification",
    "deletion"
  )
  
  stats_sub <- stats_dt[event == event_key]
  
  stats_sub[, sig_label := fifelse(
    p_adj_BH< 0.001, "***",
    fifelse(
      p_adj_BH < 0.01, "**",
      fifelse(p_adj_BH < 0.05, "*", "")
    )
  )]
  
  stats_sub <- stats_sub[
    ,
    .(
      cancer_type,
      manual_type = kaufman_class,
      selection_class,
      sig_label,
      p_adj_BH,
      odds_ratio
    )
  ]
  
  stats_sub <- stats_sub[manual_type %in% c("OG", "TSG")]
  
  plot_dt <- merge(
    plot_dt,
    stats_sub,
    by = c("cancer_type", "manual_type", "selection_class"),
    all.x = TRUE
  )
  
  cancer_order <- plot_dt[
    ,
    .(total = sum(N)),
    by = cancer_type
  ][order(total)]$cancer_type
  
  plot_dt[, cancer_type := factor(cancer_type, levels = cancer_order)]
  
  p <- ggplot(
    plot_dt,
    aes(x = cancer_type, y = N, fill = manual_type)
  ) +
    geom_col(position = position_dodge(width = 0.8)) +
    geom_text(
      aes(label = sig_label),
      position = position_dodge(width = 0.8),
      hjust = -0.15,
      size = 5
    ) +
    coord_flip() +
    facet_wrap(~ selection_class, nrow = 1) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(event_name, ": selected overlaps by cancer type and manual class"),
      subtitle = "Stars indicate Fisher enrichment over all assayable small-scale bins",
      x = "Cancer type",
      y = "Number of Kaufman loci with selected overlap",
      fill = "Kaufman locus class",
      caption = "* BH-adjusted p < 0.05; ** < 0.01; *** < 0.001"
    )
  
  ggsave(out_file, p, width = 9, height = 5)
  return(p)
}

plot_top_loci <- function(overlaps_dt, event_name, out_file, top_n = 20) {
  
  plot_dt <- unique(
    overlaps_dt[, .(locus_id, manual_type, genes, cancer_type, selection_class)]
  )
  
  plot_dt <- plot_dt[
    ,
    .N,
    by = .(locus_id, manual_type, genes, selection_class)
  ]
  
  plot_dt[
    ,
    locus_label := paste0(
      locus_id,
      " | ",
      substr(gsub("\\[|\\]|'", "", genes), 1, 40)
    )
  ]
  
  top_loci <- plot_dt[
    ,
    .(total = sum(N)),
    by = .(locus_id, locus_label)
  ][order(-total)][1:min(top_n, .N)]
  
  plot_dt <- plot_dt[locus_id %in% top_loci$locus_id]
  plot_dt[, locus_label := factor(locus_label, levels = rev(top_loci$locus_label))]
  
  p <- ggplot(
    plot_dt,
    aes(x = locus_label, y = N, fill = selection_class)
  ) +
    geom_col(position = "stack") +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(event_name, ": top manual loci by selected overlap"),
      x = "Manual locus",
      y = "Number of cancer-type selected overlaps",
      fill = "Small-scale annotation"
    )
  
  ggsave(out_file, p, width = 9, height = 7)
  return(p)
}

plot_locus_heatmap <- function(overlaps_dt, event_name, out_file) {
  
  heat_dt <- unique(
    overlaps_dt[, .(locus_id, locus_n, manual_type, cancer_type, selection_class)]
  )
  
  heat_dt <- heat_dt[
    ,
    .(
      status = paste(sort(unique(selection_class)), collapse = " + ")
    ),
    by = .(locus_id, locus_n, manual_type, cancer_type)
  ]
  
  heat_dt[
    ,
    locus_label := paste0(
      sprintf("%03d", locus_n),
      " | ",
      manual_type,
      " | ",
      locus_id
    )
  ]
  
  heat_dt[, locus_label := factor(locus_label, levels = rev(unique(heat_dt[order(locus_n)]$locus_label)))]
  
  p <- ggplot(
    heat_dt,
    aes(x = cancer_type, y = locus_label, fill = status)
  ) +
    geom_tile(color = "grey85") +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6)
    ) +
    labs(
      title = paste0(event_name, ": locus-level selected overlaps"),
      x = "Cancer type",
      y = "Manual locus",
      fill = "Overlap status"
    )
  
  ggsave(out_file, p, width = 9, height = 12)
  return(p)
}

# -------------------------
# Generate plots
# -------------------------

# p1_ampl <- plot_locus_counts(
#   ampl_overlaps,
#   "Amplification",
#   file.path(out_dir, "plot_amplification_locus_counts.pdf")
# )

# p2_ampl <- plot_by_cancer_type(
#   ampl_overlaps,
#   "Amplification",
#   file.path(out_dir, "plot_amplification_by_cancer_type.pdf")
# )

ggplot(ampl_overlaps) +
  geom_histogram(aes(x = frac_manual_locus_covered))
ggplot(del_overlaps) +
  geom_histogram(aes(x = frac_manual_locus_covered))

p2_ampl_OG_TSG_sig <- plot_by_cancer_type_OG_TSG_with_p(
  overlaps_dt = ampl_overlaps,
  stats_dt = overlap_stats_by_cancer_type,
  event_name = "Amplification",
  out_file = file.path(out_dir, "plot_amplification_by_cancer_type_OG_TSG_with_significance.pdf")
)

p2_del_OG_TSG_sig <- plot_by_cancer_type_OG_TSG_with_p(
  overlaps_dt = del_overlaps,
  stats_dt = overlap_stats_by_cancer_type,
  event_name = "Deletion",
  out_file = file.path(out_dir, "plot_deletion_by_cancer_type_OG_TSG_with_significance.pdf")
)

# p3_ampl <- plot_top_loci(
#   ampl_overlaps,
#   "Amplification",
#   file.path(out_dir, "plot_amplification_top_loci.pdf")
# )
# 
# p4_ampl <- plot_locus_heatmap(
#   ampl_overlaps,
#   "Amplification",
#   file.path(out_dir, "plot_amplification_locus_heatmap.pdf")
# )
# 
# p1_del <- plot_locus_counts(
#   del_overlaps,
#   "Deletion",
#   file.path(out_dir, "plot_deletion_locus_counts.pdf")
# )
# 
# p2_del <- plot_by_cancer_type(
#   del_overlaps,
#   "Deletion",
#   file.path(out_dir, "plot_deletion_by_cancer_type.pdf")
# )
# 
# p3_del <- plot_top_loci(
#   del_overlaps,
#   "Deletion",
#   file.path(out_dir, "plot_deletion_top_loci.pdf")
# )
# 
# p4_del <- plot_locus_heatmap(
#   del_overlaps,
#   "Deletion",
#   file.path(out_dir, "plot_deletion_locus_heatmap.pdf")
# )


# ============================================================
# Require at least X% of each Kaufman locus to be covered
# by selected focal bins
# ============================================================
plot_by_cancer_type_OG_TSG_min_coverage <- function(
    coverage_dt,
    event_name,
    min_kaufman_coverage,
    out_file
) {
  
  plot_dt <- coverage_dt[
    frac_kaufman_locus_covered >= min_kaufman_coverage
  ]
  
  plot_dt <- unique(
    plot_dt[, .(
      locus_id,
      manual_type,
      cancer_type,
      selection_class
    )]
  )
  
  plot_dt <- plot_dt[
    ,
    .N,
    by = .(cancer_type, manual_type, selection_class)
  ]
  
  cancer_order <- plot_dt[
    ,
    .(total = sum(N)),
    by = cancer_type
  ][order(total)]$cancer_type
  
  plot_dt[, cancer_type := factor(cancer_type, levels = cancer_order)]
  
  p <- ggplot(
    plot_dt,
    aes(x = cancer_type, y = N, fill = manual_type)
  ) +
    geom_col(position = "dodge") +
    coord_flip() +
    facet_wrap(~ selection_class, nrow = 1) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(
        event_name,
        ": Kaufman loci with ≥",
        min_kaufman_coverage * 100,
        "% coverage by selected focal bins"
      ),
      x = "Cancer type",
      y = paste0(
        "Number of Kaufman loci ≥",
        min_kaufman_coverage * 100,
        "% covered"
      ),
      fill = "Kaufman locus class"
    )
  
  ggsave(out_file, p, width = 9, height = 5)
  return(p)
}

min_kaufman_coverage <- 0.50
# use 0.60 for 60%
# min_kaufman_coverage <- 0.60

compute_kaufman_locus_coverage <- function(
    manual_gr,
    selected_gr,
    event_name
) {
  
  hits <- findOverlaps(
    query = manual_gr,
    subject = selected_gr,
    ignore.strand = TRUE
  )
  
  if (length(hits) == 0) {
    return(data.table())
  }
  
  q <- queryHits(hits)
  s <- subjectHits(hits)
  
  overlap_gr <- pintersect(manual_gr[q], selected_gr[s])
  
  overlap_dt <- data.table(
    locus_i = q,
    
    event = event_name,
    locus_id = mcols(manual_gr)$locus_id[q],
    locus_n = mcols(manual_gr)$locus_n[q],
    manual_type = mcols(manual_gr)$manual_type[q],
    genes = mcols(manual_gr)$genes[q],
    
    chr = as.character(seqnames(manual_gr[q])),
    manual_start = start(manual_gr[q]),
    manual_end = end(manual_gr[q]),
    manual_width_bp = width(manual_gr[q]),
    
    cancer_type = mcols(selected_gr)$cancer_type[s],
    selection_class = mcols(selected_gr)$selection_class[s],
    annot_final = mcols(selected_gr)$annot_final[s],
    small_bin_uid = mcols(selected_gr)$bin_uid[s],
    
    overlap_start = start(overlap_gr),
    overlap_end = end(overlap_gr),
    overlap_width_bp = width(overlap_gr)
  )
  
  # Collapse all overlapping focal bins within each:
  # Kaufman locus × cancer type × selection class
  coverage_dt <- overlap_dt[
    ,
    {
      reduced_overlap <- reduce(IRanges(
        start = overlap_start,
        end = overlap_end
      ))
      
      .(
        covered_bp = sum(width(reduced_overlap)),
        n_overlap_rows = .N,
        n_unique_focal_bins = uniqueN(small_bin_uid),
        annotations_seen = paste(sort(unique(annot_final)), collapse = "; ")
      )
    },
    by = .(
      event,
      locus_i,
      locus_id,
      locus_n,
      manual_type,
      genes,
      chr,
      manual_start,
      manual_end,
      manual_width_bp,
      cancer_type,
      selection_class
    )
  ]
  
  coverage_dt[
    ,
    pct_kaufman_locus_covered := 100 * covered_bp / manual_width_bp
  ]
  
  coverage_dt[
    ,
    frac_kaufman_locus_covered := covered_bp / manual_width_bp
  ]
  
  coverage_dt[order(event, cancer_type, selection_class, locus_n)]
}

ampl_kaufman_coverage <- compute_kaufman_locus_coverage(
  manual_gr = manual_gr,
  selected_gr = ampl_selected_gr,
  event_name = "amplification"
)

del_kaufman_coverage <- compute_kaufman_locus_coverage(
  manual_gr = manual_gr,
  selected_gr = del_selected_gr,
  event_name = "deletion"
)

ampl_kaufman_coverage_filtered <- ampl_kaufman_coverage[
  frac_kaufman_locus_covered >= min_kaufman_coverage
]

del_kaufman_coverage_filtered <- del_kaufman_coverage[
  frac_kaufman_locus_covered >= min_kaufman_coverage
]

p_ampl_50cov <- plot_by_cancer_type_OG_TSG_min_coverage(
  coverage_dt = ampl_kaufman_coverage,
  event_name = "Amplification",
  min_kaufman_coverage = min_kaufman_coverage,
  out_file = file.path(out_dir, paste0(
    "plot_amplification_Kaufman_loci_",
    min_kaufman_coverage * 100,
    "pct_coverage_by_cancer_type_OG_TSG.pdf"
  ))
)

p_del_50cov <- plot_by_cancer_type_OG_TSG_min_coverage(
  coverage_dt = del_kaufman_coverage,
  event_name = "Deletion",
  min_kaufman_coverage = min_kaufman_coverage,
  out_file = file.path(out_dir, paste0(
    "plot_deletion_Kaufman_loci_",
    min_kaufman_coverage * 100,
    "pct_coverage_by_cancer_type_OG_TSG.pdf"
  ))
)

