# rm(list=ls())
gc(full=T)

library(data.table)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(readxl)


# ============================================================
# Input files
# ============================================================

biscut_file <- "../BISCUT_Shih/41586_2023_6266_MOESM6_ESM.xlsx"

# ampl_file <- "annotation-matrix/annotation_ampl_Small-scale.tsv"
# del_file  <- "annotation-matrix/annotation_del_Small-scale.tsv"
ampl_file <- "annotation-matrix/annotation_ampl_Mid-length.tsv"
del_file  <- "annotation-matrix/annotation_del_Mid-length.tsv"

out_dir <- "../BISCUT_Shih/biscut_vs_selected_midlength_bins"
dir.create(out_dir, showWarnings = FALSE)

# Include "Likely Positive/Negative Selection" in our annotations?
keep_likely <- TRUE

# Minimum fraction of BISCUT peak covered by one of our focal bins
coverage_cutoff <- 0.01
# coverage_cutoff <- 0.60

# ============================================================
# Helper functions
# ============================================================

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

clean_gene_label <- function(x, max_chars = 70) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "No genes"
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars), "..."), x)
}


# ============================================================
# 1. Load BISCUT and collapse to peak-level loci
# ============================================================

make_biscut_granges <- function(biscut_raw) {
  
  dt <- as.data.table(biscut_raw)
  
  required_cols <- c(
    "Chr",
    "Gene",
    "type",
    "direction",
    "negpos",
    "Peak.Start",
    "Peak.End"
  )
  
  missing_cols <- setdiff(required_cols, names(dt))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing expected BISCUT columns:\n",
      paste(missing_cols, collapse = ", "),
      "\n\nAvailable columns are:\n",
      paste(names(dt), collapse = "\n")
    )
  }
  
  dt[, chr := add_chr(Chr)]
  
  dt[, event := fifelse(
    direction == "amp", "amplification",
    fifelse(direction == "del", "deletion", NA_character_)
  )]
  
  dt[, reference_selection_class := fifelse(
    negpos == "p", "Positive",
    fifelse(negpos == "n", "Negative", NA_character_)
  )]
  
  dt[, reference_cancer_type := as.character(type)]
  
  # Harmonize BISCUT labels to match your annotation matrices.
  # Edit this mapping if needed after checking setdiff() below.
  dt[reference_cancer_type %in% c("GBM", "LGG"),
     reference_cancer_type := "GBMLGG"]
  
  dt[reference_cancer_type %in% c("COAD", "READ"),
     reference_cancer_type := "COADREAD"]
  
  dt[, peak_start := as.integer(Peak.Start)]
  dt[, peak_end := as.integer(Peak.End)]
  
  dt <- dt[
    !is.na(chr) &
      !is.na(peak_start) &
      !is.na(peak_end) &
      !is.na(event) &
      !is.na(reference_selection_class) &
      !is.na(reference_cancer_type)
  ]
  
  # Collapse repeated gene rows into unique BISCUT peak intervals.
  biscut_loci_dt <- dt[
    ,
    .(
      genes = paste(sort(unique(Gene)), collapse = ";"),
      n_genes_observed = uniqueN(Gene),
      num_genes = if ("num_genes" %in% names(dt)) {
        suppressWarnings(max(as.numeric(num_genes), na.rm = TRUE))
      } else {
        NA_real_
      },
      combined_sig = if ("combined_sig" %in% names(dt)) {
        suppressWarnings(max(as.numeric(combined_sig), na.rm = TRUE))
      } else {
        NA_real_
      },
      ks_p = if ("ks_p" %in% names(dt)) {
        suppressWarnings(min(as.numeric(ks_p), na.rm = TRUE))
      } else {
        NA_real_
      },
      log10_ks_p = if ("log10_ks_p" %in% names(dt)) {
        suppressWarnings(max(as.numeric(log10_ks_p), na.rm = TRUE))
      } else {
        NA_real_
      },
      cytoband = if ("Cytoband" %in% names(dt)) {
        paste(sort(unique(Cytoband)), collapse = ";")
      } else {
        NA_character_
      }
    ),
    by = .(
      reference_cancer_type,
      event,
      reference_selection_class,
      chr,
      peak_start,
      peak_end
    )
  ]
  
  biscut_loci_dt[, reference_id := paste0(
    "BISCUT_",
    reference_cancer_type, "_",
    event, "_",
    reference_selection_class, "_",
    chr, "_", peak_start, "_", peak_end
  )]
  
  gr <- GRanges(
    seqnames = biscut_loci_dt$chr,
    ranges = IRanges(
      start = biscut_loci_dt$peak_start,
      end = biscut_loci_dt$peak_end
    ),
    reference_id = biscut_loci_dt$reference_id,
    reference_cancer_type = biscut_loci_dt$reference_cancer_type,
    event = biscut_loci_dt$event,
    reference_selection_class = biscut_loci_dt$reference_selection_class,
    genes = biscut_loci_dt$genes,
    n_genes_observed = biscut_loci_dt$n_genes_observed,
    num_genes = biscut_loci_dt$num_genes,
    combined_sig = biscut_loci_dt$combined_sig,
    ks_p = biscut_loci_dt$ks_p,
    log10_ks_p = biscut_loci_dt$log10_ks_p,
    cytoband = biscut_loci_dt$cytoband
  )
  
  gr
}


# ============================================================
# 2. Load our amplification/deletion matrices
# ============================================================

make_all_small_granges <- function(small_raw, event_name, keep_likely = TRUE) {
  
  small <- as.data.table(small_raw)
  
  small[, chr := add_chr(chr)]
  small[, event := event_name]
  
  # In your annotation matrices, "type" is the cancer type.
  small[, cancer_type := as.character(type)]
  
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
      end = as.integer(small$end)
    ),
    bin_uid = small$bin_uid,
    binID = small$binID,
    event = small$event,
    cancer_type = small$cancer_type,
    annot_final = small$annot_final,
    selection_class = small$selection_class,
    observed = small$observed,
    prediction = small$prediction,
    residual = small$residual,
    pval.residual = small$pval.residual
  )
}

make_selected_small_granges <- function(small_raw, event_name, keep_likely = TRUE) {
  
  gr <- make_all_small_granges(
    small_raw = small_raw,
    event_name = event_name,
    keep_likely = keep_likely
  )
  
  gr[!is.na(mcols(gr)$selection_class)]
}


# ============================================================
# 3. Compute BISCUT overlaps with our selected bins
# ============================================================

overlap_biscut_with_selected <- function(
    biscut_gr,
    selected_gr,
    event_name,
    require_same_cancer_type = TRUE
) {
  
  ref <- biscut_gr[mcols(biscut_gr)$event == event_name]
  
  hits <- findOverlaps(
    query = ref,
    subject = selected_gr,
    ignore.strand = TRUE
  )
  
  if (length(hits) == 0) {
    return(data.table())
  }
  
  q <- queryHits(hits)
  s <- subjectHits(hits)
  
  overlap_gr <- pintersect(ref[q], selected_gr[s])
  
  out <- data.table(
    event = event_name,
    
    reference_id = mcols(ref)$reference_id[q],
    reference_cancer_type = mcols(ref)$reference_cancer_type[q],
    reference_selection_class = mcols(ref)$reference_selection_class[q],
    genes = mcols(ref)$genes[q],
    cytoband = mcols(ref)$cytoband[q],
    combined_sig = mcols(ref)$combined_sig[q],
    ks_p = mcols(ref)$ks_p[q],
    log10_ks_p = mcols(ref)$log10_ks_p[q],
    
    cancer_type = mcols(selected_gr)$cancer_type[s],
    selection_class = mcols(selected_gr)$selection_class[s],
    selected_bin_uid = mcols(selected_gr)$bin_uid[s],
    binID = mcols(selected_gr)$binID[s],
    
    chr = as.character(seqnames(ref[q])),
    reference_start = start(ref[q]),
    reference_end = end(ref[q]),
    reference_width_bp = width(ref[q]),
    
    selected_start = start(selected_gr[s]),
    selected_end = end(selected_gr[s]),
    selected_width_bp = width(selected_gr[s]),
    
    overlap_start = start(overlap_gr),
    overlap_end = end(overlap_gr),
    overlap_width_bp = width(overlap_gr),
    
    frac_reference_locus_covered = width(overlap_gr) / width(ref[q]),
    frac_selected_bin_covered = width(overlap_gr) / width(selected_gr[s])
  )
  
  if (require_same_cancer_type) {
    out <- out[reference_cancer_type == cancer_type]
  }
  
  out[]
}


# ============================================================
# 4. Fisher enrichment test per cancer type
# ============================================================

bin_overlaps_reference_with_min_coverage <- function(
    bins_gr,
    ref_gr,
    coverage_cutoff = 0.50
) {
  
  out <- rep(FALSE, length(bins_gr))
  
  if (length(bins_gr) == 0 || length(ref_gr) == 0) {
    return(out)
  }
  
  hits <- findOverlaps(
    query = bins_gr,
    subject = ref_gr,
    ignore.strand = TRUE
  )
  
  if (length(hits) == 0) {
    return(out)
  }
  
  q <- queryHits(hits)
  s <- subjectHits(hits)
  
  overlap_gr <- pintersect(bins_gr[q], ref_gr[s])
  
  hit_dt <- data.table(
    bin_i = q,
    ref_i = s,
    overlap_width_bp = width(overlap_gr),
    ref_width_bp = width(ref_gr[s])
  )
  
  hit_dt[, frac_reference_locus_covered := overlap_width_bp / ref_width_bp]
  
  hit_dt <- hit_dt[
    ,
    .(
      max_frac_reference_locus_covered = max(frac_reference_locus_covered)
    ),
    by = bin_i
  ]
  
  out[hit_dt$bin_i] <- hit_dt$max_frac_reference_locus_covered >= coverage_cutoff
  
  out
}

run_biscut_fisher_by_cancer_type <- function(
    all_bins_gr,
    biscut_gr,
    event_name,
    coverage_cutoff = 0.50,
    reference_selection_classes = c("Positive", "Negative"),
    selection_classes = c("Positive", "Negative")
) {
  
  cancer_types <- sort(unique(mcols(all_bins_gr)$cancer_type))
  
  results <- list()
  idx <- 1
  
  for (ct in cancer_types) {
    
    bins_ct <- all_bins_gr[mcols(all_bins_gr)$cancer_type == ct]
    
    for (ref_sel in reference_selection_classes) {
      
      ref_ct <- biscut_gr[
        mcols(biscut_gr)$event == event_name &
          mcols(biscut_gr)$reference_cancer_type == ct &
          mcols(biscut_gr)$reference_selection_class == ref_sel
      ]
      
      for (sel in selection_classes) {
        
        is_selected <- mcols(bins_ct)$selection_class == sel
        is_selected[is.na(is_selected)] <- FALSE
        
        overlaps_reference <- bin_overlaps_reference_with_min_coverage(
          bins_gr = bins_ct,
          ref_gr = ref_ct,
          coverage_cutoff = coverage_cutoff
        )
        
        tab <- table(
          selected = factor(is_selected, levels = c(FALSE, TRUE)),
          overlaps_reference = factor(overlaps_reference, levels = c(FALSE, TRUE))
        )
        
        fisher_res <- fisher.test(tab, alternative = "greater")
        
        n_selected <- sum(is_selected)
        n_selected_overlapping <- sum(is_selected & overlaps_reference)
        
        n_background <- sum(!is_selected)
        n_background_overlapping <- sum(!is_selected & overlaps_reference)
        
        results[[idx]] <- data.table(
          event = event_name,
          cancer_type = ct,
          reference_selection_class = ref_sel,
          selection_class = sel,
          coverage_cutoff = coverage_cutoff,
          
          n_reference_loci = length(ref_ct),
          n_total_bins = length(bins_ct),
          
          n_selected_bins = n_selected,
          n_selected_bins_overlapping_reference = n_selected_overlapping,
          pct_selected_bins_overlapping_reference = ifelse(
            n_selected == 0,
            NA_real_,
            100 * n_selected_overlapping / n_selected
          ),
          
          n_background_bins = n_background,
          n_background_bins_overlapping_reference = n_background_overlapping,
          pct_background_bins_overlapping_reference = ifelse(
            n_background == 0,
            NA_real_,
            100 * n_background_overlapping / n_background
          ),
          
          odds_ratio = unname(fisher_res$estimate),
          p_value = fisher_res$p.value
        )
        
        idx <- idx + 1
      }
    }
  }
  
  res <- rbindlist(results, fill = TRUE)
  res[, p_adj_BH := p.adjust(p_value, method = "BH")]
  res[]
}


# ============================================================
# 5. Summary table: number of BISCUT loci recovered
# ============================================================

summarize_biscut_overlap_by_cancer_type <- function(
    overlaps_dt,
    biscut_gr,
    event_name,
    coverage_cutoff = 0.50
) {
  
  ref_dt <- data.table(
    event = mcols(biscut_gr)$event,
    cancer_type = mcols(biscut_gr)$reference_cancer_type,
    reference_selection_class = mcols(biscut_gr)$reference_selection_class,
    reference_id = mcols(biscut_gr)$reference_id
  )
  
  ref_dt <- ref_dt[event == event_name]
  
  denom <- ref_dt[
    ,
    .(
      n_reference_loci_total = uniqueN(reference_id)
    ),
    by = .(event, cancer_type, reference_selection_class)
  ]
  
  ov <- copy(overlaps_dt)
  ov <- ov[
    frac_reference_locus_covered >= coverage_cutoff
  ]
  
  numer <- ov[
    ,
    .(
      n_reference_loci_recovered = uniqueN(reference_id)
    ),
    by = .(
      event,
      cancer_type,
      reference_selection_class,
      selection_class
    )
  ]
  
  out <- merge(
    numer,
    denom,
    by = c("event", "cancer_type", "reference_selection_class"),
    all.x = TRUE
  )
  
  out[, pct_reference_loci_recovered := 100 * n_reference_loci_recovered / n_reference_loci_total]
  
  out[order(event, cancer_type, selection_class, reference_selection_class)]
}


# ============================================================
# 6. Plot by cancer type with stars
# ============================================================

plot_biscut_by_cancer_type_with_p <- function(
    overlaps_dt,
    stats_dt,
    event_name,
    coverage_cutoff = 0.50,
    out_file
) {
  
  event_key <- ifelse(
    event_name == "Amplification",
    "amplification",
    "deletion"
  )
  
  plot_dt <- copy(overlaps_dt)
  
  plot_dt <- plot_dt[
    frac_reference_locus_covered >= coverage_cutoff
  ]
  
  plot_dt <- unique(
    plot_dt[, .(
      reference_id,
      reference_selection_class,
      cancer_type,
      selection_class
    )]
  )
  
  plot_dt <- plot_dt[
    ,
    .N,
    by = .(
      cancer_type,
      reference_selection_class,
      selection_class
    )
  ]
  
  if (nrow(plot_dt) == 0) {
    warning("No overlaps to plot for ", event_name)
    return(NULL)
  }
  
  stats_sub <- copy(stats_dt[
    event == event_key &
      abs(coverage_cutoff - get("coverage_cutoff", parent.frame())) < 1e-9
  ])
  
  stats_sub[, sig_label := fifelse(
    p_adj_BH_global < 0.001, "***",
    fifelse(
      p_adj_BH_global < 0.01, "**",
      fifelse(p_adj_BH_global < 0.05, "*", "")
    )
  )]
  
  stats_sub <- stats_sub[
    ,
    .(
      cancer_type,
      reference_selection_class,
      selection_class,
      sig_label,
      p_adj_BH_global,
      odds_ratio
    )
  ]
  
  plot_dt <- merge(
    plot_dt,
    stats_sub,
    by = c(
      "cancer_type",
      "reference_selection_class",
      "selection_class"
    ),
    all.x = TRUE
  )
  
  plot_dt[is.na(sig_label), sig_label := ""]
  
  cancer_order <- plot_dt[
    ,
    .(total = sum(N)),
    by = cancer_type
  ][order(total)]$cancer_type
  
  plot_dt[, cancer_type := factor(cancer_type, levels = cancer_order)]
  
  y_max <- max(plot_dt$N, na.rm = TRUE)
  
  p <- ggplot(
    plot_dt,
    aes(
      x = cancer_type,
      y = N,
      fill = reference_selection_class
    )
  ) +
    geom_col(position = position_dodge(width = 0.8)) +
    geom_text(
      aes(label = sig_label),
      position = position_dodge(width = 0.8),
      hjust = -0.15,
      size = 5
    ) +
    coord_flip(clip = "off") +
    facet_wrap(~ selection_class, nrow = 1) +
    theme_classic(base_size = 14) +
    theme(
      plot.margin = margin(5.5, 40, 5.5, 5.5)
    ) +
    labs(
      title = paste0(
        event_name,
        ": selected focal bins overlapping BISCUT loci"
      ),
      subtitle = paste0(
        "Only overlaps covering ≥",
        coverage_cutoff * 100,
        "% of the BISCUT peak are shown"
      ),
      x = "Cancer type",
      y = paste0(
        "Number of BISCUT loci with ≥",
        coverage_cutoff * 100,
        "% coverage"
      ),
      fill = "BISCUT class",
      caption = "* BH-adjusted p < 0.05; ** < 0.01; *** < 0.001"
    ) +
    ylim(0, y_max * 1.15)
  
  ggsave(out_file, p, width = 9, height = 5)
  
  return(p)
}


# ============================================================
# 7. Plot top recovered BISCUT loci with genes
# ============================================================

plot_top_biscut_loci <- function(
    overlaps_dt,
    event_name,
    selection_class_name = "Positive",
    reference_selection_class_name = "Positive",
    coverage_cutoff = 0.50,
    top_n = 25,
    out_file
) {
  
  event_key <- ifelse(
    event_name == "Amplification",
    "amplification",
    "deletion"
  )
  
  plot_dt <- copy(overlaps_dt)
  
  plot_dt <- plot_dt[
    event == event_key &
      selection_class == selection_class_name &
      reference_selection_class == reference_selection_class_name &
      frac_reference_locus_covered >= coverage_cutoff
  ]
  
  if (nrow(plot_dt) == 0) {
    warning("No loci to plot for ", event_name)
    return(NULL)
  }
  
  plot_dt <- plot_dt[
    ,
    .(
      max_frac_reference_locus_covered = max(frac_reference_locus_covered),
      n_cancer_types = uniqueN(cancer_type),
      cancer_types_seen = paste(sort(unique(cancer_type)), collapse = "; "),
      genes = unique(genes)[1],
      cytoband = unique(cytoband)[1],
      combined_sig = max(combined_sig, na.rm = TRUE),
      log10_ks_p = max(log10_ks_p, na.rm = TRUE),
      chr = unique(chr)[1],
      reference_start = unique(reference_start)[1],
      reference_end = unique(reference_end)[1]
    ),
    by = .(reference_id)
  ]
  
  plot_dt[, genes_clean := clean_gene_label(genes)]
  plot_dt[, locus_label := paste0(
    genes_clean,
    " | ",
    cytoband,
    " | ",
    chr,
    ":",
    format(reference_start, big.mark = ",", scientific = FALSE),
    "-",
    format(reference_end, big.mark = ",", scientific = FALSE)
  )]
  
  plot_dt <- plot_dt[order(-combined_sig, -log10_ks_p)]
  plot_dt <- plot_dt[1:min(top_n, .N)]
  
  plot_dt[, locus_label := factor(locus_label, levels = rev(locus_label))]
  
  p <- ggplot(
    plot_dt,
    aes(
      x = locus_label,
      y = max_frac_reference_locus_covered * 100
    )
  ) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 11) +
    labs(
      title = paste0(
        event_name,
        ": top recovered BISCUT loci"
      ),
      subtitle = paste0(
        "Our ",
        selection_class_name,
        " bins vs BISCUT ",
        reference_selection_class_name,
        " loci"
      ),
      x = "",
      y = "% of BISCUT peak covered"
    )
  
  ggsave(out_file, p, width = 11, height = 8)
  
  return(p)
}


# ============================================================
# 8. Run analysis
# ============================================================

message("Loading files...")

biscut_raw <- read_xlsx(biscut_file, sheet = 1, skip = 16)
ampl_raw <- fread(ampl_file)
del_raw <- fread(del_file)

message("Building GRanges...")

biscut_gr <- make_biscut_granges(biscut_raw)

ampl_all_gr <- make_all_small_granges(
  small_raw = ampl_raw,
  event_name = "amplification",
  keep_likely = keep_likely
)

del_all_gr <- make_all_small_granges(
  small_raw = del_raw,
  event_name = "deletion",
  keep_likely = keep_likely
)

ampl_selected_gr <- make_selected_small_granges(
  small_raw = ampl_raw,
  event_name = "amplification",
  keep_likely = keep_likely
)

del_selected_gr <- make_selected_small_granges(
  small_raw = del_raw,
  event_name = "deletion",
  keep_likely = keep_likely
)

message("BISCUT peak summary:")
print(table(mcols(biscut_gr)$event, mcols(biscut_gr)$reference_selection_class))

message("Cancer-type overlap between BISCUT and our matrices:")
print(intersect(
  sort(unique(mcols(biscut_gr)$reference_cancer_type)),
  sort(unique(ampl_raw$type))
))

message("BISCUT cancer types not in our matrices:")
print(setdiff(
  sort(unique(mcols(biscut_gr)$reference_cancer_type)),
  sort(unique(ampl_raw$type))
))

message("Our cancer types not in BISCUT:")
print(setdiff(
  sort(unique(ampl_raw$type)),
  sort(unique(mcols(biscut_gr)$reference_cancer_type))
))


# -------------------------
# Overlap tables
# -------------------------

message("Computing overlaps...")

biscut_ampl_overlaps <- overlap_biscut_with_selected(
  biscut_gr = biscut_gr,
  selected_gr = ampl_selected_gr,
  event_name = "amplification",
  require_same_cancer_type = TRUE
)

biscut_del_overlaps <- overlap_biscut_with_selected(
  biscut_gr = biscut_gr,
  selected_gr = del_selected_gr,
  event_name = "deletion",
  require_same_cancer_type = TRUE
)

fwrite(
  biscut_ampl_overlaps,
  file.path(out_dir, "BISCUT_amplification_overlaps_all.tsv"),
  sep = "\t"
)

fwrite(
  biscut_del_overlaps,
  file.path(out_dir, "BISCUT_deletion_overlaps_all.tsv"),
  sep = "\t"
)

fwrite(
  biscut_ampl_overlaps[frac_reference_locus_covered >= coverage_cutoff],
  file.path(out_dir, paste0(
    "BISCUT_amplification_overlaps_min_",
    coverage_cutoff * 100,
    "pct_peak_coverage.tsv"
  )),
  sep = "\t"
)

fwrite(
  biscut_del_overlaps[frac_reference_locus_covered >= coverage_cutoff],
  file.path(out_dir, paste0(
    "BISCUT_deletion_overlaps_min_",
    coverage_cutoff * 100,
    "pct_peak_coverage.tsv"
  )),
  sep = "\t"
)


# -------------------------
# Fisher statistics
# -------------------------

message("Running Fisher enrichment tests...")

biscut_stats_ampl <- run_biscut_fisher_by_cancer_type(
  all_bins_gr = ampl_all_gr,
  biscut_gr = biscut_gr,
  event_name = "amplification",
  coverage_cutoff = coverage_cutoff
)

biscut_stats_del <- run_biscut_fisher_by_cancer_type(
  all_bins_gr = del_all_gr,
  biscut_gr = biscut_gr,
  event_name = "deletion",
  coverage_cutoff = coverage_cutoff
)

biscut_overlap_stats_by_cancer_type <- rbindlist(list(
  biscut_stats_ampl,
  biscut_stats_del
), fill = TRUE)

biscut_overlap_stats_by_cancer_type[
  ,
  p_adj_BH_global := p.adjust(p_value, method = "BH")
]

fwrite(
  biscut_overlap_stats_by_cancer_type,
  file.path(out_dir, paste0(
    "BISCUT_overlap_stats_by_cancer_type_min_",
    coverage_cutoff * 100,
    "pct_peak_coverage.tsv"
  )),
  sep = "\t"
)

# Direction-matched stats: our Positive with BISCUT Positive,
# our Negative with BISCUT Negative.
biscut_direction_matched_stats <- biscut_overlap_stats_by_cancer_type[
  selection_class == reference_selection_class
]

fwrite(
  biscut_direction_matched_stats,
  file.path(out_dir, paste0(
    "BISCUT_direction_matched_stats_min_",
    coverage_cutoff * 100,
    "pct_peak_coverage.tsv"
  )),
  sep = "\t"
)


# -------------------------
# Recovery summary
# -------------------------

biscut_ampl_summary <- summarize_biscut_overlap_by_cancer_type(
  overlaps_dt = biscut_ampl_overlaps,
  biscut_gr = biscut_gr,
  event_name = "amplification",
  coverage_cutoff = coverage_cutoff
)

biscut_del_summary <- summarize_biscut_overlap_by_cancer_type(
  overlaps_dt = biscut_del_overlaps,
  biscut_gr = biscut_gr,
  event_name = "deletion",
  coverage_cutoff = coverage_cutoff
)

biscut_overlap_summary <- rbindlist(list(
  biscut_ampl_summary,
  biscut_del_summary
), fill = TRUE)

fwrite(
  biscut_overlap_summary,
  file.path(out_dir, paste0(
    "BISCUT_reference_loci_recovered_by_cancer_type_min_",
    coverage_cutoff * 100,
    "pct_peak_coverage.tsv"
  )),
  sep = "\t"
)


# -------------------------
# Plots
# -------------------------

message("Plotting...")

(p_biscut_ampl <- plot_biscut_by_cancer_type_with_p(
  overlaps_dt = biscut_ampl_overlaps,
  stats_dt = biscut_overlap_stats_by_cancer_type,
  event_name = "Amplification",
  coverage_cutoff = coverage_cutoff,
  out_file = file.path(out_dir, paste0(
    "plot_BISCUT_amplification_by_cancer_type_min_",
    coverage_cutoff * 100,
    "pct_peak_coverage_with_significance.pdf"
  ))
))

(p_biscut_del <- plot_biscut_by_cancer_type_with_p(
  overlaps_dt = biscut_del_overlaps,
  stats_dt = biscut_overlap_stats_by_cancer_type,
  event_name = "Deletion",
  coverage_cutoff = coverage_cutoff,
  out_file = file.path(out_dir, paste0(
    "plot_BISCUT_deletion_by_cancer_type_min_",
    coverage_cutoff * 100,
    "pct_peak_coverage_with_significance.pdf"
  ))
))

# Top direction-matched examples
(p_top_ampl_pos <- plot_top_biscut_loci(
  overlaps_dt = biscut_ampl_overlaps,
  event_name = "Amplification",
  selection_class_name = "Positive",
  reference_selection_class_name = "Positive",
  coverage_cutoff = coverage_cutoff,
  top_n = 25,
  out_file = file.path(out_dir, paste0(
    "plot_top_BISCUT_positive_amplification_loci_recovered_min_",
    coverage_cutoff * 100,
    "pct.pdf"
  ))
))

(p_top_del_pos <- plot_top_biscut_loci(
  overlaps_dt = biscut_del_overlaps,
  event_name = "Deletion",
  selection_class_name = "Positive",
  reference_selection_class_name = "Positive",
  coverage_cutoff = coverage_cutoff,
  top_n = 25,
  out_file = file.path(out_dir, paste0(
    "plot_top_BISCUT_positive_deletion_loci_recovered_min_",
    coverage_cutoff * 100,
    "pct.pdf"
  ))
))

(p_top_ampl_neg <- plot_top_biscut_loci(
  overlaps_dt = biscut_ampl_overlaps,
  event_name = "Amplification",
  selection_class_name = "Negative",
  reference_selection_class_name = "Negative",
  coverage_cutoff = coverage_cutoff,
  top_n = 25,
  out_file = file.path(out_dir, paste0(
    "plot_top_BISCUT_negative_amplification_loci_recovered_min_",
    coverage_cutoff * 100,
    "pct.pdf"
  ))
))

(p_top_del_neg <- plot_top_biscut_loci(
  overlaps_dt = biscut_del_overlaps,
  event_name = "Deletion",
  selection_class_name = "Negative",
  reference_selection_class_name = "Negative",
  coverage_cutoff = coverage_cutoff,
  top_n = 25,
  out_file = file.path(out_dir, paste0(
    "plot_top_BISCUT_negative_deletion_loci_recovered_min_",
    coverage_cutoff * 100,
    "pct.pdf"
  ))
))

