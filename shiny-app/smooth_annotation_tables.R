## ------------------------------------------------------------------------
## smooth_annotation_tables.R
## ------------------------------------------------------------------------
## Aggregates the 100kb-bin annotation tables into coarser genomic windows
## to shrink file size / memory footprint (to fit the Shiny free tier),
## while respecting that this is an ANNOTATION track, not a plain signal:
##
## - Continuous columns (observed, prediction) are averaged (mean) within
##   each window. This is genuine smoothing: it also reduces per-bin noise,
##   not just row count.
##
## - Categorical columns (annot_final, top1) are reduced with a "priority
##   mode" instead of a plain majority vote: if a window contains ANY class
##   other than the configured "background" classes (see BACKGROUND_CLASSES
##   below), the most frequent NON-background class wins, even if the
##   background classes are individually more numerous in that window.
##   Plain majority voting would let a single interesting, isolated 100kb
##   bin (e.g. "Positive Selection") disappear under a sea of surrounding
##   "No Detectable Force"/"Occurrence" bins once you coarsen the
##   resolution - priority-mode keeps that signal visible instead.
##
## Output: for every existing Data/annotation_<ampl|del>_<cluster>.tsv,
## writes Data/annotation_<ampl|del>_<cluster>_smoothed.tsv (new files,
## originals untouched) plus a Data/smoothing_summary.tsv with before/after
## row counts and file sizes so you can judge whether the reduction is
## enough before switching the app over to the smoothed files.
##
## Usage: adjust WINDOW_BP below if you want a different reduction factor,
## then run:  Rscript smooth_annotation_tables.R
## ------------------------------------------------------------------------

library(data.table)

# Merge every ~WINDOW_BP/100000 original 100kb bins into one output bin.
# 500,000 -> ~5x fewer rows. 1,000,000 -> ~10x fewer rows. Increase this if
# the app is still too heavy after a first pass.
WINDOW_BP <- 500000

DATA_DIR   <- "Data"
SUFFIX_OUT <- "_smoothed"

model_types <- c("Mid-length", "Small-scale", "Chromosome-level", "Arm-level", "No-cluster")
models      <- c("ampl", "del")

# Keep this in sync with NEEDED_COLUMNS in 0_LoadData.R - these are the
# only columns the app actually reads.
NEEDED_COLUMNS <- c(
  "chr", "start", "end", "binID", "type",
  "observed", "prediction",
  "top1", "annot_final"
)

# Classes treated as "uninteresting background" when reducing a categorical
# column. Leave a vector empty to fall back to plain majority voting for
# that column (used for top1, whose vocabulary isn't a fixed, known set).
BACKGROUND_CLASSES <- list(
  annot_final = c("No Detectable Force", "Occurrence"),
  top1        = character(0)
)

## ---- helpers -------------------------------------------------------------

read_raw_table <- function(path) {
  
  dt <- data.table::fread(path, sep = "\t", nThread = data.table::getDTthreads())
  
  # Same defensive header handling as 0_LoadData.R: trim whitespace, match
  # column names case-insensitively, keep only what we need (this also
  # transparently drops the stray leading row-index column some of these
  # files have, since it never matches a name in NEEDED_COLUMNS).
  data.table::setnames(dt, trimws(names(dt)))
  current <- names(dt)
  lower_current <- tolower(current)
  
  for (w in NEEDED_COLUMNS) {
    idx <- match(tolower(w), lower_current)
    if (is.na(idx)) {
      stop(
        "Column \"", w, "\" not found in ", path, ".\n",
        "Columns found: ", paste(current, collapse = ", ")
      )
    }
    if (current[idx] != w) data.table::setnames(dt, current[idx], w)
  }
  
  dt <- dt[, NEEDED_COLUMNS, with = FALSE]
  
  if (!all(startsWith(as.character(dt$chr), "chr"))) {
    dt[, chr := paste0("chr", chr)]
  }
  
  dt
}

#' Priority-mode: the most frequent value in x that is NOT in `background`,
#' if any such value is present; otherwise the plain majority value.
priority_mode <- function(x, background) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  
  tab <- sort(table(x), decreasing = TRUE)
  candidates <- names(tab)
  
  non_bg <- candidates[!(candidates %in% background)]
  if (length(non_bg)) return(non_bg[1])
  
  candidates[1]
}

smooth_one_table <- function(dt, window_bp) {
  
  dt <- data.table::copy(dt)
  dt[, window_id := start %/% window_bp]
  
  agg <- dt[, list(
    start         = min(start),
    end           = max(end),
    observed      = mean(observed, na.rm = TRUE),
    prediction    = mean(prediction, na.rm = TRUE),
    annot_final   = priority_mode(annot_final, BACKGROUND_CLASSES$annot_final),
    top1          = priority_mode(top1, BACKGROUND_CLASSES$top1),
    n_bins_merged = .N
  ), by = list(chr, type, window_id)]
  
  agg[, binID := paste0(chr, "_", window_id)]
  agg[, window_id := NULL]
  
  data.table::setcolorder(
    agg,
    c("chr", "start", "end", "binID", "type", "observed", "prediction", "top1", "annot_final", "n_bins_merged")
  )
  data.table::setorder(agg, chr, start)
  
  agg
}

## ---- main ------------------------------------------------------------

results <- list()

for (model_type in model_types) {
  for (model in models) {
    
    in_path  <- file.path(DATA_DIR, paste0("annotation_", model, "_", model_type, ".tsv"))
    out_path <- file.path(DATA_DIR, paste0("annotation_", model, "_", model_type, SUFFIX_OUT, ".tsv"))
    
    if (!file.exists(in_path)) {
      message("Skipping (not found): ", in_path)
      next
    }
    
    message("Processing ", in_path, " ...")
    
    dt <- read_raw_table(in_path)
    smoothed <- smooth_one_table(dt, WINDOW_BP)
    
    data.table::fwrite(smoothed, out_path, sep = "\t")
    
    key <- paste(model_type, model)
    results[[key]] <- list(
      in_rows  = nrow(dt),
      out_rows = nrow(smoothed),
      in_size  = file.info(in_path)$size,
      out_size = file.info(out_path)$size
    )
    
    message(sprintf(
      "  %d -> %d rows (%.1fx fewer), %.1f MB -> %.1f MB",
      nrow(dt), nrow(smoothed), nrow(dt) / nrow(smoothed),
      file.info(in_path)$size / 1e6, file.info(out_path)$size / 1e6
    ))
  }
}

if (length(results)) {
  
  summary_dt <- data.table::rbindlist(lapply(names(results), function(k) {
    r <- results[[k]]
    data.table::data.table(
      table    = k,
      in_rows  = r$in_rows,
      out_rows = r$out_rows,
      reduction_x = round(r$in_rows / r$out_rows, 1),
      in_MB    = round(r$in_size / 1e6, 1),
      out_MB   = round(r$out_size / 1e6, 1)
    )
  }))
  
  print(summary_dt)
  data.table::fwrite(summary_dt, file.path(DATA_DIR, "smoothing_summary.tsv"), sep = "\t")
  
} else {
  
  message("No source .tsv files found under '", DATA_DIR, "' matching annotation_<ampl|del>_<cluster>.tsv")
}