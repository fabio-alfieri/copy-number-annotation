## ------------------------------------------------------------------------
## Data loading
## ------------------------------------------------------------------------
## Loads the FULL (non-subsampled) annotation tables, lazily and with a
## minimal memory footprint so the app can run on the Shiny (shinyapps.io)
## free tier.
##
## File naming
## -----------
## Tables are expected at:
##   Data/annotation_<ampl|del>_<model_type>.tsv
## e.g. Data/annotation_ampl_Mid-length.tsv, Data/annotation_del_Small-scale.tsv
## No renaming is needed - just drop the files you were given into Data/
## with their current names.
##
## Cluster categories currently available: Mid-length, Small-scale,
## Chromosome-level, Arm-level, No-cluster (shown in the UI as "All CNAs" -
## see CLUSTER_LABELS below). The No-cluster tables come from a source with
## a slightly different schema and are normalised into the layout used by
## every other table by prepare_no_cluster_tables.R, which writes
## Data/annotation_<ampl|del>_No-cluster.tsv; smooth_annotation_tables.R
## then produces the _smoothed files this script actually reads.
##
## Memory optimizations
## ---------------------
## 1) Column pruning: the app only ever reads chr, start, end, binID, type,
##    observed, prediction, top1 and annot_final (see NEEDED_COLUMNS below -
##    cross-checked against every $col / [["col"]] access in the app).
##    Everything else in the source tables (labels, residual, pval.residual,
##    pred_is_correct, region_is_diploid, shap_top1, annot_3_classes,
##    annot_5_classes, is_candidate, sas_ns, is_strong_candidate, the stray
##    row-index column, ...) is dropped right after parsing and never held
##    in memory. On these files this cuts memory roughly in half.
## 2) Lazy loading: the model_type x model tables are read only the first
##    time a user selection actually needs them (`meta_list_get()` below),
##    not all upfront, and are cached in memory for the life of the R
##    process afterwards (`.table_cache`).
## 3) .rds cache, usable WITHOUT the source .tsv: on first parse of a table,
##    a compact binary copy is written next to the source file as
##    Data/annotation_<ampl|del>_<model_type>.rds. On every subsequent load
##    (including a fresh app restart) that .rds is used directly if present,
##    regardless of whether the .tsv is even there. This means that once
##    you've generated the .rds files once (e.g. by running the app
##    locally), you can deploy ONLY the small .rds files to shinyapps.io -
##    skipping the much larger raw .tsv files entirely, which helps with
##    both the free tier's bundle-size limit and cold-start memory/time.
##    (If Data/ is read-only at runtime - as the deployed app bundle is on
##    shinyapps.io - writing a NEW .rds simply fails silently and the app
##    keeps working from the .tsv; it just means you generate the .rds
##    locally first and upload it, rather than letting the deployed app
##    write it for you.)
## ------------------------------------------------------------------------

library(data.table)

model_types <- c("Mid-length", "Small-scale", "Chromosome-level", "Arm-level", "No-cluster")
lis_names   <- c("ampl", "del")

# Labels shown to the user. The INTERNAL identifier stays the one used in
# the file names ("No-cluster"); the label is only for the UI dropdown and
# for plot titles/axis labels.
CLUSTER_LABELS <- c(
  "Mid-length"       = "Mid-length",
  "Small-scale"      = "Small-scale",
  "Arm-level"        = "Arm-level",
  "Chromosome-level" = "Chromosome-level",
  "No-cluster"       = "All CNAs"
)

cluster_label <- function(x) {
  unname(ifelse(x %in% names(CLUSTER_LABELS), CLUSTER_LABELS[x], x))
}

# annot_final vocabulary, in a FIXED order: this is what guarantees that a
# given class always gets the same colour across different selections
# (the palette in landscape_plot_observed_prediction.R is indexed by
# position in this vector).
#
# The "All CNAs"/No-cluster source tables originally wrote the top1
# vocabulary ("Incorrect prediction - Positive/Negative Selection") into
# annot_final; prepare_no_cluster_tables.R recodes those to the two
# "Excess of ..." classes below, so nothing extra is needed here.
ANNOT_FINAL_CLASSES <- c(
  "Negative Selection",
  "Excess of Prediction - Likely Negative Selection",
  "Positive Selection",
  "Excess of Observation - Likely Positive Selection",
  "Occurrence",
  "No Detectable Force"
)

# The only columns the app ever reads (verified against every $col /
# [["col"]] access across app.R, ui.R, 1_dynamic_plotting_functions.R,
# landscape_plot_observed_prediction.R and plotting_helpers.R). Extend this
# if a future feature needs another column from the source tables.
NEEDED_COLUMNS <- c(
  "chr", "start", "end", "binID", "type",
  "observed", "prediction",
  "top1", "annot_final"     # the two selectable "Select annotation to plot" options
)

## ---- small reference tables: safe to load eagerly -----------------------

centromere_table <- data.table::fread(
  "Data/cytoBand.txt",
  header    = FALSE,
  col.names = c("chr", "start", "end", "cytoband", "type"),
  nThread   = data.table::getDTthreads()
)
centromere_table <- centromere_table[type == "acen"]
centromere_table <- as.data.frame(centromere_table)

load("Data/All_levels_backbonetables.RData")
backbone.100kbp <- chr_backbone_namesfixed[["0.1Mbp"]]
backbone.100kbp <- data.table::rbindlist(backbone.100kbp)
backbone.100kbp[, binID := paste0(chr, "_", bin)]
backbone.100kbp[, chr := paste0("chr", chr)]
backbone.100kbp <- as.data.frame(backbone.100kbp)

backbone.100kbp_granges <- GRanges(
  seqnames = backbone.100kbp$chr,
  ranges   = IRanges(start = backbone.100kbp$start_bin, end = backbone.100kbp$end_bin),
  binID    = backbone.100kbp$binID
)

## ---- large annotation tables: lazy-loaded + cached -----------------------

.table_cache <- new.env(parent = emptyenv())

#' Case/whitespace-insensitive column rename-in-place, keeping only the
#' requested columns (dropping everything else - this is where the memory
#' savings from column pruning actually happen). Errors clearly, listing
#' what was actually found, if a required column is genuinely missing -
#' rather than a cryptic "object 'x' not found" three call-frames away.
select_and_rename_columns <- function(dt, wanted, source_path) {
  
  data.table::setnames(dt, trimws(names(dt)))
  current <- names(dt)
  lower_current <- tolower(current)
  
  for (w in wanted) {
    idx <- match(tolower(w), lower_current)
    if (is.na(idx)) {
      stop(
        "Expected column \"", w, "\" not found in ", source_path, ".\n",
        "Columns found: ", paste(current, collapse = ", "), "\n",
        "This usually means the header and data columns are misaligned ",
        "(e.g. the file has a leading, unnamed row-index column) or the ",
        "column was renamed upstream. Check the first two lines of the ",
        "file with readLines(path, n = 2)."
      )
    }
    if (current[idx] != w) data.table::setnames(dt, current[idx], w)
  }
  
  dt[, wanted, with = FALSE]
}

read_annotation_table <- function(model_type, model) {
  
  cache_key <- paste(model_type, model, sep = "__")
  cached <- get0(cache_key, envir = .table_cache, inherits = FALSE)
  if (!is.null(cached)) return(cached)
  
  tsv_path <- file.path("Data", paste0("annotation_", model, "_", model_type, "_smoothed.tsv"))
  rds_path <- file.path("Data", paste0("annotation_", model, "_", model_type, "_smoothed.rds"))
  
  # The .rds cache is checked FIRST and is usable on its own: if it exists
  # (and there either is no .tsv to compare against, or it's at least as
  # recent as the .tsv), we never need to touch the .tsv at all - which is
  # exactly what lets you deploy a lean bundle with only .rds files.
  dt <- NULL
  if (file.exists(rds_path)) {
    rds_is_current <- !file.exists(tsv_path) ||
      file.info(rds_path)$mtime >= file.info(tsv_path)$mtime
    if (rds_is_current) {
      dt <- tryCatch(readRDS(rds_path), error = function(e) NULL)
    }
  }
  
  if (is.null(dt)) {
    
    if (!file.exists(tsv_path)) {
      stop(
        "Annotation table not found (looked for both\n  ", rds_path, "\nand\n  ", tsv_path, ")."
      )
    }
    
    # NOTE: header is intentionally left as "auto" (fread's default) rather
    # than forced to TRUE. These files are written with an unnamed leading
    # row-index column (write.table()'s default row.names = TRUE). fread's
    # automatic header detection recognises that mismatch and shifts the
    # real header into alignment (adding a guessed row-index column) - but
    # that heuristic only runs when header is "auto"; forcing header = TRUE
    # disables it and silently misreads the first DATA row as the header
    # instead.
    #
    # select = NEEDED_COLUMNS is intentionally NOT used here even though it
    # would save some parse-time memory too: fread's select matching is
    # case-sensitive and silently DROPS a requested column with a near-miss
    # name (e.g. "Chr" vs "chr") instead of erroring, which previously
    # caused a hard-to-diagnose "object not found" failure downstream. We
    # read every column, resolve names defensively, THEN prune - slightly
    # more peak memory during this one-off parse, but no silent data loss.
    dt <- data.table::fread(
      tsv_path,
      sep     = "\t",
      nThread = data.table::getDTthreads()
    )
    
    dt <- select_and_rename_columns(dt, NEEDED_COLUMNS, tsv_path)
    
    if (!all(startsWith(as.character(dt$chr), "chr"))) {
      dt[, chr := paste0("chr", chr)]
    }
    
    # factor-encode chr/type purely to get a numeric chr1..chr22 sort order
    # via setorder() below; both are converted back to character right
    # after so every downstream comparison (%in%, ==, etc.) keeps working
    # exactly as before.
    chr_levels <- paste0("chr", 1:22)
    dt[, chr := factor(chr, levels = chr_levels)]
    dt[, type := factor(type)]
    
    data.table::setorder(dt, chr, start)
    dt[, chr := as.character(chr)]
    dt[, type := as.character(type)]
    
    # Persist a fast, ALREADY-PRUNED binary copy next to the source file.
    # This is a pure optimization: if Data/ is read-only at runtime (e.g.
    # the deployed bundle on shinyapps.io), the write fails silently and
    # the app just re-parses the .tsv next time. Run the app once locally
    # to generate these .rds files, then deploy them alongside (or instead
    # of) the .tsv files.
    tryCatch(
      saveRDS(dt, rds_path, compress = FALSE),
      error = function(e) invisible(NULL)
    )
  }
  
  data.table::setDT(dt)
  # setindex() (not setkey()): setkey() physically re-sorts the table by
  # its key columns using a plain (lexicographic) sort on character
  # columns, which would undo the chr1, chr2, ..., chr22 numeric ordering
  # set above via setorder(). setindex() builds a lookup index for fast
  # filtering without touching row order.
  data.table::setindex(dt, type, chr)
  
  assign(cache_key, dt, envir = .table_cache)
  dt
}

#' Lazily fetch one (model_type x model) annotation table.
#' Only reads/parses a table on first request; subsequent requests for the
#' same table within the same R process are served from `.table_cache`.
meta_list_get <- function(model_type, model) {
  if (!(model_type %in% model_types)) stop("Unknown model_type: ", model_type)
  if (!(model %in% lis_names)) stop("Unknown model: ", model)
  read_annotation_table(model_type, model)
}