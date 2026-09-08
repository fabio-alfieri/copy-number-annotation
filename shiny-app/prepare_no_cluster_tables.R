#!/usr/bin/env Rscript
## ------------------------------------------------------------------------
## prepare_no_cluster_tables.R
## ------------------------------------------------------------------------
## Converte i due file "All CNAs" (ampl / del) nello STESSO formato e con lo
## STESSO schema di naming delle altre tabelle di annotazione, cioe':
##
##   Data/annotation_ampl_No-cluster.tsv
##   Data/annotation_del_No-cluster.tsv
##
## Differenze note fra i due formati di input:
##   - i file "All CNAs" hanno due colonne in piu' (annot_3_classes,
##     annot_5_classes) che le altre tabelle non hanno -> vengono scartate;
##   - annot_final usa nei file "All CNAs" il vocabolario di top1
##     ("Incorrect prediction - Positive/Negative Selection") invece di
##     quello delle altre tabelle ("Excess of Observation/Prediction -
##     Likely ... Selection") -> viene ricodificata (vedi
##     ANNOT_FINAL_RECODE), con verifica sul segno di observed - prediction;
##   - alcuni file sono stati scritti con write.table(row.names = TRUE),
##     quindi hanno una colonna indice iniziale SENZA nome -> fread la
##     riconosce (header = "auto") e qui viene semplicemente scartata perche'
##     non corrisponde a nessun nome dello schema canonico.
## Il resto delle colonne e' identico (stessi nomi, stesso ordine).
##
## Lo schema canonico viene letto, quando possibile, direttamente
## dall'header di una tabella gia' esistente (REFERENCE_TABLE), cosi' se
## quello schema cambia in futuro non c'e' nulla da aggiornare qui.
##
## Uso:
##   Rscript prepare_no_cluster_tables.R
##
## Poi, per generare le versioni smoothed effettivamente lette dall'app:
##   Rscript smooth_annotation_tables.R No-cluster
## ------------------------------------------------------------------------

library(data.table)

## ---- configurazione ------------------------------------------------------

DATA_DIR     <- "Data"
CLUSTER_NAME <- "No-cluster"          # nome interno usato dall'app
# (etichetta mostrata: "All CNAs")

# Se l'auto-rilevamento non trova i file (o ne trova piu' di uno), indica qui
# i percorsi a mano, es.:
#   INPUT_FILES <- list(ampl = "Data/All_CNAs_ampl.tsv",
#                       del  = "Data/All_CNAs_del.tsv")
INPUT_FILES <- list(ampl = NULL, del = NULL)

# Tabella di riferimento da cui leggere lo schema canonico delle colonne.
REFERENCE_TABLE <- file.path(DATA_DIR, "annotation_ampl_Mid-length.tsv")

# Schema usato solo se REFERENCE_TABLE non esiste (e' l'header che hai
# incollato per le tabelle "normali").
FALLBACK_COLUMNS <- c(
  "chr", "start", "end", "binID", "type", "labels",
  "observed", "prediction", "residual", "pval.residual",
  "pred_is_correct", "region_is_diploid",
  "top1", "shap_top1",
  "is_candidate", "sas_ns", "is_strong_candidate",
  "annot_final"
)

# Colonne che l'app legge davvero (deve restare allineato a NEEDED_COLUMNS
# in 0_LoadData.R e in smooth_annotation_tables.R): se una di queste manca,
# lo script si ferma subito invece di produrre un file inutilizzabile.
NEEDED_COLUMNS <- c(
  "chr", "start", "end", "binID", "type",
  "observed", "prediction",
  "top1", "annot_final"
)

# --- ricodifica di annot_final -------------------------------------------
# Le tabelle "All CNAs" usano in annot_final il vocabolario di top1
# ("Incorrect prediction - ..."), mentre le altre tabelle usano la forma
# "Excess of ... - Likely ...". E' lo stesso concetto: predizione errata,
# con la direzione data dal segno di (observed - prediction).
#   observed > prediction  -> eccesso di OSSERVATO  -> likely positive
#   observed < prediction  -> eccesso di PREDETTO   -> likely negative
# NB: la ricodifica riguarda SOLO annot_final. top1 resta invariata, perche'
# li' la forma "Incorrect prediction - ..." e' corretta anche nelle tabelle
# standard (verificato sulle righe di esempio delle tabelle gia' in uso).
ANNOT_FINAL_RECODE <- c(
  "Incorrect prediction - Positive Selection" = "Excess of Observation - Likely Positive Selection",
  "Incorrect prediction - Negative Selection" = "Excess of Prediction - Likely Negative Selection"
)

# Vocabolario atteso in uscita: dopo la ricodifica non dovrebbe restare
# nulla fuori da questa lista (stesse 6 classi delle altre tabelle).
CANONICAL_ANNOT_FINAL <- c(
  "Negative Selection",
  "Excess of Prediction - Likely Negative Selection",
  "Positive Selection",
  "Excess of Observation - Likely Positive Selection",
  "Occurrence",
  "No Detectable Force"
)

# Pattern per riconoscere automaticamente i file di input.
INPUT_PATTERN  <- "(all[^a-z0-9]*cna|no[^a-z0-9]*cluster)"
MODEL_PATTERNS <- list(ampl = "ampl", del = "del")

## ---- helper --------------------------------------------------------------

output_path <- function(model) {
  file.path(DATA_DIR, paste0("annotation_", model, "_", CLUSTER_NAME, ".tsv"))
}

#' Trova il file di input per un modello, o restituisce quello forzato a mano.
find_input_file <- function(model) {
  
  if (!is.null(INPUT_FILES[[model]])) {
    p <- INPUT_FILES[[model]]
    if (!file.exists(p)) stop("File indicato in INPUT_FILES non trovato: ", p)
    return(p)
  }
  
  files <- list.files(DATA_DIR, pattern = "\\.tsv$", full.names = FALSE)
  if (!length(files)) stop("Nessun .tsv in '", DATA_DIR, "'.")
  
  nm      <- tolower(files)
  out_nm  <- tolower(basename(output_path(model)))
  
  keep <- grepl(INPUT_PATTERN, nm) &
    grepl(MODEL_PATTERNS[[model]], nm) &
    !grepl("_smoothed", nm, fixed = TRUE) &
    nm != out_nm
  
  hits <- files[keep]
  
  if (length(hits) == 0) {
    stop(
      "Nessun file 'All CNAs' trovato per il modello '", model, "' in '", DATA_DIR, "'.\n",
      "File presenti: ", paste(files, collapse = ", "), "\n",
      "Indica il percorso a mano nella lista INPUT_FILES in cima allo script."
    )
  }
  
  if (length(hits) > 1) {
    stop(
      "Piu' di un candidato per il modello '", model, "': ",
      paste(hits, collapse = ", "), "\n",
      "Indica quale usare nella lista INPUT_FILES in cima allo script."
    )
  }
  
  file.path(DATA_DIR, hits)
}

#' Schema canonico: header della tabella di riferimento, altrimenti fallback.
canonical_columns <- function() {
  if (file.exists(REFERENCE_TABLE)) {
    hdr <- names(data.table::fread(REFERENCE_TABLE, nrows = 0))
    hdr <- trimws(hdr)
    # scarta l'eventuale colonna indice iniziale senza nome (fread la
    # nomina V1 quando l'header e' disallineato)
    hdr <- hdr[nzchar(hdr) & hdr != "V1"]
    if (length(hdr)) {
      message("Schema canonico letto da ", REFERENCE_TABLE, " (", length(hdr), " colonne).")
      return(hdr)
    }
  }
  message("REFERENCE_TABLE non disponibile: uso FALLBACK_COLUMNS.")
  FALLBACK_COLUMNS
}

#' Rinomina case-insensitive + selezione delle sole colonne richieste.
#' Stessa logica difensiva di 0_LoadData.R::select_and_rename_columns().
align_columns <- function(dt, wanted, source_path, required) {
  
  data.table::setnames(dt, trimws(names(dt)))
  current       <- names(dt)
  lower_current <- tolower(current)
  
  found <- character(0)
  for (w in wanted) {
    idx <- match(tolower(w), lower_current)
    if (is.na(idx)) next
    if (current[idx] != w) data.table::setnames(dt, current[idx], w)
    found <- c(found, w)
  }
  
  missing_required <- setdiff(required, found)
  if (length(missing_required)) {
    stop(
      "Colonne obbligatorie mancanti in ", source_path, ": ",
      paste(missing_required, collapse = ", "), "\n",
      "Colonne trovate: ", paste(current, collapse = ", ")
    )
  }
  
  missing_optional <- setdiff(wanted, found)
  if (length(missing_optional)) {
    message("  ATTENZIONE - colonne dello schema assenti in input (verranno omesse): ",
            paste(missing_optional, collapse = ", "))
  }
  
  dropped <- setdiff(current, c(found, ""))
  if (length(dropped)) {
    message("  Colonne scartate (non nello schema canonico): ",
            paste(dropped, collapse = ", "))
  }
  
  dt[, found, with = FALSE]
}

## ---- conversione ---------------------------------------------------------

CANONICAL <- canonical_columns()

stopifnot(all(NEEDED_COLUMNS %in% CANONICAL))

summary_rows <- list()

for (model in c("ampl", "del")) {
  
  in_path  <- find_input_file(model)
  out_path <- output_path(model)
  
  message("\n[", model, "] ", in_path, "  ->  ", out_path)
  
  # header = "auto" (default), NON forzato a TRUE: e' cio' che permette a
  # fread di riconoscere l'eventuale colonna indice iniziale senza nome.
  dt <- data.table::fread(in_path, sep = "\t", nThread = data.table::getDTthreads())
  
  n_in <- nrow(dt)
  
  dt <- align_columns(dt, CANONICAL, in_path, NEEDED_COLUMNS)
  
  # normalizzazione minima, identica a quella fatta a valle dall'app
  if (!all(startsWith(as.character(dt$chr), "chr"))) {
    dt[, chr := paste0("chr", chr)]
  }
  
  # --- controlli di sanita' ------------------------------------------------
  n_bad_coord <- sum(dt$start >= dt$end, na.rm = TRUE)
  if (n_bad_coord) message("  ATTENZIONE - ", n_bad_coord, " righe con start >= end.")
  
  n_na <- sum(is.na(dt$observed) | is.na(dt$prediction))
  if (n_na) message("  ATTENZIONE - ", n_na, " righe con observed/prediction NA.")
  
  extra_chr <- setdiff(unique(dt$chr), paste0("chr", 1:22))
  if (length(extra_chr)) {
    message("  NOTA - cromosomi fuori da chr1-chr22 presenti (l'app li ignora): ",
            paste(extra_chr, collapse = ", "))
  }
  
  # --- ricodifica di annot_final -------------------------------------------
  dt[, annot_final := as.character(annot_final)]
  
  to_recode <- dt$annot_final %in% names(ANNOT_FINAL_RECODE)
  n_recoded <- sum(to_recode)
  
  if (n_recoded) {
    
    # Verifica della direzione PRIMA di ricodificare: la classe attesa in
    # base al segno di (observed - prediction) deve coincidere con quella
    # della mappa. Se non coincide, la mappa e' sbagliata (o il file usa
    # una convenzione diversa) e lo script si ferma.
    src      <- dt$annot_final[to_recode]
    mapped   <- unname(ANNOT_FINAL_RECODE[src])
    by_sign  <- ifelse(
      dt$observed[to_recode] >= dt$prediction[to_recode],
      "Excess of Observation - Likely Positive Selection",
      "Excess of Prediction - Likely Negative Selection"
    )
    
    n_mismatch <- sum(mapped != by_sign, na.rm = TRUE)
    if (n_mismatch) {
      bad <- head(which(mapped != by_sign), 5)
      stop(
        "Ricodifica di annot_final incoerente in ", in_path, ": ", n_mismatch,
        " righe su ", n_recoded, " hanno un segno di (observed - prediction) ",
        "opposto alla classe attesa.\n",
        "Esempi (annot_final originale -> atteso dalla mappa / atteso dal segno):\n",
        paste0("  ", src[bad], " -> ", mapped[bad], " / ", by_sign[bad], collapse = "\n"), "\n",
        "Controlla ANNOT_FINAL_RECODE in cima allo script prima di procedere."
      )
    }
    
    dt[to_recode, annot_final := unname(ANNOT_FINAL_RECODE[annot_final])]
    message("  Ricodificate ", n_recoded, " righe di annot_final (",
            round(100 * n_recoded / nrow(dt), 1), "%), direzione verificata sul segno del residuo.")
  }
  
  classes_annot <- sort(unique(dt$annot_final))
  message("  Classi in annot_final (", length(classes_annot), "): ",
          paste(classes_annot, collapse = " | "))
  
  unexpected <- setdiff(classes_annot, CANONICAL_ANNOT_FINAL)
  if (length(unexpected)) {
    message("  ATTENZIONE - classi fuori dal vocabolario delle altre tabelle: ",
            paste(unexpected, collapse = " | "), "\n",
            "  Vanno aggiunte a ANNOT_FINAL_RECODE (se sono sinonimi) oppure ",
            "alla palette e alla UI dell'app (se sono classi nuove).")
  }
  message("  Tipi tumorali (", length(unique(dt$type)), "): ",
          paste(sort(unique(dt$type)), collapse = ", "))
  
  data.table::setorder(dt, chr, start)
  
  data.table::fwrite(dt, out_path, sep = "\t")
  
  message("  Scritte ", nrow(dt), " righe / ", ncol(dt), " colonne (input: ",
          n_in, " righe), ", round(file.info(out_path)$size / 1e6, 1), " MB.")
  
  summary_rows[[model]] <- data.table::data.table(
    model         = model,
    input         = basename(in_path),
    output        = basename(out_path),
    rows          = nrow(dt),
    cols          = ncol(dt),
    out_MB        = round(file.info(out_path)$size / 1e6, 1),
    annot_classes = paste(classes_annot, collapse = " | ")
  )
}

summary_dt <- data.table::rbindlist(summary_rows)
print(summary_dt[, .(model, input, output, rows, cols, out_MB)])

data.table::fwrite(summary_dt, file.path(DATA_DIR, "no_cluster_conversion_summary.tsv"), sep = "\t")

message("\nFatto. Passo successivo:  Rscript smooth_annotation_tables.R ", CLUSTER_NAME)