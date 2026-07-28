# merge objects
rm(list=ls())
gc(full=T)
packages <- c(
  "stringr", "reshape2", "tidyr",
  "ggplot2","fastcluster", "tidyverse"
)
installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)

model.outputs <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/TTS", full.names = TRUE)

##### extract cts_label from filename (token "TTS") and keep only SAME cancer-type pairs
extract_cts_label <- function(path){
  fname.parts <- strsplit(basename(path), "_")[[1]]
  tts.idx <- which(fname.parts == "TTS")
  length.idx <- which(fname.parts == "Length")
  if(length(tts.idx) == 0 || length(length.idx) == 0 || length.idx <= tts.idx + 1){
    return(NA_character_)
  }
  cts.tokens <- fname.parts[(tts.idx + 1):(length.idx - 1)]
  if(length(cts.tokens) != 2) return(NA_character_)
  if(cts.tokens[1] != cts.tokens[2]) return(NA_character_) # drop cross-type pairs
  paste(cts.tokens, collapse = "_")
}

cts_labels_by_file <- sapply(model.outputs, extract_cts_label)
model.outputs <- model.outputs[!is.na(cts_labels_by_file)]
cts_labels_by_file <- cts_labels_by_file[!is.na(cts_labels_by_file)]
names(cts_labels_by_file) <- NULL

#patterns <- c('Arm-level','Chromosome-level','Mid-length','Small-scale','no_cluster')
patterns <- "Mid-length"

shap.outputs <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/SHAP/TTS/", full.names = TRUE)
shap.outputs <- grep("::", shap.outputs, value = TRUE)

for(i in patterns){
  print(i)
  
  ##### subset to this Length pattern, Amplification model only, same-cancer-type files only
  pattern.idx <- grep(i, model.outputs)
  ampl.idx <- grep('ampl', model.outputs)
  keep.idx <- intersect(pattern.idx, ampl.idx)
  
  pattern.files <- model.outputs[keep.idx]
  pattern.cts_labels <- cts_labels_by_file[keep.idx]
  
  if(length(pattern.files) == 0){
    warning(paste("No matching amplification files (same cancer-type pair) found for pattern:", i))
    next
  }
  
  ##### loop over every same-type cts_label found for this pattern
  for(k in seq_along(pattern.files)){
    
    cts_label <- pattern.cts_labels[k]
    print(cts_label)
    
    # Feature Matrix Values
    load(pattern.files[k])
    
    ##### index through the cts_label layer, Amplification model only
    if(!("ampl_score" %in% names(Output.regressor))){
      warning(paste("ampl_score not found in", basename(pattern.files[k]), "- skipping"))
      rm(Output.regressor)
      next
    }
    model_ampl <- Output.regressor$ampl_score
    
    if(!(cts_label %in% names(model_ampl))){
      warning(paste("cts_label", cts_label, "not found in ampl_score for file", basename(pattern.files[k]), "- skipping"))
      rm(Output.regressor, model_ampl)
      next
    }
    model_ampl <- model_ampl[[cts_label]]
    
    rm(Output.regressor)
    
    models.X <- list(ampl = rbind(model_ampl[[1]]$X_train, model_ampl[[2]]$X_test))
    names(models.X) <- paste0(i,'::Amplification model::',cts_label)
    
    # SHAP values
    ##### match against three-part SHAP filenames (class::model::cts_label)
    shap.match <- shap.outputs[grepl(paste0(i,'::'), shap.outputs) &
                                 grepl('Amplification', shap.outputs) &
                                 grepl(paste0('::',cts_label,'_'), shap.outputs)]
    
    if(length(shap.match) == 0){
      warning(paste("No matching SHAP file found for", i, cts_label, "- skipping"))
      next
    }
    
    load(shap.match[1])
    shap_ampl <- models.shap.df[[1]] %>% dplyr::select(-BIAS)
    rm(models.shap.df, models.X.test, models.X.train, models.shap.values)
    
    models.shap.df <- list(ampl = shap_ampl)
    names(models.shap.df) <- paste0(i,'::Amplification model::',cts_label)
    
    # aggregate and save
    feature_and_shap <- list(models.X,models.shap.df)
    names(feature_and_shap) <- c('models.X','models.shap.df')
    
    write_rds(feature_and_shap, file = paste0("/home/ieo7429/Scrivania/results_regressor_gab/SHAP_and_featurematrix/TTS/SHAP_and_FeatureMatrix_",i,"_",cts_label,"_Ampl.rds"))
    
  }
  
}
