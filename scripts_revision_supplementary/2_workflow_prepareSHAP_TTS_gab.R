rm(list = ls())
gc(full=T)
library(reshape2)
library(tidyverse)
#Interactome INSIDER (alternative to HIPPIE)
model.outputs <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/TTS/", full.names = TRUE)
##### GAB EDITED HERE
# model.outputs <- grep("2025-08", model.outputs, value = TRUE) # selected only models from august 2025
# model.outputs <- model.outputs[-grep("NoCluster",model.outputs)]
# model.outputs <- model.outputs[c(1,2,3,4,13,14,15,16)]
for(file in model.outputs){
  
  ##### NEW: extract cts_label (the two cancer-type tokens between "TTS" and "Length")
  # e.g. ".../Output.regressor_TTS_STAD_COADREAD_Length_Mid-length_ampl_..." -> "STAD_COADREAD"
  fname.parts <- strsplit(basename(file), "_")[[1]]
  tts.idx <- which(fname.parts == "TTS")
  length.idx <- which(fname.parts == "Length")
  
  if(length(tts.idx) == 0 || length(length.idx) == 0 || length.idx <= tts.idx + 1){
    warning(paste("Could not extract cts_label from filename:", basename(file), "- skipping"))
    next
  }
  
  cts.tokens <- fname.parts[(tts.idx + 1):(length.idx - 1)]
  
  if(length(cts.tokens) != 2){
    warning(paste("Unexpected number of cancer-type tokens in filename:", basename(file), "- skipping"))
    next
  }
  
  cts_label <- paste(cts.tokens, collapse = "_")
  
  ##### NEW: keep only SAME cancer-type pairs (e.g. STAD_STAD), skip cross pairs (e.g. STAD_COADREAD)
  if(cts.tokens[1] != cts.tokens[2]){
    next
  }
  
  print(cts_label)
  
  models.shap.df <- list()
  models.shap.values <- list()
  
  models.X.test <- list()
  models.X.train <- list()
  
  name <- strsplit(gsub("_0\\.1Mbp.*$", "", basename(file)),"_")
  print(name)
  
  model <- name[[1]][length(name[[1]])]
  class <- gsub("Output.regressor_Length_|Output.regressor_Length_Location_","",gsub("_0\\.1Mbp.*$", "", basename(file)))
  class <- gsub(paste0("_",model),"",class)
  if(model == "ampl"){
    model <- "Amplification model"
  }else{
    model <- "Deletion model"
  }
  
  load(file)
  
  ##### NEW: index into the cts_label layer of Output.regressor
  if(!(cts_label %in% names(Output.regressor[[1]]))){
    warning(paste("cts_label", cts_label, "not found in Output.regressor[[1]] for file", basename(file), "- skipping"))
    next
  }
  Output.regressor.cts <- Output.regressor[[1]][[cts_label]]
  
  X_test <- Output.regressor.cts[[1]]$X_test
  X_train <- Output.regressor.cts[[1]]$X_train
  
  # NOTA per gab del futuro. In 2_runRegressor_gab.R line 12 we know that train and test are always the same
  
  models.X.test[[paste(class,model,cts_label,sep = "::")]] <- X_test
  models.X.train[[paste(class,model,cts_label,sep = "::")]] <- X_train
  
  avg.shap.df <- c()
  avg.shap.values <- c()
  
  for(iter in 1:10){
    
    # SHAP df
    
    shap.df_train <- Output.regressor.cts[[iter]]$SHAP.values_train
    shap.df_test <- Output.regressor.cts[[iter]]$SHAP.values_test
    
    labels_train <- paste(Output.regressor.cts[[iter]]$Train.labels$bin,Output.regressor.cts[[iter]]$Train.labels$Type, sep = "-")
    labels_test <- paste(Output.regressor.cts[[iter]]$Test.labels$bin,Output.regressor.cts[[iter]]$Test.labels$Type, sep = "-")
    
    shap.df_train$labels <- labels_train
    shap.df_test$labels <- labels_test
    
    #NOTA per gab del futuro: shap.df è semplicemente rbind senza sortare perchè poi la media viene fatta con group by
    #NOTA per fabio del futuro: così è molto meglio perchè tutte le misure hanno 10 replicates da fare average
    
    shap.df <- rbind(shap.df_train, shap.df_test)
    
    melted.shap.df <- melt(shap.df)
    melted.shap.df$iter <- iter
    
    avg.shap.df <- rbind(avg.shap.df,melted.shap.df)
    
    # SHAP values
    shap.value <- Output.regressor.cts[[iter]]$SHAP.imp
    shap.value$iter <- iter
    avg.shap.values <- rbind(avg.shap.values,shap.value)
    
  }
  
  # SHAP df
  avg.shap.df <- as.data.frame(avg.shap.df)
  avg.shap.df <- avg.shap.df %>% group_by(labels,variable) %>% summarise("value" = mean(value))
  avg.shap.df <- dcast(avg.shap.df, labels ~ variable, value.var="value")
  models.shap.df[[paste(class,model,cts_label,sep = "::")]] <- avg.shap.df
  
  # SHAP values
  avg.shap.values <- as.data.frame(avg.shap.values)
  avg.shap.values <- avg.shap.values %>% group_by(Feature) %>% summarise("Mean_Abs_SHAP" = mean(Mean_Abs_SHAP))
  models.shap.values[[paste(class,model,cts_label,sep = "::")]] <- avg.shap.values
  
  # major edit here, shifted this inside loop and sourcing next code
  
  save(models.shap.df,
       models.shap.values,
       models.X.train,
       models.X.test,
       file = paste0("/home/ieo7429/Scrivania/results_regressor_gab/SHAP/",class,"::",model,"::",cts_label,"_Avg_shap_values_InteractomeINSIDER.RData") # remember here it's related to last trained model. Improve
  )
  
}
