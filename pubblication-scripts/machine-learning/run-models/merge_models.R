# merge objects
rm(list=ls())
gc(full=T)

library(ggplot2)
library(tidyr)
library(stringr)
library(tidyverse)
library(fastcluster)

setwd("/home/ieo7429/Scrivania/")

model.outputs <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/", full.names = TRUE)
model.outputs <- grep("2025-08", model.outputs, value = TRUE)
patterns <- c('Arm-level','Chromosome-level','Mid-length','Small-scale','no_cluster')

shap.outputs <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/SHAP/", full.names = TRUE)
shap.outputs <- grep("::", shap.outputs, value = TRUE)

for(i in patterns){
  print(i)
  
  # Feature Matrix Values
  load(model.outputs[grep(i, model.outputs)][grep('ampl',model.outputs[grep(i, model.outputs)])])
  model_ampl <- Output.regressor; model_ampl <- model_ampl$ampl_score 
  
  load(model.outputs[grep(i, model.outputs)][grep('del',model.outputs[grep(i, model.outputs)])])
  model_del <- Output.regressor; model_del <- model_del$del_score
  
  rm(Output.regressor)
  
  models.X <- list(ampl = rbind(model_ampl[[1]]$X_train, model_ampl[[2]]$X_test),
                   del = rbind(model_del[[1]]$X_train, model_del[[2]]$X_test))
  
  names(models.X) <- c(paste0(i,'::Amplification model'),
                       paste0(i,'::Deletion model'))
  
  # SHAP values
  load(shap.outputs[grep(paste0(i,'::'), shap.outputs)][1])
  shap_ampl <- models.shap.df[[1]] %>% dplyr::select(-BIAS)
  rm(models.shap.df)
  
  load(shap.outputs[grep(paste0(i,'::'), shap.outputs)][2])
  shap_del <- models.shap.df[[1]] %>% dplyr::select(-BIAS)
  rm(models.shap.df, models.X.test, models.X.train, models.shap.values)
  
  models.shap.df <- list(ampl = shap_ampl,
                         del = shap_del)
  
  names(models.shap.df) <- c(paste0(i,'::Amplification model'),
                             paste0(i,'::Deletion model'))
  
  # aggregate and save
  feature_and_shap <- list(models.X,models.shap.df)
  names(feature_and_shap) <- c('models.X','models.shap.df')
  
  write_rds(feature_and_shap, file = paste0("dev/Data/SHAP_and_FeatureMatrix_",i,"_AmplDel.rds"))
  
}
