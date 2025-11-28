rm(list=ls())
gc(full=T)

library(ggplot2)
library(tidyr)
library(stringr)
library(tidyverse)
library(fastcluster)

setwd("/home/ieo7429/Scrivania/")

model.outputs <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/", full.names = TRUE)
model.outputs <- grep("2025-08", model.outputs, value = TRUE) # selected only models from august 2025
patterns <- c('Arm-level','Chromosome-level','Mid-length','Small-scale','NoCluster')

for(i in patterns){
  
  load(model.outputs[grep(i, model.outputs)][grep('ampl',model.outputs[grep(i, model.outputs)])])
  model_ampl <- Output.regressor; model_ampl <- model_ampl$ampl_score 
  
  load(model.outputs[grep(i, model.outputs)][grep('del',model.outputs[grep(i, model.outputs)])])
  model_del <- Output.regressor; model_del <- model_del$del_score

  # load("results_regressor_gab/InteractomeINSIDER/Output.regressor_Length_Mid-length_ampl_0.1Mbp_covThr_zero_2025-06-04 18:24:55.RData")
  # model_ampl <- Output.regressor; model_ampl <- model_ampl$ampl_score 
  # 
  # load("results_regressor_gab/InteractomeINSIDER/Output.regressor_Length_Mid-length_del_0.1Mbp_covThr_zero_2025-06-04 19:57:16.RData")
  # model_del <- Output.regressor; model_del <- model_del$del_score
  
  rm(Output.regressor)
  
  ampl_pred_list <- lapply(X = seq_along(along.with = model_ampl), FUN = function(x){
    model_output <- model_ampl[[x]]
    full_preds <- rbind(model_output$Train.labels, model_output$Test.labels)
  })
  
  del_pred_list <- lapply(X = seq_along(along.with = model_del), FUN = function(x){
    model_output <- model_del[[x]]
    full_preds <- rbind(model_output$Train.labels, model_output$Test.labels)
  })
  
  mean_preds_ampl <- Reduce(f = rbind, x = ampl_pred_list) %>%
                          group_by(bin, Type) %>%
                          summarise(
                          ampl_score = mean(ampl_score),
                          prediction = mean(prediction)) %>% 
                          ungroup()
  
  mean_preds_del <- Reduce(f = rbind, x = del_pred_list) %>%
                          group_by(bin, Type) %>%
                          summarise(
                          del_score = mean(del_score),
                          prediction = mean(prediction)) %>% 
                          ungroup()
  
  
  saveRDS(object = mean_preds_ampl, file = paste0("results_regressor_gab/InteractomeINSIDER/",i,"-pred_ampl.rds"))
  saveRDS(object = mean_preds_del, file = paste0("results_regressor_gab/InteractomeINSIDER/",i,"-pred_del.rds"))

}








