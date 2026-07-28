rm(list=ls())
gc(full=T)

packages <- c(
  "stringr", "tidyr", "ggplot2", "fastcluster", "tidyverse"
)

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

model.outputs_ext <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/chromosome_specific/", full.names = TRUE)

chr_nums <- 1:22
patterns <- c('Arm-level','Chromosome-level','Mid-length','Small-scale','NoCluster')

for(chr_num in chr_nums){

  model.outputs <- model.outputs_ext[grepl(pattern = paste0("CS", "_", chr_num, "_"), x = model.outputs_ext)]
  
  for(i in patterns){
    
    if((i == "Arm-level") & (chr_num %in% c(13,14,15,21,22))){
      next  
    }
    
    load(model.outputs[grep(i, model.outputs)][grep('ampl',model.outputs[grep(i, model.outputs)])])
    model_ampl <- Output.regressor; model_ampl <- model_ampl$ampl_score 
    
    load(model.outputs[grep(i, model.outputs)][grep('del',model.outputs[grep(i, model.outputs)])])
    model_del <- Output.regressor; model_del <- model_del$del_score
    
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
    
    pred_ampl_path <- paste0("/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/predictions_chromosome_specific/", i, "_", chr_num, "-pred_ampl.rds")
    pred_del_path <- paste0("/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/predictions_chromosome_specific/", i, "_", chr_num, "-pred_del.rds")
    
    saveRDS(object = mean_preds_ampl, file = pred_ampl_path)
    saveRDS(object = mean_preds_del, file = pred_del_path)
    
    }
}
