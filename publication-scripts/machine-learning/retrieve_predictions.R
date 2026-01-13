rm(list=ls())
gc(full=T)

wd <- 'path/to/GitHub/copy-number-annotation/'

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

model.outputs <- list.files(path = paste0(wd, "data/InteractomeINSIDER/"), full.names = TRUE)
patterns <- c('Arm-level','Chromosome-level','Mid-length','Small-scale','NoCluster')

for(i in patterns){
  
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
  
  pred_ampl_path <- paste0(wd, "data/InteractomeINSIDER/",i,"-pred_ampl.rds")
  pred_del_path <- paste0(wd, "data/InteractomeINSIDER/",i,"-pred_del.rds")
    
  saveRDS(object = mean_preds_ampl, file = pred_ampl_path)
  saveRDS(object = mean_preds_del, file = pred_del_path)

}
