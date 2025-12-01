
cts_all <- ml_table$Type

cts <- sort(unique(ml_table$Type))

for(ct in cts){
  print(ct)
  source('/home/ieo7429/Scrivania/scripts/1_hyperTuning_LOCTO_gab.R', local = T)
  
  hyper_tuning <- readRDS(file = paste0('/home/ieo7429/Scrivania/hyperparam_tuning/', 
                                        classS, 
                                        '_hyper_tuning_LOCTO_', ct,'.rds'))
  
  
  params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
  params <- as.list(params[1,-c(7:13)])
  
  source('/home/ieo7429/Scrivania/scripts/2_runRegressor_LOCTO_gab.R', local = T)
  
}