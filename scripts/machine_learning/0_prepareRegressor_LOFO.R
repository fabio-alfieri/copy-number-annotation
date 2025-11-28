
features_to_exclude <- c("Centromere_Length", "Chromosome_Length", "Centromere_Type", "distance.to.centromere", "distance.to.telomere")

for(feature_to_exclude in features_to_exclude){
  
  source('/home/ieo7429/Scrivania/scripts/1_hyperTuning_LOFO_gab.R', local = T)
  
  hyper_tuning <- readRDS(file = paste0('/home/ieo7429/Scrivania/hyperparam_tuning/', 
                                        classS, 
                                        '_hyper_tuning_LOFO_', feature_to_exclude,'.rds'))
  
  
  params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
  params <- as.list(params[1,-c(7:13)])
  
  source('/home/ieo7429/Scrivania/scripts/2_runRegressor_LOFO_gab.R', local = T)
  
}