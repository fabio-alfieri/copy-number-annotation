
features_to_keep <- c("Centromere_Length", "Chromosome_Length", "Centromere_Type", "distance.to.centromere", "distance.to.telomere")

for(feature_to_keep in features_to_keep){
  
  features_to_discard <- features_to_keep[features_to_keep != feature_to_keep]
  
  source('/home/ieo7429/Scrivania/scripts/1_hyperTuning_LAFO_PROXY_gab.R', local = T)
  
  hyper_tuning <- readRDS(file = paste0('/home/ieo7429/Scrivania/hyperparam_tuning/', 
                                        classS, 
                                        '_hyper_tuning_LAFO_PROXY_', feature_to_keep,'.rds'))
  
  
  params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
  params <- as.list(params[1,-c(7:13)])
  
  source('/home/ieo7429/Scrivania/scripts/2_runRegressor_LAFO_PROXY_gab.R', local = T)
  
}