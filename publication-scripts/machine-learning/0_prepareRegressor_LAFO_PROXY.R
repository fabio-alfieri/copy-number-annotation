hyperTuning_script <- './1_hyperTuning_LAFO_PROXY.R'
runRegressor_script <- './2_runRegressor_LAFO_PROXY.R'

features_to_keep <- c("Centromere_Length", 
                      "Chromosome_Length", 
                      "Centromere_Type", 
                      "distance.to.centromere", 
                      "distance.to.telomere")

for(feature_to_keep in features_to_keep){
  
  hyperTuning_path <- paste0('../Data/hyperparam_tuning/', classS, '_hyper_tuning_LAFO_PROXY_', feature_to_keep,'.rds')
  
  features_to_discard <- features_to_keep[features_to_keep != feature_to_keep]
  
  source(hyperTuning_script, local = T)
  
  hyper_tuning <- readRDS(file = hyperTuning_path)
  
  
  params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
  params <- as.list(params[1,-c(7:13)])
  
  source(runRegressor_script, local = T)
  
}