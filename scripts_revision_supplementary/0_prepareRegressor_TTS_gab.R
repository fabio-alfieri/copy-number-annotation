cts_all <- ml_table$Type
cts <- sort(unique(ml_table$Type))
cts_pairs <- expand.grid(cts, cts)

for(idx in 1:nrow(cts_pairs)){
  
  cts_pair <- as.character(unlist(cts_pairs[idx,]))
  
  cts_label <- paste0(cts_pair[1], "_", cts_pair[2])
  
  print(cts_label)
  
  source('/home/ieo7429/Scrivania/scripts/refactored/1_hyperTuning_TTS_gab.R', local = T)
  
  hyper_tuning <- readRDS(file = paste0('/home/ieo7429/Scrivania/hyperparam_tuning/', classS, '_hyper_tuning_TTS_', cts_label,'.rds'))
  
  params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
  params <- as.list(params[1,-c(7:13)])
  
  source('/home/ieo7429/Scrivania/scripts/refactored/2_runRegressor_TTS_gab.R', local = T)
  
}