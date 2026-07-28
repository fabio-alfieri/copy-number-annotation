chr_nums <- as.factor(unlist(lapply(X = ml_table$bin, FUN = function(x){chr_number <- strsplit(x, split = "_")[[1]][[1]]})))
chromosomes <- sort(as.numeric(levels(chr_nums)))

chromosomes_pairs <- expand.grid(chromosomes, chromosomes)
# 4_9

for(idx in 1:nrow(chromosomes_pairs)){
  
  chromosomes_pair <- as.character(unlist(chromosomes_pairs[idx,]))
  
  chromosomes_label <- paste0(chromosomes_pair[1], "_", chromosomes_pair[2])
  
  print(chromosomes_label)
  
  source('/home/ieo7429/Scrivania/scripts/refactored/1_hyperTuning_CSP_gab.R', local = T)
  
  hyper_tuning <- readRDS(file = paste0('/home/ieo7429/Scrivania/hyperparam_tuning/', classS, '_hyper_tuning_CSP_scaled_', chromosomes_label,'.rds'))
  
  params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
  params <- as.list(params[1,-c(7:13)])
  
  source('/home/ieo7429/Scrivania/scripts/refactored/2_runRegressor_CSP_gab.R', local = T)
  
}