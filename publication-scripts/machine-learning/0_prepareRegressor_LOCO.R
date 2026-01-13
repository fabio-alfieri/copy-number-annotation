hyperTuning_script <- './1_hyperTuning_LOCO.R'
runRegressor_script <- './2_runRegressor_LOCO.R'

chr_nums <- as.factor(unlist(lapply(X = ml_table$bin, FUN = function(x){
  chr_number <- strsplit(x, split = "_")[[1]][[1]]
})))

chromosomes <- sort(as.numeric(levels(chr_nums)))

for(chromosome in chromosomes){
  
hyperTuning_path <- paste0('../Data/hyperparam_tuning/', classS, '_hyper_tuning_LOCO_', chromosome,'.rds')

source(hyperTuning_script, local = T)

hyper_tuning <- readRDS(file = hyperTuning_path)

params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
params <- as.list(params[1,-c(7:13)])

source(runRegressor_script, local = T)

}