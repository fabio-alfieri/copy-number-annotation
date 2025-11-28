
chr_nums <- as.factor(unlist(lapply(X = ml_table$bin, FUN = function(x){
  chr_number <- strsplit(x, split = "_")[[1]][[1]]
})))

chromosomes <- sort(as.numeric(levels(chr_nums)))

for(chromosome in chromosomes){
  
source('/home/ieo7429/Scrivania/scripts/1_hyperTuning_LOCO_gab.R', local = T)

hyper_tuning <- readRDS(file = paste0('/home/ieo7429/Scrivania/hyperparam_tuning/', 
                                      classS, 
                                      '_hyper_tuning_LOCO_', chromosome,'.rds'))


params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
params <- as.list(params[1,-c(7:13)])

source('/home/ieo7429/Scrivania/scripts/2_runRegressor_LOCO_gab.R', local = T)

}