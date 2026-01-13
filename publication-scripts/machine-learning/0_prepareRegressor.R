hyperTuning_script <- './1_hyperTuning.R'
runRegressor_script <- './2_runRegressor.R'
hyperTuning_path <- paste0('../Data/hyperparam_tuning/', classS, '_hyper_tuning.rds')

source(hypertuning_script, local = T)

hyper_tuning <- readRDS(file = hyperTuning_path)

params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
params <- as.list(params[1,-c(7:13)])

source(runRegressor_script, local = T)
