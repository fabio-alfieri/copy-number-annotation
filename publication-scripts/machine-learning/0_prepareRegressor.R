wd <- 'path/to/GitHub/copy-number-annotation/'

hyperTuning_script <- 'publication-scripts/machine-learning/1_hyperTuning.R'
runRegressor_script <- 'publication-scripts/machine-learning/2_runRegressor.R'
hyperTuning_path <- paste0(wd, 'data/hyperparam_tuning/', classS, '_hyper_tuning.rds')

source(paste0(wd,hypertuning_script), local = T)

hyper_tuning <- readRDS(file = hyperTuning_path)

params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
params <- as.list(params[1,-c(7:13)])

source(paste0(wd,runRegressor_script), local = T)
