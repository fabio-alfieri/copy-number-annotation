wd <- 'path/to/GitHub/copy-number-annotation/'

hyperTuning_script <- 'publication-scripts/machine-learning/1_hyperTuning_OO.R'
runRegressor_script <- 'publication-scripts/machine-learning/2_runRegressor_OO.R'

source(paste0(wd, hyperTuning_script), local = T)

hyperTuning_path <- paste0(wd, 'data/hyperparam_tuning/', classS, '_hyper_tuning_OO.rds')

hyper_tuning <- readRDS(file = hyperTuning_path)

params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
params <- as.list(params[1,-c(7:13)])

source(paste0(wd, runRegressor_script), local = T)
