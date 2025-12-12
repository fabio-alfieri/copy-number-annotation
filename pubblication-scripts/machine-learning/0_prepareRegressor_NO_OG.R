hyperTuning_script <- './1_hyperTuning_NO_OG.R'
runRegressor_script <- './2_runRegressor_NO_OG.R'

source(hyperTuning_script, local = T)

hyperTuning_path <- paste0('../Data/hyperparam_tuning/', classS, '_hyper_tuning_NO_OG.rds')

hyper_tuning <- readRDS(file = hyperTuning_path)

params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
params <- as.list(params[1,-c(7:13)])

source(runRegressor_script, local = T)
