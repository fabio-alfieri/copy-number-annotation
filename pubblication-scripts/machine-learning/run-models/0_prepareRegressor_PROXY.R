# Run hyperparameter tuning 
source('/home/ieo7429/Scrivania/scripts/1_hyperTuning_PROXY_gab.R', local = T)

hyper_tuning <- readRDS(file = paste0('/home/ieo7429/Scrivania/hyperparam_tuning/', 
                                      classS, 
                                      '_hyper_tuning_PROXY.rds'))

# Parameters for running regression model
#params <- hyper_tuning$results[order(hyper_tuning$results$correlation, decreasing = T),] #### MAJOR EDIT GAB, since caret changed, now we sort by RMSE increasing (decreasing = F)
#params <- as.list(params[1,-c(7:13)])

params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
params <- as.list(params[1,-c(7:13)])

# Run regression model
source('/home/ieo7429/Scrivania/scripts/2_runRegressor_PROXY_gab.R', local = T)
