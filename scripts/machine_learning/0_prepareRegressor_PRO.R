# Run hyperparameter tuning

train_lengths <- c(1150,1200,1250,1300,1350,1400)
test_lengths <- c(460,480,500,520,540,560)

for (idx in seq(train_lengths)) {

ml_table <- ml_table_backup
train_length <- train_lengths[idx]
test_length <- test_lengths[idx]
  
source('/home/ieo7429/Scrivania/scripts/1_hyperTuning_PRO_gab.R', local = T)

hyper_tuning <- readRDS(file = paste0('/home/ieo7429/Scrivania/hyperparam_tuning/', 
                                      classS,'_hyper_tuning_PRO_', train_length, "_", test_length, "_", num_of_skipped_chrs, '.rds'))

# Parameters for running regression model
#params <- hyper_tuning$results[order(hyper_tuning$results$correlation, decreasing = T),] #### MAJOR EDIT GAB, since caret changed, now we sort by RMSE increasing (decreasing = F)
#params <- as.list(params[1,-c(7:13)])

params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
params <- as.list(params[1,-c(7:13)])

# Run regression model
source('/home/ieo7429/Scrivania/scripts/2_runRegressor_PRO_gab.R', local = T)

}