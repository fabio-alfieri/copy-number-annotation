hyperTuning_script <- './1_hyperTuning_PRO.R'
runRegressor_script <- './2_runRegressor_PRO.R'

train_lengths <- c(1150,1200,1250,1300,1350,1400)
test_lengths <- c(460,480,500,520,540,560)

for (idx in seq(train_lengths)) {

ml_table <- ml_table_backup
train_length <- train_lengths[idx]
test_length <- test_lengths[idx]
  
source(hyperTuning_script, local = T)

hyperTuning_path <- paste0('../Data/hyperparam_tuning/', classS,'_hyper_tuning_PRO_', train_length, "_", test_length, "_", num_of_skipped_chrs, '.rds')

hyper_tuning <- readRDS(file = hyperTuning_path)

params <- hyper_tuning$results[order(hyper_tuning$results$RMSE, decreasing = F),]
params <- as.list(params[1,-c(7:13)])

source(runRegressor_script, local = T)

}