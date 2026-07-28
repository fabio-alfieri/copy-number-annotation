X_data <- ml_table[,!colnames(ml_table) %in% c('bin','Type','ampl_score', 'del_score')]
features <- colnames(X_data)
y_data <- ml_table[,label]

if(length(unique(cts_pair)) > 1){
  ct_idxs_train <- which(cts_all == cts_pair[1])
  ct_idxs_test <- which(cts_all == cts_pair[2])
  print(unique(ml_table[ct_idxs_train,]$Type))
  print(unique(ml_table[ct_idxs_test,]$Type))
} else if (length(unique(cts_pair)) == 1){
  ct_idxs <- which(cts_all == cts_pair[1])
  X_data <- X_data[ct_idxs,]
  y_data <- y_data[ct_idxs]
  set.seed(122)
  ct_idxs_train <- createDataPartition(y_data, p = 0.8, list = FALSE)
  ct_idxs_test <- -ct_idxs_train
  print(unique(ml_table[ct_idxs,][ct_idxs_train,]$Type))
  print(unique(ml_table[ct_idxs,][ct_idxs_test,]$Type))
}

# LOCO split

X_train <- X_data[ct_idxs_train, ]
X_test <- X_data[ct_idxs_test, ]
y_train <- y_data[ct_idxs_train]
y_test <- y_data[ct_idxs_test]

# Define the parameter grid
param_grid <- expand.grid(nrounds = c(300),
                          max_depth = c(3, 4, 5),
                          eta = c(0.025, 0.05, 0.1, 0.15, 0.25),
                          gamma = c(0.2, 0.3, 0.4),
                          colsample_bytree = c(1),
                          min_child_weight = c(1), # Additional XGBoost parameter
                          subsample = c(0.8))

# Define cross-validation: To summarize performance during cross-validation
correlation_summary <- function(data, lev = NULL, model = NULL) {
  cor_coef <- cor(data$obs, data$pred)
  res <- if (is.numeric(cor_coef)) cor_coef else numeric(1)
  names(res) <- "correlation"
  res
}

# Set up Cross-validation
ctrl <- trainControl(method = "cv", # Use k-fold CV
                     number = 5, # Perform 5 folds
                     verboseIter = TRUE, # Print progress during training
                     summaryFunction = defaultSummary) ### editing here, instead of correlation_summary, using defaultSummary

# Hyperparameter tuning
xgb_tune <- caret::train(
  x = X_train,
  y = y_train, 
  method = "xgbTree", 
  trControl = ctrl,
  tuneGrid = param_grid,
  metric = "RMSE", # Use RMSE as the optimization metric
  maximize = FALSE, # Since we want to minimize RMSE
  nthread = 15
)

# Print the best parameters
print(xgb_tune$bestTune)

saveRDS(xgb_tune, file = paste0('/home/ieo7429/Scrivania/hyperparam_tuning/', classS, '_hyper_tuning_TTS_', cts_label, '.rds')) # remember here it's related to last trained model. Improve

rm(ctrl, correlation_summary, param_grid, xgb_tune, X_data,y_data,X_train,X_test,y_train,y_test)
