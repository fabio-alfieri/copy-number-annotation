wd <- 'path/to/GitHub/copy-number-annotation/'

# Extract features
X_data <- ml_table[,!colnames(ml_table) %in% c('bin','Type','ampl_score', 'del_score', feature_to_exclude)]
features <- colnames(X_data)
y_data <- ml_table[,label]

# 80/20 split
set.seed(122)
split_indices <- createDataPartition(y_data, p = 0.8, list = FALSE)
X_train <- X_data[split_indices, ]
X_test <- X_data[-split_indices, ]
y_train <- y_data[split_indices]
y_test <- y_data[-split_indices]


# Define the parameter grid
param_grid <- expand.grid(nrounds = c(300),
                          max_depth = c(3, 4, 5),
                          eta = c(0.025, 0.05, 0.1, 0.15, 0.25),
                          gamma = c(0.2, 0.3, 0.4),
                          colsample_bytree = c(1),
                          min_child_weight = c(1),
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
  metric = "RMSE",
  maximize = FALSE,
  nthread = 15
)

# Print the best parameters
print(xgb_tune$bestTune)

hyperTuning_path <- paste0(wd, 'data/hyperparam_tuning/', classS, '_hyper_tuning_LOFO_', feature_to_exclude,'.rds')

saveRDS(xgb_tune, file = hyperTuning_path)

rm(ctrl, correlation_summary, param_grid, xgb_tune, X_data,y_data,X_train,X_test,y_train,y_test)
