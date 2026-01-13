X_data <- ml_table[,!colnames(ml_table) %in% c('bin','Type','ampl_score', 'del_score')]
features <- colnames(X_data)
y_data <- ml_table[,label]

set.seed(122)
split_indices <- createDataPartition(y_data, p = 0.8, list = FALSE)
X_train <- X_data[split_indices, ]
X_test <- X_data[-split_indices, ]
y_train <- y_data[split_indices]
y_test <- y_data[-split_indices]

param_grid <- expand.grid(nrounds = c(300),
                          max_depth = c(3, 4, 5),
                          eta = c(0.025, 0.05, 0.1, 0.15, 0.25),
                          gamma = c(0.2, 0.3, 0.4),
                          colsample_bytree = c(1),
                          min_child_weight = c(1),
                          subsample = c(0.8))

correlation_summary <- function(data, lev = NULL, model = NULL) {
  cor_coef <- cor(data$obs, data$pred)
  res <- if (is.numeric(cor_coef)) cor_coef else numeric(1)
  names(res) <- "correlation"
  res
}

ctrl <- trainControl(method = "cv",
                     number = 5,
                     verboseIter = TRUE,
                     summaryFunction = defaultSummary)

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

print(xgb_tune$bestTune)

hyperTuning_path <- paste0(wd, 'data/hyperparam_tuning/', classS, '_hyper_tuning.rds')

saveRDS(xgb_tune, file = hyperTuning_path)

rm(ctrl, correlation_summary, param_grid, xgb_tune, X_data,y_data,X_train,X_test,y_train,y_test)
