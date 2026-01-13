# 1. Initialization
nrounds <- 10 # Number of iterations for model averaging
cores <- 15
max_nrounds <- 300

X_data <- ml_table[,!colnames(ml_table) %in% c('bin','Type','ampl_score', 'del_score')]
features <- colnames(X_data)


# 80/20 train and test set (Use the seed in hypertuning to obtain the same train and test) 
set.seed(122)
train <- createDataPartition(y = ml_table[,label], 
                             times = 1, # number of partitions to create
                             p = 0.8, # percentage of the data goes to training
                             list = FALSE)

dftrain <- ml_table[train,]
dftest <- ml_table[-train,]

# 2. Convert data to XGBoost format
xgb_train_all = xgb.DMatrix(data = data.matrix(dftrain[,features]),
                            label = dftrain[,label],nthread = 6)
xgb_test_all = xgb.DMatrix(data = data.matrix(dftest[,features]),
                           label = dftest[,label],nthread = 6)

# 3. Model averaging loop (10 iterations, for each split the train set into train and test, evaluate the performance)
Output.regressor <- list()

for(j in 1:nrounds){
  
  set.seed(sample(1:1000, size = 1)) # set random seed to ensure reproducibility of random number generations
  
  train.train <- sample(1:nrow(dftrain), size = round(0.8*nrow(dftrain)))
  
  dftrain.dftrain <- dftrain[train.train,] # train
  dftest.dftrain <- dftrain[-train.train,] # validation
  
  xgb_train = xgb.DMatrix(data = data.matrix(dftrain.dftrain[,features]), 
                          label = dftrain.dftrain[,label],nthread = 6)
  
  xgb_test = xgb.DMatrix(data = data.matrix(dftest.dftrain[,features]), 
                         label = dftest.dftrain[,label],nthread = 6)
  
  # Create a watchlist to monitor performance on both the training and test sets
  watchlist <- list(train = xgb_train, test = xgb_test)
  
  # Train the XGBoost model with watchlist
  xgb_model <- xgb.train(data = xgb_train, params = params, 
                         watchlist = watchlist,
                         nrounds = max_nrounds, nthread = cores, verbose = 0)
  
  # Extract evaluation metrics from the watchlist
  evaluation <- xgb_model$evaluation_log
  best_iter <- min(evaluation[evaluation$test_rmse == min(evaluation$test_rmse),]$iter)
  
  # Model after watching, with best iteration
  xgb_model <- xgb.train(data = xgb_train, params = params,
                         watchlist = watchlist, verbose = 0,
                         nrounds = best_iter, nthread = cores)
  
  # Feature importance
  imp <- xgb.importance(feature_names = colnames(xgb_train),
                        model = xgb_model)
  
  # The original test data
  pred_train <- predict(xgb_model, xgb_train_all)
  pred_test <- predict(xgb_model, xgb_test_all)
  
  y_train <- dftrain[,label]
  y_test <- dftest[,label]
  
  save.data_train <- dftrain[,c("bin","Type",label)]
  save.data_train$prediction <- pred_train
  
  save.data_test <- dftest[,c("bin","Type",label)]
  save.data_test$prediction <- pred_test
  
  spearman_corr <- cor(save.data_test[,label], save.data_test$prediction, method = "spearman")
  pearson_corr <- cor(save.data_test[,label], save.data_test$prediction, method = "pearson")
  
  r2 <- caret::R2(save.data_test$prediction, save.data_test[,label])
  rmse <- caret::RMSE(save.data_test$prediction, save.data_test[,label])
  
  # Convert the data to a DMatrix for SHAP calculation
  dmat_train <- xgb.DMatrix(data = as.matrix(dftrain[,features]))
  dmat_test <- xgb.DMatrix(data = as.matrix(dftest[,features]))
  
  # Predict with SHAP contributions
  shap_values_train <- predict(xgb_model, dmat_train, predcontrib = TRUE)
  shap_values_test <- predict(xgb_model, dmat_test, predcontrib = TRUE)
  
  # Convert SHAP values to a data frame for easier handling
  shap_df_train <- as.data.frame(shap_values_train)
  shap_df_test <- as.data.frame(shap_values_test)
  
  # Calculate the mean absolute SHAP value for each feature to determine importance
  mean_abs_shap <- colMeans(abs(shap_df_test))
  
  # Create a data frame for feature importance based on SHAP
  shap_importance <- data.frame(
    Feature = colnames(shap_df_test),
    Mean_Abs_SHAP = mean_abs_shap
  )
  
  
  Output.regressor.sub <- list("Model" = xgb_model,
                               "Performance" = list("Spearman" = spearman_corr, 
                                                    "Pearson" = pearson_corr,
                                                    "R2" = r2,
                                                    "RMSE" = rmse),
                               "X_train" = dftrain,
                               "X_test" = dftest,
                               "y_train" = y_train,
                               "y_test" = y_test, 
                               "Train pred" = pred_train,
                               "Test_pred" = pred_test,
                               "Train.labels" = save.data_train,
                               "Test.labels" = save.data_test,
                               "SHAP.values_train" = shap_df_train,
                               "SHAP.values_test" = shap_df_test,
                               "Feature.imp" = imp,
                               "SHAP.imp" = shap_importance)
  
  Output.regressor[[label]][[j]] <- Output.regressor.sub
}

output_regressor_path <- paste0("../Data/results_regressor/",
                                folder,"/",
                                "Output.regressor_",
                                clus.group,"_",
                                classS,"_",
                                ampl_or_del,"_",
                                level,"_",
                                bincov,"_",
                                Sys.time(),".RData")
  
save(Output.regressor, 
     file = output_regressor_path)


