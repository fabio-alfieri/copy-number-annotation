wd <- 'path/to/GitHub/copy-number-annotation/'

extract_indexes <- function(ml_table, train_length, test_length){
  
  bool_vec <- c(rep(T, train_length), rep(F, test_length)) # 50 20
  
  df <- as.data.frame(do.call(rbind, 
                              lapply(X = ml_table$bin, 
                                     FUN = function(x){strsplit(x, split = "_")[[1]]
                                     })))
  print(nrow(df))
  colnames(df) <- c("chr", "bin")
  df$chr <- as.numeric(df$chr)
  df$bin <- as.numeric(df$bin)
  df$Type <- ml_table$Type
  
  coord_sorted_by_tt <- unique(df[order(df$Type, as.numeric(df$chr), as.numeric(df$bin), decreasing = F),])
  grouped <- coord_sorted_by_tt %>% group_by(Type, chr) %>% summarise(numBin = n())
  
  df <- df[order(df$Type, as.numeric(df$chr), as.numeric(df$bin), decreasing = F),]
  
  df$train <- unlist(apply(X = grouped, MARGIN = 1, function(x){
    numBin = as.numeric(x["numBin"])
    TF_seq <- rep(bool_vec, times = numBin / length(bool_vec))
    final_vec <- c(TF_seq, rep(T, (numBin - length(TF_seq))))
    
    
    }
  )
)
  
  tbl <- df %>%
    mutate(train = factor(train, levels = c(FALSE, TRUE))) %>%
    group_by(chr, train) %>%
    summarise(n = n(), .groups = "drop") %>%
    tidyr::complete(chr, train, fill = list(n = 0)) %>%
    tidyr::pivot_wider(names_from = train, values_from = n)
  
  skipped_chrs <- tbl[which(tbl$`FALSE` == 0),]$chr; print(skipped_chrs)
  num_of_skipped_chrs <<- length(skipped_chrs); print(num_of_skipped_chrs)
  
  df <- df[!(df$chr %in% skipped_chrs),]
  print(nrow(df))
  
  df$test <- !(df$train)
  
  df$bin <- paste0(df$chr, "_", df$bin); df$chr <- NULL
  
  return(df)
}

train_test_indexes <- extract_indexes(ml_table, train_length, test_length)

ml_table <- merge(ml_table, train_test_indexes, by = c("bin", "Type"))

X_data <- ml_table[,!colnames(ml_table) %in% c('bin','Type','ampl_score', 'del_score')]
X_train <- X_data[X_data$train == T,!colnames(X_data) %in% c('train','test')]
X_test <- X_data[X_data$test == T,!colnames(X_data) %in% c('train','test')]

features <- colnames(X_train)[colnames(X_train) == colnames(X_test)]

y_data <- ml_table[,c("train", "test", label)]
y_train <- y_data[y_data$train == T, !colnames(y_data) %in% c('train','test')]
y_test <- y_data[y_data$test == T, !colnames(y_data) %in% c('train','test')]

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
                     summaryFunction = defaultSummary)

# Hyperparameter tuning
xgb_tune <- caret::train(
  x = X_train,
  y = y_train, 
  method = "xgbTree", 
  trControl = ctrl,
  tuneGrid = param_grid,
  metric = "RMSE",
  maximize = FALSE,
  nthread = 40
)

# Print the best parameters
print(xgb_tune$bestTune)

hyperTuning_path <- paste0(wd, 'data/hyperparam_tuning/', classS,'_hyper_tuning_PRO_', train_length, "_", test_length, "_", num_of_skipped_chrs, '.rds')

saveRDS(xgb_tune, file = hyperTuning_path)

rm(ctrl, correlation_summary, param_grid, xgb_tune, X_data,y_data,X_train,X_test,y_train,y_test)