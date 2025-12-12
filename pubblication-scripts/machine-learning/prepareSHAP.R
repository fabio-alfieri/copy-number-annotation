rm(list = ls())
gc(full=T)

packages <- c("reshape2", "tidyverse")

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

interactome_insider_path <- "../Data/InteractomeINSIDER/"
model.outputs <- list.files(path = interactome_insider_path, full.names = TRUE)

for(file in model.outputs){
  
  models.shap.df <- list()
  models.shap.values <- list()
  
  models.X.test <- list()
  models.X.train <- list()
  
  name <- strsplit(gsub("_0\\.1Mbp.*$", "", basename(file)),"_")
  print(name)
  model <- name[[1]][length(name[[1]])]
  class <- gsub("Output.regressor_Length_|Output.regressor_Length_Location_","",gsub("_0\\.1Mbp.*$", "", basename(file)))
  class <- gsub(paste0("_",model),"",class)
  if(model == "ampl"){
    model <- "Amplification model"
  }else{
    model <- "Deletion model"
  }
  
  load(file)
  
  X_test <- Output.regressor[[1]][[1]]$X_test
  X_train <- Output.regressor[[1]][[1]]$X_train
  
  models.X.test[[paste(class,model,sep = "::")]] <- X_test
  models.X.train[[paste(class,model,sep = "::")]] <- X_train
  
  avg.shap.df <- c()
  avg.shap.values <- c()
  
  for(iter in 1:10){
    
    # SHAP df
    
    shap.df_train <- Output.regressor[[1]][[iter]]$SHAP.values_train
    shap.df_test <- Output.regressor[[1]][[iter]]$SHAP.values_test
    
    labels_train <- paste(Output.regressor[[1]][[iter]]$Train.labels$bin,Output.regressor[[1]][[iter]]$Train.labels$Type, sep = "-")
    labels_test <- paste(Output.regressor[[1]][[iter]]$Test.labels$bin,Output.regressor[[1]][[iter]]$Test.labels$Type, sep = "-")
    
    shap.df_train$labels <- labels_train
    shap.df_test$labels <- labels_test
    
    shap.df <- rbind(shap.df_train, shap.df_test)
    
    melted.shap.df <- melt(shap.df)
    melted.shap.df$iter <- iter
    
    avg.shap.df <- rbind(avg.shap.df,melted.shap.df)
    
    # SHAP values
    shap.value <- Output.regressor[[1]][[iter]]$SHAP.imp
    shap.value$iter <- iter
    avg.shap.values <- rbind(avg.shap.values,shap.value)
    
  }
  
  # SHAP df
  avg.shap.df <- as.data.frame(avg.shap.df)
  avg.shap.df <- avg.shap.df %>% group_by(labels,variable) %>% summarise("value" = mean(value))
  avg.shap.df <- dcast(avg.shap.df, labels ~ variable, value.var="value")
  models.shap.df[[paste(class,model,sep = "::")]] <- avg.shap.df
  
  # SHAP values
  avg.shap.values <- as.data.frame(avg.shap.values)
  avg.shap.values <- avg.shap.values %>% group_by(Feature) %>% summarise("Mean_Abs_SHAP" = mean(Mean_Abs_SHAP))
  models.shap.values[[paste(class,model,sep = "::")]] <- avg.shap.values
  
  models_shap_path <- paste0("../Data/SHAP/",class,"::",model,"_Avg_shap_values_InteractomeINSIDER.RData")
  
  save(models.shap.df,
        models.shap.values,
        models.X.train,
        models.X.test,
        file = 
  )
  
}
