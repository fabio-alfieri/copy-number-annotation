rm(list = ls())
gc(full=T)

library(reshape2)
library(tidyverse)

#Interactome INSIDER (alternative to HIPPIE)
model.outputs <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/", full.names = TRUE)

##### GAB EDITED HERE

model.outputs <- grep("2025-08", model.outputs, value = TRUE) # selected only models from august 2025

# model.outputs <- model.outputs[-grep("NoCluster",model.outputs)]
# model.outputs <- model.outputs[c(1,2,3,4,13,14,15,16)]

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
  
  # NOTA per gab del futuro. In 2_runRegressor_gab.R line 12 we know that train and test are always the same

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
    
    #NOTA per gab del futuro: shap.df è semplicemente rbind senza sortare perchè poi la media viene fatta con group by
    #NOTA per fabio del futuro: così è molto meglio perchè tutte le misure hanno 10 replicates da fare average

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
 
  # major edit here, shifted this inside loop and sourcing next code
  
  save(models.shap.df,
        models.shap.values,
        models.X.train,
        models.X.test,
        file = paste0("/home/ieo7429/Scrivania/results_regressor_gab/SHAP/",class,"::",model,"_Avg_shap_values_InteractomeINSIDER.RData") # remember here it's related to last trained model. Improve
  )
  
  # source("/home/ieo7429/Scrivania/scripts/3_workflow_SHAP_gabR.R") 
   
}

# save(models.shap.df,
#      models.shap.values,
#      models.X.train,
#      models.X.test,
#      file = paste0("/home/ieo7429/Scrivania/results_regressor_gab/SHAP/",class, "_Avg_shap_values_InteractomeINSIDER.RData") # remember here it's related to last trained model. Improve
# )

if(F){
  #HIPPIE (alternative to Interactome INSIDER)
  model.outputs <- list.files(path = "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/XGBoost/results_regressor/HIPPIE",
                              full.names = TRUE)
  model.outputs <- model.outputs[-grep("NoCluster",model.outputs)]
  
  models.shap.df <- list()
  models.shap.values <- list()
  models.X.test <- list()
  
  for(file in model.outputs){
    name <- strsplit(gsub("_0\\.1Mbp.*$", "", basename(file)),"_")
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
    models.X.test[[paste(class,model,sep = "::")]] <- X_test
    
    avg.shap.df <- c()
    avg.shap.values <- c()
    
    for(iter in 1:10){
      # SHAP df
      shap.df <- Output.regressor[[1]][[iter]]$SHAP.values
      labels <- paste(Output.regressor[[1]][[iter]]$Test.labels$bin,
                      Output.regressor[[1]][[iter]]$Test.labels$Type, sep = "-")
      shap.df$labels <- labels
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
    avg.shap.df <- dcast(avg.shap.df, labels ~ variable, value.var = "value")
    models.shap.df[[paste(class,model,sep = "::")]] <- avg.shap.df
    
    # SHAP values
    avg.shap.values <- as.data.frame(avg.shap.values)
    avg.shap.values <- avg.shap.values %>% group_by(Feature) %>% summarise("Mean_Abs_SHAP" = mean(Mean_Abs_SHAP))
    models.shap.values[[paste(class,model,sep = "::")]] <- avg.shap.values
    
    
  }
  
  save(models.shap.df,
       models.shap.values,
       models.X.test,
       file = "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/SHAP/Avg_shap_values_HIPPIE.RData"
  )
}
