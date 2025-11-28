rm(list = ls())
gc(full=T)

library(reshape2)
library(tidyverse)

#Interactome INSIDER (alternative to HIPPIE)
model.outputs <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/", full.names = TRUE)
model.outputs <- model.outputs[grep('Output.regressor',model.outputs)]
model.outputs <- model.outputs[grep('2025-08',model.outputs)]

#model.outputs <- model.outputs[-grep('Centromere',model.outputs)]
#model.outputs <- model.outputs[-grep('Telomere',model.outputs)]
#model.outputs <- model.outputs[-grep('Interstitial',model.outputs)]

for(file in model.outputs){
  
  models.Gain <- list()
  
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
  
  model.Gain <- data.frame()
  
  for(iter in 1:10){
    
    # SHAP df
    
    tmp <- Output.regressor[[1]][[iter]]$Feature.imp
    tmp$iter <- iter
    model.Gain <- rbind(model.Gain, tmp)
    
  }
  
  # average Gain
  models.Gain[[paste(class,model,sep = "::")]] <- model.Gain %>% group_by(Feature) %>% summarise(mean(Gain, na.rm = TRUE))
  
  save(models.Gain,
       file = paste0("/home/ieo7429/Scrivania/results_regressor_gab/SHAP/",class,"::",model,"_Avg_Gain.RData") # remember here it's related to last trained model. Improve
  )
  
  
}


gain.files <- grep('_Avg_Gain', list.files('/home/ieo7429/Scrivania/results_regressor_gab/SHAP/', full.names = T), value = T)

gains <- data.frame()
for(file in gain.files){
  load(file)
  df <- models.Gain[[1]]
  tmp <- do.call(rbind, str_split(names(models.Gain), pattern = '::'))
  colnames(tmp) <- c('model','scna')
  df <- cbind(df, tmp)
  gains <- rbind(gains, df)
}
colnames(gains)[2] <- 'Gain'

gains %>% ggplot(aes(x = Gain)) +
  geom_histogram(binwidth = .01) +
  facet_wrap(~scna)

gains.total <- gains %>% group_by(scna, Feature) %>% summarize(no_zero = sum(Gain > 0))
important_features <- levels(factor(gains.total[gains.total$no_zero >= 4,]$Feature))

gains %>% filter(Feature %in% important_features) %>%
  ggplot(aes(x = Gain, y = Feature)) +
  geom_col(aes(fill = model)) +
  facet_wrap(~scna)

small_features <- c(grep('Length_',important_features, value = T), 'genes.bin', 'Ess.distance_pancancer', 'mutations_norm', 'partners.trans')

pdf('dev/imgs/feature_Gain_big_features.pdf', width = 10, height = 6)
gains %>% filter(Feature %in% important_features[!important_features %in% small_features]) %>%
  ggplot(aes(x = Gain, y = Feature)) +
  geom_col(aes(fill = model)) +
  facet_wrap(~scna)
dev.off()

pdf('dev/imgs/feature_Gain_small_features.pdf', width = 10, height = 6)
gains %>% filter(Feature %in% small_features) %>%
  ggplot(aes(x = Gain, y = Feature)) +
  geom_col(aes(fill = model)) +
  facet_wrap(~scna)
dev.off()

write.table("dev/Data/gain_table.tsv", x = gains, sep = "\t", quote = F, row.names = T, col.names = T)

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
