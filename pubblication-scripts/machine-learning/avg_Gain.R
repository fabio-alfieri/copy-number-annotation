rm(list=ls())
gc(full=T)

packages <- c("reshape2", "tidyverse")

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

#Interactome INSIDER (alternative to HIPPIE)
model.outputs <- list.files(path = "../Data/InteractomeINSIDER/", full.names = TRUE)
model.outputs <- model.outputs[grep('Output.regressor',model.outputs)]

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
       file = paste0("../Data/SHAP/",class,"::",model,"_Avg_Gain.RData") # remember here it's related to last trained model. Improve
  )
  
}


gain.files <- grep('_Avg_Gain', list.files('../Data/SHAP/', full.names = T), value = T)

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

pdf('../Data/plots/feature_Gain_big_features.pdf', width = 10, height = 6)
gains %>% filter(Feature %in% important_features[!important_features %in% small_features]) %>%
  ggplot(aes(x = Gain, y = Feature)) +
  geom_col(aes(fill = model)) +
  facet_wrap(~scna)
dev.off()

pdf('../Data/plots/feature_Gain_small_features.pdf', width = 10, height = 6)
gains %>% filter(Feature %in% small_features) %>%
  ggplot(aes(x = Gain, y = Feature)) +
  geom_col(aes(fill = model)) +
  facet_wrap(~scna)
dev.off()

write.table("../Data/gain/gain_table.tsv", x = gains, sep = "\t", quote = F, row.names = T, col.names = T)
