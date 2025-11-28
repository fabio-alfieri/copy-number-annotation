# Directionality
# rm(list = ls())
# gc()

# Color codes
background <- "#e1e5f2"
lines <- "#012a4a"
bars <- "#a4ac86" # Feature importance plots/bar plots

# Data - Interactome INSIDER
load(paste0("/home/ieo7429/Scrivania/results_regressor_gab/SHAP/",class,"::",model,"_Avg_shap_values_InteractomeINSIDER.RData")) # remember here it's related to last trained model. Improve
models <- names(models.shap.df)

# Outputs
direction.plots <- list()
direction.m <- c()

for(model in models){
  print(model)
  source('/home/ieo7429/Scrivania/scripts/8_SHAP_directionality_gab.R', local = T, verbose = F)
}

direction.m <- as.data.frame(direction.m)
colnames(direction.m) <- c("Model","Feature","Pearson","Pearson_p.val", "Spearman", "Spearman_p.val")
direction.m$Pearson <- as.numeric(direction.m$Pearson)
direction.m$Spearman <- as.numeric(direction.m$Spearman)
direction.m$Pearson_p.val <- as.numeric(direction.m$Pearson_p.val)
direction.m$Spearman_p.val <- as.numeric(direction.m$Spearman_p.val)

save(direction.m, 
     direction.plots, 
     file = paste0("/home/ieo7429/Scrivania/results_regressor_gab/SHAP/", class, "_SHAP_directions_InteractomeINSIDER.RData")) # remember here it's related to last trained model. Improve

if(F){
  # Data - HIPPIE
  load( "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/SHAP/Avg_shap_values_HIPPIE.RData")
  models <- names(models.shap.df)
  
  # Outputs
  direction.plots <- list()
  direction.m <- c()
  
  for(model in models){
    
    source('./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/8_SHAP_directionality.R', local = T, verbose = F)
    
  }
  
  direction.m <- as.data.frame(direction.m)
  colnames(direction.m) <- c("Model","Feature","Pearson","Spearman")
  direction.m$Pearson <- as.numeric(direction.m$Pearson)
  direction.m$Spearman <- as.numeric(direction.m$Spearman)
  
  save(direction.m, direction.plots, file = "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/SHAP/SHAP_directions_HIPPIE.RData")
}
