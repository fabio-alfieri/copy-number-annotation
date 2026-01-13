# Directionality

# Color codes
background <- "#e1e5f2"
lines <- "#012a4a"
bars <- "#a4ac86" # Feature importance plots/bar plots

wd <- 'path/to/GitHub/copy-number-annotation/'

# Data - Interactome INSIDER
shap_path <- paste0(wd, "data/SHAP/")
avg_shap_path <- paste0(shap_path,class,"::",model,"_Avg_shap_values_InteractomeINSIDER.RData")

directionality_script <- paste0('publication-scripts/machine-learning/8_SHAP_directionality.R')

load(avg_shap_path) 
models <- names(models.shap.df)

# Outputs
direction.plots <- list()
direction.m <- c()

for(model in models){
  print(model)
  source(directionality_script, local = T, verbose = F)
}

direction.m <- as.data.frame(direction.m)
colnames(direction.m) <- c("Model","Feature","Pearson","Pearson_p.val", "Spearman", "Spearman_p.val")
direction.m$Pearson <- as.numeric(direction.m$Pearson)
direction.m$Spearman <- as.numeric(direction.m$Spearman)
direction.m$Pearson_p.val <- as.numeric(direction.m$Pearson_p.val)
direction.m$Spearman_p.val <- as.numeric(direction.m$Spearman_p.val)

directons_path <- paste0(wd, "data/SHAP/", class, "_SHAP_directions_InteractomeINSIDER.RData")

save(direction.m, 
     direction.plots, 
     file = directons_path)
