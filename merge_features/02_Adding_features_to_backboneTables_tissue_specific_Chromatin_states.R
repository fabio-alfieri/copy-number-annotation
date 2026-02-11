
#Preparing ML tables with all features together - Part II (Tissue-specific features)

rm(list = ls())
gc()

# Tissue specific features (from GS)
# 1. Chromatin States

# Set working directory 
setwd("/mnt/fabiogokce/gokce") #If you're working on workstation

#Required libraries
library(data.table)
library(dplyr)

#NOTE
#GTEx data was prepared for all the available TCGA match. To get the data for the cohorts that we are interested, use the cohorts below. 
#Chromatin states were prepared for the cohorts common between Roadmap and cohorts of interest (the vector below)
#If we don't have available tissue specific data for the corresponding cohort, use backbone table with non tissue-specific features

cohorts <- c("BRCA","LUAD","LUSC","CESC","THCA","HNSC","PAAD","COADREAD","GBMLGG",
             "SKCM","BLCA","PCPG","PRAD","KIRC","MESO","TGCT","KIRP","SARC",    
             "LIHC","ESCA","STAD","UCS","OV") #Cohorts of interest - 23 cohorts

#Features: 1. Chromatin states (Available data for 14 tissues)
load("./Data/Chromatin_states_all_levels_tissues_all_methods.RData") #Abundance and length method, together with normalized by chromosome coverage! 
load("./Data/Backbone_tables_with_non_tissue_specific_features.RData") #Backbone tables with non-tissue specific features

levels <- names(Chromatin.States)
backbone_tables_w_all_features_GS <- list()

for(level in levels){
  level.dataset <- Chromatin.States[[level]]
  level.backbone.features <- backbone_tables_w_features_non_tissue_specific[[level]]

  level.dataset.new <- list() #New list to store outputs from the loop
  
  for(tissue in cohorts){ #For each tissue, prepare the matrix, if available, and merge it with non tissue-specific features
    if(tissue %in% names(level.dataset)){ #Check if tissue-specific data (Chromatin states) is available for the tissue
      tissue.df <- do.call(cbind, level.dataset[[tissue]])
      tissue.df$bin <- rownames(tissue.df)
      tissue.df <- merge(level.backbone.features, tissue.df, by = "bin")}
    else{tissue.df <- level.backbone.features} #If tissue-specific data is not available, assign tissue.df to the backbone table with non tissue-specific features
    level.dataset.new[[tissue]] <- tissue.df}
  backbone_tables_w_all_features_GS[[level]] <- level.dataset.new}

rm(list = setdiff(ls(),c("cohorts","backbone_tables_w_all_features_GS")))

save(backbone_tables_w_all_features_GS,
     file = "./Data/Backbone_tables_with_all_features_GS_except_GTEx.RData")

