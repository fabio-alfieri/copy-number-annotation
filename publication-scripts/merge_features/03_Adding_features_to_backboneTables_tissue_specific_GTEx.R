
# Preparing ML tables with all features together - Part II (Tissue-specific features)

rm(list = ls())
gc()

# Tissue specific features (from GS)
# 2. GTEx mRNA and protein levels

# Set working directory 
setwd("/mnt/fabiogokce/gokce")

# Required libraries
library(maditr) # for dcast

# PART B
# GTEx data (Calculated based on Interactome INSIDER interactions)

load("./Data/GTEx_rna_protein_level_differentdatasets_per_bins.RData") # Results
load("./Data/GTEx_rna_protein_level_differentdatasets_per_bins_Tissue_specific.RData") # Results.TS

# Cohorts that GTEx data is available - Skip this part (update on September 18, 2024)
#gtex.cohorts <- unique(Results$`1Mbp`$eGTEX.rna$Cohorts) # Since the data prepared for the cohorts common in all the different datasets

# Dataset to be used
dataset <- "scaled.GTEx.v8.TPM"
gtex.cohorts <- unique(Results$`1Mbp`$scaled.GTEx.v8.TPM$Cohorts)

# Feature data with the features previously added 
load("./Data/Backbone_tables_with_all_features_GS_except_GTEx.RData")

levels <- names(backbone_tables_w_all_features_GS)

backbone_tables_w_all_features_GS_complete <- list()

for(level in levels){
  cohorts <- names(backbone_tables_w_all_features_GS[[level]])
  for(cohort in cohorts){
    
    m <- backbone_tables_w_all_features_GS[[level]][[cohort]]
    
    if(cohort %in% gtex.cohorts){
      
      # GTEx data (Some rows could be missing since for bins without any gene, there is no GTEx feature, assign NA to those bins)
      gtex <- Results[[level]][[dataset]]
      gtex <- gtex[gtex$Cohorts == cohort,]
      gtex <- gtex %>% reshape2::dcast(bin ~ GeneSet, value.var = "Median")
      colnames(gtex)[which(colnames(gtex) == "genes")] <- "genes.bin"
      colnames(gtex)[which(colnames(gtex) == "all.int.trans")] <- "all.int.trans_IntINSIDER"
      colnames(gtex)[which(colnames(gtex) == "ppis.trans")] <- "ppis.trans_IntINSIDER"
      # GTEx data (Some rows could be missing since for bins without any gene, there is no GTEx feature, assign NA to those bins)
      gtexTS <- Results.TS[[level]][[dataset]]
      gtexTS <- gtexTS[gtexTS$Cohorts == cohort,]
      gtexTS <- gtexTS %>% reshape2::dcast(bin ~ GeneSet, value.var = "Median")
      colnames(gtexTS)[which(colnames(gtexTS) == "genes")] <- "genes.bin"
      colnames(gtexTS)[2:5] <- paste(colnames(gtexTS)[2:5],"TS",sep = "_")
      colnames(gtexTS)[which(colnames(gtexTS) == "all.int.trans_TS")] <- "all.int.trans_TS_IntINSIDER"
      colnames(gtexTS)[which(colnames(gtexTS) == "ppis.trans_TS")] <- "ppis.trans_TS_IntINSIDER"
      
      m <- merge(m, gtex, by = "bin", all.x = T)
      m <- merge(m, gtexTS, by = "bin", all.x = T)
      
      backbone_tables_w_all_features_GS_complete[[level]][[cohort]] <- m
      
    }
    else{
      backbone_tables_w_all_features_GS_complete[[level]][[cohort]] <- m
    }
    
  
    }
}

save(backbone_tables_w_all_features_GS_complete, file = "./Data/Backbone_tables_with_all_features_GS.RData")

rm(list = ls())
gc()

# PART B
# GTEx data (Calculated based on HIPPIE interactions -  no cutoff)

load("./Data/GTEx_rna_protein_level_differentdatasets_per_bins_HIPPIE.RData") # Results
load("./Data/GTEx_rna_protein_level_differentdatasets_per_bins_Tissue_specific_HIPPIE.RData") # Results.TS

# Dataset to be used
dataset <- "scaled.GTEx.v8.TPM"
gtex.cohorts <- unique(Results$`1Mbp`$scaled.GTEx.v8.TPM$Cohorts)

# Feature data with the features previously added 
load("./Data/Backbone_tables_with_all_features_GS.RData")

levels <- names(backbone_tables_w_all_features_GS_complete)

for(level in levels){
  cohorts <- names(backbone_tables_w_all_features_GS_complete[[level]])
  for(cohort in cohorts){
    
    m <- backbone_tables_w_all_features_GS_complete[[level]][[cohort]]
    
    if(cohort %in% gtex.cohorts){
      
      # GTEx data (Some rows could be missing since for bins without any gene, there is no GTEx feature, assign NA to those bins)
      gtex <- Results[[level]][[dataset]]
      gtex <- gtex[gtex$Cohorts == cohort,]
      gtex <- gtex %>% reshape2::dcast(bin ~ GeneSet, value.var = "Median")
      colnames(gtex)[2:3] <- paste(colnames(gtex)[2:3],"HIPPE",sep = "_")
      # GTEx data (Some rows could be missing since for bins without any gene, there is no GTEx feature, assign NA to those bins)
      gtexTS <- Results.TS[[level]][[dataset]]
      gtexTS <- gtexTS[gtexTS$Cohorts == cohort,]
      gtexTS <- gtexTS %>% reshape2::dcast(bin ~ GeneSet, value.var = "Median")
      colnames(gtexTS)[2:3] <- paste(colnames(gtexTS)[2:3],"HIPPIE_TS",sep = "_")
      
      m <- merge(m, gtex, by = "bin", all.x = T)
      m <- merge(m, gtexTS, by = "bin", all.x = T)
      
      backbone_tables_w_all_features_GS_complete[[level]][[cohort]] <- m
      
    }
    else{
      backbone_tables_w_all_features_GS_complete[[level]][[cohort]] <- m
    }
    
    
  }
}


# PART C
# GTEx data (Calculated based on HIPPIE interactions -  cut off 0.73)

rm(list = setdiff(ls(),"backbone_tables_w_all_features_GS_complete"))
gc()

levels <- names(backbone_tables_w_all_features_GS_complete)

load("./Data/GTEx_rna_protein_level_differentdatasets_per_bins_HIPPIE_cutoff_073.RData") # Results

# Dataset to be used
dataset <- "scaled.GTEx.v8.TPM"
gtex.cohorts <- unique(Results$`1Mbp`$scaled.GTEx.v8.TPM$Cohorts)

for(level in levels){
  cohorts <- names(backbone_tables_w_all_features_GS_complete[[level]])
  for(cohort in cohorts){
    
    m <- backbone_tables_w_all_features_GS_complete[[level]][[cohort]]
    
    if(cohort %in% gtex.cohorts){
      
      # GTEx data (Some rows could be missing since for bins without any gene, there is no GTEx feature, assign NA to those bins)
      gtex <- Results[[level]][[dataset]]
      gtex <- gtex[gtex$Cohorts == cohort,]
      gtex <- gtex %>% reshape2::dcast(bin ~ GeneSet, value.var = "Median")
      colnames(gtex)[2:3] <- paste(colnames(gtex)[2:3],"HIPPE_073",sep = "_")
      
      m <- merge(m, gtex, by = "bin", all.x = T)
      
      backbone_tables_w_all_features_GS_complete[[level]][[cohort]] <- m
      
    }
    else{
      backbone_tables_w_all_features_GS_complete[[level]][[cohort]] <- m
    }
    
    
  }
}

save(backbone_tables_w_all_features_GS_complete, file = "./Data/Backbone_tables_with_all_features_GS.RData")
