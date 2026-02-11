
# Date: September 18, 2024
# Preparing feature matrix for ML models

rm(list = ls())

# Set working directory
setwd("/mnt/fabiogokce/gokce")

# Required libraries
library(maditr)

# A. Availability of features
load("./Data/Backbone_tables_with_all_features_GS.RData")
load("./Data/Backbone_tables_with_features_FA.RData") # No 3Mbp, that is why it is not in the common levels, and I did not run the ML Models

# For an example focal level 
#data <- Features.FA[["1Mbp"]]
data <- backbone_tables_w_all_features_GS_updated[["1Mbp"]]

features.df <- c()
#Available features for different cancer types
for(type in names(data)){
  features.df <- rbind(features.df,
                       data.frame("Features" = colnames(data[[type]]),
                                  "Type" = rep(type,ncol(data[[type]]))))}
features <- setdiff(unique(features.df$Features),c("key","bin","chr","start_bin","end_bin","genes"))
features.info.gs <- c()
features.l <- list()
# Find the cancer types where all the features are available
for(feature in features){
  cohorts <- features.df[features.df$Features == feature,]$Type
  features.info.gs <- rbind(features.info.gs,
                            c(feature, length(cohorts), paste(cohorts,collapse = "::")))
  features.l[[feature]] <- cohorts
}
common.types <- Reduce(intersect, features.l)

# Note: Features from FA: Tissue-specific features are only for 6 cancer types, there is no tissue specificity for haploinsufficient genes
# Features from GS: Common types = 11

# B. Prepare feature lists

rm(list = setdiff(ls(),"common.types"))

# 1. Non-tissue-specific features (for 55 levels)
# Paralogs and Ohnologs
# Ensembl Biomart and CORUM related features
# Distance to fragile sites - 2 options: "dist.to.closest.FGS","Dist.FGS.weighted" 

default.cols <- c("mean.GC.content","Density.complex.proteins","total_n_partners.trans","total_n_PPIs.trans",
                  "Density.Ohnologs","total_n_ohnologs.mmpaper_trans","total_n_paralogs_trans") # Some features' names are the one after normalisation

dist.FGS <- list("Closest.FGS" = c("dist.to.closest.FGS"),
                 "Weighted.FGS" = c("Dist.FGS.weighted"))

# 2. Tissue-specific features (for 31 levels)
# 2.1. Chromatin states
states <- paste0("E",seq(1,25))

chromatin.states.groups <- list("Abundance" = paste0("Abundance_Counts.",states),
                                #"Abundance.norm" = paste0("Abundance_Normalised.by.chr.coverage.",states),
                                "Length" = paste0("Length_Counts.",states)
                                #,"Length.norm" = paste0("Length_Normalised.by.chr.coverage.",states),
                                #"Binary" = paste0("Abundance_binary.",states)
)

# 2.2. GTEx mRNA and protein levels

gtex.groups <- list("GTEx.levels" = c("genes.bin","partners.trans","all.int.trans"),
                    "GTEx.levels.TS" = c("genes.bin_TS","partners.trans_TS","all.int.trans_TS"))

# 2.3. Cancer drivers

cancer.driver.density <- list("Cancer.driver.density" = c("density.OG","density.TSG"),
                              "Cancer.driver.weighted.density" = c("Weighted.density.OG","Weighted.density.TSG"))
cancer.driver.dist <- list("Closest.cancer.driver" = c("dist.to.closest.OG","dist.to.closest.TSG"),
                           "Weighted.cancer.driver" = c("Dist.OG.weighted","Dist.TSG.weighted"))

# 3. Features from 

features.fa <- list("Features.FA.pancancer" = c("mutations_norm", "distance.to.centromere","distance.to.telomere",
                                                "Ess.distance_pancancer","ESSscore_pancancer", "HAPLOscore_pancancer"))


# C. Prepare feature tables - all features together
feature.combs <- list()
for(name1 in names(dist.FGS)){
  for(name2 in names(chromatin.states.groups)){
    for(name3 in names(gtex.groups)){
      for(name5 in names(cancer.driver.density)){
        for(name6 in names(cancer.driver.dist)){
          for(name4 in names(features.fa)){
            all.combn <- c(default.cols,
                           dist.FGS[[name1]],
                           chromatin.states.groups[[name2]],
                           gtex.groups[[name3]],
                           cancer.driver.density[[name5]],
                           cancer.driver.dist[[name6]],
                           features.fa[[name4]])
            feature.combs[[paste(name1,name2,name3,name5,name6,name4,sep = "_")]] <- all.combn
          }
          
        }
      }
      
    }
  }
}


load("./Data/Backbone_tables_with_all_features_GS.RData")
load("./Data/Backbone_tables_with_features_FA.RData")

common.levels <- intersect(names(backbone_tables_w_all_features_GS_updated),
                           names(Features.FA))
common.cohorts <- common.types# All features available for those

Feature.tables <- list()

for(level in common.levels){
  level.m <- c()
  for(cohort in common.cohorts){
    gs.features.m <- backbone_tables_w_all_features_GS_updated[[level]][[cohort]]
    fa.features.m <- Features.FA[[level]][[cohort]]
    fa.features.m <- fa.features.m[,colnames(fa.features.m) %in% c("bin","mutations_norm", "distance.to.centromere","distance.to.telomere",
                                                                   "Ess.distance_pancancer","ESSscore_pancancer", "HAPLOscore_pancancer")]
    features.all <- merge(gs.features.m,fa.features.m, by = "bin")
    features.all$Type <- cohort
    level.m <- rbind(level.m,features.all)
  }
  # Normalisation of two features
  level.m$Density.complex.proteins <- level.m$total.complex.subunit / level.m$n_genes
  level.m$Density.Ohnologs <- level.m$n_ohnolog_bin / level.m$n_genes
  
  for(name in names(feature.combs)){
    selected.features <- feature.combs[[name]]
    level.m.new <- level.m[,colnames(level.m) %in% c("bin",selected.features,"Type")]
    Feature.tables[[level]][[name]] <- level.m.new
  }
}

save(Feature.tables, file = "./Data/Feature_tables_different_feature_combn_updated_cancer_drivers_partI.RData")

# D. Prepare feature tables - all features together - for HIPPIE features

rm(list = setdiff(ls(),"common.types"))

# Feature tables with all combinations
load("./Data/Feature_tables_different_feature_combn_updated_cancer_drivers_partI.RData") # Feature.tables

# GTEx data for HIPPIE
load("./Data/GTEx_rna_protein_level_differentdatasets_per_bins_HIPPIE.RData") # Results
load("./Data/GTEx_rna_protein_level_differentdatasets_per_bins_Tissue_specific_HIPPIE.RData") # Results.TS
# Tissue-general PPIs number from HIPPIE
load("./Data/Backbone_tables_with_non_tissue_specific_features_CORUM_HIPPIE.RData") # backbone_tables_w_features

Feature.tables.hippie <- list()
for(level in names(Feature.tables)){
  for(combin in names(Feature.tables[[level]])){
    # ML table with final features
    df <- Feature.tables[[level]][[combin]]
    # Remove features based on Interactome Insider
    df <- df[,!colnames(df) %in% c("total_n_PPIs.trans","all.int.trans","all.int.trans_TS"),]
    
    # Find if gtex feature is tissue-specific or not 
    gtex.fea <- strsplit(combin,"_")[[1]][3]
    # Prepararing the data for features based on HIPPIE
    if(gtex.fea == "GTEx.levels"){
      hippie.fea <- Results[[level]][["scaled.GTEx.v8.TPM"]]
      hippie.data <- c()
      for(cohort in unique(df$Type)){
        gtex <- hippie.fea[hippie.fea$Cohorts == cohort,]
        gtex <- gtex %>% dcast(bin ~ GeneSet, value.var = "Median")
        gtex$Type <- cohort
        colnames(gtex)[which(colnames(gtex) == "all.int.trans")] <- "all.int.trans.HIPPIE"
        hippie.data <- rbind(hippie.data, gtex[,c("bin","all.int.trans.HIPPIE","Type")])
      }
    }
    else{
      hippie.fea <- Results.TS[[level]][["scaled.GTEx.v8.TPM"]]
      hippie.data <- c()
      for(cohort in unique(df$Type)){
        gtex <- hippie.fea[hippie.fea$Cohorts == cohort,]
        gtex <- gtex %>% dcast(bin ~ GeneSet, value.var = "Median")
        gtex$Type <- cohort
        colnames(gtex)[which(colnames(gtex) == "all.int.trans")] <- "all.int.trans.TS.HIPPIE"
        hippie.data <- rbind(hippie.data, gtex[,c("bin","all.int.trans.TS.HIPPIE","Type")])
      }
      
    }
    # Tissue-general number of PPIs based on HIPPIE
    n_ppis_hippie <- backbone_tables_w_features[[level]]
    hippie.data <- merge(hippie.data, n_ppis_hippie[,c("bin","total_n_PPIs.trans-HIPPIE")], by = "bin", all.x = T)
    # Add HIPPIE features to the final feature table
    df <- merge(df,hippie.data, by = c("bin","Type"), all.x = T)
    combin <- paste("HIPPIE",combin,sep = "_")
    Feature.tables.hippie[[level]][[combin]] <- df
  }
}

# Joining the two lists
Final.feature.tables <- list()
for(level in names(Feature.tables)){
  list1 <- Feature.tables[[level]]
  list2 <- Feature.tables.hippie[[level]]
  joined <- c(list1, list2)
  Final.feature.tables[[level]] <- joined
}

save(Final.feature.tables, 
     file = "./Data/Feature_tables_different_feature_combn_updated_cancer_drivers_and_HIPPIE.RData")
