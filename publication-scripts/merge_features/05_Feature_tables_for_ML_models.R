
# Date: November 12, 2024 (Updated)
# Preparing feature matrix for ML models

rm(list = ls())
gc()

# For 11 cancer types that we have all the features available

# Features datasets

backbone_path <- "./Data/Backbone_tables_with_all_features_GS.RData"
backbone_path_FA <- "./Data/Backbone_tables_with_features_FA.RData"
outpath_HIPPIE <- "./Data/Feature_tables_for_ML_models_11cohorts_HIPPIE.RData"
outpath_IntINSIDER <- "./Data/Feature_tables_for_ML_models_11cohorts_IntINSIDER.RData"

load(backbone_path)
load(backbone_path_FA)

common.levels <- intersect(names(backbone_tables_w_all_features_GS_updated),names(Features.FA))
common.cohorts <- c("BRCA","COADREAD","ESCA","GBMLGG","KIRC","KIRP","LUSC","LUAD","OV","PAAD","STAD") # All features available for those


# A. Prepare feature lists - HIPPIE

# 1. Non-tissue-specific features
# Paralogs and Ohnologs
# Ensembl Biomart, HIPPIE and CORUM related features
# Distance to fragile sites

non.specific <- c("n_genes","mean.GC.content", # Ensembl
                  "total.complex.subunit","total_n_partners.trans", # CORUM
                  "total_n_PPIs.trans_HIPPIE_073", # HIPPIE
                  "n_ohnolog_bin",
                  "total_n_ohnologs.mmpaper_trans","total_n_paralogs_trans",
                  "dist.to.closest.FGS"
                  ) # Some features' names will be modified after normalisation (divided by number of n_genes)

# 2. Tissue-specific features (for 31 levels)
# 2.1. Chromatin states
states <- paste0("E",seq(1,25))
chromatin.states.groups <- paste0("Length_Counts.",states)

# 2.2. GTEx mRNA and protein levels

gtex.groups <- c("genes.bin","partners.trans","all.int.trans_HIPPE_073")

# 2.3. Cancer drivers

cancer.drivers <- c("density.OG","density.TSG","dist.to.closest.OG","dist.to.closest.TSG") # Density actually means the number of drivers

# 3. Features from 

features.fa <- c("mutations_norm", "distance.to.centromere","distance.to.telomere",
                 "Ess.distance_pancancer","ESSscore_pancancer", "HAPLOscore_pancancer")


Feature.tables.for.ML <- list()

for(level in common.levels){
  level.m <- c()
  for(cohort in common.cohorts){
    gs.features.m <- backbone_tables_w_all_features_GS_updated[[level]][[cohort]]
    gs.features.m <- gs.features.m[,colnames(gs.features.m) %in% c("bin",non.specific,
                                      chromatin.states.groups,
                                      gtex.groups,
                                      cancer.drivers)]
    fa.features.m <- Features.FA[[level]][[cohort]]
    fa.features.m <- fa.features.m[,colnames(fa.features.m) %in% c("bin",features.fa)]
    
    features.all <- merge(gs.features.m,fa.features.m, by = "bin")
    features.all$Type <- cohort
    level.m <- rbind(level.m,features.all)
  }
  
  # Normalisation of some features
  level.m$Density.complex.proteins <- level.m$total.complex.subunit / level.m$n_genes
  level.m$Density.Ohnologs <- level.m$n_ohnolog_bin / level.m$n_genes
  level.m$density.OG <- level.m$density.OG / level.m$n_genes
  level.m$density.TSG <- level.m$density.TSG / level.m$n_genes
  
  level.m <- level.m[, !colnames(level.m) %in% c("total.complex.subunit","n_ohnolog_bin","n_genes")]

  Feature.tables.for.ML[[level]] <- level.m
}

save(Feature.tables.for.ML, file = outpath_HIPPIE)

# A. Prepare feature lists - Interactome INSIDER

# 1. Non-tissue-specific features
# Paralogs and Ohnologs
# Ensembl Biomart, HIPPIE and CORUM related features
# Distance to fragile sites

non.specific <- c("n_genes","mean.GC.content", # Ensembl
                  "total.complex.subunit","total_n_partners.trans", # CORUM
                  "total_n_PPIs.trans_IntINSIDER", # Interactome INSIDER
                  "n_ohnolog_bin",
                  "total_n_ohnologs.mmpaper_trans","total_n_paralogs_trans",
                  "dist.to.closest.FGS"
) # Some features' names will be modified after normalisation (divided by number of n_genes)

# 2. Tissue-specific features (for 31 levels)
# 2.1. Chromatin states
states <- paste0("E",seq(1,25))
chromatin.states.groups <- paste0("Length_Counts.",states)

# 2.2. GTEx mRNA and protein levels

gtex.groups <- c("genes.bin","partners.trans","all.int.trans_IntINSIDER")

# 2.3. Cancer drivers

cancer.drivers <- c("density.OG","density.TSG","dist.to.closest.OG","dist.to.closest.TSG") # Density actually means the number of drivers

# 3. Features from 

features.fa <- c("mutations_norm", "distance.to.centromere","distance.to.telomere",
                 "Ess.distance_pancancer","ESSscore_pancancer", "HAPLOscore_pancancer")


Feature.tables.for.ML <- list()

for(level in common.levels){
  level.m <- c()
  for(cohort in common.cohorts){
    gs.features.m <- backbone_tables_w_all_features_GS_updated[[level]][[cohort]]
    gs.features.m <- gs.features.m[,colnames(gs.features.m) %in% c("bin",non.specific,
                                                                   chromatin.states.groups,
                                                                   gtex.groups,
                                                                   cancer.drivers)]
    fa.features.m <- Features.FA[[level]][[cohort]]
    fa.features.m <- fa.features.m[,colnames(fa.features.m) %in% c("bin",features.fa)]
    
    features.all <- merge(gs.features.m,fa.features.m, by = "bin")
    features.all$Type <- cohort
    level.m <- rbind(level.m,features.all)
  }
  
  # Normalisation of some features
  level.m$Density.complex.proteins <- level.m$total.complex.subunit / level.m$n_genes
  level.m$Density.Ohnologs <- level.m$n_ohnolog_bin / level.m$n_genes
  level.m$density.OG <- level.m$density.OG / level.m$n_genes
  level.m$density.TSG <- level.m$density.TSG / level.m$n_genes
  
  level.m <- level.m[, !colnames(level.m) %in% c("total.complex.subunit","n_ohnolog_bin","n_genes")]
  
  Feature.tables.for.ML[[level]] <- level.m
}

save(Feature.tables.for.ML, file = outpath_IntINSIDER)

