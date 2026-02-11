
# Date: September 18, 2024

# Preparing ML tables with all features together - Part III (Tissue-specific features)

rm(list = ls())
gc()

# Tissue specific features (from GS)
# 3. Cancer drivers

cancer_drivers_path <- "./Data/Distance_to_closest_cancer-drivers_intogen.RData"
cancer_drivers_weighted_path <- "./Data/Distance_to_cancer-drivers_intogen_weighted.RData"
cancer_drivers_density_path <- "./Data/Density_cancer-drivers_intogen.RData"
backbone_path <- "./Data/Backbone_tables_with_all_features_GS.RData"

# Cancer drivers
# Distance
load(cancer_drivers_path)
load(cancer_drivers_weighted_path)

# Density and weighted density
load(cancer_drivers_density_path)

intogen.cohorts <- names(Density.drivers$`0.1Mbp`)

# Feature data with the features previously added
load(backbone_path)
levels <- names(backbone_tables_w_all_features_GS_complete)

backbone_tables_w_all_features_GS_updated <- list()

for(level in levels){
  cohorts <- names(backbone_tables_w_all_features_GS_complete[[level]])
  for(cohort in cohorts){
    
    m <- backbone_tables_w_all_features_GS_complete[[level]][[cohort]]
    
    if(cohort %in% intogen.cohorts){
      
      # Density and weighted density
      feature1 <- Density.drivers[[level]][[cohort]]
      
      m <- merge(m, feature1, by = "bin", all.x = T)
      
      if(!level %in% c("Chromosome","Arm")){
        # Distance to cancer drivers I
        feature2 <- merge(Dist.oncogenes[[level]][[cohort]], Dist.tumoursupp[[level]][[cohort]], by = "bin", all = T)
        feature2 <- feature2[,c("bin","dist.to.closest.OG","dist.to.closest.TSG")]
        
        # Distance to cancer drivers II (Weighted)
        feature3 <- merge(Dist.oncogenes.weighted[[level]][[cohort]], Dist.tumoursupp.weighted[[level]][[cohort]], by = "bin", all = T)
        feature2 <- merge(feature2, feature3, by = "bin", all = T)
        
        m <- merge(m, feature2, by = "bin", all.x = T)
        
      }
      backbone_tables_w_all_features_GS_updated[[level]][[cohort]] <- m
      
    }
    else{
      backbone_tables_w_all_features_GS_updated[[level]][[cohort]] <- m
    }
    
    
  }
}

save(backbone_tables_w_all_features_GS_updated, file = backbone_path)
