# rm(list=ls())
gc(full=T)

setwd('/mnt/fabiogokce/fabiohd/')

chr_backbone_genes_cancer <- readRDS("ml_models/data/misc/chr_backbone_genes_cancer.rds")
mut_cna_features <- readRDS("ml_models/data/mut_cna_features.rds")

mut_cna_features_mod <- list()
for(i in names(mut_cna_features)){
  chr_backbone_genes_cancer[[i]][,c('chr','bin')] <- do.call(rbind, str_split(chr_backbone_genes_cancer[[i]]$bin, '_'))
  chr_backbone_genes_cancer[[i]][,c('start_bin','end_bin','chr','bin')] <- apply(
    chr_backbone_genes_cancer[[i]][,c('start_bin','end_bin','chr','bin')], 2, as.numeric
  )
  common_columns <- c(1:6,7,14,21,28,35,36)
  LUAD_table <- chr_backbone_genes_cancer[[i]][,c(common_columns, grep("LUAD", colnames(chr_backbone_genes_cancer[[i]])))]
  LUSC_table <- chr_backbone_genes_cancer[[i]][,c(common_columns, grep("LUSC", colnames(chr_backbone_genes_cancer[[i]])))]
  COADREAD_table <- chr_backbone_genes_cancer[[i]][,c(common_columns, grep("COADREAD", colnames(chr_backbone_genes_cancer[[i]])))]
  GBMLGG_table <- chr_backbone_genes_cancer[[i]][,c(common_columns, grep("GBMLGG", colnames(chr_backbone_genes_cancer[[i]])))]
  BRCA_table <- chr_backbone_genes_cancer[[i]][,c(common_columns, grep("BRCA", colnames(chr_backbone_genes_cancer[[i]])))]
  PAAD_table <- chr_backbone_genes_cancer[[i]][,c(common_columns, grep("PAAD", colnames(chr_backbone_genes_cancer[[i]])))]
  LUAD_table$type <- 'LUAD'
  LUSC_table$type <- 'LUSC'
  COADREAD_table$type <- 'COADREAD'
  GBMLGG_table$type <- 'GBMLGG'
  BRCA_table$type <- 'BRCA'
  PAAD_table$type <-'PAAD'
  
  colnames(LUAD_table)[13:17] <- c("OGscore_TS","OGscore.potency_TS","TSGscore_TS",           
                                   "TSGscore.potency_TS","essential_TS")
  colnames(LUSC_table)[13:17] <- c("OGscore_TS","OGscore.potency_TS","TSGscore_TS",           
                                   "TSGscore.potency_TS","essential_TS")
  colnames(COADREAD_table)[13:17] <- c("OGscore_TS","OGscore.potency_TS","TSGscore_TS",           
                                   "TSGscore.potency_TS","essential_TS")
  colnames(GBMLGG_table)[13:17] <- c("OGscore_TS","OGscore.potency_TS","TSGscore_TS",           
                                   "TSGscore.potency_TS","essential_TS")
  colnames(BRCA_table)[13:17] <- c("OGscore_TS","OGscore.potency_TS","TSGscore_TS",           
                                   "TSGscore.potency_TS","essential_TS")
  colnames(PAAD_table)[13:17] <- c("OGscore_TS","OGscore.potency_TS","TSGscore_TS",           
                                   "TSGscore.potency_TS","essential_TS")
  
  chr_backbone_genes_cancer_ts <- rbind(LUAD_table,LUSC_table,COADREAD_table,
                                        BRCA_table,PAAD_table,GBMLGG_table)
  
  mut_cna_features_filt <- as.data.frame(do.call(rbind, mut_cna_features[[i]]) %>% 
                                           filter(type == 'BRCA' |type == 'LUSC' |type == 'LUAD' |
                                                    type == 'PAAD' | type == 'COADREAD' | type == 'GBMLGG'))
  mut_cna_features_filt[,c('start_bin','end_bin','chr','bin')] <- apply(
    mut_cna_features_filt[,c('start_bin','end_bin','chr','bin')], 2, as.numeric
  )
  
  merged <- full_join(mut_cna_features_filt,
                      chr_backbone_genes_cancer_ts,
            by = c('chr','bin','start_bin','end_bin','type'))
  
  mut_cna_features_mod[[i]] <- merged
}
rm(PAAD_table,BRCA_table,COADREAD_table,GBMLGG_table,LUSC_table,LUAD_table,
   merged,mut_cna_features_filt,chr_backbone_genes_cancer_ts,i,common_columns)

saveRDS(mut_cna_features_mod, file = 'ml_models/data/ml_merged.rds')

