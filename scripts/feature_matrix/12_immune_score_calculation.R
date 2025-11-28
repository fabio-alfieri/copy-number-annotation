# rm(list=ls())
gc(full=T)

library(tidyverse)
library(rstatix)
library(reshape2)
library(parallel)
library(ggplot2)

for(tt in c(
  "BRCA",
  "LUAD",
  "LUSC",
  "CESC",
  "THCA",
  "HNSC",
  "PAAD",
  "COADREAD",
  "GBMLGG",
  # # #
  "SKCM",
  "BLCA",
  "PCPG",
  "PRAD",
  "KIRC",
  "KIRP",
  "MESO",
  "TGCT",
  "SARC",
  "LIHC",
  "ESCA",
  "STAD",
  "UCS",
  "OV"
  )){
  print(tt)
  
  # Calculate an immune score to each patient
  tpm <- readRDS(file = paste0('/mnt/fabiogokce/fabiohd/mutation_compensation/data/TCGA_tpm/',tt,'_tpm.rds.gz'))
  tpm <- tpm[,c(2,4,8)]
  
  tpm <- tpm[!duplicated(paste0(tpm$patient_id,tpm$`Approved symbol`)),]
  
  colnames(tpm)[2] <- 'symbol'
  
  genes = c("CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK")
  tpm <- tpm[tpm$symbol %in% genes,]
  
  tpm_wide <- tpm %>% pivot_wider(names_from = symbol, values_from = tpm_counts)
  
  tpm_wide[,-1] <- apply(tpm_wide[,-1], 2, rank)
  tpm_wide$Immune_Score <- rank(rowSums(tpm_wide[,-1]))
  
  # Load backbone tables for specific bin size
  load("/mnt/fabiogokce/fabiohd/ml_models/ml_tables_from_0_prepare_tables_240329.RData")
  backbone <- ML.tables[['1Mbp']][[tt]][,c(1,4,5,6)]
  rm(ML.tables)
  backbone$bin <- do.call(rbind, strsplit(backbone$bin, '_'))[,2]
  backbone$bin <- as.numeric(backbone$bin)
  
  # Load patient specific CNAs
  patient_cna <- readRDS(paste0('/mnt/fabiogokce/fabiohd/US_intern/altered_bins/',tt,'_altered_bins_1Mbp.rds'))
  
  # load CIBERSORTx data
  if(tt == 'COADREAD'){
    ciberx1 <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/US_intern/data/CIBERSORTX_TCGA/',
                                       'COAD','_LM22_matrix_CIBERSORTX.tsv'), sep = '\t', header = T)
    ciberx2 <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/US_intern/data/CIBERSORTX_TCGA/',
                                       'READ','_LM22_matrix_CIBERSORTX.tsv'), sep = '\t', header = T)
    ciberx <- rbind(ciberx1,ciberx2)
    rm(ciberx1, ciberx2)
  }else if(tt =='GBMLGG'){
    ciberx1 <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/US_intern/data/CIBERSORTX_TCGA/',
                                        'GBM','_LM22_matrix_CIBERSORTX.tsv'), sep = '\t', header = T)
    ciberx2 <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/US_intern/data/CIBERSORTX_TCGA/',
                                        'LGG','_LM22_matrix_CIBERSORTX.tsv'), sep = '\t', header = T)
    ciberx <- rbind(ciberx1,ciberx2)
    rm(ciberx1, ciberx2)
  }else{
    ciberx <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/US_intern/data/CIBERSORTX_TCGA/',
                                       tt,'_LM22_matrix_CIBERSORTX.tsv'), sep = '\t', header = T)
  }
  ciberx <- ciberx[as.numeric(do.call(rbind, str_split(ciberx$Mixture, '-'))[,4]) <= 9,] # remove normal patients
  ciberx$Mixture <- substring(ciberx$Mixture, 1,12)
  
  # run analysis on immune_score
  statistics_immune_score <- mclapply(1:22, mc.cores = 22, function(chromosome){
    
    patient_ampl_tmp <- patient_cna[[tt]][['ampl']][[paste0('chr',chromosome)]]
    patient_del_tmp <- patient_cna[[tt]][['del']][[paste0('chr',chromosome)]]
    backbone_tmp <- backbone %>% filter(chr == chromosome) 
    
    stat_bin <- data.frame()
    for(i in backbone_tmp$bin){
      print(i)
      
      if(length(patient_ampl_tmp[[i]]) == 0){
        condition_ampl <- 'zero'
      }else{
        dimension.a <- nrow(patient_ampl_tmp[[i]])
        if(dimension.a <= 6){
          condition_ampl <- 'zero'
        }else(
          condition_ampl <- 'yes'
        )
      }
      if(length(patient_del_tmp[[i]]) == 0){
        condition_del <- 'zero'
      }else{
        dimension.d <- nrow(patient_del_tmp[[i]])
        if(dimension.d <= 6){
          condition_del <- 'zero'
        }else(
          condition_del <- 'yes'
        )
      }
      
      
      if(condition_ampl == 'zero' & condition_del == 'yes'){
        tpm_wide$cna_status <- ifelse(tpm_wide$patient_id %in% patient_del_tmp[[i]]$patient_id, -1, 0)
        ciberx$cna_status <- ifelse(ciberx$Mixture %in% patient_del_tmp[[i]]$patient_id, -1, 0)
        
        stat_bin <- rbind(stat_bin,
                          cbind(chr = chromosome,
                                bin = i,
                                kruskall.IS = kruskal.test(Immune_Score ~ as.factor(cna_status), data = tpm_wide)$p.value,
                                wil.ampl.IS = 1,
                                wil.del.IS = wilcox_test(tpm_wide, Immune_Score ~ cna_status)$p,
                                
                                kruskall.Dendritic_cells = kruskal.test(Dendritic_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Dendritic_cells = 1,
                                wil.del.Dendritic_cells = wilcox_test(ciberx, Dendritic_cells ~ cna_status)$p,
                                
                                kruskall.Macrophages = kruskal.test(Macrophages ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Macrophages = 1,
                                wil.del.Macrophages = wilcox_test(ciberx, Macrophages ~ cna_status)$p,
                                
                                kruskall.B_cells = kruskal.test(B_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.B_cells = 1,
                                wil.del.B_cells = wilcox_test(ciberx, B_cells ~ cna_status)$p,
                                
                                kruskall.CD4_T_cells = kruskal.test(CD4_T_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.CD4_T_cells = 1,
                                wil.del.CD4_T_cells = wilcox_test(ciberx, CD4_T_cells ~ cna_status)$p,
                                
                                kruskall.Mast_cells = kruskal.test(Mast_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Mast_cells = 1,
                                wil.del.Mast_cells = wilcox_test(ciberx, Mast_cells ~ cna_status)$p,
                                
                                kruskall.NK_cells = kruskal.test(NK_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.NK_cells = 1,
                                wil.del.NK_cells = wilcox_test(ciberx, NK_cells ~ cna_status)$p,
                                
                                kruskall.T_cells = kruskal.test(T_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.T_cells = 1,
                                wil.del.T_cells = wilcox_test(ciberx, T_cells ~ cna_status)$p,
                                
                                kruskall.Lymphocytes = kruskal.test(Lymphocytes ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Lymphocytes = 1,
                                wil.del.Lymphocytes = wilcox_test(ciberx, Lymphocytes ~ cna_status)$p,
                                
                                kruskall.Non_plasma_B_cells = kruskal.test(Non_plasma_B_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Non_plasma_B_cells = 1,
                                wil.del.Non_plasma_B_cells = wilcox_test(ciberx, Non_plasma_B_cells ~ cna_status)$p
                          ))
        
      }else if(condition_ampl == 'yes' & condition_del == 'zero'){
        tpm_wide$cna_status <- ifelse(tpm_wide$patient_id %in% patient_ampl_tmp[[i]]$patient_id, 1, 0)
        ciberx$cna_status <- ifelse(ciberx$Mixture %in% patient_ampl_tmp[[i]]$patient_id, 1, 0)
        
        stat_bin <- rbind(stat_bin,
                          cbind(chr = chromosome,
                                bin = i,
                                kruskall.IS = kruskal.test(Immune_Score ~ as.factor(cna_status), data = tpm_wide)$p.value,
                                wil.ampl.IS = wilcox_test(tpm_wide, Immune_Score ~ cna_status)$p,
                                wil.del.IS = 1,
                                
                                kruskall.Dendritic_cells = kruskal.test(Dendritic_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Dendritic_cells = wilcox_test(ciberx, Dendritic_cells ~ cna_status)$p,
                                wil.del.Dendritic_cells = 1,
                                
                                kruskall.Macrophages = kruskal.test(Macrophages ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Macrophages = wilcox_test(ciberx, Macrophages ~ cna_status)$p,
                                wil.del.Macrophages = 1,
                                
                                kruskall.B_cells = kruskal.test(B_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.B_cells = wilcox_test(ciberx, B_cells ~ cna_status)$p,
                                wil.del.B_cells = 1,
                                
                                kruskall.CD4_T_cells = kruskal.test(CD4_T_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.CD4_T_cells = wilcox_test(ciberx, CD4_T_cells ~ cna_status)$p,
                                wil.del.CD4_T_cells = 1,
                                
                                kruskall.Mast_cells = kruskal.test(Mast_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Mast_cells = wilcox_test(ciberx, Mast_cells ~ cna_status)$p,
                                wil.del.Mast_cells = 1,
                                
                                kruskall.NK_cells = kruskal.test(NK_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.NK_cells = wilcox_test(ciberx, NK_cells ~ cna_status)$p,
                                wil.del.NK_cells = 1,
                                
                                kruskall.T_cells = kruskal.test(T_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.T_cells = wilcox_test(ciberx, T_cells ~ cna_status)$p,
                                wil.del.T_cells = 1,
                                
                                kruskall.Lymphocytes = kruskal.test(Lymphocytes ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Lymphocytes = wilcox_test(ciberx, Lymphocytes ~ cna_status)$p,
                                wil.del.Lymphocytes = 1,
                                
                                kruskall.Non_plasma_B_cells = kruskal.test(Non_plasma_B_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Non_plasma_B_cells = wilcox_test(ciberx, Non_plasma_B_cells ~ cna_status)$p,
                                wil.del.Non_plasma_B_cells = 1))
        
      }else if(condition_ampl == 'zero' & condition_del == 'zero'){
        
        stat_bin <- rbind(stat_bin,
                          cbind(chr = chromosome,
                                bin = i,
                                kruskall.IS = 1,
                                wil.ampl.IS = 1,
                                wil.del.IS = 1,
                                
                                kruskall.Dendritic_cells = 1,
                                wil.ampl.Dendritic_cells = 1,
                                wil.del.Dendritic_cells = 1,
                                
                                kruskall.Macrophages = 1,
                                wil.ampl.Macrophages = 1,
                                wil.del.Macrophages = 1,
                                
                                kruskall.B_cells = 1,
                                wil.ampl.B_cells = 1,
                                wil.del.B_cells = 1,
                                
                                kruskall.CD4_T_cells = 1,
                                wil.ampl.CD4_T_cells = 1,
                                wil.del.CD4_T_cells = 1,
                                
                                kruskall.Mast_cells = 1,
                                wil.ampl.Mast_cells = 1,
                                wil.del.Mast_cells = 1,
                                
                                kruskall.NK_cells = 1,
                                wil.ampl.NK_cells = 1,
                                wil.del.NK_cells = 1,
                                
                                kruskall.T_cells = 1,
                                wil.ampl.T_cells = 1,
                                wil.del.T_cells = 1,
                                
                                kruskall.Lymphocytes = 1,
                                wil.ampl.Lymphocytes = 1,
                                wil.del.Lymphocytes = 1,
                                
                                kruskall.Non_plasma_B_cells = 1,
                                wil.ampl.Non_plasma_B_cells = 1,
                                wil.del.Non_plasma_B_cells = 1))
        
      }else{
        tpm_wide$cna_status <- apply(cbind(ifelse(tpm_wide$patient_id %in% patient_ampl_tmp[[i]]$patient_id, 1, 0),
                                           ifelse(tpm_wide$patient_id %in% patient_del_tmp[[i]]$patient_id, -1, 0)), 1, sum)
        ciberx$cna_status <- apply(cbind(ifelse(ciberx$Mixture %in% patient_ampl_tmp[[i]]$patient_id, 1, 0),
                                         ifelse(ciberx$Mixture %in% patient_del_tmp[[i]]$patient_id, -1, 0)), 1, sum)
        
        stat_bin <- rbind(stat_bin,
                          cbind(chr = chromosome,
                                bin = i,
                                kruskall.IS = kruskal.test(Immune_Score ~ as.factor(cna_status), data = tpm_wide)$p.value,
                                wil.ampl.IS = wilcox_test(tpm_wide[tpm_wide$cna_status != -1,], Immune_Score ~ cna_status)$p,
                                wil.del.IS = wilcox_test(tpm_wide[tpm_wide$cna_status != 1,], Immune_Score ~ cna_status)$p,
                                
                                kruskall.Dendritic_cells = kruskal.test(Dendritic_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Dendritic_cells = wilcox_test(ciberx[ciberx$cna_status != -1,], Dendritic_cells ~ cna_status)$p,
                                wil.del.Dendritic_cells = wilcox_test(ciberx[ciberx$cna_status != 1,], Dendritic_cells ~ cna_status)$p,
                                
                                kruskall.Macrophages = kruskal.test(Macrophages ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Macrophages = wilcox_test(ciberx[ciberx$cna_status != -1,], Macrophages ~ cna_status)$p,
                                wil.del.Macrophages = wilcox_test(ciberx[ciberx$cna_status != 1,], Macrophages ~ cna_status)$p,
                                
                                kruskall.B_cells = kruskal.test(B_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.B_cells = wilcox_test(ciberx[ciberx$cna_status != -1,], B_cells ~ cna_status)$p,
                                wil.del.B_cells = wilcox_test(ciberx[ciberx$cna_status != 1,], B_cells ~ cna_status)$p,
                                
                                kruskall.CD4_T_cells = kruskal.test(CD4_T_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.CD4_T_cells = wilcox_test(ciberx[ciberx$cna_status != -1,], CD4_T_cells ~ cna_status)$p,
                                wil.del.CD4_T_cells = wilcox_test(ciberx[ciberx$cna_status != 1,], CD4_T_cells ~ cna_status)$p,
                                
                                kruskall.Mast_cells = kruskal.test(Mast_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Mast_cells = wilcox_test(ciberx[ciberx$cna_status != -1,], Mast_cells ~ cna_status)$p,
                                wil.del.Mast_cells = wilcox_test(ciberx[ciberx$cna_status != 1,], Mast_cells ~ cna_status)$p,
                                
                                kruskall.NK_cells = kruskal.test(NK_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.NK_cells = wilcox_test(ciberx[ciberx$cna_status != -1,], NK_cells ~ cna_status)$p,
                                wil.del.NK_cells = wilcox_test(ciberx[ciberx$cna_status != 1,], NK_cells ~ cna_status)$p,
                                
                                kruskall.T_cells = kruskal.test(T_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.T_cells = wilcox_test(ciberx[ciberx$cna_status != -1,], T_cells ~ cna_status)$p,
                                wil.del.T_cells = wilcox_test(ciberx[ciberx$cna_status != 1,], T_cells ~ cna_status)$p,
                                
                                kruskall.Lymphocytes = kruskal.test(Lymphocytes ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Lymphocytes = wilcox_test(ciberx[ciberx$cna_status != -1,], Lymphocytes ~ cna_status)$p,
                                wil.del.Lymphocytes = wilcox_test(ciberx[ciberx$cna_status != 1,], Lymphocytes ~ cna_status)$p,
                                
                                kruskall.Non_plasma_B_cells = kruskal.test(Non_plasma_B_cells ~ as.factor(cna_status), data = ciberx)$p.value,
                                wil.ampl.Non_plasma_B_cells = wilcox_test(ciberx[ciberx$cna_status != -1,], Non_plasma_B_cells ~ cna_status)$p,
                                wil.del.Non_plasma_B_cells = wilcox_test(ciberx[ciberx$cna_status != 1,], Non_plasma_B_cells ~ cna_status)$p))
      }
    }
    
    return(stat_bin)
    
  })
  
  saveRDS(statistics_immune_score, file = paste0('/mnt/fabiogokce/fabiohd/US_intern/results_df/',
                                                 tt,'_1Mbp_ImmuneScore_CIBERSORTx.rds'))
}


ggplot(tpm_wide, aes(x = as.factor(cna_status), y = Immune_Score)) +
  geom_boxplot()
ggplot(ciberx, aes(x = as.factor(cna_status), y = Macrophages)) +
  geom_boxplot()
