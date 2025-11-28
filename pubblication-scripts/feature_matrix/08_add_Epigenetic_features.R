# rm(list=ls())
gc(full=T)

library(tidyverse)

load("/mnt/fabiogokce/fabiohd/ml_models/data/misc/ChromatinStates_all_levels_tissues_avgcounts&amp;scaled_avgcount.RData")
epi_features <- ChrStates.tissues
rm(ChrStates.tissues)

ml_merged <- readRDS("/mnt/fabiogokce/fabiohd/ml_models/data/ml_merged_tel.rds")

ml_merged_epi <- list()
for(i in names(ml_merged)){
  print(i)
  
  epi_merged_counts <- data.frame()
  for(n in names(epi_features[[i]])){
    epi_features[[i]][[n]]
    
    epi_features[[i]][[n]]$Counts$type <- n
    
    epi_features[[i]][[n]]$Counts <- cbind(
      as.data.frame(do.call(rbind,str_split(rownames(epi_features[[i]][[n]]$Counts),'_'))),
      epi_features[[i]][[n]]$Counts)
    colnames(epi_features[[i]][[n]]$Counts)[1:3] <- c('chr','start_bin','end_bin')
    epi_features[[i]][[n]]$Counts$chr <- parse_number(epi_features[[i]][[n]]$Counts$chr)
    
    epi_merged_counts <- rbind(epi_merged_counts, epi_features[[i]][[n]]$Counts)
  }
  
  epi_merged_counts[,2:3] <- apply(epi_merged_counts[,2:3], 2, as.numeric)
  
  ml_merged_epi[[i]] <- left_join( ml_merged[[i]], epi_merged_counts)
}

saveRDS(ml_merged_epi, "/mnt/fabiogokce/fabiohd/ml_models/data/ml_merged_tel_epi.rds")

# try to merge epigenetic features according to ChromatinState_Info table

ml_merged_epi2 <- list()
for(i in names(ml_merged_epi)){
  active_tss <- ml_merged_epi[[i]][,c('E1')]
  promoter <- apply(ml_merged_epi[[i]][,c('E2','E3','E4')], 1, sum)
  transcribed <- apply(ml_merged_epi[[i]][,c('E5','E6','E7')], 1, sum)
  weak_trascribed <- ml_merged_epi[[i]][,c('E8')]
  transcribed_regulatory_regions <- apply(ml_merged_epi[[i]][,c('E9','E10','E11','E12')], 1, sum)
  active_enhancer <- apply(ml_merged_epi[[i]][,c('E13','E14','E15')], 1, sum)
  weak_enhancer <- apply(ml_merged_epi[[i]][,c('E16','E17','E18')], 1, sum)
  primary_DNAase <- ml_merged_epi[[i]][,c('E19')]
  heterochromatin <- ml_merged_epi[[i]][,c('E21')]
  poised_promoter <- ml_merged_epi[[i]][,c('E22')]
  bivalent_promoter <- ml_merged_epi[[i]][,c('E23')]
  repressed_polycomb <- ml_merged_epi[[i]][,c('E24')]
  quiescent <- ml_merged_epi[[i]][,c('E25')]
  
  
  ml_merged_epi2[[i]] <- cbind(ml_merged[[i]], as.data.frame(cbind(active_tss,promoter,transcribed,weak_trascribed,transcribed_regulatory_regions,
        active_enhancer,weak_enhancer,primary_DNAase,heterochromatin,poised_promoter,bivalent_promoter,
        repressed_polycomb,quiescent)))
}

saveRDS(ml_merged_epi2, "/mnt/fabiogokce/fabiohd/ml_models/data/ml_merged_tel_epi2.rds")


ml_merged_epi3 <- list()
for(i in names(ml_merged_epi)){
  promoters_transcribed <- apply(ml_merged_epi[[i]][,c('E1','E2','E3','E4','E5','E6','E7',
                                                       'E8','E9','E10','E11','E12')], 1, sum)
  enhancer <- apply(ml_merged_epi[[i]][,c('E13','E14','E15','E16','E17','E18','E19')], 1, sum)
  other_promoters <- apply(ml_merged_epi[[i]][,c('E22','E23','E24')], 1, sum)
  hetero_quiesc <- apply(ml_merged_epi[[i]][,c('E20','E21','E25')], 1, sum)
  
  ml_merged_epi3[[i]] <- cbind(ml_merged[[i]], as.data.frame(cbind(promoters_transcribed,enhancer,
                                                                   other_promoters,hetero_quiesc)))
}

saveRDS(ml_merged_epi3, "/mnt/fabiogokce/fabiohd/ml_models/data/ml_merged_tel_epi3.rds")

ml_merged_epi4 <- list()
for(i in names(ml_merged_epi)){
  promoters_transcribed <- apply(ml_merged_epi[[i]][,c('E1','E2','E3','E4','E5','E6','E7',
                                                       'E8','E9','E10','E11','E12')], 1, sum)
  enhancer <- apply(ml_merged_epi[[i]][,c('E13','E14','E15','E16','E17')], 1, sum)
  dnase <- apply(ml_merged_epi[[i]][,c('E18','E19')], 1, sum)
  
  other_promoters <- apply(ml_merged_epi[[i]][,c('E22','E23','E24')], 1, sum)
  hetero_quiesc <- apply(ml_merged_epi[[i]][,c('E20','E21','E25')], 1, sum)
  
  ml_merged_epi4[[i]] <- cbind(ml_merged[[i]], as.data.frame(cbind(promoters_transcribed,enhancer,dnase,
                                                                   other_promoters,hetero_quiesc)))
}

saveRDS(ml_merged_epi4, "/mnt/fabiogokce/fabiohd/ml_models/data/ml_merged_tel_epi4.rds")

ml_merged_epi5 <- list()
for(i in names(ml_merged_epi)){
  promoters_transcribed <- apply(ml_merged_epi[[i]][,c('E1','E2','E3','E4','E5','E6','E7',
                                                       'E8','E9','E10','E11','E12')], 1, sum)
  enhancer <- apply(ml_merged_epi[[i]][,c('E13','E14','E15','E16','E17')], 1, sum)
  dnase <- apply(ml_merged_epi[[i]][,c('E18','E19')], 1, sum)
  
  other_promoters <- apply(ml_merged_epi[[i]][,c('E22','E23','E24')], 1, sum)
  znf <- ml_merged_epi[[i]][,'E20']
  heterochr <- ml_merged_epi[[i]][,'E21']
  quiesc <- ml_merged_epi[[i]][,'E25']
  
  ml_merged_epi5[[i]] <- cbind(ml_merged[[i]], as.data.frame(cbind(promoters_transcribed,enhancer,dnase,
                                                                   other_promoters,znf,heterochr,quiesc)))
}

saveRDS(ml_merged_epi5, "/mnt/fabiogokce/fabiohd/ml_models/data/ml_merged_tel_epi5.rds")

