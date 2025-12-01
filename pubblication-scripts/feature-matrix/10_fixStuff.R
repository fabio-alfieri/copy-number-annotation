# rm(list=setdiff(ls(), "ml_tables"))
# gc(full=T)

# ml_tables <- readRDS("/mnt/fabiogokce/fabiohd/ml_models/ml_tables_final.rds")
cytoband <- read.table(file = '/mnt/fabiogokce/fabiohd/ml_models/data/misc/cytoBand.txt',
                       header = T)


# Chromosome ----
variables_chr <- list()
for(tumor_type in names(ml_tables[['Chromosome']])){
  # tumor_type <- 'BRCA'
  for(chr in 1:22){
    # chr <- 1
    x <- ml_tables[['0.1Mbp']][[tumor_type]]
    x <- x[x$chr == chr,]
    mutations_raw <- sum(x$mutations_raw)
    mutations_norm <- sum(x$mutations_raw)/sum(x$total.Transcript.length, na.rm = T)
    n_OG <- sum(x$n_OG, na.rm = T)
    OGscore_pancancer <- n_OG/sum(x$n_genes, na.rm = T)
    OGscore.potency_pancancer <- sum(x$OGscore.potency_pancancer, na.rm = T)
    n_TSG <- sum(x$n_TSG, na.rm = T)
    TSGscore_pancancer <- n_TSG/sum(x$n_genes, na.rm = T)
    TSGscore.potency_pancancer <- sum(x$TSGscore.potency_pancancer, na.rm = T)
    n_ESS <- sum(x$n_ESS, na.rm = T)
    ESSscore_pancancer <- n_ESS/sum(x$n_genes, na.rm = T)
    n_HAPLO <- sum(x$n_HAPLO, na.rm = T)
    HAPLOscore_pancancer <- n_HAPLO/sum(x$n_genes, na.rm = T)
    
    variables_chr[[tumor_type]] <- as.data.frame(rbind(variables_chr[[tumor_type]], 
                                     cbind(chr = chr,
                                           mutations_raw, mutations_norm,
                                           n_OG, OGscore_pancancer, OGscore.potency_pancancer,
                                           n_TSG, TSGscore_pancancer, TSGscore.potency_pancancer,
                                           n_ESS, ESSscore_pancancer,
                                           n_HAPLO, HAPLOscore_pancancer)))
  }
  
  if(tumor_type == 'BRCA' | tumor_type == 'LUAD' | tumor_type == 'LUSC' | tumor_type == 'PAAD' | 
     tumor_type == 'COADREAD' | tumor_type == 'GBMLGG'){
    temp <- data.frame()
    for(chr in 1:22){
      # chr <- 1
      x <- ml_tables[['0.1Mbp']][[tumor_type]]
      x <- x[x$chr == chr,]
      
      n_OG_TS <- sum(x$n_OG, na.rm = T)
      OGscore_TS <- n_OG/sum(x$n_genes, na.rm = T)
      OGscore.potency_TS <- sum(x$OGscore.potency_TS, na.rm = T)
      n_TSG_TS <- sum(x$n_TSG, na.rm = T)
      TSGscore_TS <- n_TSG/sum(x$n_genes, na.rm = T)
      TSGscore.potency_TS <- sum(x$TSGscore.potency_TS, na.rm = T)
      n_ESS_TS <- sum(x$n_ESS, na.rm = T)
      ESSscore_TS <- n_ESS/sum(x$n_genes, na.rm = T)
      
      temp <- rbind(temp, 
                    cbind(chr = chr,
                          n_OG_TS, OGscore_TS, OGscore.potency_TS,
                          n_TSG_TS, TSGscore_TS, TSGscore.potency_TS,
                          n_ESS_TS, ESSscore_TS))
      
    }
    variables_chr[[tumor_type]] <- full_join(as.data.frame(variables_chr[[tumor_type]]), temp, by = 'chr')
  }
}


# Arm ----
variables_arm <- list()
for(tumor_type in names(ml_tables[['Arm']])){
  # tumor_type <- 'BRCA'
  for(arm in cytoband$arm){
    if(arm == '13p' | arm == '14p' | arm == '15p'){next}
    x <- ml_tables[['0.1Mbp']][[tumor_type]]
    x <- x[x$chr == parse_number(cytoband[cytoband$arm == arm,]$chr),]
    x <- x[x$start_bin >= cytoband[cytoband$arm == arm,]$start & x$start_bin <= cytoband[cytoband$arm == arm,]$end,]
    mutations_raw <- sum(x$mutations_raw)
    mutations_norm <- sum(x$mutations_raw)/sum(x$total.Transcript.length, na.rm = T)
    n_OG <- sum(x$n_OG, na.rm = T)
    OGscore_pancancer <- n_OG/sum(x$n_genes, na.rm = T)
    OGscore.potency_pancancer <- sum(x$OGscore.potency_pancancer, na.rm = T)
    n_TSG <- sum(x$n_TSG, na.rm = T)
    TSGscore_pancancer <- n_TSG/sum(x$n_genes, na.rm = T)
    TSGscore.potency_pancancer <- sum(x$TSGscore.potency_pancancer, na.rm = T)
    n_ESS <- sum(x$n_ESS, na.rm = T)
    ESSscore_pancancer <- n_ESS/sum(x$n_genes, na.rm = T)
    n_HAPLO <- sum(x$n_HAPLO, na.rm = T)
    HAPLOscore_pancancer <- n_HAPLO/sum(x$n_genes, na.rm = T)
    
    variables_arm[[tumor_type]] <- as.data.frame(rbind(variables_arm[[tumor_type]], 
                                                       cbind(arm = arm,
                                                             mutations_raw, mutations_norm,
                                                             n_OG, OGscore_pancancer, OGscore.potency_pancancer,
                                                             n_TSG, TSGscore_pancancer, TSGscore.potency_pancancer,
                                                             n_ESS, ESSscore_pancancer,
                                                             n_HAPLO, HAPLOscore_pancancer)))
  }
  
  if(tumor_type == 'BRCA' | tumor_type == 'LUAD' | tumor_type == 'LUSC' | tumor_type == 'PAAD' | 
     tumor_type == 'COADREAD' | tumor_type == 'GBMLGG'){
    temp <- data.frame()
    for(arm in cytoband$arm){
      if(arm == '13p' | arm == '14p' | arm == '15p'){next}
      x <- ml_tables[['0.1Mbp']][[tumor_type]]
      x <- x[x$chr == parse_number(cytoband[cytoband$arm == arm,]$chr),]
      x <- x[x$start_bin >= cytoband[cytoband$arm == arm,]$start & x$start_bin <= cytoband[cytoband$arm == arm,]$end,]
      
      n_OG_TS <- sum(x$n_OG, na.rm = T)
      OGscore_TS <- n_OG/sum(x$n_genes, na.rm = T)
      OGscore.potency_TS <- sum(x$OGscore.potency_TS, na.rm = T)
      n_TSG_TS <- sum(x$n_TSG, na.rm = T)
      TSGscore_TS <- n_TSG/sum(x$n_genes, na.rm = T)
      TSGscore.potency_TS <- sum(x$TSGscore.potency_TS, na.rm = T)
      n_ESS_TS <- sum(x$n_ESS, na.rm = T)
      ESSscore_TS <- n_ESS/sum(x$n_genes, na.rm = T)
      
      temp <- rbind(temp, 
                    cbind(arm = arm,
                          n_OG_TS, OGscore_TS, OGscore.potency_TS,
                          n_TSG_TS, TSGscore_TS, TSGscore.potency_TS,
                          n_ESS_TS, ESSscore_TS))
      
    }
    variables_arm[[tumor_type]] <- full_join(as.data.frame(variables_arm[[tumor_type]]), temp, by = 'arm')
  }
}

for(tt in names(ml_tables[['Arm']])){
  x <- ml_tables[['Arm']][[tt]]
  x <- x[,-c((ncol(x)-7):ncol(x))]
  colnames(variables_arm[[tt]])[1] <- 'bin'
  ml_tables[['Arm']][[tt]] <- full_join(x, variables_arm[[tt]], by = 'bin')
}

for(tt in names(ml_tables[['Chromosome']])){
  x <- ml_tables[['Chromosome']][[tt]]
  x <- x[,-c((ncol(x)-7):ncol(x))]
  colnames(variables_chr[[tt]])[1] <- 'chr'
  x$chr <- as.numeric(x$chr)
  variables_chr[[tt]][,1] <- as.numeric(variables_chr[[tt]][,1])
  ml_tables[['Chromosome']][[tt]] <- full_join(x, variables_chr[[tt]], by = 'chr')
}

saveRDS(ml_tables, file = '/mnt/fabiogokce/fabiohd/ml_models/ml_tables_240318.rds')
