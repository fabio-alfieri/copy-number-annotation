
#Preparing ML tables with all features together - Part I (Non tissue-specific features)

# Non-tissue specific features (from GS)
# 1. Ohnologs and paralogs
# 2. Ensembl and CORUM related features
# 3. Distance to fragile sites

rm(list = ls())
gc()

# Set working directory 

setwd("/mnt/fabiogokce/gokce")

#Required libraries
library(dplyr)

#1.Non Tissue-specific features

#Most of the non Tissue-specific features are gene based. Therefore, we need genes located on the corresponding bin.
#Genes per bins for focal CNA
load("./Data/All_levels_genes_per_bins.RData")

#All backbone tables for focal CNAs and aneuploidies
load("./Data/All_levels_backbonetables.RData")

#Features: 2. Ensembl and CORUM related features
load("./Data/Ensembl_CORUM_features_tissue_general.RData")

levels <- names(chr_backbone_namesfixed)
backbone_tables_w_features <- list()

for(level in levels){
  backbone <- do.call(rbind,chr_backbone_namesfixed[[level]])
  if(!level %in% c("Arm","Chromosome")){backbone$bin <- paste(backbone$chr,backbone$bin,sep = "_")} #To match the bin names between backbone and genes_on_bins tables
  features.df <- c()
  for(bin in backbone$bin){
    bin.genes <- genes_on_bins[[level]][[bin]]
    #Some bins do not have genes, and those are not included in genes_on_bins data
    if(length(bin.genes) == 0){row <- c(bin,"",length(bin.genes),rep(NA,29))}
    else{
      df.genes <- gene.table[gene.table$Gene.name %in% bin.genes,] #Gene table only has protein coding genes!
      
      #Find complex partners in trans
      partners <- unique(unlist(lapply(df.genes$Partners, function(x) strsplit(x,"::")[[1]])))
      partners <- partners[complete.cases(partners)]
      partners.trans <- setdiff(partners,bin.genes)
      
      #Find interacting proteins in trans (Interactome INSIDER)
      ppis_IntINSIDER <- unique(unlist(lapply(df.genes$PPIs_IntINSIDER, function(x) strsplit(x,"::")[[1]])))
      ppis_IntINSIDER <- ppis_IntINSIDER[complete.cases(ppis_IntINSIDER)]
      ppis.trans_IntINSIDER <- setdiff(ppis_IntINSIDER,bin.genes)
      
      #Find interacting proteins in trans based on HIPPIE (no-cutoff)
      ppis_HIPPIE <- unique(unlist(lapply(df.genes$PPIs_HIPPIE, function(x) strsplit(x,"::")[[1]])))
      ppis_HIPPIE <- ppis_HIPPIE[complete.cases(ppis_HIPPIE)]
      ppis.trans_HIPPIE <- setdiff(ppis_HIPPIE,bin.genes)
      #Find interacting proteins in trans based on HIPPIE (.63)
      ppis_HIPPIE_063 <- unique(unlist(lapply(df.genes$PPIs_063_HIPPIE, function(x) strsplit(x,"::")[[1]])))
      ppis_HIPPIE_063 <- ppis_HIPPIE_063[complete.cases(ppis_HIPPIE_063)]
      ppis.trans_HIPPIE_063 <- setdiff(ppis_HIPPIE_063,bin.genes)
      #Find interacting proteins in trans based on HIPPIE (.73)
      ppis_HIPPIE_073 <- unique(unlist(lapply(df.genes$PPIs_073_HIPPIE, function(x) strsplit(x,"::")[[1]])))
      ppis_HIPPIE_073 <- ppis_HIPPIE_073[complete.cases(ppis_HIPPIE_073)]
      ppis.trans_HIPPIE_073 <- setdiff(ppis_HIPPIE_073,bin.genes)
      
      
      #Create bin row with features
      row <- c(bin,
               paste(bin.genes,collapse = "::"), #Genes on the bin
               length(bin.genes), #Number of genes on the bin
               #Gene length
               mean(df.genes$Gene.length,na.rm = T),
               sum(df.genes$Gene.length,na.rm = T),
               #Gene transcript length
               mean(df.genes$Transcript.length..including.UTRs.and.CDS.,na.rm = T),
               sum(df.genes$Transcript.length..including.UTRs.and.CDS.,na.rm = T),
               #Coding region length
               mean(df.genes$cds_length,na.rm = T),
               sum(df.genes$cds_length,na.rm = T),
               #GC content
               mean(df.genes$Gene...GC.content,na.rm = T),
               sum(df.genes$Gene...GC.content,na.rm = T),
               #Number of protein complex subunit on the bin
               length(which(df.genes$ComplexSubunit == "True")),
               #Partners and number of partners
               paste(partners, collapse = "::"),length(partners), 
               #Partners in trans and number of partners in trans
               paste(partners.trans, collapse = "::"),length(partners.trans),
               #Interacting proteins and number of them (Interactome INSIDER)
               paste(ppis_IntINSIDER, collapse = "::"),length(ppis_IntINSIDER),
               #Interacting proteins in trans and number of them (Interactome INSIDER)
               paste(ppis.trans_IntINSIDER, collapse = "::"), length(ppis.trans_IntINSIDER),
               
               # Interaction from HIPPIE
               #Interacting proteins and number of them (no cut-off)
               paste(ppis_HIPPIE, collapse = "::"),length(ppis_HIPPIE),
               #Interacting proteins in trans and number of them (no cut-off)
               paste(ppis.trans_HIPPIE, collapse = "::"), length(ppis.trans_HIPPIE),
               #Interacting proteins and number of them (.63)
               paste(ppis_HIPPIE_063, collapse = "::"),length(ppis_HIPPIE_063),
               #Interacting proteins in trans and number of them (.63)
               paste(ppis.trans_HIPPIE_063, collapse = "::"), length(ppis.trans_HIPPIE_063),
               #Interacting proteins and number of them (.73)
               paste(ppis_HIPPIE_073, collapse = "::"),length(ppis_HIPPIE_073),
               #Interacting proteins in trans and number of them (.73)
               paste(ppis.trans_HIPPIE_073, collapse = "::"), length(ppis.trans_HIPPIE_073))
      
      }
    features.df <- rbind(features.df,row)}
  features.df <- as.data.frame(features.df)
  colnames(features.df) <- c("bin",
                             "genes",
                             "n_genes",
                             "mean.Gene.lenght","total.Gene.lenght",
                             "mean.Transcript.lenght","total.Transcript.length",
                             "mean.CDS.lenght","total.CDS.lenght",
                             "mean.GC.content","total.GC.content",
                             "total.complex.subunit",
                             "Partners","total_n_partners",
                             "Partners.trans","total_n_partners.trans",
                             "PPIs_IntINSIDER","total_n_PPIs_IntINSIDER",
                             "PPIs.trans_IntINSIDER","total_n_PPIs.trans_IntINSIDER",
                             "PPIs_HIPPIE","total_n_PPIs_HIPPIE",
                             "PPIs.trans_HIPPIE","total_n_PPIs.trans_HIPPIE",
                             "PPIs_HIPPIE_063","total_n_PPIs_HIPPIE_063",
                             "PPIs.trans_HIPPIE_063","total_n_PPIs.trans_HIPPIE_063",
                             "PPIs_HIPPIE_073","total_n_PPIs_HIPPIE_073",
                             "PPIs.trans_HIPPIE_073","total_n_PPIs.trans_HIPPIE_073")
  features.df[,c(3:12,14,16,18,20,22,24,26,28,30,32)] <- sapply(features.df[,c(3:12,14,16,18,20,22,24,26,28,30,32)], as.numeric)
  backbone <- merge(backbone,features.df, by="bin")
  backbone_tables_w_features[[level]] <- backbone}

save(backbone_tables_w_features, file = "./Data/Backbone_tables_with_non_tissue_specific_features.RData")

#Features: 1. Ohnologs and Paralogs

rm(list = ls())
gc()

# Ohnologs
load("./Data/Backbone_tables_with_non_tissue_specific_features.RData")
load("./Data/OhnologList.RData")
backbone_tables_w_features_added <- list()
for(level in names(backbone_tables_w_features)){
  df <- backbone_tables_w_features[[level]]
  features.df <- c() 
  for(bin in df$bin){
    genes <- strsplit(df[df$bin == bin,"genes"],"::")[[1]]
    #Some bins do not have genes, and those are not included in genes_on_bins data
    if(length(genes) == 0){row <- c(bin,rep(NA,5))}
    else{
      #Ohnologs
      #1.Makino and McLysaght paper
      ohnolog <- as.data.frame(ohnologs[["mmpaper"]])
      df.ohnolog <- ohnolog[ohnolog$Symbol %in% genes,]
      n_ohnolog_bin <- nrow(df.ohnolog) #Number of ohnologs in bin
      #All ohnologs
      if(n_ohnolog_bin != 0){
        #Ohnologs located on outside of the bin
        df.ohnolog$n_ohnolog.mmpaper_trans <- unlist(lapply(df.ohnolog$ohnologs_mmpaper, function(x)
          length(setdiff(strsplit(x,"::")[[1]],genes)))) #Find ohnologs of genes other than genes on the bin
        mean_n_ohnologs.mmpaper <- mean(df.ohnolog$n_ohnolog_mmpaper,na.rm = TRUE)
        total_n_ohnologs.mmpaper <- sum(df.ohnolog$n_ohnolog_mmpaper,na.rm = TRUE)
        mean_n_ohnologs.mmpaper_trans <- mean(df.ohnolog$n_ohnolog.mmpaper_trans,na.rm = TRUE)
        total_n_ohnologs.mmpaper_trans <- sum(df.ohnolog$n_ohnolog.mmpaper_trans,na.rm = TRUE)}
      else{
        mean_n_ohnologs.mmpaper <- NA
        total_n_ohnologs.mmpaper <- NA
        mean_n_ohnologs.mmpaper_trans <- NA
        total_n_ohnologs.mmpaper_trans <- NA}
      #All information together for the corresponding bin
      row <- c(bin,n_ohnolog_bin,
               mean_n_ohnologs.mmpaper, total_n_ohnologs.mmpaper, mean_n_ohnologs.mmpaper_trans, total_n_ohnologs.mmpaper_trans)}
    features.df <- rbind(features.df,row)}
  features.df <- as.data.frame(features.df)
  colnames(features.df) <- c("bin","n_ohnolog_bin",
                             "mean_n_ohnologs.mmpaper", "total_n_ohnologs.mmpaper", "mean_n_ohnologs.mmpaper_trans", "total_n_ohnologs.mmpaper_trans")
  features.df[,2:6] <- sapply(features.df[,2:6], as.numeric)
  df <- merge(df,features.df, by="bin")
  backbone_tables_w_features_added[[level]] <- df}
save(backbone_tables_w_features_added, file = "./Data/Backbone_tables_with_non_tissue_specific_features.RData")

# Paralogs

rm(list = ls())
load("./Data/Backbone_tables_with_non_tissue_specific_features.RData")
load("./Data/ParalogList.RData")
paralogs$n_paralogs <- as.numeric(paralogs$n_paralogs)
backbone_tables_w_features_added_new <- list()
for(level in names(backbone_tables_w_features_added)){
  df <- backbone_tables_w_features_added[[level]]
  features.df <- c() 
  for(bin in df$bin){
    genes <- strsplit(df[df$bin == bin,"genes"],"::")[[1]]
    #Some bins do not have genes, and those are not included in genes_on_bins data
    if(length(genes) == 0){row <- c(bin,rep(NA,5))}
    else{
      #Paralogs
      df.paralog <- paralogs[paralogs$symbol %in% genes,]
      n_paralog_bin <- nrow(df.paralog)
      #All ohnologs
      if(n_paralog_bin != 0){
        #Paralogs located on outside of the bin
        df.paralog$n_paralog_trans <- unlist(lapply(df.paralog$paralogs, function(x)
          length(setdiff(strsplit(x,"::")[[1]],genes))))
        mean_n_paralogs <- mean(df.paralog$n_paralogs,na.rm = TRUE)
        total_n_paralogs <- sum(df.paralog$n_paralogs,na.rm = TRUE)
        mean_n_paralogs_trans <- mean(df.paralog$n_paralog_trans,na.rm = TRUE)
        total_n_paralogs_trans <- sum(df.paralog$n_paralog_trans,na.rm = TRUE)}
      else{
        mean_n_paralogs <- NA
        total_n_paralogs <- NA
        mean_n_paralogs_trans <- NA
        total_n_paralogs_trans <- NA}
      #All information together for the corresponding bin
      row <- c(bin,n_paralog_bin,
               mean_n_paralogs, total_n_paralogs, mean_n_paralogs_trans, total_n_paralogs_trans)}
    features.df <- rbind(features.df,row)}
  features.df <- as.data.frame(features.df)
  colnames(features.df) <- c("bin","n_paralogs_bin",
                             "mean_n_paralogs", "total_n_paralogs", "mean_n_paralogs_trans", "total_n_paralogs_trans")
  features.df[,2:6] <- sapply(features.df[,2:6], as.numeric)
  df <- merge(df,features.df, by="bin")
  backbone_tables_w_features_added_new[[level]] <- df}
save(backbone_tables_w_features_added_new, file = "./Data/Backbone_tables_with_non_tissue_specific_features.RData")

#Features: 3. Distance to fragile sites

rm(list = ls())
load("./Data/Backbone_tables_with_non_tissue_specific_features.RData")
load("./Data/Distance_to_closest_fragile_site.RData")
load("./Data/FeatureOptimization_Fragile_sites_weighted_distance_scores.RData")

backbone_tables_w_features_non_tissue_specific <- list()
for(level in names(backbone_tables_w_features_added_new)){
  df <- backbone_tables_w_features_added_new[[level]]
  if(level %in% c("Arm","Chromosome")){backbone_tables_w_features_non_tissue_specific[[level]] <- df}
  else{
    fragile.sites <- Dist.fragile.sites[[level]]
    fragile.sites.weighted <- Weighted.dist.scores[[level]]
    fragile.sites <- merge(fragile.sites,fragile.sites.weighted,by="bin")
    df <- merge(df,fragile.sites[,c("bin","dist.to.closest.FGS","Dist.FGS.weighted")], by = "bin", all.x = "TRUE")
    df$dist.to.closest.FGS.scaled <- as.numeric(scale(df$dist.to.closest.FGS))
    backbone_tables_w_features_non_tissue_specific[[level]] <- df}
}

save(backbone_tables_w_features_non_tissue_specific,
     file =  "./Data/Backbone_tables_with_non_tissue_specific_features.RData") 
