# Tumour suppressors and oncogenes

rm(list = ls())

packages <- c(
  "dplyr", "openxlsx", "ggplot2", "parallel"
)

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

wd <- 'path/to/GitHub/copy-number-annotation/'
setwd(wd)

cancer_drivers_path <- "feature_matrix_data/2023-05-31_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"
cohorts_path <- "feature_matrix_data/2023-05-31_IntOGen-Cohorts/cohorts.tsv"
tcga_cohort_path <- "plots/Number_drivers_TCGA_cohorts.pdf"
tcga_cohort_xlsx <- "feature_matrix_data/IntOGen_tcga_cohorts.xlsx"
ensembl_path <- "feature_matrix_data/hg19_14022024.txt"
ccds_path <- "feature_matrix_data/CCDS.current_hg19.txt"
drivers_list_path <- "feature_matrix_data/IntOGen_cancer_drivers_tissue-specific.RData"
genes_per_bins_path <- "feature_matrix_data/All_levels_genes_per_bins.RData"
density_drivers_path <- "feature_matrix_data/Density_cancer-drivers_intogen.RData"
backbone_path <- "feature_matrix_data/All_levels_backbonetables.RData"
distance_weighted_path <- "feature_matrix_data/Distance_to_cancer-drivers_intogen_weighted.RData"

# Cancer drivers from intOGen (release date 2023.05.31)
cancer.drivers <- read.delim(cancer_drivers_path)

# 1. Some info about the cohorts and drivers
unique(cancer.drivers$CANCER_TYPE)
unique(cancer.drivers$ROLE) # "Act", "LoF", "ambiguous"
info <- cancer.drivers %>% group_by(COHORT) %>% summarise(n = n())

cohorts <- read.delim(cohorts_path) 
setdiff(cohorts$COHORT,info$COHORT)
diff <- cohorts[cohorts$COHORT %in% setdiff(cohorts$COHORT,info$COHORT),] # Note: there are 12 cohorts that are not in the drivers list because there is no driver detected in those cohorts
merged.df <- merge(info, cohorts, by = "COHORT")
merged.df <- merged.df[,c(1,2,15,3:14)]
cor(merged.df$n,merged.df$DRIVERS)

# TCGA cohorts plot
cohorts.TCGA <- cohorts[cohorts$SOURCE == "TCGA/PanCatAtlas, phs000178",]
cohorts.TCGA <- cohorts.TCGA[order(cohorts.TCGA$DRIVERS,decreasing = T),]
#write.xlsx(cohorts.TCGA, file = "./Data/IntOGen_tcga_cohorts.xlsx") #Mapping column was added manually
cancer.drivers <- cancer.drivers[cancer.drivers$COHORT %in% cohorts.TCGA$COHORT,]
info <- cancer.drivers %>% group_by(CANCER_TYPE,ROLE) %>% summarise(n = n())
info$CANCER_TYPE <- factor(info$CANCER_TYPE, levels = cohorts.TCGA$CANCER)

p <- ggplot(info, aes(x=CANCER_TYPE, y=n, fill=ROLE)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Number of genes")

pdf(tcga_cohort_path,width = 9, height = 4)
p
dev.off()

# 2. List of tissue-specific tumour suppressors (LoF) and oncogenes (Act)
rm(list = ls())

# 2.1. Adding gene coordinates
#Data from Ensembl Biomart (downloaded on 14 February 2024 - Genome assembly GRCh37.p13)
#List of genes with the largest transcript (autosomal protein coding genes)
hg19 <- read.delim(ensembl_path) # 215404
genes <- unique(hg19$Gene.name)
hg19 <- hg19[hg19$Gene.type == "protein_coding",] # 160238
hg19 <- hg19[hg19$Chromosome.scaffold.name %in% as.character(seq(1,22,1)),] # 141152 
hg19 <- hg19[order(hg19$Transcript.length..including.UTRs.and.CDS.,decreasing = T),]
hg19 <- hg19[!duplicated(hg19$Gene.name),] #19349 - The unique gene set with the longest transcript
hg19$Gene.length <- hg19$Gene.end..bp. - hg19$Gene.start..bp.
hg19 <- hg19[,c(1,9,3,4,5,8,10,11,14,15)]

#Data from CCDC - Keep the one with the longest CCDS
ccds <- read.delim(ccds_path) #27809 genes
ccds <- ccds[ccds$chromosome %in% as.character(seq(1,22,1)),] #26473 genes
ccds$CCDS.ID <- gsub("\\..*","",ccds$ccds_id)
ccds <- ccds[ccds$ccds_status == "Public",]
ccds[,c(8,9)] <- sapply(ccds[,c(8,9)],as.numeric)
#Keep the gene with the longest coding region
ccds$cds_length <- ccds$cds_to - ccds$cds_from
ccds <- ccds[order(ccds$cds_length, decreasing = T),]
ccds <- ccds[!duplicated(ccds$gene),] # 17559 - the longest cdc
ccds <- ccds[,c(1:5,8,9,12,13)]
#Merging two tables
gene.table <- merge(hg19,ccds,by.x="Gene.name",by.y = "gene",all = TRUE) #19621 genes (union)
gene.table <- gene.table[,c(1,3,4,5)]
colnames(gene.table) <- c("SYMBOL","chrom","start","end")

rm(list = setdiff(ls(),"gene.table"))

cancer.drivers <- read.delim(cancer_drivers_path) #4360
cancer.drivers <- merge(cancer.drivers,gene.table, by = "SYMBOL") #4048
# 574 drivers out of 619 left after merging. 

# 2.2. Tissue-specific driver lists
mapping <- read.xlsx(tcga_cohort_xlsx)
cancer.drivers <- cancer.drivers[cancer.drivers$COHORT %in% mapping$COHORT,]

cancer.types <- unique(mapping$Mapping) # Based on the naming in our previous data

drivers.list <- list()
for(type in cancer.types){
  CANCER_TYPE <- mapping[mapping$Mapping == type, "CANCER"]
  m <- cancer.drivers[cancer.drivers$CANCER_TYPE %in% CANCER_TYPE,]
  cohort.l <- list()
  for(role in unique(m$ROLE)){
    genes <- m[m$ROLE == role, c("SYMBOL","chrom","start","end","QVALUE_COMBINATION")]
    genes <- genes[!duplicated(genes$SYMBOL),]
    cohort.l[[role]] <- genes
  }
  drivers.list[[type]] <- cohort.l
}

save(drivers.list, file = drivers_list_path)

# 3. Calculating features
rm(list = ls())

# 3.1. Density of cancer drivers (Number of genes) -  Weighted scores
load(drivers_list_path) # Tissue-specific cancer drivers final list
load(genes_per_bins_path) # Genes on the bins

# Function for density calculation 

weighted.density <- function(genes){
  ogs.bin <- ogs[ogs$SYMBOL %in% genes,]
  if(nrow(ogs.bin) == 0){
    density.ogs <- 0
    Weighted.density.ogs <- 0}
  else{
    density.ogs <- nrow(ogs.bin)
    Weighted.density.ogs <- mean(ogs.bin$weight)}
  tsg.bin <- tsg[tsg$SYMBOL %in% genes,]
  if(nrow(tsg.bin) == 0){
    density.tsg <- 0
    Weighted.density.tsg <- 0}
  else{
    density.tsg <- nrow(tsg.bin)
    Weighted.density.tsg <- mean(tsg.bin$weight)}
  result <- c(density.ogs, Weighted.density.ogs,density.tsg,Weighted.density.tsg)
  return(result)
}

Density.drivers <- list()
for(level in names(genes_on_bins)){
  m <- genes_on_bins[[level]]
  for(tissue in names(drivers.list)){
    ogs <- drivers.list[[tissue]]$Act
    ogs$weight <- -log10(ogs$QVALUE_COMBINATION)
    tsg <- drivers.list[[tissue]]$LoF
    tsg$weight <- -log10(tsg$QVALUE_COMBINATION)
    res.df <- mclapply(m, weighted.density, mc.cores = 20)
    final.df <- do.call(rbind,res.df)
    final.df <- as.data.frame(final.df)
    colnames(final.df) <- c("density.OG", "Weighted.density.OG","density.TSG","Weighted.density.TSG")
    final.df$bin <- rownames(final.df)
    Density.drivers[[level]][[tissue]] <- final.df
  }
  print(level)
}

save(Density.drivers,
     file = density_drivers_path)

# 3.2. Distance to closest TSG/OG
load(drivers_list_path)
load(backbone_path)

# Oncogenes
Dist.oncogenes <- list()
for(cohort in names(drivers.list)){
  driver.sites <- drivers.list[[cohort]]$Act
  driver.sites$chrom <- paste0("chr",driver.sites$chrom)
  for(level in names(chr_backbone_namesfixed)){
    if(level %in% c("Arm","Chromosome"))next
    coordinates <- do.call(rbind,chr_backbone_namesfixed[[level]])
    coordinates$bin <- paste(coordinates$chr,coordinates$bin,sep = "_")
    coordinates[,3:4] <- sapply(coordinates[,3:4], as.numeric)
    df <- c()
    for(i in 1:nrow(coordinates)){
      chr <- paste0("chr",as.character(coordinates$chr)[i])
      start <- coordinates$start_bin[i]
      end <- coordinates$end_bin[i]
      #Distance to the closest oncogene
      driver.sites.chr <- driver.sites[driver.sites$chrom == chr,]
      if(nrow(driver.sites.chr) == 0)next
      df.bin <- c()
      for(j in 1:nrow(driver.sites.chr)){
        fs.start <- driver.sites.chr$start[j]
        fs.end <- driver.sites.chr$end[j]
        dist.to.driver <- ifelse(fs.end < start, start - fs.end, #If driver is upstream of the bin
                                      ifelse(fs.start > end, fs.start - end, 0 #If driver is downstream of the bin
                                      ))
        row <- c(chr,coordinates$bin[i],start,end,fs.start,fs.end,dist.to.driver)
        df.bin <- rbind(df.bin,row)}
      #Find the closest bin
      df.bin <- as.data.frame(df.bin)
      df.bin$V7 <- as.numeric(df.bin$V7)
      df.bin <- df.bin[order(df.bin$V7, decreasing = F),]
      df <- rbind(df,df.bin[1,])}
    df <- as.data.frame(df)
    colnames(df) <- c("chr","bin","start_bin","end_bin","start_OG","end_OG","dist.to.closest.OG")
    df[,3:7] <- sapply(df[,3:7], as.numeric)
    Dist.oncogenes[[level]][[cohort]] <- df
  }
  print(cohort)
  
}

# Tumour suppressors
Dist.tumoursupp <- list()
for(cohort in names(drivers.list)){
  driver.sites <- drivers.list[[cohort]]$LoF
  driver.sites$chrom <- paste0("chr",driver.sites$chrom)
  for(level in names(chr_backbone_namesfixed)){
    if(level %in% c("Arm","Chromosome"))next
    coordinates <- do.call(rbind,chr_backbone_namesfixed[[level]])
    coordinates$bin <- paste(coordinates$chr,coordinates$bin,sep = "_")
    coordinates[,3:4] <- sapply(coordinates[,3:4], as.numeric)
    df <- c()
    for(i in 1:nrow(coordinates)){
      chr <- paste0("chr",as.character(coordinates$chr)[i])
      start <- coordinates$start_bin[i]
      end <- coordinates$end_bin[i]
      #Distance to the closest oncogene
      driver.sites.chr <- driver.sites[driver.sites$chrom == chr,]
      if(nrow(driver.sites.chr) == 0)next
      df.bin <- c()
      for(j in 1:nrow(driver.sites.chr)){
        fs.start <- driver.sites.chr$start[j]
        fs.end <- driver.sites.chr$end[j]
        dist.to.driver <- ifelse(fs.end < start, start - fs.end, #If driver is upstream of the bin
                                 ifelse(fs.start > end, fs.start - end, 0 #If driver is downstream of the bin
                                 ))
        row <- c(chr,coordinates$bin[i],start,end,fs.start,fs.end,dist.to.driver)
        df.bin <- rbind(df.bin,row)}
      #Find the closest bin
      df.bin <- as.data.frame(df.bin)
      df.bin$V7 <- as.numeric(df.bin$V7)
      df.bin <- df.bin[order(df.bin$V7, decreasing = F),]
      df <- rbind(df,df.bin[1,])}
    df <- as.data.frame(df)
    colnames(df) <- c("chr","bin","start_bin","end_bin","start_TSG","end_TSG","dist.to.closest.TSG")
    df[,3:7] <- sapply(df[,3:7], as.numeric)
    Dist.tumoursupp[[level]][[cohort]] <- df
  }
  print(cohort)
  
}

save(Dist.oncogenes,Dist.tumoursupp,
     file = "./Data/Distance_to_closest_cancer-drivers_intogen.RData")

# 3.3. Weighted distance score

rm(list = ls())

load(drivers_list_path)
load(backbone_path)

# Oncogenes
Dist.oncogenes.weighted <- list()
for(cohort in names(drivers.list)){
  driver.sites <- drivers.list[[cohort]]$Act
  driver.sites$chrom <- paste0("chr",driver.sites$chrom)
  for(level in names(chr_backbone_namesfixed)){
    if(level %in% c("Arm","Chromosome"))next
    coordinates <- do.call(rbind,chr_backbone_namesfixed[[level]])
    coordinates$bin <- paste(coordinates$chr,coordinates$bin,sep = "_")
    coordinates[,3:4] <- sapply(coordinates[,3:4], as.numeric)
    df <- c()
    for(i in 1:nrow(coordinates)){
      chr <- paste0("chr",as.character(coordinates$chr)[i])
      start <- coordinates$start_bin[i]
      end <- coordinates$end_bin[i]
      #Distance to the closest oncogene
      driver.sites.chr <- driver.sites[driver.sites$chrom == chr,]
      if(nrow(driver.sites.chr) == 0)next
      df.bin <- c()
      for(j in 1:nrow(driver.sites.chr)){
        fs.start <- driver.sites.chr$start[j]
        fs.end <- driver.sites.chr$end[j]
        dist.to.driver <- ifelse(fs.end < start, start - fs.end, #If driver is upstream of the bin
                                 ifelse(fs.start > end, fs.start - end, 0 #If driver is downstream of the bin
                                 ))
        row <- c(chr,coordinates$bin[i],start,end,fs.start,fs.end,dist.to.driver)
        df.bin <- rbind(df.bin,row)}
      df <- rbind(df,df.bin)}
    df <- as.data.frame(df)
    colnames(df) <- c("chr","bin","start_bin","end_bin","start_OG","end_OG","dist.to.OG")
    df[,3:7] <- sapply(df[,3:7], as.numeric)
    # Weighted scores
    #Giving weight based on the distance
    df$weight <- ifelse(df$dist.to.OG > 10000000, 0,
                        ifelse(df$dist.to.OG <= 10000000 & df$dist.to.OG > 8000000, 2,
                               ifelse(df$dist.to.OG <= 8000000 & df$dist.to.OG > 6000000, 4,
                                      ifelse(df$dist.to.OG <= 6000000 & df$dist.to.OG > 4000000, 6,
                                             ifelse(df$dist.to.OG <= 4000000 & df$dist.to.OG > 2000000, 8,10)))))
    df.weight <- df %>% group_by(bin) %>% summarise(Dist.OG.weighted = sum(weight))
    
    Dist.oncogenes.weighted[[level]][[cohort]] <- df.weight
  }
  print(cohort)
  
}


# Tumour suppressors
Dist.tumoursupp.weighted <- list()
for(cohort in names(drivers.list)){
  driver.sites <- drivers.list[[cohort]]$LoF
  driver.sites$chrom <- paste0("chr",driver.sites$chrom)
  for(level in names(chr_backbone_namesfixed)){
    if(level %in% c("Arm","Chromosome"))next
    coordinates <- do.call(rbind,chr_backbone_namesfixed[[level]])
    coordinates$bin <- paste(coordinates$chr,coordinates$bin,sep = "_")
    coordinates[,3:4] <- sapply(coordinates[,3:4], as.numeric)
    df <- c()
    for(i in 1:nrow(coordinates)){
      chr <- paste0("chr",as.character(coordinates$chr)[i])
      start <- coordinates$start_bin[i]
      end <- coordinates$end_bin[i]
      #Distance to the closest oncogene
      driver.sites.chr <- driver.sites[driver.sites$chrom == chr,]
      if(nrow(driver.sites.chr) == 0)next
      df.bin <- c()
      for(j in 1:nrow(driver.sites.chr)){
        fs.start <- driver.sites.chr$start[j]
        fs.end <- driver.sites.chr$end[j]
        dist.to.driver <- ifelse(fs.end < start, start - fs.end, #If driver is upstream of the bin
                                 ifelse(fs.start > end, fs.start - end, 0 #If driver is downstream of the bin
                                 ))
        row <- c(chr,coordinates$bin[i],start,end,fs.start,fs.end,dist.to.driver)
        df.bin <- rbind(df.bin,row)}
      df <- rbind(df,df.bin)}
    df <- as.data.frame(df)
    colnames(df) <- c("chr","bin","start_bin","end_bin","start_TSG","end_TSG","dist.to.TSG")
    df[,3:7] <- sapply(df[,3:7], as.numeric)
    # Weighted scores
    #Giving weight based on the distance
    df$weight <- ifelse(df$dist.to.TSG > 10000000, 0,
                        ifelse(df$dist.to.TSG <= 10000000 & df$dist.to.TSG > 8000000, 2,
                               ifelse(df$dist.to.TSG <= 8000000 & df$dist.to.TSG > 6000000, 4,
                                      ifelse(df$dist.to.TSG <= 6000000 & df$dist.to.TSG > 4000000, 6,
                                             ifelse(df$dist.to.TSG <= 4000000 & df$dist.to.TSG > 2000000, 8,10)))))
    df.weight <- df %>% group_by(bin) %>% summarise(Dist.TSG.weighted = sum(weight))
    
    Dist.tumoursupp.weighted[[level]][[cohort]] <- df.weight
  }
  print(cohort)
  
}

save(Dist.oncogenes.weighted,Dist.tumoursupp.weighted,
     file = distance_weighted_path)
