
#Preprocessing of feature data
#Chromatin states from RoadmapEpigenome

# Set working directory 

setwd("/mnt/fabiogokce/gokce") #If working on workstation

#Required libraries
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(maditr)
library(reshape2)
library(factoextra)

#Annotation table for chromatin states was prepared based on the table on RoadmapEpigenome
#Path to the annotation file "./Data/ChromatinStates_Info.xlsx"

#Downloading chromatin state files from RoadmapEpigenome
#Metadata for Reprocessed data from 127 Consolidated Epigenomes (111 Roadmap + 16 ENCODE) and 
#Unconsolidated Epigenomes were download from RoadmapEpigenome (https://egg2.wustl.edu/roadmap/web_portal/meta.html)
#Then chromatin state files were downloaded for the cells retrieved from the metadata

# metadata <- read.xlsx("./Downloads/Roadmap.metadata.qc.jul2013.xlsx", startRow = 4, sheet = 1, colNames = F)
# cells <- unique(metadata$X2)
# files <- paste(cells,"25_imputed12marks_segments.bed.gz",sep = "_")
# for(file in files){
#   url <- paste0("https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/",file)
#   dir <- paste0("./Downloads/ChmmModels/",gsub(".gz","",basename(url)))
#   download.file(url,dir)}

#Optional part I
#Getting general info about chromatin states - Their distribution through the genome
# bed.files <- list.files("./Downloads/ChmmModels", recursive = F, full.names = T)
# final.df <- c()
# for(file in bed.files){
#   chmm <- read.delim(file, header = F)
#   # Calculate the length of each state
#   chmm$length <- chmm$V3-chmm$V2
#   # Total length of scanned region for the normalization because this is different for different cells
#   total <- sum(chmm$length)
#   cell <- strsplit(basename(file),"_")[[1]][1]
#   # Coverage of each chromatin state in the corresponding cell
#   cov <- chmm %>% group_by(V4) %>% summarise("Total.length"=sum(length))
#   # Coverage per percentage (normalised)
#   cov$Per <- (cov$Total.length / total)*100
#   cov$cell <- cell
#   final.df <- rbind(final.df,cov)}
# 
# avg.cov <- final.df %>% group_by(V4) %>% summarise("avg.coverage"=mean(Per))
# avg.cov <- avg.cov[order(avg.cov$avg.coverage,decreasing = T),]
# avg.cov$V4 <- factor(avg.cov$V4, levels = avg.cov$V4)
# 
# barplot <- ggplot(avg.cov, aes(x=V4,y=avg.coverage)) +
#   geom_bar(stat = "identity", fill = "#0C475B", color = "#0C475B") +
#   theme_classic2() + xlab("Chromatin states") +
#   ylab("Average coverage\nCovarage: Total number of bases divided by total number of bases") +
#   labs(title = "Average covarage of chromatin states across cells") +
#   geom_text(aes(label = round(avg.coverage, digits = 2), x = V4, y = avg.coverage),
#             position = position_dodge(width = 0.8), vjust = -0.6)
# 
# pdf("./Plots/Features/Chromatin_states/AvgCovarege_chromatinstates_acrosscells.pdf",
#     width = 10, height = 6)
# barplot
# dev.off()

# Optional Part II - Date: July 31, 2024
# Visualization in IGV
bed.files <- list.files("./Downloads/ChmmModels", recursive = F, full.names = T)
# files we are interested in
#Retrieving TCGA names for each cell
#roadmap_metadata_with_tcga_projects.tsv (Table from Richard)
roadmap.tcga.map <- read.delim("./Downloads/roadmap_metadata_with_tcga_projects.tsv",header = F) #there is no E001 but it is stem cell no TCGA matching. 
mapped.tcga.cohorts <- setdiff(unique(roadmap.tcga.map$V12),c("stem_cell","muscle","FAT","vascular","small_intestine","placenta","spleen",""))
roadmap.tcga.map <- roadmap.tcga.map[roadmap.tcga.map$V12 %in% mapped.tcga.cohorts,c(1,12)]
colnames(roadmap.tcga.map) <- c("Cell.ID","TCGA.name")
#TCGA cohort naming conversion so the cohort names will match with the bin data
roadmap.tcga.map$TCGA.nameII <- ifelse(roadmap.tcga.map$TCGA.name == "COAD; READ","COADREAD",
                                       ifelse(roadmap.tcga.map$TCGA.name == "GBM; LGG","GBMLGG",roadmap.tcga.map$TCGA.name))
roadmap.tcga.map <- separate_rows(roadmap.tcga.map,3, sep = "; ")
roadmap.tcga.map <- as.data.frame(roadmap.tcga.map)
roadmap.tcga.map <- roadmap.tcga.map[roadmap.tcga.map$TCGA.nameII %in% 
                                       c("BRCA","COADREAD","PAAD","GBMLGG","LUSC","LUAD"),]
cells <- unique(roadmap.tcga.map$Cell.ID) # same cells for LUAD and LUSC

for(file in bed.files){
  cell <- strsplit(basename(file),"_")[[1]][1]
  if(!cell %in% cells)next
  chmm <- read.delim(file, header = F)
  for(state in unique(chmm$V4)){
    chmm.state <- chmm[chmm$V4 == state,]
    write.table(chmm.state[,1:3], 
                file = paste0("./IGV_Chmm_inputs/",cell,"_",state,".bed"), 
                col.names = F, row.names = F, sep = "\t", quote = F)
  }
  
}
  

#Method I: For each cell: Calculate the total number of bases for each state in bins at different levels
#Retrieving TCGA names for each cell
#roadmap_metadata_with_tcga_projects.tsv (Table from Richard)
roadmap.tcga.map <- read.delim("./Downloads/roadmap_metadata_with_tcga_projects.tsv",header = F) #there is no E001 but it is stem cell no TCGA matching. 
mapped.tcga.cohorts <- setdiff(unique(roadmap.tcga.map$V12),c("stem_cell","muscle","FAT","vascular","small_intestine","placenta","spleen",""))
roadmap.tcga.map <- roadmap.tcga.map[roadmap.tcga.map$V12 %in% mapped.tcga.cohorts,c(1,12)]
colnames(roadmap.tcga.map) <- c("Cell.ID","TCGA.name")
#TCGA cohort naming conversion so the cohort names will match with the bin data
roadmap.tcga.map$TCGA.nameII <- ifelse(roadmap.tcga.map$TCGA.name == "COAD; READ","COADREAD",
                                       ifelse(roadmap.tcga.map$TCGA.name == "GBM; LGG","GBMLGG",roadmap.tcga.map$TCGA.name))
roadmap.tcga.map <- separate_rows(roadmap.tcga.map,3, sep = "; ")
roadmap.tcga.map <- as.data.frame(roadmap.tcga.map)

#Backbone dataset: Start and end positions of bins at all levels 
load("./Data/All_levels_backbonetables.RData")
cohorts <- c("BRCA","LUAD","LUSC","CESC","THCA","HNSC","PAAD","COADREAD","GBMLGG",
             "SKCM","BLCA","PCPG","PRAD","KIRC","MESO","TGCT","KIRP","SARC",    
             "LIHC","ESCA","STAD","UCS","OV") # 23 cohorts
cells <- unique(roadmap.tcga.map[roadmap.tcga.map$TCGA.nameII %in% cohorts,]$Cell.ID)
#levels <- names(chr_backbone_namesfixed) # To run for all the levels
levels <- c("2Mbp", "4Mbp","6Mbp","8Mbp","10Mbp", "12Mbp","14Mbp","16Mbp","18Mbp","20Mbp","22Mbp","24Mbp",
            "26Mbp","28Mbp","30Mbp","32Mbp","34Mbp","36Mbp","38Mbp","40Mbp","42Mbp","44Mbp", "46Mbp", "48Mbp",
            "1Mbp","Chromosome","Arm","0.5Mbp","0.1Mbp","0.25Mbp","3Mbp") #Levels we are interested

#Function to calculate length for each chromatin state in bp
chmm.state.length <- function(coordinates.l){
  # Number of bases for states in the ith bin
  chr <- coordinates.l$chr
  bin <- coordinates.l$bin
  s <- coordinates.l$s
  e <- coordinates.l$e
  chmm.sub <- chmm[chmm[chmm$V1 == chr & (chmm$V3 >= s & chmm$V2 <= s),"index"]:
                     chmm[chmm$V1 == chr & (chmm$V3 >= e & chmm$V2 <= e),"index"],]
  chmm.sub$length <- ifelse(chmm.sub$V3 >= s & chmm.sub$V2 <= s, chmm.sub$V3-s,
                            ifelse(chmm.sub$V3 >= e & chmm.sub$V2 <= e, e-chmm.sub$V2, chmm.sub$V3 - chmm.sub$V2))
  #chmm.sub$bin <- paste(chr,s,e,sep = "_")
  chmm.sub$bin <- paste(gsub("chr","",chr),bin,sep = "_") # May 16, 2024
  return(chmm.sub)
}

ChrStates.bins <- list()

for(level in levels){
  # Bin, start and end locations
  coordinates <- do.call(rbind,chr_backbone_namesfixed[[level]])
  
  #Make a list for each bin to call function
  coord.list <- list()
  for(i in 1:nrow(coordinates)){
    coord.list[[i]] <- list("chr" = paste0("chr",coordinates$chr[i]),
                            "bin" = coordinates$bin[i],
                            "s" = as.numeric(coordinates$start_bin[i]),
                            "e" = as.numeric(coordinates$end_bin[i]))} 
  
  level.list <- list()
  for(cell in cells){
    cell.file <- paste0("./Downloads/ChmmModels/",cell,"_25_imputed12marks_segments.bed")
    chmm <- read.delim(cell.file, header = F)
    chmm$index <- as.numeric(rownames(chmm)) # Function will use this from the global environment
    #Call the function 
    res <- mclapply(coord.list, chmm.state.length, mc.cores = 20)
    chrstate.m <- do.call(rbind,res)
    # Sum the number of bases for each state and for each bin
    chrstate.m <- chrstate.m %>% group_by(bin, V4) %>% summarise(total_bases=sum(length))
    level.list[[cell]] <- chrstate.m}
  print(level)
  ChrStates.bins[[level]] <- level.list
  #save(ChrStates.bins,file = "./Data/ChromatinStates_all_levels_cells.RData")
  save(ChrStates.bins,file = "./Data/ChromatinStates_all_levels_cells_binnamesfixed.RData")
}

#Fixing bin names for Arm and chromosome levels
load("./Data/ChromatinStates_all_levels_cells_binnamesfixed.RData")

for(level in c("Arm","Chromosome")){
  for(cell in names(ChrStates.bins[[level]])){
    df <- ChrStates.bins[[level]][[cell]]
    if(level == "Arm"){df$bin <- unlist(lapply(df$bin, function(x) strsplit(x,"_")[[1]][2]))}
    else{df$bin <- unlist(lapply(df$bin, function(x) paste(strsplit(x,"_")[[1]][2],strsplit(x,"_")[[1]][3],sep = "_")))}
    ChrStates.bins[[level]][[cell]] <- df}}

save(ChrStates.bins, file = "./Data/ChromatinStates_all_levels_cells_binnamesfixed.RData")

rm(list = ls())

#Method II: For each cell: Count the abundance for each state in bins at different levels
#Retrieving TCGA names for each cell
#roadmap_metadata_with_tcga_projects.tsv (Table from Richard)
roadmap.tcga.map <- read.delim("./Downloads/roadmap_metadata_with_tcga_projects.tsv",header = F) #there is no E001 but it is stem cell no TCGA matching. 
mapped.tcga.cohorts <- setdiff(unique(roadmap.tcga.map$V12),c("stem_cell","muscle","FAT","vascular","small_intestine","placenta","spleen",""))
roadmap.tcga.map <- roadmap.tcga.map[roadmap.tcga.map$V12 %in% mapped.tcga.cohorts,c(1,12)]
colnames(roadmap.tcga.map) <- c("Cell.ID","TCGA.name")
#TCGA cohort naming conversion so the cohort names will match with the bin data
roadmap.tcga.map$TCGA.nameII <- ifelse(roadmap.tcga.map$TCGA.name == "COAD; READ","COADREAD",
                                       ifelse(roadmap.tcga.map$TCGA.name == "GBM; LGG","GBMLGG",roadmap.tcga.map$TCGA.name))
roadmap.tcga.map <- separate_rows(roadmap.tcga.map,3, sep = "; ")
roadmap.tcga.map <- as.data.frame(roadmap.tcga.map)

#Backbone dataset: Start and end positions of bins at all levels 
load("./Data/All_levels_backbonetables.RData")
cohorts <- c("BRCA","LUAD","LUSC","CESC","THCA","HNSC","PAAD","COADREAD","GBMLGG",
             "SKCM","BLCA","PCPG","PRAD","KIRC","MESO","TGCT","KIRP","SARC",    
             "LIHC","ESCA","STAD","UCS","OV") # 23 cohorts
cells <- unique(roadmap.tcga.map[roadmap.tcga.map$TCGA.nameII %in% cohorts,]$Cell.ID)
#levels <- names(chr_backbone_namesfixed) # To run for all the levels
levels <- c("2Mbp", "4Mbp","6Mbp","8Mbp","10Mbp", "12Mbp","14Mbp","16Mbp","18Mbp","20Mbp","22Mbp","24Mbp",
            "26Mbp","28Mbp","30Mbp","32Mbp","34Mbp","36Mbp","38Mbp","40Mbp","42Mbp","44Mbp", "46Mbp", "48Mbp",
            "1Mbp","Chromosome","Arm","0.5Mbp","0.1Mbp","0.25Mbp","3Mbp") #Levels we are interested

#Function to calculate length for each chromatin state in bp
chmm.state.length <- function(coordinates.l){
  # Number of bases for states in the ith bin
  chr <- coordinates.l$chr
  bin <- coordinates.l$bin
  s <- coordinates.l$s
  e <- coordinates.l$e
  chmm.sub <- chmm[chmm[chmm$V1 == chr & (chmm$V3 >= s & chmm$V2 <= s),"index"]:
                     chmm[chmm$V1 == chr & (chmm$V3 >= e & chmm$V2 <= e),"index"],]
  chmm.sub <- as.data.frame(table(chmm.sub$V4))
  chmm.sub$bin <- paste(gsub("chr","",chr),bin,sep = "_") # May 16, 2024
  return(chmm.sub)
}

ChrStates.bins <- list()

for(level in levels){
  # Bin, start and end locations
  coordinates <- do.call(rbind,chr_backbone_namesfixed[[level]])
  
  #Make a list for each bin to call function
  coord.list <- list()
  for(i in 1:nrow(coordinates)){
    coord.list[[i]] <- list("chr" = paste0("chr",coordinates$chr[i]),
                            "bin" = coordinates$bin[i],
                            "s" = as.numeric(coordinates$start_bin[i]),
                            "e" = as.numeric(coordinates$end_bin[i]))} 
  
  level.list <- list()
  for(cell in cells){
    cell.file <- paste0("./Downloads/ChmmModels/",cell,"_25_imputed12marks_segments.bed")
    chmm <- read.delim(cell.file, header = F)
    chmm$index <- as.numeric(rownames(chmm)) # Function will use this from the global environment
    #Call the function 
    res <- mclapply(coord.list, chmm.state.length, mc.cores = 20)
    chrstate.m <- do.call(rbind,res)
    colnames(chrstate.m) <- c("Chmm.state","Abundance","bin")
    level.list[[cell]] <- chrstate.m}
  print(level)
  ChrStates.bins[[level]] <- level.list
  save(ChrStates.bins,file = "./Data/ChromatinStates_all_levels_cells_abundance_methodII.RData")
}

#Fixing bin names for Arm and chromosome levels
load("./Data/ChromatinStates_all_levels_cells_abundance_methodII.RData")

for(level in c("Arm","Chromosome")){
  for(cell in names(ChrStates.bins[[level]])){
    df <- ChrStates.bins[[level]][[cell]]
    if(level == "Arm"){df$bin <- unlist(lapply(df$bin, function(x) strsplit(x,"_")[[1]][2]))}
    else{df$bin <- unlist(lapply(df$bin, function(x) paste(strsplit(x,"_")[[1]][2],strsplit(x,"_")[[1]][3],sep = "_")))}
    ChrStates.bins[[level]][[cell]] <- df}}

save(ChrStates.bins, file = "./Data/ChromatinStates_all_levels_cells_abundance_methodII.RData")

#Optional part II - Visualization - PCA and correlations between chromatin states - Method I
#PCA for cells and correlation between different states at different levels
#PCA - 1Mbp 
load("./Data/ChromatinStates_all_levels_cells_binnamesfixed.RData")
level <- ChrStates.bins[["1Mbp"]]
m <- c()
for(cell in names(level)){
  df <- level[[cell]]
  df$column <- paste(df$bin,df$V4,sep = "_")
  df$cell <- cell
  m <- rbind(m,df)}
df <- m %>% reshape2::dcast(cell~column, value.var = "total_bases")
df <- as.data.frame(df)
#Table for Cell-Tissue name matching (Supplementary column)
roadmap.tcga.map <- read.delim("./Downloads/roadmap_metadata_with_tcga_projects.tsv",header = F) #there is no E001 but it is stem cell no TCGA matching. 
mapped.tcga.cohorts <- setdiff(unique(roadmap.tcga.map$V12),c("stem_cell","muscle","FAT","vascular","small_intestine","placenta","spleen",""))
roadmap.tcga.map <- roadmap.tcga.map[roadmap.tcga.map$V12 %in% mapped.tcga.cohorts,c(1,12)]
colnames(roadmap.tcga.map) <- c("Cell.ID","TCGA.name")

df <- merge(roadmap.tcga.map,df,by.x = "Cell.ID",by.y = "cell")
rownames(df) <- df$Cell.ID
df <- df[,-1]
df[is.na(df)] <- 0
df.pca <- df[,-1]
var.zero.columns <- which(apply(df.pca, 2, var)==0)
df.pca <- df.pca[,-var.zero.columns]
df.pca <- df.pca[,-1]
res.pca <- prcomp(df.pca, scale = TRUE)
groups <- as.factor(df$TCGA.name)
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#073b4c","#118ab2","#06d6a0","#ef476f",
                         "#ae2012","#a5be00","#ee6c4d","#f08080",
                         "#000814","#7b2cbf","#eeef20","#0466c8",
                         "#d68c45","#5e3023"),
             legend.title = "TCGA cohorts",
             repel = TRUE,
             title = "1Mbp - All bins")

#Comparing the number of bases across chromosomes
rm(list = ls())
load("./Data/All_levels_backbonetables.RData") # For chromosome length
load("./Data/ChromatinStates_all_levels_cells_binnamesfixed.RData")
#Chromosome level
coordinates <- do.call(rbind,chr_backbone_namesfixed[["Chromosome"]])
coordinates$length <- coordinates$end_bin - coordinates$start_bin
cell.l <- ChrStates.bins[["Chromosome"]] #Number of bases for chromosomes
final.m <- c()
for(cell in names(cell.l)){
  df <- cell.l[[cell]]
  df$chr <- unlist(lapply(df$bin, function(x) strsplit(x,"_")[[1]][1]))
  df <- merge(df, coordinates, by = "chr")
  for(chmm in unique(df$V4)){
    m <- df[df$V4 == chmm,]
    res <- cor.test(m$total_bases, m$length)
    final.m <- rbind(final.m, c(cell, chmm, as.numeric(res$estimate),res$p.value))}}
final.m <- as.data.frame(final.m)
colnames(final.m) <- c("Cell","ChrState","Rho","pvalue")
final.m[,3:4] <- sapply(final.m[,3:4], as.numeric)

chr.plot <- ggplot(final.m, aes(x = ChrState, y = Rho)) + geom_boxplot() +
  geom_point(size = .5) + theme_bw() + 
  labs(title = "Chromosome level") + xlab("Chromatin states") +
  ylab("Correlation between number of bases of Chr.State\nand chromosome length within a cell")
pdf("./Plots/Features/Chromatin_states/Correlation_btw_ChromatinStates_chromosomeLength.pdf", 
    width = 9, height = 5)
chr.plot
dev.off()

#Optional part II - Visualization - PCA and correlations between chromatin states - Method II
#PCA for cells and correlation between different states at different levels
#PCA - 1Mbp 
load("./Data/ChromatinStates_all_levels_cells_abundance_methodII.RData")
level <- ChrStates.bins[["1Mbp"]]
m <- c()
for(cell in names(level)){
  df <- level[[cell]]
  df$column <- paste(df$bin,df$Chmm.state,sep = "_")
  df$cell <- cell
  m <- rbind(m,df)}
df <- m %>% reshape2::dcast(cell~column, value.var = "Abundance")
df <- as.data.frame(df)
#Table for Cell-Tissue name matching (Supplementary column)
roadmap.tcga.map <- read.delim("./Downloads/roadmap_metadata_with_tcga_projects.tsv",header = F) #there is no E001 but it is stem cell no TCGA matching. 
mapped.tcga.cohorts <- setdiff(unique(roadmap.tcga.map$V12),c("stem_cell","muscle","FAT","vascular","small_intestine","placenta","spleen",""))
roadmap.tcga.map <- roadmap.tcga.map[roadmap.tcga.map$V12 %in% mapped.tcga.cohorts,c(1,12)]
colnames(roadmap.tcga.map) <- c("Cell.ID","TCGA.name")

df <- merge(roadmap.tcga.map,df,by.x = "Cell.ID",by.y = "cell")
rownames(df) <- df$Cell.ID
df <- df[,-1]
df[is.na(df)] <- 0
df.pca <- df[,-1]
var.zero.columns <- which(apply(df.pca, 2, var)==0)
df.pca <- df.pca[,-var.zero.columns]
res.pca <- prcomp(df.pca, scale = TRUE)
groups <- as.factor(df$TCGA.name)
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#073b4c","#118ab2","#06d6a0","#ef476f",
                         "#ae2012","#a5be00","#ee6c4d","#f08080",
                         "#000814","#7b2cbf","#eeef20","#0466c8",
                         "#d68c45","#5e3023"),
             legend.title = "TCGA cohorts",
             repel = TRUE,
             title = "1Mbp - All bins")

#Comparing the number of bases across chromosomes
rm(list = ls())
load("./Data/All_levels_backbonetables.RData") # For chromosome length
load("./Data/ChromatinStates_all_levels_cells_abundance_methodII.RData")
#Chromosome level
coordinates <- do.call(rbind,chr_backbone_namesfixed[["Chromosome"]])
coordinates$length <- coordinates$end_bin - coordinates$start_bin
cell.l <- ChrStates.bins[["Chromosome"]] #Number of bases for chromosomes
final.m <- c()
for(cell in names(cell.l)){
  df <- cell.l[[cell]]
  df$chr <- unlist(lapply(df$bin, function(x) strsplit(x,"_")[[1]][1]))
  df <- merge(df, coordinates, by = "chr")
  for(chmm in unique(df$Chmm.state)){
    m <- df[df$Chmm.state == chmm,]
    res <- cor.test(m$Abundance, m$length)
    final.m <- rbind(final.m, c(cell, chmm, as.numeric(res$estimate),res$p.value))}}
final.m <- as.data.frame(final.m)
colnames(final.m) <- c("Cell","ChrState","Rho","pvalue")
final.m[,3:4] <- sapply(final.m[,3:4], as.numeric)

chr.plot <- ggplot(final.m, aes(x = ChrState, y = Rho)) + geom_boxplot() +
  geom_point(size = .5) + theme_bw() + 
  labs(title = "Chromosome level") + xlab("Chromatin states") +
  ylab("Correlation between abundance of Chr.State\nand chromosome length within a cell")
pdf("./Plots/Features/Chromatin_states/Correlation_btw_ChromatinStates_abundance_chromosomeLength.pdf", 
    width = 9, height = 5)
chr.plot
dev.off()

#----------------------------Tissue-level

rm(list = ls())

#METHOD I
#Chromatin states per tissue: Averaging the number of bases across cells, 
#and then scaling within the tissue table for Cell-Tissue name matching
roadmap.tcga.map <- read.delim("./Downloads/roadmap_metadata_with_tcga_projects.tsv",header = F) #there is no E001 but it is stem cell no TCGA matching. 
mapped.tcga.cohorts <- setdiff(unique(roadmap.tcga.map$V12),c("stem_cell","muscle","FAT","vascular","small_intestine","placenta","spleen",""))
roadmap.tcga.map <- roadmap.tcga.map[roadmap.tcga.map$V12 %in% mapped.tcga.cohorts,c(1,12)]
colnames(roadmap.tcga.map) <- c("Cell.ID","TCGA.name")
#TCGA cohort naming conversion so the cohort names will match with the bin data
roadmap.tcga.map$TCGA.nameII <- ifelse(roadmap.tcga.map$TCGA.name == "COAD; READ","COADREAD",
                                       ifelse(roadmap.tcga.map$TCGA.name == "GBM; LGG","GBMLGG",roadmap.tcga.map$TCGA.name))
roadmap.tcga.map <- separate_rows(roadmap.tcga.map,3, sep = "; ")
roadmap.tcga.map <- as.data.frame(roadmap.tcga.map)
cohorts <- c("BRCA","LUAD","LUSC","CESC","THCA","HNSC","PAAD","COADREAD","GBMLGG",
             "SKCM","BLCA","PCPG","PRAD","KIRC","MESO","TGCT","KIRP","SARC",    
             "LIHC","ESCA","STAD","UCS","OV")
tissues <- unique(intersect(roadmap.tcga.map$TCGA.nameII,cohorts))

load("./Data/ChromatinStates_all_levels_cells_binnamesfixed.RData") #Chromatin states per cells
levels <- names(ChrStates.bins)
ChrStates.tissues <- list()
for(level in levels){
  level.set <- list()
  level.data <- ChrStates.bins[[level]]
  for(tissue in tissues){
    tissue.set <- list()
    cells <- roadmap.tcga.map[roadmap.tcga.map$TCGA.nameII == tissue,"Cell.ID"] #Find cells mapping to the tissue
    tissue.df <- c()
    tissue.df.norm <- c()
    for(cell in cells){
      df.cell <- level.data[[cell]]
      df.cell$Cell <- cell
      tissue.df <- rbind(tissue.df,df.cell)
      #Normalization by the total identified bases per chromosome
      df.cell$chr <- unlist(lapply(df.cell$bin, function(x) strsplit(x,"_")[[1]][1]))
      bases.per.chr <- df.cell %>% group_by(chr) %>% summarise("Sum" = sum(total_bases))
      df.cell <- merge(df.cell, bases.per.chr, by = "chr")
      df.cell$Per.per.chr <- (df.cell$total_bases / df.cell$Sum)*100
      tissue.df.norm <- rbind(tissue.df.norm,df.cell)}
    #Take average of the cells
    df <- tissue.df %>% group_by(bin,V4) %>% summarise("Mean"=mean(total_bases))
    df <- df %>% dcast(bin~V4, value.var = "Mean")
    df <- as.data.frame(df)
    rownames(df) <- df$bin
    df <- df[,-1]
    df[is.na(df)] <- 0
    m.scaled <- as.data.frame(apply(df, 2,  function(x) scale(x))) #within column scale
    rownames(m.scaled) <- rownames(df)
    tissue.set[["Scaled"]] <- m.scaled
    tissue.set[["Counts"]] <- df
    #Normalized counts
    df.norm <- tissue.df.norm %>% group_by(bin,V4) %>% summarise("Mean"=mean(Per.per.chr))
    df.norm <- df.norm %>% reshape2::dcast(bin~V4, value.var = "Mean")
    df.norm <- as.data.frame(df.norm)
    rownames(df.norm) <- df.norm$bin
    df.norm <- df.norm[,-1]
    df.norm[is.na(df.norm)] <- 0
    tissue.set[["Normalised.by.chr.coverage"]] <- df.norm
    level.set[[tissue]] <- tissue.set}
  ChrStates.tissues[[level]] <- level.set}
save(ChrStates.tissues,file = "./Data/ChromatinStates_all_levels_tissues_avgcounts&scaled_avgcount.RData")

#METHOD II

rm(list = ls())

#Chromatin states per tissue: Averaging the abundance of chromatin states across cells, 
#and then scaling within the tissue table for Cell-Tissue name matching
roadmap.tcga.map <- read.delim("./Downloads/roadmap_metadata_with_tcga_projects.tsv",header = F) #there is no E001 but it is stem cell no TCGA matching. 
mapped.tcga.cohorts <- setdiff(unique(roadmap.tcga.map$V12),c("stem_cell","muscle","FAT","vascular","small_intestine","placenta","spleen",""))
roadmap.tcga.map <- roadmap.tcga.map[roadmap.tcga.map$V12 %in% mapped.tcga.cohorts,c(1,12)]
colnames(roadmap.tcga.map) <- c("Cell.ID","TCGA.name")
#TCGA cohort naming conversion so the cohort names will match with the bin data
roadmap.tcga.map$TCGA.nameII <- ifelse(roadmap.tcga.map$TCGA.name == "COAD; READ","COADREAD",
                                       ifelse(roadmap.tcga.map$TCGA.name == "GBM; LGG","GBMLGG",roadmap.tcga.map$TCGA.name))
roadmap.tcga.map <- separate_rows(roadmap.tcga.map,3, sep = "; ")
roadmap.tcga.map <- as.data.frame(roadmap.tcga.map)
cohorts <- c("BRCA","LUAD","LUSC","CESC","THCA","HNSC","PAAD","COADREAD","GBMLGG",
             "SKCM","BLCA","PCPG","PRAD","KIRC","MESO","TGCT","KIRP","SARC",    
             "LIHC","ESCA","STAD","UCS","OV")
tissues <- unique(intersect(roadmap.tcga.map$TCGA.nameII,cohorts))

load("./Data/ChromatinStates_all_levels_cells_abundance_methodII.RData")#Chromatin states per cells
levels <- names(ChrStates.bins)
ChrStates.tissues <- list()
for(level in levels){
  level.set <- list()
  level.data <- ChrStates.bins[[level]]
  for(tissue in tissues){
    tissue.set <- list()
    cells <- roadmap.tcga.map[roadmap.tcga.map$TCGA.nameII == tissue,"Cell.ID"] #Find cells mapping to the tissue
    tissue.df <- c()
    tissue.df.norm <- c()
    for(cell in cells){
      df.cell <- level.data[[cell]]
      df.cell$Cell <- cell
      tissue.df <- rbind(tissue.df,df.cell)
      #Normalization by the total identified abundance per chromosome
      df.cell$chr <- unlist(lapply(df.cell$bin, function(x) strsplit(x,"_")[[1]][1]))
      bases.per.chr <- df.cell %>% group_by(chr) %>% summarise("Sum" = sum(Abundance))
      df.cell <- merge(df.cell, bases.per.chr, by = "chr")
      df.cell$Per.per.chr <- (df.cell$Abundance / df.cell$Sum)*100
      tissue.df.norm <- rbind(tissue.df.norm,df.cell)}
    #Take average of the cells
    df <- tissue.df %>% group_by(bin,Chmm.state) %>% summarise("Mean"=mean(Abundance))
    df <- df %>% dcast(bin~Chmm.state, value.var = "Mean")
    df <- as.data.frame(df)
    rownames(df) <- df$bin
    df <- df[,-1]
    df[is.na(df)] <- 0
    m.scaled <- as.data.frame(apply(df, 2,  function(x) scale(x))) #within column scale
    rownames(m.scaled) <- rownames(df)
    tissue.set[["Scaled"]] <- m.scaled
    tissue.set[["Counts"]] <- df
    #Normalized counts
    df.norm <- tissue.df.norm %>% group_by(bin,Chmm.state) %>% summarise("Mean"=mean(Per.per.chr))
    df.norm <- df.norm %>% reshape2::dcast(bin~Chmm.state, value.var = "Mean")
    df.norm <- as.data.frame(df.norm)
    rownames(df.norm) <- df.norm$bin
    df.norm <- df.norm[,-1]
    df.norm[is.na(df.norm)] <- 0
    tissue.set[["Normalised.by.chr.coverage"]] <- df.norm
    level.set[[tissue]] <- tissue.set}
  ChrStates.tissues[[level]] <- level.set}
save(ChrStates.tissues,file = "./Data/ChromatinStates_all_levels_tissues_avgcounts&scaled_avgcount_methodII.RData")

#-------------------------------------

#Correlations between different states
#Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)}
#load("./Data/ChromatinStates_all_levels_tissues_avgcounts&scaled_avgcount.RData")
load("./Data/ChromatinStates_all_levels_tissues_avgcounts&scaled_avgcount_methodII.RData")
tissues <- names(ChrStates.tissues[[1]]) # 14 tissues
levels <- names(ChrStates.tissues) # 31 levels
correlation.plots<-list()
for(level in levels){
  level.data <- ChrStates.tissues[[level]]
  level.plots <- list()
  for(tissue in tissues){
    df <- level.data[[tissue]]$Counts
    #df <- level.data[[tissue]]$Normalised.by.chr.coverage
    df <- round(cor(df, method = "spearman"),2)
    df <- get_upper_tri(df)
    melted_df <- reshape2::melt(df,na.rm = TRUE)
    p<-ggplot(data = melted_df, aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile(color="white") + theme_minimal() +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0, limit = c(-1,1), space = "Lab", 
                           name="Spearman\nCorrelation") +
      theme(axis.text.x = element_text(angle = 90)) +
      labs(title = tissue) + xlab("") + ylab("")
    level.plots[[tissue]] <- p}
  all <- ggarrange(plotlist = level.plots, common.legend = TRUE, legend = "bottom")
  correlation.plots[[level]] <- all}

# pdf("./Plots/Features/Chromatin_states/Correlations_among_chrstates_length_0.1mbp.pdf",width = 20,height = 15)
# correlation.plots[["0.1Mbp"]]
# dev.off()
# 
# pdf("./Plots/Features/Chromatin_states/Correlations_among_chrstates_length_10mbp.pdf",width = 20,height = 15)
# correlation.plots[["10Mbp"]]
# dev.off()
# 
# pdf("./Plots/Features/Chromatin_states/Correlations_among_chrstates_length_Armp.pdf",width = 20,height = 15)
# correlation.plots[["Arm"]]
# dev.off()

pdf("./Plots/Features/Chromatin_states/Correlations_among_chrstates_abundance_0.1mbp.pdf",width = 20,height = 15)
correlation.plots[["0.1Mbp"]]
dev.off()

pdf("./Plots/Features/Chromatin_states/Correlations_among_chrstates_abundance_10mbp.pdf",width = 20,height = 15)
correlation.plots[["10Mbp"]]
dev.off()

pdf("./Plots/Features/Chromatin_states/Correlations_among_chrstates_abundance_Armp.pdf",width = 20,height = 15)
correlation.plots[["Arm"]]
dev.off()


# Date: May 31, 2024
#Prepare one final matrix with different calculations

rm(list = ls())
library(rlist)

assign("Length",get(load("./Data/ChromatinStates_all_levels_tissues_avgcounts&scaled_avgcount.RData")))
assign("Abundance",get(load("./Data/ChromatinStates_all_levels_tissues_avgcounts&scaled_avgcount_methodII.RData")))
rm("ChrStates.tissues")

data.list <- list("Abundance" = Abundance,
                  "Length" = Length)

Chromatin.States <- list()

levels <- names(Abundance)
tissues <- names(Abundance$`2Mbp`)

for(level in levels){
  for(tissue in tissues){
    tissue.l <- list()
    for(name in names(data.list)){
      res.l <- data.list[[name]][[level]][[tissue]]
      new.names <- paste(name,names(res.l),sep = "_")
      names(res.l) <- new.names
      tissue.l <- append(tissue.l,res.l)}
    #Abundance binary matrix
    abundance.df <- tissue.l$Abundance_Counts
    abundance.df[abundance.df > 0] <- 1
    tissue.l[["Abundance_binary"]] <- abundance.df
    Chromatin.States[[level]][[tissue]] <- tissue.l}}

# Counts, scaled, normalised for length and abundance methods (with an binary abundance matrix)
save(Chromatin.States, file = "./Data/Chromatin_states_all_levels_tissues_all_methods.RData")
