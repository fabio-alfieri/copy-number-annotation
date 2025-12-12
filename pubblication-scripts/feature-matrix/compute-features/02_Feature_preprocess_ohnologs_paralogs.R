
#Pre-processing of feature data
#Ohnologs and paralogs
#Tissue-general

# Set working directory 
setwd("/mnt/fabiogokce/gokce")

#Required libraries
library(readxl)
library(dplyr)

ohnologs <- list()

#Ohnologs from different datasets
#1. From Makino and Mclysaght paper
makino_mclysaght <- read_excel("./Downloads/Ohnologs/st01.xls", sheet = 5, skip = 1, col_types = "text")
makino_mclysaght <- makino_mclysaght[-c(9058:9060),]
#The ohnologs are provided with their ensembl gene IDs, and ensembl release 52 was used on the paper.
#However, that version is not available anymore. 
#Therefore, Ensembl release 54 was used to convert IDs to symbols. 
ensembl_release56 <- read.delim("./Downloads/ensembl_release56.txt")
makino_mclysaght <- merge(makino_mclysaght,ensembl_release56,by.x = "Ohnolog1",by.y="Ensembl.Gene.ID")
colnames(makino_mclysaght)[13] <- "Associated.Gene.Name1"
makino_mclysaght <- merge(makino_mclysaght,ensembl_release56,by.x = "Ohnolog2",by.y="Ensembl.Gene.ID")
colnames(makino_mclysaght)[14] <- "Associated.Gene.Name2"
ohno1 <- makino_mclysaght %>% group_by(Associated.Gene.Name1) %>% 
  summarise("ohnologs_mmpaper" = paste(Associated.Gene.Name2,collapse = "::"),
            "n_ohnolog_mmpaper" = n())
ohno1 <- as.data.frame(ohno1)
colnames(ohno1)[1] <- "Symbol" 
ohno2 <- makino_mclysagohno1ohno2 <- makino_mclysaght %>% group_by(Associated.Gene.Name2) %>% 
  summarise("ohnologs_mmpaper" = paste(Associated.Gene.Name1,collapse = "::"),
            "n_ohnolog_mmpaper" = n())
ohno2 <- as.data.frame(ohno2)
colnames(ohno2)[1] <- "Symbol" 
ohno <- rbind(ohno1,ohno2)
ohno <- ohno %>% group_by(Symbol) %>% 
  summarise("ohnologs_mmpaper" = paste(ohnologs_mmpaper,collapse = "::"),
            "n_ohnolog_mmpaper" = sum(n_ohnolog_mmpaper)) #7290
#Add the data to the final list
ohnologs[["mmpaper"]] <- ohno
#2. From Ohnolog database (different cutoffs)
files <- list.files(path = "./Downloads/Ohnologs", pattern = "^hsapiens.Pairs*", full.names = T, recursive = F)
for(file in files){
  cutoff <- strsplit(basename(file),"\\.")[[1]][3]
  ohnodatabase <- read.delim(file)
  ohno1 <- ohnodatabase %>% group_by(Symbol1) %>%
    summarise("ohnologs" = paste(Symbol2,collapse = "::"),
              "n_ohnolog" = n())
  ohno1 <- as.data.frame(ohno1)
  colnames(ohno1)[1] <- "Symbol"
  ohno2 <- ohnodatabase %>% group_by(Symbol2) %>%
    summarise("ohnologs" = paste(Symbol1,collapse = "::"),
              "n_ohnolog" = n())
  ohno2 <- as.data.frame(ohno2)
  colnames(ohno2)[1] <- "Symbol"
  ohno <- rbind(ohno1,ohno2)
  ohno <- ohno %>% group_by(Symbol) %>% 
    summarise("ohnologs" = paste(ohnologs,collapse = "::"),
              "n_ohnolog" = sum(n_ohnolog))
  ohnologs[[paste("OhnoDatabase",cutoff,sep = "_")]] <- ohno}
save(ohnologs, file = "./Data/OhnologList.RData")

#Paralogs
#1. Colm' data
#paralogs <- read.csv("./Downloads/all_paralog_genes_min20.csv") #Paralogs
paralog.features <- read.csv("./Downloads/Table_S8.csv") #Detailed features for paralogs
#Colm's cut-off for sequence identity (reciprocal, at least 20%)
paralog.features <- paralog.features[paralog.features$min_sequence_identity >= 0.20,3:5]
paralog.genes <- union(paralog.features$A1, paralog.features$A2)
paralogs <- c()
for(gene in paralog.genes){
  df <- paralog.features[paralog.features$A1 == gene | paralog.features$A2 == gene,]
  par <- setdiff(union(df$A1,df$A2),gene)
  paralogs <- rbind(paralogs,c(gene,paste(par,collapse = "::"),length(par)))
}
paralogs <- as.data.frame(paralogs)
colnames(paralogs) <- c("symbol","paralogs","n_paralogs")
save(paralogs, file = "./Data/ParalogList.RData")

