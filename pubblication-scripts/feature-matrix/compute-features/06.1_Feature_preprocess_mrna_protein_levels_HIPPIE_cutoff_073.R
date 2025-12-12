
# Feature optimization
# GTEx mRNA for PPIs from HIPPIE (with confidence cutoff)
# Date: October 5, 2024

rm(list = ls())

# Set working directory 
setwd("/mnt/fabiogokce/gokce")

# Required libraries
library(parallel)
library(stringr)


#GTEx data - In the final feature set, I am using scaled relative tissue medians, thus only calculate features for that datasets
load("./Data/GTEx_different_versions_TCGA_cohort_map.RData")
common.tissues <- setdiff(colnames(datasets$scaled.GTEx.v8.TPM),c("Gene.stable.ID","Gene.name"))
datasets <- datasets[5]

# Analysis
# 1. For all genes 

# Functions
# Function to find gene groups for a bin - Only for changed gene lists
gene.lists.per.bin <- function(bin){
  gene.set <- list("ppis.trans" = as.character(strsplit(backbone[backbone$bin == bin, "PPIs_073.trans-HIPPIE"],"::")[[1]]), #PPIs in trans
                   "all.int.trans" = union(as.character(strsplit(backbone[backbone$bin == bin, "Partners.trans"],"::")[[1]]),
                                           as.character(strsplit(backbone[backbone$bin == bin, "PPIs_073.trans-HIPPIE"],"::")[[1]])) #All interactions
  ) #For each bin, define a new gene sets
  names(gene.set) <- paste(bin, names(gene.set), sep = "_")
  return(gene.set)}

# Function to calculate mean/median expression values
expression.scores <- function(genes){
  df <- exp.data[exp.data$Gene.name %in% genes,]
  if(nrow(df) == 0){ #In the case of no partner, PPIs or no covarage by expression data
    res <- data.frame("Cohorts" = colnames(df)[-c(1,2)],
                      "Mean" = rep(NA, length(colnames(df)[-c(1,2)])),
                      "Median" = rep(NA, length(colnames(df)[-c(1,2)])))}
  else{
    res <- data.frame("Cohorts" = colnames(df)[-c(1,2)],
                      "Mean" = colMeans(df[,-c(1,2)], na.rm = TRUE),
                      "Median" = apply(df[,-c(1,2)],2,median, na.rm = TRUE))}
  return(res)}

# Backbone tables
load("./Data/Backbone_tables_with_non_tissue_specific_features_CORUM_HIPPIE.RData")
# Columns needed
selected.columns <- c("bin","chr","start_bin","end_bin",
                      "genes","n_genes","Partners.trans","PPIs_073.trans-HIPPIE")

levels <- names(backbone_tables_w_features)


Results <- list() # Do not run this when running additional part
for(level in levels){
  backbone <- backbone_tables_w_features[[level]]
  backbone <- backbone[, colnames(backbone) %in% selected.columns]
  # Remove bins without any gene
  backbone <- backbone[backbone$n_genes != 0,]
  # Find genes in each group for bins and make one final list
  gene.sets <- mclapply(backbone$bin, gene.lists.per.bin, mc.cores = 30)
  gene.sets.merged <- list()
  for(i in 1:length(gene.sets)){gene.sets.merged <- append(gene.sets.merged, gene.sets[[i]])}
  
  for(dataset in names(datasets)){
    exp.data <- datasets[[dataset]]
    exp.data  <- exp.data[,colnames(exp.data) %in% c("Gene.stable.ID","Gene.name",common.tissues)]
    res.df <- mclapply(gene.sets.merged, expression.scores, mc.cores = 30)
    final.df <- do.call(rbind,res.df)
    
    if(level == "Arm"){
      final.df$bin <- unlist(lapply(rownames(final.df), function(x) strsplit(x,"_")[[1]][1]))
      final.df$GeneSet <- unlist(lapply(rownames(final.df), function(x) strsplit(x,"_")[[1]][2]))
    }
    else{
      final.df$bin <- unlist(lapply(rownames(final.df), function(x) paste(strsplit(x,"_")[[1]][1],strsplit(x,"_")[[1]][2],sep = "_")))
      final.df$GeneSet <- unlist(lapply(rownames(final.df), function(x) strsplit(x,"_")[[1]][3]))
    }
    
    final.df$GeneSet <- str_remove(final.df$GeneSet, paste0(".",paste(common.tissues,collapse = "|.")))
    final.df$GeneSet <- gsub("\\.\\d+","",final.df$GeneSet)
    final.df$GTExData <- dataset
    Results[[level]][[dataset]] <- final.df
    print(dataset)
  }
  
  print(level)}

save(Results, file = "./Data/GTEx_rna_protein_level_differentdatasets_per_bins_HIPPIE_cutoff_073.RData")


