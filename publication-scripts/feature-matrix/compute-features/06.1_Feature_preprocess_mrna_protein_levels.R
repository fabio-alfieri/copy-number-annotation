
# Feature optimization
# GTEx mRNA and protein levels
# Date: April 15, 2024

rm(list = ls())

packages <- c(
  "stringr", "parallel"
)

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

gtex_path <- "./Data/GTEx_different_versions_TCGA_cohort_map.RData"
gtex_additional_path <- "./Data/GTEx_rna_protein_level_differentdatasets_per_bins.RData" 
backbone_table <- "./Data/Backbone_tables_with_non_tissue_specific_features.RData"
gtex_output_path <- "./Data/GTEx_rna_protein_level_differentdatasets_per_bins.RData"
gtex_v8_path <- "./Data/Tissue_specific_genes_based_on_GTEx_v8.RData"
tissue_specific_proteins_path <- "./Data/Tissue_specific_proteins_based_on_eGTEx_TSscore.RData"
gtex_tissue_specific <- "./Data/GTEx_rna_protein_level_differentdatasets_per_bins_Tissue_specific.RData"


#GTEx data
load(gtex_path)
common.tissues <- setdiff(colnames(datasets$eGTEX.rna),c("Gene.stable.ID","Gene.name"))
for(name in names(datasets)[-1]){
  common.tissues <- intersect(common.tissues,
                              setdiff(colnames(datasets[[name]]),c("Gene.stable.ID","Gene.name")))}

# Analysis
# 1. For all genes & proteins

# Functions
# Function to find gene groups for a bin
gene.lists.per.bin <- function(bin){
  gene.set <- list("genes" = as.character(strsplit(backbone[backbone$bin == bin, "genes"],"::")[[1]]), #Genes on the bin
                   "partners.trans" = as.character(strsplit(backbone[backbone$bin == bin, "Partners.trans"],"::")[[1]]), #Partners in trans
                   "ppis.trans" = as.character(strsplit(backbone[backbone$bin == bin, "PPIs.trans"],"::")[[1]]), #PPIs in trans
                   "all.int.trans" = union(as.character(strsplit(backbone[backbone$bin == bin, "Partners.trans"],"::")[[1]]),
                                           as.character(strsplit(backbone[backbone$bin == bin, "PPIs.trans"],"::")[[1]])) #All interactions
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
load(backbone_table)
# Columns needed
selected.columns <- c("bin","chr","start_bin","end_bin",
                      "genes","n_genes","Partners.trans","PPIs.trans")

levels <- names(backbone_tables_w_features_non_tissue_specific)

#-----Additional part start-----
# This part was added on September 16, 2024 to calculate the expression-related features for the other available datasets in GTEx v6 and v8
load(gtex_additional_path) # Results
datasets <- datasets[3:5]
common.tissues <- setdiff(colnames(datasets$GTEx.v6.RPKM),c("Gene.stable.ID","Gene.name"))
for(name in names(datasets)[-1]){
  common.tissues <- intersect(common.tissues,
                              setdiff(colnames(datasets[[name]]),c("Gene.stable.ID","Gene.name")))}
#-----Additional part end-----

Results <- list() # Do not run this when running additional part
for(level in levels){
  backbone <- backbone_tables_w_features_non_tissue_specific[[level]]
  backbone <- backbone[, colnames(backbone) %in% selected.columns]
  # Remove bins without any gene
  backbone <- backbone[backbone$n_genes != 0,]
  # Find genes in each group for bins and make one final list
  gene.sets <- mclapply(backbone$bin, gene.lists.per.bin, mc.cores = 20)
  gene.sets.merged <- list()
  for(i in 1:length(gene.sets)){gene.sets.merged <- append(gene.sets.merged, gene.sets[[i]])}
  
  for(dataset in names(datasets)){
    exp.data <- datasets[[dataset]]
    exp.data  <- exp.data[,colnames(exp.data) %in% c("Gene.stable.ID","Gene.name",common.tissues)]
    res.df <- mclapply(gene.sets.merged, expression.scores, mc.cores = 20)
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

save(Results, file = gtex_output_path)


# 2. For Tissue-specific genes & proteins

rm(list = setdiff(ls(),c("backbone_tables_w_features_non_tissue_specific",
                         "datasets",
                         "common.tissues",
                         "levels",
                         "selected.columns")))
                                        
load(gtex_v8_path) #Tissue-specific genes defined based on RNA level data
load(tissue_specific_proteins_path) #Tissue-specific proteins defined based on eGTEx paper TS score

#-----Additional part start-----
# This part was added on September 16, 2024 to calculate the expression-related features for the other available datasets in GTEx v6 and v8
load(gtex_tissue_specific)
#-----Additional part end-----

# Functions
# Function to find gene groups for a bin
gene.lists.per.bin <- function(bin){
  gene.set <- list("genes" = intersect(tissue.genes,as.character(strsplit(backbone[backbone$bin == bin, "genes"],"::")[[1]])), #Genes on the bin
                   "partners.trans" = intersect(tissue.genes,as.character(strsplit(backbone[backbone$bin == bin, "Partners.trans"],"::")[[1]])), #Partners in trans
                   "ppis.trans" = intersect(tissue.genes,as.character(strsplit(backbone[backbone$bin == bin, "PPIs.trans"],"::")[[1]])), #PPIs in trans
                   "all.int.trans" = intersect(tissue.genes,union(as.character(strsplit(backbone[backbone$bin == bin, "Partners.trans"],"::")[[1]]),
                                                                  as.character(strsplit(backbone[backbone$bin == bin, "PPIs.trans"],"::")[[1]]))) #All interactions
  ) #For each bin, define a new gene sets
  names(gene.set) <- paste(bin, names(gene.set), sep = "_")
  return(gene.set)}

# Function to calculate mean/median expression values
expression.scores <- function(genes){
  df <- exp.data[exp.data$Gene.name %in% genes,]
  
  if(nrow(df) == 0){ #In the case of no partner, PPIs or no covarage by expression data
    res <- data.frame("Cohorts" = tissue,
                      "Mean" = NA,
                      "Median" = NA)}
  else{
    res <- data.frame("Cohorts" = tissue,
                      "Mean" = mean(df[,tissue],na.rm = TRUE),
                      "Median" = median(df[,tissue],na.rm = TRUE))}
  return(res)}

Results.TS <- list() # Do not run this when running additional part

for(level in levels){
  backbone <- backbone_tables_w_features_non_tissue_specific[[level]]
  backbone <- backbone[, colnames(backbone) %in% selected.columns]
  # Remove bins without any gene
  backbone <- backbone[backbone$n_genes != 0,]
  
  for(dataset in names(datasets)){
    if(dataset == "eGTEX.protein"){TS.list <- tissue.specific.proteins}
    else{TS.list <- tissue.specific.genes}
    
    tissue.final  <- c()
    for(tissue in common.tissues){
      
      #Tissue-specific genes
      tissue.genes <- TS.list[[tissue]]$Gene.name
      #Gene sets for each bin
      gene.sets <- mclapply(backbone$bin, gene.lists.per.bin)
      gene.sets.merged <- list()
      for(i in 1:length(gene.sets)){gene.sets.merged <- append(gene.sets.merged, gene.sets[[i]])}
      
      exp.data <- datasets[[dataset]]
      exp.data  <- exp.data[,colnames(exp.data) %in% c("Gene.stable.ID","Gene.name",tissue)]
      res.df <- mclapply(gene.sets.merged, expression.scores, mc.cores = 20)
      final.df <- do.call(rbind,res.df)
      
      if(level == "Arm"){
        final.df$bin <- unlist(lapply(rownames(final.df), function(x) strsplit(x,"_")[[1]][1]))
        final.df$GeneSet <- unlist(lapply(rownames(final.df), function(x) strsplit(x,"_")[[1]][2]))
      }
      else{
        final.df$bin <- unlist(lapply(rownames(final.df), function(x) paste(strsplit(x,"_")[[1]][1],strsplit(x,"_")[[1]][2],sep = "_")))
        final.df$GeneSet <- unlist(lapply(rownames(final.df), function(x) strsplit(x,"_")[[1]][3]))
      }

      final.df$GeneSet <- gsub("\\.\\d+","",final.df$GeneSet)
      final.df$GTExData <- dataset
      tissue.final <- rbind(tissue.final,final.df)
    }
    Results.TS[[level]][[dataset]] <- tissue.final
  }
  print(level)
}

save(Results.TS, file = gtex_tissue_specific)
