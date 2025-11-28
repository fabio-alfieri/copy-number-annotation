# rm(list=ls())
gc(full=T)

library(stringr)
library(dplyr)

genes <- read.delim(file = '/mnt/fabiogokce/fabiohd/ml_models/data/misc/hg19.ncbiRefSeq.gtf', header = F)
# from https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/

if(F){
  file_path <- "/mnt/fabiogokce/fabiohd/ml_models/data/misc/genes.txt"
  num_cycles <- nrow(genes)
  file_conn <- file(file_path, open = "w")
  
  for (i in 1:num_cycles) {
    line_to_write <- paste(i, str_split(genes$V9[i], ';')[[1]][1])
    writeLines(line_to_write, file_conn)
  }
  
  close(file_conn)
}

gene_names <- read.table('/mnt/fabiogokce/fabiohd/ml_models/data/misc/genes.txt')
genes$V9 <- NULL

genes$V9 <- gene_names$V3

genes <- genes[genes$V1 %in% paste0('chr', 1:22),]
rm(gene_names)

gene_lengths <- mclapply(levels(factor(genes$V9)), mc.cores = 48, function(i){
  tmp <- genes[genes$V9 == i,]
  return(cbind(i, start = min(tmp$V4), end = max(tmp$V5), chr = levels(factor(tmp$V1))))
})

gene_length <- do.call(rbind, gene_lengths)
rm(gene_lengths, genes)

colnames(gene_length) <- c('Gene', 'Start', 'End', 'Chromosome')
df_genes <- as.data.frame(gene_length)
rm(gene_length)

head(df_genes)
ccds <- read.delim('/mnt/fabiogokce/fabiohd/mutation_compensation/data/misc/CCDS.current_hg19.txt')
colnames(ccds)[3] <- 'Gene'

df_genes <- left_join(df_genes, ccds[,c(3,6)])
df_genes_backup <- df_genes

df_genes <- df_genes[df_genes$ccds_status == 'Public',]
df_genes <- df_genes[!is.na(df_genes$ccds_status),]
df_genes$ccds_status <- NULL

df_genes <- df_genes[!duplicated(df_genes),]
rm(ccds)

df_genes$chr <- parse_number(df_genes$chr)

# Function to assign bins to genes
# assign_bins <- function(gene_data, bin_data) {
#   overlapping_bins <- mclapply(1:nrow(gene_data), function(i) {
#     overlapping_bins <- which(
#       gene_data$chr[i] == bin_data$chr &
#         gene_data$Start[i] <= bin_data$end_bin & gene_data$End[i] >= bin_data$start_bin
#     )
#     if (length(overlapping_bins) > 0) {
#       overlap_lengths <- numeric(length(overlapping_bins))
#       # length_coding <- sum(overlapping_bins$End-overlapping_bins$Start)
#       for (j in seq_along(overlapping_bins)) {
#         overlap_lengths[j] <- min(gene_data$End[i], bin_data$end_bin[overlapping_bins[j]]) -
#           max(gene_data$Start[i], bin_data$start_bin[overlapping_bins[j]])
#       }
#       return(cbind(bin = bin_data$bin[overlapping_bins[which.max(overlap_lengths)]],
#                    overlap_length = overlap_lengths[j]))
#     } else {
#       return(cbind(bin = "Not in Bin",
#                    overlap_length = 0))
#     }
#   }, mc.cores = 30)
#   gene_data <- cbind(gene_data, do.call(rbind,overlapping_bins))
#   return(gene_data)
# }

assign_bins <- function(gene_data, bin_data) {
  overlapping_bins <- mclapply(1:nrow(gene_data), function(i) {
    overlapping_bins <- which(
      gene_data$chr[i] == bin_data$chr &
        gene_data$Start[i] <= bin_data$end_bin & gene_data$End[i] >= bin_data$start_bin
    )
    if (length(overlapping_bins) > 0) {
      overlap_lengths <- numeric(length(overlapping_bins))
      for (j in seq_along(overlapping_bins)) {
        overlap_start <- max(gene_data$Start[i], bin_data$start_bin[overlapping_bins[j]])
        overlap_end <- min(gene_data$End[i], bin_data$end_bin[overlapping_bins[j]])
        overlap_length <- max(0, overlap_end - overlap_start + 1)  # Length of overlap
        coding_overlap <- min(gene_data$End[i], bin_data$end_bin[overlapping_bins[j]]) -
          max(gene_data$Start[i], bin_data$start_bin[overlapping_bins[j]])
        overlap_lengths[j] <- coding_overlap
      }
      return(cbind(bin = bin_data$bin[overlapping_bins[which.max(overlap_lengths)]],
                   coding_length = sum(overlap_lengths)))
    } else {
      return(cbind(bin = "Not in Bin",
                   coding_length = 0))
    }
  }, mc.cores = 30)
  gene_data <- cbind(gene_data, do.call(rbind, overlapping_bins))
  return(gene_data)
}

# Load the ML table 
genes_per_bin <- list()
chr_backbone_genes <- list()

PATH <- '/mnt/fabiogokce/fabiohd/ml_models/'
mut_cna_features <- readRDS(paste0(PATH,"data/misc/mut_cna_features.rds"))

for(segment_length in c(0.1,0.25,0.5,1:50)){
  print(segment_length)
  
  chr_backbone <- readRDS(paste0(PATH,"/data/misc/chr_backbone_wSmallerBins.rds"))[[paste0(segment_length,'Mbp')]]
  chr_backbone <- as.data.frame(do.call(rbind, chr_backbone))
  colnames(chr_backbone)[3:4] <- c('start_bin', 'end_bin')
  chr_backbone[,1:4] <- apply(chr_backbone[,1:4], 2, as.numeric)
  
  df_genes[,2:3] <- apply(df_genes[,2:3], 2, as.numeric)
  colnames(df_genes)[4] <- 'chr'
  
  # Assign bins to genes
  result_df <- assign_bins(df_genes, chr_backbone)
  
  result_df$bin <- paste0(result_df$chr,'_',result_df$bin)
  
  coding_length <- data.frame()
  for(bin in levels(factor(result_df$bin))){
    coding_length <- rbind(
      coding_length,
      cbind(bin = bin,
            coding_length = sum(as.numeric(result_df[result_df$bin == bin,]$coding_length)),
            n_genes = length(result_df[result_df$bin == bin,]$coding_length))
    )
  }
  
  chr_backbone$bin <- paste0(chr_backbone$chr,'_',chr_backbone$bin)
  chr_backbone <- full_join(chr_backbone, coding_length, by = 'bin')
  chr_backbone[is.na(chr_backbone)] <- 0
  chr_backbone[,c(1,3:6)] <- apply(chr_backbone[,c(1,3:6)], 2, as.numeric)
  
  chr_backbone$coding_length <- chr_backbone$coding_length/(chr_backbone$end_bin-chr_backbone$start_bin)
  chr_backbone$coding_length <- ifelse(chr_backbone$coding_length>1,1,chr_backbone$coding_length)
  # Group genes separating them using a comma3
  grouped_genes <- split(result_df$Gene, result_df$bin)
  
  genes_per_bin[[paste0(segment_length,'Mbp')]] <- grouped_genes
  chr_backbone_genes[[paste0(segment_length,'Mbp')]] <- chr_backbone
  
}

saveRDS(genes_per_bin, file = paste0(PATH,'data/genes_per_bin.rds'))
saveRDS(chr_backbone_genes, file = paste0(PATH,'/data/misc/chr_backbone_genes.rds'))
