# rm(list=ls())
gc(full=T)

library(stringr)
library(dplyr)

genes <- read.delim(file = '/home/ieo5099/mountHD/buffering_events/data/hg19.ncbiRefSeq.gtf', header = F)
# from https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/

if(F){
  file_path <- "/home/ieo5099/mountHD/buffering_events/data/genes.txt"
  num_cycles <- nrow(genes)
  file_conn <- file(file_path, open = "w")
  
  for (i in 1:num_cycles) {
    line_to_write <- paste(i, str_split(genes$V9[i], ';')[[1]][1])
    writeLines(line_to_write, file_conn)
  }
  
  close(file_conn)
  cat("File writing completed.\n")
}

gene_names <- read.table('/home/ieo5099/mountHD/buffering_events/data/genes.txt')
genes$V9 <- NULL

genes$V9 <- gene_names$V3

genes <- genes[genes$V1 %in% paste0('chr', 1:22),]
rm(gene_names)

gene_lengths <- mclapply(levels(factor(genes$V9)), mc.cores = 35, function(i){
  tmp <- genes[genes$V9 == i,]
  return(cbind(i, start = min(tmp$V4), end = max(tmp$V5), chr = levels(factor(tmp$V1))))
})

gene_length <- do.call(rbind, gene_lengths)
rm(gene_lengths, genes)

colnames(gene_length) <- c('Gene', 'Start', 'End', 'Chromosome')
df_genes <- as.data.frame(gene_length)
rm(gene_length)

head(df_genes)
ccds <- read.delim('/home/ieo5099/mountHD/mutation_compensation/data/misc/CCDS.current_hg19.txt')
colnames(ccds)[3] <- 'Gene'

df_genes <- left_join(df_genes, ccds[,c(3,6)])
df_genes_backup <- df_genes

df_genes <- df_genes[df_genes$ccds_status == 'Public',]
df_genes <- df_genes[!is.na(df_genes$ccds_status),]
df_genes$ccds_status <- NULL

df_genes <- df_genes[!duplicated(df_genes),]
rm(ccds)


# Function to assign bins to genes
assign_bins <- function(gene_data, bin_data) {
  overlapping_bins <- mclapply(1:nrow(gene_data), function(i) {
    overlapping_bins <- which(
      gene_data$Chromosome[i] == bin_data$Chromosome &
        gene_data$Start[i] <= bin_data$Bin_End & gene_data$End[i] >= bin_data$Bin_Start
    )
    if (length(overlapping_bins) > 0) {
      overlap_lengths <- numeric(length(overlapping_bins))
      for (j in seq_along(overlapping_bins)) {
        overlap_lengths[j] <- min(gene_data$End[i], bin_data$Bin_End[overlapping_bins[j]]) -
          max(gene_data$Start[i], bin_data$Bin_Start[overlapping_bins[j]])
      }
      return(bin_data$Bin[overlapping_bins[which.max(overlap_lengths)]])
    } else {
      return("Not in Bin")
    }
  }, mc.cores = 30)
  gene_data$Bin <- unlist(overlapping_bins)
  return(gene_data)
}


# Load the ML table 
genes_per_bin <- list()
df_ml <- list()

PATH <- '/home/ieo5099/mountHD/ml_models/'
 
for(segment_length in c(
  1
  ,2,3,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,
  38,40,42,44,46,48
                        )){
  print(segment_length)
  # segment_length <- 4
  df_bin = read.table(paste0(PATH,'data/seg_split/Level_',segment_length,'Mbp.txt'), header = T)
  colnames(df_bin)[c(1,2,3,4)] <- c('Bin', 'Bin_Start', 'Bin_End', 'Chromosome')
  df_bin[,4] <- paste0('chr',df_bin[,4])
  
  df_genes[,2:3] <- apply(df_genes[,2:3], 2, as.numeric)
  
  
  # Assign bins to genes
  result_df <- assign_bins(df_genes, df_bin)
  result_df$Bin <- as.factor(result_df$Bin)
  
  
  # Group genes separating them using a comma3
  grouped_genes <- split(result_df$Gene, result_df$Bin)
  result <- data.frame(
    Bin = names(grouped_genes),
    Grouped_Genes = sapply(grouped_genes, paste, collapse = ", ")
  )
  
  genes_per_bin[[segment_length]] <- grouped_genes
  
  # Calculate haploinsufficient gene density ----
  pLI_scores <- readxl::read_xlsx("/home/ieo5099/mountHD/mutation_compensation/data/misc/pLI_scores.xlsx", sheet = 2)
  pLI_scores_intolerant <- pLI_scores[pLI_scores$pLI >= 0.2,]
  
  haplo_genes <- data.frame()
  for(bin in names(grouped_genes)){
    haplo_genes <- rbind(haplo_genes, 
                         cbind(bin, sum(grouped_genes[[bin]] %in% pLI_scores_intolerant$gene)))
  }
  
  rm(pLI_scores, pLI_scores_intolerant)
  
  # Calculate the number of OGs and TSGs ----
  TSGs <- readxl::read_xlsx("/home/ieo5099/mountHD/mutation_compensation/data/misc/Davoli2013_TUSON.xlsx", sheet = 1)
  OGs <- readxl::read_xlsx("/home/ieo5099/mountHD/mutation_compensation/data/misc/Davoli2013_TUSON.xlsx", sheet = 2)
  
  tsg_genes <- data.frame()
  for(bin in names(grouped_genes)){
    tsg_genes <- rbind(tsg_genes, 
                       cbind(bin, 
                             pancancer = sum(grouped_genes[[bin]] %in% 
                                               TSGs[TSGs$`PAN-Cancer_q-value` <= 0.1,]$Gene),
                             BRCA = sum(grouped_genes[[bin]] %in% 
                                          TSGs[TSGs$`Breast_q-value` <= 0.1,]$Gene),
                             COADREAD = sum(grouped_genes[[bin]] %in% 
                                              TSGs[TSGs$`Colon.Rectum_q-value` <= 0.1,]$Gene),
                             PAAD = sum(grouped_genes[[bin]] %in% 
                                          TSGs[TSGs$`Pancreas_p-value` <= 0.1,]$Gene),
                             GBMLGG = sum(grouped_genes[[bin]] %in% 
                                            TSGs[TSGs$`Glioblastoma_q-value` <= 0.1,]$Gene),
                             LUSC = sum(grouped_genes[[bin]] %in% 
                                          TSGs[TSGs$`Lung.Squamous.Cell_q-value` <= 0.1,]$Gene),
                             LUAD = sum(grouped_genes[[bin]] %in% 
                                          TSGs[TSGs$`Lung.Adenocarcinoma_q-value` <= 0.1,]$Gene)
                       ))
  }
  
  og_genes <- data.frame()
  for(bin in names(grouped_genes)){
    og_genes <- rbind(og_genes, 
                      cbind(bin, 
                            pancancer = sum(grouped_genes[[bin]] %in% 
                                              OGs[OGs$`PAN-Cancer_q-value` <= 0.1,]$Gene),
                            BRCA = sum(grouped_genes[[bin]] %in% 
                                         OGs[OGs$`Breast_q-value` <= 0.1,]$Gene),
                            COADREAD = sum(grouped_genes[[bin]] %in% 
                                             OGs[OGs$`Colon.Rectum_q-value` <= 0.1,]$Gene),
                            PAAD = sum(grouped_genes[[bin]] %in% 
                                         OGs[OGs$`Pancreas_p-value` <= 0.1,]$Gene),
                            GBMLGG = sum(grouped_genes[[bin]] %in% 
                                           OGs[OGs$`Glioblastoma_q-value` <= 0.1,]$Gene),
                            LUSC = sum(grouped_genes[[bin]] %in% 
                                         OGs[OGs$`Lung.Squamous.Cell_q-value` <= 0.1,]$Gene),
                            LUAD = sum(grouped_genes[[bin]] %in% 
                                         OGs[OGs$`Lung.Adenocarcinoma_q-value` <= 0.1,]$Gene)
                      ))
  }
  
  
  tsg_genes_potency <- data.frame()
  for(bin in names(grouped_genes)){
    tsg_genes_potency <- rbind(tsg_genes_potency, 
                               cbind(bin, 
                                     pancancer = sum(-log10(TSGs[TSGs$`PAN-Cancer_q-value` <= 0.1 &
                                                                   TSGs$Gene %in% grouped_genes[[bin]],]$`PAN-Cancer_q-value`)),
                                     BRCA = sum(-log10(TSGs[TSGs$`Breast_q-value` <= 0.1 &
                                                              TSGs$Gene %in% grouped_genes[[bin]],]$`Breast_q-value`)),
                                     COADREAD = sum(-log10(TSGs[TSGs$`Colon.Rectum_q-value` <= 0.1 &
                                                                  TSGs$Gene %in% grouped_genes[[bin]],]$`Colon.Rectum_q-value`)),
                                     PAAD = sum(-log10(TSGs[TSGs$`Pancreas_q-value` <= 0.1 &
                                                              TSGs$Gene %in% grouped_genes[[bin]],]$`Pancreas_q-value`)),
                                     GBMLGG = sum(-log10(TSGs[TSGs$`Glioblastoma_q-value` <= 0.1 &
                                                                TSGs$Gene %in% grouped_genes[[bin]],]$`Glioblastoma_q-value`)),
                                     LUSC = sum(-log10(TSGs[TSGs$`Lung.Squamous.Cell_q-value` <= 0.1 &
                                                              TSGs$Gene %in% grouped_genes[[bin]],]$`Lung.Squamous.Cell_q-value`)),
                                     LUAD = sum(-log10(TSGs[TSGs$`Lung.Adenocarcinoma_q-value` <= 0.1 &
                                                              TSGs$Gene %in% grouped_genes[[bin]],]$`Lung.Adenocarcinoma_q-value`))
                               ))
  }
  
  og_genes_potency <- data.frame()
  for(bin in names(grouped_genes)){
    og_genes_potency <- rbind(og_genes_potency, 
                              cbind(bin, 
                                    pancancer = sum(-log10(OGs[OGs$`PAN-Cancer_q-value` <= 0.1 &
                                                                 OGs$Gene %in% grouped_genes[[bin]],]$`PAN-Cancer_q-value`)),
                                    BRCA = sum(-log10(OGs[OGs$`Breast_q-value` <= 0.1 &
                                                            TSGs$Gene %in% grouped_genes[[bin]],]$`Breast_q-value`)),
                                    COADREAD = sum(-log10(OGs[OGs$`Colon.Rectum_q-value` <= 0.1 &
                                                                OGs$Gene %in% grouped_genes[[bin]],]$`Colon.Rectum_q-value`)),
                                    PAAD = sum(-log10(OGs[OGs$`Pancreas_q-value` <= 0.1 &
                                                            OGs$Gene %in% grouped_genes[[bin]],]$`Pancreas_q-value`)),
                                    GBMLGG = sum(-log10(OGs[OGs$`Glioblastoma_q-value` <= 0.1 &
                                                              OGs$Gene %in% grouped_genes[[bin]],]$`Glioblastoma_q-value`)),
                                    LUSC = sum(-log10(OGs[OGs$`Lung.Squamous.Cell_q-value` <= 0.1 &
                                                            OGs$Gene %in% grouped_genes[[bin]],]$`Lung.Squamous.Cell_q-value`)),
                                    LUAD = sum(-log10(OGs[OGs$`Lung.Adenocarcinoma_q-value` <= 0.1 &
                                                            OGs$Gene %in% grouped_genes[[bin]],]$`Lung.Adenocarcinoma_q-value`))
                              ))
  }
  
  rm(OGs, TSGs)
  
  
  # Essential (common and tissue-specific) ----
  common_essential <- read.table('/home/ieo5099/mountHD/mutation_compensation/data/misc/CRISPR_common_essentials.csv', header = T)
  load("/home/ieo5099/mountHD/mutation_compensation/data/misc/mean_crispr_effect.RData")
  
  essential_genes <- data.frame()
  for(bin in names(grouped_genes)){
    essential_genes <- rbind(essential_genes, 
                             cbind(bin, 
                                   essential_pancancer = sum(common_essential$gene %in% grouped_genes[[bin]]),
                                   essential_BRCA = sum(mean_crispr_effect[order(
                                     mean_crispr_effect$`Breast Cancer`, decreasing = F),][1:250,]$genes %in% grouped_genes[[bin]]),
                                   essential_COADREAD = sum(mean_crispr_effect[order(
                                     mean_crispr_effect$`Colon/Colorectal Cancer`, decreasing = F),][1:250,]$genes %in% grouped_genes[[bin]]),
                                   essential_PAAD = sum(mean_crispr_effect[order(
                                     mean_crispr_effect$`Pancreatic Cancer`, decreasing = F),][1:250,]$genes %in% grouped_genes[[bin]]),
                                   essential_GBMLGG = sum(mean_crispr_effect[order(
                                     mean_crispr_effect$`Brain Cancer`, decreasing = F),][1:250,]$genes %in% grouped_genes[[bin]]),
                                   essential_LUSC = sum(mean_crispr_effect[order(
                                     mean_crispr_effect$`Lung Cancer`, decreasing = F),][1:250,]$genes %in% grouped_genes[[bin]]),
                                   essential_LUAD = sum(mean_crispr_effect[order(
                                     mean_crispr_effect$`Lung Cancer`, decreasing = F),][1:250,]$genes %in% grouped_genes[[bin]]))
    )
  }
  rm(common_essential)
  
  
  # add features to df_bin (ML table) ----
  
  # gene_count
  gene_count2 <- as.data.frame(sapply(grouped_genes, length))
  gene_count2$Bin <- rownames(gene_count2)
  rownames(gene_count2) <- NULL
  colnames(gene_count2)[1] <- 'gene_count2'
  
  df_bin2 <- full_join(gene_count2, df_bin, by = 'Bin')
  
  df_bin2 <- df_bin2[!is.na(df_bin2$gene_count2),]
  df_bin2 <- df_bin2[!df_bin2$Bin == 'Not in Bin',]
  
  # OGs
  og_genes <- og_genes[!og_genes$bin == 'Not in Bin',]
  colnames(og_genes)[1] <- 'Bin'
  colnames(og_genes)[2:8] <- paste0('OGscore_',colnames(og_genes)[2:8])
  df_bin2 <- full_join(og_genes, df_bin2)
  
  og_genes_potency <- og_genes_potency[!og_genes_potency$bin == 'Not in Bin',]
  colnames(og_genes_potency)[1] <- 'Bin'
  colnames(og_genes_potency)[2:8] <- paste0('OGscore.potency_',colnames(og_genes_potency)[2:8])
  df_bin2 <- full_join(og_genes_potency, df_bin2)
  
  # TSGs
  tsg_genes <- tsg_genes[!tsg_genes$bin == 'Not in Bin',]
  colnames(tsg_genes)[1] <- 'Bin'
  colnames(tsg_genes)[2:8] <- paste0('TSGscore_',colnames(tsg_genes)[2:8])
  df_bin2 <- full_join(tsg_genes, df_bin2)
  
  tsg_genes_potency <- tsg_genes_potency[!tsg_genes_potency$bin == 'Not in Bin',]
  colnames(tsg_genes_potency)[1] <- 'Bin'
  colnames(tsg_genes_potency)[2:8] <- paste0('TSGscore.potency_',colnames(tsg_genes_potency)[2:8])
  df_bin2 <- full_join(tsg_genes_potency, df_bin2)
  
  # haploinsufficient genes
  head(haplo_genes)
  colnames(haplo_genes) <- c('Bin', 'haploinsufficient_genes')
  df_bin2 <- left_join(df_bin2, haplo_genes, by = 'Bin')
  df_bin2$haploinsufficient_genes <- as.integer(df_bin2$haploinsufficient_genes)
  
  df_bin2$haploinsuff_density <- df_bin2$haploinsufficient_genes/df_bin2$gene_count2
  
  # essential genes
  head(essential_genes)
  colnames(essential_genes)[1] <- 'Bin'
  df_bin2 <- left_join(df_bin2, essential_genes)
  
  # df_bin2[,c(2:30,61:68)] <- apply(df_bin2[,c(2:30,61:68)], 2, as.numeric)
  
  
  colnames(df_bin2)
  df_bin2$Chromosome <- parse_number(df_bin2$Chromosome)
  
  df_ml[[segment_length]] <- df_bin2
}

saveRDS(genes_per_bin, '/home/ieo5099/mountHD/ml_models/data/genes_per_bin.rds')
saveRDS(df_ml, '/home/ieo5099/mountHD/ml_models/data/ml_tables.rds')
