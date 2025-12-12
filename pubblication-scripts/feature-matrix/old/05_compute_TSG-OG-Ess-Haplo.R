# rm(list=ls())
gc(full=T)

library(tidyverse)
  
setwd('/mnt/fabiogokce/fabiohd/')

chr_backbone_genes <- readRDS("ml_models/data/misc/chr_backbone_genes.rds")
genes_per_bin <- readRDS("ml_models/data/genes_per_bin.rds")

for(segment_length in c(0.1,0.25,0.5,1:50)){
  print(segment_length)
  
  grouped_genes <- genes_per_bin[[paste0(segment_length,'Mbp')]]
  
  # Calculate haploinsufficient gene density ----
  pLI_scores <- readxl::read_xlsx("mutation_compensation/data/misc/pLI_scores.xlsx", sheet = 2)
  pLI_scores_intolerant <- pLI_scores[pLI_scores$pLI >= 0.2,]
  
  haplo_genes <- data.frame()
  for(bin in names(grouped_genes)){
    haplo_genes <- rbind(haplo_genes, 
                         cbind(bin, sum(grouped_genes[[bin]] %in% pLI_scores_intolerant$gene)))
  }
  
  rm(pLI_scores, pLI_scores_intolerant)
  
  
  # Calculate the number of OGs and TSGs ----
  TSGs <- readxl::read_xlsx("mutation_compensation/data/misc/Davoli2013_TUSON.xlsx", sheet = 1)
  OGs <- readxl::read_xlsx("mutation_compensation/data/misc/Davoli2013_TUSON.xlsx", sheet = 2)
  
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
  common_essential <- read.table('mutation_compensation/data/misc/CRISPR_common_essentials.csv', header = T)
  load("mutation_compensation/data/misc/mean_crispr_effect.RData")
  
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
  
  
  # OGs
  og_genes <- og_genes[!og_genes$bin == 'Not in Bin',]
  colnames(og_genes)[2:8] <- paste0('OGscore_',colnames(og_genes)[2:8])
  chr_backbone_genes[[paste0(segment_length,'Mbp')]] <- 
    full_join(chr_backbone_genes[[paste0(segment_length,'Mbp')]], og_genes)
  
  og_genes_potency <- og_genes_potency[!og_genes_potency$bin == 'Not in Bin',]
  colnames(og_genes_potency)[2:8] <- paste0('OGscore.potency_',colnames(og_genes_potency)[2:8])
  chr_backbone_genes[[paste0(segment_length,'Mbp')]] <- 
    full_join(chr_backbone_genes[[paste0(segment_length,'Mbp')]], og_genes_potency)
  
  # TSGs
  tsg_genes <- tsg_genes[!tsg_genes$bin == 'Not in Bin',]
  colnames(tsg_genes)[2:8] <- paste0('TSGscore_',colnames(tsg_genes)[2:8])
  chr_backbone_genes[[paste0(segment_length,'Mbp')]] <- 
    full_join(chr_backbone_genes[[paste0(segment_length,'Mbp')]], tsg_genes)
  
  tsg_genes_potency <- tsg_genes_potency[!tsg_genes_potency$bin == 'Not in Bin',]
  colnames(tsg_genes_potency)[2:8] <- paste0('TSGscore.potency_',colnames(tsg_genes_potency)[2:8])
  chr_backbone_genes[[paste0(segment_length,'Mbp')]] <- 
    full_join(chr_backbone_genes[[paste0(segment_length,'Mbp')]], tsg_genes_potency)
  
  # haploinsufficient genes
  colnames(haplo_genes) <- c('bin', 'haploinsufficient_genes')
  chr_backbone_genes[[paste0(segment_length,'Mbp')]] <- 
    full_join(chr_backbone_genes[[paste0(segment_length,'Mbp')]], haplo_genes)
  
  # essential genes
  colnames(essential_genes)[1] <- 'bin'
  chr_backbone_genes[[paste0(segment_length,'Mbp')]] <- 
    full_join(chr_backbone_genes[[paste0(segment_length,'Mbp')]], essential_genes)
}

saveRDS(chr_backbone_genes, file = 'ml_models/data/misc/chr_backbone_genes_cancer.rds')
