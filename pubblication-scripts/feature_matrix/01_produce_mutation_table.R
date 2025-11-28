# BIN-LEVEL ANALYSIS
#
# This script produces as output mutation score and amplification frequency for
# each chromosome (for each tumor type)

# rm(list=ls())
gc(full=T)

suppressMessages({
  library(readxl)
  library(tibble)
  library(ggplot2)
  library(tidyr)
  library(readr)
  library(stringr)
  library(data.table)
  library(dplyr)
  library(ggExtra)
  library(crayon)
  library(parallel)
})


# setwd("../")
setwd("/mnt/fabiogokce/fabiohd/mutation_compensation/")

tumor_types <- c(
  "BRCA",
  "LUAD",
  "LUSC",
  "CESC",
  "THCA",
  "HNSC",
  "PAAD",
  "COADREAD",
  "GBMLGG",
  # # #
  "SKCM",
  "BLCA",
  "PCPG",
  "PRAD",
  "KIRC",
  "MESO",
  "TGCT",
  "KIRP",
  "SARC",
  "LIHC",
  "ESCA",
  "STAD",
  "UCS",
  "OV"
)


# 5 10 20 30
segment_cutoff <- 20
cat("\n\n segment cutoff:", segment_cutoff, "\n\n")
stringent_mutations <- F

cores <- 10

fixed_bin_length <- 500000

# for (tumor_type in tumor_types) {
mclapply(tumor_types, mc.cores = cores, function(tumor_type){
  ## 1st Loop: load sCNAs and SNVs files ----
  cat(" \n\n\n >>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<< \n")
  cat(" >     ----------", paste0(tumor_type, " ---------- \n"))
  
  # Load SNVs file
  cat("\n > Load", tumor_type, "Mutation File \n")
  snv <-
    read.csv(file = paste0("data/FireBrowse_SNVs/", tumor_type, "_mutations.csv"))
  snv <-
    snv[!duplicated(snv[, c(6, 7, 8, 9, 17)]),] # remove duplicated mutations in the same patient
  # cat("  done!")
  
  # Load CNAs file
  cat("\n > Load", tumor_type, "SCNA File \n")
  scna <-
    read.delim(
      paste0(
        "data/FireBrowse_CNAs/",
        tumor_type,
        "_hg19_FireBrowse_totalCopyNumber_classifiedALL.tsv"
      ),
      header = T
    )
  # keep patients with both snv and scna data
  snv_filt <-
    snv[snv$patient_id %in% levels(factor(scna$patient_id)),]
  # keep patients with both snv and scna data
  scna <-
    scna[scna$patient_id %in% snv$patient_id,]
  
  common_patients <-
    levels(as.factor(snv[snv$patient_id %in% scna$patient_id,]$patient_id))
  # keep only non copy-neutral segments (amplified or deleted)
  scna <-
    scna[scna$Segment_Mean >= as.numeric(segment_cutoff)/100 | 
           scna$Segment_Mean <= -as.numeric(segment_cutoff)/100,]
  
  rm(snv)
  gc(full = T)
  
  
  ##  initialize table for chromosome and arm amplification frequencies ----
  cat("\n\n >> Start", tumor_type, "chromosome Loop \n")
  
  chr_df <- data.frame()
  ## >> 2nd Loop: chromosome << ----
  for (chr in 1:22) {
    gc(full = TRUE)
    # print(chr)
    
    ## 2nd loop: load chromosome region lengths ----
    chr_info <-
      read.table("data/misc/chr_info_h19.txt", header = TRUE)
    chr_arms <-
      read.table(file = "data/misc/cytoBand.txt", header = T)
    chr_arms[, 2:3] <-
      apply(chr_arms[, 2:3] / fixed_bin_length, 2, as.integer)
    
    all_patients <-
      levels(as.factor(as.character(scna$patient_id)))
    
    ## 2nd loop: create temporary tables for SNVs and CNAs ----
    temp_cna <- scna[scna$Chromosome == chr,]
    temp_snv <- snv_filt[snv_filt$Chromosome == chr,]
    
    ## 2nd loop: define binning parameters ----
    n_bins <-
      as.integer(chr_info[chr_info$Chromosome == paste0("chr", chr),]$Length /
                   fixed_bin_length)
    length_bin <-
      as.integer(chr_info[chr_info$Chromosome == paste0("chr", chr),]$Length /
                   n_bins)
    
    ## 2nd loop: compute start and length for each CNA according to bins ----
    end <- length_bin
    start_bin <- NULL
    temp_cna_length <- NULL
    i <- 1
    
    # this function assign each CNA event to segments 
    whichbin <- function(data, end, i) {
      if (any(data$Start < end) == TRUE) {
        temp2 <- data[data$Start < end,]
        temp2 <-
          cbind(
            temp2,
            start_bin = rep(i, nrow(temp2)),
            bin_length = (temp2$End - temp2$Start) / fixed_bin_length
          )
      }
      return(temp2)
    }
    
    # apply the whichbin() function
    for (i in 1:n_bins) {
      if (nrow(temp_cna) == 0) {
        break
      }
      if (any(temp_cna$Start < end)) {
        temp_cna_length <-
          rbind(temp_cna_length, whichbin(temp_cna, end, i))
      }
      temp_cna <- temp_cna[!temp_cna$Start < end,]
      end <- end + length_bin
    }
    
    # keep CNA events higher that 0,5 Mbp (half of bin size)
    temp_cna_length <-
      temp_cna_length[temp_cna_length$length > fixed_bin_length/2,]
    amplified_patients <-
      levels(factor(temp_cna_length$patient_id))
    
    ## 2nd loop: load segmented chromosome/gene structure ----
    bin_gene <- read_rds(file = '/mnt/fabiogokce/fabiohd/ml_models/data/misc/chr_backbone_wSmallerBins.rds')
    bin_gene <- bin_gene[[paste0(fixed_bin_length/1000000,'Mbp')]][[chr]]
    
    # bin_gene <-
    #   read.table(
    #     paste0(
    #       "data/ChromosomeGeneStructure/chr_",
    #       chr,
    #       "_binSize_",
    #       fixed_bin_length,
    #       ".txt"
    #     )
    #   )
    
    ## 2nd loop: compute amplification and deletion frequencies ----
    chr_bins <- data.frame(bins = 1:n_bins)
    
    # temp_cna_length <- temp_cna_length_backup
    
    for (i in 1:nrow(temp_cna_length)) {
      if (temp_cna_length$bin_length[i] > 0.5) {
        y <-
          temp_cna_length$bin_length[i] - as.integer(temp_cna_length$bin_length[i])
        if (y >= 0.5) {
          x <-
            data.frame(
              bins = temp_cna_length$start_bin[i]:as.integer(
                temp_cna_length$start_bin[i] + temp_cna_length$bin_length[i]
              ),
              c(rep(
                1, as.integer(temp_cna_length$bin_length[i] + 1)
              ))
            )
        } else{
          x <-
            data.frame(
              bins = temp_cna_length$start_bin[i]:as.integer(
                temp_cna_length$start_bin[i] + temp_cna_length$bin_length[i] - 1
              ),
              c(rep(
                1, as.integer(temp_cna_length$bin_length[i])
              ))
            )
        }
      } else{
        next
      }
      
      chr_bins <- full_join(chr_bins, x, by = "bins")
      
      colnames(chr_bins) <-
        c(colnames(chr_bins)[-length(colnames(chr_bins))], paste0("CNA", i))
    }
    rm(i)
    rm(x)
    rm(y)
    
    chr_bins <- t(as.matrix(chr_bins))
    chr_bins <- as.data.frame(chr_bins)
    chr_bins <- chr_bins[-1,]
    
    # CNA shorter than half of the bin size (0.5 Mbp in this case) were removed from calculation
    chr_bins <-
      cbind(
        Segment_Mean = temp_cna_length$Segment_Mean,
        classified = temp_cna_length$classified,
        patient = temp_cna_length$patient_id,
        chr_bins
      )
    
    # chr_bins[is.na(chr_bins)] <- 0
    
    cna_freq_ampl <-
      as.numeric(colSums(chr_bins[chr_bins$Segment_Mean >= as.numeric(segment_cutoff)/100, 4:ncol(chr_bins)], na.rm = TRUE) /
                   length(common_patients))
    cna_freq_del <-
      as.numeric(colSums(chr_bins[chr_bins$Segment_Mean <= -as.numeric(segment_cutoff)/100, 4:ncol(chr_bins)], na.rm = TRUE) /
                   length(common_patients))
    cna_freq_total <-
      as.numeric(colSums(chr_bins[, 4:ncol(chr_bins)], na.rm = TRUE) / length(common_patients))
    
    
    unamplified_patients <-
        all_patients[!(all_patients %in% chr_bins$patient)]
    
    ## 2nd loop: compute the mutation score in copy-neutral regions ----
    # (normalized for the number of patients (1), the coding region (2) and thus log10 (3))
    chr_bins[is.na(chr_bins)] <- 0
    
    # write_tsv(chr_bins, file = paste0("/home/ieo5099/mountHD/ml_models/data/CNA_0-1_tables/",
                                      # tumor_type,"_chr",chr,"_",fixed_bin_length,"bins.tsv"))
    
    
    chr_bins_pt <- data.frame()
    for (pt in levels(factor(chr_bins$patient))) {
      # colSums(chr_bins[chr_bins$patient == pt,][,-c(1:3)])
      chr_bins_pt <-
        rbind(chr_bins_pt, c(pt, as.numeric(colSums(chr_bins[chr_bins$patient == pt,][,-c(1:3)]))))
    }
    colnames(chr_bins_pt) <- c("patients", 1:n_bins)
    
    x <- NULL
    end <- length_bin
    start_bin <- 0
    mutations_raw <- NULL
    mutations_norm <- NULL
    n_patients <- NULL
    
    mut_withinCNA <- F
    for (i in 1:n_bins) {
      if(mut_withinCNA == T){
        n_pts <- length(all_patients)
        mut <- temp_snv[str_sub(temp_snv$Tumor_Sample_Barcode, end = 12) %in% all_patients,]
        mut <-
          mut[mut$Start_Position >= start_bin &
                mut$End_Position < end,]
        mut_a <- as.data.frame(table(mut$patient_id))
      }else{
        n_pts <-
          length(c(chr_bins_pt[chr_bins_pt[, i + 1] == 0,]$patient, unamplified_patients))
        mut <-
          temp_snv[str_sub(temp_snv$Tumor_Sample_Barcode, end = 12) %in%
                     c(chr_bins_pt[chr_bins_pt[, i +1] == 0,]$patient,
                       unamplified_patients),]
        mut <-
          mut[mut$Start_Position >= start_bin &
                mut$End_Position < end,]
        mut_a <- as.data.frame(table(mut$patient_id))
      }
      
      if (dim(mut_a)[1] != 0) {
        colnames(mut_a) <- c("patient_id", "mutations_raw")
        x <- sum(mut_a$mutations_raw , na.rm = T)
        mutations_raw <- c(mutations_raw, x)
        mutations_norm <- c(mutations_norm, x / n_pts)
      } else{
        mutations_raw <- c(mutations_raw, 0)
        mutations_norm <- c(mutations_norm, 0)
      }
      n_patients <- c(n_patients, n_pts)
      start_bin <- start_bin + length_bin
      end <- end + length_bin
    }
    
    bin_gene_mut <-
      cbind(bin_gene, mutations_raw, mutations_norm, n_patients)
    bin_gene_mut <- as.data.frame(bin_gene_mut)
    
    # end of the 2nd loop: write chromosome table ----
    
    chr_df <- rbind(chr_df, bin_gene_mut)
  }
  write.table(
    chr_df,
    file = paste0('/mnt/fabiogokce/fabiohd/ml_models/data/mutation_tables/',
                  tumor_type,
                  "_",
                  fixed_bin_length/1000000,
                  "Mbp.txt"
    ), 
    col.names = T,
    row.names = F, 
    quote = F
  )
  
})


