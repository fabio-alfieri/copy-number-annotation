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
setwd("~/mountHD/mutation_compensation/")

tumor_types <- c(
  "BRCA",
  "LUAD",
  "LUSC",
  "HNSC",
  "PAAD",
  "COADREAD",
  "GBMLGG"
)

mutation_types <- c(
  "all_mutations"
)

# 5 10 20 30
segment_cutoff <- 20
cat("\n\n segment cutoff:", segment_cutoff, "\n\n")
stringent_mutations <- F

cores <- 10
# cores <- length(tumor_types)-8

# produce HAPLOINSUFFICIENT_GHIS score table
# ghis <- as.data.frame(readxl::read_xlsx("data/misc/nar-03716-met-n-2014-File006.xlsx", sheet = 3))
# colnames(ghis) <- c("ENSEMBL", "GHIS")
# library(org.Hs.eg.db)
# genes <- as.vector(ghis[,1])
# annots <- select(org.Hs.eg.db, keys=genes,
#                  columns="SYMBOL", keytype="ENSEMBL")
# result <- merge(ghis, annots, by.x="ENSEMBL", by.y="ENSEMBL")
# colnames(result) <- c("ENSEMBL", "GHIS", "Hugo_Symbol")
# write.table(result, file = "data/misc/ghis_scores.tsv", sep = "\t", quote = F, row.names = F)

fixed_bin_length <- 1000000

for(fixed_bin_length in c(1,2
                          ,4,6,8,10,12,14,16,18,20,22,24,26,
                          28,30,32,34,36,38,40,42,44,46,48,50
)*1000000){
  for (mutation_type in mutation_types) {
    
    mclapply(tumor_types, mc.cores = cores, function(tumor_type){
      
      # Load SNV file
      snv <-
        read.csv(file = paste0("data/FireBrowse_SNVs/", tumor_type, "_mutations.csv"))
      snv <-
        snv[!duplicated(snv[, c(6, 7, 8, 9, 17)]),] # remove duplicated mutations in the same patient
      
      # Load CNA file
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
      scna$length <- scna$length+1
      
      rm(snv)
      gc(full = T)
      
      
      ##  initialize table for chromosome and arm amplification frequencies ----
      for (chr in 1:22) {
        gc(full = TRUE)
        
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
        amplified_patients <-
          levels(factor(temp_cna_length$patient_id))
        
        ## 2nd loop: load segmented chromosome/gene structure ----
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

        bin_gene_new <- 
          read.table(
            paste0(
              '/home/ieo5099/mountHD/ml_models/backbone_tables/Level_',
              fixed_bin_length/1000000,'Mbp_gendist.txt'
            ), header = T
          )
        colnames(bin_gene_new)[2:3] <- c('bin_start', 'bin_end')
        bin_gene_new <- bin_gene_new[bin_gene_new$chr == chr,]
        
        ## 2nd loop: compute chromosome and arm amplification frequencies ----
        temp_cna_length_backup <- temp_cna_length
        temp_cna_length <-
          temp_cna_length %>% filter(Segment_Mean >= as.numeric(segment_cutoff)/100)
        
        altered_bins <- lapply(1:nrow(bin_gene_new), function(i){
          # print(i)
          
          start <- bin_gene_new[i,]$bin_start
          end <- bin_gene_new[i,]$bin_end
          bin_id <- bin_gene_new[i,]$bin
          
          in_bin <- temp_cna_length[!(temp_cna_length$Start < start & 
                                        temp_cna_length$End < start),]
          in_bin <- in_bin[!(in_bin$End > end &
                               in_bin$Start > end),]
          
          all_bin <- in_bin[in_bin$Start <= start & in_bin$End >= end,]
          partial_bin <- in_bin[!(in_bin$Start <= start & in_bin$End >= end),]
          
          if(nrow(all_bin)==0 & nrow(partial_bin)==0){return(NA)}
          
          small_bin <- partial_bin[partial_bin$length < 1000000,]
          left_bin <- partial_bin[partial_bin$End > start & partial_bin$End < end,]
          right_bin <- partial_bin[partial_bin$Start > start & partial_bin$End > end,]
          
          if(nrow(all_bin)>0){
            all <- as.data.frame(cbind(patient_id = all_bin$patient_id, 
                                       bin = (bin_gene_new[i,]$bin_end-bin_gene_new[i,]$bin_start),
                                       bin_id = bin_id))
          }else{
            all <- data.frame()
          }
          if(nrow(small_bin)>0){
            small <- as.data.frame(cbind(patient_id = small_bin$patient_id, 
                                         bin = small_bin$End-small_bin$Start,
                                         bin_id = bin_id))
            all <- rbind(all, small)
          }
          if(nrow(left_bin)>0){
            left <- as.data.frame(cbind(patient_id = left_bin$patient_id, 
                                        bin = left_bin$End-start,
                                        bin_id = bin_id))
            all <- rbind(all, left)
          }
          if(nrow(right_bin)>0){
            right <- as.data.frame(cbind(patient_id = right_bin$patient_id, 
                                         bin = end-right_bin$Start,
                                         bin_id = bin_id))
            all <- rbind(all, right)
          }
          
          return(all)
        })
        
        
        scores_ampl <- data.frame()
        for(i in 1:length(altered_bins)){
          
          if(any(is.na(altered_bins[[i]]))){
            score <- 0
          }else{
            score <- sum(as.numeric(altered_bins[[i]][,2]))/
              (length(common_patients)*(bin_gene_new[i,]$bin_end-bin_gene_new[i,]$bin_start))
          }
          
          scores_ampl <- rbind(scores_ampl, 
                          cbind(bin = bin_gene_new[i,]$bin,
                                ampl_score = score))
        }
        
        
        ## 2nd loop: compute amplification and deletion frequencies ----
        chr_bins <- data.frame(bins = 1:n_bins)
        
        temp_cna_length <-
          temp_cna_length_backup %>% filter(Segment_Mean <= -as.numeric(segment_cutoff)/100)
        
        altered_bins <- lapply(1:nrow(bin_gene_new), function(i){
          # print(i)
          
          start <- bin_gene_new[i,]$bin_start
          end <- bin_gene_new[i,]$bin_end
          bin_id <- bin_gene_new[i,]$bin
          
          in_bin <- temp_cna_length[!(temp_cna_length$Start < start & 
                                        temp_cna_length$End < start),]
          in_bin <- in_bin[!(in_bin$End > end &
                               in_bin$Start > end),]
          
          all_bin <- in_bin[in_bin$Start <= start & in_bin$End >= end,]
          partial_bin <- in_bin[!(in_bin$Start <= start & in_bin$End >= end),]
          
          if(nrow(all_bin)==0 & nrow(partial_bin)==0){return(NA)}
          
          small_bin <- partial_bin[partial_bin$length < 1000000,]
          left_bin <- partial_bin[partial_bin$End > start & partial_bin$End < end,]
          right_bin <- partial_bin[partial_bin$Start > start & partial_bin$End > end,]
          
          if(nrow(all_bin)>0){
            all <- as.data.frame(cbind(patient_id = all_bin$patient_id, 
                                       bin = (bin_gene_new[i,]$bin_end-bin_gene_new[i,]$bin_start),
                                       bin_id = bin_id))
          }else{
            all <- data.frame()
          }
          if(nrow(small_bin)>0){
            small <- as.data.frame(cbind(patient_id = small_bin$patient_id, 
                                         bin = small_bin$End-small_bin$Start,
                                         bin_id = bin_id))
            all <- rbind(all, small)
          }
          if(nrow(left_bin)>0){
            left <- as.data.frame(cbind(patient_id = left_bin$patient_id, 
                                        bin = left_bin$End-start,
                                        bin_id = bin_id))
            all <- rbind(all, left)
          }
          if(nrow(right_bin)>0){
            right <- as.data.frame(cbind(patient_id = right_bin$patient_id, 
                                         bin = end-right_bin$Start,
                                         bin_id = bin_id))
            all <- rbind(all, right)
          }
          
          return(all)
        })
        
        
        scores_del <- data.frame()
        for(i in 1:length(altered_bins)){
          
          if(any(is.na(altered_bins[[i]]))){
            score <- 0
          }else{
            score <- sum(as.numeric(altered_bins[[i]][,2]))/
              (length(common_patients)*(bin_gene_new[i,]$bin_end-bin_gene_new[i,]$bin_start))
          }
          
          scores_del <- rbind(scores_del, 
                              cbind(bin = bin_gene_new[i,]$bin,
                               del_score = score))
        }
        
        scores <- full_join(scores_del, scores_ampl, by = 'bin')
        
        scores <- left_join(bin_gene_new, scores, by = 'bin')
        
        write.table(scores, file = paste0("/home/ieo5099/mountHD/ml_models/data/chr_split/",
                                          tumor_type,"_chr",chr,"_",fixed_bin_length/1000000,"Mbp.txt"),
                  quote = F, row.names = F, col.names = T)
        
      }
      
      # end of the 1st loop: write chromosome/arm amplification frequencies ----
      cat("\n\n> >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Analysis on", tumor_type, " ended.\n\n")
    })
    # }
  }
}
