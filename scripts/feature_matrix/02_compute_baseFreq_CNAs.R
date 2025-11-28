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

focal <- T

cores <- length(tumor_types)

for(fixed_bin_length in c(#0.1,0.25,0.5,
                          1,2,3,4#,6,8,10,12,14,16,18,20,22,24,26,
                          #28,30,32,34,36,38,40,42,44,46,48
                          )*1000000){
    
  mclapply(tumor_types, mc.cores = cores, function(tumor_type){
    
    for_immune_score <- list()
    for_immune_score[[tumor_type]] <- list()
    for_immune_score[[tumor_type]]$ampl <- list()
    for_immune_score[[tumor_type]]$del <- list()
    
      ## 1st Loop: load sCNAs and SNVs files ----
      snv <-
        read.csv(file = paste0("data/FireBrowse_SNVs/", tumor_type, "_mutations.csv"))
      snv <-
        snv[!duplicated(snv[, c(6, 7, 8, 9, 17)]),] # remove duplicated mutations in the same patient
      # cat("  done!")
      
      # Load CNAs file
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
      
      scna$length <- scna$length+1
      
      if(focal == T){
        scna <- scna[scna$classified == 'focal',]
        cytoband <- read.table(file = '/mnt/fabiogokce/fabiohd/mutation_compensation/data/misc/chr_info_h19.txt',
                               header = T)
      }
      
      rm(snv)
      gc(full = T)
      
      chr_scores <- data.frame()
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
        
        if(focal == T){
          filter_chr <- 0.35*(cytoband[cytoband$Chromosome == paste0('chr',chr),]$Length/2)
          temp_cna <- temp_cna[temp_cna$length <= filter_chr,]
        }
        
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
        colnames(bin_gene)[3:4] <- c('start_bin', 'end_bin')
        bin_gene <- as.data.frame(apply(bin_gene, 2, as.numeric))
        
        ## 2nd loop: compute chromosome and arm amplification frequencies ----
        temp_cna_length_backup <- temp_cna_length
        temp_cna_length <-
          temp_cna_length %>% filter(Segment_Mean >= as.numeric(segment_cutoff)/100)
        
        # rm(altered_bins)
        altered_bins <- lapply(1:nrow(bin_gene), function(i){
          # print(i)
          
          start <- as.numeric(bin_gene[i,]$start_bin)
          end <- as.numeric(bin_gene[i,]$end_bin)
          
          in_bin <- temp_cna_length[!(temp_cna_length$Start < start & 
                                        temp_cna_length$End < start),]
          in_bin <- in_bin[!(in_bin$End > end &
                               in_bin$Start > end),]
          
          all_bin <- in_bin[in_bin$Start <= start & in_bin$End >= end,]
          partial_bin <- in_bin[!(in_bin$Start <= start & in_bin$End >= end),]
          
          if(nrow(all_bin)==0 & nrow(partial_bin)==0){return(NA)}
          
          small_bin <- partial_bin[partial_bin$length < end-start,]
          left_bin <- partial_bin[partial_bin$End > start & partial_bin$End < end,]
          right_bin <- partial_bin[partial_bin$Start > start & partial_bin$End > end,]
          
          if(nrow(all_bin)>0){
            all <- as.data.frame(cbind(patient_id = all_bin$patient_id, 
                                       length = end-start,
                                       sg = 2*2^all_bin$Segment_Mean))
          }else{
            all <- data.frame()
          }
          if(nrow(small_bin)>0){
            small <- as.data.frame(cbind(patient_id = small_bin$patient_id,
                                         length = small_bin$End-small_bin$Start,
                                         sg = 2*2^small_bin$Segment_Mean))
            all <- rbind(all, small)
          }
          if(nrow(left_bin)>0){
            left <- as.data.frame(cbind(patient_id = left_bin$patient_id, 
                                        length = left_bin$End-start,
                                        sg = 2*2^left_bin$Segment_Mean))
            all <- rbind(all, left)
          }
          if(nrow(right_bin)>0){
            right <- as.data.frame(cbind(patient_id = right_bin$patient_id, 
                                         length = end-right_bin$Start,
                                         sg = 2*2^right_bin$Segment_Mean))
            all <- rbind(all, right)
          }
          
          return(all)
        })
        altered_bins <- lapply(altered_bins, function(df){
          df <- as.data.frame(df)
          df <- df[as.numeric(df$length) >= as.numeric(fixed_bin_length/2),]
          return(df)
        })
        
        for_immune_score[[tumor_type]][['ampl']][[paste0('chr',chr)]] <- altered_bins
        
        scores_ampl <- data.frame()
        for(i in 1:length(altered_bins)){

          if(is.logical(altered_bins
                       [[i]])){
            score <- 0
            score_sg <- 2
          }else{
            score <- sum(as.numeric(altered_bins[[i]][,2]))/
              (length(common_patients)*(bin_gene[i,]$end_bin-bin_gene[i,]$start_bin))

            patients <- names(table(altered_bins[[i]]$patient_id)[table(altered_bins[[i]]$patient_id) > 1])
            for(pt in patients){
              length <-sum(as.numeric(altered_bins[[i]][altered_bins[[i]]$patient_id == pt,]$length))
              sg <- weighted.mean(as.numeric(altered_bins[[i]][altered_bins[[i]]$patient_id == pt,]$sg),
                                  as.numeric(altered_bins[[i]][altered_bins[[i]]$patient_id == pt,]$length),
                                  na.rm = T)
              altered_bins[[i]] <- altered_bins[[i]][altered_bins[[i]]$patient_id != pt,]
              altered_bins[[i]] <- rbind(altered_bins[[i]],
                                         cbind(
                                           patient_id = pt,
                                           length = length,
                                           sg = sg
                                         ))
            }

            bin_length <- (bin_gene[bin_gene$bin == i,]$end_bin -
                             bin_gene[bin_gene$bin == i,]$start_bin)
            if(any(as.numeric(altered_bins[[i]]$length) < bin_length)){
              altered_bin_ok <- altered_bins[[i]][!as.numeric(altered_bins[[i]]$length) < bin_length,]
              altered_bin_no <- altered_bins[[i]][as.numeric(altered_bins[[i]]$length) < bin_length,]

              altered_bin_no_corrected <- data.frame()
              for(n in 1:nrow(altered_bin_no)){
                x <- altered_bin_no[n,]
                x <- rbind(x,
                           cbind(
                             patient_id = x$patient_id,
                             length = (bin_length-as.numeric(x$length)),
                             sg = 2
                           ))
                altered_bin_no_corrected <- rbind(altered_bin_no_corrected,
                                                  cbind(
                                                    patient_id = x$patient_id[1],
                                                    length = bin_length,
                                                    sg = weighted.mean(as.numeric(x$sg),
                                                                       as.numeric(x$length),
                                                                       na.rm = T)
                                                  ))
              }

              altered_bins[[i]] <- rbind(altered_bin_ok,
                                         altered_bin_no_corrected)
            }


            score_sg <- mean(c(as.numeric(altered_bins[[i]]$sg),rep(2,length(common_patients)-length(levels(factor(altered_bins[[i]]$patient_id))))))
          }

          scores_ampl <- rbind(scores_ampl,
                               cbind(bin = paste0(chr,'_',i),
                                     start_bin = bin_gene[i,]$start_bin,
                                     end_bin = bin_gene[i,]$end_bin,
                                     chr = chr,
                                     ampl_score = score,
                                     ampl_score_sg = score_sg))
        }

        bin_gene$bin <- paste0(chr,'_',bin_gene$bin)
        bin_gene$gene_count <- NULL
        bin_gene$length_coding <- NULL
        bin_gene$length_perc <- NULL
        bin_gene$gene_id <- NULL
        
        
        ## 2nd loop: compute amplification and deletion frequencies ----
        chr_bins <- data.frame(bins = 1:n_bins)
        
        temp_cna_length <-
          temp_cna_length_backup %>% filter(Segment_Mean <= -as.numeric(segment_cutoff)/100)
        
        rm(altered_bins)
        altered_bins <- lapply(1:nrow(bin_gene), function(i){
          
          start <- bin_gene[i,]$start_bin
          end <- bin_gene[i,]$end_bin
          
          in_bin <- temp_cna_length[!(temp_cna_length$Start < start & 
                                        temp_cna_length$End < start),]
          in_bin <- in_bin[!(in_bin$End > end &
                               in_bin$Start > end),]
          
          all_bin <- in_bin[in_bin$Start <= start & in_bin$End >= end,]
          partial_bin <- in_bin[!(in_bin$Start <= start & in_bin$End >= end),]
          
          if(nrow(all_bin)==0 & nrow(partial_bin)==0){return(NA)}
          
          small_bin <- partial_bin[partial_bin$length < end-start,]
          left_bin <- partial_bin[partial_bin$End > start & partial_bin$End < end,]
          right_bin <- partial_bin[partial_bin$Start > start & partial_bin$End > end,]
          
          if(nrow(all_bin)>0){
            all <- as.data.frame(cbind(patient_id = all_bin$patient_id, 
                                       length = end-start,
                                       sg = 2*2^all_bin$Segment_Mean))
          }else{
            all <- data.frame()
          }
          if(nrow(small_bin)>0){
            small <- as.data.frame(cbind(patient_id = small_bin$patient_id,
                                         length = small_bin$End-small_bin$Start,
                                         sg = 2*2^small_bin$Segment_Mean))
            all <- rbind(all, small)
          }
          if(nrow(left_bin)>0){
            left <- as.data.frame(cbind(patient_id = left_bin$patient_id, 
                                        length = left_bin$End-start,
                                        sg = 2*2^left_bin$Segment_Mean))
            all <- rbind(all, left)
          }
          if(nrow(right_bin)>0){
            right <- as.data.frame(cbind(patient_id = right_bin$patient_id, 
                                         length = end-right_bin$Start,
                                         sg = 2*2^right_bin$Segment_Mean))
            all <- rbind(all, right)
          }
          # colnames(all)[ncol(all)] <- paste0('bin_',i)
          
          return(all)
        })
        
        altered_bins <- lapply(altered_bins, function(df){
          df <- as.data.frame(df)
          df <- df[as.numeric(df$length) >= as.numeric(fixed_bin_length/2),]
          return(df)
        })
        
        for_immune_score[[tumor_type]][['del']][[paste0('chr',chr)]] <- altered_bins
        
        
        scores_del <- data.frame()
        for(i in 1:length(altered_bins)){

          if(is.logical(altered_bins
                       [[i]])){
            score <- 0
            score_sg <- 2
          }else{
            score <- sum(as.numeric(altered_bins[[i]][,2]))/
              (length(common_patients)*(bin_gene[i,]$end_bin-bin_gene[i,]$start_bin))

            patients <- names(table(altered_bins[[i]]$patient_id)[table(altered_bins[[i]]$patient_id) > 1])
            for(pt in patients){
              # altered_bins[[i]][altered_bins[[i]]$patient_id == pt,]
              length <-sum(as.numeric(altered_bins[[i]][altered_bins[[i]]$patient_id == pt,]$length))
              sg <- weighted.mean(as.numeric(altered_bins[[i]][altered_bins[[i]]$patient_id == pt,]$sg),
                                  as.numeric(altered_bins[[i]][altered_bins[[i]]$patient_id == pt,]$length),
                                  na.rm = T)
              altered_bins[[i]] <- altered_bins[[i]][altered_bins[[i]]$patient_id != pt,]
              altered_bins[[i]] <- rbind(altered_bins[[i]],
                                         cbind(
                                           patient_id = pt,
                                           length = length,
                                           sg = sg
                                         ))
            }

            bin_length <- (bin_gene[bin_gene$bin == i,]$end_bin -
                             bin_gene[bin_gene$bin == i,]$start_bin)
            if(any(as.numeric(altered_bins[[i]]$length) < bin_length)){
              altered_bin_ok <- altered_bins[[i]][!as.numeric(altered_bins[[i]]$length) < bin_length,]
              altered_bin_no <- altered_bins[[i]][as.numeric(altered_bins[[i]]$length) < bin_length,]

              altered_bin_no_corrected <- data.frame()
              for(n in 1:nrow(altered_bin_no)){
                x <- altered_bin_no[n,]
                x <- rbind(x,
                           cbind(
                             patient_id = x$patient_id,
                             length = (bin_length-as.numeric(x$length)),
                             sg = 2
                           ))
                altered_bin_no_corrected <- rbind(altered_bin_no_corrected,
                                                  cbind(
                                                    patient_id = x$patient_id[1],
                                                    length = bin_length,
                                                    sg = weighted.mean(as.numeric(x$sg),
                                                                       as.numeric(x$length),
                                                                       na.rm = T)
                                                  ))
              }

              altered_bins[[i]] <- rbind(altered_bin_ok,
                                         altered_bin_no_corrected)
            }


            score_sg <- mean(c(as.numeric(altered_bins[[i]]$sg),rep(2,length(common_patients)-length(levels(factor(altered_bins[[i]]$patient_id))))))
          }

          scores_del <- rbind(scores_del,
                              cbind(bin = paste0(chr,'_',i),
                                    start_bin = bin_gene[i,]$start_bin,
                                    end_bin = bin_gene[i,]$end_bin,
                                    chr = chr,
                                    del_score = score,
                                    del_score_sg = score_sg))
        }

        scores_del[,2:3] <- apply(scores_del[,2:3], 2, as.numeric)
        scores_ampl[,2:3] <- apply(scores_ampl[,2:3], 2, as.numeric)

        scores <- full_join(scores_del, scores_ampl, by = c('bin', 'start_bin', 'end_bin', 'chr'))

        chr_scores <- rbind(chr_scores, scores)
        
      }
      
      # saveRDS(for_immune_score, file = 
      #           paste0('/mnt/fabiogokce/fabiohd/US_intern/altered_bins/',tumor_type,'_altered_bins_',fixed_bin_length/1000000,
      #                  'Mbp.rds'))
      
      write.table(chr_scores, file = paste0("/mnt/fabiogokce/fabiohd/ml_models/data/cna_tables_focal/",
                                        tumor_type,'_',fixed_bin_length/1000000,"Mbp.txt"),
                  quote = F, row.names = F, col.names = T)
    })
}

cat(
  "\n\n OUTPUT of the script: \n \t (1) raw tables path: results/tables/01_binLevelAnalysis/ \n"
)

