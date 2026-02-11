# rm(list=ls())
gc(full=T)

wd <- 'path/to/GitHub/copy-number-annotation/'
setwd(wd)

ml_merged_path <- "feature_matrix_data/ml_merged.rds"
chr_info_path <- 'feature_matrix_data/cytoBand.txt'
chr_coverage_path <- 'feature_matrix_data/averageChrCoverage.txt'
ml_merged_tel_centr_path <- "feature_matrix_data/ml_merged_tel.rds"
ampl_freq_path <- 'plots/00_ampl_freq.pdf'
del_freq_path <- 'plots/00_del_freq.pdf'
all_freq_path <- 'plots/00_all_freq.pdf'

ml_merged <- readRDS(ml_merged_path)

# add genomic locations ----
chr_info <- read.table(chr_info_path, header = T)
chr_info_seq <- read.table(chr_coverage_path, header = T)
# chr_info_seq <- chr_info_seq[!duplicated(chr_info_seq$chr),]

bins <- c(0.1,0.25,0.5,1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)

ml_merged_tel_centr <- list()
for(bin in bins){
  print(bin)
  df <- ml_merged[[paste0(bin,'Mbp')]]
  # colnames(df)[2:3] <- c('start', 'end')
  # df$chr <- as.integer(do.call(rbind, str_split(df$bin, pattern = '_'))[,1])
  
  final <- data.frame()
  for(chr in 1:22){
    df_chr <- df[df$chr == chr,]
    
    distance.to.telomere1 <- df_chr$end_bin - chr_info_seq[chr_info_seq$chr == chr,]$start
    distance.to.telomere2 <- -(df_chr$start_bin - chr_info_seq[chr_info_seq$chr == chr,]$end)
    distance.to.telomere <- apply(cbind(distance.to.telomere1, distance.to.telomere2), 1, min)
    
    distance.to.centromere1 <- chr_info[chr_info$arm == paste0(chr,'p'),]$end - df_chr$end_bin
    distance.to.centromere2 <- df_chr$end_bin - chr_info[chr_info$arm == paste0(chr,'q'),]$start
    
    distance.to.centromere <- apply(cbind(distance.to.centromere1, distance.to.centromere2), 1, max)
    
    tmp <- as.data.frame(cbind(bin = df_chr$bin, 
                 distance.to.centromere, 
                 distance.to.telomere, 
                 start_bin = df_chr$start_bin, 
                 end_bin = df_chr$end_bin, 
                 chr = chr))
    colnames(tmp)[1] <- 'bin'
    tmp <- tmp[!duplicated(tmp),]
    
    final <- rbind(final, tmp)
  }
  
  x <- full_join(df, final, relationship = "many-to-many")
  x <- x[!duplicated(x),]

  x$distance.to.centromere <- ifelse(x$distance.to.centromere < 0, 0, x$distance.to.centromere)
  x$distance.to.telomere <- ifelse(x$distance.to.telomere < 0, 0, x$distance.to.telomere)
  
  ml_merged_tel_centr[[paste0(bin,'Mbp')]] <- x
}

saveRDS(ml_merged_tel_centr, ml_merged_tel_centr_path)


if(F){
  # control the tables making plots ----
  bins <- c(1
            #,2,3,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50
  )
  
  bin <- 1
  for(bin in bins){
    print(bin)
    
    df_ml <- ml_merged_tel_centr[[paste0(bin,'Mbp')]]
    
    df_ml_tts <- df_ml[df_ml$type == 'BRCA',]
    # hist(df_ml_tts$cna_freq_ampl, breaks = 100)
    # df_ml_tts$binz <- do.call(rbind, str_split(df_ml_tts$bin, pattern = '_'))[,2]
    
    plots_ampl <- list()
    plots_del <- list()
    plots_all <- list()
    
    for(chr in 1:22){
      df_ml_tts_chr <- df_ml_tts[df_ml_tts$chr == chr,]
      
      plots_ampl[[chr]] <- ggplot(df_ml_tts_chr) + 
        geom_col(aes(y = as.numeric(ampl_score), x = as.numeric(start_bin))) +
        geom_line(aes(y = distance.to.telomere/(max(distance.to.telomere)/max(as.numeric(ampl_score))), x = as.numeric(start_bin)), group = 1, col = 'red') +
        geom_line(aes(y = distance.to.centromere/(max(distance.to.centromere)/max(as.numeric(ampl_score))), x = as.numeric(start_bin)), group = 1, col = 'blue') +
        ggtitle(paste0('Chromosome ', chr)) +
        xlab('Genomc position') +
        ylab('Amplification Frequency') +
        theme_classic() 
      plots_del[[chr]] <- ggplot(df_ml_tts_chr) + 
        geom_col(aes(y = as.numeric(del_score), x = as.numeric(start_bin))) +
        geom_line(aes(y = distance.to.telomere/(max(distance.to.telomere)/max(as.numeric(del_score))), x = as.numeric(start_bin)), group = 1, col = 'red') +
        geom_line(aes(y = distance.to.centromere/(max(distance.to.centromere)/max(as.numeric(del_score))), x = as.numeric(start_bin)), group = 1, col = 'blue') +
        ggtitle(paste0('Chromosome ', chr)) +
        xlab('Genomc position') +
        ylab('Deletion Frequency') +
        theme_classic()
    }
  }
  
  pdf(file = ampl_freq_path, width = 20, height = 20)
  ggpubr::ggarrange(plotlist = plots_ampl)
  dev.off()
  
  pdf(file = del_freq_path, width = 20, height = 20)
  ggarrange(plotlist = plots_del)
  dev.off()
  
  pdf(file = all_freq_path, width = 20, height = 20)
  ggarrange(plotlist = plots_all)
  dev.off()
}
