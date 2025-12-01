# rm(list=ls())
gc(full=T)

ml_merged <- readRDS("/mnt/fabiogokce/fabiohd/ml_models/data/ml_merged_tel.rds")

ml_merged_tel_centr <- list()
for(bin in names(ml_merged)){
  print(bin)
  df <- ml_merged[[bin]]
  # colnames(df)[2:3] <- c('start', 'end')
  # df$chr <- as.integer(do.call(rbind, str_split(df$bin, pattern = '_'))[,1])
  
  # PANCACER OG/TSG/Ess
  distances <- data.frame()
  for(chr in 1:22){
    df_chr <- df[df$chr == chr,]
    df_chr[is.na(df_chr)] <- 0
    df_chr <- df_chr[df_chr$type == 'BRCA',]
    
    dist <- data.frame(row.names = df_chr$bin)
    for(i in df_chr[df_chr$OGscore_pancancer >= 1,]$bin){
      dist <- cbind(dist, abs(df_chr$bin - i))
    }
    OG.distance_pancancer <- abs(as.numeric((apply(dist, 1, min))))
    rm(dist)
    
    dist <- data.frame(row.names = df_chr$bin)
    for(i in df_chr[df_chr$TSGscore_pancancer >= 1,]$bin){
      dist <- cbind(dist, abs(df_chr$bin - i))
    }
    TSG.distance_pancancer <-  abs(as.numeric((apply(dist, 1, min))))
    rm(dist)
    
    dist <- data.frame(row.names = df_chr$bin)
    for(i in df_chr[df_chr$essential_pancancer >= 1,]$bin){
      dist <- cbind(dist, abs(df_chr$bin - i))
    }
    Ess.distance_pancancer <-  abs(as.numeric((apply(dist, 1, min))))
    rm(dist)
    
    distances <- rbind(distances, as.data.frame(cbind(bin = df_chr$bin, 
                                     chr = chr,
                                     OG.distance_pancancer, 
                                     TSG.distance_pancancer, 
                                     Ess.distance_pancancer)))
    rm(dist)
  }
  
  df <- full_join(df, distances, by = c('bin', 'chr'))
  rm(distances, Ess.distance_pancancer, OG.distance_pancancer, 
     TSG.distance_pancancer, df_chr, i, chr)
  
  distances_TS <- data.frame()
  for(tumor_type in levels(factor(df$type))){
    df_ts <- df[df$type == tumor_type,]
    df_ts[is.na(df_ts)] <- 0
    
    for(chr in 1:22){
      df_ts_chr <- df_ts[df_ts$chr == chr,]
      colnames(df_ts_chr)[21:25]
      df_ts_chr[,21:25] <- apply(df_ts_chr[,21:25], 2, as.numeric)
      
      if(any(df_ts_chr$OGscore_TS >= 1)){
        dist <- data.frame(row.names = df_ts_chr$bin)
        for(i in df_ts_chr[df_ts_chr$OGscore_TS >= 1,]$bin){
          dist <- cbind(dist, abs(df_ts_chr$bin - i))
        }
        OG.distance_TS <- abs(as.numeric((apply(dist, 1, min))))
        rm(dist)
      }else{
        OG.distance_TS <- rep(NA, length(df_ts_chr$bin))
      }
      
      if(any(df_ts_chr$TSGscore_TS >= 1)){
        dist <- data.frame(row.names = df_ts_chr$bin)
        for(i in df_ts_chr[df_ts_chr$TSGscore_TS >= 1,]$bin){
          dist <- cbind(dist, abs(df_ts_chr$bin - i))
        }
        TSG.distance_TS <- abs(as.numeric((apply(dist, 1, min))))
        rm(dist)
      }else{
        TSG.distance_TS <- rep(NA, length(df_ts_chr$bin))
      }
      
      if(any(df_ts_chr$essential_TS >= 1)){
        dist <- data.frame(row.names = df_ts_chr$bin)
        for(i in df_ts_chr[df_ts_chr$essential_TS >= 1,]$bin){
          dist <- cbind(dist, abs(df_ts_chr$bin - i))
        }
        Ess.distance_TS <- abs(as.numeric((apply(dist, 1, min))))
        rm(dist)
      }else{
        Ess.distance_TS <- rep(NA, length(df_ts_chr$bin))
      }
      
      distances_TS <- rbind(distances_TS, cbind(chr = chr,
                                                bin = df_ts_chr$bin,
                                                type = tumor_type,
                                                OG.distance_TS, 
                                                TSG.distance_TS,
                                                Ess.distance_TS))
      
    }
  }
  
  distances_TS[is.na(distances_TS$OG.distance_TS),'OG.distance_TS'] <-
    max(as.numeric(distances_TS$OG.distance_TS), na.rm = T)
  distances_TS[is.na(distances_TS$TSG.distance_TS),'TSG.distance_TS'] <-
    max(as.numeric(distances_TS$TSG.distance_TS), na.rm = T)
  distances_TS[is.na(distances_TS$Ess.distance_TS),'Ess.distance_TS'] <-
    max(as.numeric(distances_TS$Ess.distance_TS), na.rm = T)
  
  distances_TS[,c(1:2,4:6)] <- apply(distances_TS[,c(1:2,4:6)], 2, as.numeric)
  
  df <- full_join(df, distances_TS, by = c('bin', 'chr', 'type'))
  rm(distances_TS, Ess.distance_TS, OG.distance_TS, 
     TSG.distance_TS, tumor_type, chr, i, df_ts_chr,df_ts)
  
  ml_merged_tel_centr[[bin]] <- df
}

saveRDS(ml_merged_tel_centr, "/mnt/fabiogokce/fabiohd/ml_models/data/ml_merged_tel_cg.rds")

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
        geom_line(aes(y = distance.to.centromere/(max(distance.to.centromere)/max(as.numeric(del_score))), x = as.numeric(start_bin)), group = 1, col = 'grey') +
        geom_line(aes(y = (OG.distance_pancancer), x = as.numeric(start_bin)), group = 1, col = 'red') +
        geom_line(aes(y = (TSG.distance_pancancer), x = as.numeric(start_bin)), group = 1, col = 'blue') +
        geom_line(aes(y = (n_genes), x = as.numeric(start_bin)), group = 1, col = 'black') +
        geom_line(aes(y = Ess.distance_pancancer, x = as.numeric(start_bin)), group = 1, col = 'green') +
        ggtitle(paste0('Chromosome ', chr)) +
        xlab('Genomic position') +
        ylab('Scores') +
        theme_classic() 
      # plots_del[[chr]] <- ggplot(df_ml_tts_chr) + 
      #   geom_col(aes(y = as.numeric(del_score), x = as.numeric(start_bin))) +
      #   geom_line(aes(y = distance.to.telomere/(max(distance.to.telomere)/max(as.numeric(del_score))), x = as.numeric(start_bin)), group = 1, col = 'red') +
      #   geom_line(aes(y = distance.to.centromere/(max(distance.to.centromere)/max(as.numeric(del_score))), x = as.numeric(start_bin)), group = 1, col = 'blue') +
      #   ggtitle(paste0('Chromosome ', chr)) +
      #   xlab('Genomc position') +
      #   ylab('Deletion Frequency') +
      #   theme_classic()
    }
  }
  
  pdf(file = '/home/ieo5099/Desktop_linux/aneuploidy_determinants/data/tables_for_ml/00_ampl_freq.pdf', width = 20, height = 20)
  ggpubr::ggarrange(plotlist = plots_ampl)
  dev.off()
  
  pdf(file = '/home/ieo5099/Desktop_linux/aneuploidy_determinants/data/tables_for_ml/00_del_freq.pdf', width = 20, height = 20)
  ggarrange(plotlist = plots_del)
  dev.off()
  
  pdf(file = '/home/ieo5099/Desktop_linux/aneuploidy_determinants/data/tables_for_ml/00_all_freq.pdf', width = 20, height = 20)
  ggarrange(plotlist = plots_all)
  dev.off()
}
