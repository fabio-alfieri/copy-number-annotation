# rm(list=ls())
gc(full=T)

load("~/mountHD/ml_models/data/FinalStates_Averaged_Scaled&amp;counts_9Feb24.RData")
ml_tables <- readRDS("~/mountHD/ml_models/data/ml_tables.rds")

df_merged <- list()
for(bin_size in c(1,2,3,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)){
  for(tumor_type in c('LUSC','LUAD','COADREAD','PAAD','BRCA','GBMLGG')){
    df <- chromStates.tissues[[tumor_type]][[paste0(bin_size,'Mbp')]]$Counts
    chr_bin <- as.data.frame(do.call(rbind, str_split(rownames(df), '_')))
    df <- cbind(chr_bin, df)
    rownames(df) <- NULL
    colnames(df)[1:3] <- c('Chromosome','Bin_Start','Bin_End')
    df$Chromosome <- parse_number(df$Chromosome)
    df$type <- tumor_type
    df[,1:3] <- apply(df[,1:3], 2, as.numeric)
    
    feature_table <- ml_tables[[bin_size]] %>% filter(type == tumor_type)
    feature_table <- left_join(feature_table, df, by = c("Chromosome","Bin_Start","Bin_End","type"))
    
    colnames(feature_table)
    feature_table_filt <- feature_table[,c(1,31,32,33,41,36,37,
                                           34,35,40,42,43,44,
                                           9,2,23,16)]
    feature_table_filt <- cbind(feature_table_filt, feature_table[,grep(tumor_type, colnames(feature_table))])
    
    df_merged[[tumor_type]][[paste0(bin_size,'Mbp')]] <- feature_table_filt
  }
}
