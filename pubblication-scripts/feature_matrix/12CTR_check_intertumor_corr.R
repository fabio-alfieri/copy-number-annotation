# rm(list=ls())
gc(full=T)

# check how immune scores and amplification frequencies are correlated each other across different cell lines ----
load("/mnt/fabiogokce/fabiohd/ml_models/ml_tables_from_0_prepare_tables_240329.RData")  

i <- '1Mbp'
ml_table <- list()
ml_table.subset <- list()
cor.df <- data.frame()
for(tt in c('LUSC','LUAD','BRCA','PAAD','SKCM','COADREAD','OV')){
  immune <- readRDS(file = paste0('/mnt/fabiogokce/fabiohd/US_intern/results_df/',
                                  tt,'_',i,'_ImmuneScore_CIBERSORTx.rds'))
  immune.df <- do.call(rbind, immune)
  colnames(immune.df)[grep("wil.ampl",colnames(immune.df))]
  colnames(immune.df)[grep("wil.del",colnames(immune.df))]
  rm(immune)
  
  immune.df$bin <- paste0(immune.df$chr,'_', immune.df$bin)
  
  ml_table[[tt]] <- full_join(ML.tables$`1Mbp`$LUSC[,c('ampl_score','del_score','bin')], 
                              immune.df[,c(1,2,grep('kruskall', colnames(immune.df)))],
                              by = 'bin')
  rm(immune.df)
  colnames(ml_table[[tt]])[-c(3:4)] <- paste0(colnames(ml_table[[tt]])[-c(3:4)], '__', tt)
}

cor.df <- full_join(ml_table$LUSC, ml_table$LUAD)
cor.df <- full_join(cor.df, ml_table$BRCA)
cor.df <- full_join(cor.df, ml_table$PAAD)
cor.df <- full_join(cor.df, ml_table$SKCM)
cor.df <- full_join(cor.df, ml_table$COADREAD)
cor.df <- full_join(cor.df, ml_table$OV)

cor.table <- cor(cor.df[,-c(3,4)])
ComplexHeatmap::Heatmap(cor.table, row_order = 1:nrow(cor.table),
                        column_order = 1:ncol(cor.table))

