# rm(list=ls())
gc(full=T)

load("/mnt/fabiogokce/fabiohd/ml_models/ml_tables_from_0_prepare_tables_240329.RData")

tumor_types <- names(ML.tables$`0.1Mbp`)

ml_tables <- list()
for(i in names(ML.tables)[1:28]){
  for(tt in names(ML.tables[[i]])){
    cna_focal <- read.table(paste0('/mnt/fabiogokce/fabiohd/ml_models/data/cna_tables_focal/',
                                   tt,'_',i,'.txt'), header = T)
    cna_focal <- cna_focal[,c(1,5,7)]
    colnames(cna_focal) <- c('bin', 'del_score_focal', 'ampl_score_focal')
    
    ml_tables[[i]][[tt]] <- full_join(ML.tables[[i]][[tt]], cna_focal, by = 'bin')
  }
}
  
ML.tables <- ml_tables
save(ML.tables, file = "/mnt/fabiogokce/fabiohd/ml_models/ml_tables_from_11_focalCNAscore.RData")
