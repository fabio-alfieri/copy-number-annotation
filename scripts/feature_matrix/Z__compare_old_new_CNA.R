# rm(list=ls())
gc(full=T)

library(ggpubr)

# tumor_type <- 'LUSC'
# chr <- 1
fixed_bin_length <- 1000000

p <- list()
for(tumor_type in c('LUSC','LUAD','COADREAD','BRCA','PAAD')){
  final <- data.frame()
  rows <- 1
  for(chr in 1:22){
    new_freq <- read.table(paste0("/home/ieo5099/mountHD/buffering_events/data/new_CNA_computation/",
                                  tumor_type,"_chr",chr,"_",fixed_bin_length,"bins.tsv"), header = T)
    old_freq <- read.table(paste0("/home/ieo5099/mountHD/mutation_compensation/results/tables/01_binLevelAnalysis/all_mutations_vs_CN/",
                                  tumor_type,"_chr",chr,"_",as.integer(fixed_bin_length),"BIN_table.txt"), header = T)
    old_freq$bin <- paste0(chr,'_',old_freq$bin)
    merged <- right_join(old_freq, new_freq, by = 'bin')
    merged$gene_id <- NULL
    
    merged$counts <- rows:(rows+nrow(merged)-1)
    rows <- rows+nrow(merged)
    
    final <- rbind(final, merged)
  }
  
  final$gene_count <- NULL  
  final$length_coding <- NULL
  final$mutations_norm <- NULL
  final$length_perc <- NULL
  final$cna_freq_total <- NULL
  final$cna_freq_ampl <- NULL
  final$cna_freq_del <- NULL
  final$n_patients <- NULL
  final$mutations_raw <- NULL
  
  write.table(final, file = paste0('/home/ieo5099/mountHD/ml_models/',tumor_type,'_newCNAcalc.tsv'), 
              sep = '\t', quote = F, col.names = T, row.names = F)
  
  # p[[tumor_type]][[1]] <- ggplot(final, aes(x = counts)) +
  #   geom_line(aes(y = ampl_score)) +
  #   geom_line(aes(y = cna_freq_ampl, color = 'red')) +
  #   ggtitle(paste(tumor_type, fixed_bin_length, ' - Amplifications (red is the old freq)')) +
  #   theme_classic() +
  #   theme(legend.position = "None")
  # 
  # p[[tumor_type]][[2]] <- ggplot(final, aes(x = cna_freq_ampl, y = ampl_score)) +
  #   geom_point() +
  #   theme_classic() +
  #   ggtitle(paste(tumor_type, fixed_bin_length, ' - Amplifications'))
  # 
  # p[[tumor_type]][[3]] <- ggplot(final, aes(x = counts)) +
  #   geom_line(aes(y = del_score)) +
  #   geom_line(aes(y = cna_freq_del, color = 'red')) +
  #   ggtitle(paste(tumor_type, fixed_bin_length, ' - Deletions (red is the old freq)')) +
  #   theme_classic() +
  #   theme(legend.position = "None")
  # 
  # p[[tumor_type]][[4]] <- ggplot(final, aes(x = cna_freq_del, y = del_score)) +
  #   geom_point() +
  #   theme_classic() +
  #   ggtitle(paste(tumor_type, fixed_bin_length, ' - Deletions'))
}

pdf(file = '/home/ieo5099/mountHD/ml_models/compare_old_new.pdf', width = 16)
for(tumor_type in c('LUSC','LUAD','COADREAD','BRCA','PAAD')){
  print(ggarrange(plotlist = p[[tumor_type]], 
                  nrow = 2, ncol = 2, 
                  widths = c(3,1)))
}
dev.off()
