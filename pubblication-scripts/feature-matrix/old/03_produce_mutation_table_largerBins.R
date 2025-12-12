# rm(list=ls())
gc(full=T)

add.col <- function(df, segment_length) {
  # df <- df[df$gene_count != 0, ]
  df <- cbind(df,
              resize = rep(
                1:round(dim(df)[1] / segment_length + 0.5),
                each = round(dim(df)[1] / round(
                  dim(df)[1] / segment_length + 0.5
                ) + 0.5)
              )[1:nrow(df)])
  return(df)
}
merge.bins <- function(df) {
  df2 <- data.frame()
  for (i in levels(factor(df$resize))) {
    # df <- df[!(df$gene_count <= 1 & df$length_perc <= 0.05), ]
    start <- 
      df[df$resize == i, ]$start_bin[1]
    end <-
      df[df$resize == i, ]$end_bin[length(df[df$resize == i, ]$start_bin)]
    mutations_raw <- 
      sum(as.numeric(df[df$resize == i, ]$mutations_raw))
    mutations_norm <-
      mean(as.numeric(df[df$resize == i, ]$mutations_norm))
    
    df2 <- rbind.data.frame(
      df2,
      c(
        i,
        start,
        end,
        mutations_raw,
        mutations_norm
      ),
      stringsAsFactors = FALSE
    )
  }
  colnames(df2) <- c(
    "bin",
    "start_bin",
    "end_bin",
    "mutations_raw",
    "mutations_norm"
  )
  return(df2)
}

library(parallel)

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

all_features <-list()

for(segment_length in c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30
                        ,32,34,36,38,40,42,44,46,48
                        )){
  print(segment_length)
  all_features[[paste0(segment_length,'Mbp')]] <- list()
  
  all_features[[paste0(segment_length,'Mbp')]] <- mclapply(tumor_types, mc.cores = 23, function(tumor_type){
    chr_backbone <- list()
    
    bin_gene_mut <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/ml_models/data/mutation_tables/',
                                             tumor_type,'_1Mbp.txt'), sep = ' ')
    colnames(bin_gene_mut) <- bin_gene_mut[1,]
    bin_gene_mut <- bin_gene_mut[-1,]
    
    merge_chr <- data.frame()
    for(chr in 1:22){
      bin_gene_chr <- bin_gene_mut[bin_gene_mut$chr == chr,1:6]
      bin_gene_sg <- add.col(bin_gene_chr, segment_length)
      merge_chr <- rbind(merge_chr, cbind(chr = as.character(chr),
                                          merge.bins(bin_gene_sg)))
    }
    
    bin_gene_cna <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/ml_models/data/cna_tables/',
                                             tumor_type,'_',segment_length,'Mbp.txt'), sep = ' ')
    colnames(bin_gene_cna) <- bin_gene_cna[1,]
    bin_gene_cna <- bin_gene_cna[-1,]
    
    return(cbind(full_join(bin_gene_cna[,-1], merge_chr),
                 type = tumor_type)
           )
  })
  
  names(all_features[[paste0(segment_length,'Mbp')]]) <- tumor_types
  # for(tumor_type in tumor_types){
    # all_features[[segment_length]][[tumor_type]] <- full_join(bin_gene_cna[,-1], merge_chr) 
  # }
}


for(tumor_type in tumor_types){
  bin_gene_mut <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/ml_models/data/mutation_tables/',
                                           tumor_type,'_1Mbp.txt'), sep = ' ')
  colnames(bin_gene_mut) <- bin_gene_mut[1,]
  bin_gene_mut <- bin_gene_mut[-1,]
  
  bin_gene_cna <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/ml_models/data/cna_tables/',
                                           tumor_type,'_1Mbp.txt'), sep = ' ')
  colnames(bin_gene_cna) <- bin_gene_cna[1,]
  bin_gene_cna <- bin_gene_cna[-1,]
  
  all_features[['1Mbp']][[tumor_type]] <- cbind(full_join(bin_gene_cna[,-1], bin_gene_mut[,-7]),
                                                type = tumor_type)
}
names(all_features[['1Mbp']]) <- tumor_types

mut_cna_features <- all_features

saveRDS(mut_cna_features, file = '/mnt/fabiogokce/fabiohd/ml_models/data/mut_cna_features.rds')


# add 0.1, 0.25 and 0.5 Mbp ----
mut_cna_features <- readRDS("/mnt/fabiogokce/fabiohd/ml_models/data/misc/mut_cna_features.rds")

for(segment_length in c(0.1,0.25,0.5)){
  mut_cna_features[[paste0(segment_length,'Mbp')]] <- list()
  for(tumor_type in tumor_types){
    bin_gene_mut <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/ml_models/data/mutation_tables/',
                                             tumor_type,'_',segment_length,'Mbp.txt'), sep = ' ', header = T)
    bin_gene_cna <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/ml_models/data/cna_tables/',
                                             tumor_type,'_',segment_length,'Mbp.txt'), sep = ' ', header = T)
    bin_gene_cna <- bin_gene_cna[,-1]
    
    mut_cna_features[[paste0(segment_length,'Mbp')]][[tumor_type]] <- cbind(full_join(bin_gene_mut, bin_gene_cna), type = tumor_type)
  }
}

saveRDS(mut_cna_features, file = '/mnt/fabiogokce/fabiohd/ml_models/data/mut_cna_features.rds')
