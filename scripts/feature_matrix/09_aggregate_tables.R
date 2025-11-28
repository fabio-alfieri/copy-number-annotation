# merge gokce and fabio's tables

# rm(list=ls())
gc(full=T)

load("/mnt/fabiogokce/gokce/Data/Backbone_tables_with_all_features_GS.RData")
ml_tables_gc <- backbone_tables_w_all_features_GS_v2
rm(backbone_tables_w_all_features_GS_v2)

# add mutation and aneuploidy frequencies ----
mut_cna_features <- readRDS("/mnt/fabiogokce/fabiohd/ml_models/data/mut_cna_features.rds")
for(i in names(mut_cna_features)){
  for(tt in names(mut_cna_features[[i]])){
    mut_cna_features[[i]][[tt]][,1:10] <- apply(mut_cna_features[[i]][[tt]][,1:10], 2, as.numeric)
    mut_cna_features[[i]][[tt]]$bin <- paste0(mut_cna_features[[i]][[tt]]$chr, '_', mut_cna_features[[i]][[tt]]$bin)
    ml_tables_gc[[i]][[tt]][,c('start_bin','end_bin','chr')] <- apply(ml_tables_gc[[i]][[tt]][,c('start_bin','end_bin','chr')], 2, as.numeric)
    print(any(is.na(ml_tables_gc[[i]][[tt]][,3:5])))
    
    ml_tables_gc[[i]][[tt]] <- full_join(ml_tables_gc[[i]][[tt]], 
                                         mut_cna_features[[i]][[tt]], 
                                         by = c('start_bin','end_bin','chr','bin'))
  }
}
rm(mut_cna_features)

# add pancacer driver gene and distances (drivers + tel/cent) ----
ml_merged_tel_cg <- readRDS("/mnt/fabiogokce/fabiohd/ml_models/data/ml_merged_tel_cg.rds")
for(i in names(ml_merged_tel_cg)){
  to_merge <- ml_merged_tel_cg[[i]]
  to_merge$n_patients <- NULL
  to_merge <- to_merge[-c(4:13,20:24,30:32)]
  to_merge <- unique(to_merge)
  to_merge$bin <- NULL
  
  to_merge_new <- data.frame()
  for(chr in 1:22){
    to_merge_chr <- to_merge[to_merge$chr == chr,]
    to_merge_chr$distance.to.centromere_scaled <- as.numeric(scale(to_merge_chr$distance.to.centromere))
    to_merge_chr$distance.to.telomere_scaled <- as.numeric(scale(to_merge_chr$distance.to.telomere))
    # to_merge_chr$OGscore_pancancer_scaled <- as.numeric(scale(as.numeric(to_merge_chr$OGscore_pancancer)))
    # to_merge_chr$TSGscore_pancancer_scaled <- as.numeric(scale(as.numeric(to_merge_chr$TSGscore_pancancer)))
    # to_merge_chr$Ess.distance_pancancer_ls(scaled <- as.numeric(scale(as.numeric(to_merge_chr$Ess.distance_pancancer)))
    to_merge_new <- rbind(to_merge_new, to_merge_chr)
  }
  
  to_merge_new[,3:9] <- apply(to_merge_new[,3:9], 2, as.numeric)
  
  for(tt in names(ml_tables_gc[[i]])){
    ml_tables_gc[[i]][[tt]] <- full_join(ml_tables_gc[[i]][[tt]], to_merge_new)
    ml_tables_gc[[i]][[tt]]$n_OG <- ml_tables_gc[[i]][[tt]]$OGscore_pancancer
    ml_tables_gc[[i]][[tt]]$OGscore_pancancer <- 
      ml_tables_gc[[i]][[tt]]$n_OG/ml_tables_gc[[i]][[tt]]$n_genes
    ml_tables_gc[[i]][[tt]]$n_TSG <- ml_tables_gc[[i]][[tt]]$TSGscore_pancancer
    ml_tables_gc[[i]][[tt]]$TSGscore_pancancer <- 
      ml_tables_gc[[i]][[tt]]$n_TSG/ml_tables_gc[[i]][[tt]]$n_genes
    ml_tables_gc[[i]][[tt]]$n_ESS <- ml_tables_gc[[i]][[tt]]$essential_pancancer
    ml_tables_gc[[i]][[tt]]$ESSscore_pancancer <- 
      ml_tables_gc[[i]][[tt]]$n_ESS/ml_tables_gc[[i]][[tt]]$n_genes
    ml_tables_gc[[i]][[tt]]$n_HAPLO <- ml_tables_gc[[i]][[tt]]$haploinsufficient_genes
    ml_tables_gc[[i]][[tt]]$haploinsufficient_genes <- NULL
    ml_tables_gc[[i]][[tt]]$HAPLOscore_pancancer <- 
      ml_tables_gc[[i]][[tt]]$n_HAPLO/ml_tables_gc[[i]][[tt]]$n_genes
  }
}
rm(to_merge_chr,to_merge,to_merge_new,chr,i,tt,ml_merged_tel_cg)

# add pancancer driver gene and distances (drivers + tel/cent) ----
ml_merged_tel_cg <- readRDS("/mnt/fabiogokce/fabiohd/ml_models/data/ml_merged_tel_cg.rds")
for(i in names(ml_merged_tel_cg)){
  for(tt in levels(factor(ml_merged_tel_cg[[i]]$type))){
    df <- ml_merged_tel_cg[[i]][ml_merged_tel_cg[[i]]$type == tt,
                                c('start_bin', 'end_bin', 'chr',
                                  'OGscore_TS', 'OGscore.potency_TS',
                                  'TSGscore_TS', 'TSGscore.potency_TS',
                                  'essential_TS','OG.distance_TS','TSG.distance_TS',
                                  'Ess.distance_TS')]
    df <- as.data.frame(apply(df, 2, as.numeric))
    df$n_OG_TS <- df$OGscore_TS
    df$OGscore_TS <- df$n_OG_TS
    df$n_TSG_TS <- df$TSGscore_TS
    df$TSGscore_TS <- df$df_n_TSG_TS
    df$n_ESS_TS <- df$essential_TS
    df$ESSscore_TS <- df$essential_TS
    df$essential_TS <- NULL
  
    ml_tables_gc[[i]][[tt]] <- full_join(ml_tables_gc[[i]][[tt]], df)
    
    ml_tables_gc[[i]][[tt]]$OGscore_TS <- 
      ml_tables_gc[[i]][[tt]]$OGscore_TS/ml_tables_gc[[i]][[tt]]$n_genes
    ml_tables_gc[[i]][[tt]]$TSGscore_TS <- 
      ml_tables_gc[[i]][[tt]]$n_TSG_TS/ml_tables_gc[[i]][[tt]]$n_genes
    ml_tables_gc[[i]][[tt]]$ESSscore_TS <- 
      ml_tables_gc[[i]][[tt]]$n_ESS_TS/ml_tables_gc[[i]][[tt]]$n_genes
    
  }
}

saveRDS(ml_tables_gc, file = '/mnt/fabiogokce/fabiohd/ml_models/ml_tables.rds')
# rm(list=ls())

# add chromosome and arm level info ----
ml_merged <- readRDS(file = '/mnt/fabiogokce/fabiohd/ml_models/ml_tables.rds')

arm <- read.table(file = '/mnt/fabiogokce/fabiohd/ml_models/data/misc/arm_Davoli2013.tsv', header = T)
chr <- read.table(file = '/mnt/fabiogokce/fabiohd/ml_models/data/misc/chr_Davoli2013.tsv', header = T)
arm$bin <- arm$chrarm
arm$chr <- NULL
arm$arm <- NULL
arm$chrarm <- NULL

for(i in c('Arm', 'Chromosome')){
  for(tt in names(ml_merged[[i]])){
    if(i == 'Chromosome'){
      ml_merged[[i]][[tt]] <- full_join(ml_merged[[i]][[tt]], chr, by = 'chr')
    }else{
      ml_merged[[i]][[tt]] <- full_join(ml_merged[[i]][[tt]], arm, by = 'bin')
    }
  }
}
rm(arm, chr, i, tt)

for(tt in names(ml_merged[['Arm']])){
  logical <- any(tt == 'SKCM' | tt == 'BLCA' | tt == 'PCPG' | tt == 'PRAD' |
        tt == 'KIRC' | tt == 'MESO' | tt == 'TGCT' | tt == 'KIRP' | tt == 'OV' |
        tt == 'SARC' | tt == 'LIHC' | tt == 'ESCA' | tt == 'STAD' | tt == 'UCS')
  if(logical){
    arm <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/mutation_compensation/results/tables/02_produceStatistics/',
                                    tt,'_chrArm.tsv'), header = T)
    arm$ampl <- NULL
    colnames(arm) <- c('bin', 'mutations_norm')
    ml_merged[['Arm']][[tt]] <- full_join(ml_merged[['Arm']][[tt]], arm, by = 'bin')
    
    arm$bin <- parse_number(arm$bin)
    arm <- aggregate(bin ~ mutations_norm, arm, mean)
    arm$mutations_norm <- as.numeric(scale(arm$mutations_norm))
    ml_merged[['Chromosome']][[tt]] <- full_join(ml_merged[['Chromosome']][[tt]], chr, by = 'chr')
  }else{
    chr <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/mutation_compensation/results/tables/02_produceStatistics/',
                                    tt,'_chr.tsv'), header = T)
    chr$ampl <- NULL
    colnames(chr)[2] <- 'mutations_norm'
    chr$chr <- as.numeric(chr$chr)
    chr$mutations_norm <- as.numeric(scale(chr$mutations_norm))
    ml_merged[['Chromosome']][[tt]] <- full_join(ml_merged[['Chromosome']][[tt]], chr, by = 'chr')
    
    arm <- read.table(file = paste0('/mnt/fabiogokce/fabiohd/mutation_compensation/results/tables/02_produceStatistics/',
                                    tt,'_chrArm.tsv'), header = T)
    arm$ampl <- NULL
    colnames(arm) <- c('bin', 'mutations_norm')
    arm$mutations_norm <- as.numeric(scale(arm$mutations_norm))
    colnames(arm)[1] <- 'bin'
    ml_merged[['Arm']][[tt]] <- full_join(ml_merged[['Arm']][[tt]], arm, by = 'bin')
  }
}

saveRDS(ml_merged, file = '/mnt/fabiogokce/fabiohd/ml_models/ml_tables_final.rds')

