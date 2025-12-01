rm(list=ls())
gc(full=T)

library(stringr)
library(parallel)
library(reshape2)
library(dplyr)
library(ggplot2)
library(Matrix)
library(caret)
library(caTools)
library(dplyr)
require(xgboost)
require(tidyverse)

setwd("/mnt/ghost/fabiogokce/gokce")

# Note: HIPPIE and Interactome INSIDER parts are alternative to each other

# Interactome INSIDER
#-----
output.files <- list.files(path = "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/MLTables_updated/InteractomeINSIDER",
                           full.names = TRUE, recursive = FALSE, pattern = "ampl_or_del_0.1Mbp_covThr_zero")

output.files <- output.files[2]
#output.files <- output.files[-3] # remove smoothed file
#output.files <- output.files[-3] # remove smoothed file
#output.files <- output.files[2] # select only clustered mid and small
#output.files <- output.files[1] # remove location

name <- "ampl_or_del_0.1Mbp_covThr_zero"
level <- "0.1Mbp"
bincov <- "covThr_zero"
folder <- "InteractomeINSIDER"
scale <- FALSE
cancer_type <- 'pancancer'
#-----

if(F){
  #HIPPIE
  #-----
  output.files <- list.files(path = "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/MLTables_updated/HIPPIE",
                             full.names = TRUE, recursive = FALSE, pattern = "ampl_or_del_0.1Mbp_covThr_zero")
  output.files <- output.files[grep("Length|Length_Location",output.files)]
  output.files <- output.files[-3] # remove smoothed file if exist
  
  name <- "ampl_or_del_0.1Mbp_covThr_zero"
  level <- "0.1Mbp"
  bincov <- "covThr_zero"
  folder <- "HIPPIE"
  #-----
}

chromosome <- read.table(file = 'Codes/Codes-CNAs/MethodII/fabio/Data/chromosome.tsv', header = T)
centromeres <- read.table(file = 'Codes/Codes-CNAs/MethodII/fabio/Data/centomere.tsv', header = T)

tmp <- full_join(chromosome, centromeres, by = 'chr') %>% dplyr::select(chr, Centromere_Length, Chromosome_Length, Centromere_Type)
tmp$chr <- as.character(tmp$chr)

for(file in output.files){
  load(file)
  clus.group <- gsub("_ampl_or_del_0.1Mbp_covThr_zero.RData","",basename(file))
  
  ############## GAB MADE A MAJOR EDIT HERE ##############
  
  if (clus.group == "Length") {
    
    clusters <- names(ML.Tables[["ampl_or_del"]])[1:2] # select only Mid-Length models
    
  } else if (clus.group == "Length_Location") {
    
    clusters <- names(ML.Tables[["ampl_or_del"]])[1:2] # select only clustered mid and small [3:9] #[1:2] # select only arm and chr
    
  } else if (clus.group == "NoCluster") {
    
    ML.Tables$ampl_or_del <- list(no_cluster = ML.Tables$ampl_or_del) # since for no cluster ML.Tables$ampl_or_del is the data.frame, I put it into a list.
    
    clusters <- "no_cluster"
    
  }
  
  
  for(classS in clusters){
    ml_table_backup <- ML.Tables[["ampl_or_del"]][[classS]]
    ml_table_backup$bin <- as.character(ml_table_backup$bin)
    print(classS)

    ml_table_backup$chr <- do.call(rbind, str_split(ml_table_backup$bin,'_'))[,1]
    ml_table_backup <- left_join(ml_table_backup, tmp, by = 'chr')
    
    if(scale == T){
      ml_table_backup$mean.GC.content <- ml_table_backup %>% select(mean.GC.content) %>% scale(center = 0) %>% as.vector()
      ml_table_backup$Centromere_Length <- ml_table_backup %>% select(Centromere_Length) %>% scale(center = 0) %>% as.vector()
      ml_table_backup$Chromosome_Length <- ml_table_backup %>% select(Chromosome_Length) %>% scale(center = 0) %>% as.vector()
      
      ml_table_backup <- ml_table_backup %>% group_by(chr) %>%
        mutate(Length_Counts.E1 = scale(Length_Counts.E1, center = 0),
               Length_Counts.E2 = scale(Length_Counts.E2, center = 0),
               Length_Counts.E3 = scale(Length_Counts.E3, center = 0),
               Length_Counts.E4 = scale(Length_Counts.E4, center = 0),
               Length_Counts.E5 = scale(Length_Counts.E5, center = 0),
               Length_Counts.E6 = scale(Length_Counts.E6, center = 0),
               Length_Counts.E7 = scale(Length_Counts.E7, center = 0),
               Length_Counts.E8 = scale(Length_Counts.E8, center = 0),
               Length_Counts.E9 = scale(Length_Counts.E9, center = 0),
               Length_Counts.E10 = scale(Length_Counts.E10, center = 0),
               Length_Counts.E11 = scale(Length_Counts.E11, center = 0),
               Length_Counts.E12 = scale(Length_Counts.E12, center = 0),
               Length_Counts.E13 = scale(Length_Counts.E13, center = 0),
               Length_Counts.E14 = scale(Length_Counts.E14, center = 0),
               Length_Counts.E15 = scale(Length_Counts.E15, center = 0),
               Length_Counts.E16 = scale(Length_Counts.E16, center = 0),
               Length_Counts.E17 = scale(Length_Counts.E17, center = 0),
               Length_Counts.E18 = scale(Length_Counts.E18, center = 0),
               Length_Counts.E19 = scale(Length_Counts.E19, center = 0),
               Length_Counts.E20 = scale(Length_Counts.E20, center = 0),
               Length_Counts.E21 = scale(Length_Counts.E21, center = 0),
               Length_Counts.E22 = scale(Length_Counts.E22, center = 0),
               Length_Counts.E23 = scale(Length_Counts.E23, center = 0),
               Length_Counts.E24 = scale(Length_Counts.E24, center = 0),
               Length_Counts.E25 = scale(Length_Counts.E25, center = 0),
               dist.to.closest.FGS = scale(dist.to.closest.FGS, center = 0),
               dist.to.closest.OG = scale(dist.to.closest.OG, center = 0),
               dist.to.closest.TSG = scale(dist.to.closest.TSG, center = 0),
               distance.to.centromere = scale(distance.to.centromere, center = 0),
               distance.to.telomere = scale(distance.to.telomere, center = 0),
               Ess.distance_pancancer = scale(Ess.distance_pancancer, center = 0)
               ) %>% ungroup() 
      
      ml_table_backup <- as.data.frame(apply(ml_table_backup, 2, as.numeric))
    }
    
    ml_table_backup$chr <- NULL
    
    if(cancer_type != 'pancancer'){
      ml_table_backup <- ml_table_backup %>% filter(Type == cancer_type)
    }

    source('/home/ieo7429/Scrivania/scripts/run_models_gab.R', local = T, verbose = F)
    
  }
  
}

