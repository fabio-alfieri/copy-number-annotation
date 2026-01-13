rm(list=ls())
gc(full=T)

packages <- c(
  "stringr", "parallel", "reshape2", "dplyr",
  "ggplot2", "Matrix", "caret",
  "caTools", "xgboost", "tidyverse"
)

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

wd <- 'path/to/GitHub/copy-number-annotation/'

datasets_path <- paste0(wd, 'data/InteractomeINSIDER/')
chromosome_path <- paste0(wd, 'data/centromeres_and_chromosomes/chromosome.tsv')
centromere_path <- paste0(wd, 'data/centromeres_and_chromosomes/centomere.tsv')
run_models_script <- paste0(wd, 'publication-scripts/machine-learning/run_models.R')

output.files <- list.files(path = datasets_path, full.names = TRUE, recursive = FALSE, pattern = "ampl_or_del_0.1Mbp_covThr_zero")

name <- "ampl_or_del_0.1Mbp_covThr_zero"
level <- "0.1Mbp"
bincov <- "covThr_zero"
folder <- "InteractomeINSIDER"

chromosome <- read.table(file = chromosome_path, header = T)
centromeres <- read.table(file = centromere_path, header = T)

tmp <- full_join(chromosome, centromeres, by = 'chr') %>% 
  dplyr::select(chr, Centromere_Length, Chromosome_Length, Centromere_Type)

tmp$chr <- as.character(tmp$chr)

for(file in output.files){
  load(file)
  
  clus.group <- gsub("_ampl_or_del_0.1Mbp_covThr_zero.RData","",basename(file))
  
  if (clus.group == "Length") {
    
    clusters <- names(ML.Tables[["ampl_or_del"]])
    
  } else if (clus.group == "Length_Location") {
    
    clusters <- names(ML.Tables[["ampl_or_del"]])
    
  } else if (clus.group == "NoCluster") {
    
    ML.Tables$ampl_or_del <- list(no_cluster = ML.Tables$ampl_or_del)
    
    clusters <- "no_cluster"
    
  }
  
  for(classS in clusters){
    
    ml_table_backup <- ML.Tables[["ampl_or_del"]][[classS]]
    ml_table_backup$bin <- as.character(ml_table_backup$bin)

    ml_table_backup$chr <- do.call(rbind, str_split(ml_table_backup$bin,'_'))[,1]
    ml_table_backup <- left_join(ml_table_backup, tmp, by = 'chr')
    ml_table_backup$chr <- NULL
    
    source(run_models_script, local = T, verbose = F)
    
  }
  
}

