rm(list=ls())
gc(full=T)

library(ggplot2)
library(tidyr)
library(stringr)
library(tidyverse)

setwd("/home/ieo7429/Scrivania/")

model.outputs <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/", full.names = TRUE)
model.outputs <- grep(pattern = "OO", x = model.outputs, value = T)

patterns <- c("Arm-level",
              "Chromosome-level",
              "Mid-length",
              "Small-scale",
              "NoCluster", 
              "ampl", "del", 
              "LOCO", 
              as.character(0:2000),
              "LOCTO",
              c("BRCA", "COADREAD", "ESCA", 
                "GBMLGG", "KIRC", "KIRP", 
                "LUAD", "LUSC", "OV", 
                "PAAD", "STAD"),
              "NSS", 
              "OO", 
              "PRO", 
              "new", 
              "XL", 
              "PROXY", 
              "LOFO",
              "LAFO",
              "Chromosome", 
              "Centromere", 
              "Centromere", 
              "distance.to.telomere", 
              "distance.to.centromere", 
              "Length",
              "NO",
              "OG",
              "Type")

avg_performances <- lapply(X = model.outputs, 
       FUN = function(model_path){
  
  load(model_path)
        
  parts <- strsplit(x = model_path, split = "_")[[1]]
  model_name <- paste(parts[parts %in% patterns], collapse = "_")
  
  if(("LOCO" %in% parts) & !("NSS" %in% parts)){
  chr <- as.numeric(parts[5])
  model <- Output.regressor[[1]][[chr]]
  } else if(("LOCO" %in% parts) & ("NSS" %in% parts)) {
  chr <- as.numeric(parts[6])
  model <- Output.regressor[[1]][[chr]]
  } else if("LOCTO" %in% parts){
  ct <- parts[5]
  model <- Output.regressor[[1]][[ct]]
  } else {
  model <- Output.regressor[[1]]
  }
  
  r2_array <- unlist(lapply(model, function(x){x$Performance$R2}))
  rmse_array <- unlist(lapply(model, function(x){x$Performance$RMSE}))
  pearson_array <- unlist(lapply(model, function(x){x$Performance$Pearson}))
  spearman_array <- unlist(lapply(model, function(x){x$Performance$Spearman}))
  
  out.list <- list(avg_r2 = mean(r2_array),
                   avg_rmse = mean(rmse_array),
                   avg_pearson = mean(pearson_array),
                   avg_spearman = mean(spearman_array)
                   )
  
  out.list <- list(out.list)
  
  names(out.list) <- model_name
  print(out.list[[1]]$avg_r2)
  
  print(paste0("processed ", model_name))
  
  return(out.list)
  
})

avg_performance_table <- do.call(what = rbind, 
                                 lapply(X = avg_performances, 
                                        FUN = function(sublist){
                                          
                                          row.name <- names(sublist)
                                          vals <- unlist(sublist[[1]])
                                          row.df <- t(setNames(data.frame(vals), row.name))
                                          
                                        }))

to.discard <- c("PRO_500_200_Length_Mid-length_ampl", "PRO_50_20_Length_Mid-length_ampl", "PRO_5_2_Length_Mid-length_ampl", "PRO_250_100_Length_Mid-length_ampl", "PRO_1000_400_Length_Mid-length_ampl", "NSS_Length_Mid-length_ampl")
avg_performance_table <- avg_performance_table[!(rownames(avg_performance_table) %in% to.discard),]

write.table(x = avg_performance_table, 
            file = "results_regressor_gab/avg_performances_OO.tsv",
            quote = F, sep = "\t", row.names = T, col.names = T)










