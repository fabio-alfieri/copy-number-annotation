wd <- 'path/to/GitHub/copy-number-annotation/'

rm(list=ls())
gc(full=T)

packages <- c("stringr", "tidyr", "ggplot2", "tidyverse")

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

model.outputs <- list.files(path = paste0(wd, "data/InteractomeINSIDER/"), full.names = TRUE)

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

avg_performance_table <- avg_performance_table[!(rownames(avg_performance_table) %in% to.discard),]

avg_performance_path <- paste0(wd, "data/results_regressor/avg_performances_complete.tsv")

write.table(x = avg_performance_table, 
            file = avg_performance_path,
            quote = F, sep = "\t", row.names = T, col.names = T)










