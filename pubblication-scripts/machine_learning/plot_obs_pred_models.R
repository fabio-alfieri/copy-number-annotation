rm(list=ls())
gc(full=T)

library(ggplot2)
library(tidyr)
library(stringr)
library(tidyverse)

setwd("/home/ieo7429/Scrivania/")

model.outputs <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/", full.names = TRUE)
model.outputs <- grep("Output.regressor", model.outputs, value = TRUE)
model.outputs <- model.outputs <- model.outputs[c(93,91,92,63,17,54,5)]


patterns <- c("Arm-level","Chromosome-level","Mid-length","Small-scale","NoCluster", "ampl", "del", "LOCO", as.character(1:22), "NSS", "OO", "PRO", "new", "XL")

for(model_path in model.outputs){

load(model_path)

parts <- strsplit(x = model_path, split = "_")[[1]]
model_name <- paste(parts[parts %in% patterns], collapse = "_")

if(("LOCO" %in% parts) & !("NSS" %in% parts)){
  chr <- as.numeric(parts[5])
  model <- Output.regressor[[1]][[chr]]
} else if(("LOCO" %in% parts) & ("NSS" %in% parts)) {
  chr <- as.numeric(parts[6])
  model <- Output.regressor[[1]][[chr]]
} else {
  model <- Output.regressor[[1]]
}

#if("LOCO" %in% parts){
#  chr <- as.numeric(parts[5])
#  model <- Output.regressor[[1]][[chr]]
#} else {
#  model <- Output.regressor[[1]]
#}

  obs_pred <- model[[1]]$Test.labels
  obs_pred <- obs_pred[order(as.numeric(unlist(lapply(X = obs_pred$bin, FUN = function(x){strsplit(x, split = "_")[[1]][[2]]}))), decreasing = F),]
  obs <- obs_pred[,3]
  pred <- obs_pred[,4]
  
  df <- data.frame(
    index = 1:length(obs),
    obs = obs,
    pred = pred
  )
  
  df_long <- tidyr::pivot_longer(df, cols = c(obs, pred), names_to = "type", values_to = "value")
  
  p <- ggplot(df_long, aes(x = index, y = value, color = type)) +
    geom_line(size = 1) +
    labs(
      x = "Index",
      y = "Value",
      title = paste0("Observed vs Predicted", chr)
    ) +
    theme_minimal()
  
  ggsave(filename = paste0(model_name, ".pdf"), plot = p, width = 30, height = 8)

}




