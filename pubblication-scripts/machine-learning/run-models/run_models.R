# rm(list=ls())

mode <- "ONLY_OCCURRENCE"

if(mode == "LOCO"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_LOCO_gab.R'
} else if(mode == "CLASSIC"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_gab.R'
} else if(mode == "NO_SEGMENT_SPECIFIC"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_NSS_gab.R'
} else if(mode == "ONLY_OCCURRENCE"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_OO_gab.R'
} else if(mode == "PROXIMITY"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_PRO_gab.R'
} else if(mode == "LOCO_NSS"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_LOCO_NSS_gab.R'
} else if(mode == "LOFO"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_LOFO_gab.R'
} else if(mode == "PROXY"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_PROXY_gab.R'
} else if(mode == "LOFO_PROXY"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_LOFO_PROXY_gab.R'
} else if(mode == "NO_OG"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_NO_OG_gab.R'
} else if(mode == "LOCTO"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_LOCTO_gab.R'
} else if(mode == "LAFO_PROXY"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_LAFO_PROXY_gab.R'
}


require(xgboost)
require(tidyverse)
library(Matrix)
library(caret)
library(caTools)
library(dplyr)
library(ggplot2)

for(ampl_or_del in c("ampl", "del")){
  ml_table <- ml_table_backup
  
  if(ampl_or_del == 'ampl'){
    label <- 'ampl_score'
  }else if(ampl_or_del == 'del'){
    label <- 'del_score'
  }
  
  cat('Running Regressor for', ampl_or_del, '\n')
  cat('Running ', script, '\n')
  
  source(script, local = T, verbose = F)

}

cat('Done! \n')


