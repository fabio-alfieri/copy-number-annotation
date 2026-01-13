
mode <- "CLASSIC"

if(mode == "LOCO"){
  script <- './0_prepareRegressor_LOCO.R'
} else if(mode == "CLASSIC"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor.R'
} else if(mode == "NO_SEGMENT_SPECIFIC"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_NSS.R'
} else if(mode == "ONLY_OCCURRENCE"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_OO.R'
} else if(mode == "PROXIMITY"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_PRO.R'
} else if(mode == "LOCO_NSS"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_LOCO_NSS.R'
} else if(mode == "LOFO"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_LOFO.R'
} else if(mode == "PROXY"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_PROXY.R'
} else if(mode == "LOFO_PROXY"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_LOFO_PROXY.R'
} else if(mode == "NO_OG"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_NO_OG.R'
} else if(mode == "LOCTO"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_LOCTO.R'
} else if(mode == "LAFO_PROXY"){
  script <- '/home/ieo7429/Scrivania/scripts/0_prepareRegressor_LAFO_PROXY.R'
}

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


