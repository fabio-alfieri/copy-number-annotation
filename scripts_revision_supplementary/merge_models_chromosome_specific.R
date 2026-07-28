# merge objects
rm(list=ls())
gc(full=T)

packages <- c(
  "stringr", "reshape2", "tidyr",
  "ggplot2","fastcluster", "tidyverse"
)

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

model.outputs_ext <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/chromosome_specific/", full.names = TRUE)
shap.outputs_ext <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/SHAP/chromosome_specific/avg_shap/", full.names = TRUE)
shap.outputs_ext <- grep("::", shap.outputs_ext, value = TRUE)

chr_nums <- 1:22
patterns <- c('Arm-level','Chromosome-level','Mid-length','Small-scale','no_cluster')

for(chr_num in chr_nums){
  
  model.outputs <- model.outputs_ext[grepl(pattern = paste0("CS", "_", chr_num, "_"), x = model.outputs_ext)]
  shap.outputs <- shap.outputs_ext[grepl(pattern = paste0("CS", "_", chr_num, "_"), x = shap.outputs_ext)]

  for(i in patterns){
    
    if((i == "Arm-level") & (chr_num %in% c(13,14,15,21,22))){
      next  
    }
    
    print(i)
    
    # Feature Matrix Values
    load(model.outputs[grep(i, model.outputs)][grep('ampl',model.outputs[grep(i, model.outputs)])])
    model_ampl <- Output.regressor; model_ampl <- model_ampl$ampl_score 
    
    load(model.outputs[grep(i, model.outputs)][grep('del',model.outputs[grep(i, model.outputs)])])
    model_del <- Output.regressor; model_del <- model_del$del_score
    
    rm(Output.regressor)
    
    models.X <- list(ampl = rbind(model_ampl[[1]]$X_train, model_ampl[[2]]$X_test),
                     del = rbind(model_del[[1]]$X_train, model_del[[2]]$X_test))
    
    names(models.X) <- c(paste0(i,'::Amplification model'),
                         paste0(i,'::Deletion model'))
    
    # SHAP values
    load(shap.outputs[grep(paste0(i,'::'), shap.outputs)][1])
    shap_ampl <- models.shap.df[[1]] %>% dplyr::select(-BIAS)
    rm(models.shap.df)
    
    load(shap.outputs[grep(paste0(i,'::'), shap.outputs)][2])
    shap_del <- models.shap.df[[1]] %>% dplyr::select(-BIAS)
    rm(models.shap.df, models.X.test, models.X.train, models.shap.values)
    
    models.shap.df <- list(ampl = shap_ampl,
                           del = shap_del)
    
    names(models.shap.df) <- c(paste0(i,'::Amplification model'),
                               paste0(i,'::Deletion model'))
    
    # aggregate and save
    feature_and_shap <- list(models.X,models.shap.df)
    names(feature_and_shap) <- c('models.X','models.shap.df')
    
    write_rds(feature_and_shap, file = paste0("/home/ieo7429/Scrivania/results_regressor_gab/SHAP_and_featurematrix/chromosome_specific/", i, "_", chr_num, "_AmplDel.rds"))
    
    }
  }
