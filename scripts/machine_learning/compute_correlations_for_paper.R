# compute correlations across models
rm(list=ls()) # clean variables
gc(full=T) # garbage collector

#######################  import libraries #######################

library(ggplot2)
library(tidyr)
library(stringr)
library(tidyverse)
library(disttree)
library(fitdistrplus)
library(colorspace)
library(dplyr)
library(ggsignif)
library(patchwork)

#######################

#######################  set working directory ####################### 

setwd("/home/ieo7429/Scrivania/")

####################### 

#######################  load helpers ####################### 

source("scripts/get_residuals.R")
source("scripts/compute_ratios.R")
source("scripts/change_distance.R")
source("scripts/aggregate_chromatin_states.R")
source("scripts/define_annotation_rules.R")

####################### 

#######################  read model's output ####################### 

load("dev/Data/All_levels_backbonetables.RData")

model_classes <- c("Chromosome-level", "Arm-level", "Mid-length", "Small-scale", "no_cluster")
model_classes <- c("Mid-length")

cor.final <- data.frame()

for(model_class in model_classes){
  
  shap_and_feature_path <- paste0("dev/Data/SHAP_and_FeatureMatrix_",model_class,"_AmplDel.rds")
  pred_ampl_path <- paste0("results_regressor_gab/InteractomeINSIDER/",model_class,"-pred_ampl.rds")
  pred_del_path <- paste0("results_regressor_gab/InteractomeINSIDER/",model_class,"-pred_del.rds")
  
  df <- readRDS(file = shap_and_feature_path)
  pred_ampl <- readRDS(file = pred_ampl_path) 
  pred_del <- readRDS(file = pred_del_path)
  
  #######################  manipulate data and define variables ####################### 
  
  res.ampl <- get_residuals(pred_ampl); res.ampl$labels <- paste0(res.ampl$bin, "-", res.ampl$Type) # load ampl residuals
  res.del <- get_residuals(pred_del); res.del$labels <- paste0(res.del$bin, "-", res.del$Type) # load del residuals
  
  bool.plot <- F
  
  output <- list() # define output list
  
  model_types <- c("ampl", "del") # define model types
  
  for (i in model_types) { # iterate over model types
  
  # i <- model_types[1]
  
  message(paste("Now processing", model_class, i, "model"))
  
  if(i == 'ampl'){ # if model is amplification 
    
    shap <- df$models.shap.df[[grep('Amplification',names(df$models.shap.df))]] # take shap values
    values <- df$models.X[[grep('Amplification',names(df$models.X))]] # take actual values
    res <- res.ampl[,c("labels", "observed", "prediction", "residual")] # take residuals 
    score_type <- "ampl_score" # define the target variable
    pred <- pred_ampl # take predictions
    
  } else if(i == 'del'){ # if model is deletion, do the same
    
    shap <- df$models.shap.df[[grep('Deletion',names(df$models.shap.df))]]
    values <- df$models.X[[grep('Deletion',names(df$models.X))]]
    res <- res.del[,c("labels", "observed", "prediction", "residual")]
    score_type <- "del_score"
    pred <- pred_del
    
  }
  
  values$labels <- paste0(values$bin,'-',values$Type)
  values <- change_distance(values)
  
  # cols.to.discard <- c("Chromosome_Length","Centromere_Length","Centromere_Type", # select cols to discard
  #                      "BIAS","density.OG","density.TSG","labels","ESSscore_pancancer",
  #                      "HAPLOscore_pancancer","Density.complex.proteins", "Density.Ohnologs")
  
  col.sum.zero <- colnames(shap %>% dplyr::select(-labels))[apply(shap %>% dplyr::select(-labels), 2, sum) == 0]
  non.segment.features <- c("Chromosome_Length","Centromere_Length","Centromere_Type")
  cols.to.discard <- c("labels","BIAS", col.sum.zero) # include non.segment.features if the model is != Arm- or Chr-level
  
  cols.to.keep <- colnames(shap)[!colnames(shap) %in% cols.to.discard] # extract cols to keep
  cols.to.keep <- c(cols.to.keep,"labels")
  
  shap_agg <- shap[, cols.to.keep] # subset for cols to keep 
  values_agg <- values[, cols.to.keep] # subset for cols to keep
  values_agg <- values_agg[match(shap_agg$labels, values_agg$labels), ]
  
  actual_shap_mat <- shap_agg[,-ncol(shap_agg)] # remove labels (last column)
  actual_values_mat <- values_agg[,-ncol(values_agg)] # remove labels (last column)
  
  distfit_output <- distfit(y = res$residual) # fit distribution of residuals
  mu <- distfit_output$par["mu"] # take mean
  sigma <- distfit_output$par["sigma"] # take standard deviation
  pvals <- 2 * pmin(pnorm(res$residual, mean = mu, sd = sigma), 1 - pnorm(res$residual, mean = mu, sd = sigma)) # take p-values
  res$pval.residual <- pvals # assign p-values
  pval_thr <- 0.6 # select threshold to discriminate "correct" vs "incorrect" prediction
  res$pred_is_correct <- ifelse(res$pval.residual > pval_thr, T, F) # binarize correct prediction
  
  ####################### 
  
  good.segments <- res %>% filter(pval.residual >= pval_thr)
  
  actual_shap_mat <- cbind(actual_shap_mat, labels = shap_agg$labels) %>% 
    as_tibble() %>% 
    filter(labels %in% good.segments$labels) %>% 
    dplyr::select(-labels)
  
  actual_values_mat <- cbind(actual_values_mat, labels = shap_agg$labels) %>% 
    as_tibble() %>% 
    filter(labels %in% good.segments$labels) %>% 
    dplyr::select(-labels)
  
  common_cols <- intersect(colnames(actual_shap_mat), colnames(actual_values_mat)) 
  
  cor_results <- lapply(common_cols, function(col) {
    corP <- cor.test(actual_shap_mat[[col]], actual_values_mat[[col]], method = "pearson", na.action = "na.omit")
    corS <- cor.test(actual_shap_mat[[col]], actual_values_mat[[col]], method = "spearman", na.action = "na.omit", exact = F)
    
    cor.out.t <- cbind(estimate.P = corP$estimate, pval.P = corP$p.value, estimate.S = corS$estimate, pval.S = corS$p.value)
    
    plot_df <- data.frame(
      SHAP = actual_shap_mat[[col]],
      Feature = actual_values_mat[[col]]
    ) %>%
      filter(is.finite(SHAP) & is.finite(Feature))
    
    cor.out.p <- ggplot(plot_df, aes(x = Feature, y = SHAP)) +
      geom_point(alpha = 0.6, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "darkred") +
      theme_minimal() +
      labs(
        title = paste("SHAP vs Feature:", col),
        x = "Feature Value",
        y = "SHAP Value",
        caption = paste0(
          "Pearson: r = ", round(corP$estimate, 2), ", p = ", signif(corP$p.value, 3), "\n",
          "Spearman: rho = ", round(corS$estimate, 2), ", p = ", signif(corS$p.value, 3)
        )
      ) +
      geom_density_2d(color="black") +
      theme(
        plot.caption = element_text(hjust = 0, face = "italic"),
        plot.title = element_text(face = "bold")
      )
    
    
    outlist <- list(table = cor.out.t, plot = cor.out.p)
    
    return(outlist)
    
  })
  
  cor_tables <- lapply(cor_results, function(x) x$table)
  
  cor_df <- data.frame(Column = common_cols, do.call(rbind, cor_tables))
  cor_df$Column <- factor(cor_df$Column, levels = cor_df[order(cor_df$estimate.S, decreasing = T),'Column'])
  
  features_occurrence <- c("mean.GC.content", "Length_Counts.E1", "Length_Counts.E2",
                            "Length_Counts.E3", "Length_Counts.E4", "Length_Counts.E5",
                            "Length_Counts.E6", "Length_Counts.E7", "Length_Counts.E8",
                            "Length_Counts.E9", "Length_Counts.E10", "Length_Counts.E11",
                            "Length_Counts.E12", "Length_Counts.E13", "Length_Counts.E14",
                            "Length_Counts.E15", "Length_Counts.E16", "Length_Counts.E17",
                            "Length_Counts.E18", "Length_Counts.E19", "Length_Counts.E20",
                            "Length_Counts.E21", "Length_Counts.E22", "Length_Counts.E23",
                            "Length_Counts.E24", "Length_Counts.E25", "distance.to.centromere",
                            "distance.to.telomere", "dist.to.closest.FGS", "Centromere_Length",
                            "Chromosome_Length", "Centromere_Type")
  
  cor_df$feature_type <- ifelse(cor_df$Column %in% features_occurrence, "Occurrence", "Selection")
  
  cor_df$signif <- ifelse(cor_df$pval.S < 0.001, "***",
                          ifelse(cor_df$pval.S < 0.01, "**",
                                 ifelse(cor_df$pval.S < 0.05, "*", "")))
  
  cor_df$pval.P_fdr <- p.adjust(cor_df$pval.P, method = 'bonferroni') # ho provato a vedere se facendo multiple testing correction perdo significatività su OG, ma purtroppo no! :(
  cor_df$model <- model_class
  cor_df$ampl_del <- i
  
  cor.final <- rbind(cor.final, cor_df)
  
  write.table(file = "/home/ieo7429/Scrivania/results_regressor_gab/cor_final.tsv", x = cor.final, sep = "\t", quote = F, row.names = F, col.names = T)
  
  #cor_df$is_signif_both <- ifelse(cor_df$estimate.P > 0, 'pos', 'neg')
  
  # pS <- ggplot(cor_df, aes(x = Column, y = estimate.S)) +
  #   geom_col(aes(fill = feature_type)) +
  #   geom_text(aes(
  #     label = paste0(round(estimate.S, 2), signif)
  #   )) +
  #   theme_minimal() +
  #   labs(
  #     title = "Correlation Between SHAP and Feature Matrix",
  #     y = "Spearman Correlation"
  #   ) +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # 
  # 
  # cor_df$Column <- factor(cor_df$Column, levels = cor_df[order(cor_df$estimate.P, decreasing = T),'Column'])
  # 
  # pP <- ggplot(cor_df, aes(x = Column, y = estimate.P)) +
  #   geom_col(aes(fill = feature_type)) +
  #   geom_text(aes(label = paste0(round(estimate.P, 2), signif), vjust = -0.5)) +
  #   theme_minimal() +
  #   labs(title = "Correlation Between SHAP and Feature Matrix",
  #        y = "Pearson Correlation") +
  #   theme(axis.text = element_text(angle = 90))
  # 
  # ggsave(
  #   filename = paste0(model_class,'::',ifelse(i=='ampl','AmplificationModel','DeletionModel'), "_pancancer_correlations_SHAP_actualValues_Pearson.pdf"),
  #   plot = pP, device = "pdf",
  #   width = 15, height = 10
  # )
  # ggsave(
  #   filename = paste0(ifelse(i=='ampl','AmplificationModel','DeletionModel'), "_pancancer_correlations_SHAP_actualValues_Spearman.pdf"),
  #   plot = pS, device = "pdf",
  #   width = 15, height = 10
  # )
  }
}
