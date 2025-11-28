rm(list=ls())
gc(full=T)

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
source("scripts/define_annotation_rules_correlation_based.R")
source("scripts/do_pancancer_barplots_and_boxplots.R")
source("scripts/do_tumor_specific_barplots_and_boxplots.R")
source("scripts/do_do_tumor_specific_barplots_and_boxplots.R")
source("scripts/do_top1_annotation.R")
source("scripts/compare_ratio_based_top1.R")
source("scripts/do_plot_grouped_by_chr_type.R")
source("scripts/do_plot_grouped_by_chr.R")
source("scripts/do_plot_grouped_by_type.R")
source("scripts/plotting_helpers.R")
source("scripts/landscape_plot_observed_prediction.R")
source("scripts/do_plot_residuals.R")

####################### 

#######################  read model's output ####################### 

model_classes <- c("Mid-length", 'Chromosome-level', 'Arm-level', 'Small-scale', 'no_cluster')

for (model_class in model_classes) {

shap_and_feature_path <- paste0("dev/Data/SHAP_and_FeatureMatrix_",model_class,"_AmplDel.rds")
pred_ampl_path <- paste0("results_regressor_gab/InteractomeINSIDER/",model_class,"-pred_ampl.rds")
pred_del_path <- paste0("results_regressor_gab/InteractomeINSIDER/",model_class,"-pred_del.rds")
correlations_table_path <- paste0("results_regressor_gab/final_verdict.tsv")

df <- readRDS(file = shap_and_feature_path)
pred_ampl <- readRDS(file = pred_ampl_path) 
pred_del <- readRDS(file = pred_del_path)
correlations_table <- read.table(file = correlations_table_path, header = T, sep = "\t")

####################### 

setwd("dev/imgs/")

#######################  manipulate data and define variables ####################### 

res.ampl <- get_residuals(pred_ampl); res.ampl$labels <- paste0(res.ampl$bin, "-", res.ampl$Type) # load ampl residuals
res.del <- get_residuals(pred_del); res.del$labels <- paste0(res.del$bin, "-", res.del$Type) # load del residuals

output <- list()

model_types <- c("ampl", "del")

for (i in model_types) { # iterate over model types

message(paste0("Now processing ", model_class, " ", i, " model"))

if(i == 'ampl'){
  
  shap <- df$models.shap.df[[grep('Amplification',names(df$models.shap.df))]]
  values <- df$models.X[[grep('Amplification',names(df$models.X))]]
  res <- res.ampl[,c("labels", "observed", "prediction", "residual")] 
  score_type <- "ampl_score"
  pred <- pred_ampl
  
} else if(i == 'del'){
  
  shap <- df$models.shap.df[[grep('Deletion',names(df$models.shap.df))]]
  values <- df$models.X[[grep('Deletion',names(df$models.X))]]
  res <- res.del[,c("labels", "observed", "prediction", "residual")]
  score_type <- "del_score"
  pred <- pred_del
  
}

values$labels <- paste0(values$bin,'-',values$Type)
values <- change_distance(values)

col.sum.zero <- colnames(shap %>% dplyr::select(-labels))[apply(shap %>% dplyr::select(-labels), 2, sum) == 0]
non.segment.features <- c("Chromosome_Length","Centromere_Length","Centromere_Type")

if(model_class %in% c('Mid-length','Small-scale')){cols.to.discard <- c("labels","BIAS", col.sum.zero, non.segment.features)
}else{cols.to.discard <- c("labels","BIAS")}

cols.to.keep <- colnames(shap)[!colnames(shap) %in% cols.to.discard]
cols.to.keep <- c(cols.to.keep,"labels")

shap_agg <- shap[, cols.to.keep] 
values_agg <- values[, cols.to.keep]
values_agg <- values_agg[match(shap_agg$labels, values_agg$labels), ]

actual_shap_mat <- shap_agg[,-ncol(shap_agg)]
actual_values_mat <- values_agg[,-ncol(values_agg)]

distfit_output <- distfit(y = res$residual)
mu <- distfit_output$par["mu"]
sigma <- distfit_output$par["sigma"]
pvals <- 2 * pmin(pnorm(res$residual, mean = mu, sd = sigma), 1 - pnorm(res$residual, mean = mu, sd = sigma)) # take p-values
res$pval.residual <- pvals
pval_thr <- 0.6
res$pred_is_correct <- ifelse(res$pval.residual > pval_thr, T, F)

labels_to_keep <- res[res$pred_is_correct == T, "labels"]
values_agg_filt <- values_agg[values_agg$labels %in% labels_to_keep$labels,]
shap_agg_filt <- shap_agg[shap_agg$labels %in% labels_to_keep$labels,]

vars <- colnames(shap_agg_filt)[-ncol(shap_agg_filt)]
method <- "pearson"

for (var in vars) {
print(var)
labels_pos <- shap_agg_filt[shap_agg_filt[,var] > 0, "labels"]
shap_agg_pos <- shap_agg_filt[shap_agg_filt$labels %in% labels_pos,]
values_agg_pos <- values_agg_filt[values_agg_filt$labels %in% labels_pos,]

labels_neg <- shap_agg_filt$labels[!shap_agg_filt$labels %in% labels_pos]
shap_agg_neg <- shap_agg_filt[shap_agg_filt$labels %in% labels_neg,]
values_agg_neg <- values_agg_filt[values_agg_filt$labels %in% labels_neg,]

plot_with_cor <- function(shap, vals, var, method = "pearson") {

  df <- data.frame(x = shap, y = vals)
  df <- df[is.finite(df$x) & is.finite(df$y), ]
  
  if (nrow(df) == 0) {
    stop("Dopo il filtraggio, non ci sono dati validi per il plot.")
  }
  
  cor_res <- cor.test(df$x, df$y, method = method, exact = FALSE)
  
  corr_est <- round(cor_res$estimate, 3)
  p_val <- signif(cor_res$p.value, 3)
  label_text <- paste0(method, ": r = ", corr_est, ", p = ", p_val)
  
  if (length(unique(sign(df$x))) > 1) {
    to.attach <- ""
  } else {
    to.attach <- gsub(pattern = "1", replacement = "", 
                    x = as.character(unique(sign(df$x))))
  }
  
  
  plot_title <- paste0(to.attach, var)
  
  g <- ggplot(df, aes(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    theme_minimal() +
    labs(
      x = "SHAP", 
      y = "Value", 
      title = plot_title,
      subtitle = label_text
    )
  
  return(g)
}

pos_bool <- ((length(shap_agg_pos[,var]) > 0) & ((length(values_agg_pos[,var]) > 0)))
neg_bool <- ((length(shap_agg_neg[,var]) > 0) & ((length(values_agg_neg[,var]) > 0)))

if (pos_bool){
  p_pos <- invisible(plot_with_cor(shap_agg_pos[,var], values_agg_pos[,var], var = var, method = method))
}

if (neg_bool){
  p_neg <- invisible(plot_with_cor(shap_agg_neg[,var], values_agg_neg[,var], var = var, method = method))
}

if(pos_bool & neg_bool){
p_sep <- invisible(gridExtra::grid.arrange(p_neg, p_pos, ncol = 2))

ggsave(
  filename = paste0("p_sep_", model_class, "_", i, "_", var, "_", method, ".pdf"),
  plot = p_sep,
  width = 8, height = 6, units = "in"
)

}

p_tot <- invisible(plot_with_cor(shap_agg_filt[,var], values_agg_filt[,var], var = var, method = method))

ggsave(
  filename = paste0("p_tot_", model_class, "_", i, "_", var, "_", method, ".pdf"),
  plot = p_tot,
  width = 8, height = 6, units = "in"
)

  }
}

setwd("../../")

}



