rm(list=ls())
gc(full=T)

packages <- c(
  "stringr", "tidyr", "reshape2", "dplyr",
  "ggplot2", "colorspace", "ggsignif",
  "disttree", "fitdistrplus", "tidyverse",
  "patchwork"
)

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

#######################  load helpers ####################### 

source("./get_residuals.R")
source("./change_distance.R")
source("./define_annotation_rules_correlation_based.R")
source("./do_pancancer_barplots_and_boxplots.R")
source("./do_tumor_specific_barplots_and_boxplots.R")
source("./do_do_tumor_specific_barplots_and_boxplots.R")
source("./do_top1_annotation.R")
source("./do_plot_grouped_by_chr_type.R")
source("./do_plot_grouped_by_chr.R")
source("./do_plot_grouped_by_type.R")
source("./plotting_helpers.R")
source("./landscape_plot_observed_prediction.R")
source("./do_plot_residuals.R")

####################### 

#######################  read model's output ####################### 

model_classes <- c('Mid-length', 'Arm-level', 'Chromosome-level', 'Small-scale')

for (model_class in model_classes) {
  
shap_and_feature_path <- paste0("../Data/merged/SHAP_and_FeatureMatrix_",model_class,"_AmplDel.rds")
pred_ampl_path <- paste0("../Data/InteractomeINSIDER/",model_class,"-pred_ampl.rds")
pred_del_path <- paste0("../Data/InteractomeINSIDER/",model_class,"-pred_del.rds")
correlations_table_path <- paste0("../Data/results_regressor/final_verdict.tsv")

df <- readRDS(file = shap_and_feature_path)
pred_ampl <- readRDS(file = pred_ampl_path) 
pred_del <- readRDS(file = pred_del_path)
correlations_table <- read.table(file = correlations_table_path, header = T, sep = "\t")

#######################  manipulate data and define variables ####################### 

res.ampl <- get_residuals(pred_ampl); res.ampl$labels <- paste0(res.ampl$bin, "-", res.ampl$Type) # load ampl residuals
res.del <- get_residuals(pred_del); res.del$labels <- paste0(res.del$bin, "-", res.del$Type)      # load del residuals

bool.plot <- F

output <- list() # define output list

model_types <- c("ampl", "del") # define model types

for (i in model_types) { # iterate over model types

message(paste0("Now processing ", model_class, " ", i, " model"))

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

col.sum.zero <- colnames(shap %>% dplyr::select(-labels))[apply(shap %>% dplyr::select(-labels), 2, sum) == 0]
non.segment.features <- c("Chromosome_Length","Centromere_Length","Centromere_Type")

if(model_class %in% c('Mid-length','Small-scale')){
  
  cols.to.discard <- c("labels","BIAS", col.sum.zero, non.segment.features) # include non.segment.features if the model is != Arm- or Chr-level

  }else{
  
  cols.to.discard <- c("labels","BIAS")
  
  }

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


pval_thr <- 0.6                                                                              # select threshold to discriminate "correct" vs "incorrect" prediction
diploidy_threshold <- median(res$prediction)                                                 # select threshold to discriminate "frequently diploid" vs "frequently non diploid"

ns_thresholds <- data.frame(model_class = c("Mid-length", "Mid-length",
                                            "Arm-level", "Arm-level",
                                            "Chromosome-level", "Chromosome-level",
                                            "Small-scale", "Small-scale",
                                            "no_cluster", "no_cluster"),
                            model_type = rep(c("ampl", "del"),5),
                            thresholds = c(0.02, 0.025, 0.002, 0.010, 0.01, 0, 0.0004, 0, 0.04, 0.04))

ns_threshold <- ns_thresholds[(ns_thresholds$model_class == model_class & ns_thresholds$model_type == i),"thresholds"]

print(ns_threshold)

res$pred_is_correct <- ifelse(res$pval.residual > pval_thr, T, F)

res$region_is_diploid <- ifelse(res$observed > diploidy_threshold, F, T)

####################### 

####################### examine residuals ####################### 

if (bool.plot) {
  
  residuals_plots <- do_plot_residuals(res, save = F)
  
}

####################### 

#######################  define ruling system for annotation ####################### 

define_annotation_rules_correlation_based(i, correlations_table)

####################### 

#######################  top1 annotation #######################

message(paste0("Performing annotation based on most important feature"))

res_annot <- t(apply(X = actual_shap_mat, 
                     MARGIN = 1, 
                     FUN = do_top1_annotation))

res_annot <- cbind(res_annot, shap_agg$labels)
colnames(res_annot) <- c(paste0("top",1:num_of_top_features), paste0("shap_top",1:num_of_top_features), "labels")
res_annot <- as.data.frame(res_annot)

merged_res_annot <-  merge(res, res_annot, on = "labels") # merge with residuals

merged_res_annot$type <- do.call(rbind,str_split(merged_res_annot$labels,'-'))[,2] # take tumor type

merged_res_annot$is_occurrence <- ifelse(merged_res_annot$top1 %in% features_occurrence, TRUE, FALSE) 
merged_res_annot$is_selection <- ifelse(merged_res_annot$top1 %in% features_selection, TRUE, FALSE) 

merged_res_annot$is_positive <- ifelse(merged_res_annot$top1 %in% features_positive_selection, TRUE, FALSE)
merged_res_annot$is_negative <- ifelse(merged_res_annot$top1 %in% features_negative_selection, TRUE, FALSE)

merged_res_annot$is_relaxed_positive_selection <- ifelse(merged_res_annot$top1 %in% features_relaxed_positive_selection, TRUE, FALSE) # remove this category from here
merged_res_annot$is_relaxed_negative_selection <- ifelse(merged_res_annot$top1 %in% features_relaxed_negative_selection, TRUE, FALSE)

merged_res_annot$annot <- ifelse(merged_res_annot$is_occurrence, "Occurrence",
                                 ifelse(merged_res_annot$is_positive, "Positive Selection",
                                        ifelse(merged_res_annot$is_negative, "Negative Selection",
                                               ifelse(merged_res_annot$is_relaxed_positive_selection, "Relaxed Positive Selection",
                                                      ifelse(merged_res_annot$is_relaxed_negative_selection, "Relaxed Negative Selection",
                                                             "Unknown")
                                               )
                                        )
                                 )
)

merged_res_annot$annot_3_classes <- ifelse(merged_res_annot$pred_is_correct,
                                           merged_res_annot$annot, 
                                           ifelse(merged_res_annot$observed > merged_res_annot$prediction, 
                                                  'Positive Selection',
                                                  'Negative Selection'))

merged_res_annot$annot_5_classes <- ifelse(merged_res_annot$pred_is_correct, 
                                           merged_res_annot$annot, 
                                           ifelse(merged_res_annot$observed > merged_res_annot$prediction, 
                                                  'Incorrect prediction - Positive Selection',
                                                  'Incorrect prediction - Negative Selection'))

merged_res_annot$top1 <- ifelse(!(merged_res_annot$pred_is_correct), merged_res_annot$annot_5_classes, merged_res_annot$top1)

merged_res_annot$top1 <- ifelse(merged_res_annot$top1 %in% features_occurrence, 
                                gsub(x = merged_res_annot$top1, pattern = "-", replacement = ""), 
                                merged_res_annot$top1)

merged_res_annot$annot <- NULL
merged_res_annot$is_occurrence <- NULL
merged_res_annot$is_selection <- NULL
merged_res_annot$is_positive <- NULL
merged_res_annot$is_negative <- NULL
merged_res_annot$is_relaxed_positive_selection <- NULL
merged_res_annot$is_relaxed_negative_selection <- NULL

#######################

compute_feats_from_shap_mat <- function(row){
  
  features <- colnames(actual_shap_mat);
  signs <- sign(row); signs <- gsub(pattern = "1", replacement = "", signs)
  features <- paste0(signs, features)
  
  return(features)
  
}

compute_sas_ns <- function(idx){
  
  feat_row <- feats[idx,]
  shap_vals_row <- actual_shap_mat[idx,]
  
  idx_ns <- which(feat_row %in% features_negative_selection)

  shap_vals_row_ns <- shap_vals_row[idx_ns]
  sas_ns <- sum(abs(shap_vals_row_ns))
  
  return(sas_ns)
  
}

feats <- t(apply(X = actual_shap_mat, 
                 MARGIN = 1, 
                 FUN = compute_feats_from_shap_mat))

labels_true_negative_sel <- merged_res_annot[(merged_res_annot$pred_is_correct == F &
                                                  merged_res_annot$annot_5_classes == 'Incorrect prediction - Negative Selection'),]$labels

labels_candidate_negative_sel <- merged_res_annot[(merged_res_annot$pred_is_correct == T &
                                                     merged_res_annot$region_is_diploid == T &
                                                       merged_res_annot$annot_5_classes == "Negative Selection"),]$labels


merged_res_annot$is_candidate <- ifelse(merged_res_annot$labels %in% labels_candidate_negative_sel, T, F)

sas_ns <- sapply(X = 1:nrow(actual_shap_mat), FUN = compute_sas_ns, simplify = T)
names(sas_ns) <- shap_agg$labels

sas_ns_tns <- sas_ns[labels_true_negative_sel]
sas_ns_cns <- sas_ns[labels_candidate_negative_sel]

plot_df_combined_sas_ns <- bind_rows(
  data.frame(col = sas_ns_tns) %>%
    mutate(Group = "True Negative"),
  data.frame(col = sas_ns_cns) %>%
    mutate(Group = "Candidate Negative")
)

hist_sas_ns <- ggplot(plot_df_combined_sas_ns, aes(x = col, fill = Group)) +
  geom_vline(xintercept = ns_threshold) +
  geom_histogram(bins = 75, color = "white", alpha = 0.5, position = "identity") +
  labs(
    x = "sum(|SHAP_negative_selection|)",
    y = "Frequency",
    title = paste0(model_class, "_", i, "_", diploidy_threshold
                   ),
    fill = "Group"
  ) +
  scale_fill_manual(values = c(
    "True Negative" = "blue",
    "Candidate Negative" = "red"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

filename_sas_ns <- paste0("..Data/plots/hist_sas_ns_", model_class, "_", i, "_"
                          , diploidy_threshold
                          , "_comparison.pdf")

if(!(file.exists(filename_sas_ns))){
  ggsave(filename = filename_sas_ns, plot = hist_sas_ns, width = 8, height = 6, units = "in")
  }

if(all(names(sas_ns) == merged_res_annot$labels)){
  merged_res_annot$sas_ns <- sas_ns
  merged_res_annot$is_strong_candidate <- ifelse((merged_res_annot$is_candidate == T & 
                                                    merged_res_annot$sas_ns > ns_threshold), T, F)
}


merged_res_annot$annot_final <- NA

bool1 <- (merged_res_annot$pred_is_correct == FALSE) # first case: divide correct and incorrect

bool2 <- (merged_res_annot$pred_is_correct == TRUE & # second case: correct predictions, not diploid, annotated for negative --> error
            merged_res_annot$region_is_diploid == FALSE & 
            merged_res_annot$annot_5_classes == "Negative Selection")

bool3 <- (merged_res_annot$pred_is_correct == TRUE & # third case: correct predictions, not diploid, not annotated for negative --> occurrence or positive
            merged_res_annot$region_is_diploid == FALSE & 
            merged_res_annot$annot_5_classes != "Negative Selection")

bool4 <- (merged_res_annot$pred_is_correct == TRUE & # fourth case: correct predictions, diploid, not annotated for negative --> error
           merged_res_annot$region_is_diploid == TRUE  & 
           merged_res_annot$is_candidate == FALSE)

bool5 <- (merged_res_annot$pred_is_correct == TRUE & # fifth case: correct predictions, diploid, annotated for negative with high sum abs --> negative
            merged_res_annot$region_is_diploid == TRUE  & 
            merged_res_annot$is_candidate == TRUE & 
            merged_res_annot$is_strong_candidate == TRUE)

bool6 <- (merged_res_annot$pred_is_correct == TRUE & # sizth case: correct predictions, diploid, annotated for negative with low sum abs --> error
            merged_res_annot$region_is_diploid == TRUE  & 
            merged_res_annot$is_candidate == TRUE & 
            merged_res_annot$is_strong_candidate == FALSE)

bool7 <- (merged_res_annot$annot_5_classes %in%  # seventh case: anything that is relaxed selection becomes no detectable forces
            c("Relaxed Positive Selection", "Relaxed Negative Selection"))


if(as.logical(sum(bool1))){
  merged_res_annot[bool1,]$annot_final <- merged_res_annot[bool1,"annot_5_classes"]
}
 
if(as.logical(sum(bool2))){
  merged_res_annot[bool2,]$annot_final <- "No Detectable Force"
}
  
if(as.logical(sum(bool3))){
  merged_res_annot[bool3,]$annot_final <- merged_res_annot[bool3, "annot_5_classes"] 
}  

if(as.logical(sum(bool4))){
  merged_res_annot[bool4,]$annot_final <- "No Detectable Force"
}  

if(as.logical(sum(bool5))){
  merged_res_annot[bool5,]$annot_final <- merged_res_annot[bool5, "annot_5_classes"]
}

if(as.logical(sum(bool6))){
  merged_res_annot[bool6,]$annot_final <- "No Detectable Force"
}

if(as.logical(sum(bool7))){
  merged_res_annot[bool7,]$annot_final <- "No Detectable Force"
}

merged_res_annot$chr <- sapply(X = strsplit(x = merged_res_annot$labels, split = "[_-]"), FUN = `[`, 1)

filename_annot <- paste0("../Data/annotation/merged_res_annot_", model_class, "_", i, "_", pval_thr, "_", diploidy_threshold, "_", ns_threshold, ".tsv")
write.table(x = merged_res_annot, file = filename_annot, quote = F, sep = "\t", row.names = F, col.names = T)
      











    
    
    