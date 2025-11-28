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

load("dev/Data/All_levels_backbonetables.RData")

model_class <- 'Mid-length'
# 'Chromosome-level' 'Arm-level' 'Mid-length' 'Small-scale' 'no_cluster'

# shap_and_feature_path <- "dev/Data/SHAP_and_FeatureMatrix_Mid-length_AmplDel.rds"
# pred_ampl_path <- "results_regressor_gab/InteractomeINSIDER/pred_ampl.rds"
# pred_del_path <- "results_regressor_gab/InteractomeINSIDER/pred_del.rds"

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

bool.plot <- F

output <- list() # define output list

model_types <- c("ampl", "del") # define model types

# for (i in model_types) { # iterate over model types

i <- model_types[1]

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

# cols.to.discard <- c("Chromosome_Length","Centromere_Length","Centromere_Type", # select cols to discard
#                      "BIAS","density.OG","density.TSG","labels","ESSscore_pancancer",
#                      "HAPLOscore_pancancer","Density.complex.proteins", "Density.Ohnologs")

col.sum.zero <- colnames(shap %>% dplyr::select(-labels))[apply(shap %>% dplyr::select(-labels), 2, sum) == 0]
non.segment.features <- c("Chromosome_Length","Centromere_Length","Centromere_Type")

if(model_class %in% c('Mid-length','Small-scale')){cols.to.discard <- c("labels","BIAS", col.sum.zero, non.segment.features) # include non.segment.features if the model is != Arm- or Chr-level
}else{cols.to.discard <- c("labels","BIAS")}

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

####################### examine residuals ####################### 

if (bool.plot) {
  
  residuals_plots <- do_plot_residuals(res, save = F)

}

####################### 

#######################  define ruling system for annotation ####################### 

# define_annotation_rules(i)
define_annotation_rules_correlation_based(i, correlations_table)

####################### 

#######################  definition of ratios ####################### 

message(paste0("Computing ratios"))

features <- colnames(actual_shap_mat) # take all features
labels <- shap_agg$labels # take labels

selection_to_occurrence_ratio <- apply(X = actual_shap_mat, MARGIN = 1, FUN = compute_ratio, features_selection, features_occurrence)
log2_selection_to_occurrence_ratio <- log2(selection_to_occurrence_ratio) # log2(ratio)

positive_to_negative_ratio <- apply(X = actual_shap_mat, MARGIN = 1, FUN = compute_ratio, features_positive_selection, features_negative_selection)
log2_positive_to_negative_ratio <- log2(positive_to_negative_ratio) # log2(ratio)

res_ratio <- as.data.frame(cbind(selection_to_occurrence_ratio, log2_selection_to_occurrence_ratio, 
                                 positive_to_negative_ratio, log2_positive_to_negative_ratio))


res_ratio <- as.data.frame(apply(X = res_ratio, MARGIN = 2, as.numeric))

res_ratio$labels <- labels # add labels
                           # Labels are in the same order, because ratio derive from actual_shap_mat, which derives from shap_agg

merged_res_ratio <-  merge(res, res_ratio, on = "labels") # merge with residuals

merged_res_ratio$type <- do.call(rbind,str_split(merged_res_ratio$labels,'-'))[,2] # take tumor type

merged_res_ratio$is_selection <- ifelse(merged_res_ratio$log2_selection_to_occurrence_ratio > 0, TRUE, FALSE) 
# binarize selection based on log2(ratio)
# if log2(ratio) > 0, selection outweights occurrence
# if log2(ratio) < 0, occurrence outweights selection

merged_res_ratio$is_positive <- ifelse(((merged_res_ratio$log2_selection_to_occurrence_ratio > 0) & 
                                          (merged_res_ratio$log2_positive_to_negative_ratio > 0)), TRUE, 
                                       ifelse(((merged_res_ratio$log2_selection_to_occurrence_ratio > 0) & 
                                                 (merged_res_ratio$log2_positive_to_negative_ratio < 0)), FALSE, 
                                              "NO_SELECTION"))

# binarize positive selection based on log2(ratio) and presence of selection
# if log2(ratio) > 0, positive selection outweights negative selection 
# if log2(ratio) < 0, negative selection outweights positive selection

####################### 

#######################  ratio based annotation ####################### 

message(paste0("Annotating based on ratios"))

merged_res_ratio$annot <- ifelse(merged_res_ratio$is_selection == TRUE & merged_res_ratio$is_positive == TRUE, 'Positive Selection', 
                                 ifelse(merged_res_ratio$is_selection == TRUE & merged_res_ratio$is_positive == FALSE,  'Negative Selection', 
                                        'Occurrence'))

# annotation based on 3 classes, meaning that incorrect prediction flow into postive and negative selection
merged_res_ratio$annot_3_classes <- ifelse(merged_res_ratio$pred_is_correct,
                                           merged_res_ratio$annot, 
                                           ifelse(merged_res_ratio$observed > merged_res_ratio$prediction, 
                                           'Positive Selection',
                                           'Negative Selection'))

# annotation based on 5 classes, meaning that incorrect prediction do not flow into postive and negative selection
merged_res_ratio$annot_5_classes <- ifelse(merged_res_ratio$pred_is_correct, 
                                           merged_res_ratio$annot, 
                                           ifelse(merged_res_ratio$observed > merged_res_ratio$prediction, 
                                                  'Incorrect prediction - Positive Selection',
                                                  'Incorrect prediction - Negative Selection'))

merged_res_ratio$annot <- NULL # remove because it's not the official one
merged_res_ratio$is_selection <- NULL
merged_res_ratio$is_positive <- NULL
merged_res_ratio$pred_is_correct <- NULL

table(merged_res_ratio$annot_3_classes, useNA = "ifany")
table(merged_res_ratio$annot_5_classes, useNA = "ifany")

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
                                       # ifelse(merged_res_annot$top1 %in% features_negative_selection, FALSE, 
                                       #        "NO_SELECTION"))

merged_res_annot$is_negative <- ifelse(merged_res_annot$top1 %in% features_negative_selection, TRUE, FALSE)

# merged_res_annot$is_relaxed_selection <- ifelse(merged_res_annot$is_selection == TRUE & 
#                                                   merged_res_annot$is_positive == "NO_SELECTION", TRUE, FALSE)

merged_res_annot$is_relaxed_positive_selection <- ifelse(merged_res_annot$top1 %in% features_relaxed_positive_selection, TRUE, FALSE)
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

# merged_res_annot$annot <- ifelse(merged_res_annot$is_selection == TRUE & merged_res_annot$is_positive == TRUE, 'Positive Selection', 
#                                  ifelse(merged_res_annot$is_selection == TRUE & merged_res_annot$is_positive == FALSE,  'Negative Selection',
#                                         ifelse(merged_res_annot$is_relaxed_selection == TRUE,  'Relaxed Selection',
#                                         'Occurrence')))

table(merged_res_annot$annot, useNA = "ifany")

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
merged_res_annot$pred_is_correct <- NULL

#######################

#######################  explore annotation results ####################### 

message(paste0("Plotting annotations"))

merged_res_annot$annot_5_classes <- ifelse(merged_res_annot$annot_5_classes == 'Relaxed Negative Selection' | merged_res_annot$annot_5_classes == 'Relaxed Positive Selection', 
                                           'Relaxed Selection', merged_res_annot$annot_5_classes)

save(merged_res_annot, file = paste0('/home/ieo7429/Scrivania/dev/Data/merged_res_annot_',model_class,'_',i,'.RData'))
save(merged_res_ratio, file = paste0('/home/ieo7429/Scrivania/dev/Data/merged_res_ratio_',model_class,'_',i,'.RData'))

if (bool.plot) {
  
  pancancer_plots_ratio <- do_pancancer_barplots_and_boxplots(input_df = merged_res_ratio, save = F)
  tumor_specific_plots_ratio <- do_tumor_specific_barplots_and_boxplots(input_df = merged_res_ratio, what.to.plot = c("3","5"), save = F)
  
  pancancer_plots_annot <- do_pancancer_barplots_and_boxplots(input_df = merged_res_annot, save = F)
  tumor_specific_plots_annot <- do_tumor_specific_barplots_and_boxplots(input_df = merged_res_annot, what.to.plot = c("3","5", "detailed"), save = F)
  
} 

####################### 

####################### explore top1 annotation #######################

message(paste0("Compare the two annotations"))

res_ratio_with_annot <- merge(x = merged_res_ratio, y = merged_res_annot, by = "labels", suffixes = c("_ratio", "_annot"))
res_ratio_with_annot <- res_ratio_with_annot[!(duplicated(as.list(res_ratio_with_annot)))]

colnames(res_ratio_with_annot) <- c("labels", "observed", "prediction", "residual", "pval.residual", "selection_to_occurrence_ratio",    
                                    "log2_selection_to_occurrence_ratio", "positive_to_negative_ratio", "log2_positive_to_negative_ratio",  
                                    "type", "annot_3_classes_ratio", "annot_5_classes_ratio", "top1", "shap_top1", "annot_3_classes_annot", "annot_5_classes_annot")

res_ratio_with_annot <- res_ratio_with_annot[, c(1, 10, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16)]

cont_table_long <- as.data.frame(table(res_ratio_with_annot$annot_5_classes_ratio, res_ratio_with_annot$annot_5_classes_annot))
cont_table_wide <- as.data.frame(pivot_wider(data = cont_table_long, names_from = Var2, values_from = Freq))
rownames(cont_table_wide) <- cont_table_wide$Var1
cont_table_wide$Var1 <- NULL
cont_table_wide <- cont_table_wide[-c(1,2),-c(1,2,6,7)]
chisq_res <- chisq.test(cont_table_wide)
print(chisq_res)

if (bool.plot) {
  
  compare_ratio_based_top1(input_df = res_ratio_with_annot, save = F)
  
}

####################### 

####################### tissue specificity ####################### 

message(paste0("tissue specific"))

res_ratio_with_annot$chr <- sapply(X = strsplit(x = res_ratio_with_annot$labels, split = "[_-]"), FUN = `[`, 1)

tot_type_bins <- res_ratio_with_annot %>% group_by(type) %>% summarise(count_bins = n())

tot_chrom_bins <- res_ratio_with_annot %>% group_by(chr) %>% summarise(count_bins = n())
tot_chrom_bins$chr <- as.integer(tot_chrom_bins$chr)
tot_chrom_bins <- tot_chrom_bins[order(tot_chrom_bins$chr, decreasing = F),]

tot_chrom_bins_per_type <- res_ratio_with_annot %>% group_by(type, chr) %>% summarise(count_bins = n())
tot_chrom_bins_per_type$chr <- as.integer(tot_chrom_bins_per_type$chr)
tot_chrom_bins_per_type <- tot_chrom_bins_per_type[order(tot_chrom_bins_per_type$type, tot_chrom_bins_per_type$chr, decreasing = F),]

distinct_colors <- c(
  "#DC6D56", "#EAA051", "#E3D85B", "#B5E067", "#6FDD8E", "#51D4C3", "#59C5F9", "#77A7F7", "#A48EEC", "#D583E4",
  "#ED82B7", "#EC9B9B", "#D97F66", "#DEA863", "#D9CF63", "#ACD26A", "#6AD48E", "#52CFC2", "#59BFF0", "#709DEB",
  "#9C85E0", "#C87BD6", "#E97AB4", "#E78F96", "#D66C56", "#DB975B", "#D3C45E", "#9FCA68", "#63CD8B", "#4BC7BD",
  "#53B7E6", "#6993E0", "#8C7BD4", "#B572C9", "#E071AF", "#DF857F", "#CC634F", "#D08B58", "#C7B55B", "#92BD63",
  "#57C080", "#41B9B3", "#49A9DB", "#5F89D4", "#7A72C8", "#9F69BE", "#D466A9", "#D27A6C", "#BE5946", "#C38052",
  "#BAA453", "#84AC5D", "#4BAF77", "#37A89E", "#3E98C9", "#537CC1", "#6D65B5", "#905CAC", "#CB5A9F", "#C06D61",
  "#AF4F3D", "#B47447", "#AB9848", "#769F51", "#409F6D", "#2A9993", "#3087BC", "#456CAE", "#5D579F", "#7E4D95",
  "#C05493", "#A95D54", "#9B4231", "#9F663A", "#967A3A", "#6A843F", "#36865C", "#217F81", "#256C9E", "#3B548F",
  "#52407D", "#733871", "#B94A86", "#964C3E", "#863224", "#8B4F2E", "#816131", "#5B682F", "#2D6A4B", "#17646E",
  "#1A537D", "#2E4070", "#443160", "#622A56", "#A64079", "#803F2A", "#712414", "#763E22", "#6C4E24", "#4B531F",
  "#1E553A", "#094F5B", "#0C3F6A", "#1E2E5C", "#311F4C", "#4C1A42"
)

cols.to.use <- c("top1", "annot_5_classes_annot", "annot_3_classes_annot",
                         "annot_5_classes_ratio", "annot_3_classes_ratio")

for (col in cols.to.use) {
  
grouped_by_type <- do_plot_grouped_by_type(input_df = res_ratio_with_annot, 
                        tot_type_bins = tot_type_bins, 
                        col.to.group = col, 
                        distinct_colors = distinct_colors, 
                        save = FALSE)

grouped_by_chr <- do_plot_grouped_by_chr(input_df = res_ratio_with_annot, 
                        tot_chrom_bins = tot_chrom_bins, 
                        col.to.group = col, 
                        distinct_colors = distinct_colors, 
                        save = FALSE)

grouped_by_chr_type <- do_plot_grouped_by_chr_type(input_df = res_ratio_with_annot, 
                            tot_chrom_bins_per_type = tot_chrom_bins_per_type, 
                            col.to.group = col, 
                            distinct_colors = distinct_colors, 
                            save = FALSE)

}

grouped_by_chr_type <- do_plot_grouped_by_chr_type(input_df = res_ratio_with_annot, 
                                                   tot_chrom_bins_per_type = tot_chrom_bins_per_type, 
                                                   col.to.group = "annot_5_classes_ratio", 
                                                   distinct_colors = distinct_colors, 
                                                   save = FALSE)

grouped_by_type <- do_plot_grouped_by_type(input_df = res_ratio_with_annot, 
                                           tot_type_bins = tot_type_bins, 
                                           col.to.group = 'annot_5_classes_annot', 
                                           distinct_colors = distinct_colors, 
                                           save = FALSE)

pdf(paste0('plots_annot_final/',model_class,'_',i,'_barplot_grouped_by_type.pdf'), width = 10)
grouped_by_type$plot
dev.off()


output_table <- grouped_by_chr_type$data
output_table$count_bins <- NULL

write.table(x = output_table, file = paste0("../Data/", i, "_data_grouped_by_chr_tt.tsv"), quote = F, sep = "\t", row.names = F, col.names = T)

#######################

####################### fabio's correlation plots ####################### 

good.segments <- merged_res_ratio %>% filter(pval.residual >= pval_thr)

actual_shap_mat <- cbind(actual_shap_mat, labels = merged_res_ratio$labels) %>% as_tibble() %>% filter(labels %in% good.segments$labels) %>% dplyr::select(-labels)
actual_values_mat <- cbind(actual_values_mat, labels = merged_res_ratio$labels) %>% as_tibble() %>% filter(labels %in% good.segments$labels) %>% dplyr::select(-labels)

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
cor_df$feature_type <- ifelse(cor_df$Column %in% features_occurrence, "Occurrence", "Selection")
cor_df$signif <- ifelse(cor_df$pval.S < 0.001, "***",
                        ifelse(cor_df$pval.S < 0.01, "**",
                               ifelse(cor_df$pval.S < 0.05, "*", "")))

cor_df$pval.P_fdr <- p.adjust(cor_df$pval.P, method = 'bonferroni') # ho provato a vedere se facendo multiple testing correction perdo significatività su OG, ma purtroppo no! :(

#cor_df$is_signif_both <- ifelse(cor_df$estimate.P > 0, 'pos', 'neg')

pS <- ggplot(cor_df, aes(x = Column, y = estimate.S)) +
  geom_col(aes(fill = feature_type)) +
  geom_text(aes(
    label = paste0(round(estimate.S, 2), signif)
  )) +
  theme_minimal() +
  labs(
    title = "Correlation Between SHAP and Feature Matrix",
    y = "Spearman Correlation"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


cor_df$Column <- factor(cor_df$Column, levels = cor_df[order(cor_df$estimate.P, decreasing = T),'Column'])

pP <- ggplot(cor_df, aes(x = Column, y = estimate.P)) +
  geom_col(aes(fill = feature_type)) +
  geom_text(aes(label = paste0(round(estimate.P, 2), signif), vjust = -0.5)) +
  theme_minimal() +
  labs(title = "Correlation Between SHAP and Feature Matrix",
       y = "Pearson Correlation") +
  theme(axis.text = element_text(angle = 90))

ggsave(
  filename = paste0(model_type,'::',ifelse(i=='ampl','AmplificationModel','DeletionModel'), "_pancancer_correlations_SHAP_actualValues_Pearson.pdf"),
  plot = pP, device = "pdf",
  width = 15, height = 10
)
ggsave(
  filename = paste0(ifelse(i=='ampl','AmplificationModel','DeletionModel'), "_pancancer_correlations_SHAP_actualValues_Spearman.pdf"),
  plot = pS, device = "pdf",
  width = 15, height = 10
)

####################### 

####################### landscape #######################

res_ratio_with_annot$binID <- unlist(lapply(res_ratio_with_annot$labels, function(x){strsplit(x, "-")[[1]][1]}))
# res_ratio_with_annot <- res_ratio_with_annot[, c(1,18,
#                                                  17,2,3:16)]

backbone.100kb <- do.call(rbind, chr_backbone_namesfixed$`0.1Mbp`)
backbone.100kb$binID <- paste0(backbone.100kb$chr, "_", backbone.100kb$bin)

backbone.100kb <- backbone.100kb %>% 
                    dplyr::select(binID = binID, 
                    start = start_bin, 
                    end = end_bin,
                    chr = chr)

res_ratio_with_annot_with_backbone <- merge(x = res_ratio_with_annot, y = backbone.100kb, by = c("binID","chr"))
# res_ratio_with_annot_with_backbone <- res_ratio_with_annot_with_backbone[, c(1,2,4,3,19,20,5:18)]

res_ratio_with_annot_with_backbone <- res_ratio_with_annot_with_backbone[order(as.integer(res_ratio_with_annot_with_backbone$chr),
                                                                               as.integer(res_ratio_with_annot_with_backbone$start),
                                                                               decreasing = FALSE),]

res_ratio_with_annot_with_backbone$chr <- paste0("chr",res_ratio_with_annot_with_backbone$chr)
res_ratio_with_annot_with_backbone <- res_ratio_with_annot_with_backbone[,c(2,20,18,19,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)]

write.table(x = res_ratio_with_annot_with_backbone, 
            file = paste0("/home/ieo7429/Scrivania/dev/Data/res_ratio_with_annot_with_backbone_",
                                                                  i, "_", model_class, ".tsv"),
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)


# tt <- 'BRCA'
for(tt in levels(factor(res_ratio_with_annot_with_backbone$type))){
  print(tt)
  p <- res_ratio_with_annot_with_backbone %>% dplyr::filter(type == tt) %>%
    landscape_plot_observed_prediction(c(paste0('chr',1:22)), 
                                       type_mask = tt, 
                                       model_mask = i, # 'ampl' or 'del'
                                       plot_observed = TRUE, 
                                       plot_predicted = TRUE,
                                       annot_to_plot = 'annot_5_classes_annot', # 'annot_5_classes_ratio' or 'annot_5_classes_annot',
                                       annot_to_plot_ticks = "all",
                                       annot_to_plot_kde = "all",
                                       annot_to_plot_meta = "log2_selection_to_occurrence_ratio", # 'log2_selection_to_occurrence_ratio' or 'log2_positive_to_negative_ratio'
                                       make.interactive = F)
  pdf(file = paste0('../imgs/plots_annot_final/landscape_plots/',i,'_landscape_',tt,'.pdf'), width = 25, height = 15)
  print(p)
  dev.off()
}

#######################



