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
source("scripts/landscape_plot_observed_only.R")
source("scripts/do_plot_residuals.R")

####################### 

#######################  read model's output ####################### 

model_classes <- c('Arm-level') 

#for (model_class in model_classes) {

model_class <- model_classes[1]

shap_and_feature_path <- paste0("dev/Data/SHAP_and_FeatureMatrix_",model_class,"_AmplDel.rds")
pred_ampl_path <- paste0("results_regressor_gab/InteractomeINSIDER/",model_class,"-pred_ampl.rds")
pred_del_path <- paste0("results_regressor_gab/InteractomeINSIDER/",model_class,"-pred_del.rds")
correlations_table_path <- paste0("results_regressor_gab/final_verdict.tsv")

df <- readRDS(file = shap_and_feature_path)
pred_ampl <- readRDS(file = pred_ampl_path) 
pred_del <- readRDS(file = pred_del_path)
correlations_table <- read.table(file = correlations_table_path, header = T, sep = "\t")

####################### 

#######################  manipulate data and define variables ####################### 

res.ampl <- get_residuals(pred_ampl); res.ampl$labels <- paste0(res.ampl$bin, "-", res.ampl$Type) # load ampl residuals
res.del <- get_residuals(pred_del); res.del$labels <- paste0(res.del$bin, "-", res.del$Type)      # load del residuals

bool.plot <- F

output <- list() # define output list

model_types <- c("ampl", "del") # define model types

#for (i in model_types) { # iterate over model types

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

mutations_norm <- values[,c("mutations_norm", "labels")]

annot_df_list <- c(ampl = "dev/Data/merged_res_annot_Arm-level_ampl_0.6_0.0205180818215013_0.002.tsv",
                   del = "dev/Data/merged_res_annot_Arm-level_del_0.6_0.0201767960563302_0.01.tsv")

annot_df <- read.table(file = annot_df_list[i], header = T, sep = "\t")
annot_df <- annot_df[,c("labels","observed", "prediction", "top1", "shap_top1", "annot_final")]
annot_df_with_mut <- merge(mutations_norm, annot_df, by = "labels")

annot_df_with_mut$shap_top1 <- ifelse(annot_df_with_mut$annot_final %in% c("Positive Selection","Negative Selection"), annot_df_with_mut$shap_top1, 
                                      ifelse(annot_df_with_mut$annot_final %in% c("Incorrect prediction - Positive Selection","Incorrect prediction - Negative Selection"), -0.01, 0))
annot_df_with_mut$annot_final <- ifelse(annot_df_with_mut$annot_final %in% c("Positive Selection","Negative Selection"), annot_df_with_mut$top1, annot_df_with_mut$annot_final)

position_and_tt <- as.data.frame(do.call(rbind, lapply(X = annot_df_with_mut$labels, FUN = function(x){
  parts <- strsplit(x, split="-")[[1]]
  binID <- parts[[1]]
  
  chr <- strsplit(parts, split = "_")[[1]][[1]]
  prog <- strsplit(parts, split = "_")[[1]][[2]]
  tt <- parts[[2]]
  
  return(c(x, chr, prog, tt))
  
  })))


colnames(position_and_tt) <- c("labels", "chr", "prog", "tt")

position_and_tt$chr <- as.integer(position_and_tt$chr)
position_and_tt$prog <- as.integer(position_and_tt$prog)

plotting_df <- merge(annot_df_with_mut, position_and_tt, by = "labels")
plotting_df <- plotting_df[order(plotting_df$chr, plotting_df$prog, decreasing = F),]
plotting_df <- plotting_df %>% filter(tt == "BRCA") %>% mutate(pos = row_number())
plotting_df$chr <- paste0("chr", plotting_df$chr)

plotting_df <- plotting_df %>%
  mutate(annot_final = recode(annot_final,
                              `-dist.to.closest.TSG` = "Proximity to closest TSG",
                              `dist.to.closest.TSG` = "Proximity to closest TSG",
                              `dist.to.closest.OG` = "Proximity to closest OG",
                              `-dist.to.closest.OG` = "Proximity to closest OG",
                              `Incorrect prediction - Positive Selection` = "Positive Selection - Excess of Observation",
                              `Incorrect prediction - Negative Selection` = "Negative Selection - Excess of Prediction"))

plotting_df$binID <- paste0(plotting_df$chr, "_", plotting_df$prog)

chr_bounds <- plotting_df %>%
  mutate(chr = as.character(chr)) %>%
  group_by(chr) %>%
  summarise(start = min(pos), end = max(pos), .groups = "drop") %>%
  mutate(chr_num = readr::parse_number(chr)) %>%
  arrange(chr_num) %>%
  mutate(fill = rep(c("white", "#e7deed"), length.out = n()))


library(ggplot2)
library(dplyr)
library(ggnewscale)  # for new_scale_fill()

library(ggnewscale)  # for new_scale_fill()

# Define your palette: 8 distinct colors + grey for "No Detectable Force"
custom_palette <- c(
  "Proximity to closest TSG" = "#D95F02",               # muted orange-red
  "Proximity to closest OG" = "#7570B3",                # dusty blue
  "all.int.trans" = "#1B9E77",                          # calmer green
  "genes.bin" = "#A6761D",                              # earthy mustard
  "mutations_norm" = "#00FF00",                         # ⚡ flashy neon magenta
  "Negative Selection - Excess of Prediction" = "#666666",  # dark grey-brown
  "Positive Selection - Excess of Observation" = "#E78AC3"  # pastel pink
)

shap_shift <- 0  # You can change this value to shift more or less

p <- ggplot() +
  geom_rect(
    data = chr_bounds,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = fill),
    alpha = 0.3
  ) +
  scale_fill_identity() +
  new_scale_fill() +
  geom_bar(
    data = plotting_df,
    aes(x = pos, y = shap_top1 - shap_shift, fill = annot_final),  # SHIFT APPLIED HERE
    stat = "identity"
  ) +
  scale_fill_manual(values = custom_palette) +
  
  geom_hline(yintercept = 0 - shap_shift, linetype = "dashed", color = "grey", linewidth = 0.2) +  # ADJUSTED Y-INTERCEPT
  
  scale_x_continuous(
    breaks = chr_bounds %>% mutate(center = (start + end)/2) %>% pull(center),
    labels = chr_bounds$chr,
    expand = c(0.005, 0)
  ) +
  
  labs(
    title = "Segment Annotation (based on SHAP values)",
    subtitle = "[WHOLE GENOME] [BRCA]",
    x = "Genomic Position",
    y = "Contribution to SCNA frequency (Mid-length)"
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.line.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  )

landscape_plot <- landscape_plot_observed_prediction(
  filtered_landscape = plotting_df,
  genome_mask = unique(plotting_df$chr),
  type_mask = "BRCA",
  model_mask = i  # "ampl" or "del"
)

# Add your custom layers to the landscape_plot
combined_plot <- landscape_plot +
  # Add chromosome background
  geom_rect(
    data = chr_bounds,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = fill),
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  scale_fill_identity() +
  new_scale_fill() +
  
  # Add SHAP bar plot layer on top - SHIFT APPLIED HERE
  geom_bar(
    data = plotting_df,
    aes(x = pos, y = shap_top1 - shap_shift, fill = annot_final),  # SHIFT APPLIED HERE
    stat = "identity",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = custom_palette) +
  
  # Optional aesthetics
  labs(
    title = "Segment Annotation (Observed + SHAP)",
    subtitle = "[WHOLE GENOME] [BRCA]",
    x = "Genomic Position",
    y = "Contribution to SCNA frequency (Mid-length)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.line.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave(
  filename = paste0("arm_level_", i, "_BRCA_shap.pdf"),
  plot = combined_plot,
  device = "pdf",
  width = 12,
  height = 8
)






