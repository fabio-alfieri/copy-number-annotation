rm(list=ls())
gc(full=T)

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

load("/home/ieo7429/Scrivania/dev/Data/hMEC_frequencies_day2-day28.RData")
hMEC <- to_save
rm(to_save)
hMEC_day2 <- hMEC$day2$gain[,c("chr", "hMEC_150nM")]; colnames(hMEC_day2) <- c("armID", "hMEC_150nM_day2")
hMEC_day28 <- hMEC$day28$gain[,c("chr", "hMEC_150nM")]; colnames(hMEC_day28) <- c("armID", "hMEC_150nM_day28")
hMEC <- merge(hMEC_day2, hMEC_day28, by = "armID")

load("/home/ieo7429/Scrivania/dev/Data/hPDE_frequencies_day2-day28.RData")
hPDE <- to_save
rm(to_save)
hPDE_day2 <- hPDE$day2$gain[,c("chr", "hPDE_75nM")]; colnames(hPDE_day2) <- c("armID", "hPDE_75nM_day2")
hPDE_day28 <- hPDE$day28$gain[,c("chr", "hPDE_75nM")]; colnames(hPDE_day28) <- c("armID", "hPDE_75nM_day28")
hPDE <- merge(hPDE_day2, hPDE_day28, by = "armID")

match_100kbp_arm <- read.table(file = "/home/ieo7429/Scrivania/dev/Data/match_bins_arms_clean.tsv", header = T)
colnames(match_100kbp_arm) <- c("binID", "armID")

setwd("/home/ieo7429/Scrivania/")

source("scripts/get_residuals.R")

model_classes <- c('Mid-length', 'Arm-level', 'Chromosome-level', 'Small-scale', 'no_cluster')

annots <- c("dev/Data/new_stuff/merged_res_annot_Mid-length_ampl_0.6_0.0632314186543226_0.02.tsv",
            "dev/Data/new_stuff/merged_res_annot_Arm-level_ampl_0.6_0.0205180818215013_0.002.tsv",
            "dev/Data/new_stuff/merged_res_annot_Chromosome-level_ampl_0.6_0.00890466617420316_0.01.tsv",
            "dev/Data/new_stuff/merged_res_annot_Small-scale_ampl_0.6_0.00773112666793168_4e-04.tsv",
            "dev/Data/new_stuff/merged_res_annot_no_cluster_ampl_0.6_0.0974942784756422_0.04.tsv")
names(annots) <- model_classes

for(model_class in model_classes){

#model_class <- model_classes[1]

annot_path <- annots[model_class]

annot <- read.delim(file = annot_path); annot$type <- NULL
binID_and_type <- do.call(rbind, strsplit(annot$labels, split = "-"))
colnames(binID_and_type) <- c("binID", "type")

annot <- cbind(binID_and_type, annot)

annot$annot_final <- ifelse(annot$annot_final == "Incorrect prediction - Positive Selection", "Positive Selection",
                            ifelse(annot$annot_final == "Incorrect prediction - Negative Selection", "Negative Selection", annot$annot_final)
                            )

annot_with_arms <- merge(x = annot, y = match_100kbp_arm, by = "binID")
annot_with_arms <- annot_with_arms[,c(1,2,3,4,5,6,17,19)]

annot_with_arms_BRCA <- annot_with_arms[annot_with_arms$type == "BRCA",c("armID", "binID", "type", "labels", "prediction", "observed", "residual", "annot_final")]
annot_with_arms_PAAD <- annot_with_arms[annot_with_arms$type == "PAAD",c("armID", "binID", "type", "labels", "prediction", "observed", "residual", "annot_final")]


aggregated_BRCA <- data.frame(annot_with_arms_BRCA %>%
  group_by(armID) %>%
  summarise(median_obs = median(observed),
            median_pred = median(prediction),
            freq = list((table(annot_final) / sum((table(annot_final))))
            )) %>%
  unnest_wider(freq, names_sep = "_"))

aggregated_BRCA[is.na(aggregated_BRCA)] <- 0

aggregated_PAAD <- data.frame(annot_with_arms_PAAD %>%
                                group_by(armID) %>%
                                summarise(median_obs = median(observed),
                                          median_pred = median(prediction),
                                          freq = list((table(annot_final) / sum((table(annot_final))))
                                          )) %>%
                                unnest_wider(freq, names_sep = "_"))

aggregated_PAAD[is.na(aggregated_PAAD)] <- 0

aggregated_BRCA$freq_No.Detectable.Force <- NULL
aggregated_BRCA$freq_Positive.Selection <- NULL
aggregated_BRCA$freq_Negative.Selection <- NULL
aggregated_BRCA$freq_Occurrence <- NULL

aggregated_PAAD$freq_No.Detectable.Force <- NULL
aggregated_PAAD$freq_Positive.Selection <- NULL
aggregated_PAAD$freq_Negative.Selection <- NULL
aggregated_PAAD$freq_Occurrence <- NULL



setwd("dev/imgs/plots_annot_final")

library(ggplot2)
library(dplyr)
library(tidyr)

plot_with_spearman <- function(df, 
                               xvar, 
                               yvar, 
                               title=NULL, 
                               labelvar = NULL) {
  
  pearson_test <- cor.test(df[[xvar]], df[[yvar]], method = "pearson", exact = FALSE)
  pearson_r <- pearson_test$estimate
  pearson_p <- pearson_test$p.value
  
  pearson_label <- paste0("pearson R = ", round(pearson_r, 2),
                           ", p = ", signif(pearson_p, 3))
  
  p <- ggplot(df, aes_string(x = xvar, y = yvar)) +
    geom_point(color = "steelblue", size = 2) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    annotate("text", 
             x = min(df[[xvar]], na.rm = TRUE), 
             y = max(df[[yvar]], na.rm = TRUE), 
             label = pearson_label, 
             hjust = 0, vjust = 1, size = 5, color = "black") +
    labs(
      x = xvar,
      y = yvar,
      title = title
    ) +
    theme_minimal(base_size = 14)
  
  if(!is.null(labelvar) && labelvar %in% colnames(df)) {
    p <- p + geom_text(aes_string(label = labelvar), 
                       hjust = 0, vjust = -0.5, size = 3, color = "darkred")
  }
  
  return(p)
}

plot_list <- list(
  
 
  list(df = data.frame(
    BRCA_observed = aggregated_BRCA$median_obs,
    PAAD_observed = aggregated_PAAD$median_obs,
    armID = aggregated_BRCA$armID),
    x = "BRCA_observed", 
    y = "PAAD_observed", 
    title = "Observed median comparison", 
    label = "armID"),
  
  list(df = data.frame(
    BRCA_predicted = aggregated_BRCA$median_pred,
    PAAD_predicted = aggregated_PAAD$median_pred,
    armID = aggregated_BRCA$armID),
    x = "BRCA_predicted", 
    y = "PAAD_predicted", 
    title = "Predicted median comparison", 
    label = "armID")
)

all_plots <- lapply(plot_list, function(p) {
  plot_with_spearman(p$df, 
                     xvar = p$x, 
                     yvar = p$y, 
                     title = p$title, 
                     labelvar = p$label)
})

for (i in seq_along(all_plots)) {
  plt <- all_plots[[i]]
  safe_title <- gsub("[^A-Za-z0-9_]", "_", plot_list[[i]]$title)

  filename <- paste0(safe_title, "_", model_class, ".pdf")
  
  ggsave(filename, plt, width = 8, height = 6, dpi = 300)
  print(paste("Saved plot:", filename))
}

setwd("../../../")

}

df <- merge(hMEC[,c("armID", "hMEC_150nM_day28")], hPDE[,c("armID", "hPDE_75nM_day28")], by = "armID")

p <- plot_with_spearman(df, 
                   xvar = "hMEC_150nM_day28", 
                   yvar = "hPDE_75nM_day28", 
                   title = "Cell Line comparison", 
                   labelvar = "armID")


ggsave("cell_line_plot.pdf", p, width = 8, height = 6, dpi = 300)

