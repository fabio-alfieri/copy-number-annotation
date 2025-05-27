rm(list=ls())
gc(full=T)

library(ggplot2)
library(tidyr)
library(stringr)
library(tidyverse)

df <- readRDS(file = 'dev/Data/SHAP_and_FeatureMatrix_Mid-length_AmplDel.rds')

source('dev/residuals_plot.R')
res.ampl <- plot_residuals(read_rds('dev/Data/pred_ampl.rds'))
res.del <- plot_residuals(read_rds('dev/Data/pred_del.rds'))

res.filt <- list()
quantile_filt <- 0.95
res.filt[['ampl']] <- res.ampl %>% filter(residual <= as.numeric(quantile(res.ampl$residual, prob = quantile_filt)))
res.filt[['del']] <- res.del %>% filter(residual <= as.numeric(quantile(res.del$residual, prob = quantile_filt)))

tt <- 'BRCA'

output <- list()

i <- "ampl"

midlength <- list()
if(i == 'ampl'){
  midlength[['shap']] <- df$models.shap.df$`Mid-length::Amplification model`
  midlength[['values']] <- df$models.X.test$`Mid-length::Amplification model`
}else if(i == 'del'){
  midlength[['shap']] <- df$models.shap.df$`Mid-length::Deletion model`
  midlength[['values']] <- df$models.X.test$`Mid-length::Deletion model`
}

# Extract matrices
shap <- midlength$shap
values <- midlength$values
values$labels <- paste0(values$bin,'-',values$Type)

values <- left_join(res.filt[[i]] %>% select(bin, Type, residual),
                    values)
shap <- left_join(values %>% select(labels), shap)

# Aggregate Chromatin state into categories in both shap and values
aggregate_chromatin_states <- function(df) {
  df$promoters <- df$Length_Counts.E1 + df$Length_Counts.E2 + df$Length_Counts.E3 +
    df$Length_Counts.E4 + df$Length_Counts.E22 + df$Length_Counts.E23
  df$transcribed <- df$Length_Counts.E5 + df$Length_Counts.E6 + df$Length_Counts.E7 +
    df$Length_Counts.E8 + df$Length_Counts.E9 + df$Length_Counts.E10 +
    df$Length_Counts.E11 + df$Length_Counts.E12
  df$enhancers <- df$Length_Counts.E13 + df$Length_Counts.E14 + df$Length_Counts.E15 +
    df$Length_Counts.E16 + df$Length_Counts.E17 + df$Length_Counts.E18
  df$repressed <- df$Length_Counts.E19 + df$Length_Counts.E20 + df$Length_Counts.E24 +
    df$Length_Counts.E25
  df$accessible <- df$Length_Counts.E21
  return(df)
}
change_distance <- function(df) {
  distances <- c("dist.to.closest.OG", "dist.to.closest.TSG", "dist.to.closest.FGS",
                 "distance.to.centromere", "distance.to.telomere", "Ess.distance_pancancer")
  df[,distances] <- -1*df[,distances]
  return(df)
}

shap <- aggregate_chromatin_states(shap)
values <- aggregate_chromatin_states(values)
# shap <- change_distance(shap)
values <- change_distance(values)

# Subset to relevant columns
columns <- c("dist.to.closest.OG", "dist.to.closest.TSG", "dist.to.closest.FGS",
             "distance.to.centromere", "distance.to.telomere", "mutations_norm",
             "genes.bin", "promoters", "transcribed", "Ess.distance_pancancer",
             "enhancers", "repressed","accessible",
             # "mean.GC.content","total_n_partners.trans","total_n_PPIs.trans",
             # "total_n_ohnologs.mmpaper_trans", "total_n_paralogs_trans",
             # "all.int.trans","partners.trans",
             # "HAPLOscore_pancancer","Density.complex.proteins","Density.Ohnologs",
             'labels')
shap_agg <- shap[, columns]
values_agg <- values[, columns]


num_of_top_features <- 5
actual_values_mat <- values_agg[,-ncol(values_agg)]
actual_shap_mat <- shap_agg[,-ncol(shap_agg)] 

res <- t(apply(X = actual_shap_mat, 
      MARGIN = 1, 
      FUN = function(row){
        idxs <- order(unlist(abs(row)), na.last = TRUE)[1:num_of_top_features]
        cluster_identity <- paste0(idxs, collapse = "")
        top_features <- colnames(actual_shap_mat)[idxs]
        return(c(cluster_identity, top_features))
      }))

# unique_elems <- unique(as.vector(res))
# palette <- rainbow(n = length(unique_elems))
# ComplexHeatmap::Heatmap(res[1:100,], col = palette)

actual_shap_mat$cluster_id <- res[,1]

actual_shap_mat_medians <- actual_shap_mat %>%
  group_by(cluster_id) %>%
  summarise(across(where(is.numeric), \(x) median(x, na.rm = TRUE)))

# ComplexHeatmap::Heatmap(cor(t(actual_shap_mat_medians[1:100,-1])))

cluster_id_first_round <- actual_shap_mat_medians$cluster_id
actual_shap_mat_medians$cluster_id <- NULL
rownames(actual_shap_mat_medians) <- cluster_id_first_round

dist_matrix <- dist(actual_shap_mat_medians, method = "euclidean")
hc <- hclust(dist_matrix, method = "complete")

num_of_clusters = 16
new_cluster_ids <- cutree(hc, k = num_of_clusters)
actual_shap_mat_medians$new_cluster_id <- new_cluster_ids[rownames(actual_shap_mat_medians)]
actual_shap_mat_medians$cluster_id <- cluster_id_first_round

# dend <- as.dendrogram(hc)
# height_thr <- 0.05
# collapsed_dend <- cut(dend, h = height_thr)$upper
# plot(collapsed_dend, main = paste("Dendrogram Cut at Height", height_thr))

# per_cluster_purity <- actual_shap_mat_medians %>%
#  group_by(new_cluster_id) %>%
#  summarise(
#    total = n(),
#    most_common_label_count = max(table(cluster_id)),
#    purity = most_common_label_count / total
#  )

actual_shap_mat <- actual_shap_mat %>%
  left_join(
    actual_shap_mat_medians %>% select(cluster_id, new_cluster_id),
    by = "cluster_id"
  )

actual_values_mat <- data.frame(
  apply(X = actual_values_mat,
        MARGIN = 2,
        FUN = scale)
  )

actual_values_mat$cluster_id <- actual_shap_mat$cluster_id
actual_values_mat$new_cluster_id <- actual_shap_mat$new_cluster_id

df_long <- actual_values_mat %>%
  select(-cluster_id) %>%
  pivot_longer(
    cols = -new_cluster_id,
    names_to = "feature",
    values_to = "value"
  )

df_medians <- df_long %>%
  group_by(new_cluster_id, feature) %>%
  summarise(median_value = median(value, na.rm = TRUE), .groups = "drop")

plots <- df_medians %>%
  group_split(new_cluster_id) %>%
  lapply(function(sub_df) {
    cluster_label <- unique(sub_df$new_cluster_id)
    
    ggplot(sub_df, aes(x = feature, y = median_value)) +
      geom_col(aes(fill = feature)) +
      labs(title = paste("Cluster", cluster_label),
           x = "Feature", y = "Median Value") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.position = "none")
  })

do.call(gridExtra::grid.arrange, c(plots, ncol = 4))




                            