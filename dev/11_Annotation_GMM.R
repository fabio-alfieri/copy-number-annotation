rm(list=ls())
gc(full=T)

# SHAP value annotation
library(ggplot2)
library(tidyr)
library(stringr)
library(tidyverse)

# load("~/Desktop/fabio/results_regressor_unscaled/SHAP/Avg_shap_values_InteractomeINSIDER.RData")
# 
# df <- list(models.X.test = list(`Mid-length::Amplification model` = models.X.test$`Mid-length::Amplification model`,
#                                `Mid-length::Deletion model` = models.X.test$`Mid-length::Deletion model`),
#           models.shap.df = list(`Mid-length::Amplification model`= models.shap.df$`Mid-length::Amplification model`,
#                                 `Mid-length::Deletion model` = models.shap.df$`Mid-length::Deletion model`))
# write_rds(df, file = 'dev/Data/SHAP_and_FeatureMatrix_Mid-length_AmplDel.rds')
df <- readRDS(file = 'dev/Data/SHAP_and_FeatureMatrix_Mid-length_AmplDel.rds')

output <- list()
for(i in c('ampl','del')){
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
  shap <- change_distance(shap)
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
  
  # K-means clustering on SHAP aggregated values
  # set.seed(123)
  # k <- 50
  # km <- kmeans(shap_agg %>% select(-labels), centers = k)
  # shap_agg$cluster <- km$cluster
  
  library(mclust)
  set.seed(1234)
  gmm <- Mclust(shap_agg %>% select(-labels), G = 25)
  
  shap_agg$cluster <- predict(gmm)$classification
  table(shap_agg$cluster)
  
  # Merge cluster info into values_agg
  values_agg <- left_join(shap_agg %>% select('labels','cluster'), 
                          values_agg,
                          by = 'labels')
  
  values_agg[,-c(1:2)] <- apply(values_agg[,-c(1:2)], 2, scale)
  
  # Summarize value contributions per cluster
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  value_long <- pivot_longer(values_agg, cols = -c(labels, cluster), names_to = "feature", values_to = "value")
  
  summary_df <- value_long %>%
    group_by(cluster, feature) %>%
    summarize(mean_value = median(value, na.rm = TRUE)) %>%
    ungroup()
  
  # Plot feature contribution by cluster
  library(scico)
  ggplot(summary_df, aes(x = factor(cluster), y = mean_value, fill = feature)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Cluster", y = "Mean Feature Value", title = "Region Clustering by SHAP: Feature Contributions") +
    theme_minimal() +
    scale_fill_viridis_d(option = "H")
  # scale_fill_scico_d(palette = "vik")
  
  top_features <- summary_df %>%
    group_by(cluster) %>%
    top_n(2, wt = mean_value) %>%
    arrange(cluster, desc(mean_value)) %>%
    summarise(label = paste(feature, collapse = ", "))
  # View(top_features)
  
  aggregated <- top_features %>%
    group_by(label) %>%
    summarise(
      n_clusters = n(),
      clusters = paste(sort(unique(cluster)), collapse = ", ")
    ) %>%
    arrange(desc(n_clusters))
  # View(aggregated)
  
  if(F){
    library(ggplot2)
    library(Rtsne)
    # Scale the data for PCA/t-SNE
    shap_scaled <- shap_agg %>% select(-labels) %>% scale()
    cluster_labels <- shap_agg$cluster
    # PCA
    pca_res <- prcomp(shap_scaled, center = TRUE, scale. = TRUE)
    pca_df <- as.data.frame(pca_res$x[, 1:2])
    pca_df$cluster <- as.factor(cluster_labels)
    
    ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
      geom_point(alpha = 0.7) +
      labs(title = "PCA of Regions Based on SHAP Values") +
      theme_minimal()
  }
  
  # try to plot regions
  toplot.plot <- full_join(values_agg %>% select(labels,cluster),values %>% select(bin,,Type,labels,ampl_score,del_score)) %>% select(-labels) %>%
    filter(Type == 'BRCA') %>% separate(bin, sep = '_', into = c('chr','bin'))
  toplot.plot <- toplot.plot %>%
    arrange(as.numeric(as.character(chr)), as.numeric(bin)) %>%
    mutate(pos = row_number())
  
  
  aggregated$clusters.aggregated <- rownames(aggregated)
  cluster_mapping <- aggregated %>%
    separate_rows(clusters, sep = ",\\s*") %>%  # splits by comma and optional space
    mutate(clusters = as.integer(clusters)) %>%
    select(original_cluster = clusters, clusters.aggregated)
  toplot.plot <- toplot.plot %>%
    left_join(cluster_mapping, by = c("cluster" = "original_cluster")) %>%
    mutate(cluster = clusters.aggregated) %>%
    select(-clusters.aggregated)
  toplot.plot$cluster <- as.numeric(toplot.plot$cluster)
  
  # Compute chromosome ranges for background shading
  chr_bounds <- toplot.plot %>%
    group_by(chr) %>%
    summarize(start = min(pos), end = max(pos)) %>%
    mutate(fill = rep(c("lightgrey", "darkgrey"), length.out = n()))
  
  base_plot <- ggplot() +
    geom_rect(data = chr_bounds,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.3) +
    scale_fill_manual(values = c("#eeeeee", "#cccccc")) +
    geom_line(data = toplot.plot, aes(x = pos, y = ampl_score), color = "red") +
    geom_line(data = toplot.plot, aes(x = pos, y = -del_score), color = "blue") +
    geom_hline(yintercept = 0, colour = 'grey', linetype = 'dashed') +
    labs(x = "Genomic Position", y = "Amplification Score (Mid-length)") +
    theme_classic()
  
  # Add cluster lines above plot
  cluster_ticks <- toplot.plot %>%
    mutate(cluster_ymin = max(toplot.plot$ampl_score) * 1.07 + cluster * 0.015,
           cluster_ymax = cluster_ymin + 0.01)
  
  p.final <- base_plot +
    geom_segment(data = cluster_ticks,
                 aes(x = pos, xend = pos, y = cluster_ymin, yend = cluster_ymax, col = as.character(cluster)),
                 size = 0.3) +
    scale_colour_viridis_d(option = "H")
  
  output[[i]]$p.final <- p.final
  output[[i]]$toplot <- toplot.plot
  output[[i]]$aggregated <- aggregated
  output[[i]]$top_features <- top_features
}

# Visualize plots with annotatiion for ampl or del
output$ampl$p.final
output$del$p.final

# Visualize aggregated fatures
output$ampl$aggregated
output$del$aggregated




