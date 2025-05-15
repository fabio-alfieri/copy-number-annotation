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
  pivot_longer(values_agg, cols = -c('labels','cluster')) %>% 
    ggplot() + geom_histogram(aes(x = value)) + facet_wrap(~name)
  
  # Summarize value contributions per cluster
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  value_long <- pivot_longer(values_agg, cols = -c(labels, cluster), names_to = "feature", values_to = "value")
  
  summary_df <- value_long %>%
    group_by(cluster, feature) %>%
    summarize(mean_value = median(value, na.rm = TRUE)) %>%
    ungroup()
  
  summary_df <- summary_df %>%
    group_by(feature) %>%
    mutate(mean_value = ifelse(is.na(mean_value), min(mean_value, na.rm = TRUE), mean_value)) %>%
    ungroup()
  
  clusters <- unique(summary_df$cluster)
  features <- unique(summary_df$feature)
  
  mat <- matrix(data = NA, nrow = 25, ncol = 13)
  rownames(mat) <- clusters; colnames(mat) <- features
  
  for (j in seq_along(1:nrow(summary_df))) {
    clus <- as.integer(summary_df[j, "cluster"])
    feat <- as.character(summary_df[j, "feature"])
    val <- as.numeric(summary_df[j, "mean_value"])
    mat[as.character(clus), feat] <- val
  }
  
  mat_scaled <- apply(mat, 2, scale)
  rownames(mat_scaled) <- rownames(mat)
  colnames(mat_scaled) <- colnames(mat)
  
  mat_scaled_long <- mat_scaled %>%
    as.data.frame() %>%
    mutate(cluster = rownames(mat_scaled)) %>%
    pivot_longer(
      cols = -cluster,                  # All columns except "cluster"
      names_to = "feature",
      values_to = "value"
    )
  
  # OLD don't consider here
  if(F){
    library(ComplexHeatmap)
    set.seed(123)
    # ComplexHeatmap::Heatmap(mat, cluster_columns = F)
    ht <- Heatmap(mat_scaled, cluster_columns = F)
    # ht
    draw_ht <- draw(ht)
    library(dendextend)
    row_hclust <- as.hclust(row_dend(draw_ht))
    
    k <- 2
    feature_clusters <- cutree(row_hclust, k = k)
    
    cluster_df <- tibble(
      feature = names(feature_clusters),
      cluster = factor(feature_clusters)
    )
    
    # Add matrix values
    mat_df <- as.data.frame(mat_scaled) %>%
      rownames_to_column("feature")
    
    # Combine and compute per-cluster means
    cluster_annotation <- cluster_df %>%
      left_join(mat_df, by = "feature") %>%
      group_by(cluster) %>%
      summarise(across(-feature, mean, .names = "mean_{.col}"))  # mean abundance per cluster
    
    # Convert to long format and scale within each feature
    annot_long <- cluster_annotation %>%
      pivot_longer(-cluster, names_to = "feature", values_to = "value") %>%
      group_by(feature) %>%
      mutate(scaled_value = scale(value)[,1])  # Z-score per feature
    
    # Test 1
    (top_features_per_cluster <- annot_long %>%
        group_by(cluster) %>%
        arrange(desc(value)) %>%
        slice_head(n = 3) %>%
        summarise(top_features = paste(feature, collapse = ", ")))
    # Test 2
    (top_differential_features <- annot_long %>%
        group_by(feature) %>%
        mutate(
          mean_others = (sum(value) - value) / (n() - 1),
          delta = abs(value - mean_others),
          direction = ifelse(value > mean_others, "HIGH", "LOW")
        ) %>%
        ungroup() %>%
        group_by(cluster) %>%
        arrange(desc(delta)) %>%
        slice_head(n = 3) %>%
        mutate(annotated_feature = paste0(feature, " [", direction, "]")) %>%
        summarise(top_features = paste(annotated_feature, collapse = ", ")))
    
    # Add annotation to the Heatmap
    row_ha <- rowAnnotation(
      cluster = cluster_df$cluster,
      annotation_name_side = "top"
    )
    
    Heatmap(mat_scaled, cluster_rows = row_hclust, right_annotation = row_ha,
            cluster_columns = F)
    
  }
  source('dev/12_hclust_nested.R')
  
  # build the heatmap with all the side annotations
  clusters <- as.data.frame(cbind(k2 = cluster_assignments$k2,
                                  k4 = cluster_assignments$k4,
                                  k8 = cluster_assignments$k8,
                                  k16 = cluster_assignments$k16))
  clusters$cluster <- as.integer(rownames(clusters))
  
  # values_agg <- full_join(clusters, values_agg, by = 'cluster') %>% filter(Type == 'BRCA')
  
  toplot.plot <- full_join(values_agg %>% select(labels,cluster),values %>% 
              select(bin,,Type,labels,ampl_score,del_score)) %>% 
    select(-labels) %>% full_join(clusters, by = 'cluster') %>%
    filter(Type == 'BRCA') %>% separate(bin, sep = '_', into = c('chr','bin'))
  
  
  # Plot feature contribution by cluster
  # library(scico)
  # region_clustering <- ggplot(mat_scaled_long, aes(x = factor(cluster), y = value, fill = feature)) +
  #    geom_bar(stat = "identity", position = "dodge") +
  #    labs(x = "Cluster", y = "Mean Feature Value", title = "Region Clustering by SHAP: Feature Contributions") +
  #    theme_minimal() +
  #    scale_fill_viridis_d(option = "H")
  #  scale_fill_scico_d(palette = "vik")
  # 
  # ggplot(mat_scaled_long, aes(x = feature, y = value, fill = as.factor(cluster))) +
  #   geom_bar(stat = "identity", position = "dodge") +
  #   labs(x = "Cluster", y = "Mean Feature Value", title = "Region Clustering by SHAP: Feature Contributions") +
  #   theme_minimal() +
  #   scale_fill_viridis_d(option = "H")
  # save(summary_df, file = 'summary_df.RData')
  # 
  # top_features <- mat_scaled_long %>%
  #   group_by(cluster) %>%
  #   top_n(2, wt = value) %>%
  #   arrange(cluster, desc(value)) %>%
  #   summarise(label = paste(feature, collapse = ", "))
  # # View(top_features)
  # 
  # aggregated <- top_features %>%
  #   group_by(label) %>%
  #   summarise(
  #     n_clusters = n(),
  #     clusters = paste(sort(unique(cluster)), collapse = ", ")
  #   ) %>%
  #   arrange(desc(n_clusters))
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
  
  # LANDSCAPE PLOT 
  # toplot.plot <- full_join(values_agg %>% select(labels,cluster),values %>% select(bin,,Type,labels,ampl_score,del_score)) %>% select(-labels) %>%
  #   filter(Type == 'BRCA') %>% separate(bin, sep = '_', into = c('chr','bin'))
  toplot.plot <- toplot.plot %>%
    arrange(as.numeric(as.character(chr)), as.numeric(bin)) %>%
    mutate(pos = row_number())
  
  # aggregated$clusters.aggregated <- rownames(aggregated)
  # cluster_mapping <- aggregated %>%
  #   separate_rows(clusters, sep = ",\\s*") %>%  # splits by comma and optional space
  #   mutate(clusters = as.integer(clusters)) %>%
  #   select(original_cluster = clusters, clusters.aggregated)
  # toplot.plot <- toplot.plot %>%
  #   left_join(cluster_mapping, by = c("cluster" = "original_cluster")) %>%
  #   mutate(cluster = clusters.aggregated) %>%
  #   select(-clusters.aggregated)
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
    mutate(cluster_ymin2 = max(toplot.plot$ampl_score) * 1.07 + k2 * 0.015, 
           cluster_ymax2 = cluster_ymin2 + 0.01) %>%
    mutate(cluster_ymin4 = max(toplot.plot$ampl_score) * 1.07 + (k4+4) * 0.015,
           cluster_ymax4 = cluster_ymin4 + 0.01 ) %>%
    mutate(cluster_ymin8 = max(toplot.plot$ampl_score) * 1.07 + (k8+10) * 0.015,
           cluster_ymax8 = cluster_ymin8 + 0.01) %>%
    mutate(cluster_ymin16 = max(toplot.plot$ampl_score) * 1.07 + (k16+20) * 0.015,
           cluster_ymax16 = cluster_ymin16 + 0.01)
  
  (p.final <- base_plot +
    geom_segment(data = cluster_ticks,
                 aes(x = pos, xend = pos, y = cluster_ymin2, yend = cluster_ymax2, col = as.character(k2)),
                 size = 0.3)+
      geom_segment(data = cluster_ticks,
                   aes(x = pos, xend = pos, y = cluster_ymin4, yend = cluster_ymax4, col = as.character(k4)),
                   size = 0.3)+
      geom_segment(data = cluster_ticks,
                   aes(x = pos, xend = pos, y = cluster_ymin8, yend = cluster_ymax8, col = as.character(k8)),
                   size = 0.3)+
      geom_segment(data = cluster_ticks,
                   aes(x = pos, xend = pos, y = cluster_ymin16, yend = cluster_ymax16, col = as.character(k16)),
                   size = 0.3) +
    scale_colour_viridis_d(option = "H"))
  
  output[[i]]$p.final <- p.final
  output[[i]]$toplot <- toplot.plot
  output[[i]]$aggregated <- annotations_list
}

# Visualize plots with annotation for ampl or del
output$ampl$p.final
output$ampl$aggregated

output$del$p.final
output$del$aggregated


# Explore Clusters
ggplot(output$ampl$toplot) +
  geom_boxplot(aes(x = as.factor(k4), y = ampl_score))
ggplot(output$ampl$toplot) +
  geom_boxplot(aes(x = as.factor(k8), y = ampl_score))
ggplot(output$ampl$toplot) +
  geom_boxplot(aes(x = as.factor(k16), y = ampl_score))

output$ampl$toplot %>% ggplot(aes(x = ampl_score, color = as.factor(k8), fill = as.factor(k8))) +
  geom_density(alpha = 0.4) +
  labs(
    x    = "Amplification score",
    y    = "Density",
    color = "Cluster (k = 8)",
    fill  = "Cluster (k = 8)"
  ) +
  theme_minimal() +
  xlim(-0.2,0.75)


output$ampl$toplot %>% group_by(k8) %>% 
  summarise(mean.ampl = mean(ampl_score),mean.del = mean(del_score)) %>%
  pivot_longer(cols = c(mean.ampl, mean.del)) %>%
  ggplot(aes(x = factor(k8), y = value, fill = factor(name))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  ggtitle('Amplification Model Clusters')
output$ampl$toplot %>% group_by(k16) %>% 
  summarise(mean.ampl = mean(ampl_score),mean.del = mean(del_score)) %>%
  pivot_longer(cols = c(mean.ampl, mean.del)) %>%
  ggplot(aes(x = factor(k16), y = value, fill = factor(name))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  ggtitle('Amplification Model Clusters')

output$del$toplot %>% group_by(k8) %>% 
  summarise(mean.ampl = mean(ampl_score),mean.del = mean(del_score)) %>%
  pivot_longer(cols = c(mean.ampl, mean.del)) %>%
  ggplot(aes(x = factor(k8), y = value, fill = factor(name))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  ggtitle('Amplification Model Clusters')
output$del$toplot %>% group_by(k16) %>% 
  summarise(mean.ampl = mean(ampl_score),mean.del = mean(del_score)) %>%
  pivot_longer(cols = c(mean.ampl, mean.del)) %>%
  ggplot(aes(x = factor(k16), y = value, fill = factor(name))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  ggtitle('Amplification Model Clusters')



# Data
# 1.  shap values
# 2.  feature matrix (toplot.plot)
# 3.  clustering + annotation
        # list(cluster_explained = aggregated,
        #      landscape = toplot.plot)
# 4.  hclust




