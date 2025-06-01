rm(list=ls())
gc(full=T)

library(ggplot2)
library(tidyr)
library(stringr)
library(tidyverse)

setwd("/Users/gabry/OneDrive/Desktop/shiny_app/")

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

i <- "del"

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

values <- left_join(res.filt[[i]] %>% select(bin, Type, residual), values)
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
             "enhancers", "repressed","accessible", 'labels')

shap_agg <- shap[, columns]
values_agg <- values[, columns]

num_of_top_features <- 5
actual_values_mat <- values_agg[,-ncol(values_agg)]
actual_shap_mat <- shap_agg[,-ncol(shap_agg)] 

res <- t(apply(X = actual_shap_mat, 
      MARGIN = 1, 
      FUN = function(row){
        
        idxs <- order(unlist(abs(row)))[1:num_of_top_features]
        cluster_identity <- paste0(idxs, collapse = "")
        
        signs <- sign(row[idxs]); signs <- gsub(pattern = "1", replacement = "", signs)
        top_features <- colnames(actual_shap_mat)[idxs]
        top_features <- paste0(signs, top_features)
        
        return(c(cluster_identity, top_features))
      
        }))

actual_shap_mat$cluster_id <- res[,1]

##########################################################################################

# HYPOTHESIS: Counter intuitive clusters actually make sense to be found, since the prediction happens to be 
              # in WHOLE GENOME PAN CANCER context, and we do not have a filtering step, based on p-value or log2FC
              # What I mean is that if bin3 of cancer X has a frequency of amplification of 0.1, and it belongs to cluster j
              # probably this would not be a high signal result, but may still be informative, 
              # pointing towards a cluster full of noisy regions
              # just randomly together, while some high signal regions may happen in more organized clusters

              # In other words, not all predicted regions will be tumor drivers, otherwise the WHOLE GENOME 
              # would be important for tumorigenesis. What I expect is that if we were the model, and tried to predict the amplification / deletion
              # frequency in a non important region, we would still be able to assign a value, which would probably be meaningless.

num_clusters <- length(table(res[,1]))

cumulative_members_sum <- unlist(lapply(X = seq(from = 1, to = num_clusters, by = 15), 
                          FUN = function(x){
                            sum(sort(table(res[,1]), decreasing = T)[1:x])
                            }))

plot(seq(from = 1, to = num_clusters, by = 15), cumulative_members_sum); 
  abline(v = 700, col = "red"); abline(h = 30000, col = "blue");
  abline(v = 3100, col = "red"); abline(h = 45000, col = "blue")

# What I do is basically count the number of regions for each initial cluster, just by merging regions with top k features similar
# Something interesting appears in the plot, since the vast majority of regions (30k - 45k, meaning from 50% to 75%) is
# explained by the minority of clusters (700 - 3100, meaning 5% - 22%)
# Are the clusters above the blue line the background noise????

putative_informative_clusters <- names(sort(table(res[,1]), decreasing = T))[1:700]
putative_noisy_clusters <- names(sort(table(res[,1]), decreasing = T))[701:num_clusters]

actual_shap_mat_new <- actual_shap_mat
actual_values_mat_new <- actual_values_mat

actual_shap_mat_new[actual_shap_mat_new$cluster_id %in% putative_informative_clusters, "putative_cluster_id"] <- 1
actual_shap_mat_new[actual_shap_mat_new$cluster_id %in% putative_noisy_clusters, "putative_cluster_id"] <- 2

actual_values_mat_new <- data.frame(
  apply(X = actual_values_mat_new,
        MARGIN = 2,
        FUN = scale)
)

actual_values_mat_new$cluster_id <- actual_shap_mat_new$cluster_id

# put the other option if you want to see the two clusters
actual_values_mat_new$putative_cluster_id <- 1 #actual_shap_mat$putative_cluster_id 


df_long_new <- actual_values_mat_new %>%
  select(-cluster_id) %>%
  pivot_longer(
    cols = -putative_cluster_id,
    names_to = "feature",
    values_to = "value"
  )

df_medians_new <- df_long_new %>%
  group_by(putative_cluster_id, feature) %>%
  summarise(median_value = median(value, na.rm = TRUE), .groups = "drop")

y_range_new <- c(-4.414780, 5.194962) # same range as below

plots <- df_medians_new %>%
  group_split(putative_cluster_id) %>%
  lapply(function(sub_df) {
    cluster_label <- unique(sub_df$putative_cluster_id)
    
    ggplot(sub_df, aes(x = feature, y = median_value)) +
      geom_col(aes(fill = feature)) +
      labs(title = paste("Cluster", cluster_label),
           x = "Feature", y = "Median Value") +
      coord_cartesian(ylim = y_range) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.position = "none")
  })

do.call(gridExtra::grid.arrange, c(plots, ncol = 1))

# It's clear that just the first approach does not separate driver events clearly, but still this
# (the two plots are very similar). Leveraging that we can estimate a genomic pancancer feature shape,
# to help interpreting the clusters. I would have expected two very different plots, but still it's not enough.
# So just reasoning, we can assume that shape to be a noisy shape, with mix up of several different regions 

##########################################################################################

actual_shap_mat_medians <- actual_shap_mat %>%
  group_by(cluster_id) %>%
  summarise(across(where(is.numeric), \(x) median(x, na.rm = TRUE)))

cluster_id_first_round <- actual_shap_mat_medians$cluster_id
actual_shap_mat_medians$cluster_id <- NULL
rownames(actual_shap_mat_medians) <- cluster_id_first_round

dist_matrix <- dist(actual_shap_mat_medians, method = "euclidean")
hc <- hclust(dist_matrix, method = "complete")

num_of_clusters = 16
new_cluster_ids <- cutree(hc, k = num_of_clusters)
actual_shap_mat_medians$new_cluster_id <- new_cluster_ids[rownames(actual_shap_mat_medians)]
actual_shap_mat_medians$cluster_id <- cluster_id_first_round

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

y_range <- range(df_medians$median_value, na.rm = TRUE)

plots <- df_medians %>%
  group_split(new_cluster_id) %>%
  lapply(function(sub_df) {
    cluster_label <- unique(sub_df$new_cluster_id)
    
    ggplot(sub_df, aes(x = feature, y = median_value)) +
      geom_col(aes(fill = feature)) +
      labs(title = paste("Cluster", cluster_label),
           x = "Feature", y = "Median Value") +
      coord_cartesian(ylim = y_range) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.position = "none")
  })

do.call(gridExtra::grid.arrange, c(plots, ncol = 4))

titles <- c("NOISE","NOISE","NOISE","NOISE",
            "NOISE","OG[LOW]; CENT[LOW]", "NOISE", "NOISE",
            "FGS[LOW]; OG[LOW]; ESS[LOW]","FGS[LOW]; OG[LOW]","D","OG[LOW]; TSG[HI]",
            "TRANSCRIBED[HI]","G","ENH[HI]","ESS[LOW]")

plots_with_titles <- Map(function(plot, title) {
  plot + labs(title = title)
}, plots, titles)

do.call(gridExtra::grid.arrange, c(plots_with_titles, ncol = 4))


