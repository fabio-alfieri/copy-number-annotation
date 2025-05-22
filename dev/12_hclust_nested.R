library(ComplexHeatmap)
library(dendextend)
library(tibble)
library(dplyr)
library(tidyr)

# 1. compute your row dendrogram once
set.seed(123)
ht <- Heatmap(mat_scaled, cluster_columns = FALSE)
drawn <- draw(ht)
row_h <- as.hclust(row_dend(drawn))

# 2. prepare storage
ks <- c(2, 4, 8, 16)
cluster_assignments <- list()
annotations_list    <- list()

# 3. loop over each k, now taking the top 2 features per cluster
for(k in ks) {
  # a) cut tree
  cl <- cutree(row_h, k = k)
  cluster_assignments[[paste0("k",k)]] <- cl
  
  # b) compute per-cluster top-2 differential features
  df_cl <- tibble(feature = names(cl), cluster = cl) %>%
    left_join(
      as.data.frame(mat_scaled) %>% rownames_to_column("feature"),
      by = "feature"
    )
  
  annot_k <- df_cl %>%
    group_by(cluster) %>% summarise(across(-feature, mean, .names = "{col}")) %>%
    pivot_longer(-cluster, names_to = 'feature', values_to = 'value') %>%
    group_by(feature) %>%
    mutate(
      mean_others = (sum(value) - value) / (n() - 1),
      delta = abs(value - mean_others),
      direction = ifelse(value > mean_others, "HIGH", "LOW")
    ) %>% # filter(direction == 'HIGH') %>%
    ungroup() %>%
    group_by(cluster) %>%
    arrange(desc(delta)) %>%
    slice_head(n = 4) %>%
    mutate(annotated_feature = paste0(feature, " [", direction, "]")) %>%
    summarise(top_features = paste(annotated_feature, collapse = ", "))
  
  annotations_list[[paste0("k",k)]] <- annot_k
}

# 4. build hierarchical path for each feature
features <- names(cluster_assignments$k2)
hierarchy <- tibble(feature = features) %>%
  rowwise() %>%
  mutate(
    path = paste0(
      map_chr(ks, function(k) {
        ktag <- paste0("k", k)
        # which cluster this feature falls into at level k
        clu  <- cluster_assignments[[ktag]][feature]
        # pull the correct column (`top_features`)
        annotations_list[[ktag]] %>%
          filter(cluster == clu) %>%
          pull(top_features)
      }),
      collapse = " > "
    )
  ) %>%
  ungroup()

# 5. annotate heatmap
row_ha <- rowAnnotation(
  split_hierarchy = hierarchy$path,
  annotation_name_side = "top"
)

f <- Heatmap(mat_scaled,
        cluster_rows     = row_h,
        right_annotation = row_ha,
        cluster_columns  = FALSE)

print(f)
