pattern_shap <- "Mid-length"
pattern_predictions <- "Mid-length(_[0-9]+)?-pred_ampl"

model.outputs_SHAP <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/SHAP_and_featurematrix/chromosome_specific", pattern = pattern_shap, recursive = T, full.names = T)
#model.outputs_SHAP <- model.outputs_SHAP[-c(1,25,26)]

types <- lapply(X = model.outputs_SHAP, FUN = function(x){strsplit(x = x, split = "/")}[[1]][7])
models <- gsub(pattern = "_AmplDel.rds", replacement = "", x = lapply(X = model.outputs_SHAP, FUN = function(x){strsplit(x = x, split = "/")}[[1]][8]))
models <- gsub(pattern = "SHAP_and_FeatureMatrix_", replacement = "", x = models)
#models[1] <- paste0(models[1], "_WG")
#models[length(models)] <- paste0(models[length(models)], "_EDITED")

############ SEPARATE CHROMOSOME ANALYSIS

model.outputs_SHAP_chromosome_specific <- model.outputs_SHAP[which(types == "chromosome_specific")]

# models[-1] drops the WG model (index 1); chromosome-specific models start at index 2
rankings_ampl <- lapply(X = seq_along(model.outputs_SHAP_chromosome_specific), 
                        FUN = function(idx){
                          shap_and_features <- readRDS(model.outputs_SHAP_chromosome_specific[idx])
                          shap_ampl <- shap_and_features$models.shap.df[[1]]
                          imp_ampl <- apply(abs(shap_ampl[, -1]), 2, sum)
                          ranking_ampl <- data.frame(feature = names(imp_ampl), importance = imp_ampl, rank = rank(-imp_ampl, ties.method = "min"))
                          ranking_ampl <- ranking_ampl[order(ranking_ampl$rank), ]
                          ranking_ampl$model <- models[-1][idx]
                          return(ranking_ampl)
                        }
)

rankings_del <- lapply(X = seq_along(model.outputs_SHAP_chromosome_specific), 
                       FUN = function(idx){
                         shap_and_features <- readRDS(model.outputs_SHAP_chromosome_specific[idx])
                         shap_del <- shap_and_features$models.shap.df[[2]]
                         imp_del <- apply(abs(shap_del[, -1]), 2, sum)
                         ranking_del <- data.frame(feature = names(imp_del), importance = imp_del, rank = rank(-imp_del, ties.method = "min"))
                         ranking_del <- ranking_del[order(ranking_del$rank), ]
                         ranking_del$model <- models[-1][idx]
                         return(ranking_del)
                       }
)

library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(grid)          # required by save_pheatmap_pair (grid.newpage, pushViewport, etc.)
library(ComplexHeatmap) # required by save_pheatmap_pair (draw())

rankings_df_ampl <- bind_rows(rankings_ampl)
rankings_df_del <- bind_rows(rankings_del)

feature_order_ampl <- rankings_df_ampl %>% group_by(feature) %>% summarise(mean_rank = mean(rank, na.rm = TRUE)) %>% arrange(mean_rank) %>% pull(feature)
feature_order_del <- rankings_df_del %>% group_by(feature) %>% summarise(mean_rank = mean(rank, na.rm = TRUE)) %>% arrange(mean_rank) %>% pull(feature)

rankings_df_ampl$feature <- factor(rankings_df_ampl$feature, levels = feature_order_ampl)
rankings_df_del$feature <- factor(rankings_df_del$feature, levels = feature_order_del)

rank_matrix_ampl <- rankings_df_ampl %>% dplyr::select(feature, model, rank) %>% tidyr::pivot_wider(names_from = model, values_from = rank) %>% tibble::column_to_rownames("feature") %>% as.matrix()
rank_matrix_del <- rankings_df_del %>% dplyr::select(feature, model, rank) %>% tidyr::pivot_wider(names_from = model, values_from = rank) %>% tibble::column_to_rownames("feature") %>% as.matrix()

sd_rank_ampl <- apply(X = rank_matrix_ampl, MARGIN = 1, FUN = sd)
sd_rank_del <- apply(X = rank_matrix_del, MARGIN = 1, FUN = sd)

rankings_df_ampl$sd_rank <- sd_rank_ampl[as.character(rankings_df_ampl$feature)]
rankings_df_del$sd_rank <- sd_rank_del[as.character(rankings_df_del$feature)]

plot_df_ampl <- rankings_df_ampl %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(importance, na.rm = TRUE),
    median_importance = median(importance, na.rm = TRUE),
    sd_rank = first(sd_rank)
  )

plot_df_del <- rankings_df_del %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(importance, na.rm = TRUE),
    median_importance = median(importance, na.rm = TRUE),
    sd_rank = first(sd_rank)
  )

df_combined <- bind_rows(
  data.frame(Feature = names(sd_rank_ampl), SD = sd_rank_ampl, Type = "Amplification"),
  data.frame(Feature = names(sd_rank_del),  SD = sd_rank_del,  Type = "Deletion")
)

cor_ampl <- cor.test(plot_df_ampl$median_importance, plot_df_ampl$sd_rank, method = "pearson")
r_ampl <- round(cor_ampl$estimate, 3)
p_ampl <- signif(cor_ampl$p.value, 3)

plot_df_ampl_lm <- plot_df_ampl %>% filter(median_importance > 0)

cor_del <- cor.test(plot_df_del$median_importance, plot_df_del$sd_rank, method = "pearson")
r_del <- round(cor_del$estimate, 3)
p_del <- signif(cor_del$p.value, 3)

plot_df_del_lm <- plot_df_del %>% filter(median_importance > 0)

pca_ampl <- prcomp(x = t(rank_matrix_ampl), center = TRUE, scale. = TRUE)

pca_df_ampl <- data.frame(
  model = rownames(pca_ampl$x),
  PC1   = pca_ampl$x[, 1],
  PC2   = pca_ampl$x[, 2]
)

var_explained <- summary(pca_ampl)$importance[2, ] * 100

ranking_plot_ampl_chr <- ggplot(rankings_df_ampl,
                                aes(x = feature, y = rank, fill = model)) +
  geom_col(position = position_dodge(width = 0.9)) +
  coord_flip() +
  scale_y_reverse() +
  labs(x = "Feature", y = "Rank", fill = "Model") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8), legend.position = "right")

ranking_plot_del_chr <- ggplot(rankings_df_del,
                               aes(x = feature, y = rank, fill = model)) +
  geom_col(position = position_dodge(width = 0.9)) +
  coord_flip() +
  scale_y_reverse() +
  labs(x = "Feature", y = "Rank", fill = "Model") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8), legend.position = "right")

rank_heatmap_ampl_chr <- pheatmap(
  rank_matrix_ampl,
  cluster_rows = FALSE, cluster_cols = TRUE, scale = "none",
  color = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))(40)
)

rank_heatmap_del_chr <- pheatmap(
  rank_matrix_del,
  cluster_rows = FALSE, cluster_cols = TRUE, scale = "none",
  color = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))(40)
)

sd_plot_chr <- ggplot(df_combined, aes(x = reorder(Feature, SD), y = SD, fill = Type)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Variability of Feature Importance Ranks", x = "Feature", y = "SD of Rank") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

library(ggplot2)
library(ggrepel)
library(scales)

scatterplot_ampl_chr <- ggplot(
  plot_df_ampl,
  aes(x = median_importance, y = sd_rank)
) +
  geom_point(
    size = 3,
    alpha = 0.8,
    colour = "#2C7FB8"
  ) +
  geom_text_repel(
    data = plot_df_ampl_lm,
    aes(label = feature),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.color = "grey60"
  ) +
  geom_smooth(
    data = plot_df_ampl_lm,
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = "#D7301F",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.05,
    vjust = 1.2,
    size = 4,
    label = sprintf(
      "Pearson r = %.2f\np = %.2g",
      r_ampl,
      p_ampl
    )
  ) +
  scale_x_log10(
    name = "Median feature importance (log-scale)",
    breaks = scales::breaks_log(n = 6),
    labels = scales::label_number()
  ) +
  scale_y_continuous(
    name = "Stability of feature across chromosomes"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    plot.margin = margin(10, 15, 10, 10)
  )




scatterplot_del_chr <- ggplot(plot_df_del, aes(x = median_importance, y = sd_rank)) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(data = plot_df_del_lm, aes(label = feature), size = 3, max.overlaps = Inf) +
  geom_smooth(data = plot_df_del_lm, method = "lm", formula = y ~ x, se = TRUE, color = "red") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("Pearson r = ", r_del, "\np = ", p_del)) +
  scale_x_log10() +
  theme_bw()

pca_plot_ampl <- ggplot(pca_df_ampl, aes(x = PC1, y = PC2, label = model)) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(
    title = "PCA of Feature Importance Ranks across Models (Amplification)",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

############ SHARED HELPER: per-chromosome rankings from a single SHAP matrix

# Defined once; reused by both WG and EDITED sections.
# `label_prefix` becomes the model column value, e.g. "Mid-length_WG_chrX"
compute_rankings_by_chr <- function(shap_df, chr_vec, label_prefix) {
  lapply(chr_vec, function(chr) {
    rows       <- extract_chr(shap_df$labels) == chr
    sub_shap   <- shap_df[rows, , drop = FALSE]
    feat_mat   <- sub_shap[, colnames(sub_shap) != "labels", drop = FALSE]
    imp        <- colSums(abs(feat_mat))
    data.frame(
      feature    = names(imp),
      importance = imp,
      rank       = rank(-imp, ties.method = "min"),
      model      = paste0(label_prefix, "_chr", chr),
      row.names  = NULL
    ) %>% arrange(rank)
  })
}

extract_chr <- function(labels) sub("_.*", "", labels)

############ WHOLE GENOME ANALYSIS

shap_and_features_WG <- readRDS(model.outputs_SHAP[1])

shap_ampl_WG <- shap_and_features_WG$models.shap.df[[1]]
shap_del_WG  <- shap_and_features_WG$models.shap.df[[2]]

feat_ampl_WG <- shap_and_features_WG$models.X[[1]][, -c(3, 4)]
feat_del_WG  <- shap_and_features_WG$models.X[[2]][, -c(3, 4)]

feat_ampl_WG$labels <- paste0(feat_ampl_WG$bin, "-", feat_ampl_WG$Type)
feat_ampl_WG <- feat_ampl_WG[, -c(1, 2)]
feat_ampl_WG <- feat_ampl_WG[, colnames(shap_ampl_WG)]

feat_del_WG$labels <- paste0(feat_del_WG$bin, "-", feat_del_WG$Type)
feat_del_WG <- feat_del_WG[, -c(1, 2)]
feat_del_WG <- feat_del_WG[, colnames(shap_del_WG)]

stopifnot(all(colnames(shap_ampl_WG) == colnames(feat_ampl_WG)))
stopifnot(all(colnames(shap_del_WG)  == colnames(feat_del_WG)))

chromosomes_WG <- sort(unique(extract_chr(shap_ampl_WG$labels)))

rankings_ampl_WG <- compute_rankings_by_chr(shap_ampl_WG, chromosomes_WG, paste0(pattern_shap, "_WG"))
rankings_del_WG  <- compute_rankings_by_chr(shap_del_WG,  chromosomes_WG, paste0(pattern_shap, "_WG"))

rankings_df_ampl_WG <- bind_rows(rankings_ampl_WG)
rankings_df_del_WG  <- bind_rows(rankings_del_WG)

feature_order_ampl_WG <- rankings_df_ampl_WG %>%
  group_by(feature) %>% summarise(mean_rank = mean(rank, na.rm = TRUE)) %>%
  arrange(mean_rank) %>% pull(feature)

feature_order_del_WG <- rankings_df_del_WG %>%
  group_by(feature) %>% summarise(mean_rank = mean(rank, na.rm = TRUE)) %>%
  arrange(mean_rank) %>% pull(feature)

rankings_df_ampl_WG$feature <- factor(rankings_df_ampl_WG$feature, levels = feature_order_ampl_WG)
rankings_df_del_WG$feature  <- factor(rankings_df_del_WG$feature,  levels = feature_order_del_WG)

rank_matrix_ampl_WG <- rankings_df_ampl_WG %>%
  dplyr::select(feature, model, rank) %>%
  tidyr::pivot_wider(names_from = model, values_from = rank) %>%
  tibble::column_to_rownames("feature") %>% as.matrix()

rank_matrix_del_WG <- rankings_df_del_WG %>%
  dplyr::select(feature, model, rank) %>%
  tidyr::pivot_wider(names_from = model, values_from = rank) %>%
  tibble::column_to_rownames("feature") %>% as.matrix()

sd_rank_ampl_WG <- apply(rank_matrix_ampl_WG, 1, sd, na.rm = TRUE)
sd_rank_del_WG  <- apply(rank_matrix_del_WG,  1, sd, na.rm = TRUE)

rankings_df_ampl_WG$sd_rank <- sd_rank_ampl_WG[as.character(rankings_df_ampl_WG$feature)]
rankings_df_del_WG$sd_rank  <- sd_rank_del_WG[as.character(rankings_df_del_WG$feature)]

plot_df_ampl_WG <- rankings_df_ampl_WG %>%
  group_by(feature) %>%
  summarise(mean_importance = mean(importance, na.rm = TRUE),
            median_importance = median(importance, na.rm = TRUE),
            sd_rank = first(sd_rank))

plot_df_del_WG <- rankings_df_del_WG %>%
  group_by(feature) %>%
  summarise(mean_importance = mean(importance, na.rm = TRUE),
            median_importance = median(importance, na.rm = TRUE),
            sd_rank = first(sd_rank))

df_combined_WG <- bind_rows(
  data.frame(Feature = names(sd_rank_ampl_WG), SD = sd_rank_ampl_WG, Type = "Amplification"),
  data.frame(Feature = names(sd_rank_del_WG),  SD = sd_rank_del_WG,  Type = "Deletion")
)

cor_ampl_WG <- cor.test(plot_df_ampl_WG$median_importance, plot_df_ampl_WG$sd_rank, method = "pearson")
r_ampl_WG   <- round(cor_ampl_WG$estimate, 3)
p_ampl_WG   <- signif(cor_ampl_WG$p.value, 3)

cor_del_WG  <- cor.test(plot_df_del_WG$median_importance, plot_df_del_WG$sd_rank, method = "pearson")
r_del_WG    <- round(cor_del_WG$estimate, 3)
p_del_WG    <- signif(cor_del_WG$p.value, 3)

plot_df_ampl_WG_lm <- plot_df_ampl_WG %>% filter(median_importance > 0)
plot_df_del_WG_lm  <- plot_df_del_WG  %>% filter(median_importance > 0)

ranking_plot_ampl_WG <- ggplot(rankings_df_ampl_WG, aes(x = feature, y = rank, fill = model)) +
  geom_col(position = position_dodge(width = 0.9)) +
  coord_flip() + scale_y_reverse() +
  labs(title = "WG Model – Feature Ranks per Chromosome (Amplification)",
       x = "Feature", y = "Rank", fill = "Chromosome") +
  theme_bw() + theme(axis.text.y = element_text(size = 8), legend.position = "right")

ranking_plot_del_WG <- ggplot(rankings_df_del_WG, aes(x = feature, y = rank, fill = model)) +
  geom_col(position = position_dodge(width = 0.9)) +
  coord_flip() + scale_y_reverse() +
  labs(title = "WG Model – Feature Ranks per Chromosome (Deletion)",
       x = "Feature", y = "Rank", fill = "Chromosome") +
  theme_bw() + theme(axis.text.y = element_text(size = 8), legend.position = "right")

rank_heatmap_ampl_WG <- pheatmap(
  rank_matrix_ampl_WG, cluster_rows = FALSE, cluster_cols = TRUE, scale = "none",
  main  = "WG Model – Feature Rank Heatmap (Amplification)",
  color = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))(40)
)

rank_heatmap_del_WG <- pheatmap(
  rank_matrix_del_WG, cluster_rows = FALSE, cluster_cols = TRUE, scale = "none",
  main  = "WG Model – Feature Rank Heatmap (Deletion)",
  color = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))(40)
)

sd_plot_WG <- ggplot(df_combined_WG, aes(x = reorder(Feature, SD), y = SD, fill = Type)) +
  geom_col(position = "dodge") + coord_flip() +
  labs(title = "WG Model – Variability of Feature Importance Ranks",
       x = "Feature", y = "SD of Rank") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

scatterplot_ampl_WG <- ggplot(plot_df_ampl_WG, aes(x = median_importance, y = sd_rank)) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(data = plot_df_ampl_WG_lm, aes(label = feature), size = 3, max.overlaps = Inf) +
  geom_smooth(data = plot_df_ampl_WG_lm, method = "lm", formula = y ~ x, se = TRUE, color = "red") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("Pearson r = ", r_ampl_WG, "\np = ", p_ampl_WG)) +
  scale_x_log10() +
  labs(title = "WG Model – Median Importance vs Rank SD (Amplification)",
       x = "Median SHAP Importance (log10)", y = "SD of Rank") +
  theme_bw()

scatterplot_del_WG <- ggplot(plot_df_del_WG, aes(x = median_importance, y = sd_rank)) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(data = plot_df_del_WG_lm, aes(label = feature), size = 3, max.overlaps = Inf) +
  geom_smooth(data = plot_df_del_WG_lm, method = "lm", formula = y ~ x, se = TRUE, color = "red") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("Pearson r = ", r_del_WG, "\np = ", p_del_WG)) +
  scale_x_log10() +
  labs(title = "WG Model – Median Importance vs Rank SD (Deletion)",
       x = "Median SHAP Importance (log10)", y = "SD of Rank") +
  theme_bw()

############ EDITED ANALYSIS

shap_and_features_EDITED <- readRDS(model.outputs_SHAP[length(model.outputs_SHAP)])

shap_ampl_EDITED <- shap_and_features_EDITED$models.shap.df[[1]]
shap_del_EDITED  <- shap_and_features_EDITED$models.shap.df[[2]]

feat_ampl_EDITED <- shap_and_features_EDITED$models.X[[1]][, -c(3, 4)]
feat_del_EDITED  <- shap_and_features_EDITED$models.X[[2]][, -c(3, 4)]

feat_ampl_EDITED$labels <- paste0(feat_ampl_EDITED$bin, "-", feat_ampl_EDITED$Type)
feat_ampl_EDITED <- feat_ampl_EDITED[, -c(1, 2)]
feat_ampl_EDITED <- feat_ampl_EDITED[, colnames(shap_ampl_EDITED)]

feat_del_EDITED$labels <- paste0(feat_del_EDITED$bin, "-", feat_del_EDITED$Type)
feat_del_EDITED <- feat_del_EDITED[, -c(1, 2)]
feat_del_EDITED <- feat_del_EDITED[, colnames(shap_del_EDITED)]

stopifnot(all(colnames(shap_ampl_EDITED) == colnames(feat_ampl_EDITED)))
stopifnot(all(colnames(shap_del_EDITED)  == colnames(feat_del_EDITED)))

# extract_chr is already defined above; chromosomes derived from EDITED SHAP labels
chromosomes_EDITED <- sort(unique(extract_chr(shap_ampl_EDITED$labels)))

# Reuse the shared helper — only the label prefix differs from WG
rankings_ampl_EDITED <- compute_rankings_by_chr(shap_ampl_EDITED, chromosomes_EDITED, paste0(pattern_shap, "_EDITED"))
rankings_del_EDITED  <- compute_rankings_by_chr(shap_del_EDITED,  chromosomes_EDITED, paste0(pattern_shap, "_EDITED"))

rankings_df_ampl_EDITED <- bind_rows(rankings_ampl_EDITED)
rankings_df_del_EDITED  <- bind_rows(rankings_del_EDITED)

feature_order_ampl_EDITED <- rankings_df_ampl_EDITED %>%
  group_by(feature) %>% summarise(mean_rank = mean(rank, na.rm = TRUE)) %>%
  arrange(mean_rank) %>% pull(feature)

feature_order_del_EDITED <- rankings_df_del_EDITED %>%
  group_by(feature) %>% summarise(mean_rank = mean(rank, na.rm = TRUE)) %>%
  arrange(mean_rank) %>% pull(feature)

rankings_df_ampl_EDITED$feature <- factor(rankings_df_ampl_EDITED$feature, levels = feature_order_ampl_EDITED)
rankings_df_del_EDITED$feature  <- factor(rankings_df_del_EDITED$feature,  levels = feature_order_del_EDITED)

rank_matrix_ampl_EDITED <- rankings_df_ampl_EDITED %>%
  dplyr::select(feature, model, rank) %>%
  tidyr::pivot_wider(names_from = model, values_from = rank) %>%
  tibble::column_to_rownames("feature") %>% as.matrix()

rank_matrix_del_EDITED <- rankings_df_del_EDITED %>%
  dplyr::select(feature, model, rank) %>%
  tidyr::pivot_wider(names_from = model, values_from = rank) %>%
  tibble::column_to_rownames("feature") %>% as.matrix()

sd_rank_ampl_EDITED <- apply(rank_matrix_ampl_EDITED, 1, sd, na.rm = TRUE)
sd_rank_del_EDITED  <- apply(rank_matrix_del_EDITED,  1, sd, na.rm = TRUE)

rankings_df_ampl_EDITED$sd_rank <- sd_rank_ampl_EDITED[as.character(rankings_df_ampl_EDITED$feature)]
rankings_df_del_EDITED$sd_rank  <- sd_rank_del_EDITED[as.character(rankings_df_del_EDITED$feature)]

plot_df_ampl_EDITED <- rankings_df_ampl_EDITED %>%
  group_by(feature) %>%
  summarise(mean_importance = mean(importance, na.rm = TRUE),
            median_importance = median(importance, na.rm = TRUE),
            sd_rank = first(sd_rank))

plot_df_del_EDITED <- rankings_df_del_EDITED %>%
  group_by(feature) %>%
  summarise(mean_importance = mean(importance, na.rm = TRUE),
            median_importance = median(importance, na.rm = TRUE),
            sd_rank = first(sd_rank))

df_combined_EDITED <- bind_rows(
  data.frame(Feature = names(sd_rank_ampl_EDITED), SD = sd_rank_ampl_EDITED, Type = "Amplification"),
  data.frame(Feature = names(sd_rank_del_EDITED),  SD = sd_rank_del_EDITED,  Type = "Deletion")
)

cor_ampl_EDITED <- cor.test(plot_df_ampl_EDITED$median_importance, plot_df_ampl_EDITED$sd_rank, method = "pearson")
r_ampl_EDITED   <- round(cor_ampl_EDITED$estimate, 3)
p_ampl_EDITED   <- signif(cor_ampl_EDITED$p.value, 3)

cor_del_EDITED  <- cor.test(plot_df_del_EDITED$median_importance, plot_df_del_EDITED$sd_rank, method = "pearson")
r_del_EDITED    <- round(cor_del_EDITED$estimate, 3)
p_del_EDITED    <- signif(cor_del_EDITED$p.value, 3)

plot_df_ampl_EDITED_lm <- plot_df_ampl_EDITED %>% filter(median_importance > 0)
plot_df_del_EDITED_lm  <- plot_df_del_EDITED  %>% filter(median_importance > 0)

ranking_plot_ampl_EDITED <- ggplot(rankings_df_ampl_EDITED, aes(x = feature, y = rank, fill = model)) +
  geom_col(position = position_dodge(width = 0.9)) +
  coord_flip() + scale_y_reverse() +
  labs(title = "EDITED Model – Feature Ranks per Chromosome (Amplification)",
       x = "Feature", y = "Rank", fill = "Chromosome") +
  theme_bw() + theme(axis.text.y = element_text(size = 8), legend.position = "right")

ranking_plot_del_EDITED <- ggplot(rankings_df_del_EDITED, aes(x = feature, y = rank, fill = model)) +
  geom_col(position = position_dodge(width = 0.9)) +
  coord_flip() + scale_y_reverse() +
  labs(title = "EDITED Model – Feature Ranks per Chromosome (Deletion)",
       x = "Feature", y = "Rank", fill = "Chromosome") +
  theme_bw() + theme(axis.text.y = element_text(size = 8), legend.position = "right")

rank_heatmap_ampl_EDITED <- pheatmap(
  rank_matrix_ampl_EDITED, cluster_rows = FALSE, cluster_cols = TRUE, scale = "none",
  main  = "EDITED Model – Feature Rank Heatmap (Amplification)",
  color = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))(40)
)

rank_heatmap_del_EDITED <- pheatmap(
  rank_matrix_del_EDITED, cluster_rows = FALSE, cluster_cols = TRUE, scale = "none",
  main  = "EDITED Model – Feature Rank Heatmap (Deletion)",
  color = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))(40)
)

sd_plot_EDITED <- ggplot(df_combined_EDITED, aes(x = reorder(Feature, SD), y = SD, fill = Type)) +
  geom_col(position = "dodge") + coord_flip() +
  labs(title = "EDITED Model – Variability of Feature Importance Ranks",
       x = "Feature", y = "SD of Rank") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

scatterplot_ampl_EDITED <- ggplot(plot_df_ampl_EDITED, aes(x = median_importance, y = sd_rank)) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(data = plot_df_ampl_EDITED_lm, aes(label = feature), size = 3, max.overlaps = Inf) +
  geom_smooth(data = plot_df_ampl_EDITED_lm, method = "lm", formula = y ~ x, se = TRUE, color = "red") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("Pearson r = ", r_ampl_EDITED, "\np = ", p_ampl_EDITED)) +
  scale_x_log10() +
  labs(title = "EDITED Model – Median Importance vs Rank SD (Amplification)",
       x = "Median SHAP Importance (log10)", y = "SD of Rank") +
  theme_bw()

scatterplot_del_EDITED <- ggplot(plot_df_del_EDITED, aes(x = median_importance, y = sd_rank)) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(data = plot_df_del_EDITED_lm, aes(label = feature), size = 3, max.overlaps = Inf) +
  geom_smooth(data = plot_df_del_EDITED_lm, method = "lm", formula = y ~ x, se = TRUE, color = "red") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("Pearson r = ", r_del_EDITED, "\np = ", p_del_EDITED)) +
  scale_x_log10() +
  labs(title = "EDITED Model – Median Importance vs Rank SD (Deletion)",
       x = "Median SHAP Importance (log10)", y = "SD of Rank") +
  theme_bw()

############ SAVE ALL PAIRED COMPARISON PLOTS

library(patchwork)

output_dir <- "/home/ieo7429/Scrivania/results_regressor_gab/comparison_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Chr-specific vs WG ----

ggsave(file.path(output_dir, "scatterplot_ampl_chr_vs_WG.pdf"),
       scatterplot_ampl_chr + scatterplot_ampl_WG + plot_layout(ncol = 2),
       width = 16, height = 7)

ggsave(file.path(output_dir, "scatterplot_del_chr_vs_WG.pdf"),
       scatterplot_del_chr + scatterplot_del_WG + plot_layout(ncol = 2),
       width = 16, height = 7)

ggsave(file.path(output_dir, "ranking_plot_ampl_chr_vs_WG.pdf"),
       ranking_plot_ampl_chr + ranking_plot_ampl_WG + plot_layout(ncol = 2),
       width = 22, height = 12)

ggsave(file.path(output_dir, "ranking_plot_del_chr_vs_WG.pdf"),
       ranking_plot_del_chr + ranking_plot_del_WG + plot_layout(ncol = 2),
       width = 22, height = 12)

ggsave(file.path(output_dir, "sd_plot_chr_vs_WG.pdf"),
       sd_plot_chr + sd_plot_WG + plot_layout(ncol = 2),
       width = 16, height = 7)

# ---- Chr-specific vs EDITED ----

ggsave(file.path(output_dir, "scatterplot_ampl_chr_vs_EDITED.pdf"),
       scatterplot_ampl_chr + scatterplot_ampl_EDITED + plot_layout(ncol = 2),
       width = 16, height = 7)

ggsave(file.path(output_dir, "scatterplot_del_chr_vs_EDITED.pdf"),
       scatterplot_del_chr + scatterplot_del_EDITED + plot_layout(ncol = 2),
       width = 16, height = 7)

ggsave(file.path(output_dir, "ranking_plot_ampl_chr_vs_EDITED.pdf"),
       ranking_plot_ampl_chr + ranking_plot_ampl_EDITED + plot_layout(ncol = 2),
       width = 22, height = 12)

ggsave(file.path(output_dir, "ranking_plot_del_chr_vs_EDITED.pdf"),
       ranking_plot_del_chr + ranking_plot_del_EDITED + plot_layout(ncol = 2),
       width = 22, height = 12)

ggsave(file.path(output_dir, "sd_plot_chr_vs_EDITED.pdf"),
       sd_plot_chr + sd_plot_EDITED + plot_layout(ncol = 2),
       width = 16, height = 7)

# ---- WG vs EDITED ----

ggsave(file.path(output_dir, "scatterplot_ampl_WG_vs_EDITED.pdf"),
       scatterplot_ampl_WG + scatterplot_ampl_EDITED + plot_layout(ncol = 2),
       width = 16, height = 7)

ggsave(file.path(output_dir, "scatterplot_del_WG_vs_EDITED.pdf"),
       scatterplot_del_WG + scatterplot_del_EDITED + plot_layout(ncol = 2),
       width = 16, height = 7)

ggsave(file.path(output_dir, "ranking_plot_ampl_WG_vs_EDITED.pdf"),
       ranking_plot_ampl_WG + ranking_plot_ampl_EDITED + plot_layout(ncol = 2),
       width = 22, height = 12)

ggsave(file.path(output_dir, "ranking_plot_del_WG_vs_EDITED.pdf"),
       ranking_plot_del_WG + ranking_plot_del_EDITED + plot_layout(ncol = 2),
       width = 22, height = 12)

ggsave(file.path(output_dir, "sd_plot_WG_vs_EDITED.pdf"),
       sd_plot_WG + sd_plot_EDITED + plot_layout(ncol = 2),
       width = 16, height = 7)

# ---- Heatmap pairs (pheatmap objects require pdf + grid layout) ----

save_pheatmap_pair <- function(hm1, hm2, filename, width = 22, height = 12,
                               label1 = "Chr-specific", label2 = "WG") {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  ComplexHeatmap::draw(hm1, newpage = FALSE)
  grid::grid.text(label1, y = grid::unit(1, "npc") + grid::unit(0.2, "cm"),
                  gp = grid::gpar(fontsize = 14, fontface = "bold"))
  grid::popViewport()
  
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  ComplexHeatmap::draw(hm2, newpage = FALSE)
  grid::grid.text(label2, y = grid::unit(1, "npc") + grid::unit(0.2, "cm"),
                  gp = grid::gpar(fontsize = 14, fontface = "bold"))
  grid::popViewport()
  
  dev.off()
}

# Chr vs WG
save_pheatmap_pair(rank_heatmap_ampl_chr, rank_heatmap_ampl_WG,
                   file.path(output_dir, "rank_heatmap_ampl_chr_vs_WG.pdf"),
                   label1 = "Chr-specific", label2 = "WG")

save_pheatmap_pair(rank_heatmap_del_chr, rank_heatmap_del_WG,
                   file.path(output_dir, "rank_heatmap_del_chr_vs_WG.pdf"),
                   label1 = "Chr-specific", label2 = "WG")

# Chr vs EDITED
save_pheatmap_pair(rank_heatmap_ampl_chr, rank_heatmap_ampl_EDITED,
                   file.path(output_dir, "rank_heatmap_ampl_chr_vs_EDITED.pdf"),
                   label1 = "Chr-specific", label2 = "EDITED")

save_pheatmap_pair(rank_heatmap_del_chr, rank_heatmap_del_EDITED,
                   file.path(output_dir, "rank_heatmap_del_chr_vs_EDITED.pdf"),
                   label1 = "Chr-specific", label2 = "EDITED")

# WG vs EDITED
save_pheatmap_pair(rank_heatmap_ampl_WG, rank_heatmap_ampl_EDITED,
                   file.path(output_dir, "rank_heatmap_ampl_WG_vs_EDITED.pdf"),
                   label1 = "WG", label2 = "EDITED")

save_pheatmap_pair(rank_heatmap_del_WG, rank_heatmap_del_EDITED,
                   file.path(output_dir, "rank_heatmap_del_WG_vs_EDITED.pdf"),
                   label1 = "WG", label2 = "EDITED")

###############################################################################
## COMPARISON: CHROMOSOME-SPECIFIC vs WHOLE-GENOME SHAP IMPORTANCE
###############################################################################

library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

## ---- 1. Whole-genome importance per feature (unchanged) ----

feature_importance_ampl <- colSums(abs(shap_ampl_WG[, -1]), na.rm = TRUE)
feature_importance_del  <- colSums(abs(shap_del_WG[, -1]),  na.rm = TRUE)

wg_df_ampl <- data.frame(
  feature       = names(feature_importance_ampl),
  WG_importance = as.numeric(feature_importance_ampl),
  stringsAsFactors = FALSE
)

wg_df_del <- data.frame(
  feature       = names(feature_importance_del),
  WG_importance = as.numeric(feature_importance_del),
  stringsAsFactors = FALSE
)

## ---- 2. Median chromosome-specific importance per feature ----

chr_median_ampl <- plot_df_ampl %>%
  transmute(feature = as.character(feature), chr_median_importance = median_importance)

chr_median_del <- plot_df_del %>%
  transmute(feature = as.character(feature), chr_median_importance = median_importance)

## ---- 3. One row per feature: WG importance vs median chromosome importance ----

comparison_ampl <- chr_median_ampl %>%
  left_join(wg_df_ampl, by = "feature") %>%
  mutate(
    log10_WG          = log10(WG_importance + 1),
    log10_chr_median  = log10(chr_median_importance + 1)
  )

comparison_del <- chr_median_del %>%
  left_join(wg_df_del, by = "feature") %>%
  mutate(
    log10_WG          = log10(WG_importance + 1),
    log10_chr_median  = log10(chr_median_importance + 1)
  )

if (anyNA(comparison_ampl$WG_importance)) {
  warning("Amplification: features with no WG match -> ",
          paste(comparison_ampl$feature[is.na(comparison_ampl$WG_importance)], collapse = ", "))
}
if (anyNA(comparison_del$WG_importance)) {
  warning("Deletion: features with no WG match -> ",
          paste(comparison_del$feature[is.na(comparison_del$WG_importance)], collapse = ", "))
}

## ---- 4. Correlations -- n = number of features, not features x chromosomes ----

pearson_cor_ampl  <- cor(comparison_ampl$log10_WG, comparison_ampl$log10_chr_median,
                         method = "pearson", use = "complete.obs")
spearman_cor_ampl <- cor(comparison_ampl$WG_importance, comparison_ampl$chr_median_importance,
                         method = "spearman", use = "complete.obs")

pearson_cor_del  <- cor(comparison_del$log10_WG, comparison_del$log10_chr_median,
                        method = "pearson", use = "complete.obs")
spearman_cor_del <- cor(comparison_del$WG_importance, comparison_del$chr_median_importance,
                        method = "spearman", use = "complete.obs")

cat("Amplification — Pearson (log10+1):", round(pearson_cor_ampl, 3),
    "| Spearman:", round(spearman_cor_ampl, 3),
    "| n features:", nrow(comparison_ampl), "\n")
cat("Deletion — Pearson (log10+1):", round(pearson_cor_del, 3),
    "| Spearman:", round(spearman_cor_del, 3),
    "| n features:", nrow(comparison_del), "\n")

## ---- 5. Main plot: one labeled point per feature ----

comparison_scatter_ampl <- ggplot(comparison_ampl, aes(x = log10_WG, y = log10_chr_median)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, alpha = 0.85, color = "steelblue") +
  ggrepel::geom_text_repel(aes(label = feature), size = 3, max.overlaps = Inf) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 3.5,
           label = paste0("Pearson r (log10) = ", round(pearson_cor_ampl, 3),
                          "\nSpearman rho = ", round(spearman_cor_ampl, 3),
                          "\nn = ", nrow(comparison_ampl), " features")) +
  theme_bw(base_size = 12) +
  labs(title = "Amplification: WG importance vs median chromosome-specific importance",
       subtitle = "One point per feature (median taken across chromosomes)",
       x = "log10(WG importance + 1)",
       y = "log10(Median chromosome-specific importance + 1)")

comparison_scatter_del <- ggplot(comparison_del, aes(x = log10_WG, y = log10_chr_median)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, alpha = 0.85, color = "indianred") +
  ggrepel::geom_text_repel(aes(label = feature), size = 3, max.overlaps = Inf) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 3.5,
           label = paste0("Pearson r (log10) = ", round(pearson_cor_del, 3),
                          "\nSpearman rho = ", round(spearman_cor_del, 3),
                          "\nn = ", nrow(comparison_del), " features")) +
  theme_bw(base_size = 12) +
  labs(title = "Deletion: WG importance vs median chromosome-specific importance",
       subtitle = "One point per feature (median taken across chromosomes)",
       x = "log10(WG importance + 1)",
       y = "log10(Median chromosome-specific importance + 1)")

print(comparison_scatter_ampl)
print(comparison_scatter_del)

## ---- 6. FEATURE NAME -> PUBLICATION NAME LOOKUP ----

feature_lookup <- data.frame(
  feature = c(
    "mean.GC.content","total_n_partners.trans","total_n_PPIs.trans",
    "total_n_ohnologs.mmpaper_trans","total_n_paralogs_trans","dist.to.closest.FGS",
    "Length_Counts.E1","Length_Counts.E10","Length_Counts.E11","Length_Counts.E12",
    "Length_Counts.E13","Length_Counts.E14","Length_Counts.E15","Length_Counts.E16",
    "Length_Counts.E17","Length_Counts.E18","Length_Counts.E19","Length_Counts.E2",
    "Length_Counts.E20","Length_Counts.E21","Length_Counts.E22","Length_Counts.E23",
    "Length_Counts.E24","Length_Counts.E25","Length_Counts.E3","Length_Counts.E4",
    "Length_Counts.E5","Length_Counts.E6","Length_Counts.E7","Length_Counts.E8",
    "Length_Counts.E9","all.int.trans","genes.bin","partners.trans","density.OG",
    "density.TSG","dist.to.closest.OG","dist.to.closest.TSG","mutations_norm",
    "distance.to.centromere","distance.to.telomere","Ess.distance_pancancer",
    "ESSscore_pancancer","HAPLOscore_pancancer","Density.complex.proteins",
    "Density.Ohnologs","Centromere_Length","Chromosome_Length","Centromere_type"
  ),
  publication_name = c(
    "GC content","Number complex partner","Number PPIs","Number ohnologs",
    "Number paralogs","Distance to FGS","Active TSS",
    "Transcribed 5' preferential and Enh","Transcribed 3' preferential and Enh",
    "Transcribed and Weak Enhancer","Active Enhancer 1","Active Enhancer 2",
    "Active Enhancer Flank","Weak Enhancer 1","Weak Enhancer 2",
    "Primary H3K27ac possible Enhancer","Primary DNase","Promoter Upstream TSS",
    "ZNF genes & repeats","Heterochromatin","Poised Promoter","Bivalent Promoter",
    "Repressed Polycomb","Quiescent/Low","Promoter Downstream TSS 1",
    "Promoter Downstream TSS 2","Transcribed - 5' preferential",
    "Strong transcription","Transcribed - 3' preferential","Weak transcription",
    "Transcribed & regulatory (Prom/Enh)","Normal tissue all interactions expression",
    "Normal tissue expression","Normal tissue partner expression","OG density",
    "TSG density","Distance to OG","Distance to TSG","Mutation score",
    "Distance to centromere","Distance to telomere","Distance to essential gene",
    "Essentiality score","Haploinsufficiency score","Complex protein density",
    "Onholog density","Centromere length","Chromosome length","Centromere type"
  ),
  stringsAsFactors = FALSE
)

## ---- 7. Faceted diagnostic plots, using publication names for labels ----
## Same idea as before, but instead of collapsing to the median you see each
## chromosome's own value against WG. Point labels now use the publication
## name from the table above; any feature not found in the lookup falls
## back to its raw feature name so nothing is dropped silently.

facet_data_ampl <- rankings_df_ampl %>%
  mutate(feature = as.character(feature)) %>%
  left_join(wg_df_ampl, by = "feature") %>%
  left_join(feature_lookup, by = "feature") %>%
  mutate(
    publication_name = ifelse(is.na(publication_name), feature, publication_name),
    chromosome = str_extract(model, "(?<=_)[0-9XYxy]+$"),  # adjust regex if your model names differ
    log10_WG   = log10(WG_importance + 1),
    log10_chr  = log10(importance + 1)
  )

facet_data_del <- rankings_df_del %>%
  mutate(feature = as.character(feature)) %>%
  left_join(wg_df_del, by = "feature") %>%
  left_join(feature_lookup, by = "feature") %>%
  mutate(
    publication_name = ifelse(is.na(publication_name), feature, publication_name),
    chromosome = str_extract(model, "(?<=_)[0-9XYxy]+$"),
    log10_WG   = log10(WG_importance + 1),
    log10_chr  = log10(importance + 1)
  )

chr_levels_ampl <- str_sort(unique(facet_data_ampl$chromosome), numeric = TRUE)
chr_levels_del  <- str_sort(unique(facet_data_del$chromosome),  numeric = TRUE)
facet_data_ampl$chromosome <- factor(facet_data_ampl$chromosome, levels = chr_levels_ampl)
facet_data_del$chromosome  <- factor(facet_data_del$chromosome,  levels = chr_levels_del)

facet_cor_ampl <- facet_data_ampl %>%
  group_by(chromosome) %>%
  summarise(spearman = cor(WG_importance, importance, method = "spearman", use = "complete.obs"),
            .groups = "drop") %>%
  mutate(label = paste0("rho = ", round(spearman, 2)))

facet_cor_del <- facet_data_del %>%
  group_by(chromosome) %>%
  summarise(spearman = cor(WG_importance, importance, method = "spearman", use = "complete.obs"),
            .groups = "drop") %>%
  mutate(label = paste0("rho = ", round(spearman, 2)))

facet_plot_ampl <- ggplot(facet_data_ampl, aes(x = log10_WG, y = log10_chr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 1.6, alpha = 0.7, color = "steelblue") +
  ggrepel::geom_text_repel(aes(label = publication_name), size = 2, max.overlaps = 8) +
  geom_text(data = facet_cor_ampl, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.3, size = 2.6) +
  facet_wrap(~ chromosome, scales = "free") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 8)) +
  labs(title = "Amplification: WG vs per-chromosome importance, by chromosome",
       x = "log10(WG importance + 1)", y = "log10(Chromosome-specific importance + 1)")

facet_plot_del <- ggplot(facet_data_del, aes(x = log10_WG, y = log10_chr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 1.6, alpha = 0.7, color = "indianred") +
  ggrepel::geom_text_repel(aes(label = publication_name), size = 2, max.overlaps = 8) +
  geom_text(data = facet_cor_del, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.3, size = 2.6) +
  facet_wrap(~ chromosome, scales = "free") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 8)) +
  labs(title = "Deletion: WG vs per-chromosome importance, by chromosome",
       x = "log10(WG importance + 1)", y = "log10(Chromosome-specific importance + 1)")

print(facet_plot_ampl)
print(facet_plot_del)

## ---- 8. Save everything alongside the rest of the comparison plots ----

ggsave(file.path(output_dir, "scatterplot_ampl_WG_vs_chrMedian.pdf"), comparison_scatter_ampl, width = 8, height = 7)
ggsave(file.path(output_dir, "scatterplot_del_WG_vs_chrMedian.pdf"),  comparison_scatter_del,  width = 8, height = 7)
ggsave(file.path(output_dir, "facet_scatterplot_ampl_WG_vs_chr_pubnames.pdf"), facet_plot_ampl, width = 18, height = 14)
ggsave(file.path(output_dir, "facet_scatterplot_del_WG_vs_chr_pubnames.pdf"),  facet_plot_del,  width = 18, height = 14)

## ---- 9. Histogram of per-chromosome correlations (Amplification only) ----

hist_cor_ampl <- ggplot(facet_cor_ampl, aes(x = spearman)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "steelblue", color = "white", alpha = 0.85) +
  geom_vline(xintercept = spearman_cor_ampl, linetype = "dashed", color = "black") +
  annotate("text", x = spearman_cor_ampl, y = Inf, hjust = -0.1, vjust = 1.5, size = 3.2,
           label = paste0("median-based rho = ", round(spearman_cor_ampl, 3))) +
  theme_bw(base_size = 12) +
  labs(title = "Amplification: distribution of per-chromosome Spearman correlations",
       subtitle = "Each bar = one chromosome's WG vs chromosome-specific importance correlation",
       x = "Spearman rho (per chromosome)",
       y = "Count of chromosomes")

print(hist_cor_ampl)