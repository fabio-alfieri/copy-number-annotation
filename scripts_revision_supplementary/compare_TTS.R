pattern_shap <- "Mid-length"

##### NEW: point directly at the TTS subfolder (files are Ampl only, no more type subfolders to split on)
model.outputs_SHAP <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/SHAP_and_featurematrix/TTS", 
                                 pattern = pattern_shap, recursive = F, full.names = TRUE)

##### NEW: model name = pattern + cts_label, extracted from filename
# e.g. "SHAP_and_FeatureMatrix_Mid-length_STAD_STAD_Ampl.rds" -> "Mid-length_STAD_STAD"
models <- gsub(pattern = "_Ampl\\.rds$", replacement = "", x = basename(model.outputs_SHAP))
models <- gsub(pattern = "SHAP_and_FeatureMatrix_", replacement = "", x = models)

############ AMPLIFICATION-ONLY ANALYSIS (per cancer-type-pair model)

rankings_ampl <- lapply(X = seq_along(model.outputs_SHAP), 
                        FUN = function(idx){
                          shap_and_features <- readRDS(model.outputs_SHAP[idx])
                          shap_ampl <- shap_and_features$models.shap.df[[1]]
                          imp_ampl <- apply(abs(shap_ampl[, -1]), 2, sum)
                          ranking_ampl <- data.frame(feature = names(imp_ampl), importance = imp_ampl, rank = rank(-imp_ampl, ties.method = "min"))
                          ranking_ampl <- ranking_ampl[order(ranking_ampl$rank), ]
                          ranking_ampl$model <- models[idx]
                          return(ranking_ampl)
                        }
)

library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(grid)          # required by save_pheatmap_pair (grid.newpage, pushViewport, etc.)
library(ComplexHeatmap) # required by save_pheatmap_pair (draw())

rankings_df_ampl <- bind_rows(rankings_ampl)

feature_order_ampl <- rankings_df_ampl %>% group_by(feature) %>% summarise(mean_rank = mean(rank, na.rm = TRUE)) %>% arrange(mean_rank) %>% pull(feature)

rankings_df_ampl$feature <- factor(rankings_df_ampl$feature, levels = feature_order_ampl)

rank_matrix_ampl <- rankings_df_ampl %>% dplyr::select(feature, model, rank) %>% tidyr::pivot_wider(names_from = model, values_from = rank) %>% tibble::column_to_rownames("feature") %>% as.matrix()

sd_rank_ampl <- apply(X = rank_matrix_ampl, MARGIN = 1, FUN = sd)

rankings_df_ampl$sd_rank <- sd_rank_ampl[as.character(rankings_df_ampl$feature)]

plot_df_ampl <- rankings_df_ampl %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(importance, na.rm = TRUE),
    median_importance = median(importance, na.rm = TRUE),
    sd_rank = first(sd_rank)
  )

##### NEW: single-type df_combined since Deletion no longer exists
df_combined <- data.frame(Feature = names(sd_rank_ampl), SD = sd_rank_ampl, Type = "Amplification")

cor_ampl <- cor.test(plot_df_ampl$median_importance, plot_df_ampl$sd_rank, method = "pearson")
r_ampl <- round(cor_ampl$estimate, 3)
p_ampl <- signif(cor_ampl$p.value, 3)

plot_df_ampl_lm <- plot_df_ampl %>% filter(median_importance > 0)

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

rank_heatmap_ampl_chr <- pheatmap(
  rank_matrix_ampl,
  cluster_rows = FALSE, cluster_cols = TRUE, scale = "none",
  color = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))(40)
)

sd_plot_chr <- ggplot(df_combined, aes(x = reorder(Feature, SD), y = SD, fill = Type)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Variability of Feature Importance Ranks", x = "Feature", y = "SD of Rank") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

scatterplot_ampl_chr <- ggplot(plot_df_ampl, aes(x = median_importance, y = sd_rank)) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(data = plot_df_ampl_lm, aes(label = feature), size = 3, max.overlaps = Inf) +
  geom_smooth(data = plot_df_ampl_lm, method = "lm", formula = y ~ x, se = TRUE, color = "red") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("Pearson r = ", r_ampl, "\np = ", p_ampl)) +
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

##### NEW: save all plots
out_dir <- "/home/ieo7429/Scrivania/results_regressor_gab/SHAP_and_featurematrix/TTS/plots/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(filename = paste0(out_dir, "ranking_plot_ampl_", pattern_shap, ".pdf"),
       plot = ranking_plot_ampl_chr, width = 10, height = 8)

# pheatmap objects aren't ggplot objects; wrap the gtable for ggsave, or save directly via pdf()
ggsave(filename = paste0(out_dir, "rank_heatmap_ampl_", pattern_shap, ".pdf"),
       plot = rank_heatmap_ampl_chr$gtable, width = 8, height = 8)

ggsave(filename = paste0(out_dir, "sd_plot_", pattern_shap, ".pdf"),
       plot = sd_plot_chr, width = 8, height = 8)

ggsave(filename = paste0(out_dir, "scatterplot_ampl_", pattern_shap, ".pdf"),
       plot = scatterplot_ampl_chr, width = 7, height = 6)

ggsave(filename = paste0(out_dir, "pca_plot_ampl_", pattern_shap, ".pdf"),
       plot = pca_plot_ampl, width = 7, height = 6)

###############################################################################
## COMPARISON: TUMOR-TYPE-SPECIFIC vs WHOLE-GENOME SHAP IMPORTANCE (Ampl only)
##
## Mirrors the "chromosome-specific vs WG" comparison from the chromosome
## script, adapted for tumor-type-pair models. As before: correlate WG
## importance against the MEDIAN importance across tumor-type models per
## feature (plot_df_ampl$median_importance, already computed above) rather
## than against every tumor-type value individually, to avoid pseudo-
## replication in the correlation.
###############################################################################

library(stringr)

## ---- 1. Load the WG reference file (same one used in the chromosome-specific
##         script: first file matching pattern_shap in the ORIGINAL, non-TTS path) ----

model.outputs_SHAP_WG_source <- list.files(
  path = "/home/ieo7429/Scrivania/results_regressor_gab/SHAP_and_featurematrix",
  pattern = pattern_shap, recursive = TRUE, full.names = TRUE
)

shap_and_features_WG <- readRDS(model.outputs_SHAP_WG_source[2])
shap_ampl_WG <- shap_and_features_WG$models.shap.df[[1]]

## ---- 2. Whole-genome importance per feature ----

feature_importance_ampl <- colSums(abs(shap_ampl_WG[, -1]), na.rm = TRUE)

wg_df_ampl <- data.frame(
  feature       = names(feature_importance_ampl),
  WG_importance = as.numeric(feature_importance_ampl),
  stringsAsFactors = FALSE
)

## ---- 3. Median tumor-type-specific importance per feature ----
## (already computed above as plot_df_ampl$median_importance)

tt_median_ampl <- plot_df_ampl %>%
  transmute(feature = as.character(feature), tt_median_importance = median_importance)

## ---- 4. One row per feature: WG importance vs median tumor-type importance ----

comparison_ampl <- tt_median_ampl %>%
  left_join(wg_df_ampl, by = "feature") %>%
  mutate(
    log10_WG        = log10(WG_importance + 1),
    log10_tt_median = log10(tt_median_importance + 1)
  )

if (anyNA(comparison_ampl$WG_importance)) {
  warning("Amplification: features with no WG match -> ",
          paste(comparison_ampl$feature[is.na(comparison_ampl$WG_importance)], collapse = ", "))
}

## ---- 5. Correlations -- n = number of features, not features x tumor-type models ----

pearson_cor_ampl  <- cor(comparison_ampl$log10_WG, comparison_ampl$log10_tt_median,
                         method = "pearson", use = "complete.obs")
spearman_cor_ampl <- cor(comparison_ampl$WG_importance, comparison_ampl$tt_median_importance,
                         method = "spearman", use = "complete.obs")

cat("Amplification (tumor-type) — Pearson (log10+1):", round(pearson_cor_ampl, 3),
    "| Spearman:", round(spearman_cor_ampl, 3),
    "| n features:", nrow(comparison_ampl), "\n")

## ---- 6. Main plot: one labeled point per feature ----

comparison_scatter_ampl <- ggplot(comparison_ampl, aes(x = log10_WG, y = log10_tt_median)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, alpha = 0.85, color = "steelblue") +
  ggrepel::geom_text_repel(aes(label = feature), size = 3, max.overlaps = Inf) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 3.5,
           label = paste0("Pearson r (log10) = ", round(pearson_cor_ampl, 3),
                          "\nSpearman rho = ", round(spearman_cor_ampl, 3),
                          "\nn = ", nrow(comparison_ampl), " features")) +
  theme_bw(base_size = 12) +
  labs(title = "Amplification: WG importance vs median tumor-type-specific importance",
       subtitle = "One point per feature (median taken across tumor-type models)",
       x = "log10(WG importance + 1)",
       y = "log10(Median tumor-type-specific importance + 1)")

print(comparison_scatter_ampl)

## ---- 7. Optional diagnostic: every tumor-type model individually, faceted ----
## `model` is already the tumor-type-pair label here, so no regex extraction
## step is needed (unlike the chromosome script, which had to pull the
## chromosome out of the model name).

facet_data_ampl <- rankings_df_ampl %>%
  mutate(feature = as.character(feature)) %>%
  left_join(wg_df_ampl, by = "feature") %>%
  mutate(
    log10_WG  = log10(WG_importance + 1),
    log10_tt  = log10(importance + 1)
  )

facet_cor_ampl <- facet_data_ampl %>%
  group_by(model) %>%
  summarise(spearman = cor(WG_importance, importance, method = "spearman", use = "complete.obs"),
            .groups = "drop") %>%
  mutate(label = paste0("rho = ", round(spearman, 2)))

facet_plot_ampl <- ggplot(facet_data_ampl, aes(x = log10_WG, y = log10_tt)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 1.6, alpha = 0.7, color = "steelblue") +
  ggrepel::geom_text_repel(aes(label = feature), size = 2, max.overlaps = 8) +
  geom_text(data = facet_cor_ampl, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.3, size = 2.6) +
  facet_wrap(~ model, scales = "free") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 8)) +
  labs(title = "Amplification: WG vs per-tumor-type importance, by tumor-type model",
       x = "log10(WG importance + 1)", y = "log10(Tumor-type-specific importance + 1)")

print(facet_plot_ampl)

## ---- 8. Histogram of per-tumor-type correlations ----
## Same idea as the chromosome-script histogram: one Spearman rho per
## tumor-type model (from facet_cor_ampl), showing how much agreement
## with WG importance varies across tumor-type pairs.

hist_cor_ampl <- ggplot(facet_cor_ampl, aes(x = spearman)) +
  geom_histogram(bins = 5, boundary = 0, fill = "steelblue", color = "white", alpha = 0.85) +
  geom_vline(xintercept = spearman_cor_ampl, linetype = "dashed", color = "black") +
  annotate("text", x = spearman_cor_ampl, y = Inf, hjust = -0.1, vjust = 1.5, size = 3.2,
           label = paste0("median-based rho = ", round(spearman_cor_ampl, 3))) +
  theme_bw(base_size = 12) +
  labs(title = "Amplification: distribution of per-tumor-type Spearman correlations",
       subtitle = "Each bar = one tumor-type model's WG vs tumor-type-specific importance correlation",
       x = "Spearman rho (per tumor-type model)",
       y = "Count of tumor-type models")

print(hist_cor_ampl)

## ---- 9. Save everything ----

ggsave(filename = paste0(out_dir, "scatterplot_ampl_WG_vs_ttMedian_", pattern_shap, ".pdf"),
       plot = comparison_scatter_ampl, width = 8, height = 7)

ggsave(filename = paste0(out_dir, "facet_scatterplot_ampl_WG_vs_tt_", pattern_shap, ".pdf"),
       plot = facet_plot_ampl, width = 18, height = 14)

ggsave(filename = paste0(out_dir, "histogram_ampl_tt_correlations_", pattern_shap, ".pdf"),
       plot = hist_cor_ampl, width = 7, height = 5)
