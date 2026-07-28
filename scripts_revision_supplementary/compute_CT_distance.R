rm(list=ls())
gc(full=T)

library(stringr)
library(parallel)
library(reshape2)
library(dplyr)
library(ggplot2)
library(Matrix)
library(caret)
library(caTools)
library(dplyr)
require(xgboost)
require(tidyverse)
library(ggrepel)

setwd("/mnt/ghost/fabiogokce/gokce")

# ── Output directory for all plots ────────────────────────────────────────────
output_dir <- "/home/ieo7429/Scrivania/results_regressor_gab/transfer_similarity_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Note: HIPPIE and Interactome INSIDER parts are alternative to each other

# Interactome INSIDER
#-----
output.files <- list.files(path = "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/MLTables_updated/InteractomeINSIDER",
                           full.names = TRUE, 
                           recursive = FALSE, 
                           pattern = "ampl_or_del_0.1Mbp_covThr_zero")

output.files <- output.files[1] # remove smoothed file
file <- output.files[1]
classS <- "Mid-length"

name <- "ampl_or_del_0.1Mbp_covThr_zero"
level <- "0.1Mbp"
bincov <- "covThr_zero"
folder <- "InteractomeINSIDER"
scale <- FALSE
cancer_type <- 'pancancer'


#-----

chromosome <- read.table(file = 'Codes/Codes-CNAs/MethodII/fabio/Data/chromosome.tsv', header = T)
centromeres <- read.table(file = 'Codes/Codes-CNAs/MethodII/fabio/Data/centomere.tsv', header = T)

tmp <- full_join(chromosome, centromeres, by = 'chr') %>% dplyr::select(chr, Centromere_Length, Chromosome_Length, Centromere_Type)
tmp$chr <- as.character(tmp$chr)

load(file)

ml_table_backup <- ML.Tables[["ampl_or_del"]][[classS]]
ml_table_backup$bin <- as.character(ml_table_backup$bin)
print(classS)

ml_table_backup$chr <- do.call(rbind, str_split(ml_table_backup$bin,'_'))[,1]
ml_table_backup <- left_join(ml_table_backup, tmp, by = 'chr')
ml_table_backup$chr <- NULL

# ── 0. Re-attach chr if dropped ───────────────────────────────────────────────
ml_table_backup <- ml_table_backup %>%
  mutate(chr = str_split_fixed(bin, "_", 2)[, 1])

# ── 1. Per-chromosome pairwise correlations between cancer types ───────────────
cor_per_chr <- ml_table_backup %>%
  dplyr::select(bin, chr, Type, ampl_score) %>%
  group_by(chr) %>%
  group_split() %>%
  lapply(function(df) {
    ch <- unique(df$chr)
    wide <- df %>% dplyr::select(bin, Type, ampl_score) %>% pivot_wider(names_from = Type, values_from = ampl_score) %>% column_to_rownames("bin")          # rows = bins (genomic positions aligned!)
    wide <- wide[rowSums(is.na(wide)) == 0, colSums(is.na(wide)) == 0]
    if (ncol(wide) < 2 || nrow(wide) < 3) return(NULL)   # need ≥3 bins to correlate
    cor(wide, method = "pearson") %>%    # swap "pearson" → "spearman" if preferred
      as.data.frame() %>%
      rownames_to_column("train_type") %>%
      pivot_longer(-train_type, names_to = "test_type", values_to = "cor") %>%
      mutate(chr = ch, n_bins = nrow(wide))
  }) %>%
  bind_rows()

# ── 2. Aggregate across chromosomes per pair ──────────────────────────────────
# Weighted mean: chromosomes with more bins contribute more
cor_agg <- cor_per_chr %>%
  filter(train_type != test_type) %>%
  group_by(train_type, test_type) %>%
  summarise(
    weighted_mean_cor = weighted.mean(cor, w = n_bins, na.rm = TRUE),
    mean_cor          = mean(cor, na.rm = TRUE),
    median_cor        = median(cor, na.rm = TRUE),
    n_chr             = n(),
    .groups = "drop"
  )

# ── 3. Load & parse performance table ─────────────────────────────────────────
perf <- read.table(
  "/home/ieo7429/Scrivania/results_regressor_gab/avg_performances_complete_TTS.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)

suffix <- paste0("_Length_", classS, "_ampl")

perf <- perf %>%
  rownames_to_column("pair_id") %>%
  mutate(
    pair_stripped = str_remove(pair_id, fixed(suffix)),
    train_type    = str_split_fixed(pair_stripped, "_", 2)[, 1],
    test_type     = str_split_fixed(pair_stripped, "_", 2)[, 2]
  )

# ── 4. Join ────────────────────────────────────────────────────────────────────
perf_cor <- perf %>%
  left_join(cor_agg, by = c("train_type", "test_type")) %>%
  filter(!is.na(weighted_mean_cor))

# ── 4.5. Collapse symmetric pairs (e.g. BRCA-OV and OV-BRCA) ──────────────────
perf_cor_avg <- perf_cor %>%
  mutate(
    type_a = pmin(train_type, test_type),
    type_b = pmax(train_type, test_type)
  ) %>%
  group_by(type_a, type_b) %>%
  summarise(
    avg_r2            = mean(avg_r2, na.rm = TRUE),
    weighted_mean_cor = mean(weighted_mean_cor, na.rm = TRUE),
    mean_cor          = mean(mean_cor, na.rm = TRUE),
    median_cor        = mean(median_cor, na.rm = TRUE),
    n_directions      = n(),     # 1 if only one direction was present, 2 if both
    .groups = "drop"
  ) %>%
  rename(train_type = type_a, test_type = type_b)

# ── 5. Correlation-of-correlations ────────────────────────────────────────────
res_spearman <- cor.test(perf_cor_avg$weighted_mean_cor, perf_cor_avg$avg_r2, method = "spearman", exact = F)
res_pearson  <- cor.test(perf_cor_avg$weighted_mean_cor, perf_cor_avg$avg_r2, method = "pearson")

cat(sprintf("Spearman rho = %.3f (p = %.4f) | Pearson r = %.3f (p = %.4f)  |  n = %d pairs\n",
            res_spearman$estimate, res_spearman$p.value,
            res_pearson$estimate,  res_pearson$p.value,
            nrow(perf_cor_avg)))

# ── 6. Plot: similarity vs R² ──────────────────────────────────────────────────
ann <- sprintf("Pearson r = %.3f, p = %.4f",
               res_spearman$estimate, res_spearman$p.value,
               res_pearson$estimate,  res_pearson$p.value)

p_similarity_vs_r2 <- ggplot(perf_cor_avg, aes(x = weighted_mean_cor, y = avg_r2)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1), se = TRUE, color = "grey60", linetype = "dashed") +
  geom_point(aes(color = (train_type == test_type)), size = 3) +
  geom_text_repel(aes(label = paste0(train_type, " - ", test_type)), size = 2, max.overlaps = 25) +
  scale_color_manual(values = c("FALSE" = "#3266ad", "TRUE" = "#d94f3d"), name = NULL) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = ann, size = 4, fontface = "italic") +
  labs(
    title = paste0("Amplification profile similarity (chromosome-wise) vs Average Test R²"),
    x     = "Weighted mean Spearman correlation across chromosomes",
    y     = expression(Average ~ R^2)
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

print(p_similarity_vs_r2)

ggsave(
  filename = file.path(output_dir, "similarity_vs_avgR2.pdf"),
  plot     = p_similarity_vs_r2,
  width    = 10,
  height   = 8
)

# ── 7. Average transfer performance by tumor type ────────────────────────────
transfer_by_type <- perf_cor_avg %>%
  pivot_longer(cols = c(train_type, test_type), names_to = "role", values_to = "cancer_type") %>%
  group_by(cancer_type) %>%
  summarise(
    avg_transfer_r2 = mean(avg_r2, na.rm = TRUE),
    n_pairs         = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(avg_transfer_r2))

print(transfer_by_type)

p_transfer_by_type <- ggplot(transfer_by_type, aes(x = reorder(cancer_type, avg_transfer_r2), y = avg_transfer_r2)) +
  geom_col(fill = "#3266ad") +
  geom_text(aes(label = sprintf("%.3f", avg_transfer_r2)), hjust = -0.15, size = 3) +
  coord_flip() +
  labs(
    title = paste0("Average cross-type transfer R² by tumor type  |  ", classS),
    x = NULL,
    y = expression(Average ~ transfer ~ R^2)
  ) +
  theme_bw(base_size = 13)

print(p_transfer_by_type)

ggsave(
  filename = file.path(output_dir, "avg_transfer_r2_by_type.pdf"),
  plot     = p_transfer_by_type,
  width    = 12,
  height   = 8
)

# ── Optional: heatmap of per-chromosome correlations for one pair ─────────────
# e.g. inspect KIRP vs ESCA across chromosomes
p_kirp_esca_perchr <- cor_per_chr %>%
  filter(train_type == "KIRP", test_type == "ESCA") %>%
  mutate(chr = factor(chr, levels = as.character(1:22))) %>%
  ggplot(aes(x = chr, y = cor, fill = cor)) +
  geom_col() +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850", midpoint = 0) +
  labs(title = "KIRP vs ESCA — per-chromosome ampl_score correlation",
       x = "Chromosome", y = "Pearson r") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_kirp_esca_perchr)

ggsave(
  filename = file.path(output_dir, "KIRP_vs_ESCA_perchr_correlation.pdf"),
  plot     = p_kirp_esca_perchr,
  width    = 10,
  height   = 6
)
