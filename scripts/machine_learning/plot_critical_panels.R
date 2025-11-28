library(ggplot2)

model_performances <- read.delim(file = "/home/ieo7429/Scrivania/results_regressor_gab/avg_performances_OO.tsv")

#rownames(model_performances) <- gsub(pattern = "_Length_Mid-length_ampl", replacement = "", x = rownames(model_performances))


PRO_idxs <- grep(pattern = "PRO_[0-9]+_[0-9]+_[0-9]+", x = rownames(model_performances))
LOCO_idxs <- grep(pattern = "LOCO_[0-9]+", x = rownames(model_performances))
LOCO_NSS_idxs <- grep(pattern = "LOCO_NSS_[1-9]+", x = rownames(model_performances))
LOFO_idxs <- grep(pattern = "^LOFO_(?!PROXY)", x = rownames(model_performances), perl = T)
PROXY_idxs <- grep(pattern = "^PROXY$",  x = rownames(model_performances))
LOFO_PROXY_idxs <- grep(pattern = "LOFO_PROXY_",  x = rownames(model_performances))
NSS_new_idxs <- grep(pattern = "NSS_new",  x = rownames(model_performances))
NO_OG_idxs <- grep(pattern = "NO_OG",  x = rownames(model_performances))
LOTO_idxs <- grep(pattern = "LOCTO_",  x = rownames(model_performances))
LAFO_PROXY_idxs <- grep(pattern = "LAFO_PROXY_",  x = rownames(model_performances))
OO_idxs <- grep(pattern = "OO",  x = rownames(model_performances))

PRO_performances <- model_performances[PRO_idxs,]
PRO_performances <- PRO_performances[order(unlist(as.numeric(lapply(X = rownames(PRO_performances), FUN = function(x){strsplit(x, split = "_")[[1]][[2]]}))), decreasing = F),]
PRO_performances <- PRO_performances[c(1:13),]

LOCO_performances <- model_performances[LOCO_idxs,]
LOCO_performances <- LOCO_performances[order(unlist(as.numeric(lapply(X = rownames(LOCO_performances), FUN = function(x){strsplit(x, split = "_")[[1]][[2]]}))), decreasing = F),]

LOTO_performances <- model_performances[LOTO_idxs,]

LOCO_NSS_performances <- model_performances[LOCO_NSS_idxs,]
LOFO_performances <- model_performances[LOFO_idxs,]
PROXY_performances <- model_performances[PROXY_idxs,]
LOFO_PROXY_performances <- model_performances[LOFO_PROXY_idxs,]
NSS_new_performances <- model_performances[NSS_new_idxs,]
NO_OG_perforamances <- model_performances[NO_OG_idxs,]
LAFO_PROXY_performances <- model_performances[LAFO_PROXY_idxs,]
OO_performances <- model_performances[OO_idxs,]


##################### PRO PLOTS

df_long <- PRO_performances %>%
  rownames_to_column("model") %>%
  mutate(model = factor(model, levels = rownames(PRO_performances))) %>%
  pivot_longer(
    cols = -model,
    names_to = "metric",
    values_to = "value"
  )

p1 <- ggplot(df_long[df_long$metric == "avg_r2",], aes(x = model, y = value, group = 1)) +
  geom_line() +
  geom_point(size = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Model",
    y = "Value",
    title = "Performance Metrics by Model"
  )

ggsave(filename = "dev/imgs/critical_panel/perf_decay.pdf", plot = p1, width = 20, height = 10)


##################### LOCO PLOTS

library(tidyverse)

df_stats_loco <- LOCO_performances %>%
  rownames_to_column("model") %>%
  pivot_longer(-model, names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  summarise(
    mean = mean(value),
    sd   = sd(value),
    ymin = mean - sd,
    ymax = mean + sd
  ) %>%
  mutate(
    y_min_limit = case_when(
      metric == "avg_r2" ~ 0,
      metric == "avg_rmse" ~ 0,
      metric %in% c("avg_pearson", "avg_spearman") ~ -1,
      TRUE ~ NA_real_
    ),
    y_max_limit = case_when(
      metric == "avg_r2" ~ 1,
      metric == "avg_rmse" ~ max(ymax) * 1.1,     # add a bit of headroom
      metric %in% c("avg_pearson", "avg_spearman") ~ 1,
      TRUE ~ NA_real_
    )
  )

p2 <- ggplot(df_stats_loco[df_stats_loco$metric == "avg_r2",], aes(x = "value", y = mean)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  geom_blank(aes(y = y_min_limit)) +
  geom_blank(aes(y = y_max_limit)) +
  
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = "",
    y = "Mean value",
    title = "Average LOCO Performance Metrics with Error Bars"
  )

ggsave(filename = "dev/imgs/critical_panel/avg_LOCO.pdf", plot = p2, width = 5, height = 8)

LOCO_performances$chr <- rownames(LOCO_performances)

avg_loco <- as.data.frame(t(df_stats_loco[,2]))
rownames(avg_loco) <- "Average LOCO"
colnames(avg_loco) <- df_stats_loco$metric
avg_loco

avg_loco$chr <- "Average_LOCO"
avg_r2_error_loco <- df_stats_loco[df_stats_loco$metric == "avg_r2",]

p5 <- ggplot(rbind(avg_loco, LOCO_performances), aes(x = factor(chr, levels = chr), y = avg_r2, group = 1)) + 
  geom_bar(stat = "identity", fill = "steelblue") +  # Bar plot for avg_r2
  geom_point(size = 2) +  # Points on top of the bars
  geom_errorbar(
    aes(
      ymin = avg_r2 - avg_r2_error_loco$sd, 
      ymax = avg_r2 + avg_r2_error_loco$sd
    ),
    data = subset(rbind(avg_loco, LOCO_performances), chr == "Average_LOCO"),  # Only apply error bars to Average_LOCO
    width = 0.2, 
    color = "black"
  ) +
  theme_bw() +  # Use theme_bw for a clean plot
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  labs(
    x = "LOCO",
    y = "Average R2",
    title = "Average R2 by LOCO Model"
  )

ggsave(filename = "dev/imgs/critical_panel/all_LOCOs.pdf", plot = p5, width = 20, height = 10)


##################### LOTO PLOTS

library(tidyverse)

df_stats_loto <- LOTO_performances %>%
  rownames_to_column("model") %>%
  pivot_longer(-model, names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  summarise(
    mean = mean(value),
    sd   = sd(value),
    ymin = mean - sd,
    ymax = mean + sd
  ) %>%
  mutate(
    y_min_limit = case_when(
      metric == "avg_r2" ~ 0,
      metric == "avg_rmse" ~ 0,
      metric %in% c("avg_pearson", "avg_spearman") ~ -1,
      TRUE ~ NA_real_
    ),
    y_max_limit = case_when(
      metric == "avg_r2" ~ 1,
      metric == "avg_rmse" ~ max(ymax) * 1.1,     # add a bit of headroom
      metric %in% c("avg_pearson", "avg_spearman") ~ 1,
      TRUE ~ NA_real_
    )
  )

p3 <- ggplot(df_stats_loto[df_stats_loto$metric == "avg_r2",], aes(x = "value", y = mean)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  geom_blank(aes(y = y_min_limit)) +
  geom_blank(aes(y = y_max_limit)) +
  
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = "",
    y = "Mean value",
    title = "Average LOTO Performance Metrics with Error Bars"
  )

ggsave(filename = "dev/imgs/critical_panel/avg_LOTO.pdf", plot = p3, width = 5, height = 8)

LOTO_performances$type <- unlist(lapply(X = rownames(LOTO_performances), FUN = function(x){strsplit(x, split = "_")[[1]][[2]]}))

avg_loto <- as.data.frame(t(df_stats_loto[,2]))
rownames(avg_loto) <- "Average LOTO"
colnames(avg_loto) <- df_stats_loto$metric
avg_loto

avg_loto$type <- "Average_LOTO"
avg_r2_error_loto <- df_stats_loto[df_stats_loto$metric == "avg_r2",]

p6 <- ggplot(rbind(avg_loto, LOTO_performances), aes(x = factor(type, levels = type), y = avg_r2, group = 1)) + 
  geom_bar(stat = "identity", fill = "steelblue") +  # Bar plot for avg_r2
  geom_point(size = 2) +  # Points on top of the bars
  geom_errorbar(
    aes(
      ymin = avg_r2 - avg_r2_error_loto$sd, 
      ymax = avg_r2 + avg_r2_error_loto$sd
    ),
    data = subset(rbind(avg_loto, LOTO_performances), type == "Average_LOTO"),  # Only apply error bars to Average_LOCO
    width = 0.2, 
    color = "black"
  ) +
  theme_bw() +  # Use theme_bw for a clean plot
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  labs(
    x = "LOTO",
    y = "Average R2",
    title = "Average R2 by LOTO Model"
  )

ggsave(filename = "dev/imgs/critical_panel/all_LOTOs.pdf", plot = p6, width = 20, height = 10)



##################### ALL TOGETHER

all_together_df <- rbind(model_performances[6,],
                         NSS_new_performances,
                         PROXY_performances,
                         LOFO_performances,
                         LOFO_PROXY_performances)

all_together_df <- all_together_df[c(1,2,3,4,9,5,10,6,11,7,12,8,13),]

plot_df <- all_together_df %>%
  rownames_to_column("model") %>%
  mutate(model = factor(model, levels = rownames(all_together_df))) %>%
  pivot_longer(
    cols = -model,
    names_to = "metric",
    values_to = "value"
  )

hlines <- plot_df[c(1),c(2,3)]

p4 <- ggplot(plot_df[plot_df$metric == "avg_r2",], aes(x = model, y = value)) +
  geom_col(fill = "steelblue")+
  geom_hline(data = hlines, aes(yintercept = value), linetype = "dotted") + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(size = 12)
  ) +
  labs(
    x = "Model",
    y = "Value",
    title = "Performance Metrics (All Models Together)"
  )

ggsave(filename = "dev/imgs/critical_panel/all_models_performances.pdf", plot = p4, width = 10, height = 7)

#################################

OO_performances_long <- OO_performances %>%
  tibble::rownames_to_column("model") %>%      # keep rownames
  dplyr::select(model, avg_r2) %>%                    # keep only R2 column
  pivot_longer(cols = avg_r2,
               names_to = "metric",
               values_to = "value")


df <- OO_performances_long[OO_performances_long$metric == "avg_r2", ]

p8 <- ggplot(df, aes(x = value, y = model)) +
  geom_col(fill = "steelblue") +
  theme_bw() +
  labs(
    x = "Average Test R2",
    y = "",
    title = "Model with only occurrence features"
  )

ggsave(filename = "dev/imgs/critical_panel/OO_performances.pdf", plot = p8, width = 10, height = 7)











