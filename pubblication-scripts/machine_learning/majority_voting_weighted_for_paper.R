# compute correlations across models
rm(list=ls()) # clean variables
gc(full=T) # garbage collector

#######################  import libraries #######################

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
library(ggplot2)
library(ggtext)
library(dplyr)

#######################

#######################  set working directory ####################### 

setwd("/home/ieo7429/Scrivania/")

#######################

#######################  import needed tables ####################### 

avg_performances <- read.table(file = "results_regressor_gab/avg_performances_table.tsv", header = T, sep = "\t")
correlations <- read.table(file = "results_regressor_gab/cor_final.tsv", header = T, sep = "\t")

avg_performances[,c("model", "ampl_del")] <- do.call(rbind, strsplit(rownames(avg_performances), split = "_")) # take class of CNA and model type
avg_performances[avg_performances$model == "NoCluster",]$model <- "no_cluster" # fix NoCluster
rownames(avg_performances) <- NULL # remove rownames for readibility

correlations_with_avg_performances <- merge(x = correlations, # join tables to include performances
                                            y = avg_performances, 
                                            on = c("model", "ampl_del"))

# remove useless columns
correlations_with_avg_performances$avg_rmse <- NULL 
correlations_with_avg_performances$avg_pearson <- NULL
correlations_with_avg_performances$avg_spearman <- NULL
correlations_with_avg_performances$signif <- NULL
correlations_with_avg_performances$pval.P_fdr <- NULL

# add sign of correlation (pearson), and weighted correlations
correlations_with_avg_performances$cor.sign <- ifelse(sign(correlations_with_avg_performances$estimate.P) == 1, "Positive", "Negative")
correlations_with_avg_performances$weighted_cor.P <- correlations_with_avg_performances$estimate.P * correlations_with_avg_performances$avg_r2
correlations_with_avg_performances$weighted_cor.S <- correlations_with_avg_performances$estimate.S * correlations_with_avg_performances$avg_r2

# tot contributors are the number of model (for each feature, for each type) contributing with a correlation
# so we take the total number for feature x and model y, meaning that that is the maximum
# for R2 we sum, summing the "cumulative variance" explained
tot_contributors <- correlations_with_avg_performances %>%
  group_by(Column, ampl_del) %>%
  summarise(num_of_contributors = n(),
            sum_of_r2_contributors = sum(avg_r2)) %>% 
  ungroup() 

# tot_contributors <- tot_contributors[tot_contributors$Column %in% tot_contributors[tot_contributors$num_of_contributors >= 4,]$Column,]

# merge 
correlations_with_avg_performances_and_tot_contributors <- merge(x = correlations_with_avg_performances, 
                                                                 y = tot_contributors, 
                                                                 on = c("Column", "ampl_del"))

# features <- c('')

corP_imp_features <- correlations_with_avg_performances_and_tot_contributors %>% 
  filter(Column %in% important_features) %>%
  ggplot(aes(y = Column, x = weighted_cor.P, fill = model)) +
  geom_bar(stat = 'identity'#, position = "dodge"
  ) +
  facet_wrap(~ampl_del)  +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

corP_selection <- correlations_with_avg_performances_and_tot_contributors %>% 
  filter(feature_type == "Selection") %>%
  ggplot(aes(y = Column, x = weighted_cor.P, fill = model)) +
  geom_bar(stat = 'identity'#, position = "dodge"
           ) +
  facet_wrap(~ampl_del)  +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
corP_occurrence <- correlations_with_avg_performances_and_tot_contributors %>% 
  filter(feature_type == "Occurrence") %>%
  ggplot(aes(y = Column, x = weighted_cor.P, fill = model)) +
  geom_bar(stat = 'identity'#, position = "dodge"
  ) +
  facet_wrap(~ampl_del)  +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
corS_selection <- correlations_with_avg_performances_and_tot_contributors %>% 
  filter(feature_type == "Selection"
  ) %>%
  ggplot(aes(y = Column, x = weighted_cor.S, fill = model)) +
  geom_bar(stat = 'identity'#, position = "dodge"
  ) +
  facet_wrap(~ampl_del)  +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

pdf('dev/imgs/ampl_del_weigthed_corP_selectionFeatures.pdf', width = 12, height = 8)
corP_imp_features
dev.off()

pdf('dev/imgs/ampl_del_weigthed_corP_selectionFeatures.pdf', width = 12, height = 4)
corP_selection
dev.off()
pdf('dev/imgs/ampl_del_weigthed_corP_occurrenceFeatures.pdf', width = 12, height = 8)
corP_occurrence
dev.off()
pdf('dev/imgs/ampl_del_weigthed_corS_selectionFeatures.pdf', width = 12, height = 4)
corS_selection
dev.off()




# define either the percentage of models (for each feature, for each type) contributing with a certain kind of correlation
# in case of R2 define as the percentage of variance explained by the models for a certain kind of correlation over the total
# let's call R2 as "trust" we have in the models. In total contributors we define the maximum trust we can give them (max 5)
# we take the total trust and then we compute the percentage divided by positive / negative correlation
voting_output <- correlations_with_avg_performances_and_tot_contributors %>% 
  group_by(Column, ampl_del, cor.sign) %>% 
  summarise(perc_of_voting = n() / dplyr::first(num_of_contributors) * 100, 
            perc_of_r2 = sum(avg_r2) / dplyr::first(sum_of_r2_contributors) * 100
            )

voting_output <- merge(unique(correlations_with_avg_performances_and_tot_contributors[,c("Column", "feature_type")]), 
                      voting_output, 
                      on = "Column")

# the two metrics correlate
cor.test(voting_output$perc_of_voting, voting_output$perc_of_r2)
plot(voting_output$perc_of_voting, voting_output$perc_of_r2)

# final verdict is taking the max
final_verdict <- voting_output %>% group_by(Column, ampl_del) %>% 
  slice_max(perc_of_voting, n = 1, with_ties = F) %>% 
  ungroup() %>% 
  dplyr::select(Column, feature_type, ampl_del, cor.sign)%>% 
  pivot_wider(id_cols = c(Column, feature_type), names_from = ampl_del, values_from = cor.sign)

# final_verdict <- merge(x = final_verdict, y = unique(correlations[,c("Column", "feature_type")]), on = "Column")

colnames(final_verdict) <- c("feature", "feature_type", "ampl", "del")

write.table(x = final_verdict, file = "results_regressor_gab/final_verdict.tsv", sep = "\t", 
             quote = F, row.names = F, col.names = T)

feats <- c("dist.to.closest.TSG", "dist.to.closest.OG",
           "all.int.trans", "partners.trans", 
           "total_n_partners.trans", "total_n_PPIs.trans", 
           "mutations_norm", "Ess.distance_pancancer",
           "total_n_paralogs_trans", "genes.bin",
           "total_n_ohnologs.mmpaper_trans"
)


heatmap_data <- final_verdict %>%
  filter(feature %in% feats) %>%
  pivot_longer(cols = c(ampl, del), names_to = "Type", values_to = "Sign")

ggplot(heatmap_data, aes(x = Type, y = feature, fill = Sign)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c("Positive" = "steelblue", "Negative" = "tomato", "NA" = "grey80"),
    na.value = "grey80"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_markdown(size = 8),  # use element_markdown for bold
    panel.grid = element_blank()
  )  +  # apply the conditional labels
  labs(
    x = "",
    y = "Column",
    fill = "Sign",
    title = "Amplification vs Deletion Heatmap"
  )


######################

tmp <- voting_output %>% group_by(Column, ampl_del) %>% 
  slice_max(perc_of_voting, n = 1, with_ties = F) %>% 
  ungroup() %>%
  mutate(
    perc_of_voting = if_else(cor.sign == "Negative", -perc_of_voting, perc_of_voting),
    perc_of_r2     = if_else(cor.sign == "Negative", -perc_of_r2, perc_of_r2)
  )  %>% 
  pivot_wider(id_cols = Column, names_from = ampl_del, values_from = c(perc_of_voting, perc_of_r2))

colnames(tmp) <- c("feature", "ampl_model_voting", "del_model_voting", "ampl_model_r2", "del_model_r2")

# sel_feats <- c("Ess.distance_pancancer", "dist.to.closest.TSG", "dist.to.closest.OG", "all.int.trans", "genes.bin", "mutation_score")

heatmap_data <- tmp %>%
  filter(feature %in% feats) %>%
  pivot_longer(cols = c(ampl_model_voting, del_model_voting, 
                        # ampl_model_r2, del_model_r2
                        ), names_to = "Type", values_to = "Sign") %>%
  dplyr::select(feature, Type, Sign)

ggplot(heatmap_data, aes(x = Type, y = feature, fill = Sign)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "tomato",      # negative values
    mid = "grey80",      # zero
    high = "steelblue",  # positive values
    midpoint = 0,        # center the scale at 0
    na.value = "grey80"  # for NA values
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  labs(
    x = "",
    y = "Column",
    fill = "Sign",
    title = "Amplification vs Deletion Heatmap"
  )

#######################

library(pheatmap)

mat <- heatmap_data %>%
  tidyr::pivot_wider(names_from = Type, values_from = Sign) %>%
  tibble::column_to_rownames("feature") %>%
  as.matrix()

pheatmap(
  mat,
  color = colorRampPalette(c("tomato", "grey80", "steelblue"))(100),
  clustering_method = "complete",   # hierarchical clustering
  cluster_rows = T,
  cluster_cols = F,
  main = "Amplification vs Deletion Heatmap"
)

