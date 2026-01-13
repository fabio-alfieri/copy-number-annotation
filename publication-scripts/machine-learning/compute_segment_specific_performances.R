rm(list=ls())
gc(full=T)

wd <- 'path/to/GitHub/copy-number-annotation/'

packages <- c(
  "stringr", "parallel", "reshape2", "dplyr",
  "ggplot2", "tidyr", "caret", "colorspace",
  "fitdistrplus", "disttree", "tidyverse",
  "colorspace", "ggsignif", "patchwork",
  "GenomicRanges"
)

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

#######################

#######################  load helpers ####################### 
wd.ml <- paste0(wd, "publication-scripts/machine-learning/")

source(paste0(wd.ml, "get_residuals.R"))
source(paste0(wd.ml, "compute_ratios.R"))
source(paste0(wd.ml, "aggregate_chromatin_states.R"))
source(paste0(wd.ml, "define_annotation_rules.R"))
source(paste0(wd.ml, "compare_ratio_based_top1.R"))
source(paste0(wd.ml, "change_distance.R"))
source(paste0(wd.ml, "define_annotation_rules_correlation_based.R"))
source(paste0(wd.ml, "do_pancancer_barplots_and_boxplots.R"))
source(paste0(wd.ml, "do_tumor_specific_barplots_and_boxplots.R"))
source(paste0(wd.ml, "do_do_tumor_specific_barplots_and_boxplots.R"))
source(paste0(wd.ml, "do_top1_annotation.R"))
source(paste0(wd.ml, "do_plot_grouped_by_chr_type.R"))
source(paste0(wd.ml, "do_plot_grouped_by_chr.R"))
source(paste0(wd.ml, "do_plot_grouped_by_type.R"))
source(paste0(wd.ml, "plotting_helpers.R"))
source(paste0(wd.ml, "landscape_plot_observed_prediction.R"))
source(paste0(wd.ml, "do_plot_residuals.R"))

####################### 

#######################  read model's output ####################### 

output <- list() # define output list

model_classes <- c("Mid-length", "Arm-level", 'Chromosome-level', 'Small-scale', 'no_cluster')

for (model_class in model_classes) {

shap_and_feature_path <- paste0(wd, "data/merged/SHAP_and_FeatureMatrix_",model_class,"_AmplDel.rds")
pred_ampl_path <- paste0(wd, "data/InteractomeINSIDER/",model_class,"-pred_ampl.rds")
pred_del_path <- paste0(wd, "data/InteractomeINSIDER/",model_class,"-pred_del.rds")

df <- readRDS(file = shap_and_feature_path)
pred_ampl <- readRDS(file = pred_ampl_path) 
pred_del <- readRDS(file = pred_del_path)

####################### 

#######################  manipulate data and define variables ####################### 

res.ampl <- get_residuals(pred_ampl); res.ampl$labels <- paste0(res.ampl$bin, "-", res.ampl$Type) # load ampl residuals
res.del <- get_residuals(pred_del); res.del$labels <- paste0(res.del$bin, "-", res.del$Type)      # load del residuals

bool.plot <- F

model_types <- c("ampl", "del") # define model types

for (i in model_types) { # iterate over model types

message(paste0("Now processing ", model_class, " ", i, " model"))

if(i == 'ampl'){ # if model is amplification 
  
  shap <- df$models.shap.df[[grep('Amplification',names(df$models.shap.df))]] # take shap values
  values <- df$models.X[[grep('Amplification',names(df$models.X))]] # take actual values
  res <- res.ampl[,c("labels", "observed", "prediction", "residual")] # take residuals 
  score_type <- "ampl_score" # define the target variable
  pred <- pred_ampl # take predictions
  
} else if(i == 'del'){ # if model is deletion, do the same
  
  shap <- df$models.shap.df[[grep('Deletion',names(df$models.shap.df))]]
  values <- df$models.X[[grep('Deletion',names(df$models.X))]]
  res <- res.del[,c("labels", "observed", "prediction", "residual")]
  score_type <- "del_score"
  pred <- pred_del
  
}

values$labels <- paste0(values$bin,'-',values$Type)
values <- change_distance(values)

col.sum.zero <- colnames(shap %>% dplyr::select(-labels))[apply(shap %>% dplyr::select(-labels), 2, sum) == 0]
non.segment.features <- c("Chromosome_Length","Centromere_Length","Centromere_Type")

if(model_class %in% c('Mid-length','Small-scale')){
  
  cols.to.discard <- c("labels","BIAS", col.sum.zero, non.segment.features) # include non.segment.features if the model is != Arm- or Chr-level
  
}else{
  
  cols.to.discard <- c("labels","BIAS")
  
}

cols.to.keep <- colnames(shap)[!colnames(shap) %in% cols.to.discard] # extract cols to keep
cols.to.keep <- c(cols.to.keep,"labels")

shap_agg <- shap[, cols.to.keep] # subset for cols to keep 
values_agg <- values[, cols.to.keep] # subset for cols to keep
values_agg <- values_agg[match(shap_agg$labels, values_agg$labels), ]

actual_shap_mat <- shap_agg[,-ncol(shap_agg)] # remove labels (last column)
actual_values_mat <- values_agg[,-ncol(values_agg)] # remove labels (last column)

distfit_output <- distfit(y = res$residual) # fit distribution of residuals
mu <- distfit_output$par["mu"] # take mean
sigma <- distfit_output$par["sigma"] # take standard deviation
pvals <- 2 * pmin(pnorm(res$residual, mean = mu, sd = sigma), 1 - pnorm(res$residual, mean = mu, sd = sigma)) # take p-values
res$pval.residual <- pvals # assign p-values

pval_thr <- 0.6                                                                              # select threshold to discriminate "correct" vs "incorrect" prediction
diploidy_threshold <- median(res$prediction)                                                 # select threshold to discriminate "frequently diploid" vs "frequently non diploid"

ns_thresholds <- data.frame(model_class = c("Mid-length", "Mid-length",
                                            "Arm-level", "Arm-level",
                                            "Chromosome-level", "Chromosome-level",
                                            "Small-scale", "Small-scale",
                                            "no_cluster", "no_cluster"),
                            model_type = rep(c("ampl", "del"),5),
                            thresholds = c(0.02, 0.025, 0.002, 0.010, 0.01, 0, 0.0004, 0, 0.04, 0.04))

ns_threshold <- ns_thresholds[(ns_thresholds$model_class == model_class & ns_thresholds$model_type == i),"thresholds"]

res$pred_is_correct <- ifelse(res$pval.residual > pval_thr, T, F)

res$region_is_diploid <- ifelse(res$observed > diploidy_threshold, F, T)

res$binID <- gsub(pattern = "-(.*)", replacement = "", x = res$labels)

res <- res[,c("labels", "binID", "observed", "prediction", "pred_is_correct")]

####################### 

####################### compute overlap between genes and bins ####################### 

backbone <- do.call(rbind, chr_backbone_namesfixed$`0.1Mbp`)
colnames(backbone) <- c("chr", "bin", "start", "end")
backbone$binID <- paste0(backbone$chr, "_", backbone$bin)
backbone$bin <- NULL

gene.table.coords <- gene.table %>% 
                      dplyr::select(
                        Gene.name, 
                        chromosome,
                        cds_from,
                        cds_to
                      )

to.discard <- apply(X = gene.table.coords, 
                    MARGIN = 1, 
                    FUN = function(row){
                      return(any(is.na(row)))
                      })

gene.table.coords.filtered <- gene.table.coords[!to.discard, ]
colnames(gene.table.coords.filtered) <- c("name", "chr", "start", "end")
gene.table.coords.filtered <- gene.table.coords.filtered[,c(2,3,4,1)]

backbone_gr <- GRanges(backbone)
gene.table_gr <- GRanges(gene.table.coords.filtered)

hits <- findOverlaps(query = backbone_gr, subject = gene.table_gr)

overlapping_binIDs <- mcols(backbone_gr[queryHits(hits)])$binID
overlapping_gene_names <- mcols(gene.table_gr[subjectHits(hits)])$name

overlap_df <- data.frame(binID = overlapping_binIDs, 
                         gene_name = overlapping_gene_names)

overlap_df <- overlap_df %>% group_by(binID) %>%
  summarise(
    num_of_genes = n(),
    genes = paste0(gene_name, collapse = "::")
  )

non_overlapping_binIDs <- mcols(backbone_gr[-unique(queryHits(hits))])$binID
non_overlapping_gene_num_of_genes <- rep(0,length(non_overlapping_binIDs))
non_overlapping_genes <- rep(NA,length(non_overlapping_binIDs))

non_overlap_df <- data.frame(binID = non_overlapping_binIDs, 
                             num_of_genes = non_overlapping_gene_num_of_genes,
                             genes = non_overlapping_genes)


global_overlap_df <- rbind(overlap_df, non_overlap_df)
  
####################### 

####################### compute segment specific performances ####################### 

performances_with_num_genes <- merge(res, global_overlap_df, by = "binID")
performances_with_num_genes$has_genes <- performances_with_num_genes$num_of_genes > 0

genes_free_segments <- performances_with_num_genes[performances_with_num_genes$has_genes == F,]
genes_populated_segments <- performances_with_num_genes[performances_with_num_genes$has_genes == T,]

gene_free_test <- cor.test(genes_free_segments$observed, genes_free_segments$prediction, method = "spearman", exact = F)
gene_populated_test <- cor.test(genes_populated_segments$observed, genes_populated_segments$prediction, method = "spearman", exact = F)
global_performances <- cor.test(res$observed, res$prediction, method = "spearman", exact = F)

p1_name <- paste0(wd, "data/plots/", model_class, "_", i, "_genes_free_scatterplot.pdf")
pdf(p1_name)
plot(genes_free_segments$observed,
     genes_free_segments$prediction,
     main = "Genes-free regions Scatterplot",
     xlab = "Observed",
     ylab = "Predicted")
dev.off()

p2_name <- paste0(wd, "data/plots/", model_class, "_", i, "_genes_populated_scatterplot.pdf")
pdf(p2_name)
plot(genes_populated_segments$observed,
     genes_populated_segments$prediction,
     main = "Genes-populated regions Scatterplot",
     xlab = "Observed",
     ylab = "Predicted")
dev.off()

p3_name <- paste0(wd, "data/plots/", model_class, "_", i, "_all_regions_scatterplot.pdf")
pdf(p3_name)
plot(res$observed,
     res$prediction,
     main = "All Regions Scatterplot",
     xlab = "Observed",
     ylab = "Predicted")
dev.off()

row <- c(gene_free_test$estimate,
         gene_populated_test$estimate,
         global_performances$estimate)

names(row) <- c("gene_free_all",
                "gene_populated_all",
                "overall_performances")
row <- list(row)
names(row) <- paste0(model_class, "-", i)

output <- c(output, row)

####################### 

 }
}

output <- do.call(rbind, output)

write.table(x = output, 
            file = paste0(wd, "data/segment_specific/segment_specific_performances.tsv"), 
            quote = F, sep = "\t", row.names = T, col.names = T)

df_long <- output %>%
  pivot_longer(
    cols = c(gene_free_all, gene_populated_all, overall_performances),
    names_to = "metric",
    values_to = "value"
  )

p <- ggplot(df_long, aes(x = category, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    x = "Category",
    y = "Value",
    fill = "Metric",
    title = "Performance metrics per model category"
  )

ggsave(paste0(wd, "data/plots/CNA_performance_grouped_barplot.pdf"), plot = p,
       width = 10, height = 6, dpi = 300)



