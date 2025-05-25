rm(list=ls())
gc(full=T)

setwd('/Users/gabry/OneDrive/Desktop/shiny_app/') # Gab
# setwd('/Users/ieo5099/Desktop/copy-number-annotation/') # Fab

source('dev/0_LoadData.R')
source('dev/1_dynamic_plotting_functions.R')
source('dev/1bis_parse_input_data_FA.R') # overwrite parse_input_data()

processed_data <- parse_input_data(shap.list = shap.list,
                                   out_annot_list = out_annot_list,
                                   chr_backbone_namesfixed = chr_backbone_namesfixed, 
                                   centromere_table = centromere_table,
                                   clustering_depth = 4)


shap.list <- processed_data$shap.list
out_annot_list_processed <- processed_data$out_annot_list_processed
backbone.100kb <- processed_data$backbone.100kb
centromere_table <- processed_data$centromere_table

shap_ampl_clean <- shap.list$ampl
shap_ampl_clean_BRCA <- shap_ampl_clean[shap_ampl_clean$type == "BRCA",]

rownames_ampl_BRCA <- shap_ampl_clean_BRCA$binID

shap_ampl_clean_BRCA$binID <- NULL; shap_ampl_clean_BRCA$chr <- NULL; 
shap_ampl_clean_BRCA$type <- NULL; shap_ampl_clean_BRCA$labels <- NULL

shap_ampl_clean_BRCA <- t(shap_ampl_clean_BRCA)
colnames(shap_ampl_clean_BRCA) <- rownames_ampl_BRCA

cor.mat_ampl_BRCA <- cor(shap_ampl_clean_BRCA[1:13,1:100])

rand_idxs <- sample(x = 1:ncol(shap_ampl_clean_BRCA), size = 100)
cor.mat_ampl_random_BRCA <- cor(shap_ampl_clean_BRCA[1:13,rand_idxs])

shap_del_clean <- shap.list$del

shap_del_clean$binID <- NULL; shap_del_clean$chr <- NULL; 
shap_del_clean$type <- NULL; shap_del_clean$labels <- NULL




