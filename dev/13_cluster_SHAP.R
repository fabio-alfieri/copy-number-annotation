# Signature dei cluster (?)
rm(list=ls())
gc(full=T)

# setwd('/Users/gabry/OneDrive/Desktop/shiny_app/') # Gab
setwd('/Users/ieo5099/Desktop/copy-number-annotation/') # Fab

shap_sum_abs_cluster <- list()
for(clustering_depth in 1:4){
  
  source('dev/0_LoadData.R')
  source('dev/1_dynamic_plotting_functions.R')
  source('dev/1bis_parse_input_data_FA.R') # overwrite parse_input_data()
  
  processed_data <- parse_input_data(shap.list = shap.list,
                                     out_annot_list = out_annot_list,
                                     chr_backbone_namesfixed = chr_backbone_namesfixed, 
                                     centromere_table = centromere_table,
                                     clustering_depth = clustering_depth)
  
  # these files will be saved after
  
  shap.list <- processed_data$shap.list
  out_annot_list_processed <- processed_data$out_annot_list_processed
  backbone.100kb <- processed_data$backbone.100kb
  centromere_table <- processed_data$centromere_table
  
  k <- paste0('k',2^clustering_depth)
  top_k <- paste0('top_',k)
  shap_sum_abs_cluster[[clustering_depth]] <- list()
  
  for(i in c('ampl','del')){
    shap_annotation <- full_join(shap.list[[i]], out_annot_list_processed[[i]])
    
    shap_sum_abs_cluster[[clustering_depth]][[i]] <- shap_annotation %>% group_by(.data[[k]], .data[[top_k]]) %>% 
      summarize(across(
      -c(binID, pos, ampl, del, chr, type, labels),
      ~ sum(.),
      .names = "{.col}"))
  }
}


pdf(file = 'dev/Plots/cluster_SHAP.pdf', height = 12, width = 12)
for(clustering_depth in 1:4){
  for(i in 1:2){
    top_k <- paste0('top_k',2^clustering_depth)
    print(pivot_longer(shap_sum_abs_cluster[[clustering_depth]][[i]], cols = -c(1:2)) %>%
      ggplot(aes(x = reorder(name, value), y = value, fill = name)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = paste('Clustering depth',clustering_depth, '-', ifelse(i == 1, 'Ampl', 'Del')),
           x = "Feature",
           y = "Sum of Absolute SHAP Values") +
      theme_minimal() +
      theme(legend.position = "none") +
      facet_wrap(~.data[[top_k]]))
  }
}
dev.off()
