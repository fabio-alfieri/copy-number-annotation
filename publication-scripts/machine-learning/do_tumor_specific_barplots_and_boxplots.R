do_tumor_specific_barplots_and_boxplots <- function(input_df, save = TRUE) {
  
  cancer_types <- unique(input_df$type)
  boxplots <- list()
  barplots <- list()
  
  input_df$annot_to_plot <- input_df$annot_final
  num <- "final strategy"
  
    for (ct in cancer_types) {
      
      plot_data <- input_df %>% filter(type == ct)
      
      ct_spec_plots <- do_do_tumor_specific_barplots_and_boxplots(ct = ct,
                                                 input_df = plot_data, 
                                                 save = save)
      
      p_boxplot <- ct_spec_plots$p_boxplot
      p_barplot <- ct_spec_plots$p_barplot
      
      boxplots[[paste0(num, "_", ct)]] <- p_boxplot
      barplots[[paste0(num, "_", ct)]] <- p_barplot
      
    }

  outlist <- list(boxplots = boxplots,
                  barplots = barplots)
  
  return(outlist)
  
}
