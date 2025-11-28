do_tumor_specific_barplots_and_boxplots <- function(input_df, save = TRUE) {
  
# prefix <- strsplit(x = deparse(substitute(input_df)), split = "_")[[1]][3]  
  
#  comparisons_5_classes <- list(
#    c("Negative Selection", "Positive Selection"),
#    c("Negative Selection", 'Incorrect prediction - Positive Selection'),
#    c("Negative Selection", 'Incorrect prediction - Negative Selection'),
#    c('Incorrect prediction - Negative Selection', "Positive Selection"),
#    c('Incorrect prediction - Positive Selection', "Positive Selection"),
#    c('Incorrect prediction - Negative Selection', 'Incorrect prediction - Positive Selection')
#  )
  
#  comparisons_3_classes <- list(
#    c("Negative Selection", "Occurrence"),
#    c("Negative Selection", "Positive Selection"),
#    c("Occurrence", "Positive Selection")
#  )
  
  cancer_types <- unique(input_df$type)
  boxplots <- list()
  barplots <- list()
  
  input_df$annot_to_plot <- input_df$annot_final
  num <- "final strategy"
  
#  for (num in what.to.plot) {
    
#    if (num == "3") {
#      input_df$annot_to_plot <- input_df$annot_3_classes
#      comparisons_to_plot <- comparisons_3_classes
#    } else if (num == "5") {
#      input_df$annot_to_plot <- input_df$annot_5_classes
#      comparisons_to_plot <- comparisons_5_classes
#    } else if (num == "detailed") {
#      input_df$annot_to_plot <- input_df$top1
#      comparisons_to_plot <- list()
#    }
    
    for (ct in cancer_types) {
      
      plot_data <- input_df %>% filter(type == ct)
      
      ct_spec_plots <- do_do_tumor_specific_barplots_and_boxplots(ct = ct,
#                                                 comparisons_to_plot = comparisons_to_plot, 
                                                 input_df = plot_data, 
#                                                 num = num, 
#                                                 prefix = prefix, 
                                                 save = save)
      
      p_boxplot <- ct_spec_plots$p_boxplot
      p_barplot <- ct_spec_plots$p_barplot
      
      boxplots[[paste0(num, "_", ct)]] <- p_boxplot
      barplots[[paste0(num, "_", ct)]] <- p_barplot
      
    }
#  }
    
  outlist <- list(boxplots = boxplots,
                  barplots = barplots)
  
  return(outlist)
  
}
