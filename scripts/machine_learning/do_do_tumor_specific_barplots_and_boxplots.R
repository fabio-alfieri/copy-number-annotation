do_do_tumor_specific_barplots_and_boxplots <- function(ct,
                                                       comparisons_to_plot, 
                                                       input_df, 
#                                                       num, prefix, 
                                                       save = TRUE){

# x_axis_label <- num

#  fill_colors_detailed <- c(
#    "Negative Selection" = "purple",
#    "Occurrence" = "gray",
#    "Positive Selection" = "green"
#  )
  
#  if (num == "detailed") {
    
#    p_boxplot <- ggplot(input_df, aes(x = annot_to_plot, y = observed, fill = annot_3_classes)) +
#      geom_boxplot(color = "black", outlier.size = 0.8) +
#      scale_fill_manual(values = fill_colors_detailed) +
#      geom_signif(
#        comparisons = comparisons_to_plot,
#        test = function(x, y) wilcox.test(x, y, exact = FALSE),
#        step_increase = 0.1,
#        textsize = 2.5,
#        map_signif_level = TRUE
#      ) +
#      ggtitle(ct) +
#      xlab(x_axis_label) +
#      ylab(score_type) +
#      theme_classic(base_size = 8) +
#      theme(
#        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
#        axis.title.x = element_text(size = 8, face = "bold"),
#        axis.title.y = element_text(size = 8, face = "bold"),
#        axis.text.x = element_text(size = 5, angle = 70, hjust = 1),
#        axis.text.y = element_text(size = 7),
#        strip.text = element_text(size = 7)
#      )
#    
#    p_barplot <- ggplot(input_df, aes(x = annot_to_plot, fill = annot_3_classes)) +
#      geom_bar(color = "black") +
#      scale_fill_manual(values = fill_colors_detailed) +
#      ggtitle(ct) +
#      xlab(x_axis_label) +
#      ylab("# of regions") +
#      theme_classic(base_size = 8) +
#      theme(
#        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
#        axis.title.x = element_text(size = 8, face = "bold"),
#        axis.title.y = element_text(size = 8, face = "bold"),
#        axis.text.x = element_text(size = 5, angle = 70, hjust = 1),
#        axis.text.y = element_text(size = 7),
#        legend.position = "right"
#      )
#    
#  } else {
    
    p_boxplot <- ggplot(input_df, aes(x = annot_to_plot, y = observed)) +
      geom_boxplot(fill = "skyblue", color = "black", outlier.size = 0.8) +
#      geom_signif(
#        comparisons = comparisons_to_plot,
#        test = function(x, y) wilcox.test(x, y, exact = FALSE), 
#        step_increase = 0.1,
#        textsize = 2.5,
#        map_signif_level = TRUE
#      ) +
      ggtitle(ct) +
      xlab("final strategy") +
      ylab(score_type) +
      theme_classic(base_size = 8) +
      theme(
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 8, face = "bold"),
        axis.title.y = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 5, angle = 70, hjust = 1),
        axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 7)
      )
    
    p_barplot <- ggplot(input_df, aes(x = annot_to_plot)) +
      geom_bar(fill = "steelblue", color = "black") +
      ggtitle(ct) +
      xlab("final_strategy") +
      ylab("# of regions") +
      theme_classic(base_size = 8) +
      theme(
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 8, face = "bold"),
        axis.title.y = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 5, angle = 70, hjust = 1),
        axis.text.y = element_text(size = 7),
        legend.position = "none"
      )
#  }
  
  if (save) {
    
    ggsave(filename = paste0(model_class, "_", i, "_", pval_thr, "_", diploidy_threshold, "_", ns_threshold, "_", ct, "_boxplot_final.pdf"), plot = p_boxplot, width = 14, height = 6, dpi = 300)
    ggsave(filename = paste0(model_class, "_", i, "_", pval_thr, "_", diploidy_threshold, "_", ns_threshold, "_", ct, "_barplot_final.pdf"), plot = p_barplot, width = 14, height = 6, dpi = 300)
    
  }
  
  outlist <- list(p_boxplot = p_boxplot,
                  p_barplot = p_barplot)
  
  return(outlist)
    
}
