do_pancancer_barplots_and_boxplots <- function(input_df, save = TRUE){
  
  clean_theme <- theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid = element_blank()
    )
  
  pancancer_boxplot_3 <- ggplot(input_df, aes(x = annot_final, y = observed)) +
    geom_boxplot(fill = "skyblue", color = "black", outlier.size = 0.8) +
    ggtitle("Pan-cancer observed (final strategy)") +
    ylab(score_type) +
    xlab(NULL) +
    clean_theme
  
  pancancer_barplot_3 <- input_df %>%
    ggplot(aes(x = annot_final)) +
    geom_bar(fill = "steelblue") +
    ggtitle("Pan-cancer count (final strategy)") +
    ylab("Count") +
    xlab(NULL) +
    clean_theme

  outlist <- list(pancancer_boxplot_3 = pancancer_boxplot_3,
                  pancancer_barplot_3 = pancancer_barplot_3)
                  
  
  if (save) {
    
    ggsave(filename = paste0("../Data/plots/", model_class, "_", i, "_", pval_thr, "_", diploidy_threshold, "_", ns_threshold, "_pancancer_boxplot_final.pdf"), plot = pancancer_boxplot_3, width = 14, height = 6, dpi = 300)
    ggsave(filename = paste0("../Data/plots/", model_class, "_", i, "_", pval_thr, "_", diploidy_threshold, "_", ns_threshold, "_pancancer_barplot_final.pdf"), plot = pancancer_barplot_3, width = 14, height = 6, dpi = 300)

  }

  return(outlist)
}
  