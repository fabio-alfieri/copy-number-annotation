do_pancancer_barplots_and_boxplots <- function(input_df, save = TRUE){
  
  # prefix <- strsplit(x = deparse(substitute(input_df)), split = "_")[[1]][3]
  
  clean_theme <- theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid = element_blank()
    )
  
# comparisons_5_classes <- list(
#   c("Negative Selection", "Positive Selection"),
#    c("Negative Selection", 'Incorrect prediction - Positive Selection'),
#    c("Negative Selection", 'Incorrect prediction - Negative Selection'),
#    c('Incorrect prediction - Negative Selection', "Positive Selection"),
#    c('Incorrect prediction - Positive Selection', "Positive Selection"),
#    c('Incorrect prediction - Negative Selection', 'Incorrect prediction - Positive Selection'),
#    c('Negative Selection', 'Relaxed Selection'),
#    c('Positive Selection', 'Relaxed Selection')
#  )
  
#  comparisons_3_classes <- list(
#    c("Negative Selection", "Occurrence"),
#    c("Negative Selection", "Positive Selection"),
#    c("Occurrence", "Positive Selection")
#  )
  
  pancancer_boxplot_3 <- ggplot(input_df, aes(x = annot_final, y = observed)) +
    geom_boxplot(fill = "skyblue", color = "black", outlier.size = 0.8) +
#    geom_signif(
#      comparisons = comparisons_3_classes,
#      test = "wilcox.test",
#      step_increase = 0.1,
#      textsize = 2.5,
#      map_signif_level = TRUE
#    ) +
    ggtitle("Pan-cancer observed (final strategy)") +
    ylab(score_type) +
    xlab(NULL) +
    clean_theme
  
  
#  pancancer_boxplot_5 <- ggplot(input_df, aes(x = annot_5_classes, y = observed)) +
#    geom_boxplot(fill = "lightgreen", color = "black", outlier.size = 0.8) +
#    geom_signif(
#      comparisons = comparisons_5_classes,
#      test = "wilcox.test",
#      step_increase = 0.1,
#      textsize = 2.5,
#      map_signif_level = TRUE
#    ) +
#    ggtitle("Pan-cancer observed (5 classes)") +
#    ylab(score_type) +
#    xlab(NULL) +
#    clean_theme
  

  pancancer_barplot_3 <- input_df %>%
    ggplot(aes(x = annot_final)) +
    geom_bar(fill = "steelblue") +
    ggtitle("Pan-cancer count (final strategy)") +
    ylab("Count") +
    xlab(NULL) +
    clean_theme

#  pancancer_barplot_5 <- input_df %>%
#    ggplot(aes(x = annot_5_classes)) +
#    geom_bar(fill = "darkseagreen") +
#    ggtitle("Pan-cancer count (5 classes)") +
#    ylab("Count") +
#    xlab(NULL) +
#    clean_theme
  
#  outlist <- list(pancancer_boxplot_3 = pancancer_boxplot_3,
#                  pancancer_boxplot_5 = pancancer_boxplot_5,
#                  pancancer_barplot_3 = pancancer_barplot_3,
#                  pancancer_barplot_5 = pancancer_barplot_5)
  
  outlist <- list(pancancer_boxplot_3 = pancancer_boxplot_3,
                  pancancer_barplot_3 = pancancer_barplot_3)
                  
  
  if (save) {
    
    ggsave(filename = paste0(model_class, "_", i, "_", pval_thr, "_", diploidy_threshold, "_", ns_threshold, "_pancancer_boxplot_final.pdf"), plot = pancancer_boxplot_3, width = 14, height = 6, dpi = 300)
#    ggsave(filename = paste0(i, "_", prefix, "_pancancer_boxplot_5.pdf"), plot = pancancer_boxplot_5, width = 5, height = 4, dpi = 300)
    ggsave(filename = paste0(model_class, "_", i, "_", pval_thr, "_", diploidy_threshold, "_", ns_threshold, "_pancancer_barplot_final.pdf"), plot = pancancer_barplot_3, width = 14, height = 6, dpi = 300)
#    ggsave(filename = paste0(i, "_", prefix, "_pancancer_barplot_5.pdf"), plot = pancancer_barplot_5, width = 5, height = 4, dpi = 300)
    
  }
  
#   if (prefix == "annot") {
    
#     pancancer_boxplot_detailed <- input_df %>%
#      ggplot(aes(x = top1, y = observed, fill = annot_3_classes)) +
#      geom_boxplot(outlier.size = 0.7) +
#      scale_fill_manual(values = c("Positive Selection" = "green", "Negative Selection" = "purple", "Occurrence" = "grey"), name = "Feature Type") +
#      labs(
#        title = paste0("Observed Value ", i, " over top1 by Pos/Neg Feature Type"),
#        x = "top1",
#        y = score_type
#      ) +
#      clean_theme
    
#    pancancer_barplot_detailed <- input_df %>%
#      ggplot(aes(x = top1, fill = annot_3_classes)) +
#      geom_bar(color = "black") +
#      scale_fill_manual(values = c("Positive Selection" = "green", "Negative Selection" = "purple", "Occurrence" = "grey"), name = "Feature Type") +
#      ggtitle("Pan-cancer count (5 classes)") +
#      ylab("Count") +
#      xlab(NULL) +
#      clean_theme
    
#    outlist[["pancancer_boxplot_detailed"]] <- pancancer_boxplot_detailed
#    outlist[["pancancer_barplot_detailed"]] <- pancancer_barplot_detailed
    
#    if (save) {
      
#      ggsave(filename = paste0(i, "_", prefix, "_pancancer_boxplot_detailed.pdf"), plot = pancancer_boxplot_detailed, width = 15, height = 10, dpi = 300)
#      ggsave(filename = paste0(i, "_", prefix, "_pancancer_barplot_detailed.pdf"), plot = pancancer_barplot_detailed, width = 15, height = 10, dpi = 300)
      
#    }
    
#  }
  
  return(outlist)
}
  