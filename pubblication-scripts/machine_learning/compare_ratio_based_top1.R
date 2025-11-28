compare_ratio_based_top1 <- function(input_df, save){
  
  clean_base_theme <- theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 9, face = "bold"),
      axis.title.y = element_text(size = 9, face = "bold"),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8)
    )
  
  log2_pos_neg_over_top1 <<- input_df %>%
    filter(is.finite(log2_positive_to_negative_ratio)) %>%
    ggplot(aes(x = top1, y = log2_positive_to_negative_ratio, fill = annot_3_classes_annot)) +
    geom_boxplot(outlier.size = 0.7) +
    scale_fill_manual(values = c("Positive Selection" = "green", "Negative Selection" = "purple", "Occurrence" = "grey"), name = "Feature Type") +
    labs(
      title = "log2(pos_sel/neg_sel) Ratio over top1 by Feature Type",
      x = "top1",
      y = "log2(pos_sel/neg_sel)"
    ) +
    clean_base_theme
  
  shap_summary_pos_neg <- input_df %>%
    group_by(annot_5_classes_ratio, annot_5_classes_annot) %>%
    summarise(sum_abs_shap = sum(abs(as.numeric(shap_top1)), na.rm = TRUE)) %>%
    ungroup()
  
  importance_contributions <<- ggplot(shap_summary_pos_neg, aes(x = annot_5_classes_annot, y = sum_abs_shap, fill = annot_5_classes_annot)) +
    geom_col() +
    scale_fill_manual(values = c("Positive Selection" = "green",
                                 "Incorrect prediction - Positive Selection" = "lightgreen",
                                 "Negative Selection" = "purple", 
                                 "Incorrect prediction - Negative Selection" = "thistle",
                                 "Occurrence" = "grey")) +
    facet_wrap(~ annot_5_classes_ratio, scales = "free_x", nrow = 1, labeller = label_wrap_gen(width = 25)) +
    labs(
      fill = "Feature Type",
      x = "Feature Type",
      y = "Sum(abs(SHAP top1))",
      title = "SHAP Importance by Annotation (5 classes) and Pos/Neg Feature Type"
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title = element_text(size = 9, face = "bold"),
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "grey90"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
  
  if (save) {
    
    ggsave(paste0(i, "_log2_pos_neg_sel_ratio_over_top1_by_type_of_feature.pdf"), plot = log2_pos_neg_over_top1, width = 10, height = 6)
    ggsave(paste0(i, "_sum_abs_SHAP_top1_by_annot_5_classes_and_pos_neg_type_of_feature.pdf"), plot = importance_contributions, width = 12, height = 6)
    
  }
}
