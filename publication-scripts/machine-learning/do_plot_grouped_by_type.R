wd <- 'path/to/GitHub/copy-number-annotation/'

do_plot_grouped_by_type <- function(input_df, tot_type_bins, col.to.group, distinct_colors, save = TRUE) {
  
  num_per_type <- input_df %>% group_by(type, !!sym(col.to.group)) %>% summarise(count_annots = n())
  
  annot_per_type <- merge(x = tot_type_bins, y = num_per_type, by = "type")
  
  annot_per_type$type <- as.factor(annot_per_type$type)
  annot_per_type <- annot_per_type[order(annot_per_type$type, decreasing = F),]
  annot_per_type <- annot_per_type %>% group_by(type, !!sym(col.to.group)) %>% summarise(count_bins = unique(count_bins), perc_annot = (count_annots / count_bins)*100)
  
  perc_grouped_by_type <- ggplot(annot_per_type, aes(x = type, y = perc_annot, fill = !!sym(col.to.group))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = distinct_colors) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0("Annotation Distribution by Tumor Type [", i, "]"),
      subtitle = paste0(col.to.group, ": Percentage of genomic bins annotated per category per tumor type"),
      x = "Tumor Type",
      y = "Percentage of Annotated Bins",
      fill = col.to.group
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(margin = margin(b = 10))
    )
  
  if (save) {
    
    if ("top1" %in% col.to.group) {
      filename <- paste0(wd, "data/plots/", i, "_annot_grouped_by_tt_detailed.pdf")
    } else {
      filename <- paste0(wd, "data/plots/", i, "_", model_class, "_", pval_thr, "_", diploidy_threshold, "_", ns_threshold, "_grouped_by_tt.pdf")
    }
    
    ggsave(filename = filename, plot = perc_grouped_by_type, width = 12, height = 6, dpi = 300)
    
  }
  
  outlist <- list(data = annot_per_type,
                  plot = perc_grouped_by_type)
  
  return(outlist)
  
}
