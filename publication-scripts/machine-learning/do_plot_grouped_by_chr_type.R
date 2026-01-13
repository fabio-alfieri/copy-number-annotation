do_plot_grouped_by_chr_type <- function(input_df, tot_chrom_bins_per_type, col.to.group, distinct_colors, save = TRUE){
  
  num_per_chr_per_type <- input_df %>% group_by(type, chr, !!sym(col.to.group)) %>% summarise(count_annots = n())
  
  annot_per_chr_per_type <- merge(x = tot_chrom_bins_per_type, y = num_per_chr_per_type, by = c("type", "chr"))
  
  annot_per_chr_per_type <- annot_per_chr_per_type %>% group_by(type, chr, !!sym(col.to.group)) %>% summarise(count_bins = unique(count_bins), perc_annot = (count_annots / count_bins)*100)
  
  perc_grouped_by_chr_type <- ggplot(annot_per_chr_per_type, aes(x = chr, y = perc_annot, fill = !!sym(col.to.group))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = distinct_colors) +
    facet_wrap(~ type, scales = "free_y") +
    scale_x_discrete(drop = FALSE) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0("Annotation Distribution by chromosome and type [", i, "]"),
      subtitle = paste0(col.to.group, ": Percentage of genomic bins annotated per category per chromosome per type"),
      x = "Chromosome",
      y = "Percentage of Annotated Bins",
      fill = col.to.group
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 9, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(margin = margin(b = 10))
    )
  
  if (save) {
    
    if ("top1" %in% col.to.group) {
      filename <- paste0("../Data/plots/", i, "_annot_grouped_by_chr_type_detailed.pdf")
    } else {
      filename <- paste0("../Data/plots/", i, "_", model_class, "_", pval_thr, "_", diploidy_threshold, "_", ns_threshold, "_grouped_by_chr_tt.pdf")
    }
    
    ggsave(filename = filename, plot = perc_grouped_by_chr_type, width = 12, height = 6, dpi = 300)
    
  }
  
  outlist <- list(data = annot_per_chr_per_type,
                  plot = perc_grouped_by_chr_type)
  
  return(outlist)
  
}
