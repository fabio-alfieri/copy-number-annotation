landscape_plot_observed_prediction <- function(filtered_landscape,
                                         genome_mask,
                                         type_mask,
                                         model_mask) {
  
  # Basic input checks
  if (!model_mask %in% c("ampl", "del")) {
    stop("Invalid model selected. Use 'ampl' or 'del'.")
  }
  
  # Build readable labels
  track_label <- "Observed"
  genome_label <- if (length(genome_mask) == 22) "WHOLE GENOME" else paste(genome_mask, collapse = ", ")
  
  # Add genomic position index
  filtered_landscape <- filtered_landscape %>% mutate(pos = row_number())
  
  # Get chromosome boundaries (assumes get_chr_bounds is sourced or defined)
  chr_bounds <- get_chr_bounds(filtered_landscape = filtered_landscape)
  chr_to_plot <- unique(filtered_landscape$chr)
  
  # Start base plot
  base_plot <- plot_base_layer(chr_bounds = chr_bounds, filtered_landscape = filtered_landscape)
  
  # Add only observed layer (assumes plot_obs_layer is sourced or defined)
  base_plot <- plot_obs_layer(base_plot = base_plot, 
                              chr_to_plot = chr_to_plot, 
                              filtered_landscape = filtered_landscape, 
                              model = model_mask)
  
  base_plot <- plot_pred_layer(base_plot = base_plot, 
                              chr_to_plot = chr_to_plot, 
                              filtered_landscape = filtered_landscape, 
                              model = model_mask)
  
  # Add reference line and labels
  base_plot <- base_plot +
    geom_hline(yintercept = 0, 
               linetype = "dashed", 
               color = "grey", 
               linewidth = 0.2) +
    labs(
      title = "Observed SCNA Landscape",
      subtitle = paste0("[", genome_label, "] [", type_mask, "] [", track_label, "]"),
      x = "Genomic Position",
      y = "SCNA frequency (Mid-length)"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.line.y = element_blank()
    ) +
    scale_x_continuous(
      breaks = chr_bounds %>% mutate(center = (start + end)/2) %>% pull(center),
      labels = chr_bounds$chr,
      expand = c(0.005, 0)
    )
  
  return(base_plot)
}
