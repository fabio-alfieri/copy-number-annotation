landscape_plot_observed_prediction <- function(filtered_landscape,
                                               genome_mask, type_mask, model_mask,
                                               annot_to_plot,
                                               annot_to_plot_kde = "all",
                                               make.interactive = TRUE) {
  
  
  valid_input <- c("ampl", "del")
  if (!all(model_mask %in% c(valid_input))) stop("Invalid model selected. Use 'ampl' and/or 'del'.")
  
  valid_tracks <- c("Observed", "Predicted")
  track_mask <- valid_tracks
  if (length(track_mask) == 0) { track_mask <- NA }
  if (length(track_mask) > 1) track_mask <- paste(track_mask, collapse = ", ")
  
  if (length(genome_mask) == 22) genome_mask <- "WHOLE GENOME"
  if (length(genome_mask) > 1) genome_mask <- paste(genome_mask, collapse = ", ")
  
  unique_classes <- c("Negative Selection", "Incorrect prediction - Negative Selection", "Positive Selection", "Incorrect prediction - Positive Selection",  "Occurrence", "No Detectable Force")
  n_of_classes <- length(unique_classes)
  
  annot_to_plot_ticks <- factor(unique_classes)
  annot_to_plot_kde <- factor(unique_classes)
  
  title <- ifelse(make.interactive, "", "Segment Annotation (based on SHAP values)")
  subtitle <- ifelse(make.interactive, "", paste0("[", genome_mask, "] [", type_mask, "] [", track_mask, "]"))
  
  filtered_landscape <- filtered_landscape %>% mutate(pos = row_number())
  
  chr_bounds <- get_chr_bounds(filtered_landscape = filtered_landscape)
  
  base_plot <- plot_base_layer(chr_bounds = chr_bounds, filtered_landscape = filtered_landscape)
  
  chr_to_plot <- unique(filtered_landscape$chr)
  
  base_plot <- plot_obs_layer(base_plot = base_plot, 
                              chr_to_plot = chr_to_plot, 
                              filtered_landscape = filtered_landscape, 
                              model = model_mask)

  
  base_plot <- plot_pred_layer(base_plot = base_plot, 
                               chr_to_plot = chr_to_plot, 
                               filtered_landscape = filtered_landscape, 
                               model = model_mask)

  base_plot <- base_plot +
    geom_hline(yintercept = 0, 
               linetype = "dashed", 
               color = "grey", 
               linewidth = 0.2) +
    labs(title = title, 
         subtitle = subtitle, 
         x = "Genomic Position", 
         y = "SCNA frequency (Mid-length)") +
    theme_classic() + 
    theme(legend.position = "none")
  
  all_colors <- c(
    "#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe",
    "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000",
    "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080",
    "#FFFFFF", "#000000", "#bcf60c", "#ff4500", "#6495ED",
    "#ff69b4", "#7fff00", "#ff1493", "#00ced1", "#dda0dd"
  )
  
  all_modes <- c("obs", "pred")
  
  upper_limit <- ceiling(max(filtered_landscape$observed) * 10) / 10
  lower_limit <- floor(min(filtered_landscape$prediction) * 10) / 10
  
  color_palette_ticks <- all_colors[1:n_of_classes]
  names(color_palette_ticks) <- annot_to_plot_ticks
  
  ticksize <- 0.1
  linewidth <- 0.1
  

  layers_kde <- lapply(X = annot_to_plot_kde, FUN = function(x){
    y <- (as.integer(x) %% 2) + 1
    list(mode = ifelse(n_of_classes > 50, all_modes[y], "obs"),
         name = annot_to_plot_kde[x])
  })
    
  for (lay in layers_kde) {
    base_plot <- add_density_layer(base_plot = base_plot, 
                                   input_df = filtered_landscape, 
                                   name_annot = lay$name,
                                   mode = lay$mode,  
                                   clustering_col = annot_to_plot, 
                                   unique_classes = unique_classes,
                                   linewidth = linewidth, 
                                   color_palette_ticks = color_palette_ticks, 
                                   lower_limit = lower_limit, 
                                   upper_limit = upper_limit)
  }
  
  base_plot <- base_plot +
    scale_fill_manual(values = color_palette_ticks) +
    scale_colour_manual(values = color_palette_ticks)
  
  base_plot <- base_plot + ggnewscale::new_scale_fill()
  
  base_plot <- base_plot +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
  
  y_breaks <- pretty(c(lower_limit, upper_limit))
  
  actual_limits <- ggplot_build(base_plot)$layout$panel_scales_y[[1]]$range$range
  low <- actual_limits[1] - 0.1
  high <- actual_limits[2] +  0.1
  
  base_plot <- base_plot +
    geom_segment(
      aes(x = -Inf, xend = -Inf,
          y = (min(y_breaks)),
          yend = (max(y_breaks))),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 1
    ) +
    scale_x_continuous(
      breaks = chr_bounds %>% 
        mutate(center = (start + end)/2) %>% 
        pull(center),
      labels = chr_bounds$chr,
      expand = c(0.005, 0)
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = function(x){abs(x)},
      expand = c(0, 0)
    ) +
    coord_cartesian(ylim = c(low, high)) +
    theme(
      axis.line.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  
  if (make.interactive) {
    base_plot <- girafe(
      ggobj = base_plot,
      fonts = list(sans = "Roboto"),
      width_svg  = 10,
      height_svg = 6,
      options  = list(
        opts_tooltip(
          delay_mouseover = 0,
          delay_mouseout  = 0,
          offx            = 10,
          offy            = -10
        )
      )
    )
  } else {
    base_plot <- base_plot +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.box.just = "center",
        legend.margin = margin(t = 5, b = 5),
        legend.text = element_text(size = 8)
      )
  }  
  
  return(base_plot)
  
}






