landscape_plot_observed_prediction <- function(filtered_landscape,
                                               cluster_mask, genome_mask, type_mask, model_mask,
                                               plot_observed = TRUE, 
                                               plot_predicted = TRUE,
                                               annot_to_plot,
                                               annot_to_plot_ticks = "all",
                                               annot_to_plot_kde = "all",
                                               make.interactive = TRUE) {
  # make.interactive controls what is returned:
  #   TRUE     -> a single girafe (interactive) object   [on-screen / HTML export]
  #   FALSE    -> a single ggplot (static, with legend)  [PDF export]
  #   "both"   -> list(interactive = <girafe>, static = <ggplot>)
  #               Use this so the (expensive) data prep and layer-building
  #               below runs ONCE per user selection, and both the screen
  #               view and every download format are built from the same
  #               result instead of re-filtering/re-plotting from scratch
  #               for each output.
  
  valid_input <- c("ampl", "del")
  if (!all(model_mask %in% valid_input)) stop("Invalid model selected. Use 'ampl' and/or 'del'.")
  
  valid_tracks <- c("Observed", "Predicted")
  track_mask <- valid_tracks[c(plot_observed, plot_predicted)]
  if (length(track_mask) == 0) track_mask <- NA
  if (length(track_mask) > 1) track_mask <- paste(track_mask, collapse = ", ")
  
  if (length(genome_mask) == 22) genome_mask <- "WHOLE GENOME"
  if (length(genome_mask) > 1) genome_mask <- paste(genome_mask, collapse = ", ")
  
  if (annot_to_plot == "annot_final"){
  
    # ANNOT_FINAL_CLASSES (0_LoadData.R) FIRST, so a given class always
    # keeps the same colour across selections. union() rather than the
    # fixed list alone: should an unforeseen class ever reach the data, it
    # still gets a palette entry instead of silently dropping out of the
    # plot with a missing fill.
    unique_classes <- union(
      ANNOT_FINAL_CLASSES,
      unique(as.character(filtered_landscape[[annot_to_plot]]))
    )
  
  } else {
  
    unique_classes <- unique(filtered_landscape[,annot_to_plot])

  }
  
  n_of_classes <- length(unique_classes)
  
  if (identical(annot_to_plot_ticks, "all")) {
    annot_to_plot_ticks <- unique_classes
  } else if (!isFALSE(annot_to_plot_ticks)) {
    annot_to_plot_ticks <- as.character(annot_to_plot_ticks)
    missing_ticks <- setdiff(annot_to_plot_ticks, unique_classes)
    if (length(missing_ticks)) warning("Some annot_to_plot_ticks not found: ", paste(missing_ticks, collapse = ", "))
    annot_to_plot_ticks <- intersect(annot_to_plot_ticks, unique_classes)
  }
  
  if (identical(annot_to_plot_kde, "all")) {
    annot_to_plot_kde <- unique_classes
  } else if (!isFALSE(annot_to_plot_kde)) {
    annot_to_plot_kde <- as.character(annot_to_plot_kde)
    missing_kde <- setdiff(annot_to_plot_kde, unique_classes)
    if (length(missing_kde)) warning("Some annot_to_plot_kde not found: ", paste(missing_kde, collapse = ", "))
    annot_to_plot_kde <- intersect(annot_to_plot_kde, unique_classes)
  }
  
  filtered_landscape <- filtered_landscape %>% dplyr::mutate(pos = dplyr::row_number())
  chr_bounds <- get_chr_bounds(filtered_landscape)
  base_plot <- plot_base_layer(chr_bounds, filtered_landscape)
  chr_to_plot <- unique(filtered_landscape$chr)
  
  if (plot_observed) base_plot <- plot_obs_layer(base_plot, chr_to_plot, filtered_landscape, model_mask)
  if (plot_predicted) base_plot <- plot_pred_layer(base_plot, chr_to_plot, filtered_landscape, model_mask)
  
  base_plot <- base_plot +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.2) +
    labs(x = "Genomic Position", y = paste0("CNA frequency (", cluster_label(cluster_mask), ")")) +
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
  upper_limit <- ceiling(max(filtered_landscape$observed, na.rm = TRUE) * 10) / 10
  lower_limit <- floor(min(filtered_landscape$prediction, na.rm = TRUE) * 10) / 10
  color_palette_all <- setNames(all_colors[seq_len(n_of_classes)], unique_classes)
  
  ticksize <- 0.1
  linewidth <- 0.1
  
  build_layers <- function(selected_names) {
    lapply(selected_names, function(name) {
      idx <- match(name, unique_classes)
      y <- (idx %% 2) + 1
      list(mode = ifelse(n_of_classes > 10, all_modes[y], "obs"), name = name)
    })
  }
  
  if (!isFALSE(annot_to_plot_ticks) && length(annot_to_plot_ticks)) {
    for (lay in build_layers(annot_to_plot_ticks)) {
      base_plot <- add_segment_layer(base_plot, filtered_landscape, lay$name, lay$mode, annot_to_plot, unique_classes,
                                     ticksize, color_palette_all, lower_limit, upper_limit)
    }
  }
  
  if (!isFALSE(annot_to_plot_kde) && length(annot_to_plot_kde)) {
    for (lay in build_layers(annot_to_plot_kde)) {
      base_plot <- add_density_layer(base_plot, filtered_landscape, lay$name, lay$mode, annot_to_plot, unique_classes,
                                     linewidth, color_palette_all, lower_limit, upper_limit)
    }
  }
  
  base_plot <- base_plot +
    scale_fill_manual(values = color_palette_all) +
    scale_colour_manual(values = color_palette_all) +
    ggnewscale::new_scale_fill()
  
  base_plot <- base_plot +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
  
  y_breaks <- pretty(c(lower_limit, upper_limit))
  actual_limits <- ggplot_build(base_plot)$layout$panel_scales_y[[1]]$range$range
  low <- actual_limits[1] - 0.1
  high <- actual_limits[2] + 0.1
  
  base_plot <- base_plot +
    geom_segment(aes(x = -Inf, xend = -Inf, y = min(y_breaks), yend = max(y_breaks)),
                 inherit.aes = FALSE, color = "black", linewidth = 1) +
    scale_x_continuous(breaks = chr_bounds %>% dplyr::mutate(center = (start + end)/2) %>% dplyr::pull(center),
                       labels = chr_bounds$chr, expand = c(0.005, 0)) +
    scale_y_continuous(breaks = y_breaks, labels = function(x) abs(x), expand = c(0, 0)) +
    coord_cartesian(ylim = c(low, high)) +
    theme(axis.line.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  if ("is_centromere" %in% colnames(filtered_landscape)) {
    centromere_positions <- filtered_landscape %>%
      dplyr::filter(is_centromere == TRUE) %>%
      dplyr::pull(pos)
    
    if (length(centromere_positions) > 0) {

      y_max <- max(filtered_landscape$prediction, na.rm = TRUE)
      
      base_plot <- base_plot +
        geom_segment(
          data = data.frame(x = centromere_positions),
          aes(x = x, xend = x, y = 0, yend = y_max),
          linetype = "dashed",
          color = "grey",
          linewidth = 0.3
        )
    }
  }
  
  # ---- finalize as an interactive girafe widget (screen / HTML export) ----
  finalize_interactive <- function(p) {
    
    p <- p + labs(title = "", subtitle = "")
    
    g <- girafe(
      ggobj    = p,
      fonts    = list(sans = "Roboto"),
      width_svg  = 10,
      height_svg = 6
    )
    
    g <- girafe_options(
      g,
      opts_tooltip(css = tooltip_css, delay_mouseover = 0, delay_mouseout = 0, offx = 10, offy = -10),
      opts_hover(css = hover_css),
      opts_toolbar(saveaspng = FALSE)
    )
    
    g
  }
  
  # ---- finalize as a static ggplot with a full legend (PDF / print) ----
  # Sized/styled so that a larger export canvas (see ggsave() call in
  # app.R) is used mostly for the PLOT panel: the legend keys/text are kept
  # compact so they don't grow proportionally with the export size.
  finalize_static <- function(p) {
    
    title <- "Segment Annotation (based on SHAP values)"
    subtitle <- paste0("[", genome_mask, "] [", type_mask, "] [", track_mask, "]")
    
    p <- p +
      labs(title = title, subtitle = subtitle) +
      guides(
        fill = guide_legend(
          title = "Selection Classes",
          title.position = "top",
          nrow = 2,
          byrow = TRUE,
          override.aes = list(shape = 15, size = 2.2, linetype = 0)
        ),
        colour = "none",
        linetype = "none",
        shape = "none"
      ) +
      geom_blank(aes(linetype = if (plot_observed) "Observed" else NULL)) +
      geom_blank(aes(linetype = if (plot_predicted) "Predicted" else NULL)) +
      scale_linetype_manual(
        name = "Data Type",
        values = c("Observed" = "solid", "Predicted" = "dashed"),
        guide = guide_legend(
          title.position = "top",
          override.aes = list(color = "black", fill = NA)
        )
      ) +
      theme(
        plot.title    = element_text(size = 15, face = "bold"),
        plot.subtitle = element_text(size = 10.5),
        axis.title    = element_text(size = 11.5),
        axis.text     = element_text(size = 8.5),
        legend.position   = "bottom",
        legend.direction  = "horizontal",
        legend.box        = "vertical",
        legend.box.just   = "left",
        legend.margin     = margin(t = 2, b = 2),
        legend.spacing.y  = unit(0.05, "cm"),
        legend.spacing.x  = unit(0.15, "cm"),
        legend.key.size   = unit(0.35, "cm"),
        legend.title      = element_text(size = 7.5, face = "bold"),
        legend.text       = element_text(size = 6.5),
        plot.margin       = margin(t = 6, r = 12, b = 2, l = 6)
      )
    
    p
  }
  
  if (identical(make.interactive, "both")) {
    return(list(interactive = finalize_interactive(base_plot), static = finalize_static(base_plot)))
  } else if (isTRUE(make.interactive)) {
    return(finalize_interactive(base_plot))
  } else {
    return(finalize_static(base_plot))
  }
}
