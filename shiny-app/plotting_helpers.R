get_chr_bounds <- function(filtered_landscape){
  
  chr_bounds <- filtered_landscape %>%
    group_by(chr) %>%
    summarize(start = min(pos), end = max(pos), .groups = "drop") %>%
    mutate(chr_num = readr::parse_number(chr)) %>%
    arrange(chr_num) %>%
    mutate(fill = rep(c("white", "#e7deed"), length.out = n())) %>%
    dplyr::select(-chr_num)
  
}

plot_base_layer <- function(chr_bounds, filtered_landscape){
  
  base_plot <- ggplot() +
    geom_rect(
      data = chr_bounds,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = fill),
      alpha = 0.3
    ) +
    scale_fill_identity() +
    new_scale_fill()
  
  return(base_plot)
}

plot_obs_layer <- function(base_plot, chr_to_plot, filtered_landscape, model){
  
  if (model == "ampl") {
    bg_obs <- "#FF0000"; fg_obs <- get_contrast(bg_obs)
    model_text <- "Amplification"
  } else if (model == "del") {
    bg_obs <- "#0000FF"; fg_obs <- get_contrast(bg_obs)
    model_text <- "Deletion"
  }
  
  for (chr in chr_to_plot) {
    chr_data <- filtered_landscape[filtered_landscape$chr == chr, ]
    chr_data <- chr_data %>% 
      mutate(
        coord = paste0(chr_data$chr, ":", chr_data$start, "-", chr_data$end),
        data_id = binID
      )
    base_plot <- base_plot +
      geom_line(data = chr_data, aes(x = pos, y = observed), color = bg_obs) +
      geom_point_interactive(
        data = chr_data,
        aes(x = pos, y = observed, 
            tooltip = sprintf(
              "<div style='background:%s; 
                  color:%s; 
                  padding:4px; 
                  border-radius:0px; 
                  border:none; 
                  outline:none; 
                  box-shadow:none;'>
                  Coordinates: %s<br>
                  Observed %s Frequency: %s</div>",
              bg_obs, fg_obs, coord, model_text, round(observed,3)
            ), 
            data_id = data_id),
        size = 3, color = "transparent"
      )
  }
  
  return(base_plot)
}

plot_pred_layer <- function(base_plot, chr_to_plot, filtered_landscape, model){
  
  if (model == "ampl") {
    bg_pred <- "#FF5257"; fg_pred <- get_contrast(bg_pred)
    model_text <- "Amplification"
  } else if (model == "del") {
    bg_pred <- "#1671FF"; fg_pred <- get_contrast(bg_pred)
    model_text <- "Deletion"
  }
  
  for (chr in chr_to_plot) {
    chr_data <- filtered_landscape[filtered_landscape$chr == chr, ]
    chr_data <- chr_data %>% 
      mutate(
        coord = paste0(chr_data$chr, ":", chr_data$start, "-", chr_data$end),
        data_id = binID
      )
    base_plot <- base_plot +
      geom_line(data = chr_data, aes(x = pos, y = prediction), color = bg_pred, alpha = 0.45) +
      geom_point_interactive(
        data = chr_data,
        aes(x = pos, y = prediction, 
            tooltip = sprintf(
              "<div style='background:%s; 
                  color:%s; 
                  padding:4px; 
                  border-radius:0px; 
                  border:none; 
                  outline:none; 
                  box-shadow:none;'>
                  Coordinates: %s<br>
                  Predicted %s Frequency: %s</div>",
              bg_pred, fg_pred, coord, model_text, round(prediction,3)
            ), 
            data_id = data_id),
        size = 3, color = "transparent"
      )
  }
  
  return(base_plot)
}

get_contrast <- function(hexcol) {
  rgb <- col2rgb(hexcol) / 255
  lum <- 0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]
  ifelse(lum > 0.5, "#000000", "#FFFFFF")
}

add_segment_layer <- function(base_plot, input_df, 
                              name_annot, mode, 
                              clustering_col, unique_classes, 
                              ticksize, color_palette_ticks, 
                              lower_limit, upper_limit,
                              available_height = 0.3) {

  generate_tick_df <- function(input_df, 
                               name_annot, mode, 
                               clustering_col, unique_classes,
                               color_palette_ticks, 
                               lower_limit, upper_limit,
                               available_height) {

    unique_vec <- as.character(unique_classes)
    n_classes <- max(length(unique_vec), 1)

    spacing <- available_height / n_classes
    height  <- min(0.015, spacing * 0.4)  # keep ticks small relative to spacing

    idx <- match(as.character(name_annot), unique_vec)
    if (is.na(idx)) idx <- 1L

    if (mode == "obs") {
      cluster_ymid <- upper_limit + (idx * spacing)
    } else {
      cluster_ymid <- lower_limit - (idx * spacing)
    }

    df <- input_df %>%
      filter(.data[[clustering_col]] == name_annot)

    if (nrow(df) == 0) return(NULL)

    df <- df %>%
      mutate(
        coord = paste0(chr, ":", start, "-", end),
        tooltip = sprintf(
          "<div style='background:%s; color:%s; padding:4px;'>Coords: %s<br>%s</div>",
          color_palette_ticks[as.character(name_annot)],
          get_contrast(color_palette_ticks[as.character(name_annot)]),
          coord,
          .data[[clustering_col]]
        ),
        cluster_ymid = cluster_ymid,
        cluster_ymin = cluster_ymid - height,
        cluster_ymax = cluster_ymid + height,
        cluster = factor(name_annot, levels = names(color_palette_ticks))
      ) %>%
      ungroup()

    return(df)
  }

  cluster_ticks <- generate_tick_df(input_df = input_df, 
                                    name_annot = name_annot, 
                                    mode = mode, 
                                    clustering_col = clustering_col,
                                    unique_classes = unique_classes,
                                    color_palette_ticks = color_palette_ticks, 
                                    lower_limit = lower_limit, upper_limit = upper_limit,
                                    available_height = available_height)

  if (is.null(cluster_ticks)) return(base_plot)

  base_plot <- base_plot +
    geom_rect_interactive(
      data = cluster_ticks,
      aes(
        xmin = pos - 1,
        xmax = pos + 1,
        ymin = cluster_ymin,
        ymax = cluster_ymax,
        tooltip = tooltip,
        data_id = paste0(.data[[clustering_col]], "_", binID)
      ),
      color = NA,
      linewidth = ticksize * 10,
      fill = NA
    ) +
    geom_rect(
      data = cluster_ticks,
      aes(
        xmin = pos - 0.5,
        xmax = pos + 0.5,
        ymin = cluster_ymin,
        ymax = cluster_ymax,
        fill = cluster
      ),
      color = NA,
      linewidth = ticksize
    )

  return(base_plot)
}


add_density_layer <- function(base_plot, input_df, 
                              name_annot, mode, 
                              clustering_col, unique_classes, 
                              linewidth, color_palette_ticks, 
                              lower_limit, upper_limit,
                              available_height = 0.3) {

  generate_density_df <- function(input_df, 
                                  name_annot, mode, 
                                  clustering_col, unique_classes, 
                                  color_palette_ticks, 
                                  lower_limit, upper_limit,
                                  available_height) {

    unique_vec <- as.character(unique_classes)
    n_classes <- max(length(unique_vec), 1)

    spacing <- available_height / n_classes
    height  <- min(0.015, spacing * 0.4)

    idx <- match(as.character(name_annot), unique_vec)
    if (is.na(idx)) idx <- 1L

    if (mode == "obs") {
      cluster_ymid <- upper_limit + (idx * spacing)
    } else {
      cluster_ymid <- lower_limit - (idx * spacing)
    }

    df <- input_df %>%
      filter(.data[[clustering_col]] == name_annot)

    if (nrow(df) == 0) return(NULL)

    df <- df %>%
      mutate(
        coord = paste0(chr, ":", start, "-", end),
        tooltip = sprintf(
          "<div style='background:%s; color:%s; padding:4px;'>%s</div>",
          color_palette_ticks[as.character(name_annot)],
          get_contrast(color_palette_ticks[as.character(name_annot)]),
          .data[[clustering_col]]
        ),
        cluster_ymid = cluster_ymid,
        cluster_ymin = cluster_ymid - height,
        cluster_ymax = cluster_ymid + height,
        cluster = factor(name_annot, levels = names(color_palette_ticks))
      ) %>%
      ungroup()

    return(df)
  }

  densities_input <- generate_density_df(input_df = input_df,
                                         name_annot = name_annot,
                                         mode = mode,
                                         clustering_col = clustering_col,
                                         unique_classes = unique_classes,
                                         color_palette_ticks = color_palette_ticks,
                                         lower_limit = lower_limit, upper_limit = upper_limit,
                                         available_height = available_height)

  if (is.null(densities_input)) return(base_plot)

  dens <- density(densities_input$pos, bw = 5)
  ymin <- unique(densities_input$cluster_ymin)
  ymax <- unique(densities_input$cluster_ymax)
  cluster <- unique(densities_input$cluster)
  tooltip <- unique(densities_input$tooltip)

  scaled_y <- (dens$y / max(dens$y)) * (ymax - ymin) * 0.9 + ymin

  dens_df <- data.frame(
    pos = dens$x,
    ymin = ymin,
    y = scaled_y,
    cluster = cluster,
    tooltip = tooltip
  )

  dens_df$cluster <- factor(dens_df$cluster, levels = names(color_palette_ticks))

  density_layer <- geom_ribbon_interactive(
    data = dens_df,
    aes(
      x = pos,
      ymin = ymin,
      ymax = y,
      fill = cluster,
      tooltip = tooltip,
      data_id = cluster
    ),
    linetype = "blank",
    alpha = 0.65
  )

  segment_layer <- geom_segment_interactive(
    data = dens_df,
    mapping = aes(
      y = ymin,
      yend = ymin,
      colour = cluster,
      tooltip = tooltip,
      data_id = pos
    ),
    x = min(input_df$pos),
    xend = max(input_df$pos),
    linewidth = linewidth
  )

  base_plot <- base_plot +
    density_layer +
    segment_layer

  return(base_plot)
}



add_segment_gradient_layer <- function(base_plot, input_df, 
                                       value_col,             # column name as string
                                       mode = "obs", 
                                       clustering_col,        # column name for class labels (string)
                                       unique_classes,        # vector of unique classes
                                       ticksize = 0.1,
                                       lower_limit, upper_limit,
                                       available_height = 0.3,
                                       mid_point = 0) {
  
  unique_vec <- as.character(unique_classes)
  n_classes <- max(length(unique_vec), 1)
  
  spacing <- available_height / n_classes
  height  <- min(0.015, spacing * 0.4)
  
  if (value_col == "shap_top1") {
    input_df <- input_df %>%
      filter(!annot_final %in% c(
        "Incorrect prediction - Positive Selection", 
        "Incorrect prediction - Negative Selection"
      ))
  } else if (value_col %in% c("log2_selection_to_occurrence_ratio", 
                              "log2_positive_to_negative_ratio")) {
    input_df <- input_df %>%
      filter(!annot_final %in% c(
        "Incorrect prediction - Positive Selection", 
        "Incorrect prediction - Negative Selection"
      ))
  }
  
  cluster_ticks <- input_df %>%
    mutate(
      value_numeric = as.numeric(.data[[value_col]]),
      cluster_idx = match(as.character(.data[[clustering_col]]), unique_vec),
      cluster_idx = ifelse(is.na(cluster_idx), 1L, cluster_idx),
      cluster_ymid = if (mode == "obs") {
        upper_limit + cluster_idx * spacing
      } else {
        lower_limit - cluster_idx * spacing
      },
      cluster_ymin = cluster_ymid - height,
      cluster_ymax = cluster_ymid + height,
      coord = paste0(chr, ":", start, "-", end),
      tooltip = sprintf(
        "<div style='padding:4px;'>Coords: %s<br>%s: %.3f</div>",
        coord,
        value_col,
        value_numeric
      )
    ) %>%
    filter(!is.na(value_numeric))
  
  if (nrow(cluster_ticks) == 0) return(base_plot)
  
  base_plot <- base_plot + ggnewscale::new_scale_fill()
  
  base_plot <- base_plot +
    geom_rect_interactive(
      data = cluster_ticks,
      aes(
        xmin = pos - 0.5,
        xmax = pos + 0.5,
        ymin = cluster_ymin,
        ymax = cluster_ymax,
        fill = value_numeric,
        tooltip = tooltip,
        data_id = paste0("grad_", binID)
      ),
      color = NA,
      linewidth = ticksize
    ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = mid_point,
      name = value_col
    )
  
  return(base_plot)
}




