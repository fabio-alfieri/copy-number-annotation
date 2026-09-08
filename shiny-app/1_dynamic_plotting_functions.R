parse_input_coord <- function(input, chrom_sizes){
  
  # chrom_sizes: named numeric vector, one entry per chromosome, giving the
  # maximum coordinate seen in the LOADED backbone/annotation files (not an
  # external genome assembly) - e.g. chrom_sizes["chr1"]. This is what lets
  # us validate/clamp a user-entered region against the data that is
  # actually available, without any extra genome-annotation package.
  
  raw_input <- input
  
  # Be forgiving of common typing mistakes: stray whitespace, thousands
  # separators ("1,000,000"), and "Chr"/"CHR" instead of "chr".
  input <- tolower(trimws(gsub("[,\\s]+", "", input, perl = TRUE)))
  
  valid_pattern <- "^chr(?:[1-9]|1[0-9]|2[0-2]):[0-9]+-[0-9]+$"
  
  if (!grepl(valid_pattern, input, ignore.case = FALSE)) {
    stop(
      "Invalid coordinate format: \"", raw_input, "\".\n",
      "Use the format chr:start-end, e.g. chr1:1-5000000 ",
      "(chr1-chr22 only, with integer start and end, start < end)."
    )
  }
  
  elems <- strsplit(input, "[:-]")[[1]]
  chr   <- elems[1]
  start <- suppressWarnings(as.numeric(elems[2]))
  end   <- suppressWarnings(as.numeric(elems[3]))
  
  if (is.na(start) || is.na(end)) {
    stop("Invalid coordinates: START and END must be integers. E.g. chr1:1-5000000")
  }
  
  if (start >= end) {
    stop("START (", start, ") must be smaller than END (", end, ").")
  }
  
  if (!(chr %in% names(chrom_sizes)) || is.na(chrom_sizes[[chr]])) {
    stop(
      "Chromosome \"", chr, "\" not found in the loaded data.\n",
      "Available chromosomes: ", paste(sort(names(chrom_sizes)), collapse = ", ")
    )
  }
  
  chrom_size <- chrom_sizes[[chr]]
  clamped <- FALSE
  
  if (start < 1) {
    start <- 1
    clamped <- TRUE
  }
  
  if (end > chrom_size) {
    end <- chrom_size
    clamped <- TRUE
  }
  
  if (start >= end) {
    stop(
      "The requested interval is entirely outside the boundaries of ", chr,
      " in the loaded data (maximum size: ", chrom_size, ")."
    )
  }
  
  input_granges <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges   = IRanges::IRanges(start = start, end = end)
  )
  
  # Flag on the object (rather than a side-effect / UI call here) so the
  # caller (filter_df -> app.R) can decide how to surface the "we trimmed
  # your region to the chromosome's actual boundaries" message to the user.
  attr(input_granges, "clamped") <- clamped
  
  return(input_granges)
}
parse_input_chr <- function(input){
  
  valid_input <- "^chr(?:[1-9]|1[0-9]|2[0-2])$"
  
  is_valid <- all(grepl(pattern = valid_input, x = input, ignore.case = F))
  
  if (is_valid) {
    
    return(input)
    
  } else {
    
    stop("INVALID chromosome specification. \nENTER a valid chromosome specification. \nExample: chr[1-22]")
    
  }
  
}
parse_input_type <- function(input){
  
  accepted_types <- c("STAD", "GBMLGG", "COADREAD",
                      "KIRP", "KIRC", "OV", "ESCA",
                      "LUAD", "LUSC", "PAAD", "BRCA")
  
  if (length(input) > 1) {
    stop(paste("Only 1 cancer type can be selected at a time. \n\nSelect among: \n", 
               paste(accepted_types, collapse = ", ")))
  }
  
  if (!(input %in% accepted_types)) {
    stop(paste("This research only supports some cancer types. \n\nSelect among: \n", 
               paste(accepted_types, collapse = ", ")))
  }
  
  return(input)
}
parse_input_model <- function(input){
  
  valid_input <- c("ampl","del")
  is_valid <- input %in% valid_input
  
  if (is_valid) {
    
    not_selected <- valid_input[valid_input != input]
    
    outlist <- list(selected = input,
                    not_selected = not_selected)
    
    return(outlist)
    
  } else {
    
    stop("Invalid model selected. \n Models are either \"ampl\" or \"del\"")
    
  }
}


parse_input_cluster <- function(input){
  # Keep in sync with model_types in 0_LoadData.R. "No-cluster" is the
  # internal id of the category shown to the user as "All CNAs".
  valid_input <- c("Mid-length", "Small-scale", "Arm-level", "Chromosome-level", "No-cluster")
  is_valid <- input %in% valid_input
  
  if (is_valid) {
    
    not_selected <- valid_input[valid_input != input]
    
    outlist <- list(selected = input,
                    not_selected = not_selected)
    
    return(outlist)
    
  } else {
    
    stop("Invalid cluster selected. \n Models are either \"Mid-length\", \"Small-scale\", \"Arm-level\", \"Chromosome-level\" or \"No-cluster\"")
    
  }  
}

parse_annot_to_plot <- function(clustering_depth, input){
  
  if (isFALSE(input)) {
    return(input)
  }
  
  valid_input <- seq_len(clustering_depth)
  
  if ((length(input) == 1) && (input == "all")) {
    return(valid_input)
  } 
  
  is_valid <- all(as.integer(input) %in% valid_input)
  
  if (is_valid) { 
    return(sort(as.integer(input)))
  } else {
    stop("Invalid clusters selected")
  }
}

filter_df <- function(input_obj, 
                      backbone_granges,
                      cluster_input = NULL,
                      type_input = NULL, 
                      model_input = NULL,
                      chr_input = NULL, 
                      coord_input = NULL){
  
  # model filtering policy:
  # either "ampl" or "del" must be specified
  
  cluster_mask <- parse_input_cluster(cluster_input)
  
  model_mask <- parse_input_model(model_input)
  
  if (length(model_mask$selected) > 1) {
    stop("SHAP list or landscape df MUST be explicitly filtered")       
  }
  
  if (is.function(input_obj)) {
    
    # lazy-loading path (see 0_LoadData.R::meta_list_get): the requested
    # table is read/parsed on first use and cached in memory afterwards,
    # so only the model_type x model combinations the user actually
    # selects are ever loaded.
    df_input <- input_obj(cluster_mask$selected, model_mask$selected)
    
  } else if (is.list(input_obj)) {
    
    df_input <- input_obj[[cluster_mask$selected]][[model_mask$selected]]
    
  } else {
    
    stop("input_obj MUST either be a lazy-loading function or a nested list") 
    
  }
  
  # Type filtering policy:
  # Only one cancer type to specify
  # Must be in the 11 cancer types
  #
  # NOTE: filtering by type FIRST (before the region filter below) shrinks
  # the table before the coordinate-overlap / chromosome membership check
  # runs, which matters once tables hold every cancer type at full size
  # rather than a subsample.
  
  type_mask <- parse_input_type(type_input)
  
  type_filtered_df <- df_input[df_input$type %in% type_mask, ]
  
  # Region filtering policy:
  # if both chromosome and coordinate level are NULL --> plot whole genome
  # if both are defined --> plot the most detailed one (coordinates)
  #
  # Chromosome sizes used to validate/clamp a manually-typed region are
  # derived from the LOADED backbone bins themselves (max bin end per
  # chromosome) rather than from an external genome-assembly package, per
  # the "usa i confini nei file di input" requirement.
  
  region_clamped <- FALSE
  
  # Overlap the requested region against the FILTERED TABLE'S OWN start/end
  # coordinates directly (not against backbone_granges' binID scheme).
  # binID matching against the original 100kb backbone breaks the moment
  # the table's bins don't use that exact scheme any more - e.g. after
  # smoothing/downsampling into coarser windows with their own binID
  # ("chr1_0", "chr1_1", ...), which never matches the backbone's
  # ("chr1_36", ...) and silently returned zero rows. Overlapping on actual
  # genomic coordinates works regardless of bin size/naming.
  overlap_by_coord <- function(region_granges, df) {
    row_granges <- GenomicRanges::GRanges(
      seqnames = df$chr,
      ranges   = IRanges::IRanges(start = df$start, end = df$end)
    )
    hits <- GenomicRanges::findOverlaps(query = region_granges, subject = row_granges)
    df[unique(S4Vectors::subjectHits(hits)), ]
  }
  
  if ((!is.null(chr_input)) && (!is.null(coord_input))) {
    
    chrom_sizes <- tapply(GenomicRanges::end(backbone_granges),
                          as.character(GenomicRanges::seqnames(backbone_granges)),
                          max)
    
    genome_mask <- parse_input_coord(coord_input, chrom_sizes)
    region_clamped <- isTRUE(attr(genome_mask, "clamped"))
    
    coord_filtered_df <- overlap_by_coord(genome_mask, type_filtered_df)
    
  } else if ((is.null(chr_input)) && (is.null(coord_input))) {
    
    genome_mask <- unique(type_filtered_df$chr)
    coord_filtered_df <- type_filtered_df
    
  } else if ((!is.null(chr_input)) && (is.null(coord_input))){
    
    genome_mask <- parse_input_chr(chr_input)
    coord_filtered_df <- type_filtered_df[type_filtered_df$chr %in% genome_mask, ]
    
  } else {
    
    chrom_sizes <- tapply(GenomicRanges::end(backbone_granges),
                          as.character(GenomicRanges::seqnames(backbone_granges)),
                          max)
    
    genome_mask <- parse_input_coord(coord_input, chrom_sizes)
    region_clamped <- isTRUE(attr(genome_mask, "clamped"))
    
    coord_filtered_df <- overlap_by_coord(genome_mask, type_filtered_df)
    
  }
  
  if (nrow(coord_filtered_df) == 0) {
    stop(
      "No bins found for the requested region/selection. ",
      "Try a wider interval or check the chromosome selection."
    )
  }
  
  outlist <- list(final_df = as.data.frame(coord_filtered_df),
                  model_mask = model_mask$selected,
                  type_mask = type_mask,
                  genome_mask = as.character(genome_mask),
                  region_clamped = region_clamped)
  
  return(outlist)
  
}









































































########################################################################## OBSOLETE

landscape_plot_interactive_prediction <- function(filtered_landscape,
                                                  backbone.100kb,
                                                  genome_mask, type_mask, model_mask,
                                                  plot_observed = TRUE, plot_predicted = TRUE,
                                                  annot_to_plot_ticks = "all",
                                                  annot_to_plot_kde = "all") {
  
  get_chr_bounds <- function(filtered_landscape){
    
    chr_bounds <- filtered_landscape %>%
      group_by(chr) %>%
      summarize(start = min(pos), end = max(pos), .groups = "drop") %>%
      mutate(chr_num = readr::parse_number(chr)) %>%
      arrange(chr_num) %>%
      mutate(fill = rep(c("white", "#e7deed"), length.out = n())) %>%
      select(-chr_num)
    
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
  
  plot_obs_layer <- function(base_plot, chr_to_plot, filtered_landscape, backbone.100kb, model){
    
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
        rowwise() %>% 
        mutate(
          coord = as.character(backbone.100kb[mcols(backbone.100kb)$binID == binID][1]),
          data_id = binID
        )
      base_plot <- base_plot +
        geom_line(data = chr_data, aes(x = pos, y = obs), color = bg_obs) +
        geom_point_interactive(
          data = chr_data,
          aes(x = pos, y = obs, 
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
                bg_obs, fg_obs, coord, model_text, round(obs,3)
              ), 
              data_id = data_id),
          size = 3, color = "transparent"
        )
    }
    
    return(base_plot)
  }
  
  plot_pred_layer <- function(base_plot, chr_to_plot, filtered_landscape, backbone.100kb, model){
    
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
        rowwise() %>% 
        mutate(
          coord = as.character(backbone.100kb[mcols(backbone.100kb)$binID == binID][1]),
          data_id = binID
        )
      base_plot <- base_plot +
        geom_line(data = chr_data, aes(x = pos, y = pred), color = bg_pred, alpha = 0.45) +
        geom_point_interactive(
          data = chr_data,
          aes(x = pos, y = pred, 
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
                bg_pred, fg_pred, coord, model_text, round(obs,3)
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
  
  add_density_layer <- function(base_plot, input_df, 
                                name_annot,mode, 
                                clustering_col, top_clustering_col, clustering_depth, 
                                backbone.100kb, 
                                linewidth, color_palette_ticks, 
                                lower_limit, upper_limit) {
    
    
    generate_density_df <- function(input_df, 
                                    name_annot,mode, 
                                    clustering_col, top_clustering_col, clustering_depth, 
                                    backbone.100kb,
                                    color_palette_ticks, 
                                    lower_limit, upper_limit) {
      
      height <- 0.015
      bg <- color_palette_ticks[as.character(name_annot)]
      fg <- get_contrast(bg)
      
      df <- input_df %>%
        filter(.data[[clustering_col]] == name_annot)
      
      if (nrow(df) == 0) return(NULL)
      
      if (mode == "obs") {
        start <- upper_limit
        df <- df %>%
          rowwise() %>%
          mutate(
            tooltip = sprintf(
              "<div style='background:%s; color:%s; padding:4px;'>%s</div>",
              bg, fg, .data[[top_clustering_col]]
            ),
            cluster_ymid = (round(start, 1) + 0.1) + ((.data[[clustering_col]] / (clustering_depth / 3.8)) * 0.1),
            cluster_ymin = cluster_ymid - height,
            cluster_ymax = cluster_ymid + height
          ) %>%
          ungroup()
      } else {
        start <- lower_limit - 0.05
        df <- df %>%
          rowwise() %>%
          mutate(
            coord = as.character(backbone.100kb[mcols(backbone.100kb)$binID == binID][1]),
            tooltip = sprintf(
              "<div style='background:%s; color:%s; padding:4px;'>%s</div>",
              bg, fg, .data[[top_clustering_col]]
            ),
            cluster_ymid = (round(start, 1) - 0.1) - ((.data[[clustering_col]] / (clustering_depth / 3.8)) * 0.1),
            cluster_ymin = cluster_ymid - height,
            cluster_ymax = cluster_ymid + height
          ) %>%
          ungroup()
        
      }
      df$cluster <- factor(name_annot, levels = names(color_palette_ticks))
      return(df)
      
    }
    
    
    densities_input <- generate_density_df(input_df = input_df, 
                                           name_annot = name_annot, 
                                           mode = mode, 
                                           clustering_col = clustering_col, 
                                           top_clustering_col = top_clustering_col, 
                                           clustering_depth = clustering_depth, 
                                           backbone.100kb =  backbone.100kb, 
                                           color_palette_ticks = color_palette_ticks, 
                                           lower_limit = lower_limit, upper_limit = upper_limit)
    
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
    
  }
  
  
  add_segment_layer <- function(base_plot, input_df, 
                                name_annot,mode, 
                                clustering_col, top_clustering_col, clustering_depth, 
                                backbone.100kb, 
                                ticksize, color_palette_ticks, 
                                lower_limit, upper_limit) {
    
    
    generate_tick_df <- function(input_df, 
                                 name_annot,mode, 
                                 clustering_col, top_clustering_col, clustering_depth, 
                                 backbone.100kb,
                                 color_palette_ticks, 
                                 lower_limit, upper_limit) {
      
      height <- 0.015
      bg <- color_palette_ticks[as.character(name_annot)]
      fg <- get_contrast(bg)
      
      df <- input_df %>%
        filter(.data[[clustering_col]] == name_annot)
      
      if (nrow(df) == 0) return(NULL)
      
      if (mode == "obs") {
        start <- upper_limit
        df <- df %>%
          rowwise() %>%
          mutate(
            coord = as.character(backbone.100kb[mcols(backbone.100kb)$binID == binID][1]),
            tooltip = sprintf(
              "<div style='background:%s; color:%s; padding:4px;'>Coords: %s<br>%s</div>",
              bg, fg, coord, .data[[top_clustering_col]]
            ),
            cluster_ymid = (round(start, 1) + 0.1) + ((.data[[clustering_col]] / (clustering_depth / 3.8)) * 0.1),
            cluster_ymin = cluster_ymid - height,
            cluster_ymax = cluster_ymid + height
          ) %>%
          ungroup()
      } else {
        start <- lower_limit - 0.05
        df <- df %>%
          rowwise() %>%
          mutate(
            coord = as.character(backbone.100kb[mcols(backbone.100kb)$binID == binID][1]),
            tooltip = sprintf(
              "<div style='background:%s; color:%s; padding:4px;'>Coords: %s<br>%s</div>",
              bg, fg, coord, .data[[top_clustering_col]]
            ),
            cluster_ymid = (round(start, 1) - 0.1) - ((.data[[clustering_col]] / (clustering_depth / 3.8)) * 0.1),
            cluster_ymin = cluster_ymid - height,
            cluster_ymax = cluster_ymid + height
          ) %>%
          ungroup()
        
      }
      
      df$cluster <- factor(name_annot, levels = names(color_palette_ticks))
      return(df)
      
    }
    
    
    
    
    cluster_ticks <- generate_tick_df(input_df = input_df, 
                                      name_annot = name_annot, 
                                      mode = mode, 
                                      clustering_col = clustering_col, 
                                      top_clustering_col = top_clustering_col, 
                                      clustering_depth = clustering_depth, 
                                      backbone.100kb =  backbone.100kb, 
                                      color_palette_ticks = color_palette_ticks, 
                                      lower_limit = lower_limit, upper_limit = upper_limit)
    
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
          data_id = paste0(.data[[top_clustering_col]], "_", binID)
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
  }
  
  valid_input <- c("ampl", "del")
  valid_tracks <- c("Observed", "Predicted")
  track_mask <- valid_tracks[c(plot_observed, plot_predicted)]
  
  if (length(track_mask) == 0) {
    model_mask <- NA
  }
  if (!all(model_mask %in% c(valid_input,NA))) stop("Invalid model selected. Use 'ampl' and/or 'del'.")
  if (length(genome_mask) == 22) genome_mask <- "WHOLE GENOME"
  if (length(genome_mask) > 1) genome_mask <- paste(genome_mask, collapse = ", ")
  if (length(track_mask) > 1) track_mask <- paste(track_mask, collapse = ", ")
  
  title    <- "Segment Annotation (based on SHAP values)"; title <- ""
  subtitle <- paste0("[", genome_mask, "] [", type_mask, "] [", track_mask, "]"); subtitle <- ""
  
  filtered_landscape <- filtered_landscape %>% mutate(pos = row_number())
  
  chr_bounds <- get_chr_bounds(filtered_landscape = filtered_landscape)
  
  message("Plotting base layer...")
  base_plot <- plot_base_layer(chr_bounds = chr_bounds, filtered_landscape = filtered_landscape)
  
  chr_to_plot <- unique(filtered_landscape$chr)
  
  if (plot_observed) {
    message("Plotting Observed track...")
    base_plot <- plot_obs_layer(base_plot = base_plot, 
                                chr_to_plot = chr_to_plot, 
                                filtered_landscape = filtered_landscape, 
                                backbone.100kb = backbone.100kb, 
                                model = model_mask)
  }
  
  if (plot_predicted) {
    message("Plotting Prediction track...")
    base_plot <- plot_pred_layer(base_plot = base_plot, 
                                 chr_to_plot = chr_to_plot, 
                                 filtered_landscape = filtered_landscape, 
                                 backbone.100kb = backbone.100kb, 
                                 model = model_mask)
  }
  
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
  
  all_colors <-c( "#412336", "#cc0000", "#0000cc", "#007700", 
                  "#800080", "#ff8000", "#999900", "#00aaaa", 
                  "#ff66cc", "#8b4513", "#cc00cc", "#3399cc", 
                  "#FFC0CB", "#000000", "#FF0000", "#EE82EE") # madonna che bello così simmetrico
  
  all_modes <- c("obs", "pred")
  
  upper_limit <- ceiling(max(filtered_landscape$obs) * 10) / 10
  lower_limit <- 0
  
  clustering_col <- grep(pattern = "^k\\d{1,2}$", x = colnames(filtered_landscape), value = T)
  top_clustering_col <- paste0("top_",clustering_col)
  clustering_depth <- as.integer(gsub(pattern = "k", x = clustering_col, replacement = ""))
  
  annot_to_plot_ticks <- parse_annot_to_plot(clustering_depth = clustering_depth, input = annot_to_plot_ticks)
  annot_to_plot_kde <- parse_annot_to_plot(clustering_depth = clustering_depth, input = annot_to_plot_kde)
  
  color_palette_ticks <- all_colors[1:clustering_depth]
  names(color_palette_ticks) <- levels(factor(1:clustering_depth))
  
  ticksize <- 0.1
  linewidth <- 0.1
  
  if (!isFALSE(annot_to_plot_ticks)) {
    message("Plotting Ticks layers....")
    layers_ticks <- lapply(X = annot_to_plot_ticks, FUN = function(x){
      y <- (x %% 2) + 1
      return(
        list(mode = all_modes[y],
             name = names(color_palette_ticks)[x]
        )
      )
    }
    )
    
    for (lay in layers_ticks) {
      
      base_plot <- add_segment_layer(base_plot = base_plot, 
                                     input_df = filtered_landscape, 
                                     name_annot = lay$name,
                                     mode = lay$mode, 
                                     clustering_col = clustering_col,
                                     top_clustering_col = top_clustering_col,
                                     clustering_depth = clustering_depth,
                                     backbone.100kb = backbone.100kb,
                                     ticksize = ticksize, 
                                     color_palette_ticks = color_palette_ticks, 
                                     lower_limit = lower_limit, upper_limit = upper_limit
      )
    }
  }
  
  if (!isFALSE(annot_to_plot_kde)) {
    message("Plotting KDE layers....")
    layers_kde <- lapply(X = annot_to_plot_kde, FUN = function(x){
      y <- (x %% 2) + 1
      return(
        list(mode = all_modes[y],
             name = names(color_palette_ticks)[x]
        )
      )
    }
    )
    
    for (lay in layers_kde) {
      
      base_plot <- add_density_layer(base_plot = base_plot, 
                                     input_df = filtered_landscape, 
                                     name_annot = lay$name,
                                     mode = lay$mode, 
                                     clustering_col = clustering_col,
                                     top_clustering_col = top_clustering_col,
                                     clustering_depth = clustering_depth,
                                     backbone.100kb = backbone.100kb,
                                     linewidth = linewidth, 
                                     color_palette_ticks = color_palette_ticks, 
                                     lower_limit = lower_limit, upper_limit = upper_limit
      )
    }
  }
  
  message("Almost done!")                                     
  
  y_breaks <- pretty(c(lower_limit, upper_limit))
  
  base_plot <- base_plot +
    geom_segment(
      aes(x = -Inf, xend = -Inf,
          y = (min(y_breaks)-0.003),
          yend = (max(y_breaks)+0.001)),
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
    coord_cartesian(ylim = c(-0.6, 1.2)) +
    theme(
      axis.line.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) + 
    scale_fill_manual(values = color_palette_ticks) +
    scale_colour_manual(values = color_palette_ticks)
  
  message("Making the plot interactive...")
  p <- girafe(
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
  message("The plot is ready! Enjoy!")
  return(p)
}

prepare_shap_to_plot <- function(filtered_shap_ampl, filtered_shap_del){
  
  to_clean <- c("labels", "type", "chr", "binID")
  mask <- !(colnames(filtered_shap_ampl) %in% to_clean)
  
  filtered_shap_abs_sum_ampl <- data.frame(value = apply(X = filtered_shap_ampl[,mask], MARGIN = 2, FUN = function(x){sum(abs(x))}))
  filtered_shap_abs_sum_del <- data.frame(value = apply(X = filtered_shap_del[,mask], MARGIN = 2, FUN = function(x){sum(abs(x))}))
  
  filtered_shap_abs_sum_ampl$feature <- rownames(filtered_shap_abs_sum_ampl); filtered_shap_abs_sum_ampl$color <- rainbow(nrow(filtered_shap_abs_sum_ampl))
  filtered_shap_abs_sum_del$feature <- rownames(filtered_shap_abs_sum_del); filtered_shap_abs_sum_del$color <- rainbow(nrow(filtered_shap_abs_sum_del))
  
  outlist <- list(filtered_shap_abs_sum_ampl = filtered_shap_abs_sum_ampl,
                  filtered_shap_abs_sum_del = filtered_shap_abs_sum_del)
  
  return(outlist)
}

barplot_shap <- function(shap.abs.sum, genome_mask, type_mask, model_mask){
  
  if (length(genome_mask) == 22) {
    genome_mask <- "WHOLE GENOME"
  }
  
  if (length(genome_mask) > 1) {
    genome_mask <- paste(genome_mask, collapse = ", ")
  } 
  
  title <- "SHAP Value Contribution per Feature"
  subtitle <- paste0(paste0(" [", genome_mask, "] "), 
                     paste0(" [", type_mask, "] "), 
                     paste0(" [", model_mask, "] "))
  
  ggplot(shap.abs.sum, aes(x = reorder(feature, value), y = value, fill = feature)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = title,
         subtitle = subtitle,
         x = "Feature",
         y = "Sum of Absolute SHAP Values") +
    theme_minimal() +
    theme(legend.position = "none")
  
}

prepare_landscape_to_plot <- function(model_input, shap_plotting_list, 
                                      filtered_shap_output, filtered_landscape, pred_list){
  
  if (model_input$selected == "ampl") {
    
    filtered_shap_abs_sum <- shap_plotting_list$filtered_shap_abs_sum_ampl
    genome_mask <- filtered_shap_output$genome_mask
    model_mask <- filtered_shap_output$model_mask
    type_mask <- filtered_shap_output$type_mask
    
    pred_df <- pred_list[[model_input$selected]]; pred_df$ampl_score <- NULL
    
    filtered_landscape <- merge(x = filtered_landscape, by.x = c("type", "binID"),
                                y = pred_df,            by.y = c("Type", "bin"),
                                sort = F
    )
    
    clustering_col <- grep(pattern = "^k\\d{1,2}$", x = colnames(filtered_landscape), value = T)
    top_clustering_col <- paste0("top_",clustering_col)
    clustering_depth <- as.integer(gsub(pattern = "k", x = clustering_col, replacement = ""))
    
    filtered_landscape <- filtered_landscape %>% 
      dplyr::select(type, binID, all_of(clustering_col), 
                    chr, ampl, pos, 
                    all_of(top_clustering_col), prediction)
    
    filtered_landscape$obs <- filtered_landscape$ampl
    filtered_landscape$pred <- filtered_landscape$prediction
    filtered_landscape$ampl <- NULL; filtered_landscape$prediction <- NULL
    
    
  } else {
    
    filtered_shap_abs_sum <- shap_plotting_list$filtered_shap_abs_sum_del
    genome_mask <- filtered_shap_output$genome_mask
    model_mask <- filtered_shap_output$model_mask
    type_mask <- filtered_shap_output$type_mask
    
    pred_df <- pred_list[[model_input$selected]]; pred_df$del_score <- NULL
    
    filtered_landscape <- merge(x = filtered_landscape, by.x = c("type", "binID"),
                                y = pred_df,            by.y = c("Type", "bin"),
                                sort = F
    )
    
    clustering_col <- grep(pattern = "^k\\d{1,2}$", x = colnames(filtered_landscape), value = T)
    top_clustering_col <- paste0("top_",clustering_col)
    clustering_depth <- as.integer(gsub(pattern = "k", x = clustering_col, replacement = ""))
    
    filtered_landscape <- filtered_landscape %>% 
      dplyr::select(type, binID, all_of(clustering_col), 
                    chr, del, pos, 
                    all_of(top_clustering_col), prediction)
    
    filtered_landscape$obs <- filtered_landscape$del
    filtered_landscape$pred <- filtered_landscape$prediction
    filtered_landscape$del <- NULL; filtered_landscape$prediction <- NULL
    
  }
  
  outlist <- list(filtered_landscape = filtered_landscape,
                  genome_mask = genome_mask,
                  model_mask = model_mask,
                  type_mask = type_mask)
  
  return(outlist)
}