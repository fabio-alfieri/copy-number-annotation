if(F){parse_input_data <- function(shap.list, toplot.plot, clusters_explained, 
                                   chr_backbone_namesfixed, centromere_table){
  
  shap_clean_list <- list()
  
  shap_cols_to_discard <- c("Chromosome_Length", "Centromere_Length", "Centromere_Type", 'BIAS')
  shap_cols_to_keep <- c("labels", "type", "chr", "bin",
                         "dist.to.closest.OG", "dist.to.closest.TSG", "dist.to.closest.FGS",
                         "distance.to.centromere", "distance.to.telomere", "mutations_norm",
                         "genes.bin", "Length_Counts.E17", "Length_Counts.E19", "Ess.distance_pancancer",
                         "Length_Counts.E25", "Length_Counts.E1")
  
  for (idx in seq_along(along.with = shap.list)) {
    
    name <- names(shap.list)[idx]
    shap.df <- shap.list[[idx]]
    
    shap.df <- shap.df[,apply(shap.df, 2, function(x){all(x!=0)})]
    
    shap.df$type <- do.call(rbind, str_split(shap.df$labels, pattern = '-'))[,2]
    shap.df$chr <- do.call(rbind, str_split(shap.df$labels, '_'))[,1]
    shap.df$bin <- do.call(rbind, str_split(do.call(rbind, str_split(shap.df$labels, '_'))[,2], '-'))[,1]
    
    shap.df <- shap.df %>% select(-all_of(shap_cols_to_discard))
    shap.df <- shap.df %>% select(all_of(shap_cols_to_keep))
    
    
    shap.df$binID <- paste0(shap.df$chr, "_", shap.df$bin)
    shap.df$chr <- paste0("chr", shap.df$chr); shap.df$bin <- NULL
    
    shap_clean_list[[name]] <- shap.df
    
  }
  
  toplot.plot$binID <- paste0(toplot.plot$chr, "_", toplot.plot$bin)
  toplot.plot$chr <- paste0("chr", toplot.plot$chr); toplot.plot$bin <- NULL
  colnames(toplot.plot)[c(3,4)] <- c("ampl", "del")
  toplot.plot$clusters10 <- NULL; toplot.plot$clusters15 <- NULL
  
  toplot.plot$order <- seq_along(1:nrow(toplot.plot))
  
  toplot.plot <- merge(x = toplot.plot, 
                       y = clusters_explained, 
                       by = "clusters20", 
                       all.x = TRUE)
  
  toplot.plot <- toplot.plot[order(toplot.plot$order), ]; toplot.plot$order <- NULL
  
  toplot.plot[which(is.na(toplot.plot$reason)),]$reason <- "Unknown"
  toplot.plot[which(toplot.plot$reason == "Unknown"),]$clusters20 <- max(toplot.plot[which(toplot.plot$reason == "Unknown"),]$clusters20)
  
  to_flip <- c(1,7,9,10)
  clusters_to_flip <- unique(toplot.plot$clusters20)[to_flip]
  toplot.plot[toplot.plot$clusters20 %in% clusters_to_flip, ]$clusters20 <- 
    -toplot.plot[toplot.plot$clusters20 %in% clusters_to_flip, ]$clusters20
  
  positive_clusters <- sort(unique(toplot.plot[sign(toplot.plot$clusters20) == 1, ]$clusters20))
  negative_clusters <- sort(unique(toplot.plot[sign(toplot.plot$clusters20) == -1, ]$clusters20))
  
  lapply(X = seq_along(positive_clusters), FUN = function(idx){
    curr_cluster <- positive_clusters[idx]
    toplot.plot[toplot.plot$clusters20 == curr_cluster, ]$clusters20 <<- idx
  })
  
  lapply(X = seq_along(negative_clusters), FUN = function(idx){
    curr_cluster <- negative_clusters[idx]
    toplot.plot[toplot.plot$clusters20 == curr_cluster, ]$clusters20 <<- -idx
  })
  
  
  backbone.100kb <- chr_backbone_namesfixed$`0.1Mbp`; backbone.100kb <- dplyr::bind_rows(backbone.100kb)
  backbone.100kb$binID <- paste0(backbone.100kb$chr, "_", backbone.100kb$bin)
  backbone.100kb$chr <- paste0("chr", backbone.100kb$chr); backbone.100kb$bin <- NULL
  
  backbone.100kb <- GenomicRanges::GRanges(seqnames = backbone.100kb$chr, 
                                           ranges = IRanges(start = backbone.100kb$start_bin, end = backbone.100kb$end_bin),
                                           binID = backbone.100kb$binID)
  
  centromere_table$Centromere_Length <- NULL
  centromere_table$Centromere_Type <- NULL
  centromere_table$Centromere <- NULL
  centromere_table$chr <- paste0("chr", centromere_table$chr)
  centromere_table$midpoint <- (centromere_table$end - centromere_table$start) / 2
  centromere_table$midpoint_start <- centromere_table$midpoint
  centromere_table$midpoint_end <- centromere_table$midpoint
  
  centromere_gr <- GRanges(seqnames = centromere_table$chr, 
                           ranges = IRanges(start = centromere_table$midpoint_start, 
                                            end = centromere_table$midpoint_end))
  
  hits <- findOverlaps(query = centromere_gr, subject = backbone.100kb)
  centromere_table_out <- data.frame(backbone.100kb[subjectHits(hits)])
  centromere_table_out$width <- NULL; centromere_table_out$strand <- NULL
  colnames(centromere_table_out) <- c("chr", "start", "end", "binID")
  
  
  outlist <- list(shap.list = shap_clean_list,
                  toplot.plot = toplot.plot,
                  backbone.100kb = backbone.100kb, 
                  centromere_table = centromere_table_out)
  
  return(outlist)
}} # old parse_input_data()

parse_input_coord <- function(input){
  
  chrom_sizes <- seqlengths(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)[1:22]
  valid_input <- "^chr(?:[1-9]|1[0-9]|2[0-2]):\\d+-\\d+$"
  
  is_valid <- grepl(pattern = valid_input, x = input, ignore.case = F)
  
  if (is_valid) {
    
    elems <- strsplit(x = input, split = "[:-]")[[1]]
    chr <- elems[1]; start <- as.numeric(elems[2]); end <- as.numeric(elems[3])
    corresponding_chr_size <- chrom_sizes[chr]

    if (end > corresponding_chr_size) {
      
      stop("END coordinate > than chromosome size. \nENTER a valid END coordinate")
      
    }
    
  } else {
    
    stop("INVALID coordinate PATTERN. \nENTER a valid PATTERN. \nExample: chr:start-end")
    
  }
  
  input_granges <- GenomicRanges::GRanges(seqnames = chr, 
                                          ranges = IRanges(start = start, 
                                                           end = end))
  
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
  
  if ((length(input) == 1) && (input == "all")) {
    
    outlist <- list(selected = c("ampl","del"),
                    not_selected = NA)
    
    return(outlist)
  } 
  
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
parse_clustering_depth <- function(input){
  
  valid_input <- c(1,2,3,4)
  valid_input_chr <- c("1","2","3","4")
  
  is_valid <- ((input %in% valid_input) || (input %in% valid_input_chr))
  
  if (is_valid) {
    
    outlist <- list(selected = as.integer(input),
                    not_selected = valid_input[!(valid_input %in% as.integer(input))])
    
    return(outlist)
    
  } else {
    
    stop("Invalid clustering depth selected. Values must be 1, 2, 3 or 4")
    
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

filter_df <- function(input_obj, backbone_granges, type_input = NULL, model_input = NULL, chr_input = NULL, coord_input = NULL){
  
  # model filtering policy:
  # either "ampl" or "del" must be specified

  model_mask <- parse_input_model(model_input)
  
  if (any(class(input_obj) == "list")) {
    
    if (length(model_mask$selected) == 2) {
      stop("SHAP list must be explicitly filtered")       
    }
    
    df_input <- input_obj[[model_mask$selected]]

  } else if (any(class(input_obj) == "data.frame")) {
    
    df_input <- input_obj
      
  } else {
    
    stop("input_obj MUST either be a \"data.frame\" or a \"list\"") 
  
  }

  # Region filtering policy:
  # if both chromosome and coordinate level are NULL --> plot whole genome
  # if both are defined --> plot the most detailed one (coordinates)
  
  if ((!is.null(chr_input)) && (!is.null(coord_input))) {
    
    genome_mask <- parse_input_coord(coord_input)
    
    hits <- GenomicRanges::findOverlaps(query = genome_mask, 
                                        subject = backbone_granges)
    
    overlapping_bins <- values(backbone_granges[subjectHits(hits)])
    
    coord_filtered_df <- df_input[df_input$binID %in% overlapping_bins$binID,]
    
  } else if ((is.null(chr_input)) && (is.null(coord_input))) {

    genome_mask <- unique(df_input$chr)
    coord_filtered_df <- df_input[df_input$chr %in% genome_mask,]
    
  } else if ((!is.null(chr_input)) && (is.null(coord_input))){
    
    genome_mask <- parse_input_chr(chr_input)
    coord_filtered_df <- df_input[df_input$chr %in% genome_mask,]
    
  } else {
    
    genome_mask <- parse_input_coord(coord_input)
    hits <- GenomicRanges::findOverlaps(query = genome_mask, 
                                        subject = backbone_granges)
    
    overlapping_bins <- values(backbone_granges[subjectHits(hits)])
    coord_filtered_df <- df_input[df_input$binID %in% overlapping_bins$binID,]
    
  }
  
  # Type filtering policy:
  # Only one cancer type to specify
  # Must be in the 11 cancer types
  
  type_mask <- parse_input_type(type_input)
  
  type_filtered_df <- coord_filtered_df[coord_filtered_df$type == type_mask,]
  
  outlist <- list(final_df = type_filtered_df,
                  model_mask = model_mask$selected,
                  type_mask = type_mask,
                  genome_mask = as.character(genome_mask))

  return(outlist)
  
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

landscape_plot_interactive <- function(filtered_landscape,
                                       backbone.100kb,
                                       genome_mask, type_mask, model_mask,
                                       plot_ampl = TRUE, plot_del = TRUE,
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
  plot_ampl_layer <- function(base_plot, chr_to_plot, filtered_landscape, backbone.100kb){
    
    bg_ampl <- "#FF0000"; fg_ampl <- get_contrast(bg_ampl)
    
    for (chr in chr_to_plot) {
      chr_data <- filtered_landscape[filtered_landscape$chr == chr, ]
      chr_data <- chr_data %>% 
        rowwise() %>% 
        mutate(
        coord = as.character(backbone.100kb[mcols(backbone.100kb)$binID == binID][1]),
        data_id = binID
      )
      base_plot <- base_plot +
        geom_line(data = chr_data, aes(x = pos, y = ampl), color = "red") +
        geom_point_interactive(
          data = chr_data,
          aes(x = pos, y = ampl, 
              tooltip = sprintf(
                "<div style='background:%s; 
                  color:%s; 
                  padding:4px; 
                  border-radius:0px; 
                  border:none; 
                  outline:none; 
                  box-shadow:none;'>
                  Coords: %s<br>
                  Ampl: %s</div>",
                bg_ampl, fg_ampl, coord, round(ampl,3)
              ), 
              data_id = data_id),
          size = 3, color = "transparent"
        )
    }
    return(base_plot)
  }
  plot_del_layer <- function(base_plot, chr_to_plot, filtered_landscape, backbone.100kb){
    
    bg_del <- "#0000FF"; fg_del <- get_contrast(bg_del)
    for (chr in chr_to_plot) {
      chr_data <- filtered_landscape[filtered_landscape$chr == chr, ]
      chr_data <- chr_data %>% 
        rowwise() %>% 
        mutate(
        coord = as.character(backbone.100kb[mcols(backbone.100kb)$binID == binID][1]),
        data_id = binID
      )
      base_plot <- base_plot +
        geom_line(data = chr_data, aes(x = pos, y = -del), color = "blue") +
        geom_point_interactive(
          data = chr_data,
          aes(x = pos, y = -del, 
              tooltip = sprintf(
                "<div style='background:%s; 
                  color:%s; 
                  padding:4px; 
                  border-radius:0px; 
                  border:none; 
                  outline:none; 
                  box-shadow:none;'>
                  Coords: %s<br>
                  Del: %s</div>",
                bg_del, fg_del, coord, round(del,3)
              ), data_id = data_id),
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
                                linewidth, color_palette_ticks) {
    
    
    generate_density_df <- function(input_df, 
                                    name_annot,mode, 
                                    clustering_col, top_clustering_col, clustering_depth, 
                                    backbone.100kb,
                                    color_palette_ticks) {
      
      height <- 0.015
      bg <- color_palette_ticks[as.character(name_annot)]
      fg <- get_contrast(bg)
      
      df <- input_df %>%
        filter(.data[[clustering_col]] == name_annot)
      
      if (nrow(df) == 0) return(NULL)
      
      if (mode == "ampl") {
        start <- max(input_df$ampl)
        df <- df %>%
          rowwise() %>%
          mutate(
            tooltip = sprintf(
              "<div style='background:%s; color:%s; padding:4px;'>%s</div>",
              bg, fg, .data[[top_clustering_col]]
            ),
            cluster_ymid = (round(start, 1) + 0.05) + ((.data[[clustering_col]] / (clustering_depth / 3.8)) * 0.1),
            cluster_ymin = cluster_ymid - height,
            cluster_ymax = cluster_ymid + height
          ) %>%
          ungroup()
      } else {
        start <- min(-input_df$del) - 0.05
        df <- df %>%
          rowwise() %>%
          mutate(
            coord = as.character(backbone.100kb[mcols(backbone.100kb)$binID == binID][1]),
            tooltip = sprintf(
              "<div style='background:%s; color:%s; padding:4px;'>%s</div>",
              bg, fg, .data[[top_clustering_col]]
            ),
            cluster_ymid = (round(start, 1) - 0.05) - ((.data[[clustering_col]] / (clustering_depth / 3.8)) * 0.1),
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
                                           color_palette_ticks = color_palette_ticks)
    
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
      alpha = 0.4
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
                                ticksize, color_palette_ticks) {
    
    
    generate_tick_df <- function(input_df, 
                                 name_annot,mode, 
                                 clustering_col, top_clustering_col, clustering_depth, 
                                 backbone.100kb,
                                 color_palette_ticks) {
      
      height <- 0.015
      bg <- color_palette_ticks[as.character(name_annot)]
      fg <- get_contrast(bg)
      
      df <- input_df %>%
        filter(.data[[clustering_col]] == name_annot)
      
      if (nrow(df) == 0) return(NULL)
      
      if (mode == "ampl") {
        start <- max(input_df$ampl)
        df <- df %>%
          rowwise() %>%
          mutate(
            coord = as.character(backbone.100kb[mcols(backbone.100kb)$binID == binID][1]),
            tooltip = sprintf(
              "<div style='background:%s; color:%s; padding:4px;'>Coords: %s<br>%s</div>",
              bg, fg, coord, .data[[top_clustering_col]]
            ),
            cluster_ymid = (round(start, 1) + 0.05) + ((.data[[clustering_col]] / (clustering_depth / 3.8)) * 0.1),
            cluster_ymin = cluster_ymid - height,
            cluster_ymax = cluster_ymid + height
          ) %>%
          ungroup()
      } else {
        start <- min(-input_df$del) - 0.05
        df <- df %>%
          rowwise() %>%
          mutate(
            coord = as.character(backbone.100kb[mcols(backbone.100kb)$binID == binID][1]),
            tooltip = sprintf(
              "<div style='background:%s; color:%s; padding:4px;'>Coords: %s<br>%s</div>",
              bg, fg, coord, .data[[top_clustering_col]]
            ),
            cluster_ymid = (round(start, 1) - 0.05) - ((.data[[clustering_col]] / (clustering_depth / 3.8)) * 0.1),
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
                                      color_palette_ticks = color_palette_ticks)
    
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
  
  if (!all(model_mask %in% valid_input)) stop("Invalid model selected. Use 'ampl' and/or 'del'.")
  if (length(genome_mask) == 22) genome_mask <- "WHOLE GENOME"
  if (length(genome_mask) > 1) genome_mask <- paste(genome_mask, collapse = ", ")
  if (length(model_mask) > 1) model_mask <- paste(model_mask, collapse = ", ")
  
  title    <- "Segment Annotation (based on SHAP values)"
  subtitle <- paste0("[", genome_mask, "] [", type_mask, "] [", model_mask, "]")
  
  filtered_landscape <- filtered_landscape %>% mutate(pos = row_number())
  
  chr_bounds <- get_chr_bounds(filtered_landscape = filtered_landscape)
  
  base_plot <- plot_base_layer(chr_bounds = chr_bounds, filtered_landscape = filtered_landscape)
  
  chr_to_plot <- unique(filtered_landscape$chr)
  
  if (plot_ampl) {
    base_plot <- plot_ampl_layer(base_plot = base_plot, 
                                  chr_to_plot = chr_to_plot, 
                                  filtered_landscape = filtered_landscape, 
                                  backbone.100kb = backbone.100kb)
    }
  
  if (plot_del) {
    base_plot <- plot_del_layer(base_plot = base_plot, 
                   chr_to_plot = chr_to_plot, 
                   filtered_landscape = filtered_landscape, 
                   backbone.100kb = backbone.100kb)
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
  
  all_modes <- c("ampl", "del")
  
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
                                   color_palette_ticks = color_palette_ticks
     )
   }
  }
  
  if (!isFALSE(annot_to_plot_kde)) {
    
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
                                   color_palette_ticks = color_palette_ticks
    )
   }
  }
                                            
  upper_limit <- ceiling(max(filtered_landscape$ampl) * 10) / 10
  lower_limit <- floor(min(-filtered_landscape$del) * 10) / 10
  sym_limit   <- min(abs(upper_limit), abs(lower_limit))
  y_breaks    <- pretty(c(-sym_limit, sym_limit))
  
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
    coord_cartesian(ylim = c(-1.2, 1.2)) +
    theme(
      axis.line.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) + 
    scale_fill_manual(values = color_palette_ticks) +
    scale_colour_manual(values = color_palette_ticks)
  
  
  girafe(
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
}

