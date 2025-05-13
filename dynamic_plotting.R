setwd('/Users/gabry/OneDrive/Desktop/shiny_app/')

library(tidyverse)
library(ggplot2)
library(IRanges) 
library(stringr)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)

shap.list <- readRDS("dev/Data/shap_Mid-length_AmplDel.rds")
toplot.plot <- readRDS("dev/Data/landscape_plot.rds")
clusters_explained <- readRDS("dev/Data/clusters_explained.rds")
centromere_table <- read.table("dev/Data/centomere.tsv", header = T)
load("dev/Data/All_levels_backbonetables.RData")

parse_input_data <- function(shap.list, 
                             toplot.plot,
                             clusters_explained,
                             chr_backbone_namesfixed,
                             centromere_table){
  
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
}

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

filter_df <- function(input_obj, 
                      backbone_granges, 
                      type_input = NULL, 
                      model_input = NULL, 
                      chr_input = NULL,
                      coord_input = NULL){
  
  # model filtering policy:
  # either "ampl" or "del" must be specified
  
  model_mask <- parse_input_model(model_input)
  
  if (class(input_obj) == "list"){
    
    df_input <- input_obj[[model_mask$selected]]

  } else if (class(input_obj) == "data.frame") {
    
    df_input <- input_obj[,colnames(input_obj) != model_mask$not_selected]
      
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

generate_tick_df <- function(input_df, name_annot, mode){
  
  
  if (mode == "ampl") {
    output_df <- input_df %>%
      filter(reason == name_annot) %>%
      mutate(cluster_ymin = max(input_df$ampl) * 1.07 + clusters20 * 0.015,
             cluster_ymax = cluster_ymin + 0.03)
  } else {
    
    output_df <- input_df %>%
      filter(reason == name_annot) %>%
      mutate(cluster_ymax = min(-input_df$del) * 1.07 + clusters20 * 0.015,
             cluster_ymin = cluster_ymax - 0.03)
    
  }
  
  return(output_df)
}

landscape_plot <- function(filtered_landscape_ampl, 
                           filtered_landscape_del, 
                           genome_mask, 
                           type_mask, 
                           model_mask,
                           plot_ampl = TRUE, 
                           plot_del = TRUE,
                           plot_unknown = TRUE, 
                           plot_essential = TRUE, 
                           plot_accessible = TRUE,
                           plot_hiexpr = TRUE, 
                           plot_og_centr_lowmu = TRUE, 
                           plot_active = TRUE,
                           plot_tsg_centr_tel_lowmu = TRUE, 
                           plot_fgs = TRUE,
                           plot_acc_enh_prom_trx_rep_lowexp_himu = TRUE,
                           plot_tsg_fgs_tel = TRUE, 
                           plot_og = TRUE, 
                           plot_rep = TRUE){ 
  
  valid_input <- c("ampl","del")
  is_valid <- all(model_mask %in% valid_input)
  
  if (length(genome_mask) == 22) {
    genome_mask <- "WHOLE GENOME"
  }
  
  if (length(genome_mask) > 1) {
    genome_mask <- paste(genome_mask, collapse = ", ")
  } 
  
  if (is_valid) {
    
    if (length(model_mask) > 1) {
      
      model_mask <- paste(model_mask, collapse = ", ")
      
    }
  } else {
    
    stop("Invalid model selected. \n Models are either \"ampl\" or \"del\"")
  
  }
  
  title <- "Segment Annotation (based on SHAP values)"
  subtitle <- paste0(paste0(" [", genome_mask, "] "), 
                     paste0(" [", type_mask, "] "), 
                     paste0(" [", model_mask, "] "))
  
  filtered_landscape_ampl <- filtered_landscape_ampl %>% mutate(pos = 1:n())
  filtered_landscape_del <- filtered_landscape_del %>% mutate(pos = 1:n())
  
  chr_bounds <- filtered_landscape_ampl %>%
    group_by(chr) %>%
    summarize(start = min(pos), end = max(pos), .groups = "drop") %>%
    mutate(chr_num = readr::parse_number(chr)) %>%
    arrange(chr_num) %>%
    mutate(fill = rep(c("white", "#e7deed"), length.out = n())) %>%
    select(-chr_num)
  
  color_palette_background <- c("white", "#e7deed")
  
  base_plot <- ggplot() +
    geom_rect(data = chr_bounds,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.3) +
    scale_fill_manual(values = color_palette_background)
  
  if (plot_ampl) {
    
    chr_to_plot <- unique(filtered_landscape_ampl$chr) 
    
    for (chr in chr_to_plot) {
      base_plot <- base_plot + geom_line(data = filtered_landscape_ampl[filtered_landscape_ampl$chr == chr,], 
                                         aes(x = pos, y = ampl), 
                                         color = "red")
    }
  }
  
  if (plot_del) {
    
    chr_to_plot <- unique(filtered_landscape_ampl$chr) 
    
    for (chr in chr_to_plot) {
      base_plot <- base_plot + geom_line(data = filtered_landscape_del[filtered_landscape_del$chr == chr,], 
                                         aes(x = pos, y = -del), 
                                         color = "blue")
    }
  }
  
  base_plot <- base_plot +
    geom_hline(yintercept = 0, colour = 'grey', linetype = 'dashed', linewidth = 0.2) +
    labs(title = title, 
         subtitle = subtitle, 
         x = "Genomic Position",
         y = "SCNA frequency (Mid-length)") +
    theme_classic() +
    theme(legend.position = 'none')
  
  color_palette_ticks <- c(
    "Unknown" = "#666666",
    "Essential" = "#cc0000",
    "Accessible / Low Expression / High Mu" = "#0000cc",
    "High Expression" = "#007700",
    "OG / Centromere / Low Mu" = "#800080",
    "Enhancer / Promoters / Transcribed = ACTIVE" = "#ff8000",
    "TSG / Centromere / Telomere / Low Mu" = "#999900",
    "FGS" = "#00aaaa",
    "Accessible / Enhancers / Promoters / Transcribed / Repressed / Low Expression / High Mu (?)" = "#ff66cc",
    "TSG / FGS / Telomere" = "#8b4513",
    "OG" = "#cc00cc",
    "Repressed" = "#3399cc"
  )
  
  
  ticksize <- 0.1
  
  if (plot_unknown) {
    
    cluster_ticks_unknown <- generate_tick_df(filtered_landscape_ampl, 
                                              name_annot = names(color_palette_ticks)[1], 
                                              mode = "ampl")
    
    base_plot <- base_plot +
      geom_segment(data = cluster_ticks_unknown,
                   aes(x = pos, xend = pos, 
                       y = cluster_ymin, 
                       yend = cluster_ymax,
                       color = reason),
                   linewidth = ticksize)
  }
  if (plot_essential) {
    
    cluster_ticks_essential <- generate_tick_df(filtered_landscape_ampl, 
                                                name_annot = names(color_palette_ticks)[2], 
                                                mode = "ampl")
    
    base_plot <- base_plot +
      geom_segment(data = cluster_ticks_essential,
                   aes(x = pos, xend = pos, 
                       y = cluster_ymin, 
                       yend = cluster_ymax,
                       color = reason),
                   linewidth = ticksize)
  }
  if (plot_accessible) {
    
    cluster_ticks_accessible <- generate_tick_df(filtered_landscape_ampl, 
                                                 name_annot = names(color_palette_ticks)[3], 
                                                 mode = "ampl")
    
    base_plot <- base_plot +
      geom_segment(data = cluster_ticks_accessible,
                   aes(x = pos, xend = pos, 
                       y = cluster_ymin, 
                       yend = cluster_ymax,
                       color = reason),
                   linewidth = ticksize)
  }
  if (plot_hiexpr) {
    
    cluster_ticks_hiexpr <- generate_tick_df(filtered_landscape_ampl, 
                                             name_annot = names(color_palette_ticks)[4], 
                                             mode = "ampl")
    
    base_plot <- base_plot +
      geom_segment(data = cluster_ticks_hiexpr,
                   aes(x = pos, xend = pos, 
                       y = cluster_ymin, 
                       yend = cluster_ymax,
                       color = reason),
                   linewidth = ticksize)
  }
  if (plot_og_centr_lowmu) {
    
    cluster_ticks_og_centr_lowmu <- generate_tick_df(filtered_landscape_del, 
                                                     name_annot = names(color_palette_ticks)[5], 
                                                     mode = "del")
    

    base_plot <- base_plot +
      geom_segment(data = cluster_ticks_og_centr_lowmu,
                   aes(x = pos, xend = pos, 
                       y = cluster_ymin, 
                       yend = cluster_ymax,
                       color = reason),
                   linewidth = ticksize)
  }
  if (plot_active) {
    
    cluster_ticks_active <- generate_tick_df(filtered_landscape_ampl, 
                                             name_annot = names(color_palette_ticks)[6], 
                                             mode = "ampl")
    
    base_plot <- base_plot +
      geom_segment(data = cluster_ticks_active,
                   aes(x = pos, xend = pos, 
                       y = cluster_ymin, 
                       yend = cluster_ymax,
                       color = reason),
                   linewidth = ticksize)
  }
  if (plot_tsg_centr_tel_lowmu) {
    
    cluster_ticks_tsg_centr_tel_lowmu <- generate_tick_df(filtered_landscape_del, 
                                                          name_annot = names(color_palette_ticks)[7], 
                                                          mode = "del")
    
    base_plot <- base_plot +
      geom_segment(data = cluster_ticks_tsg_centr_tel_lowmu,
                   aes(x = pos, xend = pos, 
                       y = cluster_ymin, 
                       yend = cluster_ymax,
                       color = reason),
                   linewidth = ticksize)
  }
  if (plot_fgs) {
    
    cluster_ticks_fgs <- generate_tick_df(filtered_landscape_del, 
                                          name_annot = names(color_palette_ticks)[8], 
                                          mode = "del")
    
    base_plot <- base_plot +
      geom_segment(data = cluster_ticks_fgs,
                   aes(x = pos, xend = pos, 
                       y = cluster_ymin, 
                       yend = cluster_ymax,
                       color = reason),
                   linewidth = ticksize)
  }
  if (plot_acc_enh_prom_trx_rep_lowexp_himu) {
    
    cluster_ticks_acc_enh_prom_trx_rep_lowexp_himu <- generate_tick_df(filtered_landscape_ampl, 
                                                                       name_annot = names(color_palette_ticks)[9], 
                                                                       mode = "ampl")
    
    base_plot <- base_plot +
      geom_segment(data = cluster_ticks_acc_enh_prom_trx_rep_lowexp_himu,
                   aes(x = pos, xend = pos, 
                       y = cluster_ymin, 
                       yend = cluster_ymax,
                       color = reason),
                   linewidth = ticksize)
  }
  if (plot_tsg_fgs_tel) {
    
    cluster_ticks_tsg_fgs_tel <- generate_tick_df(filtered_landscape_del, 
                                                  name_annot = names(color_palette_ticks)[10], 
                                                  mode = "del")
    
    base_plot <- base_plot +
      geom_segment(data = cluster_ticks_tsg_fgs_tel,
                   aes(x = pos, xend = pos, 
                       y = cluster_ymin, 
                       yend = cluster_ymax,
                       color = reason),
                   linewidth = ticksize)
  }
  if (plot_og) {
    
    cluster_ticks_og <- generate_tick_df(filtered_landscape_ampl, 
                                         name_annot = names(color_palette_ticks)[11], 
                                         mode = "ampl")
    
    base_plot <- base_plot +
      geom_segment(data = cluster_ticks_og,
                   aes(x = pos, xend = pos, 
                       y = cluster_ymin, 
                       yend = cluster_ymax,
                       color = reason),
                   linewidth = ticksize)
  }
  if (plot_rep) {
    
    cluster_ticks_rep <- generate_tick_df(filtered_landscape_ampl, 
                                          name_annot = names(color_palette_ticks)[12], 
                                          mode = "ampl")
    
    base_plot <- base_plot +
      geom_segment(data = cluster_ticks_rep,
                   aes(x = pos, xend = pos, 
                       y = cluster_ymin, 
                       yend = cluster_ymax,
                       color = reason),
                   linewidth = ticksize)
  }
  
  upper_limit <- ceiling(max(filtered_landscape_ampl$ampl) * 10) / 10
  lower_limit <- floor(min(-filtered_landscape_del$del) * 10) / 10
  
  sym_limit <- min(abs(upper_limit), abs(lower_limit))
  
  y_breaks <- pretty(c(-sym_limit, sym_limit))
  
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
      breaks = chr_bounds %>% mutate(center = (start + end)/2) %>% pull(center),
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
    scale_color_manual(values = color_palette_ticks)
  
  base_plot
}

processed_data <- parse_input_data(shap.list = shap.list, 
                                   toplot.plot = toplot.plot,
                                   clusters_explained = clusters_explained,
                                   chr_backbone_namesfixed = chr_backbone_namesfixed, 
                                   centromere_table = centromere_table)

shap.list <- processed_data$shap.list; 
toplot.plot <- processed_data$toplot.plot; 
backbone.100kb <- processed_data$backbone.100kb
centromere_table <- processed_data$centromere_table

type_input <- "BRCA"; 
model_input_ampl <- "ampl"; model_input_del <- "del"
coord_input <- NULL; chr_input <- NULL

filtered_shap_output_ampl <- filter_df(input_obj = shap.list, backbone_granges = backbone.100kb, 
                                       type_input = type_input, model_input = model_input_ampl, 
                                       chr_input = chr_input, coord_input = coord_input)

filtered_shap_output_del <- filter_df(input_obj = shap.list, backbone_granges = backbone.100kb, 
                                      type_input = type_input, model_input = model_input_del, 
                                      chr_input = chr_input, coord_input = coord_input)


filtered_landscape_output_ampl <- filter_df(input_obj = toplot.plot, backbone_granges = backbone.100kb,
                                            type_input = type_input, model_input = model_input_ampl, 
                                            chr_input = chr_input, coord_input = coord_input)

filtered_landscape_output_del <- filter_df(input_obj = toplot.plot, backbone_granges = backbone.100kb,
                                           type_input = type_input, model_input = model_input_del,  
                                           chr_input = chr_input, coord_input = coord_input)


filtered_shap_ampl <- filtered_shap_output_ampl$final_df; filtered_shap_del <- filtered_shap_output_del$final_df
filtered_landscape_ampl <- filtered_landscape_output_ampl$final_df; filtered_landscape_del <- filtered_landscape_output_del$final_df

genome_mask_ampl <- filtered_shap_output_ampl$genome_mask; genome_mask_del <- filtered_shap_output_del$genome_mask
model_mask_ampl <- filtered_shap_output_ampl$model_mask; model_mask_del <- filtered_shap_output_del$model_mask
type_mask_ampl <- filtered_shap_output_ampl$type_mask; type_mask_del <- filtered_shap_output_del$type_mask

shap_plotting_list <- prepare_shap_to_plot(filtered_shap_ampl = filtered_shap_ampl, 
                                           filtered_shap_del = filtered_shap_del)

filtered_shap_abs_sum_ampl <- shap_plotting_list$filtered_shap_abs_sum_ampl
filtered_shap_abs_sum_del <- shap_plotting_list$filtered_shap_abs_sum_del

barplot_shap(shap.abs.sum = filtered_shap_abs_sum_ampl, genome_mask = genome_mask_ampl, type_mask = type_mask_ampl, model_mask = model_mask_ampl)
barplot_shap(shap.abs.sum = filtered_shap_abs_sum_del, genome_mask = genome_mask_del, type_mask = type_mask_del, model_mask = model_mask_del)

landscape_plot(filtered_landscape_ampl = filtered_landscape_ampl, 
               filtered_landscape_del = filtered_landscape_del, 
               genome_mask = genome_mask_ampl, 
               type_mask = type_mask_ampl, 
               model_mask = c("ampl","del"),
               plot_ampl = TRUE, 
               plot_del = TRUE,
               plot_unknown = FALSE, 
               plot_essential = TRUE, 
               plot_accessible = TRUE, 
               plot_hiexpr = TRUE, 
               plot_og_centr_lowmu = TRUE, 
               plot_active = TRUE, 
               plot_tsg_centr_tel_lowmu = TRUE ,
               plot_fgs = TRUE, 
               plot_acc_enh_prom_trx_rep_lowexp_himu = TRUE, 
               plot_tsg_fgs_tel = TRUE, 
               plot_og = TRUE, 
               plot_rep = TRUE)

