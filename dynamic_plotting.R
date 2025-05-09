setwd('/Users/gabry/OneDrive/Desktop/copy-number-annotation-main/') # set this folder as directory

# DO NOT TOUCH
##############################################

# useful notes/data for shiny app
library(tidyverse)
library(ggplot2)
library(IRanges)

load('dev/Data/shap.RData')

# remove features that are always 0 across columns
shap.df <- shap.df[,apply(shap.df, 2, function(x){all(x!=0)})]

shap.df$type <- do.call(rbind, str_split(shap.df$labels, pattern = '-'))[,2]
shap.df$chr <- do.call(rbind, str_split(shap.df$labels, '_'))[,1]
shap.df$bin <- do.call(rbind, str_split(do.call(rbind, str_split(shap.df$labels, '_'))[,2], '-'))[,1]

# remove all non-segment features and useless columns
shap.df <- shap.df %>% select(-c("Chromosome_Length", "Centromere_Length", "Centromere_Type", 'BIAS'))
shap.df <- shap.df %>% select("labels", "type", "chr", "bin",
                              "dist.to.closest.OG", "dist.to.closest.TSG", "dist.to.closest.FGS",
                              "distance.to.centromere", "distance.to.telomere", "mutations_norm",
                              "genes.bin", "Length_Counts.E17", "Length_Counts.E19", "Ess.distance_pancancer",
                              "Length_Counts.E25", "Length_Counts.E1")

View(shap.df) # this is the table you should use for plotting sum(abs(SHAP)) for each point

saveRDS(object = shap.df, file = "dev/Data/shap_df_clean_NOT_DEF.rds")

##############################################

shap.df <- readRDS("dev/Data/shap_df_clean_NOT_DEF.rds")
toplot.plot <- readRDS('dev/Data/landscape_plot.rds')
load("dev/Data/All_levels_backbonetables.RData")

parse_input_data <- function(shap.df, toplot.plot, chr_backbone_namesfixed){
  
  shap.df$binID <- paste0(shap.df$chr, "_", shap.df$bin)
  shap.df$chr <- paste0("chr", shap.df$chr); shap.df$bin <- NULL

  shap.list <- list(ampl = shap.df, del = shap.df) # CHANGE THIS

  toplot.plot$binID <- paste0(toplot.plot$chr, "_", toplot.plot$bin)
  toplot.plot$chr <- paste0("chr", toplot.plot$chr); toplot.plot$bin <- NULL
  colnames(toplot.plot)[c(3,4)] <- c("ampl", "del")
  
  
  backbone.100kb <- chr_backbone_namesfixed$`0.1Mbp`; backbone.100kb <- dplyr::bind_rows(backbone.100kb)
  backbone.100kb$binID <- paste0(backbone.100kb$chr, "_", backbone.100kb$bin)
  backbone.100kb$chr <- paste0("chr", backbone.100kb$chr); backbone.100kb$bin <- NULL
  
  backbone.100kb <- GenomicRanges::GRanges(seqnames = backbone.100kb$chr, 
                                           ranges = IRanges(start = backbone.100kb$start_bin, end = backbone.100kb$end_bin),
                                           binID = backbone.100kb$binID)
  
  outlist <- list(shap.list = shap.list,
                  toplot.plot = toplot.plot,
                  backbone.100kb = backbone.100kb)
  
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

filter_df <- function(input_obj, backbone_granges, type_input = NULL, model_input = NULL, chr_input = NULL, coord_input = NULL){
  
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
landscape_plot <- function(filtered_landscape_ampl, filtered_landscape_del, genome_mask, type_mask, model_mask){ 
  
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
      plot_ampl <- TRUE
      plot_del <- TRUE
      
    } else if (length(model_mask) == 1) {
      
      if (model_mask == "ampl") {
        
        plot_ampl <- TRUE
        plot_del <- FALSE
        
      } else {
        
        plot_ampl <- FALSE
        plot_del <- TRUE
        
      }
    }
    
  } else {
    
    stop("Invalid model selected. \n Models are either \"ampl\" or \"del\"")
  
  }
  
  title <- "SHAP Value Contribution per Feature"
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
    mutate(fill = rep(c("lightgrey", "darkgrey"), length.out = n())) %>%
    select(-chr_num)
  
  base_plot <- ggplot() +
    geom_rect(data = chr_bounds,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.3) +
    scale_fill_manual(values = c("#eeeeee", "#cccccc"))
  
  if (plot_ampl) {
    base_plot <- base_plot + geom_line(data = filtered_landscape_ampl, aes(x = pos, y = ampl), color = "red")
  }
  
  if (plot_del) {
    base_plot <- base_plot + geom_line(data = filtered_landscape_del, aes(x = pos, y = -del), color = "blue")
  }
  
  base_plot <- base_plot +
    geom_hline(yintercept = 0, colour = 'grey', linetype = 'dashed') +
    labs(title = title, 
         subtitle = subtitle, 
         x = "Genomic Position", 
         y = "Amplification Score (Mid-length)") +
    theme_classic() +
    theme(legend.position = 'none')
  
  cluster_ticks <- filtered_landscape_ampl %>%
    mutate(cluster_ymin = max(filtered_landscape_ampl$ampl) * 1.07 + clusters15 * 0.015,
           cluster_ymax = cluster_ymin + 0.01)
  
  base_plot +
    geom_segment(data = cluster_ticks,
                 aes(x = pos, xend = pos, y = cluster_ymin, yend = cluster_ymax),
                 color = "black", size = 0.3)
  
}

processed_data <- parse_input_data(shap.df = shap.df, toplot.plot = toplot.plot, chr_backbone_namesfixed = chr_backbone_namesfixed)

shap.list <- processed_data$shap.list; toplot.plot <- processed_data$toplot.plot; backbone.100kb <- processed_data$backbone.100kb

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

to_clean <- c("labels", "type", "chr", "binID"); mask <- !(colnames(filtered_shap_ampl) %in% to_clean)

filtered_shap_abs_sum_ampl <- data.frame(value = apply(X = filtered_shap_ampl[,mask], 
                                                       MARGIN = 2,
                                                       FUN = function(x){sum(abs(x))}))

filtered_shap_abs_sum_del <- data.frame(value = apply(X = filtered_shap_del[,mask], 
                                                       MARGIN = 2,
                                                       FUN = function(x){sum(abs(x))}))

filtered_shap_abs_sum_ampl$feature <- rownames(filtered_shap_abs_sum_ampl)
filtered_shap_abs_sum_ampl$color <- rainbow(nrow(filtered_shap_abs_sum_ampl))

filtered_shap_abs_sum_del$feature <- rownames(filtered_shap_abs_sum_del)
filtered_shap_abs_sum_del$color <- rainbow(nrow(filtered_shap_abs_sum_del))

genome_mask_ampl <- filtered_shap_output_ampl$genome_mask; genome_mask_del <- filtered_shap_output_del$genome_mask
model_mask_ampl <- filtered_shap_output_ampl$model_mask; model_mask_del <- filtered_shap_output_del$model_mask
type_mask_ampl <- filtered_shap_output_ampl$type_mask; type_mask_del <- filtered_shap_output_del$type_mask

barplot_shap(shap.abs.sum = filtered_shap_abs_sum_ampl, 
             genome_mask = genome_mask_ampl, 
             type_mask = type_mask_ampl, 
             model_mask = model_mask_ampl)

barplot_shap(shap.abs.sum = filtered_shap_abs_sum_del, 
             genome_mask = genome_mask_del, 
             type_mask = type_mask_del, 
             model_mask = model_mask_del)

landscape_plot(filtered_landscape_ampl = filtered_landscape_ampl, 
               filtered_landscape_del = filtered_landscape_del, 
               genome_mask = genome_mask_ampl, 
               type_mask = type_mask_ampl, 
               model_mask = c("ampl","del"))





















