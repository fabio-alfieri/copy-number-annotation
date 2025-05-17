parse_input_data <- function(shap.list, 
                             toplot.plot, 
                             clusters_explained, 
                             chr_backbone_namesfixed, 
                             centromere_table, 
                             clustering_depth){
  
  clustering_depth <- parse_clustering_depth(clustering_depth)
  
  selected_depth <- clustering_depth$selected
  not_selected_depth <- clustering_depth$not_selected
  
  # SHAP DF
  shap_clean_list <- list()
  
  # Chromatin features with similar function are aggregated
  aggregate_chromatin_states <- function(df) {
    df$promoters <- df$Length_Counts.E1 + df$Length_Counts.E2 + df$Length_Counts.E3 +
      df$Length_Counts.E4 + df$Length_Counts.E22 + df$Length_Counts.E23
    df$transcribed <- df$Length_Counts.E5 + df$Length_Counts.E6 + df$Length_Counts.E7 +
      df$Length_Counts.E8 + df$Length_Counts.E9 + df$Length_Counts.E10 +
      df$Length_Counts.E11 + df$Length_Counts.E12
    df$enhancers <- df$Length_Counts.E13 + df$Length_Counts.E14 + df$Length_Counts.E15 +
      df$Length_Counts.E16 + df$Length_Counts.E17 + df$Length_Counts.E18
    df$repressed <- df$Length_Counts.E19 + df$Length_Counts.E20 + df$Length_Counts.E24 +
      df$Length_Counts.E25
    df$accessible <- df$Length_Counts.E21
    return(df)
  }
  
  # Each "distance" feature is multiplied by -1
  change_distance <- function(df) {
    distances <- c("dist.to.closest.OG", "dist.to.closest.TSG", "dist.to.closest.FGS",
                   "distance.to.centromere", "distance.to.telomere", "Ess.distance_pancancer")
    df[,distances] <- -1*df[,distances]
    return(df)
  }
  
  shap_cols_to_discard <- c("Chromosome_Length", "Centromere_Length", "Centromere_Type", 'BIAS')
  shap_cols_to_keep <- c(#"labels", "type", "chr", "bin",
                         # "dist.to.closest.OG", "dist.to.closest.TSG", "dist.to.closest.FGS",
                         # "distance.to.centromere", "distance.to.telomere", "mutations_norm",
                         # "genes.bin", "Length_Counts.E17", "Length_Counts.E19", "Ess.distance_pancancer",
                         # "Length_Counts.E25", "Length_Counts.E1"
                         "dist.to.closest.OG", "dist.to.closest.TSG", "dist.to.closest.FGS",
                         "distance.to.centromere", "distance.to.telomere", "mutations_norm",
                         "genes.bin", "promoters", "transcribed", "Ess.distance_pancancer",
                         "enhancers", "repressed","accessible",
                         "labels")
  
  for (idx in seq_along(along.with = shap.list)) {
    
    name <- names(shap.list)[idx]
    shap.df <- shap.list[[idx]]
    
    shap.df <- shap.df[,apply(shap.df, 2, function(x){all(x!=0)})]
    
    shap.df <- aggregate_chromatin_states(shap.df)
    shap.df <- change_distance(shap.df)
    
    shap.df <- shap.df %>% select(-all_of(shap_cols_to_discard))
    shap.df <- shap.df %>% select(all_of(shap_cols_to_keep))
    
    shap.df$type <- do.call(rbind, str_split(shap.df$labels, pattern = '-'))[,2]
    shap.df$chr <- do.call(rbind, str_split(shap.df$labels, '_'))[,1]
    shap.df$bin <- do.call(rbind, str_split(do.call(rbind, str_split(shap.df$labels, '_'))[,2], '-'))[,1]
    
    shap.df$binID <- paste0(shap.df$chr, "_", shap.df$bin)
    shap.df$chr <- paste0("chr", shap.df$chr); shap.df$bin <- NULL
    
    shap_clean_list[[name]] <- shap.df
    
  }
  
  # LANDSCAPE DF
  toplot.plot$binID <- paste0(toplot.plot$chr, "_", toplot.plot$bin)
  toplot.plot$chr <- paste0("chr", toplot.plot$chr); toplot.plot$bin <- NULL; toplot.plot$cluster <- NULL
  colnames(toplot.plot)[str_detect(colnames(toplot.plot), pattern = 'ampl|del')] <- c("ampl", "del")

  toplot.plot$order <- seq_along(1:nrow(toplot.plot))
  
  selected_depth_code <- paste0("k",2**selected_depth)
  not_selected_depth_code <- paste0("k",2**not_selected_depth)
  
  top_selected_depth_code <- paste0("top_", selected_depth_code)

  colnames(clusters_explained[[selected_depth_code]]) <- c(selected_depth_code, top_selected_depth_code) 
  
  toplot.plot <- left_join(toplot.plot, clusters_explained[[selected_depth_code]], by = selected_depth_code)
  toplot.plot[,not_selected_depth_code] <- NULL
  
  colnames(toplot.plot)[str_detect(colnames(toplot.plot), pattern = 'Type')] <- 'type'
  
  toplot.plot <- toplot.plot[order(toplot.plot$order), ]; toplot.plot$order <- NULL

  # BACKBONE DF
  backbone.100kb <- chr_backbone_namesfixed$`0.1Mbp`; backbone.100kb <- dplyr::bind_rows(backbone.100kb)
  backbone.100kb$binID <- paste0(backbone.100kb$chr, "_", backbone.100kb$bin)
  backbone.100kb$chr <- paste0("chr", backbone.100kb$chr); backbone.100kb$bin <- NULL
  
  backbone.100kb <- GenomicRanges::GRanges(seqnames = backbone.100kb$chr, 
                                           ranges = IRanges(start = backbone.100kb$start_bin, end = backbone.100kb$end_bin),
                                           binID = backbone.100kb$binID)
  
  
  backbone.500kb <- chr_backbone_namesfixed$`0.5Mbp`; backbone.500kb <- dplyr::bind_rows(backbone.500kb)
  backbone.500kb$binID <- paste0(backbone.500kb$chr, "_", backbone.500kb$bin)
  backbone.500kb$chr <- paste0("chr", backbone.500kb$chr); backbone.500kb$bin <- NULL
  
  backbone.500kb <- GenomicRanges::GRanges(seqnames = backbone.500kb$chr, 
                                           ranges = IRanges(start = as.numeric(backbone.500kb$start_bin), 
                                                            end = as.numeric(backbone.500kb$end_bin)),
                                           binID = backbone.500kb$binID)
  
  
  # CENTROMERE DF
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
                  backbone.500kb = backbone.500kb,
                  centromere_table = centromere_table_out)
  
  return(outlist)
}
