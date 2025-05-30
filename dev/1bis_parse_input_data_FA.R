
process_out_annot_list <- function(out_annot_list){
  
  model_types <- c("ampl", "del")
  
  outlist <- list()
  
  for (idx in seq_along(along.with = model_types)){
    
    model <- model_types[idx]
    
    source <- out_annot_list[[model]]
    source$p.final <- NULL
    annot <- source$aggregated
    clustering_depths <- names(annot)
    final_toplot <- source$toplot
    
    for (clustering_depth in clustering_depths){
      
      corresponding_col <- match(clustering_depth, colnames(final_toplot))
      curr_annot <- annot[[clustering_depth]]
      
      new_colname <- paste0("top_", clustering_depth)
      
      final_toplot <- merge(x = final_toplot, 
                            y = curr_annot, 
                            by.x = clustering_depth,
                            by.y = "cluster")
      
      final_toplot[[new_colname]] <- final_toplot$top_features; final_toplot$top_features <- NULL
      
    }
    outlist[[model]] <- final_toplot
  }
  
  return(outlist)
}
process_toplot.plot <- function(toplot.plot, clustering_depth){
  
  selected_depth <- clustering_depth$selected
  not_selected_depth <- clustering_depth$not_selected
  
  toplot.plot <- toplot.plot[order(as.integer(toplot.plot$chr),
                                   as.integer(toplot.plot$bin)),]
  
  toplot.plot$binID <- paste0(toplot.plot$chr, "_", toplot.plot$bin)
  toplot.plot$chr <- paste0("chr", toplot.plot$chr); toplot.plot$bin <- NULL; toplot.plot$cluster <- NULL
  colnames(toplot.plot)[str_detect(colnames(toplot.plot), pattern = 'ampl|del')] <- c("ampl", "del")

  toplot.plot$order <- seq_along(1:nrow(toplot.plot))
  
  selected_depth_code <- paste0("k",2**selected_depth)
  not_selected_depth_code <- paste0("k",2**not_selected_depth)
  
  top_selected_depth_code <- paste0("top_", selected_depth_code)
  top_not_selected_depth_code <- paste0("top_", not_selected_depth_code)
  
  # colnames(clusters_explained[[selected_depth_code]]) <- c(selected_depth_code, top_selected_depth_code) 
  # toplot.plot <- left_join(toplot.plot, clusters_explained[[selected_depth_code]], by = selected_depth_code)
  
  toplot.plot[,not_selected_depth_code] <- NULL
  toplot.plot[,top_not_selected_depth_code] <- NULL
  
  colnames(toplot.plot)[str_detect(colnames(toplot.plot), pattern = 'Type')] <- 'type'
  
  toplot.plot <- toplot.plot[order(toplot.plot$order), ]; toplot.plot$order <- NULL
  
  return(toplot.plot)
}
process_shap_list <- function(shap.list){
  
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
  return(shap_clean_list)
}
process_backbone <- function(chr_backbone_namesfixed){
  
  backbone.100kb <- chr_backbone_namesfixed$`0.1Mbp`; backbone.100kb <- dplyr::bind_rows(backbone.100kb)
  backbone.100kb$binID <- paste0(backbone.100kb$chr, "_", backbone.100kb$bin)
  backbone.100kb$chr <- paste0("chr", backbone.100kb$chr); backbone.100kb$bin <- NULL
  
  backbone.100kb <- GenomicRanges::GRanges(seqnames = backbone.100kb$chr, 
                                           ranges = IRanges(start = backbone.100kb$start_bin, end = backbone.100kb$end_bin),
                                           binID = backbone.100kb$binID)
  
  
  return(backbone.100kb)
  
}
process_centromere_table <- function(centromere_table, backbone.100kb){
  
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
  
  return(centromere_table_out)
  
}



parse_input_data <- function(shap.list, 
                             out_annot_list, 
                             chr_backbone_namesfixed, 
                             centromere_table, 
                             clustering_depth){
  
  clustering_depth <- parse_clustering_depth(clustering_depth)
  
  shap_clean_list <- process_shap_list(shap.list = shap.list)
  
  backbone.100kb <- process_backbone(chr_backbone_namesfixed = chr_backbone_namesfixed)
  
  centromere_table_out <- process_centromere_table(centromere_table = centromere_table, 
                                                   backbone.100kb = backbone.100kb)
  
  out_annot_list_processed <- lapply(X = process_out_annot_list(out_annot_list), 
                                     FUN = function(x){
                                       process_toplot.plot(x, clustering_depth = clustering_depth)
                                     })
  
  outlist <- list(shap.list = shap_clean_list,
                  out_annot_list_processed = out_annot_list_processed,
                  backbone.100kb = backbone.100kb,
                  centromere_table = centromere_table_out)
  
  return(outlist)
}
