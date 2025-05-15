rm(list=ls())
gc(full=T)

shap.list <- readRDS("dev/Data/shap_Mid-length_AmplDel.rds")
output.annotation <- readRDS("dev/Data/output_annotation.rds") # NEW
toplot.plot <- output.annotation$ampl$toplot # updated
clusters_explained <- output.annotation$ampl$aggregated # updated
centromere_table <- read.table("dev/Data/centomere.tsv", header = T)
load("dev/Data/All_levels_backbonetables.RData")

parse_input_data <- function(shap.list, toplot.plot, clusters_explained, chr_backbone_namesfixed, centromere_table){
  
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
  toplot.plot$chr <- paste0("chr", toplot.plot$chr); toplot.plot$bin <- NULL
  colnames(toplot.plot)[str_detect(colnames(toplot.plot), pattern = 'ampl|del')] <- c("ampl", "del")
  # toplot.plot$clusters10 <- NULL; toplot.plot$clusters15 <- NULL
  
  toplot.plot$order <- seq_along(1:nrow(toplot.plot))
  
  # Qui bisogna capire come fare per 
  
  colnames(clusters_explained$k2) <- c('k2','top_k2')
  colnames(clusters_explained$k4) <- c('k4','top_k4')
  colnames(clusters_explained$k8) <- c('k8','top_k8')
  colnames(clusters_explained$k16)<- c('k16','top_k16')
  
  toplot.plot <- left_join(
    left_join(
      left_join(left_join(toplot.plot, 
                          clusters_explained$k2, by = 'k2'),
                clusters_explained$k4, by = 'k4'), 
      clusters_explained$k8, by = 'k8'),
    clusters_explained$k16, by = 'k16')
  
  # toplot.plot <- merge(x = toplot.plot, 
  #                      y = clusters_explained, 
  #                      by = "clusters20", 
  #                      all.x = TRUE)
  
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
  
  # BACKBONE DF
  backbone.100kb <- chr_backbone_namesfixed$`0.1Mbp`; backbone.100kb <- dplyr::bind_rows(backbone.100kb)
  backbone.100kb$binID <- paste0(backbone.100kb$chr, "_", backbone.100kb$bin)
  backbone.100kb$chr <- paste0("chr", backbone.100kb$chr); backbone.100kb$bin <- NULL
  
  backbone.100kb <- GenomicRanges::GRanges(seqnames = backbone.100kb$chr, 
                                           ranges = IRanges(start = backbone.100kb$start_bin, end = backbone.100kb$end_bin),
                                           binID = backbone.100kb$binID)
  
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
                  centromere_table = centromere_table_out)
  
  return(outlist)
}
