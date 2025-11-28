
process_centromere_table <- function(centromere_table, backbone.100kb){
  
  centromere_table$cytoband <- NULL
  centromere_table$type <- NULL
  centromere_table <- centromere_table %>% filter(!(chr %in% c("chrX", "chrY")))
  centromere_table <- centromere_table %>% group_by(chr) %>% summarise(start= min(start), end = max(end))
  centromere_table$midpoint <- (centromere_table$end + centromere_table$start) / 2
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

