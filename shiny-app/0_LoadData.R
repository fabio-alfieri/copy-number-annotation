rm(list=ls())
gc(full=T)

centromere_table <- read.table("Data/cytoBand.txt", header = F)
colnames(centromere_table) <- c("chr", "start", "end", "cytoband", "type")

centromere_table <- centromere_table[centromere_table$type == "acen",]

lis_names <- c("ampl", "del")
model_types <- c("Mid-length", "Small-scale", "Chromosome-level", "Arm-level")

meta_list <- list()
for (model_type in model_types) {

  name_ampl <- paste0("Data/res_ratio_with_annot_with_backbone_ampl_sub_new_", model_type, ".tsv")
  name_del <- paste0("Data/res_ratio_with_annot_with_backbone_del_sub_new_", model_type, ".tsv")
  
  gc()
  res_ratio_with_annot_with_backbone_ampl <- read.table(name_ampl, header = T, sep = "\t")
  res_ratio_with_annot_with_backbone_del <- read.table(name_del, header = T, sep = "\t")
  gc()
  
  if(!(all(startsWith(x = as.character(res_ratio_with_annot_with_backbone_ampl$chr), prefix = "chr")))){
    res_ratio_with_annot_with_backbone_ampl$chr <- paste0("chr", res_ratio_with_annot_with_backbone_ampl$chr)
  }
  
  if(!(all(startsWith(x = as.character(res_ratio_with_annot_with_backbone_del$chr), prefix = "chr")))){
    res_ratio_with_annot_with_backbone_del$chr <- paste0("chr", res_ratio_with_annot_with_backbone_del$chr)
    
  }
  
  res_ratio_with_annot_with_backbone_ampl <- res_ratio_with_annot_with_backbone_ampl[order(as.integer(gsub(x = res_ratio_with_annot_with_backbone_ampl$chr, 
                                                               pattern = "chr", 
                                                               replacement = "")), 
                                               res_ratio_with_annot_with_backbone_ampl$start, 
                                               decreasing = F),]
  
  res_ratio_with_annot_with_backbone_del <- res_ratio_with_annot_with_backbone_del[order(as.integer(gsub(x = res_ratio_with_annot_with_backbone_del$chr, 
                                                               pattern = "chr", 
                                                               replacement = "")), 
                                               res_ratio_with_annot_with_backbone_del$start, 
                                               decreasing = F),]
  
  dfs_list <- list(ampl = res_ratio_with_annot_with_backbone_ampl,
                   del = res_ratio_with_annot_with_backbone_del)
  
  meta_list[[model_type]] <- dfs_list
  
}

load("Data/All_levels_backbonetables.RData")
backbone.100kbp <- chr_backbone_namesfixed[["0.1Mbp"]]
backbone.100kbp <- do.call(rbind, backbone.100kbp)
backbone.100kbp$binID <- paste0(backbone.100kbp$chr, "_", backbone.100kbp$bin)
backbone.100kbp$chr <- paste0("chr", backbone.100kbp$chr)

backbone.100kbp_granges <- GRanges(seqnames = backbone.100kbp$chr, 
                                   ranges = IRanges(start = backbone.100kbp$start_bin, 
                                                    end = backbone.100kbp$end_bin), 
                                   binID = backbone.100kbp$binID)

