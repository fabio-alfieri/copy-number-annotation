rm(list=ls())
gc(full=T)

wd <- 'path/to/GitHub/copy-number-annotation/'

packages <- c("stringr", "reshape2", "tidyr", "tidyverse")

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

annots <- list.files(paste0(wd, "data/annotation/"), pattern = "res_ratio_with")

patterns <- c("Arm-level",
              "Chromosome-level",
              "Mid-length",
              "Small-scale",
              "no", "cluster", 
              "ampl", "del"
              )

annot_statistics <- lapply(X = annots, 
                           FUN = function(annot_path){
                             
                             annot <- read.table(annot_path, header = T, sep = "\t")
                             
                             parts <- strsplit(x = annot_path, split = "[_.]")[[1]]
                             annot_name <- paste(parts[parts %in% patterns], collapse = "_")
                             print(annot_name)
                             
                             tot_type_bins <- annot %>% group_by(type) %>% summarise(count_bins = n())
                             num_per_type <- annot %>% group_by(type, annot_final) %>% summarise(count_annots = n())
                             annot_per_type <- merge(x = tot_type_bins, y = num_per_type, by = "type")
                               
                             annot_per_type$type <- as.factor(annot_per_type$type)
                             annot_per_type <- annot_per_type[order(annot_per_type$type, decreasing = F),]
                             annot_per_type <- annot_per_type %>% 
                               group_by(type, annot_final) %>% 
                               summarise(
                                 count_bins = unique(count_bins),
                                 perc_annot = (count_annots / count_bins) * 100
                               )
                             
                             annot_per_type$count_bins <- NULL
                             colnames(annot_per_type) <- c("type","annot_final",paste0(annot_name, "_perc"))
                             
                             return(annot_per_type)
                             
                           })


summary_statistics <- Reduce(function(x, y) merge(x, y, by = c("type","annot_final"), all = TRUE), annot_statistics)
summary_statistics[is.na(summary_statistics)] <- 0

write.table(x = summary_statistics, file = paste0(wd, "data/annotation/summary_statistics.tsv"), quote = F, col.names = T)

rownames(summary_statistics) <- paste0(summary_statistics$type, "_", summary_statistics$annot_final)
summary_statistics$type <- NULL
summary_statistics$annot_final <- NULL


pdf(paste0(wd, "data/plots/heatmap.pdf"), width = 20, height = 20) 
pheatmap(summary_statistics, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c("white", "red"))(50))
dev.off()
