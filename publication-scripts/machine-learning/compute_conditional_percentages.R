rm(list=ls())
gc(full=T)


wd <- 'path/to/GitHub/copy-number-annotation/'

packages <- c("stringr", "ggplot2", "tidyr", "tidyverse")

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

annot_statistics <- lapply(X = annots[c(3,8)], 
                           FUN = function(annot_path){
                             
                             annot <- read.table(annot_path, header = T, sep = "\t")
                             annot <- annot[annot$type == "BRCA",]
                             
                             parts <- strsplit(x = annot_path, split = "[_.]")[[1]]
                             annot_name <- paste(parts[parts %in% patterns], collapse = "_")
                             print(annot_name)
                             
                             perc <- (nrow(annot[annot$region_is_diploid == T,]) / nrow(annot)) * 100
                             
                             return(perc)
                             
                           })


summary_statistics <- Reduce(function(x, y) merge(x, y, by = c("region_is_diploid","annot_final"), all = TRUE), annot_statistics)
summary_statistics[is.na(summary_statistics)] <- 0

summary_statistics$avg <- (summary_statistics$`ampl_Mid-length_perc` + summary_statistics$`del_Mid-length_perc`) / 2

write.table(x = summary_statistics, file = paste0(wd, "data/annotation/summary_statistics_conditional.tsv"), quote = F, col.names = T)

