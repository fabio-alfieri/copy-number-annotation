rm(list=ls())
gc(full=T)

packages <- c("stringr", "ggplot2", "tidyr", "tidyverse")

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

paths <- list.files("../Data/annotation/", pattern = "merged_res_annot")

annotations_list_names <- c("Arm-level_ampl",
                            "Arm-level_del",
                            "Mid-length_ampl",
                            "Mid-length_del",
                            "Chromosome-level_ampl",
                            "Chromosome-level_del",
                            "Small-scale_ampl",
                            "Small-scale_del",
                            "no_cluster_ampl",
                            "no_cluster_del")

annotations_list <- lapply(X = paths, FUN = function(x){
  
  table <- read.table(file = x, 
                      header = T,
                      sep = "\t")
  
})

names(annotations_list) <- annotations_list_names

annotations_counts <- lapply(X = annotations_list, FUN = function(x){table(x$annot_final)})
annotations_counts <- as.matrix(bind_rows(annotations_counts))
annotations_counts[is.na(annotations_counts)] <- 0
rownames(annotations_counts) <- annotations_list_names

write.table(x = annotations_counts, 
            file = "../Data/annotation/annotations_counts.tsv", 
            quote = F, sep = "\t", row.names = T, col.names = T)



