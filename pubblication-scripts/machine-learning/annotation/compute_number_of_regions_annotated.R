paths <- c("dev/Data/merged_res_annot_Arm-level_ampl_0.6_0.0205180818215013_0.002.tsv",
           "dev/Data/merged_res_annot_Arm-level_del_0.6_0.0201767960563302_0.01.tsv",
           "dev/Data/merged_res_annot_Mid-length_ampl_0.6_0.0632314186543226_0.02.tsv",
           "dev/Data/merged_res_annot_Mid-length_del_0.6_0.0869447428733111_0.025.tsv",
           "dev/Data/merged_res_annot_Chromosome-level_ampl_0.6_0.00890466617420316_0.01.tsv",
           "dev/Data/merged_res_annot_Chromosome-level_del_0.6_0.00725054279901087_0.tsv",
           "dev/Data/merged_res_annot_Small-scale_ampl_0.6_0.00773112666793168_4e-04.tsv",
           "dev/Data/merged_res_annot_Small-scale_del_0.6_0.00519783566705883_0.tsv",
           "dev/Data/merged_res_annot_no_cluster_ampl_0.6_0.0974942784756422_0.04.tsv",
           "dev/Data/merged_res_annot_no_cluster_del_0.6_0.121279673650861_0.04.tsv")

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
            file = "dev/Data/annotations_counts.tsv", 
            quote = F, sep = "\t", row.names = T, col.names = T)



