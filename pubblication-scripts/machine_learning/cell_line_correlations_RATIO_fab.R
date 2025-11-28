# compare with Cell Line data
rm(list=ls())
gc(full= T)

# remotes::install_github("fabio-alfieri/karyotapR")
library(dplyr)
library(karyotapR)

model_class <- 'no_cluster' # Chromosome-level

# annotation table
merged_res_ratio.files <- grep('merged_res_ratio', list.files('/home/ieo7429/Scrivania/dev/Data/', full.names = T), value = T)
load(grep(model_class, merged_res_ratio.files, value = T))

# backbone table
load("/home/ieo7429/Scrivania/dev/Data/All_levels_backbonetables.RData")

####################### 

#######################  calculate arm specific scores for Selection and Occurrence annotations   ####################### 

merged_res_ratio$binID <- unlist(lapply(merged_res_ratio$labels,
                                        function(x){strsplit(x, "-")[[1]][1]}))

backbone.100kb <- do.call(rbind, chr_backbone_namesfixed$`0.1Mbp`)
backbone.100kb$binID <- paste0(backbone.100kb$chr, "_", backbone.100kb$bin)

backbone.100kb <- backbone.100kb %>% 
  dplyr::select(binID = binID, 
                start = start_bin, 
                end = end_bin,
                chr = chr)

# add p and q arm to backbone table
centromere <- read.table('/home/ieo7429/Scrivania/dev/Data/cell line data/hg19.centromere.bed', header = T)
# from https://github.com/brentp/poverlap/blob/master/data/hg19.centromere.bed
# I noticed that all centromeres are 3Mbp
# depends on the assembly. hg19 has fixed centromeres

centromere <- centromere %>% mutate(chr = gsub("chr", "", chrom)) %>% dplyr::filter(!chr %in% c('X','Y'))
centromere$chrom <- NULL
centromere <- apply(centromere, 2, as.numeric) %>% as.data.frame()

backbone_annotated <- backbone.100kb %>%
  left_join(centromere, by = c("chr" = "chr")) %>%
  mutate(
    arm = case_when(
      end.x < start.y ~ "p",              # before centromere start
      start.x > end.y ~ "q"              # after centromere end
    )
  ) %>% dplyr::select(binID, start = start.x, end = end.x, chr, arm)
backbone_annotated$arm <- paste0(backbone_annotated$chr,backbone_annotated$arm)

merged_res_ratio_with_backbone <- merge(x = merged_res_ratio, 
                                        y = backbone_annotated, by = "binID")

merged_res_ratio_with_backbone <- merged_res_ratio_with_backbone[order(as.integer(merged_res_ratio_with_backbone$chr),
                                                                       as.integer(merged_res_ratio_with_backbone$start),
                                                                       decreasing = FALSE),]

merged_res_ratio_with_backbone_only_centromeric_regions <- merged_res_ratio_with_backbone[grep('NA', merged_res_ratio_with_backbone$arm),]

merged_res_ratio_with_backbone <- merged_res_ratio_with_backbone[!merged_res_ratio_with_backbone$arm %in% 
                                                                   levels(factor(grep('NA', merged_res_ratio_with_backbone$arm, value = T))),]

merged_res_ratio_with_backbone <- merged_res_ratio_with_backbone[!merged_res_ratio_with_backbone$annot_3_classes %in% c('Relaxed Positive Selection', 'Relaxed Negative Selection'),]


#scores <- list()
#for(tt in levels(factor(merged_res_ratio_with_backbone$type))){
#  # calculate scores for arm and chromosome level
#  arm.scores <- merged_res_ratio_with_backbone %>% filter(type == tt) %>%
#    group_by(arm) %>%
#    summarise(median_sel_to_occ_ratio = median(selection_to_occurrence_ratio, na.rm = T),
#              median_pos_to_neg_ratio = median(positive_to_negative_ratio, na.rm = T),
#              mean_sel_to_occ_ratio = mean(selection_to_occurrence_ratio, na.rm = T),
#              mean_pos_to_neg_ratio = mean(positive_to_negative_ratio, na.rm = T), .groups = "drop_last") %>%
#    pivot_longer(cols = c(median_sel_to_occ_ratio, median_pos_to_neg_ratio, mean_sel_to_occ_ratio, mean_pos_to_neg_ratio))
#  
#  arm.scores <- arm.scores %>% mutate(name = factor(name))
#
#  classes <- levels(factor(arm.scores$name))
#  arms    <- sort(unique(arm.scores$value))
#  
#  # full <- merge(
#  #   expand.grid(arm = arms, name = classes, stringsAsFactors = FALSE),
#  #   arm.scores,
#  #   by = c("arm", "name"),
#  #   all.x = TRUE
#  # )
#  
#  # Replace missing values with 0 as requested
#  # for (col in c("n", "total", "freq")) {
#  #   full[[col]][is.na(full[[col]])] <- 0
#  # }
#  
#  split_tables <- split(arm.scores[order(arm.scores$arm), ], arm.scores$name)
#  
#  scores[[tt]] <- split_tables
#  
#}

scores <- list()
for(tt in levels(factor(merged_res_ratio_with_backbone$type))){
  arm.scores <- merged_res_ratio_with_backbone %>% filter(type == tt) %>%
    group_by(arm, annot_3_classes) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(
      total = sum(n),
      freq = n / total
    ) %>%
    arrange(arm, desc(freq))
  
  arm.scores <- arm.scores %>% mutate(annot_3_classes = factor(annot_3_classes))
  
  classes <- levels(factor(arm.scores$annot_3_classes))
  arms    <- sort(unique(arm.scores$arm))
  
  full <- merge(
    expand.grid(arm = arms, annot_3_classes = classes, stringsAsFactors = FALSE),
    arm.scores,
    by = c("arm", "annot_3_classes"),
    all.x = TRUE
  )
  
  for (col in c("n", "total", "freq")) {
    full[[col]][is.na(full[[col]])] <- 0
  }
  
  split_tables <- split(full[order(full$arm), ], full$annot_3_classes)
  
  scores[[tt]] <- split_tables
  
}

####################### 

# combine sub-tables (cancer types) for one class
combine_class_tables <- function(class_list) {
  cleaned <- lapply(names(class_list), function(ct) {
    df <- class_list[[ct]][, c("arm", "freq")]
    names(df)[2] <- ct
    df
  })
  combined <- Reduce(function(x, y) merge(x, y, by = "arm", all = TRUE), cleaned)
  return(combined)
}

combined_tables <- lapply(scores, combine_class_tables)
names(combined_tables) <- names(scores)

score.df <- do.call(rbind, combined_tables)
score.df$type <- do.call(rbind, strsplit(rownames(score.df), split = "\\."))[,1]

make_matrix <- function(df, selection) {
  df %>%
    dplyr::select(arm, type, !!sym(selection)) %>%
    pivot_wider(names_from = type, values_from = !!sym(selection), values_fill = 0) %>%
    column_to_rownames("arm")
}

selection_types <- setdiff(colnames(score.df), c("arm", "type"))

cor_mats.spearman <- lapply(selection_types, function(sel) {
  mat <- make_matrix(score.df, sel)
  cor(as.matrix(mat), use = "pairwise.complete.obs", method = "spearman")
})
names(cor_mats.spearman) <- selection_types

pheatmap::pheatmap(cor_mats.spearman[[1]], main = 'Positive Selection Cor Matrix (Spearman - Arm-level)')
pheatmap::pheatmap(cor_mats.spearman[[2]], main = 'Negative Selection Cor Matrix (Spearman - Arm-level)')


####################### 

#######################  get frequency of cell lines  ####################### 

scores.brca <- score.df %>% filter(type == 'BRCA')
tcga <- read.delim(file = '/home/ieo7429/Scrivania/dev/Data/cell line data/TCGA-Average.Arm.Amp.Frequency.txt')
tcga <- tcga[grep('Arm',tcga$X),]
tcga$arm <- str_replace_all(tcga$X, 'Arm','chr')

scores.brca$arm <- paste0('chr',scores.brca$arm)

tmp <- merge(scores.brca, tcga %>% dplyr::select(arm,BRCA))

cor.test(tmp$`Negative Selection`, tmp$BRCA)
cor.test(tmp$Occurrence, tmp$BRCA)
cor.test(tmp$`Positive Selection`, tmp$BRCA)

# cor.test(tmp$mean_pos_to_neg_ratio, tmp$BRCA)
# cor.test(tmp$mean_sel_to_occ_ratio, tmp$BRCA)
# cor.test(tmp$median_pos_to_neg_ratio, tmp$BRCA)
# cor.test(tmp$median_sel_to_occ_ratio, tmp$BRCA)

day28 <- read.delim(file = '/home/ieo7429/Scrivania/dev/Data/cell line data/day28-hMEC_gain_freq.tsv')
day2 <- read.delim(file = '/home/ieo7429/Scrivania/dev/Data/cell line data/day2-hMEC_gain_freq.tsv')
colnames(day28)[4] <- 'arm'
colnames(day2)[4] <- 'arm'

scores.day28 <- merge(tmp, day28)
scores.day2 <-  merge(tmp, day2)

cor.test(scores.day28$BRCA, scores.day28$hMEC_75nM)
cor.test(scores.day2$BRCA, scores.day2$hMEC_75nM)

cor.test(scores.day2$`Positive Selection`, scores.day2$hMEC_150nM)
cor.test(scores.day2$`Negative Selection`, scores.day2$hMEC_150nM)
cor.test(scores.day2$Occurrence, scores.day2$hMEC_150nM)

cor.test(scores.day28$`Positive Selection`, scores.day28$hMEC_150nM)
cor.test(scores.day28$`Negative Selection`, scores.day28$hMEC_150nM)

# PAAD
scores.paad <- score.df %>% filter(type == 'PAAD')
tcga <- read.delim(file = '/home/ieo7429/Scrivania/dev/Data/cell line data/TCGA-Average.Arm.Amp.Frequency.txt')
tcga <- tcga[grep('Arm',tcga$X),]
tcga$arm <- str_replace_all(tcga$X, 'Arm','chr')

scores.paad$arm <- paste0('chr',scores.paad$arm)

tmp <- merge(scores.paad, tcga %>% dplyr::select(arm,PAAD))

cor.test(tmp$`Positive Selection`, tmp$PAAD)
cor.test(tmp$`Negative Selection`, tmp$PAAD)
cor.test(tmp$Occurrence, tmp$PAAD)

day28 <- read.delim(file = '/home/ieo7429/Scrivania/dev/Data/cell line data/day28-hPDE_gain_freq.tsv')
day2 <- read.delim(file = '/home/ieo7429/Scrivania/dev/Data/cell line data/day2-hPDE_gain_freq.tsv')
colnames(day28)[4] <- 'arm'
colnames(day2)[4] <- 'arm'

scores.day28 <- merge(scores.paad, day28)
scores.day2 <-  merge(scores.paad, day2)

cor.test(scores.day2$`Positive Selection`, scores.day2$hPDE_75nM)
cor.test(scores.day2$`Negative Selection`, scores.day2$hPDE_75nM)
cor.test(scores.day2$Occurrence, scores.day2$hPDE_75nM)

cor.test(scores.day28$`Positive Selection`, scores.day28$hPDE_150nM)
cor.test(scores.day28$`Positive Selection`, scores.day28$hPDE_75nM)
