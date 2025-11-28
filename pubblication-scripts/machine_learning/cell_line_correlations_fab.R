# compare with Cell Line data
rm(list=ls())
gc(full= T)

library(dplyr)
# remotes::install_github("fabio-alfieri/karyotapR")
library(karyotapR)
library(tidyverse)

model_class <- 'Arm-level' # Chromosome-level

# annotation table
merged_res_annot.files <- grep('merged_res_annot', list.files('/home/ieo7429/Scrivania/dev/Data/', 
                                                     full.names = T), value = T)

# load(grep('del', grep(model_class, merged_res_annot.files, value = T), value = T))
# merged_res_annot_del <- merged_res_annot

# load(grep('ampl', grep(model_class, merged_res_annot.files, value = T), value = T))
merged_res_annot <- read_tsv(grep('ampl_0.6', grep(model_class, merged_res_annot.files, value = T), value = T))
# backbone table
load("/home/ieo7429/Scrivania/dev/Data/All_levels_backbonetables.RData")


####################### 

#######################  calculate arm specific scores for Selection and Occurrence annotations   ####################### 

merged_res_annot$binID <- unlist(lapply(merged_res_annot$labels,
                                        function(x){strsplit(x, "-")[[1]][1]}))
# merged_res_annot_del$binID <- unlist(lapply(merged_res_annot_del$labels,
#                                         function(x){strsplit(x, "-")[[1]][1]}))

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
centromere <- centromere %>% mutate(chr = gsub("chr", "", chrom)) %>% 
  dplyr::filter(!chr %in% c('X','Y'))
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

merged_res_annot_with_backbone <- merge(x = merged_res_annot, 
                                        y = backbone_annotated, by = "binID")

merged_res_annot_with_backbone <- merged_res_annot_with_backbone[order(as.integer(merged_res_annot_with_backbone$chr),
                                                                       as.integer(merged_res_annot_with_backbone$start),
                                                                       decreasing = FALSE),]

merged_res_annot_with_backbone_only_centromeric_regions <- merged_res_annot_with_backbone[grep('NA', merged_res_annot_with_backbone$arm),]

merged_res_annot_with_backbone <- merged_res_annot_with_backbone[!merged_res_annot_with_backbone$arm %in% 
                                                                   levels(factor(grep('NA', merged_res_annot_with_backbone$arm, value = T))),]

# merged_res_annot_with_backbone <- merged_res_annot_with_backbone[!merged_res_annot_with_backbone$annot_3_classes %in% c('Relaxed Positive Selection', 'Relaxed Negative Selection'),]


scores <- list()
for(tt in levels(factor(merged_res_annot_with_backbone$type))){
  # calculate scores for arm and chromosome level
  arm.scores <- merged_res_annot_with_backbone %>% filter(type == tt) %>%
    # group_by(arm) %>% summarise(mean = mean(prediction, na.rm = T)) %>%
    group_by(arm, annot_3_classes) %>%
    summarise(n = n(), max = max(prediction, na.rm = T), .groups = "drop_last") %>%
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
  
  # Replace missing values with 0 as requested
  for (col in c("n", "total", "freq")) {
    full[[col]][is.na(full[[col]])] <- 0
  }
  
  split_tables <- split(full[order(full$arm), ], full$annot_3_classes)
  
  scores[[tt]] <- split_tables
  
}

# scores <- list()
# for(tt in levels(factor(merged_res_annot_with_backbone$type))){
#   # calculate scores for arm and chromosome level
#   arm.scores <- merged_res_annot_with_backbone %>% filter(type == tt) %>%
#     group_by(arm, annot_3_classes) %>%
#     summarise(n = n(), .groups = "drop_last") %>%
#     mutate(
#       total = sum(n),
#       freq = n / total
#     ) %>%
#     arrange(arm, desc(freq))
# 
#   arm.scores <- arm.scores %>% mutate(annot_3_classes = factor(annot_3_classes))
# 
#   classes <- levels(factor(arm.scores$annot_3_classes))
#   arms    <- sort(unique(arm.scores$arm))
# 
#   full <- merge(
#     expand.grid(arm = arms, annot_3_classes = classes, stringsAsFactors = FALSE),
#     arm.scores,
#     by = c("arm", "annot_3_classes"),
#     all.x = TRUE
#   )
# 
#   # Replace missing values with 0 as requested
#   for (col in c("n", "total", "freq")) {
#     full[[col]][is.na(full[[col]])] <- 0
#   }
# 
#   split_tables <- split(full[order(full$arm), ], full$annot_3_classes)
# 
#   scores[[tt]] <- split_tables
# 
# }

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

score.df <- bind_rows(combined_tables)
# score.df <- do.call(rbind, combined_tables)
# score.df$type <- do.call(rbind, strsplit(rownames(score.df), split = "\\."))[,1]

score.df$type <- rep(names(combined_tables),each = 34)

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

pheatmap::pheatmap(cor_mats.spearman[["Positive Selection"]], main = 'Positive Selection Cor Matrix (Spearman - Arm-level)')
pheatmap::pheatmap(cor_mats.spearman[["Negative Selection"]], main = 'Negative Selection Cor Matrix (Spearman - Arm-level)')


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

####################### 

#######################  Correlate ML scores with Cell Line RATIOS  ####################### 

load(file = '/home/ieo7429/Scrivania/dev/Data/cell line data/cellLine_ratios/arm_ratio_gain_day28_hMEC.RData')
day28 <- freq_ratio_gain_tidy; rm(freq_ratio_gain_tidy)
load(file = '/home/ieo7429/Scrivania/dev/Data/cell line data/cellLine_ratios/arm_ratio_gain_day2_hMEC.RData')
day2 <- freq_ratio_gain_tidy; rm(freq_ratio_gain_tidy)

cellLines.ratio <- rbind(day28,day2) #%>% filter(!arm %in% c('chr20p','chr20q'))
cellLines.ratio <- merge(tmp, cellLines.ratio)
cellLines.ratio$pos_neg_selection <- cellLines.ratio$`Positive Selection`-cellLines.ratio$`Negative Selection`

save(cellLines.ratio, file = 'cellLine_ratios/cellLine.ratio_hMEC_gain.RData')

ggplot(cellLines.ratio, aes(x = ratio, y = pos_neg_selection)) +
  geom_point() +
  facet_wrap(~day+comparison) +
  geom_smooth(method = 'lm') +
  geom_text(aes(label = arm))
ggplot(cellLines.ratio, aes(x = ratio, y = `Positive Selection`)) +
  geom_point() +
  facet_wrap(~day+comparison, scale = 'free') +
  geom_smooth(method = 'lm') +
  geom_text(aes(label = arm))
ggplot(cellLines.ratio, aes(x = ratio, y = `Negative Selection`)) +
  geom_point() +
  facet_wrap(~day+comparison, scale = 'free') +
  geom_smooth(method = 'lm') +
  geom_text(aes(label = arm))
ggplot(cellLines.ratio, aes(x = ratio, y = Occurrence+`Negative Selection`+`Relaxed Negative Selection`)) +
  geom_point() +
  facet_wrap(~day+comparison) +
  geom_smooth(method = 'lm') +
  geom_text(aes(label = arm))
ggplot(cellLines.ratio, aes(x = BRCA, y = ratio)) +
  geom_point() +
  facet_wrap(~day+comparison) +
  geom_smooth(method = 'lm') +
  geom_text(aes(label = arm))

x <- cellLines.ratio %>% filter(day == 'day28', comparison == '75nM_vs_ctr')
cor.test(x = x$ratio, y = x$`Positive Selection`)
cor.test(x = x$ratio, y = x$pos_neg_selection)
cor.test(x = x$ratio, y = x$BRCA)

library(betareg)
# Normalizzare ratio (se >1)
x$ratio_scaled <- (x$ratio - min(x$ratio) + 0.001) / (max(x$ratio) - min(x$ratio) + 0.002)

model <- betareg(ratio_scaled ~ `Positive Selection`, data = x)
summary(model)
model <- betareg(ratio_scaled ~ `Negative Selection`+Occurrence, data = x)
summary(model)
model <- betareg(ratio_scaled ~ Occurrence, data = x)
summary(model)


#######################  get frequency of cell lines  ####################### 

scores.brca <- score.df %>% filter(type == 'PAAD')
tcga <- read.delim(file = '/home/ieo7429/Scrivania/dev/Data/cell line data/TCGA-Average.Arm.Amp.Frequency.txt')
tcga <- tcga[grep('Arm',tcga$X),]
tcga$arm <- stringr::str_replace_all(tcga$X, 'Arm','chr')

scores.brca$arm <- paste0('chr',scores.brca$arm)

tmp <- merge(scores.brca, tcga %>% dplyr::select(arm,PAAD))

cor.test(tmp$`Negative Selection`, tmp$PAAD)
cor.test(tmp$Occurrence, tmp$PAAD)
cor.test(tmp$`Positive Selection`, tmp$PAAD)

load(file = '/home/ieo7429/Scrivania/dev/Data/cell line data/cellLine_ratios/arm_ratio_gain_day28_hPDE.RData')
day28 <- freq_ratio_gain_tidy; rm(freq_ratio_gain_tidy)
load(file = '/home/ieo7429/Scrivania/dev/Data/cell line data/cellLine_ratios/arm_ratio_gain_day2_hPDE.RData')
day2 <- freq_ratio_gain_tidy; rm(freq_ratio_gain_tidy)

cellLines.ratio <- rbind(day28,day2) #%>% filter(!arm %in% c('chr20p','chr20q'))
cellLines.ratio <- merge(tmp, cellLines.ratio)
cellLines.ratio$pos_neg_selection <- cellLines.ratio$`Positive Selection`-cellLines.ratio$`Negative Selection`
# cellLines.ratio <- cellLines.ratio[!cellLines.ratio$ratio == 'Inf',]

save(cellLines.ratio, file = 'cellLine_ratios/cellLine.ratio_hPDE_gain.RData')


cellLines.ratio %>% filter(day == 'day28', comparison == '75nM_vs_ctr') %>%
  ggplot(aes(x = ratio, y = PAAD)) +
  geom_point() +
  geom_text(aes(label = arm)) +
  geom_smooth(method = 'lm')

ggplot(cellLines.ratio, aes(x = ratio, y = pos_neg_selection)) +
  geom_point() +
  facet_wrap(~day+comparison) +
  geom_smooth(method = 'lm') +
  geom_text(aes(label = arm))
ggplot(cellLines.ratio, aes(x = ratio, y = `Positive Selection`)) +
  geom_point() +
  facet_wrap(~day+comparison, scale = 'free') +
  geom_smooth(method = 'lm') +
  geom_text(aes(label = arm))

x <- cellLines.ratio %>% filter(day == 'day28', comparison == '75nM_vs_ctr')
cor.test(x = x$ratio, y = x$`Positive Selection`)
cor.test(x = x$ratio, y = x$`Negative Selection`)
cor.test(x = x$ratio, y = x$pos_neg_selection)


x <- cellLines.ratio %>% filter(day == 'day2', comparison == '75nM_vs_ctr')
cor.test(x = x$ratio, y = x$`Negative Selection`)
cor.test(x = x$ratio, y = x$pos_neg_selection)

library(betareg)
# Normalizzare ratio (se >1)
x$ratio_scaled <- (x$ratio - min(x$ratio) + 0.001) / (max(x$ratio) - min(x$ratio) + 0.002)
x$positive.selection <- x$`Positive Selection`+x$`Relaxed Positive Selection`
x$negative.selection <- x$`Negative Selection`+x$`Relaxed Negative Selection`

model <- betareg(ratio_scaled ~ `Positive Selection`, data = x)
summary(model)
model <- betareg(ratio_scaled ~ positive.selection, data = x)
summary(model)

model <- betareg(ratio_scaled ~ `Negative Selection`, data = x)
summary(model)
model <- betareg(ratio_scaled ~ negative.selection, data = x)
summary(model)

model <- betareg(ratio_scaled ~ `Negative Selection`+Occurrence, data = x)
summary(model)
model <- betareg(ratio_scaled ~ Occurrence, data = x)
summary(model)

ggplot(x, aes(x = positive.selection, y = ratio_scaled)) +
  # geom_jitter(width = 0.01, height = 0, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point()

ggplot(x, aes(x = as.factor(positive.selection == 0), y = ratio_scaled)) +
  geom_boxplot() +
  geom_point()

wilcox.test(x$ratio_scaled ~ as.factor(x$positive.selection == 0), exact = F)


######################## Correlate ML scores with Cell Line FREQUENCIES  ####################### 
if(F){
  day28 <- read.delim(file = '/home/ieo7429/Scrivania/dev/Data/cell line data/day28-hMEC_gain_freq_filtTotalReads.tsv')
  day2 <- read.delim(file = '/home/ieo7429/Scrivania/dev/Data/cell line data/day2-hMEC_gain_freq_filtTotalReads.tsv')
  colnames(day28)[4] <- 'arm'
  colnames(day2)[4] <- 'arm'
  
  
  scores.day28 <- merge(scores.brca, day28)
  scores.day2 <-  merge(scores.brca, day2)
  
  scores.day28 <- scores.day28 %>% filter(!arm %in% c('chr20p','chr20q'))
  scores.day2 <- scores.day2 %>% filter(!arm %in% c('chr20p','chr20q'))
  
  cor.test(scores.day28$`Positive Selection`, scores.day28$hMEC_150nM)
  cor.test(scores.day28$`Positive Selection`, scores.day28$hMEC_75nM)
  cor.test(scores.day28$Occurrence, scores.day28$hMEC_150nM)
  cor.test(scores.day28$`Negative Selection`, scores.day28$hMEC_150nM)
  
  cor.test(scores.day2$`Positive Selection`, scores.day2$hMEC_150nM)
  cor.test(scores.day2$`Negative Selection`, scores.day2$hMEC_150nM)
  cor.test(scores.day2$Occurrence, scores.day2$hMEC_150nM)
  
  
  scores.day28$pos_neg <- scores.day28$`Positive Selection`-scores.day28$`Negative Selection`
  cor.test(scores.day28$pos_neg, scores.day28$hMEC_150nM)
  cor.test(scores.day28$pos_neg, scores.day28$hMEC_75nM)
  
  scores.day28 %>% #filter(!arm %in% c('chr20')) %>% 
    ggplot(aes(x = pos_neg, y = hMEC_150nM, labels = arm)) +
    geom_point() +
    geom_smooth(method = 'lm')
  
  
  
  
  
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
  
  cor.test(scores.day28$`Positive Selection`, scores.day28$hPDE_150nM)
  cor.test(scores.day28$`Positive Selection`, scores.day28$hPDE_75nM)
  
  cor.test(scores.day2$`Positive Selection`, scores.day2$hPDE_75nM)
  cor.test(scores.day2$`Negative Selection`, scores.day2$hPDE_75nM)
  cor.test(scores.day2$Occurrence, scores.day2$hPDE_75nM)
  
}
