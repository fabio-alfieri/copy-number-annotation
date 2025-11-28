rm(list = ls())
gc()

load("~/Scrivania/dev/Data/All_levels_backbonetables.RData")

cnames <- c("chr", "binID", "start", "end")

backbone_arm <- do.call(rbind, chr_backbone_namesfixed$Arm)
colnames(backbone_arm) <- cnames
backbone_arm_gr <- GRanges(backbone_arm)

backbone_100kbp <- do.call(rbind, chr_backbone_namesfixed$`0.1Mbp`)
colnames(backbone_100kbp) <- cnames
backbone_100kbp$binID <- paste0(backbone_100kbp$chr, "_", backbone_100kbp$binID)
backbone_100kbp_gr <- GRanges(backbone_100kbp)

hits <- findOverlaps(query = backbone_100kbp_gr, backbone_arm_gr)
query <- queryHits(hits)
subject <- subjectHits(hits)

match_100kbp_arm <- data.frame(
                      binID_100kbp = mcols(backbone_100kbp_gr[query])$binID,
                      binID_arm = mcols(backbone_arm_gr[subject])$binID,
                      idx = 1:length(mcols(backbone_arm_gr[subject])$binID)
                    )

overlaps <- data.frame(pintersect(backbone_100kbp_gr[query], backbone_arm_gr[subject]))
overlaps$strand <- NULL; overlaps$hit <- NULL
overlaps$idx <- 1:nrow(overlaps)

colnames(overlaps) <- c("chr", "start", "end", "width", "binID", "idx")

duplicated_idxs <- which(duplicated(data.frame(pintersect(backbone_100kbp_gr[query], backbone_arm_gr[subject]))$binID))
duplicated_bins <- mcols(pintersect(backbone_100kbp_gr[query], backbone_arm_gr[subject]))[duplicated_idxs,"binID"]

to.remove <- overlaps[overlaps$binID %in% duplicated_bins,c("binID", "idx", "width")] %>%
  group_by(binID) %>%
  slice_min(width, n = 1) %>%
  pull(idx)

match_100kbp_arm_clean <- match_100kbp_arm[!(match_100kbp_arm$idx %in% to.remove),]
match_100kbp_arm_clean$idx <- NULL

write.table(x = match_100kbp_arm_clean, file = "dev/Data/match_bins_arms_clean.tsv", quote = F, row.names = F, col.names = T)





