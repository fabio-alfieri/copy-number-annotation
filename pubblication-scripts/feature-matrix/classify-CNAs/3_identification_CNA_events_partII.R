
# Based on fraction (total segment length/ bin length) - Decide CNA event
# If both amlification and deletion fraction are higher than coverage threshold, assign the event with a higher fraction

# Find the total CNA score
scores_final$Total.CNA <- scores_final$Ampl.fraction + scores_final$Del.fraction
 
# Based on coverage threshold, assign the event
scores_final$Ampl <- ifelse(scores_final$Ampl.fraction > covThr,1,0)
scores_final$Del <- ifelse(scores_final$Del.fraction > covThr,-1,0)
scores_final$CNA <- ifelse(scores_final$Total.CNA > covThr,1,0)

# Assign event based on amplification and deletion frequencies
scores_final$Score <- ifelse(scores_final$Ampl == 1 & scores_final$Del == 0, 1,
                             ifelse(scores_final$Ampl == 0 & scores_final$Del == -1, -1,
                                    ifelse(scores_final$Ampl == 0 & scores_final$Del == 0, 0,
                                           ifelse(scores_final$Ampl.fraction > scores_final$Del.fraction,1,-1))))

# Assign event based on the Total CNA
scores_final$Total.Score <- ifelse(scores_final$CNA == 1, 1, 0)

# Final Patient-bin matrix (Amplification and deletion scores)
final.df <- scores_final[,colnames(scores_final) %in% c("Sample","bin","Score")]
final.df <- final.df %>% dcast(Sample ~ bin, value.var = "Score")
final.df <- as.data.frame(final.df)

bin_order <- paste(coordinates$chr,coordinates$bin, sep = "_")
bin_order <- intersect(bin_order, colnames(final.df))
final.df <- final.df[,c("Sample",bin_order)]

# Final Patient-bin matrix (Total CNA scores)
final.df.total <- scores_final[,colnames(scores_final) %in% c("Sample","bin","Total.Score")]
final.df.total <- final.df.total %>% dcast(Sample ~ bin, value.var = "Total.Score")
final.df.total <- as.data.frame(final.df.total)

bin_order <- paste(coordinates$chr,coordinates$bin, sep = "_")
bin_order <- intersect(bin_order, colnames(final.df.total))
final.df.total <- final.df.total[,c("Sample",bin_order)]

Results.data[[paste0("ArmFraction_",armf)]][[paste0("SegMean_",sg.mean)]][[level]][[sg.class]][[tumor_type]] <- final.df

Results.data.totalCNA[[paste0("ArmFraction_",armf)]][[paste0("SegMean_",sg.mean)]][[level]][[sg.class]][[tumor_type]] <- final.df.total
