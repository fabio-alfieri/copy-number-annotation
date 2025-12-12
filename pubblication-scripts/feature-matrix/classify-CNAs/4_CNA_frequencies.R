
# if the genome coverage is same for all the samples, we can calculate frequencies by dividing amp/del to the total number of samples
# which the length(unique(samples)) of segment file.

# Calculate amplification and deletion frequencies for each bin by using patient matrix

# Number of patients for each bin - scna_w_armfraction for the specific tumor_type
segmentfile <- paste0("./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Segments_mapped_to_bins/",
                      cohort,"_",level,".RData")
load(segmentfile)
countsamples <- function(bin){
  samplesize <- length(unique(bin$Sample))
  return(samplesize)
}
samplesize <- mclapply(res, countsamples, mc.cores = 15)
samplesize <- do.call(rbind,samplesize)
samplesize <- as.data.frame(samplesize)
samplesize$bin <- rownames(samplesize)

rm("res")
gc()

# Prepare the matrix
rownames(final.df) <- final.df$Sample
final.df <- final.df[,-1]

# Learn the number of bin
num.bin <- ncol(final.df)

if(event == "TotalCNA"){
  # Counting the events
  events <- final.df %>%  summarise_all(~ sum(.x == 1, na.rm = TRUE))
  events <- reshape2::melt(events)
  colnames(events)[2] <- "Total.CNA"
  
  all <- merge(events, samplesize, by.x = "variable", by.y = "bin", all.x = T)
  all$Total.CNA.Freq <- all$Total.CNA / all$V1
  
  #all$class <- segclass # mute this if running for NoCluster
  all$level <- level
  all$type <- cohort
  
  All.classes <- rbind(All.classes,all)
}else{
  # Counting the events (amplification and deletions separately)
  gains <- final.df %>%  summarise_all(~ sum(.x == 1, na.rm = TRUE))
  gains <- reshape2::melt(gains)
  colnames(gains)[2] <- "amp"
  
  losses <- final.df %>%  summarise_all(~ sum(.x == -1, na.rm = TRUE))
  losses <- reshape2::melt(losses)
  colnames(losses)[2] <- "del"
  
  # Merging all the counts data
  all <- merge(gains,losses,by = "variable")
  all <- merge(all, samplesize, by.x = "variable", by.y = "bin", all.x = T)
  all$amp.freq <- all$amp / all$V1
  all$del.freq <- all$del / all$V1
  
  #all$class <- segclass # mute this if running for NoCluster
  all$level <- level
  all$type <- cohort
  
  All.classes <- rbind(All.classes,all)
  
}










