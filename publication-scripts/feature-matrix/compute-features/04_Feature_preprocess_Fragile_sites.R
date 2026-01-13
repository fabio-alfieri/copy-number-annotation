
#Fragile sites

# Set working directory 
setwd("/mnt/fabiogokce/gokce") #If working on workstation

#Required libraries
library(openxlsx)

fragile.sites <- read.xlsx("./Downloads/fragile_sites_li_2020.xlsx", startRow = 3)

load("./Data/All_levels_backbonetables.RData")
Dist.fragile.sites <- list()
for(level in names(chr_backbone_namesfixed)){
  if(level %in% c("Arm","Chromosome"))next
  coordinates <- do.call(rbind,chr_backbone_namesfixed[[level]])
  coordinates$bin <- paste(coordinates$chr,coordinates$bin,sep = "_")
  coordinates[,3:4] <- sapply(coordinates[,3:4], as.numeric)
  df <- c()
  for(i in 1:nrow(coordinates)){
    chr <- paste0("chr",as.character(coordinates$chr)[i])
    start <- coordinates$start_bin[i]
    end <- coordinates$end_bin[i]
    #Distance to the closest fragile site
    fragile.sites.chr <- fragile.sites[fragile.sites$chrom == chr,]
    if(nrow(fragile.sites.chr) == 0)next
    df.bin <- c()
    for(j in 1:nrow(fragile.sites.chr)){
      fs.start <- fragile.sites.chr$start[j]
      fs.end <- fragile.sites.chr$end[j]
      dist.to.closest.FGS <- ifelse(fs.end < start, start - fs.end, #If FGS is upstream of the bin
                                   ifelse(fs.start > end, fs.start - end, 0 #If FGS is downstream of the bin
                                          ))
      row <- c(chr,coordinates$bin[i],start,end,fs.start,fs.end,dist.to.closest.FGS)
      df.bin <- rbind(df.bin,row)}
    #Find the closest bin
    df.bin <- as.data.frame(df.bin)
    df.bin$V7 <- as.numeric(df.bin$V7)
    df.bin <- df.bin[order(df.bin$V7, decreasing = F),]
    df <- rbind(df,df.bin[1,])}
  df <- as.data.frame(df)
  colnames(df) <- c("chr","bin","start_bin","end_bin","start_FGS","end_FGS","dist.to.closest.FGS")
  df[,3:7] <- sapply(df[,3:7], as.numeric)
  Dist.fragile.sites[[level]] <- df
  
}

save(Dist.fragile.sites,
     file = "./Data/Distance_to_closest_fragile_site.RData")

# Feature optimization

#Feature optimization
#Fragile sites

# Set working directory 
setwd("/mnt/fabiogokce/gokce") #If working on workstation

#Required libraries
library(openxlsx)

fragile.sites <- read.xlsx("./Downloads/fragile_sites_li_2020.xlsx", startRow = 3)

load("/mnt/fabiogokce/gokce/Data/All_levels_backbonetables.RData")
Dist.fragile.sites <- list()

for(level in names(chr_backbone_namesfixed)){
  if(level %in% c("Arm","Chromosome"))next
  coordinates <- do.call(rbind,chr_backbone_namesfixed[[level]])
  coordinates$bin <- paste(coordinates$chr,coordinates$bin,sep = "_")
  coordinates[,3:4] <- sapply(coordinates[,3:4], as.numeric)
  df <- c()
  for(i in 1:nrow(coordinates)){
    chr <- paste0("chr",as.character(coordinates$chr)[i])
    start <- coordinates$start_bin[i]
    end <- coordinates$end_bin[i]
    #Distance to the closest fragile site
    fragile.sites.chr <- fragile.sites[fragile.sites$chrom == chr,]
    if(nrow(fragile.sites.chr) == 0)next
    df.bin <- c()
    for(j in 1:nrow(fragile.sites.chr)){
      fs.start <- fragile.sites.chr$start[j]
      fs.end <- fragile.sites.chr$end[j]
      dist.to.closest.FGS <- ifelse(fs.end < start, start - fs.end, #If FGS is upstream of the bin
                                    ifelse(fs.start > end, fs.start - end, 0 #If FGS is downstream of the bin
                                    ))
      row <- c(chr,coordinates$bin[i],start,end,fs.start,fs.end,dist.to.closest.FGS)
      df.bin <- rbind(df.bin,row)}
    df <- rbind(df,df.bin)}
  df <- as.data.frame(df)
  colnames(df) <- c("chr","bin","start_bin","end_bin","start_FGS","end_FGS","dist.to.FGS")
  df[,3:7] <- sapply(df[,3:7], as.numeric)
  Dist.fragile.sites[[level]] <- df}

rm(list = setdiff(ls(),"Dist.fragile.sites"))

#Weighted distance scores to the fragile sites
Weighted.dist.scores <- list()
for(level in names(Dist.fragile.sites)){
  df <- Dist.fragile.sites[[level]]
  #Giving weight based on the distance
  df$weight <- ifelse(df$dist.to.FGS > 10000000, 0,
                      ifelse(df$dist.to.FGS <= 10000000 & df$dist.to.FGS > 8000000, 2,
                             ifelse(df$dist.to.FGS <= 8000000 & df$dist.to.FGS > 6000000, 4,
                                    ifelse(df$dist.to.FGS <= 6000000 & df$dist.to.FGS > 4000000, 6,
                                           ifelse(df$dist.to.FGS <= 4000000 & df$dist.to.FGS > 2000000, 8,10)))))
  df.weight <- df %>% group_by(bin) %>% summarise(Dist.FGS.weighted = sum(weight))
  Weighted.dist.scores[[level]] <- df.weight}

save(Weighted.dist.scores,
     file = "./Data/FeatureOptimization_Fragile_sites_weighted_distance_scores.RData")



