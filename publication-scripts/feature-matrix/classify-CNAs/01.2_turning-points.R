
# Learning location cutoff

# Date: October 27, 2024
# Location of segments

rm(list=ls())

centromere_path <- "./Data/hg19_centromere.txt"
telomere_path <- "./Data/hg19_telomere.txt"
chr_arms_path <- "./Data/ChrArmCoverage.txt"
segments_path <- "./Data/Mapping_SCNAs/temp4"
dist_centr_tel_outpath <- "./Data/Segments_distance_to_telomere_centromere.RData"
distribution_path <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/Data_and_plots_dist_centromere_telomere_armratio.RData"
centromere_turningpoints_outpath <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/Centromere-turningpoints.RData"
centromere_turningpoints_plot <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/TurningPoints_centromere.pdf"
telomere_turningpoints_outpath <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/Telomere-turningpoints.RData"
telomere_turningpoints_plot <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/TurningPoints_telomere.pdf"
dist_breakpoints_cent_tel <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/Distance_breakpoints_cent_tel.pdf"

# Centromere and telomere positions (UCSC)
hg19_centromere <- read.table(file = centromere_path,skip = 2)
hg19_centromere <- hg19_centromere[!hg19_centromere$V1 %in% c("chrX","chrY"),]
colnames(hg19_centromere) <- c("Chromosome","start.pos","end.pos")

hg19_telomere <- read.table(file = telomere_path,skip = 2)
hg19_telomere <- hg19_telomere[!hg19_telomere$V1 %in% c("chrX","chrY"),]
colnames(hg19_telomere) <- c("Chromosome","start.pos","end.pos")
# In the UCSC table, chr17 is missing, add the telomere position for that (10000 distance from both ends)
chr_arms <- read.delim(chr_arms_path)
chr17 <- data.frame("Chromosome" = c("chr17","chr17"),
                    "start.pos" = c(0,81195210-10000),
                    "end.pos" = c(10000,81195210))
hg19_telomere <- rbind(hg19_telomere,chr17)

rm(list = setdiff(ls(),c("hg19_centromere","hg19_telomere")))

# Segments files with arm location
segment.files <- list.files(path = segments_path, full.names = T, recursive = F)

Dist.centromere <- list()
for(file in segment.files){
  load(file)
  tumor_type <- strsplit(basename(file),"_")[[1]][1]
  for(i in 1:nrow(scna_w_armfraction)){
    chr <- paste0("chr",as.character(scna_w_armfraction$Chromosome)[i])
    start <- scna_w_armfraction$Start[i]
    end <- scna_w_armfraction$End[i]
    #Centromere position on the chromosome
    cent.start <- hg19_centromere[hg19_centromere$Chromosome == chr,"start.pos"]
    cent.end <- hg19_centromere[hg19_centromere$Chromosome == chr,"end.pos"]
    dist.to.centromere <- ifelse(cent.end < start, start - cent.end, #If Centromere is upstream of the bin
                                 ifelse(cent.start > end, cent.start - end, 0 #If centromere is downstream of the bin
                                 ))
    scna_w_armfraction$Dist.centromere[i] <- dist.to.centromere
  }
  Dist.centromere[[tumor_type]] <- scna_w_armfraction
  print(tumor_type)
}

Dist.telomere.centromere <- list()
for(tumor_type in names(Dist.centromere)){
  scna_w_armfraction <- Dist.centromere[[tumor_type]]
  for(i in 1:nrow(scna_w_armfraction)){
    chr <- paste0("chr",as.character(scna_w_armfraction$Chromosome)[i])
    start <- scna_w_armfraction$Start[i]
    end <- scna_w_armfraction$End[i]
    #Telomere position for the chromosome
    telomeres.chr <- hg19_telomere[hg19_telomere == chr,]
    # Calculate distance for each telemore
    distances <- c()
    for(j in 1:nrow(telomeres.chr)){
      tel.start <- telomeres.chr$start.pos[j]
      tel.end <- telomeres.chr$end.pos[j]
      dist.to.tel <- ifelse(tel.end < start, start - tel.end, #If FGS is upstream of the bin
                            ifelse(tel.start > end, tel.start - end, 0 #If FGS is downstream of the bin
                            ))
      distances <- c(distances,dist.to.tel)}
    #Find the closest bin
    dist.to.closest.telomere <- as.numeric(min(distances))
    scna_w_armfraction$Dist.closest.telomere[i] <- dist.to.closest.telomere
  }
  Dist.telomere.centromere[[tumor_type]] <- scna_w_armfraction
  print(tumor_type)
}

save(Dist.telomere.centromere, file = dist_centr_tel_outpath)

# # Distribution of absolute distances and arm-ratio
# 
# rm(list = ls())
# 
# # Required libraries
# library(ggplot2)
# 
# load("./Data/Segments_distance_to_telomere_centromere.RData")
# # To calculate the distance ratio
# chr_arms <- read.delim("./Data/ChrArmCoverage.txt")
# chr_arms$armlength <- (chr_arms$end_bin - chr_arms$start_bin)
# chr_arms <- chr_arms[,5:6]
# colnames(chr_arms) <- c("Arm","armlength")
# 
# plots <- list()
# Data <- list()
# for(cohort in names(Dist.telomere.centromere)){
#   segments <- Dist.telomere.centromere[[cohort]]
#   segments$ArmII <- str_replace_all(segments$Arm, "[:digit:]", "")
#   segments.new <- c()
#   # Add the distance ratio (distance to centromere/telomere divided by arm length)
#   for(arm in c("p","q","p-q")){
#     if(arm %in% c("p","q")){
#       m <- segments[segments$ArmII == arm,]
#       m <- merge(m,chr_arms,by="Arm")
#       m$ArmRatio.centromere <- m$Dist.centromere / m$armlength
#       m$ArmRatio.telomere <- m$Dist.closest.telomere / m$armlength
#       segments.new <- rbind(segments.new,m)}
#     else{
#       m <- segments[segments$ArmII == arm,]
#       # Find the arm where segment covers more so closer to the telomere
#       m$Whicharm <- colnames(m[,c("ArmFraction.p","ArmFraction.q")])[apply(m[,c("ArmFraction.p","ArmFraction.q")],1, which.max)]
#       m$Arm.new <- paste0(m$Chromosome,unlist(lapply(m$Whicharm, function(x) strsplit(x,"\\.")[[1]][2])))
#       m <- merge(m,chr_arms,by.x = "Arm.new",by.y="Arm")
#       m$ArmRatio.centromere <- NA
#       m$ArmRatio.telomere <- m$Dist.closest.telomere / m$armlength
#       segments.new <- rbind(segments.new,m[,!colnames(m) %in% c("Arm.new","Whicharm")])}
#     print(arm)
#   }
#   Data[[cohort]] <- segments.new
# 
#   # Plotting the distance and distance-ratio distribution -  Density plot
#   plot.data <- data.frame("DistanceRatio" = c(segments.new$ArmRatio.centromere,segments.new$ArmRatio.telomere),
#                           "Distance" = c(segments.new$Dist.centromere,segments.new$Dist.closest.telomere),
#                           "Chromosome" = c(segments.new$Chromosome,segments.new$Chromosome),
#                           "Name" = c(rep("Centromere",length(segments.new$ArmRatio.centromere)),
#                                      rep("Telomere",length(segments.new$ArmRatio.telomere))))
#   hist.dist <- ggplot(plot.data, aes(Distance, color = Name, fill = Name)) + geom_histogram(position = "identity", alpha = .5) +
#     scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
#     facet_wrap(~Chromosome) + theme_bw() + labs(title = cohort)
#   
#   density.dist <- ggplot(plot.data, aes(x = Distance, color = Name)) + 
#     geom_density() + scale_color_manual(values=c("#E69F00", "#56B4E9")) +
#     theme_bw() + labs(title = cohort) + facet_wrap(~Chromosome)
#   
#   hist.ratio <- ggplot(plot.data, aes(DistanceRatio, color = Name, fill = Name)) + geom_histogram(position = "identity", alpha = .5) +
#     scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
#     facet_wrap(~Chromosome) + theme_bw() + labs(title = cohort)
#   
#   density.ratio <- ggplot(plot.data, aes(x = DistanceRatio, color = Name)) + 
#     geom_density() + scale_color_manual(values=c("#E69F00", "#56B4E9")) +
#     theme_bw() + labs(title = cohort) + facet_wrap(~Chromosome)
#   
#   plots[[cohort]][["Histogram"]][["Absolute"]] <- hist.dist
#   plots[[cohort]][["Histogram"]][["Ratio"]] <- hist.ratio
#   plots[[cohort]][["Density"]][["Absolute"]] <- density.dist
#   plots[[cohort]][["Density"]][["Ratio"]] <- density.ratio
# 
#     
# }
# 
# save(Data,
#      plots, 
#      file = "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/Data_and_plots_dist_centromere_telomere_armratio.RData")



# Find distribution peak to define cutoff
load(distribution_path)
rm(list = setdiff(ls(),"Data"))

# Centromere
stat <- c()
plots <- list()

for(cohort in names(Data)){
  segments <- Data[[cohort]]
  # Separately for each chromosome
  chr.plots <- list()
  #chr.plots.ratio <- list()
  for(chr in unique(segments$Chromosome)){
    segments.chr <- segments[segments$Chromosome == chr,]
    # Centromere - Absolute distance
    Da = density(segments.chr$Dist.centromere)
    DeltaY = diff(Da$y)
    Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
    peaks <- density(segments.chr$Dist.centromere, na.rm = T)$x[Turns]
    p <- ggplot(segments.chr, aes(Dist.centromere)) + 
      geom_density() + geom_vline(xintercept = peaks[1]) +
      theme_bw() + xlab("Distance to centromere") + ylab("density") + labs(title = paste0("chr",chr))
    chr.plots[[as.character(chr)]] <- p
    stat <- rbind(stat,c(cohort,chr,peaks[1],peaks[2],"Absolute"))
    
    # # Centromere - Arm ratio
    # Da = density(segments.chr$ArmRatio.centromere,na.rm = T)
    # DeltaY = diff(Da$y)
    # Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
    # peaks <- density(segments.chr$ArmRatio.centromere, na.rm = T)$x[Turns]
    # p <- ggplot(segments.chr, aes(ArmRatio.centromere)) + 
    #   geom_density() + geom_vline(xintercept = peaks[1]) +
    #   theme_bw() + xlab("Arm ratio - centromere") + ylab("density")
    # chr.plots.ratio[[as.character(chr)]] <- p
    # stat <- rbind(stat,c(cohort,chr,peaks[1],peaks[2],"Ratio"))
  }
  cohort.plot <- ggarrange(plotlist = chr.plots)
  cohort.plot <- annotate_figure(cohort.plot, top = text_grob(cohort,face = "bold", size = 10))
  plots[["Centromere"]][["Absolute"]][[cohort]] <- cohort.plot
  
  # cohort.plot.ratio <- ggarrange(plotlist = chr.plots.ratio)
  # plots[["Centromere"]][["Ratio"]][[cohort]] <- cohort.plot.ratio

}
stat <- as.data.frame(stat)
colnames(stat) <- c("Cohort","Chromosome","TP1","TP2","Value")
stat[,3:4] <- sapply(stat[,3:4],as.numeric)
save(stat,plots, file = centromere_turningpoints_outpath)

# Telomere
rm(list = setdiff(ls(),"Data"))
stat <- c()
plots <- list()

for(cohort in names(Data)){
  segments <- Data[[cohort]]
  # Separately for each chromosome
  chr.plots <- list()
  #chr.plots.ratio <- list()
  for(chr in unique(segments$Chromosome)){
    segments.chr <- segments[segments$Chromosome == chr,]
    # Telomere - Absolute distance
    Da = density(segments.chr$Dist.closest.telomere)
    DeltaY = diff(Da$y)
    Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
    peaks <- density(segments.chr$Dist.closest.telomere, na.rm = T)$x[Turns]
    p <- ggplot(segments.chr, aes(Dist.closest.telomere)) + 
      geom_density() + geom_vline(xintercept = peaks[1]) +
      theme_bw() + xlab("Distance to telomere") + ylab("density") + labs(title = paste0("chr",chr))
    chr.plots[[as.character(chr)]] <- p
    stat <- rbind(stat,c(cohort,chr,peaks[1],peaks[2],"Absolute"))
    
    # # Telomere - Arm ratio
    # Da = density(segments.chr$ArmRatio.telomere,na.rm = T)
    # DeltaY = diff(Da$y)
    # Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
    # peaks <- density(segments.chr$ArmRatio.telomere, na.rm = T)$x[Turns]
    # p <- ggplot(segments.chr, aes(ArmRatio.telomere)) + 
    #   geom_density() + geom_vline(xintercept = peaks[1]) +
    #   theme_bw() + xlab("Arm ratio - telomere") + ylab("density")
    # chr.plots.ratio[[as.character(chr)]] <- p
    # stat <- rbind(stat,c(cohort,chr,peaks[1],peaks[2],"Ratio"))
  }
  cohort.plot <- ggarrange(plotlist = chr.plots)
  cohort.plot <- annotate_figure(cohort.plot, top = text_grob(cohort,face = "bold", size = 10))
  plots[["Telomere"]][["Absolute"]][[cohort]] <- cohort.plot
  
  # cohort.plot.ratio <- ggarrange(plotlist = chr.plots.ratio)
  # plots[["Telomere"]][["Ratio"]][[cohort]] <- cohort.plot.ratio
  
}
stat <- as.data.frame(stat)
colnames(stat) <- c("Cohort","Chromosome","TP1","TP2","Value")
stat[,3:4] <- sapply(stat[,3:4],as.numeric)
save(stat,plots,file = telomere_turningpoints_outpath)

# Date: October 28, 2024
# Visualization of peak points

rm(list=ls())

setwd("/mnt/fabiogokce/gokce")

# 1. Centromere
load(centromere_turningpoints_outpath)
pdf(centromere_turningpoints_plot,width = 15, height = 9)
plots$Centromere$Absolute$BRCA
plots$Centromere$Absolute$COADREAD
plots$Centromere$Absolute$ESCA
plots$Centromere$Absolute$GBMLGG
plots$Centromere$Absolute$KIRC
plots$Centromere$Absolute$KIRP
plots$Centromere$Absolute$LUAD
plots$Centromere$Absolute$LUSC
plots$Centromere$Absolute$OV
plots$Centromere$Absolute$PAAD
plots$Centromere$Absolute$STAD
dev.off()

stat$Chromosome <- factor(stat$Chromosome,
                          levels = as.character(seq(1,22)))
m <- stat[stat$Value == "Absolute",1:4]
m <- melt(m)
m$Mb <- m$value / 1000000


p1<-ggplot(m[m$variable == "TP1",], aes(x = Chromosome, y = Mb)) + 
  geom_boxplot() +
  geom_point() +
  theme_bw() + labs(title = "Distance to centromere")

# 2. Telomere
load(telomere_turningpoints_outpath)
pdf(telomere_turningpoints_plot,width = 15, height = 9)
plots$Telomere$Absolute$BRCA
plots$Telomere$Absolute$COADREAD
plots$Telomere$Absolute$ESCA
plots$Telomere$Absolute$GBMLGG
plots$Telomere$Absolute$KIRC
plots$Telomere$Absolute$KIRP
plots$Telomere$Absolute$LUAD
plots$Telomere$Absolute$LUSC
plots$Telomere$Absolute$OV
plots$Telomere$Absolute$PAAD
plots$Telomere$Absolute$STAD
dev.off()

stat$Chromosome <- factor(stat$Chromosome,
                          levels = as.character(seq(1,22)))
m <- stat[stat$Value == "Absolute",1:4]
m <- melt(m)
m$Mb <- m$value / 1000000

p2<-ggplot(m[m$variable == "TP1",], aes(x = Chromosome, y = Mb)) + 
  geom_boxplot() +
  geom_point() +
  theme_bw() + labs(title = "Distance to telomere")

p <- ggarrange(p1,p2, nrow = 2)
pdf(dist_breakpoints_cent_tel, width = 10)
p
dev.off()
# # Decision: Mean of first turning point for each chromosome (among the set of cohorts)
# 
# # Date: October 29, 2024
# # Considering the closest segments to centromere and telomere for the distribution (for each arm separately)
# # Absolute distances
# 
# rm(list = ls())
# 
# # Required libraries
# library(ggplot2)
# library(dplyr)
# 
# load("./Data/Segments_distance_to_telomere_centromere.RData")
# 
# 
# plots <- list()
# stat <- c()
# for(cohort in names(Dist.telomere.centromere)){
#   segments <- Dist.telomere.centromere[[cohort]]
#   segments$Arm <- ifelse(segments$Arm == "p-q", paste0(segments$Chromosome,segments$Arm), segments$Arm)
#   
#   # Closest segments in each patient
#   data.centromere <- segments %>% group_by(Sample,Arm) %>% 
#     summarise("Min" = min(Dist.centromere))
#   data.centromere$Chromosome <- as.character(unlist(lapply(data.centromere$Arm, function(x) gsub("-","",gsub("[^0-9.-]", "", x)))))
#   data.centromere$Name <- "Dist.centromere"
#   data.telomere <- segments %>% group_by(Sample,Arm) %>% 
#     summarise("Min" = min(Dist.closest.telomere))
#   data.telomere$Chromosome <- as.character(unlist(lapply(data.telomere$Arm, function(x) gsub("-","",gsub("[^0-9.-]", "", x)))))
#   data.telomere$Name <- "Dist.telomere"
#   
#   data <- rbind(data.centromere,data.telomere)
#   
# 
#   hist.dist <- ggplot(data, aes(Min, color = Name, fill = Name)) + geom_histogram(position = "identity", alpha = .5) +
#     scale_color_manual(values=c("#E69F00", "#56B4E9")) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
#     facet_wrap(~Chromosome) + theme_bw() + labs(title = cohort)
#   
#   density.dist <- ggplot(data, aes(x = Min, color = Name)) + 
#     geom_density() + scale_color_manual(values=c("#E69F00", "#56B4E9")) +
#     theme_bw() + labs(title = cohort) + facet_wrap(~Chromosome, scales = "free")
#   
#   
#   plots[[cohort]][["Histogram"]][["Absolute"]] <- hist.dist
#   plots[[cohort]][["Density"]][["Absolute"]] <- density.dist
#   
#   # Turning points
#   for(chr in unique(data$Chromosome)){
#     #Centromere
#     segments.chr <- data.centromere[data.centromere$Chromosome == chr,]
#     Da = density(segments.chr$Min)
#     DeltaY = diff(Da$y)
#     Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
#     peaks <- density(segments.chr$Min, na.rm = T)$x[Turns]
#     stat <- rbind(stat,c(cohort,chr,peaks[1],peaks[2],"Centromere"))
#     #Telomere
#     segments.chr <- data.telomere[data.telomere$Chromosome == chr,]
#     Da = density(segments.chr$Min)
#     DeltaY = diff(Da$y)
#     Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
#     peaks <- density(segments.chr$Min, na.rm = T)$x[Turns]
#     stat <- rbind(stat,c(cohort,chr,peaks[1],peaks[2],"Telomere"))
#     }
# 
# }
# 
# stat <- as.data.frame(stat)
# colnames(stat) <- c("Cohort","Chromosome","TP1","TP2","Group")
# stat[,3:4] <- sapply(stat[,3:4],as.numeric)



