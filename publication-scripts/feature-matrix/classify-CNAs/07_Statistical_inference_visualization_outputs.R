
# Visualization of output files and statistical analysis on the data
# Outputs from Pipeline.R
# Arm fraction cutoff is 90%

# Clean the environment
rm(list = ls())
gc()

packages <- c(
  "openxlsx", "dplyr", "ggplot2", "ggpubr", "mclust", "ComplexHeatmap"
)

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

# A. Statistical Analysis: The cases in segments files when they are grouped by length and centromere-telomere location
#-----
cluster_2_path <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/ClusterII/armfraction_0.9.RData"
cluster2_outpath <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/Segment_cluster_cases.xlsx"
final_cluster_path <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/Final_clustering/armfraction_0.9.RData"
final_cluster_outpath <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/Percentage_of_CNAsegments_cluster1_across_cohorts.pdf"
final_cluster_outpath_2 <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/Percentage_of_CNAsegments_cluster1_across_cohorts_perchr.pdf"
percentage_CNAsegments_across_cohorts_ampl_del_plot <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/Percentage_of_CNAsegments_cluster1_across_cohorts_amp_del.pdf"
percentage_CNAsegments_across_cohorts_perchr_ampl_del_plot <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/Percentage_of_CNAsegments_cluster1_across_cohorts_perchr_amp_del.pdf"
heatmap_per_segmentcluster <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/Heatmap_scaled_per_segmentclusters.pdf"
freq_of_top5_path <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/Freq_of_top5_in_segmentclusters_chromosomes.pdf"
backbone_path <- "./Data/All_levels_backbonetables.RData"
totalCNA_biosamples <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/CNA_scores_11cohorts_interstitial_cov_zero__TotalCNA_divided_binsamples.pdf"
CNA_scores_cohorts <- "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/CNA_scores_11cohorts_interstitial_cov_zero_Arm_Chr_levels_TotalCNA_divided_binsamples.pdf"


load(cluster_2_path) # Output.List.new
all.data <- c()
for(cohort in names(Output.List.new)){
  scna <- Output.List.new[[cohort]]
  # If a segment spans on both arms, change "Centromere-bounded" to Centromere-spanning"
  scna$Centromere.class <- ifelse(scna$ArmII == "p-q", "Centromere-spanning",scna$Centromere.class)
  # Remove segments under 1kb
  scna <- scna[scna$CLUSTER1 != "Not classified",]
  scna$cohort <- cohort
  all.data <- rbind(all.data, scna)
}

stat <- all.data %>% group_by(CLUSTER1,Centromere.class,Telomere.class) %>% summarise("Freq" = n())
stat.cohort <- all.data %>% group_by(CLUSTER1,Centromere.class,Telomere.class,cohort) %>% summarise("Freq" = n())
write.xlsx(stat, file = cluster2_outpath)
#-----

# B. Segments after the clustering - Final clusters and number of segments in each cluster
# B.1. Cluster I
# B.1.1 ignore the chromosome information
#-----
  
load(final_cluster_path) # where arm fraction to separate mid-length and arm/chromosome-level segments is 90%

Info <- c() # Ignore the Info loaded with the datasets - make new one where you remove the diploid segments!
for(cohort in names(Cluster.data)){
  m <- Cluster.data[[cohort]]
  m$event <- ifelse(m$Segment_Mean > 0.2, "Amplification",
                    ifelse(m$Segment_Mean < -0.2, "Deletion","Diploid"))
  numbers <- m %>% group_by(Chromosome,Final.cluster,event) %>% summarise("Number" = n())
  numbers <- as.data.frame(numbers)
  numbers$cohort <- cohort
  Info <- rbind(Info,numbers)
}

Info$CLUSTER1 <- unlist(lapply(Info$Final.cluster, function(x) strsplit(as.character(x),"_")[[1]][1]))

# Ignore chromosome information for this section
data <- Info[Info$event != "Diploid",]
data <- data %>% group_by(cohort,CLUSTER1) %>% summarise("Total" = sum(Number))
allsegments <- data %>% group_by(cohort) %>% summarise("Allsegments" = sum(Total))
#allsegments <- Info %>% group_by(cohort) %>% summarise("Allsegments" = sum(Number)) # the order doesn't change if you count all segments or changed ones
data <- merge(data, allsegments, by = "cohort")
data$Percentage <- (data$Total / data$Allsegments)*100

# Plotting the Percentage
plots <- list()
for(clus in unique(data$CLUSTER1)){
  m <- data[data$CLUSTER1 == clus,]
  m <- m[order(m$Percentage, decreasing = T),]
  m$cohort <- factor(m$cohort, levels = m$cohort)
  p <- ggplot(m, aes(x = cohort, y = Percentage)) +  geom_bar(stat="identity") +
    theme_bw() + labs(title = clus) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")
  plots[[clus]] <- p
}
final.p <- ggarrange(plotlist = plots)
final.p <- annotate_figure(final.p,top = text_grob("Percentage: Number of amp and del segments in a class / Total number of amp and del segments\nArm fraction - 0.9",
                                            color = "red", face = "bold", size = 10))
pdf(final_cluster_outpath,
    width = 13, height = 8)
final.p
dev.off()
#-----

# B.1.2 Chromosome information: To understand which chromosome(s) drives the overall pattern for each cohort
# Cluster I - Deeper understanding while focusing on chromosome
#-----
rm(list = setdiff(ls(),"Info"))
data <- Info[Info$event != "Diploid",]
data <- data %>% group_by(cohort,CLUSTER1,Chromosome) %>% summarise("Total" = sum(Number))
allsegments <- data %>% group_by(cohort, Chromosome) %>% summarise("Allsegments" = sum(Total))
data <- merge(data, allsegments, by = c("cohort","Chromosome"))
data$Percentage <- (data$Total / data$Allsegments)*100

# Percentage
plots <- list()
for(clus in unique(data$CLUSTER1)){
  cohort.plots <- list()
  for(cohort in as.character(unique(data$cohort))){
    m <- data[data$CLUSTER1 == clus & data$cohort == cohort,]

    m <- m[order(m$Percentage, decreasing = T),]
    m$Chromosome <- factor(m$Chromosome, levels = m$Chromosome)
    p <- ggplot(m, aes(x = Chromosome, y = Percentage)) +  geom_bar(stat="identity") +
      theme_bw() + labs(title = cohort) + xlab("Chromosomes")
    cohort.plots[[cohort]] <- p
  }
  clus.p <- ggarrange(plotlist = cohort.plots)
  clus.p <- annotate_figure(clus.p,
                            top = text_grob(clus,color = "red",face = "bold",size = 10))
  plots[[clus]] <- clus.p
}

pdf(final_cluster_outpath_2,
    width = 18, height = 7)
plots$`Arm-level`
plots$`Chromosome-level`
plots$`Mid-length`
plots$`Small-scale`
dev.off()
#-----

# B.2 Visualization of segment clusters but separate amplification and deletion cases this time
# B.2.1 Ignore chromosome information for this section
#-----
rm(list = setdiff(ls(),"Info"))

data <- Info[Info$event != "Diploid",]
data <- data %>% group_by(cohort,CLUSTER1,event) %>% summarise("Total" = sum(Number))
allsegments <- data %>% group_by(cohort) %>% summarise("Allsegments" = sum(Total))
data <- merge(data, allsegments, by = "cohort")
data$Percentage <- (data$Total / data$Allsegments)*100

# Plotting the Percentage
plots <- list()
for(clus in unique(data$CLUSTER1)){
    m <- data[data$CLUSTER1 == clus,]
    # Order
    order.m <- m %>% group_by(cohort) %>% summarise("total" = sum(Percentage))
    order.m <- order.m[order(order.m$total, decreasing = T),]
    m$cohort <- factor(m$cohort, levels = order.m$cohort)
    p <- ggplot(m, aes(x = cohort, y = Percentage, fill = event)) +  
      geom_bar(stat="identity", position = position_dodge()) +
      scale_fill_manual(values = c("#780000","#003049")) +
      theme_bw() + labs(title = clus) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")
  plots[[clus]] <- p
}
final.p <- ggarrange(plotlist = plots, common.legend = TRUE)
final.p <- annotate_figure(final.p,top = text_grob("Percentage: Number of amp(del) segments in a class / Total number of amp and del segments\nArm fraction - 0.9",
                                                   color = "red", face = "bold", size = 10))

pdf(percentage_CNAsegments_across_cohorts_ampl_del_plot,
    width = 13, height = 8)
final.p
dev.off()
#-----

# B.2.2 Chromosome information: To understand which chromosome(s) drives the overall pattern for each cohort
#-----
rm(list = setdiff(ls(),"Info"))
data <- Info[Info$event != "Diploid",]
data <- data %>% group_by(cohort,CLUSTER1,Chromosome,event) %>% summarise("Total" = sum(Number))
allsegments <- data %>% group_by(cohort, Chromosome) %>% summarise("Allsegments" = sum(Total))
data <- merge(data, allsegments, by = c("cohort","Chromosome"))
data$Percentage <- (data$Total / data$Allsegments)*100

# Percentage
plots <- list()
for(clus in unique(data$CLUSTER1)){
  cohort.plots <- list()
  for(cohort in as.character(unique(data$cohort))){
    m <- data[data$CLUSTER1 == clus & data$cohort == cohort,]
    # Order chromosomes
    order.m <- m %>% group_by(cohort,Chromosome) %>% summarise("total" = sum(Percentage))
    order.m <- order.m[order(order.m$total, decreasing = T),]
    m$Chromosome <- factor(m$Chromosome, levels = order.m$Chromosome)
    p <- ggplot(m, aes(x = Chromosome, y = Percentage, fill = event)) + 
      geom_bar(stat="identity", position = position_dodge()) +
      scale_fill_manual(values = c("#780000","#003049")) +
      theme_bw() + labs(title = cohort) + xlab("Chromosomes")
    cohort.plots[[cohort]] <- p
  }
  clus.p <- ggarrange(plotlist = cohort.plots, common.legend = TRUE)
  clus.p <- annotate_figure(clus.p,
                            top = text_grob(clus,color = "red",face = "bold",size = 10))
  plots[[clus]] <- clus.p
}

pdf(percentage_CNAsegments_across_cohorts_perchr_ampl_del_plot,
    width = 18, height = 7)
plots$`Arm-level`
plots$`Chromosome-level`
plots$`Mid-length`
plots$`Small-scale`
dev.off()
#-----

# A.3 Clustering chromosomes based on the number of each cluster
# A.3.1 Heatmap
#----
rm(list = setdiff(ls(),"Info"))
data <- Info[Info$event != "Diploid",]
data <- data %>% group_by(cohort,CLUSTER1,Chromosome) %>% summarise("Total" = sum(Number))
allsegments <- data %>% group_by(cohort, Chromosome) %>% summarise("Allsegments" = sum(Total))
data <- merge(data, allsegments, by = c("cohort","Chromosome"))
data$Percentage <- (data$Total / data$Allsegments)*100
data$sample <- paste(data$cohort, data$Chromosome, sep = "_")
counts <- dcast(data, sample ~ CLUSTER1, value.var = "Percentage")
#counts[is.na(counts)] <- 0
rownames(counts) <- counts$sample
counts <- counts[,-1]
counts.scaled <- as.data.frame(apply(counts, 2, scale))
rownames(counts.scaled) <- rownames(counts)

# Heatmap
m <- as.matrix(counts.scaled)
ht <- Heatmap(m, rect_gp = gpar(col = "white", lwd = 0.5),
        row_names_gp = gpar(fontsize = 3), 
        name = "Scaled percentage",
        heatmap_width = unit(50,"mm"),
        heatmap_height = unit(300,"mm"))

pdf(heatmap_per_segmentcluster, height = 15)
ht
dev.off()
#-----

# A.3.2 Frequency of chromosomes pop up in top5 of each cluster
#-----
rm(list = setdiff(ls(),"Info"))
data <- Info[Info$event != "Diploid",]
data <- data %>% group_by(cohort,CLUSTER1,Chromosome) %>% summarise("Total" = sum(Number))
allsegments <- data %>% group_by(cohort, Chromosome) %>% summarise("Allsegments" = sum(Total))
data <- merge(data, allsegments, by = c("cohort","Chromosome"))
data$Percentage <- (data$Total / data$Allsegments)*100

stat <- c()
for(cohort in unique(data$cohort)){
  for(clus in unique(data$CLUSTER1)){
    m <- data[data$cohort == cohort & data$CLUSTER1 == clus,]
    m <- m[order(m$Percentage, decreasing = T),]
    stat <- rbind(stat, m[1:5,])
  }
}

freq.m <- stat %>% group_by(Chromosome,CLUSTER1) %>% summarise("Freq" = n(),
                                                               "Cohorts" = paste(cohort,collapse = ";"))
freq.m$Chromosome <- as.character(freq.m$Chromosome)
freq.m$Chromosome <- factor(freq.m$Chromosome, levels = as.character(seq(1,22)))
freq.m$label <- ifelse(freq.m$Freq <= 2,freq.m$Cohorts, freq.m$Freq )

p<- ggplot(freq.m, aes(x=Chromosome, y=Freq, fill=CLUSTER1, label = label)) +
  geom_bar(stat="identity", position=position_dodge2(preserve = "single")) +
  geom_text(position = position_dodge2(width = 0.9, preserve = "single"), 
            vjust=0.5, angle = 90, hjust = -0.5) + 
  scale_fill_manual(values = c("#c8553d","#f28f3b","#ffd5c2","#588b8b")) + 
  theme_classic()

pdf(freq_of_top5_path, 
    width = 15, height = 6)
p
dev.off()
#------


# B. Amplification, deletion or Total CNA frequencies - visualization
# B.1. Total CNA - divided by bin-samples
#-----
# Date: November 26, 2024
# Visualization of CNA frequencies

rm(list = ls())
gc()

# Color codes
vline <- "#d6ccc2"

# Bins
load(backbone_path)

# CNA scores
output.files <- list.files(path = "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/CNA_frequencies",
                           full.names = TRUE, recursive = FALSE, pattern = "divided_binsamples")

cohorts <- c("BRCA","COADREAD","ESCA","GBMLGG","KIRC","KIRP","LUAD","LUSC","OV","PAAD","STAD")

# Data preparation
plot.datasets <- list()
for(file in output.files){
  covThr <- paste(strsplit(basename(file),"_")[[1]][3],strsplit(basename(file),"_")[[1]][4],sep = "_")
  load(file)
  
  for(level in names(All.class.data)){
    # Bin coordinates
    coordinates <- do.call(rbind,chr_backbone_namesfixed[[level]])
    coordinates$bin <- paste(coordinates$chr, coordinates$bin, sep = "_")
    
    for(cohort in cohorts){
      cna.data <- coordinates
      for(seg.class in names(All.class.data[[level]])){
        m <- All.class.data[[level]][[seg.class]]
        m.cohort <- m[m$type == cohort,]
        colnames(m.cohort)[4] <- seg.class
        
        # Skip this part if the output files with "divided_binsamples"
        #m.cohort$Total.sample <- m.cohort$amp + m.cohort$del + m.cohort$diploid
        #m.cohort <- m.cohort[m.cohort$Total.sample >= 10,]
        # Skip part ends
        
        # Amplification scores
        cna.data <- merge(cna.data,m.cohort[,c(1,4)], by = "bin", by.y = "variable", all.x = TRUE)
      }
      plot.datasets[[covThr]][[level]][[cohort]] <- cna.data
    }}
  
}

# Plots - all genome
plots <- list()
#seg.classes <- c("Small-scale_Interstitial","Mid-length_Interstitial")
seg.classes <- c("Arm-level","Chromosome-level")

for(covThr in names(plot.datasets)){
  for(level in names(plot.datasets[[covThr]])){
    for(cohort in cohorts){
      cohort.p <- list()
      m <- plot.datasets[[covThr]][[level]][[cohort]]
      
      # Order m based on bin
      m$bin.num <- as.numeric(unlist(lapply(as.character(m$bin), function(x) strsplit(x,"_")[[1]][2])))
      m_ordered <- c()
      x_breaks <- c()
      for(chr in seq(1,22)){
        sub <- m[m$chr == chr,]
        sub <- sub[order(sub$bin.num, decreasing = F),]
        m_ordered <- rbind(m_ordered,sub)
        x_breaks <- c(x_breaks,as.character(sub[1,]$bin))
      }
      
      m_ordered$bin <- factor(m_ordered$bin, levels = unique(m_ordered$bin))
      melted <- reshape2::melt(m_ordered[,c(1,5:13)])
      melted <- melted[melted$variable %in% seg.classes,]
      
      p <- ggplot(melted, aes(x = bin, y = value, group=variable)) + 
        geom_line(aes(color = variable), na.rm = TRUE) + theme_classic() + 
        scale_color_manual(values = c("#003049","#a7c957")) +
        ylab("Total CNA frequency") + xlab("Genomic coordinates (Bins)") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_x_discrete(breaks = x_breaks) + 
        geom_vline(xintercept = x_breaks, color = vline) +
        labs(title = paste(level, cohort,sep = "-")) +
        theme(legend.position="top")
      
      plots[[covThr]][[level]][[cohort]] <- p
      
    }
    
  }
  
}

pdf(totalCNA_biosamples,width = 20, height = 7)
plots$covThr_zero$`0.1Mbp`$BRCA
plots$covThr_zero$`0.1Mbp`$COADREAD
plots$covThr_zero$`0.1Mbp`$ESCA
plots$covThr_zero$`0.1Mbp`$GBMLGG
plots$covThr_zero$`0.1Mbp`$KIRC
plots$covThr_zero$`0.1Mbp`$KIRP
plots$covThr_zero$`0.1Mbp`$LUAD
plots$covThr_zero$`0.1Mbp`$LUSC
plots$covThr_zero$`0.1Mbp`$OV
plots$covThr_zero$`0.1Mbp`$PAAD
plots$covThr_zero$`0.1Mbp`$STAD
dev.off()

pdf(CNA_scores_cohorts,width = 20, height = 7)
plots$covThr_zero$`0.1Mbp`$BRCA
plots$covThr_zero$`0.1Mbp`$COADREAD
plots$covThr_zero$`0.1Mbp`$ESCA
plots$covThr_zero$`0.1Mbp`$GBMLGG
plots$covThr_zero$`0.1Mbp`$KIRC
plots$covThr_zero$`0.1Mbp`$KIRP
plots$covThr_zero$`0.1Mbp`$LUAD
plots$covThr_zero$`0.1Mbp`$LUSC
plots$covThr_zero$`0.1Mbp`$OV
plots$covThr_zero$`0.1Mbp`$PAAD
plots$covThr_zero$`0.1Mbp`$STAD
dev.off()

# Plots - chr9

nth_element <- function(vector, starting_position, n) { 
  vector[seq(starting_position, length(vector), n)] 
}
plots <- list()
seg.classes <- c("Small-scale_Interstitial","Mid-length_Interstitial")

for(covThr in names(plot.datasets)){
  for(level in names(plot.datasets[[covThr]])){
    for(cohort in cohorts){
      m <- plot.datasets[[covThr]][[level]][[cohort]]
      
      # Order m based on bin
      m$bin.num <- as.numeric(unlist(lapply(as.character(m$bin), function(x) strsplit(x,"_")[[1]][2])))
      m_ordered <- c()
      for(chr in seq(1,22)){
        sub <- m[m$chr == chr,]
        sub <- sub[order(sub$bin.num, decreasing = F),]
        m_ordered <- rbind(m_ordered,sub)
      }
      
      m_ordered$bin <- factor(m_ordered$bin, levels = unique(m_ordered$bin))
      # Filter: Chromosome
      m_ordered <- m_ordered[m_ordered$chr == 9,]
      
      
      melted <- reshape2::melt(m_ordered[,c(1,5:13)])
      melted <- melted[melted$variable %in% seg.classes,]
      
      x_breaks <- unique(melted$bin)
      x_breaks <- nth_element(x_breaks, 1, 100)
      
      p <- ggplot(melted, aes(x = bin, y = value, group=variable)) + 
        geom_line(aes(color = variable), na.rm = TRUE) + theme_classic() + 
        scale_color_manual(values = c("#003049","#a7c957")) +
        ylab("Total CNA frequency") + xlab("Genomic coordinates (Bins)") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_x_discrete(breaks = x_breaks) + 
        labs(title = paste(level, cohort,sep = "-")) +
        theme(legend.position="top")
      
      plots[[covThr]][[level]][[cohort]] <- p
      
    }
    
  }
  
}

pdf("./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/CNA_scores_11cohorts_interstitial_cov_zero__TotalCNA_chr9_divided_binsamples.pdf",width = 20, height = 7)
plots$covThr_zero$`0.1Mbp`$BRCA
plots$covThr_zero$`0.1Mbp`$COADREAD
plots$covThr_zero$`0.1Mbp`$ESCA
plots$covThr_zero$`0.1Mbp`$GBMLGG
plots$covThr_zero$`0.1Mbp`$KIRC
plots$covThr_zero$`0.1Mbp`$KIRP
plots$covThr_zero$`0.1Mbp`$LUAD
plots$covThr_zero$`0.1Mbp`$LUSC
plots$covThr_zero$`0.1Mbp`$OV
plots$covThr_zero$`0.1Mbp`$PAAD
plots$covThr_zero$`0.1Mbp`$STAD
dev.off()

#-----

# B.2. Turning points (Peaks) for CNA distribution - Without smoothing
#-----
rm(list = ls())
gc()

load(backbone_path)

# CNA scores
output.files <- list.files(path = "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/CNA_frequencies",
                           full.names = TRUE, recursive = FALSE, pattern = "divided_binsamples")

cohorts <- c("BRCA","COADREAD","ESCA","GBMLGG","KIRC","KIRP","LUAD","LUSC","OV","PAAD","STAD")

# Data preparation
plot.datasets <- list()
for(file in output.files){
  covThr <- paste(strsplit(basename(file),"_")[[1]][3],strsplit(basename(file),"_")[[1]][4],sep = "_")
  load(file)
  
  for(level in names(All.class.data)){
    # Bin coordinates
    coordinates <- do.call(rbind,chr_backbone_namesfixed[[level]])
    coordinates$bin <- paste(coordinates$chr, coordinates$bin, sep = "_")
    
    for(cohort in cohorts){
      cna.data <- coordinates
      for(seg.class in names(All.class.data[[level]])){
        m <- All.class.data[[level]][[seg.class]]
        m.cohort <- m[m$type == cohort,]
        colnames(m.cohort)[4] <- seg.class
        
        # Skip this part if the output files with "divided_binsamples"
        #m.cohort$Total.sample <- m.cohort$amp + m.cohort$del + m.cohort$diploid
        #m.cohort <- m.cohort[m.cohort$Total.sample >= 10,]
        # Skip part ends
        
        # Amplification scores
        cna.data <- merge(cna.data,m.cohort[,c(1,4)], by = "bin", by.y = "variable", all.x = TRUE)
      }
      plot.datasets[[covThr]][[level]][[cohort]] <- cna.data
    }}
  
}

# Functions to find peaks
find_peaks <- function(x) {
  diff_sign <- sign(diff(x))
  peak_indices <- which(diff_sign[-length(diff_sign)] > 0 & diff_sign[-1] < 0)
  return(peak_indices + 1) # Adjust for shifted index
}

peak.datasets <- list()
seg.classes <- c("Small-scale_Interstitial","Mid-length_Interstitial")

for(covThr in names(plot.datasets)){
  for(level in names(plot.datasets[[covThr]])){
    for(cohort in cohorts){
      m <- plot.datasets[[covThr]][[level]][[cohort]]
      
      # Order m based on bin
      m$bin.num <- as.numeric(unlist(lapply(as.character(m$bin), function(x) strsplit(x,"_")[[1]][2])))
      m_ordered <- c()
      for(chr in seq(1,22)){
        sub <- m[m$chr == chr,]
        sub <- sub[order(sub$bin.num, decreasing = F),]
        m_ordered <- rbind(m_ordered,sub)
      }
      
      m_ordered$bin <- factor(m_ordered$bin, levels = unique(m_ordered$bin))
      
      all.peaks <- c()
      for(chr in seq(1,22)){
        # Filter: Chromosome
        m_ordered.chr <- m_ordered[m_ordered$chr == chr,]
        m_ordered.chr <- m_ordered.chr[,c("bin","chr","start_bin","end_bin",seg.classes)]
        # Find peaks for Small-scale interstitial
        peak_indices <- find_peaks(m_ordered.chr$`Small-scale_Interstitial`)
        m_ordered.chr$peak.small <- FALSE
        m_ordered.chr$peak.small[peak_indices] <- TRUE
        
        # Find peaks for Mid-length interstitial
        peak_indices <- find_peaks(m_ordered.chr$`Mid-length_Interstitial`)
        m_ordered.chr$peak.mid <- FALSE
        m_ordered.chr$peak.mid[peak_indices] <- TRUE
        
        all.peaks <- rbind(all.peaks,m_ordered.chr)
      }
      
      peak.datasets[[covThr]][[level]][[cohort]] <- all.peaks
    }
    
  }
  
}
#-----


# to be checked
#------
p1 <- ggplot(Info, aes(x = Cohort, y = Size)) +  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Number of segments") + xlab("") +
  labs(title = paste0("Arm fraction - ",armf)) + facet_wrap(~Cluster, scales = "free")

  sum <- Info %>% group_by(Cluster) %>% summarise("Total" = sum(Size))
  sum$Per <- round((sum$Total / sum(sum$Total)) * 100, digits = 2)
  sum$Text <- paste0(sum$Total, "\n(",sum$Per,")")
  p2 <- ggplot(sum, aes(x = Cluster, y = Total)) +  geom_bar(stat="identity") +
    theme_bw() + ylab("Number of segments") + xlab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_text(aes(label=Text), vjust=0)
  p <- ggarrange(p1, p2)


pdf("./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/NumSegments_armfraction_90percent.pdf", width = 18, height = 10)
num.segments$`0.9`
dev.off()


# Date: November 14, 2024
# Visualization of CNA frequencies

rm(list = ls())
gc()

# Color codes
vline <- "#d6ccc2"

# Bins
load(backbone_path)

# CNA scores
# output.files <- list.files(path = "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/CNA_frequencies",
#                            full.names = TRUE, recursive = FALSE)
output.files <- list.files(path = "./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/CNA_frequencies",
                           full.names = TRUE, recursive = FALSE, pattern = "divided_binsamples")

cohorts <- c("BRCA","COADREAD","ESCA","GBMLGG","KIRC","KIRP","LUAD","LUSC","OV","PAAD","STAD")

# Data preparation
plot.datasets <- list()
for(file in output.files){
  covThr <- paste(strsplit(basename(file),"_")[[1]][3],strsplit(basename(file),"_")[[1]][4],sep = "_")
  load(file)
  
  for(level in names(All.class.data)){
    # Bin coordinates
    coordinates <- do.call(rbind,chr_backbone_namesfixed[[level]])
    coordinates$bin <- paste(coordinates$chr, coordinates$bin, sep = "_")
    
    for(cohort in cohorts){
      cna.data <- coordinates
      cna.data.del <- coordinates
      for(seg.class in names(All.class.data[[level]])){
        m <- All.class.data[[level]][[seg.class]]
        m.cohort <- m[m$type == cohort,]
        
        # Skip this part if the output files with "divided_binsamples"
        #m.cohort$Total.sample <- m.cohort$amp + m.cohort$del + m.cohort$diploid
        #m.cohort <- m.cohort[m.cohort$Total.sample >= 10,]
        # Skip part ends
        
        # Amplification scores
        m.cohort.amp <- m.cohort[,c("variable","amp.freq")]
        colnames(m.cohort.amp)[2] <- seg.class
        cna.data <- merge(cna.data,m.cohort.amp, by.x = "bin", by.y = "variable", all.x = TRUE)
        
        # Deletion scores
        m.cohort.del <- m.cohort[,c("variable","del.freq")]
        colnames(m.cohort.del)[2] <- seg.class
        cna.data.del <- merge(cna.data.del,m.cohort.del, by.x = "bin", by.y = "variable", all.x = TRUE)
      }
      plot.datasets[[covThr]][[level]][[cohort]][["Amp"]] <- cna.data
      plot.datasets[[covThr]][[level]][[cohort]][["Del"]] <- cna.data.del
    }}
  
}


# Plots 
seg.classes <- c("Small-scale_Interstitial","Mid-length_Interstitial")
plots <- list()

for(covThr in names(plot.datasets)){
  for(level in names(plot.datasets[[covThr]])){
    for(cohort in cohorts){
      cohort.p <- list()
      for(name in c("Amp","Del")){
        m <- plot.datasets[[covThr]][[level]][[cohort]][[name]]
        
        # Order m based on bin
        m$bin.num <- as.numeric(unlist(lapply(as.character(m$bin), function(x) strsplit(x,"_")[[1]][2])))
        m_ordered <- c()
        x_breaks <- c()
        for(chr in seq(1,22)){
          sub <- m[m$chr == chr,]
          sub <- sub[order(sub$bin.num, decreasing = F),]
          m_ordered <- rbind(m_ordered,sub)
          x_breaks <- c(x_breaks,as.character(sub[1,]$bin))
        }
        
        m_ordered$bin <- factor(m_ordered$bin, levels = unique(m_ordered$bin))
        melted <- melt(m_ordered[,c(1,5:6)])
        #melted <- melt(m_ordered[,c(1,5:11)])
        melted <- melted[melted$variable %in% seg.classes,]
        
        p <- ggplot(melted, aes(x = bin, y = value, group=variable)) + 
          geom_line(aes(color = variable), na.rm = TRUE) + theme_classic() + 
          scale_color_manual(values = c("#003049","#d9d9d9")) +
          ylab(paste(name,"frequency", sep = " ")) + xlab("Genomic coordinates (Bins)") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          scale_x_discrete(breaks = x_breaks) + 
          geom_vline(xintercept = x_breaks, color = vline) +
          labs(title = paste(level, cohort, name, sep = "-")) +
          theme(legend.position="top")
        cohort.p[[name]] <- p
        
      }
      final <- ggarrange(plotlist = cohort.p, common.legend = TRUE, nrow = 2)
      plots[[covThr]][[level]][[cohort]] <- final
    }
    
  }
  
}


pdf("./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/CNA_scores_11cohorts_interstitial_cov_0.5.pdf",width = 20, height = 7)
plots$covThr_0.5.RData$`0.1Mbp`$BRCA
plots$covThr_0.5.RData$`0.1Mbp`$COADREAD
plots$covThr_0.5.RData$`0.1Mbp`$ESCA
plots$covThr_0.5.RData$`0.1Mbp`$GBMLGG
plots$covThr_0.5.RData$`0.1Mbp`$KIRC
plots$covThr_0.5.RData$`0.1Mbp`$KIRP
plots$covThr_0.5.RData$`0.1Mbp`$LUAD
plots$covThr_0.5.RData$`0.1Mbp`$LUSC
plots$covThr_0.5.RData$`0.1Mbp`$OV
plots$covThr_0.5.RData$`0.1Mbp`$PAAD
plots$covThr_0.5.RData$`0.1Mbp`$STAD
dev.off()

pdf("./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Plots/CNA_scores_11cohorts_interstitial_cov_zero_divided_binsamples.pdf",width = 20, height = 7)
plots$covThr_zero$`0.1Mbp`$BRCA
plots$covThr_zero$`0.1Mbp`$COADREAD
plots$covThr_zero$`0.1Mbp`$ESCA
plots$covThr_zero$`0.1Mbp`$GBMLGG
plots$covThr_zero$`0.1Mbp`$KIRC
plots$covThr_zero$`0.1Mbp`$KIRP
plots$covThr_zero$`0.1Mbp`$LUAD
plots$covThr_zero$`0.1Mbp`$LUSC
plots$covThr_zero$`0.1Mbp`$OV
plots$covThr_zero$`0.1Mbp`$PAAD
plots$covThr_zero$`0.1Mbp`$STAD
dev.off()

# Date: November 19, 2024 - continue
# Amplification and deletion frequencies 

rm(list = setdiff(ls(),"plot.datasets"))
gc()

#load("./Data/Feature_tables_for_ML_models_11cohorts.RData")

# Plots - Focusing on small regions
seg.classes <- c("Small-scale_Interstitial","Mid-length_Interstitial")
plots <- list()

for(covThr in names(plot.datasets)){
  for(level in names(plot.datasets[[covThr]])){
    for(cohort in cohorts){
      cohort.p <- list()
      for(name in c("Amp","Del")){
        m <- plot.datasets[[covThr]][[level]][[cohort]][[name]]
        
        # Order m based on bin
        m$bin.num <- as.numeric(unlist(lapply(as.character(m$bin), function(x) strsplit(x,"_")[[1]][2])))
        m_ordered <- c()
        x_breaks <- c()
        for(chr in seq(1,22)){
          sub <- m[m$chr == chr,]
          sub <- sub[order(sub$bin.num, decreasing = F),]
          m_ordered <- rbind(m_ordered,sub)
          x_breaks <- c(x_breaks,as.character(sub[1,]$bin))
        }
        
        m_ordered$bin <- factor(m_ordered$bin, levels = unique(m_ordered$bin))
        melted <- melt(m_ordered[,c(1,5:11)])
        melted <- melted[melted$variable %in% seg.classes,]
        
        # Focusing on very small region
        selection.df <- m_ordered[m_ordered$chr == "7",]
        selection <- paste0("7_",seq(1095,1110))
        p <- ggplot(melted[melted$bin %in% selection,], aes(x = bin, y = value, group=variable)) +
          geom_line(aes(color = variable), na.rm = TRUE) + theme_classic() +
          ylab(paste(name,"frequency", sep = " ")) + xlab("Genomic coordinates (Bins)") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          scale_x_discrete(breaks = x_breaks) +
          geom_vline(xintercept = x_breaks, color = vline) +
          labs(title = paste(level, cohort, name, sep = "-")) +
          theme(legend.position="top")
        cohort.p[[name]] <- p
        
      }
      final <- ggarrange(plotlist = cohort.p, common.legend = TRUE, nrow = 2)
      plots[[covThr]][[level]][[cohort]] <- final
    }
    
  }
  
}
#------
