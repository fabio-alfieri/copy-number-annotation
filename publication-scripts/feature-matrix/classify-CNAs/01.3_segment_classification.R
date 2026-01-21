
final_clustering_outpath <- paste0("./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/Final_clustering/","armfraction_",armf,".RData")

Cluster.data <- list()
Info <- c()
for(cohort in names(Output.List.new)){
  scna <- Output.List.new[[cohort]]
  # If a segment spans on both arms, change "Centromere-bounded" to Centromere-spanning"
  scna$Centromere.class <- ifelse(scna$ArmII == "p-q", "Centromere-spanning",scna$Centromere.class)
  # Remove segments under 1kb
  scna <- scna[scna$CLUSTER1 != "Not classified",]
  #info <- scna %>% group_by(Centromere.class,Telomere.class,CLUSTER1) %>% summarise(Freq = n())
  scna.final <- c()
  for(cl1 in c("Small-scale","Mid-length","Arm-level","Chromosome-level")){
    if(cl1 %in% c("Small-scale","Mid-length")){
      cluster.m <- scna[scna$CLUSTER1 == cl1,]
      cluster.m$CLUSTER2 <- ifelse(cluster.m$Centromere.class == "Centromere-bounded" & cluster.m$Telomere.class == "Interstitial","Centromere-bounded",
                                   ifelse(cluster.m$Centromere.class == "Centromere-spanning" & cluster.m$Telomere.class == "Interstitial","Centromere-spanning",
                                          ifelse(cluster.m$Centromere.class == "Centromere-bounded" & cluster.m$Telomere.class == "Telomere-bounded","Arm-level",
                                                 ifelse(cluster.m$Centromere.class == "Interstitial" & cluster.m$Telomere.class == "Telomere-bounded","Telomere-bounded","Interstitial"))))
      
      cluster.m$Final.cluster <- ifelse(cluster.m$CLUSTER2 == "Arm-level" & cluster.m$Chromosome %in% c("13","14","15","21","22"), "Chromosome-level", 
                                        ifelse(cluster.m$CLUSTER2 == "Arm-level" & !cluster.m$Chromosome %in% c("13","14","15","21","22"), "Arm-level", paste(cluster.m$CLUSTER1,cluster.m$CLUSTER2,sep = "_")))
      scna.final <- rbind(scna.final,cluster.m[,colnames(cluster.m) != "CLUSTER2"])
    }
    else{
      cluster.m <- scna[scna$CLUSTER1 == cl1,]
      cluster.m$Final.cluster <- cl1
      scna.final <- rbind(scna.final,cluster.m)
    }
  }
  Cluster.data[[cohort]] <- scna.final
  m <- as.data.frame(table(scna.final$Final.cluster))
  colnames(m) <- c("Cluster","Size")
  m$Cohort <- cohort
  Info <- rbind(Info,m)
}
save(Cluster.data,Info,
     file = final_clustering_outpath)


# To visualize cluster in IGV (Needs to be modified)
# for(cluster1 in unique(data.new$CLUSTER1)){
#   for(cluster2 in unique(data$CLUSTER2)){
#     m <- data.new[data.new$CLUSTER1 == cluster1 & data.new$CLUSTER2 == cluster2,c(2,3,4,5,7)]
#     colnames(m) <- c("ID","chrom","loc.start","loc.end","seg.mean")
#     write.table(m, 
#                 file = paste0("./IGV/Cluster_","_",cohort,"_",paste(cluster1,cluster2,sep = "_"),".seg"), 
#                 quote = F, row.names = F,col.names = T,sep = "\t")}
# }
