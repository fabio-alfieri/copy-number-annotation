
# Classification of segments - Location based

outpath <- paste0("./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/ClusterII/","armfraction_",as.character(armf),".RData")

Output.List.new <- list()

for(cohort in cohorts){
  # Centromere/telomere distances
  m <- Data[[cohort]] 
  m <- m[,c(2:5,16,17)]
  
  # Data from the first classification
  scna <- Output.List[[cohort]]
  
  # Merge two data
  scna <- merge(scna, m, by = c("Sample","Chromosome","Start","End"))
  scna.new <- c()
  for(chr in unique(scna$Chromosome)){
    cent.mean <- as.numeric(centromere[centromere$Chromosome == as.character(chr),"Mean"])
    tel.mean <- as.numeric(telomere[telomere$Chromosome == as.character(chr),"Mean"])
    
    scna.chr <- scna[scna$Chromosome == chr,]
    scna.chr$Centromere.class <- ifelse(scna.chr$Dist.centromere <= cent.mean, "Centromere-bounded","Interstitial")
    scna.chr$Telomere.class <- ifelse(scna.chr$Dist.closest.telomere <= tel.mean, "Telomere-bounded","Interstitial")
    scna.new <- rbind(scna.new,scna.chr)
  }
  Output.List.new[[cohort]] <- scna.new
}

save(Output.List.new,
     file = outpath)
