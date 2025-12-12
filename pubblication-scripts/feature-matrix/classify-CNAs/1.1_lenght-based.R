
# Classification of segments - Length based

Output.List <- list() # Output file

for(tumor_type in cohorts){
  file <- paste0("./Data/Mapping_SCNAs/temp4/",tumor_type,"_armfraction_scna_segmentmean_nofilter.RData")
  # Load scna_w_armfraction for the specific tumor_type
  load(file)
  scna_w_armfraction$ArmII <- str_replace_all(scna_w_armfraction$Arm, "[:digit:]", "")
  
  # Based on length
  scna_w_armfraction$Lengthclass <- ifelse(scna_w_armfraction$length <= 3000000 & scna_w_armfraction$length >= 1000, "Focal",
                                           ifelse(scna_w_armfraction$length > 3000000, "Large","Not classified"))
  
  # within p-arm
  within.p <- scna_w_armfraction[scna_w_armfraction$ArmII == "p",]
  within.p$Class <- ifelse(within.p$ArmFraction.p >= armf, "Arm-level",
                           ifelse(within.p$ArmFraction.p < armf & within.p$ArmFraction.p >= 0.20, "Mid-length", "Small-scale"))
  # within q-arm
  within.q <- scna_w_armfraction[scna_w_armfraction$ArmII == "q",]
  within.q$Class <- ifelse(within.q$ArmFraction.q >= armf, "Arm-level",
                           ifelse(within.q$ArmFraction.q < armf & within.q$ArmFraction.q >= 0.20, "Mid-length", "Small-scale"))
  
  # Both arms
  rest <- scna_w_armfraction[scna_w_armfraction$Arm == "p-q",]
  rest$ArmFraction <- ifelse(rest$ArmFraction.p > rest$ArmFraction.q, rest$ArmFraction.p, rest$ArmFraction.q)
  rest$Class <- ifelse(rest$ArmFraction >= armf, "Arm-level",
                       ifelse(rest$ArmFraction < armf & rest$ArmFraction >= 0.20, "Mid-length", "Small-scale"))
  # Merge all data
  scna_w_armfraction.new <- rbind(within.p,within.q)
  scna_w_armfraction.new <- rbind(scna_w_armfraction.new, rest[,colnames(rest) != "ArmFraction"])
  
  # Chromosome-level segments
  # Grouping chromosomes
  m.acro <- scna_w_armfraction.new[scna_w_armfraction.new$Chromosome %in% c(13,14,15,21,22),]
  m.meta <- scna_w_armfraction.new[!scna_w_armfraction.new$Chromosome %in% c(13,14,15,21,22),]
  
  #For (sub)metacentric chromosomes
  m.meta$ClassII <- ifelse(m.meta$ArmFraction.p >= armf & m.meta$ArmFraction.q >= armf, "Chromosome-level", m.meta$Class)
  
  #For acrocentric chromosomes
  m.acro$ClassII <- ifelse(m.acro$ArmFraction.q >= armf, "Chromosome-level",  m.acro$Class)
  
  scna_w_armfraction.new <- rbind(m.acro,m.meta)
  scna_w_armfraction.new$CLUSTER1 <- ifelse(scna_w_armfraction.new$ClassII == "Chromosome-level", "Chromosome-level",
                                            ifelse(scna_w_armfraction.new$ClassII == "Arm-level", "Arm-level",
                                                   ifelse(scna_w_armfraction.new$ClassII == "Small-scale" & scna_w_armfraction.new$Lengthclass == "Focal", "Small-scale",
                                                          ifelse(scna_w_armfraction.new$Lengthclass == "Not classified","Not classified","Mid-length"))))
  scna_w_armfraction.new <- scna_w_armfraction.new[,c(1:4,6,13:16,20)]

  Output.List[[tumor_type]] <- scna_w_armfraction.new
}

save(Output.List,
     file = paste0("./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Data/ClusterI","armfraction_",as.character(armf),".RData"))