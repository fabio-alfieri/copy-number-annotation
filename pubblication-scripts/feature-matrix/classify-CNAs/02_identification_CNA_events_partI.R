
# Identification of CNA events 

# Functions
# To find bins within the corresponding bin or intersecting with the bin
scna.int.bins <- function(coordinates.l){
  # Number of bases for states in the ith bin
  chr <- gsub("chr","",coordinates.l$chr)
  bin <- coordinates.l$bin
  s <- coordinates.l$s
  e <- coordinates.l$e
  bin_length <- e - s
  
  scna_w_bin_position <- c()
  # SCNAs within the bin
  within <- scna.m[scna.m$Chromosome == chr & (scna.m$Start >= s & scna.m$End <= e),]
  if(nrow(within) > 0){
    within$Length.in.bin <- within$End - within$Start
    scna_w_bin_position <- rbind(scna_w_bin_position,within)}
  # SCNAs cover the bin but longer
  cover <- scna.m[scna.m$Chromosome == chr & (scna.m$Start < s & scna.m$End > e),]
  if(nrow(cover) > 0){
    cover$Length.in.bin <- e - s # will be bin length
    scna_w_bin_position <- rbind(scna_w_bin_position,cover)}
  # SCNAs start on the right side of the bin but end inside the bin
  right_start <- scna.m[scna.m$Chromosome == chr & (scna.m$Start < s & scna.m$End <= e & scna.m$End > s),]
  if(nrow(right_start) > 0){
    right_start$Length.in.bin <- right_start$End - s
    scna_w_bin_position <- rbind(scna_w_bin_position,right_start)
  }
  # SCNAs start inside the bin but end outside
  left_end <- scna.m[scna.m$Chromosome == chr & (scna.m$Start >= s & scna.m$Start < e & scna.m$End > e),]
  if(nrow(left_end) > 0){
    left_end$Length.in.bin <- e - left_end$Start
    scna_w_bin_position <- rbind(scna_w_bin_position,left_end)
  }
  
  if(length(scna_w_bin_position) > 0){
    scna_w_bin_position$binlength <- bin_length
    scna_w_bin_position$bin <- paste(chr,bin,sep = "_")}
  
  return(scna_w_bin_position)
} 

# Calculate amplification and deletion score for each bin (fraction of the bin)
# Variables needed: segMean
patient_matrix <- function(df){ # df is bin mappend segments
  
  # Filter for the SCNA class
  df <- df[df$Final.cluster == segment_class,]
  
  if(nrow(df) > 0){
    bin.l <- df$binlength[1]
    bin <- df$bin[1]
    patients <- unique(df$Sample)
    patient.df <- c()
    for(patient in patients){
      df_ampl <- df[df$Sample == patient & df$Segment_Mean > segMean,] # Amplified segments
      total.ampl.l <- sum(df_ampl$Length.in.bin) # Total amplified length in a bin
      total.ampl.f <- total.ampl.l / bin.l
      
      df_del <- df[df$Sample == patient & df$Segment_Mean < -segMean,] # Deleted segments
      total.del.l <- sum(df_del$Length.in.bin) # Total deleted length in a bin
      total.del.f <- total.del.l / bin.l
      
      df_dip <- df[df$Sample == patient & abs(df$Segment_Mean) <= segMean,] # Diploid segments
      total.dip.l <- sum(df_dip$Length.in.bin) # Total deleted length in a bin
      total.dip.f <- total.dip.l / bin.l
      
      scores_patient <- c(bin, patient, total.ampl.f, #ampl, 
                          total.del.f, #del, 
                          total.dip.f #,dip
      )
      patient.df <- rbind(patient.df,scores_patient)
    }

    patient.df <- as.data.frame(patient.df)
    colnames(patient.df) <- c("bin","Sample","Ampl.fraction","Del.fraction","Dip.fraction")
    patient.df[,3:5] <- sapply(patient.df[,3:5], as.numeric)
    return(patient.df)
  }
}

Mapping_info <- list()

for(level in levels){
  # Bin, start and end locations
  coordinates <- do.call(rbind,chr_backbone_namesfixed[[level]])
  #Make a list for each bin to call function
  coord.list <- list()
  for(i in 1:nrow(coordinates)){
    bin <- paste(coordinates$chr[i],
                 coordinates$bin[i],
                 sep = "_")
    coord.list[[bin]] <- list("chr" = paste0("chr",coordinates$chr[i]),
                              "bin" = coordinates$bin[i],
                              "s" = as.numeric(coordinates$start_bin[i]),
                              "e" = as.numeric(coordinates$end_bin[i]))}
  
  for(tumor_type in tumor_types){
    
    scna.m <- Cluster.data[[tumor_type]]
    
    # PICK ONLY ONE OF THEM
    # If consider segment clusters based on length and location
    #segment_class_list <- unique(scna.m$Final.cluster)
    # If consider segment clusters based on length only
    scna.m$Final.cluster <- unlist(lapply(scna.m$Final.cluster, function(x) strsplit(x,"_")[[1]][1]))
    segment_class_list <- unique(scna.m$Final.cluster)
    
    
    # Find segments in each bin
    # Call the function to find scnas within/intersecting with the corresponding bin and to calculate their length 
    res <- mclapply(coord.list, scna.int.bins, mc.cores = 15)
    
    # Find bins for which no segment mapped
    no_segment_bins <- which(lengths(res) == 0)
    # Keep this info
    Mapping_info[[level]][[tumor_type]] <- names(res)[no_segment_bins] # Before selecting segment clusters
    
    res <- res[-no_segment_bins]
    
    for(segment_class in segment_class_list){
      
      # Patient-based matrix
      scores_final <- mclapply(res, patient_matrix, mc.cores = 20)
      
      scores_final <- as.data.frame(do.call(rbind,scores_final))
      
      # PICK ONLY ONE OF THEM
      # If consider segment clusters based on length and location
      # save(scores_final, 
      #      file = paste0("./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/CNA_fractions_per_sample_per_bin/Length_Location/",
      #                    tumor_type,"_",level,"_class_",segment_class,"_Armfraction_",armf,"_segMean_",segMean,".RData"))
      #If consider segment clusters based on length only
      save(scores_final, 
           file = paste0("./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/CNA_fractions_per_sample_per_bin/Length/",
                         tumor_type,"_",level,"_class_",segment_class,"_Armfraction_",armf,"_segMean_",segMean,".RData"))
      
      # Remove from memory and trigger garbage collection
      rm(scores_final)
      gc()
      
    }
    
    rm(res)
    gc()
  }
}