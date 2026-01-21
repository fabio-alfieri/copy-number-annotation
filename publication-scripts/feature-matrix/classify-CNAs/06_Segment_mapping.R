# Functions
# To find bins within the corresponding bin or intersecting with the bin
scna.int.bins <- function(coordinates.l){
  # Number of bases for states in the ith bin
  chr <- as.numeric(gsub("chr","",coordinates.l$chr))
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

  # Find segments in each bin
  # Call the function to find scnas within/intersecting with the corresponding bin and to calculate their length 
  res <- mclapply(coord.list, scna.int.bins, mc.cores = 15)
  outpath <- paste0("./Codes/Codes-CNAs/MethodII/Parameter_tuning_segments/Segments_mapped_to_bins/", tumor_type,"_",level,".RData")
  save(res,
       file = outpath)
      
  # Remove from memory and trigger garbage collection
  rm(res)
  gc()
}
