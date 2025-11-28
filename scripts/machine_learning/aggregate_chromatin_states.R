aggregate_chromatin_states <- function(df, col.name) {
  
  promoters <- c("Length_Counts.E1", "Length_Counts.E2", "Length_Counts.E3", 
                 "Length_Counts.E4", "Length_Counts.E22", "Length_Counts.E23")
  
  transcribed <- c("Length_Counts.E5", "Length_Counts.E6", "Length_Counts.E7", 
                   "Length_Counts.E8", "Length_Counts.E9", "Length_Counts.E10",
                   "Length_Counts.E11", "Length_Counts.E12")
  
  enhancers <- c("Length_Counts.E13", "Length_Counts.E14", "Length_Counts.E15", 
                 "Length_Counts.E16", "Length_Counts.E17", "Length_Counts.E18")
  
  repressed <- c("Length_Counts.E19", "Length_Counts.E20", "Length_Counts.E24", 
                 "Length_Counts.E25")
  
  accessible <- c("Length_Counts.E21")
  
  df[df[,col.name] %in% promoters,][,col.name] <- "promoters"
  df[df[,col.name] %in% transcribed,][,col.name] <- "transcribed"
  df[df[,col.name] %in% enhancers,][,col.name] <- "enhancers"
  df[df[,col.name] %in% repressed,][,col.name] <- "repressed"
  df[df[,col.name] %in% accessible,][,col.name] <- "accessible"
  
  return(df)
}
