define_annotation_rules <- function(model_type, correlations_table){
  
  if(model_type == 'ampl'){ # define ruling system for amplification
    
    # occurrence is every occurrence related feature regardless of the sign,
    # -dist.to.closest.OG --> not wanting to amplify an OG, means it's not worth it
    # dist.to.closest.TSG --> not wanting to amplify a TSG, means it's not worth it
    features_occurrence <<- c("mean.GC.content",        "-mean.GC.content",
                              "Length_Counts.E1",       "-Length_Counts.E1", 
                              "Length_Counts.E2",       "-Length_Counts.E2",
                              "Length_Counts.E3",       "-Length_Counts.E3",
                              "Length_Counts.E4",       "-Length_Counts.E4",
                              "Length_Counts.E5",       "-Length_Counts.E5",
                              "Length_Counts.E6",       "-Length_Counts.E6",
                              "Length_Counts.E7",       "-Length_Counts.E7",
                              "Length_Counts.E8",       "-Length_Counts.E8",
                              "Length_Counts.E9",       "-Length_Counts.E9",
                              "Length_Counts.E10",      "-Length_Counts.E10",
                              "Length_Counts.E11",      "-Length_Counts.E11",
                              "Length_Counts.E12",      "-Length_Counts.E12",
                              "Length_Counts.E13",      "-Length_Counts.E13",
                              "Length_Counts.E14",      "-Length_Counts.E14",
                              "Length_Counts.E15",      "-Length_Counts.E15",
                              "Length_Counts.E16",      "-Length_Counts.E16",
                              "Length_Counts.E17",      "-Length_Counts.E17",
                              "Length_Counts.E18",      "-Length_Counts.E18",
                              "Length_Counts.E19",      "-Length_Counts.E19",
                              "Length_Counts.E20",      "-Length_Counts.E20",
                              "Length_Counts.E21",      "-Length_Counts.E21",
                              "Length_Counts.E22",      "-Length_Counts.E22",
                              "Length_Counts.E23",      "-Length_Counts.E23",
                              "Length_Counts.E24",      "-Length_Counts.E24",
                              "Length_Counts.E25",      "-Length_Counts.E25",
                              "distance.to.centromere", "-distance.to.centromere",
                              "distance.to.telomere",   "-distance.to.telomere",
                              "dist.to.closest.FGS",    "-dist.to.closest.FGS"
    )
    
    # positive selection is every positive selection related feature, 
    features_positive_selection <<- c("dist.to.closest.OG",    
                                      "partners.trans", 
                                      "total_n_partners.trans", 
                                      "total_n_PPIs.trans",   
                                      "total_n_ohnologs.mmpaper_trans",
                                      "total_n_paralogs_trans", 
                                      "all.int.trans", 
                                      "genes.bin",              
                                      "mutations_norm",         
                                      "Ess.distance_pancancer")
    
    # negative selection is every negative selection related feature, 
    features_negative_selection <<- c("-dist.to.closest.TSG" 
                                      )
    
    features_relaxed_selection <<- c("dist.to.closest.TSG",
                                     "-dist.to.closest.OG",
                                     "-all.int.trans",
                                     "-partners.trans", 
                                     "-total_n_partners.trans", 
                                     "-total_n_PPIs.trans", 
                                     "-mutations_norm", 
                                     "-Ess.distance_pancancer",
                                     "-total_n_paralogs_trans",
                                     "-genes.bin",
                                     "-total_n_ohnologs.mmpaper_trans"
                                     )
    
    features_selection <<- c(features_positive_selection, features_negative_selection, features_relaxed_selection)
    
  } else if(model_type == 'del'){
    
    # occurrence is every occurrence related feature regardless of the sign,
    # dist.to.closest.OG --> wanting to amplify an OG, means it's not worth it
    # -dist.to.closest.TSG --> not wanting to amplify a TSG, means it's not worth it
    features_occurrence <<- c("mean.GC.content",        "-mean.GC.content",
                              "Length_Counts.E1",       "-Length_Counts.E1", 
                              "Length_Counts.E2",       "-Length_Counts.E2",
                              "Length_Counts.E3",       "-Length_Counts.E3",
                              "Length_Counts.E4",       "-Length_Counts.E4",
                              "Length_Counts.E5",       "-Length_Counts.E5",
                              "Length_Counts.E6",       "-Length_Counts.E6",
                              "Length_Counts.E7",       "-Length_Counts.E7",
                              "Length_Counts.E8",       "-Length_Counts.E8",
                              "Length_Counts.E9",       "-Length_Counts.E9",
                              "Length_Counts.E10",      "-Length_Counts.E10",
                              "Length_Counts.E11",      "-Length_Counts.E11",
                              "Length_Counts.E12",      "-Length_Counts.E12",
                              "Length_Counts.E13",      "-Length_Counts.E13",
                              "Length_Counts.E14",      "-Length_Counts.E14",
                              "Length_Counts.E15",      "-Length_Counts.E15",
                              "Length_Counts.E16",      "-Length_Counts.E16",
                              "Length_Counts.E17",      "-Length_Counts.E17",
                              "Length_Counts.E18",      "-Length_Counts.E18",
                              "Length_Counts.E19",      "-Length_Counts.E19",
                              "Length_Counts.E20",      "-Length_Counts.E20",
                              "Length_Counts.E21",      "-Length_Counts.E21",
                              "Length_Counts.E22",      "-Length_Counts.E22",
                              "Length_Counts.E23",      "-Length_Counts.E23",
                              "Length_Counts.E24",      "-Length_Counts.E24",
                              "Length_Counts.E25",      "-Length_Counts.E25",
                              "distance.to.centromere", "-distance.to.centromere",
                              "distance.to.telomere",   "-distance.to.telomere",
                              "dist.to.closest.FGS",    "-dist.to.closest.FGS"
    )
    
    # positive selection is every positive selection related feature, 
    features_positive_selection <<- c("dist.to.closest.TSG",
                                      "genes.bin",
                                      "all.int.trans",
                                      "partners.trans",
                                      "total_n_partners.trans", 
                                      "total_n_PPIs.trans",     
                                      "total_n_ohnologs.mmpaper_trans",
                                      "total_n_paralogs_trans"
                                      )
    
    # negative selection is every negative selection related feature, 
    features_negative_selection <<- c("-dist.to.closest.OG",    
                                      "-Ess.distance_pancancer",
                                      "-mutations_norm"
                                      )
    
    features_relaxed_selection <<- c("-dist.to.closest.TSG",
                                     "Ess.distance_pancancer",
                                     "dist.to.closest.OG",
                                     "-all.int.trans", 
                                     "-genes.bin", 
                                     "mutations_norm",
                                     "-partners.trans",
                                     "-total_n_partners.trans", 
                                     "-total_n_PPIs.trans",     
                                     "-total_n_ohnologs.mmpaper_trans",
                                     "-total_n_paralogs_trans"
                                     )
    
    
    features_selection <<- c(features_positive_selection, features_negative_selection)
    
  }
  
  num_of_top_features <<- 1 # define rule for top1 annotation
  
}
