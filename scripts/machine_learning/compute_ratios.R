compute_selection_occurrence_ratio <- function(row){ # iterate over rows
  
  signs <- sign(row); signs <- gsub(pattern = "1", replacement = "", signs) # take signs as strings
  sign_aware_features <- paste0(signs, features) # annotate features sign-aware
  names(row) <- sign_aware_features
  
  sel <- sum(abs(row[names(row) %in% features_selection]))+0.01  # selection weight 
  occ <- sum(abs(row[names(row) %in% features_occurrence]))+0.01 # occurrence weight
  
  ratio <- sel / occ # ratio
  
  return(ratio)
} # generalize

compute_positive_negative_ratio <- function(row){ # iterate over rows
  
  signs <- sign(row); signs <- gsub(pattern = "1", replacement = "", signs) # take signs as strings
  sign_aware_features <- paste0(signs, features) # annotate features sign-aware
  names(row) <- sign_aware_features
  
  pos_sel <- sum(abs(row[names(row) %in% features_positive_selection]))+0.01 # positive selection weight
  neg_sel <- sum(abs(row[names(row) %in% features_negative_selection]))+0.01 # negative selection weight
  
  ratio <- pos_sel / neg_sel # ratio
  
  return(ratio)
}

compute_ratio <- function(row, top_features, bottom_features){
  
  signs <- sign(row); signs <- gsub(pattern = "1", replacement = "", signs) # take signs as strings
  sign_aware_features <- paste0(signs, features) # annotate features sign-aware
  names(row) <- sign_aware_features
  
  top <- sum(abs(row[names(row) %in% top_features]))+0.01
  bottom <- sum(abs(row[names(row) %in% bottom_features]))+0.01
  
  ratio <- top / bottom # ratio
  
  return(ratio)
  
}



