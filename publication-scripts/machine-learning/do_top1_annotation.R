do_top1_annotation <- function(row){
  
  idxs <- order(unlist(abs(row)), decreasing = TRUE)[1:num_of_top_features] # idx of top features
  top_shap_vals <- row[idxs]
  top_features <- colnames(actual_shap_mat)[idxs];
  signs <- sign(top_shap_vals); signs <- gsub(pattern = "1", replacement = "", signs) # sign of top features
  
  top_features <- paste0(signs, top_features)
  
  return(c(top_features, top_shap_vals))
  
}
