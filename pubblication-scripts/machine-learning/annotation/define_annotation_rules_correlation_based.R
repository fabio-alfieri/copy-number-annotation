define_annotation_rules_correlation_based <- function(model_type, correlations_table){
  
  features_occurrence <<- correlations_table[correlations_table$feature_type == "Occurrence", "feature"]
  features_occurrence <<- c(features_occurrence, paste0("-", features_occurrence))
  
  features_selection <- correlations_table[correlations_table$feature_type == "Selection", c("feature", model_type)]
  features_selection <- features_selection[!is.na(features_selection[,model_type,]),]
  
  features_selection$what.to.paste <- ifelse(features_selection[,model_type] == "Positive", "", "-")
  features_selection$feature_name_selection <- paste0(features_selection$what.to.paste, features_selection$feature)
  
  features_selection$what.not.to.paste <- ifelse(features_selection[,model_type] == "Positive", "-", "")
  features_selection$feature_name_relaxed_selection <- paste0(features_selection$what.not.to.paste, features_selection$feature)
  
  features_positive_selection <<- features_selection[features_selection[,model_type] == "Positive", "feature_name_selection"]
  features_relaxed_positive_selection <<- features_selection[features_selection[,model_type] == "Positive", "feature_name_relaxed_selection"]
  
  features_negative_selection <<- features_selection[features_selection[,model_type] == "Negative", "feature_name_selection"]
  features_relaxed_negative_selection <<- features_selection[features_selection[,model_type] == "Negative", "feature_name_relaxed_selection"]
  
  features_selection <- features_selection$feature
  features_selection <<- c(features_selection, paste0("-", features_selection))
  
  
  num_of_top_features <<- 1 # define rule for top1 annotation
  
}
