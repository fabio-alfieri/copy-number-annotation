generate_signature_matrix <- function(rds_path, model_name = "Mid-length::Amplification model") {

  lis <- readRDS(rds_path)
  feature_ampl <- lis$models.X.test[[model_name]]
  shap_ampl <- lis$models.shap.df[[model_name]]
  
  feature_ampl <- feature_ampl[, !colnames(feature_ampl) %in% c("bin", "Type", "ampl_score", "del_score")]
  labels <- shap_ampl$labels; shap_ampl$BIAS <- NULL; shap_ampl$labels <- NULL
  
  major_drivers <- apply(shap_ampl[,-1], 1, function(x) colnames(shap_ampl)[which.max(abs(x))])
  
  df <- as.data.frame(table(major_drivers))
  p <- ggplot(df, aes(x = major_drivers, y = Freq)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Major Drivers", y = "Frequency", title = "Barplot of Major Drivers")
  print(p)
  
  major_driver_int <- as.integer(factor(major_drivers))
  feature_ampl$driven_by <- major_driver_int
  
  mat_to_scale <- feature_ampl[, !colnames(feature_ampl) %in% "driven_by"]
  mat_scaled <- as.data.frame(apply(mat_to_scale, 2, scale))
  mat_scaled$driven_by <- feature_ampl$driven_by
  
  signature_matrix <- mat_scaled %>%
    group_by(driven_by) %>%
    summarise(across(everything(), ~ median(.x, na.rm = TRUE)))
  
  signature_matrix[is.na(signature_matrix)] <- 0
  
  heatmap_mat <- as.matrix(signature_matrix[,-1])
  rownames(heatmap_mat) <- signature_matrix$driven_by
  ht <- ComplexHeatmap::Heatmap(heatmap_mat, name = "Median z-score", show_row_names = TRUE, show_column_names = TRUE)
  ComplexHeatmap::draw(ht)
  
  return(signature_matrix)
}
