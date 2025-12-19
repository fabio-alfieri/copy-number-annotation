packages <- c("dplyr", "ggplot2", "ggiraph", "scales")

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

get_residuals <- function(prediction_df){

  observed_col <- grep(pattern = "*_score", x = colnames(prediction_df), value = T)
  is_ampl <- grepl("ampl", x = observed_col)
  
  prediction_df$observed <- prediction_df[[observed_col]]
  prediction_df[[observed_col]] <- NULL

  prediction_df <- prediction_df %>%
    mutate(
      residual = observed - prediction,
      )

  return(prediction_df)
}