library(ggplot2)
library(ggiraph)
library(dplyr)
library(scales)

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