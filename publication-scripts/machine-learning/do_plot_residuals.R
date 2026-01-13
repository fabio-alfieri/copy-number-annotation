wd <- 'path/to/GitHub/copy-number-annotation/'

do_plot_residuals <- function(input_df, save = TRUE){
  
df <- input_df %>%
  mutate(prediction_status = case_when(
    pred_is_correct == TRUE ~ "Correct Prediction",
    pred_is_correct == FALSE & residual < 0 ~ "Incorrect Prediction - Negative Selection",
    pred_is_correct == FALSE & residual > 0 ~ "Incorrect Prediction - Positive Selection"
  ))

df$type <- unlist(lapply(X = df$labels, FUN = function(x){strsplit(x, split = "-")[[1]][2]}))

pancancer_residuals <- ggplot(df, aes(x = residual, fill = prediction_status)) +
  geom_histogram(bins = 100, color = "black", position = "identity", alpha = 0.7) +
  scale_fill_manual(
    values = c(
      "Correct Prediction" = "steelblue",
      "Incorrect Prediction - Negative Selection" = "orange",
      "Incorrect Prediction - Positive Selection" = "red"
    )
  ) +
  theme_minimal() +
  labs(
    title = "Residuals (Pan-Cancer)",
    x = "Residual",
    y = "Count",
    fill = "Prediction Status"
  )

tumor_specific_residuals <- ggplot(df, aes(x = residual, fill = prediction_status)) +
  geom_histogram(bins = 50, color = "black", position = "identity", alpha = 0.7) +
  scale_fill_manual(
    values = c(
      "Correct Prediction" = "steelblue",
      "Incorrect Prediction - Negative Selection" = "orange",
      "Incorrect Prediction - Positive Selection" = "red"
    )
  ) +
  facet_wrap(~ type, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Residuals by Cancer Type",
    x = "Residual",
    y = "Count",
    fill = "Prediction Status"
  )

if (save) {
  
  ggsave(filename = paste0(wd, "data/plots/", i, "_pancancer_residuals.pdf"), plot = pancancer_residuals, width = 15, height = 14, dpi = 300)
  ggsave(filename = paste0(wd, "data/plots/", i, "_tumor_specific_residuals.pdf"), plot = tumor_specific_residuals , width = 15, height = 14, dpi = 300)
  
}

outlist <- list(pancancer_residuals = pancancer_residuals,
                tumor_specific_residuals = tumor_specific_residuals)

return(outlist)

}
