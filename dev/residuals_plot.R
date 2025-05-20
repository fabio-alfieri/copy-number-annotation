library(ggplot2)
library(ggiraph)
library(dplyr)
library(scales)

plot_residuals <- function(prediction_df){

  get_contrast <- function(hexcol) {
    rgb <- col2rgb(hexcol) / 255
    lum <- 0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]
    ifelse(lum > 0.5, "#000000", "#FFFFFF")
}
  
  observed_col <- grep(pattern = "*_score", x = colnames(prediction_df), value = T)
  is_ampl <- grepl("ampl", x = observed_col)
  
  prediction_df$observed <- prediction_df[[observed_col]]
  prediction_df[[observed_col]] <- NULL
  
  palette_low <- "beige"
  palette_high <- ifelse(is_ampl, "red", "blue")
  
  gradient_fn <- col_numeric(palette = c(palette_low, 
                                         palette_high), 
                             domain = NULL)
  
  prediction_df <- prediction_df %>%
    mutate(
      uid         = row_number(),
      ideal       = observed,
      residual    = abs(prediction - ideal),
      fill_color  = gradient_fn(residual),
      text_color  = get_contrast(fill_color),
      alpha_val   = rescale(residual, to = c(0.2, 1)),
      tooltip_html = paste0(
        "<div style='padding:5px; background-color:", fill_color, ";",
        "color:", text_color, "; border-radius:4px;'>",
        "<strong>Prediction:</strong> ", round(prediction, 3), "<br/>",
        "<strong>Ideal:</strong> ", round(ideal, 3), "<br/>",
        "<strong>Residual:</strong> ", round(residual, 3),
        "</div>"
      )
    )
  
  p <- ggplot(prediction_df) +
    geom_line(aes(x = observed, y = ideal),
              color = "black", size = 0.8) +
    geom_point_interactive(
      aes(x = observed, y = prediction,
          tooltip = tooltip_html,
          data_id = uid,
          fill = residual,
          alpha = alpha_val),
      shape     = 21,
      color     = "black",
      size      = 3.5,
      stroke    = 0.3,
      hover_css = paste(
        "fill-opacity:1;",
        "stroke-opacity:1;",
        "stroke-width:1px;",
        "r:6px;"
      )
    ) +
    geom_segment_interactive(
      aes(x = observed, xend = observed,
          y = ideal, yend = prediction,
          data_id = uid),
      color     = "black",
      size      = 0.3,
      alpha     = 0,
      hover_css = paste(
        "stroke-opacity:1;",
        "stroke-width:0.3px;"
      )
    ) +
    scale_fill_gradient(
      low = palette_low, high = palette_high,
      name = "Absolute Residual"
    ) +
    guides(alpha = "none") +
    labs(title = "Interactive Residual Plot with Gradient Fill & Transparency",
         x = paste0("Observed", ifelse(is_ampl, "Amplification", "Deletion"), "Score"),
         y = paste0("Predicted", ifelse(is_ampl, "Amplification", "Deletion"), "Score")) +
    theme_minimal()
  
  girafe(
    ggobj = p,
    options = list(
      opts_hover(css     = "cursor:pointer;"),
      opts_hover_inv(css = "opacity:0.05;"),
      opts_tooltip(css   = "box-shadow:2px 2px 6px rgba(0,0,0,0.3); font-size:12px;")
    )
  )
}





  