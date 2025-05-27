rm(list=ls())
gc(full=T)

setwd('/Users/gabry/OneDrive/Desktop/shiny_app/') # Gab
# setwd('/Users/ieo5099/Desktop/copy-number-annotation/') # Fab

source('dev/0_LoadData.R')
source('dev/1_dynamic_plotting_functions.R')
source('dev/1bis_parse_input_data_FA.R') # overwrite parse_input_data()

processed_data <- parse_input_data(shap.list = shap.list,
                                   out_annot_list = out_annot_list,
                                   chr_backbone_namesfixed = chr_backbone_namesfixed, 
                                   centromere_table = centromere_table,
                                   clustering_depth = 4)

# these files will be saved after

shap.list <- processed_data$shap.list
out_annot_list_processed <- processed_data$out_annot_list_processed
backbone.100kb <- processed_data$backbone.100kb
centromere_table <- processed_data$centromere_table

# this will be replaced by user input

plot_shap <- TRUE
plot_landscape <- TRUE

type_input <- "BRCA"
model_input <- parse_input_model("ampl")
coord_input <- NULL
chr_input <- NULL

filtered_shap_output_ampl <- filter_df(input_obj = shap.list, backbone_granges = backbone.100kb, 
                                       type_input = type_input, model_input = "ampl", 
                                       chr_input = chr_input, coord_input = coord_input)

filtered_shap_output_del <- filter_df(input_obj = shap.list, backbone_granges = backbone.100kb, 
                                      type_input = type_input, model_input = "del", 
                                      chr_input = chr_input, coord_input = coord_input)


filtered_landscape_output_ampl <- filter_df(input_obj = out_annot_list_processed, backbone_granges = backbone.100kb,
                                       type_input = type_input, model_input = "ampl", chr_input = chr_input, coord_input = coord_input)

filtered_landscape_output_del <- filter_df(input_obj = out_annot_list_processed, backbone_granges = backbone.100kb,
                                        type_input = type_input, model_input = "del", chr_input = chr_input, coord_input = coord_input)


filtered_shap_ampl <- filtered_shap_output_ampl$final_df; filtered_shap_del <- filtered_shap_output_del$final_df
filtered_landscape_ampl <- filtered_landscape_output_ampl$final_df; filtered_landscape_del <- filtered_landscape_output_del$final_df

shap_plotting_list <- prepare_shap_to_plot(filtered_shap_ampl = filtered_shap_ampl, 
                                           filtered_shap_del = filtered_shap_del)

if (model_input$selected == "ampl") {
  
  filtered_shap_abs_sum <- shap_plotting_list$filtered_shap_abs_sum_ampl
  genome_mask <- filtered_shap_output_ampl$genome_mask
  model_mask <- filtered_shap_output_ampl$model_mask
  type_mask <- filtered_shap_output_ampl$type_mask
  filtered_landscape <- filtered_landscape_ampl
  
} else {
  
  filtered_shap_abs_sum <- shap_plotting_list$filtered_shap_abs_sum_del
  genome_mask <- filtered_shap_output_del$genome_mask
  model_mask <- filtered_shap_output_del$model_mask
  type_mask <- filtered_shap_output_del$type_mask
  filtered_landscape <- filtered_landscape_del
  
}

if (plot_shap) {
  barplot_shap(shap.abs.sum = filtered_shap_abs_sum, 
             genome_mask = genome_mask, 
             type_mask = type_mask, 
             model_mask = model_mask)
}

if (plot_landscape) {
  
p <- landscape_plot_interactive(filtered_landscape = filtered_landscape, 
                           genome_mask = genome_mask, 
                           type_mask = type_mask,
                           plot_ampl = TRUE, 
                           plot_del = TRUE, 
                           annot_to_plot_ticks = "all", 
                           annot_to_plot_kde = "all",
                           backbone.100kb = backbone.100kb)


tooltip_css <- "
    [class^='tooltip_svg_'] {
      background: transparent !important;
      color: #fafafa;
      padding: 2px 6px;
      border-radius: 2px;
      font-family: 'Roboto', sans-serif;
      font-size: 12px;
      line-height: 1.2;
      box-shadow: 0 0px 0px rgba(0,0,0,0);
      opacity: 0;
      transform: translateY(2px);
      transition: opacity 0.15s ease, transform 0.15s ease;
      pointer-events: auto;
      white-space: nowrap;
      border: none !important;
      z-index: 1000;
    }
    [class^='tooltip_svg_']:hover {
      opacity: 0.9;
      transform: translateY(0);
    }
    [class^='tooltip_svg_']::after { display: none !important; }
    [class^='tooltip_svg_'] strong { font-weight: bold; }
    @media (max-width: 600px) {
      [class^='tooltip_svg_'] {
        font-size: 11px;
        padding: 2px 4px;
      }
    }
  "

p2 <- htmlwidgets::prependContent(p, tags$head(tags$style(HTML(tooltip_css))))
p2

# htmlwidgets::saveWidget(p2, 
#                         paste0("dev/Plots/landscape_", type_mask, ".html"), 
#                         selfcontained = TRUE)

}




