library(shiny)
library(ggiraph)
library(htmlwidgets)
library(ggplot2)
library(dplyr)

# Source data loading and helper functions
source('dev/0_LoadData.R')
source('dev/1_dynamic_plotting_functions.R')
source('dev/1bis_parse_input_data_FA.R')

# Preprocess data once at startup
data_processed <- parse_input_data(
  shap.list = shap.list,
  out_annot_list = out_annot_list,
  chr_backbone_namesfixed = chr_backbone_namesfixed,
  centromere_table = centromere_table,
  clustering_depth = 4
)
# Extract processed objects
shap.list <- data_processed$shap.list
out_annot_list_processed <- data_processed$out_annot_list_processed
backbone.100kb <- data_processed$backbone.100kb
centromere_table <- data_processed$centromere_table

# Define CSS for tooltips
tooltip_css <- "
  background: transparent !important;
  color: #fafafa;
  padding: 2px 6px;
  border-radius: 2px;
  font-family: 'Roboto', sans-serif;
  font-size: 12px;
  line-height: 1.2;
  box-shadow: 0 0px 0px rgba(0,0,0,0);
  transform: translateY(2px);
  transition: opacity 0.15s ease, transform 0.15s ease;
  white-space: nowrap;
  border: none !important;
  z-index: 1000;
"
# Define CSS for hover state
hover_css <- "
  opacity: 0.9 !important;
  transform: translateY(0) !important;
"

# Define UI: allow selection of cancer type
ui <- fluidPage(
  titlePanel("SCNA Landscape Plot"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "type_input",
        label = "Select cancer type:",
        choices = c("STAD","GBMLGG","COADREAD","KIRP","KIRC",
                    "OV","ESCA","LUAD","LUSC","PAAD","BRCA"),
        selected = "BRCA"
      )
    ),
    mainPanel(
      girafeOutput("landscape_plot", width = "100%", height = "600px")
    )
  )
)

# Server logic: dynamic type input
server <- function(input, output, session) {
  output$landscape_plot <- renderGirafe({
    # Filter based on user-selected type and ampl
    shap_res <- filter_df(
      input_obj = shap.list,
      backbone_granges = backbone.100kb,
      type_input = input$type_input,
      model_input = "ampl",
      chr_input = NULL,
      coord_input = NULL
    )
    land_res <- filter_df(
      input_obj = out_annot_list_processed,
      backbone_granges = backbone.100kb,
      type_input = input$type_input,
      model_input = "ampl",
      chr_input = NULL,
      coord_input = NULL
    )
    df_land <- land_res$final_df
    req(nrow(df_land) > 0)
    
    # Generate interactive ggiraph object
    girafe_obj <- landscape_plot_interactive(
      filtered_landscape = df_land,
      backbone.100kb      = backbone.100kb,
      genome_mask         = shap_res$genome_mask,
      type_mask           = shap_res$type_mask,
      plot_ampl           = FALSE,
      plot_del            = FALSE,
      annot_to_plot_ticks = FALSE,
      annot_to_plot_kde   = "all"
    )
    # Apply tooltip and hover CSS without pointer-events
    girafe_options(
      girafe_obj,
      opts_tooltip(
        css = tooltip_css,
        delay_mouseover = 0,
        delay_mouseout = 0,
        offx = 10,
        offy = -10
      ),
      opts_hover(css = hover_css)
    )
  })
}

# Launch app
shinyApp(ui, server)
