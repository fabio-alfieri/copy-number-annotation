# app.R

library(shiny)
library(ggiraph)
library(htmlwidgets)
library(ggplot2)
library(dplyr)
library(shinycssloaders)

setwd("./")

# Load helper functions
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
shap.list                <- data_processed$shap.list
out_annot_list_processed <- data_processed$out_annot_list_processed
backbone.100kb           <- data_processed$backbone.100kb
centromere_table         <- data_processed$centromere_table

# CSS styles
tooltip_css <- "
  background: transparent !important;
  color: #fafafa;
  padding: 2px 6px;
  border-radius: 2px;
  font-family: 'Roboto', sans-serif;
  font-size: 12px;
  line-height: 1.2;
  box-shadow: none !important;
  transform: translateY(2px);
  white-space: nowrap;
"
hover_css <- "opacity: 0.9 !important; transform: translateY(0) !important;"

# UI
ui <- fluidPage(
  titlePanel("SCNA Landscape Plot"),
  sidebarLayout(
    sidebarPanel(
      selectInput("type_input", "Select cancer type:",
                  choices = c("STAD","GBMLGG","COADREAD","KIRP","KIRC",
                              "OV","ESCA","LUAD","LUSC","PAAD","BRCA"),
                  selected = "BRCA"),
      
      selectInput("chr_input", "Select chromosome:",
                  choices = c("Whole Genome" = "whole_genome", paste0("chr", 1:22)),
                  selected = "whole_genome"),
      
      textInput("genomic_coords", 
                "Genomic coordinates (chr:start-end):", 
                ""),  # Must be "" in UI, will convert to NULL in server
      
      checkboxInput("plot_ampl", "Show amplification track", TRUE),
      checkboxInput("plot_del",  "Show deletion track",    TRUE),
      
      checkboxInput("enable_ticks", "Enable annotation ticks", TRUE),
      conditionalPanel(
        condition = "input.enable_ticks",
        textInput("annot_ticks_input", 
                  "Ticks clusters ('all' or numbers comma-separated):", 
                  "all")
      ),
      
      checkboxInput("enable_kde", "Enable KDE layers", TRUE),
      conditionalPanel(
        condition = "input.enable_kde",
        textInput("annot_kde_input", 
                  "KDE clusters ('all' or numbers comma-separated):", 
                  "all")
      ),
      
      actionButton("go", "Go!", class = "btn-primary")
    ),
    mainPanel(
      withSpinner(
        girafeOutput("landscape_plot", width = "100%", height = "600px"),
        type = 6
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  plot_reactive <- eventReactive(input$go, {
    withProgress(message = "Building SCNA landscape…", value = 0, {
      incProgress(0.1, detail = "Filtering data")
      
      selected_chr <- if (input$chr_input == "whole_genome") NULL else input$chr_input
      input_coord  <- if (input$genomic_coords == "") NULL else input$genomic_coords
      
      shap_res <- filter_df(
        shap.list,
        backbone.100kb,
        input$type_input,
        "ampl",
        selected_chr,
        input_coord
      )
      
      land_res <- filter_df(
        out_annot_list_processed,
        backbone.100kb,
        input$type_input,
        "ampl",
        selected_chr,
        input_coord
      )
      
      df_land <- land_res$final_df
      req(nrow(df_land) > 0)
      
      ticks_val <- if (!input$enable_ticks) FALSE else {
        txt <- tolower(input$annot_ticks_input)
        if (txt == "all") "all" else as.numeric(strsplit(txt, ",")[[1]])
      }
      
      kde_val <- if (!input$enable_kde) FALSE else {
        txt <- tolower(input$annot_kde_input)
        if (txt == "all") "all" else as.numeric(strsplit(txt, ",")[[1]])
      }
      
      incProgress(0.6, detail = "Generating plot")
      gir_widget <- landscape_plot_interactive(
        filtered_landscape  = df_land,
        backbone.100kb      = backbone.100kb,
        genome_mask         = shap_res$genome_mask,
        type_mask           = shap_res$type_mask,
        plot_ampl           = input$plot_ampl,
        plot_del            = input$plot_del,
        annot_to_plot_ticks = ticks_val,
        annot_to_plot_kde   = kde_val
      )
      
      incProgress(0.2, detail = "Configuring interactivity")
      gir_widget <- girafe_options(
        gir_widget,
        opts_tooltip(
          css             = tooltip_css,
          delay_mouseover = 0,
          delay_mouseout  = 0,
          offx            = 10,
          offy            = -10
        ),
        opts_hover(css = hover_css)
      )
      
      incProgress(1)
      showNotification("Landscape ready!", type = "message", duration = 2)
      gir_widget
    })
  }, ignoreNULL = FALSE)
  
  output$landscape_plot <- renderGirafe({
    plot_reactive()
  })
}

shinyApp(ui, server)
