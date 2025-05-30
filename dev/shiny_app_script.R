library(shiny)
library(ggiraph)
library(htmlwidgets)
library(ggplot2)
library(dplyr)
library(shinycssloaders)

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
tooltip_css <- "[class^='tooltip_svg_'] {
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
}"
hover_css <- "opacity: 0.9 !important; transform: translateY(0) !important;"

# UI
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      html, body {
        height: 100%;
        margin: 0;
        padding: 0;
        overflow: hidden;
        font-family: 'Segoe UI', sans-serif;
      }

      .container-fluid {
        padding: 0 !important;
        margin: 0 !important;
        height: 100vh;
        display: flex;
        flex-direction: column;
      }

      .title-panel {
        flex: 0 0 auto;
        padding: 10px 0;
        font-size: 22px;
        font-weight: bold;
        background-color: #f0f4fa;
        border-bottom: 1px solid #dbe9ff;
        text-align: center;
        width: 100vw;
        margin: 0;
        box-sizing: border-box;
      }

      #flex-container {
        display: flex;
        flex-grow: 1;
        overflow: hidden;
      }

      #sidebar-container {
        width: 300px;
        overflow-y: auto;
        flex-shrink: 0;
        transition: width 0.3s ease;
        display: flex;
        flex-direction: column;
        background-color: #ffffff;
        border-right: 1px solid #ffffff;
        padding: 10px;
        align-items: center;
      }

      #sidebar-container.collapsed {
        width: 50px !important;
        padding: 10px 5px;
        display: flex !important;
        flex-direction: column;
        align-items: center;
        justify-content: flex-start;
      }

      #sidebar-container.collapsed .well,
      #sidebar-container.collapsed .form-group,
      #sidebar-container.collapsed .checkbox,
      #sidebar-container.collapsed .shiny-input-container {
        display: none !important;
      }

      #sidebar-container.collapsed .btn {
        width: 100%;
        margin: 0 0 10px 0;
        font-size: 18px;
        padding: 6px 0;
        text-align: center;
      }

      #sidebar-container .btn {
        width: auto;
        margin: 0 auto;
      }

      #main-panel-container {
        flex-grow: 1;
        padding-left: 20px;
        overflow: hidden;
        display: flex;
        flex-direction: column;
      }

      #main-panel-container.expanded {
        margin-left: 0;
      }

      #landscape_plot {
        flex-grow: 1;
        min-height: 0;
      }

      .form-control, .selectize-input {
        border-radius: 6px;
        border-color: #cfe2ff;
      }

      .btn-primary {
        background-color: #2a78e8;
        border-color: #2a78e8;
      }

      .btn-primary:hover {
        background-color: #155ab3;
        border-color: #155ab3;
      }

      .checkbox input[type='checkbox'] {
        accent-color: #2a78e8;
      }

      .shiny-input-container label {
        font-weight: 600;
        font-size: 14px;
        margin-bottom: 4px;
      }
    ")),
    tags$script(HTML(
      "$(document).on('shiny:connected', function(){
         $('#toggle_sidebar').on('click', function(){
           $('#sidebar-container').toggleClass('collapsed');
           $('#main-panel-container').toggleClass('expanded');
         });
       });"
    ))
  ),
  
  div(class = "title-panel", "SCNA Landscape Plot"),
  
  div(id = "flex-container",
      div(id = "sidebar-container", class = "sidebar-panel",
          actionButton("toggle_sidebar", "☰", class = "btn btn-secondary", style = "margin-bottom:10px; width:100%; font-size:18px;"),
          sidebarPanel(width = 12,
                       selectInput("type_input", "Select cancer type:",
                                   choices = c("STAD","GBMLGG","COADREAD","KIRP","KIRC","OV","ESCA","LUAD","LUSC","PAAD","BRCA"),
                                   selected = "BRCA"),
                       selectInput("model_input", "Select model:",
                                   choices = c("Amplification" = "ampl", "Deletion" = "del"),
                                   selected = "ampl"),
                       selectInput("chr_input", "Select chromosome(s):",
                                   choices = paste0("chr", 1:22),
                                   selected = paste0("chr", 1:22),
                                   multiple = TRUE,
                                   selectize = TRUE),
                       textInput("genomic_coords", "Genomic coordinates (chr:start-end):", ""),
                       checkboxInput("plot_observed", "Show observed track", TRUE),
                       checkboxInput("plot_predicted", "Show predicted track", TRUE),
                       checkboxInput("enable_ticks", "Enable annotation ticks", TRUE),
                       conditionalPanel(
                         condition = "input.enable_ticks",
                         textInput("annot_ticks_input", "Ticks clusters ('all' or numbers comma-separated):", "all")
                       ),
                       checkboxInput("enable_kde", "Enable KDE layers", TRUE),
                       conditionalPanel(
                         condition = "input.enable_kde",
                         textInput("annot_kde_input", "KDE clusters ('all' or numbers comma-separated):", "all")
                       ),
                       actionButton("go", "Go!", class = "btn-primary")
          )
      ),
      div(id = "main-panel-container",
          girafeOutput("landscape_plot", width = "100%", height = "calc(100vh - 70px)")
      )
  )
)

# Server
server <- function(input, output, session) {
  plot_reactive <- eventReactive(input$go, {
    withProgress(message = "Building SCNA landscape…", value = 0, {
      incProgress(0.1, detail = "Filtering data")
      sel_chr   <- input$chr_input
      sel_coord <- if (input$genomic_coords == "") NULL else input$genomic_coords
      
      land_out <- filter_df(
        input_obj        = out_annot_list_processed,
        backbone_granges = backbone.100kb,
        type_input       = input$type_input,
        model_input      = input$model_input,
        chr_input        = sel_chr,
        coord_input      = sel_coord
      )
      df_land <- land_out$final_df
      req(nrow(df_land) > 0)
      
      shap_amp <- filter_df(shap.list, backbone.100kb, input$type_input, 'ampl', sel_chr, sel_coord)$final_df
      shap_del <- filter_df(shap.list, backbone.100kb, input$type_input, 'del',  sel_chr, sel_coord)$final_df
      shap_plot_list <- prepare_shap_to_plot(shap_amp, shap_del)
      
      outl <- prepare_landscape_to_plot(
        model_input         = parse_input_model(input$model_input),
        shap_plotting_list  = shap_plot_list,
        filtered_shap_output = land_out,
        filtered_landscape  = df_land,
        pred_list           = pred_list
      )
      fland <- outl$filtered_landscape
      gmask <- outl$genome_mask
      tmask <- outl$type_mask
      mmask <- outl$model_mask
      
      ticks_val <- if (!input$enable_ticks) FALSE else { txt <- tolower(input$annot_ticks_input); if (txt == "all") "all" else as.numeric(strsplit(txt, ",")[[1]]) }
      kde_val   <- if (!input$enable_kde) FALSE else { txt <- tolower(input$annot_kde_input);   if (txt == "all") "all" else as.numeric(strsplit(txt, ",")[[1]]) }
      
      incProgress(0.6, detail = "Generating plot")
      p <- landscape_plot_interactive_prediction(
        filtered_landscape  = fland,
        backbone.100kb      = backbone.100kb,
        genome_mask         = gmask,
        type_mask           = tmask,
        model_mask          = mmask,
        plot_observed       = input$plot_observed,
        plot_predicted      = input$plot_predicted,
        annot_to_plot_ticks = ticks_val,
        annot_to_plot_kde   = kde_val
      )
      
      incProgress(0.2, detail = "Configuring interactivity")
      p <- girafe_options(
        p,
        opts_tooltip(css = tooltip_css, delay_mouseover = 0, delay_mouseout = 0, offx = 10, offy = -10),
        opts_hover(css = hover_css)
      )
      
      incProgress(1)
      showNotification("Landscape ready!", type = "message", duration = 2)
      p
    })
  }, ignoreNULL = FALSE)
  
  output$landscape_plot <- renderGirafe({
    plot_reactive()
  })
}

# Run the app
shinyApp(ui = ui, server = server)
