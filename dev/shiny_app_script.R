library(shiny)
library(ggiraph)
library(htmlwidgets)
library(ggplot2)
library(dplyr)
library(shinycssloaders)

# Load helper functions
source('dev/0_LoadData.R')  # loads shap.list, out_annot_list, etc.
source('dev/1_dynamic_plotting_functions.R')  # contains landscape_plot_interactive_prediction
source('dev/1bis_parse_input_data_FA.R')    # parse_input_data(), filter_df(), etc.

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
tooltip_css <- "[class^='tooltip_svg_'] { background: transparent !important; color: #fafafa; padding: 2px 6px; border-radius: 2px; font-family: 'Roboto', sans-serif; font-size: 12px; line-height: 1.2; box-shadow: none !important; transform: translateY(2px); white-space: nowrap; }"
hover_css   <- "opacity: 0.9 !important; transform: translateY(0) !important;"

# UI definition
ui <- fluidPage(
  tags$head(
    tags$style(HTML(
      "html, body {height:100%; margin:0; padding:0; overflow:hidden;}\
       .container-fluid {height:100vh; display:flex; flex-direction:column;}\
       .title-panel {flex:0 0 auto; padding-bottom:5px;}\
       #flex-container {display:flex; flex-grow:1; overflow:hidden;}\
       #sidebar-container {width:280px; overflow-y:auto; flex-shrink:0; transition:width 0.3s ease; display:flex; flex-direction:column;}\
       #sidebar-container.collapsed {width:40px !important;}\
       #sidebar-container.collapsed .well {display:none;}\
       #sidebar-container .btn {width:100%;}\
       #main-panel-container {flex-grow:1; padding-left:20px; overflow:hidden; display:flex; flex-direction:column;}\
       #main-panel-container.expanded {margin-left:0;}\
       #landscape_plot {flex-grow:1; min-height:0;}"
    )),
    tags$script(HTML(
      "$(document).on('shiny:connected', function(){\n         $('#toggle_sidebar').on('click', function(){\n           $('#sidebar-container').toggleClass('collapsed');\n           $('#main-panel-container').toggleClass('expanded');\n         });\n       });"
    ))
  ),
  titlePanel("SCNA Landscape Plot"),
  div(id = "flex-container",
      div(id = "sidebar-container", class = "sidebar-panel",
          actionButton("toggle_sidebar", "☰", class = "btn btn-secondary", style = "margin-bottom:10px; width:100%;"),
          sidebarPanel(width = 12,
                       selectInput("type_input", "Select cancer type:",
                                   choices = c("STAD","GBMLGG","COADREAD","KIRP","KIRC","OV","ESCA","LUAD","LUSC","PAAD","BRCA"),
                                   selected = "BRCA"),
                       selectInput("model_input", "Select model:",
                                   choices = c("ampl","del"),
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

# Server logic
server <- function(input, output, session) {
  plot_reactive <- eventReactive(input$go, {
    withProgress(message = "Building SCNA landscape…", value = 0, {
      # Progress: Filtering
      incProgress(0.1, detail = "Filtering data")
      sel_chr   <- input$chr_input
      sel_coord <- if (input$genomic_coords == "") NULL else input$genomic_coords
      
      # Filter landscape data for selected model
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
      
      # Filter SHAP for both models for ticks/KDE
      shap_amp <- filter_df(shap.list, backbone.100kb, input$type_input, 'ampl', sel_chr, sel_coord)$final_df
      shap_del <- filter_df(shap.list, backbone.100kb, input$type_input, 'del',  sel_chr, sel_coord)$final_df
      shap_plot_list <- prepare_shap_to_plot(shap_amp, shap_del)
      
      # Prepare landscape with obs/pred and masks
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
      
      # Progress: Building
      # Parse tick/KDE cluster selections
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
  
  # Render the interactive plot
  output$landscape_plot <- renderGirafe({
    plot_reactive()
  })
}

# Launch the application
shinyApp(ui = ui, server = server)
