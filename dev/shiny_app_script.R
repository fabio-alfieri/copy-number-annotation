library(shiny)
library(ggiraph)
library(htmlwidgets)
library(ggplot2)
library(dplyr)
library(shinycssloaders)

setwd("/Users/gabry/OneDrive/Desktop/shiny_app/")

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

# CSS styles for tooltips and hover
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


ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      html, body {
        height: 100%;
        margin: 0; padding: 0;
        overflow: hidden; /* no vertical scroll on page */
      }
      .container-fluid {
        height: 100vh; /* full viewport height */
        display: flex;
        flex-direction: column;
      }
      .title-panel {
        flex: 0 0 auto;
        padding-bottom: 5px;
      }
      #flex-container {
        display: flex;
        flex-grow: 1;
        overflow: hidden;
      }
      #sidebar-container {
        width: 280px;
        overflow-y: auto; /* scroll inside sidebar if needed */
        flex-shrink: 0;
        transition: width 0.3s ease;
        display: flex;
        flex-direction: column;
      }
      #sidebar-container.collapsed {
        width: 40px !important;
      }
      #sidebar-container.collapsed .well {
        display: none;
      }
      #sidebar-container .btn {
        width: 100%;
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
      /* girafeOutput fills all vertical space */
      #landscape_plot {
        flex-grow: 1;
        min-height: 0;
      }
    ")),
    # JS to toggle sidebar class on button click
    tags$script(HTML("
      $(document).on('shiny:connected', function() {
        $('#toggle_sidebar').on('click', function() {
          $('#sidebar-container').toggleClass('collapsed');
          $('#main-panel-container').toggleClass('expanded');
        });
      });
    "))
  ),
  
  titlePanel("SCNA Landscape Plot"),
  
  div(
    id = "flex-container",
    div(
      id = "sidebar-container",
      class = "sidebar-panel",
      actionButton("toggle_sidebar", "☰", class = "btn btn-secondary", style = "margin-bottom: 10px; width: 100%;"),
      sidebarPanel(
        width = 12,
        selectInput("type_input", "Select cancer type:",
                    choices = c("STAD","GBMLGG","COADREAD","KIRP","KIRC",
                                "OV","ESCA","LUAD","LUSC","PAAD","BRCA"),
                    selected = "BRCA"),
        
        selectInput("chr_input", "Select chromosome(s):",
                    choices = paste0("chr", 1:22),
                    selected = paste0("chr", 1:22),
                    multiple = TRUE,
                    selectize = TRUE),
        
        textInput("genomic_coords", 
                  "Genomic coordinates (chr:start-end):", 
                  ""),
        
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
      )
    ),
    
    div(
      id = "main-panel-container",
      # Here's the key height adjustment:
      girafeOutput("landscape_plot", width = "100%", height = "calc(100vh - 70px)")
      # 70px chosen as approximate height of title panel + padding; adjust if needed
    )
  )
)

server <- function(input, output, session) {
  
  plot_reactive <- eventReactive(input$go, {
    withProgress(message = "Building SCNA landscape…", value = 0, {
      incProgress(0.1, detail = "Filtering data")
      
      selected_chr <- input$chr_input
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
