library(shiny)
library(ggiraph)
library(glue)

palette_custom <- list(
  background = "#f0f4fa",
  primary = "#2a78e8",
  primary_dark = "#155ab3",
  border = "#dbe9ff"
)

tooltip_css <- "
[class^='tooltip_svg_'],
[class^='tooltip_svg_'] * {
  background: rgba(42, 120, 232, 0.9) !important; /* blue bg instead of black */
  color: #ffffff !important;
  padding: 2px 6px !important;
  border: none !important;
  box-shadow: none !important;
  border-radius: 4px !important;
  font-family: 'Roboto', sans-serif !important;
  font-size: 12px !important;
  line-height: 1.2 !important;
  white-space: nowrap !important;
}
[class^='tooltip_svg_'] rect { fill: none !important; stroke: none !important; }
[class^='tooltip_svg_'] { border: none !important; box-shadow: none !important; background: none !important; }
"

hover_css <- "opacity: 0.9 !important; transform: translateY(0) !important;"

ui <- fluidPage(
  tags$head(
    tags$style(HTML(glue::glue("
      html, body {{height:100%; margin:0; padding:0; overflow:hidden; font-family:'Segoe UI', sans-serif; background-color:{palette_custom$background};}}
      .container-fluid {{padding:0 !important; margin:0 !important; height:100vh; display:flex; flex-direction:column;}}
      .title-panel {{flex:0 0 auto; padding:10px 0; font-size:22px; font-weight:bold; background-color:{palette_custom$background};
                      border-bottom:1px solid {palette_custom$border}; text-align:center; width:100vw; margin:0; box-sizing:border-box; color:{palette_custom$primary_dark};}}
      #flex-container {{display:flex; flex-grow:1; overflow:hidden;}}
      #sidebar-container {{width:300px; overflow-y:auto; flex-shrink:0; transition: width 0.3s ease; display:flex; flex-direction:column;
                           background-color:{palette_custom$background}; border-right:1px solid {palette_custom$border}; padding:10px; align-items:center;}}
      #sidebar-container.collapsed {{width:50px !important; padding:10px 5px; display:flex !important; flex-direction:column; align-items:center; justify-content:flex-start;}}
      #sidebar-container.collapsed .well,
      #sidebar-container.collapsed .form-group,
      #sidebar-container.collapsed .checkbox,
      #sidebar-container.collapsed .shiny-input-container {{display:none !important;}}
      #sidebar-container.collapsed .btn {{width:100%; margin:0 0 10px 0; font-size:18px; padding:6px 0; text-align:center;}}
      #sidebar-container .btn {{width:auto; margin:0 auto;}}
      #sidebar-container.collapsed img {{display: none !important;}}
      #main-panel-container {{flex-grow:1; padding-left:20px; overflow:hidden; display:flex; flex-direction:column;}}
      #main-panel-container.expanded {{margin-left:0;}}
      #landscape_plot {{flex-grow:1; min-height:0;}}
      .form-control, .selectize-input {{border-radius:6px; border-color:{palette_custom$border};}}
      .btn-primary {{background-color:{palette_custom$primary}; border-color:{palette_custom$primary}; color:white;}}
      .btn-primary:hover {{background-color:{palette_custom$primary_dark}; border-color:{palette_custom$primary_dark}; color:white;}}
      .checkbox input[type='checkbox'] {{accent-color:{palette_custom$primary};}}
      .shiny-input-container label {{font-weight:600; font-size:14px; margin-bottom:4px; color:{palette_custom$primary_dark};}}
    "))),
    tags$script(HTML("
      $(document).on('shiny:connected', function(){
         $('#toggle_sidebar').on('click', function(){
           $('#sidebar-container').toggleClass('collapsed');
           $('#main-panel-container').toggleClass('expanded');
         });
       });
    "))
  ),
  
  div(class = "title-panel", "PENNE: Prediction & Explanation of Non neutral copy Number Events"),
  
  div(id = "flex-container",
      div(id = "sidebar-container", class = "sidebar-panel",
          actionButton("toggle_sidebar", "☰", class = "btn btn-secondary", 
                       style = "margin-bottom:10px; width:100%; font-size:18px;"),
          tags$img(src = "penne.png", height = "180px", style = "margin-bottom:15px; display:block; margin-left:auto; margin-right:auto;"),
          sidebarPanel(width = 12,
                       selectInput("type_input", "Select cancer type:",
                                   choices = c("STAD","GBMLGG","COADREAD","KIRP","KIRC","OV",
                                               "ESCA","LUAD","LUSC","PAAD","BRCA"),
                                   selected = "BRCA"),
                       
                       selectInput("model_input", "Select model:",
                                   choices = c("Amplification" = "ampl", "Deletion" = "del"),
                                   selected = "ampl"),
                       
                       selectInput(
                         "cluster_input", 
                         "Select CNA cluster:", 
                         # The value sent to the server is the INTERNAL id
                         # used in the file names ("No-cluster"); only the
                         # displayed label is "All CNAs". Keep in sync with
                         # CLUSTER_LABELS in 0_LoadData.R.
                         choices = c("Mid-length", "Small-scale", "Arm-level", "Chromosome-level",
                                     "All CNAs" = "No-cluster"),
                         selected = "Mid-length"
                       ),
                       
                       selectInput("chr_input", "Select chromosome(s):",
                                   choices = paste0("chr", 1:22),
                                   selected = paste0("chr", 1:22),
                                   multiple = TRUE,
                                   selectize = TRUE),
                       
                       div(
                         class = "shiny-input-container",
                         tags$label("Select annotation to plot:"),
                         
                         tags$select(
                           id = "annot_to_plot",
                           class = "form-control",
                           tags$option(value = "annot_final", "6 class annotation"),
                           tags$option(value = "top1", "Top 1 Feature")
                         )
                       ),
                       
                       div(
                         class = "shiny-input-container",
                         tags$label("Select genomic coordinates (chr:start-end):"),
                         
                         tags$input(
                           id = "genomic_coords",
                           type = "text",
                           class = "form-control",
                           value = "",
                           placeholder = "e.g. chr1:1-5000000"
                         ),
                         tags$small(
                           style = "color:#6b7280; display:block; margin-top:4px;",
                           "Optional: if filled in, it takes priority over the chromosomes selected above. Leave empty to use the chromosome selection."
                         )
                       ),
                       checkboxInput("plot_observed", "Show observed track", TRUE),
                       checkboxInput("plot_predicted", "Show predicted track", TRUE),
                       checkboxInput("enable_ticks", "Show annotation ticks", FALSE),
                       conditionalPanel(
                         condition = "input.enable_ticks",
                         selectizeInput(
                           inputId = "annot_ticks_input",   # ← unchanged
                           label = "Annotation classes to plot:",
                           choices = c(
                             "Negative Selection", 
                             "Excess of Prediction - Likely Negative Selection", 
                             "Positive Selection", 
                             "Excess of Observation - Likely Positive Selection", 
                             "Occurrence", 
                             "No Detectable Force"
                           ),
                           selected = c(
                             "Negative Selection", 
                             "Excess of Prediction - Likely Negative Selection", 
                             "Positive Selection", 
                             "Excess of Observation - Likely Positive Selection", 
                             "Occurrence", 
                             "No Detectable Force"
                           ),
                           multiple = TRUE
                         )
                       ),
                       checkboxInput("enable_kde", "Show annotation density", TRUE),
                       conditionalPanel(
                         condition = "input.enable_kde",
                         selectizeInput(
                           inputId = "annot_kde_input",
                           label = "Annotation classes to plot:",
                           choices = c(
                             "Negative Selection", 
                             "Excess of Prediction - Likely Negative Selection", 
                             "Positive Selection", 
                             "Excess of Observation - Likely Positive Selection", 
                             "Occurrence", 
                             "No Detectable Force"
                           ),
                           selected = c(
                             "Negative Selection", 
                             "Excess of Prediction - Likely Negative Selection", 
                             "Positive Selection", 
                             "Excess of Observation - Likely Positive Selection", 
                             "Occurrence", 
                             "No Detectable Force"
                           ),
                           multiple = TRUE
                         )
                       ),
                       
                       # ---- Each button on its own row ----
                       div(style = "margin-top:10px;",
                           actionButton("go", "Go!", class = "btn-primary")
                       ),
                       div(style = "margin-top:10px;",
                           downloadButton("download_html", "Download HTML", class = "btn btn-secondary")
                       ),
                       div(style = "margin-top:10px;",
                           downloadButton("download_pdf", "Download PDF", class = "btn btn-secondary")
                       ),
                       div(
                         style = "margin-top:10px;",
                         title = "Download the annotation (BED) for the current selection only",
                         downloadButton("download_annotation", "Download Annotation", class = "btn btn-secondary")
                       )
          )
      ),
      div(id = "main-panel-container",
          girafeOutput("landscape_plot", width = "100%", height = "calc(100vh - 70px)")
      )
  )
)