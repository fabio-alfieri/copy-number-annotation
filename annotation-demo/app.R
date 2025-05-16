library(shiny)
library(htmlwidgets)
library(htmltools)
library(ggnewscale)

# source('dev/MAIN_dynamic_plotting.R')

ui <- fluidPage(
    # 1) CSS to make the plot container fill the window
    tags$head(
        tags$style(HTML("
      html, body { height: 100%; margin: 0; padding: 0; }
      #plot_container {
        width: 100%;
        height: calc(100vh - 120px); /* adjust 120px to match title + input row height */
      }
      #plot_container .html-widget {
        width: 100% !important;
        height: 100% !important;
      }
    "))
    ),
    
    titlePanel("Copy Number Annotation Browser - CATCHY NAME NEEDED!"),
    
    fluidRow(
        column(width = 2,
               selectInput(
                   inputId  = "type_input",
                   label    = "Tumor type:",
                   choices  = c("BRCA", "COADREAD", "ESCA", 
                                "GBMLGG", "KIRC", "KIRP", 
                                "LUAD", "LUSC", "OV", "PAAD", "STAD"),
                   selected = "BRCA"
               ),
               selectInput(
                   inputId  = "model_input",
                   label    = "Model type:",
                   choices  = c("Amplification", "Deletions"),
                   selected = "Amplification"
               ),
               textInput(
                   inputId     = "coord_input",
                   label       = "Genomic coordinates:",
                   value       = "chr1:1-1000000",
                   placeholder = "e.g. chr1:100000-200000"
               )
        ),
        
        column(width = 10,
               # the resizable container for your widget
               div(id = "plot_container",
                   uiOutput("distPlot")
               )
        )
    )
)

server <- function(input, output) {
    
    widget_reactive <- reactive({
        req(input$type_input, input$model_input, input$coord_input)
    
        htmlwidgets::prependContent(
          p,
          tags$head(tags$style(HTML(tooltip_css)))
        )
    })
    
    output$distPlot <- renderUI({
        widget_reactive()
    })
    
    # output$distPlot <- renderUI({
    #     p2
    # })
}

# Run the application 
shinyApp(ui = ui, server = server)

shinylive::export(appdir = "annotation-demo/", destdir = "docs/")
