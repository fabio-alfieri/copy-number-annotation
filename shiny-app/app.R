set.seed(123)

library(shiny)
library(ggiraph)
library(htmlwidgets)
library(ggplot2)
library(dplyr)
library(shinycssloaders)
library(zip)
library(glue)
library(tidyverse)
library(IRanges) 
library(stringr)
library(readr)
library(htmltools)
library(ggnewscale)
library(GenomicRanges)


source('0_LoadData.R')
source('1_dynamic_plotting_functions.R')
source('1bis_parse_input_data_FA.R')
source('plotting_helpers.R')
source('landscape_plot_observed_prediction.R')
source('ui.R')

centromere_table_out <- process_centromere_table(
  centromere_table = centromere_table, 
  backbone.100kb = backbone.100kbp_granges
)

server <- function(input, output, session) {
  
  current_plot <- reactiveVal(NULL)
  
  plot_reactive <- eventReactive(input$go, {
    withProgress(message = "Building CNA landscape…", value = 0, {
      
      incProgress(0.1, detail = "Filtering data")
      
      sel_coord <- if (input$genomic_coords == "") NULL else input$genomic_coords
      
      land_out <- filter_df(
        input_obj        = meta_list,
        backbone_granges = backbone.100kbp_granges,
        cluster_input    = input$cluster_input,
        type_input       = input$type_input,
        model_input      = input$model_input,
        chr_input        = input$chr_input,
        coord_input      = sel_coord
      )
      
      df_land <- land_out$final_df
      df_land$is_centromere <- df_land$binID %in% centromere_table_out$binID
      
      if (nrow(df_land) == 0) {
        showNotification("No data available with current filters.", type = "error", duration = 4)
        return(NULL)
      }
      
      annot_col <- "annot_final"
      
      ticks_val <- if (!input$enable_ticks) FALSE else {
        txt <- trimws(input$annot_ticks_input)
        if (identical(tolower(txt), "all")) "all" else unlist(strsplit(txt, "\\s*,\\s*"))
      }
      
      kde_val <- if (!input$enable_kde) FALSE else {
        txt <- trimws(input$annot_kde_input)
        if (identical(tolower(txt), "all")) "all" else unlist(strsplit(txt, "\\s*,\\s*"))
      }
      
      incProgress(0.6, detail = "Generating plot")
      
      p <- landscape_plot_observed_prediction(
        filtered_landscape  = df_land,
        cluster_mask        = input$cluster_input,
        genome_mask         = input$chr_input,
        type_mask           = input$type_input,
        model_mask          = input$model_input,
        plot_observed       = input$plot_observed,
        plot_predicted      = input$plot_predicted,
        annot_to_plot       = annot_col,
        annot_to_plot_ticks = ticks_val,
        annot_to_plot_kde   = kde_val,
        make.interactive    = TRUE
      )
      
      incProgress(0.2, detail = "Configuring interactivity")
      
      if (inherits(p, "gg")) { 
        p <- girafe_options(
          girafe(ggobj = p, fonts = list(sans = "Roboto"), width_svg = 10, height_svg = 6),
          opts_tooltip(css = tooltip_css, delay_mouseover = 0, delay_mouseout = 0, offx = 10, offy = -10),
          opts_hover(css = hover_css),
          opts_toolbar(saveaspng = FALSE)
        )
      }
      
      incProgress(1)
      showNotification("Landscape ready!", type = "message", duration = 2)
      
      p
    })
  }, ignoreNULL = FALSE)
  
  output$landscape_plot <- renderGirafe({
    p <- plot_reactive()
    
    on.exit({
      current_plot(NULL) 
      gc()
    })
    
    p
  })
  
  output$download_html <- downloadHandler(
    filename = function() {
      paste0("CNA_landscape_", Sys.Date(), "_", input$cluster_input, "_", input$type_input, "_", input$model_input, ".zip")
    },
    content = function(file) {
      td <- tempfile("CNA_export_")
      dir.create(td)
      on.exit(unlink(td, recursive = TRUE), add = TRUE)
      
      html_file <- file.path(td, "CNA_landscape.html")
      widget_obj <- plot_reactive()
      if (is.null(widget_obj)) stop("No plot to export.")
      
      saveWidget(widget_obj, html_file, selfcontained = FALSE)
      
      assets_dir <- sub("\\.html$", "_files", html_file)
      files_to_zip <- c("CNA_landscape.html")
      if (dir.exists(assets_dir)) {
        files_to_zip <- c(files_to_zip, basename(assets_dir))
      }
      
      oldwd <- setwd(td)
      on.exit(setwd(oldwd), add = TRUE)
      zipfile <- file.path(td, "CNA_landscape_bundle.zip")
      zip::zip(zipfile = zipfile, files = files_to_zip, recurse = TRUE)
      
      file.copy(zipfile, file, overwrite = TRUE)
    }
  )
  
  output$download_pdf <- downloadHandler(
    filename = function() {
      paste0("CNA_landscape_", Sys.Date(), "_", input$cluster_input, "_", input$type_input, "_", input$model_input, ".pdf")
    },
    content = function(file) {
      sel_coord <- if (input$genomic_coords == "") NULL else input$genomic_coords
      
      land_out <- filter_df(
        input_obj        = meta_list,
        backbone_granges = backbone.100kbp_granges,
        cluster_input    = input$cluster_input,
        type_input       = input$type_input,
        model_input      = input$model_input,
        chr_input        = input$chr_input,
        coord_input      = sel_coord
      )
      
      df_land <- land_out$final_df
      df_land$is_centromere <- df_land$binID %in% centromere_table_out$binID
      
      if (nrow(df_land) == 0) stop("No data to export.")
      
      annot_col <- "annot_final"
      
      ticks_val <- if (!input$enable_ticks) FALSE else {
        txt <- trimws(input$annot_ticks_input)
        if (identical(tolower(txt), "all")) "all" else unlist(strsplit(txt, "\\s*,\\s*"))
      }
      
      kde_val <- if (!input$enable_kde) FALSE else {
        txt <- trimws(input$annot_kde_input)
        if (identical(tolower(txt), "all")) "all" else unlist(strsplit(txt, "\\s*,\\s*"))
      }
      
      p <- landscape_plot_observed_prediction(
        filtered_landscape  = df_land,
        cluster_mask        = input$cluster_input,
        genome_mask         = input$chr_input,
        type_mask           = input$type_input,
        model_mask          = input$model_input,
        plot_observed       = input$plot_observed,
        plot_predicted      = input$plot_predicted,
        annot_to_plot       = annot_col,
        annot_to_plot_ticks = ticks_val,
        annot_to_plot_kde   = kde_val,
        make.interactive    = FALSE
      )
      
      if (!inherits(p, "ggplot")) {
        stop("PDF export failed: Plot is not a ggplot object.")
      }
      
      ggsave(file, p, device = "pdf", width = 10, height = 6)
    }
  )
}

shinyApp(ui = ui, server = server)

