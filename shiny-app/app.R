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
library(data.table)


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

## ---------------------------------------------------------------------
## Filename helper
## ---------------------------------------------------------------------
## Every exported file (PDF, HTML bundle, BED) is named after the CURRENT
## selection, so a BRCA / chr1 export always looks like:
##   CNA_landscape_BRCA_ampl_Mid-length_chr1_20260903.pdf
## regardless of which download button was used.

sanitize_component <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_-]+", "", x)
  x
}

build_export_basename <- function(type_input, model_input, cluster_input, chr_input, genomic_coords) {
  
  region_part <- if (!is.null(genomic_coords) && nzchar(genomic_coords)) {
    sanitize_component(gsub("[:-]", "_", genomic_coords))
  } else if (length(chr_input) == 22) {
    "WHOLEGENOME"
  } else if (length(chr_input) == 1) {
    sanitize_component(chr_input)
  } else {
    paste0(length(chr_input), "chroms")
  }
  
  paste(
    "CNA_landscape",
    sanitize_component(type_input),
    sanitize_component(model_input),
    sanitize_component(cluster_input),
    region_part,
    format(Sys.Date(), "%Y%m%d"),
    sep = "_"
  )
}

server <- function(input, output, session) {
  
  observeEvent(input$annot_to_plot, {
    
    sel_coord <- trimws(input$genomic_coords)
    sel_coord <- if (identical(sel_coord, "")) NULL else sel_coord
    
    land_out <- tryCatch(
      filter_df(
        input_obj        = meta_list_get,
        backbone_granges = backbone.100kbp_granges,
        cluster_input    = input$cluster_input,
        type_input       = input$type_input,
        model_input      = input$model_input,
        chr_input        = input$chr_input,
        coord_input      = sel_coord
      ),
      error = function(e) NULL
    )
    
    if (is.null(land_out)) return(NULL)
    
    df_tmp <- land_out$final_df
    if (nrow(df_tmp) == 0) return(NULL)
    
    annot <- input$annot_to_plot
    
    # Same union() logic as landscape_plot_observed_prediction(): the fixed
    # vocabulary keeps the ordering stable, anything unexpected in the data
    # is still offered to the user rather than being silently unselectable.
    classes <- if (annot == "annot_final") {
      union(ANNOT_FINAL_CLASSES, sort(unique(as.character(df_tmp[[annot]]))))
    } else {
      sort(unique(df_tmp[[annot]]))
    }
    
    updateSelectizeInput(
      session, "annot_ticks_input",
      choices = classes,
      selected = classes
    )
    
    updateSelectizeInput(
      session, "annot_kde_input",
      choices = classes,
      selected = classes
    )
  })
  
  ## ---------------------------------------------------------------------
  ## 1) Filter the data ONCE per "Go" click.
  ## ---------------------------------------------------------------------
  filtered_data <- eventReactive(input$go, {
    withProgress(message = "Filtering data…", value = 0.2, {
      
      sel_coord <- trimws(input$genomic_coords)
      sel_coord <- if (identical(sel_coord, "")) NULL else sel_coord
      
      # A manually-typed region can be malformed (bad pattern, start > end,
      # a chromosome outside chr1-chr22, etc.) - filter_df()/parse_input_coord()
      # raise a clear, user-facing error message via stop() in that case;
      # we surface it as a notification instead of letting it crash the
      # reactive chain, so the user can just fix the text box and retry.
      land_out <- tryCatch(
        filter_df(
          input_obj        = meta_list_get,
          backbone_granges = backbone.100kbp_granges,
          cluster_input    = input$cluster_input,
          type_input       = input$type_input,
          model_input      = input$model_input,
          chr_input        = input$chr_input,
          coord_input      = sel_coord
        ),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 8)
          NULL
        }
      )
      
      if (is.null(land_out)) return(NULL)
      
      if (isTRUE(land_out$region_clamped)) {
        showNotification(
          "The requested interval went beyond the chromosome boundaries in the loaded data: it was trimmed to the available limits.",
          type = "warning", duration = 6
        )
      }
      
      df_land <- land_out$final_df
      df_land$is_centromere <- df_land$binID %in% centromere_table_out$binID
      
      if (nrow(df_land) == 0) {
        showNotification("No data available with current filters.", type = "error", duration = 4)
        return(NULL)
      }
      
      list(
        df_land        = df_land,
        annot_col      = input$annot_to_plot,
        ticks_val      = if (!input$enable_ticks) FALSE else {
          txt <- trimws(input$annot_ticks_input)
          if (identical(tolower(txt), "all")) "all" else unlist(strsplit(txt, "\\s*,\\s*"))
        },
        kde_val        = if (!input$enable_kde) FALSE else {
          txt <- trimws(input$annot_kde_input)
          if (identical(tolower(txt), "all")) "all" else unlist(strsplit(txt, "\\s*,\\s*"))
        },
        cluster_input  = input$cluster_input,
        type_input     = input$type_input,
        model_input    = input$model_input,
        chr_input      = input$chr_input,
        genomic_coords = trimws(input$genomic_coords),
        plot_observed  = input$plot_observed,
        plot_predicted = input$plot_predicted
      )
    })
  }, ignoreNULL = FALSE)
  
  ## ---------------------------------------------------------------------
  ## 2) Build the plot ONCE per "Go" click (both interactive & static
  ##    versions come out of a single call), reused by the screen render
  ##    AND by every download handler below - nothing gets re-filtered or
  ##    re-plotted from scratch on download anymore.
  ## ---------------------------------------------------------------------
  built_plots <- reactive({
    fd <- filtered_data()
    req(fd)
    
    withProgress(message = "Building CNA landscape…", value = 0.5, {
      
      out <- landscape_plot_observed_prediction(
        filtered_landscape  = fd$df_land,
        cluster_mask        = fd$cluster_input,
        genome_mask         = fd$chr_input,
        type_mask           = fd$type_input,
        model_mask          = fd$model_input,
        plot_observed       = fd$plot_observed,
        plot_predicted      = fd$plot_predicted,
        annot_to_plot       = fd$annot_col,
        annot_to_plot_ticks = fd$ticks_val,
        annot_to_plot_kde   = fd$kde_val,
        make.interactive    = "both"
      )
      
      incProgress(1)
      out
    })
  })
  
  ## Selection-consistent basename shared by every downloadHandler below.
  current_basename <- function() {
    fd <- filtered_data()
    req(fd)
    build_export_basename(fd$type_input, fd$model_input, fd$cluster_input, fd$chr_input, fd$genomic_coords)
  }
  
  output$landscape_plot <- renderGirafe({
    bp <- built_plots()
    req(bp)
    showNotification("Landscape ready!", type = "message", duration = 2)
    bp$interactive
  })
  
  output$download_html <- downloadHandler(
    filename = function() paste0(current_basename(), ".zip"),
    content = function(file) {
      td <- tempfile("CNA_export_")
      dir.create(td)
      on.exit(unlink(td, recursive = TRUE), add = TRUE)
      
      bp <- built_plots()
      if (is.null(bp)) stop("No plot to export.")
      
      base_name <- current_basename()
      html_file <- file.path(td, paste0(base_name, ".html"))
      
      saveWidget(bp$interactive, html_file, selfcontained = FALSE)
      
      assets_dir <- sub("\\.html$", "_files", html_file)
      files_to_zip <- basename(html_file)
      if (dir.exists(assets_dir)) {
        files_to_zip <- c(files_to_zip, basename(assets_dir))
      }
      
      oldwd <- setwd(td)
      on.exit(setwd(oldwd), add = TRUE)
      zipfile <- file.path(td, paste0(base_name, "_bundle.zip"))
      zip::zip(zipfile = zipfile, files = files_to_zip, recurse = TRUE)
      
      file.copy(zipfile, file, overwrite = TRUE)
    }
  )
  
  output$download_pdf <- downloadHandler(
    filename = function() paste0(current_basename(), ".pdf"),
    content = function(file) {
      bp <- built_plots()
      if (is.null(bp) || !inherits(bp$static, "ggplot")) {
        stop("PDF export failed: Plot is not a ggplot object.")
      }
      
      # Bigger canvas, but the legend was made compact in
      # landscape_plot_observed_prediction() so the extra size goes to the
      # plot panel rather than growing the legend proportionally.
      ggsave(file, bp$static, device = "pdf", width = 16, height = 8, units = "in", dpi = 300)
    }
  )
  
  output$download_annotation <- downloadHandler(
    filename = function() paste0(current_basename(), ".bed"),
    content = function(file) {
      fd <- filtered_data()
      req(fd)
      
      df <- fd$df_land
      annot_col <- fd$annot_col
      
      if (!(annot_col %in% colnames(df))) {
        stop("Selected annotation column not found in the filtered data.")
      }
      
      # Only chr, start, end and the selected annotation - nothing else.
      bed_df <- data.frame(
        chrom = as.character(df$chr),
        start = pmax(as.integer(df$start) - 1L, 0L),   # BED is 0-based, half-open
        end   = as.integer(df$end),
        annot = as.character(df[[annot_col]]),
        stringsAsFactors = FALSE
      )
      
      chr_order <- paste0("chr", 1:22)
      bed_df <- bed_df[order(match(bed_df$chrom, chr_order), bed_df$start), ]
      
      readr::write_tsv(bed_df, file, col_names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)