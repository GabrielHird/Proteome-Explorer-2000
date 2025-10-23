##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

# Server ------------------------------------------------------------------

function(input, output, session) {
  

# Identifications ---------------------------------------------------------
  
  output$barplot <- renderPlot({
    req(input$bar_color, input$bar_sort)
    p <- Barplot_IDs(qf, i = "Proteins", color = input$bar_color, sort = input$bar_sort)
    print(p)
  })
  
  output$download_barplot <- downloadHandler(
    filename = function() {
      paste0(results_folder,"barplot_", Sys.Date(), ".png")
    },
    content = function(file) {
      png(file, width = 920, height = 550)
      p <- Barplot_IDs(qf, i = "Proteins", color = input$bar_color, sort = input$bar_sort)
      print(p)
      dev.off()
    }
  )

# Missing values ----------------------------------------------------------

  output$missMap <- renderPlot({
    p <- miss_map(qf, i = "Proteins")
    print(p)
  })


# Normalization ----------------------------------------------------------

  output$pre_norm_dens <- renderPlotly({
    p <- Density_plotter(qf, i = "Peptides")
    print(p)
  })
  
  output$post_norm_dens <- renderPlotly({
    p <- Density_plotter(qf, i = "PeptidesProccessed")
    print(p)
  })
  
  output$pre_norm_box <- renderPlot({
    p <- Norm_boxplotter(qf, i = "Peptides")
    print(p)
  })
  
  output$post_norm_box <- renderPlot({
    p <- Norm_boxplotter(qf, i = "PeptidesProccessed")
    print(p)
  })


# Protein -----------------------------------------------------------------

  output$aggregationPlot <- renderPlot({
    req(input$gene)
    p <- Aggregation_plotter(qf, input$gene)
    print(p)
  })
  
  output$prot_boxplot <- renderPlot({
    req(input$gene)
    p <- Boxplot_plotter(qf,input$gene)
    print(p)
  })
  
  output$volcanoPlotProt <- renderPlotly({
    p <- Volcano_plotter_highlight(qf, input$gene)
    ggplotly(p, tooltip = "text")
  })
  
  output$protein_info <- renderUI({
    
    # Get the protein data frame
    protein_data <- as.data.frame(rowData(qf[["Results"]]))
    
    # Identify the selected protein based on the input "gene"
    selected_gene <- input$gene
    protein_row <- protein_data[protein_data$Gene == selected_gene, ]
    
    # Use the UniProt ID(s) from the protein_row to query for function text.
    # Here we assume the UniProt accessions are stored in 'Protein.Group'
    ids <- protein_row$Protein.Group
    
    # Build a query to retrieve the function annotation using the "cc_function" field.
    query <- list("accession_id" = ids)
    function_df <- query_uniprot(query, 
                                 columns = c("accession", "cc_function"), 
                                 show_progress = FALSE)
    
    # Get the function text; if multiple entries exist, take the first.
    function_text <- if (nrow(function_df) > 0) {
      function_df[["Function [CC]"]][1]
    } else {
      "Not available"
    }
    
    # Build a custom HTML table string including the function text row.
    info_table <- paste0(
      "<table style='border-collapse: collapse; width: 100%;'>",
      "<tr><td><b>Protein:</b></td><td>", protein_row$Gene, "</td></tr>",
      "<tr><td><b>Full name:</b></td><td>", protein_row$First.Protein.Description, "</td></tr>",
      "<tr><td><b>Uniprot ID(s):</b></td><td>", protein_row$Protein.Group, "</td></tr>",
      "<tr><td><b>EntrezID:</b></td><td>", protein_row$ENTREZID, "</td></tr>",
      "<tr><td><b>Number of peptides:</b></td><td>", protein_row$Number.Peptides, "</td></tr>",
      "<tr><td><b>Significant:</b></td><td>",
      "<div style='display:inline-block; padding:4px 8px; border-radius:5px; border: 2px solid ",
      ifelse(tolower(as.character(protein_row$signif)) == "true", "green", "red"),
      "; color:",
      ifelse(tolower(as.character(protein_row$signif)) == "true", "green", "red"),
      ";'>",
      protein_row$signif,
      "</div>",
      "</td></tr>",
      "<tr><td><b>LogFC:</b></td><td>", protein_row$logFC, "</td></tr>",
      "<tr><td><b>Q-value:</b></td><td>", protein_row$qval, "</td></tr>",
      "<tr><td><b>Function:</b></td><td>", function_text, "</td></tr>",
      "</table>"
    )
    
    HTML(info_table)
  })


# Stats table -------------------------------------------------------------

  df <- as.data.frame(rowData(qf[["Results"]]))
  numeric_cols <- names(df)[sapply(df, is.numeric)]
  numeric_cols <- numeric_cols[numeric_cols != "Number.Peptides"]
  
  output$results_table <- renderDT({
    datatable(
      df,
      rownames = FALSE,
      class = 'display compact cell-border stripe',
      extensions = c("Scroller", "FixedHeader"),
      options = list(
        deferRender = TRUE,
        scrollY = 900,
        scrollX = TRUE,
        scrollCollapse = TRUE,
        scroller = TRUE,
        autoWidth = FALSE,
        columnDefs = list(
          list(width = '150px', targets = "_all")
        )
      )
    ) %>% formatRound(columns = numeric_cols, digits = 4)
  })
    

# Stats results -----------------------------------------------------------
  
  output$volcanoPlot <- renderPlot({
    p <- Volcano_plotter(qf, Title = "Volcano Plot")
    print(p)
  })
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste0(results_folder, "volcano_", Sys.Date(), ".png")
    },
    content = function(file) {
      png(file, width = 443, height = 610)
      p <- Volcano_plotter(qf, Title = "Volcano Plot")
      print(p)
      dev.off()
    }
  )


# Heatmap -----------------------------------------------------------------
  
  annotation_list <- reactive({
    req(input$selected_annotations)
    setNames(as.list(input$selected_annotations), input$selected_annotations)
  })
  
  # Full heatmap
  default_heatmap <- eventReactive(input$updateHeatmap, {
    make_heatmap(qf, 
                 i = "Results", 
                 annot_list = annotation_list(), 
                 only_significant = FALSE, 
                 n = NULL)
  }, ignoreNULL = FALSE)
  
  output$defaultHeatmapPlot <- renderPlot({
    default_heatmap()
  })
  
  # Stats heatmap
  custom_heatmap <- reactive({
    req(input$n_value)
    make_heatmap(qf, 
                 i = "Results", 
                 annot_list = annotation_list(), 
                 only_significant = TRUE, 
                 n = input$n_value)
  })
  
  output$customHeatmapPlot <- renderPlot({
    custom_heatmap()
  })
  
  observe({
    req(qf)
    slider_max <- nrow(rowData(qf[["Results"]]))
    updateSliderInput(
      session, 
      "n_value", 
      max = slider_max, 
      value = min(50, slider_max)
    )
  })
  
  # Download full heatmap
  
  output$download_heatmap_full <- downloadHandler(
    filename = function() {
      paste0("Heatmap_full_", Sys.Date(), ".png")
    },
    content = function(file) {
      png(file, width = 600, height = 600)
      print(default_heatmap()) 
      dev.off()
    }
  )
  
  # Download stats heatmap
  output$download_heatmap_sig <- downloadHandler(
    filename = function() {
      paste0("Heatmap_sig_", Sys.Date(), ".png")
    },
    content = function(file) {
      png(file, width = 600, height = 600)
      print(custom_heatmap())
      dev.off()
    }
  )

# PCA ---------------------------------------------------------------------
  
  output$pcaPlot <- renderPlot({
    p <- PCA_plotter(qf, i = "Proteins", color = "group")
    print(p)
  })
  
  # Download
  output$download_pca <- downloadHandler(
    filename = function() {
      paste0("PCA", Sys.Date(), ".png")
    },
    content = function(file) {
      png(file, width = 530, height = 360)
      p <- PCA_plotter(qf, i = "Proteins", color = "group", label = "none")
      print(p)
      dev.off()
    }
  )
    


}
