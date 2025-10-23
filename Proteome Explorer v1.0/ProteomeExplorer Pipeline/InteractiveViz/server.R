##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

# Server ------------------------------------------------------------------

function(input, output, session) {

  pipeline_dir <- normalizePath("..", winslash = "/", mustWork = TRUE)
  repo_root <- normalizePath(file.path(pipeline_dir, "..", ".."), winslash = "/", mustWork = TRUE)

  default_config <- list(
    Project_name = "PExA Lungcancer v6",
    Analysis_run = "EBP_v1_detect",
    First_pass = TRUE,
    Run_normalizer = TRUE,
    Run_DEA = TRUE,
    Run_GSEA = TRUE,
    run_pathfind = FALSE,
    save_plots = TRUE,
    export_data = TRUE,
    subtype = FALSE,
    Run_Interactive = FALSE,
    Treatment_group = "Cancer",
    Reference_group = "Normal",
    Import_qval = 0.01,
    Import_pg_qval = 0.01,
    Group_treshold = 0.75,
    Global_treshold = 0.6,
    Treshold_level = "Group",
    peptide_level = TRUE,
    Norm_method = "CycLoess",
    Batch_corr = "eigenms",
    Agg_method = "robustsummary",
    Impute = FALSE,
    impute_MAR = NA,
    impute_MNAR = NA,
    impute_method = NA,
    DEA_method = "msqrob",
    Study_design = "unpaired",
    alpha = 0.05,
    logFC_cutoff = NA,
    Org_db = "org.Hs.eg.db",
    Boxplot_prot = c("CSTA"),
    Agg_prot = c("CSTA"),
    Group_colors = c("Cancer" = "#FABC3C", "Normal" = "#006D77"),
    heatmap_annot = c(
      "Group" = "group",
      "Batch" = "Extraction_batch",
      "Study" = "Study",
      "Membrane" = "PExA_membrane"
    ),
    Expression_lvl_color = c(OVER = "#CB4335", UNDER = "#2E86C1"),
    DIA_nn_path = "./Data/DIA-NN output/diann_report_EBP.tsv",
    FASTA_path = "./Data/FASTA file/EBP.fasta",
    Sample_data_path = "./Data/Sample Metadata/EBP_sample_data.xlsx",
    Your_name = "Gabriel Hirdman",
    Lab_name = "Lindstedt Lab",
    Contact = "gabriel.hirdman@med.lu.se"
  )

  parse_comma_list <- function(x) {
    items <- trimws(unlist(strsplit(x %||% "", ",")))
    items <- items[nzchar(items)]
    if (length(items) == 0) character() else unique(items)
  }

  parse_named_lines <- function(x) {
    lines <- strsplit(x %||% "", "\n", fixed = TRUE)[[1]]
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
    if (length(lines) == 0) {
      return(character())
    }
    split_lines <- strsplit(lines, "=", fixed = TRUE)
    vals <- vapply(split_lines, function(parts) {
      parts <- trimws(parts)
      if (length(parts) < 2) "" else parts[2]
    }, character(1))
    names(vals) <- vapply(split_lines, function(parts) {
      parts <- trimws(parts)
      parts[1]
    }, character(1))
    vals <- vals[nzchar(vals) & nzchar(names(vals))]
    vals
  }

  empty_to_na <- function(x) {
    if (is.null(x) || x == "") {
      NA
    } else {
      x
    }
  }

  safe_as_numeric <- function(x) {
    x <- trimws(x %||% "")
    if (!nzchar(x)) {
      return(NA_real_)
    }
    suppressWarnings(as.numeric(x))
  }

  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }

  analysis_data <- reactiveValues(
    qf = NULL,
    Norm_method = default_config$Norm_method,
    results_folder = NULL,
    config = default_config,
    logs = character(),
    status = "Awaiting analysis",
    running = FALSE
  )

  append_log <- function(text) {
    analysis_data$logs <- c(analysis_data$logs, paste0(format(Sys.time(), "%H:%M:%S"), " - ", text))
  }

  observeEvent(input$reset_defaults, {
    updateTextInput(session, "project_name", value = default_config$Project_name)
    updateTextInput(session, "analysis_run", value = default_config$Analysis_run)
    updateCheckboxInput(session, "first_pass", value = default_config$First_pass)
    updateCheckboxInput(session, "run_normalizer", value = default_config$Run_normalizer)
    updateCheckboxInput(session, "run_dea", value = default_config$Run_DEA)
    updateCheckboxInput(session, "run_gsea", value = default_config$Run_GSEA)
    updateCheckboxInput(session, "run_pathfind", value = default_config$run_pathfind)
    updateCheckboxInput(session, "save_plots", value = default_config$save_plots)
    updateCheckboxInput(session, "export_data", value = default_config$export_data)
    updateCheckboxInput(session, "subtype", value = default_config$subtype)
    updateCheckboxInput(session, "run_interactive", value = default_config$Run_Interactive)
    updateTextInput(session, "treatment_group", value = default_config$Treatment_group)
    updateTextInput(session, "reference_group", value = default_config$Reference_group)
    updateNumericInput(session, "import_qval", value = default_config$Import_qval)
    updateNumericInput(session, "import_pg_qval", value = default_config$Import_pg_qval)
    updateNumericInput(session, "group_treshold", value = default_config$Group_treshold)
    updateNumericInput(session, "global_treshold", value = default_config$Global_treshold)
    updateSelectInput(session, "treshold_level", selected = default_config$Treshold_level)
    updateCheckboxInput(session, "peptide_level", value = default_config$peptide_level)
    updateSelectInput(session, "norm_method", selected = default_config$Norm_method)
    updateSelectInput(session, "batch_corr", selected = default_config$Batch_corr)
    updateSelectInput(session, "agg_method", selected = default_config$Agg_method)
    updateCheckboxInput(session, "impute", value = default_config$Impute)
    updateTextInput(session, "impute_mar", value = ifelse(is.na(default_config$impute_MAR), "", default_config$impute_MAR))
    updateTextInput(session, "impute_mnar", value = ifelse(is.na(default_config$impute_MNAR), "", default_config$impute_MNAR))
    updateTextInput(session, "impute_method", value = ifelse(is.na(default_config$impute_method), "", default_config$impute_method))
    updateSelectInput(session, "dea_method", selected = default_config$DEA_method)
    updateSelectInput(session, "study_design", selected = default_config$Study_design)
    updateNumericInput(session, "alpha", value = default_config$alpha)
    updateTextInput(session, "logfc_cutoff", value = ifelse(is.na(default_config$logFC_cutoff), "", default_config$logFC_cutoff))
    updateTextInput(session, "org_db", value = default_config$Org_db)
    updateTextInput(session, "boxplot_prot", value = paste(default_config$Boxplot_prot, collapse = ", "))
    updateTextInput(session, "agg_prot", value = paste(default_config$Agg_prot, collapse = ", "))
    updateTextAreaInput(session, "group_colors", value = paste(sprintf("%s = %s", names(default_config$Group_colors), default_config$Group_colors), collapse = "\n"))
    updateTextAreaInput(session, "expression_lvl_color", value = paste(sprintf("%s = %s", names(default_config$Expression_lvl_color), default_config$Expression_lvl_color), collapse = "\n"))
    updateTextAreaInput(session, "heatmap_annot", value = paste(sprintf("%s = %s", names(default_config$heatmap_annot), default_config$heatmap_annot), collapse = "\n"))
    updateTextInput(session, "dia_path", value = default_config$DIA_nn_path)
    updateTextInput(session, "fasta_path", value = default_config$FASTA_path)
    updateTextInput(session, "sample_path", value = default_config$Sample_data_path)
    updateTextInput(session, "your_name", value = default_config$Your_name)
    updateTextInput(session, "lab_name", value = default_config$Lab_name)
    updateTextInput(session, "contact", value = default_config$Contact)
  })

  observeEvent(input$run_pipeline, {
    if (analysis_data$running) {
      return()
    }

    config <- list(
      Project_name = input$project_name,
      Analysis_run = input$analysis_run,
      First_pass = isTRUE(input$first_pass),
      Run_normalizer = isTRUE(input$run_normalizer),
      Run_DEA = isTRUE(input$run_dea),
      Run_GSEA = isTRUE(input$run_gsea),
      run_pathfind = isTRUE(input$run_pathfind),
      save_plots = isTRUE(input$save_plots),
      export_data = isTRUE(input$export_data),
      subtype = isTRUE(input$subtype),
      Run_Interactive = FALSE,
      Treatment_group = input$treatment_group,
      Reference_group = input$reference_group,
      Import_qval = input$import_qval,
      Import_pg_qval = input$import_pg_qval,
      Group_treshold = input$group_treshold,
      Global_treshold = input$global_treshold,
      Treshold_level = input$treshold_level,
      peptide_level = isTRUE(input$peptide_level),
      Norm_method = input$norm_method,
      Batch_corr = input$batch_corr,
      Agg_method = input$agg_method,
      Impute = isTRUE(input$impute),
      impute_MAR = empty_to_na(input$impute_mar),
      impute_MNAR = empty_to_na(input$impute_mnar),
      impute_method = empty_to_na(input$impute_method),
      DEA_method = input$dea_method,
      Study_design = input$study_design,
      alpha = input$alpha,
      logFC_cutoff = safe_as_numeric(input$logfc_cutoff),
      Org_db = input$org_db,
      Boxplot_prot = parse_comma_list(input$boxplot_prot),
      Agg_prot = parse_comma_list(input$agg_prot),
      Group_colors = parse_named_lines(input$group_colors),
      heatmap_annot = parse_named_lines(input$heatmap_annot),
      Expression_lvl_color = parse_named_lines(input$expression_lvl_color),
      DIA_nn_path = input$dia_path,
      FASTA_path = input$fasta_path,
      Sample_data_path = input$sample_path,
      Your_name = input$your_name,
      Lab_name = input$lab_name,
      Contact = input$contact
    )

    if (length(config$Group_colors) == 0) {
      showNotification("Please provide at least one group color (format: Group = #HEX).", type = "error")
      return()
    }

    if (length(config$heatmap_annot) == 0) {
      showNotification("Please provide at least one heatmap annotation mapping.", type = "error")
      return()
    }

    analysis_data$running <- TRUE
    analysis_data$status <- "Running analysis..."
    analysis_data$logs <- character()
    append_log("Pipeline started")

    withProgress(message = "Executing pipeline", value = 0, {
      pipeline_env <- new.env(parent = globalenv())
      for (nm in names(config)) {
        assign(nm, config[[nm]], envir = pipeline_env)
      }

      options(pex_pipeline_dir = pipeline_dir)
      old_wd <- getwd()
      on.exit({
        setwd(old_wd)
        options(pex_pipeline_dir = NULL)
      }, add = TRUE)
      setwd(repo_root)

      result <- tryCatch({
        withCallingHandlers({
          sys.source(file.path(pipeline_dir, "1 - Initiation.R"), envir = pipeline_env)
        },
        message = function(m) {
          append_log(m$message)
          invokeRestart("muffleMessage")
        },
        warning = function(w) {
          append_log(paste0("Warning: ", w$message))
          invokeRestart("muffleWarning")
        })
        TRUE
      }, error = function(e) {
        append_log(paste0("Error: ", e$message))
        showNotification(paste("Pipeline failed:", e$message), type = "error")
        FALSE
      })

      if (isTRUE(result) && exists("qf", envir = pipeline_env)) {
        analysis_data$qf <- pipeline_env$qf
        analysis_data$Norm_method <- config$Norm_method
        analysis_data$results_folder <- pipeline_env$results_folder %||% NULL
        analysis_data$config <- config
        analysis_data$status <- "Analysis completed"
        append_log("Pipeline completed successfully")
      } else if (!isTRUE(result)) {
        analysis_data$status <- "Pipeline failed"
      } else {
        analysis_data$status <- "Pipeline finished without returning results"
        showNotification("Pipeline finished but no results were produced.", type = "warning")
      }
    })

    analysis_data$running <- FALSE
  })

  output$pipeline_status <- renderUI({
    status_text <- analysis_data$status
    status_class <- if (grepl("error|fail", tolower(status_text))) {
      "text-danger"
    } else if (grepl("complete", tolower(status_text))) {
      "text-success"
    } else {
      "text-muted"
    }
    div(class = status_class, strong(status_text))
  })

  output$pipeline_log <- renderText({
    if (length(analysis_data$logs) == 0) {
      "Logs will appear here when the pipeline runs."
    } else {
      paste(analysis_data$logs, collapse = "\n")
    }
  })

  output$norm_method_display <- renderText({
    analysis_data$Norm_method %||% "Not available"
  })

  output$bar_color_input <- renderUI({
    req(analysis_data$qf)
    df <- as.data.frame(colData(analysis_data$qf))
    selected <- if ("shortname" %in% names(df)) "shortname" else names(df)[1]
    varSelectInput(
      inputId = "bar_color",
      label = "Color by:",
      data = df,
      selected = selected
    )
  })

  output$bar_sort_input <- renderUI({
    req(analysis_data$qf)
    df <- as.data.frame(colData(analysis_data$qf))
    selected <- if ("Nr_prot" %in% names(df)) "Nr_prot" else names(df)[1]
    varSelectInput(
      inputId = "bar_sort",
      label = "Sort by (low to high):",
      data = df,
      selected = selected
    )
  })

  output$gene_selector <- renderUI({
    req(analysis_data$qf)
    req(!is.null(analysis_data$qf[["Results"]]))
    genes <- as.character(as.data.frame(rowData(analysis_data$qf[["Results"]]))[["Gene"]])
    genes <- genes[nzchar(genes)]
    selectInput(
      inputId = "gene",
      label = "Select gene:",
      choices = sort(unique(genes)),
      selected = if (length(genes)) genes[1] else NULL
    )
  })

  output$heatmap_annotations <- renderUI({
    req(analysis_data$qf)
    df <- as.data.frame(colData(analysis_data$qf))
    choices <- names(df)
    selected <- intersect(c("group", "Extraction_batch"), choices)
    if (length(selected) == 0) {
      selected <- head(choices, 2)
    }
    checkboxGroupInput(
      inputId = "selected_annotations",
      label = "Select annotation columns:",
      choices = choices,
      selected = selected
    )
  })


# Identifications ---------------------------------------------------------

  output$barplot <- renderPlot({
    req(analysis_data$qf, input$bar_color, input$bar_sort)
    p <- Barplot_IDs(analysis_data$qf, i = "Proteins", color = input$bar_color, sort = input$bar_sort)
    print(p)
  })

  output$download_barplot <- downloadHandler(
    filename = function() {
      req(analysis_data$results_folder)
      paste0(analysis_data$results_folder, "barplot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(analysis_data$qf)
      png(file, width = 920, height = 550)
      p <- Barplot_IDs(analysis_data$qf, i = "Proteins", color = input$bar_color, sort = input$bar_sort)
      print(p)
      dev.off()
    }
  )

# Missing values ----------------------------------------------------------

  output$missMap <- renderPlot({
    req(analysis_data$qf)
    p <- miss_map(analysis_data$qf, i = "Proteins")
    print(p)
  })


# Normalization ----------------------------------------------------------

  output$pre_norm_dens <- renderPlotly({
    req(analysis_data$qf)
    p <- Density_plotter(analysis_data$qf, i = "Peptides")
    print(p)
  })

  output$post_norm_dens <- renderPlotly({
    req(analysis_data$qf)
    p <- Density_plotter(analysis_data$qf, i = "PeptidesProccessed")
    print(p)
  })

  output$pre_norm_box <- renderPlot({
    req(analysis_data$qf)
    p <- Norm_boxplotter(analysis_data$qf, i = "Peptides")
    print(p)
  })

  output$post_norm_box <- renderPlot({
    req(analysis_data$qf)
    p <- Norm_boxplotter(analysis_data$qf, i = "PeptidesProccessed")
    print(p)
  })


# Protein -----------------------------------------------------------------

  output$aggregationPlot <- renderPlot({
    req(analysis_data$qf, input$gene)
    p <- Aggregation_plotter(analysis_data$qf, input$gene)
    print(p)
  })

  output$prot_boxplot <- renderPlot({
    req(analysis_data$qf, input$gene)
    p <- Boxplot_plotter(analysis_data$qf,input$gene)
    print(p)
  })

  output$volcanoPlotProt <- renderPlotly({
    req(analysis_data$qf, input$gene)
    p <- Volcano_plotter_highlight(analysis_data$qf, input$gene)
    ggplotly(p, tooltip = "text")
  })

  output$protein_info <- renderUI({
    req(analysis_data$qf, input$gene)
    req(!is.null(analysis_data$qf[["Results"]]))

    protein_data <- as.data.frame(rowData(analysis_data$qf[["Results"]]))
    selected_gene <- input$gene
    protein_row <- protein_data[protein_data$Gene == selected_gene, ]

    if (nrow(protein_row) == 0) {
      return(HTML("<p>No protein information available.</p>"))
    }

    ids <- protein_row$Protein.Group
    query <- list("accession_id" = ids)
    function_df <- query_uniprot(query,
                                 columns = c("accession", "cc_function"),
                                 show_progress = FALSE)

    function_text <- if (nrow(function_df) > 0) {
      function_df[["Function [CC]"]][1]
    } else {
      "Not available"
    }

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

  results_table_data <- reactive({
    req(analysis_data$qf)
    req(!is.null(analysis_data$qf[["Results"]]))
    as.data.frame(rowData(analysis_data$qf[["Results"]]))
  })

  output$results_table <- renderDT({
    df <- results_table_data()
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    numeric_cols <- numeric_cols[numeric_cols != "Number.Peptides"]

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
    req(analysis_data$qf)
    p <- Volcano_plotter(analysis_data$qf, Title = "Volcano Plot")
    print(p)
  })

  output$download_volcano <- downloadHandler(
    filename = function() {
      req(analysis_data$results_folder)
      paste0(analysis_data$results_folder, "volcano_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(analysis_data$qf)
      png(file, width = 443, height = 610)
      p <- Volcano_plotter(analysis_data$qf, Title = "Volcano Plot")
      print(p)
      dev.off()
    }
  )


# Heatmap -----------------------------------------------------------------

  annotation_list <- reactive({
    req(input$selected_annotations)
    setNames(as.list(input$selected_annotations), input$selected_annotations)
  })

  default_heatmap <- eventReactive(input$updateHeatmap, {
    req(analysis_data$qf)
    make_heatmap(analysis_data$qf,
                 i = "Results",
                 annot_list = annotation_list(),
                 only_significant = FALSE,
                 n = NULL)
  }, ignoreNULL = FALSE)

  output$defaultHeatmapPlot <- renderPlot({
    req(default_heatmap())
    default_heatmap()
  })

  custom_heatmap <- reactive({
    req(analysis_data$qf, input$n_value)
    make_heatmap(analysis_data$qf,
                 i = "Results",
                 annot_list = annotation_list(),
                 only_significant = TRUE,
                 n = input$n_value)
  })

  output$customHeatmapPlot <- renderPlot({
    req(custom_heatmap())
    custom_heatmap()
  })

  observe({
    req(analysis_data$qf)
    req(!is.null(analysis_data$qf[["Results"]]))
    slider_max <- nrow(rowData(analysis_data$qf[["Results"]]))
    updateSliderInput(
      session,
      "n_value",
      max = slider_max,
      value = min(50, slider_max)
    )
  })

  output$download_heatmap_full <- downloadHandler(
    filename = function() {
      paste0("Heatmap_full_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(default_heatmap())
      png(file, width = 600, height = 600)
      print(default_heatmap())
      dev.off()
    }
  )

  output$download_heatmap_sig <- downloadHandler(
    filename = function() {
      paste0("Heatmap_sig_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(custom_heatmap())
      png(file, width = 600, height = 600)
      print(custom_heatmap())
      dev.off()
    }
  )

# PCA ---------------------------------------------------------------------

  output$pcaPlot <- renderPlot({
    req(analysis_data$qf)
    p <- PCA_plotter(analysis_data$qf, i = "Proteins", color = "group")
    print(p)
  })

  output$download_pca <- downloadHandler(
    filename = function() {
      paste0("PCA", Sys.Date(), ".png")
    },
    content = function(file) {
      req(analysis_data$qf)
      png(file, width = 530, height = 360)
      p <- PCA_plotter(analysis_data$qf, i = "Proteins", color = "group", label = "none")
      print(p)
      dev.off()
    }
  )

}
