##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

# Server ------------------------------------------------------------------

function(input, output, session) {

  repo_root <- get0(
    "repo_root",
    ifnotfound = normalizePath("..", winslash = "/", mustWork = FALSE),
    inherits = TRUE
  )
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = FALSE)
  if (!dir.exists(repo_root)) {
    repo_root <- normalizePath(".", winslash = "/", mustWork = FALSE)
  }
  if (!dir.exists(file.path(repo_root, "R"))) {
    stop(
      "The Proteome Explorer 'R' directory could not be found at ",
      repo_root,
      ". Launch the app from the repository root so the shared scripts are available."
    )
  }

  pipeline_dir <- get0(
    "pipeline_dir",
    ifnotfound = file.path(repo_root, "pipeline"),
    inherits = TRUE
  )
  pipeline_dir <- normalizePath(pipeline_dir, winslash = "/", mustWork = FALSE)

  repo_file <- function(...) {
    path <- file.path(repo_root, ...)
    if (!file.exists(path)) {
      stop(
        sprintf(
          "Expected to find '%s' inside the repository root (%s) but it was missing.",
          file.path(...),
          repo_root
        )
      )
    }
    path
  }

  sys.source(repo_file("R", "targets_config.R"), envir = environment())
  sys.source(repo_file("R", "targets_pipeline.R"), envir = environment())

  default_config <- list(
    Project_name = "PExA Lungcancer v6",
    Analysis_run = "EBP_v1_detect",
    First_pass = TRUE,
    Run_normalizer = TRUE,
    Run_DEA = TRUE,
    Run_GSEA = TRUE,
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

  build_status_badge <- function(text, variant) {
    badge_fn <- NULL
    if (exists("bs_badge", mode = "function")) {
      badge_fn <- get("bs_badge", mode = "function")
    } else if (exists("bs_badge", envir = asNamespace("bslib"), inherits = FALSE)) {
      badge_fn <- get("bs_badge", envir = asNamespace("bslib"))
    }

    classes <- paste0("bg-", variant)

    if (!is.null(badge_fn)) {
      args <- list(text, class = classes)
      if ("bs_theme" %in% names(formals(badge_fn))) {
        args$bs_theme <- bslib::bs_theme()
      }
      return(do.call(badge_fn, args))
    }

    htmltools::span(text, class = paste("badge", classes))
  }

  resolve_path <- function(path) {
    path <- trimws(path %||% "")
    if (!nzchar(path)) {
      return("")
    }
    expanded <- path.expand(path)
    if (.Platform$OS.type == "windows") {
      is_abs <- grepl("^[A-Za-z]:[\\/]|^\\\\", expanded)
    } else {
      is_abs <- substr(expanded, 1, 1) == "/"
    }
    target <- if (is_abs) expanded else file.path(repo_root, expanded)
    normalizePath(target, winslash = "/", mustWork = FALSE)
  }

  analysis_data <- reactiveValues(
    qf = NULL,
    Norm_method = default_config$Norm_method,
    results_folder = NULL,
    config = default_config,
    logs = character(),
    status = "Awaiting analysis",
    running = FALSE,
    plot_functions = list(),
    pipeline_graph = NULL
  )

  flush_pipeline_outputs <- function(scroll = FALSE) {
    try(shiny::flushReact(), silent = TRUE)
    if (isTRUE(scroll)) {
      session$sendCustomMessage("pex-scroll-log", list())
    }
  }

  set_pipeline_status <- function(text, immediate = TRUE) {
    analysis_data$status <- text
    if (isTRUE(immediate)) {
      flush_pipeline_outputs()
    }
  }

  append_log <- function(text) {
    if (length(text) == 0) {
      return(invisible(NULL))
    }
    text <- text[!is.na(text)]
    text <- text[nzchar(text)]
    if (!length(text)) {
      return(invisible(NULL))
    }
    timestamp <- format(Sys.time(), "%H:%M:%S")
    analysis_data$logs <- c(analysis_data$logs, paste(timestamp, "-", text))
    flush_pipeline_outputs(scroll = TRUE)
    invisible(NULL)
  }

  root_volume <- if (.Platform$OS.type == "windows") "C:/" else "/"
  volumes <- c(
    "Workspace" = repo_root,
    "Pipeline" = pipeline_dir,
    "Home" = path.expand("~"),
    "Root" = root_volume
  )
  volumes <- volumes[dir.exists(volumes)]
  if (!length(volumes)) {
    volumes <- c("Workspace" = repo_root)
  }

  shinyFiles::shinyFileChoose(
    input,
    id = "dia_path_browse",
    session = session,
    roots = volumes,
    filetypes = c("", "tsv", "csv", "txt")
  )
  shinyFiles::shinyFileChoose(
    input,
    id = "fasta_path_browse",
    session = session,
    roots = volumes,
    filetypes = c("", "fasta", "fa")
  )
  shinyFiles::shinyFileChoose(
    input,
    id = "sample_path_browse",
    session = session,
    roots = volumes,
    filetypes = c("", "xlsx", "xls", "csv")
  )

  observeEvent(input$dia_path_browse, {
    files <- shinyFiles::parseFilePaths(volumes, input$dia_path_browse)
    if (nrow(files) > 0) {
      path <- normalizePath(files$datapath[1], winslash = "/", mustWork = FALSE)
      updateTextInput(session, "dia_path", value = path)
    }
  })

  observeEvent(input$fasta_path_browse, {
    files <- shinyFiles::parseFilePaths(volumes, input$fasta_path_browse)
    if (nrow(files) > 0) {
      path <- normalizePath(files$datapath[1], winslash = "/", mustWork = FALSE)
      updateTextInput(session, "fasta_path", value = path)
    }
  })

  observeEvent(input$sample_path_browse, {
    files <- shinyFiles::parseFilePaths(volumes, input$sample_path_browse)
    if (nrow(files) > 0) {
      path <- normalizePath(files$datapath[1], winslash = "/", mustWork = FALSE)
      updateTextInput(session, "sample_path", value = path)
    }
  })

  observeEvent(input$reset_defaults, {
    updateTextInput(session, "project_name", value = default_config$Project_name)
    updateTextInput(session, "analysis_run", value = default_config$Analysis_run)
    updateCheckboxInput(session, "first_pass", value = default_config$First_pass)
    updateCheckboxInput(session, "run_normalizer", value = default_config$Run_normalizer)
    updateCheckboxInput(session, "run_dea", value = default_config$Run_DEA)
    updateCheckboxInput(session, "run_gsea", value = default_config$Run_GSEA)
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

  output$pipeline_status <- renderUI({
    status <- analysis_data$status
    badge_status <- dplyr::case_when(
      status == "Awaiting analysis" ~ "secondary",
      status == "Running analysis..." ~ "primary",
      startsWith(status, "Running ") ~ "primary",
      status == "Analysis completed" ~ "success",
      status == "Pipeline failed" ~ "danger",
      status == "Pipeline finished without returning results" ~ "warning",
      status == "File validation failed" ~ "danger",
      status == "Pipeline finished but no results were produced." ~ "warning",
      TRUE ~ "secondary"
    )
    tagList(
      strong("Status:"),
      build_status_badge(status, badge_status)
    )
  })

  output$pipeline_log <- renderText({
    if (!length(analysis_data$logs)) {
      "Pipeline logs will appear here."
    } else {
      paste(analysis_data$logs, collapse = "\n")
    }
  })

  observeEvent(input$run_pipeline, {

    if (analysis_data$running) {
      showNotification("Pipeline is already running.", type = "message")
      return()
    }

    config <- list(
      Project_name = input$project_name,
      Analysis_run = input$analysis_run,
      First_pass = input$first_pass,
      Run_normalizer = input$run_normalizer,
      Run_DEA = input$run_dea,
      Run_GSEA = input$run_gsea,
      save_plots = input$save_plots,
      export_data = input$export_data,
      subtype = input$subtype,
      Run_Interactive = input$run_interactive,
      Treatment_group = input$treatment_group,
      Reference_group = input$reference_group,
      Import_qval = input$import_qval,
      Import_pg_qval = input$import_pg_qval,
      Group_treshold = input$group_treshold,
      Global_treshold = input$global_treshold,
      Treshold_level = input$treshold_level,
      peptide_level = input$peptide_level,
      Norm_method = input$norm_method,
      Batch_corr = input$batch_corr,
      Agg_method = input$agg_method,
      Impute = input$impute,
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
      DIA_nn_path = resolve_path(input$dia_path),
      FASTA_path = resolve_path(input$fasta_path),
      Sample_data_path = resolve_path(input$sample_path),
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

    required_paths <- c(
      "DIA-NN report" = config$DIA_nn_path,
      "FASTA file" = config$FASTA_path,
      "Sample metadata" = config$Sample_data_path
    )

    missing_paths <- names(required_paths)[!nzchar(unname(required_paths)) | !file.exists(unname(required_paths))]
    if (length(missing_paths) > 0) {
      msg <- paste0("Missing required file(s): ", paste(missing_paths, collapse = ", "), ".")
      showNotification(msg, type = "error")
      set_pipeline_status("File validation failed")
      append_log(msg)
      return()
    }

    step_labels <- c(
      import = "Import data",
      preprocess = "Data preprocessing",
      dea = "Differential expression analysis",
      postprocess = "Post-processing",
      pathway = "Pathway analysis",
      export = "Export and visualization",
      report = "Generate report"
    )

    analysis_data$plot_functions <- list()
    analysis_data$running <- TRUE
    set_pipeline_status("Running analysis...")
    analysis_data$logs <- character()
    append_log("Pipeline started")

    analysis_data$pipeline_graph <- tryCatch({
      targets::tar_visnetwork(callr_function = NULL)
    }, error = function(e) {
      append_log(paste0("Warning: unable to render targets graph - ", e$message))
      NULL
    })

    withProgress(message = "Executing pipeline", value = 0, {
      captured_output <- character()
      prepared_config <- NULL
      pipeline_state <- NULL
      setProgress(value = 0, detail = "Preparing configuration...")

      result <- tryCatch({
        local_config <- prepare_pipeline_config(config = config, repo_root = repo_root)
        local_state <- NULL

        captured_output <<- capture.output({
          withCallingHandlers({
            local_state <<- run_full_pipeline(
              local_config,
              callbacks = list(
                on_step_start = function(step, index, total) {
                  label <- step_labels[[step]] %||% step
                  append_log(sprintf("Starting %s (%d/%d)", label, index, total))
                  set_pipeline_status(sprintf("Running %s (%d/%d)", label, index, total))
                  setProgress(
                    value = (index - 1) / total,
                    detail = sprintf("Running %s...", label)
                  )
                },
                on_step_complete = function(step, index, total) {
                  label <- step_labels[[step]] %||% step
                  append_log(sprintf("Completed %s (%d/%d)", label, index, total))
                  setProgress(
                    value = index / total,
                    detail = sprintf("Completed %s", label)
                  )
                }
              )
            )
          },
          message = function(m) {
            append_log(m$message)
            invokeRestart("muffleMessage")
          },
          warning = function(w) {
            append_log(paste0("Warning: ", w$message))
            invokeRestart("muffleWarning")
          })
        }, type = "output")

        prepared_config <<- local_config
        pipeline_state <<- local_state
        TRUE
      }, error = function(e) {
        append_log(paste0("Error: ", e$message))
        setProgress(value = 1, detail = "Pipeline failed")
        showNotification(paste("Pipeline failed:", e$message), type = "error")
        FALSE
      })

      if (length(captured_output)) {
        append_log(captured_output)
      }

      if (isTRUE(result) && !is.null(pipeline_state$qf)) {
        analysis_data$qf <- pipeline_state$qf
        analysis_data$Norm_method <- prepared_config$Norm_method
        analysis_data$results_folder <- prepared_config$results_folder
        analysis_data$config <- prepared_config
        set_pipeline_status("Analysis completed")
        setProgress(value = 1, detail = "Pipeline completed successfully")
        append_log("Pipeline completed successfully")

        required_functions <- c(
          "Barplot_IDs",
          "miss_map",
          "Density_plotter",
          "Norm_boxplotter",
          "Aggregation_plotter",
          "Boxplot_plotter",
          "Volcano_plotter",
          "Volcano_plotter_highlight",
          "make_heatmap",
          "PCA_plotter"
        )

        analysis_data$plot_functions <- lapply(setNames(required_functions, required_functions), function(fn) {
          pipeline_state[[fn]] %||% NULL
        })
      } else if (!isTRUE(result)) {
        set_pipeline_status("Pipeline failed")
        setProgress(value = 1, detail = "Pipeline failed")
      } else {
        set_pipeline_status("Pipeline finished without returning results")
        showNotification("Pipeline finished but no results were produced.", type = "warning")
        setProgress(value = 1, detail = "Pipeline finished without results")
      }
    })

    analysis_data$running <- FALSE

    if (analysis_data$status == "Analysis completed" && isTRUE(config$Run_Interactive)) {
      append_log("Launching interactive dashboard...")
      shiny::runApp(file.path(repo_root, "app"))
    }
  })

  output$pipeline_network <- visNetwork::renderVisNetwork({
    graph <- analysis_data$pipeline_graph
    if (is.null(graph)) {
      nodes <- data.frame(
        id = "awaiting_run",
        label = "Launch the pipeline to view the targets graph.",
        stringsAsFactors = FALSE
      )
      edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)
      return(visNetwork::visNetwork(nodes, edges))
    }
    graph
  })

  observeEvent(analysis_data$plot_functions, {
    if (!length(analysis_data$plot_functions)) {
      return()
    }

    output$norm_plots <- renderUI({
      if (!length(analysis_data$plot_functions)) {
        return(helpText("Run the pipeline to see normalisation diagnostics."))
      }
      norm_fns <- c("Density_plotter", "Norm_boxplotter")
      missing <- setdiff(norm_fns, names(analysis_data$plot_functions))
      if (length(missing)) {
        return(helpText(paste("Missing plotting functions:", paste(missing, collapse = ", "))))
      }

      tagList(norm_cards)
    })

    output$agg_plot <- renderPlot({
      fn <- analysis_data$plot_functions$Aggregation_plotter
      if (is.null(fn)) {
        return(NULL)
      }
      fn()
    })

    output$volcano_plot <- renderPlot({
      fn <- analysis_data$plot_functions$Volcano_plotter
      if (is.null(fn)) {
        return(NULL)
      }
      fn()
    })

    output$volcano_plot_highlight <- renderPlot({
      fn <- analysis_data$plot_functions$Volcano_plotter_highlight
      if (is.null(fn)) {
        return(NULL)
      }
      fn()
    })

    output$heatmap <- renderPlot({
      fn <- analysis_data$plot_functions$make_heatmap
      if (is.null(fn)) {
        return(NULL)
      }
      fn()
    })

    output$pca_plot <- renderPlot({
      fn <- analysis_data$plot_functions$PCA_plotter
      if (is.null(fn)) {
        return(NULL)
      }
      fn()
    })

    output$pre_norm_dens <- renderPlotly({
      fn <- analysis_data$plot_functions$Density_plotter
      if (is.null(fn)) {
        return(NULL)
      }
      fn(type = "raw")
    })

    output$post_norm_dens <- renderPlotly({
      fn <- analysis_data$plot_functions$Density_plotter
      if (is.null(fn)) {
        return(NULL)
      }
      fn(type = "normalized")
    })

    output$pre_norm_box <- renderPlot({
      fn <- analysis_data$plot_functions$Norm_boxplotter
      if (is.null(fn)) {
        return(NULL)
      }
      fn(type = "raw")
    })

    output$post_norm_box <- renderPlot({
      fn <- analysis_data$plot_functions$Norm_boxplotter
      if (is.null(fn)) {
        return(NULL)
      }
      fn(type = "normalized")
    })

    output$boxplot_plot <- renderPlot({
      fn <- analysis_data$plot_functions$Boxplot_plotter
      if (is.null(fn)) {
        return(NULL)
      }
      fn()
    })

    output$barplot_ids <- renderPlot({
      fn <- analysis_data$plot_functions$Barplot_IDs
      if (is.null(fn)) {
        return(NULL)
      }
      fn()
    })

    output$missingness_heatmap <- renderPlot({
      fn <- analysis_data$plot_functions$miss_map
      if (is.null(fn)) {
        return(NULL)
      }
      fn()
    })
  })

  output$analysis_summary <- renderDT({
    if (is.null(analysis_data$qf)) {
      return(NULL)
    }
    qf <- analysis_data$qf

    norm_method <- analysis_data$Norm_method

    stats <- data.frame(
      Metric = c(
        "Number of assays",
        "Number of features",
        "Number of proteins",
        "Normalization method"
      ),
      Value = c(
        length(BiocParallel::bpparam("SerialParam")),
        nrow(qf),
        length(unique(qf$protein_group)),
        norm_method
      ),
      stringsAsFactors = FALSE
    )

    datatable(
      stats,
      options = list(pageLength = 5, dom = "tip"),
      rownames = FALSE,
      caption = "Analysis summary"
    )
  })

}
