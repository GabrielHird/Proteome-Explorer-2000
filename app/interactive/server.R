##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

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
      "The Proteome Explorer repository root could not be located. ",
      "Launch the app with runApp('app') from the project root."
    )
  }

  sys.source(file.path(repo_root, "R", "targets_config.R"), envir = environment())
  sys.source(file.path(repo_root, "R", "targets_pipeline.R"), envir = environment())
  sys.source(file.path(repo_root, "R", "app_targets.R"), envir = environment())

  default_config <- pipeline_default_config()
  runtime_config <- reactiveVal(app_read_runtime_config(repo_root))
  run_log <- reactiveVal(character())
  meta_data <- reactiveVal(data.frame())
  results_listing <- reactiveVal(data.frame())
  outdated_targets <- reactiveVal(character())

  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }

  to_named_string <- function(x) {
    if (is.null(x)) {
      return("")
    }
    if (length(x) == 0) {
      return("")
    }
    if (is.null(names(x)) || all(names(x) == "")) {
      return(paste(x, collapse = ", "))
    }
    paste(sprintf("%s = %s", names(x), x), collapse = "\n")
  }

  parse_comma_list <- function(x) {
    vals <- trimws(unlist(strsplit(x %||% "", ",")))
    vals <- vals[nzchar(vals)]
    if (length(vals)) vals else character()
  }

  parse_named_lines <- function(x) {
    lines <- strsplit(x %||% "", "\n", fixed = TRUE)[[1]]
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
    if (!length(lines)) {
      return(character())
    }
    split <- strsplit(lines, "=", fixed = TRUE)
    values <- vapply(split, function(parts) {
      parts <- trimws(parts)
      if (length(parts) < 2) "" else parts[2]
    }, character(1))
    names(values) <- vapply(split, function(parts) {
      parts <- trimws(parts)
      parts[1]
    }, character(1))
    values <- values[nzchar(values) & nzchar(names(values))]
    values
  }

  empty_to_na <- function(x) {
    x <- trimws(x %||% "")
    if (!nzchar(x)) {
      return(NA)
    }
    x
  }

  safe_numeric <- function(x) {
    if (is.null(x)) {
      return(NA_real_)
    }
    if (is.numeric(x)) {
      return(as.numeric(x))
    }
    x <- trimws(as.character(x))
    if (!nzchar(x)) {
      return(NA_real_)
    }
    suppressWarnings(as.numeric(x))
  }

  update_inputs_from_config <- function(cfg) {
    updateTextInput(session, "project_name", value = cfg$Project_name)
    updateTextInput(session, "analysis_run", value = cfg$Analysis_run)
    updateTextInput(session, "dia_nn_path", value = cfg$DIA_nn_path)
    updateTextInput(session, "fasta_path", value = cfg$FASTA_path)
    updateTextInput(session, "sample_data_path", value = cfg$Sample_data_path)
    updateNumericInput(session, "import_qval", value = cfg$Import_qval)
    updateNumericInput(session, "import_pg_qval", value = cfg$Import_pg_qval)
    updateNumericInput(session, "group_threshold", value = cfg$Group_treshold)
    updateNumericInput(session, "global_threshold", value = cfg$Global_treshold)
    updateSelectInput(session, "threshold_level", selected = cfg$Treshold_level)
    updateSelectInput(session, "agg_method", selected = cfg$Agg_method)
    updateCheckboxInput(session, "peptide_level", value = isTRUE(cfg$peptide_level))
    updateCheckboxInput(session, "run_normalizer", value = isTRUE(cfg$Run_normalizer))
    updateSelectInput(session, "norm_method", selected = cfg$Norm_method)
    updateSelectInput(session, "batch_corr", selected = cfg$Batch_corr)
    updateCheckboxInput(session, "impute", value = isTRUE(cfg$Impute))
    updateTextInput(session, "impute_method", value = cfg$impute_method %||% "")
    updateNumericInput(session, "impute_mar", value = cfg$impute_MAR)
    updateNumericInput(session, "impute_mnar", value = cfg$impute_MNAR)
    updateCheckboxInput(session, "run_dea", value = isTRUE(cfg$Run_DEA))
    updateSelectInput(session, "dea_method", selected = cfg$DEA_method)
    updateSelectInput(session, "study_design", selected = cfg$Study_design)
    updateNumericInput(session, "alpha", value = cfg$alpha)
    updateNumericInput(session, "logfc_cutoff", value = cfg$logFC_cutoff)
    updateCheckboxInput(session, "run_gsea", value = isTRUE(cfg$Run_GSEA))
    updateTextInput(session, "org_db", value = cfg$Org_db)
    updateCheckboxInput(session, "save_plots", value = isTRUE(cfg$save_plots))
    updateCheckboxInput(session, "export_data", value = isTRUE(cfg$export_data))
    updateTextInput(session, "treatment_group", value = cfg$Treatment_group)
    updateTextInput(session, "reference_group", value = cfg$Reference_group)
    updateTextAreaInput(session, "boxplot_prot", value = paste(cfg$Boxplot_prot, collapse = ", "))
    updateTextAreaInput(session, "agg_prot", value = paste(cfg$Agg_prot, collapse = ", "))
    updateTextAreaInput(session, "group_colors", value = to_named_string(cfg$Group_colors))
    updateTextAreaInput(session, "heatmap_annotations", value = to_named_string(cfg$heatmap_annot))
    updateTextAreaInput(session, "expression_levels", value = to_named_string(cfg$Expression_lvl_color))
    updateTextInput(session, "your_name", value = cfg$Your_name)
    updateTextInput(session, "lab_name", value = cfg$Lab_name)
    updateTextInput(session, "contact", value = cfg$Contact)
  }

  observe({
    cfg <- runtime_config()
    update_inputs_from_config(cfg)
    outdated_targets(tryCatch(app_targets_outdated(repo_root), error = function(e) character()))
    meta_data(tryCatch(app_targets_meta(repo_root), error = function(e) data.frame()))
    res_info <- app_collect_results(cfg, repo_root)
    files <- res_info$files
    if (!length(files)) {
      results_listing(data.frame())
    } else {
      rel_paths <- files
      repo_prefix <- paste0(repo_root, "/")
      rel_paths[startsWith(rel_paths, repo_prefix)] <- substring(rel_paths[startsWith(rel_paths, repo_prefix)], nchar(repo_prefix) + 1)
      results_listing(data.frame(
        name = basename(files),
        path = rel_paths,
        stringsAsFactors = FALSE
      ))
    }
  })

  observeEvent(input$reset_defaults, {
    runtime_config(default_config)
    update_inputs_from_config(default_config)
    showNotification("Settings reset to defaults.", type = "message")
  })

  observeEvent(input$run_pipeline, {
    overrides <- list(
      Project_name = trimws(input$project_name %||% default_config$Project_name),
      Analysis_run = trimws(input$analysis_run %||% default_config$Analysis_run),
      DIA_nn_path = trimws(input$dia_nn_path %||% default_config$DIA_nn_path),
      FASTA_path = trimws(input$fasta_path %||% default_config$FASTA_path),
      Sample_data_path = trimws(input$sample_data_path %||% default_config$Sample_data_path),
      Import_qval = safe_numeric(input$import_qval),
      Import_pg_qval = safe_numeric(input$import_pg_qval),
      Group_treshold = safe_numeric(input$group_threshold),
      Global_treshold = safe_numeric(input$global_threshold),
      Treshold_level = input$threshold_level %||% default_config$Treshold_level,
      Agg_method = input$agg_method %||% default_config$Agg_method,
      peptide_level = isTRUE(input$peptide_level),
      Run_normalizer = isTRUE(input$run_normalizer),
      Norm_method = input$norm_method %||% default_config$Norm_method,
      Batch_corr = input$batch_corr %||% default_config$Batch_corr,
      Impute = isTRUE(input$impute),
      impute_method = empty_to_na(input$impute_method),
      impute_MAR = safe_numeric(input$impute_mar),
      impute_MNAR = safe_numeric(input$impute_mnar),
      Run_DEA = isTRUE(input$run_dea),
      DEA_method = input$dea_method %||% default_config$DEA_method,
      Study_design = input$study_design %||% default_config$Study_design,
      alpha = safe_numeric(input$alpha),
      logFC_cutoff = safe_numeric(input$logfc_cutoff),
      Run_GSEA = isTRUE(input$run_gsea),
      Org_db = trimws(input$org_db %||% default_config$Org_db),
      save_plots = isTRUE(input$save_plots),
      export_data = isTRUE(input$export_data),
      Treatment_group = trimws(input$treatment_group %||% default_config$Treatment_group),
      Reference_group = trimws(input$reference_group %||% default_config$Reference_group),
      Boxplot_prot = parse_comma_list(input$boxplot_prot),
      Agg_prot = parse_comma_list(input$agg_prot),
      Group_colors = parse_named_lines(input$group_colors),
      heatmap_annot = parse_named_lines(input$heatmap_annotations),
      Expression_lvl_color = parse_named_lines(input$expression_levels),
      Your_name = trimws(input$your_name %||% default_config$Your_name),
      Lab_name = trimws(input$lab_name %||% default_config$Lab_name),
      Contact = trimws(input$contact %||% default_config$Contact)
    )

    combined <- modifyList(default_config, overrides, keep.null = TRUE)
    runtime_config(combined)
    app_write_runtime_config(combined, repo_root)

    showNotification("Starting {targets} pipeline run...", type = "message")

    log_messages <- character()
    result <- tryCatch({
      withCallingHandlers({
        app_run_targets_pipeline(repo_root)
      },
      message = function(m) {
        log_messages <<- c(log_messages, m$message)
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        log_messages <<- c(log_messages, paste("Warning:", w$message))
        invokeRestart("muffleWarning")
      })
      TRUE
    }, error = function(e) {
      log_messages <<- c(log_messages, paste("Error:", e$message))
      showNotification(paste("Pipeline run failed:", e$message), type = "error")
      FALSE
    })

    run_log(log_messages)

    if (isTRUE(result)) {
      showNotification("Pipeline run completed.", type = "message")
    }

    outdated_targets(tryCatch(app_targets_outdated(repo_root), error = function(e) character()))
    meta_data(tryCatch(app_targets_meta(repo_root), error = function(e) data.frame()))
    cfg <- runtime_config()
    res_info <- app_collect_results(cfg, repo_root)
    files <- res_info$files
    if (!length(files)) {
      results_listing(data.frame())
    } else {
      rel_paths <- files
      repo_prefix <- paste0(repo_root, "/")
      rel_paths[startsWith(rel_paths, repo_prefix)] <- substring(rel_paths[startsWith(rel_paths, repo_prefix)], nchar(repo_prefix) + 1)
      results_listing(data.frame(
        name = basename(files),
        path = rel_paths,
        stringsAsFactors = FALSE
      ))
    }
  })

  output$outdated_table <- DT::renderDT({
    targets <- outdated_targets()
    if (!length(targets)) {
      return(DT::datatable(data.frame(Status = "All targets up to date"), options = list(dom = "t"), rownames = FALSE))
    }
    DT::datatable(
      data.frame(Target = targets, stringsAsFactors = FALSE),
      options = list(pageLength = 10),
      rownames = FALSE
    )
  })

  output$meta_table <- DT::renderDT({
    meta <- meta_data()
    if (!nrow(meta)) {
      return(DT::datatable(data.frame(Message = "No completed runs yet."), options = list(dom = "t"), rownames = FALSE))
    }
    meta$seconds <- round(meta$seconds, 2)
    meta$warnings <- vapply(meta$warnings, function(x) {
      if (is.null(x) || !length(x)) "" else paste(x, collapse = "\n")
    }, character(1))
    meta$error <- vapply(meta$error, function(x) {
      if (is.null(x) || !length(x)) "" else paste(x, collapse = "\n")
    }, character(1))
    meta$started <- as.character(meta$started)
    meta$finished <- as.character(meta$finished)
    DT::datatable(meta, options = list(pageLength = 10))
  })

  output$run_log <- renderText({
    log <- run_log()
    if (!length(log)) {
      "Pipeline output will appear here after the next run."
    } else {
      paste(log, collapse = "\n")
    }
  })

  output$report_link <- renderUI({
    cfg <- runtime_config()
    res_info <- app_collect_results(cfg, repo_root)
    report_file <- tryCatch({
      old_dir <- setwd(repo_root)
      on.exit(setwd(old_dir))
      targets::tar_read(report_file, store = app_targets_store(repo_root))
    }, error = function(e) NA_character_)
    if (is.null(report_file) || (length(report_file) && is.na(report_file))) {
      report_file <- res_info$report_file
    }
    if (!is.na(report_file) && file.exists(report_file)) {
      rel <- report_file
      repo_prefix <- paste0(repo_root, "/")
      if (startsWith(rel, repo_prefix)) {
        rel <- substring(rel, nchar(repo_prefix) + 1)
      }
      tags$p(tags$a(href = rel, rel = "noopener", target = "_blank", "Open generated report"))
    } else {
      tags$p("Report has not been generated yet.")
    }
  })

  output$results_table <- DT::renderDT({
    listing <- results_listing()
    if (!nrow(listing)) {
      return(DT::datatable(data.frame(Message = "No exported files detected."), options = list(dom = "t"), rownames = FALSE))
    }
    DT::datatable(listing, options = list(pageLength = 15))
  })
}
