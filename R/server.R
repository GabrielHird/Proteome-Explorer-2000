server <- function(input, output, session) {
  transient_alert()
  repo_root <- project_repo_root()
  defaults <- project_default_config()

  runtime_config <- reactiveVal(project_runtime_config(repo_root))
  run_log <- reactiveVal(character())
  meta_data <- reactiveVal(data.frame())
  progress_data <- reactiveVal(data.frame())
  results_table_data <- reactiveVal(data.frame())
  outdated_targets <- reactiveVal(character())

  control_set()

  refresh_results <- function(cfg) {
    outdated_targets(results_outdated(repo_root))
    meta_data(results_meta(repo_root))
    progress_data(results_progress(repo_root))
    listing <- results_collect(cfg, repo_root)
    files <- listing$files
    if (!length(files)) {
      results_table_data(data.frame())
    } else {
      rel_paths <- files
      repo_prefix <- paste0(repo_root, "/")
      rel_paths[startsWith(rel_paths, repo_prefix)] <- substring(rel_paths[startsWith(rel_paths, repo_prefix)], nchar(repo_prefix) + 1)
      results_table_data(data.frame(
        name = basename(files),
        path = rel_paths,
        stringsAsFactors = FALSE
      ))
    }
  }

  observe({
    cfg <- runtime_config()
    project_update_inputs(session, cfg)
    refresh_results(cfg)
  })

  observeEvent(input$reset_defaults, {
    runtime_config(defaults)
    project_write_config(defaults, repo_root)
    run_log(character())
    refresh_results(defaults)
    showNotification("Settings reset to defaults.", type = "message")
    control_set()
  })

  observeEvent(input$run_pipeline, {
    overrides <- project_collect_overrides(input, defaults)
    combined <- project_merge_config(defaults, overrides)
    runtime_config(combined)
    project_write_config(combined, repo_root)

    showNotification("Starting {targets} pipeline run...", type = "message")
    process_set_running(TRUE)
    control_set()

    result <- process_run_pipeline(repo_root)
    run_log(result$log)

    if (isTRUE(result$success)) {
      showNotification("Pipeline run completed.", type = "message")
    } else {
      showNotification(paste("Pipeline run failed:", result$message), type = "error")
    }

    process_set_running(FALSE)
    control_set()
    refresh_results(runtime_config())
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

  output$progress_table <- DT::renderDT({
    progress <- progress_data()
    if (!nrow(progress)) {
      progress <- data.frame(Message = "No progress information available.")
      return(DT::datatable(progress, options = list(dom = "t"), rownames = FALSE))
    }
    progress$time <- round(progress$time, 2)
    progress$progress <- round(progress$progress, 2)
    DT::datatable(progress, options = list(pageLength = 10))
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
    log_text(run_log())
  })

  output$report_link <- renderUI({
    cfg <- runtime_config()
    report_file <- results_report_path(repo_root, cfg)
    if (!is.null(report_file) && length(report_file) && !is.na(report_file) && file.exists(report_file)) {
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
    listing <- results_table_data()
    if (!nrow(listing)) {
      return(DT::datatable(data.frame(Message = "No exported files detected."), options = list(dom = "t"), rownames = FALSE))
    }
    DT::datatable(listing, options = list(pageLength = 15))
  })
}
