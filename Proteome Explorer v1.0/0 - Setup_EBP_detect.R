##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

suppressPackageStartupMessages({
  library(shiny)
})

pipeline_dir <- normalizePath(file.path("Proteome Explorer v1.0", "ProteomeExplorer Pipeline"), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

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

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

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

landing_panel <- div(
  class = "text-center",
  style = "margin-top: 80px;",
  h1("Proteome Explorer Pipeline"),
  p(class = "lead", "Launch the analysis pipeline with your dataset and configuration."),
  actionButton("start_setup", "Launch Pipeline Setup", class = "btn btn-primary btn-lg")
)

configuration_panel <- fluidRow(
  column(
    width = 4,
    wellPanel(
      h3("Pipeline control"),
      actionButton("run_pipeline", "Run Pipeline", class = "btn btn-success btn-lg btn-block"),
      br(),
      actionButton("reset_defaults", "Reset defaults", class = "btn btn-link"),
      hr(),
      uiOutput("pipeline_status"),
      hr(),
      strong("Log"),
      verbatimTextOutput("pipeline_log", placeholder = TRUE)
    )
  ),
  column(
    width = 8,
    tabsetPanel(
      id = "config_tabs",
      tabPanel(
        "Project",
        textInput("project_name", "Project name", value = default_config$Project_name),
        textInput("analysis_run", "Analysis run", value = default_config$Analysis_run),
        textInput("your_name", "Your name", value = default_config$Your_name),
        textInput("lab_name", "Lab name", value = default_config$Lab_name),
        textInput("contact", "Contact", value = default_config$Contact)
      ),
      tabPanel(
        "Data",
        textInput("dia_path", "DIA-NN report", value = default_config$DIA_nn_path),
        textInput("fasta_path", "FASTA file", value = default_config$FASTA_path),
        textInput("sample_path", "Sample metadata", value = default_config$Sample_data_path)
      ),
      tabPanel(
        "Run options",
        checkboxInput("first_pass", "First pass", value = default_config$First_pass),
        checkboxInput("run_normalizer", "Run normalizer", value = default_config$Run_normalizer),
        checkboxInput("run_dea", "Run DEA", value = default_config$Run_DEA),
        checkboxInput("run_gsea", "Run GSEA", value = default_config$Run_GSEA),
        checkboxInput("run_pathfind", "Run pathfindR", value = default_config$run_pathfind),
        checkboxInput("save_plots", "Save plots", value = default_config$save_plots),
        checkboxInput("export_data", "Export data", value = default_config$export_data),
        checkboxInput("subtype", "Subtype analysis", value = default_config$subtype),
        checkboxInput("run_interactive", "Launch interactive dashboard after run", value = default_config$Run_Interactive)
      ),
      tabPanel(
        "Groups",
        textInput("treatment_group", "Treatment group", value = default_config$Treatment_group),
        textInput("reference_group", "Reference group", value = default_config$Reference_group)
      ),
      tabPanel(
        "Filtering",
        numericInput("import_qval", "Import q-value", value = default_config$Import_qval, min = 0, step = 0.001),
        numericInput("import_pg_qval", "Import protein group q-value", value = default_config$Import_pg_qval, min = 0, step = 0.001),
        numericInput("group_treshold", "Group threshold", value = default_config$Group_treshold, min = 0, max = 1, step = 0.05),
        numericInput("global_treshold", "Global threshold", value = default_config$Global_treshold, min = 0, max = 1, step = 0.05),
        selectInput("treshold_level", "Threshold level", choices = c("Group", "Global"), selected = default_config$Treshold_level),
        checkboxInput("peptide_level", "Peptide level filtering", value = default_config$peptide_level)
      ),
      tabPanel(
        "Normalization",
        selectInput("norm_method", "Normalization method", choices = c("CycLoess", "mean", "median", "Quantile", "RLR", "GI", "log2", "VSN"), selected = default_config$Norm_method),
        selectInput("batch_corr", "Batch correction", choices = c("eigenms", "sva", "comBat", "none"), selected = default_config$Batch_corr),
        selectInput("agg_method", "Aggregation method", choices = c("maxlfq", "robustsummary", "medianpolish"), selected = default_config$Agg_method)
      ),
      tabPanel(
        "Imputation",
        checkboxInput("impute", "Enable imputation", value = default_config$Impute),
        textInput("impute_mar", "MAR method", value = ifelse(is.na(default_config$impute_MAR), "", default_config$impute_MAR)),
        textInput("impute_mnar", "MNAR method", value = ifelse(is.na(default_config$impute_MNAR), "", default_config$impute_MNAR)),
        textInput("impute_method", "Combined method", value = ifelse(is.na(default_config$impute_method), "", default_config$impute_method))
      ),
      tabPanel(
        "Statistics",
        selectInput("dea_method", "DEA method", choices = c("proDA", "limma", "DEqMS", "msqrob", "msempire"), selected = default_config$DEA_method),
        selectInput("study_design", "Study design", choices = c("unpaired", "paired"), selected = default_config$Study_design),
        numericInput("alpha", "Alpha (FDR)", value = default_config$alpha, min = 0, max = 1, step = 0.01),
        textInput("logfc_cutoff", "logFC cutoff", value = ifelse(is.na(default_config$logFC_cutoff), "", default_config$logFC_cutoff))
      ),
      tabPanel(
        "Pathway & plotting",
        textInput("org_db", "Organism database", value = default_config$Org_db),
        textInput("boxplot_prot", "Boxplot proteins (comma separated)", value = paste(default_config$Boxplot_prot, collapse = ", ")), 
        textInput("agg_prot", "Aggregation proteins (comma separated)", value = paste(default_config$Agg_prot, collapse = ", ")), 
        textAreaInput("group_colors", "Group colors (Group = #HEX per line)", value = paste(sprintf("%s = %s", names(default_config$Group_colors), default_config$Group_colors), collapse = "\n"), rows = 3),
        textAreaInput("expression_lvl_color", "Expression level colors (Name = #HEX)", value = paste(sprintf("%s = %s", names(default_config$Expression_lvl_color), default_config$Expression_lvl_color), collapse = "\n"), rows = 2),
        textAreaInput("heatmap_annot", "Heatmap annotations (Label = column per line)", value = paste(sprintf("%s = %s", names(default_config$heatmap_annot), default_config$heatmap_annot), collapse = "\n"), rows = 4)
      )
    )
  )
)

ui <- fluidPage(
  titlePanel("Pipeline Launcher"),
  conditionalPanel("input.start_setup == 0", landing_panel),
  conditionalPanel("input.start_setup > 0", configuration_panel)
)

server <- function(input, output, session) {
  analysis_state <- reactiveValues(
    status = "Awaiting pipeline launch",
    logs = character(),
    running = FALSE
  )

  append_log <- function(text) {
    analysis_state$logs <- c(analysis_state$logs, paste0(format(Sys.time(), "%H:%M:%S"), " - ", text))
  }

  output$pipeline_status <- renderUI({
    status_text <- analysis_state$status
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
    if (length(analysis_state$logs) == 0) {
      "Logs will appear here after the pipeline starts."
    } else {
      paste(analysis_state$logs, collapse = "\n")
    }
  })

  observeEvent(input$reset_defaults, {
    updateTextInput(session, "project_name", value = default_config$Project_name)
    updateTextInput(session, "analysis_run", value = default_config$Analysis_run)
    updateTextInput(session, "your_name", value = default_config$Your_name)
    updateTextInput(session, "lab_name", value = default_config$Lab_name)
    updateTextInput(session, "contact", value = default_config$Contact)
    updateTextInput(session, "dia_path", value = default_config$DIA_nn_path)
    updateTextInput(session, "fasta_path", value = default_config$FASTA_path)
    updateTextInput(session, "sample_path", value = default_config$Sample_data_path)
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
  })

  observeEvent(input$run_pipeline, {
    req(input$start_setup > 0)

    if (analysis_state$running) {
      showNotification("Pipeline is already running.", type = "warning")
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
      Run_Interactive = isTRUE(input$run_interactive),
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
      showNotification("Please provide at least one heatmap annotation (Label = column).", type = "error")
      return()
    }

    analysis_state$running <- TRUE
    analysis_state$status <- "Running analysis..."
    analysis_state$logs <- character()
    append_log("Pipeline started")

    on.exit({
      analysis_state$running <- FALSE
    }, add = TRUE)

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

      if (isTRUE(result)) {
        if (!is.null(pipeline_env$results_folder)) {
          append_log(paste0("Results saved to ", pipeline_env$results_folder))
        }
        analysis_state$status <- "Analysis completed"
        append_log("Pipeline completed successfully")
      } else {
        analysis_state$status <- "Pipeline failed"
      }
    })
  })
}

app <- shinyApp(ui = ui, server = server)

if (interactive()) {
  runApp(app)
} else {
  app
}
