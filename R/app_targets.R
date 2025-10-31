# Helper utilities for coordinating the Shiny app with the {targets}
# pipeline. These functions make it easy to persist UI driven settings,
# translate them into stage specific configurations, and run the pipeline
# while collecting useful metadata for display.

app_repo_root <- function() {
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidate <- cwd

  while (!identical(candidate, dirname(candidate))) {
    if (dir.exists(file.path(candidate, "R")) &&
        file.exists(file.path(candidate, "_targets.R"))) {
      return(candidate)
    }
    candidate <- dirname(candidate)
  }

  stop(
    "Unable to determine the repository root. Launch the app from inside ",
    "the Proteome Explorer repository."
  )
}

app_runtime_config_path <- function(repo_root = app_repo_root()) {
  normalizePath(
    file.path(repo_root, "config", "ui_runtime.yml"),
    winslash = "/",
    mustWork = FALSE
  )
}

app_ensure_runtime_config <- function(repo_root = app_repo_root()) {
  path <- app_runtime_config_path(repo_root)
  if (!file.exists(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    yaml::write_yaml(pipeline_default_config(), path)
  }
  path
}

app_read_runtime_config <- function(repo_root = app_repo_root()) {
  path <- app_ensure_runtime_config(repo_root)
  cfg <- yaml::read_yaml(path)
  modifyList(pipeline_default_config(), cfg, keep.null = TRUE)
}

app_write_runtime_config <- function(config, repo_root = app_repo_root()) {
  defaults <- pipeline_default_config()
  combined <- modifyList(defaults, config, keep.null = TRUE)
  path <- app_ensure_runtime_config(repo_root)
  yaml::write_yaml(combined, path)
  invisible(combined)
}

app_stage_fields <- function(stage) {
  common <- c(
    "Project_name",
    "Analysis_run",
    "First_pass",
    "Run_Interactive",
    "subtype",
    "Treatment_group",
    "Reference_group"
  )
  data_inputs <- c("DIA_nn_path", "FASTA_path", "Sample_data_path")
  filtering <- c(
    "Import_qval",
    "Import_pg_qval",
    "Group_treshold",
    "Global_treshold",
    "Treshold_level",
    "peptide_level",
    "Agg_method"
  )
  normalization <- c(
    "Run_normalizer",
    "Norm_method",
    "Batch_corr",
    "Impute",
    "impute_MAR",
    "impute_MNAR",
    "impute_method"
  )
  dea <- c(
    "Run_DEA",
    "DEA_method",
    "Study_design",
    "alpha",
    "logFC_cutoff"
  )
  visuals <- c(
    "save_plots",
    "Boxplot_prot",
    "Agg_prot",
    "Group_colors",
    "heatmap_annot",
    "Expression_lvl_color"
  )
  pathway <- c("Run_GSEA", "Org_db")
  export <- c("export_data", "save_plots")
  report <- c("Your_name", "Lab_name", "Contact")

  switch(
    stage,
    import = unique(c(common, data_inputs, filtering)),
    preprocess = unique(c(common, normalization, filtering, visuals)),
    dea = unique(c(common, dea, filtering, normalization)),
    postprocess = unique(c(common, visuals, dea)),
    pathway = unique(c(common, pathway, dea, visuals)),
    export = unique(c(common, export, visuals, dea, pathway)),
    report = unique(c(common, report, export, visuals, pathway, dea)),
    stop(sprintf("Unknown pipeline stage '%s'", stage))
  )
}

app_build_stage_config <- function(config,
                                   stage,
                                   repo_root = app_repo_root(),
                                   create_dirs = TRUE) {
  fields <- app_stage_fields(stage)
  overrides <- config[intersect(names(config), fields)]
  prepare_pipeline_config(
    config = overrides,
    repo_root = repo_root,
    create_dirs = create_dirs
  )
}

app_targets_store <- function(repo_root = app_repo_root()) {
  normalizePath(file.path(repo_root, "_targets"), winslash = "/", mustWork = FALSE)
}

app_run_targets_pipeline <- function(repo_root = app_repo_root(), names = NULL) {
  old_dir <- setwd(repo_root)
  on.exit(setwd(old_dir))

  targets::tar_make(names = names)
}

app_targets_outdated <- function(repo_root = app_repo_root()) {
  old_dir <- setwd(repo_root)
  on.exit(setwd(old_dir))
  targets::tar_outdated(store = app_targets_store(repo_root))
}

app_targets_meta <- function(repo_root = app_repo_root()) {
  old_dir <- setwd(repo_root)
  on.exit(setwd(old_dir))
  targets::tar_meta(
    fields = c("name", "seconds", "warnings", "error", "started", "finished"),
    store = app_targets_store(repo_root)
  )
}

app_targets_progress <- function(repo_root = app_repo_root()) {
  old_dir <- setwd(repo_root)
  on.exit(setwd(old_dir))
  progress <- targets::tar_progress(store = app_targets_store(repo_root))
  if (is.null(progress) || !nrow(progress)) {
    progress <- data.frame(
      name = character(),
      time = numeric(),
      progress = numeric(),
      stringsAsFactors = FALSE
    )
  }
  progress
}

app_collect_results <- function(config, repo_root = app_repo_root()) {
  cfg <- prepare_pipeline_config(
    config = config,
    repo_root = repo_root,
    create_dirs = FALSE
  )
  results_dir <- cfg$results_folder
  report_file <- file.path(cfg$results_folder, "Proteome Explorer report.docx")

  files <- character()
  if (dir.exists(results_dir)) {
    files <- list.files(results_dir, recursive = TRUE, full.names = TRUE)
    files <- normalizePath(files, winslash = "/", mustWork = FALSE)
  }

  list(
    results_dir = results_dir,
    report_file = if (file.exists(report_file)) normalizePath(report_file, winslash = "/", mustWork = FALSE) else NA_character_,
    files = files
  )
}
