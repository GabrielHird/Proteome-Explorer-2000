# Helper functions for configuring the Proteome Explorer pipeline when
# running under {targets} or via the Shiny interface.

pipeline_required_packages <- function() {
  c(
    "yaml",
    "dplyr",
    "filenamer",
    "data.table",
    "BiocParallel",
    "conflicted",
    "readr",
    "readxl",
    "openxlsx",
    "iq",
    "Biobase",
    "SummarizedExperiment",
    "msdap",
    "QFeatures",
    "stringr",
    "NormalyzerDE",
    "imputeLCMD",
    "ProteoMM",
    "sva",
    "hexbin",
    "proDA",
    "limma",
    "DEqMS",
    "msqrob2",
    "msEmpiRe",
    "ComplexHeatmap",
    "ggrepel",
    "uwot",
    "mixOmics",
    "clusterProfiler",
    "org.Hs.eg.db",
    "enrichplot",
    "GSVA",
    "msigdbr",
    "shiny",
    "bslib",
    "plotly",
    "patchwork",
    "naniar",
    "corrplot",
    "FactoMineR",
    "factoextra",
    "DT",
    "queryup",
    "bsicons",
    "RColorBrewer",
    "forcats",
    "paletteer",
    "scales",
    "GGally",
    "hrbrthemes",
    "ggplot2",
    "tibble",
    "tidyr",
    "purrr",
    "Biostrings",
    "magrittr",
    "visNetwork",
    "targets"
  )
}

pipeline_default_config <- function() {
  list(
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
    Group_colors = c(Cancer = "#FABC3C", Normal = "#006D77"),
    heatmap_annot = c(
      Group = "group",
      Batch = "Extraction_batch",
      Study = "Study",
      Membrane = "PExA_membrane"
    ),
    Expression_lvl_color = c(OVER = "#CB4335", UNDER = "#2E86C1"),
    DIA_nn_path = "./Data/DIA-NN output/diann_report_EBP.tsv",
    FASTA_path = "./Data/FASTA file/EBP.fasta",
    Sample_data_path = "./Data/Sample Metadata/EBP_sample_data.xlsx",
    Your_name = "Gabriel Hirdman",
    Lab_name = "Lindstedt Lab",
    Contact = "gabriel.hirdman@med.lu.se"
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

resolve_config_path <- function(path, repo_root) {
  if (!length(path) || is.null(path)) {
    return("")
  }
  path <- path[[1]]
  if (!nzchar(path)) {
    return("")
  }
  expanded <- path.expand(path)
  is_abs <- if (.Platform$OS.type == "windows") {
    grepl("^[A-Za-z]:[\\/]|^\\\\", expanded)
  } else {
    substr(expanded, 1L, 1L) == "/"
  }
  target <- if (is_abs) expanded else file.path(repo_root, expanded)
  normalizePath(target, winslash = "/", mustWork = FALSE)
}

as_character_vec <- function(x) {
  if (is.null(x)) {
    return(character())
  }
  if (is.list(x)) {
    x <- unlist(x, recursive = TRUE, use.names = TRUE)
  }
  as.character(x)
}

as_named_vec <- function(x) {
  if (is.null(x)) {
    return(character())
  }
  if (is.list(x)) {
    x <- unlist(x, recursive = TRUE, use.names = TRUE)
  }
  x <- as.character(x)
  if (is.null(names(x))) {
    names(x) <- rep("", length(x))
  }
  x
}

prepare_pipeline_config <- function(config = NULL,
                                     config_path = "config/pipeline.yml",
                                     repo_root = NULL,
                                     create_dirs = TRUE) {
  repo_root <- repo_root %||% normalizePath(".", winslash = "/", mustWork = TRUE)
  defaults <- pipeline_default_config()

  overrides <- list()
  if (!is.null(config)) {
    overrides <- config
  } else if (!is.null(config_path) && file.exists(config_path)) {
    overrides <- yaml::read_yaml(config_path)
  }

  combined <- modifyList(defaults, overrides, keep.null = TRUE)

  combined$Boxplot_prot <- as_character_vec(combined$Boxplot_prot)
  combined$Agg_prot <- as_character_vec(combined$Agg_prot)
  combined$Group_colors <- as_named_vec(combined$Group_colors)
  combined$heatmap_annot <- as_named_vec(combined$heatmap_annot)
  combined$Expression_lvl_color <- as_named_vec(combined$Expression_lvl_color)

  combined$repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  combined$pipeline_dir <- normalizePath(
    file.path(combined$repo_root, "pipeline"),
    winslash = "/",
    mustWork = TRUE
  )

  combined$DIA_nn_path <- resolve_config_path(combined$DIA_nn_path, combined$repo_root)
  combined$FASTA_path <- resolve_config_path(combined$FASTA_path, combined$repo_root)
  combined$Sample_data_path <- resolve_config_path(combined$Sample_data_path, combined$repo_root)

  saved_data_folder <- file.path(combined$repo_root, "Data", "Saved_data", combined$Analysis_run)
  results_folder <- file.path(combined$repo_root, "Results", "Proteome Explorer", combined$Analysis_run)
  results_data_folder <- file.path(results_folder, "Excel files")

  if (isTRUE(create_dirs)) {
    dir.create(saved_data_folder, recursive = TRUE, showWarnings = FALSE)
    dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
    dir.create(results_data_folder, recursive = TRUE, showWarnings = FALSE)
  }

  must_work_dirs <- isTRUE(create_dirs)

  combined$saved_data_folder <- paste0(
    normalizePath(saved_data_folder, winslash = "/", mustWork = must_work_dirs),
    "/"
  )
  combined$results_folder <- paste0(
    normalizePath(results_folder, winslash = "/", mustWork = must_work_dirs),
    "/"
  )
  combined$results_data_folder <- paste0(
    normalizePath(results_data_folder, winslash = "/", mustWork = must_work_dirs),
    "/"
  )

  combined
}
