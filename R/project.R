project_repo_root <- function() {
  app_repo_root()
}

project_default_config <- function() {
  pipeline_default_config()
}

project_runtime_config <- function(repo_root = project_repo_root()) {
  app_read_runtime_config(repo_root)
}

project_write_config <- function(config, repo_root = project_repo_root()) {
  app_write_runtime_config(config, repo_root)
}

project_update_inputs <- function(session, cfg) {
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

project_collect_overrides <- function(input, defaults) {
  list(
    Project_name = trimws(input$project_name %||% defaults$Project_name),
    Analysis_run = trimws(input$analysis_run %||% defaults$Analysis_run),
    DIA_nn_path = trimws(input$dia_nn_path %||% defaults$DIA_nn_path),
    FASTA_path = trimws(input$fasta_path %||% defaults$FASTA_path),
    Sample_data_path = trimws(input$sample_data_path %||% defaults$Sample_data_path),
    Import_qval = safe_numeric(input$import_qval),
    Import_pg_qval = safe_numeric(input$import_pg_qval),
    Group_treshold = safe_numeric(input$group_threshold),
    Global_treshold = safe_numeric(input$global_threshold),
    Treshold_level = input$threshold_level %||% defaults$Treshold_level,
    Agg_method = input$agg_method %||% defaults$Agg_method,
    peptide_level = isTRUE(input$peptide_level),
    Run_normalizer = isTRUE(input$run_normalizer),
    Norm_method = input$norm_method %||% defaults$Norm_method,
    Batch_corr = input$batch_corr %||% defaults$Batch_corr,
    Impute = isTRUE(input$impute),
    impute_method = empty_to_na(input$impute_method),
    impute_MAR = safe_numeric(input$impute_mar),
    impute_MNAR = safe_numeric(input$impute_mnar),
    Run_DEA = isTRUE(input$run_dea),
    DEA_method = input$dea_method %||% defaults$DEA_method,
    Study_design = input$study_design %||% defaults$Study_design,
    alpha = safe_numeric(input$alpha),
    logFC_cutoff = safe_numeric(input$logfc_cutoff),
    Run_GSEA = isTRUE(input$run_gsea),
    Org_db = trimws(input$org_db %||% defaults$Org_db),
    save_plots = isTRUE(input$save_plots),
    export_data = isTRUE(input$export_data),
    Treatment_group = trimws(input$treatment_group %||% defaults$Treatment_group),
    Reference_group = trimws(input$reference_group %||% defaults$Reference_group),
    Boxplot_prot = parse_comma_list(input$boxplot_prot),
    Agg_prot = parse_comma_list(input$agg_prot),
    Group_colors = parse_named_lines(input$group_colors),
    heatmap_annot = parse_named_lines(input$heatmap_annotations),
    Expression_lvl_color = parse_named_lines(input$expression_levels),
    Your_name = trimws(input$your_name %||% defaults$Your_name),
    Lab_name = trimws(input$lab_name %||% defaults$Lab_name),
    Contact = trimws(input$contact %||% defaults$Contact)
  )
}

project_merge_config <- function(defaults, overrides) {
  modifyList(defaults, overrides, keep.null = TRUE)
}
