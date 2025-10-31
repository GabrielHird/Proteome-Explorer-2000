ui <- bslib::page_sidebar(
  title = "Proteome Explorer Targets Dashboard",
  sidebar = bslib::sidebar(
    width = 350,
    shinyjs::useShinyjs(),
    div(
      id = "run_settings",
      h4("Run settings"),
      textInput("project_name", "Project name"),
      textInput("analysis_run", "Analysis run ID"),
      textInput("treatment_group", "Treatment group"),
      textInput("reference_group", "Reference group"),
      textInput("dia_nn_path", "DIA-NN report path"),
      textInput("fasta_path", "FASTA path"),
      textInput("sample_data_path", "Sample metadata path"),
      numericInput("import_qval", "Precursor q-value filter", value = 0.01, min = 0, max = 1, step = 0.001),
      numericInput("import_pg_qval", "Protein group q-value filter", value = 0.01, min = 0, max = 1, step = 0.001),
      numericInput("group_threshold", "Group completeness threshold", value = 0.75, min = 0, max = 1, step = 0.01),
      numericInput("global_threshold", "Global completeness threshold", value = 0.6, min = 0, max = 1, step = 0.01),
      selectInput("threshold_level", "Filter level", choices = c("Group", "Global")),
      selectInput("agg_method", "Aggregation method", choices = c("maxlfq", "robustsummary", "medianpolish")),
      checkboxInput("peptide_level", "Keep peptide level assays", value = TRUE),
      checkboxInput("run_normalizer", "Run normalization", value = TRUE),
      selectInput("norm_method", "Normalization method", choices = c("CycLoess", "Quantile", "GlobalVsN", "Loess")),
      selectInput("batch_corr", "Batch correction", choices = c("none", "eigenms", "sva")),
      checkboxInput("impute", "Perform imputation", value = FALSE),
      textInput("impute_method", "Imputation method override"),
      numericInput("impute_mar", "MAR fraction", value = NA, min = 0, max = 1, step = 0.01),
      numericInput("impute_mnar", "MNAR fraction", value = NA, min = 0, max = 1, step = 0.01),
      checkboxInput("run_dea", "Run differential analysis", value = TRUE),
      selectInput("dea_method", "DEA method", choices = c("msqrob", "msqrob2", "proDA", "limma", "DEqMS", "msEmpiRe")),
      selectInput("study_design", "Study design", choices = c("unpaired", "paired")),
      numericInput("alpha", "Significance level", value = 0.05, min = 0, max = 1, step = 0.001),
      numericInput("logfc_cutoff", "log2 FC cutoff", value = NA, step = 0.1),
      checkboxInput("run_gsea", "Run pathway analysis", value = TRUE),
      textInput("org_db", "Organism database", value = "org.Hs.eg.db"),
      checkboxInput("save_plots", "Save plots", value = TRUE),
      checkboxInput("export_data", "Export data tables", value = TRUE),
      textAreaInput("boxplot_prot", "Boxplot proteins", rows = 2),
      textAreaInput("agg_prot", "Aggregation proteins", rows = 2),
      textAreaInput("group_colors", "Group colours (name = value per line)", rows = 3),
      textAreaInput("heatmap_annotations", "Heatmap annotations (name = column per line)", rows = 3),
      textAreaInput("expression_levels", "Expression level colours (name = value per line)", rows = 3),
      textInput("your_name", "Report author"),
      textInput("lab_name", "Laboratory"),
      textInput("contact", "Contact")
    ),
    actionButton("run_pipeline", "Run pipeline", class = "btn-primary w-100 mb-2"),
    actionButton("reset_defaults", "Reset to defaults", class = "btn-secondary w-100")
  ),
  bslib::navset_tab(
    bslib::nav_panel(
      "Status",
      bslib::layout_columns(
        bslib::card(
          bslib::card_header("Outdated targets"),
          DT::DTOutput("outdated_table")
        ),
        bslib::card(
          bslib::card_header("Pipeline progress"),
          DT::DTOutput("progress_table")
        ),
        bslib::card(
          bslib::card_header("Last run summary"),
          DT::DTOutput("meta_table")
        )
      ),
      bslib::card(
        bslib::card_header("Pipeline log"),
        verbatimTextOutput("run_log", placeholder = TRUE)
      )
    ),
    bslib::nav_panel(
      "Results",
      bslib::layout_columns(
        bslib::card(
          bslib::card_header("Report"),
          uiOutput("report_link")
        ),
        bslib::card(
          bslib::card_header("Exported files"),
          DT::DTOutput("results_table")
        )
      )
    )
  )
)
