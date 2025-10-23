##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

norm_cards <- list(
  card(
    full_screen = TRUE,
    card_header("Peptides before normalisation"),
    plotlyOutput("pre_norm_dens")
  ),
  card(
    full_screen = TRUE,
    card_header("Peptides after normalisation"),
    plotlyOutput("post_norm_dens")
  ),
  card(
    full_screen = TRUE,
    card_header("Peptides before normalisation"),
    plotOutput("pre_norm_box")
  ),
  card(
    full_screen = TRUE,
    card_header("Peptides after normalisation"),
    plotOutput("post_norm_box")
  )
)


# UI ----------------------------------------------------------------------

page_navbar(

  # Setup
  title = "Proteome Explorer Interactive",
  theme = bs_theme(bootswatch = "spacelab"),
  tags$head(
    tags$style(HTML("
    /* Force a fixed table layout */
    .dataTables_wrapper table {
      table-layout: fixed;
    }
    /* Set a uniform, smaller width for all header and data cells */
    .dataTables_wrapper table th,
    .dataTables_wrapper table td {
      width: 100px;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      padding: 2px 4px;
      font-size: 12px;
    }
    /* Expand multi-line inputs for easier reading */
    .shiny-textarea-input textarea {
      min-height: 140px;
      font-family: 'Fira Code', 'Source Code Pro', monospace;
      font-size: 13px;
      line-height: 1.35;
      white-space: pre;
    }
    .pe-data-path input.form-control {
      font-family: 'Fira Code', 'Source Code Pro', monospace;
      font-size: 13px;
    }
  "))
  ),


  nav_panel(
    "Run Analysis",
    layout_sidebar(
      sidebar = sidebar(
        title = "Launch Analysis",
        actionButton(
          inputId = "run_pipeline",
          label = tagList(bs_icon("play-fill"), "Run analysis"),
          class = "btn-primary"
        ),
        actionButton(
          inputId = "reset_defaults",
          label = tagList(bs_icon("arrow-counterclockwise"), "Reset defaults"),
          class = "btn-default"
        ),
        hr(),
        uiOutput("pipeline_status"),
        hr(),
        verbatimTextOutput("pipeline_log", placeholder = TRUE)
      ),
      layout_columns(
        col_widths = c(12),
        card(
          card_header("Project information"),
          layout_columns(
            col_widths = c(6, 6),
            textInput("project_name", "Project name", value = "PExA Lungcancer v6"),
            textInput("analysis_run", "Analysis run", value = "EBP_v1_detect")
          ),
          layout_columns(
            col_widths = c(6,6),
            textInput("your_name", "Your name", value = "Gabriel Hirdman"),
            textInput("lab_name", "Lab name", value = "Lindstedt Lab")
          ),
          textInput("contact", "Contact", value = "gabriel.hirdman@med.lu.se")
        ),
        card(
          card_header("Run settings"),
          layout_columns(
            col_widths = c(4,4,4),
            checkboxInput("first_pass", "First pass", value = TRUE),
            checkboxInput("run_normalizer", "Run normalizer", value = TRUE),
            checkboxInput("run_dea", "Run DEA", value = TRUE)
          ),
          layout_columns(
            col_widths = c(4,4,4),
            checkboxInput("run_gsea", "Run GSEA", value = TRUE),
            checkboxInput("save_plots", "Save plots", value = TRUE)
          ),
          layout_columns(
            col_widths = c(4,4,4),
            checkboxInput("export_data", "Export data", value = TRUE),
            checkboxInput("subtype", "Subtype analysis", value = FALSE),
            checkboxInput("run_interactive", "Launch dashboard after run", value = FALSE)
          )
        ),
        card(
          card_header("Group setup"),
          layout_columns(
            col_widths = c(6,6),
            textInput("treatment_group", "Treatment group", value = "Cancer"),
            textInput("reference_group", "Reference group", value = "Normal")
          )
        ),
        card(
          card_header("Filtering"),
          layout_columns(
            col_widths = c(6,6),
            numericInput("import_qval", "Import q-value", value = 0.01, min = 0, step = 0.001),
            numericInput("import_pg_qval", "Import protein group q-value", value = 0.01, min = 0, step = 0.001)
          ),
          layout_columns(
            col_widths = c(4,4,4),
            numericInput("group_treshold", "Group threshold", value = 0.75, min = 0, max = 1, step = 0.05),
            numericInput("global_treshold", "Global threshold", value = 0.6, min = 0, max = 1, step = 0.05),
            selectInput("treshold_level", "Threshold level", choices = c("Group", "Global"), selected = "Group")
          ),
          checkboxInput("peptide_level", "Peptide level filtering", value = TRUE)
        ),
        card(
          card_header("Normalization"),
          layout_columns(
            col_widths = c(6,6),
            selectInput(
              "norm_method",
              "Normalization method",
              choices = c("CycLoess", "mean", "median", "Quantile", "RLR", "GI", "log2", "VSN"),
              selected = "CycLoess"
            ),
            selectInput(
              "batch_corr",
              "Batch correction",
              choices = c("eigenms", "sva", "comBat", "none"),
              selected = "eigenms"
            )
          )
        ),
        card(
          card_header("Aggregation"),
          selectInput(
            "agg_method",
            "Aggregation method",
            choices = c("maxlfq", "robustsummary", "medianpolish"),
            selected = "robustsummary"
          )
        ),
        card(
          card_header("Imputation"),
          layout_columns(
            col_widths = c(4,4,4),
            checkboxInput("impute", "Enable imputation", value = FALSE),
            textInput("impute_mar", "MAR method", value = ""),
            textInput("impute_mnar", "MNAR method", value = "")
          ),
          textInput("impute_method", "Combined imputation method", value = "")
        ),
        card(
          card_header("DEA analysis"),
          layout_columns(
            col_widths = c(6,6),
            selectInput("dea_method", "DEA method", choices = c("proDA", "limma", "DEqMS", "msqrob", "msempire"), selected = "msqrob"),
            selectInput("study_design", "Study design", choices = c("unpaired", "paired"), selected = "unpaired")
          ),
          layout_columns(
            col_widths = c(6,6),
            numericInput("alpha", "Alpha (FDR)", value = 0.05, min = 0, max = 1, step = 0.01),
            textInput("logfc_cutoff", "logFC cutoff", value = "")
          )
        ),
        card(
          card_header("Pathway"),
          textInput("org_db", "Organism database", value = "org.Hs.eg.db")
        ),
        card(
          card_header("Points of interest"),
          layout_columns(
            col_widths = c(6,6),
            textInput("boxplot_prot", "Boxplot proteins (comma separated)", value = "CSTA"),
            textInput("agg_prot", "Aggregation proteins (comma separated)", value = "CSTA")
          )
        ),
        card(
          card_header("Colors"),
          textAreaInput(
            "group_colors",
            "Group colors (one per line: Group = #HEX)",
            value = "Cancer = #FABC3C\nNormal = #006D77",
            rows = 6,
            width = "100%"
          ),
          textAreaInput(
            "expression_lvl_color",
            "Expression level colors (Name = #HEX)",
            value = "OVER = #CB4335\nUNDER = #2E86C1",
            rows = 4,
            width = "100%"
          )
        ),
        card(
          card_header("Heatmap annotations"),
          textAreaInput(
            "heatmap_annot",
            "Annotation mappings (Label = column)",
            value = "Group = group\nBatch = Extraction_batch\nStudy = Study\nMembrane = PExA_membrane",
            rows = 6,
            width = "100%"
          )
        ),
        card(
          card_header("Data paths"),
          div(
            class = "mb-3",
            style = "display: flex; gap: 0.5rem; align-items: flex-end;",
            div(class = "pe-data-path", style = "flex: 1;", textInput("dia_path", "DIA-NN report path", value = "./Data/DIA-NN output/diann_report_EBP.tsv", width = "100%")),
            shinyFiles::shinyFilesButton("dia_path_browse", label = "Browse…", title = "Select DIA-NN report", multiple = FALSE, class = "btn btn-default")
          ),
          div(
            class = "mb-3",
            style = "display: flex; gap: 0.5rem; align-items: flex-end;",
            div(class = "pe-data-path", style = "flex: 1;", textInput("fasta_path", "FASTA path", value = "./Data/FASTA file/EBP.fasta", width = "100%")),
            shinyFiles::shinyFilesButton("fasta_path_browse", label = "Browse…", title = "Select FASTA file", multiple = FALSE, class = "btn btn-default")
          ),
          div(
            class = "mb-3",
            style = "display: flex; gap: 0.5rem; align-items: flex-end;",
            div(class = "pe-data-path", style = "flex: 1;", textInput("sample_path", "Sample metadata path", value = "./Data/Sample Metadata/EBP_sample_data.xlsx", width = "100%")),
            shinyFiles::shinyFilesButton("sample_path_browse", label = "Browse…", title = "Select sample metadata", multiple = FALSE, class = "btn btn-default")
          )
        )
      )
    )
  ),

  nav_panel(
    "View Analysis Results",
    navset_tab(
      nav_panel(
        "Identifications",
        layout_sidebar(
          sidebar = sidebar(
            title = "Identifications Controls",
            uiOutput("bar_color_input"),
            uiOutput("bar_sort_input"),
            actionButton("updateBar", "Update Plot")
          ),
          card(
            full_screen = TRUE,
            card_header("Identifications Barplot"),
            div(
              style = "display: inline-block; margin-bottom: 10px;",
              downloadButton("download_barplot", "Download Plot (.png)", class = "btn btn-default", style = "width: auto;")
            ),
            plotOutput("barplot")
          )
        )
      ),

      nav_panel(
        "Missing Values",
        card(
          full_screen = TRUE,
          card_header("Missing Values Map"),
          div(
            style = "display: inline-block; margin-bottom: 10px;",
            downloadButton("download_missMap", "Download Plot (.png)", class = "btn btn-default", style = "width: auto;")
          ),
          plotOutput("missMap")
        )
      ),

      nav_panel(
        "Normalisation",
        layout_columns(
          fill = FALSE,
          col_widths = c(2, 10),
          value_box(
            title = "Normalization method:",
            value = textOutput("norm_method_display"),
            showcase = bs_icon("clipboard-check-fill"),
            showcase_layout = c("left center"),
            theme = "primary",
            max_height = "200px"
          ),
          tagList(
            layout_columns(
              col_widths = c(6,6),
              norm_cards[[1]],
              norm_cards[[2]]
            ),
            layout_columns(
              col_widths = c(6,6),
              norm_cards[[3]],
              norm_cards[[4]]
            )
          )
        )
      ),

      nav_panel(
        "Protein",
        layout_sidebar(
          sidebar = sidebar(
            title = "Select protein",
            uiOutput("gene_selector"),
            actionButton("updateAgg", "Update Plot")
          ),
          layout_columns(
            col_widths = c(6,6),
            layout_columns(
              col_widths = c(12,12),
              row_heights = c(1,1),
              layout_columns(
                col_widths = c(6,6),
                card(
                  full_screen = TRUE,
                  card_header("Boxplot"),
                  plotOutput("prot_boxplot")
                ),
                card(
                  full_screen = TRUE,
                  card_header("Protein Details"),
                  uiOutput("protein_info")
                )
              ),
              card(
                full_screen = TRUE,
                card_header("Protein Aggregation Plot"),
                plotOutput("aggregationPlot")
              )
            ),
            card(
              full_screen = TRUE,
              card_header("Volcano Plot"),
              plotlyOutput("volcanoPlotProt")
            )
          )
        )
      ),

      nav_panel(
        "Statistics Table",
        DTOutput("results_table")
      ),

      nav_panel(
        "Heatmap",
        layout_sidebar(
          sidebar = sidebar(
            title = "Heatmap Annotation Controls",
            uiOutput("heatmap_annotations"),
            actionButton("updateHeatmap", "Update Heatmap")
          ),
          card(
            full_screen = TRUE,
            card_header(
              "Heatmap",
              tabsetPanel(
                id = "plot_tab",
                type = "pills",
                selected = "custom",
                tabPanel("All Identifications", value = "default"),
                tabPanel("Significant", value = "custom")
              )
            ),
            card_body(
              conditionalPanel(
                condition = "input.plot_tab == 'default'",
                div(
                  style = "float: right; margin-bottom: 10px;",
                  downloadButton("download_heatmap_full", "Download Plot (.png)", class = "btn btn-default", style = "width: auto;")
                ),
                plotOutput("defaultHeatmapPlot", height = "1000px")
              ),
              conditionalPanel(
                condition = "input.plot_tab == 'custom'",
                sliderInput(
                  inputId = "n_value",
                  label = "Set number for n",
                  min = 1,
                  max = 100,
                  value = 50
                ),
                div(
                  style = "float: right; margin-bottom: 10px;",
                  downloadButton("download_heatmap_sig", "Download Plot (.png)", class = "btn btn-default", style = "width: auto;")
                ),
                plotOutput("customHeatmapPlot", height = "1000px")
              )
            )
          )
        )
      ),

      nav_panel(
        "Volcano",
        card(
          full_screen = TRUE,
          card_header("Volcano Plot"),
          div(
            style = "display: inline-block; margin-bottom: 10px;",
            downloadButton("download_volcano", "Download Plot (.png)", class = "btn btn-default", style = "width: auto;")
          ),
          plotOutput("volcanoPlot")
        )
      ),

      nav_panel(
        "PCA",
        card(
          full_screen = TRUE,
          card_header("PCA Plot"),
          div(
            style = "display: inline-block; margin-bottom: 10px;",
            downloadButton("download_pca", "Download Plot (.png)", class = "btn btn-default", style = "width: auto;")
          ),
          plotOutput("pcaPlot")
        )
      )
    )
  )
)
