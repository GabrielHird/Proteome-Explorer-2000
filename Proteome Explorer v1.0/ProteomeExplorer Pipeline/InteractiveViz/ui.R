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
  "))
  ),
  

# Identifications ---------------------------------------------------------

  nav_panel(
    "Identifications",
    layout_sidebar(
      sidebar = sidebar(
        title = "Identifications Controls",
        varSelectInput(
          inputId = "bar_color",
          label = "Color by:",
          data = as.data.frame(colData(qf)),
          selected = as.name("shortname")
        ),
        varSelectInput(
          inputId = "bar_sort",
          label = "Sort by (low to high):",
          data = as.data.frame(colData(qf)),
          selected = as.name("Nr_prot")
        ),
        actionButton("updateBar", "Update Plot")
      ),
      card(
        full_screen = TRUE,
        card_header("Identifications Barplot"),
        div(style = "display: inline-block; margin-bottom: 10px;",
            downloadButton("download_barplot", "Download Plot (.png)", class = "btn btn-default", style = "width: auto;")
        ),
        plotOutput("barplot")
      )
    )
  ),
  

# Missing values ----------------------------------------------------------

  nav_panel(
    "Missing Values",
    card(
      full_screen = TRUE,
      card_header("Missing Values Map"),
      div(style = "display: inline-block; margin-bottom: 10px;",
          downloadButton("download_missMap", "Download Plot (.png)", class = "btn btn-default", style = "width: auto;")
      ),
      plotOutput("missMap")
    )
  ),
  


# Normalization ----------------------------------------------------------

  nav_panel(
    "Normalisation",
    layout_columns(
      fill = FALSE,
      col_widths = c(2, 10),
      value_box(
        title = "Normalization method:",
        value = Norm_method,
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
      
# Protein -----------------------------------------------------------------

  nav_panel(
    "Protein",
    layout_sidebar(
      sidebar = sidebar(
        title = "Select protein",
        selectInput(
          inputId = "gene",
          label = "Select gene:",
          choices = as.character(as.data.frame(rowData(qf[["Results"]]))[["Gene"]]),
          selected = "ALB"
        ),
        actionButton("updateAgg", "Update Plot")
      ),
      layout_columns(
        col_widths = c(6,6),
        
        # Left column
        layout_columns(
          col_widths = c(12,12),
          row_heights = c(1,1),
          
          # First row: boxplot and table side-by-side
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
          
          # Second row: Protein Aggregation Plot
          card(
            full_screen = TRUE,
            card_header("Protein Aggregation Plot"),
            plotOutput("aggregationPlot")
          )
        ),
        
        # Right column remains the Volcano Plot
        card(
          full_screen = TRUE,
          card_header("Volcano Plot"),
          plotlyOutput("volcanoPlotProt")
        )
      )
    )
  ),


# Stats table -------------------------------------------------------------

  nav_panel(
    "Statistics Table",
    DTOutput("results_table")
  ),
  
# Heatmap -----------------------------------------------------------

  nav_panel(
    "Heatmap",
    layout_sidebar(
      sidebar = sidebar(
        title = "Heatmap Annotation Controls",
        checkboxGroupInput(
          inputId = "selected_annotations",
          label = "Select annotation columns:",
          choices = names(as.data.frame(colData(qf))),
          selected = c("group", "Extraction_batch")
        ),
        actionButton("updateHeatmap", "Update Heatmap")
      ),
      card(
        full_screen = TRUE,
        card_header(
          "Heatmap",
          tabsetPanel(
            id = "plot_tab",    # used to track which tab is active
            type = "pills",     # gives you the pill-style tabs
            selected = "custom",  # make "custom" the default active tab
            tabPanel("All Identifications", value = "default"),
            tabPanel("Significant", value = "custom")
          )
        ),
        card_body(
          conditionalPanel(
            condition = "input.plot_tab == 'default'",
            div(style = "float: right; margin-bottom: 10px;",
                downloadButton("download_heatmap_full", "Download Plot (.png)", class = "btn btn-default", style = "width: auto;")
            ),
            plotOutput("defaultHeatmapPlot", height = "1000px")
          ),
          conditionalPanel(
            condition = "input.plot_tab == 'custom'",
            # Use sliderInput with a reactive max in the server (if needed)
            sliderInput(
              inputId = "n_value", 
              label = "Set number for n", 
              min = 1, 
              max = 100,    # placeholder max; update it reactively in the server if necessary
              value = 50
            ),
            div(style = "float: right; margin-bottom: 10px;",
                downloadButton("download_heatmap_sig", "Download Plot (.png)", class = "btn btn-default", style = "width: auto;")
            ),
            plotOutput("customHeatmapPlot", height = "1000px"))
          )
        )
      )
    ),


# Stats results -----------------------------------------------------------

  nav_panel(
    "Volcano",
    card(
      full_screen = TRUE,
      card_header("Volcano Plot"),
      div(style = "display: inline-block; margin-bottom: 10px;",
          downloadButton("download_volcano", "Download Plot (.png)", class = "btn btn-default", style = "width: auto;")
      ),
      plotOutput("volcanoPlot")
    )
  ),
  

# PCA ---------------------------------------------------------------------

  nav_panel(
    "PCA",
    card(
      full_screen = TRUE,
      card_header("PCA Plot"),
      div(style = "display: inline-block; margin-bottom: 10px;",
          downloadButton("download_pca", "Download Plot (.png)", class = "btn btn-default", style = "width: auto;")
      ),
      plotOutput("pcaPlot")
    )
  )
)
