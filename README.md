# Proteome Explorer 2000

This repository contains the Proteome Explorer 2000 analysis pipeline and interactive Shiny dashboard.

## Project structure

```
app/
  app.R                # Shiny entry point. Calls into the interactive UI/server definitions.
  interactive/
    server.R           # Dashboard server logic.
    ui.R               # Dashboard layout and inputs.

pipeline/
  01_initiation.R      # Main orchestrator that loads helpers and runs each pipeline stage.
  02_import.R          # Import DIA-NN output and metadata.
  03_data_preprocessing.R
  04_dea_analysis.R
  05_subtype_dea_analysis.R
  06_postprocessing.R
  07_pathway_analysis.R
  08_export_plots.R
  final_print_report.R # Report generation helpers.
  publication_plots/   # Specialised publication plotting scripts.

R/
  pipeline_functions.R # Shared helper functions sourced by the pipeline scripts.
```

## Running the interactive app

Run the Shiny app from the repository root:

```r
shiny::runApp("app")
```

`app/app.R` automatically resolves the repository root and pipeline directory, so it can be launched from anywhere as long as the working directory points to the repository.

## Executing the pipeline from R

You can also source the pipeline directly:

```r
source("pipeline/01_initiation.R", local = TRUE)
```

The script expects the working directory to be the repository root (this is handled automatically when launching the Shiny app).
