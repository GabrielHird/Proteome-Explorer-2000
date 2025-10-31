# Proteome Explorer 2000

This repository contains the Proteome Explorer 2000 analysis pipeline and interactive Shiny dashboard.

## Project structure

```
app/
  app.R                # Shiny entry point. Calls into the interactive UI/server definitions.
  interactive/
    server.R           # Dashboard server logic.
    ui.R               # Dashboard layout and inputs.

config/
  pipeline.yml         # Default configuration used by the {targets} plan.

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
  targets_config.R     # Utilities for building pipeline configuration objects.
  targets_pipeline.R   # Functional wrappers around the legacy scripts.

_targets.R             # {targets} plan wiring the functional pipeline together.
```

## Running the interactive app

Run the Shiny app from the repository root:

```r
shiny::runApp("app")
```

`app/app.R` automatically resolves the repository root and pipeline directory, so it can be launched from anywhere as long as the working directory points to the repository.

## Executing the pipeline from R

### Functional API

The numbered scripts under `pipeline/` are now wrapped by a functional API. To execute the full workflow from R you can write:

```r
cfg <- prepare_pipeline_config()
state <- run_full_pipeline(cfg)
```

The returned `state` list contains the `qf` object, downstream results, and plotting helpers. You can customise the run by editing `config/pipeline.yml` or by passing an overrides list:

```r
cfg <- prepare_pipeline_config(config = list(Analysis_run = "my_run"))
state <- run_full_pipeline(cfg)
```

### {targets} workflow

A declarative pipeline is available via `{targets}`. From the project root run:

```r
targets::tar_make()
```

The main artefacts are exposed as individual targets (`import_state`, `preprocessed_state`, …, `report_file`). Adjust `config/pipeline.yml` before calling `tar_make()` to change parameters.
