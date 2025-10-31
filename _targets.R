library(targets)

source("R/targets_config.R")
source("R/targets_pipeline.R")
source("R/app_targets.R")

tar_option_set(packages = pipeline_required_packages())

list(
  tar_target(ui_runtime_config, app_read_runtime_config()),
  tar_target(import_config, app_build_stage_config(ui_runtime_config, "import")),
  tar_target(preprocess_config, app_build_stage_config(ui_runtime_config, "preprocess")),
  tar_target(dea_config, app_build_stage_config(ui_runtime_config, "dea")),
  tar_target(postprocess_config, app_build_stage_config(ui_runtime_config, "postprocess")),
  tar_target(pathway_config, app_build_stage_config(ui_runtime_config, "pathway")),
  tar_target(export_config, app_build_stage_config(ui_runtime_config, "export")),
  tar_target(report_config, app_build_stage_config(ui_runtime_config, "report")),
  tar_target(import_state, run_import(import_config)),
  tar_target(preprocessed_state, run_preprocess(preprocess_config, import_state)),
  tar_target(dea_state, run_dea(dea_config, preprocessed_state)),
  tar_target(postprocess_state, run_postprocess(postprocess_config, dea_state)),
  tar_target(pathway_state, run_pathway(pathway_config, postprocess_state)),
  tar_target(export_state, run_export(export_config, pathway_state)),
  tar_target(report_state, run_report(report_config, export_state)),
  tar_target(report_file, report_state$report_file, format = "file")
)
