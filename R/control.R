control_set <- function() {
  if (process_running()) {
    control_running()
  } else {
    control_stopped()
  }
}

control_running <- function() {
  session <- getDefaultReactiveDomain()
  shinyjs::disable("run_settings")
  shinyjs::disable("run_pipeline")
  updateActionButton(session, "run_pipeline", label = "Running...")
}

control_stopped <- function() {
  session <- getDefaultReactiveDomain()
  shinyjs::enable("run_settings")
  shinyjs::enable("run_pipeline")
  updateActionButton(session, "run_pipeline", label = "Run pipeline")
}
