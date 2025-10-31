.process_state <- new.env(parent = emptyenv())
.process_state$running <- FALSE

process_status <- function() {
  list(running = process_running())
}

process_running <- function() {
  isTRUE(.process_state$running)
}

process_set_running <- function(value) {
  .process_state$running <- isTRUE(value)
  invisible(process_status())
}

process_run_pipeline <- function(repo_root, names = NULL) {
  log_messages <- character()
  result <- tryCatch({
    withCallingHandlers({
      app_run_targets_pipeline(repo_root, names = names)
    },
    message = function(m) {
      log_messages <<- c(log_messages, m$message)
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      log_messages <<- c(log_messages, paste("Warning:", w$message))
      invokeRestart("muffleWarning")
    })
    list(success = TRUE, message = NULL)
  }, error = function(e) {
    list(success = FALSE, message = conditionMessage(e))
  })
  list(success = result$success, message = result$message, log = log_messages)
}

process_cancel <- function() {
  shiny::showNotification(
    "Pipeline cancellation is not available in this configuration.",
    type = "warning"
  )
}
