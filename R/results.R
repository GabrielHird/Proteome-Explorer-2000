results_outdated <- function(repo_root) {
  tryCatch(app_targets_outdated(repo_root), error = function(e) character())
}

results_meta <- function(repo_root) {
  tryCatch(app_targets_meta(repo_root), error = function(e) data.frame())
}

results_progress <- function(repo_root) {
  tryCatch(app_targets_progress(repo_root), error = function(e) data.frame())
}

results_collect <- function(config, repo_root) {
  app_collect_results(config, repo_root)
}

results_report_path <- function(repo_root, config = NULL) {
  if (!is.null(config)) {
    info <- results_collect(config, repo_root)
    if (!is.na(info$report_file) && file.exists(info$report_file)) {
      return(info$report_file)
    }
  }
  tryCatch({
    old_dir <- setwd(repo_root)
    on.exit(setwd(old_dir))
    targets::tar_read(report_file, store = app_targets_store(repo_root))
  }, error = function(e) NA_character_)
}
