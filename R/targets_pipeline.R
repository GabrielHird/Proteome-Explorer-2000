# Functions that wrap the legacy pipeline scripts so they can be executed
# as pure functions within a {targets} workflow or from other R contexts.

load_pipeline_packages <- function(packages = NULL) {
  if (is.null(packages)) {
    packages <- if (exists("pipeline_required_packages", mode = "function")) {
      pipeline_required_packages()
    } else {
      character()
    }
  }

  packages <- unique(packages)

  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Required package '%s' is not installed.", pkg))
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }

  if ("conflicted" %in% loadedNamespaces()) {
    conflicted::conflicts_prefer(
      QFeatures::longForm,
      dplyr::slice,
      dplyr::filter,
      dplyr::rename,
      dplyr::setdiff,
      dplyr::select,
      base::as.factor,
      base::unname
    )
  }

  invisible(packages)
}

create_pipeline_helpers_env <- function(repo_root,
                                        packages = NULL,
                                        parent_env = NULL) {
  if (is.null(parent_env)) {
    parent_env <- parent.frame()
    if (identical(parent_env, emptyenv())) {
      parent_env <- .GlobalEnv
    }
  }

  loaded_packages <- load_pipeline_packages(packages)
  helpers <- new.env(parent = parent_env)
  attr(helpers, "pex_loaded_packages") <- loaded_packages
  sys.source(file.path(repo_root, "R", "pipeline_functions.R"), envir = helpers)
  helpers
}

run_pipeline_script <- function(script,
                                config,
                                state = list(),
                                helpers_env = NULL) {
  if (is.null(helpers_env)) {
    helpers_env <- create_pipeline_helpers_env(config$repo_root)
  }

  runtime_env <- new.env(parent = helpers_env)
  list2env(config, envir = runtime_env)
  if (length(state)) {
    list2env(state, envir = runtime_env)
  }

  script_path <- file.path(config$pipeline_dir, script)
  if (!file.exists(script_path)) {
    stop(sprintf("Pipeline script '%s' was not found at %s", script, script_path))
  }

  sys.source(script_path, envir = runtime_env)

  all_names <- ls(runtime_env, all.names = TRUE)
  keep <- setdiff(all_names, ".__S3MethodsTable__.")
  mget(keep, envir = runtime_env, inherits = FALSE)
}

run_import <- function(config, state = list(), helpers_env = NULL) {
  run_pipeline_script("02_import.R", config, state, helpers_env)
}

run_preprocess <- function(config, state, helpers_env = NULL) {
  run_pipeline_script("03_data_preprocessing.R", config, state, helpers_env)
}

run_dea <- function(config, state, helpers_env = NULL) {
  script <- if (isTRUE(config$subtype)) {
    "05_subtype_dea_analysis.R"
  } else {
    "04_dea_analysis.R"
  }
  run_pipeline_script(script, config, state, helpers_env)
}

run_postprocess <- function(config, state, helpers_env = NULL) {
  run_pipeline_script("06_postprocessing.R", config, state, helpers_env)
}

run_pathway <- function(config, state, helpers_env = NULL) {
  run_pipeline_script("07_pathway_analysis.R", config, state, helpers_env)
}

run_export <- function(config, state, helpers_env = NULL) {
  run_pipeline_script("08_export_plots.R", config, state, helpers_env)
}

run_report <- function(config, state, helpers_env = NULL) {
  run_pipeline_script("final_print_report.R", config, state, helpers_env)
}

run_full_pipeline <- function(config,
                              steps = c("import",
                                        "preprocess",
                                        "dea",
                                        "postprocess",
                                        "pathway",
                                        "export",
                                        "report"),
                              callbacks = list()) {
  helpers_env <- create_pipeline_helpers_env(config$repo_root)
  state <- list()
  total_steps <- length(steps)

  globals_to_expose <- intersect("Group_colors", names(config))
  original_globals <- list()
  missing_sentinel <- new.env(parent = emptyenv())

  if (length(globals_to_expose)) {
    for (nm in globals_to_expose) {
      if (exists(nm, envir = .GlobalEnv, inherits = FALSE)) {
        original_globals[[nm]] <- get(nm, envir = .GlobalEnv, inherits = FALSE)
      } else {
        original_globals[[nm]] <- missing_sentinel
      }
      assign(nm, config[[nm]], envir = .GlobalEnv)
    }

    on.exit({
      for (nm in names(original_globals)) {
        if (identical(original_globals[[nm]], missing_sentinel)) {
          rm(list = nm, envir = .GlobalEnv)
        } else {
          assign(nm, original_globals[[nm]], envir = .GlobalEnv)
        }
      }
    }, add = TRUE)
  }

  step_functions <- list(
    import = run_import,
    preprocess = run_preprocess,
    dea = run_dea,
    postprocess = run_postprocess,
    pathway = run_pathway,
    export = run_export,
    report = run_report
  )

  for (i in seq_along(steps)) {
    step <- steps[[i]]
    fn <- step_functions[[step]]
    if (is.null(fn)) {
      stop(sprintf("Unknown pipeline step '%s'", step))
    }
    if (!is.null(callbacks$on_step_start)) {
      callbacks$on_step_start(step = step, index = i, total = total_steps)
    }
    state <- fn(config, state, helpers_env = helpers_env)
    if (!is.null(callbacks$on_step_complete)) {
      callbacks$on_step_complete(step = step, index = i, total = total_steps)
    }
  }

  state
}
