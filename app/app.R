##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(plotly)
  library(DT)
  library(bsicons)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

get_app_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_arg_idx <- grep(file_arg, cmd_args)
  if (length(file_arg_idx) > 0) {
    file_path <- sub(file_arg, "", cmd_args[file_arg_idx[1]])
    return(dirname(normalizePath(file_path, winslash = "/", mustWork = TRUE)))
  }
  this_file <- NULL
  for (i in rev(seq_len(sys.nframe()))) {
    frame <- sys.frame(i)
    if (!is.null(frame$ofile)) {
      this_file <- frame$ofile
      break
    }
  }
  if (!is.null(this_file)) {
    return(dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

app_dir <- get_app_dir()
repo_root <- normalizePath(file.path(app_dir, ".."), winslash = "/", mustWork = TRUE)
pipeline_dir <- normalizePath(file.path(repo_root, "pipeline"), winslash = "/", mustWork = TRUE)
interactive_dir <- normalizePath(file.path(app_dir, "interactive"), winslash = "/", mustWork = TRUE)

ui_env <- new.env(parent = globalenv())
ui <- source(file.path(interactive_dir, "ui.R"), local = ui_env)$value

server_env <- new.env(parent = globalenv())
server_env$repo_root <- repo_root
server_env$pipeline_dir <- pipeline_dir
server <- source(file.path(interactive_dir, "server.R"), local = server_env)$value

shinyApp(ui = ui, server = server)
