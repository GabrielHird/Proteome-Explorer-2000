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

cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
repo_candidate <- normalizePath(file.path(cwd, ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_candidate, "pipeline"))) {
  repo_candidate <- cwd
}
repo_root <- normalizePath(repo_candidate, winslash = "/", mustWork = TRUE)
app_dir <- normalizePath(file.path(repo_root, "app"), winslash = "/", mustWork = TRUE)
interactive_dir <- normalizePath(file.path(app_dir, "interactive"), winslash = "/", mustWork = TRUE)
pipeline_dir <- normalizePath(file.path(repo_root, "pipeline"), winslash = "/", mustWork = TRUE)

ui_env <- new.env(parent = globalenv())
ui <- source(file.path(interactive_dir, "ui.R"), local = ui_env)$value

server_env <- new.env(parent = globalenv())
server_env$repo_root <- repo_root
server_env$pipeline_dir <- pipeline_dir
server <- source(file.path(interactive_dir, "server.R"), local = server_env)$value

shinyApp(ui = ui, server = server)
