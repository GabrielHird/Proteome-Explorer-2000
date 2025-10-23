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

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
pipeline_dir <- normalizePath(file.path("Proteome Explorer v1.0", "ProteomeExplorer Pipeline"), winslash = "/", mustWork = TRUE)
interactive_dir <- file.path(pipeline_dir, "InteractiveViz")

ui_env <- new.env(parent = globalenv())
ui <- source(file.path(interactive_dir, "ui.R"), local = ui_env)$value

server_env <- new.env(parent = globalenv())
server_env$repo_root <- repo_root
server_env$pipeline_dir <- pipeline_dir
server <- source(file.path(interactive_dir, "server.R"), local = server_env)$value

shinyApp(ui = ui, server = server)
