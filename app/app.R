##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(plotly)
  library(DT)
  library(bsicons)
  library(visNetwork)
  library(targets)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
parent_candidate <- normalizePath(file.path(cwd, ".."), winslash = "/", mustWork = FALSE)

potential_roots <- unique(c(parent_candidate, cwd))
repo_root <- NULL
for (candidate in potential_roots) {
  if (!length(candidate) || !dir.exists(candidate)) {
    next
  }
  candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
  if (dir.exists(file.path(candidate, "R"))) {
    repo_root <- candidate
    break
  }
}

if (is.null(repo_root)) {
  stop(
    "Unable to locate the Proteome Explorer repository root. ",
    "Ensure that the 'R' folder exists in the project directory and ",
    "launch the app with runApp('app') from the repository root."
  )
}

repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
app_dir <- normalizePath(file.path(repo_root, "app"), winslash = "/", mustWork = TRUE)
interactive_dir <- normalizePath(file.path(app_dir, "interactive"), winslash = "/", mustWork = TRUE)
ui_env <- new.env(parent = globalenv())
ui <- source(file.path(interactive_dir, "ui.R"), local = ui_env)$value

server_env <- new.env(parent = globalenv())
server_env$repo_root <- repo_root
server <- source(file.path(interactive_dir, "server.R"), local = server_env)$value

shinyApp(ui = ui, server = server)
