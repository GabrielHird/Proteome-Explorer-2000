library(testthat)

source(file.path("R", "targets_config.R"))

test_that("pipeline_required_packages enumerates {targets} packages", {
  pkgs <- pipeline_required_packages()
  expect_type(pkgs, "character")
  expect_gt(length(pkgs), 0)

  old_pkgs <- targets::tar_option_get("packages")
  on.exit(targets::tar_option_set(packages = old_pkgs), add = TRUE)

  targets::tar_option_set(packages = pkgs)
  expect_setequal(targets::tar_option_get("packages"), pkgs)
})

test_that("prepare_pipeline_config tolerates missing output directories when create_dirs = FALSE", {
  tmp_repo <- tempfile("pex_repo")
  dir.create(tmp_repo)
  on.exit(unlink(tmp_repo, recursive = TRUE), add = TRUE)

  dir.create(file.path(tmp_repo, "pipeline"), recursive = TRUE)

  saved_path <- file.path(tmp_repo, "Data", "Saved_data", pipeline_default_config()$Analysis_run)
  results_path <- file.path(tmp_repo, "Results", "Proteome Explorer", pipeline_default_config()$Analysis_run)
  results_data_path <- file.path(results_path, "Excel files")

  expect_false(dir.exists(saved_path))
  expect_false(dir.exists(results_path))
  expect_false(dir.exists(results_data_path))

  expect_no_error(
    prepare_pipeline_config(
      config = NULL,
      config_path = NULL,
      repo_root = tmp_repo,
      create_dirs = FALSE
    )
  )
})
