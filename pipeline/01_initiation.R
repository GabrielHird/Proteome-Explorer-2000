##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

# Import packages ---------------------------------------------------------

pipeline_dir <- getOption("pex_pipeline_dir", default = "./pipeline")
pipeline_dir <- normalizePath(pipeline_dir, winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(pipeline_dir, ".."), winslash = "/", mustWork = TRUE)

pipeline_env <- environment()

suppressPackageStartupMessages({

  # Setup/helper
  sys.source(file.path(repo_root, "R", "pipeline_functions.R"), envir = pipeline_env)
  library(dplyr)
  library(filenamer)
  library(data.table)
  library(BiocParallel)
  library(conflicted)

  # Importing
  library(readr)
  library(readxl)
  library(openxlsx)
  library(iq)
  library(Biobase)
  library(SummarizedExperiment)
  library(msdap)
  library(QFeatures)
  library(stringr)

  # Preprocessing
  library(NormalyzerDE)
  library(imputeLCMD)
  library(ProteoMM)
  library(sva)
  library(hexbin)

  # DEA Analysis
  library(proDA)
  library(limma)
  library(DEqMS)
  library(msqrob2)
  library(msEmpiRe)

  # Postprocessing
  library(ComplexHeatmap)
  library(ggrepel)
  library(uwot)
  library(mixOmics)

  # Pathway
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(GSVA)
  library(msigdbr)

  # InteractiveViz
  library(shiny)
  library(bslib)
  library(plotly)
  library(patchwork)
  library(naniar)
  library(corrplot)
  library(FactoMineR)
  library(factoextra)
  library(DT)
  library(queryup)
  library(bsicons)

  # Publication Plots
  library(RColorBrewer)
  library(forcats)
  library(paletteer)
  library(scales)
  library(GGally)
  library(hrbrthemes)

})

pix <- function(n) {
  n*0.0104166667
}

conflicts_prefer(QFeatures::longFormat,
                 dplyr::slice,
                 dplyr::filter,
                 dplyr::rename,
                 dplyr::setdiff,
                 dplyr::select,
                 base::as.factor,
                 base::unname)

# Folders -----------------------------------------------------------------

# Path names
saved_data_folder    <- paste0("./Data/Saved_data/",Analysis_run,"/")
results_folder       <- paste0("./Results/Proteome Explorer/", Analysis_run,"/")
results_data_folder  <- paste0(results_folder,"/Excel files/")

if(!dir.exists(paths = saved_data_folder)) {
  make_path(saved_data_folder)
}

if(!dir.exists(paths = results_folder)) {
  make_path(results_folder)
}

if(!dir.exists(paths = results_data_folder)) {
  make_path(results_data_folder)
}


# Multicore ---------------------------------------------------------------

available_cores <- parallel::detectCores()
n_cores <- max(1L, available_cores - 1L)

if (n_cores == 1L) {
  message("Parallel processing disabled; using a single core.")
}

BiocParallel::register(BiocParallel::SnowParam(workers = n_cores), default = TRUE)

# Run analysis ------------------------------------------------------------


sys.source(file.path(pipeline_dir, "02_import.R"), envir = pipeline_env)

sys.source(file.path(pipeline_dir, "03_data_preprocessing.R"), envir = pipeline_env)

if(subtype == TRUE) {
  sys.source(file.path(pipeline_dir, "05_subtype_dea_analysis.R"), envir = pipeline_env)
} else {
  sys.source(file.path(pipeline_dir, "04_dea_analysis.R"), envir = pipeline_env)
}

sys.source(file.path(pipeline_dir, "06_postprocessing.R"), envir = pipeline_env)

sys.source(file.path(pipeline_dir, "07_pathway_analysis.R"), envir = pipeline_env)

sys.source(file.path(pipeline_dir, "08_export_plots.R"), envir = pipeline_env)

# Print report
sys.source(file.path(pipeline_dir, "final_print_report.R"), envir = pipeline_env)

if(Run_Interactive == TRUE) {
  runApp(file.path(repo_root, "app"))
}
