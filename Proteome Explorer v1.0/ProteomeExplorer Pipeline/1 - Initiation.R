##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

# Import packages ---------------------------------------------------------

suppressPackageStartupMessages({
  
  # Setup/helper
  source("./Scripts/Proteome Explorer v1.0/ProteomeExplorer Pipeline/PExp - Functions.R")
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
  library(pathfindR)
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
                 base::as.factor)

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

n_cores <- parallel::detectCores() - 1
BiocParallel::register(BiocParallel::SnowParam(workers = n_cores), default = TRUE)

# Run analysis ------------------------------------------------------------


source("./Scripts/Proteome Explorer v1.0/ProteomeExplorer Pipeline/2 - Import.R")

source("./Scripts/Proteome Explorer v1.0/ProteomeExplorer Pipeline/3 - Data preprocessing.R")

if(subtype == TRUE) {
  source("./Scripts/Proteome Explorer v1.0/ProteomeExplorer Pipeline/4 - Subtype DEA Analysis.R")
} else {
  source("./Scripts/Proteome Explorer v1.0/ProteomeExplorer Pipeline/4 - DEA Analysis.R")
}

source("./Scripts/Proteome Explorer v1.0/ProteomeExplorer Pipeline/5 - Postprocessing.R")

source("./Scripts/Proteome Explorer v1.0/ProteomeExplorer Pipeline/6 - Pathway analysis.R")

source("./Scripts/Proteome Explorer v1.0/ProteomeExplorer Pipeline/7 - Export plots.R")

# Print report
source("./Scripts/Proteome Explorer v1.0/ProteomeExplorer Pipeline/Final - Print report.R")

if(Run_Interactive == TRUE) {
  runApp("./Scripts/Proteome Explorer v1.0/ProteomeExplorer Pipeline/InteractiveViz")
}


