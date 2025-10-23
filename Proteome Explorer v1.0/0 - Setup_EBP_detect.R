##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################
  
# Version control ---------------------------------------------------------

Project_name = "PExA Lungcancer v6"
Analysis_run = "EBP_v1_detect"

# Run settings ------------------------------------------------------------

First_pass      = T  # Set to false to use already created data
Run_normalizer  = T  # Only needs to rerun if filtering is changed
Run_DEA         = T  # Run DEA analysis
Run_GSEA        = T  # Run pathway analysis
run_pathfind    = F  # (Optional pathfindeR)
save_plots      = T  # Save plots to file
export_data     = T  # Export stats

subtype         = FALSE

Run_Interactive = T # Launch the shiny dashboard

# Analysis settings -------------------------------------------------------

# Group setup
Treatment_group = "Cancer"
Reference_group = "Normal"

# Filtering
Import_qval     = 0.01 # Only retain precursors with a Q.Value less than..
Import_pg_qval  = 0.01 # Only retain precursors with a PQ.Q.Value less than...

Group_treshold  = 0.75 # Keep proteins identified in at least ... of one group  
Global_treshold = 0.6  # Keep proteins identified in at least ... of all samples
Treshold_level  = "Group" # Set to "Global" or "Group" to filter at sample or group level

peptide_level   = TRUE # Set to TRUE if you want to filter at peptide level (Group only)

# Normalization
Norm_method     = "CycLoess" # One of "CycLoess", "mean", "median", "Quantile", "RLR", "GI", "log2", "VSN"
Batch_corr      = "eigenms" # Select either "eigenms", "sva", "comBat" or "none"

# Aggregation
Agg_method      = "robustsummary" # Select one of "maxlfq", "robustsummary" or "medianpolish"

# Imputation
Impute          = FALSE
impute_MAR      = NA
impute_MNAR     = NA
impute_method   = NA # If only one method selected, otherwise NA

# DEA analysis
DEA_method      = "msqrob" # One of: "proDA", "limma", "DEqMS", "msqrob", "msempire"
Study_design    = "unpaired" # One of: "unpaired" or "paired" (Only for MSqRob so far)
alpha           = 0.05    # Min Q-value (FDR corrected p-value) that a protein need to pass as significant
logFC_cutoff    = NA     # Min absolute logFC (+/-) a protein need to pass as significant, if set to "NA" -> use bootstrapping.



# Pathway
Org_db          = "org.Hs.eg.db"

# Plotting settings -------------------------------------------------------

# Point of Interests
Boxplot_prot <- c("CSTA")
Agg_prot     <- c("CSTA")

# Colors
Group_colors    = c(
  "Cancer"   = "#FABC3C",
  "Normal" = "#006D77"
) 

# Heatmap annotations
heatmap_annot <- c(
  "Group"    = "group",
  "Batch"    = "Extraction_batch",
  "Study"    = "Study",
  "Membrane" = "PExA_membrane"
)

Expression_lvl_color <- c(
  OVER = "#CB4335",
  UNDER = "#2E86C1"
)

# Setup -------------------------------------------------------------------

# Data paths
DIA_nn_path      = "./Data/DIA-NN output/diann_report_EBP.tsv"
FASTA_path       = "./Data/FASTA file/EBP.fasta"
Sample_data_path = "./Data/Sample Metadata/EBP_sample_data.xlsx"

# Experiment information
Your_name = "Gabriel Hirdman"
Lab_name  = "Lindstedt Lab"
Contact   = "gabriel.hirdman@med.lu.se" 

# Run ---------------------------------------------------------------------

pipeline_dir <- "./ProteomeExplorer Pipeline"
options(pex_pipeline_dir = pipeline_dir)
source(file.path(pipeline_dir, "1 - Initiation.R"))
