##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

# Import Raw data ------------------------------------------------------

if(First_pass == TRUE) {

print("Loading DIA-NN data and Sample Data")

Sample_data <- read_excel(Sample_data_path) %>% 
  dplyr::mutate(runCol = shortname)

diann <- fread(DIA_nn_path)
diann$Run <- Sample_data$shortname[match(diann$Run, Sample_data$sample_id)]

diann <- diann %>%
  mutate(Gene = str_trim(str_extract(Genes, "^[^;]+")))

}

# Import to QFeatures -----------------------------------------------------

if(First_pass == TRUE) {

  print("Creating the QFeatures object")
  
  # Create an QFeatures object
  qf <- readQFeaturesFromDIANN(diann,
                               runCol = "Run",
                               quantCols = "Precursor.Normalised",
                               colData = Sample_data)
  
  # Remove samples marked as "exclude = TRUE"
  qf <- subsetByColData(qf, qf$exclude == "FALSE")
  qf <- dropEmptyAssays(qf, dims = 1:2)
  
  # Count
  p_1 <- countUniqueDIA(qf)
  
  # Filter out precursors according to q.value
  qf <- filterFeatures(qf, 
                       ~ Q.Value < Import_qval & PG.Q.Value < Import_pg_qval, 
                       na.rm = TRUE, 
                       keep = TRUE)
  
  # Count
  p_2 <- countUniqueDIA(qf)
  
    print( paste0("Removed a total of ",p_1-p_2," of ",p_1," precursors with the current filter") )
  
  # Replace "0" with "NA"
    print("Replacing 0 with NA")
  qf <- zeroIsNA(qf, i = seq_along(qf))

  # Join all samples together
    print("Joining all samples")
  qf <- joinDIANNAssays(qf, name = "Precursors")
  
  # Delete single sample assays
  qf <- removeAssay(qf, i = names(qf)[ncols(qf) == 1])
  
  # Log2 values
  qf <- logTransform(qf,
                     base = 2,
                     i = "Precursors",
                     name = "PrecursorsLog")
  
  # Fix problematic Peptides
  qf <- fix_multiple_pep(qf, i = "PrecursorsLog")
  

# Summarize to peptides ---------------------------------------------------
  
  if(Agg_method == "maxlfq") {
    fun_agg <- maxLFQ_wrapper
  } else if(Agg_method == "robustsummary") {
    fun_agg <- MsCoreUtils::robustSummary
  } else if(Agg_method == "medianpolish") {
    fun_agg <- MsCoreUtils::medianPolish
  } 
  
  message("Starting Precursor aggregation...")
  qf <- aggregateFeatures(qf,
                          i = "PrecursorsLog",
                          fcol = "Stripped.Sequence",
                          fun = fun_agg,
                          name = "Peptides")
  message("Aggregation finished.")
  
}


# Tools to explore QFeatures ----------------------------------------------

# plot(qf)
# - Quantative
# assay(qf[[1]])
# - Feature data
# rowData(qf[[1]])
# - Sample data
# colData(hl)
# - Interactive
# display(qf)






