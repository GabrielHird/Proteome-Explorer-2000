##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

if(First_pass == TRUE) {
# Filter Contaminants -----------------------------------------------------

### Filtering data:
print(paste0("Nr of peptides pre filtering: ", nrow(rowData(qf[["Peptides"]]))))

# Filter out Contaminants
qf <- filterFeatures(qf,
                     i = "Peptides",
                     ~ !startsWith(Protein.Group, "Cont_"))

print(paste0("Nr of peptides post filtering: ", nrow(rowData(qf[["Peptides"]]))))


# Fix Canonical -----------------------------------------------------------

## Filter to keep only features that have a gene name
qf <- filterFeatures(qf,i = "Peptides",
                     ~ Genes != "")

qf <- filterCanonical(qf,
                      assay = "Peptides",
                      fasta = FASTA_path)


# Filter NA values -------------------------------------------------------

if(Treshold_level == "Group") {
  qf <- flagPeptidesForRemoval(qf, 
                               peptide_assay = "Peptides", 
                               group_var = "group", 
                               threshold = Group_treshold, 
                               protein_col = "Protein.Group",
                               peptide_level = peptide_level)
  
  qf <- filterFeatures(qf, ~ removePeptide == FALSE, i = "Peptides", keep = TRUE)
  
  # Filter out Peptides that are NA in all 
  qf <- filterNA(qf, i= "Peptides", pNA = 0.99)
  
}

if(Treshold_level == "Global") {
  qf <- filterNA(qf, i = "Peptides", pNA = Global_treshold)
}

print(paste0("Nr of peptides post filtering at group level: ", nrow(rowData(qf[["Peptides"]]))))


# Normalization -----------------------------------------------------------

normalizer_output <- file.path(results_folder, "normalyzer", paste0(Norm_method, "-normalized.txt"))

if(Run_normalizer == TRUE) {
# Run normalizer
# De Log
SumExpObj <- SummarizedExperiment(
  assays = list(counts = 2^assay(qf[["Peptides"]])),
  rowData = rowData(qf[["Peptides"]]),
  colData = colData(qf)
)

# Run Normalyzer for evaluation
normalyzer(jobName = "normalyzer",
           experimentObj = SumExpObj,
           sampleColName = "runCol",
           outputDir = results_folder)

rm(SumExpObj)
}

if(Run_normalizer == TRUE || file.exists(normalizer_output)) {
  if(!file.exists(normalizer_output)) {
    stop(paste0("Normalyzer output not found at ", normalizer_output, 
                ". Ensure Normalyzer completed successfully or provide the existing output."))
  }

  # Import normalized data from Normalyzer
  df_norm <- fread(normalizer_output) %>%
    column_to_rownames(var = "Stripped.Sequence") %>%
    dplyr::select(colData(qf)$runCol)

  se <- SummarizedExperiment(assays = list(counts = as.matrix(df_norm)))
  qf <- addAssay(qf, se, name = "PeptidesNorm")
  rowData(qf[["PeptidesNorm"]]) <- rowData(qf[["Peptides"]])
  qf <- addAssayLinkOneToOne(qf, from = "Peptides", to = "PeptidesNorm")

  rm(df_norm, se)
} else {
  qf <- addAssay(qf, qf[["Peptides"]], name = "PeptidesNorm")
  qf <- addAssayLinkOneToOne(qf, from = "Peptides", to = "PeptidesNorm")
}


# (Optional Imputation) ---------------------------------------------------

if(Impute == TRUE) {
  
  if(!is.na(impute_method)) {
    
    qf <- impute(qf, method = impute_method,
                 i = "PeptidesNorm",
                 name = "PeptidesImp")
    
  } else {
    
  randna <- model.Selector(assay(qf[["PeptidesNorm"]]))
  rowData(qf)[["PeptidesNorm"]]$randna <- as.logical(randna[[1]])
  
  qf <- impute(qf, method = "mixed",
               i = "PeptidesNorm",
               name = "PeptidesImp",
               randna = rowData(qf)[["PeptidesNorm"]]$randna,
               mar = impute_MAR,
               mnar = impute_MNAR)

  rm(randna)

  }

  qf <- addAssay(qf, qf[["PeptidesImp"]], name = "PeptidesProccessed")
  qf <- addAssayLinkOneToOne(qf, from = "PeptidesImp", to = "PeptidesProccessed")

}

if(Impute == FALSE) {
  qf <- addAssay(qf, qf[["PeptidesNorm"]], name = "PeptidesProccessed")
  qf <- addAssayLinkOneToOne(qf, from = "PeptidesNorm", to = "PeptidesProccessed")
}


# Optional Batch Correction ---------------------------------------

if(Batch_corr == "eigenms") {
  
  set.seed(12345)
  
  qf <- removeAssay(qf, i = "PeptidesProccessed")
  print("Running EigenMS Batch correction")
   
  if(Impute == TRUE) {
    i_name = "PeptidesImp"
  } else {
    i_name = "PeptidesNorm"
  }
  
  # Setup data
  expr_val       <- assay(qf, i = i_name)
  rowdata        <- as.data.frame(rowData(qf[[i_name]])) %>% 
    relocate(Stripped.Sequence, .before = everything() )
  coldata        <- as.data.frame(colData(qf))
  groups         <- factor(coldata$group)
  names(groups)  <- coldata$group
  
  # Normalize
  TREATS         <- groups
  expr_val_eig1  <- eig_norm1(m = expr_val,
                              treatment = groups,
                              prot.info = rowdata)
  expr_val_norm  <- eig_norm2(rv = expr_val_eig1)
  expr_norm      <- expr_val_norm$normalized
  
  # Add back to QFeatures
  expr <- expr_norm[, coldata$shortname, drop = FALSE]
  rownames(expr) <- rowData(qf[[i_name]])$Stripped.Sequence
  se   <- SummarizedExperiment(
    assays = list(counts = as.matrix(expr)),
    rowData = rowData(qf[[i_name]])
  )

  qf   <- addAssay(qf, se, name = "PeptidesProccessed")
  qf   <- filterFeatures(qf, ~ Stripped.Sequence %in%
                          base::intersect(as.data.frame(rowData(qf[["PeptidesNorm"]]))$Stripped.Sequence,
                                    as.data.frame(rowData(qf[["PeptidesProccessed"]]))$Stripped.Sequence))
  qf <- addAssayLinkOneToOne(qf, from = i_name, to = "PeptidesProccessed")

  rm(i_name, expr_val, rowdata, coldata, groups, TREATS, expr_val_eig1, expr_val_norm, se, expr, expr_norm)
  
}

if(Batch_corr == "sva") {
  
  if (!Impute) {
    stop("Error: Imputation must be performed when using SVA. Set Impute = TRUE.")
  }
  
  qf <- removeAssay(qf, i = "PeptidesProccessed")
  
  set.seed(222)
  
  # Setup models
  mod = model.matrix(~ as.factor(group), data = as.data.frame(colData(qf)) )
  mod0 <- model.matrix(~ 1, data = as.data.frame(colData(qf)))
  edata <- assay(qf[["PeptidesImp"]])
  
  # Run SVA
  n.sv = num.sv(edata, mod, method = "be")
  sva_res <- sva(edata, mod, mod0, n.sv = 3)
  
  # SVA Batch Corr
  res_corrected <- fsva(edata, mod, sva_res, method = "exact", newdat = edata)

  # Add back to QFeatures
  se <- SummarizedExperiment(assays = list(counts = as.matrix(res_corrected$db)))
  qf <- addAssay(qf, se, name = "PeptidesProccessed")
  rowData(qf[["PeptidesProccessed"]]) <- rowData(qf[["PeptidesImp"]])
  qf <- addAssayLinkOneToOne(qf, from = "PeptidesImp", to = "PeptidesProccessed")
  
}

if(Batch_corr == "comBat") {
  
  qf <- removeAssay(qf, i = "PeptidesProccessed")
  
  dat <- assay(qf[["PeptidesImp"]])                    
  pheno <- as.data.frame(colData(qf))
  batch <- pheno$Extraction_batch
  
  m_corr <- ComBat(
    dat = dat,
    batch = batch,
    par.prior = TRUE
  )
    
  # Add back to QFeatures
  se <- SummarizedExperiment(assays = list(counts = as.matrix(m_corr)))
  qf <- addAssay(qf, se, name = "PeptidesProccessed")
  rowData(qf[["PeptidesProccessed"]]) <- rowData(qf[["PeptidesImp"]])
  qf <- addAssayLinkOneToOne(qf, from = "PeptidesImp", to = "PeptidesProccessed")
  
  rm(dat, pheno, batch, m_corr, se)
  
}


# Aggregate to Proteins ---------------------------------------------------

if (!exists("fun_agg")) {
  if(Agg_method == "maxlfq") {
    fun_agg <- maxLFQ_wrapper
  } else if(Agg_method == "robustsummary") {
    fun_agg <- robustSummary_quiet
  } else if(Agg_method == "medianpolish") {
    fun_agg <- MsCoreUtils::medianPolish
  }
}

message("Starting Peptide aggregation...")
qf <- aggregateFeatures(qf,
                        i = "PeptidesProccessed",
                        fcol = "Protein.Group",
                        fun = fun_agg,
                        name = "Proteins")
message("Aggregation Peptide finished.")


# For QC plots ------------------------------------------------------------

# Counts
qf <- countUniqueFeatures(qf,
                          i = "Peptides",
                          colDataName = "Nr_pep")

qf <- countUniqueFeatures(qf,
                          i = "Peptides",
                          groupBy = "Genes",
                          colDataName = "Nr_prot")


# Save QFeatures for DEA --------------------------------------------------

save(qf, file = paste0(saved_data_folder, "Qfeatures_raw.R"))

}




