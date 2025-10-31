
# Aggregate via MaxLFQ ----------------------------------------------------

# Old
maxLFQ_wrapper <- function(x) {
  iq::maxLFQ(x)$estimate
}

robustSummary_quiet <- function(...) {
  withCallingHandlers(
    MsCoreUtils::robustSummary(...),
    warning = function(w) {
      if (grepl("rlm' failed to converge", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Create an ExpressionSet -------------------------------------------------

createExpressionSet <- function(df, Sample_data) {
  
  # Expression data: extract LFQ columns and convert to matrix
  LFQ_columns <- which(colnames(df) %in% Sample_data$sample_id)
  exprs <- as.matrix(df[, ..LFQ_columns])
  rownames(exprs) <- df$Stripped.Sequence
  
  # Phenotype data: reformat and order by sample IDs
  pData <- as.data.frame(Sample_data)
  rownames(pData) <- pData$sample_id
  pData <- pData[colnames(exprs), ]
  pData <- new("AnnotatedDataFrame", data = pData)
  
  # Feature data: remove LFQ columns and set rownames
  fData <- as.data.frame(df)
  fData <- fData[, -LFQ_columns]
  rownames(fData) <- fData$Stripped.Sequence
  fData <- new("AnnotatedDataFrame", data = fData)
  
  # Experiment data: create MIAME object with metadata
  eData <- new("MIAME",
               name = Your_name,
               lab = Lab_name,
               contact = Contact,
               title = Project_name)
  
  # Assemble and return the ExpressionSet
  eset <- ExpressionSet(
    assayData = exprs,
    phenoData = pData,
    featureData = fData,
    experimentData = eData
  )
  
  return(eset)
}


# Join QFeatures assays ---------------------------------------------------

joinDIANNAssays <- function(qf, name = "Precursors") {
  
  cat("Checking for duplicates\n")
  
  ## Check if any Precursor.Id is duplicated
  anyDups <- sapply(names(qf), function(i) {
    any(duplicated(rowData(qf)[[i]]$Precursor.Id))
  })
  
  if (any(anyDups)) {
    stop("Error: duplicated Precursor.Id in data, unable to join")
  }
  
  cat("No duplicates: updating rownames\n")
  
  ## Set rownames
  for (i in names(qf)) {
    rownames(qf[[i]]) <- rowData(qf[[i]])$Precursor.Id
  }
  
  cat("Joining assays into -> Precursors\n")
  
  ## Join
  qf <- joinAssays(qf,
                   i = names(qf),
                   name = name)
  
  cat("All done!\n")
  
  return(qf)
  
}


# Count unique ------------------------------------------------------------

countUniqueDIA <- function(qf) {
  all_sequences <- character()  # Initialize an empty character vector
  for (assay_name in names(qf)) {
    assay <- qf[[assay_name]]
    if ("Stripped.Sequence" %in% colnames(rowData(assay))) {
      sequences <- rowData(assay)$Stripped.Sequence
      all_sequences <- c(all_sequences, sequences)
    }
  }
  length(unique(all_sequences))
}


# Function to filter for the canonical protein per gene based on P --------

filterCanonical <- function(q, assay, fasta) {
  
  # 1) FASTA -> accession rank + PE (for tie-break)
  fa  <- Biostrings::readAAStringSet(fasta); hdr <- names(fa)
  rank_map <- tibble(
    acc  = str_match(hdr, "^[a-z]{2}\\|([^|]+)\\|")[,2],
    db   = substr(hdr, 1, 2),
    name = str_match(hdr, "^[a-z]{2}\\|[^|]+\\|([^\\s]+)")[,2],
    pe   = suppressWarnings(as.integer(str_match(hdr, "PE=(\\d+)")[,2])),
    frag = str_detect(hdr, "\\(Fragment\\)"),
    iso  = str_detect(name, "-\\d+$"),
    score = case_when(
      db == "sp" & !frag & !iso ~ 4L,
      db == "sp" & !frag        ~ 3L,
      db == "sp"                ~ 2L,
      TRUE                      ~ 1L
    )
  ) %>% transmute(acc, score = coalesce(score, 0L), pe = coalesce(pe, 999L))
  
  # 2) RowData -> long: (feature_id, Gene, acc)
  rd_long <- SummarizedExperiment::rowData(q[[assay]]) %>%
    as.data.frame() %>%
    rownames_to_column("feature_id") %>%
    mutate(Gene = coalesce(Genes, Gene)) %>%
    filter(!is.na(Gene), Gene != "") %>%
    separate_rows(Protein.Group, sep = ";\\s*") %>%
    mutate(acc = if_else(str_detect(Protein.Group, "^[a-z]{2}\\|"),
                         str_match(Protein.Group, "^[a-z]{2}\\|([^|]+)\\|")[,2],
                         Protein.Group))
  
  # 3) Pick winner accession per Gene (by score, then PE)
  winners <- rd_long %>%
    distinct(Gene, acc) %>%
    left_join(rank_map, by = "acc") %>%
    mutate(score = coalesce(score, 0L), pe = coalesce(pe, 999L)) %>%
    arrange(Gene, desc(score), pe) %>%
    group_by(Gene) %>% slice(1) %>% ungroup() %>%
    transmute(Gene, best_acc = acc)
  
  # 4) Keep all peptide rows that include the winner accession for their Gene
  keep_ids <- rd_long %>%
    inner_join(winners, by = c("Gene", "acc" = "best_acc")) %>%
    distinct(feature_id) %>% pull()
  
  q[keep_ids, , ]
}

# Convert QFeatures to eset -----------------------------------------------

QFeattoES <- function(qf, assayName) {
  se <- qf[[assayName]]
  colData(se) <- colData(qf)
  assayNames(se) <- "exprs"
  as(se, "ExpressionSet")
}


# Modify MSDAP function to import eset ------------------------------------
import_eset_msdap <- function(eset,
                              column_fdata_protein_id = "Protein.Group",
                              column_pdata_sample_group = "group", 
                              acquisition_mode = "dia",
                              is_log2 = TRUE) {
  reset_log()
  
  if (!acquisition_mode %in% c("dda", "dia")) {
    append_log("`acquisition_mode` parameter can only be 'dda' or 'dia'", type = "error")
  }
  
  if (!"ExpressionSet" %in% names(Biobase::classVersion(eset))) {
    append_log("This function requires a dataset of type: ExpressionSet", type = "error")
  }
  
  if (!column_pdata_sample_group %in% colnames(Biobase::pData(eset))) {
    cat("Column names in ExpressionSet pData:", paste(colnames(Biobase::pData(eset)), collapse = ", "), "\n")
    append_log(paste("ExpressionSet pData does not contain the sample information column you requested:", 
                     column_pdata_sample_group), type = "error")
  }
  
  if (!column_fdata_protein_id %in% colnames(Biobase::fData(eset))) {
    cat("Column names in ExpressionSet fData:", paste(colnames(Biobase::fData(eset)), collapse = ", "), "\n")
    append_log(paste("ExpressionSet fData does not contain the protein identifier column you requested:", 
                     column_fdata_protein_id), type = "error")
  }
  
  # Process sample data
  samples <- Biobase::pData(eset)
  samples$sample_id <- rownames(samples)
  samples <- as_tibble(samples) %>%
    mutate_all(as.character) %>%
    rename(group = !!column_pdata_sample_group)
  
  if ("exclude" %in% colnames(samples)) {
    samples <- samples %>% mutate(exclude = as.logical(exclude))
  } else {
    samples <- add_column(samples, exclude = FALSE)
  }
  
  if (!"shortname" %in% colnames(samples)) {
    samples$shortname <- samples$sample_id
  }
  
  # Process protein data
  proteins <- Biobase::fData(eset)
  proteins <- proteins[!duplicated(proteins[, column_fdata_protein_id]), ]
  proteins <- as_tibble(proteins) %>% mutate_all(as.character) %>%
    rename(protein_id = !!column_fdata_protein_id)
  
  if (!"gene_symbols_or_id" %in% colnames(proteins)) {
    proteins <- add_column(proteins, gene_symbols_or_id = proteins$protein_id)
  }
  
  if (!"fasta_headers" %in% colnames(proteins)) {
    proteins$fasta_headers <- proteins$protein_id
  }
  
  # Process peptide data: pivot the expression matrix to long format
  tib <- as_tibble(Biobase::exprs(eset)) %>%
    add_column(peptide_id = Biobase::featureNames(eset),
               protein_id = Biobase::fData(eset)[, column_fdata_protein_id]) %>%
    pivot_longer(cols = -c(peptide_id, protein_id), 
                 names_to = "sample_id", 
                 values_to = "intensity") %>%
    filter(is.finite(intensity) & intensity > 0)
  
  # Initialize default sequence columns (will be overridden if join finds a value)
  tib <- tib %>% 
    mutate(sequence_plain = peptide_id,
           sequence_modified = peptide_id)
  
  # Merge plain sequence if available
  col_plainseq <- grep("^sequence$|^seq$|(base|plain)[._-]*(sequence|seq)|(sequence|seq)[._-]*(base|plain)", 
                       colnames(Biobase::fData(eset)), ignore.case = TRUE, value = TRUE)
  if (length(col_plainseq) > 0) {
    tib <- tib %>% 
      left_join(
        as_tibble(Biobase::fData(eset)) %>%
          mutate(peptide_id = Biobase::featureNames(eset)) %>%
          select(peptide_id, plain_seq = !!sym(col_plainseq[1])),
        by = "peptide_id"
      ) %>%
      mutate(sequence_plain = coalesce(plain_seq, sequence_plain)) %>%
      select(-plain_seq)
    append_log(paste("Found plain peptide sequences in ExpressionSet fData, column:", col_plainseq[1]), type = "info")
  }
  
  # Merge modified sequence if available
  col_modseq <- grep("(modified|mod)[._-]*(sequence|seq)|(sequence|seq)[._-]*(modified|mod)", 
                     colnames(Biobase::fData(eset)), ignore.case = TRUE, value = TRUE)
  if (length(col_modseq) > 0) {
    tib <- tib %>% 
      left_join(
        as_tibble(Biobase::fData(eset)) %>%
          mutate(peptide_id = Biobase::featureNames(eset)) %>%
          select(peptide_id, mod_seq = !!sym(col_modseq[1])),
        by = "peptide_id"
      ) %>%
      mutate(sequence_modified = coalesce(mod_seq, sequence_modified)) %>%
      select(-mod_seq)
    append_log(paste("Found modified peptide sequences in ExpressionSet fData, column:", col_modseq[1]), type = "info")
  }
  
  # Merge retention time if available
  col_rt <- grep("retention[_ ]?time|rt", colnames(Biobase::fData(eset)), ignore.case = TRUE, value = TRUE)
  if (length(col_rt) > 0) {
    tib <- tib %>% 
      left_join(
        as_tibble(Biobase::fData(eset)) %>%
          mutate(peptide_id = Biobase::featureNames(eset)) %>%
          select(peptide_id, rt = !!sym(col_rt[1])),
        by = "peptide_id"
      )
    append_log(paste("Found retention time (rt) in ExpressionSet fData, column:", col_rt[1]), type = "info")
  } else {
    append_log("Warning: No retention time (rt) found in ExpressionSet. Defaulting to NA.", type = "warning")
    tib <- tib %>% mutate(rt = NA)
  }
  
  if (!is_log2) {
    tib <- tib %>% mutate(intensity = log2(intensity))
  }
  tib <- tib %>% 
    mutate(intensity = if_else(!is.na(intensity) & intensity < 0, 0, intensity),
           detect = TRUE,
           isdecoy = FALSE)
  
  # Create the dataset object
  dataset <- list(samples = samples, proteins = proteins, peptides = tib, acquisition_mode = acquisition_mode)
  
  # --- Incorporate peptide and protein count calculations ---
  dataset <- invalidate_cache(dataset)
  dataset$samples <- as_tibble(dataset$samples)
  dataset$samples <- msdap:::peptide_and_protein_counts_per_sample(dataset$peptides, dataset$samples, is_dia_dataset(dataset))
  
  # Add the empty plots list
  dataset$plots <- list()
  
  return(dataset)
}

# Plotting functions ------------------------------------------------------

plot_densities <- function(df, assay = "peptideRaw", title = "") {
  
  df[[assay]] %>%
    assay %>%
    as.data.frame() %>%
    gather(sample, intensity) %>%
    ggplot(aes(x = intensity, group = sample, color = sample)) +
    geom_density() +
    ggtitle(title)
  
}


# Filter function -----------------------------------------------------------------

flagPeptidesForRemoval <- function(qf,
                                   peptide_assay = "Peptides",
                                   group_var = "group",
                                   threshold = 0.75,
                                   protein_col = "Protein.Group",
                                   peptide_level = FALSE) {
  
  if (!requireNamespace("QFeatures", quietly = TRUE)) {
    stop("Package 'QFeatures' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!inherits(qf, "QFeatures")) {
    stop("Input 'qf' must be a QFeatures object.")
  }
  if (!peptide_assay %in% names(qf)) {
    stop("Peptide assay '", peptide_assay, "' not found in QFeatures object.")
  }
  if (!is.logical(peptide_level) || length(peptide_level) != 1) {
    stop("'peptide_level' must be a single logical value (TRUE or FALSE).")
  }
  
  quant <- assay(qf[[peptide_assay]])
  cd <- colData(qf)
  
  if (!(group_var %in% colnames(cd))) {
    stop("Group variable '", group_var, "' not found in colData.")
  }
  grp <- cd[[group_var]]
  
  if (is.null(grp)) {
    stop("Group variable '", group_var, "' resulted in NULL.")
  }
  if (anyNA(grp)) {
    warning("Group variable '", group_var, "' contains NA values. Samples with NA group will be ignored.")
    valid_samples <- !is.na(grp)
    quant <- quant[, valid_samples, drop = FALSE]
    grp <- grp[valid_samples]
    if (length(grp) == 0) {
      stop("No valid samples remaining after removing NAs from group variable.")
    }
  }
  
  grp <- factor(grp)
  groups <- levels(grp)
  
  if (length(groups) == 0) {
    stop("Group variable '", group_var, "' has no levels after processing.")
  }
  
  present_mat <- sapply(groups, function(g) {
    cols <- which(grp == g)
    if (length(cols) == 0) {
      return(rep(0, nrow(quant)))
    }
    rowMeans(!is.na(quant[, cols, drop = FALSE])) * 100
  })
  
  if (is.vector(present_mat)) {
    present_mat <- matrix(present_mat, ncol = 1)
  }
  colnames(present_mat) <- paste0("present_", groups)
  
  present75 <- apply(present_mat, 1, function(x) any(x / 100 >= threshold, na.rm = TRUE))
  
  rd <- DataFrame(rowData(qf[[peptide_assay]]))
  rd$present75 <- present75
  rd <- cbind(rd, DataFrame(present_mat))
  
  if (peptide_level) {
    message("Performing peptide-level filtering (protein grouping ignored).")
    rd$removePeptide <- !rd$present75
  } else {
    message("Performing protein-level filtering.")
    if (!(protein_col %in% colnames(rd))) {
      stop("Protein grouping column '", protein_col, "' not found in rowData. ",
           "Needed because peptide_level = FALSE.")
    }
    if (anyNA(rd[[protein_col]])) {
      warning("Protein grouping column '", protein_col, "' contains NA values. ",
              "Peptides with NA protein group might not be handled as expected.")
    }
    
    rd$protein_ok <- ave(rd$present75, rd[[protein_col]], FUN = function(x) any(x, na.rm = TRUE))
    rd$removePeptide <- !rd$protein_ok
  }
  
  rowData(qf[[peptide_assay]]) <- rd
  
  return(qf)
}
# Fixing issues -----------------------------------------------------------

fix_multiple_pep <- function(qf, i = "") {
  
  df <- as.data.frame(rowData(qf[[i]]))
  
  # Identify flagged peptides that match to multiple protein groups
  flagged_peptides <- df %>%
    group_by(Stripped.Sequence) %>%
    summarize(n_unique = n_distinct(Protein.Group),
              .groups = "drop") %>%
    filter(n_unique > 1) %>%
    pull(Stripped.Sequence)
  
  # Mutate the data frame: for flagged peptides, append a unique suffix to each occurrence
  df <- df %>%
    group_by(Stripped.Sequence) %>%
    mutate(Stripped.Sequence = if (dplyr::first(Stripped.Sequence) %in% flagged_peptides) {
      paste0(Stripped.Sequence, "_", row_number())
    } else {
      Stripped.Sequence
    }) %>%
    ungroup()
  
  
  rowData(qf[[i]]) <- DataFrame(df)
  
  return(qf)
  
}


# For Bootstrapping -------------------------------------------------------

permute_ab <- function(x, nmax = 100) {
  # Ensure x is an integer vector with only 1s and 2s.
  if (!is.integer(x)) {
    stop("x must be an integer vector of group IDs (1,2)")
  }
  if (sum(x == 1) < 2 || sum(x == 2) < 2) {
    stop("x must contain at least two of each group")
  }
  
  n1 <- sum(x == 1)
  n2 <- sum(x == 2)
  n <- length(x)
  
  # Calculate total number of combinations
  total_combs <- choose(n, n1)
  
  if (total_combs <= nmax) {
    # Generate all combinations if the total is manageable.
    all_combs <- combn(n, n1)  # each column is one possible combination
  } else {
    # Otherwise, randomly sample nmax unique combinations.
    all_combs <- matrix(NA, nrow = n1, ncol = nmax)
    seen <- character(0)
    count <- 0
    while (count < nmax) {
      candidate <- sort(sample(n, n1))
      key <- paste(candidate, collapse = "-")
      if (!(key %in% seen)) {
        count <- count + 1
        seen <- c(seen, key)
        all_combs[, count] <- candidate
      }
    }
  }
  
  # Construct the permutation matrix.
  # For each combination, positions specified in the combination get indices 1:n1 (group1)
  # and the rest get indices n1+1:n (group2).
  perm_matrix <- matrix(0L, nrow = ncol(all_combs), ncol = n)
  for (i in seq_len(ncol(all_combs))) {
    p <- integer(n)
    pos1 <- all_combs[, i]            # positions chosen for group 1
    p[pos1] <- seq_len(n1)            # assign group1 indices
    p[-pos1] <- seq(n1 + 1, n)         # assign group2 indices
    perm_matrix[i, ] <- p
  }
  
  return(perm_matrix)
}

bootstrap_log2FC_threshold <- function(qf, 
                                       Reference_group, 
                                       Treatment_group, 
                                       nmax = 1000, 
                                       probs = 0.95, 
                                       seed = 123, 
                                       return_fc_matrix = FALSE) {
  # Extract the SummarizedExperiment from your QFeatures object.
  se <- qf[["Proteins"]]
  colData(se) <- colData(qf)
  
  # Ensure a 'group' column exists and relevel so that the reference group comes first.
  if (!"group" %in% colnames(colData(se))) {
    stop("The SummarizedExperiment must have a 'group' column in its colData.")
  }
  se$group <- relevel(factor(se$group), ref = Reference_group)
  
  # Verify that there are exactly two groups.
  if (length(levels(se$group)) != 2) {
    stop("Exactly two groups are required for this permutation approach.")
  }
  
  # Create an integer vector for groups: 1 for Reference, 2 for Treatment.
  group_int <- ifelse(se$group == Reference_group, 1L, 2L)
  
  # Get the assay matrix (assumed to be log2-transformed).
  x <- assay(se)
  
  # Get sample names; assume the column order in x corresponds to the order of group_int.
  sample_names <- colnames(x)
  
  # Count samples in each group.
  n1 <- sum(group_int == 1)
  n2 <- sum(group_int == 2)
  
  # Generate the permutation matrix using our robust permute_ab function.
  set.seed(seed)
  perm_matrix <- permute_ab(as.integer(group_int), nmax = nmax)
  nperm <- nrow(perm_matrix)
  
  # Allocate a matrix to hold fold-change values for each protein across permutations.
  fc_matrix <- matrix(NA, nrow = nrow(x), ncol = nperm)
  
  # For each permutation, reassign sample labels and compute fold changes.
  for (i in 1:nperm) {
    perm_indices <- perm_matrix[i, ]  # permutation indices for this iteration
    
    # Reorder sample names according to the permutation.
    samples_permuted <- sample_names[perm_indices]
    
    # Split samples: first n1 are group1 (Reference), remaining n2 are group2 (Treatment).
    sample_id_cond1 <- samples_permuted[1:n1]
    sample_id_cond2 <- samples_permuted[(n1 + 1):(n1 + n2)]
    
    # Compute fold change for each protein as the difference in row means.
    x1 <- rowMeans(x[, sample_id_cond1, drop = FALSE], na.rm = TRUE)
    x2 <- rowMeans(x[, sample_id_cond2, drop = FALSE], na.rm = TRUE)
    fc_matrix[, i] <- x1 - x2
  }
  
  # Determine the threshold as the maximum absolute quantile from the distribution.
  threshold <- max(abs(quantile(fc_matrix, probs = c(1 - probs, probs), na.rm = TRUE)))
  
  if (return_fc_matrix) {
    return(list(threshold = threshold, fc_matrix = fc_matrix))
  } else {
    return(threshold)
  }
}

