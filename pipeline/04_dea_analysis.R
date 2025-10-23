
##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

# Load QFeatures for DEA --------------------------------------------------

if(Run_DEA == TRUE & First_pass == FALSE) {
  print("Loading previously saved QFeatures object")
  load(paste0(saved_data_folder, "Qfeatures_raw.R"))
}


# Bootstrapping log2 FC cutoff --------------------------------------------

if(Run_DEA == TRUE) {
  if (is.na(logFC_cutoff)) {
   logFC_cutoff <- bootstrap_log2FC_threshold(qf, 
                                  Reference_group = Reference_group, 
                                  Treatment_group = Treatment_group, 
                                  nmax = 1000, 
                                  probs = 0.95, 
                                  seed = 123, 
                                  return_fc_matrix = FALSE)
   logFC_cutoff <- round(logFC_cutoff,2)
  }

  save(logFC_cutoff, file = paste0(saved_data_folder, "logFC_cutoff.R") )
} else {
  load(file= paste0(saved_data_folder, "logFC_cutoff.R"))
}

# ProDA -------------------------------------------------------------------

if(DEA_method == "proDA" & Run_DEA == TRUE) {
  
  # Prep data
  se <- qf[["Proteins"]]
  colData(se) <- colData(qf)
  se$group <- relevel(factor(se$group), ref = Reference_group)
  cont <- paste0("group", Treatment_group)
  
  # Prep the fit
  fit <- proDA(data = se,
               design = ~ group)
  
  # Extract results
  results <- test_diff(fit, cont, sort_by = "pval")
  
  # Standardize
  results <- results %>%
    as_tibble() %>%
    dplyr::select(name, pval, adj_pval, diff) %>% 
    dplyr::rename(
      Protein.Group = name,
      pval = pval,
      qval = adj_pval,
      logFC = diff
    )
  
  # Cleanup
  rm(fit,se)
  
}

# Limma -------------------------------------------------------------------

if(DEA_method == "limma" & Run_DEA == TRUE) {
  
  # Prep data
  se <- qf[["Proteins"]]
  colData(se) <- colData(qf)
  cont <- paste0("se$group", Treatment_group)
  
  # Prepare Data
  se$group <- relevel(factor(se$group), ref = Reference_group)
  model_design <- model.matrix(~ se$group)
  
  # Prep the fit
  fit <- se %>% 
    assay() %>% 
    lmFit(design = model_design)
  fit <- eBayes(fit = fit,
                trend = TRUE)
  
  # Extract results
  results <- topTable(fit = fit,
                      coef = cont,
                      adjust.method = "BH",
                      number = Inf) %>% 
    rownames_to_column("Protein")
  
  # Standardize
  results <- results %>%
    as_tibble() %>% 
    dplyr::select(Protein, logFC, P.Value, adj.P.Val) %>% 
    dplyr::rename(
      Protein.Group = Protein,
      pval = P.Value,
      qval = adj.P.Val,
      logFC = logFC
    )
  
  # Cleanup
  rm(fit,se,model_design)
}

# DEqMS -------------------------------------------------------------------

if(DEA_method == "DEqMS" & Run_DEA == TRUE) {
  
  # Prep data
  se <- qf[["Proteins"]]
  colData(se) <- colData(qf)
  cont <- paste0("group", Treatment_group,"-","group",Reference_group)
  
  # Prepare Data
  group <- relevel(factor(se$group), ref = Reference_group)
  model_design <- model.matrix(~0 + group)
  
  # Prep the fit
  fit1 <- se %>% 
    assay() %>% 
    lmFit(design = model_design)
  cont <- makeContrasts(cont, levels = model_design)
  fit2 = contrasts.fit(fit1, contrasts = cont)
  fit3 <- eBayes(fit2)
  fit3$count <- rowData(se)[rownames(fit3$coefficients), ".n"]
  fit4 = spectraCounteBayes(fit3)
  
  # Extract results
  results = outputResult(fit4, coef_col = 1)
  
  # Standardize
  results <- results %>% 
    rownames_to_column(var = "Protein.Group") %>% 
    as_tibble () %>% 
    dplyr::select(Protein.Group, logFC, sca.P.Value, sca.adj.pval ) %>% 
    dplyr::rename(
      pval = sca.P.Value,
      qval = sca.adj.pval,
      logFC = logFC
    )
  
  # Cleanup
  rm(se, group, model_design, fit1, cont, fit2, fit3, fit4)
}

# MSqRob ------------------------------------------------------------------

if(DEA_method == "msqrob" & Run_DEA == TRUE) {
  
  # Prepare Data
  
  colData(qf)$group    <- relevel(factor(colData(qf)$group), ref = Reference_group)
  colData(qf)$Bio_repl <- factor(colData(qf)$Bio_repl)
  
  if(Study_design == "unpaired") {
    formula_model <- ~ group
  } else if (Study_design == "paired") {
    formula_model <- ~ group + (1 | Bio_repl)
  }
  
  df <- msqrob2::msqrob(object = qf,
                        i = "Proteins",
                        formula = formula_model,
                        overwrite = TRUE)
  
  cont1 <- paste0("group",Treatment_group,"=0")
  cont2 <- paste0("group",Treatment_group)
  
  L <- makeContrast(cont1, parameterNames = c(cont2))
  
  # Prep the fit 
  df <- msqrob2::hypothesisTest(object = df,
                                i = "Proteins",
                                contrast = L,
                                overwrite = TRUE)
  
  # Standardize
  results <- rowData(df[["Proteins"]]) %>% 
    as_tibble() %>% 
    dplyr::select(Protein.Group, ends_with(".logFC"), ends_with(".pval"), ends_with(".adjPval")) %>% 
    dplyr::rename(
      pval = ends_with(".pval"),
      qval = ends_with(".adjPval"),
      logFC = ends_with(".logFC")
    )
  
  # Cleanup
  rm(df, L)
  
}

# MS-EmpiRe ---------------------------------------------------------------

if(DEA_method == "msempire" & Run_DEA == TRUE) {
  
  # Setup data
  es <- QFeattoES(qf, assayName = "PeptidesProccessed")
  msdap <- import_eset_msdap(es)
  msdap = setup_contrasts(msdap, contrast_list = list(c(Reference_group, Treatment_group)) )
  msdap = dea(msdap, dea_algorithm = "msempire")
  
  # Standardize
  results <- msdap$de_proteins %>% 
    as_tibble() %>% 
    dplyr::select(protein_id, pvalue, qvalue, foldchange.log2) %>% 
    dplyr::rename(
      Protein.Group = protein_id,
      pval = pvalue,
      qval = qvalue,
      logFC = foldchange.log2
    )
  rm(es, msdap)
  
}

# Final cleanup -----------------------------------------------------------

if(Run_DEA == TRUE) {
  
r <- as.data.frame(rowData(qf[["Proteins"]]))
r <- clusterProfiler::bitr(r$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = Org_db)
r <- r[!duplicated(r$SYMBOL), ]
  
results <- results %>% 
  inner_join(dplyr::select(as.data.frame(rowData(qf[["Proteins"]])),
                           Protein.Group,
                           Protein.Names,
                           Genes,
                           Gene,
                           First.Protein.Description,
                           .n),
             by = "Protein.Group") %>% 
  left_join(r, by = c("Gene" = "SYMBOL")) %>% 
  dplyr::rename(Number.Peptides = .n) %>% 
  relocate(logFC, pval, qval, .after = last_col()) %>% 
  mutate(
         log10_qval = -1*log10(qval),
         logFC = round(logFC, 3),
         qval  = round(qval, 5),
         signif = if_else(qval < 0.05 & abs(logFC) > logFC_cutoff, TRUE, FALSE),
         diff_prot = case_when(
           signif == "TRUE" & logFC >= 0 ~ "OVER",
           signif == "TRUE" & logFC <= 0 ~ "UNDER",
           signif == "FALSE" ~ "FALSE"),
         label = case_when(
           qval <= 0.001 & signif == "TRUE" ~ "***",
           qval <= 0.01 & signif == "TRUE" ~ "**",
           qval <= 0.05 & signif == "TRUE" ~ "*",
           TRUE ~ "NS") ) %>% 
  relocate(log10_qval, .after = qval) 

rm(r)

# Save the data to qf
qf <- addAssay(qf, qf[["Proteins"]], name = "Results")
qf <- addAssayLinkOneToOne(qf, from = "Proteins", to = "Results")
rownames(results) <- results$Protein.Group
if(!all(rownames(qf[["Results"]]) %in% rownames(results))){
  stop("Some feature IDs in the current data are missing in the new feature data.")
}
rowData(qf[["Results"]]) <- results[rownames(qf[["Results"]]), ]
save(qf, file = paste0(saved_data_folder, "DEA_results.R"))

rm(results)

}


