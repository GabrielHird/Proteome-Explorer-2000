##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

# Special Subtype Analysis

# Update colData ----------------------------------------------------------

qf$subtype_contrast <- dplyr::case_when(
  qf$LUAD == 1 ~ "LUAD",
  qf$SCC  == 1 ~ "SCC",
  .default = NA_character_
)

qf <- qf[, !is.na(qf$subtype_contrast), ]


# Run Msqrob --------------------------------------------------------------

colData(qf)$subtype_contrast    <- relevel(factor(colData(qf)$subtype_contrast), ref = Reference_group)
colData(qf)$Bio_repl <- factor(colData(qf)$Bio_repl)

if(Study_design == "unpaired") {
  formula_model <- ~ subtype_contrast
} else if (Study_design == "paired") {
  formula_model <- ~ subtype_contrast + (1 | Bio_repl)
}

df <- msqrob2::msqrob(object = qf,
                      i = "Proteins",
                      formula = formula_model,
                      overwrite = TRUE)

cont1 <- paste0("subtype_contrast",Treatment_group,"=0")
cont2 <- paste0("subtype_contrast",Treatment_group)

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

# Final Cleanup -----------------------------------------------------------

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