
report_lines <- c(
  "############################################################",
  "#  Gabriel Hirdman Proteome Explorer v1.0® 2025               #",
  "############################################################",
  "",
  "== Version Control ==",
  paste("Project Name:", Project_name),
  paste("Analysis Run:", Analysis_run),
  "",
  "== Importation ==",
  paste("Treatment Group:", Treatment_group),
  paste("Reference Group:", Reference_group),
  paste("Importing qval:", Import_qval),
  paste("Importing pg qval:", Import_pg_qval),
  "",
  "== Preprocessing ==",
  ## NOTE the comma *before* the if() and the corrected variable name
  if (Treshold_level == "Global") {
    paste("Global filtering at:", Global_treshold, "at peptide level =", peptide_level)
  } else if (Treshold_level == "Group") {
    paste("Group filtering at:", Group_treshold, "at peptide level =", peptide_level)
  } else {
    NULL
  },
  paste("Normalization Method:", Norm_method),
  paste("Batch Correction:", Batch_corr),
  paste("Aggregation method:", Agg_method),
  ## another comma before the if()
  if (Impute) {
    paste("Multiple imputation with MAR:", impute_MAR, "and MNAR:", impute_MNAR)
  } else {
    NULL
  },
  "",
  "== DEA Analysis ==",
  paste("DEA Method: ", DEA_method),
  paste("Study design: ", Study_design),
  paste("Alpha (min Q-value):", alpha),
  paste("LogFC cutoff:", logFC_cutoff),
  if (is.na(logFC_cutoff)) {
    "Warning, Bootstrapped LogFC cutoffs"
  } else {
    NULL
  },
  paste("Organism Database:", Org_db),
  "",
  "== Plotting Settings ==",
  "Group Colors:"
)

# Append group colors
for (group in names(Group_colors)) {
  report_lines <- c(report_lines, paste("  ", group, ": ", Group_colors[group], sep = ""))
}

report_lines <- c(report_lines, "Expression Level Colors:")
for (lvl in names(Expression_lvl_color)) {
  report_lines <- c(report_lines, paste("  ", lvl, ": ", Expression_lvl_color[lvl], sep = ""))
}

report_lines <- c(
  report_lines,
  "",
  "== Setup ==",
  paste("DIA-NN Path:", DIA_nn_path),
  paste("FASTA Path:", FASTA_path),
  paste("Sample Data Path:", Sample_data_path),
  paste("Experimenter:", Your_name),
  paste("Lab Name:", Lab_name),
  paste("Contact:", Contact)
)

# -------------------------------------------------
# Write the report content to a text file
# -------------------------------------------------

report_file <- paste0(results_folder,"analysis_report.txt")
writeLines(report_lines, con = report_file)
cat("Analysis report successfully created in", report_file, "\n")
