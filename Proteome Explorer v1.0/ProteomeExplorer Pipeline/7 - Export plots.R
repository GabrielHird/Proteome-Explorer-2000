##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

# Save plots to file

if(save_plots == TRUE) {

# Normalisation
p <- Norm_boxplotter(qf, i = "PeptidesProccessed")
ggsave(paste0(results_folder,"boxplott.png"), p, width = 800/72, height = 600/72, dpi = 72)

# Volcano
p <- Volcano_plotter(qf)
p <- p +
  scale_x_continuous(limits = c(-0.75, 0.75)) +
  scale_y_continuous(limits = c(0, 12))
ggsave(paste0(results_folder,"volcano.png"), p, width = 440/72, height = 600/72, dpi = 72)

# Heatmap
png(paste0(results_folder, "Heatmap.png"),width=600/72,height=600/72,units="in",res=1200)
print(make_heatmap(qf, i = "Results", annot_list = c(group = "group"), only_significant = TRUE))
dev.off()

# PCA
p <- PCA_plotter(qf, i = "Proteins",color = "group", label = "none")
ggsave(paste0(results_folder,"PCA.png"), p, width = 530/72, height = 360/72, dpi = 72)

sign_counts <- as.data.frame(rowData(qf[["Results"]])) %>%
  dplyr::filter(diff_prot %in% c("OVER", "UNDER")) %>%
  dplyr::count(diff_prot) %>% 
  dplyr::mutate(TOTAL = sum(n))

}

# Save Excel

if(export_data == TRUE) {
  
  # Stats
  df1 <- as.data.frame(rowData(qf[["Results"]]))
  openxlsx::write.xlsx(df1, file = paste0(results_data_folder,"Stats_",Analysis_run,".xlsx"),
                       keepNA = FALSE, rowNames = FALSE)
  
  # Sample data
  df2 <- as.data.frame(colData(qf)) 
  openxlsx::write.xlsx(df2, file = paste0(results_data_folder,"Samples_",Analysis_run,".xlsx"),
                       keepNA = FALSE, rowNames = FALSE)
  
  # ExpressionMatrix
  df3 <- assay(qf[["Results"]]) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Protein.Group") %>% 
    left_join(dplyr::select(df1, Protein.Group, Gene), by = "Protein.Group") %>% 
    relocate(Gene, .after = Protein.Group)
  openxlsx::write.xlsx(df3, file = paste0(results_data_folder,"Protein_expression_",Analysis_run,".xlsx"),
                       keepNA = FALSE, rowNames = FALSE)
  
  # OmicLearn
  df4 <- assay(qf[["Results"]]) %>% 
    t() %>% 
    as.data.frame() %>%
    rownames_to_column(var = "shortname") %>%
    left_join(dplyr::select(as.data.frame(colData(qf)), shortname, group), by = "shortname") %>% 
    dplyr::rename(`_group` = group) %>% 
    dplyr::select(-shortname)
  openxlsx::write.xlsx(df4, file = paste0(results_data_folder,"OmicLearn",Analysis_run,".xlsx"),
                       keepNA = FALSE, rowNames = FALSE)
  
  # GOAT
  df5 <- qf[["Results"]] %>% 
    rowData() %>% 
    as.data.frame() %>% 
    dplyr::select(Gene, ENTREZID, qval, logFC, signif) %>% 
    rename(
      gene = ENTREZID,
      symbol = Gene,
      pvalue = qval,
      effectsize = logFC
    ) %>% 
    dplyr::relocate(gene, .before = symbol)
  rownames(df5) <- NULL
  openxlsx::write.xlsx(df5, file = paste0(results_data_folder,"GOAT_",Analysis_run,".xlsx"),
                       keepNA = FALSE, rowNames = FALSE)
  
  rm(df1, df2, df3, df4, df5)
  
  # Pathway analysis
  f <- purrr::possibly(~ .x@result, otherwise = NULL)
  pr <- purrr::compact(list(
    "GO Upregulated"   = f(Pathway_analysis[["GO_results_over"]]),
    "GO Downregulated" = f(Pathway_analysis[["GO_results_under"]]),
    "GSEA GO Terms"    = f(Pathway_analysis[["GSEA_res"]]),
    "GSEA KEGG"        = f(Pathway_analysis[["GSEA_KEGG"]])
  ))
  
  if(length(pr) > 0) {
    openxlsx::write.xlsx(pr,
                         file = paste0(results_data_folder, "Pathway_",Analysis_run,".xlsx"),
                         keepNA = FALSE, rowNames = FALSE)
  }
  rm(pr, f)
  
}

