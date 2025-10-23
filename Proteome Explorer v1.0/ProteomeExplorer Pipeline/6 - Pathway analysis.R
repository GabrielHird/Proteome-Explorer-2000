
##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################


# GSEA analysis -----------------------------------------------------------

Pathway_analysis <- list()

# GSEA --------------------------------------------------------------------

# Data setup

if(Run_GSEA == TRUE) {
  
ranked_data_creator <- function(qf) {
  
  gsea <- qf[["Results"]] %>%
    rowData() %>%
    as.data.frame() %>% 
    dplyr::select(ENTREZID, logFC) %>% 
    na.omit() %>% 
    arrange(desc(logFC)) %>% 
    dplyr::distinct(ENTREZID, .keep_all = TRUE) %>% 
    mutate(named_difference = setNames(logFC, ENTREZID) ) %>% 
    pull(named_difference)
  
}

Pathway_analysis$Ranked_data <- ranked_data_creator(qf)

# Analysis
GSEA_analyser <- function(df, ont = "") {
  
  gseGO(geneList = df,
        ont = "BP",
        keyType = "ENTREZID",
        minGSSize = 2,
        maxGSSize = 1000,
        pvalueCutoff = 0.9,
        verbose = TRUE,
        OrgDb = Org_db,
        pAdjustMethod = "fdr")
}

Pathway_analysis$GSEA_res <- GSEA_analyser(Pathway_analysis$Ranked_data, ont = "BP")


# GSEA_KEGG ---------------------------------------------------------------

KEGG_GSEA <- function(df) {

   gseKEGG(geneList = df,
           organism = "hsa",
           minGSSize = 2,
           maxGSSize = 1000,
           pvalueCutoff = 0.9,
           verbose = TRUE,
           keyType = "kegg",
           pAdjustMethod = "fdr")
}

Pathway_analysis$GSEA_KEGG <- KEGG_GSEA(Pathway_analysis$Ranked_data)

# GSVA KEGG ---------------------------------------------------------------

# Load ExpressionSet
gsva_data <- as.matrix(assay(qf[["Results"]]))
entrez <- as.data.frame(rowData(qf[["Results"]]))
rownames(gsva_data) <- entrez$ENTREZID
rownames(gsva_data) <- make.unique(trimws(rownames(gsva_data)), sep = "_")
smeta <- as.data.frame(colData(qf))

# Load GeneSetCollection
KEGG_db <- download_KEGG("hsa")
KEGG_gs <- split(KEGG_db$KEGGPATHID2EXTID$to, KEGG_db$KEGGPATHID2EXTID$from)

# All
GSVA_param <- gsvaParam(gsva_data,
                        KEGG_gs,
                        minSize = 2,
                        maxSize = 500)

GSVA_res <- gsva(GSVA_param)



# # By group ----------------------------------------------------------------
# 
# # Order samples by group
# ord <- order(smeta$group)
# GSVA_res_ord <- GSVA_res_filt[, ord]
# ann_col_ord  <- ann_col[ord, , drop = FALSE]
# 
# # Heatmap
# pheatmap(GSVA_res_ord,
#          scale = "row",
#          clustering_distance_rows = "euclidean",
#          clustering_method = "complete",
#          cluster_cols = FALSE,    # <- don’t cluster columns
#          annotation_col = ann_col_ord,
#          annotation_colors = list(Group = Group_colors),
#          show_rownames = TRUE,
#          show_colnames = FALSE)
# 
# 
# # Clean
# rm(entrez, KEGG_db)

# GO Terms ----------------------------------------------------------------

# Data setup

go_over_rep_creator <- function(df, type) {
  
  go <- qf[["Results"]] %>% 
    rowData() %>%
    as.data.frame() %>% 
    dplyr::select(ENTREZID, diff_prot) %>% 
    na.omit() %>% 
    dplyr::distinct(ENTREZID, .keep_all = TRUE) %>% 
    filter(diff_prot == type) %>% 
    pull(ENTREZID)
  
}

Pathway_analysis$GO_setup_OVER  <- go_over_rep_creator(qf, type = "OVER")
Pathway_analysis$GO_setup_UNDER <- go_over_rep_creator(qf, type = "UNDER")

# Analysis

go_over_rep <- function(gene_list) {
  
  ego <- enrichGO(gene          = gene_list,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
}

Pathway_analysis$GO_results_over  <- go_over_rep(Pathway_analysis$GO_setup_OVER)
Pathway_analysis$GO_results_under <- go_over_rep(Pathway_analysis$GO_setup_UNDER)

# Analysis (Universe = PExA)

# go_over_rep <- function(gene_list) {
#   
#   ego <- enrichGO(gene          = gene_list,
#                   OrgDb         = org.Hs.eg.db,
#                   universe      = rowData(qf[["Results"]])[["ENTREZID"]],
#                   ont           = "BP",
#                   pAdjustMethod = "BH",
#                   pvalueCutoff  = 0.01,
#                   qvalueCutoff  = 0.05,
#                   readable      = TRUE)
# }
# 
# Pathway_analysis$GO_res_universe_over  <- go_over_rep(Pathway_analysis$GO_setup_UNDER)

# String DB ---------------------------------------------------------------

# Output for loading into Cytoscape

Pathway_analysis$string_table <- as.data.frame(rowData(qf[["Results"]])) %>% 
  dplyr::filter(signif == TRUE) %>%
  dplyr::select(ENTREZID, logFC)
rownames(Pathway_analysis$string_table) <- NULL

write.csv(Pathway_analysis$string_table, file = paste0(results_data_folder, "stringdb_signif.csv"),
          row.names = FALSE)

# Save --------------------------------------------------------------------

save(Pathway_analysis, file = paste0(saved_data_folder, "GSEA.R"))

}

# Load --------------------------------------------------------------------
if(Run_GSEA == FALSE) {
  if(file.exists(paste0(saved_data_folder, "GSEA.R"))) {
    load(paste0(saved_data_folder, "GSEA.R"))
  } else {
    print("No GSEA data found")
  }
}

# PathfinderR -------------------------------------------------------------
if(run_pathfind == TRUE ) {
  
  # Data setup
  df <- qf[["Results"]] %>%
    rowData() %>%
    as.data.frame() %>% 
    dplyr::select(Gene, logFC, qval) %>% 
    na.omit() %>% 
    dplyr::rename(Gene.symbol = Gene,
                  adj.P.Val = qval) 
  
  output_df <- run_pathfindR(df,  p_val_threshold = 0.9, gene_sets = "GO-BP")
  enrichment_chart(output_df)
  
}

rm(gsva_data, GSVA_param, GSVA_res, KEGG_db, KEGG_gs)

