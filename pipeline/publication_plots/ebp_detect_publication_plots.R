##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

# Publication Level Plots

# PCA Plot Batches --------------------------------------------------------

pca_batch_plot <- function(qf,i,group, title) {
  
  # Data
  X <- t(assay(qf[[i]], i = 1, withDimnames = TRUE))
  y <- factor( colData(qf)[[group]] )
  # Colors
  n.grp <- length(levels(y))
  pal <- brewer.pal(12, "Set3")[1:n.grp]
  # PCA
  res.pca <- mixOmics::pca(X)
  # Plot
  plotIndiv(res.pca,
            group = y,
            col.per.group = pal,
            style = "lattice",
            ellipse = TRUE,
            ind.names = FALSE,
            pch = 19,
            cex = 1,
            legend = TRUE,
            legend.title = "Batch",
            title = title)
  # Save
  dev.print(
    device = pdf,
    file = paste0(results_folder, title, ".pdf"),
    width = pix(490),
    height = pix(420)
  )
  
}

# Plot
pca_batch_plot(qf, i = "Peptides", group = "Extraction_batch", title = "Pre Batch correction")
pca_batch_plot(qf, i = "PeptidesProccessed", group = "Extraction_batch", title = "Post Batch correction")


# PCA All -----------------------------------------------------------------

pca_all <- function(qf,i, title) {
  
  # Data
  X <- t(assay(qf[[i]], i = 1, withDimnames = TRUE))
  y <- factor( colData(qf)$group)
  # PCA
  res.pca <- mixOmics::pca(X)
  # Plot
  plotIndiv(res.pca,
            group = y,
            col.per.group = Group_colors,
            style = "lattice",
            ellipse = TRUE,
            ind.names = FALSE,
            pch = 20,
            cex = 1.5,
            legend = TRUE,
            legend.title = "Group",
            title = title)
  # Save
  dev.print(
    device = pdf,
    file = paste0(results_folder, title, ".pdf"),
    width = pix(500),
    height = pix(420)
  )
  
}

pca_all(qf, i = "Results", title = "Global PCA 2")

# Volcano -----------------------------------------------------------------

plot_vol <- Volcano_plotter(qf) +
  theme_classic() +
  theme(legend.position="none", 
        plot.margin = margin(t = 10, r = 30,b = 10, l = 30),
        axis.title.y = element_text(vjust=5),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white") ) 

ggsave(
  filename = paste0(results_folder, "EBP Volcano",".pdf"),
  plot     = plot_vol,
  width    = pix(460),
  height   = pix(570),
  device   = cairo_pdf 
)


# Heatmap -----------------------------------------------------------------

Heatmap <- make_heatmap(qf, i = "Results",
                  annot_list = c(Group = "group"), 
                  col_list = list(Group = Group_colors),
                  only_significant = TRUE,
                  title = "none")

Cairo::CairoPDF(paste0(results_folder, "EBP Sign Heatmap.pdf"),
                height = pix(670),
                width = pix(520) )
draw(Heatmap)
dev.off()

# GO OVER REP -------------------------------------------------------------

# Plot
Go_plot <- function(df, color) {
  
  bp <- pairwise_termsim(df)
  bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
  
  mutate(bp2, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore") +
    scale_fill_continuous(low = color, high = color, guide = FALSE)
  
}

plot_go_over  <- Go_plot(Pathway_analysis$GO_results_over, color = "#e06663ff")
# Clean up under:
Pathway_analysis$GO_results_under <- dropGO(Pathway_analysis$GO_results_under, term = c("GO:0002455", "GO:0002460") )
plot_go_under <- Go_plot(Pathway_analysis$GO_results_under, color = "#327ebaff")

# Save
ggsave(
  paste0(results_folder, "EBP GO ORA OVER",".pdf"),
  plot = plot_go_over,
  width    = pix(530),
  height   = pix(470),
  device   = cairo_pdf 
)

ggsave(
  paste0(results_folder, "EBP GO ORA UNDER",".pdf"),
  plot = plot_go_under,
  width    = pix(530),
  height   = pix(470),
  device   = cairo_pdf 
)


# KEGG --------------------------------------------------------------------

KEGG_plot <- ggplot(Pathway_analysis$GSEA_KEGG, showCategory=10, aes(NES, fct_reorder(Description, NES),
                                      fill=qvalue)) +
                geom_col() +
                scale_fill_continuous(
                  name = "Q-value",
                  low  = "#e06663", 
                  high = "#997795ff"
                ) +
                scale_y_discrete(labels = label_wrap(width = 20)) +
                theme_classic(base_size = 15) + 
                xlab("Normalized Enrichment Score") +
                ylab(NULL)

ggsave(
  paste0(results_folder, "EBP KEGG",".pdf"),
  plot = KEGG_plot,
  width    = pix(630),
  height   = pix(607),
  device   = cairo_pdf 
)

# KEGG plot

get_id <- function(strings) {
  idx <- match(strings,
               Pathway_analysis$GSEA_KEGG@result$Description)
  as.list(idx)
}

KEGG_pathway_plot <- gseaplot2(Pathway_analysis$GSEA_KEGG,
                               geneSetID = get_id(c("IL-17 signaling pathway", "MAPK signaling pathway", "Necroptosis")), 
                               color= c("#E495A5", "#86B875", "#7DB0DD"),
                               subplots = 1)
ggsave(
  paste0(results_folder, "EBP KEGG pathways",".pdf"),
  plot = KEGG_pathway_plot,
  width    = pix(630),
  height   = pix(265),
  device   = cairo_pdf 
)

# KEGG Cnetplot
Pathway_analysis$GSEA_KEGG_redeable <- setReadable(Pathway_analysis$GSEA_KEGG, "org.Hs.eg.db", "ENTREZID")
KEGG_cnetplot <- cnetplot(Pathway_analysis$GSEA_KEGG_redeable, foldChange = Pathway_analysis$Ranked_data)  +
  scale_color_gradient2(low = "#327ebaff", mid = "white", high = "#e06663ff")

ggsave(
  paste0(results_folder, "EBP KEGG pathways Cnet",".pdf"),
  plot = KEGG_cnetplot,
  width    = pix(600),
  height   = pix(500),
  device   = cairo_pdf 
)
  
# library("pathview")
# hsa04010 <- pathview(gene.data  = Pathway_analysis$Ranked_data,
#                      pathway.id = "hsa05132",
#                      species    = "hsa",
#                      limit      = list(gene=max(abs(Pathway_analysis$Ranked_data)), cpd=1))


# Profile Plots -----------------------------------------------------------

# Lists:
Proteosome_list = c( "UBB", "PSMA4", "PSMB1", "PSMA1", "PSMB3", "YOD1")
Neutrophil_list = c("ELANE", "MUC5B", "LYZ", "PRTN3")
Cornified_list = c("CSTA", "IVL", "EVPL")

Proteosome_lineplot <- Line_profile_plotter(qf,
                                            select_prots = Proteosome_list,
                                            title = "Proteosome")

Neutrophil_lineplot <- Line_profile_plotter(qf,
                                            select_prots = Neutrophil_list,
                                            title = "Neutrophil mediated cytotoxicity")

Cornified_lineplot <- Line_profile_plotter(qf,
                                           select_prots = Cornified_list,
                                           title = "Cornified envelope")

ggsave(
  paste0(results_folder, "Proteosome lineplot",".pdf"),
  plot = Proteosome_lineplot,
  width    = pix(540),
  height   = pix(320),
  device   = cairo_pdf 
)

ggsave(
  paste0(results_folder, "Neutrophil lineplot",".pdf"),
  plot = Neutrophil_lineplot,
  width    = pix(540),
  height   = pix(320),
  device   = cairo_pdf 
)

ggsave(
  paste0(results_folder, "Cornified lineplot",".pdf"),
  plot = Cornified_lineplot,
  width    = pix(540),
  height   = pix(320),
  device   = cairo_pdf 
)
