
##############################################
#Gabriel Hirdman Proteome Explorer v1.0® 2025#
##############################################

if(Run_DEA == FALSE) {
  print("Loading previously analyzed DEA results")
  load(paste0(saved_data_folder, "DEA_results.R") ) 
}

# Protein IDs -------------------------------------------------------------

Barplot_IDs <- function(qf, i = "", color = "", sort = "") {
  
  df <- countUniqueFeatures(qf, i, colDataName = "Counts") %>% 
    colData(df) %>% 
    as_tibble() %>% 
    ggplot(aes(x = reorder(shortname, .data[[sort]]),
               y = Counts,
               fill = .data[[color]])) + 
     geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank() )
  
}

# Missing values ---------------------------------------------------------

miss_map <- function(qf, i = "") {
  vis_miss(as.data.frame(assay(qf[["Proteins"]])),
           cluster = TRUE)
}

# Normalization -----------------------------------------------------------------

Density_plotter <- function(qf, i = "") {
  p <- qf[[i]] %>%
    QFeatures::longForm() %>%
    as_tibble() %>% 
    ggplot(aes(x = value, group = colname, color = colname)) +
    geom_density() +
    ggtitle(i) +
    theme_bw() +
    theme(legend.position = "none")
  return(p)
}

Norm_boxplotter <- function(qf, i = "") {
  p <- qf[[i]] %>%
    QFeatures::longForm() %>%
    as_tibble() %>% 
    ggplot(aes(x = colname, y=value)) +
    geom_boxplot() +
    ggtitle(i) +
    theme_bw() +
    theme(legend.position = "none")
  return(p)
}

# Intensity Aggregation ---------------------------------------------------

Aggregation_plotter <- function(qf, gene) {
  
  # Get the Protein.Group
  row_data <- as.data.frame(rowData(qf[["Results"]]))
  pg <- row_data$Protein.Group[row_data$Gene == gene]
  if(length(pg) == 0) {
    stop("Gene not found in the row data.")
  }
  pg <- pg[1]
  
  # Get data
  feature_data <- subsetByFeature(qf, pg)
  df <- feature_data[,,c("PeptidesNorm", "Proteins")] %>%
    QFeatures::longForm(colvars = c("group")) %>%
    as_tibble() %>% 
    mutate(assay_order = factor(
      assay,
      levels = c("PeptidesNorm", "Proteins"),
      labels = c("Peptides", "Proteins")
    ))
  
  # Plot
  p <- ggplot(df, aes(x = colname, y = value, colour = assay)) +
    geom_point(size = 3, ) +
    geom_line(aes(group = rowname)) +
    facet_grid(~ assay_order + group, scales = "free_x", space = "free_x") +
    labs(x = "Sample", y = "Abundance") +
    ggtitle(paste(gene, "abundance profiles")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "none")
  
  return(p)
}

# Heatmap -----------------------------------------------------------------

make_heatmap <- function(qf,
                         i,
                         annot_list = NULL,
                         col_list = NULL,
                         only_significant = FALSE,
                         n = 10,
                         title = "") {
  
  uit <- title 
  
  if (!(i %in% names(qf))) stop("Assay 'i' does not exist.")
  
  Annot <- NULL
  if (!is.null(annot_list) && length(annot_list)) {
    missing_cols <- setdiff(annot_list, colnames(colData(qf)))
    if(length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))
    
    annot_args <- lapply(annot_list, function(x) colData(qf)[[x]])
    names(annot_args) <- names(annot_list)
    
    if ("group" %in% names(annot_list)) {
      group_values <- unique(as.character(colData(qf)[["group"]]))
      group_values <- group_values[!is.na(group_values) & nzchar(group_values)]
      if (length(group_values)) {
        group_colors <- tryCatch(get("Group_colors", envir = .GlobalEnv), error = function(...) NULL)
        available <- NULL
        if (!is.null(group_colors)) {
          group_colors <- stats::setNames(as.character(group_colors), names(group_colors))
          available <- group_colors[names(group_colors) %in% group_values]
        } else if (!is.null(col_list) && !is.null(col_list$group)) {
          available <- col_list$group
        }
        if (is.null(available)) {
          available <- grDevices::rainbow(length(group_values))
          names(available) <- group_values
        } else {
          available <- as.character(available)
          if (is.null(names(available)) || any(!nzchar(names(available)))) {
            names(available) <- group_values[seq_along(available)]
          }
        }
        missing_levels <- setdiff(group_values, names(available))
        if (length(missing_levels)) {
          fallback <- grDevices::rainbow(length(missing_levels))
          names(fallback) <- missing_levels
          available <- c(available, fallback)
        }
        if (is.null(col_list)) col_list <- list(group = available) else col_list$group <- available
      }
    }
    
    annot_args$col <- col_list
    Annot <- do.call(HeatmapAnnotation, annot_args)
  }
  
  rdf <- as.data.frame(rowData(qf[[i]])) 
  
  et <- uit # Effective Title
  
  if (only_significant && sum(rdf$signif, na.rm = TRUE) >= 5) {
    qf <- filterFeatures(qf, i = i, ~ signif == TRUE )
    et <- "Heatmap of significant DEPs"
  } else if (only_significant) { # and < 5 significant features (or none were significant)
    top_n <- rownames(slice_max(rdf, log10_qval, n = n))
    qf <- subsetByFeature(qf, top_n) # qf is now filtered
    et <- paste("Heatmap of top ", n," proteins by Q-value")
  } else if (uit == "none") { # only_significant is FALSE
    et <- NULL 
  } else { # only_significant is FALSE AND uit was not "none"
    et <- "Heatmap of all identified proteins" # Original behavior
  }
  
  fct <- if (!is.null(uit) && uit == "none") NULL else et
  
  mat <- assay(qf[[i]])
  mat <- t(scale(t(mat)))
  
  # color
  suppressPackageStartupMessages(library(circlize))
  col_fun = colorRamp2(c(-2, 0, 2), c("#327ebaff", "white", "#e06663ff") )
  col_fun(seq(-3,3))
  
  ComplexHeatmap::Heatmap(
    matrix = mat,
    col = col_fun,
    clustering_distance_columns = "spearman",
    clustering_distance_rows = "spearman",
    clustering_method_rows =  "average",
    clustering_method_columns = "average",
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = Annot,
    column_title = fct,
    heatmap_legend_param = list(title = NULL, legend_height = unit(6, "cm"))
  )
}




# Volcano -----------------------------------------------------------------

# All
Volcano_plotter <- function(qf, Title = "") {
  
  top_genes <- rowData(qf[["Results"]]) %>% 
    as.data.frame() %>% 
    dplyr::filter(!grepl("\\|", Gene)) %>%
    arrange(qval) %>%
    slice(1:15) %>%
    pull(Gene)
  
  plot <- rowData(qf[["Results"]]) %>%
    as_tibble() %>% 
    ggplot(aes(logFC,
               y = log10_qval,
               colour = diff_prot,
               text = paste("Gene:", Gene,
                            "<br>logFC:", logFC,
                            "<br>-q-value:", qval))) +
      geom_point(size = 3) +
      geom_text_repel(aes(label = ifelse(Gene %in% top_genes, Gene, NA)), #Only text at top 25 proteins and not unnamed proteins
                      max.overlaps = 100,
                      size = 4,
                      colour="black" ) +
      scale_colour_manual(values = (c(Expression_lvl_color, "FALSE"= "#A5A6A6"))) +
      xlab("log2 foldchange") + ylab("-log10 FDR adjusted p-value") +
      theme_light() +
      theme(legend.position="none", 
            plot.margin = margin(t = 10, r = 30,b = 10, l = 30),
            axis.title.y = element_text(vjust=5),
            panel.background = element_rect(fill = "white", color = "white"),
            plot.background = element_rect(fill = "white", color = "white") ) +
      ggtitle(Title) +
      geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +  #The lines for the volcano
      geom_vline(xintercept = logFC_cutoff, linetype = 2, alpha = 0.5) +
      geom_vline(xintercept = -logFC_cutoff, linetype = 2, alpha = 0.5)
}

# Highlight
Volcano_plotter_highlight <- function(data, gene) {
  
  df <- as_tibble(rowData(data[["Results"]]))
  
  ggplot(df, aes(logFC,
                 y = log10_qval,
                 text = paste("Gene:", Gene,
                              "<br>logFC:", logFC,
                              "<br>q-value:", qval))) +
    
    geom_point(aes(colour = "#A5A6A6"), size = 3) +
    geom_point(data = df %>% filter(Gene == gene),
               aes(logFC, y = log10_qval),
               colour = "firebrick",
               size = 4) +
    geom_text_repel(data = df %>% filter(Gene == gene),
                    aes(label = Gene),
                    size = 4,
                    colour = "black",
                    nudge_y = 0.08,
                    nudge_x = 0.1) +
    scale_colour_identity() +
    xlab("log2 foldchange") +
    ylab("-log10 FDR adjusted p-value") +
    theme_light() +
    theme(legend.position = "none",
          plot.margin = margin(t = 10, r = 30, b = 10, l = 30),
          axis.title.y = element_text(vjust = 5),
          panel.background = element_rect(fill = "white", color = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
    geom_vline(xintercept = logFC_cutoff, linetype = 2, alpha = 0.5) +
    geom_vline(xintercept = -logFC_cutoff, linetype = 2, alpha = 0.5)
}

# Boxplot -----------------------------------------------------------------

Boxplot_plotter <- function(qf, Boxplot_prot = "") {
  
  plots <- list()
  
  for (prot in Boxplot_prot) {
    df <- QFeatures::longForm(qf[, , "Results"],
                                colvars = c("group"),
                                rowvars = c("Protein.Names", "Gene", "First.Protein.Description")) %>% 
      as_tibble() %>% 
      dplyr::filter(Gene == prot)
    
    if (nrow(df) > 0) {
      
      # Setup
      Title <- unique(df$Gene)
      Subtitle <- unique(df$First.Protein.Description)
      
      # Plot
      p <- ggplot(df, aes(x = group, y = value, fill = group)) +
        geom_boxplot() +
        scale_fill_manual(values = Group_colors) +
        theme_bw() +
        labs(y = expression("LFQ intensity Log"[2]),
             title = Title,
             subtitle = Subtitle) +
        theme(legend.position = "none",
              axis.title.x = element_blank(),
              axis.text.x  = element_text(size = 15, angle = 10, vjust = 0.5, hjust = 0.5))
      
      # Save
      plots[[prot]] <- p
      
    } else {
      warning(paste("No data found for protein:", prot))
    }
  }
  return(plots)
}

Boxplot_plotter(qf, Boxplot_prot = "ELANE")

BatchBoxplot_plotter <- function(qf, Boxplot_prot = "", subtitle = "") {
  
  plots <- list()
  
  for (prot in Boxplot_prot) {
    df <- QFeatures::longForm(qf[, , "Proteins"],
                                colvars = c("group", "Extraction_batch"),
                                rowvars = c("Protein.Names", "Gene", "First.Protein.Description")) %>% 
      as_tibble() %>% 
      filter(Gene == prot)
    
    if (nrow(df) > 0) {
      
      # Z-scale
      df <- df %>%
        mutate(zvalue = as.numeric(scale(value)))
      
      # Setup
      Title <- unique(df$Gene)
      Subtitle <- unique(df$First.Protein.Description)
      
      p <- ggplot(df, aes(x = Extraction_batch, y = zvalue, fill = group)) +
        geom_boxplot() +
        theme_bw() +
        labs(y = expression("LFQ intensity Log"[2]),
             title = paste(Title, subtitle),
             subtitle = Subtitle) +
        scale_y_continuous(limits = c(-3.5, 3.5)) +
        scale_fill_manual(values = Group_colors) +
        theme(legend.position = "none",
              axis.title.x = element_blank(),
              axis.text.x  = element_text(size = 15, angle = 10, vjust = 0.5, hjust = 0.5))
      
      # Save the plot
      plots[[prot]] <- p
      
    } else {
      warning(paste("No data found for protein:", prot))
    }
  }
  return(plots)
}

# Correlation -------------------------------------------------------------

Corr_plotter <- function(qf, i = "") {
  df <- qf[["Proteins"]] %>% 
    assay() %>% 
    as.data.frame() %>% 
    cor(method = "pearson",
        use = "pairwise.complete.obs") %>% 
    corrplot(method = "color",
             col = NULL,
             type = "upper",
             # addCoef.col = "white",
             diag = FALSE,
             tl.col = "black",
             tl.srt = 45,
             outline = TRUE)
  return(df)
}

# PCA plot ----------------------------------------------------------------

PCA_plotter <- function(qf, i = "", color = "", label = "all") {
  
  df <- qf[["Proteins"]] %>% 
        filterNA() %>% 
        assay() %>% 
        t()
        
  res <- PCA(df, graph = FALSE)
  
  # add background:
  p <- factoextra::fviz_pca_ind(
    res,
    pointshape = 19,
     pointsize = 3,
     habillage = as.factor(colData(qf)[[color]]),
     mean.point = FALSE,
     label = label,
     repel = TRUE) +
    theme_bw()  
  
    return(p)
}


# UMAP --------------------------------------------------------------------

UMAP_plotter <- function(qf, 
                         i = "Proteins", 
                         color = "group", 
                         umap_config = umap.defaults, 
                         dataset_name = "UMAP Plot") {
  
  # Extract and prepare the assay data from the QFeatures object
  umap_data <- qf[[i]] %>% 
    assay() %>% 
    t()
  
  # Check for missing values and perform simple mean imputation if needed
  if (any(is.na(umap_data))) {
    mean_impute <- function(x) { 
      mean_val <- mean(x, na.rm = TRUE)
      x[is.na(x)] <- mean_val 
      return(x)
    }
    umap_data <- as.data.frame(lapply(as.data.frame(umap_data), mean_impute))
    umap_data <- umap_data[, colSums(is.na(umap_data)) == 0, drop = FALSE]
  }
  
  # Scale the data
  scaled_data <- scale(umap_data)
  
  umap_res <- umap::umap(scale(umap_data), config = umap_config)
  
  # Prepare data for plotting
  umap_plot_df <- data.frame(
    UMAP1 = umap_res$layout[, 1],
    UMAP2 = umap_res$layout[, 2],
    Group  = factor(colData(qf)[[color]])
  )
  
  # Generate the UMAP plot using ggplot2
  p <- ggplot(umap_plot_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
    geom_point(alpha = 0.8) +
    labs(title = dataset_name, x = "UMAP 1", y = "UMAP 2", color = color) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(p)
}


# Profile Line Plots ------------------------------------------------------

Line_profile_plotter <- function(qf,
                                 i = "Results",
                                 select_prots = Proteosome_list,
                                 title = "",
                                 colors = c("Normal" = "grey",
                                            "Cancer" = "#006D77") ) {
                                   
  # Data setup                                 
  df <- QFeatures::longForm(qf[, , i],
                             colvars = c("group"),
                             rowvars = c("Gene")) %>%
   as_tibble() %>%
   dplyr::filter(Gene %in% select_prots ) %>% 
   select(-assay, -primary, -rowname) %>%
   pivot_wider(names_from = Gene, values_from = value)
  
  cols <- match(select_prots, names(df))
  cols <- cols[!is.na(cols)]

  # Plot
  df_plot <- df %>%
    ggparcoord(
      columns = which(names(df) %in% select_prots),
      scale = "center",
      groupColumn = "group",
      order = cols,
      showPoints = FALSE,
      title = title,
      alphaLines = 0.3
    ) +
    scale_color_manual(values = colors) +
    theme_ipsum()+
    theme(
      legend.position = "Default",
      plot.title = element_text(face = "plain", size=12)
    ) +
    xlab("") + 
    ylab("")
}


