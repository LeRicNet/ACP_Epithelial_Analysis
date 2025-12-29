# ==============================================================================
# R/figures/figure_generators.R
# ==============================================================================
# Functions to generate publication-ready figures
# Each function generates one figure panel or composite figure
# ==============================================================================

#' Generate All Main Figures
#'
#' Master function to regenerate all main manuscript figures
#'
#' @param seurat_obj Seurat object with classifications
#' @param trajectory_results Trajectory analysis results
#' @param config Configuration list
#' @return List of figure paths
#' @export
generate_all_figures <- function(seurat_obj, trajectory_results, config) {
  
  message("\n========================================")
  message("Generating All Figures")
  message("========================================\n")
  
  figure_paths <- list()
  
  # Main figures
  figure_paths$fig1 <- generate_figure_1(seurat_obj, config)
  figure_paths$fig2 <- generate_figure_2(seurat_obj, config)
  figure_paths$fig3 <- generate_figure_3(trajectory_results, config)
  figure_paths$fig4 <- generate_figure_4(seurat_obj, trajectory_results, config)
  
  # Supplementary figures
  figure_paths$supp_classification <- generate_supp_classification(seurat_obj, config)
  figure_paths$supp_trajectory <- generate_supp_trajectory(trajectory_results, config)
  
  message("\nAll figures generated!")
  return(figure_paths)
}

# ==============================================================================
# FIGURE 1: Epithelial Subtype Classification
# ==============================================================================

#' Generate Figure 1: Classification Overview
#'
#' Panel A: UMAP colored by subtype
#' Panel B: Module score heatmap
#' Panel C: Subtype proportions
#' Panel D: Key marker expression
#'
#' @param seurat_obj Seurat object
#' @param config Configuration list
#' @return Path to saved figure
#' @export
generate_figure_1 <- function(seurat_obj, config) {
  
  message("Generating Figure 1: Classification Overview...")
  
  library(ggplot2)
  library(patchwork)
  
  colors <- get_colors(config, "epithelial_subtypes")
  
  # Panel A: UMAP
  p_umap <- Seurat::DimPlot(seurat_obj, 
                            group.by = "module_score_subtype",
                            cols = colors,
                            pt.size = 0.5) +
    labs(title = "Epithelial Subtypes", color = "Subtype") +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Panel B: Module score heatmap
  score_cols <- grep("^score_", colnames(seurat_obj@meta.data), value = TRUE)
  score_cols <- score_cols[!grepl("Proliferation", score_cols)]
  
  if (length(score_cols) > 0) {
    score_means <- seurat_obj@meta.data %>%
      dplyr::group_by(module_score_subtype) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(score_cols), mean, na.rm = TRUE)) %>%
      tidyr::pivot_longer(-module_score_subtype, names_to = "Signature", values_to = "Score")
    
    score_means$Signature <- gsub("^score_", "", score_means$Signature)
    
    p_heatmap <- ggplot(score_means, 
                        aes(x = Signature, y = module_score_subtype, fill = Score)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%.2f", Score)), size = 2.5) +
      scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Mean Module Scores", x = "", y = "")
  } else {
    p_heatmap <- ggplot() + theme_void() + 
      annotate("text", x = 0.5, y = 0.5, label = "No scores available")
  }
  
  # Panel C: Subtype proportions
  subtype_counts <- as.data.frame(table(seurat_obj$module_score_subtype))
  colnames(subtype_counts) <- c("Subtype", "Count")
  subtype_counts$Percentage <- subtype_counts$Count / sum(subtype_counts$Count) * 100
  
  p_bar <- ggplot(subtype_counts, aes(x = reorder(Subtype, -Count), y = Count, fill = Subtype)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)), vjust = -0.3, size = 3) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Subtype Distribution", x = "", y = "Number of Cells")
  
  # Panel D: Key markers
  markers <- config$markers$keratins[1:4]  # KRT5, KRT8, KRT14, KRT17
  markers <- markers[markers %in% rownames(seurat_obj)]
  
  if (length(markers) > 0) {
    p_markers <- Seurat::VlnPlot(seurat_obj, 
                                  features = markers,
                                  group.by = "module_score_subtype",
                                  cols = colors,
                                  pt.size = 0) +
      plot_layout(ncol = 2)
  } else {
    p_markers <- ggplot() + theme_void()
  }
  
  # Combine panels
  fig1 <- (p_umap | p_heatmap) / (p_bar | p_markers) +
    plot_annotation(title = "Figure 1: Epithelial Subtype Classification",
                    tag_levels = "A")
  
  # Save
  path <- save_figure(fig1, "figure_1_classification", config, 
                      subdir = "main", width = 14, height = 12)
  
  return(path)
}

# ==============================================================================
# FIGURE 2: Transit-Amplifying Population Characterization
# ==============================================================================

#' Generate Figure 2: Transit-Amplifying Cells
#'
#' Panel A: Proliferation score distribution
#' Panel B: Cell cycle phase distribution
#' Panel C: UMI counts by subtype
#' Panel D: Key proliferation markers
#'
#' @param seurat_obj Seurat object
#' @param config Configuration list
#' @return Path to saved figure
#' @export
generate_figure_2 <- function(seurat_obj, config) {
  
  message("Generating Figure 2: Transit-Amplifying Characterization...")
  
  library(ggplot2)
  library(patchwork)
  
  colors <- get_colors(config, "epithelial_subtypes")
  
  # Panel A: Proliferation score
  if ("score_Proliferation" %in% colnames(seurat_obj@meta.data)) {
    p_prolif <- ggplot(seurat_obj@meta.data, 
                       aes(x = module_score_subtype, y = score_Proliferation, 
                           fill = module_score_subtype)) +
      geom_violin(scale = "width") +
      geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(title = "Proliferation Score", x = "", y = "Score")
  } else {
    p_prolif <- ggplot() + theme_void() + 
      annotate("text", x = 0.5, y = 0.5, label = "Proliferation score not available")
  }
  
  # Panel B: Cell cycle phases
  if ("Phase" %in% colnames(seurat_obj@meta.data)) {
    phase_data <- seurat_obj@meta.data %>%
      dplyr::count(module_score_subtype, Phase) %>%
      dplyr::group_by(module_score_subtype) %>%
      dplyr::mutate(Proportion = n / sum(n))
    
    cycle_colors <- get_colors(config, "cell_cycle")
    
    p_cycle <- ggplot(phase_data, 
                      aes(x = module_score_subtype, y = Proportion, fill = Phase)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = cycle_colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Cell Cycle Distribution", x = "", y = "Proportion")
  } else {
    p_cycle <- ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "Cell cycle not scored")
  }
  
  # Panel C: UMI counts
  umi_col <- intersect(c("nCount_RNA", "nCount_Spatial"), colnames(seurat_obj@meta.data))[1]
  
  if (!is.na(umi_col)) {
    p_umi <- ggplot(seurat_obj@meta.data, 
                    aes(x = module_score_subtype, y = .data[[umi_col]], 
                        fill = module_score_subtype)) +
      geom_violin(scale = "width") +
      geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
      scale_fill_manual(values = colors) +
      scale_y_log10() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(title = "Transcriptional Activity", x = "", y = "UMI Counts (log10)")
  } else {
    p_umi <- ggplot() + theme_void()
  }
  
  # Panel D: Proliferation markers
  prolif_markers <- c("MKI67", "TOP2A", "PCNA", "STMN1")
  prolif_markers <- prolif_markers[prolif_markers %in% rownames(seurat_obj)]
  
  if (length(prolif_markers) > 0) {
    p_markers <- Seurat::VlnPlot(seurat_obj,
                                  features = prolif_markers[1:min(4, length(prolif_markers))],
                                  group.by = "module_score_subtype",
                                  cols = colors,
                                  pt.size = 0) +
      plot_layout(ncol = 2)
  } else {
    p_markers <- ggplot() + theme_void()
  }
  
  # Combine
  fig2 <- (p_prolif | p_cycle) / (p_umi | p_markers) +
    plot_annotation(title = "Figure 2: Transit-Amplifying Population",
                    tag_levels = "A")
  
  path <- save_figure(fig2, "figure_2_transit_amplifying", config,
                      subdir = "main", width = 14, height = 12)
  
  return(path)
}

# ==============================================================================
# FIGURE 3: Trajectory Analysis
# ==============================================================================

#' Generate Figure 3: Trajectory Analysis
#'
#' Panel A: Trajectory plot colored by subtype
#' Panel B: Trajectory plot colored by pseudotime
#' Panel C: Pseudotime violin by subtype
#' Panel D: Mean pseudotime bar plot
#'
#' @param trajectory_results Trajectory analysis results
#' @param config Configuration list
#' @return Path to saved figure
#' @export
generate_figure_3 <- function(trajectory_results, config) {
  
  message("Generating Figure 3: Trajectory Analysis...")
  
  library(ggplot2)
  library(patchwork)
  library(monocle3)
  
  colors <- get_colors(config, "epithelial_subtypes")
  
  # Handle per-sample results or single cds
  if (inherits(trajectory_results, "cell_data_set")) {
    cds <- trajectory_results
  } else if (is.list(trajectory_results) && !is.null(trajectory_results$pooled)) {
    cds <- trajectory_results$pooled$cds
  } else {
    # Try to find first valid cds
    for (res in trajectory_results) {
      if (!res$skipped && !is.null(res$cds)) {
        cds <- res$cds
        break
      }
    }
  }
  
  if (!exists("cds")) {
    warning("No valid trajectory results found")
    return(NULL)
  }
  
  # Panel A: Trajectory by subtype
  p_traj_subtype <- monocle3::plot_cells(cds,
                                          color_cells_by = "subtype",
                                          show_trajectory_graph = TRUE,
                                          label_groups_by_cluster = FALSE,
                                          label_cell_groups = FALSE,
                                          cell_size = 0.5) +
    scale_color_manual(values = colors) +
    labs(title = "Trajectory by Subtype") +
    theme(legend.position = "right")
  
  # Panel B: Trajectory by pseudotime
  p_traj_pt <- monocle3::plot_cells(cds,
                                     color_cells_by = "pseudotime",
                                     show_trajectory_graph = TRUE,
                                     label_groups_by_cluster = FALSE,
                                     label_cell_groups = FALSE,
                                     cell_size = 0.5) +
    scale_color_viridis_c() +
    labs(title = "Trajectory by Pseudotime")
  
  # Panel C: Pseudotime violin
  pseudotime_df <- data.frame(
    subtype = SummarizedExperiment::colData(cds)$subtype,
    pseudotime = monocle3::pseudotime(cds)
  )
  
  p_violin <- ggplot(pseudotime_df, aes(x = reorder(subtype, pseudotime), 
                                         y = pseudotime, fill = subtype)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Pseudotime Distribution", x = "", y = "Pseudotime")
  
  # Panel D: Mean pseudotime
  mean_pt <- pseudotime_df %>%
    dplyr::group_by(subtype) %>%
    dplyr::summarise(
      mean_pt = mean(pseudotime, na.rm = TRUE),
      se_pt = sd(pseudotime, na.rm = TRUE) / sqrt(dplyr::n())
    )
  
  p_bar <- ggplot(mean_pt, aes(x = reorder(subtype, mean_pt), y = mean_pt, fill = subtype)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mean_pt - se_pt, ymax = mean_pt + se_pt), width = 0.2) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Mean Pseudotime", x = "", y = "Pseudotime")
  
  # Combine
  fig3 <- (p_traj_subtype | p_traj_pt) / (p_violin | p_bar) +
    plot_annotation(title = "Figure 3: Pseudotemporal Trajectory Analysis",
                    tag_levels = "A")
  
  path <- save_figure(fig3, "figure_3_trajectory", config,
                      subdir = "main", width = 14, height = 12)
  
  return(path)
}

# ==============================================================================
# FIGURE 4: Gene Expression Dynamics
# ==============================================================================

#' Generate Figure 4: Gene Expression Along Trajectory
#'
#' Shows expression of key markers along pseudotime
#'
#' @param seurat_obj Seurat object
#' @param trajectory_results Trajectory results
#' @param config Configuration list
#' @return Path to saved figure
#' @export
generate_figure_4 <- function(seurat_obj, trajectory_results, config) {
  
  message("Generating Figure 4: Gene Expression Dynamics...")
  
  library(ggplot2)
  library(patchwork)
  
  colors <- get_colors(config, "epithelial_subtypes")
  
  # Get pseudotime from trajectory
  if (inherits(trajectory_results, "cell_data_set")) {
    cds <- trajectory_results
  } else if (is.list(trajectory_results)) {
    # Find valid cds
    for (res in trajectory_results) {
      if (!is.null(res$cds)) {
        cds <- res$cds
        break
      }
    }
  }
  
  pseudotime <- monocle3::pseudotime(cds)
  subtypes <- SummarizedExperiment::colData(cds)$subtype
  
  # Match cells between cds and seurat
  common_cells <- intersect(names(pseudotime), colnames(seurat_obj))
  
  if (length(common_cells) == 0) {
    warning("No common cells between trajectory and Seurat object")
    return(NULL)
  }
  
  # Genes to plot
  gene_groups <- list(
    "Basal Keratins" = c("KRT5", "KRT14"),
    "Intermediate Keratins" = c("KRT8", "KRT19"),
    "Cell Cycle" = c("MKI67", "CDKN1A"),
    "Whorl Markers" = c("SHH", "CTNNB1")
  )
  
  plots <- list()
  
  for (group_name in names(gene_groups)) {
    genes <- gene_groups[[group_name]]
    genes <- genes[genes %in% rownames(seurat_obj)]
    
    if (length(genes) == 0) next
    
    # Build expression data
    expr_data <- data.frame(
      cell = common_cells,
      pseudotime = pseudotime[common_cells],
      subtype = subtypes[common_cells]
    )
    
    for (gene in genes) {
      expr_data[[gene]] <- as.numeric(seurat_obj@assays$RNA@data[gene, common_cells])
    }
    
    # Reshape for plotting
    expr_long <- tidyr::pivot_longer(expr_data, cols = dplyr::all_of(genes),
                                      names_to = "Gene", values_to = "Expression")
    
    p <- ggplot(expr_long, aes(x = pseudotime, y = Expression, color = subtype)) +
      geom_point(size = 0.3, alpha = 0.3) +
      geom_smooth(aes(group = 1), method = "loess", color = "black", se = FALSE) +
      facet_wrap(~Gene, scales = "free_y", ncol = 2) +
      scale_color_manual(values = colors) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(title = group_name, x = "Pseudotime", y = "Expression")
    
    plots[[group_name]] <- p
  }
  
  # Combine
  if (length(plots) >= 4) {
    fig4 <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]]) +
      plot_annotation(title = "Figure 4: Gene Expression Dynamics Along Trajectory",
                      tag_levels = "A")
  } else {
    fig4 <- wrap_plots(plots) +
      plot_annotation(title = "Figure 4: Gene Expression Dynamics")
  }
  
  path <- save_figure(fig4, "figure_4_expression_dynamics", config,
                      subdir = "main", width = 14, height = 12)
  
  return(path)
}

# ==============================================================================
# SUPPLEMENTARY FIGURES
# ==============================================================================

#' Generate Supplementary Classification Figure
#'
#' @param seurat_obj Seurat object
#' @param config Configuration list
#' @return Path to saved figure
#' @export
generate_supp_classification <- function(seurat_obj, config) {
  
  message("Generating Supplementary: Classification Diagnostics...")
  
  library(ggplot2)
  library(patchwork)
  
  colors <- get_colors(config, "epithelial_subtypes")
  
  # Get all classification diagnostic plots if available
  class_plots <- attr(seurat_obj, "classification_plots")
  
  if (!is.null(class_plots)) {
    supp <- wrap_plots(class_plots[1:min(4, length(class_plots))]) +
      plot_annotation(title = "Supplementary: Classification Diagnostics",
                      tag_levels = "A")
  } else {
    # Generate basic comparison plots
    
    # Module score distributions
    score_cols <- grep("^score_", colnames(seurat_obj@meta.data), value = TRUE)
    
    if (length(score_cols) > 0) {
      score_data <- seurat_obj@meta.data[, c(score_cols, "module_score_subtype")]
      score_long <- tidyr::pivot_longer(score_data, cols = -module_score_subtype,
                                         names_to = "Signature", values_to = "Score")
      score_long$Signature <- gsub("^score_", "", score_long$Signature)
      
      p1 <- ggplot(score_long, aes(x = Score, fill = Signature)) +
        geom_density(alpha = 0.5) +
        facet_wrap(~Signature) +
        theme_minimal() +
        labs(title = "Module Score Distributions")
      
      p2 <- ggplot(score_long, aes(x = module_score_subtype, y = Score, 
                                    fill = module_score_subtype)) +
        geom_boxplot() +
        facet_wrap(~Signature, scales = "free_y") +
        scale_fill_manual(values = colors) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none") +
        labs(title = "Scores by Subtype")
      
      supp <- p1 / p2 +
        plot_annotation(title = "Supplementary: Classification Diagnostics")
    } else {
      supp <- ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No score data available")
    }
  }
  
  path <- save_figure(supp, "supp_classification", config,
                      subdir = "supplementary", width = 14, height = 14)
  
  return(path)
}

#' Generate Supplementary Trajectory Figure
#'
#' Per-sample trajectory analysis results
#'
#' @param trajectory_results Per-sample trajectory results
#' @param config Configuration list
#' @return Path to saved figure
#' @export
generate_supp_trajectory <- function(trajectory_results, config) {
  
  message("Generating Supplementary: Per-Sample Trajectory...")
  
  library(ggplot2)
  library(patchwork)
  
  colors <- get_colors(config, "epithelial_subtypes")
  
  # Check if we have per-sample results
  if (!is.list(trajectory_results) || inherits(trajectory_results, "cell_data_set")) {
    message("  No per-sample results available")
    return(NULL)
  }
  
  # Get concordance info
  concordance <- attr(trajectory_results, "concordance")
  
  plots <- list()
  
  # Plot 1: Rank matrix heatmap
  if (!is.null(concordance$rank_matrix)) {
    rank_long <- reshape2::melt(concordance$rank_matrix)
    colnames(rank_long) <- c("Subtype", "Sample", "Rank")
    
    plots$rank_heatmap <- ggplot(rank_long, aes(x = Sample, y = Subtype, fill = Rank)) +
      geom_tile() +
      geom_text(aes(label = round(Rank, 1)), size = 3) +
      scale_fill_viridis_c() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Subtype Rank by Sample")
  }
  
  # Plot 2: Branch point counts
  branch_data <- data.frame(
    Sample = names(trajectory_results),
    BranchPoints = sapply(trajectory_results, function(x) {
      if (x$skipped) return(NA)
      x$n_branch_points
    })
  )
  
  plots$branch_points <- ggplot(branch_data, aes(x = Sample, y = BranchPoints, fill = Sample)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Branch Points per Sample", y = "Number of Branch Points")
  
  # Plot 3: Per-sample pseudotime distributions
  pt_data <- do.call(rbind, lapply(names(trajectory_results), function(s) {
    res <- trajectory_results[[s]]
    if (res$skipped) return(NULL)
    
    data.frame(
      Sample = s,
      Subtype = names(res$pseudotime_by_subtype),
      MeanPseudotime = as.numeric(res$pseudotime_by_subtype)
    )
  }))
  
  if (!is.null(pt_data)) {
    plots$pt_by_sample <- ggplot(pt_data, aes(x = Subtype, y = MeanPseudotime, 
                                               fill = Subtype)) +
      geom_bar(stat = "identity") +
      facet_wrap(~Sample) +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(title = "Mean Pseudotime by Sample")
  }
  
  # Combine
  supp <- wrap_plots(plots) +
    plot_annotation(title = sprintf("Supplementary: Per-Sample Trajectory Analysis\nMean concordance ρ = %.3f",
                                    concordance$mean_correlation),
                    tag_levels = "A")
  
  path <- save_figure(supp, "supp_trajectory_per_sample", config,
                      subdir = "supplementary", width = 14, height = 10)
  
  return(path)
}
