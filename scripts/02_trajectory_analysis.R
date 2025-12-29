#!/usr/bin/env Rscript
# ==============================================================================
# scripts/02_trajectory_analysis.R
# ==============================================================================
# Multi-Method Trajectory Analysis with Consensus and Per-Sample Validation
#
# This script performs trajectory inference using multiple methods:
#   1. Monocle3 - graph-based trajectory
#   2. Slingshot - principal curves through clusters
#   3. Diffusion pseudotime (optional)
#
# Input:
#   - Annotated Seurat object from 01_cell_type_annotation.R
#
# Output:
#   - Trajectory results with consensus pseudotime
#   - Per-sample validation
#   - Diagnostic plots
#
# Usage:
#   Rscript scripts/02_trajectory_analysis.R
#   Rscript scripts/02_trajectory_analysis.R --config path/to/config.yaml
#
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL

if (length(args) > 0) {
  if (args[1] == "--config" && length(args) > 1) {
    config_path <- args[2]
  } else if (!startsWith(args[1], "--")) {
    config_path <- args[1]
  }
}

# Find project root
find_project_root <- function() {
  if (file.exists("config/config.yaml")) return(getwd())
  if (file.exists("../config/config.yaml")) return(normalizePath(".."))
  if (requireNamespace("here", quietly = TRUE)) {
    root <- here::here()
    if (file.exists(file.path(root, "config/config.yaml"))) return(root)
  }
  current <- getwd()
  for (i in 1:5) {
    if (file.exists(file.path(current, "config/config.yaml"))) return(current)
    parent <- dirname(current)
    if (parent == current) break
    current <- parent
  }
  return(NULL)
}

project_root <- find_project_root()
if (is.null(project_root)) {
  stop("Could not find project root directory")
}
setwd(project_root)

message("\n")
message("================================================================")
message("  Step 2: Multi-Method Trajectory Analysis")
message("================================================================")
message(paste("  Started:", Sys.time()))
message(paste("  Working directory:", getwd()))
message("\n")

# Load packages
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Source utility functions
source("R/utils/config.R")
source("R/trajectory/multi_method_trajectory.R")

# Load configuration
config <- load_config(config_path)
print_config_summary(config)

# Set seed
set.seed(config$reproducibility$seed)

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("\n========================================")
message("Loading Annotated Data")
message("========================================\n")

# Load from Step 1
input_path <- file.path(get_path(config, config$paths$objects_dir), "01_seurat_annotated.rds")

if (!file.exists(input_path)) {
  stop(paste("Annotated object not found:", input_path,
             "\nPlease run 01_cell_type_annotation.R first"))
}

message(paste("Loading:", input_path))
seurat_obj <- readRDS(input_path)

message(sprintf("  Cells: %d", ncol(seurat_obj)))
message(sprintf("  Genes: %d", nrow(seurat_obj)))

# Check for classification
if (!"module_score_subtype" %in% colnames(seurat_obj@meta.data)) {
  stop("module_score_subtype not found. Run classification first.")
}

message("\nSubtype distribution:")
print(table(seurat_obj$module_score_subtype))

# ==============================================================================
# FILTER TO EPITHELIAL CELLS
# ==============================================================================

message("\n========================================")
message("Filtering to Epithelial Cells")
message("========================================\n")

# Remove unassigned and non-epithelial cells for trajectory
exclude_subtypes <- c("Unassigned", "Non-epithelial", "Unknown")
epithelial_mask <- !seurat_obj$module_score_subtype %in% exclude_subtypes &
  !is.na(seurat_obj$module_score_subtype)

n_epithelial <- sum(epithelial_mask)
message(sprintf("Epithelial cells for trajectory: %d (%.1f%%)",
                n_epithelial, n_epithelial / ncol(seurat_obj) * 100))

if (n_epithelial < 100) {
  stop("Insufficient epithelial cells for trajectory analysis (< 100)")
}

seurat_epi <- subset(seurat_obj, cells = colnames(seurat_obj)[epithelial_mask])

message("\nFiltered subtype distribution:")
print(table(seurat_epi$module_score_subtype))

# ==============================================================================
# CHECK/CREATE DIMENSIONALITY REDUCTION
# ==============================================================================

message("\n========================================")
message("Checking Dimensionality Reduction")
message("========================================\n")

# Check for existing reductions
available_reductions <- names(seurat_epi@reductions)
message(sprintf("Available reductions: %s",
                ifelse(length(available_reductions) > 0,
                       paste(available_reductions, collapse = ", "),
                       "none")))

# Run PCA if not present
if (!"pca" %in% tolower(available_reductions)) {
  message("\nRunning PCA...")
  seurat_epi <- NormalizeData(seurat_epi, verbose = FALSE)
  seurat_epi <- FindVariableFeatures(seurat_epi, verbose = FALSE)
  seurat_epi <- ScaleData(seurat_epi, verbose = FALSE)
  seurat_epi <- RunPCA(seurat_epi, npcs = 50, verbose = FALSE)
  message("  PCA complete")
}

# Run UMAP if not present
if (!"umap" %in% tolower(available_reductions)) {
  message("\nRunning UMAP...")
  seurat_epi <- RunUMAP(seurat_epi, dims = 1:30, verbose = FALSE)
  message("  UMAP complete")
}

# ==============================================================================
# RUN MULTI-METHOD TRAJECTORY ANALYSIS
# ==============================================================================

message("\n========================================")
message("Running Multi-Method Trajectory Analysis")
message("========================================\n")

# Get methods from config
trajectory_methods <- c("monocle3", "slingshot")
if (!is.null(config$trajectory$methods)) {
  trajectory_methods <- config$trajectory$methods
}

message(sprintf("Methods: %s", paste(trajectory_methods, collapse = ", ")))

# Run pooled analysis
pooled_results <- run_multi_method_trajectory(
  seurat_epi,
  config = config,
  methods = trajectory_methods,
  subtype_column = "module_score_subtype"
)

# Update seurat object with pseudotime
seurat_epi <- pooled_results$seurat_obj

message("\nPooled analysis complete!")
message(sprintf("  Methods used: %s", paste(pooled_results$methods_used, collapse = ", ")))

# ==============================================================================
# VALIDATE AGAINST EXPECTED ORDER
# ==============================================================================

message("\n========================================")
message("Validating Trajectory Order")
message("========================================\n")

expected_order <- c("Basal-like", "Transit-Amplifying", "Intermediate", "Specialized")
if (!is.null(config$trajectory$expected_order)) {
  expected_order <- config$trajectory$expected_order
}

validation_pooled <- validate_trajectory_order(pooled_results, expected_order)
print(validation_pooled)

# ==============================================================================
# PER-SAMPLE ANALYSIS
# ==============================================================================

per_sample_enabled <- TRUE
if (!is.null(config$trajectory$per_sample$enabled)) {
  per_sample_enabled <- config$trajectory$per_sample$enabled
}

per_sample_results <- NULL

if (per_sample_enabled && "sample" %in% colnames(seurat_epi@meta.data)) {

  message("\n========================================")
  message("Running Per-Sample Analysis")
  message("========================================\n")

  per_sample_results <- run_per_sample_multi_trajectory(
    seurat_epi,
    config = config,
    methods = trajectory_methods,
    sample_column = "sample",
    subtype_column = "module_score_subtype"
  )

  # Validate per-sample
  validation_per_sample <- validate_trajectory_order(per_sample_results, expected_order)

  message("\nPer-sample validation:")
  print(validation_per_sample)

} else {
  message("\nSkipping per-sample analysis (disabled or no sample column)")
}

# ==============================================================================
# GENERATE DIAGNOSTIC PLOTS
# ==============================================================================

message("\n========================================")
message("Generating Diagnostic Plots")
message("========================================\n")

fig_dir <- get_path(config, config$paths$figures_supp_dir)
ensure_dir(fig_dir)

# Get colors
colors <- get_colors(config, "epithelial_subtypes")

# Plot 1: Multi-method comparison
message("Creating method comparison plot...")
comparison_plot <- tryCatch({
  plot_trajectory_comparison(pooled_results, config)
}, error = function(e) {
  message(sprintf("  Could not create comparison plot: %s", e$message))
  NULL
})

if (!is.null(comparison_plot)) {
  ggsave(file.path(fig_dir, "02_trajectory_method_comparison.pdf"),
         comparison_plot, width = 14, height = 10)
  ggsave(file.path(fig_dir, "02_trajectory_method_comparison.png"),
         comparison_plot, width = 14, height = 10, dpi = 300)
  message("  Saved method comparison plot")
}

# Plot 2: Pseudotime on UMAP
message("Creating UMAP pseudotime plots...")

umap_plots <- list()

# Consensus pseudotime
umap_plots$consensus <- FeaturePlot(seurat_epi, features = "pseudotime_consensus",
                                    pt.size = 0.5) +
  scale_color_viridis_c() +
  labs(title = "Consensus Pseudotime", color = "Pseudotime")

# By subtype
umap_plots$subtype <- DimPlot(seurat_epi, group.by = "module_score_subtype",
                              cols = colors, pt.size = 0.5) +
  labs(title = "Epithelial Subtypes")

# Confidence
if ("pseudotime_confidence" %in% colnames(seurat_epi@meta.data)) {
  umap_plots$confidence <- FeaturePlot(seurat_epi, features = "pseudotime_confidence",
                                       pt.size = 0.5) +
    scale_color_viridis_c(option = "plasma") +
    labs(title = "Pseudotime Confidence", color = "Confidence")
}

# Per-method pseudotime
for (method in pooled_results$methods_used) {
  pt_col <- paste0("pseudotime_", method)
  if (pt_col %in% colnames(seurat_epi@meta.data)) {
    umap_plots[[method]] <- FeaturePlot(seurat_epi, features = pt_col, pt.size = 0.5) +
      scale_color_viridis_c() +
      labs(title = paste(method, "Pseudotime"), color = "Pseudotime")
  }
}

umap_combined <- wrap_plots(umap_plots, ncol = 2) +
  plot_annotation(title = "Trajectory Analysis: Pseudotime on UMAP")

ggsave(file.path(fig_dir, "02_trajectory_umap.pdf"), umap_combined, width = 12, height = 12)
ggsave(file.path(fig_dir, "02_trajectory_umap.png"), umap_combined, width = 12, height = 12, dpi = 300)
message("  Saved UMAP plots")

# Plot 3: Pseudotime by subtype
message("Creating pseudotime distribution plots...")

pt_violin <- ggplot(seurat_epi@meta.data,
                    aes(x = reorder(module_score_subtype, pseudotime_consensus),
                        y = pseudotime_consensus, fill = module_score_subtype)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Consensus Pseudotime Distribution",
       x = "Subtype (ordered by mean pseudotime)",
       y = "Pseudotime")

ggsave(file.path(fig_dir, "02_trajectory_pseudotime_violin.pdf"), pt_violin, width = 10, height = 6)
ggsave(file.path(fig_dir, "02_trajectory_pseudotime_violin.png"), pt_violin, width = 10, height = 6, dpi = 300)
message("  Saved pseudotime violin plot")

# Plot 4: Per-sample concordance (if available)
if (!is.null(per_sample_results)) {
  message("Creating per-sample concordance plot...")

  concordance <- attr(per_sample_results, "concordance")

  if (!is.null(concordance$rank_matrix)) {
    rank_long <- reshape2::melt(concordance$rank_matrix)
    colnames(rank_long) <- c("Subtype", "Sample", "Rank")

    concordance_plot <- ggplot(rank_long, aes(x = Sample, y = Subtype, fill = Rank)) +
      geom_tile() +
      geom_text(aes(label = round(Rank, 1)), size = 3) +
      scale_fill_viridis_c() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = sprintf("Per-Sample Subtype Ranking (Mean ρ = %.3f)",
                           concordance$mean_correlation),
           x = "Sample", y = "Subtype", fill = "Rank")

    ggsave(file.path(fig_dir, "02_trajectory_sample_concordance.pdf"),
           concordance_plot, width = 10, height = 8)
    ggsave(file.path(fig_dir, "02_trajectory_sample_concordance.png"),
           concordance_plot, width = 10, height = 8, dpi = 300)
    message("  Saved concordance plot")
  }
}

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

message("\n========================================")
message("Saving Outputs")
message("========================================\n")

objects_dir <- get_path(config, config$paths$objects_dir)
tables_dir <- get_path(config, config$paths$tables_dir)
ensure_dir(objects_dir)
ensure_dir(tables_dir)

# Save trajectory results
trajectory_output <- list(
  pooled = pooled_results,
  per_sample = per_sample_results,
  validation_pooled = validation_pooled,
  methods_used = pooled_results$methods_used,
  method_correlations = pooled_results$method_correlations
)

saveRDS(trajectory_output, file.path(objects_dir, "02_trajectory_results.rds"))
message(sprintf("Saved trajectory results: %s",
                file.path(objects_dir, "02_trajectory_results.rds")))

# Save updated seurat object with pseudotime
saveRDS(seurat_epi, file.path(objects_dir, "02_seurat_with_pseudotime.rds"))
message(sprintf("Saved Seurat object: %s",
                file.path(objects_dir, "02_seurat_with_pseudotime.rds")))

# Save pseudotime summary table
pt_summary <- seurat_epi@meta.data %>%
  group_by(module_score_subtype) %>%
  summarise(
    n_cells = n(),
    mean_pseudotime = mean(pseudotime_consensus, na.rm = TRUE),
    sd_pseudotime = sd(pseudotime_consensus, na.rm = TRUE),
    median_pseudotime = median(pseudotime_consensus, na.rm = TRUE),
    mean_confidence = mean(pseudotime_confidence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_pseudotime)

write.csv(pt_summary, file.path(tables_dir, "02_pseudotime_summary.csv"), row.names = FALSE)
message(sprintf("Saved pseudotime summary: %s",
                file.path(tables_dir, "02_pseudotime_summary.csv")))

# Save validation results
write.csv(validation_pooled, file.path(tables_dir, "02_trajectory_validation.csv"), row.names = FALSE)
message(sprintf("Saved validation results: %s",
                file.path(tables_dir, "02_trajectory_validation.csv")))

# Save method correlations
if (!is.null(pooled_results$method_correlations)) {
  write.csv(as.data.frame(pooled_results$method_correlations),
            file.path(tables_dir, "02_method_correlations.csv"))
  message(sprintf("Saved method correlations: %s",
                  file.path(tables_dir, "02_method_correlations.csv")))
}

# ==============================================================================
# STATISTICAL EVALUATION
# ==============================================================================

message("\n========================================")
message("Statistical Evaluation")
message("========================================\n")

# Define differentiation markers for validation
# Adjust these based on your tissue/system
diff_markers <- list(
  # Early/progenitor markers (expected to DECREASE along differentiation)
  early = c("SOX9", "KRT5", "KRT14", "TP63", "ITGA6", "ITGB1", "MKI67", "TOP2A"),
  # Late/mature markers (expected to INCREASE along differentiation)
  late = c("KRT20", "MUC2", "FABP1", "SI", "ALPI", "VIL1", "TFF3", "PIGR")
)

# Run statistical evaluation
trajectory_stats <- tryCatch({
  evaluate_trajectory_statistics(
    seurat_obj = seurat_epi,
    results = pooled_results,
    config = config,
    differentiation_markers = diff_markers,
    expected_order = expected_order
  )
}, error = function(e) {
  message(sprintf("  Could not complete statistical evaluation: %s", e$message))
  NULL
})

if (!is.null(trajectory_stats)) {
  # Save statistics as RDS
  stats_path <- file.path(objects_dir, "02_trajectory_statistics.rds")
  saveRDS(trajectory_stats, stats_path)
  message(sprintf("Saved trajectory statistics: %s", stats_path))

  # Generate and save markdown report
  report_path <- file.path(tables_dir, "02_trajectory_evaluation_report.md")
  tryCatch({
    generate_trajectory_report(trajectory_stats, report_path)
  }, error = function(e) {
    message(sprintf("  Could not generate report: %s", e$message))
  })

  # Save marker validation as CSV
  if (!is.null(trajectory_stats$marker_validation)) {
    marker_path <- file.path(tables_dir, "02_trajectory_marker_validation.csv")
    write.csv(trajectory_stats$marker_validation, marker_path, row.names = FALSE)
    message(sprintf("Saved marker validation: %s", marker_path))
  }

  # Save pseudotime by subtype
  if (!is.null(trajectory_stats$mean_pseudotime_by_subtype)) {
    pt_path <- file.path(tables_dir, "02_pseudotime_by_subtype.csv")
    write.csv(trajectory_stats$mean_pseudotime_by_subtype, pt_path, row.names = FALSE)
    message(sprintf("Saved pseudotime by subtype: %s", pt_path))
  }

  # Save cycling along trajectory
  if (!is.null(trajectory_stats$cycling_along_trajectory)) {
    cycling_path <- file.path(tables_dir, "02_cycling_along_trajectory.csv")
    write.csv(trajectory_stats$cycling_along_trajectory, cycling_path, row.names = FALSE)
    message(sprintf("Saved cycling along trajectory: %s", cycling_path))
  }

  # Add stats to output object
  trajectory_output$statistics <- trajectory_stats
  saveRDS(trajectory_output, file.path(objects_dir, "02_trajectory_results.rds"))
}

# Plot cell cycle vs trajectory (if cell cycle data available)
if ("cc_consensus" %in% colnames(seurat_epi@meta.data)) {
  message("\nCreating cell cycle vs trajectory plot...")

  cc_trajectory_plot <- tryCatch({
    plot_cell_cycle_trajectory(seurat_epi, config)
  }, error = function(e) {
    message(sprintf("  Could not generate plot: %s", e$message))
    NULL
  })

  if (!is.null(cc_trajectory_plot)) {
    cc_plot_path <- file.path(fig_dir, "02_cell_cycle_trajectory")
    ggsave(paste0(cc_plot_path, ".pdf"), cc_trajectory_plot, width = 14, height = 10)
    ggsave(paste0(cc_plot_path, ".png"), cc_trajectory_plot, width = 14, height = 10, dpi = 300)
    message(sprintf("Saved: %s.pdf/.png", cc_plot_path))
  }
}

# ==============================================================================
# COMPLETION
# ==============================================================================

message("\n")
message("================================================================")
message("  Trajectory Analysis Complete!")
message("================================================================")
message(sprintf("  Finished: %s", Sys.time()))
message(sprintf("  Epithelial cells analyzed: %d", ncol(seurat_epi)))
message(sprintf("  Methods used: %s", paste(pooled_results$methods_used, collapse = ", ")))

if (length(pooled_results$methods_used) > 1) {
  cor_vals <- pooled_results$method_correlations[upper.tri(pooled_results$method_correlations)]
  message(sprintf("  Method correlation: %.3f", mean(cor_vals, na.rm = TRUE)))
}

message(sprintf("  Validation: %s",
                ifelse(validation_pooled$concordant[1], "CONCORDANT", "NOT CONCORDANT")))

# Report quality score if available
if (!is.null(trajectory_stats) && !is.null(trajectory_stats$quality_score)) {
  message(sprintf("  Quality score: %.2f (%s)",
                  trajectory_stats$quality_score$overall,
                  trajectory_stats$quality_score$grade))
}

# Report cell cycle confounding if available
if (!is.null(trajectory_stats) && !is.null(trajectory_stats$cell_cycle_correlation)) {
  message(sprintf("  Cell cycle confounding: %s (ρ = %.3f)",
                  trajectory_stats$cell_cycle_correlation$interpretation,
                  trajectory_stats$cell_cycle_correlation$rho))
}

message("\nOutputs:")
message(sprintf("  Objects: %s", objects_dir))
message(sprintf("  Tables: %s", tables_dir))
message(sprintf("  Figures: %s", fig_dir))

message("\nKey output files:")
message("  - 02_trajectory_results.rds (full results with statistics)")
message("  - 02_trajectory_evaluation_report.md (statistical report)")
message("  - 02_cell_cycle_trajectory.pdf (cell cycle vs trajectory)")

message("\nNext step: Run 03_figure_generation.R (or generate_figures.R)")
message("================================================================\n")
