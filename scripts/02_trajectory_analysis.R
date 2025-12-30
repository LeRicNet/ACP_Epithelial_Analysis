#!/usr/bin/env Rscript
# ==============================================================================
# scripts/02_trajectory_analysis.R
# ==============================================================================
# Multi-Method Trajectory Analysis with HPC Parallelization and
# Comprehensive Statistical Evaluation
#
# This script performs trajectory inference using multiple methods:
#   1. Monocle3 - graph-based trajectory
#   2. Slingshot - principal curves through clusters
#   3. Diffusion pseudotime (optional)
#
# Statistical Framework:
#   - Per-method quality metrics
#   - Per-sample heterogeneity testing
#   - Bootstrap confidence intervals
#   - Permutation significance testing
#
# HPC Features:
#   - Automatic SLURM detection
#   - Parallel per-sample analysis
#   - Parallelized bootstrap/permutation tests
#
# Input:
#   - Annotated Seurat object from 01_cell_type_annotation.R
#
# Output:
#   - Trajectory results with consensus pseudotime
#   - Per-sample validation with heterogeneity statistics
#   - Bootstrap confidence intervals
#   - Permutation test results
#   - Comprehensive diagnostic plots
#
# Usage:
#   Rscript scripts/02_trajectory_analysis.R
#   Rscript scripts/02_trajectory_analysis.R --config path/to/config.yaml
#   Rscript scripts/02_trajectory_analysis.R --cores 16
#   Rscript scripts/02_trajectory_analysis.R --no-bootstrap  # Skip bootstrap (faster)
#
# ==============================================================================

# ==============================================================================
# SETUP AND ARGUMENT PARSING
# ==============================================================================

# Increase memory limit for parallel processing
options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
n_cores_override <- NULL
skip_bootstrap <- FALSE
skip_permutation <- FALSE
n_bootstrap <- 100
n_permutations <- 1000

i <- 1
while (i <= length(args)) {
  if (args[i] == "--config" && i < length(args)) {
    config_path <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--cores" && i < length(args)) {
    n_cores_override <- as.integer(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--no-bootstrap") {
    skip_bootstrap <- TRUE
    i <- i + 1
  } else if (args[i] == "--no-permutation") {
    skip_permutation <- TRUE
    i <- i + 1
  } else if (args[i] == "--n-bootstrap" && i < length(args)) {
    n_bootstrap <- as.integer(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--n-permutations" && i < length(args)) {
    n_permutations <- as.integer(args[i + 1])
    i <- i + 2
  } else if (!startsWith(args[i], "--")) {
    config_path <- args[i]
    i <- i + 1
  } else {
    i <- i + 1
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
message("  With HPC Parallelization & Comprehensive Statistics")
message("================================================================")
message(paste("  Started:", Sys.time()))
message(paste("  Working directory:", getwd()))
message("\n")

# Detect HPC environment
slurm_job_id <- Sys.getenv("SLURM_JOB_ID", unset = NA)
slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)
pbs_job_id <- Sys.getenv("PBS_JOBID", unset = NA)

if (!is.na(slurm_job_id)) {
  message(sprintf("  Detected SLURM environment (Job ID: %s)", slurm_job_id))
  if (!is.na(slurm_cpus)) {
    message(sprintf("  SLURM CPUs allocated: %s", slurm_cpus))
  }
} else if (!is.na(pbs_job_id)) {
  message(sprintf("  Detected PBS environment (Job ID: %s)", pbs_job_id))
}

# Load packages
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Check for parallel packages
parallel_available <- list(
  BiocParallel = requireNamespace("BiocParallel", quietly = TRUE),
  future = requireNamespace("future", quietly = TRUE),
  future.apply = requireNamespace("future.apply", quietly = TRUE)
)

message("\nParallel packages available:")
for (pkg in names(parallel_available)) {
  status <- ifelse(parallel_available[[pkg]], "OK", "MISSING")
  message(sprintf("  %s: %s", pkg, status))
}

# Source utility functions
source("R/utils/config.R")
source("R/trajectory/trajectory_analysis.R")

# Load configuration
config <- load_config(config_path)
print_config_summary(config)

# Determine number of cores
if (!is.null(n_cores_override)) {
  n_cores <- n_cores_override
} else if (!is.na(slurm_cpus)) {
  n_cores <- as.integer(slurm_cpus)
} else if (!is.null(config$reproducibility$n_cores)) {
  n_cores <- config$reproducibility$n_cores
} else {
  n_cores <- max(1, parallel::detectCores() - 1)
}

message(sprintf("\nUsing %d cores for parallel processing", n_cores))

# Set seed
set.seed(config$reproducibility$seed)

# Setup parallel backend
if (n_cores > 1 && parallel_available$BiocParallel && parallel_available$future) {
  backend <- setup_parallel_backend(n_cores = n_cores)
} else {
  message("Running in serial mode")
  backend <- list(n_cores = 1, type = "serial")
}

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

# Update object for Seurat version compatibility (fixes VisiumV1 slot issues)
if (inherits(seurat_obj, "Seurat")) {
  seurat_obj <- UpdateSeuratObject(seurat_obj)
  message("  Seurat object updated for version compatibility")
}

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
message(sprintf("Bootstrap CI: %s (n=%d)", ifelse(skip_bootstrap, "NO", "YES"), n_bootstrap))

# Run pooled analysis
pooled_results <- run_multi_method_trajectory(
  seurat_epi,
  config = config,
  methods = trajectory_methods,
  subtype_column = "module_score_subtype",
  n_cores = n_cores,
  bootstrap_ci = !skip_bootstrap,
  n_bootstrap = n_bootstrap
)

# Update seurat object with pseudotime
seurat_epi <- pooled_results$seurat_obj

message("\nPooled analysis complete!")
message(sprintf("  Methods used: %s", paste(pooled_results$methods_used, collapse = ", ")))

# Report method quality scores
if (!is.null(pooled_results$method_quality)) {
  message("\nPer-method quality scores:")
  for (method in names(pooled_results$method_quality)) {
    mq <- pooled_results$method_quality[[method]]
    message(sprintf("  %s: %.3f", method, mq$overall_quality))
  }
}

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
# PERMUTATION TEST FOR TRAJECTORY SIGNIFICANCE
# ==============================================================================

permutation_results <- NULL

if (!skip_permutation) {
  message("\n========================================")
  message("Permutation Test for Trajectory Significance")
  message("========================================\n")

  permutation_results <- permutation_test_trajectory(
    seurat_epi,
    subtype_column = "module_score_subtype",
    expected_order = expected_order,
    n_permutations = n_permutations,
    parallel = n_cores > 1
  )

  message(sprintf("\nTrajectory significance: %s", permutation_results$interpretation))

  # Also test method concordance if multiple methods
  if (length(pooled_results$methods_used) > 1) {
    method_concordance_perm <- permutation_test_method_concordance(
      pooled_results$pseudotime_matrix,
      n_permutations = n_permutations,
      parallel = n_cores > 1
    )

    pooled_results$method_concordance_permutation <- method_concordance_perm
    message(sprintf("\nMethod concordance: %s (p=%.4f)",
                    ifelse(method_concordance_perm$significant, "Significant", "Not significant"),
                    method_concordance_perm$p_value))
  }
}

# ==============================================================================
# PER-SAMPLE ANALYSIS (PARALLELIZED)
# ==============================================================================

per_sample_enabled <- TRUE
if (!is.null(config$trajectory$per_sample$enabled)) {
  per_sample_enabled <- config$trajectory$per_sample$enabled
}

per_sample_results <- NULL

if (per_sample_enabled && "sample" %in% colnames(seurat_epi@meta.data)) {

  message("\n========================================")
  message("Running Per-Sample Analysis (Parallelized)")
  message("========================================\n")

  per_sample_results <- run_per_sample_multi_trajectory(
    seurat_epi,
    config = config,
    methods = trajectory_methods,
    sample_column = "sample",
    subtype_column = "module_score_subtype",
    parallel = n_cores > 1
  )

  # Validate per-sample
  validation_per_sample <- validate_trajectory_order(per_sample_results, expected_order)

  message("\nPer-sample validation:")
  print(validation_per_sample)

  # Get heterogeneity statistics
  heterogeneity <- attr(per_sample_results, "heterogeneity")
  if (!is.null(heterogeneity)) {
    message("\nSample heterogeneity summary:")
    message(sprintf("  Heterogeneous: %s",
                    ifelse(heterogeneity$summary$is_heterogeneous, "YES", "NO")))
    if (length(heterogeneity$summary$reasons) > 0) {
      for (reason in heterogeneity$summary$reasons) {
        message(sprintf("    - %s", reason))
      }
    }
  }

} else {
  message("\nSkipping per-sample analysis (disabled or no sample column)")
}

# ==============================================================================
# COMPREHENSIVE STATISTICAL EVALUATION
# ==============================================================================

message("\n========================================")
message("Comprehensive Statistical Evaluation")
message("========================================\n")

# Define differentiation markers for validation
diff_markers <- list(
  early = c("SOX9", "KRT5", "KRT14", "TP63", "ITGA6", "ITGB1", "MKI67", "TOP2A"),
  late = c("KRT20", "MUC2", "FABP1", "SI", "ALPI", "VIL1", "TFF3", "PIGR")
)

# Run comprehensive evaluation
trajectory_stats <- evaluate_trajectory_statistics(
  seurat_obj = seurat_epi,
  results = pooled_results,
  config = config,
  differentiation_markers = diff_markers,
  expected_order = expected_order,
  run_permutation = !skip_permutation,
  n_permutations = n_permutations
)

# Add permutation results if run separately
if (!is.null(permutation_results)) {
  trajectory_stats$permutation_test <- permutation_results
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

# Bootstrap CI width (if available)
if ("pseudotime_ci_width" %in% colnames(seurat_epi@meta.data)) {
  umap_plots$ci_width <- FeaturePlot(seurat_epi, features = "pseudotime_ci_width",
                                     pt.size = 0.5) +
    scale_color_viridis_c(option = "magma") +
    labs(title = "Bootstrap CI Width", color = "CI Width")
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

# Plot 3: Pseudotime by subtype with CI
message("Creating pseudotime distribution plots...")

pt_data <- data.frame(
  pseudotime = seurat_epi$pseudotime_consensus,
  subtype = seurat_epi$module_score_subtype,
  confidence = seurat_epi$pseudotime_confidence
)

if ("pseudotime_ci_width" %in% colnames(seurat_epi@meta.data)) {
  pt_data$ci_width <- seurat_epi$pseudotime_ci_width
}

pt_violin <- ggplot(pt_data,
                    aes(x = reorder(subtype, pseudotime),
                        y = pseudotime, fill = subtype)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Consensus Pseudotime Distribution",
       subtitle = sprintf("Quality Score: %.2f (%s)",
                          trajectory_stats$quality_score$overall,
                          trajectory_stats$quality_score$grade),
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

# Plot 5: Permutation test results (if available)
if (!is.null(permutation_results)) {
  message("Creating permutation test plot...")

  perm_data <- data.frame(correlation = permutation_results$permutation_correlations)

  perm_plot <- ggplot(perm_data, aes(x = correlation)) +
    geom_histogram(bins = 50, fill = "gray70", color = "white") +
    geom_vline(xintercept = permutation_results$observed_correlation,
               color = "red", linewidth = 1.5, linetype = "solid") +
    annotate("text",
             x = permutation_results$observed_correlation,
             y = Inf,
             label = sprintf("Observed\nρ = %.3f", permutation_results$observed_correlation),
             vjust = 2, hjust = -0.1, color = "red", fontface = "bold") +
    theme_minimal() +
    labs(title = "Permutation Test for Trajectory Significance",
         subtitle = sprintf("p = %.4f, Effect size (d) = %.2f | %s",
                            permutation_results$p_value,
                            permutation_results$effect_size,
                            permutation_results$interpretation),
         x = "Correlation with Expected Order",
         y = "Permutation Count")

  ggsave(file.path(fig_dir, "02_trajectory_permutation_test.pdf"), perm_plot, width = 10, height = 6)
  ggsave(file.path(fig_dir, "02_trajectory_permutation_test.png"), perm_plot, width = 10, height = 6, dpi = 300)
  message("  Saved permutation test plot")
}

# Plot 6: Cell cycle vs trajectory
cc_col <- ifelse("cc_consensus" %in% colnames(seurat_epi@meta.data), "cc_consensus", "Phase")
if (cc_col %in% colnames(seurat_epi@meta.data)) {
  message("Creating cell cycle vs trajectory plot...")

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
    message("  Saved cell cycle trajectory plot")
  }
}

# Plot 7: Method quality comparison
if (!is.null(pooled_results$method_quality) && length(pooled_results$method_quality) > 1) {
  message("Creating method quality comparison plot...")

  quality_df <- do.call(rbind, lapply(names(pooled_results$method_quality), function(m) {
    mq <- pooled_results$method_quality[[m]]
    data.frame(
      method = m,
      metric = names(mq$quality_components),
      value = as.numeric(mq$quality_components)
    )
  }))

  quality_plot <- ggplot(quality_df, aes(x = method, y = value, fill = metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    labs(title = "Method Quality Components",
         x = "Method", y = "Score", fill = "Component") +
    ylim(0, 1)

  ggsave(file.path(fig_dir, "02_method_quality_comparison.pdf"), quality_plot, width = 10, height = 6)
  ggsave(file.path(fig_dir, "02_method_quality_comparison.png"), quality_plot, width = 10, height = 6, dpi = 300)
  message("  Saved method quality plot")
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
  validation_per_sample = if (!is.null(per_sample_results)) validation_per_sample else NULL,
  methods_used = pooled_results$methods_used,
  method_correlations = pooled_results$method_correlations,
  method_quality = pooled_results$method_quality,
  permutation_test = permutation_results,
  statistics = trajectory_stats,
  config_used = list(
    n_cores = n_cores,
    n_bootstrap = n_bootstrap,
    n_permutations = n_permutations,
    methods = trajectory_methods
  )
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
write.csv(validation_pooled, file.path(tables_dir, "02_trajectory_validation_pooled.csv"), row.names = FALSE)
message(sprintf("Saved validation results: %s",
                file.path(tables_dir, "02_trajectory_validation_pooled.csv")))

if (!is.null(per_sample_results)) {
  write.csv(validation_per_sample, file.path(tables_dir, "02_trajectory_validation_per_sample.csv"), row.names = FALSE)
  message(sprintf("Saved per-sample validation: %s",
                  file.path(tables_dir, "02_trajectory_validation_per_sample.csv")))
}

# Save method correlations
if (!is.null(pooled_results$method_correlations)) {
  write.csv(as.data.frame(pooled_results$method_correlations),
            file.path(tables_dir, "02_method_correlations.csv"))
  message(sprintf("Saved method correlations: %s",
                  file.path(tables_dir, "02_method_correlations.csv")))
}

# Save method quality scores
if (!is.null(pooled_results$method_quality)) {
  quality_summary <- do.call(rbind, lapply(names(pooled_results$method_quality), function(m) {
    mq <- pooled_results$method_quality[[m]]
    data.frame(
      method = m,
      overall_quality = mq$overall_quality,
      separation = ifelse(is.null(mq$separation), NA, mq$separation),
      continuity = ifelse(is.null(mq$continuity), NA, mq$continuity)
    )
  }))
  write.csv(quality_summary, file.path(tables_dir, "02_method_quality_scores.csv"), row.names = FALSE)
  message(sprintf("Saved method quality scores: %s",
                  file.path(tables_dir, "02_method_quality_scores.csv")))
}

# Save permutation test results
if (!is.null(permutation_results)) {
  perm_summary <- data.frame(
    observed_correlation = permutation_results$observed_correlation,
    p_value = permutation_results$p_value,
    p_value_two_tailed = permutation_results$p_value_two_tailed,
    effect_size = permutation_results$effect_size,
    n_permutations = permutation_results$n_permutations,
    significant = permutation_results$significant,
    interpretation = permutation_results$interpretation
  )
  write.csv(perm_summary, file.path(tables_dir, "02_permutation_test_results.csv"), row.names = FALSE)
  message(sprintf("Saved permutation test results: %s",
                  file.path(tables_dir, "02_permutation_test_results.csv")))
}

# Save bootstrap summary (if available)
if (!is.null(pooled_results$bootstrap)) {
  boot_summary <- data.frame(
    cell = colnames(seurat_epi),
    pseudotime = seurat_epi$pseudotime_consensus,
    ci_lower = pooled_results$bootstrap$ci_lower,
    ci_upper = pooled_results$bootstrap$ci_upper,
    ci_width = pooled_results$bootstrap$ci_width,
    boot_se = pooled_results$bootstrap$boot_se
  )
  write.csv(boot_summary, file.path(tables_dir, "02_bootstrap_confidence_intervals.csv"), row.names = FALSE)
  message(sprintf("Saved bootstrap CIs: %s",
                  file.path(tables_dir, "02_bootstrap_confidence_intervals.csv")))
}

# Save marker validation
if (!is.null(trajectory_stats$marker_validation)) {
  write.csv(trajectory_stats$marker_validation,
            file.path(tables_dir, "02_marker_validation.csv"), row.names = FALSE)
  message(sprintf("Saved marker validation: %s",
                  file.path(tables_dir, "02_marker_validation.csv")))
}

# Save heterogeneity statistics (if available)
if (!is.null(per_sample_results)) {
  heterogeneity <- attr(per_sample_results, "heterogeneity")
  if (!is.null(heterogeneity)) {
    saveRDS(heterogeneity, file.path(objects_dir, "02_sample_heterogeneity.rds"))
    message(sprintf("Saved heterogeneity statistics: %s",
                    file.path(objects_dir, "02_sample_heterogeneity.rds")))
  }
}

# Generate and save markdown report
report_path <- file.path(tables_dir, "02_trajectory_evaluation_report.md")
generate_trajectory_report(trajectory_stats, report_path)

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
  message(sprintf("  Method correlation: %.3f", trajectory_stats$mean_correlation))
}

message(sprintf("  Validation: %s",
                ifelse(validation_pooled$concordant[1], "CONCORDANT", "NOT CONCORDANT")))

# Report quality score
message(sprintf("  Quality score: %.2f (%s)",
                trajectory_stats$quality_score$overall,
                trajectory_stats$quality_score$grade))

# Report permutation test
if (!is.null(permutation_results)) {
  message(sprintf("  Permutation test: p=%.4f (%s)",
                  permutation_results$p_value,
                  permutation_results$interpretation))
}

# Report cell cycle confounding
if (!is.null(trajectory_stats$cell_cycle_correlation)) {
  message(sprintf("  Cell cycle confounding: %s (ρ = %.3f)",
                  trajectory_stats$cell_cycle_correlation$interpretation,
                  trajectory_stats$cell_cycle_correlation$rho))
}

# Report bootstrap if available
if (!is.null(pooled_results$bootstrap)) {
  message(sprintf("  Bootstrap CI mean width: %.3f",
                  mean(pooled_results$bootstrap$ci_width, na.rm = TRUE)))
}

# Report sample heterogeneity
if (!is.null(per_sample_results)) {
  heterogeneity <- attr(per_sample_results, "heterogeneity")
  if (!is.null(heterogeneity)) {
    message(sprintf("  Sample heterogeneity: %s",
                    ifelse(heterogeneity$summary$is_heterogeneous, "DETECTED", "NOT DETECTED")))
  }
}

message("\nOutputs:")
message(sprintf("  Objects: %s", objects_dir))
message(sprintf("  Tables: %s", tables_dir))
message(sprintf("  Figures: %s", fig_dir))

message("\nKey output files:")
message("  - 02_trajectory_results.rds (full results)")
message("  - 02_trajectory_evaluation_report.md (statistical report)")
message("  - 02_permutation_test_results.csv (significance test)")
message("  - 02_bootstrap_confidence_intervals.csv (CIs)")
message("  - 02_method_quality_scores.csv (per-method evaluation)")

message("\nNext step: Run 03_figure_generation.R (or generate_figures.R)")
message("================================================================\n")
