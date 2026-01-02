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
#   Rscript scripts/02_trajectory_analysis.R                           # Default: spatial
#   Rscript scripts/02_trajectory_analysis.R --dataset spatial         # Visium data
#   Rscript scripts/02_trajectory_analysis.R --dataset snrnaseq        # snRNA-seq GSE data
#   Rscript scripts/02_trajectory_analysis.R --dataset acp_scn         # ACP scRNA-seq data
#   Rscript scripts/02_trajectory_analysis.R --dataset merged          # Merged dataset
#   Rscript scripts/02_trajectory_analysis.R --dataset cellxgene       # CELLxGENE keratinocytes
#   Rscript scripts/02_trajectory_analysis.R --config path/to/config.yaml
#   Rscript scripts/02_trajectory_analysis.R --cores 16
#   Rscript scripts/02_trajectory_analysis.R --no-bootstrap            # Skip bootstrap (faster)
#
# ==============================================================================

# ==============================================================================
# SETUP AND ARGUMENT PARSING
# ==============================================================================

# Increase memory limit for parallel processing (20 GiB for large datasets)
options(future.globals.maxSize = 20 * 1024^3)  # 20 GiB

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
dataset_type <- "spatial"
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
  } else if (args[i] == "--dataset" && i < length(args)) {
    dataset_type <- tolower(args[i + 1])
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

# Validate dataset argument
valid_datasets <- c("spatial", "snrnaseq", "acp_scn", "merged", "cellxgene")
if (!dataset_type %in% valid_datasets) {
  stop(sprintf("Invalid --dataset value '%s'. Use one of: %s",
               dataset_type, paste(valid_datasets, collapse = ", ")))
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
  stop("Could not find project root directory.\n",
       "Please run this script from the project root or scripts/ directory.")
}
setwd(project_root)

message("\n")
message("================================================================")
message("  Step 2: Multi-Method Trajectory Analysis")
message(sprintf("  Dataset type: %s", toupper(dataset_type)))
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

# Load from Step 1 - with dataset-specific suffix
input_path <- file.path(get_path(config, config$paths$objects_dir),
                        sprintf("01_seurat_annotated_%s.rds", dataset_type))

if (!file.exists(input_path)) {
  stop(paste("Annotated object not found:", input_path,
             sprintf("\nPlease run: Rscript scripts/01_cell_type_annotation.R --dataset %s", dataset_type)))
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
# FILTER TO EPITHELIAL CELLS (if needed)
# ==============================================================================

message("\n========================================")
message("Preparing Epithelial Cells")
message("========================================\n")

# Check if we have is_epithelial_cell column
if ("is_epithelial_cell" %in% colnames(seurat_obj@meta.data)) {
  n_epi <- sum(seurat_obj$is_epithelial_cell, na.rm = TRUE)
  message(sprintf("Found %d epithelial cells (of %d total)", n_epi, ncol(seurat_obj)))

  # If object contains non-epithelial cells, filter
  if (n_epi < ncol(seurat_obj)) {
    message("Filtering to epithelial cells only...")
    seurat_epi <- subset(seurat_obj, is_epithelial_cell == TRUE)
  } else {
    seurat_epi <- seurat_obj
  }
} else {
  # Assume all cells are epithelial (from step 1 filtering)
  message("No is_epithelial_cell column - assuming all cells are epithelial")
  seurat_epi <- seurat_obj
}

message(sprintf("Proceeding with %d epithelial cells", ncol(seurat_epi)))

# ==============================================================================
# PREPARE DIMENSIONALITY REDUCTION
# ==============================================================================

message("\n========================================")
message("Checking Dimensionality Reduction")
message("========================================\n")

available_reductions <- names(seurat_epi@reductions)
message(sprintf("Available reductions: %s",
                ifelse(length(available_reductions) > 0,
                       paste(available_reductions, collapse = ", "),
                       "none")))

# Check if dataset is very large - if so, run preprocessing in serial mode
n_cells <- ncol(seurat_epi)
n_genes <- nrow(seurat_epi)
object_size_gb <- as.numeric(object.size(seurat_epi)) / 1024^3
is_large_dataset <- n_cells > 50000 || object_size_gb > 2

if (is_large_dataset) {
  message(sprintf("\nLarge dataset detected (%d cells, %.1f GB)", n_cells, object_size_gb))
  message("Running preprocessing in serial mode to avoid memory issues...")
  # Temporarily disable parallel processing for preprocessing
  if (requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    future::plan("sequential")
  }
}

# Run PCA if not present
if (!"pca" %in% tolower(available_reductions)) {
  message("\nRunning PCA...")

  # Check if data needs normalization
  needs_normalization <- tryCatch({
    available_layers <- Layers(seurat_epi[["RNA"]])
    if (!"data" %in% available_layers) {
      TRUE
    } else {
      # Check if data layer is same as counts (not normalized)
      counts_sample <- LayerData(seurat_epi, layer = "counts")[1:min(10, nrow(seurat_epi)), 1:min(10, ncol(seurat_epi))]
      data_sample <- LayerData(seurat_epi, layer = "data")[1:min(10, nrow(seurat_epi)), 1:min(10, ncol(seurat_epi))]
      all(counts_sample == data_sample)
    }
  }, error = function(e) TRUE)

  if (needs_normalization) {
    message("  Normalizing data...")
    seurat_epi <- NormalizeData(seurat_epi, verbose = FALSE)
  }

  message("  Finding variable features...")
  seurat_epi <- FindVariableFeatures(seurat_epi, verbose = FALSE)

  message("  Scaling data...")
  seurat_epi <- ScaleData(seurat_epi, verbose = FALSE)

  message("  Running PCA...")
  seurat_epi <- RunPCA(seurat_epi, npcs = 50, verbose = FALSE)
  message("  PCA complete")
}

# Run UMAP if not present
if (!"umap" %in% tolower(available_reductions)) {
  message("\nRunning UMAP...")
  seurat_epi <- RunUMAP(seurat_epi, dims = 1:30, verbose = FALSE)
  message("  UMAP complete")
}

# Restore parallel plan if we changed it
if (is_large_dataset && requireNamespace("future", quietly = TRUE)) {
  future::plan(old_plan)
  message("Restored parallel processing for trajectory analysis")
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

# ==============================================================================
# PER-SAMPLE ANALYSIS
# ==============================================================================

message("\n========================================")
message("Per-Sample Trajectory Analysis")
message("========================================\n")

# Skip per-sample analysis for cellxgene reference dataset
per_sample_results <- NULL
if (dataset_type == "cellxgene") {
  message("Skipping per-sample analysis for CELLxGENE reference dataset")
} else {
  # Find sample column
  sample_columns <- c("sample", "Sample", "sample_id", "Sample_ID", "orig.ident",
                      "patient", "Patient", "donor", "Donor")
  sample_col <- NULL
  for (col in sample_columns) {
    if (col %in% colnames(seurat_epi@meta.data)) {
      sample_col <- col
      break
    }
  }

  if (!is.null(sample_col)) {
    n_samples <- length(unique(seurat_epi@meta.data[[sample_col]]))
    message(sprintf("Found sample column: %s (%d samples)", sample_col, n_samples))

    if (n_samples >= 2) {
      per_sample_results <- tryCatch({
        run_per_sample_multi_trajectory(
          seurat_epi,
          config = config,
          methods = trajectory_methods,
          sample_column = sample_col,
          subtype_column = "module_score_subtype",
          parallel = (n_cores > 1)
        )
      }, error = function(e) {
        warning(sprintf("Per-sample analysis failed: %s", e$message))
        NULL
      })
    } else {
      message("Only 1 sample - skipping per-sample analysis")
    }
  } else {
    message("No sample column found - skipping per-sample analysis")
  }
}

# ==============================================================================
# VALIDATE TRAJECTORY
# ==============================================================================

message("\n========================================")
message("Validating Trajectory")
message("========================================\n")

# Get expected order from config
expected_order <- NULL
if (!is.null(config$trajectory$expected_order)) {
  expected_order <- config$trajectory$expected_order
}

# Validate pooled trajectory
validation_pooled <- tryCatch({
  validate_trajectory_order(
    seurat_epi,
    pooled_results,
    expected_order = expected_order
  )
}, error = function(e) {
  message(sprintf("  Validation warning: %s", e$message))
  NULL
})

# Validate per-sample if available
validation_per_sample <- NULL
if (!is.null(per_sample_results)) {
  validation_per_sample <- tryCatch({
    validate_trajectory_order(
      seurat_epi,
      per_sample_results,
      expected_order = expected_order,
      per_sample = TRUE
    )
  }, error = function(e) {
    message(sprintf("  Per-sample validation warning: %s", e$message))
    NULL
  })
}

# ==============================================================================
# COMPREHENSIVE STATISTICS
# ==============================================================================

message("\n========================================")
message("Computing Trajectory Statistics")
message("========================================\n")

trajectory_stats <- evaluate_trajectory_statistics(
  seurat_epi,
  pooled_results,
  config = config,
  run_permutation = !skip_permutation,
  n_permutations = n_permutations
)

# Run permutation test separately if requested
permutation_results <- NULL
if (!skip_permutation) {
  message("\nRunning permutation test...")
  permutation_results <- trajectory_stats$permutation_test
}

# ==============================================================================
# CELLXGENE-SPECIFIC ANALYSIS: Compare module_score_subtype vs cell_type
# ==============================================================================

cellxgene_comparison <- NULL
if (dataset_type == "cellxgene" && "cell_type" %in% colnames(seurat_epi@meta.data)) {
  message("\n========================================")
  message("CELLxGENE-Specific Analysis")
  message("========================================\n")

  message("Comparing pseudotime across classification schemes:")
  message("  1. module_score_subtype (our classification)")
  message("  2. cell_type (original CELLxGENE annotation)")

  # Pseudotime summary by module_score_subtype
  message("\n--- Pseudotime by Module Score Subtype ---")
  pt_by_subtype <- seurat_epi@meta.data %>%
    group_by(module_score_subtype) %>%
    summarise(
      n_cells = n(),
      pct_total = n() / nrow(seurat_epi@meta.data) * 100,
      mean_pseudotime = mean(pseudotime_consensus, na.rm = TRUE),
      sd_pseudotime = sd(pseudotime_consensus, na.rm = TRUE),
      median_pseudotime = median(pseudotime_consensus, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_pseudotime)

  print(pt_by_subtype)

  # Pseudotime summary by original cell_type
  message("\n--- Pseudotime by Original Cell Type ---")
  pt_by_celltype <- seurat_epi@meta.data %>%
    group_by(cell_type) %>%
    summarise(
      n_cells = n(),
      pct_total = n() / nrow(seurat_epi@meta.data) * 100,
      mean_pseudotime = mean(pseudotime_consensus, na.rm = TRUE),
      sd_pseudotime = sd(pseudotime_consensus, na.rm = TRUE),
      median_pseudotime = median(pseudotime_consensus, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_pseudotime)

  print(pt_by_celltype)

  # Cross-tabulation of classifications
  message("\n--- Cross-tabulation: Module Score Subtype vs Cell Type ---")
  cross_tab <- table(
    ModuleScore = seurat_epi$module_score_subtype,
    CellType = seurat_epi$cell_type
  )
  print(cross_tab)

  # Proportions within each cell_type
  message("\n--- Module Score Subtype proportions within each Cell Type ---")
  cross_prop <- prop.table(cross_tab, margin = 2) * 100
  print(round(cross_prop, 1))

  # Statistical test for association
  chi_test <- tryCatch({
    chisq.test(cross_tab)
  }, warning = function(w) {
    chisq.test(cross_tab, simulate.p.value = TRUE, B = 10000)
  }, error = function(e) NULL)

  if (!is.null(chi_test)) {
    message(sprintf("\nChi-square test (module_score_subtype × cell_type): X² = %.2f, p = %.2e",
                    chi_test$statistic, chi_test$p.value))
  }

  # Kruskal-Wallis test for pseudotime differences by cell_type
  kw_celltype <- tryCatch({
    kruskal.test(pseudotime_consensus ~ cell_type, data = seurat_epi@meta.data)
  }, error = function(e) NULL)

  if (!is.null(kw_celltype)) {
    message(sprintf("Kruskal-Wallis (pseudotime ~ cell_type): H = %.2f, p = %.2e",
                    kw_celltype$statistic, kw_celltype$p.value))
  }

  # Store results
  cellxgene_comparison <- list(
    pt_by_subtype = pt_by_subtype,
    pt_by_celltype = pt_by_celltype,
    cross_tabulation = cross_tab,
    cross_proportions = cross_prop,
    chi_square_test = chi_test,
    kruskal_wallis_celltype = kw_celltype
  )
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
  ggsave(file.path(fig_dir, sprintf("02_trajectory_method_comparison_%s.pdf", dataset_type)),
         comparison_plot, width = 14, height = 10)
  ggsave(file.path(fig_dir, sprintf("02_trajectory_method_comparison_%s.png", dataset_type)),
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
  plot_annotation(title = sprintf("Trajectory Analysis: Pseudotime on UMAP (%s)", toupper(dataset_type)))

ggsave(file.path(fig_dir, sprintf("02_trajectory_umap_%s.pdf", dataset_type)), umap_combined, width = 12, height = 12)
ggsave(file.path(fig_dir, sprintf("02_trajectory_umap_%s.png", dataset_type)), umap_combined, width = 12, height = 12, dpi = 300)
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

ggsave(file.path(fig_dir, sprintf("02_trajectory_pseudotime_violin_%s.pdf", dataset_type)), pt_violin, width = 10, height = 6)
ggsave(file.path(fig_dir, sprintf("02_trajectory_pseudotime_violin_%s.png", dataset_type)), pt_violin, width = 10, height = 6, dpi = 300)
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

    ggsave(file.path(fig_dir, sprintf("02_trajectory_sample_concordance_%s.pdf", dataset_type)),
           concordance_plot, width = 10, height = 8)
    ggsave(file.path(fig_dir, sprintf("02_trajectory_sample_concordance_%s.png", dataset_type)),
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

  ggsave(file.path(fig_dir, sprintf("02_trajectory_permutation_test_%s.pdf", dataset_type)), perm_plot, width = 10, height = 6)
  ggsave(file.path(fig_dir, sprintf("02_trajectory_permutation_test_%s.png", dataset_type)), perm_plot, width = 10, height = 6, dpi = 300)
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
    cc_plot_path <- file.path(fig_dir, sprintf("02_cell_cycle_trajectory_%s", dataset_type))
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

  ggsave(file.path(fig_dir, sprintf("02_method_quality_comparison_%s.pdf", dataset_type)), quality_plot, width = 10, height = 6)
  ggsave(file.path(fig_dir, sprintf("02_method_quality_comparison_%s.png", dataset_type)), quality_plot, width = 10, height = 6, dpi = 300)
  message("  Saved method quality plot")
}

# ==============================================================================
# CELLXGENE-SPECIFIC PLOTS
# ==============================================================================

if (dataset_type == "cellxgene" && !is.null(cellxgene_comparison)) {
  message("\nCreating CELLxGENE-specific comparison plots...")

  # Plot A: Pseudotime by cell_type (violin)
  celltype_violin <- ggplot(seurat_epi@meta.data,
                            aes(x = reorder(cell_type, pseudotime_consensus),
                                y = pseudotime_consensus, fill = cell_type)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Pseudotime by Original Cell Type (CELLxGENE)",
         x = "Cell Type (ordered by mean pseudotime)",
         y = "Consensus Pseudotime")

  # Plot B: Side-by-side comparison of classifications
  subtype_violin <- ggplot(seurat_epi@meta.data,
                           aes(x = reorder(module_score_subtype, pseudotime_consensus),
                               y = pseudotime_consensus, fill = module_score_subtype)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Pseudotime by Module Score Subtype",
         x = "Subtype (ordered by mean pseudotime)",
         y = "Consensus Pseudotime")

  # Plot C: UMAP colored by cell_type
  umap_celltype <- DimPlot(seurat_epi, group.by = "cell_type", pt.size = 0.3) +
    labs(title = "Original Cell Type (CELLxGENE)")

  # Plot D: Confusion-style heatmap
  cross_df <- as.data.frame(as.table(cellxgene_comparison$cross_proportions))
  colnames(cross_df) <- c("ModuleScore", "CellType", "Percentage")

  heatmap_comparison <- ggplot(cross_df, aes(x = CellType, y = ModuleScore, fill = Percentage)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f", Percentage)), size = 3) +
    scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue",
                         midpoint = 50, name = "% of\nCell Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()) +
    labs(title = "Module Score Classification within Each Cell Type",
         x = "Original Cell Type", y = "Module Score Subtype")

  # Combine cellxgene comparison plots
  cellxgene_combined <- (subtype_violin | celltype_violin) /
    (umap_plots$subtype | umap_celltype) /
    heatmap_comparison +
    plot_annotation(
      title = "CELLxGENE Reference: Classification Comparison",
      subtitle = sprintf("n = %d cells | Comparing module score classification vs original cell_type annotation",
                         ncol(seurat_epi)),
      tag_levels = "A"
    )

  ggsave(file.path(fig_dir, sprintf("02_cellxgene_classification_comparison_%s.pdf", dataset_type)),
         cellxgene_combined, width = 14, height = 16)
  ggsave(file.path(fig_dir, sprintf("02_cellxgene_classification_comparison_%s.png", dataset_type)),
         cellxgene_combined, width = 14, height = 16, dpi = 300)
  message("  Saved CELLxGENE classification comparison plot")

  # Additional: Pseudotime density by cell_type
  density_celltype <- ggplot(seurat_epi@meta.data, aes(x = pseudotime_consensus, fill = cell_type)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    labs(title = "Pseudotime Distribution by Original Cell Type",
         x = "Consensus Pseudotime", y = "Density", fill = "Cell Type")

  ggsave(file.path(fig_dir, sprintf("02_cellxgene_pseudotime_density_%s.pdf", dataset_type)),
         density_celltype, width = 10, height = 6)
  ggsave(file.path(fig_dir, sprintf("02_cellxgene_pseudotime_density_%s.png", dataset_type)),
         density_celltype, width = 10, height = 6, dpi = 300)
  message("  Saved CELLxGENE pseudotime density plot")
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
  validation_per_sample = validation_per_sample,
  methods_used = pooled_results$methods_used,
  method_correlations = pooled_results$method_correlations,
  method_quality = pooled_results$method_quality,
  permutation_test = permutation_results,
  statistics = trajectory_stats,
  cellxgene_comparison = cellxgene_comparison,  # NULL for non-cellxgene datasets
  config_used = list(
    dataset_type = dataset_type,
    n_cores = n_cores,
    n_bootstrap = n_bootstrap,
    n_permutations = n_permutations,
    methods = trajectory_methods
  )
)

saveRDS(trajectory_output, file.path(objects_dir, sprintf("02_trajectory_results_%s.rds", dataset_type)))
message(sprintf("Saved trajectory results: %s",
                file.path(objects_dir, sprintf("02_trajectory_results_%s.rds", dataset_type))))

# Save updated seurat object with pseudotime
saveRDS(seurat_epi, file.path(objects_dir, sprintf("02_seurat_with_pseudotime_%s.rds", dataset_type)))
message(sprintf("Saved Seurat object: %s",
                file.path(objects_dir, sprintf("02_seurat_with_pseudotime_%s.rds", dataset_type))))

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

write.csv(pt_summary, file.path(tables_dir, sprintf("02_pseudotime_summary_%s.csv", dataset_type)), row.names = FALSE)
message(sprintf("Saved pseudotime summary: %s",
                file.path(tables_dir, sprintf("02_pseudotime_summary_%s.csv", dataset_type))))

# Save validation results
if (!is.null(validation_pooled)) {
  write.csv(validation_pooled, file.path(tables_dir, sprintf("02_trajectory_validation_pooled_%s.csv", dataset_type)), row.names = FALSE)
  message(sprintf("Saved validation results: %s",
                  file.path(tables_dir, sprintf("02_trajectory_validation_pooled_%s.csv", dataset_type))))
}

if (!is.null(validation_per_sample)) {
  write.csv(validation_per_sample, file.path(tables_dir, sprintf("02_trajectory_validation_per_sample_%s.csv", dataset_type)), row.names = FALSE)
  message(sprintf("Saved per-sample validation: %s",
                  file.path(tables_dir, sprintf("02_trajectory_validation_per_sample_%s.csv", dataset_type))))
}

# Save method correlations
if (!is.null(pooled_results$method_correlations)) {
  write.csv(as.data.frame(pooled_results$method_correlations),
            file.path(tables_dir, sprintf("02_method_correlations_%s.csv", dataset_type)))
  message(sprintf("Saved method correlations: %s",
                  file.path(tables_dir, sprintf("02_method_correlations_%s.csv", dataset_type))))
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
  write.csv(quality_summary, file.path(tables_dir, sprintf("02_method_quality_scores_%s.csv", dataset_type)), row.names = FALSE)
  message(sprintf("Saved method quality scores: %s",
                  file.path(tables_dir, sprintf("02_method_quality_scores_%s.csv", dataset_type))))
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
  write.csv(perm_summary, file.path(tables_dir, sprintf("02_permutation_test_results_%s.csv", dataset_type)), row.names = FALSE)
  message(sprintf("Saved permutation test results: %s",
                  file.path(tables_dir, sprintf("02_permutation_test_results_%s.csv", dataset_type))))
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
  write.csv(boot_summary, file.path(tables_dir, sprintf("02_bootstrap_confidence_intervals_%s.csv", dataset_type)), row.names = FALSE)
  message(sprintf("Saved bootstrap CIs: %s",
                  file.path(tables_dir, sprintf("02_bootstrap_confidence_intervals_%s.csv", dataset_type))))
}

# Save marker validation
if (!is.null(trajectory_stats$marker_validation)) {
  write.csv(trajectory_stats$marker_validation,
            file.path(tables_dir, sprintf("02_marker_validation_%s.csv", dataset_type)), row.names = FALSE)
  message(sprintf("Saved marker validation: %s",
                  file.path(tables_dir, sprintf("02_marker_validation_%s.csv", dataset_type))))
}

# Save heterogeneity statistics (if available)
if (!is.null(per_sample_results)) {
  heterogeneity <- attr(per_sample_results, "heterogeneity")
  if (!is.null(heterogeneity)) {
    saveRDS(heterogeneity, file.path(objects_dir, sprintf("02_sample_heterogeneity_%s.rds", dataset_type)))
    message(sprintf("Saved heterogeneity statistics: %s",
                    file.path(objects_dir, sprintf("02_sample_heterogeneity_%s.rds", dataset_type))))
  }
}

# Generate and save markdown report
report_path <- file.path(tables_dir, sprintf("02_trajectory_evaluation_report_%s.md", dataset_type))
generate_trajectory_report(trajectory_stats, report_path)

# ==============================================================================
# CELLXGENE-SPECIFIC OUTPUTS
# ==============================================================================

if (dataset_type == "cellxgene" && !is.null(cellxgene_comparison)) {
  message("\nSaving CELLxGENE-specific outputs...")

  # Pseudotime by cell_type
  write.csv(cellxgene_comparison$pt_by_celltype,
            file.path(tables_dir, sprintf("02_pseudotime_by_celltype_%s.csv", dataset_type)),
            row.names = FALSE)
  message(sprintf("Saved: 02_pseudotime_by_celltype_%s.csv", dataset_type))

  # Cross-tabulation (counts)
  cross_tab_df <- as.data.frame.matrix(cellxgene_comparison$cross_tabulation)
  cross_tab_df$module_score_subtype <- rownames(cross_tab_df)
  cross_tab_df <- cross_tab_df[, c("module_score_subtype", setdiff(colnames(cross_tab_df), "module_score_subtype"))]
  write.csv(cross_tab_df,
            file.path(tables_dir, sprintf("02_classification_crosstab_%s.csv", dataset_type)),
            row.names = FALSE)
  message(sprintf("Saved: 02_classification_crosstab_%s.csv", dataset_type))

  # Cross-tabulation (proportions within cell_type)
  cross_prop_df <- as.data.frame.matrix(cellxgene_comparison$cross_proportions)
  cross_prop_df$module_score_subtype <- rownames(cross_prop_df)
  cross_prop_df <- cross_prop_df[, c("module_score_subtype", setdiff(colnames(cross_prop_df), "module_score_subtype"))]
  write.csv(cross_prop_df,
            file.path(tables_dir, sprintf("02_classification_proportions_%s.csv", dataset_type)),
            row.names = FALSE)
  message(sprintf("Saved: 02_classification_proportions_%s.csv", dataset_type))

  # Statistical comparison summary
  cellxgene_stats <- data.frame(
    test = c("Chi-square (module_score × cell_type)",
             "Kruskal-Wallis (pseudotime ~ cell_type)"),
    statistic = c(
      ifelse(!is.null(cellxgene_comparison$chi_square_test),
             cellxgene_comparison$chi_square_test$statistic, NA),
      ifelse(!is.null(cellxgene_comparison$kruskal_wallis_celltype),
             cellxgene_comparison$kruskal_wallis_celltype$statistic, NA)
    ),
    p_value = c(
      ifelse(!is.null(cellxgene_comparison$chi_square_test),
             cellxgene_comparison$chi_square_test$p.value, NA),
      ifelse(!is.null(cellxgene_comparison$kruskal_wallis_celltype),
             cellxgene_comparison$kruskal_wallis_celltype$p.value, NA)
    )
  )
  write.csv(cellxgene_stats,
            file.path(tables_dir, sprintf("02_cellxgene_statistical_tests_%s.csv", dataset_type)),
            row.names = FALSE)
  message(sprintf("Saved: 02_cellxgene_statistical_tests_%s.csv", dataset_type))
}

# ==============================================================================
# COMPLETION
# ==============================================================================

message("\n")
message("================================================================")
message("  Trajectory Analysis Complete!")
message("================================================================")
message(sprintf("  Dataset type: %s", toupper(dataset_type)))
message(sprintf("  Finished: %s", Sys.time()))
message(sprintf("  Epithelial cells analyzed: %d", ncol(seurat_epi)))
message(sprintf("  Methods used: %s", paste(pooled_results$methods_used, collapse = ", ")))
message(sprintf("  Quality score: %.2f (%s)",
                trajectory_stats$quality_score$overall,
                trajectory_stats$quality_score$grade))
if (!is.null(per_sample_results)) {
  concordance <- attr(per_sample_results, "concordance")
  if (!is.null(concordance)) {
    message(sprintf("  Cross-sample concordance: ρ = %.3f", concordance$mean_correlation))
  }
}
if (dataset_type == "cellxgene" && !is.null(cellxgene_comparison)) {
  message("\n  CELLxGENE-specific analysis:")
  message(sprintf("    Original cell_type categories: %d", length(unique(seurat_epi$cell_type))))
  message(sprintf("    Module score subtypes: %d", length(unique(seurat_epi$module_score_subtype))))
}
message("\nOutputs:")
message(sprintf("  Objects: %s", objects_dir))
message(sprintf("  Tables:  %s", tables_dir))
message(sprintf("  Figures: %s", fig_dir))
message("\nKey files:")
message(sprintf("  - 02_seurat_with_pseudotime_%s.rds", dataset_type))
message(sprintf("  - 02_trajectory_results_%s.rds", dataset_type))
message(sprintf("  - 02_pseudotime_summary_%s.csv", dataset_type))
if (dataset_type == "cellxgene") {
  message(sprintf("  - 02_pseudotime_by_celltype_%s.csv", dataset_type))
  message(sprintf("  - 02_classification_crosstab_%s.csv", dataset_type))
  message(sprintf("  - 02_cellxgene_classification_comparison_%s.pdf", dataset_type))
}
message("================================================================\n")
