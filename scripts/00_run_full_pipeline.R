#!/usr/bin/env Rscript
# ==============================================================================
# scripts/00_run_full_pipeline.R
# ==============================================================================
# Master script to run the complete ACP epithelial analysis pipeline
# 
# Usage:
#   Rscript scripts/00_run_full_pipeline.R
#   Or source in RStudio: source("scripts/00_run_full_pipeline.R")
#
# This script:
#   1. Loads configuration
#   2. Loads and preprocesses data
#   3. Runs module score classification
#   4. Runs trajectory analysis
#   5. Generates all figures
#   6. Saves results
# ==============================================================================

# Set working directory to project root
if (interactive()) {
  # In RStudio, use here package
  setwd(here::here())
} else {
  # From command line, find project root
  script_dir <- dirname(sys.frame(1)$ofile)
  setwd(file.path(script_dir, ".."))
}

message("\n")
message("================================================================")
message("  ACP Epithelial Differentiation Analysis Pipeline")
message("================================================================")
message(paste("  Started:", Sys.time()))
message(paste("  Working directory:", getwd()))
message("\n")

# ==============================================================================
# SETUP
# ==============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(yaml)
})

# Source all R functions
message("Loading project functions...")
source("R/utils/config.R")
source("R/classification/module_score_classification.R")
source("R/trajectory/trajectory_analysis.R")
source("R/figures/figure_generators.R")

# Load configuration
config <- load_config()
print_config_summary(config)

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================

message("\n========================================")
message("Step 1: Loading Data")
message("========================================\n")

# Check for data file
data_path <- get_path(config, config$paths$spatial_object_with_subtypes)

if (!file.exists(data_path)) {
  # Try alternative path
  data_path <- get_path(config, config$paths$spatial_object)
}

if (!file.exists(data_path)) {
  stop(paste("Data file not found:", data_path, 
             "\nPlease update paths in config/config.yaml"))
}

message(paste("Loading data from:", data_path))
seurat_obj <- readRDS(data_path)

message(sprintf("Loaded %d cells with %d genes", ncol(seurat_obj), nrow(seurat_obj)))

# Basic QC info
if ("sample" %in% colnames(seurat_obj@meta.data)) {
  message("Samples:")
  print(table(seurat_obj$sample))
}

# ==============================================================================
# STEP 2: CLASSIFICATION
# ==============================================================================

message("\n========================================")
message("Step 2: Module Score Classification")
message("========================================\n")

# Run classification pipeline
seurat_obj <- run_classification_pipeline(seurat_obj, config, 
                                          compare_with_binary = TRUE)

# Save intermediate result
ensure_dir(get_path(config, config$paths$objects_dir))
saveRDS(seurat_obj, 
        file.path(get_path(config, config$paths$objects_dir), 
                  "seurat_classified.rds"))
message("Saved classified object")

# ==============================================================================
# STEP 3: TRAJECTORY ANALYSIS
# ==============================================================================

message("\n========================================")
message("Step 3: Trajectory Analysis")
message("========================================\n")

# Check if we have enough cells for trajectory
n_epi <- sum(seurat_obj$module_score_subtype != "Non-epithelial", na.rm = TRUE)
message(sprintf("Epithelial cells for trajectory: %d", n_epi))

if (n_epi >= 100) {
  
  # Subset to epithelial cells
  seurat_epi <- subset(seurat_obj, 
                       module_score_subtype != "Non-epithelial" & 
                       !is.na(module_score_subtype))
  
  # Run pooled trajectory
  message("\n--- Pooled Analysis ---")
  cds_pooled <- tryCatch({
    run_trajectory_analysis(seurat_epi, config, "module_score_subtype")
  }, error = function(e) {
    message(paste("ERROR in pooled trajectory:", e$message))
    NULL
  })
  
  # Run per-sample trajectory if sample column exists
  trajectory_results <- list()
  
  if (!is.null(cds_pooled)) {
    trajectory_results$pooled <- list(
      cds = cds_pooled,
      skipped = FALSE
    )
  }
  
  if ("sample" %in% colnames(seurat_epi@meta.data)) {
    message("\n--- Per-Sample Analysis ---")
    per_sample_results <- tryCatch({
      run_per_sample_trajectory(seurat_epi, config, "sample", "module_score_subtype")
    }, error = function(e) {
      message(paste("ERROR in per-sample trajectory:", e$message))
      NULL
    })
    
    if (!is.null(per_sample_results)) {
      trajectory_results <- c(trajectory_results, per_sample_results)
    }
  }
  
  # Validate ordering
  if (length(trajectory_results) > 0) {
    expected_order <- c("Basal-like", "Transit-Amplifying", "Intermediate", "Specialized")
    validation <- validate_trajectory_ordering(trajectory_results, expected_order, config)
    
    # Save validation results
    save_table(validation, "trajectory_validation", config)
  }
  
  # Save trajectory results
  saveRDS(trajectory_results,
          file.path(get_path(config, config$paths$objects_dir),
                    "trajectory_results.rds"))
  message("Saved trajectory results")
  
} else {
  message("Insufficient epithelial cells for trajectory analysis")
  trajectory_results <- NULL
}

# ==============================================================================
# STEP 4: GENERATE FIGURES
# ==============================================================================

message("\n========================================")
message("Step 4: Generating Figures")
message("========================================\n")

if (!is.null(trajectory_results)) {
  figure_paths <- generate_all_figures(seurat_obj, trajectory_results, config)
  
  message("\nGenerated figures:")
  for (name in names(figure_paths)) {
    message(sprintf("  %s: %s", name, figure_paths[[name]][1]))
  }
}

# ==============================================================================
# STEP 5: SAVE FINAL RESULTS
# ==============================================================================

message("\n========================================")
message("Step 5: Saving Results")
message("========================================\n")

# Classification summary table
class_summary <- seurat_obj@meta.data %>%
  group_by(module_score_subtype) %>%
  summarise(
    n_cells = n(),
    pct = n() / nrow(seurat_obj@meta.data) * 100,
    mean_UMI = mean(nCount_Spatial, na.rm = TRUE),
    mean_genes = mean(nFeature_Spatial, na.rm = TRUE)
  ) %>%
  arrange(desc(n_cells))

save_table(class_summary, "classification_summary", config)

# Pseudotime summary if available
if (!is.null(trajectory_results) && !is.null(trajectory_results$pooled)) {
  pt_summary <- summarize_pseudotime(trajectory_results$pooled$cds)
  
  pt_stats <- pt_summary %>%
    group_by(subtype) %>%
    summarise(
      n = n(),
      mean_pseudotime = mean(pseudotime, na.rm = TRUE),
      sd_pseudotime = sd(pseudotime, na.rm = TRUE),
      median_pseudotime = median(pseudotime, na.rm = TRUE)
    ) %>%
    arrange(mean_pseudotime)
  
  save_table(pt_stats, "pseudotime_summary", config)
}

# Save session info for reproducibility
session_info <- sessionInfo()
saveRDS(session_info, 
        file.path(get_path(config, config$paths$objects_dir), "session_info.rds"))

# ==============================================================================
# DONE
# ==============================================================================

message("\n")
message("================================================================")
message("  Pipeline Complete!")
message("================================================================")
message(paste("  Finished:", Sys.time()))
message(paste("  Results saved to:", get_path(config, "results")))
message("\n")
