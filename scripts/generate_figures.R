#!/usr/bin/env Rscript
# ==============================================================================
# scripts/generate_figures.R
# ==============================================================================
# Script to regenerate specific figures without rerunning full pipeline
#
# Usage:
#   Rscript scripts/generate_figures.R                # Generate all figures
#   Rscript scripts/generate_figures.R 1             # Generate Figure 1 only
#   Rscript scripts/generate_figures.R 1 3           # Generate Figures 1 and 3
#   Rscript scripts/generate_figures.R supp          # Generate all supplementary
# ==============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set working directory
if (interactive()) {
  setwd(here::here())
} else {
  script_dir <- dirname(sys.frame(1)$ofile)
  setwd(file.path(script_dir, ".."))
}

# Load packages
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

# Source functions
source("R/utils/config.R")
source("R/figures/figure_generators.R")

# Load config
config <- load_config()

# Load pre-computed results
message("Loading pre-computed results...")

seurat_path <- file.path(get_path(config, config$paths$objects_dir), "seurat_classified.rds")
traj_path <- file.path(get_path(config, config$paths$objects_dir), "trajectory_results.rds")

if (!file.exists(seurat_path)) {
  stop("Classified Seurat object not found. Run 00_run_full_pipeline.R first.")
}

seurat_obj <- readRDS(seurat_path)

trajectory_results <- NULL
if (file.exists(traj_path)) {
  trajectory_results <- readRDS(traj_path)
}

# Determine which figures to generate
if (length(args) == 0) {
  # Generate all
  figures_to_generate <- c("1", "2", "3", "4", "supp")
} else {
  figures_to_generate <- args
}

message(paste("Generating figures:", paste(figures_to_generate, collapse = ", ")))

# Generate requested figures
for (fig in figures_to_generate) {
  
  tryCatch({
    
    if (fig == "1") {
      message("\n--- Figure 1 ---")
      generate_figure_1(seurat_obj, config)
      
    } else if (fig == "2") {
      message("\n--- Figure 2 ---")
      generate_figure_2(seurat_obj, config)
      
    } else if (fig == "3") {
      message("\n--- Figure 3 ---")
      if (is.null(trajectory_results)) {
        warning("Trajectory results not found, skipping Figure 3")
      } else {
        generate_figure_3(trajectory_results, config)
      }
      
    } else if (fig == "4") {
      message("\n--- Figure 4 ---")
      if (is.null(trajectory_results)) {
        warning("Trajectory results not found, skipping Figure 4")
      } else {
        generate_figure_4(seurat_obj, trajectory_results, config)
      }
      
    } else if (fig == "supp" || fig == "supplementary") {
      message("\n--- Supplementary Figures ---")
      generate_supp_classification(seurat_obj, config)
      if (!is.null(trajectory_results)) {
        generate_supp_trajectory(trajectory_results, config)
      }
      
    } else {
      warning(paste("Unknown figure:", fig))
    }
    
  }, error = function(e) {
    message(paste("ERROR generating figure", fig, ":", e$message))
  })
}

message("\nFigure generation complete!")
