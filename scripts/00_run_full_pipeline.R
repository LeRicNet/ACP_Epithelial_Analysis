#!/usr/bin/env Rscript
# ==============================================================================
# 00_run_full_pipeline.R
# ==============================================================================
# Master script to run the complete ACP Epithelial Differentiation Analysis
# pipeline. Orchestrates all steps from data preparation through trajectory
# analysis and figure generation.
#
# Usage:
#   Rscript scripts/00_run_full_pipeline.R
#   Rscript scripts/00_run_full_pipeline.R --dataset merged --skip-reference
#   Rscript scripts/00_run_full_pipeline.R --help
#
# Steps:
#   0a. Download reference data (optional, one-time)
#   0.  Merge datasets (if using merged analysis)
#   1.  Cell type annotation
#   1b. Reference validation and enhanced statistics
#   2.  Trajectory analysis
#   3.  Figure generation
#
# Author: [Your Name]
# Date: 2025-12-30
# ==============================================================================

# --- Setup --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
})

# --- Command Line Arguments ---------------------------------------------------

option_list <- list(
  make_option(c("-d", "--dataset"), type = "character", default = "merged",
              help = "Dataset to analyze: spatial, snrnaseq, acp_scn, merged [default: %default]"),
  make_option(c("-c", "--config"), type = "character", default = "config/config.yaml",
              help = "Path to configuration file [default: %default]"),
  make_option(c("--cores"), type = "integer", default = NULL,
              help = "Number of cores for parallel processing [default: auto-detect]"),
  make_option(c("--skip-reference"), action = "store_true", default = FALSE,
              help = "Skip reference data download"),
  make_option(c("--skip-merge"), action = "store_true", default = FALSE,
              help = "Skip dataset merging (use existing merged file)"),
  make_option(c("--skip-validation"), action = "store_true", default = FALSE,
              help = "Skip reference validation step"),
  make_option(c("--skip-trajectory"), action = "store_true", default = FALSE,
              help = "Skip trajectory analysis"),
  make_option(c("--skip-figures"), action = "store_true", default = FALSE,
              help = "Skip figure generation"),
  make_option(c("--integrate"), action = "store_true", default = FALSE,
              help = "Use Seurat integration for batch correction when merging"),
  make_option(c("--reference-only"), action = "store_true", default = FALSE,
              help = "Only download reference data, then stop"),
  make_option(c("--dry-run"), action = "store_true", default = FALSE,
              help = "Print steps without executing"),
  make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
              help = "Print detailed progress messages")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Run the complete ACP Epithelial Differentiation Analysis pipeline")
opt <- parse_args(opt_parser)

# --- Helper Functions ---------------------------------------------------------

log_step <- function(step_num, step_name, status = "START") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  if (status == "START") {
    cat("\n")
    cat(strrep("=", 70), "\n")
    cat(sprintf("[%s] STEP %s: %s\n", timestamp, step_num, step_name))
    cat(strrep("=", 70), "\n")
  } else if (status == "DONE") {
    cat(sprintf("[%s] STEP %s: %s - COMPLETED\n", timestamp, step_num, step_name))
  } else if (status == "SKIP")
    cat(sprintf("[%s] STEP %s: %s - SKIPPED\n", timestamp, step_num, step_name))
}

run_script <- function(script_path, args = character(), dry_run = FALSE) {
  cmd <- paste("Rscript", script_path, paste(args, collapse = " "))

  if (dry_run) {
    cat("  [DRY RUN] Would execute:", cmd, "\n")
    return(0)
  }

  cat("  Executing:", cmd, "\n\n")
  result <- system(cmd)

  if (result != 0) {
    stop(sprintf("Script failed with exit code %d: %s", result, script_path))
  }

  return(result)
}

check_file_exists <- function(filepath, description) {
  if (!file.exists(filepath)) {
    warning(sprintf("%s not found: %s", description, filepath))
    return(FALSE)
  }
  return(TRUE)
}

# --- Print Configuration ------------------------------------------------------

cat("\n")
cat(strrep("#", 70), "\n")
cat("#  ACP Epithelial Differentiation Analysis - Full Pipeline\n")
cat(strrep("#", 70), "\n")
cat("\nConfiguration:\n")
cat("  Dataset:          ", opt$dataset, "\n")
cat("  Config file:      ", opt$config, "\n")
cat("  Cores:            ", ifelse(is.null(opt$cores), "auto-detect", opt$cores), "\n")
cat("  Skip reference:   ", opt$`skip-reference`, "\n
")
cat("  Skip merge:       ", opt$`skip-merge`, "\n")
cat("  Skip validation:  ", opt$`skip-validation`, "\n")
cat("  Skip trajectory:  ", opt$`skip-trajectory`, "\n")
cat("  Skip figures:     ", opt$`skip-figures`, "\n")
cat("  Integration:      ", opt$integrate, "\n")
cat("  Dry run:          ", opt$`dry-run`, "\n")
cat("\n")

# --- Determine Script Directory -----------------------------------------------

# Find scripts relative to this file or use default
script_dir <- "scripts"
if (!dir.exists(script_dir)) {
  # Try relative to current working directory
  if (dir.exists("../scripts")) {
    script_dir <- "../scripts"
  } else {
    stop("Cannot find scripts directory. Run from project root.")
  }
}

# --- Step 0a: Download Reference Data -----------------------------------------

if (!opt$`skip-reference`) {
  log_step("0a", "Download Reference Data")

  ref_script <- file.path(script_dir, "00a_download_skin_reference.R")
  ref_output <- "data/reference/skin_reference_census.rds"

  if (file.exists(ref_output) && !opt$`dry-run`) {
    cat("  Reference data already exists:", ref_output, "\n")
    cat("  Skipping download. Delete file to re-download.\n")
    log_step("0a", "Download Reference Data", "SKIP")
  } else if (check_file_exists(ref_script, "Reference download script")) {
    ref_args <- c(
      "--cell_type", "keratinocyte",
      "--max_cells", "50000",
      "--assay", "10x",
      "--output_dir", "data/reference"
    )
    run_script(ref_script, ref_args, opt$`dry-run`)
    log_step("0a", "Download Reference Data", "DONE")
  }

  if (opt$`reference-only`) {
    cat("\n--reference-only flag set. Stopping after reference download.\n")
    quit(status = 0)
  }
} else {
  log_step("0a", "Download Reference Data", "SKIP")
}

# --- Step 0: Merge Datasets ---------------------------------------------------

if (opt$dataset == "merged" && !opt$`skip-merge`) {
  log_step("0", "Merge Datasets")

  merge_script <- file.path(script_dir, "00_merge_datasets.R")
  merge_output <- "data/processed/merged_acp_gse.rds"

  if (file.exists(merge_output) && !opt$`dry-run`) {
    cat("  Merged dataset already exists:", merge_output, "\n")
    cat("  Skipping merge. Delete file to re-merge.\n")
    log_step("0", "Merge Datasets", "SKIP")
  } else if (check_file_exists(merge_script, "Merge script")) {
    merge_args <- c("--config", opt$config)
    if (opt$integrate) {
      merge_args <- c(merge_args, "--integrate")
    }
    run_script(merge_script, merge_args, opt$`dry-run`)
    log_step("0", "Merge Datasets", "DONE")
  }
} else if (opt$dataset == "merged") {
  log_step("0", "Merge Datasets", "SKIP")
}

# --- Step 1: Cell Type Annotation ---------------------------------------------

log_step("1", "Cell Type Annotation")

annotation_script <- file.path(script_dir, "01_cell_type_annotation.R")
# Try v2 first, fall back to original
if (!file.exists(annotation_script)) {
  annotation_script <- file.path(script_dir, "01_cell_type_annotation_v2.R")
}

if (check_file_exists(annotation_script, "Annotation script")) {
  annotation_args <- c(
    "--dataset", opt$dataset,
    "--config", opt$config
  )
  run_script(annotation_script, annotation_args, opt$`dry-run`)
  log_step("1", "Cell Type Annotation", "DONE")
}

# --- Step 1b: Reference Validation --------------------------------------------

if (!opt$`skip-validation`) {
  log_step("1b", "Reference Validation & Enhanced Statistics")

  validation_script <- file.path(script_dir, "01b_reference_validation_enhanced_stats.R")
  ref_file <- "data/reference/skin_reference_census.rds"

  if (check_file_exists(validation_script, "Validation script")) {
    validation_args <- c(
      "--dataset", opt$dataset,
      "--config", opt$config
    )

    # Add reference path if available
    if (file.exists(ref_file)) {
      validation_args <- c(validation_args,
                           "--reference", "custom",
                           "--ref_path", ref_file)
    } else {
      validation_args <- c(validation_args, "--reference", "signatures")
      cat("  Note: Using signature-based validation (reference file not found)\n")
    }

    run_script(validation_script, validation_args, opt$`dry-run`)
    log_step("1b", "Reference Validation & Enhanced Statistics", "DONE")
  }
} else {
  log_step("1b", "Reference Validation & Enhanced Statistics", "SKIP")
}

# --- Step 2: Trajectory Analysis ----------------------------------------------

if (!opt$`skip-trajectory`) {
  log_step("2", "Trajectory Analysis")

  trajectory_script <- file.path(script_dir, "02_trajectory_analysis.R")

  if (check_file_exists(trajectory_script, "Trajectory script")) {
    trajectory_args <- c(
      "--dataset", opt$dataset,
      "--config", opt$config
    )
    if (!is.null(opt$cores)) {
      trajectory_args <- c(trajectory_args, "--cores", opt$cores)
    }
    run_script(trajectory_script, trajectory_args, opt$`dry-run`)
    log_step("2", "Trajectory Analysis", "DONE")
  }
} else {
  log_step("2", "Trajectory Analysis", "SKIP")
}

# --- Step 3: Figure Generation ------------------------------------------------

if (!opt$`skip-figures`) {
  log_step("3", "Figure Generation")

  figures_script <- file.path(script_dir, "generate_figures.R")

  if (check_file_exists(figures_script, "Figure generation script")) {
    figures_args <- c("--config", opt$config)
    run_script(figures_script, figures_args, opt$`dry-run`)
    log_step("3", "Figure Generation", "DONE")
  } else {
    cat("  Figure generation script not found. Skipping.\n")
    log_step("3", "Figure Generation", "SKIP")
  }
} else {
  log_step("3", "Figure Generation", "SKIP")
}

# --- Summary ------------------------------------------------------------------

cat("\n")
cat(strrep("=", 70), "\n")
cat("PIPELINE COMPLETE\n")
cat(strrep("=", 70), "\n")
cat("\nOutput locations:\n")
cat("  Annotated cells:    results/objects/01_seurat_annotated_", opt$dataset, ".rds\n", sep = "")
cat("  Validated cells:    results/objects/01b_seurat_validated_", opt$dataset, ".rds\n", sep = "")
cat("  Trajectory:         results/objects/02_seurat_with_pseudotime_", opt$dataset, ".rds\n", sep = "")
cat("  Tables:             results/tables/\n")
cat("  Figures:            results/figures/main/\n")
cat("                      results/figures/supplementary/\n")
cat("\nReference data:       data/reference/\n")
cat("\n")

# Record completion time
cat("Completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
