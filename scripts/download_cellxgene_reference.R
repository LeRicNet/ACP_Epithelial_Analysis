#!/usr/bin/env Rscript
#' =============================================================================
#' CELLxGENE Reference Data: Download and Convert
#' =============================================================================
#'
#' Downloads skin keratinocyte reference data from CELLxGENE Census and converts
#' to Seurat format. Handles Python environment setup automatically.
#'
#' Usage (command line):
#'   Rscript scripts/convert_h5ad_to_seurat.R [--force]
#'
#' Usage (interactive):
#'   source("scripts/convert_h5ad_to_seurat.R")
#'   skin_ref <- get_cellxgene_reference()  # Downloads if needed, returns Seurat
#'
#' Configuration:
#'   Reads paths from config/config.yaml:
#'     paths:
#'       cellxgene_keratinocytes: "data/reference/cellxgene_skin_keratinocytes.h5ad"
#' =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(reticulate)
  library(Matrix)
  library(yaml)
})

# =============================================================================
# Python Environment Setup
# =============================================================================

#' Setup Python environment for scanpy
#'
#' Creates virtualenv 'acp_epi' if needed and installs required packages.
#' Call this once before using convert_h5ad_to_seurat().
#'
#' @param envname Name of virtualenv (default: "acp_epi")
#' @param force_install Reinstall packages even if env exists (default: FALSE)
#' @return Invisible TRUE on success
setup_python_env <- function(envname = "acp_epi", force_install = FALSE) {

  # Clear any conflicting environment variable
  Sys.unsetenv("RETICULATE_PYTHON")

  # Check if virtualenv exists
  existing_envs <- tryCatch(
    reticulate::virtualenv_list(),
    error = function(e) character(0)
  )
  env_exists <- envname %in% existing_envs

  if (!env_exists) {
    message("Creating Python virtual environment: ", envname)
    reticulate::virtualenv_create(envname)
    force_install <- TRUE
  }

  # Activate the environment
  reticulate::use_virtualenv(envname, required = TRUE)
  message("Using Python environment: ", envname)

  # Install packages if needed
  if (force_install) {
    message("Installing Python packages (this may take a few minutes)...")
    packages <- c(
      "numpy<2.4",
      "cellxgene-census",
      "scanpy",
      "anndata",
      "scipy",
      "pandas",
      "h5py"
    )
    reticulate::py_install(packages, envname = envname, pip = TRUE)
    message("Python packages installed successfully")
  }

  invisible(TRUE)
}

#' Ensure Python environment is ready
#'
#' Internal function to verify scanpy is importable
#' @return TRUE if ready, stops with error if not
ensure_python_ready <- function(envname = "acp_epi") {

  # Clear conflicting env var
  Sys.unsetenv("RETICULATE_PYTHON")

  # Try to use the environment
  tryCatch({
    reticulate::use_virtualenv(envname, required = TRUE)
  }, error = function(e) {
    stop("Python environment '", envname, "' not found. ",
         "Run setup_python_env() first.", call. = FALSE)
  })

  # Test scanpy import
  tryCatch({
    reticulate::import("scanpy")
    TRUE
  }, error = function(e) {
    stop("scanpy not available. Run setup_python_env(force_install = TRUE)",
         call. = FALSE)
  })
}

# =============================================================================
# Conversion Function
# =============================================================================

#' Convert H5AD file to Seurat object
#'
#' @param h5ad_path Path to h5ad file
#' @param assay_name Name for the assay in Seurat object (default: "RNA")
#' @param save_rds Save as RDS file? (default: FALSE)
#' @param output_path Path for RDS output (default: input path with .rds extension)
#' @param envname Python virtualenv name (default: "acp_epi")
#' @return Seurat object
#'
#' @examples
#' \dontrun{
#' # First time setup
#' setup_python_env()
#'
#' # Convert file
#' skin_ref <- convert_h5ad_to_seurat("data/reference/cellxgene_skin_keratinocytes.h5ad")
#'
#' # Convert and save as RDS
#' skin_ref <- convert_h5ad_to_seurat(
#'   "data/reference/cellxgene_skin_keratinocytes.h5ad",
#'   save_rds = TRUE
#' )
#' }
convert_h5ad_to_seurat <- function(h5ad_path,
                                   assay_name = "RNA",
                                   save_rds = FALSE,
                                   output_path = NULL,
                                   envname = "acp_epi") {

  # Validate input
  if (!file.exists(h5ad_path)) {
    stop("File not found: ", h5ad_path, call. = FALSE)
  }

  if (!grepl("\\.h5ad$", h5ad_path, ignore.case = TRUE)) {
    warning("File does not have .h5ad extension: ", h5ad_path)
  }

  # Ensure Python is ready
  ensure_python_ready(envname)

  # Import scanpy
  message("Loading h5ad file: ", h5ad_path)
  sc <- reticulate::import("scanpy")

  # Read h5ad
  adata <- sc$read_h5ad(h5ad_path)
  message("Loaded AnnData: ", adata$n_obs, " cells x ", adata$n_vars, " genes")

  # Extract counts matrix (transpose to genes x cells for Seurat)
  message("Extracting expression matrix...")
  counts <- t(adata$X)
  rownames(counts) <- adata$var_names$tolist()
  colnames(counts) <- adata$obs_names$tolist()

  # Convert to sparse matrix if not already
  if (!inherits(counts, "sparseMatrix")) {
    message("Converting to sparse matrix...")
    counts <- as(counts, "CsparseMatrix")
  }

  # Extract metadata
  message("Extracting metadata...")
  meta <- as.data.frame(adata$obs)

  # Create Seurat object
  message("Creating Seurat object...")
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    assay = assay_name
  )

  message("Created Seurat object: ", ncol(seurat_obj), " cells x ",
          nrow(seurat_obj), " genes")

  # Save if requested
  if (save_rds) {
    if (is.null(output_path)) {
      output_path <- sub("\\.h5ad$", ".rds", h5ad_path, ignore.case = TRUE)
    }
    message("Saving to: ", output_path)
    saveRDS(seurat_obj, output_path)
    message("Saved successfully")
  }

  return(seurat_obj)
}

# =============================================================================
# Configuration
# =============================================================================

#' Get configuration value
#' @param key Dot-separated path (e.g., "paths.cellxgene_keratinocytes")
#' @param config_path Path to config file
#' @return Configuration value or NULL
get_config <- function(key, config_path = "config/config.yaml") {
  if (!file.exists(config_path)) {
    return(NULL)
  }


  config <- yaml::read_yaml(config_path)
  keys <- strsplit(key, "\\.")[[1]]

  value <- config
  for (k in keys) {
    if (is.null(value[[k]])) return(NULL)
    value <- value[[k]]
  }
  return(value)
}

#' Default paths
DEFAULT_H5AD_PATH <- "data/reference/cellxgene_skin_keratinocytes.h5ad"
DEFAULT_RDS_PATH <- "data/reference/cellxgene_skin_keratinocytes.rds"
DEFAULT_DOWNLOAD_SCRIPT <- "../python/download_cellxgene_reference.py"

# =============================================================================
# Download Function
# =============================================================================

#' Download CELLxGENE skin keratinocyte reference data
#'
#' Runs the Python download script if h5ad file doesn't exist or force=TRUE.
#'
#' @param h5ad_path Output path for h5ad file (default: from config or default)
#' @param max_cells Maximum cells to download (default: 150000)
#' @param force Re-download even if file exists (default: FALSE)
#' @param script_path Path to Python download script
#' @return Path to downloaded h5ad file (invisibly)
#'
#' @examples
#' \dontrun
#' # Download if not present
#' download_cellxgene_reference()
#'
#' # Force re-download
#' download_cellxgene_reference(force = TRUE)
#' }
download_cellxgene_reference <- function(h5ad_path = NULL,
                                         max_cells = 150000,
                                         force = FALSE,
                                         script_path = DEFAULT_DOWNLOAD_SCRIPT) {


  # Get path from config if not specified
  if (is.null(h5ad_path)) {
    h5ad_path <- get_config("paths.cellxgene_keratinocytes")
    if (is.null(h5ad_path)) {
      h5ad_path <- DEFAULT_H5AD_PATH
    }
  }

  # Check if download needed
  if (file.exists(h5ad_path) && !force) {
    message("Reference file already exists: ", h5ad_path)
    message("Use force=TRUE to re-download")
    return(invisible(h5ad_path))
  }

  # Verify Python script exists
  if (!file.exists(script_path)) {
    stop("Download script not found: ", script_path, "\n",
         "Expected at: ", normalizePath(script_path, mustWork = FALSE),
         call. = FALSE)
  }

  # Create output directory
  dir.create(dirname(h5ad_path), recursive = TRUE, showWarnings = FALSE)

  # Build command
  cmd <- sprintf(
    "python %s --tissue skin --max_cells %d",
    shQuote(script_path),
    max_cells
  )

  message("Downloading CELLxGENE skin keratinocyte reference...")
  message("Command: ", cmd)
  message("This may take several minutes...")

  # Run download
  exit_code <- system(cmd)

  if (exit_code != 0) {
    stop("Download failed with exit code: ", exit_code, call. = FALSE)
  }

  # Verify file was created

  if (!file.exists(h5ad_path)) {
    stop("Download completed but file not found at: ", h5ad_path, call. = FALSE)
  }

  file_size <- file.info(h5ad_path)$size / 1e6
  message("Download complete: ", h5ad_path, " (", round(file_size, 1), " MB)")

  return(invisible(h5ad_path))
}

# =============================================================================
# Main Reference Function
# =============================================================================

#' Get CELLxGENE skin keratinocyte reference as Seurat object
#'
#' Downloads reference data if needed, converts to Seurat, and optionally caches as RDS.
#' This is the main function most users should call.
#'
#' @param force_download Re-download h5ad even if exists (default: FALSE)
#' @param force_convert Re-convert to Seurat even if RDS exists (default: FALSE
#' @param max_cells Maximum cells to download (default: 150000)
#' @param save_rds Cache the Seurat object as RDS (default: TRUE)
#' @param h5ad_path Path to h5ad file (default: from config)
#' @param rds_path Path to RDS cache (default: from config)
#' @return Seurat object with skin keratinocyte reference
#'
#' @examples
#' \dontrun{
#' # First time: downloads and converts
#' skin_ref <- get_cellxgene_reference()
#'
#' # Subsequent: loads from RDS cache
#' skin_ref <- get_cellxgene_reference()
#'
#' # Force fresh download
#' skin_ref <- get_cellxgene_reference(force_download = TRUE)
#' }
get_cellxgene_reference <- function(force_download = FALSE,
                                    force_convert = FALSE,
                                    max_cells = 150000,
                                    save_rds = TRUE,
                                    h5ad_path = NULL,
                                    rds_path = NULL) {

  # Get paths from config if not specified
  if (is.null(h5ad_path)) {
    h5ad_path <- get_config("paths.cellxgene_keratinocytes")
    if (is.null(h5ad_path)) h5ad_path <- DEFAULT_H5AD_PATH
  }

  if (is.null(rds_path)) {
    rds_path <- get_config("paths.cellxgene_keratinocytes_rds")
    if (is.null(rds_path)) {
      rds_path <- sub("\\.h5ad$", ".rds", h5ad_path, ignore.case = TRUE)
    }
  }

  # Check for cached RDS first (fastest path)
  if (file.exists(rds_path) && !force_convert && !force_download) {
    message("Loading cached reference: ", rds_path)
    return(readRDS(rds_path))
  }

  # Download h5ad if needed
  if (!file.exists(h5ad_path) || force_download) {
    download_cellxgene_reference(
      h5ad_path = h5ad_path,
      max_cells = max_cells,
      force = force_download
    )
  }

  # Convert to Seurat
  seurat_obj <- convert_h5ad_to_seurat(
    h5ad_path = h5ad_path,
    save_rds = save_rds,
    output_path = rds_path
  )

  return(seurat_obj)
}

# =============================================================================
# Command Line Interface
# =============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  # Parse arguments
  force_download <- "--force" %in% args || "-f" %in% args
  help_requested <- "--help" %in% args || "-h" %in% args

  if (help_requested) {
    cat("Usage: Rscript convert_h5ad_to_seurat.R [OPTIONS]\n")
    cat("\nDownloads CELLxGENE skin keratinocyte reference and converts to Seurat.\n")
    cat("\nOptions:\n")
    cat("  --force, -f    Re-download and re-convert even if files exist\n")
    cat("  --help, -h     Show this help message\n")
    cat("\nPaths are read from config/config.yaml:\n")
    cat("  paths:\n")
    cat("    cellxgene_keratinocytes: data/reference/cellxgene_skin_keratinocytes.h5ad\n")
    quit(status = 0)
  }

  # Setup Python environment
  setup_python_env()

  # Get reference (download + convert as needed)
  skin_ref <- get_cellxgene_reference(
    force_download = force_download,
    force_convert = force_download
  )

  # Print summary
  cat("\n")
  cat("Reference loaded successfully!\n")
  cat("Cells: ", ncol(skin_ref), "\n")
  cat("Genes: ", nrow(skin_ref), "\n")
  cat("\nCell type distribution:\n")
  print(table(skin_ref$cell_type))
  cat("\nDone!\n")
}
