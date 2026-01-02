#!/usr/bin/env Rscript
#' =============================================================================
#' CELLxGENE Reference Data: Download and Convert
#' =============================================================================
#'
#' Downloads reference data from CELLxGENE Census (via query) or directly from
#' a URL, and converts to Seurat format.
#'
#' Usage (command line):
#'   Rscript scripts/download_cellxgene_reference.R [--force] [--dataset NAME]
#'
#' Usage (interactive):
#'   source("scripts/download_cellxgene_reference.R")
#'   skin_ref <- get_cellxgene_reference()  # Query-based download
#'   oral_ref <- get_reference_from_url()   # Direct URL download
#'
#' Configuration (config/config.yaml):
#'   external:
#'     cellxgene_oral_atlas: "https://datasets.cellxgene.cziscience.com/..."
#'   paths:
#'     cellxgene_oral_atlas: "data/reference/cellxgene_oral_atlas.rds"
#' =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(reticulate)
  library(Matrix)
  library(yaml)
})

cat("Working directory:", getwd(), "\n")

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

#' Get Python executable path from virtualenv
#'
#' @param envname Name of virtualenv
#' @return Path to Python executable
get_virtualenv_python <- function(envname = "acp_epi") {
  # Get the virtualenv root directory
  venv_root <- reticulate::virtualenv_root()
  venv_path <- file.path(venv_root, envname)

  # Python executable location differs by OS
  if (.Platform$OS.type == "windows") {
    python_path <- file.path(venv_path, "Scripts", "python.exe")
  } else {
    python_path <- file.path(venv_path, "bin", "python")
  }

  if (!file.exists(python_path)) {
    stop("Python executable not found at: ", python_path, call. = FALSE)
  }

  return(python_path)
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
DEFAULT_DOWNLOAD_SCRIPT <- "python/download_cellxgene_reference.py"

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

  # Get gene names from feature_name column
  message("Extracting gene names...")
  gene_names <- adata$var$feature_name

  # Handle NULL or missing feature_name
  if (is.null(gene_names) || length(gene_names) == 0) {
    message("  feature_name not found, trying var_names...")
    gene_names <- adata$var_names$tolist()
  }

  # Convert to R character vector if needed
  if (!is.character(gene_names)) {
    gene_names <- as.character(gene_names)
  }

  message(sprintf("  Gene names extracted: %d", length(gene_names)))
  message(sprintf("  Examples: %s", paste(head(gene_names, 5), collapse = ", ")))

  # Check for and handle duplicates
  n_dups <- sum(duplicated(gene_names))
  if (n_dups > 0) {
    message(sprintf("  Found %d duplicate gene names, making unique...", n_dups))
    gene_names <- make.unique(gene_names, sep = "_dup")
  }

  # Assign gene names to matrix
  rownames(counts) <- gene_names
  colnames(counts) <- adata$obs_names$tolist()

  # Convert to sparse matrix if not already
  if (!inherits(counts, "sparseMatrix")) {
    message("Converting to sparse matrix...")
    counts <- as(counts, "CsparseMatrix")
  } else if (inherits(counts, "dgRMatrix")) {
    # Convert row-sparse to column-sparse (required by Seurat)
    message("Converting dgRMatrix to dgCMatrix...")
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

  # Normalize data (required for AddModuleScore in Seurat v5)
  message("\nNormalizing data...")
  message("  (Required for module scoring in Seurat v5)")
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

  # Verify layers
  available_layers <- Layers(seurat_obj[[assay_name]])
  message(sprintf("  Available layers: %s", paste(available_layers, collapse = ", ")))

  # Verify gene names
  message("\nVerifying gene names...")
  test_genes <- c("KRT5", "KRT14", "KRT1", "KRT10", "KRT4", "KRT13", "TP63", "GAPDH", "ACTB")
  present <- test_genes[test_genes %in% rownames(seurat_obj)]
  message(sprintf("  Key markers present: %d/%d", length(present), length(test_genes)))
  message(sprintf("  Found: %s", paste(present, collapse = ", ")))

  # Save if requested
  if (save_rds) {
    if (is.null(output_path)) {
      output_path <- sub("\\.h5ad$", ".rds", h5ad_path, ignore.case = TRUE)
    }
    message("Saving to: ", output_path)
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(seurat_obj, output_path)
    message("Saved successfully")
  }

  return(seurat_obj)
}

# =============================================================================
# DIRECT URL DOWNLOAD (NEW FUNCTIONALITY)
# =============================================================================

#' Download h5ad file directly from URL
#'
#' Downloads an h5ad file from a direct URL (e.g., CELLxGENE permanent link).
#' This is useful when you have a specific dataset URL rather than querying
#' the Census API.
#'
#' @param url URL to the h5ad file
#' @param h5ad_path Local path to save the h5ad file
#' @param force Re-download even if file exists (default: FALSE)
#' @param timeout Download timeout in seconds (default: 3600 = 1 hour)
#' @return Path to downloaded h5ad file (invisibly)
#'
#' @examples
#' \dontrun{
#' download_h5ad_from_url(
#'   url = "https://datasets.cellxgene.cziscience.com/xyz.h5ad",
#'   h5ad_path = "data/reference/my_reference.h5ad"
#' )
#' }
download_h5ad_from_url <- function(url,
                                   h5ad_path,
                                   force = FALSE,
                                   timeout = 3600) {

  # Check if download needed
  if (file.exists(h5ad_path) && !force) {
    message("File already exists: ", h5ad_path)
    message("Use force=TRUE to re-download")
    return(invisible(h5ad_path))
  }

  # Validate URL

  if (!grepl("^https?://", url)) {
    stop("Invalid URL: ", url, call. = FALSE)
  }

  # Create output directory
  dir.create(dirname(h5ad_path), recursive = TRUE, showWarnings = FALSE)

  message("=" , strrep("=", 69))
  message("Downloading h5ad from URL")
  message(strrep("=", 70))
  message("\nSource: ", url)
  message("Destination: ", h5ad_path)
  message("\nThis may take several minutes for large files...")

  # Set timeout for large file downloads
  old_timeout <- getOption("timeout")
  options(timeout = timeout)

  # Download file
  tryCatch({
    download.file(
      url = url,
      destfile = h5ad_path,
      mode = "wb",  # Binary mode for h5ad
      quiet = FALSE
    )
  }, error = function(e) {
    # Clean up partial download
    if (file.exists(h5ad_path)) {
      file.remove(h5ad_path)
    }
    stop("Download failed: ", e$message, call. = FALSE)
  }, finally = {
    # Restore original timeout
    options(timeout = old_timeout)
  })

  # Verify download
  if (!file.exists(h5ad_path)) {
    stop("Download completed but file not found!", call. = FALSE)
  }

  file_size <- file.info(h5ad_path)$size / 1e6
  message(sprintf("\nDownload complete: %s (%.1f MB)", h5ad_path, file_size))

  return(invisible(h5ad_path))
}


#' Get reference from URL (download + convert)
#'
#' Downloads an h5ad file from a URL and converts to Seurat object.
#' Reads URL and output paths from config if not specified.
#'
#' @param url URL to h5ad file (default: from config external.cellxgene_oral_atlas)
#' @param h5ad_path Local path for h5ad (default: derived from rds_path)
#' @param rds_path Path for RDS cache (default: from config paths.cellxgene_oral_atlas)
#' @param force_download Re-download even if h5ad exists (default: FALSE)
#' @param force_convert Re-convert even if RDS exists (default: FALSE)
#' @param save_rds Save converted Seurat as RDS (default: TRUE)
#' @param config_key_url Config key for URL (default: "external.cellxgene_oral_atlas")
#' @param config_key_path Config key for output path (default: "paths.cellxgene_oral_atlas")
#' @param envname Python virtualenv name (default: "acp_epi")
#' @return Seurat object
#'
#' @examples
#' \dontrun
#' # Use config defaults (oral atlas)
#' oral_ref <- get_reference_from_url()
#'
#' # Specify custom URL and path
#' custom_ref <- get_reference_from_url(
#'   url = "https://example.com/data.h5ad",
#'   rds_path = "data/reference/custom.rds"
#' )
#' }
get_reference_from_url <- function(url = NULL,
                                   h5ad_path = NULL,
                                   rds_path = NULL,
                                   force_download = FALSE,
                                   force_convert = FALSE,
                                   save_rds = TRUE,
                                   config_key_url = "external.cellxgene_oral_atlas",
                                   config_key_path = "paths.cellxgene_oral_atlas",
                                   envname = "acp_epi") {

  # Get URL from config if not specified
  if (is.null(url)) {
    url <- get_config(config_key_url)
    if (is.null(url)) {
      stop("URL not specified and not found in config at: ", config_key_url, call. = FALSE)
    }
  }

  # Get RDS path from config if not specified
  if (is.null(rds_path)) {
    rds_path <- get_config(config_key_path)
    if (is.null(rds_path)) {
      stop("Output path not specified and not found in config at: ", config_key_path, call. = FALSE)
    }
  }

  # Derive h5ad path from rds path if not specified
  if (is.null(h5ad_path)) {
    h5ad_path <- sub("\\.rds$", ".h5ad", rds_path, ignore.case = TRUE)
  }

  message("\n", strrep("=", 70))
  message("Loading Reference from URL")
  message(strrep("=", 70))
  message("URL: ", url)
  message("H5AD: ", h5ad_path)
  message("RDS: ", rds_path)

  # Check for cached RDS first (fastest path)
  if (file.exists(rds_path) && !force_convert && !force_download) {
    message("\nLoading cached reference: ", rds_path)
    return(readRDS(rds_path))
  }

  # Download h5ad if needed
  if (!file.exists(h5ad_path) || force_download) {
    download_h5ad_from_url(
      url = url,
      h5ad_path = h5ad_path,
      force = force_download
    )
  }

  # Setup Python env
  setup_python_env(envname = envname, force_install = FALSE)

  # Convert to Seurat
  seurat_obj <- convert_h5ad_to_seurat(
    h5ad_path = h5ad_path,
    save_rds = save_rds,
    output_path = rds_path,
    envname = envname
  )

  return(seurat_obj)
}


# =============================================================================
# QUERY-BASED DOWNLOAD (EXISTING FUNCTIONALITY)
# =============================================================================

#' Download CELLxGENE skin keratinocyte reference data
#'
#' Runs the Python download script using the reticulate virtualenv.
#'
#' @param h5ad_path Output path for h5ad file (default: from config or default)
#' @param max_cells Maximum cells to download (default: 150000)
#' @param force Re-download even if file exists (default: FALSE)
#' @param script_path Path to Python download script
#' @param envname Name of Python virtualenv to use (default: "acp_epi")
#' @param use_reticulate_env Use reticulate virtualenv Python (default: TRUE)
#' @return Path to downloaded h5ad file (invisibly)
download_cellxgene_reference <- function(h5ad_path = NULL,
                                         max_cells = 150000,
                                         force = FALSE,
                                         script_path = DEFAULT_DOWNLOAD_SCRIPT,
                                         envname = "acp_epi",
                                         use_reticulate_env = TRUE) {

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

  message("Downloading CELLxGENE skin keratinocyte reference...")
  message("This may take several minutes...")

  if (use_reticulate_env) {
    # Use reticulate virtualenv Python executable
    setup_python_env(envname = envname, force_install = FALSE)
    python_path <- get_virtualenv_python(envname)
    message("Using Python from virtualenv: ", python_path)

    cmd <- sprintf(
      "%s %s --tissue skin --max_cells %d",
      shQuote(python_path),
      shQuote(normalizePath(script_path)),
      max_cells
    )

    message("Command: ", cmd)
    exit_code <- system(cmd)

  } else {
    # Use reticulate to run Python code directly
    setup_python_env(envname = envname, force_install = FALSE)
    ensure_python_ready(envname)

    message("Running download via reticulate...")
    sys <- reticulate::import("sys")
    original_argv <- sys$argv
    sys$argv <- list(
      script_path,
      "--tissue", "skin",
      "--max_cells", as.character(max_cells)
    )

    exit_code <- tryCatch({
      reticulate::source_python(script_path)
      0L
    }, error = function(e) {
      message("Error running Python script: ", e$message)
      1L
    }, finally = {
      sys$argv <- original_argv
    })
  }

  if (exit_code != 0) {
    stop("Download failed with exit code: ", exit_code, call. = FALSE)
  }

  if (!file.exists(h5ad_path)) {
    stop("Download completed but file not found at: ", h5ad_path, call. = FALSE)
  }

  file_size <- file.info(h5ad_path)$size / 1e6
  message("Download complete: ", h5ad_path, " (", round(file_size, 1), " MB)")

  return(invisible(h5ad_path))
}


#' Download CELLxGENE reference using pure reticulate (no external script)
#'
#' Alternative download method that runs Python code directly through reticulate
#' without needing the external Python script file.
#'
#' @param h5ad_path Output path for h5ad file
#' @param max_cells Maximum cells to download
#' @param envname Name of Python virtualenv
#' @param force Re-download even if exists
#' @return Path to h5ad file (invisibly)
download_cellxgene_via_reticulate <- function(h5ad_path = DEFAULT_H5AD_PATH,
                                              max_cells = 150000,
                                              envname = "acp_epi",
                                              force = FALSE) {

  # Check if download needed
  if (file.exists(h5ad_path) && !force) {
    message("Reference file already exists: ", h5ad_path)
    return(invisible(h5ad_path))
  }

  # Setup and activate environment
  setup_python_env(envname = envname, force_install = FALSE)
  ensure_python_ready(envname)

  message("Downloading CELLxGENE skin keratinocyte reference via reticulate...")
  message("Max cells: ", max_cells)
  message("Output: ", h5ad_path)
  message("This may take several minutes...")

  # Create output directory
  dir.create(dirname(h5ad_path), recursive = TRUE, showWarnings = FALSE)

  # Import required Python modules
  cellxgene_census <- reticulate::import("cellxgene_census")
  np <- reticulate::import("numpy")

  message("\nOpening CELLxGENE Census...")

  # Use context manager equivalent
  census <- cellxgene_census$open_soma()

  tryCatch({
    # Step 1: Explore available cell types in skin tissue
    message("Exploring available skin cell types...")

    skin_filter <- "tissue_general == 'skin of body' and is_primary_data == True and disease == 'normal'"

    obs_explore <- cellxgene_census$get_obs(
      census,
      organism = "Homo sapiens",
      value_filter = skin_filter,
      column_names = list("cell_type", "tissue", "disease")
    )

    message(sprintf("Total skin cells in Census: %d", nrow(obs_explore)))

    if (nrow(obs_explore) > 0) {
      cell_type_counts <- table(obs_explore$cell_type)
      cell_type_counts <- sort(cell_type_counts, decreasing = TRUE)

      message("\nTop 20 cell types in skin:")
      top_types <- head(cell_type_counts, 20)
      for (i in seq_along(top_types)) {
        message(sprintf("  %s: %d", names(top_types)[i], top_types[i]))
      }

      kc_types <- names(cell_type_counts)[grepl("keratinocyte|basal|epiderm|epithelial",
                                                names(cell_type_counts),
                                                ignore.case = TRUE)]
      message("\nKeratinocyte-related types found:")
      for (ct in kc_types) {
        message(sprintf("  %s: %d", ct, cell_type_counts[ct]))
      }
    }

    # Step 2: Build the actual query for keratinocytes
    cell_type_filter <- paste0(
      "tissue_general == 'skin of body'",
      " and is_primary_data == True",
      " and disease == 'normal'"
    )

    message("\nQuerying with filter:")
    message(sprintf("  %s", cell_type_filter))

    obs_df <- cellxgene_census$get_obs(
      census,
      organism = "Homo sapiens",
      value_filter = cell_type_filter,
      column_names = list("cell_type", "tissue", "disease", "assay", "donor_id")
    )

    total_cells <- nrow(obs_df)
    message(sprintf("\nFound %d total skin cells", total_cells))

    if (total_cells == 0) {
      message("\nTrying broader query (all skin, any disease status)...")
      broad_filter <- "tissue_general == 'skin of body' and is_primary_data == True"

      obs_df <- cellxgene_census$get_obs(
        census,
        organism = "Homo sapiens",
        value_filter = broad_filter,
        column_names = list("cell_type", "tissue", "disease")
      )

      total_cells <- nrow(obs_df)
      message(sprintf("Broad query found %d cells", total_cells))

      if (total_cells == 0) {
        stop("No skin cells found in CELLxGENE Census. The API schema may have changed.",
             call. = FALSE)
      }

      cell_type_filter <- broad_filter
    }

    # Step 3: Filter for keratinocyte-related cell types
    all_cell_types <- unique(obs_df$cell_type)
    kc_pattern <- "keratinocyte|basal cell of epidermis|spinous|granular|epithelial cell of epidermis"
    kc_cell_types <- all_cell_types[grepl(kc_pattern, all_cell_types, ignore.case = TRUE)]

    if (length(kc_cell_types) == 0) {
      message("No keratinocyte-specific types found. Using all skin cells.")
      kc_cell_types <- all_cell_types
    } else {
      message(sprintf("\nFiltering for %d keratinocyte-related cell types:", length(kc_cell_types)))
      for (ct in kc_cell_types) {
        message(sprintf("  - %s", ct))
      }
    }

    kc_types_escaped <- gsub("'", "\\'", kc_cell_types, fixed = TRUE)
    kc_types_str <- paste0("'", kc_types_escaped, "'", collapse = ", ")

    final_filter <- sprintf(
      "tissue_general == 'skin of body' and is_primary_data == True and disease == 'normal' and cell_type in [%s]",
      kc_types_str
    )

    message("\nFinal query filter:")
    message(sprintf("  %s", substr(final_filter, 1, 200), "..."))

    obs_df <- cellxgene_census$get_obs(
      census,
      organism = "Homo sapiens",
      value_filter = final_filter,
      column_names = list("cell_type", "tissue", "disease")
    )

    total_cells <- nrow(obs_df)
    message(sprintf("\nKeratinocyte cells to download: %d", total_cells))

    if (total_cells == 0) {
      stop("No keratinocyte cells found after filtering!", call. = FALSE)
    }

    # Step 4: Subsample if needed
    obs_coords <- NULL
    if (total_cells > max_cells) {
      message(sprintf("Subsampling from %d to %d cells...", total_cells, max_cells))
      set.seed(42)
      sample_idx <- sample(total_cells, max_cells, replace = FALSE) - 1L
      obs_coords <- as.integer(sample_idx)
    }

    # Step 5: Download the expression data
    message("\nDownloading gene expression data...")
    message("(This may take several minutes...)")

    adata <- cellxgene_census$get_anndata(
      census,
      organism = "Homo sapiens",
      obs_value_filter = final_filter,
      obs_column_names = list(
        "cell_type",
        "cell_type_ontology_term_id",
        "tissue",
        "tissue_general",
        "disease",
        "assay",
        "donor_id",
        "dataset_id",
        "sex",
        "development_stage"
      ),
      obs_coords = obs_coords
    )

    n_cells <- as.integer(adata$n_obs)
    n_genes <- as.integer(adata$n_vars)

    message(sprintf("Downloaded: %d cells x %d genes", n_cells, n_genes))

    if (n_cells == 0) {
      stop("Downloaded 0 cells! Something went wrong with the query.", call. = FALSE)
    }

    cell_types_downloaded <- table(reticulate::py_to_r(adata$obs$cell_type))
    message("\nCell types downloaded:")
    for (ct in names(cell_types_downloaded)) {
      message(sprintf("  %s: %d", ct, cell_types_downloaded[ct]))
    }

    # Step 6: Save to h5ad
    message(sprintf("\nSaving to: %s", h5ad_path))
    adata$write_h5ad(h5ad_path)

  }, finally = {
    census$close()
  })

  if (!file.exists(h5ad_path)) {
    stop("Download completed but file not found!", call. = FALSE)
  }

  file_size <- file.info(h5ad_path)$size / 1e6
  message(sprintf("Download complete: %s (%.1f MB)", h5ad_path, file_size))

  return(invisible(h5ad_path))
}

# =============================================================================
# Main Reference Function (Query-based)
# =============================================================================

#' Get CELLxGENE skin keratinocyte reference as Seurat object
#'
#' Downloads reference data if needed, converts to Seurat, and optionally caches as RDS.
#' This is the main function most users should call for query-based downloads.
#'
#' @param force_download Re-download h5ad even if exists (default: FALSE)
#' @param force_convert Re-convert to Seurat even if RDS exists (default: FALSE)
#' @param max_cells Maximum cells to download (default: 150000)
#' @param save_rds Cache the Seurat object as RDS (default: TRUE)
#' @param h5ad_path Path to h5ad file (default: from config)
#' @param rds_path Path to RDS cache (default: from config)
#' @param use_reticulate_download Use pure reticulate download (default: TRUE)
#' @param envname Python virtualenv name (default: "acp_epi")
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
                                    rds_path = NULL,
                                    use_reticulate_download = TRUE,
                                    envname = "acp_epi") {

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
    if (use_reticulate_download) {
      download_cellxgene_via_reticulate(
        h5ad_path = h5ad_path,
        max_cells = max_cells,
        envname = envname,
        force = force_download
      )
    } else {
      download_cellxgene_reference(
        h5ad_path = h5ad_path,
        max_cells = max_cells,
        force = force_download,
        envname = envname,
        use_reticulate_env = TRUE
      )
    }
  }

  # Convert to Seurat
  seurat_obj <- convert_h5ad_to_seurat(
    h5ad_path = h5ad_path,
    save_rds = save_rds,
    output_path = rds_path,
    envname = envname
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
  use_script <- "--use-script" %in% args

  # Check for --dataset argument
  dataset_arg <- grep("^--dataset", args, value = TRUE)
  if (length(dataset_arg) > 0) {
    dataset <- sub("^--dataset[= ]?", "", dataset_arg)
  } else {
    dataset <- "skin"  # Default
  }

  if (help_requested) {
    cat("Usage: Rscript download_cellxgene_reference.R [OPTIONS]\n")
    cat("\nDownloads CELLxGENE reference data and converts to Seurat.\n")
    cat("\nOptions:\n")
    cat("  --force, -f           Re-download and re-convert even if files exist\n")
    cat("  --dataset=NAME        Dataset to download: 'skin' (query) or 'oral' (URL)\n")
    cat("  --use-script          Use external Python script (default: use reticulate)\n")
    cat("  --help, -h            Show this help message\n")
    cat("\nDatasets:\n")
    cat("  skin    Query-based download of skin keratinocytes from Census\n")
    cat("  oral    Direct URL download of oral atlas from config\n")
    cat("\nConfig paths (config/config.yaml):\n")
    cat("  external:\n")
    cat("    cellxgene_oral_atlas: <URL>\n")
    cat("  paths:\n")
    cat("    cellxgene_keratinocytes: data/reference/cellxgene_skin_keratinocytes.rds\n")
    cat("    cellxgene_oral_atlas: data/reference/cellxgene_oral_atlas.rds\n")
    quit(status = 0)
  }

  # Setup Python environment
  setup_python_env()

  if (dataset == "oral") {
    # Direct URL download
    cat("\n")
    cat("Downloading oral atlas from URL...\n")
    ref <- get_reference_from_url(
      force_download = force_download,
      force_convert = force_download
    )
  } else {
    # Query-based download (skin)
    cat("\n")
    cat("Downloading skin keratinocytes via Census query...\n")
    ref <- get_cellxgene_reference(
      force_download = force_download,
      force_convert = force_download,
      use_reticulate_download = !use_script
    )
  }

  # Print summary
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("Reference loaded successfully!\n")
  cat(strrep("=", 70), "\n")
  cat("Cells: ", ncol(ref), "\n")
  cat("Genes: ", nrow(ref), "\n")
  cat("\nExample gene names:\n")
  print(head(rownames(ref), 10))
  cat("\nCell type distribution:\n")
  print(head(sort(table(ref$cell_type), decreasing = TRUE), 20))
  cat("\nDone!\n")
}
