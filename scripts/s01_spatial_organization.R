#!/usr/bin/env Rscript
# ==============================================================================
# scripts/s01_spatial_organization.R
# ==============================================================================
# Spatial Organization Analysis for ACP Epithelial Differentiation
#
# This script performs comprehensive spatial analysis to test whether
# epithelial subtypes show spatial organization reflecting developmental
# hierarchy, using both label-free and label-transfer approaches.
#
# Analysis Tracks:
#   LABEL-FREE: CytoTRACE2, spatial expression gradients, spatial autocorrelation
#   LABEL-TRANSFER: GSE subtypes → Spatial via Harmony + KNN
#   COMPARISON: Spatial gradients vs differentiation metrics
#
# Input:
#   - Spatial Visium object (epithelial spots identified via 'cellstates' column)
#   - GSE215932 snRNA-seq object (published subtypes in 'cell_type' column)
#
# Output:
#   - Annotated spatial object with transferred labels + CytoTRACE scores
#   - Comprehensive results list
#   - Figure 1 panels
#   - Summary tables
#
# Usage:
#   Rscript scripts/s01_spatial_organization.R
#   Rscript scripts/s01_spatial_organization.R --config path/to/config.yaml
#   Rscript scripts/s01_spatial_organization.R --cores 16
#
# ==============================================================================

# ==============================================================================
# SECTION 0: SETUP & CONFIGURATION
# ==============================================================================

message("\n")
message("╔══════════════════════════════════════════════════════════════════════╗")
message("║     Section 1: Spatial Organization of ACP Epithelial Subtypes       ║")
message("╚══════════════════════════════════════════════════════════════════════╝")
message("\n")

# --- Parse command line arguments ---------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
n_cores_override <- NULL

i <- 1
while (i <= length(args)) {

  if (args[i] == "--config" && i < length(args)) {
    config_path <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--cores" && i < length(args)) {
    n_cores_override <- as.integer(args[i + 1])
    i <- i + 2
  } else {
    i <- i + 1
  }
}

# --- Find project root --------------------------------------------------------
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
message(sprintf("Project root: %s", project_root))

# --- Load required packages ---------------------------------------------------
message("\n--- Loading packages ---")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(future)
  library(future.apply)
})

# Source utility functions
source("R/utils/config.R")

# --- Load configuration -------------------------------------------------------
if (is.null(config_path)) {
  config_path <- "config/config.yaml"
}
config <- load_config(config_path)
print_config_summary(config)

# --- Setup parallelization ----------------------------------------------------
n_cores <- if (!is.null(n_cores_override)) {
  n_cores_override
} else if (!is.null(config$reproducibility$n_cores)) {
  config$reproducibility$n_cores
} else {
  # Use more cores by default - at least 8 or leave 1 for system

  max(8, parallel::detectCores() - 1)
}

message(sprintf("\n--- Setting up parallel processing: %d cores ---", n_cores))
plan(multisession, workers = n_cores)
options(future.globals.maxSize = 4 * 1024^3)  # 4 GB

# --- Create output directories ------------------------------------------------
tables_dir <- get_path(config, config$paths$tables_dir)
fig_main_dir <- get_path(config, config$paths$figures_main_dir)
fig_supp_dir <- get_path(config, config$paths$figures_supp_dir)
objects_dir <- get_path(config, config$paths$objects_dir)

ensure_dir(tables_dir)
ensure_dir(fig_main_dir)
ensure_dir(fig_supp_dir)
ensure_dir(objects_dir)

# --- Initialize results list --------------------------------------------------
results <- list(
  metadata = list(
    script = "s01_spatial_organization.R",
    timestamp = Sys.time(),
    config_path = config_path,
    n_cores = n_cores
  )
)

# ==============================================================================
# SECTION 0.5: LOAD DATA
# ==============================================================================

message("\n========================================")
message("Loading Data")
message("========================================\n")

# --- Load spatial object ------------------------------------------------------
spatial_path <- get_path(config, config$paths$spatial_object)
message(sprintf("Loading spatial object: %s", spatial_path))

if (!file.exists(spatial_path)) {
  stop("Spatial object not found: ", spatial_path)
}

spatial_obj <- readRDS(spatial_path)
spatial_obj <- UpdateSeuratObject(spatial_obj)

message(sprintf("  Total spots: %d", ncol(spatial_obj)))
message(sprintf("  Genes: %d", nrow(spatial_obj)))
message(sprintf("  Samples: %s", paste(unique(spatial_obj$sample), collapse = ", ")))

# --- Filter to epithelial spots -----------------------------------------------
message("\n--- Filtering to epithelial spots ---")

if (!"cellstates" %in% colnames(spatial_obj@meta.data)) {
  stop("Column 'cellstates' not found in spatial object metadata")
}

epithelial_spots <- spatial_obj$cellstates == "Epithelial"
message(sprintf("  Epithelial spots: %d (%.1f%%)",
                sum(epithelial_spots),
                mean(epithelial_spots) * 100))

# Store full object reference, subset to epithelial
spatial_full <- spatial_obj
spatial_epi <- subset(spatial_obj, cells = colnames(spatial_obj)[epithelial_spots])

message(sprintf("  Spots per sample:"))
print(table(spatial_epi$sample))

results$metadata$n_epithelial_spots <- ncol(spatial_epi)
results$metadata$samples <- unique(spatial_epi$sample)

# --- Load GSE snRNA-seq object ------------------------------------------------
message("\n--- Loading GSE215932 snRNA-seq object ---")

gse_path <- get_path(config, config$paths$snrnaseq_processed)
message(sprintf("Loading: %s", gse_path))

if (!file.exists(gse_path)) {
  stop("GSE object not found: ", gse_path)
}

gse_obj <- readRDS(gse_path)
gse_obj <- UpdateSeuratObject(gse_obj)

message(sprintf("  Total cells: %d", ncol(gse_obj)))
message(sprintf("  Genes: %d", nrow(gse_obj)))

# Check for cell_type column
if (!"cell_type" %in% colnames(gse_obj@meta.data)) {
  stop("Column 'cell_type' not found in GSE object metadata")
}

message("  Published subtypes (cell_type):")
print(table(gse_obj$cell_type))

# Filter GSE to epithelial cells only (those with subtype labels)
# Assuming all labeled cells are epithelial
gse_subtypes <- unique(gse_obj$cell_type)
gse_epi_labels <- gse_subtypes[grepl("^E_", gse_subtypes)]  # E_C1_PE, etc.

if (length(gse_epi_labels) == 0) {
  # Try alternative - assume all are epithelial if no E_ prefix
  gse_epi_labels <- gse_subtypes
  message("  Note: Using all cell_type labels as epithelial")
}

gse_epi <- subset(gse_obj, cells = colnames(gse_obj)[gse_obj$cell_type %in% gse_epi_labels])
message(sprintf("  GSE epithelial cells: %d", ncol(gse_epi)))

results$metadata$gse_subtypes <- gse_epi_labels
results$metadata$n_gse_epithelial <- ncol(gse_epi)


# ==============================================================================
# SECTION 1: LABEL-FREE ANALYSIS - CYTOTRACE2
# ==============================================================================

message("\n========================================")
message("Section 1: CytoTRACE2 Analysis")
message("========================================\n")

# Check if CytoTRACE2 is available
cytotrace_available <- requireNamespace("CytoTRACE2", quietly = TRUE)

if (cytotrace_available) {
  message("Running CytoTRACE2 on epithelial spots...")

  library(CytoTRACE2)

  # CytoTRACE2 expects a Seurat object or matrix
  # Run on epithelial spots only
  # Note: CytoTRACE2 should auto-detect the assay, don't pass assay parameter

  message(sprintf("  Default assay: %s", DefaultAssay(spatial_epi)))
  message(sprintf("  Available assays: %s", paste(Assays(spatial_epi), collapse = ", ")))

  # First try: Pass matrix directly to avoid Seurat assay detection issues
  cytotrace_result <- tryCatch({
    message("  Extracting counts matrix for CytoTRACE2...")
    # Extract counts as matrix (CytoTRACE2 works better with raw matrices)
    counts_mat <- as.matrix(LayerData(spatial_epi, layer = "counts"))
    message(sprintf("  Matrix dimensions: %d genes x %d spots", nrow(counts_mat), ncol(counts_mat)))

    cytotrace2(counts_mat,
               is_seurat = FALSE,
               species = "human",
               parallelize_models = (n_cores > 1),
               parallelize_smoothing = (n_cores > 1),
               ncores = n_cores,
               seed = config$reproducibility$seed)
  }, error = function(e) {
    message(sprintf("  CytoTRACE2 matrix approach error: %s", e$message))
    message("  Trying Seurat object directly without assay parameter...")

    tryCatch({
      cytotrace2(spatial_epi,
                 is_seurat = TRUE,
                 species = "human",
                 slot_type = "counts",
                 parallelize_models = FALSE,
                 parallelize_smoothing = FALSE)
    }, error = function(e2) {
      message(sprintf("  Seurat approach also failed: %s", e2$message))
      NULL
    })
  })

  if (!is.null(cytotrace_result)) {
    # CytoTRACE2 returns different formats depending on input type
    # When is_seurat = FALSE (matrix input), it returns a data.frame
    # When is_seurat = TRUE (Seurat input), it returns an annotated Seurat object

    if (is.data.frame(cytotrace_result)) {
      # Matrix input returns data.frame with cell barcodes as row names
      message("  CytoTRACE2 returned data.frame format")

      # Check expected columns
      if ("CytoTRACE2_Score" %in% colnames(cytotrace_result)) {
        # Map scores back to spatial object by cell barcode
        common_cells <- intersect(rownames(cytotrace_result), colnames(spatial_epi))
        message(sprintf("  Mapping %d scores to spatial object", length(common_cells)))

        spatial_epi$cytotrace2_score <- NA
        spatial_epi$cytotrace2_potency <- NA

        spatial_epi$cytotrace2_score[common_cells] <- cytotrace_result[common_cells, "CytoTRACE2_Score"]

        if ("CytoTRACE2_Potency" %in% colnames(cytotrace_result)) {
          spatial_epi$cytotrace2_potency[common_cells] <- cytotrace_result[common_cells, "CytoTRACE2_Potency"]
        } else if ("Potency" %in% colnames(cytotrace_result)) {
          spatial_epi$cytotrace2_potency[common_cells] <- cytotrace_result[common_cells, "Potency"]
        }

        message(sprintf("  CytoTRACE2 scores computed for %d spots",
                        sum(!is.na(spatial_epi$cytotrace2_score))))
        message(sprintf("  Score range: [%.3f, %.3f]",
                        min(spatial_epi$cytotrace2_score, na.rm = TRUE),
                        max(spatial_epi$cytotrace2_score, na.rm = TRUE)))
        if (sum(!is.na(spatial_epi$cytotrace2_potency)) > 0) {
          message("  Potency distribution:")
          print(table(spatial_epi$cytotrace2_potency, useNA = "ifany"))
        }

        results$cytotrace <- list(
          success = TRUE,
          n_scored = sum(!is.na(spatial_epi$cytotrace2_score)),
          score_range = range(spatial_epi$cytotrace2_score, na.rm = TRUE),
          potency_distribution = table(spatial_epi$cytotrace2_potency, useNA = "ifany")
        )
      } else {
        message("  Warning: CytoTRACE2 data.frame missing expected columns")
        message(sprintf("  Available columns: %s", paste(colnames(cytotrace_result), collapse = ", ")))
        results$cytotrace <- list(success = FALSE, reason = "missing_columns")
      }

    } else if (inherits(cytotrace_result, "Seurat")) {
      # Seurat input returns annotated Seurat object
      message("  CytoTRACE2 returned Seurat object format")

      if ("CytoTRACE2_Score" %in% colnames(cytotrace_result@meta.data)) {
        spatial_epi$cytotrace2_score <- cytotrace_result$CytoTRACE2_Score
        spatial_epi$cytotrace2_potency <- cytotrace_result$CytoTRACE2_Potency

        message(sprintf("  CytoTRACE2 scores computed for %d spots",
                        sum(!is.na(spatial_epi$cytotrace2_score))))
        message(sprintf("  Score range: [%.3f, %.3f]",
                        min(spatial_epi$cytotrace2_score, na.rm = TRUE),
                        max(spatial_epi$cytotrace2_score, na.rm = TRUE)))
        message("  Potency distribution:")
        print(table(spatial_epi$cytotrace2_potency))

        results$cytotrace <- list(
          success = TRUE,
          n_scored = sum(!is.na(spatial_epi$cytotrace2_score)),
          score_range = range(spatial_epi$cytotrace2_score, na.rm = TRUE),
          potency_distribution = table(spatial_epi$cytotrace2_potency)
        )
      } else {
        message("  Warning: CytoTRACE2 Seurat output missing expected columns")
        results$cytotrace <- list(success = FALSE, reason = "missing_columns")
      }

    } else {
      message(sprintf("  Warning: Unexpected CytoTRACE2 output type: %s", class(cytotrace_result)[1]))
      results$cytotrace <- list(success = FALSE, reason = "unexpected_type")
    }
  } else {
    message("  CytoTRACE2 failed - will proceed without differentiation scores")
    results$cytotrace <- list(success = FALSE, reason = "execution_error")
  }
} else {
  message("CytoTRACE2 not installed - skipping")
  message("  Install with: devtools::install_github('digitalcytometry/CytoTRACE2')")
  results$cytotrace <- list(success = FALSE, reason = "not_installed")
}


# ==============================================================================
# SECTION 2: LABEL-FREE ANALYSIS - SPATIAL EXPRESSION GRADIENTS
# ==============================================================================

message("\n========================================")
message("Section 2: Spatial Expression Gradients")
message("========================================\n")

# --- Helper function to extract spatial coordinates ---------------------------
get_spatial_coordinates <- function(seurat_obj) {
  coords <- NULL

  # Try different methods to get coordinates
  if (length(seurat_obj@images) > 0) {
    # Standard Visium format
    for (img_name in names(seurat_obj@images)) {
      img_coords <- GetTissueCoordinates(seurat_obj, image = img_name)
      if (!is.null(img_coords)) {
        img_coords$image <- img_name
        coords <- rbind(coords, img_coords)
      }
    }
  }

  # Fallback: check for coordinate columns in metadata
  if (is.null(coords) || nrow(coords) == 0) {
    meta <- seurat_obj@meta.data
    coord_cols <- c("imagerow", "imagecol", "row", "col", "x", "y")
    found_cols <- intersect(coord_cols, colnames(meta))

    if (length(found_cols) >= 2) {
      coords <- meta[, found_cols, drop = FALSE]
      rownames(coords) <- colnames(seurat_obj)
    }
  }

  if (is.null(coords) || nrow(coords) == 0) {
    stop("Could not extract spatial coordinates")
  }

  # Standardize column names
  if ("imagerow" %in% colnames(coords) && "imagecol" %in% colnames(coords)) {
    coords$x <- coords$imagecol
    coords$y <- coords$imagerow
  } else if (!"x" %in% colnames(coords)) {
    colnames(coords)[1:2] <- c("x", "y")
  }

  # Filter to cells in object
  common_cells <- intersect(rownames(coords), colnames(seurat_obj))
  coords <- coords[common_cells, ]

  return(coords)
}

# --- Get coordinates for epithelial spots -------------------------------------
message("Extracting spatial coordinates...")
coords <- get_spatial_coordinates(spatial_epi)
message(sprintf("  Coordinates for %d spots", nrow(coords)))

# Store coordinates in metadata
spatial_epi$spatial_x <- coords[colnames(spatial_epi), "x"]
spatial_epi$spatial_y <- coords[colnames(spatial_epi), "y"]

# --- Run PCA if not present ---------------------------------------------------
message("\n--- Computing PCA for expression gradients ---")

if (!"pca" %in% names(spatial_epi@reductions)) {
  message("  Running PCA...")
  spatial_epi <- NormalizeData(spatial_epi, verbose = FALSE)
  spatial_epi <- FindVariableFeatures(spatial_epi, nfeatures = 2000, verbose = FALSE)
  spatial_epi <- ScaleData(spatial_epi, verbose = FALSE)
  spatial_epi <- RunPCA(spatial_epi, npcs = 30, verbose = FALSE)
}

# Extract PC embeddings
pca_embeddings <- Embeddings(spatial_epi, "pca")[, 1:10]

# --- Compute spatial autocorrelation of PC1 -----------------------------------
message("\n--- Spatial autocorrelation of transcriptional variation ---")

# Function to compute Moran's I
compute_morans_i <- function(values, coords, k = 6) {
  n <- length(values)

  # Build k-NN weights
  dist_matrix <- as.matrix(dist(coords[, c("x", "y")]))
  W <- matrix(0, n, n)

  for (i in 1:n) {
    dists <- dist_matrix[i, ]
    dists[i] <- Inf
    nn_idx <- order(dists)[1:k]
    W[i, nn_idx] <- 1
  }

  # Row standardize
  row_sums <- rowSums(W)
  row_sums[row_sums == 0] <- 1
  W <- W / row_sums

  # Moran's I calculation
  values_centered <- values - mean(values, na.rm = TRUE)
  n_valid <- sum(!is.na(values))

  numerator <- sum(W * outer(values_centered, values_centered), na.rm = TRUE)
  denominator <- sum(values_centered^2, na.rm = TRUE)

  S0 <- sum(W, na.rm = TRUE)
  I <- (n_valid / S0) * (numerator / denominator)

  # Expected value under null
  E_I <- -1 / (n_valid - 1)

  # Permutation test for p-value
  n_perm <- 999
  I_perm <- numeric(n_perm)

  for (p in 1:n_perm) {
    perm_values <- sample(values)
    perm_centered <- perm_values - mean(perm_values, na.rm = TRUE)
    num_perm <- sum(W * outer(perm_centered, perm_centered), na.rm = TRUE)
    den_perm <- sum(perm_centered^2, na.rm = TRUE)
    I_perm[p] <- (n_valid / S0) * (num_perm / den_perm)
  }

  p_value <- mean(I_perm >= I)

  return(list(
    I = I,
    expected = E_I,
    p_value = p_value,
    interpretation = ifelse(I > E_I & p_value < 0.05, "CLUSTERED",
                            ifelse(I < E_I & p_value < 0.05, "DISPERSED", "RANDOM"))
  ))
}

# Compute per-sample spatial autocorrelation
samples <- unique(spatial_epi$sample)
gradient_results <- list()

for (samp in samples) {
  message(sprintf("  Processing sample: %s", samp))

  samp_idx <- which(spatial_epi$sample == samp)
  samp_coords <- coords[colnames(spatial_epi)[samp_idx], ]
  samp_pc1 <- pca_embeddings[samp_idx, 1]

  # Moran's I for PC1
  morans_pc1 <- compute_morans_i(samp_pc1, samp_coords, k = 6)

  # Also test CytoTRACE if available
  morans_cyto <- NULL
  if ("cytotrace2_score" %in% colnames(spatial_epi@meta.data)) {
    samp_cyto <- spatial_epi$cytotrace2_score[samp_idx]
    if (sum(!is.na(samp_cyto)) > 50) {
      morans_cyto <- compute_morans_i(samp_cyto, samp_coords, k = 6)
    }
  }

  gradient_results[[samp]] <- list(
    n_spots = length(samp_idx),
    pc1_morans_i = morans_pc1,
    cytotrace_morans_i = morans_cyto
  )

  message(sprintf("    PC1 Moran's I: %.3f (p = %.4f) - %s",
                  morans_pc1$I, morans_pc1$p_value, morans_pc1$interpretation))
  if (!is.null(morans_cyto)) {
    message(sprintf("    CytoTRACE Moran's I: %.3f (p = %.4f) - %s",
                    morans_cyto$I, morans_cyto$p_value, morans_cyto$interpretation))
  }
}

results$spatial_gradients <- gradient_results

# --- Create gradient summary table --------------------------------------------
gradient_summary <- do.call(rbind, lapply(names(gradient_results), function(samp) {
  res <- gradient_results[[samp]]
  data.frame(
    sample = samp,
    n_spots = res$n_spots,
    pc1_morans_i = res$pc1_morans_i$I,
    pc1_p_value = res$pc1_morans_i$p_value,
    pc1_interpretation = res$pc1_morans_i$interpretation,
    cytotrace_morans_i = ifelse(!is.null(res$cytotrace_morans_i),
                                res$cytotrace_morans_i$I, NA),
    cytotrace_p_value = ifelse(!is.null(res$cytotrace_morans_i),
                               res$cytotrace_morans_i$p_value, NA),
    stringsAsFactors = FALSE
  )
}))

write.csv(gradient_summary,
          file.path(tables_dir, "s01_spatial_gradient_stats.csv"),
          row.names = FALSE)
message(sprintf("\nSaved: %s", file.path(tables_dir, "s01_spatial_gradient_stats.csv")))


# ==============================================================================
# SECTION 3: LABEL TRANSFER (GSE → SPATIAL)
# ==============================================================================

message("\n========================================")
message("Section 3: Label Transfer (GSE → Spatial)")
message("========================================\n")

# --- Find common genes --------------------------------------------------------
message("Finding common genes...")
common_genes <- intersect(rownames(spatial_epi), rownames(gse_epi))
message(sprintf("  Common genes: %d", length(common_genes)))

if (length(common_genes) < 1000) {
  warning("Low gene overlap - label transfer may be unreliable")
}

# --- Prepare minimal metadata for merge ---------------------------------------
# FIX: Clean metadata to avoid factor level issues in Harmony
message("\n--- Preparing objects for Harmony integration ---")
message("  Cleaning metadata to prevent factor level conflicts...")

# Create new Seurat objects with minimal, clean metadata
# This prevents the "contrasts can be applied only to factors with 2 or more levels" error

# For spatial: extract counts and create clean object
# FIX: Seurat v5 uses 'layer' instead of 'slot'
spatial_counts <- tryCatch({
  # Try Seurat v5 syntax first
  LayerData(spatial_epi, assay = "Spatial", layer = "counts")[common_genes, ]
}, error = function(e) {
  # Fallback for different assay names or older syntax
  tryCatch({
    LayerData(spatial_epi, layer = "counts")[common_genes, ]
  }, error = function(e2) {
    # Try getting counts from default assay
    GetAssayData(spatial_epi, layer = "counts")[common_genes, ]
  })
})
spatial_meta_clean <- data.frame(
  row.names = colnames(spatial_epi),
  dataset = "Spatial",
  sample = as.character(spatial_epi$sample),
  orig_barcode = colnames(spatial_epi),
  stringsAsFactors = FALSE
)

spatial_for_merge <- CreateSeuratObject(
  counts = spatial_counts,
  meta.data = spatial_meta_clean,
  project = "Spatial"
)

# For GSE: extract counts and create clean object
# FIX: Seurat v5 uses 'layer' instead of 'slot'
gse_counts <- tryCatch({
  LayerData(gse_epi, assay = "RNA", layer = "counts")[common_genes, ]
}, error = function(e) {
  tryCatch({
    LayerData(gse_epi, layer = "counts")[common_genes, ]
  }, error = function(e2) {
    GetAssayData(gse_epi, layer = "counts")[common_genes, ]
  })
})
gse_meta_clean <- data.frame(
  row.names = colnames(gse_epi),
  dataset = "GSE",
  sample = "GSE_pooled",  # or use actual sample info if available
  orig_barcode = colnames(gse_epi),
  cell_type = as.character(gse_epi$cell_type),
  stringsAsFactors = FALSE
)

gse_for_merge <- CreateSeuratObject(
  counts = gse_counts,
  meta.data = gse_meta_clean,
  project = "GSE"
)

# Merge with clean objects
message("  Merging objects...")
merged_obj <- merge(spatial_for_merge, gse_for_merge, add.cell.ids = c("Spatial", "GSE"))
message(sprintf("  Merged object: %d cells", ncol(merged_obj)))

# Ensure dataset is a proper factor with both levels
merged_obj$dataset <- factor(merged_obj$dataset, levels = c("Spatial", "GSE"))
message(sprintf("  Dataset distribution: Spatial=%d, GSE=%d",
                sum(merged_obj$dataset == "Spatial"),
                sum(merged_obj$dataset == "GSE")))

# --- Run standard preprocessing -----------------------------------------------
message("\n--- Preprocessing merged object ---")
merged_obj <- NormalizeData(merged_obj, verbose = FALSE)
merged_obj <- FindVariableFeatures(merged_obj, nfeatures = 2000, verbose = FALSE)
merged_obj <- ScaleData(merged_obj, verbose = FALSE)
merged_obj <- RunPCA(merged_obj, npcs = 30, verbose = FALSE)

# --- Run Harmony integration --------------------------------------------------
message("\n--- Running Harmony batch correction ---")

if (!requireNamespace("harmony", quietly = TRUE)) {
  stop("harmony package required. Install with: install.packages('harmony')")
}

library(harmony)

merged_obj <- RunHarmony(merged_obj,
                         group.by.vars = "dataset",
                         dims.use = 1:30,
                         verbose = TRUE,
                         plot_convergence = FALSE)

message("  Harmony integration complete")

# --- Compute UMAP on Harmony embedding ----------------------------------------
message("  Computing UMAP...")
merged_obj <- RunUMAP(merged_obj,
                      reduction = "harmony",
                      dims = 1:30,
                      verbose = FALSE)

# --- KNN label transfer -------------------------------------------------------
message("\n--- Performing KNN label transfer ---")

# Get Harmony embeddings
harmony_embeddings <- Embeddings(merged_obj, "harmony")

# Split into spatial and GSE cells
spatial_cells <- grep("^Spatial_", colnames(merged_obj), value = TRUE)
gse_cells <- grep("^GSE_", colnames(merged_obj), value = TRUE)

spatial_embedding <- harmony_embeddings[spatial_cells, ]
gse_embedding <- harmony_embeddings[gse_cells, ]

# Get GSE labels
gse_labels <- merged_obj$cell_type[gse_cells]

# KNN transfer (k = 30)
k <- 30
message(sprintf("  Using k = %d nearest neighbors", k))

# Use RANN for fast KNN
if (!requireNamespace("RANN", quietly = TRUE)) {
  stop("RANN package required. Install with: install.packages('RANN')")
}

library(RANN)

# Find k nearest GSE neighbors for each spatial spot
nn_result <- nn2(gse_embedding, spatial_embedding, k = k)
nn_indices <- nn_result$nn.idx
nn_distances <- nn_result$nn.dists

# Weighted voting based on distance
transferred_labels <- character(nrow(spatial_embedding))
transfer_confidence <- numeric(nrow(spatial_embedding))

for (i in 1:nrow(spatial_embedding)) {
  neighbor_labels <- gse_labels[nn_indices[i, ]]
  neighbor_dists <- nn_distances[i, ]

  # Inverse distance weighting
  weights <- 1 / (neighbor_dists + 1e-6)
  weights <- weights / sum(weights)

  # Aggregate votes by label
  label_weights <- tapply(weights, neighbor_labels, sum)

  # Assign most voted label
  transferred_labels[i] <- names(which.max(label_weights))
  transfer_confidence[i] <- max(label_weights)
}

names(transferred_labels) <- spatial_cells
names(transfer_confidence) <- spatial_cells

# --- Add transferred labels to spatial object ---------------------------------
# Map back to original cell names (remove "Spatial_" prefix)
original_names <- gsub("^Spatial_", "", spatial_cells)

spatial_epi$transferred_subtype <- NA
spatial_epi$transfer_confidence <- NA

spatial_epi$transferred_subtype[original_names] <- transferred_labels
spatial_epi$transfer_confidence[original_names] <- transfer_confidence

# Flag low confidence transfers
confidence_threshold <- 0.5
spatial_epi$transfer_confident <- spatial_epi$transfer_confidence >= confidence_threshold

message("\n  Label transfer results:")
print(table(spatial_epi$transferred_subtype))
message(sprintf("\n  High confidence (>= %.1f): %d spots (%.1f%%)",
                confidence_threshold,
                sum(spatial_epi$transfer_confident, na.rm = TRUE),
                mean(spatial_epi$transfer_confident, na.rm = TRUE) * 100))

results$label_transfer <- list(
  k = k,
  n_transferred = sum(!is.na(spatial_epi$transferred_subtype)),
  confidence_threshold = confidence_threshold,
  n_confident = sum(spatial_epi$transfer_confident, na.rm = TRUE),
  pct_confident = mean(spatial_epi$transfer_confident, na.rm = TRUE) * 100,
  label_distribution = table(spatial_epi$transferred_subtype),
  mean_confidence = mean(spatial_epi$transfer_confidence, na.rm = TRUE)
)

# --- Save transfer confidence table -------------------------------------------
transfer_summary <- data.frame(
  transferred_subtype = names(table(spatial_epi$transferred_subtype)),
  n_spots = as.numeric(table(spatial_epi$transferred_subtype)),
  stringsAsFactors = FALSE
)
transfer_summary$pct_total <- transfer_summary$n_spots / sum(transfer_summary$n_spots) * 100

# Add mean confidence per subtype
conf_by_subtype <- tapply(spatial_epi$transfer_confidence,
                          spatial_epi$transferred_subtype,
                          mean, na.rm = TRUE)
transfer_summary$mean_confidence <- conf_by_subtype[transfer_summary$transferred_subtype]

write.csv(transfer_summary,
          file.path(tables_dir, "s01_label_transfer_summary.csv"),
          row.names = FALSE)
message(sprintf("\nSaved: %s", file.path(tables_dir, "s01_label_transfer_summary.csv")))


# ==============================================================================
# SECTION 4: SPATIAL CLUSTERING (RIPLEY'S K)
# ==============================================================================

message("\n========================================")
message("Section 4: Spatial Clustering Analysis")
message("========================================\n")

# --- Function to compute Ripley's K -------------------------------------------
compute_ripleys_k <- function(coords, labels, target_label, r_values = NULL, n_sim = 99) {

  # Filter to target label
  target_idx <- which(labels == target_label)
  if (length(target_idx) < 10) {
    return(list(success = FALSE, reason = "insufficient_points"))
  }

  target_coords <- coords[target_idx, c("x", "y")]
  all_coords <- coords[, c("x", "y")]

  # Define study region
  x_range <- range(all_coords$x)
  y_range <- range(all_coords$y)
  area <- diff(x_range) * diff(y_range)

  # Default r values
  if (is.null(r_values)) {
    max_r <- min(diff(x_range), diff(y_range)) / 4
    r_values <- seq(0, max_r, length.out = 50)
  }

  n_target <- nrow(target_coords)
  lambda <- n_target / area  # intensity

  # Compute observed K(r)
  dist_matrix <- as.matrix(dist(target_coords))
  diag(dist_matrix) <- Inf

  K_obs <- numeric(length(r_values))
  for (i in seq_along(r_values)) {
    r <- r_values[i]
    K_obs[i] <- sum(dist_matrix <= r) / (n_target * lambda)
  }

  # Expected K under CSR: pi * r^2
  K_csr <- pi * r_values^2

  # L transformation: L(r) = sqrt(K(r)/pi) - r
  # Under CSR, L(r) = 0
  L_obs <- sqrt(K_obs / pi) - r_values

  # Simulation envelope (random labeling)
  K_sim <- matrix(NA, n_sim, length(r_values))

  for (sim in 1:n_sim) {
    # Random permutation of labels
    sim_idx <- sample(1:nrow(all_coords), n_target)
    sim_coords <- all_coords[sim_idx, ]

    sim_dist <- as.matrix(dist(sim_coords))
    diag(sim_dist) <- Inf

    for (i in seq_along(r_values)) {
      r <- r_values[i]
      K_sim[sim, i] <- sum(sim_dist <= r) / (n_target * lambda)
    }
  }

  # Envelope
  K_lower <- apply(K_sim, 2, quantile, 0.025)
  K_upper <- apply(K_sim, 2, quantile, 0.975)

  # Significance: max deviation outside envelope
  above_envelope <- K_obs > K_upper
  below_envelope <- K_obs < K_lower

  significant_clustering <- any(above_envelope)
  significant_dispersion <- any(below_envelope)

  return(list(
    success = TRUE,
    r = r_values,
    K_observed = K_obs,
    K_expected = K_csr,
    K_lower = K_lower,
    K_upper = K_upper,
    L_observed = L_obs,
    n_points = n_target,
    area = area,
    significant_clustering = significant_clustering,
    significant_dispersion = significant_dispersion
  ))
}

# --- Compute Ripley's K for each subtype per sample ---------------------------
message("Computing Ripley's K for each transferred subtype...")

subtypes <- unique(spatial_epi$transferred_subtype)
subtypes <- subtypes[!is.na(subtypes)]

ripley_results <- list()

for (samp in samples) {
  message(sprintf("\n  Sample: %s", samp))

  samp_idx <- which(spatial_epi$sample == samp)
  samp_coords <- coords[colnames(spatial_epi)[samp_idx], ]
  samp_labels <- spatial_epi$transferred_subtype[samp_idx]

  ripley_results[[samp]] <- list()

  for (st in subtypes) {
    n_st <- sum(samp_labels == st, na.rm = TRUE)
    message(sprintf("    %s (n=%d)...", st, n_st))

    if (n_st >= 10) {
      rk <- compute_ripleys_k(samp_coords, samp_labels, st, n_sim = 99)
      ripley_results[[samp]][[st]] <- rk

      if (rk$success) {
        message(sprintf("      Clustered: %s", rk$significant_clustering))
      }
    } else {
      ripley_results[[samp]][[st]] <- list(success = FALSE, reason = "insufficient_points", n = n_st)
    }
  }
}

results$ripley <- ripley_results

# --- Create summary table -----------------------------------------------------
ripley_summary <- do.call(rbind, lapply(names(ripley_results), function(samp) {
  do.call(rbind, lapply(names(ripley_results[[samp]]), function(st) {
    res <- ripley_results[[samp]][[st]]
    data.frame(
      sample = samp,
      subtype = st,
      n_spots = ifelse(res$success, res$n_points, res$n),
      success = res$success,
      significant_clustering = ifelse(res$success, res$significant_clustering, NA),
      stringsAsFactors = FALSE
    )
  }))
}))

write.csv(ripley_summary,
          file.path(tables_dir, "s01_ripley_k_results.csv"),
          row.names = FALSE)
message(sprintf("\nSaved: %s", file.path(tables_dir, "s01_ripley_k_results.csv")))


# ==============================================================================
# SECTION 5: NEIGHBORHOOD ENRICHMENT ANALYSIS
# ==============================================================================

message("\n========================================")
message("Section 5: Neighborhood Enrichment")
message("========================================\n")

# --- Function to compute neighborhood enrichment ------------------------------
compute_neighborhood_enrichment <- function(coords, labels, k = 6) {

  unique_labels <- unique(labels[!is.na(labels)])
  n_labels <- length(unique_labels)

  # Build KNN graph
  dist_matrix <- as.matrix(dist(coords[, c("x", "y")]))
  n <- nrow(coords)

  # Count label frequencies
  label_counts <- table(labels)
  label_freqs <- label_counts / sum(label_counts)

  # Initialize enrichment matrix
  enrichment <- matrix(0, n_labels, n_labels)
  rownames(enrichment) <- colnames(enrichment) <- unique_labels

  observed_counts <- matrix(0, n_labels, n_labels)
  rownames(observed_counts) <- colnames(observed_counts) <- unique_labels

  # For each cell, count neighbor labels
  for (i in 1:n) {
    if (is.na(labels[i])) next

    dists <- dist_matrix[i, ]
    dists[i] <- Inf
    nn_idx <- order(dists)[1:k]

    neighbor_labels <- labels[nn_idx]
    neighbor_labels <- neighbor_labels[!is.na(neighbor_labels)]

    if (length(neighbor_labels) > 0) {
      for (nl in neighbor_labels) {
        observed_counts[labels[i], nl] <- observed_counts[labels[i], nl] + 1
      }
    }
  }

  # Calculate enrichment (observed / expected)
  # Expected = row_total * col_freq
  row_totals <- rowSums(observed_counts)

  for (i in 1:n_labels) {
    for (j in 1:n_labels) {
      expected <- row_totals[i] * label_freqs[unique_labels[j]]
      if (expected > 0) {
        enrichment[i, j] <- log2((observed_counts[i, j] + 1) / (expected + 1))
      }
    }
  }

  # Homotypic scores (diagonal)
  homotypic_scores <- diag(enrichment)
  names(homotypic_scores) <- unique_labels

  return(list(
    enrichment_matrix = enrichment,
    observed_counts = observed_counts,
    label_frequencies = label_freqs,
    homotypic_scores = homotypic_scores
  ))
}

# --- Compute per-sample and pooled neighborhood enrichment --------------------
message("Computing neighborhood enrichment...")

neighborhood_results <- list()

# Per-sample
for (samp in samples) {
  message(sprintf("  Sample: %s", samp))

  samp_idx <- which(spatial_epi$sample == samp)
  samp_coords <- coords[colnames(spatial_epi)[samp_idx], ]
  samp_labels <- spatial_epi$transferred_subtype[samp_idx]

  neighborhood_results[[samp]] <- compute_neighborhood_enrichment(samp_coords, samp_labels, k = 6)
}

# Pooled (using per-sample results)
# Average enrichment matrices
all_matrices <- lapply(neighborhood_results, function(x) x$enrichment_matrix)
common_labels <- Reduce(intersect, lapply(all_matrices, rownames))

if (length(common_labels) > 1) {
  pooled_enrichment <- Reduce(`+`, lapply(all_matrices, function(m) m[common_labels, common_labels])) / length(all_matrices)

  neighborhood_results$pooled <- list(
    enrichment_matrix = pooled_enrichment,
    homotypic_scores = diag(pooled_enrichment)
  )
}

results$neighborhood <- neighborhood_results

# --- Save neighborhood enrichment matrix --------------------------------------
if (!is.null(neighborhood_results$pooled)) {
  enrichment_df <- as.data.frame(neighborhood_results$pooled$enrichment_matrix)
  enrichment_df$row_label <- rownames(enrichment_df)
  enrichment_df <- enrichment_df[, c("row_label", setdiff(names(enrichment_df), "row_label"))]

  write.csv(enrichment_df,
            file.path(tables_dir, "s01_neighborhood_enrichment_matrix.csv"),
            row.names = FALSE)
  message(sprintf("Saved: %s", file.path(tables_dir, "s01_neighborhood_enrichment_matrix.csv")))
}

# --- Save homotypic scores ----------------------------------------------------
homotypic_df <- data.frame(
  subtype = names(neighborhood_results$pooled$homotypic_scores),
  homotypic_enrichment = neighborhood_results$pooled$homotypic_scores,
  stringsAsFactors = FALSE
)

# Add tissue abundance
subtype_freq <- table(spatial_epi$transferred_subtype) / sum(!is.na(spatial_epi$transferred_subtype))
homotypic_df$tissue_abundance <- subtype_freq[homotypic_df$subtype]
homotypic_df$enrichment_ratio <- 2^homotypic_df$homotypic_enrichment / homotypic_df$tissue_abundance

write.csv(homotypic_df,
          file.path(tables_dir, "s01_homotypic_scores.csv"),
          row.names = FALSE)
message(sprintf("Saved: %s", file.path(tables_dir, "s01_homotypic_scores.csv")))


# ==============================================================================
# SECTION 6: BRIDGE ANALYSIS - CYTOTRACE BY SUBTYPE
# ==============================================================================

message("\n========================================")
message("Section 6: CytoTRACE by Transferred Subtype")
message("========================================\n")

if ("cytotrace2_score" %in% colnames(spatial_epi@meta.data) &&
    sum(!is.na(spatial_epi$cytotrace2_score)) > 100) {

  message("Comparing CytoTRACE2 scores across subtypes...")

  # Kruskal-Wallis test
  kw_test <- kruskal.test(cytotrace2_score ~ transferred_subtype, data = spatial_epi@meta.data)

  message(sprintf("  Kruskal-Wallis: chi-squared = %.2f, p = %.2e",
                  kw_test$statistic, kw_test$p.value))

  # Pairwise Wilcoxon tests
  pairwise_tests <- pairwise.wilcox.test(
    spatial_epi$cytotrace2_score,
    spatial_epi$transferred_subtype,
    p.adjust.method = "BH"
  )

  # Summary by subtype
  # FIX: Use dplyr:: explicitly to avoid plyr conflict
  cyto_by_subtype <- spatial_epi@meta.data %>%
    dplyr::filter(!is.na(transferred_subtype), !is.na(cytotrace2_score)) %>%
    dplyr::group_by(transferred_subtype) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_cytotrace = mean(cytotrace2_score, na.rm = TRUE),
      sd_cytotrace = sd(cytotrace2_score, na.rm = TRUE),
      median_cytotrace = median(cytotrace2_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_cytotrace))

  message("\n  CytoTRACE2 scores by subtype (ordered by mean):")
  print(cyto_by_subtype)

  results$cytotrace_by_subtype <- list(
    kruskal_wallis = kw_test,
    pairwise = pairwise_tests,
    summary = cyto_by_subtype
  )

  write.csv(cyto_by_subtype,
            file.path(tables_dir, "s01_cytotrace_by_subtype.csv"),
            row.names = FALSE)
  message(sprintf("\nSaved: %s", file.path(tables_dir, "s01_cytotrace_by_subtype.csv")))

} else {
  message("CytoTRACE2 scores not available - skipping subtype comparison")
  results$cytotrace_by_subtype <- NULL
}


# ==============================================================================
# SECTION 7: COMPARISON - SPATIAL VS DIFFERENTIATION GRADIENTS
# ==============================================================================

message("\n========================================")
message("Section 7: Spatial vs Differentiation Gradients")
message("========================================\n")

comparison_results <- list()

# --- Per-sample correlation analysis ------------------------------------------
message("Correlating spatial position with differentiation metrics...")

for (samp in samples) {
  message(sprintf("\n  Sample: %s", samp))

  samp_idx <- which(spatial_epi$sample == samp)
  samp_data <- spatial_epi@meta.data[samp_idx, ]
  samp_coords <- coords[colnames(spatial_epi)[samp_idx], ]

  # Spatial PC1 (dominant axis of spatial variation)
  samp_spatial <- cbind(samp_coords$x, samp_coords$y)
  spatial_pca <- prcomp(samp_spatial, scale. = TRUE)
  spatial_pc1 <- spatial_pca$x[, 1]

  # Expression PC1
  expr_pc1 <- Embeddings(spatial_epi, "pca")[samp_idx, 1]

  # CytoTRACE if available
  cyto_score <- if ("cytotrace2_score" %in% colnames(samp_data)) {
    samp_data$cytotrace2_score
  } else {
    rep(NA, length(samp_idx))
  }

  # Correlations
  cor_spatial_expr <- cor(spatial_pc1, expr_pc1, method = "spearman", use = "complete.obs")
  cor_spatial_cyto <- if (sum(!is.na(cyto_score)) > 50) {
    cor(spatial_pc1, cyto_score, method = "spearman", use = "complete.obs")
  } else {
    NA
  }
  cor_expr_cyto <- if (sum(!is.na(cyto_score)) > 50) {
    cor(expr_pc1, cyto_score, method = "spearman", use = "complete.obs")
  } else {
    NA
  }

  comparison_results[[samp]] <- list(
    cor_spatial_expression = cor_spatial_expr,
    cor_spatial_cytotrace = cor_spatial_cyto,
    cor_expression_cytotrace = cor_expr_cyto
  )

  message(sprintf("    Spatial PC1 vs Expression PC1: rho = %.3f", cor_spatial_expr))
  if (!is.na(cor_spatial_cyto)) {
    message(sprintf("    Spatial PC1 vs CytoTRACE: rho = %.3f", cor_spatial_cyto))
  }
  if (!is.na(cor_expr_cyto)) {
    message(sprintf("    Expression PC1 vs CytoTRACE: rho = %.3f", cor_expr_cyto))
  }
}

results$gradient_comparison <- comparison_results

# --- Create comparison summary table ------------------------------------------
comparison_summary <- do.call(rbind, lapply(names(comparison_results), function(samp) {
  res <- comparison_results[[samp]]
  data.frame(
    sample = samp,
    cor_spatial_expression = res$cor_spatial_expression,
    cor_spatial_cytotrace = res$cor_spatial_cytotrace,
    cor_expression_cytotrace = res$cor_expression_cytotrace,
    stringsAsFactors = FALSE
  )
}))

write.csv(comparison_summary,
          file.path(tables_dir, "s01_gradient_comparison.csv"),
          row.names = FALSE)
message(sprintf("\nSaved: %s", file.path(tables_dir, "s01_gradient_comparison.csv")))


# ==============================================================================
# SECTION 8: FIGURE GENERATION
# ==============================================================================

message("\n========================================")
message("Section 8: Generating Figures")
message("========================================\n")

# Get color palette
colors <- tryCatch({
  get_colors(config, "epithelial_subtypes")
}, error = function(e) {
  # Default colors if not in config
  c("E_C1_PE" = "#E41A1C", "E_C2_WC" = "#377EB8", "E_C3_KC" = "#4DAF4A",
    "E_C4_RHCG" = "#984EA3", "E_C5_Prolif" = "#FF7F00")
})

# Ensure colors exist for all subtypes
all_subtypes <- unique(spatial_epi$transferred_subtype)
all_subtypes <- all_subtypes[!is.na(all_subtypes)]
missing_colors <- setdiff(all_subtypes, names(colors))
if (length(missing_colors) > 0) {
  extra_colors <- scales::hue_pal()(length(missing_colors))
  names(extra_colors) <- missing_colors
  colors <- c(colors, extra_colors)
}

# --- Panel A: Spatial distribution of subtypes --------------------------------
message("Creating Panel A: Spatial subtype maps...")

spatial_plots <- list()
for (samp in samples) {
  samp_idx <- which(spatial_epi$sample == samp)
  samp_coords <- coords[colnames(spatial_epi)[samp_idx], ]
  samp_labels <- spatial_epi$transferred_subtype[samp_idx]

  plot_data <- data.frame(
    x = samp_coords$x,
    y = samp_coords$y,
    subtype = samp_labels,
    confident = spatial_epi$transfer_confident[samp_idx]
  )

  p <- ggplot(plot_data, aes(x = x, y = y, color = subtype)) +
    geom_point(aes(alpha = confident), size = 0.8) +
    scale_color_manual(values = colors, na.value = "grey80") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
    coord_fixed() +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(title = samp)

  spatial_plots[[samp]] <- p
}

panel_a <- wrap_plots(spatial_plots, ncol = 2) +
  plot_annotation(title = "A. Spatial distribution of epithelial subtypes")

# --- Panel B: Ripley's K curves -----------------------------------------------
message("Creating Panel B: Ripley's K curves...")

# Combine Ripley's K data for plotting
ripley_plot_data <- do.call(rbind, lapply(names(ripley_results), function(samp) {
  if (samp == "pooled") return(NULL)

  do.call(rbind, lapply(names(ripley_results[[samp]]), function(st) {
    res <- ripley_results[[samp]][[st]]
    if (!res$success) return(NULL)

    data.frame(
      sample = samp,
      subtype = st,
      r = res$r,
      K_observed = res$K_observed,
      K_expected = res$K_expected,
      K_lower = res$K_lower,
      K_upper = res$K_upper,
      stringsAsFactors = FALSE
    )
  }))
}))

# Check if we have valid data for plotting
if (!is.null(ripley_plot_data) && nrow(ripley_plot_data) > 0 &&
    "subtype" %in% colnames(ripley_plot_data) &&
    length(unique(ripley_plot_data$subtype)) > 0) {

  # Average across samples for cleaner plot
  # FIX: Use dplyr:: explicitly to avoid plyr conflict
  ripley_avg <- ripley_plot_data %>%
    dplyr::group_by(subtype, r) %>%
    dplyr::summarise(
      K_observed = mean(K_observed, na.rm = TRUE),
      K_expected = mean(K_expected, na.rm = TRUE),
      K_lower = mean(K_lower, na.rm = TRUE),
      K_upper = mean(K_upper, na.rm = TRUE),
      n_samples = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::ungroup() %>%
    as.data.frame()  # Ensure it's a plain data.frame

  # Only plot subtypes with data
  subtypes_with_data <- unique(ripley_avg$subtype)
  message(sprintf("  Plotting Ripley's K for %d subtypes: %s",
                  length(subtypes_with_data),
                  paste(subtypes_with_data, collapse = ", ")))

  # Debug: check structure
  message(sprintf("  ripley_avg columns: %s", paste(colnames(ripley_avg), collapse = ", ")))
  message(sprintf("  ripley_avg rows: %d", nrow(ripley_avg)))

  if (length(subtypes_with_data) > 0 && nrow(ripley_avg) > 0) {
    panel_b <- ggplot(ripley_avg, aes(x = r)) +
      geom_ribbon(aes(ymin = K_lower, ymax = K_upper), alpha = 0.2) +
      geom_line(aes(y = K_expected), linetype = "dashed", color = "grey50") +
      geom_line(aes(y = K_observed, color = subtype), linewidth = 1) +
      scale_color_manual(values = colors) +
      facet_wrap(~subtype, scales = "free_y") +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(
        title = "B. Ripley's K: Spatial clustering analysis",
        subtitle = "Solid = observed, Dashed = CSR expectation, Ribbon = 95% envelope",
        x = "Distance (r)",
        y = "K(r)"
      )
  } else {
    panel_b <- ggplot() +
      theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "Insufficient data for Ripley's K", size = 5) +
      labs(title = "B. Ripley's K: Insufficient data")
  }
} else {
  message("  No valid Ripley's K data available for plotting")
  panel_b <- ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient data for Ripley's K\n(subtypes need ≥10 spots per sample)", size = 4) +
    labs(title = "B. Ripley's K: Insufficient data")
}

# --- Panel C: Neighborhood enrichment heatmap ---------------------------------
message("Creating Panel C: Neighborhood enrichment heatmap...")

if (!is.null(neighborhood_results$pooled)) {
  enrichment_mat <- neighborhood_results$pooled$enrichment_matrix

  # Convert to long format for ggplot
  enrichment_long <- reshape2::melt(enrichment_mat)
  colnames(enrichment_long) <- c("From", "To", "Enrichment")

  panel_c <- ggplot(enrichment_long, aes(x = To, y = From, fill = Enrichment)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1f", Enrichment)), size = 3) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      name = "Log2\nEnrichment"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank()
    ) +
    labs(title = "C. Neighborhood enrichment",
         subtitle = "Row = cell type, Column = neighbor type")
} else {
  panel_c <- ggplot() + theme_void() +
    labs(title = "C. Neighborhood enrichment: Data unavailable")
}

# --- Panel D: CytoTRACE by subtype --------------------------------------------
message("Creating Panel D: CytoTRACE by subtype...")

if ("cytotrace2_score" %in% colnames(spatial_epi@meta.data) &&
    sum(!is.na(spatial_epi$cytotrace2_score)) > 100) {

  cyto_data <- spatial_epi@meta.data %>%
    dplyr::filter(!is.na(transferred_subtype), !is.na(cytotrace2_score))

  # Order subtypes by mean CytoTRACE
  subtype_order <- cyto_data %>%
    dplyr::group_by(transferred_subtype) %>%
    dplyr::summarise(mean_cyto = mean(cytotrace2_score), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_cyto)) %>%
    dplyr::pull(transferred_subtype)

  cyto_data$transferred_subtype <- factor(cyto_data$transferred_subtype, levels = subtype_order)

  panel_d <- ggplot(cyto_data, aes(x = transferred_subtype, y = cytotrace2_score,
                                   fill = transferred_subtype)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = "D. Differentiation potential by subtype",
      subtitle = "Higher CytoTRACE2 = more stem-like",
      x = "",
      y = "CytoTRACE2 Score"
    )
} else {
  panel_d <- ggplot() + theme_void() +
    labs(title = "D. CytoTRACE: Not available")
}

# --- Panel E: Spatial CytoTRACE gradient --------------------------------------
message("Creating Panel E: Spatial CytoTRACE maps...")

if ("cytotrace2_score" %in% colnames(spatial_epi@meta.data)) {
  cyto_spatial_plots <- list()

  for (samp in samples[1:min(2, length(samples))]) {  # Show first 2 samples
    samp_idx <- which(spatial_epi$sample == samp)
    samp_coords <- coords[colnames(spatial_epi)[samp_idx], ]
    samp_cyto <- spatial_epi$cytotrace2_score[samp_idx]

    plot_data <- data.frame(
      x = samp_coords$x,
      y = samp_coords$y,
      cytotrace = samp_cyto
    )

    p <- ggplot(plot_data, aes(x = x, y = y, color = cytotrace)) +
      geom_point(size = 0.8) +
      scale_color_viridis_c(option = "plasma", na.value = "grey80") +
      coord_fixed() +
      theme_minimal() +
      theme(
        legend.position = "right",
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()
      ) +
      labs(title = samp, color = "CytoTRACE2")

    cyto_spatial_plots[[samp]] <- p
  }

  panel_e <- wrap_plots(cyto_spatial_plots, ncol = 2) +
    plot_annotation(title = "E. Spatial distribution of differentiation potential")
} else {
  panel_e <- ggplot() + theme_void() +
    labs(title = "E. Spatial CytoTRACE: Not available")
}

# --- Panel F: Gradient correlations -------------------------------------------
message("Creating Panel F: Gradient correlation summary...")

if (nrow(comparison_summary) > 0) {
  corr_long <- comparison_summary %>%
    tidyr::pivot_longer(
      cols = starts_with("cor_"),
      names_to = "comparison",
      values_to = "correlation"
    ) %>%
    dplyr::mutate(
      comparison = dplyr::case_when(
        comparison == "cor_spatial_expression" ~ "Spatial vs\nExpression",
        comparison == "cor_spatial_cytotrace" ~ "Spatial vs\nCytoTRACE",
        comparison == "cor_expression_cytotrace" ~ "Expression vs\nCytoTRACE",
        TRUE ~ comparison
      )
    )

  panel_f <- ggplot(corr_long, aes(x = comparison, y = correlation, fill = sample)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    labs(
      title = "F. Correlation between gradient types",
      x = "",
      y = "Spearman ρ"
    )
} else {
  panel_f <- ggplot() + theme_void() +
    labs(title = "F. Gradient correlations: Data unavailable")
}

# --- Combine panels -----------------------------------------------------------
message("\nCombining panels...")

# Top row: A (spatial maps)
# Middle row: B (Ripley's K) + C (neighborhood enrichment)
# Bottom row: D (CytoTRACE violin) + E (CytoTRACE spatial) + F (correlations)

fig1 <- (panel_a) /
  (panel_b | panel_c) /
  (panel_d | panel_e | panel_f) +
  plot_layout(heights = c(2, 1.5, 1.5)) +
  plot_annotation(
    title = "Figure 1: Spatial organization of ACP epithelial subtypes",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# Save figure
ggsave(file.path(fig_main_dir, "s01_figure1_spatial_organization.pdf"),
       fig1, width = 16, height = 18)
ggsave(file.path(fig_main_dir, "s01_figure1_spatial_organization.png"),
       fig1, width = 16, height = 18, dpi = 300)

message(sprintf("Saved: %s", file.path(fig_main_dir, "s01_figure1_spatial_organization.pdf")))


# ==============================================================================
# SECTION 9: SUMMARY & SAVE RESULTS
# ==============================================================================

message("\n========================================")
message("Section 9: Summary & Save")
message("========================================\n")

# --- Generate markdown summary ------------------------------------------------
summary_lines <- c(
  "# Section 1: Spatial Organization Analysis - Summary",
  "",
  sprintf("**Analysis date:** %s", Sys.time()),
  sprintf("**Samples:** %s", paste(samples, collapse = ", ")),
  sprintf("**Epithelial spots:** %d", ncol(spatial_epi)),
  "",
  "## Key Findings",
  "",
  "### 1. Label Transfer (GSE → Spatial)",
  sprintf("- Transferred %d labels with mean confidence %.2f",
          results$label_transfer$n_transferred,
          results$label_transfer$mean_confidence),
  sprintf("- High-confidence transfers (>%.1f): %.1f%%",
          results$label_transfer$confidence_threshold,
          results$label_transfer$pct_confident),
  "",
  "### 2. Spatial Autocorrelation",
  "| Sample | PC1 Moran's I | p-value | Interpretation |",
  "|--------|---------------|---------|----------------|"
)

for (samp in names(results$spatial_gradients)) {
  res <- results$spatial_gradients[[samp]]
  summary_lines <- c(summary_lines,
                     sprintf("| %s | %.3f | %.4f | %s |",
                             samp,
                             res$pc1_morans_i$I,
                             res$pc1_morans_i$p_value,
                             res$pc1_morans_i$interpretation))
}

summary_lines <- c(summary_lines,
                   "",
                   "### 3. Spatial Clustering (Ripley's K)",
                   "Subtypes showing significant spatial clustering:")

clustered_subtypes <- ripley_summary %>%
  dplyr::filter(significant_clustering == TRUE) %>%
  dplyr::group_by(subtype) %>%
  dplyr::summarise(n_samples = dplyr::n(), .groups = "drop")

for (i in 1:nrow(clustered_subtypes)) {
  summary_lines <- c(summary_lines,
                     sprintf("- %s: clustered in %d/%d samples",
                             clustered_subtypes$subtype[i],
                             clustered_subtypes$n_samples[i],
                             length(samples)))
}

summary_lines <- c(summary_lines,
                   "",
                   "### 4. CytoTRACE2 by Subtype")

if (!is.null(results$cytotrace_by_subtype)) {
  summary_lines <- c(summary_lines,
                     sprintf("- Kruskal-Wallis p-value: %.2e",
                             results$cytotrace_by_subtype$kruskal_wallis$p.value),
                     "",
                     "| Subtype | Mean CytoTRACE2 | Interpretation |",
                     "|---------|-----------------|----------------|")

  cyto_summary <- results$cytotrace_by_subtype$summary
  for (i in 1:nrow(cyto_summary)) {
    interp <- ifelse(cyto_summary$mean_cytotrace[i] > 0.5, "More stem-like", "More differentiated")
    summary_lines <- c(summary_lines,
                       sprintf("| %s | %.3f | %s |",
                               cyto_summary$transferred_subtype[i],
                               cyto_summary$mean_cytotrace[i],
                               interp))
  }
} else {
  summary_lines <- c(summary_lines, "CytoTRACE2 analysis not available")
}

summary_lines <- c(summary_lines,
                   "",
                   "## Output Files",
                   "",
                   "### Tables",
                   "- `s01_spatial_gradient_stats.csv`",
                   "- `s01_label_transfer_summary.csv`",
                   "- `s01_ripley_k_results.csv`",
                   "- `s01_neighborhood_enrichment_matrix.csv`",
                   "- `s01_homotypic_scores.csv`",
                   "- `s01_cytotrace_by_subtype.csv`",
                   "- `s01_gradient_comparison.csv`",
                   "",
                   "### Figures",
                   "- `s01_figure1_spatial_organization.pdf`",
                   "",
                   "### Objects",
                   "- `s01_spatial_analysis_results.rds`",
                   "- `s01_spatial_epithelial_annotated.rds`")

writeLines(summary_lines, file.path(tables_dir, "s01_analysis_summary.md"))
message(sprintf("Saved: %s", file.path(tables_dir, "s01_analysis_summary.md")))

# --- Save results object ------------------------------------------------------
saveRDS(results, file.path(objects_dir, "s01_spatial_analysis_results.rds"))
message(sprintf("Saved: %s", file.path(objects_dir, "s01_spatial_analysis_results.rds")))

# --- Save annotated spatial object --------------------------------------------
saveRDS(spatial_epi, file.path(objects_dir, "s01_spatial_epithelial_annotated.rds"))
message(sprintf("Saved: %s", file.path(objects_dir, "s01_spatial_epithelial_annotated.rds")))

# --- Clean up parallelization -------------------------------------------------
plan(sequential)

# --- Final summary ------------------------------------------------------------
message("\n")
message("╔══════════════════════════════════════════════════════════════════════╗")
message("║                    Analysis Complete                                  ║")
message("╚══════════════════════════════════════════════════════════════════════╝")
message("\n")

message("Key findings:")
message(sprintf("  • %d epithelial spots analyzed across %d samples",
                ncol(spatial_epi), length(samples)))
message(sprintf("  • Label transfer: %.1f%% high-confidence",
                results$label_transfer$pct_confident))
message(sprintf("  • Spatial autocorrelation: %d/%d samples show clustered transcriptional gradients",
                sum(sapply(results$spatial_gradients, function(x) x$pc1_morans_i$interpretation == "CLUSTERED")),
                length(samples)))

if (!is.null(results$cytotrace_by_subtype)) {
  message(sprintf("  • CytoTRACE2: Significant differences across subtypes (p = %.2e)",
                  results$cytotrace_by_subtype$kruskal_wallis$p.value))
}

message("\nScript completed successfully!")
message(sprintf("Total runtime: %.1f minutes",
                difftime(Sys.time(), results$metadata$timestamp, units = "mins")))
