#!/usr/bin/env Rscript
# ==============================================================================
# scripts/00_merge_datasets.R
# ==============================================================================
# Merge ACP scRNA-seq and GSE215932 snRNA-seq datasets for combined analysis
#
# This script:
#   1. Loads both Seurat objects
#   2. Extracts RAW COUNTS (not filtered variable features) from each
#   3. Harmonizes metadata columns
#   4. Merges on common genes
#   5. Normalizes the combined data
#   6. Optionally integrates (batch correction)
#
# Usage:
#   Rscript scripts/00_merge_datasets.R
#   Rscript scripts/00_merge_datasets.R --integrate    # With batch correction
#   Rscript scripts/00_merge_datasets.R --config path/to/config.yaml
#
# Output:
#   data/processed/merged_acp_gse.rds
#
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
do_integration <- FALSE

for (arg in args) {

  if (arg == "--integrate") {
    do_integration <- TRUE
  } else if (arg == "--config") {
    idx <- which(args == "--config")
    if (idx < length(args)) config_path <- args[idx + 1]
  } else if (!startsWith(arg, "--") && is.null(config_path)) {
    config_path <- arg
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
  stop("Could not find project root directory.")
}
setwd(project_root)

message("\n")
message("================================================================")
message("  Dataset Merge: ACP_SCN + GSE215932")
message("================================================================")
message(paste("  Started:", Sys.time()))
message(paste("  Working directory:", getwd()))
message(paste("  Integration:", ifelse(do_integration, "Yes", "No")))
message("\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# Report Seurat version
seurat_version <- packageVersion("Seurat")
message(sprintf("  Seurat version: %s", seurat_version))
is_seurat_v5 <- seurat_version >= "5.0.0"
message(sprintf("  Using Seurat v5 API: %s", is_seurat_v5))

source("R/utils/config.R")
config <- load_config(config_path)
set.seed(config$reproducibility$seed)

# ==============================================================================
# LOAD DATASETS
# ==============================================================================

message("========================================")
message("Loading Datasets")
message("========================================\n")

# --- Load ACP_SCN ---
acp_path <- if (!is.null(config$paths$acp_scn_annotated)) {
  get_path(config, config$paths$acp_scn_annotated)
} else {
  get_path(config, "data/raw/acp_scn_annotated.rds")
}

message(sprintf("Loading ACP_SCN: %s", acp_path))
acp_obj <- readRDS(acp_path)
acp_obj <- UpdateSeuratObject(acp_obj)

message(sprintf("  ACP_SCN - Cells: %d, Default assay genes: %d",
                ncol(acp_obj), nrow(acp_obj)))

# --- Load GSE/snRNA-seq ---
gse_path <- get_path(config, config$paths$snrnaseq_processed)

message(sprintf("Loading GSE215932: %s", gse_path))
gse_obj <- readRDS(gse_path)
gse_obj <- UpdateSeuratObject(gse_obj)

message(sprintf("  GSE215932 - Cells: %d, Default assay genes: %d",
                ncol(gse_obj), nrow(gse_obj)))

# ==============================================================================
# EXTRACT RAW COUNTS FROM ACP_SCN
# ==============================================================================

message("\n========================================")
message("Extracting Raw Counts from ACP_SCN")
message("========================================\n")

# The ACP object may have counts in various locations
# Priority: RNA assay counts > SCT counts > integrated counts

acp_counts <- NULL
acp_counts_source <- NULL

# Check available assays
message("Available assays in ACP_SCN:")
for (assay_name in Assays(acp_obj)) {
  assay_obj <- acp_obj[[assay_name]]
  if (inherits(assay_obj, "Assay5")) {
    layers <- Layers(assay_obj)
    message(sprintf("  %s (Assay5): %d features, layers: %s",
                    assay_name, nrow(assay_obj), paste(layers, collapse = ", ")))
  } else {
    message(sprintf("  %s (Assay): %d features", assay_name, nrow(assay_obj)))
  }
}

# Try to get counts from RNA assay first
if ("RNA" %in% Assays(acp_obj)) {
  rna_assay <- acp_obj[["RNA"]]
  rna_class <- class(rna_assay)[1]
  message(sprintf("  RNA assay class: %s", rna_class))

  # In Seurat v5, use LayerData for all assay types
  # Get available layers
  available_layers <- tryCatch({
    Layers(rna_assay)
  }, error = function(e) {
    # Fallback for older objects
    c("counts", "data")
  })
  message(sprintf("    Available layers: %s", paste(available_layers, collapse = ", ")))

  # Try counts layer first
  if ("counts" %in% available_layers) {
    counts_data <- tryCatch({
      LayerData(acp_obj, assay = "RNA", layer = "counts")
    }, error = function(e) {
      message(sprintf("    Error accessing counts layer: %s", e$message))
      NULL
    })

    if (!is.null(counts_data)) {
      counts_dims <- dim(counts_data)
      message(sprintf("    Counts layer dimensions: %d x %d", counts_dims[1], counts_dims[2]))

      if (counts_dims[1] > 0 && counts_dims[2] > 0) {
        acp_counts <- counts_data
        acp_counts_source <- "RNA assay counts layer"
        message(sprintf("    Successfully extracted %d genes from counts layer", nrow(acp_counts)))
      }
    }
  }

  # If counts not found or empty, try data layer
  if (is.null(acp_counts) || nrow(acp_counts) < 5000) {
    if ("data" %in% available_layers) {
      message("    Trying data layer...")
      data_layer <- tryCatch({
        LayerData(acp_obj, assay = "RNA", layer = "data")
      }, error = function(e) {
        message(sprintf("    Error accessing data layer: %s", e$message))
        NULL
      })

      if (!is.null(data_layer) && nrow(data_layer) > 0) {
        data_dims <- dim(data_layer)
        message(sprintf("    Data layer dimensions: %d x %d", data_dims[1], data_dims[2]))

        if (is.null(acp_counts) || data_dims[1] > nrow(acp_counts)) {
          acp_counts <- data_layer
          acp_counts_source <- "RNA assay data layer (log-normalized)"
          message(sprintf("    Using data layer with %d genes", nrow(acp_counts)))
        }
      }
    }
  }
}

# If RNA counts not found or too few genes, try SCT
if (is.null(acp_counts) || nrow(acp_counts) < 5000) {
  message("\n  Checking SCT assay...")
  if ("SCT" %in% Assays(acp_obj)) {
    sct_assay <- acp_obj[["SCT"]]
    sct_class <- class(sct_assay)[1]
    message(sprintf("    SCT assay class: %s", sct_class))

    # Get available layers
    sct_layers <- tryCatch({
      Layers(sct_assay)
    }, error = function(e) {
      c("counts", "data")
    })
    message(sprintf("    Available layers: %s", paste(sct_layers, collapse = ", ")))

    # Try counts layer
    if ("counts" %in% sct_layers) {
      sct_counts <- tryCatch({
        LayerData(acp_obj, assay = "SCT", layer = "counts")
      }, error = function(e) {
        message(sprintf("    Error accessing SCT counts: %s", e$message))
        NULL
      })

      if (!is.null(sct_counts)) {
        sct_dims <- dim(sct_counts)
        message(sprintf("    SCT counts dimensions: %d x %d", sct_dims[1], sct_dims[2]))

        if (sct_dims[1] > 0 && (is.null(acp_counts) || sct_dims[1] > nrow(acp_counts))) {
          acp_counts <- sct_counts
          acp_counts_source <- "SCT assay counts layer"
          message(sprintf("    Using SCT counts: %d genes", nrow(acp_counts)))
        }
      }
    }

    # Try data layer if counts not found
    if ((is.null(acp_counts) || nrow(acp_counts) < 5000) && "data" %in% sct_layers) {
      sct_data <- tryCatch({
        LayerData(acp_obj, assay = "SCT", layer = "data")
      }, error = function(e) NULL)

      if (!is.null(sct_data) && nrow(sct_data) > 0) {
        message(sprintf("    SCT data layer: %d genes", nrow(sct_data)))
        if (is.null(acp_counts) || nrow(sct_data) > nrow(acp_counts)) {
          acp_counts <- sct_data
          acp_counts_source <- "SCT assay data layer"
        }
      }
    }
  }
}

# Final diagnostic before error
if (is.null(acp_counts)) {
  message("\n  ERROR: Could not extract any count data!")
  message("  Detailed assay inspection:")
  for (assay_name in Assays(acp_obj)) {
    assay_obj <- acp_obj[[assay_name]]
    message(sprintf("    %s (class: %s):", assay_name, class(assay_obj)[1]))

    # Try to get layers
    assay_layers <- tryCatch(Layers(assay_obj), error = function(e) NULL)
    if (!is.null(assay_layers)) {
      message(sprintf("      Layers: %s", paste(assay_layers, collapse = ", ")))
      for (layer in assay_layers) {
        layer_data <- tryCatch({
          LayerData(acp_obj, assay = assay_name, layer = layer)
        }, error = function(e) NULL)
        if (!is.null(layer_data)) {
          message(sprintf("      %s: %d x %d", layer, nrow(layer_data), ncol(layer_data)))
        }
      }
    }
  }
  stop("Could not extract raw counts from ACP_SCN object!")
}

message(sprintf("\nExtracted ACP counts from: %s", acp_counts_source))
message(sprintf("  Genes: %d", nrow(acp_counts)))
message(sprintf("  Cells: %d", ncol(acp_counts)))

# Check if this is actually more genes than default
if (nrow(acp_counts) <= 3000) {
  warning("ACP counts matrix has <= 3000 genes - may still be variable features only!")
  message("\nAttempting to find full counts in all possible locations...")

  # Last resort: check all assays for the largest counts matrix
  best_counts <- acp_counts
  best_source <- acp_counts_source

  for (assay_name in Assays(acp_obj)) {
    assay_obj <- acp_obj[[assay_name]]

    # Try to get counts layer using Seurat v5 API
    test_counts <- tryCatch({
      assay_layers <- Layers(assay_obj)
      if ("counts" %in% assay_layers) {
        LayerData(acp_obj, assay = assay_name, layer = "counts")
      } else if ("data" %in% assay_layers) {
        LayerData(acp_obj, assay = assay_name, layer = "data")
      } else {
        NULL
      }
    }, error = function(e) NULL)

    if (!is.null(test_counts) && nrow(test_counts) > nrow(best_counts)) {
      best_counts <- test_counts
      best_source <- sprintf("%s assay", assay_name)
    }
  }

  if (nrow(best_counts) > nrow(acp_counts)) {
    acp_counts <- best_counts
    acp_counts_source <- best_source
    message(sprintf("  Found better source: %s with %d genes", best_source, nrow(best_counts)))
  } else {
    message("\n  WARNING: No full counts matrix found in ACP_SCN object.")
    message("  The merged dataset will be limited to available genes.")
    message("  For better results, provide ACP data before variable feature filtering.")
  }
}

# ==============================================================================
# EXTRACT RAW COUNTS FROM GSE
# ==============================================================================

message("\n========================================")
message("Extracting Raw Counts from GSE215932")
message("========================================\n")

gse_counts <- NULL
gse_counts_source <- NULL

# GSE should have full counts in RNA assay
if ("RNA" %in% Assays(gse_obj)) {
  rna_assay <- gse_obj[["RNA"]]
  rna_class <- class(rna_assay)[1]
  message(sprintf("  GSE RNA assay class: %s", rna_class))

  # Get available layers
  available_layers <- tryCatch({
    Layers(rna_assay)
  }, error = function(e) {
    c("counts", "data")
  })
  message(sprintf("  Available layers: %s", paste(available_layers, collapse = ", ")))

  # Try counts layer
  if ("counts" %in% available_layers) {
    gse_counts <- tryCatch({
      LayerData(gse_obj, assay = "RNA", layer = "counts")
    }, error = function(e) {
      message(sprintf("  Error accessing counts: %s", e$message))
      NULL
    })

    if (!is.null(gse_counts)) {
      gse_counts_source <- "RNA assay counts layer"
      message(sprintf("  Successfully extracted %d genes from counts layer", nrow(gse_counts)))
    }
  }

  # Fallback to data layer if counts empty
  if (is.null(gse_counts) && "data" %in% available_layers) {
    gse_counts <- tryCatch({
      LayerData(gse_obj, assay = "RNA", layer = "data")
    }, error = function(e) NULL)

    if (!is.null(gse_counts)) {
      gse_counts_source <- "RNA assay data layer (normalized)"
    }
  }
}

if (is.null(gse_counts)) {
  stop("Could not extract raw counts from GSE215932 object!")
}

message(sprintf("Extracted GSE counts from: %s", gse_counts_source))
message(sprintf("  Genes: %d", nrow(gse_counts)))
message(sprintf("  Cells: %d", ncol(gse_counts)))

# ==============================================================================
# HARMONIZE METADATA
# ==============================================================================

message("\n========================================")
message("Harmonizing Metadata")
message("========================================\n")

# --- ACP metadata ---
acp_meta <- acp_obj@meta.data

# Standardize column names
if ("sample_id" %in% colnames(acp_meta) && !"sample" %in% colnames(acp_meta)) {
  acp_meta$sample <- acp_meta$sample_id
}

# Add dataset identifier
acp_meta$dataset <- "ACP_SCN"
acp_meta$dataset_source <- "ACP scRNA-seq"

# Standardize epithelial column
if ("celltypes" %in% colnames(acp_meta)) {
  acp_meta$is_epithelial <- acp_meta$celltypes == "Epithelial"
  acp_meta$original_celltype <- acp_meta$celltypes
}

message("ACP_SCN metadata:")
message(sprintf("  Samples: %s", paste(unique(acp_meta$sample), collapse = ", ")))
message(sprintf("  Epithelial cells: %d", sum(acp_meta$is_epithelial, na.rm = TRUE)))

# --- GSE metadata ---
gse_meta <- gse_obj@meta.data

# Add dataset identifier
gse_meta$dataset <- "GSE215932"
gse_meta$dataset_source <- "GSE215932 snRNA-seq"

# Standardize epithelial - should already have is_epithelial
if (!"is_epithelial" %in% colnames(gse_meta)) {
  # Infer from cell_type if available
  if ("cell_type" %in% colnames(gse_meta)) {
    gse_meta$is_epithelial <- grepl("^E_", gse_meta$cell_type)
  }
}

# Keep original cell type annotations
if ("cell_type" %in% colnames(gse_meta)) {
  gse_meta$original_celltype <- gse_meta$cell_type
}

message("\nGSE215932 metadata:")
message(sprintf("  Samples: %s", paste(unique(gse_meta$sample), collapse = ", ")))
message(sprintf("  Epithelial cells: %d", sum(gse_meta$is_epithelial, na.rm = TRUE)))

# --- Find common metadata columns ---
common_cols <- intersect(colnames(acp_meta), colnames(gse_meta))
message(sprintf("\nCommon metadata columns: %s", paste(common_cols, collapse = ", ")))

# Essential columns to keep
essential_cols <- c("sample", "dataset", "dataset_source", "is_epithelial",
                    "original_celltype", "nCount_RNA", "nFeature_RNA")
essential_cols <- essential_cols[essential_cols %in% common_cols |
                                   essential_cols %in% colnames(acp_meta) |
                                   essential_cols %in% colnames(gse_meta)]

# ==============================================================================
# FIND COMMON GENES AND MERGE
# ==============================================================================

message("\n========================================")
message("Finding Common Genes")
message("========================================\n")

acp_genes <- rownames(acp_counts)
gse_genes <- rownames(gse_counts)

common_genes <- intersect(acp_genes, gse_genes)
acp_only <- setdiff(acp_genes, gse_genes)
gse_only <- setdiff(gse_genes, acp_genes)

message(sprintf("ACP genes: %d", length(acp_genes)))
message(sprintf("GSE genes: %d", length(gse_genes)))
message(sprintf("Common genes: %d", length(common_genes)))
message(sprintf("ACP-only genes: %d", length(acp_only)))
message(sprintf("GSE-only genes: %d", length(gse_only)))

if (length(common_genes) < 5000) {
  warning(sprintf("Only %d common genes - this may limit analysis!", length(common_genes)))
}

# Check signature gene coverage in common genes
message("\n--- Signature Gene Coverage in Common Genes ---")
all_sig_genes <- unique(unlist(config$signatures$epithelial_subtypes))
sig_in_common <- sum(all_sig_genes %in% common_genes)
message(sprintf("Signature genes in common set: %d/%d (%.1f%%)",
                sig_in_common, length(all_sig_genes),
                100 * sig_in_common / length(all_sig_genes)))

# Subset to common genes
acp_counts_common <- acp_counts[common_genes, ]
gse_counts_common <- gse_counts[common_genes, ]

message(sprintf("\nSubsetted matrices to %d common genes", length(common_genes)))

# ==============================================================================
# CREATE MERGED SEURAT OBJECT
# ==============================================================================

message("\n========================================")
message("Creating Merged Seurat Object")
message("========================================\n")
# Create separate Seurat objects with harmonized data
message("Creating ACP Seurat object...")
acp_new <- CreateSeuratObject(
  counts = acp_counts_common,
  meta.data = acp_meta[colnames(acp_counts_common), , drop = FALSE],
  project = "ACP_SCN"
)

message("Creating GSE Seurat object...")
gse_new <- CreateSeuratObject(
  counts = gse_counts_common,
  meta.data = gse_meta[colnames(gse_counts_common), , drop = FALSE],
  project = "GSE215932"
)

# Merge
message("Merging objects...")
merged_obj <- merge(
  x = acp_new,
  y = gse_new,
  add.cell.ids = c("ACP", "GSE"),
  project = "ACP_GSE_merged"
)

message(sprintf("\nMerged object:"))
message(sprintf("  Total cells: %d", ncol(merged_obj)))
message(sprintf("  Genes: %d", nrow(merged_obj)))
message(sprintf("  ACP cells: %d", sum(merged_obj$dataset == "ACP_SCN")))
message(sprintf("  GSE cells: %d", sum(merged_obj$dataset == "GSE215932")))

# Dataset breakdown
message("\nDataset breakdown:")
print(table(merged_obj$dataset))

message("\nEpithelial cells by dataset:")
print(table(merged_obj$dataset, merged_obj$is_epithelial))

# ==============================================================================
# NORMALIZE DATA
# ==============================================================================

message("\n========================================")
message("Normalizing Merged Data")
message("========================================\n")

merged_obj <- NormalizeData(merged_obj, verbose = FALSE)
message("Log-normalization complete")

merged_obj <- FindVariableFeatures(merged_obj, nfeatures = 3000, verbose = FALSE)
message("Found 3000 variable features")

merged_obj <- ScaleData(merged_obj, verbose = FALSE)
message("Scaling complete")

# ==============================================================================
# OPTIONAL: INTEGRATION (BATCH CORRECTION)
# ==============================================================================

if (do_integration) {
  message("\n========================================")
  message("Integrating Datasets (Batch Correction)")
  message("========================================\n")

  # Split by dataset
  obj_list <- SplitObject(merged_obj, split.by = "dataset")

  # Find integration anchors
  message("Finding integration anchors...")
  anchors <- FindIntegrationAnchors(
    object.list = obj_list,
    dims = 1:30,
    verbose = TRUE
  )

  # Integrate
  message("Integrating data...")
  merged_obj <- IntegrateData(
    anchorset = anchors,
    dims = 1:30,
    verbose = TRUE
  )

  # Scale integrated data
  DefaultAssay(merged_obj) <- "integrated"
  merged_obj <- ScaleData(merged_obj, verbose = FALSE)

  message("Integration complete")
  message(sprintf("  Default assay: %s", DefaultAssay(merged_obj)))
}

# ==============================================================================
# RUN PCA AND UMAP
# ==============================================================================

message("\n========================================")
message("Dimensionality Reduction")
message("========================================\n")

merged_obj <- RunPCA(merged_obj, npcs = 50, verbose = FALSE)
message("PCA complete (50 PCs)")

merged_obj <- RunUMAP(merged_obj, dims = 1:30, verbose = FALSE)
message("UMAP complete")

# ==============================================================================
# SAVE OUTPUT
# ==============================================================================

message("\n========================================")
message("Saving Output")
message("========================================\n")

# Create output directory
output_dir <- get_path(config, "data/processed")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save merged object
output_filename <- ifelse(do_integration,
                          "merged_acp_gse_integrated.rds",
                          "merged_acp_gse.rds")
output_path <- file.path(output_dir, output_filename)

saveRDS(merged_obj, output_path)
message(sprintf("Saved: %s", output_path))

# Save summary
summary_df <- data.frame(
  metric = c("Total cells", "Total genes", "ACP cells", "GSE cells",
             "ACP epithelial", "GSE epithelial", "Common genes",
             "Signature gene coverage"),
  value = c(ncol(merged_obj), nrow(merged_obj),
            sum(merged_obj$dataset == "ACP_SCN"),
            sum(merged_obj$dataset == "GSE215932"),
            sum(merged_obj$dataset == "ACP_SCN" & merged_obj$is_epithelial, na.rm = TRUE),
            sum(merged_obj$dataset == "GSE215932" & merged_obj$is_epithelial, na.rm = TRUE),
            length(common_genes),
            sprintf("%d/%d (%.1f%%)", sig_in_common, length(all_sig_genes),
                    100 * sig_in_common / length(all_sig_genes)))
)

summary_path <- file.path(output_dir, "merged_dataset_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)
message(sprintf("Saved summary: %s", summary_path))

# ==============================================================================
# COMPLETION
# ==============================================================================

message("\n")
message("================================================================")
message("  Dataset Merge Complete!")
message("================================================================")
message(sprintf("  Finished: %s", Sys.time()))
message(sprintf("  Total cells: %d", ncol(merged_obj)))
message(sprintf("  Total genes: %d", nrow(merged_obj)))
message(sprintf("  Integration: %s", ifelse(do_integration, "Yes", "No")))
message("\nTo run classification on merged data:")
message("  Rscript scripts/01_cell_type_annotation.R --dataset merged")
message("================================================================\n")
