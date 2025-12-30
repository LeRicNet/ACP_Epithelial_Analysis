#!/usr/bin/env Rscript
# ==============================================================================
# scripts/01_cell_type_annotation.R
# ==============================================================================
# Cell Type Annotation for ACP Spatial Transcriptomics or scRNA-seq Data
#
# This script performs module score-based classification of epithelial subtypes,
# with comprehensive per-sample analysis and statistical testing.
#
# Key features:
#   - Multi-dataset support (spatial, snrnaseq, acp_scn)
#   - Full gene data recovery for variable-feature filtered objects
#   - Per-sample statistical analysis with chi-square, Kruskal-Wallis tests
#   - Sample consistency metrics (Jensen-Shannon divergence)
#   - Comprehensive diagnostic visualizations
#
# Usage:
#   Rscript scripts/01_cell_type_annotation.R                          # Default: spatial
#   Rscript scripts/01_cell_type_annotation.R --dataset spatial        # Visium data
#   Rscript scripts/01_cell_type_annotation.R --dataset snrnaseq       # snRNA-seq GSE data
#   Rscript scripts/01_cell_type_annotation.R --dataset acp_scn        # ACP scRNA-seq data
#
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
dataset_type <- "spatial"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--config" && i < length(args)) {
    config_path <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--dataset" && i < length(args)) {
    dataset_type <- tolower(args[i + 1])
    i <- i + 2
  } else if (!startsWith(args[i], "--")) {
    config_path <- args[i]
    i <- i + 1
  } else {
    i <- i + 1
  }
}

valid_datasets <- c("spatial", "snrnaseq", "acp_scn", "merged")
if (!dataset_type %in% valid_datasets) {
  stop(sprintf("Invalid --dataset value. Use one of: %s", paste(valid_datasets, collapse = ", ")))
}

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
message("  Step 1: Cell Type Annotation")
message(sprintf("  Dataset type: %s", toupper(dataset_type)))
message("================================================================")
message(paste("  Started:", Sys.time()))
message(paste("  Working directory:", getwd()))
message("\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

source("R/utils/config.R")
config <- load_config(config_path)
print_config_summary(config)
set.seed(config$reproducibility$seed)

# ==============================================================================
# DATASET-SPECIFIC CONFIGURATION
# ==============================================================================

check_full_data <- FALSE  # Flag for ACP data handling

if (dataset_type == "spatial") {
  input_path <- get_path(config, config$paths$spatial_object)
  alt_path <- get_path(config, config$paths$spatial_object_with_subtypes)
  count_col <- "nCount_Spatial"
  feature_col <- "nFeature_Spatial"
  default_assay <- "Spatial"
  epi_column <- "cellstates"
  epi_values <- "Epithelial"
  epi_is_boolean <- FALSE
  message("Configured for 10x Visium spatial data")

} else if (dataset_type == "snrnaseq") {
  input_path <- get_path(config, config$paths$snrnaseq_processed)
  alt_path <- NULL
  count_col <- "nCount_RNA"
  feature_col <- "nFeature_RNA"
  default_assay <- "RNA"
  epi_column <- "is_epithelial"
  epi_values <- TRUE
  epi_is_boolean <- TRUE
  message("Configured for 10x snRNA-seq data (GSE215932)")

} else if (dataset_type == "acp_scn") {
  if (!is.null(config$paths$acp_scn_annotated)) {
    input_path <- get_path(config, config$paths$acp_scn_annotated)
  } else {
    input_path <- get_path(config, "data/raw/acp_scn_annotated.rds")
  }
  alt_path <- NULL
  count_col <- "nCount_RNA"
  feature_col <- "nFeature_RNA"
  default_assay <- "RNA"
  epi_column <- "celltypes"
  epi_values <- "Epithelial"
  epi_is_boolean <- FALSE
  check_full_data <- TRUE  # Enable full data check for ACP
  message("Configured for ACP scRNA-seq annotated data")

} else if (dataset_type == "merged") {
  # Merged ACP + GSE dataset
  if (!is.null(config$paths$merged_dataset)) {
    input_path <- get_path(config, config$paths$merged_dataset)
  } else {
    # Check for both integrated and non-integrated versions
    integrated_path <- get_path(config, "data/processed/merged_acp_gse_integrated.rds")
    standard_path <- get_path(config, "data/processed/merged_acp_gse.rds")

    if (file.exists(integrated_path)) {
      input_path <- integrated_path
      message("  Using integrated merged dataset")
    } else {
      input_path <- standard_path
    }
  }
  alt_path <- NULL
  count_col <- "nCount_RNA"
  feature_col <- "nFeature_RNA"
  default_assay <- "RNA"

  # Epithelial identification - merged object has standardized is_epithelial
  epi_column <- "is_epithelial"
  epi_values <- TRUE
  epi_is_boolean <- TRUE

  # Don't check for full data - merged object already has full counts
  check_full_data <- FALSE

  message("Configured for merged ACP + GSE dataset")
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("\n========================================")
message("Loading Data")
message("========================================\n")

if (!is.null(alt_path) && file.exists(alt_path)) {
  message("Found object with existing classifications - will compare methods")
  input_path <- alt_path
}

if (!file.exists(input_path)) {
  stop(paste("Input file not found:", input_path))
}

message(paste("Loading:", input_path))
seurat_obj <- readRDS(input_path)

if (inherits(seurat_obj, "Seurat")) {
  seurat_obj <- tryCatch({
    UpdateSeuratObject(seurat_obj)
  }, error = function(e) seurat_obj)
}

message(sprintf("  Cells: %d", ncol(seurat_obj)))
message(sprintf("  Genes: %d", nrow(seurat_obj)))

# ==============================================================================
# CHECK FOR FULL GENE DATA (ACP_SCN objects may be variable-feature filtered)
# ==============================================================================

if (check_full_data) {
  message("\n=== Checking Data Structure for Full Gene Coverage ===")

  current_genes <- nrow(seurat_obj)
  message(sprintf("  Current genes in default assay: %d", current_genes))

  if (current_genes < 5000) {
    message("  Warning: Low gene count suggests variable features only")
    full_data_found <- FALSE

    # Check Seurat v5 counts layer
    if ("RNA" %in% Assays(seurat_obj) && inherits(seurat_obj[["RNA"]], "Assay5")) {
      available_layers <- Layers(seurat_obj[["RNA"]])
      message(sprintf("  Seurat v5 - available layers: %s", paste(available_layers, collapse = ", ")))

      if ("counts" %in% available_layers) {
        counts_genes <- nrow(LayerData(seurat_obj, assay = "RNA", layer = "counts"))
        message(sprintf("  Genes in counts layer: %d", counts_genes))

        if (counts_genes > current_genes) {
          message("  → Recovering full gene data from counts layer")
          full_counts <- LayerData(seurat_obj, assay = "RNA", layer = "counts")
          seurat_obj[["RNA_full"]] <- CreateAssay5Object(counts = full_counts)
          DefaultAssay(seurat_obj) <- "RNA_full"
          seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
          default_assay <- "RNA_full"
          full_data_found <- TRUE
          message(sprintf("  Created RNA_full assay with %d genes", nrow(seurat_obj)))
        }
      }
    }

    # Check Seurat v4 counts slot
    if (!full_data_found && "RNA" %in% Assays(seurat_obj) && !inherits(seurat_obj[["RNA"]], "Assay5")) {
      counts_dim <- tryCatch(dim(GetAssayData(seurat_obj, slot = "counts", assay = "RNA")),
                             error = function(e) c(0, 0))
      if (counts_dim[1] > current_genes) {
        message(sprintf("  Found %d genes in counts slot", counts_dim[1]))
        message("  → Re-normalizing from full counts...")
        full_counts <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")
        seurat_obj <- SetAssayData(seurat_obj, slot = "data",
                                   new.data = Seurat::LogNormalize(full_counts))
        full_data_found <- TRUE
        message(sprintf("  Now using %d genes", nrow(seurat_obj)))
      }
    }

    if (!full_data_found) {
      message("\n  WARNING: Could not find full gene data.")
      message("  Module scoring will use limited gene set. Results may be affected.")
    }
  }

  message(sprintf("\n  Final: Assay=%s, Genes=%d", DefaultAssay(seurat_obj), nrow(seurat_obj)))
}

# Auto-detect assay type
if (!count_col %in% colnames(seurat_obj@meta.data)) {
  if ("nCount_Spatial" %in% colnames(seurat_obj@meta.data)) {
    count_col <- "nCount_Spatial"
    feature_col <- "nFeature_Spatial"
    default_assay <- "Spatial"
  } else if ("nCount_RNA" %in% colnames(seurat_obj@meta.data)) {
    count_col <- "nCount_RNA"
    feature_col <- "nFeature_RNA"
    default_assay <- "RNA"
  }
}

# ==============================================================================
# IDENTIFY SAMPLE COLUMN
# ==============================================================================

sample_col_candidates <- c("sample", "Sample", "sample_id", "Sample_ID",
                           "orig.ident", "patient", "Patient", "donor", "Donor")
sample_col <- NULL
for (col in sample_col_candidates) {
  if (col %in% colnames(seurat_obj@meta.data)) {
    sample_col <- col
    break
  }
}

if (!is.null(sample_col)) {
  message(sprintf("\n  Sample column identified: '%s'", sample_col))
  message(sprintf("  Samples: %s", paste(unique(seurat_obj@meta.data[[sample_col]]), collapse = ", ")))
} else {
  message("\n  No sample column found - per-sample analysis will be skipped")
}

# ==============================================================================
# EPITHELIAL CELL IDENTIFICATION
# ==============================================================================

message("\n========================================")
message("Epithelial Cell Identification")
message("========================================\n")

if (!epi_column %in% colnames(seurat_obj@meta.data)) {
  stop(paste0("Epithelial column '", epi_column, "' not found in metadata."))
}

if (epi_is_boolean) {
  epithelial_mask <- seurat_obj@meta.data[[epi_column]] == TRUE
  message(sprintf("  Using boolean column '%s'", epi_column))
} else {
  epithelial_mask <- seurat_obj@meta.data[[epi_column]] %in% epi_values
  message(sprintf("  Using character column '%s', matching: %s", epi_column, paste(epi_values, collapse = ", ")))
}

epithelial_mask[is.na(epithelial_mask)] <- FALSE

n_total <- ncol(seurat_obj)
n_epithelial <- sum(epithelial_mask)

message(sprintf("  Total cells: %d", n_total))
message(sprintf("  Epithelial cells: %d (%.1f%%)", n_epithelial, 100 * n_epithelial / n_total))

# Show per-sample breakdown of epithelial identification
if (!is.null(sample_col)) {
  message("\n  Per-sample epithelial cell counts:")
  epi_by_sample <- table(
    Sample = seurat_obj@meta.data[[sample_col]],
    Epithelial = epithelial_mask
  )
  print(epi_by_sample)
}

message(sprintf("\n  Original '%s' distribution:", epi_column))
print(table(seurat_obj@meta.data[[epi_column]], useNA = "ifany"))

seurat_obj$is_epithelial_cell <- epithelial_mask

if (n_epithelial == 0) {
  stop("No epithelial cells found!")
}

message(sprintf("\n  Subsetting to %d epithelial cells...", n_epithelial))
seurat_obj_full <- seurat_obj
seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[epithelial_mask])

message(sprintf("  After filtering - Cells: %d, Genes: %d", ncol(seurat_obj), nrow(seurat_obj)))

# Report existing annotations
if ("epithelial_subtype" %in% colnames(seurat_obj@meta.data)) {
  message("\nExisting classification (epithelial_subtype):")
  print(table(seurat_obj$epithelial_subtype))
}

if (dataset_type == "snrnaseq" && "cell_type" %in% colnames(seurat_obj@meta.data)) {
  message("\nExisting cell_type annotations (from GSE data):")
  print(table(seurat_obj$cell_type, useNA = "ifany"))
}

if (!is.null(sample_col) && sample_col %in% colnames(seurat_obj@meta.data)) {
  message("\nEpithelial cells per sample:")
  print(table(seurat_obj@meta.data[[sample_col]]))
}

# ==============================================================================
# CALCULATE MODULE SCORES
# ==============================================================================

message("\n========================================")
message("Calculating Module Scores")
message("========================================\n")

available_genes <- rownames(seurat_obj)
n_genes <- length(available_genes)
message(sprintf("Genes available in dataset: %d", n_genes))

# Gene signature coverage report
message("\n=== Gene Signature Coverage Report ===")
all_signature_genes <- unique(unlist(config$signatures$epithelial_subtypes))
n_sig_genes <- length(all_signature_genes)
n_present <- sum(all_signature_genes %in% available_genes)

message(sprintf("Total unique genes in signatures: %d", n_sig_genes))
message(sprintf("Genes present in data: %d (%.1f%%)", n_present, 100 * n_present / n_sig_genes))

coverage_pct <- n_present / n_sig_genes * 100
if (coverage_pct < 50) {
  warning(sprintf("LOW GENE COVERAGE (%.1f%%)! Results may be unreliable.", coverage_pct))
} else if (coverage_pct < 80) {
  message(sprintf("Note: Moderate gene coverage (%.1f%%).", coverage_pct))
}

optimal_nbin <- min(24, max(5, floor(n_genes / 30)))
message(sprintf("\nUsing nbin = %d for module scoring", optimal_nbin))

# Prepare signatures
signatures <- list()
message("\nPreparing epithelial subtype signatures:")
for (sig_name in names(config$signatures$epithelial_subtypes)) {
  sig_genes <- config$signatures$epithelial_subtypes[[sig_name]]
  present_genes <- intersect(sig_genes, available_genes)
  message(sprintf("  %s: %d/%d genes present", sig_name, length(present_genes), length(sig_genes)))
  if (length(present_genes) >= 3) {
    signatures[[sig_name]] <- present_genes
  } else {
    warning(sprintf("  Skipping %s: insufficient genes", sig_name))
  }
}

# Additional signatures
message("\nPreparing additional signatures:")
additional_sigs <- c("Whorl_Wnt", "Senescence", "SASP")
for (sig_name in additional_sigs) {
  if (sig_name %in% names(config$signatures)) {
    sig_genes <- config$signatures[[sig_name]]
    present_genes <- intersect(sig_genes, available_genes)
    message(sprintf("  %s: %d/%d genes present", sig_name, length(present_genes), length(sig_genes)))
    if (length(present_genes) >= 3) {
      signatures[[sig_name]] <- present_genes
    }
  }
}

# Robust module scoring function
score_module_robust <- function(seurat_obj, features, name, seed, nbin_start = 24) {
  nbin_values <- c(nbin_start, 12, 8, 5)
  nbin_values <- nbin_values[nbin_values <= nbin_start]

  for (nbin in nbin_values) {
    result <- tryCatch({
      obj <- AddModuleScore(seurat_obj, features = list(features), name = name, nbin = nbin, seed = seed)
      list(success = TRUE, obj = obj, nbin = nbin)
    }, error = function(e) list(success = FALSE, error = e$message))
    if (result$success) return(result)
  }
  return(list(success = FALSE, error = "All nbin values failed"))
}

# Calculate scores
message("\nCalculating module scores...")
for (sig_name in names(signatures)) {
  message(sprintf("  Scoring: %s", sig_name))
  result <- score_module_robust(seurat_obj, features = signatures[[sig_name]],
                                name = paste0("score_", sig_name),
                                seed = config$reproducibility$seed, nbin_start = optimal_nbin)
  if (result$success) {
    seurat_obj <- result$obj
    old_name <- paste0("score_", sig_name, "1")
    new_name <- paste0("score_", sig_name)
    if (old_name %in% colnames(seurat_obj@meta.data)) {
      seurat_obj@meta.data[[new_name]] <- seurat_obj@meta.data[[old_name]]
      seurat_obj@meta.data[[old_name]] <- NULL
    }
  } else {
    warning(sprintf("  Failed to score %s", sig_name))
  }
}

# Report score distributions
message("\nModule score summary statistics:")
score_cols <- grep("^score_", colnames(seurat_obj@meta.data), value = TRUE)

for (col in score_cols) {
  scores <- seurat_obj@meta.data[[col]]
  message(sprintf("  %s: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]",
                  gsub("^score_", "", col), mean(scores, na.rm = TRUE),
                  sd(scores, na.rm = TRUE), min(scores, na.rm = TRUE), max(scores, na.rm = TRUE)))
}

# ==============================================================================
# CELL CYCLE SCORING
# ==============================================================================

message("\n========================================")
message("Cell Cycle Scoring")
message("========================================\n")

source("R/utils/cell_cycle_scoring.R")

if ("cc_consensus" %in% colnames(seurat_obj@meta.data)) {
  message("Cell cycle consensus already present")
  print(table(seurat_obj$cc_consensus))
} else if ("Phase" %in% colnames(seurat_obj@meta.data)) {
  message("Single-method cell cycle phase already present")
  print(table(seurat_obj$Phase))
} else {
  message("Scoring cell cycle phases...")
  cc_methods <- if (!is.null(config$cell_cycle$methods)) config$cell_cycle$methods else c("seurat", "module_score")
  seurat_obj <- score_cell_cycle_consensus(seurat_obj, config = config, methods = cc_methods,
                                           seed = config$reproducibility$seed)
}

# ==============================================================================
# CLASSIFY CELLS BY MODULE SCORE
# ==============================================================================

message("\n========================================")
message("Classifying Cells")
message("========================================\n")

min_threshold <- if (!is.null(config$classification$min_score_threshold)) {
  config$classification$min_score_threshold
} else 0.1
message(sprintf("Minimum score threshold: %.2f", min_threshold))

core_signatures <- c("Basal_like", "Intermediate", "Specialized")
core_score_cols <- paste0("score_", core_signatures)
core_score_cols <- core_score_cols[core_score_cols %in% colnames(seurat_obj@meta.data)]

if (length(core_score_cols) == 0) {
  stop("No core epithelial signatures were scored.")
}

score_matrix <- seurat_obj@meta.data[, core_score_cols, drop = FALSE]
colnames(score_matrix) <- gsub("^score_", "", colnames(score_matrix))

max_scores <- apply(score_matrix, 1, max, na.rm = TRUE)
max_signature <- apply(score_matrix, 1, function(x) {
  if (all(is.na(x))) return(NA)
  names(x)[which.max(x)]
})

classification <- ifelse(max_scores >= min_threshold, max_signature, "Unassigned")
message(sprintf("\nInitial classification (threshold = %.2f):", min_threshold))
print(table(classification))

# ==============================================================================
# IDENTIFY TRANSIT-AMPLIFYING CELLS
# ==============================================================================

message("\n--- Identifying Transit-Amplifying Cells ---")

ta_config <- config$classification$transit_amplifying
prolif_threshold <- if (!is.null(ta_config$proliferation_threshold)) ta_config$proliferation_threshold else 0.15
override_threshold <- if (!is.null(ta_config$override_threshold)) ta_config$override_threshold else 0.3
require_cycling <- if (!is.null(ta_config$require_cycling)) ta_config$require_cycling else TRUE

message(sprintf("  Proliferation threshold: %.2f", prolif_threshold))
message(sprintf("  Override threshold: %.2f", override_threshold))
message(sprintf("  Require cycling: %s", require_cycling))

if ("score_Transit_Amplifying" %in% colnames(seurat_obj@meta.data)) {
  prolif_scores <- seurat_obj@meta.data$score_Transit_Amplifying
} else {
  prolif_markers <- config$markers$proliferation
  prolif_present <- intersect(prolif_markers, available_genes)
  if (length(prolif_present) >= 2) {
    result <- score_module_robust(seurat_obj, features = prolif_present,
                                  name = "score_Proliferation", seed = config$reproducibility$seed,
                                  nbin_start = optimal_nbin)
    if (result$success) {
      seurat_obj <- result$obj
      seurat_obj@meta.data$score_Proliferation <- seurat_obj@meta.data$score_Proliferation1
      seurat_obj@meta.data$score_Proliferation1 <- NULL
      prolif_scores <- seurat_obj@meta.data$score_Proliferation
    } else {
      prolif_scores <- rep(0, ncol(seurat_obj))
    }
  } else {
    prolif_scores <- rep(0, ncol(seurat_obj))
  }
}

if (require_cycling && "Phase" %in% colnames(seurat_obj@meta.data)) {
  is_cycling <- seurat_obj$Phase %in% c("S", "G2M")
  message(sprintf("  Cycling cells (S/G2M): %d (%.1f%%)", sum(is_cycling), mean(is_cycling) * 100))
} else {
  is_cycling <- rep(TRUE, ncol(seurat_obj))
}

if (!is.null(ta_config$override_on_high_proliferation) && ta_config$override_on_high_proliferation) {
  high_prolif <- prolif_scores >= override_threshold & is_cycling
  n_override <- sum(high_prolif & classification != "Transit-Amplifying")
  classification[high_prolif] <- "Transit-Amplifying"
  message(sprintf("  High proliferation override: %d cells", n_override))
}

moderate_prolif <- prolif_scores >= prolif_threshold & is_cycling
unassigned_ta <- classification == "Unassigned" & moderate_prolif
classification[unassigned_ta] <- "Transit-Amplifying"
message(sprintf("  Unassigned → Transit-Amplifying: %d cells", sum(unassigned_ta)))

# ==============================================================================
# DETECT HYBRID STATES
# ==============================================================================

message("\n--- Detecting Hybrid States ---")

hybrid_threshold <- min_threshold * 0.8
score_ratio_max <- 0.5

basal_scores <- seurat_obj@meta.data$score_Basal_like
int_scores <- seurat_obj@meta.data$score_Intermediate
spec_scores <- seurat_obj@meta.data$score_Specialized

if (!is.null(basal_scores) && !is.null(int_scores)) {
  both_high <- basal_scores >= hybrid_threshold & int_scores >= hybrid_threshold
  score_diff <- abs(basal_scores - int_scores)
  max_score <- pmax(basal_scores, int_scores)
  close_scores <- (score_diff / max_score) < score_ratio_max
  eligible <- classification %in% c("Basal_like", "Intermediate")
  basal_int_hybrid <- both_high & close_scores & eligible
  classification[basal_int_hybrid] <- "Basal-Intermediate"
  message(sprintf("  Basal-Intermediate hybrids: %d", sum(basal_int_hybrid)))
}

if (!is.null(int_scores) && !is.null(spec_scores)) {
  both_high <- int_scores >= hybrid_threshold & spec_scores >= hybrid_threshold
  score_diff <- abs(int_scores - spec_scores)
  max_score <- pmax(int_scores, spec_scores)
  close_scores <- (score_diff / max_score) < score_ratio_max
  eligible <- classification %in% c("Intermediate", "Specialized")
  int_spec_hybrid <- both_high & close_scores & eligible
  classification[int_spec_hybrid] <- "Intermediate-Specialized"
  message(sprintf("  Intermediate-Specialized hybrids: %d", sum(int_spec_hybrid)))
}

# ==============================================================================
# FINALIZE CLASSIFICATION
# ==============================================================================

message("\n--- Finalizing Classification ---")

classification <- gsub("Basal_like", "Basal-like", classification)
seurat_obj$module_score_subtype <- classification

subtype_levels <- c("Basal-like", "Transit-Amplifying", "Intermediate", "Specialized",
                    "Basal-Intermediate", "Intermediate-Specialized", "Unassigned")
subtype_levels <- subtype_levels[subtype_levels %in% unique(classification)]
seurat_obj$module_score_subtype <- factor(seurat_obj$module_score_subtype, levels = subtype_levels)

message("\nFinal classification:")
final_counts <- table(seurat_obj$module_score_subtype)
print(final_counts)

final_pct <- prop.table(final_counts) * 100
message("\nPercentages:")
for (i in seq_along(final_counts)) {
  message(sprintf("  %s: %.1f%%", names(final_counts)[i], final_pct[i]))
}

# ==============================================================================
# PER-SAMPLE STATISTICAL ANALYSIS
# ==============================================================================

message("\n========================================")
message("Per-Sample Statistical Analysis")
message("========================================\n")

per_sample_results <- NULL

if (!is.null(sample_col) && sample_col %in% colnames(seurat_obj@meta.data)) {

  samples <- seurat_obj@meta.data[[sample_col]]
  unique_samples <- unique(samples)
  n_samples <- length(unique_samples)

  message(sprintf("  Number of samples: %d", n_samples))
  message(sprintf("  Samples: %s", paste(unique_samples, collapse = ", ")))

  # -------------------------------------------------------------------------
  # Dataset-Level Analysis (for merged data)
  # -------------------------------------------------------------------------
  if ("dataset" %in% colnames(seurat_obj@meta.data) && dataset_type == "merged") {
    message("\n--- Dataset Source Breakdown (ACP vs GSE) ---")

    dataset_subtype_table <- table(
      Dataset = seurat_obj$dataset,
      Subtype = seurat_obj$module_score_subtype
    )

    message("\nCell counts by dataset source:")
    print(dataset_subtype_table)

    dataset_proportions <- prop.table(dataset_subtype_table, margin = 1) * 100
    message("\nPercentages within each dataset:")
    print(round(dataset_proportions, 1))

    # Chi-square test for dataset differences
    if (nrow(dataset_subtype_table) >= 2 && ncol(dataset_subtype_table) >= 2) {
      dataset_chi <- tryCatch({
        chisq.test(dataset_subtype_table)
      }, warning = function(w) {
        chisq.test(dataset_subtype_table, simulate.p.value = TRUE, B = 10000)
      }, error = function(e) NULL)

      if (!is.null(dataset_chi)) {
        message(sprintf("\nChi-square test (dataset × subtype): X² = %.2f, p = %.2e",
                        dataset_chi$statistic, dataset_chi$p.value))
        if (dataset_chi$p.value < 0.05) {
          message("  → Significant differences between ACP and GSE classification patterns")
        }
      }
    }
  }

  # -------------------------------------------------------------------------
  # Per-Sample Cell Counts and Proportions
  # -------------------------------------------------------------------------
  message("\n--- Per-Sample Classification Breakdown ---")

  sample_subtype_table <- table(
    Sample = seurat_obj@meta.data[[sample_col]],
    Subtype = seurat_obj$module_score_subtype
  )

  message("\nCell counts per sample:")
  print(sample_subtype_table)

  sample_proportions <- prop.table(sample_subtype_table, margin = 1) * 100
  message("\nPercentages within each sample:")
  print(round(sample_proportions, 1))

  # Summary statistics across samples
  message("\n--- Subtype Proportion Summary Across Samples ---")

  subtype_summary <- data.frame(
    Subtype = colnames(sample_proportions),
    Mean_Pct = apply(sample_proportions, 2, mean),
    SD_Pct = apply(sample_proportions, 2, sd),
    Min_Pct = apply(sample_proportions, 2, min),
    Max_Pct = apply(sample_proportions, 2, max),
    CV = apply(sample_proportions, 2, function(x) if (mean(x) > 0) sd(x) / mean(x) * 100 else NA)
  )
  rownames(subtype_summary) <- NULL

  message("\nSubtype proportion statistics (%):")
  subtype_summary_print <- subtype_summary
  subtype_summary_print[, -1] <- round(subtype_summary_print[, -1], 2)
  print(subtype_summary_print)

  # Flag highly variable subtypes
  high_cv_subtypes <- subtype_summary$Subtype[!is.na(subtype_summary$CV) & subtype_summary$CV > 50]
  if (length(high_cv_subtypes) > 0) {
    message(sprintf("\n  Note: High variability (CV > 50%%) in: %s",
                    paste(high_cv_subtypes, collapse = ", ")))
  }

  # -------------------------------------------------------------------------
  # Statistical Tests
  # -------------------------------------------------------------------------
  message("\n--- Statistical Tests ---")

  chi_test <- NULL
  if (n_samples >= 2 && ncol(sample_subtype_table) >= 2) {
    valid_table <- sample_subtype_table[rowSums(sample_subtype_table) > 0,
                                        colSums(sample_subtype_table) > 0, drop = FALSE]

    if (nrow(valid_table) >= 2 && ncol(valid_table) >= 2) {
      chi_test <- tryCatch({
        chisq.test(valid_table)
      }, warning = function(w) {
        chisq.test(valid_table, simulate.p.value = TRUE, B = 10000)
      }, error = function(e) NULL)

      if (!is.null(chi_test)) {
        message("\nChi-square test for independence (sample × subtype):")
        message(sprintf("  X-squared = %.2f, p-value = %.2e", chi_test$statistic, chi_test$p.value))
        if (chi_test$p.value < 0.05) {
          message("  → Significant heterogeneity in subtype proportions across samples (p < 0.05)")
        } else {
          message("  → No significant heterogeneity detected (p >= 0.05)")
        }
      }
    }
  }

  # -------------------------------------------------------------------------
  # Per-Sample Module Score Statistics
  # -------------------------------------------------------------------------
  message("\n--- Per-Sample Module Score Statistics ---")

  sample_score_stats <- do.call(rbind, lapply(unique_samples, function(s) {
    sample_cells <- which(seurat_obj@meta.data[[sample_col]] == s)
    sample_data <- seurat_obj@meta.data[sample_cells, score_cols, drop = FALSE]
    data.frame(
      Sample = s,
      n_cells = length(sample_cells),
      t(colMeans(sample_data, na.rm = TRUE))
    )
  }))

  message("\nMean module scores by sample:")
  sample_score_stats_print <- sample_score_stats
  numeric_cols <- sapply(sample_score_stats_print, is.numeric)
  sample_score_stats_print[, numeric_cols] <- round(sample_score_stats_print[, numeric_cols], 3)
  print(sample_score_stats_print)

  # Kruskal-Wallis tests for score differences
  message("\nKruskal-Wallis tests for score differences across samples:")

  kw_results <- list()
  for (score_col in score_cols) {
    score_name <- gsub("^score_", "", score_col)
    kw_test <- tryCatch({
      kruskal.test(seurat_obj@meta.data[[score_col]] ~ seurat_obj@meta.data[[sample_col]])
    }, error = function(e) NULL)

    if (!is.null(kw_test)) {
      kw_results[[score_name]] <- list(statistic = kw_test$statistic, p.value = kw_test$p.value)
      sig_marker <- ifelse(kw_test$p.value < 0.001, "***",
                           ifelse(kw_test$p.value < 0.01, "**",
                                  ifelse(kw_test$p.value < 0.05, "*", "")))
      message(sprintf("  %s: H = %.2f, p = %.2e %s", score_name, kw_test$statistic, kw_test$p.value, sig_marker))
    }
  }

  # -------------------------------------------------------------------------
  # Sample Consistency (Jensen-Shannon Divergence)
  # -------------------------------------------------------------------------
  message("\n--- Sample Consistency Metrics ---")

  js_divergence <- function(p, q) {
    p <- p + 1e-10
    q <- q + 1e-10
    p <- p / sum(p)
    q <- q / sum(q)
    m <- (p + q) / 2
    0.5 * sum(p * log(p / m)) + 0.5 * sum(q * log(q / m))
  }

  js_matrix <- NULL
  if (n_samples >= 2) {
    js_matrix <- matrix(0, n_samples, n_samples)
    rownames(js_matrix) <- unique_samples
    colnames(js_matrix) <- unique_samples

    for (i in 1:(n_samples - 1)) {
      for (j in (i + 1):n_samples) {
        p <- sample_subtype_table[i, ]
        q <- sample_subtype_table[j, ]
        js_matrix[i, j] <- js_divergence(p, q)
        js_matrix[j, i] <- js_matrix[i, j]
      }
    }

    message("\nJensen-Shannon divergence between samples (lower = more similar):")
    print(round(js_matrix, 4))

    upper_tri <- js_matrix[upper.tri(js_matrix)]
    mean_js <- mean(upper_tri)
    message(sprintf("\nMean pairwise JS divergence: %.4f", mean_js))

    if (mean_js < 0.1) {
      message("  → Samples are highly consistent")
    } else if (mean_js < 0.3) {
      message("  → Samples show moderate consistency")
    } else {
      message("  → Samples show substantial heterogeneity")
    }
  }

  # Store results
  per_sample_results <- list(
    sample_column = sample_col,
    n_samples = n_samples,
    sample_names = unique_samples,
    contingency_table = sample_subtype_table,
    proportions = sample_proportions,
    subtype_summary = subtype_summary,
    chi_square_test = chi_test,
    sample_score_stats = sample_score_stats,
    kruskal_wallis_tests = kw_results,
    js_divergence_matrix = js_matrix
  )

  attr(seurat_obj, "per_sample_analysis") <- per_sample_results

} else {
  message("  No sample column found - skipping per-sample analysis")
}

# ==============================================================================
# COMPARE WITH ORIGINAL CLASSIFICATION (if exists)
# ==============================================================================

if ("epithelial_subtype" %in% colnames(seurat_obj@meta.data)) {
  message("\n========================================")
  message("Comparison with Original Classification")
  message("========================================\n")

  original <- seurat_obj$epithelial_subtype
  new_class <- seurat_obj$module_score_subtype
  confusion <- table(Original = original, ModuleScore = new_class)

  common_categories <- intersect(rownames(confusion), colnames(confusion))
  if (length(common_categories) > 0) {
    agreement_count <- sum(sapply(common_categories, function(cat) {
      if (cat %in% rownames(confusion) && cat %in% colnames(confusion)) confusion[cat, cat] else 0
    }))
    agreement <- agreement_count / sum(confusion) * 100
    message(sprintf("Overall agreement: %.1f%%", agreement))
  }
}

# ==============================================================================
# GENERATE DIAGNOSTIC PLOTS
# ==============================================================================

message("\n========================================")
message("Generating Diagnostic Plots")
message("========================================\n")

fig_dir <- get_path(config, config$paths$figures_supp_dir)
ensure_dir(fig_dir)
colors <- get_colors(config, "epithelial_subtypes")

# Plot 1: Module score distributions
message("Creating score distribution plot...")
score_data <- seurat_obj@meta.data[, c(score_cols, "module_score_subtype")]
score_long <- tidyr::pivot_longer(score_data, cols = -module_score_subtype,
                                  names_to = "Signature", values_to = "Score")
score_long$Signature <- gsub("^score_", "", score_long$Signature)

p1 <- ggplot(score_long, aes(x = Score, fill = Signature)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = min_threshold, linetype = "dashed", color = "red") +
  facet_wrap(~Signature, scales = "free_y", ncol = 2) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Module Score Distributions",
       subtitle = sprintf("Red line = classification threshold (%.2f)", min_threshold),
       x = "Module Score", y = "Density")

# Plot 2: Scores by assigned subtype
message("Creating scores by subtype plot...")
p2 <- ggplot(score_long, aes(x = module_score_subtype, y = Score, fill = module_score_subtype)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~Signature, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = colors, na.value = "gray50") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Module Scores by Assigned Subtype", x = "", y = "Score")

# Plot 3: Cell counts by subtype
message("Creating cell counts plot...")
count_data <- as.data.frame(table(seurat_obj$module_score_subtype))
colnames(count_data) <- c("Subtype", "Count")
count_data$Percentage <- count_data$Count / sum(count_data$Count) * 100

p4 <- ggplot(count_data, aes(x = reorder(Subtype, -Count), y = Count, fill = Subtype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = colors, na.value = "gray50") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Cell Counts by Subtype", x = "", y = "Number of Cells")

# Plot 4: Proliferation by subtype
message("Creating proliferation analysis plot...")
prolif_col <- ifelse("score_Transit_Amplifying" %in% colnames(seurat_obj@meta.data),
                     "score_Transit_Amplifying",
                     ifelse("score_Proliferation" %in% colnames(seurat_obj@meta.data),
                            "score_Proliferation", NULL))

if (!is.null(prolif_col)) {
  p5 <- ggplot(seurat_obj@meta.data, aes(x = module_score_subtype, y = .data[[prolif_col]],
                                         fill = module_score_subtype)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
    geom_hline(yintercept = prolif_threshold, linetype = "dashed", color = "red") +
    scale_fill_manual(values = colors, na.value = "gray50") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title = "Proliferation Score by Subtype",
         subtitle = sprintf("Red line = TA threshold (%.2f)", prolif_threshold),
         x = "", y = "Proliferation Score")
} else {
  p5 <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Proliferation score not available")
}

# Plot 5: Cell cycle by subtype
message("Creating cell cycle plot...")
if ("Phase" %in% colnames(seurat_obj@meta.data) && !all(seurat_obj$Phase == "Unknown")) {
  cycle_data <- seurat_obj@meta.data %>%
    count(module_score_subtype, Phase) %>%
    group_by(module_score_subtype) %>%
    mutate(Proportion = n / sum(n))

  cycle_colors <- get_colors(config, "cell_cycle")

  p6 <- ggplot(cycle_data, aes(x = module_score_subtype, y = Proportion, fill = Phase)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cycle_colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Cell Cycle Distribution by Subtype", x = "", y = "Proportion")
} else {
  p6 <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Cell cycle not scored")
}

# Combine main diagnostic plots
message("Combining diagnostic plots...")
diagnostic_plot <- (p1 | p2) / (p4 | p5) / (p6 | plot_spacer()) +
  plot_annotation(
    title = sprintf("Cell Type Annotation Diagnostics (%s)", toupper(dataset_type)),
    subtitle = sprintf("Module score classification with TA identification | n = %d cells", ncol(seurat_obj)),
    tag_levels = "A"
  )

diagnostic_path <- file.path(fig_dir, sprintf("01_cell_type_annotation_diagnostics_%s", dataset_type))
ggsave(paste0(diagnostic_path, ".pdf"), diagnostic_plot, width = 14, height = 16)
ggsave(paste0(diagnostic_path, ".png"), diagnostic_plot, width = 14, height = 16, dpi = 300)
message(sprintf("Saved: %s.pdf/.png", diagnostic_path))

# ==============================================================================
# PER-SAMPLE DIAGNOSTIC PLOTS
# ==============================================================================

if (!is.null(per_sample_results)) {
  message("\n--- Creating Per-Sample Diagnostic Plots ---")

  # Plot A: Stacked bar of subtype proportions by sample
  sample_composition_data <- seurat_obj@meta.data %>%
    group_by(!!sym(sample_col), module_score_subtype) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(!!sym(sample_col)) %>%
    mutate(total = sum(n), percentage = n / total * 100) %>%
    ungroup()

  sample_totals <- sample_composition_data %>%
    select(!!sym(sample_col), total) %>%
    distinct() %>%
    mutate(sample_label = sprintf("%s\n(n=%d)", !!sym(sample_col), total))

  sample_composition_data <- sample_composition_data %>%
    left_join(sample_totals, by = sample_col)

  p_sample_stack <- ggplot(sample_composition_data,
                           aes(x = sample_label, y = percentage, fill = module_score_subtype)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = colors, na.value = "gray50") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right") +
    labs(title = "Subtype Composition by Sample", x = "Sample", y = "Percentage (%)", fill = "Subtype")

  # Plot B: Side-by-side comparison
  p_sample_dodge <- ggplot(sample_composition_data,
                           aes(x = module_score_subtype, y = percentage, fill = !!sym(sample_col))) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") +
    labs(title = "Subtype Proportions Across Samples", x = "Subtype", y = "Percentage (%)", fill = "Sample")

  # Plot C: Score distributions by sample
  core_score_cols_plot <- grep("^score_(Basal_like|Intermediate|Specialized|Transit_Amplifying)",
                               colnames(seurat_obj@meta.data), value = TRUE)

  if (length(core_score_cols_plot) > 0) {
    score_sample_long <- seurat_obj@meta.data[, c(sample_col, core_score_cols_plot)] %>%
      tidyr::pivot_longer(cols = -all_of(sample_col), names_to = "Signature", values_to = "Score")
    score_sample_long$Signature <- gsub("^score_", "", score_sample_long$Signature)

    p_sample_scores <- ggplot(score_sample_long, aes(x = !!sym(sample_col), y = Score, fill = !!sym(sample_col))) +
      geom_violin(scale = "width", alpha = 0.7) +
      geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
      facet_wrap(~Signature, scales = "free_y", ncol = 2) +
      geom_hline(yintercept = min_threshold, linetype = "dashed", color = "red", alpha = 0.5) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      labs(title = "Module Score Distributions by Sample",
           subtitle = sprintf("Red dashed line = threshold (%.2f)", min_threshold),
           x = "Sample", y = "Module Score")
  } else {
    p_sample_scores <- ggplot() + theme_void()
  }

  # Plot D: Heatmap of proportions
  heatmap_data <- as.data.frame(as.table(as.matrix(per_sample_results$proportions)))
  colnames(heatmap_data) <- c("Sample", "Subtype", "Percentage")

  p_heatmap <- ggplot(heatmap_data, aes(x = Subtype, y = Sample, fill = Percentage)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f", Percentage)), size = 3) +
    scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue",
                         midpoint = max(heatmap_data$Percentage) / 2, name = "% of Sample") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
    labs(title = "Subtype Proportions Heatmap", x = "Subtype", y = "Sample")

  # Plot E: Cell cycle by sample
  if ("Phase" %in% colnames(seurat_obj@meta.data) && !all(seurat_obj$Phase == "Unknown")) {
    cycle_sample_data <- seurat_obj@meta.data %>%
      group_by(!!sym(sample_col), Phase) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(!!sym(sample_col)) %>%
      mutate(proportion = n / sum(n)) %>%
      ungroup()

    p_sample_cycle <- ggplot(cycle_sample_data, aes(x = !!sym(sample_col), y = proportion, fill = Phase)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = cycle_colors) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Cell Cycle by Sample", x = "Sample", y = "Proportion")
  } else {
    p_sample_cycle <- ggplot() + theme_void()
  }

  # Combine per-sample plots
  per_sample_plot <- (p_sample_stack | p_sample_dodge) / p_sample_scores / (p_heatmap | p_sample_cycle) +
    plot_annotation(
      title = sprintf("Per-Sample Analysis (%s)", toupper(dataset_type)),
      subtitle = sprintf("n = %d cells across %d samples", ncol(seurat_obj), per_sample_results$n_samples),
      tag_levels = "A"
    )

  per_sample_path <- file.path(fig_dir, sprintf("01_per_sample_diagnostics_%s", dataset_type))
  ggsave(paste0(per_sample_path, ".pdf"), per_sample_plot, width = 16, height = 14)
  ggsave(paste0(per_sample_path, ".png"), per_sample_plot, width = 16, height = 14, dpi = 300)
  message(sprintf("Saved: %s.pdf/.png", per_sample_path))
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

# Save Seurat object
output_path <- file.path(objects_dir, sprintf("01_seurat_annotated_%s.rds", dataset_type))
saveRDS(seurat_obj, output_path)
message(sprintf("Saved: %s", output_path))

# Overall classification summary
prolif_col_name <- ifelse("score_Transit_Amplifying" %in% colnames(seurat_obj@meta.data),
                          "score_Transit_Amplifying",
                          ifelse("score_Proliferation" %in% colnames(seurat_obj@meta.data),
                                 "score_Proliferation", NA))

summary_table <- seurat_obj@meta.data %>%
  group_by(module_score_subtype) %>%
  summarise(
    n_cells = n(),
    pct_total = n() / nrow(seurat_obj@meta.data) * 100,
    mean_Basal_like = if ("score_Basal_like" %in% names(.)) mean(score_Basal_like, na.rm = TRUE) else NA,
    mean_Intermediate = if ("score_Intermediate" %in% names(.)) mean(score_Intermediate, na.rm = TRUE) else NA,
    mean_Specialized = if ("score_Specialized" %in% names(.)) mean(score_Specialized, na.rm = TRUE) else NA,
    pct_cycling = if ("Phase" %in% names(.)) mean(Phase %in% c("S", "G2M")) * 100 else NA,
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

write.csv(summary_table, file.path(tables_dir, sprintf("01_classification_summary_%s.csv", dataset_type)), row.names = FALSE)
message(sprintf("Saved: 01_classification_summary_%s.csv", dataset_type))

# Per-sample outputs
if (!is.null(per_sample_results)) {
  # Long format breakdown
  sample_breakdown <- seurat_obj@meta.data %>%
    group_by(!!sym(sample_col), module_score_subtype) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(!!sym(sample_col)) %>%
    mutate(sample_total = sum(n_cells), pct_of_sample = n_cells / sample_total * 100) %>%
    ungroup()

  write.csv(sample_breakdown, file.path(tables_dir, sprintf("01_classification_by_sample_%s.csv", dataset_type)), row.names = FALSE)
  message(sprintf("Saved: 01_classification_by_sample_%s.csv", dataset_type))

  # Wide format proportions
  sample_proportions_wide <- sample_breakdown %>%
    select(!!sym(sample_col), module_score_subtype, pct_of_sample) %>%
    tidyr::pivot_wider(names_from = module_score_subtype, values_from = pct_of_sample, values_fill = 0)

  write.csv(sample_proportions_wide, file.path(tables_dir, sprintf("01_sample_proportions_wide_%s.csv", dataset_type)), row.names = FALSE)
  message(sprintf("Saved: 01_sample_proportions_wide_%s.csv", dataset_type))

  # Per-sample score statistics
  sample_score_stats_full <- seurat_obj@meta.data %>%
    group_by(!!sym(sample_col)) %>%
    summarise(
      n_cells = n(),
      across(all_of(score_cols), list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE)), .names = "{.col}_{.fn}"),
      pct_cycling = if ("Phase" %in% names(.)) mean(Phase %in% c("S", "G2M")) * 100 else NA,
      .groups = "drop"
    )

  write.csv(sample_score_stats_full, file.path(tables_dir, sprintf("01_sample_score_statistics_%s.csv", dataset_type)), row.names = FALSE)
  message(sprintf("Saved: 01_sample_score_statistics_%s.csv", dataset_type))

  # Subtype variability summary
  write.csv(per_sample_results$subtype_summary, file.path(tables_dir, sprintf("01_subtype_variability_%s.csv", dataset_type)), row.names = FALSE)
  message(sprintf("Saved: 01_subtype_variability_%s.csv", dataset_type))

  # Statistical test results
  stats_results <- data.frame(test = character(), statistic = numeric(), p_value = numeric(), interpretation = character(), stringsAsFactors = FALSE)

  if (!is.null(per_sample_results$chi_square_test)) {
    stats_results <- rbind(stats_results, data.frame(
      test = "Chi-square (sample x subtype)",
      statistic = per_sample_results$chi_square_test$statistic,
      p_value = per_sample_results$chi_square_test$p.value,
      interpretation = ifelse(per_sample_results$chi_square_test$p.value < 0.05, "Significant heterogeneity", "No significant heterogeneity"),
      stringsAsFactors = FALSE
    ))
  }

  for (sig_name in names(per_sample_results$kruskal_wallis_tests)) {
    kw <- per_sample_results$kruskal_wallis_tests[[sig_name]]
    stats_results <- rbind(stats_results, data.frame(
      test = sprintf("Kruskal-Wallis (%s)", sig_name),
      statistic = kw$statistic,
      p_value = kw$p.value,
      interpretation = ifelse(kw$p.value < 0.05, "Significant difference", "No significant difference"),
      stringsAsFactors = FALSE
    ))
  }

  if (!is.null(per_sample_results$js_divergence_matrix)) {
    mean_js <- mean(per_sample_results$js_divergence_matrix[upper.tri(per_sample_results$js_divergence_matrix)])
    stats_results <- rbind(stats_results, data.frame(
      test = "Mean Jensen-Shannon divergence",
      statistic = mean_js,
      p_value = NA,
      interpretation = ifelse(mean_js < 0.1, "Highly consistent", ifelse(mean_js < 0.3, "Moderately consistent", "Substantial heterogeneity")),
      stringsAsFactors = FALSE
    ))
  }

  write.csv(stats_results, file.path(tables_dir, sprintf("01_statistical_tests_%s.csv", dataset_type)), row.names = FALSE)
  message(sprintf("Saved: 01_statistical_tests_%s.csv", dataset_type))

  # JS divergence matrix
  if (!is.null(per_sample_results$js_divergence_matrix)) {
    write.csv(as.data.frame(per_sample_results$js_divergence_matrix),
              file.path(tables_dir, sprintf("01_sample_similarity_matrix_%s.csv", dataset_type)), row.names = TRUE)
    message(sprintf("Saved: 01_sample_similarity_matrix_%s.csv", dataset_type))
  }
}

# Dataset comparison (for merged data)
if ("dataset" %in% colnames(seurat_obj@meta.data) && dataset_type == "merged") {
  dataset_breakdown <- seurat_obj@meta.data %>%
    group_by(dataset, module_score_subtype) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(dataset) %>%
    mutate(dataset_total = sum(n_cells), pct_of_dataset = n_cells / dataset_total * 100) %>%
    ungroup()

  write.csv(dataset_breakdown, file.path(tables_dir, sprintf("01_classification_by_dataset_%s.csv", dataset_type)), row.names = FALSE)
  message(sprintf("Saved: 01_classification_by_dataset_%s.csv", dataset_type))

  # Wide format for easy comparison
  dataset_comparison_wide <- dataset_breakdown %>%
    select(dataset, module_score_subtype, pct_of_dataset) %>%
    tidyr::pivot_wider(names_from = dataset, values_from = pct_of_dataset, values_fill = 0)

  write.csv(dataset_comparison_wide, file.path(tables_dir, sprintf("01_dataset_comparison_%s.csv", dataset_type)), row.names = FALSE)
  message(sprintf("Saved: 01_dataset_comparison_%s.csv", dataset_type))
}

# Gene signatures used
signatures_df <- do.call(rbind, lapply(names(signatures), function(sig_name) {
  data.frame(signature = sig_name, gene = signatures[[sig_name]], stringsAsFactors = FALSE)
}))
signatures_df$in_dataset <- signatures_df$gene %in% rownames(seurat_obj)

write.csv(signatures_df, file.path(tables_dir, sprintf("01_gene_signatures_used_%s.csv", dataset_type)), row.names = FALSE)
message(sprintf("Saved: 01_gene_signatures_used_%s.csv", dataset_type))

# Signature coverage summary
sig_coverage <- signatures_df %>%
  group_by(signature) %>%
  summarise(total_genes = n(), genes_present = sum(in_dataset), coverage_pct = genes_present / total_genes * 100, .groups = "drop")

write.csv(sig_coverage, file.path(tables_dir, sprintf("01_signature_coverage_%s.csv", dataset_type)), row.names = FALSE)
message(sprintf("Saved: 01_signature_coverage_%s.csv", dataset_type))

# ==============================================================================
# COMPLETION
# ==============================================================================

message("\n")
message("================================================================")
message("  Cell Type Annotation Complete!")
message("================================================================")
message(sprintf("  Dataset type: %s", toupper(dataset_type)))
message(sprintf("  Finished: %s", Sys.time()))
message(sprintf("  Total cells: %d", ncol(seurat_obj)))
message(sprintf("  Subtypes identified: %d", length(unique(seurat_obj$module_score_subtype))))
if (!is.null(per_sample_results)) {
  message(sprintf("  Samples analyzed: %d", per_sample_results$n_samples))
}
if ("dataset" %in% colnames(seurat_obj@meta.data) && dataset_type == "merged") {
  message("\n  Dataset breakdown:")
  for (ds in unique(seurat_obj$dataset)) {
    n_ds <- sum(seurat_obj$dataset == ds)
    message(sprintf("    %s: %d cells", ds, n_ds))
  }
}
message("\nOutputs:")
message(sprintf("  Object: %s", output_path))
message(sprintf("  Tables: %s", tables_dir))
message(sprintf("  Figures: %s", fig_dir))
message("\nNext step: Run 02_trajectory_analysis.R")
message("================================================================\n")
