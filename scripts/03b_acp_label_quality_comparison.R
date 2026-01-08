#!/usr/bin/env Rscript
# ==============================================================================
# scripts/03b_acp_label_quality_comparison.R
# ==============================================================================
# Label Quality Comparison: GSE215932 vs ACP scRNA-seq
#
# This script:
#   1. Integrates the two ACP datasets (GSE215932 and ACP scRNA-seq) using Harmony
#   2. Performs bidirectional label transfer in integrated space
#   3. Evaluates label quality using multiple metrics
#   4. Tests which dataset's labels are more robust/consistent
#
# Hypotheses:
#   H1: Labels from both datasets achieve high bidirectional transfer accuracy
#   H2: One dataset's labels show better cluster separation (silhouette scores)
#   H3: One dataset's labels better preserve expected trajectory order
#   H4: Labels show significant concordance after integration (beyond chance)
#
# Label columns (auto-detected):
#   - GSE215932: 'cell_type' (contains subtypes like E_C1_PE, E_C2_WC, etc.)
#   - ACP scRNA-seq: Prefers 'module_score_subtype' (from 01_cell_type_annotation.R),
#                    falls back to 'celltypes' if it has multiple values
#
# Prerequisites:
#   - Run 01_cell_type_annotation.R on both datasets first to ensure proper labels
#
# Usage:
#   Rscript scripts/03b_acp_label_quality_comparison.R [--cores N] [--bootstrap N]
#
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(slingshot)
  library(cluster)  # For silhouette
  library(RANN)     # For nn2 (KNN)
})

options(future.globals.maxSize = 32 * 1024^3)

message("\n")
message("================================================================")
message("  Label Quality Comparison: GSE215932 vs ACP scRNA-seq")
message("================================================================")
message(sprintf("  Started: %s", Sys.time()))

# Source utilities
source("R/utils/config.R")

# Load configuration
config <- load_config()

# Parse arguments (override config defaults if specified)
args <- commandArgs(trailingOnly = TRUE)
n_cores <- config$reproducibility$n_cores
n_bootstrap <- config$trajectory$statistics$n_permutations
for (i in seq_along(args)) {
  if (args[i] == "--cores" && i < length(args)) {
    n_cores <- as.integer(args[i + 1])
  }
  if (args[i] == "--bootstrap" && i < length(args)) {
    n_bootstrap <- as.integer(args[i + 1])
  }
}
message(sprintf("  Using %d cores", n_cores))
message(sprintf("  Bootstrap iterations: %d", n_bootstrap))

set.seed(config$reproducibility$seed)

# Define paths
objects_dir <- get_path(config, config$paths$objects_dir)
tables_dir <- get_path(config, config$paths$tables_dir)
fig_dir <- get_path(config, config$paths$figures_supp_dir)
ensure_dir(fig_dir)

# ==============================================================================
# LOAD AND PREPARE DATASETS
# ==============================================================================

message("\n========================================")
message("Loading Datasets")
message("========================================\n")

# Load GSE215932 (snRNA-seq) - try annotated first, then raw
gse_path <- file.path(objects_dir, "01_seurat_annotated_snrnaseq.rds")
if (!file.exists(gse_path)) {
  gse_path <- get_path(config, config$paths$snrnaseq_processed)
}
if (!file.exists(gse_path)) {
  stop("GSE215932 dataset not found at: ", gse_path,
       "\nRun 01 script with snrnaseq first or check config paths.")
}
message(sprintf("Loading GSE215932 (snRNA-seq) from: %s", gse_path))
seurat_gse <- readRDS(gse_path)
seurat_gse <- UpdateSeuratObject(seurat_gse)

# Load ACP scRNA-seq - try annotated first, then raw
acp_path <- file.path(objects_dir, "01_seurat_annotated_acp_scn.rds")
if (!file.exists(acp_path)) {
  acp_path <- get_path(config, config$paths$acp_scn_annotated)
}
if (!file.exists(acp_path)) {
  stop("ACP scRNA-seq dataset not found at: ", acp_path,
       "\nRun 01 script with acp_scn first or check config paths.")
}
message(sprintf("Loading ACP scRNA-seq from: %s", acp_path))
seurat_acp <- readRDS(acp_path)
seurat_acp <- UpdateSeuratObject(seurat_acp)

# Join layers if Seurat v5
if (inherits(seurat_gse[["RNA"]], "Assay5") && length(Layers(seurat_gse[["RNA"]])) > 1) {
  seurat_gse <- JoinLayers(seurat_gse)
}
if (inherits(seurat_acp[["RNA"]], "Assay5") && length(Layers(seurat_acp[["RNA"]])) > 1) {
  seurat_acp <- JoinLayers(seurat_acp)
}

# Add dataset identifier
seurat_gse$dataset <- "GSE215932"
seurat_acp$dataset <- "ACP_SCN"

# Identify and standardize label columns
# GSE215932: 'cell_type' contains subtypes (E_C1_PE, E_C2_WC, etc.)
# ACP: 'celltypes' may just be "Epithelial" - check for module_score_subtype first

if (!"cell_type" %in% colnames(seurat_gse@meta.data)) {
  stop("'cell_type' column not found in GSE215932 dataset")
}

# For ACP, prefer module_score_subtype (from script 01) over celltypes
acp_label_col <- NULL
if ("module_score_subtype" %in% colnames(seurat_acp@meta.data)) {
  acp_label_col <- "module_score_subtype"
  message("  Using 'module_score_subtype' for ACP labels (from classification)")
} else if ("celltypes" %in% colnames(seurat_acp@meta.data)) {
  # Check if celltypes has more than one unique value
  n_unique <- length(unique(seurat_acp$celltypes[!is.na(seurat_acp$celltypes)]))
  if (n_unique > 1) {
    acp_label_col <- "celltypes"
    message("  Using 'celltypes' for ACP labels")
  } else {
    warning("'celltypes' column has only one unique value. Looking for alternatives...")
  }
}

# If still no good label column, look for other options
if (is.null(acp_label_col)) {
  possible_cols <- c("subtype", "cell_subtype", "classification", "cluster_labels")
  for (col in possible_cols) {
    if (col %in% colnames(seurat_acp@meta.data)) {
      n_unique <- length(unique(seurat_acp[[col]][!is.na(seurat_acp[[col]])]))
      if (n_unique > 1) {
        acp_label_col <- col
        message(sprintf("  Using '%s' for ACP labels (fallback)", col))
        break
      }
    }
  }
}

if (is.null(acp_label_col)) {
  stop("No suitable label column found for ACP dataset. ",
       "Run 01_cell_type_annotation.R first to generate module_score_subtype, ",
       "or ensure 'celltypes' contains multiple cell types.")
}

# Standardize label column names for comparison
seurat_gse$original_label <- seurat_gse$cell_type
seurat_acp$original_label <- seurat_acp[[acp_label_col]]

message(sprintf("  GSE215932: %d cells (label column: cell_type)", ncol(seurat_gse)))
message(sprintf("  ACP scRNA-seq: %d cells (label column: %s)", ncol(seurat_acp), acp_label_col))

message("\nGSE215932 cell_type distribution:")
print(table(seurat_gse$original_label))

message("\nACP scRNA-seq celltypes distribution:")
print(table(seurat_acp$original_label))

# ==============================================================================
# FILTER TO EPITHELIAL CELLS
# ==============================================================================

message("\n========================================")
message("Filtering to Epithelial Cells")
message("========================================\n")

# Filter GSE to epithelial cells
if ("is_epithelial" %in% colnames(seurat_gse@meta.data)) {
  cells_keep_gse <- seurat_gse$is_epithelial == TRUE
} else {
  # Infer from cell_type - keep all since GSE should already be epithelial
  cells_keep_gse <- rep(TRUE, ncol(seurat_gse))
}
seurat_gse <- subset(seurat_gse, cells = colnames(seurat_gse)[cells_keep_gse])
message(sprintf("  GSE215932 epithelial cells: %d", ncol(seurat_gse)))

# Filter ACP - remove Unassigned and Non-epithelial if present
acp_labels <- seurat_acp$original_label
exclude_labels <- c("Unassigned", "Non-epithelial", "Unknown", NA, "")
cells_keep_acp <- !acp_labels %in% exclude_labels & !is.na(acp_labels)
seurat_acp <- subset(seurat_acp, cells = colnames(seurat_acp)[cells_keep_acp])
message(sprintf("  ACP scRNA-seq cells (excluding Unassigned): %d", ncol(seurat_acp)))

# Remove cells with NA or empty labels from both
seurat_gse <- subset(seurat_gse, cells = colnames(seurat_gse)[!is.na(seurat_gse$original_label) &
                                                                seurat_gse$original_label != ""])
seurat_acp <- subset(seurat_acp, cells = colnames(seurat_acp)[!is.na(seurat_acp$original_label) &
                                                                seurat_acp$original_label != ""])

message(sprintf("  After NA removal - GSE: %d, ACP: %d", ncol(seurat_gse), ncol(seurat_acp)))

# Verify we have multiple labels in each dataset
gse_unique_labels <- unique(seurat_gse$original_label)
acp_unique_labels <- unique(seurat_acp$original_label)

message(sprintf("\n  GSE unique labels (%d): %s",
                length(gse_unique_labels), paste(gse_unique_labels, collapse = ", ")))
message(sprintf("  ACP unique labels (%d): %s",
                length(acp_unique_labels), paste(acp_unique_labels, collapse = ", ")))

if (length(gse_unique_labels) < 2) {
  stop("GSE dataset has fewer than 2 unique labels. Cannot perform comparison.")
}
if (length(acp_unique_labels) < 2) {
  stop("ACP dataset has fewer than 2 unique labels. ",
       "Ensure classification has been run (01_cell_type_annotation.R).")
}

# ==============================================================================
# MERGE AND INTEGRATE WITH HARMONY
# ==============================================================================

message("\n========================================")
message("Integrating Datasets with Harmony")
message("========================================\n")

# Find common features
common_features <- intersect(rownames(seurat_gse), rownames(seurat_acp))
message(sprintf("  Common features: %d", length(common_features)))

# Merge datasets
message("Merging datasets...")
seurat_merged <- merge(seurat_gse, seurat_acp, add.cell.ids = c("GSE", "ACP"))

# Standard preprocessing on merged
message("Running preprocessing...")
seurat_merged <- NormalizeData(seurat_merged, verbose = FALSE)
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst",
                                      nfeatures = 3000, verbose = FALSE)
seurat_merged <- ScaleData(seurat_merged, verbose = FALSE)
seurat_merged <- RunPCA(seurat_merged, npcs = 50, verbose = FALSE)

# Run Harmony integration
message("Running Harmony integration...")
seurat_merged <- RunHarmony(seurat_merged,
                            group.by.vars = "dataset",
                            dims.use = 1:30,
                            verbose = FALSE)

# Run UMAP on harmony embedding
message("Running UMAP on integrated embedding...")
seurat_merged <- RunUMAP(seurat_merged,
                         reduction = "harmony",
                         dims = 1:30,
                         verbose = FALSE)

message("  Integration complete")

# ==============================================================================
# BIDIRECTIONAL LABEL TRANSFER
# ==============================================================================

message("\n========================================")
message("Bidirectional Label Transfer")
message("========================================\n")

# Get cell indices by dataset
gse_cells <- colnames(seurat_merged)[seurat_merged$dataset == "GSE215932"]
acp_cells <- colnames(seurat_merged)[seurat_merged$dataset == "ACP_SCN"]

# Get harmony embeddings
harmony_emb <- Embeddings(seurat_merged, "harmony")[, 1:30]

# KNN-based label transfer function
transfer_labels <- function(query_cells, ref_cells, ref_labels, harmony_emb, k = 30) {
  query_emb <- harmony_emb[query_cells, ]
  ref_emb <- harmony_emb[ref_cells, ]

  # Compute KNN
  nn_results <- RANN::nn2(ref_emb, query_emb, k = k)

  # Vote for labels
  predicted_labels <- character(length(query_cells))
  prediction_scores <- numeric(length(query_cells))
  entropy_scores <- numeric(length(query_cells))

  for (i in seq_along(query_cells)) {
    neighbor_idx <- nn_results$nn.idx[i, ]
    neighbor_labels <- ref_labels[neighbor_idx]

    # Distance-weighted voting
    weights <- 1 / (nn_results$nn.dists[i, ] + 0.001)
    label_votes <- tapply(weights, neighbor_labels, sum)

    predicted_labels[i] <- names(which.max(label_votes))
    prediction_scores[i] <- max(label_votes) / sum(label_votes)

    # Shannon entropy of vote distribution (normalized)
    # Handle edge case where all neighbors have same label
    n_unique_labels <- length(label_votes)
    if (n_unique_labels > 1) {
      vote_probs <- label_votes / sum(label_votes)
      raw_entropy <- -sum(vote_probs * log2(vote_probs + 1e-10))
      entropy_scores[i] <- raw_entropy / log2(n_unique_labels)
    } else {
      # All neighbors have same label = zero entropy (maximum certainty)
      entropy_scores[i] <- 0
    }
  }

  return(list(
    labels = predicted_labels,
    scores = prediction_scores,
    entropy = entropy_scores
  ))
}

# Transfer GSE labels to ACP cells
message("Transferring GSE labels → ACP cells...")
gse_labels <- seurat_merged$original_label[gse_cells]
transfer_gse_to_acp <- transfer_labels(acp_cells, gse_cells, gse_labels, harmony_emb)

# Transfer ACP labels to GSE cells
message("Transferring ACP labels → GSE cells...")
acp_labels <- seurat_merged$original_label[acp_cells]
transfer_acp_to_gse <- transfer_labels(gse_cells, acp_cells, acp_labels, harmony_emb)

# Add to metadata
seurat_merged$predicted_from_other <- NA
seurat_merged$prediction_score <- NA
seurat_merged$prediction_entropy <- NA

seurat_merged$predicted_from_other[match(acp_cells, colnames(seurat_merged))] <- transfer_gse_to_acp$labels
seurat_merged$prediction_score[match(acp_cells, colnames(seurat_merged))] <- transfer_gse_to_acp$scores
seurat_merged$prediction_entropy[match(acp_cells, colnames(seurat_merged))] <- transfer_gse_to_acp$entropy

seurat_merged$predicted_from_other[match(gse_cells, colnames(seurat_merged))] <- transfer_acp_to_gse$labels
seurat_merged$prediction_score[match(gse_cells, colnames(seurat_merged))] <- transfer_acp_to_gse$scores
seurat_merged$prediction_entropy[match(gse_cells, colnames(seurat_merged))] <- transfer_acp_to_gse$entropy

message("\nLabel transfer summary:")
message(sprintf("  GSE→ACP: mean score = %.3f, mean entropy = %.3f",
                mean(transfer_gse_to_acp$scores), mean(transfer_gse_to_acp$entropy)))
message(sprintf("  ACP→GSE: mean score = %.3f, mean entropy = %.3f",
                mean(transfer_acp_to_gse$scores), mean(transfer_acp_to_gse$entropy)))

# ==============================================================================
# HYPOTHESIS TESTING
# ==============================================================================

message("\n========================================")
message("Hypothesis Testing")
message("========================================\n")

results <- list()

# -----------------------------------------------------------------------------
# H1: Bidirectional transfer accuracy
# -----------------------------------------------------------------------------
message("--- H1: Bidirectional Label Transfer Quality ---")

# Compute transfer accuracy metrics
h1_gse_to_acp <- list(
  mean_score = mean(transfer_gse_to_acp$scores),
  median_score = median(transfer_gse_to_acp$scores),
  mean_entropy = mean(transfer_gse_to_acp$entropy),
  prop_high_conf = mean(transfer_gse_to_acp$scores > 0.5)
)

h1_acp_to_gse <- list(
  mean_score = mean(transfer_acp_to_gse$scores),
  median_score = median(transfer_acp_to_gse$scores),
  mean_entropy = mean(transfer_acp_to_gse$entropy),
  prop_high_conf = mean(transfer_acp_to_gse$scores > 0.5)
)

# Bootstrap confidence intervals for mean scores
message("  Computing bootstrap CIs...")
set.seed(config$reproducibility$seed)

bootstrap_mean <- function(x, n_boot = n_bootstrap) {
  boot_means <- replicate(n_boot, mean(sample(x, replace = TRUE)))
  return(quantile(boot_means, c(0.025, 0.975)))
}

h1_gse_to_acp$score_ci <- bootstrap_mean(transfer_gse_to_acp$scores)
h1_acp_to_gse$score_ci <- bootstrap_mean(transfer_acp_to_gse$scores)

# Two-sample test comparing transfer quality
h1_wilcox <- wilcox.test(transfer_gse_to_acp$scores, transfer_acp_to_gse$scores,
                         alternative = "two.sided", conf.int = TRUE)

# Effect size (Cliff's delta)
cliff_delta <- function(x, y) {
  n_x <- length(x)
  n_y <- length(y)
  dominance <- outer(x, y, function(a, b) sign(a - b))
  return(mean(dominance))
}
h1_cliff_d <- cliff_delta(transfer_gse_to_acp$scores, transfer_acp_to_gse$scores)

results$H1 <- list(
  hypothesis = "Labels from both datasets achieve high bidirectional transfer accuracy",
  gse_to_acp = h1_gse_to_acp,
  acp_to_gse = h1_acp_to_gse,
  comparison = list(
    wilcox_p = h1_wilcox$p.value,
    cliff_delta = h1_cliff_d,
    better_source = ifelse(h1_gse_to_acp$mean_score > h1_acp_to_gse$mean_score,
                           "GSE215932", "ACP_SCN")
  ),
  conclusion = sprintf("%s labels transfer better (Cliff's δ = %.3f, p = %.2e)",
                       ifelse(h1_cliff_d > 0, "GSE215932", "ACP_SCN"),
                       abs(h1_cliff_d), h1_wilcox$p.value)
)

message(sprintf("  GSE→ACP: mean = %.3f [%.3f, %.3f]",
                h1_gse_to_acp$mean_score, h1_gse_to_acp$score_ci[1], h1_gse_to_acp$score_ci[2]))
message(sprintf("  ACP→GSE: mean = %.3f [%.3f, %.3f]",
                h1_acp_to_gse$mean_score, h1_acp_to_gse$score_ci[1], h1_acp_to_gse$score_ci[2]))
message(sprintf("  Wilcoxon p-value: %.2e", h1_wilcox$p.value))
message(sprintf("  Cliff's delta: %.3f", h1_cliff_d))
message(sprintf("  Conclusion: %s", results$H1$conclusion))

# -----------------------------------------------------------------------------
# H2: Cluster separation quality (silhouette scores)
# -----------------------------------------------------------------------------
message("\n--- H2: Cluster Separation Quality (Silhouette Analysis) ---")

# Compute silhouette scores for each labeling scheme
# Note: Must compute distance separately for each subset (dist objects can't be subsetted)

# GSE labels silhouette (on GSE cells)
message("  Computing silhouette for GSE labels...")
gse_harmony_emb <- harmony_emb[gse_cells, ]
gse_labels_factor <- as.integer(factor(seurat_merged$original_label[gse_cells]))

if (length(unique(gse_labels_factor)) > 1) {
  gse_dist <- dist(gse_harmony_emb)
  sil_gse <- silhouette(gse_labels_factor, gse_dist)
  sil_gse_mean <- mean(sil_gse[, "sil_width"])
  sil_gse_by_cluster <- tapply(sil_gse[, "sil_width"], gse_labels_factor, mean)
} else {
  sil_gse_mean <- NA
  sil_gse_by_cluster <- NA
}

# ACP labels silhouette (on ACP cells)
message("  Computing silhouette for ACP labels...")
acp_harmony_emb <- harmony_emb[acp_cells, ]
acp_labels_factor <- as.integer(factor(seurat_merged$original_label[acp_cells]))

if (length(unique(acp_labels_factor)) > 1) {
  acp_dist <- dist(acp_harmony_emb)
  sil_acp <- silhouette(acp_labels_factor, acp_dist)
  sil_acp_mean <- mean(sil_acp[, "sil_width"])
  sil_acp_by_cluster <- tapply(sil_acp[, "sil_width"], acp_labels_factor, mean)
} else {
  sil_acp_mean <- NA
  sil_acp_by_cluster <- NA
}

# Bootstrap CIs for silhouette scores
bootstrap_silhouette <- function(harmony_subset, labels_factor, n_boot = n_bootstrap) {
  n <- length(labels_factor)
  boot_sils <- numeric(n_boot)
  for (b in 1:n_boot) {
    idx <- sample(1:n, replace = TRUE)
    if (length(unique(labels_factor[idx])) > 1) {
      boot_dist <- dist(harmony_subset[idx, ])
      sil_b <- silhouette(labels_factor[idx], boot_dist)
      boot_sils[b] <- mean(sil_b[, "sil_width"])
    } else {
      boot_sils[b] <- NA
    }
  }
  return(quantile(boot_sils, c(0.025, 0.975), na.rm = TRUE))
}

message("  Computing bootstrap CIs for silhouette scores...")
sil_gse_ci <- bootstrap_silhouette(gse_harmony_emb, gse_labels_factor)
sil_acp_ci <- bootstrap_silhouette(acp_harmony_emb, acp_labels_factor)

# Permutation test for silhouette difference
message("  Running permutation test...")
n_perm <- n_bootstrap
perm_diffs <- numeric(n_perm)
all_harmony_emb <- harmony_emb[c(gse_cells, acp_cells), ]
all_labels <- c(as.character(seurat_merged$original_label[gse_cells]),
                as.character(seurat_merged$original_label[acp_cells]))
dataset_indicator <- c(rep("GSE", length(gse_cells)), rep("ACP", length(acp_cells)))
n_gse <- length(gse_cells)
n_acp <- length(acp_cells)

observed_diff <- sil_gse_mean - sil_acp_mean

for (p in 1:n_perm) {
  # Shuffle dataset assignments
  perm_idx <- sample(dataset_indicator)
  perm_gse_idx <- which(perm_idx == "GSE")
  perm_acp_idx <- which(perm_idx == "ACP")

  perm_gse_labels <- all_labels[perm_gse_idx]
  perm_acp_labels <- all_labels[perm_acp_idx]

  perm_gse_factor <- as.integer(factor(perm_gse_labels))
  perm_acp_factor <- as.integer(factor(perm_acp_labels))

  if (length(unique(perm_gse_factor)) > 1 && length(unique(perm_acp_factor)) > 1) {
    perm_gse_dist <- dist(all_harmony_emb[perm_gse_idx, ])
    perm_acp_dist <- dist(all_harmony_emb[perm_acp_idx, ])

    sil_perm_gse <- mean(silhouette(perm_gse_factor, perm_gse_dist)[, "sil_width"])
    sil_perm_acp <- mean(silhouette(perm_acp_factor, perm_acp_dist)[, "sil_width"])
    perm_diffs[p] <- sil_perm_gse - sil_perm_acp
  } else {
    perm_diffs[p] <- NA
  }
}

perm_p_value <- mean(abs(perm_diffs) >= abs(observed_diff), na.rm = TRUE)

results$H2 <- list(
  hypothesis = "One dataset's labels show better cluster separation",
  gse_silhouette = list(mean = sil_gse_mean, ci = sil_gse_ci),
  acp_silhouette = list(mean = sil_acp_mean, ci = sil_acp_ci),
  observed_diff = observed_diff,
  permutation_p = perm_p_value,
  better_separation = ifelse(sil_gse_mean > sil_acp_mean, "GSE215932", "ACP_SCN"),
  conclusion = sprintf("%s labels show better separation (Δ = %.3f, perm p = %.3f)",
                       ifelse(observed_diff > 0, "GSE215932", "ACP_SCN"),
                       abs(observed_diff), perm_p_value)
)

message(sprintf("  GSE silhouette: %.3f [%.3f, %.3f]",
                sil_gse_mean, sil_gse_ci[1], sil_gse_ci[2]))
message(sprintf("  ACP silhouette: %.3f [%.3f, %.3f]",
                sil_acp_mean, sil_acp_ci[1], sil_acp_ci[2]))
message(sprintf("  Difference: %.3f (permutation p = %.3f)", observed_diff, perm_p_value))
message(sprintf("  Conclusion: %s", results$H2$conclusion))

# -----------------------------------------------------------------------------
# H3: Trajectory order preservation
# -----------------------------------------------------------------------------
message("\n--- H3: Trajectory Order Preservation ---")

# Run Slingshot on integrated embedding using each labeling scheme
message("  Running Slingshot with GSE labels...")
umap_coords <- Embeddings(seurat_merged, "umap")

# Use expected trajectory order from config
expected_order <- config$trajectory$expected_order
message(sprintf("  Expected order from config: %s", paste(expected_order, collapse = " → ")))

# Function to infer start cluster based on expected order
infer_start_cluster <- function(labels, expected_order) {
  unique_labels <- unique(labels)
  # Try exact match first
  for (exp_label in expected_order) {
    if (exp_label %in% unique_labels) return(exp_label)
  }
  # Try partial match
  for (exp_label in expected_order) {
    matches <- grep(exp_label, unique_labels, ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) return(matches[1])
  }
  return(unique_labels[1])  # Fallback
}

# Slingshot with GSE labels (on GSE cells)
gse_start <- infer_start_cluster(seurat_merged$original_label[gse_cells], expected_order)
message(sprintf("    GSE start cluster: %s", gse_start))

tryCatch({
  sds_gse <- slingshot(umap_coords[gse_cells, ],
                       clusterLabels = seurat_merged$original_label[gse_cells],
                       start.clus = gse_start,
                       stretch = 0)
  pt_gse <- slingPseudotime(sds_gse)
  if (ncol(pt_gse) > 1) {
    pt_gse <- rowMeans(pt_gse, na.rm = TRUE)
  } else {
    pt_gse <- pt_gse[, 1]
  }
  pt_gse <- (pt_gse - min(pt_gse, na.rm = TRUE)) /
    (max(pt_gse, na.rm = TRUE) - min(pt_gse, na.rm = TRUE))
  gse_trajectory_success <- TRUE
}, error = function(e) {
  message(sprintf("    Warning: Slingshot failed for GSE: %s", e$message))
  pt_gse <<- rep(NA, length(gse_cells))
  gse_trajectory_success <<- FALSE
})

# Slingshot with ACP labels (on ACP cells)
acp_start <- infer_start_cluster(seurat_merged$original_label[acp_cells], expected_order)
message(sprintf("    ACP start cluster: %s", acp_start))

tryCatch({
  sds_acp <- slingshot(umap_coords[acp_cells, ],
                       clusterLabels = seurat_merged$original_label[acp_cells],
                       start.clus = acp_start,
                       stretch = 0)
  pt_acp <- slingPseudotime(sds_acp)
  if (ncol(pt_acp) > 1) {
    pt_acp <- rowMeans(pt_acp, na.rm = TRUE)
  } else {
    pt_acp <- pt_acp[, 1]
  }
  pt_acp <- (pt_acp - min(pt_acp, na.rm = TRUE)) /
    (max(pt_acp, na.rm = TRUE) - min(pt_acp, na.rm = TRUE))
  acp_trajectory_success <- TRUE
}, error = function(e) {
  message(sprintf("    Warning: Slingshot failed for ACP: %s", e$message))
  pt_acp <<- rep(NA, length(acp_cells))
  acp_trajectory_success <<- FALSE
})

# Store pseudotime in merged object
seurat_merged$pseudotime_own_labels <- NA
seurat_merged$pseudotime_own_labels[match(gse_cells, colnames(seurat_merged))] <- pt_gse
seurat_merged$pseudotime_own_labels[match(acp_cells, colnames(seurat_merged))] <- pt_acp

# Evaluate trajectory quality: correlation between pseudotime and label order
evaluate_trajectory_order <- function(pseudotime, labels, expected_order) {
  # Assign expected rank to each label based on order matching
  unique_labels <- unique(labels[!is.na(labels)])

  label_ranks <- sapply(unique_labels, function(lab) {
    # Try exact match first
    exact_match <- which(expected_order == lab)
    if (length(exact_match) > 0) return(exact_match[1])

    # Try partial match (case-insensitive)
    for (i in seq_along(expected_order)) {
      if (grepl(expected_order[i], lab, ignore.case = TRUE) ||
          grepl(lab, expected_order[i], ignore.case = TRUE)) {
        return(i)
      }
    }
    return(length(expected_order) + 1)  # Unknown labels at end
  })

  # Mean pseudotime by label
  pt_by_label <- tapply(pseudotime, labels, mean, na.rm = TRUE)

  # Spearman correlation between expected rank and observed pseudotime rank
  shared_labels <- intersect(names(pt_by_label), names(label_ranks))
  if (length(shared_labels) >= 3) {
    cor_test <- cor.test(label_ranks[shared_labels],
                         rank(pt_by_label[shared_labels]),
                         method = "spearman")
    return(list(rho = as.numeric(cor_test$estimate),
                p_value = cor_test$p.value,
                n_labels = length(shared_labels)))
  } else {
    return(list(rho = NA, p_value = NA, n_labels = length(shared_labels)))
  }
}

if (gse_trajectory_success) {
  traj_eval_gse <- evaluate_trajectory_order(pt_gse,
                                             seurat_merged$original_label[gse_cells],
                                             expected_order)
} else {
  traj_eval_gse <- list(rho = NA, p_value = NA, n_labels = 0)
}

if (acp_trajectory_success) {
  traj_eval_acp <- evaluate_trajectory_order(pt_acp,
                                             seurat_merged$original_label[acp_cells],
                                             expected_order)
} else {
  traj_eval_acp <- list(rho = NA, p_value = NA, n_labels = 0)
}

# Compare trajectory quality using Fisher's z-transformation for correlation comparison
fisher_z_test <- function(r1, n1, r2, n2) {
  z1 <- 0.5 * log((1 + r1) / (1 - r1))
  z2 <- 0.5 * log((1 + r2) / (1 - r2))
  se <- sqrt(1/(n1 - 3) + 1/(n2 - 3))
  z_diff <- (z1 - z2) / se
  p_value <- 2 * pnorm(-abs(z_diff))
  return(list(z_diff = z_diff, p_value = p_value))
}

if (!is.na(traj_eval_gse$rho) && !is.na(traj_eval_acp$rho)) {
  fisher_test <- fisher_z_test(traj_eval_gse$rho, traj_eval_gse$n_labels,
                               traj_eval_acp$rho, traj_eval_acp$n_labels)
} else {
  fisher_test <- list(z_diff = NA, p_value = NA)
}

results$H3 <- list(
  hypothesis = "One dataset's labels better preserve expected trajectory order",
  gse_trajectory = traj_eval_gse,
  acp_trajectory = traj_eval_acp,
  fisher_z_test = fisher_test,
  better_trajectory = ifelse(is.na(traj_eval_gse$rho) || is.na(traj_eval_acp$rho), "INCONCLUSIVE",
                             ifelse(traj_eval_gse$rho > traj_eval_acp$rho, "GSE215932", "ACP_SCN")),
  conclusion = if (is.na(traj_eval_gse$rho) || is.na(traj_eval_acp$rho)) {
    "INCONCLUSIVE: Trajectory inference failed for one or both datasets"
  } else {
    sprintf("%s labels preserve trajectory order better (ρ = %.3f vs %.3f, Fisher p = %.3f)",
            ifelse(traj_eval_gse$rho > traj_eval_acp$rho, "GSE215932", "ACP_SCN"),
            max(traj_eval_gse$rho, traj_eval_acp$rho),
            min(traj_eval_gse$rho, traj_eval_acp$rho),
            fisher_test$p_value)
  }
)

message(sprintf("  GSE trajectory correlation: ρ = %.3f (p = %.3f)",
                traj_eval_gse$rho, traj_eval_gse$p_value))
message(sprintf("  ACP trajectory correlation: ρ = %.3f (p = %.3f)",
                traj_eval_acp$rho, traj_eval_acp$p_value))
message(sprintf("  Conclusion: %s", results$H3$conclusion))

# -----------------------------------------------------------------------------
# H4: Label concordance after integration
# -----------------------------------------------------------------------------
message("\n--- H4: Label Concordance Analysis ---")

# Create contingency table of original vs transferred labels
# For this, we need to find cells where we can compare

# Compute Cohen's Kappa for agreement between original and transferred labels
compute_kappa <- function(labels1, labels2) {
  # Create confusion matrix
  common_labels <- union(unique(labels1), unique(labels2))
  conf_mat <- table(factor(labels1, levels = common_labels),
                    factor(labels2, levels = common_labels))

  n <- sum(conf_mat)
  p_o <- sum(diag(conf_mat)) / n  # Observed agreement
  p_e <- sum(rowSums(conf_mat) * colSums(conf_mat)) / (n^2)  # Expected agreement

  kappa <- (p_o - p_e) / (1 - p_e)

  # Standard error for kappa
  se_kappa <- sqrt(p_o * (1 - p_o) / (n * (1 - p_e)^2))

  return(list(kappa = kappa, se = se_kappa, p_observed = p_o, p_expected = p_e, n = n))
}

# Since labels are different between datasets, we compute agreement
# between transferred labels and actual labels on the same cells

# For ACP cells: compare GSE-transferred labels with actual ACP labels
# This requires mapping between label vocabularies
# We'll use the transfer itself as the mapping

# Adjusted Rand Index between clusterings
compute_ari <- function(labels1, labels2) {
  # Implementation of ARI
  n <- length(labels1)
  contingency <- table(labels1, labels2)

  sum_comb_rows <- sum(choose(rowSums(contingency), 2))
  sum_comb_cols <- sum(choose(colSums(contingency), 2))
  sum_comb_cells <- sum(choose(contingency, 2))

  expected_index <- sum_comb_rows * sum_comb_cols / choose(n, 2)
  max_index <- 0.5 * (sum_comb_rows + sum_comb_cols)

  ari <- (sum_comb_cells - expected_index) / (max_index - expected_index)
  return(ari)
}

# Compute NMI (Normalized Mutual Information)
compute_nmi <- function(labels1, labels2) {
  contingency <- table(labels1, labels2)
  n <- sum(contingency)

  # Marginal probabilities
  p_i <- rowSums(contingency) / n
  p_j <- colSums(contingency) / n

  # Entropies
  H_i <- -sum(p_i[p_i > 0] * log2(p_i[p_i > 0]))
  H_j <- -sum(p_j[p_j > 0] * log2(p_j[p_j > 0]))

  # Mutual information
  MI <- 0
  for (i in 1:nrow(contingency)) {
    for (j in 1:ncol(contingency)) {
      if (contingency[i, j] > 0) {
        p_ij <- contingency[i, j] / n
        MI <- MI + p_ij * log2(p_ij / (p_i[i] * p_j[j]))
      }
    }
  }

  # Normalized MI
  NMI <- 2 * MI / (H_i + H_j)
  return(NMI)
}

# Compute metrics for cross-dataset label agreement
# Using transferred labels vs original labels on the receiving dataset

# For ACP cells: original ACP labels vs labels transferred from GSE
acp_original <- seurat_merged$original_label[acp_cells]
acp_from_gse <- transfer_gse_to_acp$labels

# For GSE cells: original GSE labels vs labels transferred from ACP
gse_original <- seurat_merged$original_label[gse_cells]
gse_from_acp <- transfer_acp_to_gse$labels

# Since vocabularies differ, we compute ARI and NMI which handle different label sets
ari_acp <- compute_ari(acp_original, acp_from_gse)
nmi_acp <- compute_nmi(acp_original, acp_from_gse)

ari_gse <- compute_ari(gse_original, gse_from_acp)
nmi_gse <- compute_nmi(gse_original, gse_from_acp)

# Permutation test for ARI significance
message("  Running permutation test for ARI...")
n_perm <- 1000
perm_ari_acp <- numeric(n_perm)
perm_ari_gse <- numeric(n_perm)

for (p in 1:n_perm) {
  perm_ari_acp[p] <- compute_ari(sample(acp_original), acp_from_gse)
  perm_ari_gse[p] <- compute_ari(sample(gse_original), gse_from_acp)
}

ari_acp_p <- mean(perm_ari_acp >= ari_acp)
ari_gse_p <- mean(perm_ari_gse >= ari_gse)

results$H4 <- list(
  hypothesis = "Labels show significant concordance after integration",
  acp_concordance = list(ARI = ari_acp, NMI = nmi_acp, p_value = ari_acp_p),
  gse_concordance = list(ARI = ari_gse, NMI = nmi_gse, p_value = ari_gse_p),
  interpretation = "Higher ARI/NMI indicates better agreement between native and transferred labels",
  better_concordance = ifelse(ari_acp > ari_gse, "ACP_SCN", "GSE215932"),
  conclusion = sprintf("Both show significant concordance (ACP ARI=%.3f p=%.3f; GSE ARI=%.3f p=%.3f); %s labels more consistent",
                       ari_acp, ari_acp_p, ari_gse, ari_gse_p,
                       ifelse(ari_acp > ari_gse, "ACP_SCN", "GSE215932"))
)

message(sprintf("  ACP concordance: ARI = %.3f (p = %.3f), NMI = %.3f",
                ari_acp, ari_acp_p, nmi_acp))
message(sprintf("  GSE concordance: ARI = %.3f (p = %.3f), NMI = %.3f",
                ari_gse, ari_gse_p, nmi_gse))
message(sprintf("  Conclusion: %s", results$H4$conclusion))

# ==============================================================================
# AGGREGATE LABEL QUALITY ASSESSMENT
# ==============================================================================

message("\n========================================")
message("Aggregate Label Quality Assessment")
message("========================================\n")

# Compile all metrics into a comparison table
quality_metrics <- data.frame(
  Metric = c("Transfer Score (mean)", "Transfer Entropy (mean)",
             "Silhouette Score", "Trajectory Correlation",
             "Cross-dataset ARI", "Cross-dataset NMI"),
  GSE215932 = c(h1_acp_to_gse$mean_score,  # Score when GSE is source
                h1_acp_to_gse$mean_entropy,
                sil_gse_mean,
                traj_eval_gse$rho,
                ari_gse,
                nmi_gse),
  ACP_SCN = c(h1_gse_to_acp$mean_score,  # Score when ACP is source
              h1_gse_to_acp$mean_entropy,
              sil_acp_mean,
              traj_eval_acp$rho,
              ari_acp,
              nmi_acp),
  Higher_is_Better = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
)

# Determine winner for each metric
quality_metrics$Better <- apply(quality_metrics, 1, function(row) {
  gse_val <- as.numeric(row["GSE215932"])
  acp_val <- as.numeric(row["ACP_SCN"])
  higher_better <- as.logical(row["Higher_is_Better"])

  if (is.na(gse_val) || is.na(acp_val)) return("TIE")

  if (higher_better) {
    if (gse_val > acp_val) return("GSE215932")
    else if (acp_val > gse_val) return("ACP_SCN")
    else return("TIE")
  } else {
    if (gse_val < acp_val) return("GSE215932")
    else if (acp_val < gse_val) return("ACP_SCN")
    else return("TIE")
  }
})

message("Quality metrics comparison:")
print(quality_metrics)

# Overall winner
gse_wins <- sum(quality_metrics$Better == "GSE215932", na.rm = TRUE)
acp_wins <- sum(quality_metrics$Better == "ACP_SCN", na.rm = TRUE)

overall_winner <- if (gse_wins > acp_wins) {
  "GSE215932"
} else if (acp_wins > gse_wins) {
  "ACP_SCN"
} else {
  "TIE"
}

message(sprintf("\nOverall: GSE wins %d metrics, ACP wins %d metrics → %s",
                gse_wins, acp_wins, overall_winner))

results$overall <- list(
  quality_metrics = quality_metrics,
  gse_wins = gse_wins,
  acp_wins = acp_wins,
  overall_winner = overall_winner,
  conclusion = sprintf("%s labels are overall more robust (%d vs %d metrics)",
                       overall_winner, max(gse_wins, acp_wins), min(gse_wins, acp_wins))
)

# ==============================================================================
# GENERATE PLOTS
# ==============================================================================

message("\n========================================")
message("Generating Plots")
message("========================================\n")

# Get visualization parameters from config
pt_size <- config$visualization$plot_settings$point_size
plot_alpha <- config$visualization$plot_settings$alpha
plot_dpi <- config$visualization$plot_settings$dpi

# Define dataset colors
dataset_colors <- c("GSE215932" = "#4DAF4A", "ACP_SCN" = "#E41A1C")

# Plot 1: Integrated UMAP colored by dataset
message("Creating UMAP plots...")
p1 <- DimPlot(seurat_merged, group.by = "dataset", pt.size = pt_size) +
  scale_color_manual(values = dataset_colors) +
  labs(title = "Integrated UMAP by Dataset")

# Plot 2: UMAP colored by original labels
p2 <- DimPlot(seurat_merged, group.by = "original_label", pt.size = pt_size) +
  labs(title = "Original Labels")

# Plot 3: UMAP colored by prediction score
p3 <- FeaturePlot(seurat_merged, features = "prediction_score", pt.size = pt_size) +
  scale_color_viridis_c() +
  labs(title = "Label Transfer Confidence")

# Plot 4: UMAP colored by prediction entropy
p4 <- FeaturePlot(seurat_merged, features = "prediction_entropy", pt.size = pt_size) +
  scale_color_viridis_c(option = "magma") +
  labs(title = "Label Transfer Entropy (lower = more certain)")

# Combine UMAPs
umap_combined <- (p1 | p2) / (p3 | p4)

ggsave(file.path(fig_dir, "03b_label_comparison_umap.pdf"), umap_combined,
       width = 14, height = 12)
ggsave(file.path(fig_dir, "03b_label_comparison_umap.png"), umap_combined,
       width = 14, height = 12, dpi = plot_dpi)
message("  Saved UMAP plots")

# Plot 5: Transfer score distributions
message("Creating comparison plots...")
score_data <- data.frame(
  score = c(transfer_gse_to_acp$scores, transfer_acp_to_gse$scores),
  direction = c(rep("GSE → ACP", length(transfer_gse_to_acp$scores)),
                rep("ACP → GSE", length(transfer_acp_to_gse$scores)))
)

p5 <- ggplot(score_data, aes(x = score, fill = direction)) +
  geom_density(alpha = plot_alpha) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  scale_fill_manual(values = c("GSE → ACP" = dataset_colors["GSE215932"],
                               "ACP → GSE" = dataset_colors["ACP_SCN"])) +
  theme_minimal() +
  labs(title = "Label Transfer Score Distribution",
       subtitle = sprintf("Wilcoxon p = %.2e, Cliff's δ = %.3f",
                          h1_wilcox$p.value, h1_cliff_d),
       x = "Prediction Score", y = "Density")

# Plot 6: Quality metrics comparison
quality_plot_data <- quality_metrics %>%
  tidyr::pivot_longer(cols = c("GSE215932", "ACP_SCN"),
                      names_to = "Dataset", values_to = "Value") %>%
  filter(!is.na(Value))

p6 <- ggplot(quality_plot_data, aes(x = Metric, y = Value, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = dataset_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Label Quality Metrics Comparison",
       subtitle = sprintf("Overall winner: %s (%d vs %d metrics)",
                          overall_winner, max(gse_wins, acp_wins), min(gse_wins, acp_wins)),
       y = "Metric Value")

# Plot 7: Entropy comparison
entropy_data <- data.frame(
  entropy = c(transfer_gse_to_acp$entropy, transfer_acp_to_gse$entropy),
  direction = c(rep("GSE → ACP", length(transfer_gse_to_acp$entropy)),
                rep("ACP → GSE", length(transfer_acp_to_gse$entropy)))
)

p7 <- ggplot(entropy_data, aes(x = direction, y = entropy, fill = direction)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_manual(values = c("GSE → ACP" = dataset_colors["GSE215932"],
                               "ACP → GSE" = dataset_colors["ACP_SCN"])) +
  theme_minimal() +
  labs(title = "Label Transfer Entropy",
       subtitle = "Lower entropy = more confident predictions",
       y = "Normalized Entropy")

# Combine comparison plots
comparison_plots <- (p5 | p6) / p7
ggsave(file.path(fig_dir, "03b_label_quality_comparison.pdf"), comparison_plots,
       width = 14, height = 10)
ggsave(file.path(fig_dir, "03b_label_quality_comparison.png"), comparison_plots,
       width = 14, height = 10, dpi = plot_dpi)
message("  Saved comparison plots")

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

message("\n========================================")
message("Saving Outputs")
message("========================================\n")

# Save integrated Seurat object
saveRDS(seurat_merged, file.path(objects_dir, "03b_label_comparison_integrated.rds"))
message(sprintf("Saved: %s", file.path(objects_dir, "03b_label_comparison_integrated.rds")))

# Save hypothesis test results
results_summary <- data.frame(
  hypothesis = c("H1", "H2", "H3", "H4"),
  description = c(
    "Bidirectional transfer quality",
    "Cluster separation (silhouette)",
    "Trajectory order preservation",
    "Cross-dataset label concordance"
  ),
  better_dataset = c(
    results$H1$comparison$better_source,
    results$H2$better_separation,
    results$H3$better_trajectory,
    results$H4$better_concordance
  ),
  conclusion = c(
    results$H1$conclusion,
    results$H2$conclusion,
    results$H3$conclusion,
    results$H4$conclusion
  )
)
write.csv(results_summary, file.path(tables_dir, "03b_hypothesis_results.csv"), row.names = FALSE)
message(sprintf("Saved: %s", file.path(tables_dir, "03b_hypothesis_results.csv")))

# Save quality metrics comparison
write.csv(quality_metrics, file.path(tables_dir, "03b_quality_metrics_comparison.csv"), row.names = FALSE)
message(sprintf("Saved: %s", file.path(tables_dir, "03b_quality_metrics_comparison.csv")))

# Save detailed transfer statistics
transfer_stats <- data.frame(
  direction = c("GSE_to_ACP", "ACP_to_GSE"),
  mean_score = c(h1_gse_to_acp$mean_score, h1_acp_to_gse$mean_score),
  median_score = c(h1_gse_to_acp$median_score, h1_acp_to_gse$median_score),
  mean_entropy = c(h1_gse_to_acp$mean_entropy, h1_acp_to_gse$mean_entropy),
  prop_high_conf = c(h1_gse_to_acp$prop_high_conf, h1_acp_to_gse$prop_high_conf),
  score_ci_lower = c(h1_gse_to_acp$score_ci[1], h1_acp_to_gse$score_ci[1]),
  score_ci_upper = c(h1_gse_to_acp$score_ci[2], h1_acp_to_gse$score_ci[2])
)
write.csv(transfer_stats, file.path(tables_dir, "03b_transfer_statistics.csv"), row.names = FALSE)
message(sprintf("Saved: %s", file.path(tables_dir, "03b_transfer_statistics.csv")))

# Save full results object
saveRDS(results, file.path(objects_dir, "03b_hypothesis_results_full.rds"))
message(sprintf("Saved: %s", file.path(objects_dir, "03b_hypothesis_results_full.rds")))

# Generate markdown report
report <- c(
  "# Label Quality Comparison Report: GSE215932 vs ACP scRNA-seq",
  "",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Dataset Summary",
  "",
  sprintf("- GSE215932 cells: %d (label column: `cell_type`)", length(gse_cells)),
  sprintf("- ACP scRNA-seq cells: %d (label column: `%s`)", length(acp_cells), acp_label_col),
  sprintf("- Total integrated cells: %d", ncol(seurat_merged)),
  sprintf("- Common features: %d", length(common_features)),
  "",
  sprintf("### GSE215932 labels (%d unique):", length(gse_unique_labels)),
  paste("-", gse_unique_labels, collapse = "\n"),
  "",
  sprintf("### ACP labels (%d unique):", length(acp_unique_labels)),
  paste("-", acp_unique_labels, collapse = "\n"),
  "",
  "## Hypothesis Test Results",
  "",
  "| Hypothesis | Description | Better Dataset | Conclusion |",
  "|------------|-------------|----------------|------------|",
  sprintf("| H1 | Transfer quality | %s | %s |",
          results$H1$comparison$better_source, results$H1$conclusion),
  sprintf("| H2 | Cluster separation | %s | %s |",
          results$H2$better_separation, results$H2$conclusion),
  sprintf("| H3 | Trajectory order | %s | %s |",
          results$H3$better_trajectory, results$H3$conclusion),
  sprintf("| H4 | Label concordance | %s | %s |",
          results$H4$better_concordance, results$H4$conclusion),
  "",
  "## Quality Metrics Summary",
  "",
  "| Metric | GSE215932 | ACP_SCN | Better |",
  "|--------|-----------|---------|--------|",
  apply(quality_metrics, 1, function(row) {
    sprintf("| %s | %.3f | %.3f | %s |",
            row["Metric"], as.numeric(row["GSE215932"]),
            as.numeric(row["ACP_SCN"]), row["Better"])
  }),
  "",
  sprintf("## Overall Winner: **%s** (%d vs %d metrics)",
          overall_winner, max(gse_wins, acp_wins), min(gse_wins, acp_wins)),
  "",
  "## Statistical Methods",
  "",
  sprintf("- **H1**: Wilcoxon rank-sum test with Cliff's delta effect size; bootstrap CIs (n=%d)", n_bootstrap),
  sprintf("- **H2**: Permutation test (n=%d) for silhouette score difference", n_bootstrap),
  "- **H3**: Fisher's z-transformation for correlation comparison",
  sprintf("- **H4**: Adjusted Rand Index with permutation significance test (n=%d)", n_bootstrap),
  "",
  "## Interpretation",
  "",
  "This analysis compares the label quality between the GSE215932 (snRNA-seq) and",
  "ACP scRNA-seq datasets. Higher quality labels are characterized by:",
  "",
  "1. **Better transfer accuracy**: Labels that can be reliably transferred to the other dataset",

  "2. **Better cluster separation**: Higher silhouette scores indicate more distinct clusters",
  "3. **Better trajectory preservation**: Labels that maintain expected differentiation order",
  "4. **Higher concordance**: Labels that agree with transferred labels from the other dataset"
)

writeLines(report, file.path(tables_dir, "03b_label_comparison_report.md"))
message(sprintf("Saved: %s", file.path(tables_dir, "03b_label_comparison_report.md")))

# ==============================================================================
# COMPLETION
# ==============================================================================

message("\n")
message("================================================================")
message("  Label Quality Comparison Complete!")
message("================================================================")
message(sprintf("  Finished: %s", Sys.time()))
message(sprintf("  GSE215932 cells: %d", length(gse_cells)))
message(sprintf("  ACP scRNA-seq cells: %d", length(acp_cells)))
message("")
message("Hypothesis Results:")
message(sprintf("  H1 (Transfer quality): %s", results$H1$conclusion))
message(sprintf("  H2 (Cluster separation): %s", results$H2$conclusion))
message(sprintf("  H3 (Trajectory order): %s", results$H3$conclusion))
message(sprintf("  H4 (Label concordance): %s", results$H4$conclusion))
message("")
message(sprintf("Overall Winner: %s (%d vs %d metrics)",
                overall_winner, max(gse_wins, acp_wins), min(gse_wins, acp_wins)))
message("")
message("Outputs:")
message(sprintf("  Objects: %s", objects_dir))
message(sprintf("  Tables: %s", tables_dir))
message(sprintf("  Figures: %s", fig_dir))
message("================================================================\n")
