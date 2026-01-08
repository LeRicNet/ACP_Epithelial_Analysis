#!/usr/bin/env Rscript
# ==============================================================================
# GSE cell_type Label Validation - Rigorous Analysis (Revised)
# ==============================================================================
# Tests whether GSE215932 cell_type labels correspond to transcriptionally
# distinct populations using multiple complementary approaches.
#
# Analyses performed:
#   1. KNN cross-validation with balanced accuracy
#   2. Permutation testing with full CV (proper null distribution)
#   3. Sensitivity analysis across parameters
#   4. Confusion matrix analysis (normalized by actual class)
#   5. Silhouette analysis with stratified sampling and bootstrap CIs
#   6. Marker gene validation with multiple thresholds and power consideration
#   7. Alternative embedding comparison
#   8. Cross-dataset transfer to ACP_SCN (via Harmony integration)
#
# Usage:
#   Rscript label_validation_rigorous.R [--config path/to/config.yaml]
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(caret)
  library(class)
  library(RANN)
  library(cluster)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
})

# ==============================================================================
# SETUP
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
for (i in seq_along(args)) {
  if (args[i] == "--config" && i < length(args)) config_path <- args[i + 1]
}

find_project_root <- function() {
  candidates <- c(getwd(), normalizePath(".."),
                  if (requireNamespace("here", quietly = TRUE)) here::here() else NULL)
  for (path in candidates) {
    if (!is.null(path) && file.exists(file.path(path, "config/config.yaml"))) return(path)
  }
  current <- getwd()
  for (i in 1:5) {
    if (file.exists(file.path(current, "config/config.yaml"))) return(current)
    current <- dirname(current)
  }
  return(NULL)
}

project_root <- find_project_root()
if (is.null(project_root)) stop("Could not find project root.")
setwd(project_root)

source("R/utils/config.R")
config <- load_config(config_path)
set.seed(config$reproducibility$seed)

message("\n")
message("================================================================")
message("  GSE cell_type Label Validation - Rigorous Analysis (Revised)")
message("================================================================")
message(sprintf("  Started: %s", Sys.time()))
message(sprintf("  Seed: %d", config$reproducibility$seed))

# ==============================================================================
# ANALYSIS FUNCTIONS
# ==============================================================================

#' KNN Cross-Validation with Balanced Accuracy
#'
#' @param embeddings Numeric matrix of cell embeddings
#' @param labels Character/factor vector of cell labels
#' @param k Number of neighbors
#' @param n_folds Number of CV folds
#' @param seed Random seed
#' @return List with accuracy metrics and confusion matrix
knn_cv <- function(embeddings, labels, k = 15, n_folds = 5, seed = 42) {

  valid <- !is.na(labels) & labels != "" & labels != "Unassigned"
  emb <- embeddings[valid, , drop = FALSE]
  labels <- factor(labels[valid])

  if (length(unique(labels)) < 2) {
    return(list(overall_accuracy = NA, balanced_accuracy = NA, error = "< 2 classes"))
  }

  set.seed(seed)
  folds <- createFolds(labels, k = n_folds, returnTrain = FALSE)

  all_predictions <- character(length(labels))
  all_actuals <- character(length(labels))

  fold_results <- lapply(seq_along(folds), function(i) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(labels), test_idx)

    pred <- knn(emb[train_idx, ], emb[test_idx, ], labels[train_idx], k = k)

    all_predictions[test_idx] <<- as.character(pred)
    all_actuals[test_idx] <<- as.character(labels[test_idx])

    list(
      overall_acc = mean(pred == labels[test_idx]),
      per_class_acc = sapply(levels(labels), function(cls) {
        idx <- labels[test_idx] == cls
        if (sum(idx) == 0) return(NA)
        mean(pred[idx] == labels[test_idx][idx])
      })
    )
  })

  # Aggregate confusion matrix
  conf_mat <- table(Predicted = factor(all_predictions, levels = levels(labels)),
                    Actual = factor(all_actuals, levels = levels(labels)))

  # Calculate metrics
  overall_accs <- sapply(fold_results, `[[`, "overall_acc")
  per_class_mat <- do.call(rbind, lapply(fold_results, `[[`, "per_class_acc"))
  per_class_mean <- colMeans(per_class_mat, na.rm = TRUE)

  # Balanced accuracy (mean of per-class accuracies)
  balanced_acc <- mean(per_class_mean, na.rm = TRUE)

  list(
    overall_accuracy = mean(overall_accs),
    overall_sd = sd(overall_accs),
    balanced_accuracy = balanced_acc,
    balanced_sd = sd(rowMeans(per_class_mat, na.rm = TRUE)),
    per_class_accuracy = per_class_mean,
    per_class_sd = apply(per_class_mat, 2, sd, na.rm = TRUE),
    confusion_matrix = conf_mat,
    n_cells = sum(valid),
    n_classes = length(levels(labels)),
    class_counts = table(labels)
  )
}

#' Permutation Test with Full CV (Proper Null Distribution)
#'
#' @param embeddings Numeric matrix
#' @param labels Character/factor vector
#' @param k KNN parameter
#' @param n_folds Number of CV folds
#' @param n_perm Number of permutations
#' @param seed Random seed
#' @return List with null distribution statistics
permutation_test <- function(embeddings, labels, k = 15, n_folds = 5, n_perm = 50, seed = 42) {

  valid <- !is.na(labels) & labels != "" & labels != "Unassigned"
  emb <- embeddings[valid, , drop = FALSE]
  labels <- factor(labels[valid])

  set.seed(seed)

  message(sprintf("    Running %d permutations with full %d-fold CV each...", n_perm, n_folds))

  null_overall <- numeric(n_perm)
  null_balanced <- numeric(n_perm)

  for (p in 1:n_perm) {
    if (p %% 10 == 0) message(sprintf("      Permutation %d/%d", p, n_perm))

    perm_labels <- sample(labels)

    # Full CV on permuted labels
    folds <- createFolds(perm_labels, k = n_folds, returnTrain = FALSE)

    fold_accs <- sapply(folds, function(test_idx) {
      train_idx <- setdiff(seq_along(perm_labels), test_idx)
      pred <- knn(emb[train_idx, ], emb[test_idx, ], perm_labels[train_idx], k = k)

      overall <- mean(pred == perm_labels[test_idx])
      per_class <- sapply(levels(labels), function(cls) {
        idx <- perm_labels[test_idx] == cls
        if (sum(idx) == 0) return(NA)
        mean(pred[idx] == perm_labels[test_idx][idx])
      })
      c(overall = overall, balanced = mean(per_class, na.rm = TRUE))
    })

    null_overall[p] <- mean(fold_accs["overall", ])
    null_balanced[p] <- mean(fold_accs["balanced", ])
  }

  list(
    null_overall_mean = mean(null_overall),
    null_overall_sd = sd(null_overall),
    null_overall_95_ci = quantile(null_overall, c(0.025, 0.975)),
    null_balanced_mean = mean(null_balanced),
    null_balanced_sd = sd(null_balanced),
    null_balanced_95_ci = quantile(null_balanced, c(0.025, 0.975)),
    null_overall_distribution = null_overall,
    null_balanced_distribution = null_balanced
  )
}

#' Sensitivity Analysis Across Parameters
#'
#' @param embeddings Numeric matrix
#' @param labels Character/factor vector
#' @param k_values Vector of k values to test
#' @param dim_values Vector of dimension counts to test
#' @param seed Random seed
#' @return Data frame with results for each parameter combination
sensitivity_analysis <- function(embeddings, labels,
                                 k_values = c(5, 15, 30, 50),
                                 dim_values = c(10, 20, 30, 50),
                                 seed = 42) {

  max_dims <- ncol(embeddings)
  dim_values <- dim_values[dim_values <= max_dims]

  results <- expand.grid(k = k_values, dims = dim_values)
  results$overall_accuracy <- NA
  results$balanced_accuracy <- NA

  for (i in 1:nrow(results)) {
    k <- results$k[i]
    dims <- results$dims[i]

    emb_subset <- embeddings[, 1:dims, drop = FALSE]
    cv_result <- knn_cv(emb_subset, labels, k = k, n_folds = 5, seed = seed)

    results$overall_accuracy[i] <- cv_result$overall_accuracy
    results$balanced_accuracy[i] <- cv_result$balanced_accuracy
  }

  results
}

#' Marker Gene Validation with Multiple Thresholds
#'
#' @param seurat_obj Seurat object
#' @param label_col Column name with labels
#' @param n_markers Number of top markers per cluster
#' @param min_cluster_size Subsample larger clusters to this size for power matching
#' @return List with marker analysis results
marker_validation <- function(seurat_obj, label_col, n_markers = 10, min_cluster_size = NULL) {

  Idents(seurat_obj) <- label_col

  # Power-matched analysis: subsample to smallest cluster size if requested
  if (!is.null(min_cluster_size)) {
    cluster_sizes <- table(Idents(seurat_obj))
    target_size <- min(min_cluster_size, min(cluster_sizes))

    cells_to_keep <- unlist(lapply(names(cluster_sizes), function(cl) {
      cl_cells <- WhichCells(seurat_obj, idents = cl)
      if (length(cl_cells) > target_size) {
        sample(cl_cells, target_size)
      } else {
        cl_cells
      }
    }))
    seurat_subset <- subset(seurat_obj, cells = cells_to_keep)
    Idents(seurat_subset) <- label_col
  } else {
    seurat_subset <- seurat_obj
  }

  # Find markers
  markers <- FindAllMarkers(seurat_subset,
                            only.pos = TRUE,
                            min.pct = 0.25,
                            logfc.threshold = 0.25,  # Lower threshold to capture more
                            test.use = "wilcox",
                            verbose = FALSE)

  if (nrow(markers) == 0) {
    return(list(marker_summary = NULL, error = "No significant markers found"))
  }

  # Summary per cluster at multiple thresholds
  marker_summary <- markers %>%
    group_by(cluster) %>%
    summarise(
      n_significant = n(),
      n_logfc_0.5 = sum(avg_log2FC > 0.5),
      n_logfc_1.0 = sum(avg_log2FC > 1.0),
      n_logfc_1.5 = sum(avg_log2FC > 1.5),
      top_marker = gene[which.min(p_val_adj)],
      top_marker_logfc = avg_log2FC[which.min(p_val_adj)],
      mean_logfc = mean(avg_log2FC),
      .groups = "drop"
    )

  # Top markers per cluster
  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_min(p_val_adj, n = n_markers) %>%
    ungroup()

  list(
    all_markers = markers,
    top_markers = top_markers,
    marker_summary = marker_summary,
    power_matched = !is.null(min_cluster_size),
    target_size = if (!is.null(min_cluster_size)) target_size else NA
  )
}

#' Silhouette Analysis with Stratified Sampling and Bootstrap CIs
#'
#' @param embeddings Numeric matrix
#' @param labels Character/factor vector
#' @param n_per_cluster Cells to sample per cluster
#' @param n_boot Number of bootstrap iterations
#' @param seed Random seed
#' @return List with silhouette scores and CIs
silhouette_analysis <- function(embeddings, labels, n_per_cluster = 500,
                                n_boot = 100, seed = 42) {

  valid <- !is.na(labels) & labels != "" & labels != "Unassigned"
  emb <- embeddings[valid, , drop = FALSE]
  labels_clean <- factor(labels[valid])
  label_names <- levels(labels_clean)

  if (length(unique(labels_clean)) < 2) {
    return(list(mean_silhouette = NA, error = "< 2 classes"))
  }

  # Stratified sampling function
  stratified_sample <- function(emb, labels, n_per_cluster, seed_offset = 0) {
    set.seed(seed + seed_offset)
    idx_list <- lapply(levels(labels), function(cl) {
      cl_idx <- which(labels == cl)
      if (length(cl_idx) > n_per_cluster) {
        sample(cl_idx, n_per_cluster)
      } else {
        cl_idx
      }
    })
    idx <- unlist(idx_list)
    list(emb = emb[idx, ], labels = labels[idx])
  }

  # Single silhouette computation
  compute_silhouette <- function(emb, labels) {
    labels_int <- as.integer(factor(labels))
    d <- dist(emb)
    sil <- silhouette(labels_int, d)

    overall <- mean(sil[, "sil_width"])
    per_cluster <- tapply(sil[, "sil_width"], labels, mean)

    list(overall = overall, per_cluster = per_cluster)
  }

  # Primary estimate with stratified sample
  primary_sample <- stratified_sample(emb, labels_clean, n_per_cluster, seed_offset = 0)
  primary_result <- compute_silhouette(primary_sample$emb, primary_sample$labels)

  # Bootstrap CIs
  message(sprintf("    Running %d bootstrap iterations for silhouette CIs...", n_boot))
  boot_overall <- numeric(n_boot)
  boot_per_cluster <- matrix(NA, nrow = n_boot, ncol = length(label_names))
  colnames(boot_per_cluster) <- label_names

  for (b in 1:n_boot) {
    boot_sample <- stratified_sample(emb, labels_clean, n_per_cluster, seed_offset = b)
    boot_result <- compute_silhouette(boot_sample$emb, boot_sample$labels)
    boot_overall[b] <- boot_result$overall
    boot_per_cluster[b, names(boot_result$per_cluster)] <- boot_result$per_cluster
  }

  list(
    mean_silhouette = primary_result$overall,
    silhouette_ci = quantile(boot_overall, c(0.025, 0.975)),
    per_cluster_silhouette = primary_result$per_cluster,
    per_cluster_ci_lower = apply(boot_per_cluster, 2, quantile, 0.025, na.rm = TRUE),
    per_cluster_ci_upper = apply(boot_per_cluster, 2, quantile, 0.975, na.rm = TRUE),
    n_per_cluster = n_per_cluster,
    n_boot = n_boot
  )
}

#' Analyze Confusion Matrix Patterns (Normalized by Actual Class)
#'
#' @param conf_mat Confusion matrix (table) with Predicted as rows, Actual as cols
#' @return List with misclassification analysis
analyze_confusion <- function(conf_mat) {

  # Normalize by COLUMN (actual class) to get recall/sensitivity
  # This answers: "Of cells actually in class X, what % were predicted as Y?"
  col_sums <- colSums(conf_mat)
  conf_norm <- sweep(conf_mat, 2, col_sums, "/")
  conf_norm[is.nan(conf_norm)] <- 0

  # Per-class recall (diagonal of column-normalized matrix)
  per_class_recall <- diag(conf_norm)

  # Find major misclassifications (off-diagonal > 10%)
  major_confusions <- list()
  for (i in 1:nrow(conf_norm)) {
    for (j in 1:ncol(conf_norm)) {
      if (i != j && conf_norm[i, j] > 0.10) {
        major_confusions[[length(major_confusions) + 1]] <- list(
          actual = colnames(conf_norm)[j],
          predicted_as = rownames(conf_norm)[i],
          proportion = conf_norm[i, j]
        )
      }
    }
  }

  # Identify clusters that are frequently confused with each other (bidirectional)
  confusion_pairs <- data.frame(
    cluster_a = character(),
    cluster_b = character(),
    a_as_b = numeric(),
    b_as_a = numeric(),
    bidirectional_confusion = numeric(),
    stringsAsFactors = FALSE
  )

  classes <- colnames(conf_norm)
  for (i in 1:(length(classes) - 1)) {
    for (j in (i + 1):length(classes)) {
      # % of class i predicted as class j
      i_as_j <- conf_norm[classes[j], classes[i]]
      # % of class j predicted as class i
      j_as_i <- conf_norm[classes[i], classes[j]]

      if ((i_as_j + j_as_i) > 0.15) {
        confusion_pairs <- rbind(confusion_pairs, data.frame(
          cluster_a = classes[i],
          cluster_b = classes[j],
          a_as_b = i_as_j,
          b_as_a = j_as_i,
          bidirectional_confusion = i_as_j + j_as_i
        ))
      }
    }
  }

  list(
    normalized_confusion = conf_norm,
    per_class_recall = per_class_recall,
    major_confusions = major_confusions,
    confusion_pairs = confusion_pairs
  )
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

objects_dir <- get_path(config, config$paths$objects_dir)
tables_dir <- get_path(config, config$paths$tables_dir)
fig_dir <- get_path(config, config$paths$figures_supp_dir)
ensure_dir(tables_dir)
ensure_dir(fig_dir)

results <- list()

message("\n========================================")
message("Loading GSE215932 Dataset")
message("========================================\n")

gse_path <- file.path(objects_dir, "01_seurat_annotated_snrnaseq.rds")
if (!file.exists(gse_path)) gse_path <- get_path(config, config$paths$snrnaseq_processed)
if (!file.exists(gse_path)) stop("GSE215932 not found")

seurat_gse <- readRDS(gse_path)
seurat_gse <- UpdateSeuratObject(seurat_gse)

if (!"cell_type" %in% colnames(seurat_gse@meta.data)) {
  stop("'cell_type' column not found in GSE215932")
}

# Filter to valid cells
gse_valid <- !is.na(seurat_gse$cell_type) & seurat_gse$cell_type != ""
if ("is_epithelial" %in% colnames(seurat_gse@meta.data)) {
  gse_valid <- gse_valid & seurat_gse$is_epithelial == TRUE
}
seurat_gse <- subset(seurat_gse, cells = colnames(seurat_gse)[gse_valid])

# Join layers if Seurat v5
if (inherits(seurat_gse[["RNA"]], "Assay5") && length(Layers(seurat_gse[["RNA"]])) > 1) {
  seurat_gse <- JoinLayers(seurat_gse)
}

message(sprintf("  Cells: %d", ncol(seurat_gse)))
message(sprintf("  Labels: %s", paste(unique(seurat_gse$cell_type), collapse = ", ")))
message("\n  Class distribution:")
class_counts <- table(seurat_gse$cell_type)
print(class_counts)
message(sprintf("\n  Note: Smallest cluster has %d cells; largest has %d cells (%.1fx difference)",
                min(class_counts), max(class_counts), max(class_counts)/min(class_counts)))

# Ensure PCA exists
if (!"pca" %in% Reductions(seurat_gse)) {
  message("\n  Running PCA...")
  seurat_gse <- NormalizeData(seurat_gse, verbose = FALSE)
  seurat_gse <- FindVariableFeatures(seurat_gse, nfeatures = 3000, verbose = FALSE)
  seurat_gse <- ScaleData(seurat_gse, verbose = FALSE)
  seurat_gse <- RunPCA(seurat_gse, npcs = 50, verbose = FALSE)
}

# ==============================================================================
# ANALYSIS 1: KNN CROSS-VALIDATION WITH MULTIPLE METRICS
# ==============================================================================

message("\n========================================")
message("Analysis 1: KNN Cross-Validation")
message("========================================\n")

pca_emb <- Embeddings(seurat_gse, "pca")
labels <- seurat_gse$cell_type

# Primary analysis with default parameters
message("Running primary KNN CV (k=15, dims=30)...")
cv_primary <- knn_cv(pca_emb[, 1:30], labels, k = 15, seed = config$reproducibility$seed)

message(sprintf("\n  Overall accuracy: %.1f%% ± %.1f%%",
                cv_primary$overall_accuracy * 100, cv_primary$overall_sd * 100))
message(sprintf("  Balanced accuracy: %.1f%% ± %.1f%%",
                cv_primary$balanced_accuracy * 100, cv_primary$balanced_sd * 100))
message("\n  Per-class recall (sensitivity):")
for (cls in names(cv_primary$per_class_accuracy)) {
  message(sprintf("    %s: %.1f%% (n=%d)",
                  cls, cv_primary$per_class_accuracy[cls] * 100,
                  cv_primary$class_counts[cls]))
}

results$primary_cv <- cv_primary

# ==============================================================================
# ANALYSIS 2: PERMUTATION TEST WITH FULL CV
# ==============================================================================

message("\n========================================")
message("Analysis 2: Permutation Test (Full CV, Proper Null)")
message("========================================\n")

# Note: Using 50 permutations with full 5-fold CV each for tractability
# This is computationally expensive but methodologically correct
perm_results <- permutation_test(pca_emb[, 1:30], labels, k = 15, n_folds = 5,
                                 n_perm = 50, seed = config$reproducibility$seed)

message(sprintf("\n  Null overall accuracy: %.1f%% ± %.1f%% [%.1f%%, %.1f%%]",
                perm_results$null_overall_mean * 100, perm_results$null_overall_sd * 100,
                perm_results$null_overall_95_ci[1] * 100, perm_results$null_overall_95_ci[2] * 100))
message(sprintf("  Null balanced accuracy: %.1f%% ± %.1f%% [%.1f%%, %.1f%%]",
                perm_results$null_balanced_mean * 100, perm_results$null_balanced_sd * 100,
                perm_results$null_balanced_95_ci[1] * 100, perm_results$null_balanced_95_ci[2] * 100))

# Statistical significance using proper null
p_value_overall <- mean(perm_results$null_overall_distribution >= cv_primary$overall_accuracy)
p_value_balanced <- mean(perm_results$null_balanced_distribution >= cv_primary$balanced_accuracy)
effect_size_overall <- (cv_primary$overall_accuracy - perm_results$null_overall_mean) / perm_results$null_overall_sd
effect_size_balanced <- (cv_primary$balanced_accuracy - perm_results$null_balanced_mean) / perm_results$null_balanced_sd

message(sprintf("\n  Overall: observed=%.1f%% vs null=%.1f%%, d=%.2f, p=%.4f",
                cv_primary$overall_accuracy * 100, perm_results$null_overall_mean * 100,
                effect_size_overall, p_value_overall))
message(sprintf("  Balanced: observed=%.1f%% vs null=%.1f%%, d=%.2f, p=%.4f",
                cv_primary$balanced_accuracy * 100, perm_results$null_balanced_mean * 100,
                effect_size_balanced, p_value_balanced))

if (p_value_balanced < 0.01 && effect_size_balanced > 2) {
  message("\n  → Labels capture significantly more structure than chance (p<0.01, d>2)")
} else if (p_value_balanced < 0.05) {
  message("\n  → Labels capture more structure than chance (p<0.05), but effect size is modest")
} else {
  message("\n  → WARNING: Labels may not capture substantially more structure than chance")
}

results$permutation <- perm_results
results$permutation$p_value_overall <- p_value_overall
results$permutation$p_value_balanced <- p_value_balanced
results$permutation$effect_size_overall <- effect_size_overall
results$permutation$effect_size_balanced <- effect_size_balanced

# ==============================================================================
# ANALYSIS 3: SENSITIVITY ANALYSIS
# ==============================================================================

message("\n========================================")
message("Analysis 3: Sensitivity Analysis")
message("========================================\n")

message("Testing k = {5, 15, 30, 50}, dims = {10, 20, 30, 50}...")
sensitivity <- sensitivity_analysis(pca_emb, labels,
                                    k_values = c(5, 15, 30, 50),
                                    dim_values = c(10, 20, 30, 50),
                                    seed = config$reproducibility$seed)

message("\n  Results (Overall / Balanced Accuracy %):")
sens_summary <- sensitivity %>%
  mutate(result = sprintf("%.1f / %.1f", overall_accuracy * 100, balanced_accuracy * 100)) %>%
  select(k, dims, result) %>%
  pivot_wider(names_from = dims, values_from = result, names_prefix = "dims_")
print(as.data.frame(sens_summary))

# Check stability
overall_range <- diff(range(sensitivity$overall_accuracy))
balanced_range <- diff(range(sensitivity$balanced_accuracy))
message(sprintf("\n  Overall accuracy range: %.1f%%", overall_range * 100))
message(sprintf("  Balanced accuracy range: %.1f%%", balanced_range * 100))

if (balanced_range < 0.05) {
  message("  → Results are STABLE across parameter choices")
} else if (balanced_range < 0.10) {
  message("  → Results show MODERATE sensitivity to parameters")
} else {
  message("  → WARNING: Results are SENSITIVE to parameter choices")
}

results$sensitivity <- sensitivity

# ==============================================================================
# ANALYSIS 4: CONFUSION MATRIX ANALYSIS (NORMALIZED BY ACTUAL CLASS)
# ==============================================================================

message("\n========================================")
message("Analysis 4: Confusion Matrix Analysis")
message("========================================\n")

confusion_analysis <- analyze_confusion(cv_primary$confusion_matrix)

message("Normalized confusion matrix (% of ACTUAL class predicted as each label):")
message("  Rows = Predicted, Columns = Actual")
message("  Diagonal = Recall/Sensitivity for each class")
print(round(confusion_analysis$normalized_confusion * 100, 1))

message("\n  Per-class recall (diagonal):")
for (cls in names(confusion_analysis$per_class_recall)) {
  recall <- confusion_analysis$per_class_recall[cls]
  interpretation <- if (recall >= 0.80) "good" else if (recall >= 0.50) "moderate" else "poor"
  message(sprintf("    %s: %.1f%% (%s)", cls, recall * 100, interpretation))
}

if (length(confusion_analysis$major_confusions) > 0) {
  message("\n  Major misclassifications (>10% of actual class):")
  for (mc in confusion_analysis$major_confusions) {
    message(sprintf("    %.1f%% of %s predicted as %s",
                    mc$proportion * 100, mc$actual, mc$predicted_as))
  }
}

if (nrow(confusion_analysis$confusion_pairs) > 0) {
  message("\n  Bidirectionally confused cluster pairs (>15% combined):")
  confusion_analysis$confusion_pairs <- confusion_analysis$confusion_pairs %>%
    arrange(desc(bidirectional_confusion))
  for (i in 1:nrow(confusion_analysis$confusion_pairs)) {
    cp <- confusion_analysis$confusion_pairs[i, ]
    message(sprintf("    %s <-> %s: %.1f%% + %.1f%% = %.1f%% combined",
                    cp$cluster_a, cp$cluster_b,
                    cp$a_as_b * 100, cp$b_as_a * 100,
                    cp$bidirectional_confusion * 100))
  }
}

results$confusion <- confusion_analysis

# ==============================================================================
# ANALYSIS 5: SILHOUETTE ANALYSIS WITH STRATIFIED SAMPLING AND BOOTSTRAP CIs
# ==============================================================================

message("\n========================================")
message("Analysis 5: Silhouette Analysis (Stratified + Bootstrap)")
message("========================================\n")

message("  Using stratified sampling (500 cells/cluster) with 100 bootstrap iterations...")
sil_results <- silhouette_analysis(pca_emb[, 1:30], labels,
                                   n_per_cluster = 500, n_boot = 100,
                                   seed = config$reproducibility$seed)

message(sprintf("\n  Mean silhouette: %.3f [%.3f, %.3f]",
                sil_results$mean_silhouette,
                sil_results$silhouette_ci[1], sil_results$silhouette_ci[2]))

message("\n  Per-cluster silhouette with 95% CI:")
for (cls in names(sil_results$per_cluster_silhouette)) {
  sil_val <- sil_results$per_cluster_silhouette[cls]
  ci_lo <- sil_results$per_cluster_ci_lower[cls]
  ci_hi <- sil_results$per_cluster_ci_upper[cls]
  interpretation <- if (sil_val > 0.5) "strong" else if (sil_val > 0.25) "moderate" else if (sil_val > 0) "weak" else "poor"
  message(sprintf("    %s: %.3f [%.3f, %.3f] (%s)", cls, sil_val, ci_lo, ci_hi, interpretation))
}

message("\n  Interpretation guide:")
message("    > 0.5: Strong cluster structure")
message("    0.25-0.5: Moderate structure")
message("    0-0.25: Weak structure (overlapping clusters)")
message("    < 0: Cells closer to other clusters on average")

results$silhouette <- sil_results

# ==============================================================================
# ANALYSIS 6: MARKER GENE VALIDATION WITH MULTIPLE THRESHOLDS
# ==============================================================================

message("\n========================================")
message("Analysis 6: Marker Gene Validation")
message("========================================\n")

# Run both full and power-matched analysis
message("  Finding markers (full dataset)...")
marker_results_full <- marker_validation(seurat_gse, "cell_type", n_markers = 10,
                                         min_cluster_size = NULL)

smallest_cluster <- min(class_counts)
message(sprintf("\n  Finding markers (power-matched to n=%d per cluster)...", smallest_cluster))
marker_results_matched <- marker_validation(seurat_gse, "cell_type", n_markers = 10,
                                            min_cluster_size = smallest_cluster)

message("\n  Marker summary - FULL dataset:")
if (!is.null(marker_results_full$marker_summary)) {
  print(as.data.frame(marker_results_full$marker_summary %>%
                        select(cluster, n_significant, n_logfc_0.5, n_logfc_1.0, n_logfc_1.5, top_marker, top_marker_logfc)))
}

message("\n  Marker summary - POWER-MATCHED (n=", smallest_cluster, " per cluster):")
if (!is.null(marker_results_matched$marker_summary)) {
  print(as.data.frame(marker_results_matched$marker_summary %>%
                        select(cluster, n_significant, n_logfc_0.5, n_logfc_1.0, n_logfc_1.5, top_marker, top_marker_logfc)))
}

message("\n  Note: Power-matched analysis controls for differential statistical power")
message("        due to unequal cluster sizes. Compare the two to assess bias.")

# Check if conclusions change between full and matched
if (!is.null(marker_results_full$marker_summary) && !is.null(marker_results_matched$marker_summary)) {
  full_markers <- marker_results_full$marker_summary %>% select(cluster, n_logfc_1.0) %>% rename(full = n_logfc_1.0)
  matched_markers <- marker_results_matched$marker_summary %>% select(cluster, n_logfc_1.0) %>% rename(matched = n_logfc_1.0)
  comparison <- full_markers %>% left_join(matched_markers, by = "cluster") %>%
    mutate(ratio = ifelse(matched > 0, full / matched, NA))

  message("\n  Comparison (logFC > 1.0 markers):")
  message("    Cluster       Full  Matched  Ratio")
  for (i in 1:nrow(comparison)) {
    message(sprintf("    %-12s %5d  %5d    %.1f",
                    comparison$cluster[i], comparison$full[i],
                    comparison$matched[i], comparison$ratio[i]))
  }
}

results$markers_full <- marker_results_full
results$markers_matched <- marker_results_matched

# ==============================================================================
# ANALYSIS 7: ALTERNATIVE EMBEDDING (UMAP)
# ==============================================================================

message("\n========================================")
message("Analysis 7: Alternative Embedding (UMAP)")
message("========================================\n")

if (!"umap" %in% Reductions(seurat_gse)) {
  message("Running UMAP...")
  seurat_gse <- RunUMAP(seurat_gse, dims = 1:30, verbose = FALSE)
}

umap_emb <- Embeddings(seurat_gse, "umap")

message("Running KNN CV on UMAP embedding...")
cv_umap <- knn_cv(umap_emb, labels, k = 15, seed = config$reproducibility$seed)

message(sprintf("  UMAP: overall=%.1f%%, balanced=%.1f%%",
                cv_umap$overall_accuracy * 100, cv_umap$balanced_accuracy * 100))
message(sprintf("  PCA:  overall=%.1f%%, balanced=%.1f%%",
                cv_primary$overall_accuracy * 100, cv_primary$balanced_accuracy * 100))

overall_diff <- abs(cv_umap$overall_accuracy - cv_primary$overall_accuracy)
balanced_diff <- abs(cv_umap$balanced_accuracy - cv_primary$balanced_accuracy)

if (max(overall_diff, balanced_diff) < 0.05) {
  message("  → Results CONSISTENT across embeddings (diff < 5%)")
} else {
  message(sprintf("  → Results DIFFER across embeddings (overall diff=%.1f%%, balanced diff=%.1f%%)",
                  overall_diff * 100, balanced_diff * 100))
}

results$umap_cv <- cv_umap

# ==============================================================================
# ANALYSIS 8: CROSS-DATASET TRANSFER TO ACP_SCN
# ==============================================================================

message("\n========================================")
message("Analysis 8: Cross-Dataset Transfer to ACP_SCN")
message("========================================\n")

# Load ACP_SCN dataset
acp_path <- file.path(objects_dir, "01_seurat_annotated_acp_scn.rds")
if (!file.exists(acp_path)) acp_path <- get_path(config, config$paths$acp_scn_annotated)

if (file.exists(acp_path)) {
  message(sprintf("Loading ACP_SCN: %s", acp_path))
  seurat_acp <- readRDS(acp_path)
  seurat_acp <- UpdateSeuratObject(seurat_acp)

  # Filter ACP to epithelial cells
  if ("celltypes" %in% colnames(seurat_acp@meta.data)) {
    acp_valid <- seurat_acp$celltypes == "Epithelial"
  } else {
    acp_valid <- rep(TRUE, ncol(seurat_acp))
  }
  if ("module_score_subtype" %in% colnames(seurat_acp@meta.data)) {
    acp_valid <- acp_valid & !seurat_acp$module_score_subtype %in% c("Unassigned", "Non-epithelial")
  }
  acp_valid[is.na(acp_valid)] <- FALSE
  seurat_acp <- subset(seurat_acp, cells = colnames(seurat_acp)[acp_valid])

  # Join layers if Seurat v5
  if (inherits(seurat_acp[["RNA"]], "Assay5") && length(Layers(seurat_acp[["RNA"]])) > 1) {
    seurat_acp <- JoinLayers(seurat_acp)
  }

  message(sprintf("  ACP cells after filtering: %d", ncol(seurat_acp)))

  # Add dataset identifiers
  seurat_gse$dataset <- "GSE215932"
  seurat_acp$dataset <- "ACP_SCN"

  # Find common features
  common_features <- intersect(rownames(seurat_gse), rownames(seurat_acp))
  message(sprintf("  Common features: %d", length(common_features)))

  # Merge datasets
  message("  Merging and integrating with Harmony...")
  seurat_merged <- merge(seurat_gse, seurat_acp, add.cell.ids = c("GSE", "ACP"))

  # Process and integrate
  seurat_merged <- NormalizeData(seurat_merged, verbose = FALSE)
  seurat_merged <- FindVariableFeatures(seurat_merged, nfeatures = 3000, verbose = FALSE)
  seurat_merged <- ScaleData(seurat_merged, verbose = FALSE)
  seurat_merged <- RunPCA(seurat_merged, npcs = 30, verbose = FALSE)
  seurat_merged <- RunHarmony(seurat_merged, group.by.vars = "dataset", dims.use = 1:30, verbose = FALSE)

  harmony_emb <- Embeddings(seurat_merged, "harmony")[, 1:30]
  message("  Integration complete")

  # Identify cells by dataset
  gse_cells <- colnames(seurat_merged)[seurat_merged$dataset == "GSE215932"]
  acp_cells <- colnames(seurat_merged)[seurat_merged$dataset == "ACP_SCN"]

  # Get GSE labels
  gse_labels_merged <- seurat_merged$cell_type[gse_cells]

  # Transfer labels using KNN
  message("  Transferring GSE cell_type labels to ACP cells via KNN...")

  transfer_labels_knn <- function(query_emb, ref_emb, ref_labels, k = 30) {
    nn <- nn2(ref_emb, query_emb, k = k)

    transferred <- sapply(1:nrow(query_emb), function(i) {
      neighbor_labels <- ref_labels[nn$nn.idx[i, ]]
      weights <- 1 / (nn$nn.dists[i, ] + 0.001)
      label_votes <- tapply(weights, neighbor_labels, sum)
      names(which.max(label_votes))
    })

    # Also compute confidence scores
    confidence <- sapply(1:nrow(query_emb), function(i) {
      neighbor_labels <- ref_labels[nn$nn.idx[i, ]]
      weights <- 1 / (nn$nn.dists[i, ] + 0.001)
      label_votes <- tapply(weights, neighbor_labels, sum)
      max(label_votes) / sum(label_votes)
    })

    list(labels = transferred, confidence = confidence)
  }

  transfer_result <- transfer_labels_knn(
    harmony_emb[acp_cells, ],
    harmony_emb[gse_cells, ],
    gse_labels_merged,
    k = 30
  )

  acp_transferred_labels <- transfer_result$labels
  acp_transfer_confidence <- transfer_result$confidence

  message("\n  Transfer results:")
  message("  Label distribution on ACP:")
  print(table(acp_transferred_labels))

  message(sprintf("\n  Mean transfer confidence: %.3f", mean(acp_transfer_confidence)))
  message(sprintf("  Proportion high confidence (>0.5): %.1f%%", mean(acp_transfer_confidence > 0.5) * 100))

  # -----------------------------------------------------------------------------
  # Test 8a: KNN CV on transferred labels (ACP cells in Harmony space)
  # -----------------------------------------------------------------------------
  message("\n  Test 8a: KNN CV on ACP cells with transferred labels...")

  acp_harmony_emb <- harmony_emb[acp_cells, ]
  cv_acp_transferred <- knn_cv(acp_harmony_emb, acp_transferred_labels, k = 15,
                               seed = config$reproducibility$seed)

  message(sprintf("    Overall accuracy: %.1f%% ± %.1f%%",
                  cv_acp_transferred$overall_accuracy * 100, cv_acp_transferred$overall_sd * 100))
  message(sprintf("    Balanced accuracy: %.1f%% ± %.1f%%",
                  cv_acp_transferred$balanced_accuracy * 100, cv_acp_transferred$balanced_sd * 100))

  message("\n    Per-class recall on ACP:")
  for (cls in names(cv_acp_transferred$per_class_accuracy)) {
    n_cls <- cv_acp_transferred$class_counts[cls]
    message(sprintf("      %s: %.1f%% (n=%d)", cls,
                    cv_acp_transferred$per_class_accuracy[cls] * 100, n_cls))
  }

  # -----------------------------------------------------------------------------
  # Test 8b: KNN CV on merged dataset (all cells in Harmony space)
  # -----------------------------------------------------------------------------
  message("\n  Test 8b: KNN CV on merged dataset (GSE native + ACP transferred)...")

  # Combine labels
  all_labels <- character(ncol(seurat_merged))
  all_labels[match(gse_cells, colnames(seurat_merged))] <- gse_labels_merged
  all_labels[match(acp_cells, colnames(seurat_merged))] <- acp_transferred_labels

  cv_merged <- knn_cv(harmony_emb, all_labels, k = 15, seed = config$reproducibility$seed)

  message(sprintf("    Overall accuracy: %.1f%% ± %.1f%%",
                  cv_merged$overall_accuracy * 100, cv_merged$overall_sd * 100))
  message(sprintf("    Balanced accuracy: %.1f%% ± %.1f%%",
                  cv_merged$balanced_accuracy * 100, cv_merged$balanced_sd * 100))

  # -----------------------------------------------------------------------------
  # Test 8c: Silhouette on ACP transferred labels
  # -----------------------------------------------------------------------------
  message("\n  Test 8c: Silhouette analysis on ACP with transferred labels...")

  # Determine n_per_cluster based on ACP class sizes
  acp_class_counts <- table(acp_transferred_labels)
  acp_n_per_cluster <- min(200, min(acp_class_counts))

  sil_acp <- silhouette_analysis(acp_harmony_emb, acp_transferred_labels,
                                 n_per_cluster = acp_n_per_cluster, n_boot = 50,
                                 seed = config$reproducibility$seed)

  message(sprintf("    Mean silhouette: %.3f [%.3f, %.3f]",
                  sil_acp$mean_silhouette, sil_acp$silhouette_ci[1], sil_acp$silhouette_ci[2]))

  message("\n    Per-cluster silhouette on ACP:")
  for (cls in names(sil_acp$per_cluster_silhouette)) {
    message(sprintf("      %s: %.3f [%.3f, %.3f]", cls,
                    sil_acp$per_cluster_silhouette[cls],
                    sil_acp$per_cluster_ci_lower[cls],
                    sil_acp$per_cluster_ci_upper[cls]))
  }

  # -----------------------------------------------------------------------------
  # Test 8d: Compare GSE native vs ACP transferred performance
  # -----------------------------------------------------------------------------
  message("\n  Test 8d: Cross-dataset comparison...")

  # GSE in native PCA space vs ACP in Harmony space
  # Note: Not directly comparable due to different embeddings, but informative

  # GSE in Harmony space (for fair comparison)
  gse_harmony_emb <- harmony_emb[gse_cells, ]
  cv_gse_harmony <- knn_cv(gse_harmony_emb, gse_labels_merged, k = 15,
                           seed = config$reproducibility$seed)

  comparison_df <- data.frame(
    Context = c("GSE (native PCA)", "GSE (Harmony)", "ACP (transferred)", "Merged (Harmony)"),
    Overall_Accuracy = c(cv_primary$overall_accuracy, cv_gse_harmony$overall_accuracy,
                         cv_acp_transferred$overall_accuracy, cv_merged$overall_accuracy) * 100,
    Balanced_Accuracy = c(cv_primary$balanced_accuracy, cv_gse_harmony$balanced_accuracy,
                          cv_acp_transferred$balanced_accuracy, cv_merged$balanced_accuracy) * 100,
    N_Cells = c(cv_primary$n_cells, cv_gse_harmony$n_cells,
                cv_acp_transferred$n_cells, cv_merged$n_cells)
  )

  message("\n  Cross-dataset performance comparison:")
  print(comparison_df)

  # Assess transfer quality
  gse_harmony_balanced <- cv_gse_harmony$balanced_accuracy
  acp_transferred_balanced <- cv_acp_transferred$balanced_accuracy
  transfer_drop <- (gse_harmony_balanced - acp_transferred_balanced) / gse_harmony_balanced * 100

  message(sprintf("\n  Transfer quality assessment:"))
  message(sprintf("    GSE balanced accuracy (Harmony): %.1f%%", gse_harmony_balanced * 100))
  message(sprintf("    ACP balanced accuracy (transferred): %.1f%%", acp_transferred_balanced * 100))
  message(sprintf("    Relative drop: %.1f%%", transfer_drop))

  if (transfer_drop < 10) {
    message("    → Labels transfer WELL across datasets (<10% drop)")
    transfer_quality <- "GOOD"
  } else if (transfer_drop < 25) {
    message("    → Labels transfer MODERATELY across datasets (10-25% drop)")
    transfer_quality <- "MODERATE"
  } else {
    message("    → Labels transfer POORLY across datasets (>25% drop)")
    transfer_quality <- "POOR"
  }

  # Per-class transfer assessment
  message("\n  Per-class transfer assessment:")
  message("    Cluster          GSE(Harmony)  ACP(transferred)  Drop")

  per_class_transfer <- data.frame(
    Cluster = names(cv_gse_harmony$per_class_accuracy),
    GSE_Recall = cv_gse_harmony$per_class_accuracy * 100,
    ACP_Recall = cv_acp_transferred$per_class_accuracy[names(cv_gse_harmony$per_class_accuracy)] * 100,
    stringsAsFactors = FALSE
  )
  per_class_transfer$Drop <- per_class_transfer$GSE_Recall - per_class_transfer$ACP_Recall
  per_class_transfer$Transfer_Quality <- ifelse(per_class_transfer$Drop < 10, "Good",
                                                ifelse(per_class_transfer$Drop < 25, "Moderate", "Poor"))

  for (i in 1:nrow(per_class_transfer)) {
    message(sprintf("    %-15s %6.1f%%       %6.1f%%           %+.1f%% (%s)",
                    per_class_transfer$Cluster[i],
                    per_class_transfer$GSE_Recall[i],
                    per_class_transfer$ACP_Recall[i],
                    -per_class_transfer$Drop[i],
                    per_class_transfer$Transfer_Quality[i]))
  }

  # Store results
  results$cross_dataset <- list(
    acp_n_cells = ncol(seurat_acp),
    common_features = length(common_features),
    transfer_confidence_mean = mean(acp_transfer_confidence),
    transfer_confidence_high_prop = mean(acp_transfer_confidence > 0.5),
    label_distribution_acp = table(acp_transferred_labels),
    cv_gse_harmony = cv_gse_harmony,
    cv_acp_transferred = cv_acp_transferred,
    cv_merged = cv_merged,
    silhouette_acp = sil_acp,
    comparison = comparison_df,
    per_class_transfer = per_class_transfer,
    transfer_quality = transfer_quality,
    transfer_drop = transfer_drop
  )

  # Save merged object
  seurat_merged$gse_cell_type <- all_labels
  seurat_merged$transfer_confidence <- NA
  seurat_merged$transfer_confidence[match(acp_cells, colnames(seurat_merged))] <- acp_transfer_confidence

  saveRDS(seurat_merged, file.path(objects_dir, "label_validation_merged.rds"))
  message(sprintf("\n  Saved merged object: %s", file.path(objects_dir, "label_validation_merged.rds")))

  cross_dataset_available <- TRUE

} else {
  message("  ACP_SCN dataset not found - skipping cross-dataset analysis")
  cross_dataset_available <- FALSE
  results$cross_dataset <- list(error = "ACP_SCN not found")
}

# ==============================================================================
# INTEGRATED SUMMARY
# ==============================================================================

message("\n========================================")
message("Integrated Summary")
message("========================================\n")

# Compile evidence with explicit interpretations
evidence <- data.frame(
  Analysis = c(
    "KNN CV (Overall)",
    "KNN CV (Balanced)",
    "Permutation test (balanced)",
    "Silhouette score",
    "Parameter sensitivity (balanced)",
    "Cross-embedding consistency",
    if (cross_dataset_available) "Cross-dataset transfer" else NULL
  ),
  Result = c(
    sprintf("%.1f%%", cv_primary$overall_accuracy * 100),
    sprintf("%.1f%%", cv_primary$balanced_accuracy * 100),
    sprintf("d=%.2f, p=%.4f", effect_size_balanced, p_value_balanced),
    sprintf("%.3f [%.3f, %.3f]", sil_results$mean_silhouette,
            sil_results$silhouette_ci[1], sil_results$silhouette_ci[2]),
    sprintf("%.1f%% range", balanced_range * 100),
    sprintf("%.1f%% diff", balanced_diff * 100),
    if (cross_dataset_available) sprintf("%.1f%% drop", results$cross_dataset$transfer_drop) else NULL
  ),
  Interpretation = c(
    if (cv_primary$overall_accuracy >= 0.80) "Good" else if (cv_primary$overall_accuracy >= 0.60) "Moderate" else "Poor",
    if (cv_primary$balanced_accuracy >= 0.80) "Good" else if (cv_primary$balanced_accuracy >= 0.60) "Moderate" else "Poor",
    if (p_value_balanced < 0.01 && effect_size_balanced > 2) "Strong signal" else if (p_value_balanced < 0.05) "Significant" else "Weak",
    if (sil_results$mean_silhouette > 0.25) "Good separation" else if (sil_results$mean_silhouette > 0) "Weak separation" else "Poor separation",
    if (balanced_range < 0.05) "Stable" else if (balanced_range < 0.10) "Moderate" else "Unstable",
    if (balanced_diff < 0.05) "Consistent" else "Inconsistent",
    if (cross_dataset_available) results$cross_dataset$transfer_quality else NULL
  ),
  stringsAsFactors = FALSE
)

message("Evidence summary:")
print(evidence)

# Per-cluster assessment with multiple thresholds
message("\nPer-cluster assessment:")
cluster_assessment <- data.frame(
  Cluster = names(cv_primary$per_class_accuracy),
  N_Cells = as.numeric(cv_primary$class_counts[names(cv_primary$per_class_accuracy)]),
  CV_Recall = round(cv_primary$per_class_accuracy * 100, 1),
  Silhouette = round(sil_results$per_cluster_silhouette[names(cv_primary$per_class_accuracy)], 3),
  Sil_CI_Lower = round(sil_results$per_cluster_ci_lower[names(cv_primary$per_class_accuracy)], 3),
  Sil_CI_Upper = round(sil_results$per_cluster_ci_upper[names(cv_primary$per_class_accuracy)], 3),
  Markers_logFC1_Full = if (!is.null(marker_results_full$marker_summary)) {
    marker_results_full$marker_summary$n_logfc_1.0[match(names(cv_primary$per_class_accuracy),
                                                         marker_results_full$marker_summary$cluster)]
  } else NA,
  Markers_logFC1_Matched = if (!is.null(marker_results_matched$marker_summary)) {
    marker_results_matched$marker_summary$n_logfc_1.0[match(names(cv_primary$per_class_accuracy),
                                                            marker_results_matched$marker_summary$cluster)]
  } else NA,
  stringsAsFactors = FALSE
)

# Assessment based on multiple criteria
# Note: We now explicitly acknowledge that these criteria measure different things
cluster_assessment$Assessment <- sapply(1:nrow(cluster_assessment), function(i) {
  recall <- cluster_assessment$CV_Recall[i]
  sil <- cluster_assessment$Silhouette[i]
  markers <- cluster_assessment$Markers_logFC1_Matched[i]  # Use power-matched

  # Count evidence for separability (not "validity" - different concept)
  evidence_count <- 0
  if (!is.na(recall) && recall >= 60) evidence_count <- evidence_count + 1
  if (!is.na(sil) && sil >= 0.05) evidence_count <- evidence_count + 1
  if (!is.na(markers) && markers >= 5) evidence_count <- evidence_count + 1

  # Conservative labeling
  if (evidence_count >= 3) "SEPARABLE"
  else if (evidence_count >= 2) "PARTIALLY SEPARABLE"
  else "OVERLAPPING"
})

print(cluster_assessment)

message("\n  Note on assessment criteria:")
message("    - CV Recall: Predictability from PCA embedding (threshold: 60%)")
message("    - Silhouette: Compactness vs separation in PCA space (threshold: 0.05)")
message("    - Markers: Differentially expressed genes with logFC>1 (threshold: 5, power-matched)")
message("    These measure TRANSCRIPTIONAL SEPARABILITY, not biological validity.")
message("    Low separability may indicate: transitional states, technical factors, or")
message("    labels defined by non-transcriptomic criteria (spatial, protein, etc.)")

results$evidence_summary <- evidence
results$cluster_assessment <- cluster_assessment

# ==============================================================================
# CONCLUSIONS
# ==============================================================================

message("\n========================================")
message("Conclusions")
message("========================================\n")

# Count by assessment
n_separable <- sum(cluster_assessment$Assessment == "SEPARABLE")
n_partial <- sum(cluster_assessment$Assessment == "PARTIALLY SEPARABLE")
n_overlapping <- sum(cluster_assessment$Assessment == "OVERLAPPING")

message("Cluster separability summary:")
message(sprintf("  SEPARABLE: %d clusters", n_separable))
message(sprintf("  PARTIALLY SEPARABLE: %d clusters", n_partial))
message(sprintf("  OVERLAPPING: %d clusters", n_overlapping))

# Key findings with appropriately hedged language
message("\nKey findings:")

# Overall assessment
if (effect_size_balanced > 2 && p_value_balanced < 0.01) {
  message("  1. Labels capture significantly more transcriptional structure than chance")
  message(sprintf("     (balanced accuracy: %.1f%% vs null %.1f%%, p<0.01)",
                  cv_primary$balanced_accuracy * 100, perm_results$null_balanced_mean * 100))
} else {
  message("  1. Labels capture some structure above chance, but effect size is modest")
}

# Cluster-specific findings
overlapping_clusters <- cluster_assessment$Cluster[cluster_assessment$Assessment == "OVERLAPPING"]
if (length(overlapping_clusters) > 0) {
  message(sprintf("\n  2. Clusters with substantial transcriptional overlap: %s",
                  paste(overlapping_clusters, collapse = ", ")))
  message("     Possible interpretations (not mutually exclusive):")
  message("       a) Transitional/intermediate cell states (biologically meaningful)")
  message("       b) Overclustering of a continuous population")
  message("       c) Labels defined by non-transcriptomic criteria")
  message("       d) Technical factors (batch effects, dropout)")
  message("     This analysis cannot distinguish between these explanations.")
}

# Confusion patterns
if (nrow(confusion_analysis$confusion_pairs) > 0) {
  top_confused <- confusion_analysis$confusion_pairs[1, ]
  message(sprintf("\n  3. Most confused pair: %s <-> %s (%.1f%% bidirectional)",
                  top_confused$cluster_a, top_confused$cluster_b,
                  top_confused$bidirectional_confusion * 100))
  message("     This could indicate:")
  message("       - A continuum between these states")
  message("       - Candidates for merging")
  message("       - Or biologically distinct states with similar transcriptomes")
}

# Cross-dataset transfer findings
if (cross_dataset_available) {
  message(sprintf("\n  4. Cross-dataset transfer quality: %s", results$cross_dataset$transfer_quality))
  message(sprintf("     GSE→ACP transfer drop: %.1f%% (balanced accuracy)", results$cross_dataset$transfer_drop))

  # Identify which clusters transfer well/poorly
  good_transfer <- results$cross_dataset$per_class_transfer$Cluster[
    results$cross_dataset$per_class_transfer$Transfer_Quality == "Good"]
  poor_transfer <- results$cross_dataset$per_class_transfer$Cluster[
    results$cross_dataset$per_class_transfer$Transfer_Quality == "Poor"]

  if (length(good_transfer) > 0) {
    message(sprintf("     Clusters that transfer well: %s", paste(good_transfer, collapse = ", ")))
  }
  if (length(poor_transfer) > 0) {
    message(sprintf("     Clusters that transfer poorly: %s", paste(poor_transfer, collapse = ", ")))
  }

  message("     Interpretation:")
  if (results$cross_dataset$transfer_quality == "GOOD") {
    message("       - Labels capture generalizable transcriptional patterns")
    message("       - Suitable for cross-dataset analysis")
  } else if (results$cross_dataset$transfer_quality == "MODERATE") {
    message("       - Labels partially transfer; some clusters are dataset-specific")
    message("       - Use caution in cross-dataset applications")
  } else {
    message("       - Labels may be dataset-specific or capture batch effects")
    message("       - Not recommended for cross-dataset transfer")
  }
}

# Limitations
message("\n  Limitations of this analysis:")
message("    - Tests transcriptional separability in PCA space, not biological validity")
message("    - PCA captures variance, which may not align with biological signal")
message("    - Original labeling criteria are unknown and may include non-transcriptomic data")
message("    - Low separability ≠ invalid labels; high separability ≠ valid labels")
message("    - Assessment thresholds are conventional, not empirically derived")

results$conclusions <- list(
  n_separable = n_separable,
  n_partial = n_partial,
  n_overlapping = n_overlapping,
  overlapping_clusters = overlapping_clusters,
  effect_size_balanced = effect_size_balanced,
  p_value_balanced = p_value_balanced,
  cross_dataset_quality = if (cross_dataset_available) results$cross_dataset$transfer_quality else NA
)

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

message("\n========================================")
message("Saving Outputs")
message("========================================\n")

# Summary tables
write.csv(evidence, file.path(tables_dir, "label_validation_evidence.csv"), row.names = FALSE)
write.csv(cluster_assessment, file.path(tables_dir, "label_validation_clusters.csv"), row.names = FALSE)
write.csv(sensitivity, file.path(tables_dir, "label_validation_sensitivity.csv"), row.names = FALSE)

if (!is.null(marker_results_full$marker_summary)) {
  write.csv(marker_results_full$marker_summary,
            file.path(tables_dir, "label_validation_markers_full.csv"), row.names = FALSE)
}
if (!is.null(marker_results_matched$marker_summary)) {
  write.csv(marker_results_matched$marker_summary,
            file.path(tables_dir, "label_validation_markers_matched.csv"), row.names = FALSE)
}

# Cross-dataset results
if (cross_dataset_available) {
  write.csv(results$cross_dataset$comparison,
            file.path(tables_dir, "label_validation_cross_dataset_comparison.csv"), row.names = FALSE)
  write.csv(results$cross_dataset$per_class_transfer,
            file.path(tables_dir, "label_validation_cross_dataset_per_class.csv"), row.names = FALSE)
}

message("Saved CSV tables")

# Full results object
saveRDS(results, file.path(objects_dir, "label_validation_rigorous.rds"))
message(sprintf("Saved: %s", file.path(objects_dir, "label_validation_rigorous.rds")))

# ==============================================================================
# GENERATE PLOTS
# ==============================================================================

message("\nGenerating plots...")

# Plot 1: Per-cluster metrics comparison
p1_data <- cluster_assessment %>%
  select(Cluster, CV_Recall, Silhouette, Markers_logFC1_Matched, Assessment) %>%
  rename(Markers = Markers_logFC1_Matched) %>%
  pivot_longer(cols = c(CV_Recall, Silhouette, Markers),
               names_to = "Metric", values_to = "Value")

# Normalize for visualization
p1_data <- p1_data %>%
  group_by(Metric) %>%
  mutate(Value_Scaled = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE) + 0.001)) %>%
  ungroup()

p1 <- ggplot(p1_data, aes(x = Cluster, y = Metric, fill = Value_Scaled)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Value, 1)), size = 3) +
  scale_fill_gradient2(low = "#E41A1C", mid = "#FFEDA0", high = "#4DAF4A",
                       midpoint = 0.5, name = "Scaled\nValue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Per-Cluster Separability Metrics",
       subtitle = "Values shown; colors normalized within each metric")

# Plot 2: Confusion matrix heatmap (normalized by actual class)
conf_df <- as.data.frame(as.table(confusion_analysis$normalized_confusion * 100))
colnames(conf_df) <- c("Predicted", "Actual", "Percentage")

p2 <- ggplot(conf_df, aes(x = Actual, y = Predicted, fill = Percentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.0f", Percentage)), size = 3) +
  scale_fill_gradient2(low = "white", mid = "#FEE08B", high = "#D73027",
                       midpoint = 25, name = "% of\nActual") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Confusion Matrix (Normalized by Actual Class)",
       subtitle = "Diagonal = Recall; Off-diagonal = Misclassification rate")

# Plot 3: Sensitivity analysis
p3 <- ggplot(sensitivity, aes(x = factor(dims), y = balanced_accuracy * 100,
                              color = factor(k), group = factor(k))) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Sensitivity Analysis (Balanced Accuracy)",
       subtitle = sprintf("Range across parameters: %.1f%%", balanced_range * 100),
       x = "PCA Dimensions", y = "Balanced Accuracy (%)", color = "k")

# Plot 4: Permutation null distribution
perm_df <- data.frame(acc = perm_results$null_balanced_distribution * 100)
p4 <- ggplot(perm_df, aes(x = acc)) +
  geom_histogram(bins = 15, fill = "gray70", color = "white") +
  geom_vline(xintercept = cv_primary$balanced_accuracy * 100,
             color = "#E41A1C", linewidth = 1.5) +
  geom_vline(xintercept = perm_results$null_balanced_mean * 100,
             color = "gray30", linetype = "dashed") +
  annotate("text", x = cv_primary$balanced_accuracy * 100, y = Inf,
           label = sprintf("Observed\n%.1f%%", cv_primary$balanced_accuracy * 100),
           vjust = 1.5, hjust = -0.1, color = "#E41A1C", size = 3) +
  theme_minimal() +
  labs(title = "Permutation Test (Full CV)",
       subtitle = sprintf("Effect size d=%.2f, p=%.4f (50 permutations × 5-fold CV)",
                          effect_size_balanced, p_value_balanced),
       x = "Balanced Accuracy (%)", y = "Count")

# Plot 5: Silhouette by cluster with CIs
sil_df <- data.frame(
  Cluster = names(sil_results$per_cluster_silhouette),
  Silhouette = sil_results$per_cluster_silhouette,
  CI_Lower = sil_results$per_cluster_ci_lower[names(sil_results$per_cluster_silhouette)],
  CI_Upper = sil_results$per_cluster_ci_upper[names(sil_results$per_cluster_silhouette)]
)

p5 <- ggplot(sil_df, aes(x = reorder(Cluster, Silhouette), y = Silhouette)) +
  geom_bar(stat = "identity", fill = ifelse(sil_df$Silhouette > 0, "#4DAF4A", "#E41A1C"),
           width = 0.6) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "gray50") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Silhouette Scores by Cluster",
       subtitle = "Error bars = 95% bootstrap CI; dashed line = 0.25 threshold",
       x = "", y = "Silhouette Score")

# Plot 6: Marker comparison (full vs power-matched)
if (!is.null(marker_results_full$marker_summary) && !is.null(marker_results_matched$marker_summary)) {
  marker_comp <- marker_results_full$marker_summary %>%
    select(cluster, n_logfc_1.0) %>%
    rename(Full = n_logfc_1.0) %>%
    left_join(marker_results_matched$marker_summary %>%
                select(cluster, n_logfc_1.0) %>%
                rename(Matched = n_logfc_1.0), by = "cluster") %>%
    pivot_longer(cols = c(Full, Matched), names_to = "Analysis", values_to = "N_Markers")

  p6 <- ggplot(marker_comp, aes(x = cluster, y = N_Markers, fill = Analysis)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Full" = "#377EB8", "Matched" = "#FF7F00")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Marker Genes (logFC > 1.0)",
         subtitle = sprintf("Power-matched analysis subsampled to n=%d per cluster", smallest_cluster),
         x = "", y = "Number of Markers")
} else {
  p6 <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Marker analysis not available")
}

# Plot 7: Cross-dataset comparison (if available)
if (cross_dataset_available) {
  p7 <- ggplot(results$cross_dataset$comparison,
               aes(x = Context, y = Balanced_Accuracy, fill = Context)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", Balanced_Accuracy)), vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("GSE (native PCA)" = "#4DAF4A",
                                 "GSE (Harmony)" = "#377EB8",
                                 "ACP (transferred)" = "#FF7F00",
                                 "Merged (Harmony)" = "#984EA3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Cross-Dataset Transfer Performance",
         subtitle = sprintf("Transfer quality: %s (%.1f%% drop)",
                            results$cross_dataset$transfer_quality,
                            results$cross_dataset$transfer_drop),
         x = "", y = "Balanced Accuracy (%)") +
    ylim(0, max(results$cross_dataset$comparison$Balanced_Accuracy) * 1.15)

  # Plot 8: Per-class transfer comparison
  transfer_plot_data <- results$cross_dataset$per_class_transfer %>%
    pivot_longer(cols = c(GSE_Recall, ACP_Recall),
                 names_to = "Dataset", values_to = "Recall") %>%
    mutate(Dataset = ifelse(Dataset == "GSE_Recall", "GSE (Harmony)", "ACP (transferred)"))

  p8 <- ggplot(transfer_plot_data, aes(x = Cluster, y = Recall, fill = Dataset)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("GSE (Harmony)" = "#377EB8", "ACP (transferred)" = "#FF7F00")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Per-Class Transfer Performance",
         subtitle = "Recall on GSE vs ACP after label transfer",
         x = "", y = "Recall (%)")
} else {
  p7 <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Cross-dataset analysis not available")
  p8 <- ggplot() + theme_void()
}

# Combine plots
if (cross_dataset_available) {
  combined_plot <- (p1 | p2) / (p3 | p4) / (p5 | p6) / (p7 | p8) +
    plot_annotation(
      title = "GSE cell_type Label Validation - Rigorous Analysis",
      subtitle = sprintf("n = %d GSE cells, %d ACP cells | Tests transcriptional separability and cross-dataset transfer",
                         ncol(seurat_gse), results$cross_dataset$acp_n_cells)
    )
  plot_height <- 18
} else {
  combined_plot <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
    plot_annotation(
      title = "GSE cell_type Label Validation - Rigorous Analysis",
      subtitle = sprintf("n = %d cells, %d clusters | Tests transcriptional separability, not biological validity",
                         ncol(seurat_gse), length(unique(labels)))
    )
  plot_height <- 14
}

ggsave(file.path(fig_dir, "label_validation_rigorous.pdf"), combined_plot,
       width = 14, height = plot_height)
ggsave(file.path(fig_dir, "label_validation_rigorous.png"), combined_plot,
       width = 14, height = plot_height, dpi = 300)
message(sprintf("Saved: %s", file.path(fig_dir, "label_validation_rigorous.pdf")))

# ==============================================================================
# COMPLETION
# ==============================================================================

message("\n")
message("================================================================")
message("  Label Validation Complete")
message("================================================================")
message(sprintf("  Finished: %s", Sys.time()))
message("")
message("Key results:")
message(sprintf("  Overall accuracy: %.1f%% (balanced: %.1f%%)",
                cv_primary$overall_accuracy * 100, cv_primary$balanced_accuracy * 100))
message(sprintf("  Permutation test: d=%.2f, p=%.4f (balanced accuracy)",
                effect_size_balanced, p_value_balanced))
message(sprintf("  Mean silhouette: %.3f [%.3f, %.3f]",
                sil_results$mean_silhouette, sil_results$silhouette_ci[1], sil_results$silhouette_ci[2]))
message(sprintf("  Separable clusters: %d/%d", n_separable, nrow(cluster_assessment)))
if (cross_dataset_available) {
  message(sprintf("  Cross-dataset transfer: %s (%.1f%% drop)",
                  results$cross_dataset$transfer_quality, results$cross_dataset$transfer_drop))
}
message("")
if (length(overlapping_clusters) > 0) {
  message(sprintf("  Overlapping clusters: %s", paste(overlapping_clusters, collapse = ", ")))
}
message("")
message("  Note: This analysis tests transcriptional separability in PCA space.")
message("        Low separability does not necessarily mean labels are invalid.")
message("")
message("Outputs:")
message(sprintf("  Tables: %s", tables_dir))
message(sprintf("  Figures: %s", fig_dir))
message(sprintf("  Full results: %s", file.path(objects_dir, "label_validation_rigorous.rds")))
message("================================================================\n")
