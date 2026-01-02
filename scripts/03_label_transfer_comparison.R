#!/usr/bin/env Rscript
# ==============================================================================
# scripts/03_integrated_trajectory_comparison.R
# ==============================================================================
# Integrated Trajectory Analysis: ACP vs Healthy Reference
#
# This script:
#   1. Integrates merged (ACP) and cellxgene (healthy) datasets using Harmony
#   2. Runs trajectory analysis on the integrated embedding
#   3. Transfers cell type labels in integrated space
#   4. Identifies outlier clusters in ACP relative to the healthy trajectory
#   5. Tests specific hypotheses about trajectory divergence
#
# Hypotheses:
#   H1: ACP cells can be mapped to healthy epithelial cell types
#   H2: Trajectory order is preserved in integrated space
#   H3: ACP contains outlier populations that deviate from healthy trajectory
#   H4: Outlier cells show distinct gene expression signatures
#
# Usage:
#   Rscript scripts/03_integrated_trajectory_comparison.R [--cores N]
#
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(slingshot)
})

options(future.globals.maxSize = 32 * 1024^3)

message("\n")
message("================================================================")
message("  Integrated Trajectory Analysis: ACP vs Healthy Reference")
message("================================================================")
message(sprintf("  Started: %s", Sys.time()))

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
n_cores <- 4
for (i in seq_along(args)) {
  if (args[i] == "--cores" && i < length(args)) {
    n_cores <- as.integer(args[i + 1])
  }
}
message(sprintf("  Using %d cores", n_cores))

# Source utilities
source("R/utils/config.R")

# Load configuration
config <- load_config()
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

# Load cellxgene (healthy reference)
ref_path <- file.path(objects_dir, "01_seurat_annotated_cellxgene.rds")
if (!file.exists(ref_path)) {
  stop("Reference dataset not found. Run 01 script with cellxgene first.")
}
message("Loading reference (cellxgene)...")
seurat_ref <- readRDS(ref_path)
seurat_ref <- UpdateSeuratObject(seurat_ref)

# Load merged (ACP query)
query_path <- file.path(objects_dir, "01_seurat_annotated_merged.rds")
if (!file.exists(query_path)) {
  stop("Query dataset not found. Run 01 script with merged first.")
}
message("Loading query (merged/ACP)...")
seurat_query <- readRDS(query_path)
seurat_query <- UpdateSeuratObject(seurat_query)

# Join layers if Seurat v5
if (inherits(seurat_ref[["RNA"]], "Assay5") && length(Layers(seurat_ref[["RNA"]])) > 1) {
  seurat_ref <- JoinLayers(seurat_ref)
}
if (inherits(seurat_query[["RNA"]], "Assay5") && length(Layers(seurat_query[["RNA"]])) > 1) {
  seurat_query <- JoinLayers(seurat_query)
}

# Add dataset identifier
seurat_ref$dataset <- "healthy"
seurat_query$dataset <- "ACP"

# Verify cell_type column in reference
if (!"cell_type" %in% colnames(seurat_ref@meta.data)) {
  stop("cell_type column not found in reference dataset")
}

message(sprintf("  Reference (healthy): %d cells", ncol(seurat_ref)))
message(sprintf("  Query (ACP): %d cells", ncol(seurat_query)))

message("\nReference cell_type distribution:")
print(table(seurat_ref$cell_type))

# ==============================================================================
# FILTER TO EPITHELIAL CELLS
# ==============================================================================

message("\n========================================")
message("Filtering to Epithelial Cells")
message("========================================\n")

# Filter reference to epithelial cell types only
epithelial_types <- c("basal cell of epidermis", "keratinocyte", "spinous cell of epidermis")
cells_keep_ref <- seurat_ref$cell_type %in% epithelial_types
seurat_ref <- subset(seurat_ref, cells = colnames(seurat_ref)[cells_keep_ref])
message(sprintf("  Reference epithelial cells: %d", ncol(seurat_ref)))

# Filter query to assigned epithelial cells (remove Unassigned)
if ("module_score_subtype" %in% colnames(seurat_query@meta.data)) {
  cells_keep_query <- seurat_query$module_score_subtype != "Unassigned"
  seurat_query <- subset(seurat_query, cells = colnames(seurat_query)[cells_keep_query])
}
message(sprintf("  Query epithelial cells: %d", ncol(seurat_query)))

# ==============================================================================
# MERGE AND INTEGRATE WITH HARMONY
# ==============================================================================

message("\n========================================")
message("Integrating Datasets with Harmony")
message("========================================\n")

# Find common features
common_features <- intersect(rownames(seurat_ref), rownames(seurat_query))
message(sprintf("  Common features: %d", length(common_features)))

# Merge datasets
message("Merging datasets...")
seurat_merged <- merge(seurat_ref, seurat_query, add.cell.ids = c("healthy", "ACP"))

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
# TRANSFER LABELS IN INTEGRATED SPACE
# ==============================================================================

message("\n========================================")
message("Transferring Labels in Integrated Space")
message("========================================\n")

# Use harmony-corrected PCA for label transfer
# Split back to get corrected embeddings
acp_cells <- colnames(seurat_merged)[seurat_merged$dataset == "ACP"]
healthy_cells <- colnames(seurat_merged)[seurat_merged$dataset == "healthy"]

# Get reference labels
ref_labels <- seurat_merged$cell_type[healthy_cells]
ref_labels <- ref_labels[!is.na(ref_labels)]

# Use KNN in harmony space for label transfer
message("Performing KNN-based label transfer...")
harmony_emb <- Embeddings(seurat_merged, "harmony")[, 1:30]

# For each ACP cell, find k nearest neighbors in healthy cells
k <- 30
acp_emb <- harmony_emb[acp_cells, ]
healthy_emb <- harmony_emb[healthy_cells, ]

# Compute distances
message("  Computing distances...")
nn_results <- RANN::nn2(healthy_emb, acp_emb, k = k)

# Vote for label based on neighbors
message("  Voting for labels...")
predicted_labels <- character(length(acp_cells))
prediction_scores <- numeric(length(acp_cells))

for (i in seq_along(acp_cells)) {
  neighbor_idx <- nn_results$nn.idx[i, ]
  neighbor_labels <- ref_labels[neighbor_idx]

  # Weighted voting by distance
  weights <- 1 / (nn_results$nn.dists[i, ] + 0.001)
  label_votes <- tapply(weights, neighbor_labels, sum)

  predicted_labels[i] <- names(which.max(label_votes))
  prediction_scores[i] <- max(label_votes) / sum(label_votes)
}

# Add to metadata
seurat_merged$predicted_cell_type <- NA
seurat_merged$predicted_cell_type[match(acp_cells, colnames(seurat_merged))] <- predicted_labels
seurat_merged$prediction_score <- NA
seurat_merged$prediction_score[match(acp_cells, colnames(seurat_merged))] <- prediction_scores

# For healthy cells, use original labels
seurat_merged$predicted_cell_type[match(healthy_cells, colnames(seurat_merged))] <-
  seurat_merged$cell_type[match(healthy_cells, colnames(seurat_merged))]
seurat_merged$prediction_score[match(healthy_cells, colnames(seurat_merged))] <- 1.0

message("\nTransferred label distribution (ACP only):")
print(table(predicted_labels))
message(sprintf("Mean prediction score: %.3f", mean(prediction_scores)))

# ==============================================================================
# TRAJECTORY ANALYSIS ON INTEGRATED DATA
# ==============================================================================

message("\n========================================")
message("Trajectory Analysis on Integrated Data")
message("========================================\n")

# Use Slingshot on harmony-corrected UMAP
message("Running Slingshot on integrated embedding...")

# Get UMAP coordinates
umap_coords <- Embeddings(seurat_merged, "umap")

# Define start cluster (basal cells)
cluster_labels <- seurat_merged$predicted_cell_type

# Run slingshot
sds <- slingshot(umap_coords,
                 clusterLabels = cluster_labels,
                 start.clus = "basal cell of epidermis",
                 stretch = 0)

# Extract pseudotime
pseudotime <- slingPseudotime(sds)
if (ncol(pseudotime) > 1) {
  # Average across lineages
  pt_consensus <- rowMeans(pseudotime, na.rm = TRUE)
} else {
  pt_consensus <- pseudotime[, 1]
}

# Normalize to 0-1
pt_consensus <- (pt_consensus - min(pt_consensus, na.rm = TRUE)) /
  (max(pt_consensus, na.rm = TRUE) - min(pt_consensus, na.rm = TRUE))

seurat_merged$pseudotime_integrated <- pt_consensus

message(sprintf("  Pseudotime range: [%.2f, %.2f]",
                min(pt_consensus, na.rm = TRUE), max(pt_consensus, na.rm = TRUE)))

# ==============================================================================
# HYPOTHESIS TESTING
# ==============================================================================

message("\n========================================")
message("Hypothesis Testing")
message("========================================\n")

results <- list()

# -----------------------------------------------------------------------------
# H1: ACP cells map to healthy epithelial types
# -----------------------------------------------------------------------------
message("--- H1: Label Transfer Quality ---")

acp_scores <- prediction_scores
h1_threshold <- 0.5
h1_prop_above <- mean(acp_scores > h1_threshold)

h1_test <- prop.test(sum(acp_scores > h1_threshold), length(acp_scores),
                     p = 0.5, alternative = "greater")

results$H1 <- list(
  hypothesis = "ACP cells can be mapped to healthy epithelial cell types",
  mean_score = mean(acp_scores),
  prop_above_threshold = h1_prop_above,
  p_value = h1_test$p.value,
  conclusion = ifelse(h1_test$p.value < 0.05 && h1_prop_above > 0.5,
                      "SUPPORTED", "NOT SUPPORTED")
)

message(sprintf("  Mean prediction score: %.3f", results$H1$mean_score))
message(sprintf("  Proportion > 0.5: %.1f%%", results$H1$prop_above_threshold * 100))
message(sprintf("  Conclusion: %s", results$H1$conclusion))

# -----------------------------------------------------------------------------
# H2: Trajectory order preserved in integrated space
# -----------------------------------------------------------------------------
message("\n--- H2: Trajectory Order Preservation ---")

expected_order <- c("basal cell of epidermis", "keratinocyte", "spinous cell of epidermis")

# Calculate mean pseudotime by predicted cell type for ACP cells
pt_by_type <- tapply(seurat_merged$pseudotime_integrated[seurat_merged$dataset == "ACP"],
                     seurat_merged$predicted_cell_type[seurat_merged$dataset == "ACP"],
                     mean, na.rm = TRUE)

observed_order <- names(sort(pt_by_type))
message(sprintf("  Expected: %s", paste(expected_order, collapse = " → ")))
message(sprintf("  Observed (ACP): %s", paste(observed_order, collapse = " → ")))

# Spearman correlation
shared_types <- intersect(expected_order, names(pt_by_type))
if (length(shared_types) >= 3) {
  exp_ranks <- match(shared_types, expected_order)
  obs_values <- pt_by_type[shared_types]
  h2_cor <- cor.test(exp_ranks, rank(obs_values), method = "spearman")

  results$H2 <- list(
    hypothesis = "Trajectory order is preserved in integrated space",
    expected_order = paste(expected_order, collapse = " → "),
    observed_order = paste(observed_order, collapse = " → "),
    spearman_rho = as.numeric(h2_cor$estimate),
    p_value = h2_cor$p.value,
    conclusion = ifelse(h2_cor$estimate > 0.5, "SUPPORTED", "NOT SUPPORTED")
  )
  message(sprintf("  Spearman rho: %.3f", results$H2$spearman_rho))
} else {
  results$H2 <- list(
    hypothesis = "Trajectory order preserved",
    conclusion = "INCONCLUSIVE: Insufficient shared cell types"
  )
}
message(sprintf("  Conclusion: %s", results$H2$conclusion))

# -----------------------------------------------------------------------------
# H3: Identify outlier clusters in ACP
# -----------------------------------------------------------------------------
message("\n--- H3: Outlier Detection in ACP ---")

# Cluster ACP cells in integrated space
message("  Clustering ACP cells...")
acp_subset <- subset(seurat_merged, dataset == "ACP")
acp_subset <- FindNeighbors(acp_subset, reduction = "harmony", dims = 1:30, verbose = FALSE)
acp_subset <- FindClusters(acp_subset, resolution = 0.5, verbose = FALSE)

# For each ACP cluster, calculate:
# 1. Mean distance to healthy cells in harmony space
# 2. Pseudotime distribution statistics
# 3. Deviation from expected trajectory

acp_clusters <- levels(Idents(acp_subset))
message(sprintf("  Found %d clusters in ACP", length(acp_clusters)))

outlier_stats <- data.frame(
  cluster = character(),
  n_cells = integer(),
  mean_pseudotime = numeric(),
  sd_pseudotime = numeric(),
  mean_dist_to_healthy = numeric(),
  dominant_celltype = character(),
  celltype_purity = numeric(),
  outlier_score = numeric(),
  stringsAsFactors = FALSE
)

# Get harmony embeddings
harmony_healthy <- harmony_emb[healthy_cells, ]

for (cl in acp_clusters) {
  cl_cells <- colnames(acp_subset)[Idents(acp_subset) == cl]
  cl_harmony <- harmony_emb[cl_cells, , drop = FALSE]

  # Mean distance to healthy centroid
  healthy_centroid <- colMeans(harmony_healthy)
  dists_to_centroid <- sqrt(rowSums((sweep(cl_harmony, 2, healthy_centroid))^2))
  mean_dist <- mean(dists_to_centroid)

  # Pseudotime stats
  cl_pt <- seurat_merged$pseudotime_integrated[match(cl_cells, colnames(seurat_merged))]

  # Dominant cell type and purity
  cl_types <- seurat_merged$predicted_cell_type[match(cl_cells, colnames(seurat_merged))]
  type_table <- table(cl_types)
  dominant_type <- names(which.max(type_table))
  purity <- max(type_table) / sum(type_table)

  # Outlier score: combines distance and purity (low purity + high distance = outlier)
  outlier_score <- mean_dist * (1 - purity)

  outlier_stats <- rbind(outlier_stats, data.frame(
    cluster = cl,
    n_cells = length(cl_cells),
    mean_pseudotime = mean(cl_pt, na.rm = TRUE),
    sd_pseudotime = sd(cl_pt, na.rm = TRUE),
    mean_dist_to_healthy = mean_dist,
    dominant_celltype = dominant_type,
    celltype_purity = purity,
    outlier_score = outlier_score
  ))
}

# Identify outliers (top 20% by outlier score)
outlier_threshold <- quantile(outlier_stats$outlier_score, 0.8)
outlier_stats$is_outlier <- outlier_stats$outlier_score > outlier_threshold

message("\n  Cluster statistics:")
print(outlier_stats[, c("cluster", "n_cells", "mean_pseudotime", "dominant_celltype",
                        "celltype_purity", "outlier_score", "is_outlier")])

results$H3 <- list(
  hypothesis = "ACP contains outlier populations that deviate from healthy trajectory",
  n_clusters = nrow(outlier_stats),
  n_outliers = sum(outlier_stats$is_outlier),
  outlier_clusters = outlier_stats$cluster[outlier_stats$is_outlier],
  outlier_threshold = outlier_threshold,
  conclusion = ifelse(sum(outlier_stats$is_outlier) > 0,
                      sprintf("SUPPORTED: %d outlier cluster(s) identified",
                              sum(outlier_stats$is_outlier)),
                      "NOT SUPPORTED: No outlier clusters identified")
)
message(sprintf("  Conclusion: %s", results$H3$conclusion))

# Add outlier status to merged object
acp_subset$is_outlier_cluster <- Idents(acp_subset) %in% results$H3$outlier_clusters
seurat_merged$is_outlier_cluster <- FALSE
seurat_merged$is_outlier_cluster[match(colnames(acp_subset), colnames(seurat_merged))] <-
  acp_subset$is_outlier_cluster

seurat_merged$acp_cluster <- NA
seurat_merged$acp_cluster[match(colnames(acp_subset), colnames(seurat_merged))] <-
  as.character(Idents(acp_subset))

# -----------------------------------------------------------------------------
# H4: Outlier cells show distinct gene signatures
# -----------------------------------------------------------------------------
message("\n--- H4: Outlier Gene Signatures ---")

if (sum(outlier_stats$is_outlier) > 0) {
  # Find markers for outlier vs non-outlier ACP cells
  message("  Finding outlier-specific markers...")

  Idents(acp_subset) <- acp_subset$is_outlier_cluster
  outlier_markers <- FindMarkers(acp_subset,
                                 ident.1 = TRUE,
                                 ident.2 = FALSE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.5,
                                 verbose = FALSE)

  # Top markers
  top_up <- head(rownames(outlier_markers[outlier_markers$avg_log2FC > 0, ]), 20)
  top_down <- head(rownames(outlier_markers[outlier_markers$avg_log2FC < 0, ]), 20)

  results$H4 <- list(
    hypothesis = "Outlier cells show distinct gene expression signatures",
    n_markers = nrow(outlier_markers),
    top_upregulated = top_up,
    top_downregulated = top_down,
    conclusion = ifelse(nrow(outlier_markers) > 10,
                        sprintf("SUPPORTED: %d DE genes identified", nrow(outlier_markers)),
                        "NOT SUPPORTED: Few differential genes")
  )

  message(sprintf("  DE genes: %d", nrow(outlier_markers)))
  message(sprintf("  Top upregulated: %s", paste(head(top_up, 5), collapse = ", ")))
  message(sprintf("  Top downregulated: %s", paste(head(top_down, 5), collapse = ", ")))
} else {
  results$H4 <- list(
    hypothesis = "Outlier cells show distinct gene expression",
    conclusion = "NOT TESTED: No outlier clusters identified"
  )
  outlier_markers <- NULL
}
message(sprintf("  Conclusion: %s", results$H4$conclusion))

# ==============================================================================
# PSEUDOTIME DISTRIBUTION COMPARISON
# ==============================================================================

message("\n========================================")
message("Pseudotime Distribution Analysis")
message("========================================\n")

# Compare pseudotime distributions by dataset and cell type
pt_comparison <- data.frame(
  cell = colnames(seurat_merged),
  dataset = seurat_merged$dataset,
  cell_type = seurat_merged$predicted_cell_type,
  pseudotime = seurat_merged$pseudotime_integrated
)

# Summary statistics
pt_summary <- pt_comparison %>%
  filter(!is.na(cell_type)) %>%
  group_by(dataset, cell_type) %>%
  summarise(
    n = n(),
    mean_pt = mean(pseudotime, na.rm = TRUE),
    sd_pt = sd(pseudotime, na.rm = TRUE),
    median_pt = median(pseudotime, na.rm = TRUE),
    .groups = "drop"
  )

message("Pseudotime by dataset and cell type:")
print(pt_summary)

# Statistical tests per cell type
message("\nStatistical tests (ACP vs Healthy):")
stat_tests <- list()
for (ct in epithelial_types) {
  healthy_pt <- pt_comparison$pseudotime[pt_comparison$dataset == "healthy" &
                                           pt_comparison$cell_type == ct]
  acp_pt <- pt_comparison$pseudotime[pt_comparison$dataset == "ACP" &
                                       pt_comparison$cell_type == ct]

  if (length(healthy_pt) >= 10 && length(acp_pt) >= 10) {
    ks <- ks.test(healthy_pt, acp_pt)
    wt <- wilcox.test(healthy_pt, acp_pt)

    stat_tests[[ct]] <- data.frame(
      cell_type = ct,
      n_healthy = length(healthy_pt),
      n_acp = length(acp_pt),
      mean_diff = mean(acp_pt, na.rm = TRUE) - mean(healthy_pt, na.rm = TRUE),
      ks_D = ks$statistic,
      ks_p = ks$p.value,
      wilcox_p = wt$p.value
    )
  }
}

if (length(stat_tests) > 0) {
  stat_tests_df <- do.call(rbind, stat_tests)
  stat_tests_df$ks_p_adj <- p.adjust(stat_tests_df$ks_p, method = "BH")
  stat_tests_df$wilcox_p_adj <- p.adjust(stat_tests_df$wilcox_p, method = "BH")
  print(stat_tests_df)
  results$pseudotime_tests <- stat_tests_df
}

# ==============================================================================
# GENERATE PLOTS
# ==============================================================================

message("\n========================================")
message("Generating Plots")
message("========================================\n")

# Plot 1: Integrated UMAP colored by dataset
message("Creating UMAP plots...")
p1 <- DimPlot(seurat_merged, group.by = "dataset", pt.size = 0.1) +
  scale_color_manual(values = c("healthy" = "#4DAF4A", "ACP" = "#E41A1C")) +
  labs(title = "Integrated UMAP by Dataset")

# Plot 2: UMAP colored by cell type
p2 <- DimPlot(seurat_merged, group.by = "predicted_cell_type", pt.size = 0.1) +
  labs(title = "Predicted Cell Type")

# Plot 3: UMAP colored by pseudotime
p3 <- FeaturePlot(seurat_merged, features = "pseudotime_integrated", pt.size = 0.1) +
  scale_color_viridis_c() +
  labs(title = "Integrated Pseudotime")

# Plot 4: UMAP highlighting outliers
p4 <- DimPlot(seurat_merged, group.by = "is_outlier_cluster", pt.size = 0.1) +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "#E41A1C")) +
  labs(title = "Outlier Clusters (ACP)")

# Combine UMAPs
umap_combined <- (p1 | p2) / (p3 | p4)

ggsave(file.path(fig_dir, "03_integrated_umap.pdf"), umap_combined,
       width = 14, height = 12)
ggsave(file.path(fig_dir, "03_integrated_umap.png"), umap_combined,
       width = 14, height = 12, dpi = 300)
message("  Saved UMAP plots")

# Plot 5: Pseudotime density by dataset and cell type
message("Creating pseudotime comparison plots...")
p5 <- ggplot(pt_comparison %>% filter(!is.na(cell_type)),
             aes(x = pseudotime, fill = dataset)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~cell_type, ncol = 1) +
  scale_fill_manual(values = c("healthy" = "#4DAF4A", "ACP" = "#E41A1C")) +
  theme_minimal() +
  labs(title = "Pseudotime Distribution by Cell Type",
       x = "Pseudotime", y = "Density")

# Plot 6: Pseudotime violin by cluster (ACP only)
acp_data <- pt_comparison %>%
  filter(dataset == "ACP") %>%
  left_join(data.frame(cell = colnames(seurat_merged),
                       cluster = seurat_merged$acp_cluster,
                       is_outlier = seurat_merged$is_outlier_cluster),
            by = "cell")

p6 <- ggplot(acp_data %>% filter(!is.na(cluster)),
             aes(x = reorder(cluster, pseudotime), y = pseudotime, fill = is_outlier)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "#E41A1C")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pseudotime by ACP Cluster",
       subtitle = "Red = outlier clusters",
       x = "Cluster", y = "Pseudotime")

# Plot 7: Cluster composition
p7 <- ggplot(acp_data %>% filter(!is.na(cluster)),
             aes(x = cluster, fill = cell_type)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Composition by ACP Cluster",
       x = "Cluster", y = "Proportion", fill = "Cell Type")

# Combine additional plots
additional_plots <- (p5 | (p6 / p7))
ggsave(file.path(fig_dir, "03_pseudotime_analysis.pdf"), additional_plots,
       width = 16, height = 10)
ggsave(file.path(fig_dir, "03_pseudotime_analysis.png"), additional_plots,
       width = 16, height = 10, dpi = 300)
message("  Saved pseudotime analysis plots")

# Plot 8: Outlier score visualization
p8 <- ggplot(outlier_stats, aes(x = reorder(cluster, outlier_score),
                                y = outlier_score, fill = is_outlier)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = outlier_threshold, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "#E41A1C")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Outlier Score by Cluster",
       subtitle = sprintf("Threshold = %.3f (80th percentile)", outlier_threshold),
       x = "Cluster", y = "Outlier Score")

ggsave(file.path(fig_dir, "03_outlier_scores.pdf"), p8, width = 10, height = 6)
message("  Saved outlier score plot")

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

message("\n========================================")
message("Saving Outputs")
message("========================================\n")

# Save integrated Seurat object
saveRDS(seurat_merged, file.path(objects_dir, "03_integrated_trajectory.rds"))
message(sprintf("Saved: %s", file.path(objects_dir, "03_integrated_trajectory.rds")))

# Save hypothesis test results
results_summary <- data.frame(
  hypothesis = c("H1", "H2", "H3", "H4"),
  description = c(
    "ACP maps to healthy cell types",
    "Trajectory order preserved",
    "Outlier clusters identified",
    "Outlier gene signatures"
  ),
  conclusion = c(
    results$H1$conclusion,
    results$H2$conclusion,
    results$H3$conclusion,
    results$H4$conclusion
  )
)
write.csv(results_summary, file.path(tables_dir, "03_hypothesis_results.csv"), row.names = FALSE)
message(sprintf("Saved: %s", file.path(tables_dir, "03_hypothesis_results.csv")))

# Save cluster statistics
write.csv(outlier_stats, file.path(tables_dir, "03_cluster_outlier_stats.csv"), row.names = FALSE)
message(sprintf("Saved: %s", file.path(tables_dir, "03_cluster_outlier_stats.csv")))

# Save pseudotime comparison
write.csv(pt_summary, file.path(tables_dir, "03_pseudotime_summary.csv"), row.names = FALSE)
message(sprintf("Saved: %s", file.path(tables_dir, "03_pseudotime_summary.csv")))

# Save outlier markers if available
if (!is.null(outlier_markers)) {
  outlier_markers$gene <- rownames(outlier_markers)
  write.csv(outlier_markers, file.path(tables_dir, "03_outlier_markers.csv"), row.names = FALSE)
  message(sprintf("Saved: %s", file.path(tables_dir, "03_outlier_markers.csv")))
}

# Save statistical tests
if (exists("stat_tests_df")) {
  write.csv(stat_tests_df, file.path(tables_dir, "03_pseudotime_stat_tests.csv"), row.names = FALSE)
  message(sprintf("Saved: %s", file.path(tables_dir, "03_pseudotime_stat_tests.csv")))
}

# Generate markdown report
report <- c(
  "# Integrated Trajectory Analysis Report",
  "",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Dataset Summary",
  "",
  sprintf("- Healthy reference cells: %d", sum(seurat_merged$dataset == "healthy")),
  sprintf("- ACP query cells: %d", sum(seurat_merged$dataset == "ACP")),
  sprintf("- Total integrated cells: %d", ncol(seurat_merged)),
  "",
  "## Hypothesis Test Results",
  "",
  "| Hypothesis | Description | Result |",
  "|------------|-------------|--------|",
  sprintf("| H1 | ACP maps to healthy cell types | %s |", results$H1$conclusion),
  sprintf("| H2 | Trajectory order preserved | %s |", results$H2$conclusion),
  sprintf("| H3 | Outlier clusters identified | %s |", results$H3$conclusion),
  sprintf("| H4 | Outlier gene signatures | %s |", results$H4$conclusion),
  "",
  "## Outlier Cluster Analysis",
  "",
  sprintf("- Total ACP clusters: %d", nrow(outlier_stats)),
  sprintf("- Outlier clusters: %d", sum(outlier_stats$is_outlier)),
  sprintf("- Outlier threshold (80th percentile): %.3f", outlier_threshold),
  "",
  if (sum(outlier_stats$is_outlier) > 0) {
    c("### Outlier Clusters:",
      paste("- Cluster", outlier_stats$cluster[outlier_stats$is_outlier],
            sprintf("(n=%d, score=%.3f)",
                    outlier_stats$n_cells[outlier_stats$is_outlier],
                    outlier_stats$outlier_score[outlier_stats$is_outlier])))
  } else {
    "No outlier clusters identified."
  },
  "",
  "## Interpretation",
  "",
  "The integration-first approach ensures that trajectory comparisons are made in a",
  "batch-corrected space, reducing technical artifacts. Outlier clusters represent",
  "ACP-specific cell states that deviate from the normal epithelial differentiation",
  "trajectory observed in healthy tissue."
)

writeLines(report, file.path(tables_dir, "03_integrated_trajectory_report.md"))
message(sprintf("Saved: %s", file.path(tables_dir, "03_integrated_trajectory_report.md")))

# ==============================================================================
# COMPLETION
# ==============================================================================

message("\n")
message("================================================================")
message("  Integrated Trajectory Analysis Complete!")
message("================================================================")
message(sprintf("  Finished: %s", Sys.time()))
message(sprintf("  Integrated cells: %d", ncol(seurat_merged)))
message(sprintf("  ACP clusters: %d", nrow(outlier_stats)))
message(sprintf("  Outlier clusters: %d", sum(outlier_stats$is_outlier)))
message("")
message("Hypothesis Results:")
message(sprintf("  H1 (Label transfer): %s", results$H1$conclusion))
message(sprintf("  H2 (Order preserved): %s", results$H2$conclusion))
message(sprintf("  H3 (Outliers found): %s", results$H3$conclusion))
message(sprintf("  H4 (Outlier signatures): %s", results$H4$conclusion))
message("")
message("Outputs:")
message(sprintf("  Objects: %s", objects_dir))
message(sprintf("  Tables: %s", tables_dir))
message(sprintf("  Figures: %s", fig_dir))
message("================================================================\n")
