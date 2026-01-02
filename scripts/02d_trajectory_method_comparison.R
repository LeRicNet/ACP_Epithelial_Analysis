#!/usr/bin/env Rscript
# ==============================================================================
# scripts/02d_trajectory_method_comparison.R
# ==============================================================================
# Diagnostic script to investigate disagreement between trajectory inference
# methods (Monocle3 vs Slingshot) in the CELLxGENE reference dataset.
#
# This script performs detailed comparison of:
#   1. Overall pseudotime correlation and scatter plots
#   2. Per-subtype method agreement
#   3. Per-cell_type method agreement (CELLxGENE annotations)
#   4. Lineage-specific analysis (Slingshot multi-lineage)
#   5. Spatial (UMAP) visualization of disagreement
#   6. Root cell and branch point analysis
#
# Usage:
#   Rscript scripts/02d_trajectory_method_comparison.R
#   Rscript scripts/02d_trajectory_method_comparison.R --dataset cellxgene
#
# Input:
#   - 02_seurat_with_pseudotime_{dataset}.rds
#   - 02_trajectory_results_{dataset}.rds
#
# Output:
#   - Detailed diagnostic plots
#   - Method comparison tables
#   - Recommendations for interpretation
#
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
dataset_type <- "cellxgene"  # Default to cellxgene for this diagnostic

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
message("  Trajectory Method Comparison Diagnostic")
message(sprintf("  Dataset: %s", toupper(dataset_type)))
message("================================================================")
message(paste("  Started:", Sys.time()))
message("\n")

# Load packages
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# Source utilities
source("R/utils/config.R")
config <- load_config(config_path)

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("========================================")
message("Loading Data")
message("========================================\n")

objects_dir <- get_path(config, config$paths$objects_dir)

# Load Seurat object with pseudotime
seurat_path <- file.path(objects_dir, sprintf("02_seurat_with_pseudotime_%s.rds", dataset_type))
if (!file.exists(seurat_path)) {
  stop(sprintf("Seurat object not found: %s\nRun 02_trajectory_analysis.R first.", seurat_path))
}

message(sprintf("Loading: %s", seurat_path))
seurat_obj <- readRDS(seurat_path)
message(sprintf("  Cells: %d", ncol(seurat_obj)))

# Load trajectory results
results_path <- file.path(objects_dir, sprintf("02_trajectory_results_%s.rds", dataset_type))
if (!file.exists(results_path)) {
  stop(sprintf("Trajectory results not found: %s", results_path))
}

message(sprintf("Loading: %s", results_path))
traj_results <- readRDS(results_path)

# Check what methods were used
methods_used <- traj_results$methods_used
message(sprintf("\nMethods used: %s", paste(methods_used, collapse = ", ")))

if (!all(c("monocle3", "slingshot") %in% methods_used)) {
  stop("This script requires both monocle3 and slingshot results")
}

# ==============================================================================
# EXTRACT PSEUDOTIME DATA
# ==============================================================================

message("\n========================================")
message("Extracting Pseudotime Data")
message("========================================\n")

# Get pseudotime columns
pt_monocle <- seurat_obj$pseudotime_monocle3
pt_slingshot <- seurat_obj$pseudotime_slingshot
pt_consensus <- seurat_obj$pseudotime_consensus

# Normalize to [0,1] for fair comparison
normalize_pt <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

pt_monocle_norm <- normalize_pt(pt_monocle)
pt_slingshot_norm <- normalize_pt(pt_slingshot)

# Create comparison dataframe
comparison_df <- data.frame(
  cell = colnames(seurat_obj),
  monocle3 = pt_monocle,
  slingshot = pt_slingshot,
  consensus = pt_consensus,
  monocle3_norm = pt_monocle_norm,
  slingshot_norm = pt_slingshot_norm,
  subtype = seurat_obj$module_score_subtype,
  stringsAsFactors = FALSE
)

# Add cell_type if available (CELLxGENE)
if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
  comparison_df$cell_type <- seurat_obj$cell_type
}

# Add UMAP coordinates
if ("umap" %in% tolower(names(seurat_obj@reductions))) {
  umap_coords <- Embeddings(seurat_obj, "umap")
  comparison_df$UMAP_1 <- umap_coords[, 1]
  comparison_df$UMAP_2 <- umap_coords[, 2]
}

message(sprintf("  Monocle3 pseudotime range: [%.3f, %.3f]",
                min(pt_monocle, na.rm = TRUE), max(pt_monocle, na.rm = TRUE)))
message(sprintf("  Slingshot pseudotime range: [%.3f, %.3f]",
                min(pt_slingshot, na.rm = TRUE), max(pt_slingshot, na.rm = TRUE)))

# ==============================================================================
# OVERALL CORRELATION ANALYSIS
# ==============================================================================

message("\n========================================")
message("Overall Correlation Analysis")
message("========================================\n")

# Calculate correlations
cor_pearson <- cor(pt_monocle_norm, pt_slingshot_norm, use = "pairwise.complete.obs")
cor_spearman <- cor(pt_monocle_norm, pt_slingshot_norm, use = "pairwise.complete.obs", method = "spearman")
cor_kendall <- cor(pt_monocle_norm, pt_slingshot_norm, use = "pairwise.complete.obs", method = "kendall")

message(sprintf("  Pearson correlation:  r = %.4f", cor_pearson))
message(sprintf("  Spearman correlation: ρ = %.4f", cor_spearman))
message(sprintf("  Kendall correlation:  τ = %.4f", cor_kendall))

# Check if methods are anti-correlated (reversed)
if (cor_spearman < -0.5) {
  message("\n  WARNING: Methods are anti-correlated!")
  message("  This suggests opposite trajectory directions.")
  message("  Consider flipping one method's pseudotime.")
}

# Calculate agreement metrics
comparison_df$pt_diff <- abs(comparison_df$monocle3_norm - comparison_df$slingshot_norm)
comparison_df$agree_direction <- sign(comparison_df$monocle3_norm - 0.5) ==
  sign(comparison_df$slingshot_norm - 0.5)

mean_diff <- mean(comparison_df$pt_diff, na.rm = TRUE)
pct_agree <- mean(comparison_df$agree_direction, na.rm = TRUE) * 100

message(sprintf("\n  Mean absolute difference: %.4f", mean_diff))
message(sprintf("  Direction agreement: %.1f%%", pct_agree))

# Identify highly discordant cells
comparison_df$discordant <- comparison_df$pt_diff > 0.4
n_discordant <- sum(comparison_df$discordant, na.rm = TRUE)
pct_discordant <- mean(comparison_df$discordant, na.rm = TRUE) * 100

message(sprintf("  Highly discordant cells (diff > 0.4): %d (%.1f%%)", n_discordant, pct_discordant))

# ==============================================================================
# PER-SUBTYPE ANALYSIS
# ==============================================================================

message("\n========================================")
message("Per-Subtype Correlation Analysis")
message("========================================\n")

subtype_correlations <- comparison_df %>%
  group_by(subtype) %>%
  summarise(
    n_cells = n(),
    pearson_r = cor(monocle3_norm, slingshot_norm, use = "pairwise.complete.obs"),
    spearman_rho = cor(monocle3_norm, slingshot_norm, use = "pairwise.complete.obs", method = "spearman"),
    mean_monocle = mean(monocle3_norm, na.rm = TRUE),
    mean_slingshot = mean(slingshot_norm, na.rm = TRUE),
    mean_diff = mean(pt_diff, na.rm = TRUE),
    pct_discordant = mean(discordant, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(mean_monocle)

message("Correlation by Module Score Subtype:")
print(subtype_correlations)

# Identify problematic subtypes
problematic_subtypes <- subtype_correlations %>%
  filter(abs(spearman_rho) < 0.3 | pct_discordant > 30)

if (nrow(problematic_subtypes) > 0) {
  message("\n  Subtypes with poor method agreement:")
  for (i in 1:nrow(problematic_subtypes)) {
    message(sprintf("    - %s: ρ = %.3f, %.1f%% discordant",
                    problematic_subtypes$subtype[i],
                    problematic_subtypes$spearman_rho[i],
                    problematic_subtypes$pct_discordant[i]))
  }
}

# ==============================================================================
# PER-CELL_TYPE ANALYSIS (CELLxGENE)
# ==============================================================================

celltype_correlations <- NULL
if ("cell_type" %in% colnames(comparison_df)) {
  message("\n========================================")
  message("Per-Cell Type Correlation Analysis (CELLxGENE)")
  message("========================================\n")

  celltype_correlations <- comparison_df %>%
    group_by(cell_type) %>%
    summarise(
      n_cells = n(),
      pearson_r = cor(monocle3_norm, slingshot_norm, use = "pairwise.complete.obs"),
      spearman_rho = cor(monocle3_norm, slingshot_norm, use = "pairwise.complete.obs", method = "spearman"),
      mean_monocle = mean(monocle3_norm, na.rm = TRUE),
      mean_slingshot = mean(slingshot_norm, na.rm = TRUE),
      mean_diff = mean(pt_diff, na.rm = TRUE),
      pct_discordant = mean(discordant, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    arrange(mean_monocle)

  message("Correlation by Original Cell Type:")
  print(celltype_correlations)
}

# ==============================================================================
# SLINGSHOT LINEAGE ANALYSIS
# ==============================================================================

message("\n========================================")
message("Slingshot Lineage Analysis")
message("========================================\n")

# Check if slingshot found multiple lineages
slingshot_results <- traj_results$pooled$slingshot

if (!is.null(slingshot_results)) {
  n_lineages <- slingshot_results$n_lineages
  message(sprintf("  Number of lineages detected: %d", n_lineages))

  if (n_lineages > 1) {
    message("\n  Multiple lineages detected - this may explain disagreement!")
    message("  Monocle3 creates a single graph, while Slingshot found branching.")

    # Get lineage info
    lineages <- slingshot_results$lineages
    message("\n  Lineage composition:")
    for (i in seq_along(lineages)) {
      message(sprintf("    Lineage %d: %s", i, paste(lineages[[i]], collapse = " → ")))
    }

    # If we have per-lineage pseudotime
    if (!is.null(slingshot_results$pseudotime_matrix) &&
        is.matrix(slingshot_results$pseudotime_matrix)) {
      pt_matrix <- slingshot_results$pseudotime_matrix

      message("\n  Per-lineage pseudotime coverage:")
      for (i in 1:ncol(pt_matrix)) {
        n_assigned <- sum(!is.na(pt_matrix[, i]))
        pct_assigned <- n_assigned / nrow(pt_matrix) * 100
        message(sprintf("    Lineage %d: %d cells (%.1f%%)", i, n_assigned, pct_assigned))
      }

      # Add lineage-specific pseudotime to comparison
      comparison_df$slingshot_lineage1 <- pt_matrix[, 1]
      if (ncol(pt_matrix) > 1) {
        comparison_df$slingshot_lineage2 <- pt_matrix[, 2]

        # Cells assigned to only one lineage
        lineage1_only <- !is.na(pt_matrix[, 1]) & is.na(pt_matrix[, 2])
        lineage2_only <- is.na(pt_matrix[, 1]) & !is.na(pt_matrix[, 2])
        both_lineages <- !is.na(pt_matrix[, 1]) & !is.na(pt_matrix[, 2])

        message("\n  Lineage assignment:")
        message(sprintf("    Lineage 1 only: %d cells (%.1f%%)",
                        sum(lineage1_only), mean(lineage1_only) * 100))
        message(sprintf("    Lineage 2 only: %d cells (%.1f%%)",
                        sum(lineage2_only), mean(lineage2_only) * 100))
        message(sprintf("    Both lineages: %d cells (%.1f%%)",
                        sum(both_lineages), mean(both_lineages) * 100))

        comparison_df$lineage_assignment <- case_when(
          lineage1_only ~ "Lineage 1 only",
          lineage2_only ~ "Lineage 2 only",
          both_lineages ~ "Both lineages",
          TRUE ~ "Unassigned"
        )
      }
    }
  }
}

# ==============================================================================
# MONOCLE3 GRAPH ANALYSIS
# ==============================================================================

message("\n========================================")
message("Monocle3 Graph Analysis")
message("========================================\n")

monocle_results <- traj_results$pooled$monocle3

if (!is.null(monocle_results)) {
  graph_metrics <- monocle_results$graph_metrics

  if (!is.null(graph_metrics)) {
    message(sprintf("  Graph nodes: %d", graph_metrics$n_nodes))
    message(sprintf("  Graph edges: %d", graph_metrics$n_edges))
    message(sprintf("  Branch points: %d", graph_metrics$n_branch_points))
    message(sprintf("  Tips (endpoints): %d", graph_metrics$n_tips))
    message(sprintf("  Graph diameter: %d", graph_metrics$graph_diameter))

    if (graph_metrics$n_branch_points > 0) {
      message("\n  Note: Monocle3 detected branch points but may have")
      message("  ordered cells along a single principal path.")
    }
  }
}

# ==============================================================================
# RANK-BASED ANALYSIS
# ==============================================================================

message("\n========================================")
message("Rank-Based Analysis")
message("========================================\n")

# Convert pseudotime to ranks within each method
comparison_df$rank_monocle <- rank(comparison_df$monocle3_norm, ties.method = "average")
comparison_df$rank_slingshot <- rank(comparison_df$slingshot_norm, ties.method = "average")
comparison_df$rank_diff <- abs(comparison_df$rank_monocle - comparison_df$rank_slingshot)

# Mean rank by subtype
rank_by_subtype <- comparison_df %>%
  group_by(subtype) %>%
  summarise(
    mean_rank_monocle = mean(rank_monocle),
    mean_rank_slingshot = mean(rank_slingshot),
    rank_order_monocle = rank(mean(monocle3_norm, na.rm = TRUE)),
    rank_order_slingshot = rank(mean(slingshot_norm, na.rm = TRUE)),
    .groups = "drop"
  )

message("Subtype ordering comparison:")
message("\n  By Monocle3 (early → late):")
monocle_order <- comparison_df %>%
  group_by(subtype) %>%
  summarise(mean_pt = mean(monocle3_norm, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_pt)
message(sprintf("    %s", paste(monocle_order$subtype, collapse = " → ")))

message("\n  By Slingshot (early → late):")
slingshot_order <- comparison_df %>%
  group_by(subtype) %>%
  summarise(mean_pt = mean(slingshot_norm, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_pt)
message(sprintf("    %s", paste(slingshot_order$subtype, collapse = " → ")))

# Check if orderings are similar
order_match <- all(monocle_order$subtype == slingshot_order$subtype)
message(sprintf("\n  Subtype orderings match: %s", ifelse(order_match, "YES", "NO")))

# ==============================================================================
# GENERATE DIAGNOSTIC PLOTS
# ==============================================================================

message("\n========================================")
message("Generating Diagnostic Plots")
message("========================================\n")

fig_dir <- get_path(config, config$paths$figures_supp_dir)
ensure_dir(fig_dir)

colors <- get_colors(config, "epithelial_subtypes")

# Plot 1: Scatter plot of pseudotime values
message("Creating pseudotime scatter plot...")

p1_scatter <- ggplot(comparison_df, aes(x = monocle3_norm, y = slingshot_norm, color = subtype)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid") +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(title = "Monocle3 vs Slingshot Pseudotime",
       subtitle = sprintf("Spearman ρ = %.3f, Pearson r = %.3f", cor_spearman, cor_pearson),
       x = "Monocle3 Pseudotime (normalized)",
       y = "Slingshot Pseudotime (normalized)",
       color = "Subtype") +
  coord_fixed()

# Plot 2: Scatter plot faceted by subtype
p2_scatter_facet <- ggplot(comparison_df, aes(x = monocle3_norm, y = slingshot_norm)) +
  geom_point(alpha = 0.3, size = 0.3, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~subtype, scales = "free") +
  theme_minimal() +
  labs(title = "Method Agreement by Subtype",
       x = "Monocle3 (normalized)",
       y = "Slingshot (normalized)")

# Plot 3: Distribution of pseudotime differences
p3_diff_dist <- ggplot(comparison_df, aes(x = pt_diff, fill = subtype)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "stack") +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "red") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Distribution of Pseudotime Differences",
       subtitle = "Red line = discordance threshold (0.4)",
       x = "Absolute Difference (Monocle3 - Slingshot)",
       y = "Count",
       fill = "Subtype")

# Plot 4: UMAP colored by discordance
if ("UMAP_1" %in% colnames(comparison_df)) {
  p4_umap_discord <- ggplot(comparison_df, aes(x = UMAP_1, y = UMAP_2, color = pt_diff)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_viridis_c(option = "plasma", name = "PT Difference") +
    theme_minimal() +
    labs(title = "UMAP: Method Disagreement",
         subtitle = "Yellow = high disagreement")

  # Highlight discordant cells
  p4b_umap_discordant <- ggplot(comparison_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(data = comparison_df[!comparison_df$discordant, ],
               color = "gray80", size = 0.3, alpha = 0.3) +
    geom_point(data = comparison_df[comparison_df$discordant, ],
               aes(color = subtype), size = 0.5, alpha = 0.7) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    labs(title = "UMAP: Highly Discordant Cells",
         subtitle = sprintf("n = %d cells with difference > 0.4", n_discordant),
         color = "Subtype")
}

# Plot 5: Side-by-side pseudotime on UMAP
if ("UMAP_1" %in% colnames(comparison_df)) {
  p5a_umap_monocle <- ggplot(comparison_df, aes(x = UMAP_1, y = UMAP_2, color = monocle3_norm)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_viridis_c(option = "viridis", name = "Pseudotime") +
    theme_minimal() +
    labs(title = "Monocle3 Pseudotime")

  p5b_umap_slingshot <- ggplot(comparison_df, aes(x = UMAP_1, y = UMAP_2, color = slingshot_norm)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_viridis_c(option = "viridis", name = "Pseudotime") +
    theme_minimal() +
    labs(title = "Slingshot Pseudotime")
}

# Plot 6: Violin plots comparing methods by subtype
comparison_long <- comparison_df %>%
  select(cell, subtype, monocle3_norm, slingshot_norm) %>%
  pivot_longer(cols = c(monocle3_norm, slingshot_norm),
               names_to = "method", values_to = "pseudotime") %>%
  mutate(method = gsub("_norm", "", method))

p6_violin <- ggplot(comparison_long, aes(x = subtype, y = pseudotime, fill = method)) +
  geom_violin(position = position_dodge(0.8), scale = "width", alpha = 0.7) +
  geom_boxplot(position = position_dodge(0.8), width = 0.1, fill = "white", outlier.size = 0.3) +
  scale_fill_manual(values = c("monocle3" = "#E41A1C", "slingshot" = "#377EB8")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pseudotime Distribution by Method and Subtype",
       x = "Subtype", y = "Pseudotime (normalized)", fill = "Method")

# Plot 7: Correlation heatmap by subtype
cor_matrix <- subtype_correlations %>%
  select(subtype, spearman_rho) %>%
  mutate(correlation = spearman_rho)

p7_cor_bar <- ggplot(cor_matrix, aes(x = reorder(subtype, spearman_rho),
                                     y = spearman_rho, fill = spearman_rho)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = c(-0.3, 0.3), linetype = "dashed", color = "red") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Method Correlation by Subtype",
       subtitle = "Dashed lines = weak correlation threshold (±0.3)",
       x = "Subtype", y = "Spearman Correlation", fill = "ρ")

# Plot 8: Cell type specific (if available)
p8_celltype <- NULL
if (!is.null(celltype_correlations)) {
  p8_celltype <- ggplot(celltype_correlations,
                        aes(x = reorder(cell_type, spearman_rho),
                            y = spearman_rho, fill = spearman_rho)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Method Correlation by Cell Type (CELLxGENE)",
         x = "Cell Type", y = "Spearman Correlation", fill = "ρ")
}

# Combine plots
message("Saving diagnostic plots...")

# Main comparison figure
main_comparison <- (p1_scatter | p6_violin) /
  (p3_diff_dist | p7_cor_bar) +
  plot_annotation(
    title = sprintf("Trajectory Method Comparison: %s", toupper(dataset_type)),
    subtitle = sprintf("Monocle3 vs Slingshot | Overall ρ = %.3f", cor_spearman),
    tag_levels = "A"
  )

ggsave(file.path(fig_dir, sprintf("02d_method_comparison_main_%s.pdf", dataset_type)),
       main_comparison, width = 16, height = 12)
ggsave(file.path(fig_dir, sprintf("02d_method_comparison_main_%s.png", dataset_type)),
       main_comparison, width = 16, height = 12, dpi = 300)
message("  Saved main comparison plot")

# UMAP comparison figure
if ("UMAP_1" %in% colnames(comparison_df)) {
  umap_comparison <- (p5a_umap_monocle | p5b_umap_slingshot) /
    (p4_umap_discord | p4b_umap_discordant) +
    plot_annotation(
      title = "Spatial Distribution of Method Disagreement",
      tag_levels = "A"
    )

  ggsave(file.path(fig_dir, sprintf("02d_method_comparison_umap_%s.pdf", dataset_type)),
         umap_comparison, width = 14, height = 12)
  ggsave(file.path(fig_dir, sprintf("02d_method_comparison_umap_%s.png", dataset_type)),
         umap_comparison, width = 14, height = 12, dpi = 300)
  message("  Saved UMAP comparison plot")
}

# Per-subtype scatter
ggsave(file.path(fig_dir, sprintf("02d_method_comparison_by_subtype_%s.pdf", dataset_type)),
       p2_scatter_facet, width = 14, height = 10)
ggsave(file.path(fig_dir, sprintf("02d_method_comparison_by_subtype_%s.png", dataset_type)),
       p2_scatter_facet, width = 14, height = 10, dpi = 300)
message("  Saved per-subtype scatter plot")

# Cell type specific
if (!is.null(p8_celltype)) {
  ggsave(file.path(fig_dir, sprintf("02d_method_comparison_celltype_%s.pdf", dataset_type)),
         p8_celltype, width = 10, height = 6)
  message("  Saved cell type correlation plot")
}

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

message("\n========================================")
message("Saving Results")
message("========================================\n")

tables_dir <- get_path(config, config$paths$tables_dir)
ensure_dir(tables_dir)

# Overall summary
overall_summary <- data.frame(
  metric = c("Pearson correlation", "Spearman correlation", "Kendall correlation",
             "Mean absolute difference", "Direction agreement (%)",
             "Discordant cells (%)", "N cells"),
  value = c(cor_pearson, cor_spearman, cor_kendall,
            mean_diff, pct_agree, pct_discordant, nrow(comparison_df))
)

write.csv(overall_summary,
          file.path(tables_dir, sprintf("02d_method_comparison_summary_%s.csv", dataset_type)),
          row.names = FALSE)
message(sprintf("Saved: 02d_method_comparison_summary_%s.csv", dataset_type))

# Per-subtype correlations
write.csv(subtype_correlations,
          file.path(tables_dir, sprintf("02d_correlation_by_subtype_%s.csv", dataset_type)),
          row.names = FALSE)
message(sprintf("Saved: 02d_correlation_by_subtype_%s.csv", dataset_type))

# Per-cell_type correlations (if available)
if (!is.null(celltype_correlations)) {
  write.csv(celltype_correlations,
            file.path(tables_dir, sprintf("02d_correlation_by_celltype_%s.csv", dataset_type)),
            row.names = FALSE)
  message(sprintf("Saved: 02d_correlation_by_celltype_%s.csv", dataset_type))
}

# Subtype ordering comparison
ordering_comparison <- data.frame(
  rank = 1:nrow(monocle_order),
  monocle3_order = monocle_order$subtype,
  slingshot_order = slingshot_order$subtype,
  match = monocle_order$subtype == slingshot_order$subtype
)

write.csv(ordering_comparison,
          file.path(tables_dir, sprintf("02d_subtype_ordering_%s.csv", dataset_type)),
          row.names = FALSE)
message(sprintf("Saved: 02d_subtype_ordering_%s.csv", dataset_type))

# Save full comparison data
saveRDS(comparison_df,
        file.path(objects_dir, sprintf("02d_method_comparison_data_%s.rds", dataset_type)))
message(sprintf("Saved: 02d_method_comparison_data_%s.rds", dataset_type))

# ==============================================================================
# GENERATE INTERPRETATION REPORT
# ==============================================================================

message("\n========================================")
message("Interpretation and Recommendations")
message("========================================\n")

report <- c(
  sprintf("# Trajectory Method Comparison Report: %s", toupper(dataset_type)),
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Overall Agreement",
  sprintf("- Spearman correlation: ρ = %.3f", cor_spearman),
  sprintf("- Mean absolute difference: %.3f", mean_diff),
  sprintf("- Highly discordant cells: %.1f%%", pct_discordant),
  ""
)

# Interpretation based on correlation
if (abs(cor_spearman) > 0.7) {
  report <- c(report,
              "**Interpretation: GOOD AGREEMENT**",
              "Both methods produce similar pseudotime orderings.",
              "The consensus pseudotime is reliable.",
              "")
} else if (abs(cor_spearman) > 0.4) {
  report <- c(report,
              "**Interpretation: MODERATE AGREEMENT**",
              "Methods show partial agreement. Consider:",
              "- Using method-specific pseudotime for different analyses",
              "- Investigating discordant cell populations",
              "- The consensus may smooth over biological complexity",
              "")
} else if (abs(cor_spearman) > 0.1) {
  report <- c(report,
              "**Interpretation: POOR AGREEMENT**",
              "Methods strongly disagree. This could indicate:",
              "- Multiple trajectories/branching that methods capture differently",
              "- Different starting points or trajectory directions",
              "- Noise or batch effects affecting one method more",
              "",
              "**Recommendations:**",
              "- Examine each method's results separately",
              "- Consider biological interpretation of each trajectory",
              "- Do not rely on consensus pseudotime",
              "")
} else {
  report <- c(report,
              "**Interpretation: NO AGREEMENT**",
              "Methods produce essentially uncorrelated orderings.",
              "Possible causes:",
              "- Branching trajectories (Slingshot: lineages, Monocle3: graph)",
              "- Circular or complex topology",
              "- One method may be more appropriate for this data",
              "",
              "**Recommendations:**",
              "- Analyze Slingshot lineages separately",
              "- Examine Monocle3 trajectory graph structure",
              "- Consider which method's assumptions better fit your biology",
              "")
}

# Add lineage-specific notes
if (!is.null(slingshot_results) && slingshot_results$n_lineages > 1) {
  report <- c(report,
              "## Slingshot Lineage Analysis",
              sprintf("- Detected %d lineages", slingshot_results$n_lineages),
              "- This branching may explain disagreement with Monocle3",
              "- Consider analyzing each lineage separately",
              "")
}

# Add subtype-specific notes
poor_subtypes <- subtype_correlations %>% filter(abs(spearman_rho) < 0.3)
if (nrow(poor_subtypes) > 0) {
  report <- c(report,
              "## Problematic Subtypes",
              sprintf("The following subtypes show poor method agreement (|ρ| < 0.3):"),
              paste(sprintf("- %s (ρ = %.3f)", poor_subtypes$subtype, poor_subtypes$spearman_rho),
                    collapse = "\n"),
              "")
}

report_text <- paste(report, collapse = "\n")
writeLines(report_text, file.path(tables_dir, sprintf("02d_method_comparison_report_%s.md", dataset_type)))
message(sprintf("Saved: 02d_method_comparison_report_%s.md", dataset_type))

# Print report to console
cat("\n")
cat(report_text)
cat("\n")

# ==============================================================================
# COMPLETION
# ==============================================================================

message("\n")
message("================================================================")
message("  Method Comparison Analysis Complete!")
message("================================================================")
message(sprintf("  Dataset: %s", toupper(dataset_type)))
message(sprintf("  Finished: %s", Sys.time()))
message(sprintf("  Overall Spearman ρ: %.3f", cor_spearman))
message(sprintf("  Discordant cells: %.1f%%", pct_discordant))
message("\nOutputs:")
message(sprintf("  Tables: %s", tables_dir))
message(sprintf("  Figures: %s", fig_dir))
message("================================================================\n")
