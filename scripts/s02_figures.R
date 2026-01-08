#!/usr/bin/env Rscript
# ==============================================================================
# Figure 0.2: ACP vs Healthy Epithelial Reference
# ==============================================================================
# Generates composite figure showing:
#   A: Label transfer quality (prediction score histogram)
#   B: Mean pseudotime by cell type and dataset
#   C: Pseudotime density distributions by cell type
#   D: UMAP of ACP cells with transferred labels
#   E: Pseudotime by ACP cluster (violin, outliers highlighted)
#   F: Cell type composition by ACP cluster
#
# Requires: 03_integrated_trajectory.rds from 03_integrated_trajectory_comparison.R
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  input_path <- args[1]
} else {
  input_path <- "results/objects/03_integrated_trajectory.rds"
}

# Output directories
output_dir <- dirname(dirname(input_path))
fig_main_dir <- file.path(output_dir, "figures/main")
fig_supp_dir <- file.path(output_dir, "figures/supplementary")

dir.create(fig_main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_supp_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
if (!file.exists(input_path)) {
  stop(sprintf("Input file not found: %s\nRun 03_integrated_trajectory_comparison.R first.", input_path))
}

message("Loading integrated trajectory data...")
seurat_obj <- readRDS(input_path)
message(sprintf("  Loaded %d cells", ncol(seurat_obj)))

# ==============================================================================
# THEME DEFINITION
# ==============================================================================

theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "gray30", hjust = 0),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      strip.text = element_text(face = "bold", size = base_size),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# Color palettes
dataset_colors <- c("healthy" = "#4DAF4A", "ACP" = "#E41A1C")
celltype_colors <- c(

  "basal cell of epidermis" = "#F8766D",
  "keratinocyte" = "#00BFC4",
  "spinous cell of epidermis" = "#7CAE00"
)
outlier_colors <- c("FALSE" = "#619CFF", "TRUE" = "#E41A1C")

# ==============================================================================
# PREPARE DATA
# ==============================================================================

message("Preparing data for plotting...")

# Extract metadata
meta <- seurat_obj@meta.data %>%
  mutate(
    cell_id = rownames(.),
    dataset = factor(dataset, levels = c("healthy", "ACP")),
    predicted_cell_type = factor(predicted_cell_type,
                                 levels = c("basal cell of epidermis",
                                            "keratinocyte",
                                            "spinous cell of epidermis"))
  )

# Get UMAP coordinates
umap_coords <- as.data.frame(Embeddings(seurat_obj, "umap"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
umap_coords$cell_id <- rownames(umap_coords)

# Merge
plot_data <- meta %>%
  left_join(umap_coords, by = "cell_id")

# ACP-only data
acp_data <- plot_data %>% filter(dataset == "ACP")

# Summary statistics
n_healthy <- sum(plot_data$dataset == "healthy")
n_acp <- sum(plot_data$dataset == "ACP")

message(sprintf("  Healthy cells: %d", n_healthy))
message(sprintf("  ACP cells: %d", n_acp))

# ==============================================================================
# PANEL A: Label Transfer Quality
# ==============================================================================

message("Creating Panel A: Label transfer quality...")

# Get prediction scores for ACP cells only
pred_scores <- acp_data$prediction_score[!is.na(acp_data$prediction_score)]
mean_score <- mean(pred_scores, na.rm = TRUE)
prop_above <- mean(pred_scores > 0.5, na.rm = TRUE) * 100

panel_a <- ggplot(data.frame(score = pred_scores), aes(x = score)) +
  geom_histogram(bins = 50, fill = "#619CFF", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#E41A1C", linewidth = 0.8) +
  annotate("text", x = 0.55, y = Inf, vjust = 2, hjust = 0,
           label = sprintf("Mean = %.3f\n%.1f%% > 0.5", mean_score, prop_above),
           size = 3, color = "black") +
  scale_x_continuous(limits = c(0.4, 1.0), breaks = seq(0.5, 1.0, 0.1)) +
  theme_pub() +
  labs(
    title = "A",
    subtitle = "Label transfer quality",
    x = "Prediction Score",
    y = "Count"
  )

# ==============================================================================
# PANEL B: Mean Pseudotime by Cell Type
# ==============================================================================

message("Creating Panel B: Mean pseudotime comparison...")

pt_summary <- plot_data %>%
  filter(!is.na(predicted_cell_type), !is.na(pseudotime_integrated)) %>%
  group_by(dataset, predicted_cell_type) %>%
  summarise(
    n = n(),
    mean_pt = mean(pseudotime_integrated, na.rm = TRUE),
    se_pt = sd(pseudotime_integrated, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

panel_b <- ggplot(pt_summary, aes(x = predicted_cell_type, y = mean_pt, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_pt - se_pt, ymax = mean_pt + se_pt),
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_fill_manual(values = dataset_colors, name = "Dataset") +
  scale_x_discrete(labels = function(x) gsub(" of epidermis", "\n(epidermis)", x)) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
    legend.position = c(0.85, 0.85)
  ) +
  labs(
    title = "B",
    subtitle = "Mean pseudotime by cell type",
    x = "Cell Type",
    y = "Mean Pseudotime"
  )

# ==============================================================================
# PANEL C: Pseudotime Density by Cell Type
# ==============================================================================

message("Creating Panel C: Pseudotime density distributions...")

density_data <- plot_data %>%
  filter(!is.na(predicted_cell_type), !is.na(pseudotime_integrated))

panel_c <- ggplot(density_data, aes(x = pseudotime_integrated, fill = dataset)) +
  geom_density(alpha = 0.6, color = "black", linewidth = 0.3) +
  facet_wrap(~predicted_cell_type, ncol = 1, scales = "free_y",
             labeller = labeller(predicted_cell_type = function(x) gsub(" of epidermis", "", x))) +
  scale_fill_manual(values = dataset_colors, name = "Dataset") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme_pub() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9),
    panel.spacing = unit(0.3, "lines")
  ) +
  labs(
    title = "C",
    subtitle = "Pseudotime distribution by cell type",
    x = "Pseudotime",
    y = "Density"
  )

# ==============================================================================
# PANEL D: UMAP with Transferred Labels (ACP only)
# ==============================================================================

message("Creating Panel D: UMAP of ACP cells...")

panel_d <- ggplot(acp_data %>% filter(!is.na(predicted_cell_type)),
                  aes(x = UMAP_1, y = UMAP_2, color = predicted_cell_type)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = celltype_colors,
                     name = "Transferred\nLabel",
                     labels = function(x) gsub(" of epidermis", "", x)) +
  theme_pub() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(
    title = "D",
    subtitle = "ACP cells by transferred label",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  coord_fixed()

# ==============================================================================
# PANEL E: Pseudotime by ACP Cluster
# ==============================================================================

message("Creating Panel E: Pseudotime by cluster...")

# Filter to ACP cells with cluster assignments
cluster_data <- acp_data %>%
  filter(!is.na(acp_cluster), !is.na(pseudotime_integrated)) %>%
  mutate(
    acp_cluster = factor(acp_cluster),
    is_outlier = as.character(is_outlier_cluster)
  )

# Order clusters by mean pseudotime
cluster_order <- cluster_data %>%
  group_by(acp_cluster) %>%
  summarise(mean_pt = mean(pseudotime_integrated, na.rm = TRUE)) %>%
  arrange(mean_pt) %>%
  pull(acp_cluster)

cluster_data$acp_cluster <- factor(cluster_data$acp_cluster, levels = cluster_order)

# Get outlier status per cluster for coloring
cluster_outlier <- cluster_data %>%
  group_by(acp_cluster) %>%
  summarise(is_outlier = first(is_outlier_cluster)) %>%
  mutate(is_outlier = as.character(is_outlier))

panel_e <- ggplot(cluster_data, aes(x = acp_cluster, y = pseudotime_integrated)) +
  geom_violin(aes(fill = is_outlier), scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = outlier_colors,
                    name = "Outlier",
                    labels = c("FALSE" = "No", "TRUE" = "Yes")) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = c(0.15, 0.85)
  ) +
  labs(
    title = "E",
    subtitle = "Pseudotime by ACP cluster",
    x = "Cluster",
    y = "Pseudotime"
  )

# ==============================================================================
# PANEL F: Cell Type Composition by Cluster
# ==============================================================================

message("Creating Panel F: Cell type composition...")

composition_data <- cluster_data %>%
  filter(!is.na(predicted_cell_type)) %>%
  group_by(acp_cluster, predicted_cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(acp_cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Use same cluster order as panel E
composition_data$acp_cluster <- factor(composition_data$acp_cluster, levels = cluster_order)

panel_f <- ggplot(composition_data, aes(x = acp_cluster, y = prop, fill = predicted_cell_type)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = celltype_colors,
                    name = "Cell Type",
                    labels = function(x) gsub(" of epidermis", "", x)) +
  scale_y_continuous(labels = percent_format()) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom"
  ) +
  labs(
    title = "F",
    subtitle = "Cell type composition by cluster",
    x = "Cluster",
    y = "Proportion"
  )

# ==============================================================================
# COMBINE PANELS
# ==============================================================================

message("Combining panels...")

# Layout:
# Row 1: A | B | D
# Row 2: C | E | F

# Top row
top_row <- panel_a + panel_b + panel_d +
  plot_layout(widths = c(1, 1.2, 1.3))

# Bottom row
bottom_row <- panel_c + panel_e + panel_f +
  plot_layout(widths = c(1, 1.2, 1.2))

# Combine
figure_02 <- top_row / bottom_row +
  plot_annotation(
    title = "ACP epithelial cells map to early differentiation states but show altered pseudotime distributions",
    subtitle = sprintf("Healthy reference: n = %s cells | ACP: n = %s cells | Integration: Harmony, 30 dimensions",
                       format(n_healthy, big.mark = ","),
                       format(n_acp, big.mark = ",")),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray30")
    )
  )

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

message("Saving figure...")

# Main figure
ggsave(
  file.path(fig_main_dir, "Fig02_acp_vs_healthy.pdf"),
  figure_02,
  width = 14, height = 10, units = "in"
)
ggsave(
  file.path(fig_main_dir, "Fig02_acp_vs_healthy.png"),
  figure_02,
  width = 14, height = 10, units = "in", dpi = 300
)

message(sprintf("  Saved: %s", file.path(fig_main_dir, "Fig02_acp_vs_healthy.pdf")))
message(sprintf("  Saved: %s", file.path(fig_main_dir, "Fig02_acp_vs_healthy.png")))

# ==============================================================================
# GENERATE SUMMARY STATISTICS TABLE
# ==============================================================================

message("Generating summary statistics...")

# Pseudotime comparison statistics
pt_stats <- plot_data %>%
  filter(!is.na(predicted_cell_type), !is.na(pseudotime_integrated)) %>%
  group_by(predicted_cell_type) %>%
  summarise(
    n_healthy = sum(dataset == "healthy"),
    n_acp = sum(dataset == "ACP"),
    mean_pt_healthy = mean(pseudotime_integrated[dataset == "healthy"], na.rm = TRUE),
    mean_pt_acp = mean(pseudotime_integrated[dataset == "ACP"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_diff = mean_pt_acp - mean_pt_healthy
  )

# Perform KS tests
ks_results <- plot_data %>%
  filter(!is.na(predicted_cell_type), !is.na(pseudotime_integrated)) %>%
  group_by(predicted_cell_type) %>%
  summarise(
    ks_stat = {
      healthy_pt <- pseudotime_integrated[dataset == "healthy"]
      acp_pt <- pseudotime_integrated[dataset == "ACP"]
      if (length(healthy_pt) > 5 && length(acp_pt) > 5) {
        ks.test(healthy_pt, acp_pt)$statistic
      } else {
        NA_real_
      }
    },
    ks_p = {
      healthy_pt <- pseudotime_integrated[dataset == "healthy"]
      acp_pt <- pseudotime_integrated[dataset == "ACP"]
      if (length(healthy_pt) > 5 && length(acp_pt) > 5) {
        ks.test(healthy_pt, acp_pt)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

pt_stats <- pt_stats %>%
  left_join(ks_results, by = "predicted_cell_type")

# Cluster outlier summary
outlier_summary <- cluster_data %>%
  group_by(acp_cluster, is_outlier_cluster) %>%
  summarise(
    n_cells = n(),
    mean_pt = mean(pseudotime_integrated, na.rm = TRUE),
    dominant_type = names(which.max(table(predicted_cell_type))),
    purity = max(table(predicted_cell_type)) / n(),
    .groups = "drop"
  )

# Save tables
tables_dir <- file.path(output_dir, "tables")
dir.create(tables_dir, showWarnings = FALSE)

write.csv(pt_stats, file.path(tables_dir, "Fig02_pseudotime_stats.csv"), row.names = FALSE)
write.csv(outlier_summary, file.path(tables_dir, "Fig02_cluster_summary.csv"), row.names = FALSE)

message(sprintf("  Saved: %s", file.path(tables_dir, "Fig02_pseudotime_stats.csv")))
message(sprintf("  Saved: %s", file.path(tables_dir, "Fig02_cluster_summary.csv")))

# ==============================================================================
# PRINT SUMMARY
# ==============================================================================

message("\n========================================")
message("Figure 0.2 Generation Complete")
message("========================================")
message(sprintf("  Healthy cells: %d", n_healthy))
message(sprintf("  ACP cells: %d", n_acp))
message(sprintf("  Mean prediction score: %.3f", mean_score))
message(sprintf("  Proportion > 0.5: %.1f%%", prop_above))
message(sprintf("  ACP clusters: %d", length(unique(cluster_data$acp_cluster))))
message(sprintf("  Outlier clusters: %d", sum(cluster_outlier$is_outlier == "TRUE")))
message("")
message("Pseudotime summary by cell type:")
print(as.data.frame(pt_stats))
message("")
message("Output files:")
message(sprintf("  %s", file.path(fig_main_dir, "Fig02_acp_vs_healthy.pdf")))
message(sprintf("  %s", file.path(tables_dir, "Fig02_pseudotime_stats.csv")))
message("========================================\n")
