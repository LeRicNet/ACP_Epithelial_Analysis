#!/usr/bin/env Rscript
# ==============================================================================
# Section 0.1 Figures: Label Validation
# ==============================================================================
# Generates:
#   - Main Figure 0: Streamlined 3-panel summary for main text
#   - Supplementary Figure S1: Complete 8-panel validation analysis
#
# Requires: label_validation_rigorous.rds output from label_validation_rigorous.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(gridExtra)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Paths - adjust these to match your project structure
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  results_path <- args[1]
} else {
  # Default path - modify as needed
  results_path <- "results/objects/label_validation_rigorous.rds"
}

output_dir <- dirname(dirname(results_path))
fig_main_dir <- file.path(output_dir, "figures/main")
fig_supp_dir <- file.path(output_dir, "figures/supplementary")

dir.create(fig_main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_supp_dir, recursive = TRUE, showWarnings = FALSE)

# Load results
if (!file.exists(results_path)) {
  stop(sprintf("Results file not found: %s\nRun label_validation_rigorous.R first.", results_path))
}

results <- readRDS(results_path)
message("Loaded validation results")

# ==============================================================================
# THEME AND COLOR DEFINITIONS
# ==============================================================================

# Publication-ready theme
theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, color = "gray30", hjust = 0),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      strip.text = element_text(face = "bold", size = base_size),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Color palette for clusters
cluster_colors <- c(
  "E_C1_PE" = "#4DAF4A",
  "E_C2_WC" = "#E41A1C",
  "E_C3_KC" = "#FF7F00",
  "E_C4_RHCG" = "#377EB8",
  "E_C5_Prolif" = "#984EA3"
)

# Assessment colors
assessment_colors <- c(
  "SEPARABLE" = "#4DAF4A",
  "PARTIALLY SEPARABLE" = "#FDB863",
  "OVERLAPPING" = "#E41A1C"
)

# ==============================================================================
# MAIN FIGURE (3 panels)
# ==============================================================================
# Panel A: Confusion matrix (normalized by actual class)
# Panel B: Permutation test showing observed vs null
# Panel C: Per-cluster assessment summary
# ==============================================================================

message("Generating Main Figure...")

# --- Panel A: Confusion Matrix ---
conf_norm <- results$confusion$normalized_confusion * 100
conf_df <- as.data.frame(as.table(conf_norm))
colnames(conf_df) <- c("Predicted", "Actual", "Percentage")

# Reorder factors for better visualization
cluster_order <- c("E_C1_PE", "E_C2_WC", "E_C3_KC", "E_C4_RHCG", "E_C5_Prolif")
conf_df$Predicted <- factor(conf_df$Predicted, levels = rev(cluster_order))
conf_df$Actual <- factor(conf_df$Actual, levels = cluster_order)

# Add labels for display (using original nomenclature)
label_map <- c(
  "E_C1_PE" = "Palisading\n(E_C1_PE)",
  "E_C2_WC" = "Whorl\n(E_C2_WC)",
  "E_C3_KC" = "Keratinocyte\n(E_C3_KC)",
  "E_C4_RHCG" = "RHCG+\n(E_C4_RHCG)",
  "E_C5_Prolif" = "Proliferating\n(E_C5_Prolif)"
)

panel_a <- ggplot(conf_df, aes(x = Actual, y = Predicted, fill = Percentage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f", Percentage)),
            size = 3, color = ifelse(conf_df$Percentage > 50, "white", "black")) +
  scale_fill_gradientn(
    colors = c("white", "#FEE08B", "#FC8D59", "#D73027"),
    values = scales::rescale(c(0, 20, 50, 100)),
    limits = c(0, 100),
    name = "% of\nActual"
  ) +
  scale_x_discrete(labels = label_map) +
  scale_y_discrete(labels = label_map) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm")
  ) +
  labs(
    title = "A",
    subtitle = "KNN classification confusion matrix",
    x = "Actual Label",
    y = "Predicted Label"
  ) +
  coord_fixed()

# --- Panel B: Permutation Test ---
perm_null <- results$permutation$null_balanced_distribution * 100
observed <- results$primary_cv$balanced_accuracy * 100
null_mean <- results$permutation$null_balanced_mean * 100
effect_size <- results$permutation$effect_size_balanced

perm_df <- data.frame(accuracy = perm_null)

panel_b <- ggplot(perm_df, aes(x = accuracy)) +
  geom_histogram(bins = 12, fill = "gray70", color = "white", alpha = 0.8) +
  geom_vline(xintercept = observed, color = "#D73027", linewidth = 1.2) +
  geom_vline(xintercept = null_mean, color = "gray40", linetype = "dashed", linewidth = 0.8) +
  annotate("segment", x = null_mean, xend = observed, y = 42, yend = 42,
           arrow = arrow(length = unit(0.2, "cm"), ends = "both"),
           color = "black", linewidth = 0.5) +
  annotate("text", x = (null_mean + observed) / 2, y = 44,
           label = sprintf("d = %.1f", effect_size),
           size = 3.5, fontface = "bold") +
  annotate("text", x = observed + 1.5, y = 35,
           label = sprintf("Observed\n%.1f%%", observed),
           size = 3, color = "#D73027", hjust = 0, fontface = "bold") +
  annotate("text", x = null_mean - 1, y = 25,
           label = sprintf("Null\n%.1f%%", null_mean),
           size = 3, color = "gray40", hjust = 1) +
  scale_x_continuous(limits = c(15, 60), breaks = seq(20, 60, 10)) +
  theme_pub() +
  labs(
    title = "B",
    subtitle = "Permutation test (50 × 5-fold CV)",
    x = "Balanced Accuracy (%)",
    y = "Count"
  )

# --- Panel C: Per-Cluster Assessment ---
cluster_df <- results$cluster_assessment %>%
  mutate(
    Cluster_Label = factor(label_map[Cluster], levels = rev(label_map)),
    Assessment = factor(Assessment, levels = c("SEPARABLE", "PARTIALLY SEPARABLE", "OVERLAPPING"))
  )

# Create a summary visualization
panel_c <- ggplot(cluster_df, aes(x = Cluster_Label, y = CV_Recall, fill = Assessment)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 60, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", CV_Recall)),
            hjust = -0.2, size = 3) +
  scale_fill_manual(values = assessment_colors, name = "Assessment") +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
  coord_flip() +
  theme_pub() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  labs(
    title = "C",
    subtitle = "Classification recall by cluster",
    x = "",
    y = "Recall (%)"
  )

# --- Combine Main Figure ---
main_figure <- panel_a + panel_b + panel_c +
  plot_layout(ncol = 3, widths = c(1.2, 1, 1)) +
  plot_annotation(
    title = "Validation of snRNA-seq epithelial subtype labels",
    subtitle = "n = 20,477 epithelial nuclei from 3 ACP patients (GSE215932)",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray30")
    )
  )

# Save main figure
ggsave(
  file.path(fig_main_dir, "Fig0_label_validation.pdf"),
  main_figure,
  width = 12, height = 4.5, units = "in"
)
ggsave(
  file.path(fig_main_dir, "Fig0_label_validation.png"),
  main_figure,
  width = 12, height = 4.5, units = "in", dpi = 300
)

message(sprintf("  Saved: %s", file.path(fig_main_dir, "Fig0_label_validation.pdf")))

# ==============================================================================
# SUPPLEMENTARY FIGURE (8 panels)
# ==============================================================================
# Panel A: Per-cluster metrics heatmap
# Panel B: Confusion matrix (same as main, for completeness)
# Panel C: Sensitivity analysis
# Panel D: Permutation test
# Panel E: Silhouette scores with CIs
# Panel F: Marker gene comparison (full vs power-matched)
# Panel G: Cross-dataset transfer performance
# Panel H: Per-class transfer comparison
# ==============================================================================

message("Generating Supplementary Figure...")

# --- Panel A: Per-Cluster Metrics Heatmap ---
metrics_df <- results$cluster_assessment %>%
  select(Cluster, CV_Recall, Silhouette, Markers_logFC1_Matched) %>%
  rename(Markers = Markers_logFC1_Matched) %>%
  pivot_longer(cols = c(CV_Recall, Silhouette, Markers),
               names_to = "Metric", values_to = "Value") %>%
  group_by(Metric) %>%
  mutate(
    Value_Scaled = (Value - min(Value, na.rm = TRUE)) /
      (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE) + 0.001)
  ) %>%
  ungroup() %>%
  mutate(
    Cluster = factor(Cluster, levels = cluster_order),
    Metric = factor(Metric, levels = c("CV_Recall", "Silhouette", "Markers"),
                    labels = c("CV Recall (%)", "Silhouette", "Markers (logFC>1)"))
  )

supp_a <- ggplot(metrics_df, aes(x = Cluster, y = Metric, fill = Value_Scaled)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Value, 1)), size = 2.5) +
  scale_fill_gradient2(
    low = "#E41A1C", mid = "#FFEDA0", high = "#4DAF4A",
    midpoint = 0.5, name = "Scaled\nValue"
  ) +
  scale_x_discrete(labels = function(x) gsub("E_C[0-9]_", "", x)) +
  theme_pub(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "A. Per-cluster separability metrics", x = "", y = "")

# --- Panel B: Confusion Matrix (detailed) ---
supp_b <- panel_a +
  labs(title = "B. Confusion matrix (normalized by actual class)")

# --- Panel C: Sensitivity Analysis ---
sensitivity <- results$sensitivity %>%
  mutate(k = factor(k), dims = factor(dims))

supp_c <- ggplot(sensitivity, aes(x = dims, y = balanced_accuracy * 100,
                                  color = k, group = k)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set1", name = "k") +
  theme_pub(base_size = 9) +
  labs(
    title = sprintf("C. Sensitivity analysis (range: %.1f%%)",
                    diff(range(sensitivity$balanced_accuracy)) * 100),
    x = "PCA Dimensions",
    y = "Balanced Accuracy (%)"
  )

# --- Panel D: Permutation Test ---
supp_d <- panel_b +
  labs(title = "D. Permutation null distribution")

# --- Panel E: Silhouette Scores ---
sil_df <- data.frame(
  Cluster = names(results$silhouette$per_cluster_silhouette),
  Silhouette = results$silhouette$per_cluster_silhouette,
  CI_Lower = results$silhouette$per_cluster_ci_lower[names(results$silhouette$per_cluster_silhouette)],
  CI_Upper = results$silhouette$per_cluster_ci_upper[names(results$silhouette$per_cluster_silhouette)]
) %>%
  mutate(Cluster = factor(Cluster, levels = cluster_order[order(results$silhouette$per_cluster_silhouette)]))

supp_e <- ggplot(sil_df, aes(x = reorder(Cluster, Silhouette), y = Silhouette)) +
  geom_bar(stat = "identity",
           fill = ifelse(sil_df$Silhouette > 0, "#4DAF4A", "#E41A1C"),
           width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.5, y = 0.27, label = "Moderate",
           hjust = 0, size = 2.5, color = "gray50") +
  scale_x_discrete(labels = function(x) gsub("E_C[0-9]_", "", x)) +
  coord_flip() +
  theme_pub(base_size = 9) +
  labs(
    title = sprintf("E. Silhouette scores (mean: %.3f)", results$silhouette$mean_silhouette),
    subtitle = "Error bars = 95% bootstrap CI",
    x = "",
    y = "Silhouette Score"
  )

# --- Panel F: Marker Comparison ---
if (!is.null(results$markers_full$marker_summary) &&
    !is.null(results$markers_matched$marker_summary)) {

  marker_comp <- results$markers_full$marker_summary %>%
    select(cluster, n_logfc_1.0) %>%
    rename(Full = n_logfc_1.0) %>%
    left_join(
      results$markers_matched$marker_summary %>%
        select(cluster, n_logfc_1.0) %>%
        rename(Matched = n_logfc_1.0),
      by = "cluster"
    ) %>%
    pivot_longer(cols = c(Full, Matched), names_to = "Analysis", values_to = "N_Markers") %>%
    mutate(cluster = factor(cluster, levels = cluster_order))

  supp_f <- ggplot(marker_comp, aes(x = cluster, y = N_Markers, fill = Analysis)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7,
             color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Full" = "#377EB8", "Matched" = "#FF7F00")) +
    scale_x_discrete(labels = function(x) gsub("E_C[0-9]_", "", x)) +
    theme_pub(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = c(0.85, 0.85)) +
    labs(
      title = "F. Marker genes (logFC > 1.0)",
      subtitle = "Power-matched subsampled to smallest cluster",
      x = "",
      y = "Number of Markers"
    )
} else {
  supp_f <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Marker data not available")
}

# --- Panels G & H: Cross-Dataset Transfer ---
if (!is.null(results$cross_dataset) && is.null(results$cross_dataset$error)) {

  # Panel G: Overall comparison
  comparison_df <- results$cross_dataset$comparison %>%
    mutate(Context = factor(Context, levels = c("GSE (native PCA)", "GSE (Harmony)",
                                                "ACP (transferred)", "Merged (Harmony)")))

  supp_g <- ggplot(comparison_df, aes(x = Context, y = Balanced_Accuracy, fill = Context)) +
    geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", Balanced_Accuracy)), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c(
      "GSE (native PCA)" = "#4DAF4A",
      "GSE (Harmony)" = "#377EB8",
      "ACP (transferred)" = "#FF7F00",
      "Merged (Harmony)" = "#984EA3"
    )) +
    scale_y_continuous(limits = c(0, max(comparison_df$Balanced_Accuracy) * 1.15)) +
    theme_pub(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = sprintf("G. Cross-dataset transfer (%s)", results$cross_dataset$transfer_quality),
      subtitle = sprintf("%.1f%% relative change in balanced accuracy",
                         results$cross_dataset$transfer_drop),
      x = "",
      y = "Balanced Accuracy (%)"
    )

  # Panel H: Per-class comparison
  transfer_df <- results$cross_dataset$per_class_transfer %>%
    pivot_longer(cols = c(GSE_Recall, ACP_Recall),
                 names_to = "Dataset", values_to = "Recall") %>%
    mutate(
      Dataset = factor(Dataset, levels = c("GSE_Recall", "ACP_Recall"),
                       labels = c("GSE (Harmony)", "ACP (transferred)")),
      Cluster = factor(Cluster, levels = cluster_order)
    )

  supp_h <- ggplot(transfer_df, aes(x = Cluster, y = Recall, fill = Dataset)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7,
             color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("GSE (Harmony)" = "#377EB8",
                                 "ACP (transferred)" = "#FF7F00")) +
    scale_x_discrete(labels = function(x) gsub("E_C[0-9]_", "", x)) +
    theme_pub(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = c(0.85, 0.85)
    ) +
    labs(
      title = "H. Per-class transfer performance",
      x = "",
      y = "Recall (%)"
    )

} else {
  supp_g <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Cross-dataset analysis not available")
  supp_h <- ggplot() + theme_void()
}

# --- Combine Supplementary Figure ---
supp_figure <- (supp_a | supp_b) / (supp_c | supp_d) / (supp_e | supp_f) / (supp_g | supp_h) +

  plot_annotation(
    title = "Supplementary Figure S1: Complete validation of snRNA-seq epithelial subtype labels",
    subtitle = sprintf("GSE215932: n = %d epithelial nuclei, %d clusters; ACP_SCN: n = 2,686 epithelial cells",
                       results$primary_cv$n_cells, results$primary_cv$n_classes),
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray30")
    )
  )

# Save supplementary figure
ggsave(
  file.path(fig_supp_dir, "FigS1_label_validation_complete.pdf"),
  supp_figure,
  width = 12, height = 16, units = "in"
)
ggsave(
  file.path(fig_supp_dir, "FigS1_label_validation_complete.png"),
  supp_figure,
  width = 12, height = 16, units = "in", dpi = 300
)

message(sprintf("  Saved: %s", file.path(fig_supp_dir, "FigS1_label_validation_complete.pdf")))

# ==============================================================================
# SUMMARY TABLE FOR MANUSCRIPT
# ==============================================================================

message("Generating summary table...")

summary_table <- results$cluster_assessment %>%
  select(Cluster, N_Cells, CV_Recall, Silhouette, Markers_logFC1_Matched, Assessment) %>%
  rename(
    `N Cells` = N_Cells,
    `CV Recall (%)` = CV_Recall,
    `Silhouette` = Silhouette,
    `Markers (logFC>1)` = Markers_logFC1_Matched,
    `Assessment` = Assessment
  )

write.csv(summary_table,
          file.path(output_dir, "tables/Table_S1_label_validation_summary.csv"),
          row.names = FALSE)

message(sprintf("  Saved: %s", file.path(output_dir, "tables/Table_S1_label_validation_summary.csv")))

# ==============================================================================
# COMPLETION
# ==============================================================================

message("\n========================================")
message("Figure generation complete")
message("========================================")
message(sprintf("Main figure: %s", file.path(fig_main_dir, "Fig0_label_validation.pdf")))
message(sprintf("Supplementary: %s", file.path(fig_supp_dir, "FigS1_label_validation_complete.pdf")))
message(sprintf("Summary table: %s", file.path(output_dir, "tables/Table_S1_label_validation_summary.csv")))
message("========================================\n")
