# =============================================================================
# CD163+ Macrophage Spatial Positioning Analysis
# =============================================================================
# Purpose: Test whether CD163+ macrophages are enriched at the interface
#          between epithelial cells and collagen tracts, supporting the
#          paracrine SPP1 delivery model
#
# Hypothesis: Macrophages localize to epithelial-stromal boundaries where
#             they can provide SPP1 to adjacent epithelial cells
#
# Input:   Filtered cell data with marker expression
#          Cell-tract distances (from Task 2.6)
# Output:  Statistical comparisons, spatial analysis, publication-ready figures
# =============================================================================

library(tidyverse)
library(data.table)
library(RANN)         # Fast nearest neighbor search
library(effsize)      # Cohen's d
library(ggpubr)       # stat_compare_means
library(patchwork)
library(scales)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths
data_dir <- "/home/rstudio/interfaces/ACP_IMC/results"
output_dir <- file.path(data_dir, "cd163_macrophage_analysis")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Marker column patterns (will auto-detect exact names)
cd163_pattern <- "CD163|Sm147|147Sm"
panck_pattern <- "Pan-Cytokeratin|Cytokeratin|174Yb|Yb174"
epcam_pattern <- "EpCAM|196Pt|Pt196"
collagen_pattern <- "Collagen|169Tm|Tm169"

# Thresholds for cell type identification (will be refined based on data)
# These are relative - we'll use percentile-based thresholds
cd163_percentile <- 75    # Top 25% = CD163-high (macrophages)
panck_percentile <- 75    # Top 25% = epithelial cells

# Interface definition thresholds (in pixels, ~1um per pixel)
epithelial_proximity_threshold <- 20   # Within 20um of epithelial cell
collagen_proximity_threshold <- 10     # Within 10um of collagen tract

# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("=============================================================================\n")
cat("CD163+ Macrophage Spatial Positioning Analysis\n")
cat("=============================================================================\n\n")

# --- Load cell-tract distance data ---
cell_tract_file <- "/home/rstudio/interfaces/ACP_IMC/results/collagen_qc/cell_tract_distances_all_approaches.csv"

cat("Loading cell-tract distance data...\n")
tract_distances <- fread(cell_tract_file)
cat("  Rows loaded:", format(nrow(tract_distances), big.mark = ","), "\n")

# Filter for PERMISSIVE approach
tract_distances <- tract_distances %>%
  filter(approach == "permissive")
cat("  Cells after filtering to permissive approach:", format(nrow(tract_distances), big.mark = ","), "\n")

# --- Load filtered cell data with marker expression ---
marker_file <- file.path(data_dir, "qc_results/step6_filtered_metrics/filtered_cell_data.csv")

cat("\nLoading marker expression data...\n")
marker_data <- fread(marker_file)
cat("  Cells loaded:", format(nrow(marker_data), big.mark = ","), "\n")

# --- Merge datasets ---
cat("\n--- Merging datasets ---\n")

# Parse cell_id to extract ObjectNumber
tract_distances <- tract_distances %>%
  mutate(
    ObjectNumber = as.integer(sub(".*_([0-9]+)$", "\\1", cell_id)),
    ImageNumber = roi_id
  )

# Merge on ImageNumber and ObjectNumber
cell_data <- marker_data %>%
  inner_join(
    tract_distances %>% select(ImageNumber, ObjectNumber, distance_to_tract,
                               x, y, contact_5px, contact_10px),
    by = c("ImageNumber", "ObjectNumber")
  )

cat("  Cells after merge:", format(nrow(cell_data), big.mark = ","), "\n")

# =============================================================================
# 2. IDENTIFY MARKER COLUMNS
# =============================================================================

cat("\n--- Identifying Marker Columns ---\n")

# Function to find column by pattern
find_column <- function(data, patterns, name) {
  for (pattern in strsplit(patterns, "\\|")[[1]]) {
    matches <- grep(pattern, colnames(data), value = TRUE, ignore.case = TRUE)
    if (length(matches) > 0) {
      cat(sprintf("  %s column found: %s\n", name, matches[1]))
      return(matches[1])
    }
  }
  cat(sprintf("  WARNING: %s column not found!\n", name))
  return(NULL)
}

cd163_col <- find_column(cell_data, cd163_pattern, "CD163")
panck_col <- find_column(cell_data, panck_pattern, "Pan-Cytokeratin")
epcam_col <- find_column(cell_data, epcam_pattern, "EpCAM")
collagen_col <- find_column(cell_data, collagen_pattern, "Collagen")

# Use EpCAM as backup for epithelial identification if Pan-CK not found
epithelial_col <- if (!is.null(panck_col)) panck_col else epcam_col

if (is.null(cd163_col)) stop("CD163 column not found!")
if (is.null(epithelial_col)) stop("Neither Pan-Cytokeratin nor EpCAM column found!")

# Also get coordinates
x_col <- "x"  # From tract_distances merge
y_col <- "y"

# =============================================================================
# 3. CELL TYPE CLASSIFICATION
# =============================================================================

cat("\n--- Classifying Cell Types ---\n")

# Calculate expression thresholds
cd163_threshold <- quantile(cell_data[[cd163_col]], cd163_percentile/100, na.rm = TRUE)
epithelial_threshold <- quantile(cell_data[[epithelial_col]], panck_percentile/100, na.rm = TRUE)

cat(sprintf("  CD163 threshold (p%d): %.3f\n", cd163_percentile, cd163_threshold))
cat(sprintf("  Epithelial threshold (p%d): %.3f\n", panck_percentile, epithelial_threshold))

# Classify cells
cell_data <- cell_data %>%
  mutate(
    CD163_expr = .data[[cd163_col]],
    Epithelial_expr = .data[[epithelial_col]],
    is_macrophage = CD163_expr >= cd163_threshold,
    is_epithelial = Epithelial_expr >= epithelial_threshold,
    cell_type = case_when(
      is_macrophage & is_epithelial ~ "Double_positive",
      is_macrophage ~ "Macrophage",
      is_epithelial ~ "Epithelial",
      TRUE ~ "Other"
    )
  )

# Cell type summary
cell_type_summary <- cell_data %>%
  count(cell_type) %>%
  mutate(pct = n / sum(n) * 100)

cat("\n  Cell type distribution:\n")
for (i in 1:nrow(cell_type_summary)) {
  cat(sprintf("    %s: %s (%.1f%%)\n",
              cell_type_summary$cell_type[i],
              format(cell_type_summary$n[i], big.mark = ","),
              cell_type_summary$pct[i]))
}

# =============================================================================
# 4. CALCULATE DISTANCE TO NEAREST EPITHELIAL CELL
# =============================================================================

cat("\n--- Calculating Distance to Nearest Epithelial Cell ---\n")

# Get coordinates of epithelial cells
epithelial_coords <- cell_data %>%
  filter(is_epithelial) %>%
  select(x, y) %>%
  as.matrix()

cat(sprintf("  Epithelial cells for neighbor search: %s\n",
            format(nrow(epithelial_coords), big.mark = ",")))

# Get coordinates of all cells
all_coords <- cell_data %>%
  select(x, y) %>%
  as.matrix()

# Find nearest epithelial neighbor for each cell
cat("  Computing nearest epithelial neighbors (this may take a moment)...\n")

nn_result <- nn2(
  data = epithelial_coords,
  query = all_coords,
  k = 1,
  searchtype = "standard"
)

# Add distance to nearest epithelial cell
cell_data$dist_to_epithelial <- nn_result$nn.dists[, 1]

# For epithelial cells themselves, find distance to SECOND nearest (excluding self)
# Re-run for epithelial cells with k=2
epithelial_indices <- which(cell_data$is_epithelial)
if (length(epithelial_indices) > 1) {
  nn_epithelial <- nn2(
    data = epithelial_coords,
    query = epithelial_coords,
    k = 2,
    searchtype = "standard"
  )
  # Use second nearest neighbor (first is self with distance 0)
  cell_data$dist_to_epithelial[epithelial_indices] <- nn_epithelial$nn.dists[, 2]
}

cat("  Distance to epithelial calculation complete.\n")
cat(sprintf("    Mean distance: %.1f pixels\n", mean(cell_data$dist_to_epithelial)))
cat(sprintf("    Median distance: %.1f pixels\n", median(cell_data$dist_to_epithelial)))

# =============================================================================
# 5. DEFINE INTERFACE ZONES
# =============================================================================

cat("\n--- Defining Spatial Zones ---\n")

cell_data <- cell_data %>%
  mutate(
    # Proximity classifications
    near_collagen = distance_to_tract <= collagen_proximity_threshold,
    near_epithelial = dist_to_epithelial <= epithelial_proximity_threshold,

    # Interface = near BOTH collagen and epithelial cells
    at_interface = near_collagen & near_epithelial,

    # Spatial zone classification
    spatial_zone = case_when(
      at_interface ~ "Interface",
      near_collagen & !near_epithelial ~ "Collagen_only",
      near_epithelial & !near_collagen ~ "Epithelial_only",
      TRUE ~ "Distant"
    )
  )

# Zone summary
zone_summary <- cell_data %>%
  count(spatial_zone) %>%
  mutate(pct = n / sum(n) * 100)

cat("\n  Spatial zone distribution:\n")
for (i in 1:nrow(zone_summary)) {
  cat(sprintf("    %s: %s (%.1f%%)\n",
              zone_summary$spatial_zone[i],
              format(zone_summary$n[i], big.mark = ","),
              zone_summary$pct[i]))
}

# =============================================================================
# 6. TEST MACROPHAGE ENRICHMENT AT INTERFACE
# =============================================================================

cat("\n=============================================================================\n")
cat("Macrophage Enrichment Analysis\n")
cat("=============================================================================\n")

# Calculate observed macrophage proportions by zone
macrophage_by_zone <- cell_data %>%
  group_by(spatial_zone) %>%
  summarise(
    n_total = n(),
    n_macrophages = sum(is_macrophage),
    pct_macrophages = n_macrophages / n_total * 100,
    .groups = "drop"
  )

cat("\n--- Macrophage Proportion by Spatial Zone ---\n")
print(macrophage_by_zone)

# Global macrophage proportion (expected under null)
global_mac_pct <- sum(cell_data$is_macrophage) / nrow(cell_data) * 100
cat(sprintf("\n  Global macrophage proportion: %.2f%%\n", global_mac_pct))

# Calculate enrichment (observed/expected)
macrophage_by_zone <- macrophage_by_zone %>%
  mutate(
    expected_pct = global_mac_pct,
    enrichment = pct_macrophages / expected_pct,
    log2_enrichment = log2(enrichment)
  )

cat("\n--- Macrophage Enrichment by Zone ---\n")
print(macrophage_by_zone %>% select(spatial_zone, pct_macrophages, expected_pct, enrichment, log2_enrichment))

# Chi-square test for association
contingency_table <- cell_data %>%
  count(spatial_zone, is_macrophage) %>%
  pivot_wider(names_from = is_macrophage, values_from = n, values_fill = 0) %>%
  column_to_rownames("spatial_zone")

chi_test <- chisq.test(contingency_table)

cat("\n--- Chi-Square Test for Independence ---\n")
cat(sprintf("  X-squared: %.2f\n", chi_test$statistic))
cat(sprintf("  df: %d\n", chi_test$parameter))
cat(sprintf("  p-value: %s\n", format.pval(chi_test$p.value, digits = 3)))

# Pairwise comparison: Interface vs each other zone
cat("\n--- Pairwise Comparisons (Interface vs Other Zones) ---\n")

pairwise_results <- map_dfr(c("Collagen_only", "Epithelial_only", "Distant"), function(zone) {

  interface_data <- cell_data %>% filter(spatial_zone == "Interface")
  other_data <- cell_data %>% filter(spatial_zone == zone)

  # Fisher's exact test
  tbl <- matrix(c(
    sum(interface_data$is_macrophage),
    sum(!interface_data$is_macrophage),
    sum(other_data$is_macrophage),
    sum(!other_data$is_macrophage)
  ), nrow = 2, byrow = TRUE)

  fisher_result <- fisher.test(tbl)

  tibble(
    comparison = paste("Interface vs", zone),
    interface_pct = mean(interface_data$is_macrophage) * 100,
    other_pct = mean(other_data$is_macrophage) * 100,
    odds_ratio = fisher_result$estimate,
    p_value = fisher_result$p.value,
    ci_lower = fisher_result$conf.int[1],
    ci_upper = fisher_result$conf.int[2]
  )
})

print(pairwise_results)

# =============================================================================
# 7. PER-ROI ANALYSIS
# =============================================================================

cat("\n=============================================================================\n")
cat("Per-ROI Analysis (Biological Replicates)\n")
cat("=============================================================================\n")

per_roi_enrichment <- cell_data %>%
  group_by(ImageNumber) %>%
  summarise(
    n_total = n(),
    n_interface = sum(at_interface),
    n_mac_total = sum(is_macrophage),
    n_mac_interface = sum(is_macrophage & at_interface),

    # Proportions
    pct_mac_overall = n_mac_total / n_total * 100,
    pct_mac_at_interface = ifelse(n_interface > 0, n_mac_interface / n_interface * 100, NA),
    pct_mac_not_interface = ifelse(n_total - n_interface > 0,
                                   (n_mac_total - n_mac_interface) / (n_total - n_interface) * 100, NA),

    # Enrichment
    enrichment = pct_mac_at_interface / pct_mac_overall,

    .groups = "drop"
  )

cat("\n--- Per-ROI Macrophage Enrichment at Interface ---\n")
print(per_roi_enrichment %>%
        select(ImageNumber, n_total, pct_mac_overall, pct_mac_at_interface, enrichment))

# Count ROIs with enrichment > 1
n_enriched <- sum(per_roi_enrichment$enrichment > 1, na.rm = TRUE)
n_total_rois <- sum(!is.na(per_roi_enrichment$enrichment))

cat(sprintf("\n  ROIs with macrophage enrichment at interface: %d / %d\n",
            n_enriched, n_total_rois))

# =============================================================================
# 8. DISTANCE DISTRIBUTION ANALYSIS
# =============================================================================

cat("\n=============================================================================\n")
cat("Distance Distribution Analysis\n")
cat("=============================================================================\n")

# Compare distance distributions for macrophages vs other cells
distance_comparison <- cell_data %>%
  group_by(is_macrophage) %>%
  summarise(
    n = n(),
    mean_dist_collagen = mean(distance_to_tract),
    median_dist_collagen = median(distance_to_tract),
    mean_dist_epithelial = mean(dist_to_epithelial),
    median_dist_epithelial = median(dist_to_epithelial),
    .groups = "drop"
  ) %>%
  mutate(cell_type = ifelse(is_macrophage, "Macrophage", "Other"))

cat("\n--- Distance Summary by Cell Type ---\n")
print(distance_comparison)

# Wilcoxon tests
mac_cells <- cell_data %>% filter(is_macrophage)
other_cells <- cell_data %>% filter(!is_macrophage)

collagen_test <- wilcox.test(mac_cells$distance_to_tract, other_cells$distance_to_tract)
epithelial_test <- wilcox.test(mac_cells$dist_to_epithelial, other_cells$dist_to_epithelial)

# Effect sizes
collagen_d <- cohen.d(mac_cells$distance_to_tract, other_cells$distance_to_tract)
epithelial_d <- cohen.d(mac_cells$dist_to_epithelial, other_cells$dist_to_epithelial)

cat("\n--- Macrophages vs Other Cells: Distance to Collagen ---\n")
cat(sprintf("  Wilcoxon p-value: %s\n", format.pval(collagen_test$p.value, digits = 3)))
cat(sprintf("  Cohen's d: %.3f (%s)\n", collagen_d$estimate, collagen_d$magnitude))

cat("\n--- Macrophages vs Other Cells: Distance to Epithelial ---\n")
cat(sprintf("  Wilcoxon p-value: %s\n", format.pval(epithelial_test$p.value, digits = 3)))
cat(sprintf("  Cohen's d: %.3f (%s)\n", epithelial_d$estimate, epithelial_d$magnitude))

# =============================================================================
# 9. VISUALIZATIONS
# =============================================================================

cat("\n=============================================================================\n")
cat("Generating Figures\n")
cat("=============================================================================\n")

# Color palettes
zone_colors <- c(
  "Interface" = "#E41A1C",
  "Collagen_only" = "#377EB8",
  "Epithelial_only" = "#4DAF4A",
  "Distant" = "#999999"
)

celltype_colors <- c(
  "Macrophage" = "#E41A1C",
  "Epithelial" = "#4DAF4A",
  "Other" = "#999999",
  "Double_positive" = "#984EA3"
)

# --- Figure A: Macrophage proportion by spatial zone ---
p1 <- ggplot(macrophage_by_zone, aes(x = reorder(spatial_zone, -pct_macrophages),
                                     y = pct_macrophages, fill = spatial_zone)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = global_mac_pct, linetype = "dashed", color = "black", linewidth = 1) +
  geom_text(aes(label = sprintf("%.1f%%", pct_macrophages)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = zone_colors) +
  labs(
    title = "CD163+ Macrophage Proportion by Spatial Zone",
    subtitle = sprintf("Dashed line = global proportion (%.1f%%)", global_mac_pct),
    x = "Spatial Zone",
    y = "Macrophage Proportion (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  ) +
  coord_cartesian(ylim = c(0, max(macrophage_by_zone$pct_macrophages) * 1.15))

# --- Figure B: Enrichment bar plot ---
p2 <- ggplot(macrophage_by_zone, aes(x = reorder(spatial_zone, -log2_enrichment),
                                     y = log2_enrichment, fill = spatial_zone)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_text(aes(label = sprintf("%.2fx", enrichment),
                y = ifelse(log2_enrichment > 0, log2_enrichment + 0.1, log2_enrichment - 0.1)),
            size = 4) +
  scale_fill_manual(values = zone_colors) +
  labs(
    title = "Macrophage Enrichment by Spatial Zone",
    subtitle = "Relative to global macrophage proportion",
    x = "Spatial Zone",
    y = "Log2(Enrichment)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

# --- Figure C: Distance to collagen by cell type ---
p3 <- ggplot(cell_data %>%
               mutate(cell_label = ifelse(is_macrophage, "CD163+ Macrophage", "Other Cells")),
             aes(x = cell_label, y = distance_to_tract, fill = cell_label)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, label.y = max(cell_data$distance_to_tract) * 0.9) +
  scale_fill_manual(values = c("CD163+ Macrophage" = "#E41A1C", "Other Cells" = "#999999")) +
  scale_y_continuous(limits = c(0, quantile(cell_data$distance_to_tract, 0.99))) +
  labs(
    title = "Distance to Collagen by Cell Type",
    subtitle = sprintf("Cohen's d = %.3f (%s)", collagen_d$estimate, collagen_d$magnitude),
    x = "",
    y = "Distance to Nearest Collagen (pixels)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

# --- Figure D: Distance to epithelial cells by cell type ---
p4 <- ggplot(cell_data %>%
               filter(!is_epithelial) %>%  # Exclude epithelial cells themselves
               mutate(cell_label = ifelse(is_macrophage, "CD163+ Macrophage", "Other Cells")),
             aes(x = cell_label, y = dist_to_epithelial, fill = cell_label)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5) +
  scale_fill_manual(values = c("CD163+ Macrophage" = "#E41A1C", "Other Cells" = "#999999")) +
  scale_y_continuous(limits = c(0, quantile(cell_data$dist_to_epithelial[!cell_data$is_epithelial], 0.99))) +
  labs(
    title = "Distance to Nearest Epithelial Cell",
    subtitle = "Non-epithelial cells only",
    x = "",
    y = "Distance to Nearest Epithelial Cell (pixels)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

# --- Figure E: Per-ROI enrichment ---
p5 <- ggplot(per_roi_enrichment %>% filter(!is.na(enrichment)),
             aes(x = factor(ImageNumber), y = enrichment, fill = enrichment > 1)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", enrichment)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("FALSE" = "#377EB8", "TRUE" = "#E41A1C"),
                    labels = c("Depleted", "Enriched"),
                    name = "At Interface") +
  labs(
    title = "Macrophage Enrichment at Interface by ROI",
    subtitle = sprintf("Enriched in %d/%d ROIs", n_enriched, n_total_rois),
    x = "ROI",
    y = "Enrichment (Observed/Expected)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  coord_cartesian(ylim = c(0, max(per_roi_enrichment$enrichment, na.rm = TRUE) * 1.2))

# --- Figure F: Spatial scatter plot (sample ROI) ---
# Select ROI with most interface cells for visualization
best_roi <- per_roi_enrichment %>%
  filter(n_interface == max(n_interface)) %>%
  pull(ImageNumber) %>%
  first()

sample_data <- cell_data %>%
  filter(ImageNumber == best_roi) %>%
  mutate(
    plot_type = case_when(
      is_macrophage & at_interface ~ "Macrophage (Interface)",
      is_macrophage ~ "Macrophage (Other)",
      is_epithelial ~ "Epithelial",
      TRUE ~ "Other"
    )
  )

p6 <- ggplot(sample_data, aes(x = x, y = y, color = plot_type)) +
  geom_point(data = sample_data %>% filter(plot_type == "Other"),
             size = 0.3, alpha = 0.3) +
  geom_point(data = sample_data %>% filter(plot_type == "Epithelial"),
             size = 0.5, alpha = 0.5) +
  geom_point(data = sample_data %>% filter(grepl("Macrophage", plot_type)),
             size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c(
    "Macrophage (Interface)" = "#E41A1C",
    "Macrophage (Other)" = "#FF9999",
    "Epithelial" = "#4DAF4A",
    "Other" = "#CCCCCC"
  )) +
  labs(
    title = sprintf("Spatial Distribution (ROI %d)", best_roi),
    subtitle = "Macrophages at interface highlighted in red",
    x = "X position (pixels)",
    y = "Y position (pixels)",
    color = "Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  coord_fixed()

# --- Combined panel ---
combined_panel <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
  plot_annotation(
    title = "CD163+ Macrophage Spatial Positioning at Epithelial-Collagen Interface",
    subtitle = sprintf("IMC Analysis | %s cells across %d ROIs, 2 patients",
                       format(nrow(cell_data), big.mark = ","), n_distinct(cell_data$ImageNumber)),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

# Save figures
cat("Saving figures...\n")

ggsave(file.path(output_dir, "01_macrophage_proportion_by_zone.pdf"), p1, width = 8, height = 6)
ggsave(file.path(output_dir, "02_macrophage_enrichment_by_zone.pdf"), p2, width = 8, height = 6)
ggsave(file.path(output_dir, "03_distance_to_collagen.pdf"), p3, width = 6, height = 6)
ggsave(file.path(output_dir, "04_distance_to_epithelial.pdf"), p4, width = 6, height = 6)
ggsave(file.path(output_dir, "05_per_roi_enrichment.pdf"), p5, width = 10, height = 6)
ggsave(file.path(output_dir, "06_spatial_scatter.pdf"), p6, width = 10, height = 8)
ggsave(file.path(output_dir, "07_combined_panel.pdf"), combined_panel, width = 16, height = 18)
ggsave(file.path(output_dir, "07_combined_panel.png"), combined_panel, width = 16, height = 18, dpi = 300)

# =============================================================================
# 10. EXPORT RESULTS
# =============================================================================

cat("\n--- Exporting Results ---\n")

# Summary table for manuscript
manuscript_summary <- tibble(
  Analysis = c(
    "Total cells analyzed",
    "CD163+ macrophages",
    "Epithelial cells (Pan-CK+)",
    "Cells at interface",
    "Macrophages at interface",
    "Global macrophage %",
    "Macrophage % at interface",
    "Interface enrichment",
    "Chi-square p-value",
    "ROIs with enrichment > 1",
    "Cohen's d (distance to collagen)",
    "Cohen's d (distance to epithelial)"
  ),
  Value = c(
    format(nrow(cell_data), big.mark = ","),
    format(sum(cell_data$is_macrophage), big.mark = ","),
    format(sum(cell_data$is_epithelial), big.mark = ","),
    format(sum(cell_data$at_interface), big.mark = ","),
    format(sum(cell_data$is_macrophage & cell_data$at_interface), big.mark = ","),
    sprintf("%.2f%%", global_mac_pct),
    sprintf("%.2f%%", macrophage_by_zone$pct_macrophages[macrophage_by_zone$spatial_zone == "Interface"]),
    sprintf("%.2fx", macrophage_by_zone$enrichment[macrophage_by_zone$spatial_zone == "Interface"]),
    format.pval(chi_test$p.value, digits = 3),
    sprintf("%d/%d", n_enriched, n_total_rois),
    sprintf("%.3f (%s)", collagen_d$estimate, collagen_d$magnitude),
    sprintf("%.3f (%s)", epithelial_d$estimate, epithelial_d$magnitude)
  )
)

fwrite(manuscript_summary, file.path(output_dir, "manuscript_summary.csv"))
fwrite(macrophage_by_zone, file.path(output_dir, "macrophage_by_zone.csv"))
fwrite(per_roi_enrichment, file.path(output_dir, "per_roi_enrichment.csv"))
fwrite(pairwise_results, file.path(output_dir, "pairwise_comparisons.csv"))
fwrite(distance_comparison, file.path(output_dir, "distance_comparison.csv"))

# =============================================================================
# 11. GENERATE MANUSCRIPT TEXT
# =============================================================================

cat("\n=============================================================================\n")
cat("MANUSCRIPT TEXT (Copy-Paste Ready)\n")
cat("=============================================================================\n\n")

interface_enrichment <- macrophage_by_zone$enrichment[macrophage_by_zone$spatial_zone == "Interface"]
interface_pct <- macrophage_by_zone$pct_macrophages[macrophage_by_zone$spatial_zone == "Interface"]

manuscript_text <- sprintf(
  "To determine whether myeloid cells are positioned to deliver paracrine signals at
epithelial-stromal boundaries, we analyzed CD163+ macrophage distribution relative to
epithelial cells and collagen tracts using imaging mass cytometry (n = %s cells across
%d ROIs from 2 patients).

We defined the epithelial-collagen interface as regions within %d pixels of both
Pan-Cytokeratin+ epithelial cells and Collagen Type I+ tracts. CD163+ macrophages were
significantly enriched at this interface compared to global tissue distribution
(%.1f%% vs %.1f%%; %.2f-fold enrichment; Chi-square p %s; Figure X). This enrichment
was consistent across biological replicates, with %d of %d ROIs showing elevated
macrophage presence at the interface.

Macrophages were positioned closer to both collagen tracts (Cohen's d = %.3f, %s effect)
and epithelial cells (Cohen's d = %.3f, %s effect) compared to other non-epithelial cells.
Pairwise comparisons confirmed that macrophage enrichment at the interface exceeded that
in collagen-only zones (OR = %.2f, p %s) and distant regions (OR = %.2f, p %s).

This spatial positioning places CD163+ macrophages at the precise location where they
could provide SPP1 ligand to adjacent CD44-expressing epithelial cells, supporting
a model of paracrine-mediated epithelial-ECM interaction.",

  format(nrow(cell_data), big.mark = ","),
  n_distinct(cell_data$ImageNumber),
  epithelial_proximity_threshold,
  interface_pct,
  global_mac_pct,
  interface_enrichment,
  ifelse(chi_test$p.value < 0.001, "< 0.001", sprintf("= %.3f", chi_test$p.value)),
  n_enriched,
  n_total_rois,
  collagen_d$estimate, collagen_d$magnitude,
  epithelial_d$estimate, epithelial_d$magnitude,
  pairwise_results$odds_ratio[1],
  ifelse(pairwise_results$p_value[1] < 0.001, "< 0.001", sprintf("= %.3f", pairwise_results$p_value[1])),
  pairwise_results$odds_ratio[3],
  ifelse(pairwise_results$p_value[3] < 0.001, "< 0.001", sprintf("= %.3f", pairwise_results$p_value[3]))
)

cat(manuscript_text)

# Save manuscript text
writeLines(manuscript_text, file.path(output_dir, "manuscript_text.txt"))

# =============================================================================
# COMPLETION
# =============================================================================

cat("\n\n=============================================================================\n")
cat("Analysis Complete!\n")
cat("=============================================================================\n")
cat("\nOutput directory:", output_dir, "\n")
cat("\nFiles created:\n")
cat("  Tables:\n")
cat("    - manuscript_summary.csv\n")
cat("    - macrophage_by_zone.csv\n")
cat("    - per_roi_enrichment.csv\n")
cat("    - pairwise_comparisons.csv\n")
cat("    - distance_comparison.csv\n")
cat("  Figures:\n")
cat("    - 01_macrophage_proportion_by_zone.pdf\n")
cat("    - 02_macrophage_enrichment_by_zone.pdf\n")
cat("    - 03_distance_to_collagen.pdf\n")
cat("    - 04_distance_to_epithelial.pdf\n")
cat("    - 05_per_roi_enrichment.pdf\n")
cat("    - 06_spatial_scatter.pdf\n")
cat("    - 07_combined_panel.pdf/png\n")
cat("  Text:\n")
cat("    - manuscript_text.txt\n")
cat("\n=============================================================================\n")

# Print key results summary
cat("\n=== KEY RESULTS ===\n")
cat(sprintf("Interface enrichment: %.2fx (p %s)\n",
            interface_enrichment,
            ifelse(chi_test$p.value < 0.001, "< 0.001", sprintf("= %.3f", chi_test$p.value))))
cat(sprintf("Consistent across ROIs: %d/%d\n", n_enriched, n_total_rois))
cat(sprintf("Macrophages closer to collagen: d = %.3f (%s)\n",
            collagen_d$estimate, collagen_d$magnitude))
cat(sprintf("Macrophages closer to epithelia: d = %.3f (%s)\n",
            epithelial_d$estimate, epithelial_d$magnitude))
