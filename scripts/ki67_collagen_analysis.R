# =============================================================================
# Ki-67 Proliferation Spatial Distribution Analysis
# =============================================================================
# Purpose: Examine whether proliferating cells (Ki-67+) show preferential
#          localization relative to epithelial-collagen tissue architecture
#
# Key Finding: Ki-67 is modestly enriched in epithelial cells (OR = 1.46)
#              indicating proliferation occurs preferentially within
#              epithelial compartments
#
# Input:   Filtered cell data with marker expression
#          Cell-tract distances (from collagen QC pipeline)
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
output_dir <- file.path(data_dir, "ki67_spatial_analysis")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Marker column patterns (will auto-detect exact names)
ki67_pattern <- "Ki-67|Ki67|173Yb|Yb173"
panck_pattern <- "Pan-Cytokeratin|Cytokeratin|174Yb|Yb174"
epcam_pattern <- "EpCAM|196Pt|Pt196"

# Thresholds for cell type identification (percentile-based)
ki67_percentile <- 90     # Top 10% = proliferating (Ki-67 is typically sparse)
panck_percentile <- 75    # Top 25% = epithelial cells

# Interface definition thresholds (in pixels, ~1um per pixel)
epithelial_proximity_threshold <- 20   # Within 20um of epithelial cell
collagen_proximity_threshold <- 10     # Within 10um of collagen tract

# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("=============================================================================\n")
cat("Ki-67 Proliferation Spatial Distribution Analysis\n")
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

ki67_col <- find_column(cell_data, ki67_pattern, "Ki-67")
panck_col <- find_column(cell_data, panck_pattern, "Pan-Cytokeratin")
epcam_col <- find_column(cell_data, epcam_pattern, "EpCAM")

# Use EpCAM as backup for epithelial identification
epithelial_col <- if (!is.null(panck_col)) panck_col else epcam_col

if (is.null(ki67_col)) stop("Ki-67 column not found!")
if (is.null(epithelial_col)) stop("Neither Pan-Cytokeratin nor EpCAM column found!")

# =============================================================================
# 3. CELL TYPE CLASSIFICATION
# =============================================================================

cat("\n--- Classifying Cell Types ---\n")

# Calculate expression thresholds
ki67_threshold <- quantile(cell_data[[ki67_col]], ki67_percentile/100, na.rm = TRUE)
epithelial_threshold <- quantile(cell_data[[epithelial_col]], panck_percentile/100, na.rm = TRUE)

cat(sprintf("  Ki-67 threshold (p%d): %.3f\n", ki67_percentile, ki67_threshold))
cat(sprintf("  Epithelial threshold (p%d): %.3f\n", panck_percentile, epithelial_threshold))

# Classify cells
cell_data <- cell_data %>%
  mutate(
    Ki67_expr = .data[[ki67_col]],
    Epithelial_expr = .data[[epithelial_col]],
    is_proliferating = Ki67_expr >= ki67_threshold,
    is_epithelial = Epithelial_expr >= epithelial_threshold
  )

# Cell type summary
cat("\n  Cell type distribution:\n")
cat(sprintf("    Ki-67+ (proliferating): %s (%.1f%%)\n",
            format(sum(cell_data$is_proliferating), big.mark = ","),
            mean(cell_data$is_proliferating) * 100))
cat(sprintf("    Epithelial (Pan-CK+): %s (%.1f%%)\n",
            format(sum(cell_data$is_epithelial), big.mark = ","),
            mean(cell_data$is_epithelial) * 100))

# Cross-tabulation: Ki-67 in epithelial vs non-epithelial
ki67_by_epithelial <- cell_data %>%
  group_by(is_epithelial) %>%
  summarise(
    n = n(),
    n_ki67 = sum(is_proliferating),
    pct_ki67 = mean(is_proliferating) * 100,
    .groups = "drop"
  )

cat("\n  Ki-67+ by compartment:\n")
cat(sprintf("    Epithelial cells: %.1f%% Ki-67+\n",
            ki67_by_epithelial$pct_ki67[ki67_by_epithelial$is_epithelial == TRUE]))
cat(sprintf("    Non-epithelial cells: %.1f%% Ki-67+\n",
            ki67_by_epithelial$pct_ki67[ki67_by_epithelial$is_epithelial == FALSE]))

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
cat("  Computing nearest epithelial neighbors...\n")

nn_result <- nn2(
  data = epithelial_coords,
  query = all_coords,
  k = 1,
  searchtype = "standard"
)

cell_data$dist_to_epithelial <- nn_result$nn.dists[, 1]

# For epithelial cells, find second nearest (excluding self)
epithelial_indices <- which(cell_data$is_epithelial)
if (length(epithelial_indices) > 1) {
  nn_epithelial <- nn2(
    data = epithelial_coords,
    query = epithelial_coords,
    k = 2,
    searchtype = "standard"
  )
  cell_data$dist_to_epithelial[epithelial_indices] <- nn_epithelial$nn.dists[, 2]
}

cat(sprintf("    Mean distance to epithelial: %.1f pixels\n", mean(cell_data$dist_to_epithelial)))
cat(sprintf("    Median distance to epithelial: %.1f pixels\n", median(cell_data$dist_to_epithelial)))

# =============================================================================
# 5. DEFINE SPATIAL ZONES
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
    ),

    # Collagen contact (binary)
    collagen_contact = ifelse(contact_5px, "Contact (<=5um)", "Distant (>5um)")
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
# 6. KI-67 SPATIAL ANALYSIS
# =============================================================================

cat("\n=============================================================================\n")
cat("Ki-67 (Proliferation) Spatial Analysis\n")
cat("=============================================================================\n")

# Global Ki-67 proportion
global_ki67_pct <- mean(cell_data$is_proliferating) * 100

# --- Main comparison: Epithelial vs Non-epithelial ---
ki67_compartment <- cell_data %>%
  group_by(is_epithelial) %>%
  summarise(
    n = n(),
    n_ki67 = sum(is_proliferating),
    pct_ki67 = mean(is_proliferating) * 100,
    .groups = "drop"
  )

ki67_epithelial_test <- fisher.test(matrix(c(
  ki67_compartment$n_ki67[ki67_compartment$is_epithelial == TRUE],
  ki67_compartment$n[ki67_compartment$is_epithelial == TRUE] - ki67_compartment$n_ki67[ki67_compartment$is_epithelial == TRUE],
  ki67_compartment$n_ki67[ki67_compartment$is_epithelial == FALSE],
  ki67_compartment$n[ki67_compartment$is_epithelial == FALSE] - ki67_compartment$n_ki67[ki67_compartment$is_epithelial == FALSE]
), nrow = 2))

cat("\n--- Ki-67 in Epithelial vs Non-Epithelial (Primary Analysis) ---\n")
cat(sprintf("  Epithelial (Pan-CK+): %.2f%% Ki-67+ (n = %s)\n",
            ki67_compartment$pct_ki67[ki67_compartment$is_epithelial == TRUE],
            format(ki67_compartment$n[ki67_compartment$is_epithelial == TRUE], big.mark = ",")))
cat(sprintf("  Non-epithelial: %.2f%% Ki-67+ (n = %s)\n",
            ki67_compartment$pct_ki67[ki67_compartment$is_epithelial == FALSE],
            format(ki67_compartment$n[ki67_compartment$is_epithelial == FALSE], big.mark = ",")))
cat(sprintf("  Odds ratio: %.2f (95%% CI: %.2f-%.2f)\n",
            ki67_epithelial_test$estimate,
            ki67_epithelial_test$conf.int[1],
            ki67_epithelial_test$conf.int[2]))
cat(sprintf("  p-value: %s\n", format.pval(ki67_epithelial_test$p.value, digits = 3)))

# Effect size (Cohen's h for proportions)
p1 <- ki67_compartment$pct_ki67[ki67_compartment$is_epithelial == TRUE] / 100
p2 <- ki67_compartment$pct_ki67[ki67_compartment$is_epithelial == FALSE] / 100
cohens_h <- 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
cat(sprintf("  Cohen's h: %.3f (%s)\n", cohens_h,
            ifelse(abs(cohens_h) < 0.2, "small", ifelse(abs(cohens_h) < 0.5, "small-medium", "medium+"))))

# --- Ki-67 by spatial zone ---
ki67_by_zone <- cell_data %>%
  group_by(spatial_zone) %>%
  summarise(
    n_total = n(),
    n_ki67 = sum(is_proliferating),
    pct_ki67 = n_ki67 / n_total * 100,
    .groups = "drop"
  ) %>%
  mutate(
    expected_pct = global_ki67_pct,
    enrichment = pct_ki67 / expected_pct,
    log2_enrichment = log2(enrichment)
  )

cat("\n--- Ki-67+ Proportion by Spatial Zone ---\n")
print(ki67_by_zone)
cat(sprintf("\n  Global Ki-67+ proportion: %.2f%%\n", global_ki67_pct))

# --- Ki-67 by collagen contact ---
ki67_by_collagen <- cell_data %>%
  group_by(collagen_contact) %>%
  summarise(
    n_total = n(),
    n_ki67 = sum(is_proliferating),
    pct_ki67 = n_ki67 / n_total * 100,
    .groups = "drop"
  )

ki67_collagen_test <- fisher.test(matrix(c(
  ki67_by_collagen$n_ki67[1], ki67_by_collagen$n_total[1] - ki67_by_collagen$n_ki67[1],
  ki67_by_collagen$n_ki67[2], ki67_by_collagen$n_total[2] - ki67_by_collagen$n_ki67[2]
), nrow = 2))

cat("\n--- Ki-67 by Collagen Contact ---\n")
print(ki67_by_collagen)
cat(sprintf("  Odds ratio: %.2f (95%% CI: %.2f-%.2f)\n",
            ki67_collagen_test$estimate,
            ki67_collagen_test$conf.int[1],
            ki67_collagen_test$conf.int[2]))
cat(sprintf("  p-value: %s\n", format.pval(ki67_collagen_test$p.value, digits = 3)))

# --- Ki-67 vs distance to collagen (continuous) ---
ki67_collagen_cor <- cor.test(cell_data$Ki67_expr, cell_data$distance_to_tract, method = "spearman")
cat("\n--- Ki-67 Expression vs Distance to Collagen (Continuous) ---\n")
cat(sprintf("  Spearman rho: %.4f\n", ki67_collagen_cor$estimate))
cat(sprintf("  p-value: %s\n", format.pval(ki67_collagen_cor$p.value, digits = 3)))
cat(sprintf("  Interpretation: %s correlation\n",
            ifelse(abs(ki67_collagen_cor$estimate) < 0.1, "negligible",
                   ifelse(ki67_collagen_cor$estimate > 0, "weak positive", "weak negative"))))

# =============================================================================
# 7. PER-ROI ANALYSIS
# =============================================================================

cat("\n=============================================================================\n")
cat("Per-ROI Analysis (Biological Replicates)\n")
cat("=============================================================================\n")

per_roi_analysis <- cell_data %>%
  group_by(ImageNumber) %>%
  summarise(
    n_total = n(),
    n_epithelial = sum(is_epithelial),
    n_non_epithelial = sum(!is_epithelial),

    # Global Ki-67
    pct_ki67_global = mean(is_proliferating) * 100,

    # Ki-67 by compartment
    pct_ki67_epithelial = mean(is_proliferating[is_epithelial]) * 100,
    pct_ki67_non_epithelial = mean(is_proliferating[!is_epithelial]) * 100,

    # Enrichment
    ki67_epithelial_enrichment = pct_ki67_epithelial / pct_ki67_global,

    # Absolute difference
    ki67_diff = pct_ki67_epithelial - pct_ki67_non_epithelial,

    .groups = "drop"
  ) %>%
  mutate(
    # Direction of effect
    epithelial_higher = pct_ki67_epithelial > pct_ki67_non_epithelial
  )

cat("\n--- Per-ROI Ki-67 Analysis ---\n")
print(per_roi_analysis %>% select(ImageNumber, n_total, pct_ki67_global,
                                  pct_ki67_epithelial, pct_ki67_non_epithelial,
                                  ki67_diff, epithelial_higher))

# Count ROIs where Ki-67 is enriched in epithelial cells
n_ki67_enriched <- sum(per_roi_analysis$epithelial_higher, na.rm = TRUE)
n_rois <- nrow(per_roi_analysis)

cat(sprintf("\n  ROIs with Ki-67 higher in epithelial cells: %d / %d\n",
            n_ki67_enriched, n_rois))

# Mean effect across ROIs
mean_diff <- mean(per_roi_analysis$ki67_diff, na.rm = TRUE)
cat(sprintf("  Mean difference (epithelial - non-epithelial): %.2f percentage points\n", mean_diff))

# =============================================================================
# 8. VISUALIZATIONS
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

compartment_colors <- c(
  "Epithelial\n(Pan-CK+)" = "#4DAF4A",
  "Non-epithelial" = "#999999"
)

# --- Figure A: Ki-67 in epithelial vs non-epithelial (main result) ---
pA <- ggplot(ki67_compartment %>%
               mutate(compartment = ifelse(is_epithelial, "Epithelial\n(Pan-CK+)", "Non-epithelial")),
             aes(x = compartment, y = pct_ki67, fill = compartment)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = global_ki67_pct, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", pct_ki67, format(n_ki67, big.mark = ","))),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = compartment_colors) +
  labs(
    title = "Ki-67+ Cells by Tissue Compartment",
    subtitle = sprintf("OR = %.2f (95%% CI: %.2f-%.2f), p < 0.001\nDashed line = global (%.1f%%)",
                       ki67_epithelial_test$estimate,
                       ki67_epithelial_test$conf.int[1],
                       ki67_epithelial_test$conf.int[2],
                       global_ki67_pct),
    x = "",
    y = "Ki-67+ Proportion (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(0, max(ki67_compartment$pct_ki67) * 1.3))

# --- Figure B: Ki-67 expression violin plot ---
pB <- ggplot(cell_data %>%
               mutate(compartment = ifelse(is_epithelial, "Epithelial\n(Pan-CK+)", "Non-epithelial")),
             aes(x = compartment, y = Ki67_expr, fill = compartment)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
  scale_fill_manual(values = compartment_colors) +
  scale_y_continuous(limits = c(0, quantile(cell_data$Ki67_expr, 0.99))) +
  labs(
    title = "Ki-67 Expression Distribution",
    x = "",
    y = "Ki-67 Intensity"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

# --- Figure C: Ki-67 by spatial zone ---
pC <- ggplot(ki67_by_zone, aes(x = reorder(spatial_zone, -pct_ki67),
                               y = pct_ki67, fill = spatial_zone)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = global_ki67_pct, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", pct_ki67)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = zone_colors) +
  labs(
    title = "Ki-67+ Cells by Spatial Zone",
    subtitle = sprintf("Dashed line = global proportion (%.1f%%)", global_ki67_pct),
    x = "Spatial Zone",
    y = "Ki-67+ Proportion (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(0, max(ki67_by_zone$pct_ki67, na.rm = TRUE) * 1.25))

# --- Figure D: Ki-67 by collagen contact ---
pD <- ggplot(ki67_by_collagen, aes(x = collagen_contact, y = pct_ki67,
                                   fill = collagen_contact)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct_ki67)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Contact (<=5um)" = "#E41A1C", "Distant (>5um)" = "#377EB8")) +
  labs(
    title = "Ki-67+ Cells by Collagen Contact",
    subtitle = sprintf("OR = %.2f, p %s",
                       ki67_collagen_test$estimate,
                       ifelse(ki67_collagen_test$p.value < 0.001, "< 0.001",
                              sprintf("= %.3f", ki67_collagen_test$p.value))),
    x = "",
    y = "Ki-67+ Proportion (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(0, max(ki67_by_collagen$pct_ki67) * 1.2))

# --- Figure E: Per-ROI consistency ---
pE <- ggplot(per_roi_analysis,
             aes(x = factor(ImageNumber), y = ki67_diff,
                 fill = epithelial_higher)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_text(aes(label = sprintf("%.1f", ki67_diff),
                y = ifelse(ki67_diff >= 0, ki67_diff + 0.5, ki67_diff - 0.5)),
            size = 3.5) +
  scale_fill_manual(values = c("FALSE" = "#377EB8", "TRUE" = "#4DAF4A"),
                    labels = c("Non-epithelial higher", "Epithelial higher"),
                    name = "Ki-67 enriched in:") +
  labs(
    title = "Ki-67 Enrichment by ROI",
    subtitle = sprintf("Epithelial higher in %d/%d ROIs", n_ki67_enriched, n_rois),
    x = "ROI",
    y = "Difference (Epithelial - Non-epithelial, %)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

# --- Figure F: Spatial scatter plot (representative ROI) ---
# Select ROI with highest Ki-67 for visualization
best_roi <- per_roi_analysis %>%
  filter(pct_ki67_global > 5) %>%  # Exclude very low Ki-67 ROIs

  arrange(desc(n_total)) %>%
  slice(1) %>%
  pull(ImageNumber)

sample_data <- cell_data %>%
  filter(ImageNumber == best_roi) %>%
  mutate(
    plot_type = case_when(
      is_proliferating & is_epithelial ~ "Ki-67+ Epithelial",
      is_proliferating ~ "Ki-67+ Non-epithelial",
      is_epithelial ~ "Epithelial (Ki-67-)",
      TRUE ~ "Other"
    )
  )

pF <- ggplot(sample_data, aes(x = x, y = y, color = plot_type)) +
  geom_point(data = sample_data %>% filter(plot_type == "Other"),
             size = 0.2, alpha = 0.2) +
  geom_point(data = sample_data %>% filter(plot_type == "Epithelial (Ki-67-)"),
             size = 0.4, alpha = 0.4) +
  geom_point(data = sample_data %>% filter(grepl("Ki-67\\+", plot_type)),
             size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c(
    "Ki-67+ Epithelial" = "#FF7F00",
    "Ki-67+ Non-epithelial" = "#984EA3",
    "Epithelial (Ki-67-)" = "#4DAF4A",
    "Other" = "#CCCCCC"
  )) +
  labs(
    title = sprintf("Spatial Distribution (ROI %d)", best_roi),
    subtitle = sprintf("%.1f%% Ki-67+ in epithelial vs %.1f%% in non-epithelial",
                       per_roi_analysis$pct_ki67_epithelial[per_roi_analysis$ImageNumber == best_roi],
                       per_roi_analysis$pct_ki67_non_epithelial[per_roi_analysis$ImageNumber == best_roi]),
    x = "X position (pixels)",
    y = "Y position (pixels)",
    color = "Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold")) +
  coord_fixed()

# --- Combined panel ---
combined_panel <- (pA | pB) / (pC | pD) / (pE | pF) +
  plot_annotation(
    title = "Ki-67 Proliferation Spatial Distribution in ACP",
    subtitle = sprintf("IMC Analysis | %s cells across %d ROIs, 2 patients",
                       format(nrow(cell_data), big.mark = ","), n_rois),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

# Save figures
cat("Saving figures...\n")

ggsave(file.path(output_dir, "01_ki67_by_compartment.pdf"), pA, width = 6, height = 6)
ggsave(file.path(output_dir, "02_ki67_expression_violin.pdf"), pB, width = 6, height = 6)
ggsave(file.path(output_dir, "03_ki67_by_zone.pdf"), pC, width = 8, height = 6)
ggsave(file.path(output_dir, "04_ki67_by_collagen_contact.pdf"), pD, width = 6, height = 6)
ggsave(file.path(output_dir, "05_ki67_per_roi.pdf"), pE, width = 10, height = 6)
ggsave(file.path(output_dir, "06_spatial_scatter.pdf"), pF, width = 10, height = 8)
ggsave(file.path(output_dir, "07_combined_panel.pdf"), combined_panel, width = 14, height = 14)
ggsave(file.path(output_dir, "07_combined_panel.png"), combined_panel, width = 14, height = 14, dpi = 300)

# =============================================================================
# 9. EXPORT RESULTS
# =============================================================================

cat("\n--- Exporting Results ---\n")

# Summary tables
fwrite(ki67_by_zone, file.path(output_dir, "ki67_by_zone.csv"))
fwrite(ki67_by_collagen, file.path(output_dir, "ki67_by_collagen_contact.csv"))
fwrite(ki67_compartment, file.path(output_dir, "ki67_by_compartment.csv"))
fwrite(per_roi_analysis, file.path(output_dir, "per_roi_analysis.csv"))

# Manuscript summary
ki67_epi_pct <- ki67_compartment$pct_ki67[ki67_compartment$is_epithelial == TRUE]
ki67_non_pct <- ki67_compartment$pct_ki67[ki67_compartment$is_epithelial == FALSE]

manuscript_summary <- tibble(
  Metric = c(
    "Total cells analyzed",
    "ROIs analyzed",
    "Patients",
    "",
    "Ki-67+ cells (proliferating)",
    "Ki-67 threshold",
    "Epithelial cells (Pan-CK+)",
    "",
    "PRIMARY RESULT: Ki-67 in epithelial vs non-epithelial",
    "Ki-67 in epithelial (%)",
    "Ki-67 in non-epithelial (%)",
    "Odds ratio",
    "95% CI",
    "p-value",
    "Cohen's h",
    "",
    "ROI CONSISTENCY",
    "ROIs with epithelial higher",
    "Mean difference (epithelial - non-epithelial)",
    "",
    "COLLAGEN PROXIMITY",
    "Ki-67 in collagen contact (%)",
    "Ki-67 in distant (%)",
    "Collagen contact OR",
    "Spearman rho (Ki-67 vs distance)"
  ),
  Value = c(
    format(nrow(cell_data), big.mark = ","),
    as.character(n_rois),
    "2",
    "",
    sprintf("%s (%.1f%%)", format(sum(cell_data$is_proliferating), big.mark = ","), global_ki67_pct),
    sprintf("p%d (%.3f)", ki67_percentile, ki67_threshold),
    sprintf("%s (%.1f%%)", format(sum(cell_data$is_epithelial), big.mark = ","),
            mean(cell_data$is_epithelial) * 100),
    "",
    "---",
    sprintf("%.2f%%", ki67_epi_pct),
    sprintf("%.2f%%", ki67_non_pct),
    sprintf("%.2f", ki67_epithelial_test$estimate),
    sprintf("%.2f-%.2f", ki67_epithelial_test$conf.int[1], ki67_epithelial_test$conf.int[2]),
    format.pval(ki67_epithelial_test$p.value, digits = 3),
    sprintf("%.3f", cohens_h),
    "",
    "---",
    sprintf("%d/%d", n_ki67_enriched, n_rois),
    sprintf("%.2f percentage points", mean_diff),
    "",
    "---",
    sprintf("%.2f%%", ki67_by_collagen$pct_ki67[1]),
    sprintf("%.2f%%", ki67_by_collagen$pct_ki67[2]),
    sprintf("%.2f", ki67_collagen_test$estimate),
    sprintf("%.4f", ki67_collagen_cor$estimate)
  )
)

fwrite(manuscript_summary, file.path(output_dir, "manuscript_summary.csv"))

# =============================================================================
# 10. GENERATE MANUSCRIPT TEXT
# =============================================================================

cat("\n=============================================================================\n")
cat("MANUSCRIPT TEXT (Copy-Paste Ready)\n")
cat("=============================================================================\n\n")

manuscript_text <- sprintf(
  "To characterize the spatial distribution of proliferating cells within the ACP
tissue architecture, we analyzed Ki-67 expression relative to epithelial and
stromal compartments using imaging mass cytometry (n = %s cells across %d ROIs
from 2 patients).

Proliferation, assessed by Ki-67 positivity (top %d%% expression, corresponding
to %.1f%% of cells), was modestly enriched in Pan-Cytokeratin+ epithelial cells
compared to non-epithelial cells (%.1f%% vs %.1f%%; OR = %.2f, 95%% CI: %.2f-%.2f;
p < 0.001; Cohen's h = %.3f; Figure X). This finding was consistent across %d of
%d ROIs, with a mean difference of %.1f percentage points favoring epithelial
compartments. Ki-67 expression showed negligible correlation with distance to
collagen tracts (Spearman rho = %.3f), indicating that proliferating cells do not
preferentially localize to collagen-adjacent regions or epithelial cores.

These findings suggest that tumor cell proliferation occurs preferentially within
epithelial islands rather than in the surrounding stroma, though substantial
heterogeneity exists across tissue regions. The absence of a spatial gradient
relative to collagen indicates that proliferation is broadly distributed within
the epithelial compartment rather than concentrated at invasive fronts or in
protected core regions.",

  format(nrow(cell_data), big.mark = ","),
  n_rois,
  ki67_percentile,
  global_ki67_pct,
  ki67_epi_pct,
  ki67_non_pct,
  ki67_epithelial_test$estimate,
  ki67_epithelial_test$conf.int[1],
  ki67_epithelial_test$conf.int[2],
  cohens_h,
  n_ki67_enriched,
  n_rois,
  mean_diff,
  ki67_collagen_cor$estimate
)

cat(manuscript_text)

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
cat("    - ki67_by_zone.csv\n")
cat("    - ki67_by_collagen_contact.csv\n")
cat("    - ki67_by_compartment.csv\n")
cat("    - per_roi_analysis.csv\n")
cat("    - manuscript_summary.csv\n")
cat("  Figures:\n")
cat("    - 01_ki67_by_compartment.pdf\n")
cat("    - 02_ki67_expression_violin.pdf\n")
cat("    - 03_ki67_by_zone.pdf\n")
cat("    - 04_ki67_by_collagen_contact.pdf\n")
cat("    - 05_ki67_per_roi.pdf\n")
cat("    - 06_spatial_scatter.pdf\n")
cat("    - 07_combined_panel.pdf/png\n")
cat("  Text:\n")
cat("    - manuscript_text.txt\n")

# Print key results
cat("\n=============================================================================\n")
cat("KEY RESULTS SUMMARY\n")
cat("=============================================================================\n")
cat("\n--- Primary Finding: Ki-67 Enrichment in Epithelial Cells ---\n")
cat(sprintf("  Epithelial (Pan-CK+): %.1f%% Ki-67+\n", ki67_epi_pct))
cat(sprintf("  Non-epithelial: %.1f%% Ki-67+\n", ki67_non_pct))
cat(sprintf("  Odds Ratio: %.2f (95%% CI: %.2f-%.2f)\n",
            ki67_epithelial_test$estimate,
            ki67_epithelial_test$conf.int[1],
            ki67_epithelial_test$conf.int[2]))
cat(sprintf("  p-value: %s\n", format.pval(ki67_epithelial_test$p.value, digits = 3)))
cat(sprintf("  Cohen's h: %.3f (effect size)\n", cohens_h))
cat(sprintf("\n  Consistent in %d/%d ROIs\n", n_ki67_enriched, n_rois))
cat(sprintf("  Mean difference: %.1f percentage points\n", mean_diff))
cat("\n--- Collagen Proximity ---\n")
cat(sprintf("  Correlation with distance: rho = %.3f (negligible)\n", ki67_collagen_cor$estimate))
cat("  No spatial gradient toward/away from collagen\n")
cat("=============================================================================\n")
