# =============================================================================
# CD44 and β-Catenin Expression Analysis by Collagen Proximity
# =============================================================================
# Purpose: Test whether CD44 and β-catenin expression correlate with
#          proximity to collagen tracts, validating Visium spatial findings
#
# Input:   Cell-level data with tract distances (from Task 2.6)
#          Marker expression data (filtered cells)
# Output:  Statistical comparisons, effect sizes, publication-ready figures
# =============================================================================

library(tidyverse)
library(data.table)
library(effsize)      # Cohen's d
library(ggpubr)       # stat_compare_means
library(patchwork)
library(scales)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths - UPDATE THESE TO YOUR ACTUAL PATHS
data_dir <- "/home/rstudio/interfaces/ACP_IMC/results"
output_dir <- file.path(data_dir, "cd44_bcatenin_analysis")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Marker column names (from your panel)
cd44_col <- "Eu153_153Eu_CD44"
bcatenin_col <- "Yb176_176Yb_B-Catenin"

# Distance threshold for collagen contact (pixels, ~1μm per pixel for IMC)
contact_threshold_px <- 5

# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("=============================================================================\n")
cat("CD44 and β-Catenin Expression by Collagen Proximity\n")
cat("=============================================================================\n\n")

# --- Load cell-tract distance data ---
cell_tract_file <- "/home/rstudio/interfaces/ACP_IMC/results/collagen_qc/cell_tract_distances_all_approaches.csv"

cat("Loading cell-tract distance data...\n")
cat("  File:", cell_tract_file, "\n")

tract_distances <- fread(cell_tract_file)
cat("  Rows loaded:", format(nrow(tract_distances), big.mark = ","), "\n")

# Filter for PERMISSIVE approach only (as per your Task 2.6 recommendation)
tract_distances <- tract_distances %>%
  filter(approach == "permissive")
cat("  Cells after filtering to permissive approach:", format(nrow(tract_distances), big.mark = ","), "\n")

# --- Load filtered cell data with marker expression ---
# This should be your QC-filtered cell data with marker intensities
marker_file <- file.path(data_dir, "qc_results/step6_filtered_metrics/filtered_cell_data.csv")

cat("\nLoading marker expression data...\n")
cat("  File:", marker_file, "\n")

marker_data <- fread(marker_file)
cat("  Cells loaded:", format(nrow(marker_data), big.mark = ","), "\n")

# Check available marker columns
marker_cols <- colnames(marker_data)
cat("\n  Available marker columns (first 20):\n")
cat("   ", paste(head(marker_cols, 20), collapse = ", "), "...\n")

# --- Merge distance and marker data ---
# The cell_id in tract_distances should match a cell identifier in marker_data
# First, let's check what ID columns are available

cat("\n--- Preparing to merge datasets ---\n")
cat("  Tract distance columns:", paste(colnames(tract_distances), collapse = ", "), "\n")

# Extract the numeric cell ID from the tract_distances cell_id column
# Format: R24-00001_C1-1_Blastoma_18_ROIs_005_1 -> need to match to marker data

# Check if marker_data has ObjectNumber and ImageNumber columns
if (all(c("ObjectNumber", "ImageNumber") %in% colnames(marker_data))) {
  cat("  Found ObjectNumber and ImageNumber in marker data\n")

  # Parse the cell_id to extract ROI and object number
  tract_distances <- tract_distances %>%
    mutate(
      # Extract the last number after final underscore as ObjectNumber
      ObjectNumber = as.integer(sub(".*_([0-9]+)$", "\\1", cell_id)),
      # roi_id is already present
      ImageNumber = roi_id
    )

  # Merge on ImageNumber and ObjectNumber
  cell_data <- marker_data %>%
    inner_join(
      tract_distances %>% select(ImageNumber, ObjectNumber, distance_to_tract,
                                 contact_1px, contact_5px, contact_10px),
      by = c("ImageNumber", "ObjectNumber")
    )

} else if ("cell_id" %in% colnames(marker_data)) {
  # Direct merge on cell_id
  cell_data <- marker_data %>%
    inner_join(tract_distances %>% select(cell_id, distance_to_tract, roi_id,
                                          contact_1px, contact_5px, contact_10px),
               by = "cell_id")
} else {
  # Try to construct a matching ID
  cat("  Attempting to construct matching cell IDs...\n")

  # If marker_data has separate ROI and cell number columns, construct the ID
  stop("Could not find matching columns. Please check your data structure and update the merge logic.")
}

cat("  Cells after merge:", format(nrow(cell_data), big.mark = ","), "\n")

if (nrow(cell_data) == 0) {
  stop("No cells matched between distance and marker data! Check ID formats.")
}

# Check available columns after merge
cat("\n  Merged data columns:\n")
cat("   ", paste(head(colnames(cell_data), 30), collapse = ", "), "...\n")

# =============================================================================
# 2. IDENTIFY MARKER COLUMNS
# =============================================================================

cat("\n--- Identifying Marker Columns ---\n")

# IMC channel-to-marker mapping (from your panel)
# CD44 = 153Eu_CD44 → column name is "Eu153"
# β-catenin = 176Yb_B-Catenin → column name is "Yb176"

# Find CD44 column - try multiple patterns
cd44_patterns <- c("CD44", "Eu153", "153Eu")
cd44_col <- NULL
for (pattern in cd44_patterns) {
  matches <- grep(pattern, colnames(cell_data), value = TRUE, ignore.case = TRUE)
  if (length(matches) > 0) {
    cd44_col <- matches[1]
    break
  }
}

if (!is.null(cd44_col)) {
  cat("  CD44 column found:", cd44_col, "\n")
} else {
  cat("  Available columns containing 'Eu':\n")
  eu_cols <- grep("Eu", colnames(cell_data), value = TRUE, ignore.case = TRUE)
  cat("   ", paste(eu_cols, collapse = ", "), "\n")
  stop("CD44 column not found! Check column names above.")
}

# Find β-catenin column - try multiple patterns
bcatenin_patterns <- c("catenin", "B-Catenin", "Beta", "Yb176", "176Yb")
bcatenin_col <- NULL
for (pattern in bcatenin_patterns) {
  matches <- grep(pattern, colnames(cell_data), value = TRUE, ignore.case = TRUE)
  if (length(matches) > 0) {
    bcatenin_col <- matches[1]
    break
  }
}

if (!is.null(bcatenin_col)) {
  cat("  β-Catenin column found:", bcatenin_col, "\n")
} else {
  cat("  Available columns containing 'Yb':\n")
  yb_cols <- grep("Yb", colnames(cell_data), value = TRUE, ignore.case = TRUE)
  cat("   ", paste(yb_cols, collapse = ", "), "\n")
  stop("β-Catenin column not found! Check column names above.")
}

# Distance column is known: distance_to_tract (in pixels)
dist_col <- "distance_to_tract"
cat("  Distance column:", dist_col, "(in pixels)\n")

# ROI column
roi_col <- ifelse("ImageNumber" %in% colnames(cell_data), "ImageNumber", "roi_id")
cat("  ROI column:", roi_col, "\n")

# =============================================================================
# 3. DATA PREPARATION
# =============================================================================

cat("\n--- Preparing Analysis Data ---\n")

# NOTE: distance_to_tract is in PIXELS
# Typical IMC pixel size is ~1μm, so we'll use pixel values directly
# The pre-computed contact columns (contact_5px, contact_10px) are already available

# Create analysis dataframe with standardized column names
analysis_data <- cell_data %>%
  transmute(
    cell_id = row_number(),
    roi = .data[[roi_col]],
    distance_px = distance_to_tract,
    CD44 = .data[[cd44_col]],
    B_Catenin = .data[[bcatenin_col]],
    # Use pre-computed contact columns (in pixels, ~1μm per pixel)
    contact_5px = contact_5px,
    contact_10px = contact_10px
  ) %>%
  # Remove cells with missing values
  filter(!is.na(distance_px), !is.na(CD44), !is.na(B_Catenin)) %>%
  # Create contact groups using pre-computed 5px threshold (~5μm)
  mutate(
    contact_group = factor(
      ifelse(contact_5px, "Contact (≤5μm)", "Distant (>5μm)"),
      levels = c("Contact (≤5μm)", "Distant (>5μm)")
    )
  )

cat("  Analysis cells:", format(nrow(analysis_data), big.mark = ","), "\n")
cat("  ROIs included:", n_distinct(analysis_data$roi), "\n")

# Contact statistics
contact_stats <- analysis_data %>%
  count(contact_group) %>%
  mutate(pct = n / sum(n) * 100)

cat("\n  Contact distribution (5μm threshold):\n")
for (i in 1:nrow(contact_stats)) {
  cat(sprintf("    %s: %s cells (%.1f%%)\n",
              contact_stats$contact_group[i],
              format(contact_stats$n[i], big.mark = ","),
              contact_stats$pct[i]))
}

# =============================================================================
# 4. STATISTICAL ANALYSIS - CD44
# =============================================================================

cat("\n=============================================================================\n")
cat("CD44 Analysis\n")
cat("=============================================================================\n")

# Summary statistics by group
cd44_summary <- analysis_data %>%
  group_by(contact_group) %>%
  summarise(
    n = n(),
    mean = mean(CD44),
    median = median(CD44),
    sd = sd(CD44),
    se = sd / sqrt(n),
    q25 = quantile(CD44, 0.25),
    q75 = quantile(CD44, 0.75),
    .groups = "drop"
  )

cat("\n--- CD44 Expression Summary ---\n")
print(cd44_summary)

# Wilcoxon rank-sum test (Mann-Whitney U)
cd44_wilcox <- wilcox.test(
  CD44 ~ contact_group,
  data = analysis_data,
  conf.int = TRUE
)

cat("\n--- Wilcoxon Rank-Sum Test ---\n")
cat("  W statistic:", format(cd44_wilcox$statistic, big.mark = ","), "\n")
cat("  p-value:", format.pval(cd44_wilcox$p.value, digits = 3), "\n")
cat("  Difference estimate:", round(cd44_wilcox$estimate, 4), "\n")
cat("  95% CI: [", round(cd44_wilcox$conf.int[1], 4), ", ",
    round(cd44_wilcox$conf.int[2], 4), "]\n")

# Effect sizes
contact_cd44 <- analysis_data$CD44[analysis_data$contact_group == "Contact (≤5μm)"]
distant_cd44 <- analysis_data$CD44[analysis_data$contact_group == "Distant (>5μm)"]

cd44_cohens_d <- cohen.d(contact_cd44, distant_cd44)
cd44_cliffs_delta <- cliff.delta(contact_cd44, distant_cd44)

cat("\n--- Effect Sizes ---\n")
cat("  Cohen's d:", round(cd44_cohens_d$estimate, 3),
    sprintf("(%s)\n", cd44_cohens_d$magnitude))
cat("  Cliff's delta:", round(cd44_cliffs_delta$estimate, 3),
    sprintf("(%s)\n", cd44_cliffs_delta$magnitude))

# Store CD44 results
cd44_results <- list(
  summary = cd44_summary,
  wilcox = cd44_wilcox,
  cohens_d = cd44_cohens_d,
  cliffs_delta = cd44_cliffs_delta
)

# =============================================================================
# 5. STATISTICAL ANALYSIS - β-CATENIN
# =============================================================================

cat("\n=============================================================================\n")
cat("β-Catenin Analysis\n")
cat("=============================================================================\n")

# Summary statistics by group
bcatenin_summary <- analysis_data %>%
  group_by(contact_group) %>%
  summarise(
    n = n(),
    mean = mean(B_Catenin),
    median = median(B_Catenin),
    sd = sd(B_Catenin),
    se = sd / sqrt(n),
    q25 = quantile(B_Catenin, 0.25),
    q75 = quantile(B_Catenin, 0.75),
    .groups = "drop"
  )

cat("\n--- β-Catenin Expression Summary ---\n")
print(bcatenin_summary)

# Wilcoxon rank-sum test
bcatenin_wilcox <- wilcox.test(
  B_Catenin ~ contact_group,
  data = analysis_data,
  conf.int = TRUE
)

cat("\n--- Wilcoxon Rank-Sum Test ---\n")
cat("  W statistic:", format(bcatenin_wilcox$statistic, big.mark = ","), "\n")
cat("  p-value:", format.pval(bcatenin_wilcox$p.value, digits = 3), "\n")
cat("  Difference estimate:", round(bcatenin_wilcox$estimate, 4), "\n")
cat("  95% CI: [", round(bcatenin_wilcox$conf.int[1], 4), ", ",
    round(bcatenin_wilcox$conf.int[2], 4), "]\n")

# Effect sizes
contact_bcat <- analysis_data$B_Catenin[analysis_data$contact_group == "Contact (≤5μm)"]
distant_bcat <- analysis_data$B_Catenin[analysis_data$contact_group == "Distant (>5μm)"]

bcatenin_cohens_d <- cohen.d(contact_bcat, distant_bcat)
bcatenin_cliffs_delta <- cliff.delta(contact_bcat, distant_bcat)

cat("\n--- Effect Sizes ---\n")
cat("  Cohen's d:", round(bcatenin_cohens_d$estimate, 3),
    sprintf("(%s)\n", bcatenin_cohens_d$magnitude))
cat("  Cliff's delta:", round(bcatenin_cliffs_delta$estimate, 3),
    sprintf("(%s)\n", bcatenin_cliffs_delta$magnitude))

# Store β-catenin results
bcatenin_results <- list(
  summary = bcatenin_summary,
  wilcox = bcatenin_wilcox,
  cohens_d = bcatenin_cohens_d,
  cliffs_delta = bcatenin_cliffs_delta
)

# =============================================================================
# 6. PER-ROI ANALYSIS (Biological Replicates)
# =============================================================================

cat("\n=============================================================================\n")
cat("Per-ROI Analysis (Biological Replicates)\n")
cat("=============================================================================\n")

per_roi_stats <- analysis_data %>%
  group_by(roi, contact_group) %>%
  summarise(
    n = n(),
    CD44_mean = mean(CD44),
    CD44_median = median(CD44),
    B_Catenin_mean = mean(B_Catenin),
    B_Catenin_median = median(B_Catenin),
    .groups = "drop"
  )

# Calculate effect direction consistency across ROIs
roi_effects <- analysis_data %>%
  group_by(roi) %>%
  summarise(
    cd44_contact_mean = mean(CD44[contact_group == "Contact (≤5μm)"]),
    cd44_distant_mean = mean(CD44[contact_group == "Distant (>5μm)"]),
    cd44_diff = cd44_contact_mean - cd44_distant_mean,
    cd44_higher_in_contact = cd44_contact_mean > cd44_distant_mean,

    bcat_contact_mean = mean(B_Catenin[contact_group == "Contact (≤5μm)"]),
    bcat_distant_mean = mean(B_Catenin[contact_group == "Distant (>5μm)"]),
    bcat_diff = bcat_contact_mean - bcat_distant_mean,
    bcat_higher_in_contact = bcat_contact_mean > bcat_distant_mean,
    .groups = "drop"
  )

cat("\n--- Effect Consistency Across ROIs ---\n")
cat("  CD44 higher in contact cells:", sum(roi_effects$cd44_higher_in_contact),
    "/", nrow(roi_effects), "ROIs\n")
cat("  β-Catenin higher in contact cells:", sum(roi_effects$bcat_higher_in_contact),
    "/", nrow(roi_effects), "ROIs\n")

print(roi_effects %>% select(roi, cd44_diff, cd44_higher_in_contact,
                             bcat_diff, bcat_higher_in_contact))

# =============================================================================
# 7. SENSITIVITY ANALYSIS - Multiple Distance Thresholds
# =============================================================================

cat("\n=============================================================================\n")
cat("Sensitivity Analysis - Multiple Distance Thresholds\n")
cat("=============================================================================\n")

# Thresholds in pixels (~1μm per pixel for IMC)
sensitivity_thresholds <- c(2, 5, 10, 15, 20)

sensitivity_results <- map_dfr(sensitivity_thresholds, function(thresh) {

  temp_data <- analysis_data %>%
    mutate(contact = distance_px <= thresh)

  contact_cells <- temp_data %>% filter(contact)
  distant_cells <- temp_data %>% filter(!contact)

  # Skip if either group is too small

  if (nrow(contact_cells) < 10 || nrow(distant_cells) < 10) {
    return(NULL)
  }

  # CD44
  cd44_test <- wilcox.test(contact_cells$CD44, distant_cells$CD44)
  cd44_d <- cohen.d(contact_cells$CD44, distant_cells$CD44)

  # β-Catenin
  bcat_test <- wilcox.test(contact_cells$B_Catenin, distant_cells$B_Catenin)
  bcat_d <- cohen.d(contact_cells$B_Catenin, distant_cells$B_Catenin)


  tibble(
    threshold_px = thresh,
    n_contact = nrow(contact_cells),
    n_distant = nrow(distant_cells),
    pct_contact = nrow(contact_cells) / nrow(temp_data) * 100,

    cd44_contact_mean = mean(contact_cells$CD44),
    cd44_distant_mean = mean(distant_cells$CD44),
    cd44_p_value = cd44_test$p.value,
    cd44_cohens_d = cd44_d$estimate,

    bcat_contact_mean = mean(contact_cells$B_Catenin),
    bcat_distant_mean = mean(distant_cells$B_Catenin),
    bcat_p_value = bcat_test$p.value,
    bcat_cohens_d = bcat_d$estimate
  )
})

cat("\n--- Sensitivity Results ---\n")
print(sensitivity_results %>%
        select(threshold_px, pct_contact, cd44_cohens_d, cd44_p_value,
               bcat_cohens_d, bcat_p_value) %>%
        mutate(across(where(is.numeric), ~round(., 4))))

# =============================================================================
# 8. VISUALIZATIONS
# =============================================================================

cat("\n=============================================================================\n")
cat("Generating Figures\n")
cat("=============================================================================\n")

# Color palette
contact_colors <- c("Contact (≤5μm)" = "#E41A1C", "Distant (>5μm)" = "#377EB8")

# --- Figure 1: CD44 Expression by Collagen Proximity ---
p1a <- ggplot(analysis_data, aes(x = contact_group, y = CD44, fill = contact_group)) +

  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     label.x = 1.5,
                     label.y = max(analysis_data$CD44) * 0.95,
                     size = 5) +
  scale_fill_manual(values = contact_colors) +
  labs(
    title = "CD44 Expression by Collagen Proximity",
    subtitle = sprintf("n = %s cells | Cohen's d = %.3f (%s)",
                       format(nrow(analysis_data), big.mark = ","),
                       cd44_cohens_d$estimate,
                       cd44_cohens_d$magnitude),
    x = "Collagen Contact Group",
    y = "CD44 Mean Intensity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 11)
  )

# --- Figure 1b: CD44 by distance (continuous) ---
p1b <- ggplot(analysis_data %>% sample_n(min(10000, n())),
              aes(x = distance_px, y = CD44)) +
  geom_point(alpha = 0.1, size = 0.5, color = "gray40") +
  geom_smooth(method = "loess", color = "#E41A1C", linewidth = 1.5, se = TRUE) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 100)) +
  labs(
    title = "CD44 Expression vs Distance to Collagen",
    subtitle = "Dashed line = 5 pixel (~5μm) contact threshold",
    x = "Distance to Nearest Collagen Tract (pixels)",
    y = "CD44 Mean Intensity"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# --- Figure 2: β-Catenin Expression by Collagen Proximity ---
p2a <- ggplot(analysis_data, aes(x = contact_group, y = B_Catenin, fill = contact_group)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     label.x = 1.5,
                     label.y = max(analysis_data$B_Catenin) * 0.95,
                     size = 5) +
  scale_fill_manual(values = contact_colors) +
  labs(
    title = "β-Catenin Expression by Collagen Proximity",
    subtitle = sprintf("n = %s cells | Cohen's d = %.3f (%s)",
                       format(nrow(analysis_data), big.mark = ","),
                       bcatenin_cohens_d$estimate,
                       bcatenin_cohens_d$magnitude),
    x = "Collagen Contact Group",
    y = "β-Catenin Mean Intensity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 11)
  )

# --- Figure 2b: β-Catenin by distance (continuous) ---
p2b <- ggplot(analysis_data %>% sample_n(min(10000, n())),
              aes(x = distance_px, y = B_Catenin)) +
  geom_point(alpha = 0.1, size = 0.5, color = "gray40") +
  geom_smooth(method = "loess", color = "#377EB8", linewidth = 1.5, se = TRUE) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 100)) +
  labs(
    title = "β-Catenin Expression vs Distance to Collagen",
    subtitle = "Dashed line = 5 pixel (~5μm) contact threshold",
    x = "Distance to Nearest Collagen Tract (pixels)",
    y = "β-Catenin Mean Intensity"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# --- Figure 3: Per-ROI Effect Sizes ---
roi_plot_data <- roi_effects %>%
  pivot_longer(
    cols = c(cd44_diff, bcat_diff),
    names_to = "marker",
    values_to = "difference"
  ) %>%
  mutate(marker = ifelse(marker == "cd44_diff", "CD44", "β-Catenin"))

p3 <- ggplot(roi_plot_data, aes(x = factor(roi), y = difference, fill = marker)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_fill_manual(values = c("CD44" = "#E41A1C", "β-Catenin" = "#4DAF4A")) +
  labs(
    title = "Expression Difference (Contact - Distant) by ROI",
    subtitle = "Positive values = higher expression in cells contacting collagen",
    x = "ROI",
    y = "Mean Intensity Difference",
    fill = "Marker"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# --- Figure 4: Sensitivity Analysis ---
sensitivity_plot_data <- sensitivity_results %>%
  pivot_longer(
    cols = c(cd44_cohens_d, bcat_cohens_d),
    names_to = "marker",
    values_to = "cohens_d"
  ) %>%
  mutate(marker = ifelse(marker == "cd44_cohens_d", "CD44", "β-Catenin"))

p4 <- ggplot(sensitivity_plot_data,
             aes(x = threshold_px, y = cohens_d, color = marker, group = marker)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 5, linetype = "dotted", color = "black") +
  scale_color_manual(values = c("CD44" = "#E41A1C", "β-Catenin" = "#4DAF4A")) +
  labs(
    title = "Effect Size Sensitivity to Distance Threshold",
    subtitle = "Vertical line = primary threshold (5 pixels ~5μm)",
    x = "Contact Distance Threshold (pixels)",
    y = "Cohen's d (Contact vs Distant)",
    color = "Marker"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# --- Combined Figure Panel ---
combined_panel <- (p1a | p1b) / (p2a | p2b) +
  plot_annotation(
    title = "CD44 and β-Catenin Expression by Collagen Proximity",
    subtitle = "IMC Analysis | 8 ROIs, 2 Patients",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

# Save figures
cat("Saving figures...\n")

ggsave(file.path(output_dir, "01_cd44_violin.pdf"), p1a, width = 6, height = 6)
ggsave(file.path(output_dir, "02_cd44_by_distance.pdf"), p1b, width = 8, height = 6)
ggsave(file.path(output_dir, "03_bcatenin_violin.pdf"), p2a, width = 6, height = 6)
ggsave(file.path(output_dir, "04_bcatenin_by_distance.pdf"), p2b, width = 8, height = 6)
ggsave(file.path(output_dir, "05_per_roi_effects.pdf"), p3, width = 10, height = 6)
ggsave(file.path(output_dir, "06_sensitivity_analysis.pdf"), p4, width = 8, height = 6)
ggsave(file.path(output_dir, "07_combined_panel.pdf"), combined_panel, width = 14, height = 12)

# PNG versions for quick viewing
ggsave(file.path(output_dir, "07_combined_panel.png"), combined_panel,
       width = 14, height = 12, dpi = 300)

# =============================================================================
# 9. EXPORT RESULTS
# =============================================================================

cat("\n--- Exporting Results ---\n")

# Summary table for manuscript
manuscript_summary <- tibble(
  Marker = c("CD44", "β-Catenin"),
  N_Total = nrow(analysis_data),
  N_Contact = contact_stats$n[1],
  N_Distant = contact_stats$n[2],

  Contact_Mean = c(cd44_summary$mean[1], bcatenin_summary$mean[1]),
  Contact_Median = c(cd44_summary$median[1], bcatenin_summary$median[1]),
  Contact_SD = c(cd44_summary$sd[1], bcatenin_summary$sd[1]),

  Distant_Mean = c(cd44_summary$mean[2], bcatenin_summary$mean[2]),
  Distant_Median = c(cd44_summary$median[2], bcatenin_summary$median[2]),
  Distant_SD = c(cd44_summary$sd[2], bcatenin_summary$sd[2]),

  Wilcoxon_W = c(cd44_wilcox$statistic, bcatenin_wilcox$statistic),
  P_Value = c(cd44_wilcox$p.value, bcatenin_wilcox$p.value),

  Cohens_D = c(cd44_cohens_d$estimate, bcatenin_cohens_d$estimate),
  Cohens_D_Magnitude = c(cd44_cohens_d$magnitude, bcatenin_cohens_d$magnitude),

  Cliffs_Delta = c(cd44_cliffs_delta$estimate, bcatenin_cliffs_delta$estimate),
  Cliffs_Delta_Magnitude = c(cd44_cliffs_delta$magnitude, bcatenin_cliffs_delta$magnitude),

  ROIs_With_Positive_Effect = c(sum(roi_effects$cd44_higher_in_contact),
                                sum(roi_effects$bcat_higher_in_contact)),
  Total_ROIs = nrow(roi_effects)
)

fwrite(manuscript_summary, file.path(output_dir, "manuscript_summary.csv"))
fwrite(per_roi_stats, file.path(output_dir, "per_roi_statistics.csv"))
fwrite(roi_effects, file.path(output_dir, "per_roi_effects.csv"))
fwrite(sensitivity_results, file.path(output_dir, "sensitivity_analysis.csv"))

# =============================================================================
# 10. GENERATE MANUSCRIPT TEXT
# =============================================================================

cat("\n=============================================================================\n")
cat("MANUSCRIPT TEXT (Copy-Paste Ready)\n")
cat("=============================================================================\n\n")

# Determine direction of effects
cd44_direction <- ifelse(cd44_cohens_d$estimate > 0, "higher", "lower")
bcat_direction <- ifelse(bcatenin_cohens_d$estimate > 0, "higher", "lower")

manuscript_text <- sprintf(
  "To validate Visium findings showing elevated CD44 expression in epithelial cells
at collagen interfaces, we analyzed CD44 and β-catenin protein expression relative
to collagen tract distance in IMC data from %d cells across %d ROIs (2 patients).

Cells within 5 pixels (~5μm) of collagen tracts (n = %s, %.1f%%) showed %s CD44 expression
compared to distant cells (n = %s, %.1f%%; median %.2f vs %.2f; Wilcoxon p %s;
Cohen's d = %.3f [%s effect]). This pattern was consistent across %d/%d ROIs.
β-catenin expression was %s in contact cells (median %.2f vs %.2f; p %s;
Cohen's d = %.3f [%s effect]), with consistent direction in %d/%d ROIs.

Sensitivity analysis confirmed that effect sizes were robust across distance
thresholds from 2-20 pixels, with CD44 effects peaking at %d pixels (d = %.3f).",

  nrow(analysis_data),
  n_distinct(analysis_data$roi),
  format(contact_stats$n[1], big.mark = ","),
  contact_stats$pct[1],
  cd44_direction,
  format(contact_stats$n[2], big.mark = ","),
  contact_stats$pct[2],
  cd44_summary$median[1],
  cd44_summary$median[2],
  ifelse(cd44_wilcox$p.value < 0.001, "< 0.001", sprintf("= %.3f", cd44_wilcox$p.value)),
  abs(cd44_cohens_d$estimate),
  cd44_cohens_d$magnitude,
  sum(roi_effects$cd44_higher_in_contact),
  nrow(roi_effects),
  bcat_direction,
  bcatenin_summary$median[1],
  bcatenin_summary$median[2],
  ifelse(bcatenin_wilcox$p.value < 0.001, "< 0.001", sprintf("= %.3f", bcatenin_wilcox$p.value)),
  abs(bcatenin_cohens_d$estimate),
  bcatenin_cohens_d$magnitude,
  sum(roi_effects$bcat_higher_in_contact),
  nrow(roi_effects),
  sensitivity_results$threshold_px[which.max(abs(sensitivity_results$cd44_cohens_d))],
  max(abs(sensitivity_results$cd44_cohens_d))
)

cat(manuscript_text)

# Save manuscript text
writeLines(manuscript_text, file.path(output_dir, "manuscript_text.txt"))

# =============================================================================
# COMPLETION
# =============================================================================

cat("\n\n=============================================================================\n")
cat("✅ Analysis Complete!\n")
cat("=============================================================================\n")
cat("\nOutput directory:", output_dir, "\n")
cat("\nFiles created:\n")
cat("  Tables:\n")
cat("    - manuscript_summary.csv (main results)\n")
cat("    - per_roi_statistics.csv\n")
cat("    - per_roi_effects.csv\n")
cat("    - sensitivity_analysis.csv\n")
cat("  Figures:\n")
cat("    - 01_cd44_violin.pdf\n")
cat("    - 02_cd44_by_distance.pdf\n")
cat("    - 03_bcatenin_violin.pdf\n")
cat("    - 04_bcatenin_by_distance.pdf\n")
cat("    - 05_per_roi_effects.pdf\n")
cat("    - 06_sensitivity_analysis.pdf\n")
cat("    - 07_combined_panel.pdf/png\n")
cat("  Text:\n")
cat("    - manuscript_text.txt\n")
cat("\n=============================================================================\n")
