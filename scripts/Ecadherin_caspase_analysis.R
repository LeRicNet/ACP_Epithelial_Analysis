# =============================================================================
# E-Cadherin and Cleaved Caspase 3 Spatial Analysis
# =============================================================================
# Purpose: Complement Ki-67 proliferation findings by analyzing:
#   1. E-Cadherin expression relative to collagen (epithelial integrity/EMT)
#   2. Cleaved Caspase 3 relative to collagen (apoptosis/stress)
#
# Hypotheses:
#   - E-Cadherin may be REDUCED at collagen boundaries (partial EMT)
#   - Cleaved Caspase 3 may be ELEVATED at boundaries (stress zone)
#   - Together with Ki-67, defines functional compartmentalization
#
# Input:   Filtered cell data with marker expression
#          Cell-tract distances (from collagen QC pipeline)
# Output:  Statistical comparisons, combined multi-marker visualization
# =============================================================================

library(tidyverse)
library(data.table)
library(RANN)
library(effsize)
library(ggpubr)
library(patchwork)
library(scales)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths
data_dir <- "/home/rstudio/interfaces/ACP_IMC/results"
output_dir <- file.path(data_dir, "ecadherin_caspase_analysis")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Marker column patterns
ecadherin_pattern <- "E-Cadherin|E-cadherin|Ecadherin|Cadherin|164Dy|Dy164|160Gd|Gd160"
caspase_pattern <- "Caspase|Cleaved|CC3|142Nd|Nd142"
ki67_pattern <- "Ki-67|Ki67|173Yb|Yb173"
panck_pattern <- "Pan-Cytokeratin|Cytokeratin|174Yb|Yb174"

# Thresholds
panck_percentile <- 75    # Top 25% = epithelial cells

# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("=============================================================================\n")
cat("E-Cadherin and Cleaved Caspase 3 Spatial Analysis\n")
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

# Print all column names to help identify markers
cat("\n  All marker columns:\n")
marker_cols <- colnames(marker_data)[!colnames(marker_data) %in% c("ObjectNumber", "AcquisitionId", "ImageNumber", "CentroidX", "CentroidY", "CellId")]
cat("   ", paste(marker_cols, collapse = ", "), "\n")

# --- Merge datasets ---
cat("\n--- Merging datasets ---\n")

tract_distances <- tract_distances %>%
  mutate(
    ObjectNumber = as.integer(sub(".*_([0-9]+)$", "\\1", cell_id)),
    ImageNumber = roi_id
  )

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

ecadherin_col <- find_column(cell_data, ecadherin_pattern, "E-Cadherin")
caspase_col <- find_column(cell_data, caspase_pattern, "Cleaved Caspase 3")
ki67_col <- find_column(cell_data, ki67_pattern, "Ki-67")
panck_col <- find_column(cell_data, panck_pattern, "Pan-Cytokeratin")

# Check which markers were found
markers_found <- c(
  "E-Cadherin" = !is.null(ecadherin_col),
  "Caspase" = !is.null(caspase_col),
  "Ki-67" = !is.null(ki67_col),
  "Pan-CK" = !is.null(panck_col)
)

cat("\n  Markers found:", paste(names(markers_found)[markers_found], collapse = ", "), "\n")
if (any(!markers_found)) {
  cat("  Markers NOT found:", paste(names(markers_found)[!markers_found], collapse = ", "), "\n")
}

# Stop if critical markers missing
if (is.null(panck_col)) stop("Pan-Cytokeratin column not found - cannot identify epithelial cells!")

# =============================================================================
# 3. CLASSIFY CELLS AND PREPARE DATA
# =============================================================================

cat("\n--- Preparing Analysis Data ---\n")

# Calculate epithelial threshold
epithelial_threshold <- quantile(cell_data[[panck_col]], panck_percentile/100, na.rm = TRUE)

# Classify cells
cell_data <- cell_data %>%
  mutate(
    is_epithelial = .data[[panck_col]] >= epithelial_threshold,
    collagen_contact = ifelse(contact_5px, "Contact (<=5um)", "Distant (>5um)"),
    distance_bin = cut(distance_to_tract,
                       breaks = c(0, 5, 10, 20, 50, Inf),
                       labels = c("0-5", "5-10", "10-20", "20-50", ">50"),
                       include.lowest = TRUE)
  )

# Add marker expressions if found
if (!is.null(ecadherin_col)) {
  cell_data$ECadherin_expr <- cell_data[[ecadherin_col]]
}
if (!is.null(caspase_col)) {
  cell_data$Caspase_expr <- cell_data[[caspase_col]]
}
if (!is.null(ki67_col)) {
  cell_data$Ki67_expr <- cell_data[[ki67_col]]
}

cat(sprintf("  Total cells: %s\n", format(nrow(cell_data), big.mark = ",")))
cat(sprintf("  Epithelial cells: %s (%.1f%%)\n",
            format(sum(cell_data$is_epithelial), big.mark = ","),
            mean(cell_data$is_epithelial) * 100))

# =============================================================================
# 4. E-CADHERIN ANALYSIS (within epithelial cells)
# =============================================================================

if (!is.null(ecadherin_col)) {

  cat("\n=============================================================================\n")
  cat("E-Cadherin Analysis (Epithelial Cells Only)\n")
  cat("=============================================================================\n")

  # Subset to epithelial cells
  epithelial_data <- cell_data %>% filter(is_epithelial)

  # --- E-Cadherin by collagen contact ---
  ecad_by_contact <- epithelial_data %>%
    group_by(collagen_contact) %>%
    summarise(
      n = n(),
      mean = mean(ECadherin_expr, na.rm = TRUE),
      median = median(ECadherin_expr, na.rm = TRUE),
      sd = sd(ECadherin_expr, na.rm = TRUE),
      .groups = "drop"
    )

  cat("\n--- E-Cadherin by Collagen Contact (Epithelial Cells) ---\n")
  print(ecad_by_contact)

  # Wilcoxon test
  ecad_contact <- epithelial_data %>% filter(contact_5px) %>% pull(ECadherin_expr)
  ecad_distant <- epithelial_data %>% filter(!contact_5px) %>% pull(ECadherin_expr)

  ecad_wilcox <- wilcox.test(ecad_contact, ecad_distant)
  ecad_d <- cohen.d(ecad_contact, ecad_distant)

  cat(sprintf("\n  Wilcoxon p-value: %s\n", format.pval(ecad_wilcox$p.value, digits = 3)))
  cat(sprintf("  Cohen's d: %.3f (%s)\n", ecad_d$estimate, ecad_d$magnitude))
  cat(sprintf("  Direction: E-Cadherin is %s in contact cells\n",
              ifelse(ecad_d$estimate > 0, "HIGHER", "LOWER")))

  # Spearman correlation with distance
  ecad_cor <- cor.test(epithelial_data$ECadherin_expr, epithelial_data$distance_to_tract,
                       method = "spearman")
  cat(sprintf("  Spearman rho vs distance: %.4f (p %s)\n",
              ecad_cor$estimate,
              ifelse(ecad_cor$p.value < 0.001, "< 0.001", sprintf("= %.3f", ecad_cor$p.value))))

  # --- E-Cadherin by distance bins ---
  ecad_by_distance <- epithelial_data %>%
    group_by(distance_bin) %>%
    summarise(
      n = n(),
      mean = mean(ECadherin_expr, na.rm = TRUE),
      se = sd(ECadherin_expr, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  cat("\n--- E-Cadherin by Distance Bin ---\n")
  print(ecad_by_distance)

  # Store results
  ecad_results <- list(
    by_contact = ecad_by_contact,
    wilcox = ecad_wilcox,
    cohens_d = ecad_d,
    correlation = ecad_cor,
    by_distance = ecad_by_distance
  )

} else {
  cat("\n  E-Cadherin analysis skipped (marker not found)\n")
  ecad_results <- NULL
}

# =============================================================================
# 5. CLEAVED CASPASE 3 ANALYSIS (all cells)
# =============================================================================

if (!is.null(caspase_col)) {

  cat("\n=============================================================================\n")
  cat("Cleaved Caspase 3 (Apoptosis) Analysis\n")
  cat("=============================================================================\n")

  # --- Caspase by collagen contact (all cells) ---
  casp_by_contact_all <- cell_data %>%
    group_by(collagen_contact) %>%
    summarise(
      n = n(),
      mean = mean(Caspase_expr, na.rm = TRUE),
      median = median(Caspase_expr, na.rm = TRUE),
      sd = sd(Caspase_expr, na.rm = TRUE),
      .groups = "drop"
    )

  cat("\n--- Cleaved Caspase 3 by Collagen Contact (All Cells) ---\n")
  print(casp_by_contact_all)

  # Wilcoxon test (all cells)
  casp_contact <- cell_data %>% filter(contact_5px) %>% pull(Caspase_expr)
  casp_distant <- cell_data %>% filter(!contact_5px) %>% pull(Caspase_expr)

  casp_wilcox <- wilcox.test(casp_contact, casp_distant)
  casp_d <- cohen.d(casp_contact, casp_distant)

  cat(sprintf("\n  Wilcoxon p-value: %s\n", format.pval(casp_wilcox$p.value, digits = 3)))
  cat(sprintf("  Cohen's d: %.3f (%s)\n", casp_d$estimate, casp_d$magnitude))
  cat(sprintf("  Direction: Caspase 3 is %s in contact cells\n",
              ifelse(casp_d$estimate > 0, "HIGHER", "LOWER")))

  # --- Caspase by collagen contact (epithelial only) ---
  epithelial_data <- cell_data %>% filter(is_epithelial)

  casp_by_contact_epi <- epithelial_data %>%
    group_by(collagen_contact) %>%
    summarise(
      n = n(),
      mean = mean(Caspase_expr, na.rm = TRUE),
      median = median(Caspase_expr, na.rm = TRUE),
      .groups = "drop"
    )

  cat("\n--- Cleaved Caspase 3 by Collagen Contact (Epithelial Only) ---\n")
  print(casp_by_contact_epi)

  casp_contact_epi <- epithelial_data %>% filter(contact_5px) %>% pull(Caspase_expr)
  casp_distant_epi <- epithelial_data %>% filter(!contact_5px) %>% pull(Caspase_expr)

  casp_wilcox_epi <- wilcox.test(casp_contact_epi, casp_distant_epi)
  casp_d_epi <- cohen.d(casp_contact_epi, casp_distant_epi)

  cat(sprintf("\n  Epithelial cells only:\n"))
  cat(sprintf("    Wilcoxon p-value: %s\n", format.pval(casp_wilcox_epi$p.value, digits = 3)))
  cat(sprintf("    Cohen's d: %.3f (%s)\n", casp_d_epi$estimate, casp_d_epi$magnitude))

  # Spearman correlation
  casp_cor <- cor.test(cell_data$Caspase_expr, cell_data$distance_to_tract, method = "spearman")
  cat(sprintf("  Spearman rho vs distance (all cells): %.4f\n", casp_cor$estimate))

  # --- Caspase by distance bins ---
  casp_by_distance <- cell_data %>%
    group_by(distance_bin) %>%
    summarise(
      n = n(),
      mean = mean(Caspase_expr, na.rm = TRUE),
      se = sd(Caspase_expr, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  cat("\n--- Cleaved Caspase 3 by Distance Bin ---\n")
  print(casp_by_distance)

  # Store results
  casp_results <- list(
    by_contact_all = casp_by_contact_all,
    by_contact_epi = casp_by_contact_epi,
    wilcox_all = casp_wilcox,
    wilcox_epi = casp_wilcox_epi,
    cohens_d_all = casp_d,
    cohens_d_epi = casp_d_epi,
    correlation = casp_cor,
    by_distance = casp_by_distance
  )

} else {
  cat("\n  Cleaved Caspase 3 analysis skipped (marker not found)\n")
  casp_results <- NULL
}

# =============================================================================
# 6. COMBINED MULTI-MARKER ANALYSIS
# =============================================================================

cat("\n=============================================================================\n")
cat("Combined Multi-Marker Analysis\n")
cat("=============================================================================\n")

# Build summary table of all markers vs collagen contact
marker_summary <- tibble(
  Marker = character(),
  Population = character(),
  Contact_mean = numeric(),
  Distant_mean = numeric(),
  Cohens_d = numeric(),
  Direction = character(),
  p_value = numeric()
)

# Add Ki-67 (from previous analysis, re-calculate here)
if (!is.null(ki67_col)) {
  ki67_contact <- cell_data %>% filter(contact_5px) %>% pull(Ki67_expr)
  ki67_distant <- cell_data %>% filter(!contact_5px) %>% pull(Ki67_expr)
  ki67_d <- cohen.d(ki67_contact, ki67_distant)
  ki67_wilcox <- wilcox.test(ki67_contact, ki67_distant)

  marker_summary <- marker_summary %>%
    add_row(
      Marker = "Ki-67",
      Population = "All cells",
      Contact_mean = mean(ki67_contact, na.rm = TRUE),
      Distant_mean = mean(ki67_distant, na.rm = TRUE),
      Cohens_d = ki67_d$estimate,
      Direction = ifelse(ki67_d$estimate > 0, "Higher at contact", "Lower at contact"),
      p_value = ki67_wilcox$p.value
    )
}

# Add E-Cadherin
if (!is.null(ecad_results)) {
  marker_summary <- marker_summary %>%
    add_row(
      Marker = "E-Cadherin",
      Population = "Epithelial cells",
      Contact_mean = ecad_results$by_contact$mean[ecad_results$by_contact$collagen_contact == "Contact (<=5um)"],
      Distant_mean = ecad_results$by_contact$mean[ecad_results$by_contact$collagen_contact == "Distant (>5um)"],
      Cohens_d = ecad_results$cohens_d$estimate,
      Direction = ifelse(ecad_results$cohens_d$estimate > 0, "Higher at contact", "Lower at contact"),
      p_value = ecad_results$wilcox$p.value
    )
}

# Add Caspase
if (!is.null(casp_results)) {
  marker_summary <- marker_summary %>%
    add_row(
      Marker = "Cleaved Caspase 3",
      Population = "All cells",
      Contact_mean = casp_results$by_contact_all$mean[casp_results$by_contact_all$collagen_contact == "Contact (<=5um)"],
      Distant_mean = casp_results$by_contact_all$mean[casp_results$by_contact_all$collagen_contact == "Distant (>5um)"],
      Cohens_d = casp_results$cohens_d_all$estimate,
      Direction = ifelse(casp_results$cohens_d_all$estimate > 0, "Higher at contact", "Lower at contact"),
      p_value = casp_results$wilcox_all$p.value
    )
}

cat("\n--- Summary: All Markers vs Collagen Contact ---\n")
print(marker_summary)

# =============================================================================
# 7. VISUALIZATIONS
# =============================================================================

cat("\n=============================================================================\n")
cat("Generating Figures\n")
cat("=============================================================================\n")

# Color scheme
contact_colors <- c("Contact (<=5um)" = "#E41A1C", "Distant (>5um)" = "#377EB8")

# Initialize plot list
plot_list <- list()

# --- E-Cadherin plots ---
if (!is.null(ecadherin_col)) {

  epithelial_data <- cell_data %>% filter(is_epithelial)

  # Violin plot
  p_ecad_violin <- ggplot(epithelial_data, aes(x = collagen_contact, y = ECadherin_expr,
                                               fill = collagen_contact)) +
    geom_violin(alpha = 0.7, trim = TRUE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
    scale_fill_manual(values = contact_colors) +
    scale_y_continuous(limits = c(0, quantile(epithelial_data$ECadherin_expr, 0.99, na.rm = TRUE))) +
    labs(
      title = "E-Cadherin by Collagen Contact",
      subtitle = sprintf("Epithelial cells only | Cohen's d = %.3f (%s)",
                         ecad_results$cohens_d$estimate, ecad_results$cohens_d$magnitude),
      x = "",
      y = "E-Cadherin Intensity"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

  # Distance bins plot
  p_ecad_distance <- ggplot(ecad_results$by_distance, aes(x = distance_bin, y = mean)) +
    geom_col(fill = "#4DAF4A", width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    labs(
      title = "E-Cadherin by Distance to Collagen",
      subtitle = sprintf("Spearman rho = %.3f", ecad_results$correlation$estimate),
      x = "Distance to Collagen (pixels)",
      y = "Mean E-Cadherin Intensity"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  plot_list$ecad_violin <- p_ecad_violin
  plot_list$ecad_distance <- p_ecad_distance
}

# --- Caspase plots ---
if (!is.null(caspase_col)) {

  # Violin plot (all cells)
  p_casp_violin <- ggplot(cell_data, aes(x = collagen_contact, y = Caspase_expr,
                                         fill = collagen_contact)) +
    geom_violin(alpha = 0.7, trim = TRUE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
    scale_fill_manual(values = contact_colors) +
    scale_y_continuous(limits = c(0, quantile(cell_data$Caspase_expr, 0.99, na.rm = TRUE))) +
    labs(
      title = "Cleaved Caspase 3 by Collagen Contact",
      subtitle = sprintf("All cells | Cohen's d = %.3f (%s)",
                         casp_results$cohens_d_all$estimate, casp_results$cohens_d_all$magnitude),
      x = "",
      y = "Cleaved Caspase 3 Intensity"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

  # Distance bins plot
  p_casp_distance <- ggplot(casp_results$by_distance, aes(x = distance_bin, y = mean)) +
    geom_col(fill = "#984EA3", width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    labs(
      title = "Cleaved Caspase 3 by Distance to Collagen",
      subtitle = sprintf("Spearman rho = %.3f", casp_results$correlation$estimate),
      x = "Distance to Collagen (pixels)",
      y = "Mean Caspase 3 Intensity"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  plot_list$casp_violin <- p_casp_violin
  plot_list$casp_distance <- p_casp_distance
}

# --- Combined effect size plot ---
if (nrow(marker_summary) > 0) {

  p_effects <- ggplot(marker_summary, aes(x = reorder(Marker, Cohens_d), y = Cohens_d,
                                          fill = Cohens_d > 0)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 0, linetype = "solid") +
    geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", alpha = 0.5) +
    geom_text(aes(label = sprintf("d=%.3f", Cohens_d),
                  y = ifelse(Cohens_d >= 0, Cohens_d + 0.02, Cohens_d - 0.02)),
              size = 3.5) +
    scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                      labels = c("Lower at contact", "Higher at contact"),
                      name = "") +
    labs(
      title = "Effect Sizes: Marker Expression vs Collagen Contact",
      subtitle = "Dashed lines = |d| = 0.2 (small effect threshold)",
      x = "",
      y = "Cohen's d (Contact - Distant)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom") +
    coord_flip()

  plot_list$effects <- p_effects
}

# --- Ki-67 violin for comparison (if available) ---
if (!is.null(ki67_col)) {

  p_ki67_violin <- ggplot(cell_data, aes(x = collagen_contact, y = Ki67_expr,
                                         fill = collagen_contact)) +
    geom_violin(alpha = 0.7, trim = TRUE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
    scale_fill_manual(values = contact_colors) +
    scale_y_continuous(limits = c(0, quantile(cell_data$Ki67_expr, 0.99, na.rm = TRUE))) +
    labs(
      title = "Ki-67 by Collagen Contact",
      subtitle = sprintf("All cells | Cohen's d = %.3f", ki67_d$estimate),
      x = "",
      y = "Ki-67 Intensity"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

  plot_list$ki67_violin <- p_ki67_violin
}

# --- Build combined panel ---
# Determine layout based on available plots
if (length(plot_list) >= 5) {
  combined_panel <- (plot_list$ki67_violin | plot_list$ecad_violin | plot_list$casp_violin) /
    (plot_list$ecad_distance | plot_list$casp_distance | plot_list$effects) +
    plot_annotation(
      title = "Multi-Marker Spatial Analysis: Functional Compartmentalization",
      subtitle = sprintf("IMC Analysis | %s cells across %d ROIs",
                         format(nrow(cell_data), big.mark = ","),
                         n_distinct(cell_data$ImageNumber)),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
    )
} else if (length(plot_list) >= 2) {
  combined_panel <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title = "Marker Expression vs Collagen Proximity",
      tag_levels = "A"
    )
} else {
  combined_panel <- NULL
}

# Save figures
cat("Saving figures...\n")

if (!is.null(ecadherin_col)) {
  ggsave(file.path(output_dir, "01_ecadherin_violin.pdf"), plot_list$ecad_violin, width = 6, height = 6)
  ggsave(file.path(output_dir, "02_ecadherin_by_distance.pdf"), plot_list$ecad_distance, width = 8, height = 6)
}

if (!is.null(caspase_col)) {
  ggsave(file.path(output_dir, "03_caspase_violin.pdf"), plot_list$casp_violin, width = 6, height = 6)
  ggsave(file.path(output_dir, "04_caspase_by_distance.pdf"), plot_list$casp_distance, width = 8, height = 6)
}

if (!is.null(combined_panel)) {
  ggsave(file.path(output_dir, "05_combined_panel.pdf"), combined_panel, width = 16, height = 10)
  ggsave(file.path(output_dir, "05_combined_panel.png"), combined_panel, width = 16, height = 10, dpi = 300)
}

ggsave(file.path(output_dir, "06_effect_sizes.pdf"), plot_list$effects, width = 8, height = 6)

# =============================================================================
# 8. EXPORT RESULTS
# =============================================================================

cat("\n--- Exporting Results ---\n")

fwrite(marker_summary, file.path(output_dir, "marker_summary.csv"))

if (!is.null(ecad_results)) {
  fwrite(ecad_results$by_contact, file.path(output_dir, "ecadherin_by_contact.csv"))
  fwrite(ecad_results$by_distance, file.path(output_dir, "ecadherin_by_distance.csv"))
}

if (!is.null(casp_results)) {
  fwrite(casp_results$by_contact_all, file.path(output_dir, "caspase_by_contact.csv"))
  fwrite(casp_results$by_distance, file.path(output_dir, "caspase_by_distance.csv"))
}

# =============================================================================
# 9. GENERATE MANUSCRIPT TEXT
# =============================================================================

cat("\n=============================================================================\n")
cat("MANUSCRIPT TEXT (Copy-Paste Ready)\n")
cat("=============================================================================\n\n")

# Build text based on available markers
text_parts <- c()

text_parts <- c(text_parts, sprintf(
  "To further characterize functional compartmentalization within the tumor
architecture, we analyzed additional markers relative to collagen proximity
in epithelial cells (n = %s cells across %d ROIs).",
  format(nrow(cell_data), big.mark = ","),
  n_distinct(cell_data$ImageNumber)
))

if (!is.null(ecad_results)) {
  ecad_direction <- ifelse(ecad_results$cohens_d$estimate > 0, "elevated", "reduced")
  text_parts <- c(text_parts, sprintf(
    "\n\nE-Cadherin, a marker of epithelial integrity, was %s in epithelial cells
contacting collagen compared to distant cells (mean %.2f vs %.2f; Cohen's d = %.3f,
%s; p %s). This %s suggests that epithelial cells at collagen boundaries %s
epithelial adhesion properties, %s with %s.",
    ecad_direction,
    ecad_results$by_contact$mean[1],
    ecad_results$by_contact$mean[2],
    ecad_results$cohens_d$estimate,
    ecad_results$cohens_d$magnitude,
    ifelse(ecad_results$wilcox$p.value < 0.001, "< 0.001", sprintf("= %.3f", ecad_results$wilcox$p.value)),
    ifelse(abs(ecad_results$cohens_d$estimate) < 0.2, "negligible difference", "pattern"),
    ifelse(ecad_results$cohens_d$estimate > 0, "maintain", "partially lose"),
    ifelse(ecad_results$cohens_d$estimate > 0, "inconsistent", "consistent"),
    ifelse(ecad_results$cohens_d$estimate > 0, "stable epithelial identity at boundaries",
           "partial epithelial-mesenchymal transition at the stromal interface")
  ))
}

if (!is.null(casp_results)) {
  casp_direction <- ifelse(casp_results$cohens_d_all$estimate > 0, "elevated", "reduced")
  text_parts <- c(text_parts, sprintf(
    "\n\nCleaved Caspase 3, marking apoptotic cells, was %s at collagen boundaries
(mean %.2f vs %.2f; Cohen's d = %.3f, %s; p %s). %s",
    casp_direction,
    casp_results$by_contact_all$mean[1],
    casp_results$by_contact_all$mean[2],
    casp_results$cohens_d_all$estimate,
    casp_results$cohens_d_all$magnitude,
    ifelse(casp_results$wilcox_all$p.value < 0.001, "< 0.001", sprintf("= %.3f", casp_results$wilcox_all$p.value)),
    ifelse(abs(casp_results$cohens_d_all$estimate) < 0.2,
           "This negligible effect indicates that collagen proximity does not substantially influence cell survival.",
           ifelse(casp_results$cohens_d_all$estimate > 0,
                  "This elevation suggests that collagen-adjacent regions experience increased cellular stress.",
                  "This reduction suggests that collagen-adjacent regions may be protected from apoptosis."))
  ))
}

text_parts <- c(text_parts, sprintf(
  "\n\nTogether with the Ki-67 findings (reduced proliferation at collagen boundaries,
OR = 0.85), these results support a model of functional compartmentalization:
epithelial cores serve as proliferative niches while boundary regions engage
in ECM interaction rather than active division."))

manuscript_text <- paste(text_parts, collapse = "")

cat(manuscript_text)

writeLines(manuscript_text, file.path(output_dir, "manuscript_text.txt"))

# =============================================================================
# COMPLETION
# =============================================================================

cat("\n\n=============================================================================\n")
cat("Analysis Complete!\n")
cat("=============================================================================\n")
cat("\nOutput directory:", output_dir, "\n")
cat("\nKey findings:\n")

if (!is.null(ecad_results)) {
  cat(sprintf("  E-Cadherin: d = %.3f (%s) - %s at collagen contact\n",
              ecad_results$cohens_d$estimate,
              ecad_results$cohens_d$magnitude,
              ifelse(ecad_results$cohens_d$estimate > 0, "HIGHER", "LOWER")))
}

if (!is.null(casp_results)) {
  cat(sprintf("  Cleaved Caspase 3: d = %.3f (%s) - %s at collagen contact\n",
              casp_results$cohens_d_all$estimate,
              casp_results$cohens_d_all$magnitude,
              ifelse(casp_results$cohens_d_all$estimate > 0, "HIGHER", "LOWER")))
}

cat("\n=============================================================================\n")
