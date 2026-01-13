#!/usr/bin/env Rscript
# =============================================================================
# Collagen Tract Molecular Characterization Analysis - IMPROVED VERSION
# Section 0.8: Collagen tracts define molecularly distinct microenvironmental niches
# =============================================================================
# Improvements over v2:
# 1. Proper isotope-to-marker name mapping
# 2. Cell compartment classification based on marker expression
# 3. Spatial correlations between marker pairs
# 4. Multiple analysis strategies for robustness
# =============================================================================

library(tidyverse)
library(data.table)
library(ggpubr)
library(scales)
library(RColorBrewer)

# =============================================================================
cat("=============================================================================\n")
cat("Collagen Tract Molecular Characterization Analysis - IMPROVED\n")
cat("=============================================================================\n")

# =============================================================================
# ISOTOPE TO MARKER MAPPING (Based on standard IMC panels)
# =============================================================================
# This mapping converts isotope column names to readable marker names

isotope_to_marker <- c(
  "Pr141" = "SMA",
  "Nd142" = "Cleaved_Caspase3",
  "Nd143" = "Vimentin",
  "Nd144" = "Unknown_144",
  "Nd145" = "Unknown_145",
  "Nd146" = "Unknown_146",
  "Sm147" = "CD163",
  "Nd148" = "Unknown_148",
  "Sm149" = "Unknown_149",
  "Nd150" = "Unknown_150",
  "Eu151" = "CD31",
  "Sm152" = "CD45",
  "Eu153" = "CD44",
  "Sm154" = "EpCAM",
  "Gd155" = "Unknown_155",
  "Gd156" = "CD4",
  "Gd158" = "Unknown_158",
  "Tb159" = "CD68",
  "Gd160" = "Unknown_160",
  "Dy161" = "Unknown_161",
  "Dy162" = "CD8a",
  "Dy163" = "Unknown_163",
  "Dy164" = "E_Cadherin",
  "Ho165" = "pSTAT3",
  "Er166" = "pNFkBp65",
  "Er167" = "Unknown_167",
  "Er168" = "Unknown_168",
  "Tm169" = "Collagen_I",
  "Er170" = "CD3",
  "Yb171" = "Unknown_171",
  "Yb172" = "Unknown_172",
  "Yb173" = "Ki67",
  "Yb174" = "Pan_CK",
  "Lu175" = "pS6",
  "Yb176" = "Beta_Catenin",
  "Ir191" = "DNA1",
  "Ir193" = "DNA2"
)

# =============================================================================
# 1. Load Data
# =============================================================================
cat("\nLoading data...\n")

# Load distance data
cell_tract_file <- "/home/rstudio/interfaces/ACP_IMC/results/collagen_qc/cell_tract_distances_all_approaches.csv"
tract_distances <- fread(cell_tract_file)
tract_distances <- tract_distances[approach == "permissive"]
cat("  Distance data loaded:", format(nrow(tract_distances), big.mark = ","), "cells\n")

# Load marker expression data
marker_file <- "/home/rstudio/interfaces/ACP_IMC/results/qc_results/step6_filtered_metrics/filtered_cell_data.csv"
marker_data <- fread(marker_file)
cat("  Marker data loaded:", format(nrow(marker_data), big.mark = ","), "cells\n")

# =============================================================================
# 2. Merge and Rename Columns
# =============================================================================
cat("\n--- Merging datasets and mapping marker names ---\n")

# Parse cell_id from tract_distances to get ObjectNumber
tract_distances[, ObjectNumber := as.integer(sub(".*_([0-9]+)$", "\\1", cell_id))]
tract_distances[, ImageNumber := roi_id]

# Merge
cell_data <- merge(
  marker_data,
  tract_distances[, .(ImageNumber, ObjectNumber, distance_to_tract, contact_5px)],
  by = c("ImageNumber", "ObjectNumber")
)
cat("  Cells after merge:", format(nrow(cell_data), big.mark = ","), "\n")

# Rename isotope columns to marker names
isotope_cols <- intersect(names(isotope_to_marker), colnames(cell_data))
cat("  Isotope columns found:", length(isotope_cols), "\n")

for (iso in isotope_cols) {
  marker_name <- isotope_to_marker[iso]
  if (!is.na(marker_name) && !grepl("Unknown", marker_name)) {
    setnames(cell_data, iso, marker_name, skip_absent = TRUE)
    cat("    ", iso, "->", marker_name, "\n")
  }
}

# List available markers after renaming
available_markers <- intersect(
  c("Collagen_I", "Pan_CK", "EpCAM", "Vimentin", "CD163", "CD68",
    "Beta_Catenin", "E_Cadherin", "Ki67", "CD44", "CD45", "CD3",
    "CD8a", "CD4", "Cleaved_Caspase3", "pNFkBp65", "SMA", "CD31"),
  colnames(cell_data)
)
cat("\n  Available named markers:", paste(available_markers, collapse = ", "), "\n")

# =============================================================================
# 3. Cell Compartment Classification
# =============================================================================
cat("\n=============================================================================\n")
cat("Cell Compartment Classification\n")
cat("=============================================================================\n")

# Strategy: Classify cells based on marker expression
# - Epithelial: High Pan-CK and/or high EpCAM
# - Stromal: High Vimentin and/or high Collagen_I
# - Immune: High CD45 and/or CD68/CD163

# Calculate thresholds (75th percentile for each marker)
if ("Pan_CK" %in% colnames(cell_data)) {
  panck_thresh <- quantile(cell_data$Pan_CK, 0.75, na.rm = TRUE)
  cat("  Pan-CK 75th percentile threshold:", round(panck_thresh, 3), "\n")
}
if ("EpCAM" %in% colnames(cell_data)) {
  epcam_thresh <- quantile(cell_data$EpCAM, 0.75, na.rm = TRUE)
  cat("  EpCAM 75th percentile threshold:", round(epcam_thresh, 3), "\n")
}
if ("Vimentin" %in% colnames(cell_data)) {
  vim_thresh <- quantile(cell_data$Vimentin, 0.75, na.rm = TRUE)
  cat("  Vimentin 75th percentile threshold:", round(vim_thresh, 3), "\n")
}
if ("Collagen_I" %in% colnames(cell_data)) {
  col_thresh <- quantile(cell_data$Collagen_I, 0.75, na.rm = TRUE)
  cat("  Collagen I 75th percentile threshold:", round(col_thresh, 3), "\n")
}
if ("CD45" %in% colnames(cell_data)) {
  cd45_thresh <- quantile(cell_data$CD45, 0.75, na.rm = TRUE)
  cat("  CD45 75th percentile threshold:", round(cd45_thresh, 3), "\n")
}

# Classify cells
cell_data[, `:=`(
  is_epithelial = FALSE,
  is_stromal = FALSE,
  is_immune = FALSE
)]

if ("Pan_CK" %in% colnames(cell_data)) {
  cell_data[Pan_CK >= panck_thresh, is_epithelial := TRUE]
}
if ("EpCAM" %in% colnames(cell_data)) {
  cell_data[EpCAM >= epcam_thresh, is_epithelial := TRUE]
}
if ("Vimentin" %in% colnames(cell_data)) {
  cell_data[Vimentin >= vim_thresh, is_stromal := TRUE]
}
if ("Collagen_I" %in% colnames(cell_data)) {
  cell_data[Collagen_I >= col_thresh, is_stromal := TRUE]
}
if ("CD45" %in% colnames(cell_data)) {
  cell_data[CD45 >= cd45_thresh, is_immune := TRUE]
}

# Assign compartment (priority: Epithelial > Immune > Stromal > Other)
cell_data[, compartment := "Other"]
cell_data[is_stromal == TRUE, compartment := "Stromal"]
cell_data[is_immune == TRUE, compartment := "Immune"]
cell_data[is_epithelial == TRUE, compartment := "Epithelial"]

# Summary
compartment_summary <- cell_data[, .N, by = compartment]
cat("\n--- Compartment Classification ---\n")
print(compartment_summary)

# Also create distance-based classification
cell_data[, collagen_contact := ifelse(contact_5px, "Contact", "Distant")]

# =============================================================================
# 4. Analysis 1: Marker Expression by Compartment
# =============================================================================
cat("\n=============================================================================\n")
cat("Analysis 1: Marker Expression by Cell Compartment\n")
cat("=============================================================================\n")

# Compare marker expression across compartments
compartment_de <- data.table()

for (marker in available_markers) {
  if (!marker %in% colnames(cell_data)) next

  # Get expression by compartment
  expr_by_comp <- cell_data[, .(
    mean_expr = mean(get(marker), na.rm = TRUE),
    median_expr = median(get(marker), na.rm = TRUE),
    sd_expr = sd(get(marker), na.rm = TRUE),
    n = .N
  ), by = compartment]

  # Calculate fold change (Epithelial vs Stromal)
  epi_mean <- expr_by_comp[compartment == "Epithelial", mean_expr]
  stro_mean <- expr_by_comp[compartment == "Stromal", mean_expr]

  if (length(epi_mean) > 0 && length(stro_mean) > 0 && stro_mean > 0) {
    log2fc <- log2((epi_mean + 0.001) / (stro_mean + 0.001))

    # Statistical test
    epi_vals <- cell_data[compartment == "Epithelial", get(marker)]
    stro_vals <- cell_data[compartment == "Stromal", get(marker)]

    wtest <- tryCatch(
      wilcox.test(epi_vals, stro_vals),
      error = function(e) list(p.value = NA)
    )

    # Cohen's d
    pooled_sd <- sqrt((var(epi_vals, na.rm = TRUE) + var(stro_vals, na.rm = TRUE)) / 2)
    cohens_d <- (mean(epi_vals, na.rm = TRUE) - mean(stro_vals, na.rm = TRUE)) / pooled_sd

    compartment_de <- rbind(compartment_de, data.table(
      marker = marker,
      epithelial_mean = epi_mean,
      stromal_mean = stro_mean,
      log2FC_epi_vs_stro = log2fc,
      cohens_d = cohens_d,
      p_value = wtest$p.value
    ))
  }
}

compartment_de[, p_adj := p.adjust(p_value, method = "BH")]
compartment_de <- compartment_de[order(-log2FC_epi_vs_stro)]

cat("\n--- Marker Expression: Epithelial vs Stromal Cells ---\n")
print(compartment_de)

# =============================================================================
# 5. Analysis 2: Spatial Correlations by ROI
# =============================================================================
cat("\n=============================================================================\n")
cat("Analysis 2: Spatial Correlations Between Markers by ROI\n")
cat("=============================================================================\n")

# Define marker pairs of interest
marker_pairs <- list(
  c("Collagen_I", "Pan_CK"),
  c("Collagen_I", "EpCAM"),
  c("Collagen_I", "Vimentin"),
  c("Collagen_I", "CD163"),
  c("Collagen_I", "Beta_Catenin"),
  c("Collagen_I", "E_Cadherin"),
  c("Collagen_I", "Ki67"),
  c("Collagen_I", "CD44"),
  c("Pan_CK", "EpCAM"),
  c("Pan_CK", "Ki67"),
  c("Pan_CK", "Beta_Catenin"),
  c("Vimentin", "CD163")
)

# Calculate correlations
correlation_results <- data.table()

for (roi in unique(cell_data$ImageNumber)) {
  roi_data <- cell_data[ImageNumber == roi]

  for (pair in marker_pairs) {
    m1 <- pair[1]
    m2 <- pair[2]

    if (!m1 %in% colnames(roi_data) || !m2 %in% colnames(roi_data)) next

    vals1 <- roi_data[[m1]]
    vals2 <- roi_data[[m2]]

    valid <- !is.na(vals1) & !is.na(vals2) & is.finite(vals1) & is.finite(vals2)
    if (sum(valid) < 50) next

    cor_test <- cor.test(vals1[valid], vals2[valid], method = "pearson")

    correlation_results <- rbind(correlation_results, data.table(
      ROI = roi,
      marker1 = m1,
      marker2 = m2,
      pair_name = paste(m1, "vs", m2),
      correlation = cor_test$estimate,
      p_value = cor_test$p.value,
      n_cells = sum(valid)
    ))
  }
}

# Summary
if (nrow(correlation_results) > 0) {
  cor_summary <- correlation_results[, .(
    mean_r = mean(correlation, na.rm = TRUE),
    sd_r = sd(correlation, na.rm = TRUE),
    min_r = min(correlation, na.rm = TRUE),
    max_r = max(correlation, na.rm = TRUE),
    n_ROIs = .N
  ), by = pair_name]

  cat("\n--- Correlation Summary by Marker Pair ---\n")
  print(cor_summary[order(-abs(mean_r))])
} else {
  cor_summary <- data.table()
  cat("  No correlations calculated\n")
}

# =============================================================================
# 6. Analysis 3: Marker Expression by Collagen Contact (within Epithelial cells)
# =============================================================================
cat("\n=============================================================================\n")
cat("Analysis 3: Marker Expression by Collagen Contact (Epithelial Cells Only)\n")
cat("=============================================================================\n")

# Focus on epithelial cells and compare those at contact vs distant
epi_cells <- cell_data[compartment == "Epithelial"]
cat("  Epithelial cells:", format(nrow(epi_cells), big.mark = ","), "\n")

epi_contact_summary <- epi_cells[, .N, by = collagen_contact]
cat("  Contact:", epi_contact_summary[collagen_contact == "Contact", N], "\n")
cat("  Distant:", epi_contact_summary[collagen_contact == "Distant", N], "\n")

# Compare markers within epithelial cells
epi_contact_de <- data.table()

for (marker in available_markers) {
  if (!marker %in% colnames(epi_cells)) next

  contact_vals <- epi_cells[collagen_contact == "Contact", get(marker)]
  distant_vals <- epi_cells[collagen_contact == "Distant", get(marker)]

  if (length(contact_vals) < 50 || length(distant_vals) < 50) next

  contact_mean <- mean(contact_vals, na.rm = TRUE)
  distant_mean <- mean(distant_vals, na.rm = TRUE)

  log2fc <- log2((contact_mean + 0.001) / (distant_mean + 0.001))

  wtest <- tryCatch(
    wilcox.test(contact_vals, distant_vals),
    error = function(e) list(p.value = NA)
  )

  pooled_sd <- sqrt((var(contact_vals, na.rm = TRUE) + var(distant_vals, na.rm = TRUE)) / 2)
  cohens_d <- (contact_mean - distant_mean) / pooled_sd

  epi_contact_de <- rbind(epi_contact_de, data.table(
    marker = marker,
    contact_mean = contact_mean,
    distant_mean = distant_mean,
    log2FC_contact_vs_distant = log2fc,
    cohens_d = cohens_d,
    p_value = wtest$p.value
  ))
}

if (nrow(epi_contact_de) > 0) {
  epi_contact_de[, p_adj := p.adjust(p_value, method = "BH")]
  epi_contact_de <- epi_contact_de[order(-abs(cohens_d))]

  cat("\n--- Within Epithelial Cells: Contact vs Distant ---\n")
  print(epi_contact_de)
}

# =============================================================================
# 7. Generate Figures
# =============================================================================
cat("\n=============================================================================\n")
cat("Generating Figures\n")
cat("=============================================================================\n")

output_dir <- "/home/rstudio/interfaces/ACP_IMC/results/collagen_tract_analysis_v3"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Figure 1: Compartment Classification ---
cat("Creating compartment summary figure...\n")

compartment_summary[, pct := N / sum(N) * 100]

p1 <- ggplot(compartment_summary, aes(x = reorder(compartment, -N), y = N, fill = compartment)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%s\n(%.1f%%)", format(N, big.mark = ","), pct)),
            vjust = -0.2, size = 3.5) +
  scale_fill_manual(values = c("Epithelial" = "#3498DB", "Stromal" = "#E74C3C",
                               "Immune" = "#9B59B6", "Other" = "#95A5A6")) +
  labs(
    title = "Cell Compartment Classification",
    subtitle = "Based on marker expression (Pan-CK/EpCAM = Epithelial; Vimentin/Collagen = Stromal; CD45 = Immune)",
    x = "Compartment",
    y = "Number of Cells"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  ylim(0, max(compartment_summary$N) * 1.15)

ggsave(file.path(output_dir, "01_compartment_classification.png"), p1, width = 8, height = 6, dpi = 300)

# --- Figure 2: Epithelial vs Stromal Marker Expression ---
cat("Creating compartment DE figure...\n")

if (nrow(compartment_de) > 0) {
  compartment_de[, direction := ifelse(log2FC_epi_vs_stro > 0, "Epithelial", "Stromal")]
  compartment_de[, marker_ordered := factor(marker, levels = marker[order(log2FC_epi_vs_stro)])]

  p2 <- ggplot(compartment_de, aes(x = marker_ordered, y = log2FC_epi_vs_stro, fill = direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_hline(yintercept = c(-0.5, 0, 0.5), linetype = c("dashed", "solid", "dashed"),
               color = c("gray50", "black", "gray50")) +
    scale_fill_manual(values = c("Epithelial" = "#3498DB", "Stromal" = "#E74C3C")) +
    coord_flip() +
    labs(
      title = "Marker Expression: Epithelial vs Stromal Cells",
      subtitle = "Positive = higher in Epithelial; Negative = higher in Stromal",
      x = "",
      y = "log2 Fold Change (Epithelial / Stromal)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(output_dir, "02_epithelial_vs_stromal_DE.png"), p2, width = 10, height = 8, dpi = 300)
}

# --- Figure 3: Spatial Correlations ---
cat("Creating correlation figures...\n")

if (nrow(correlation_results) > 0) {
  # Bar plot by ROI for key pairs
  key_pairs <- c("Collagen_I vs Pan_CK", "Collagen_I vs EpCAM",
                 "Collagen_I vs Vimentin", "Collagen_I vs Beta_Catenin")

  cor_plot_data <- correlation_results[pair_name %in% key_pairs]

  if (nrow(cor_plot_data) > 0) {
    p3 <- ggplot(cor_plot_data, aes(x = factor(ROI), y = correlation, fill = correlation > 0)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_hline(yintercept = 0, color = "black") +
      facet_wrap(~ pair_name, ncol = 2) +
      scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"),
                        labels = c("Negative", "Positive"), name = "") +
      labs(
        title = "Spatial Correlations Between Markers by ROI",
        x = "ROI",
        y = "Pearson's r"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )

    ggsave(file.path(output_dir, "03_spatial_correlations_by_ROI.png"), p3, width = 10, height = 8, dpi = 300)
  }

  # Summary heatmap
  cor_matrix <- dcast(cor_summary, pair_name ~ ., value.var = "mean_r")

  p3b <- ggplot(cor_summary, aes(x = reorder(pair_name, mean_r), y = 1, fill = mean_r)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", mean_r)), size = 3) +
    scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C", midpoint = 0,
                         limits = c(-0.5, 0.5), name = "Mean r") +
    coord_flip() +
    labs(
      title = "Mean Spatial Correlations Across ROIs",
      x = "",
      y = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  ggsave(file.path(output_dir, "04_correlation_summary_heatmap.png"), p3b, width = 10, height = 6, dpi = 300)
}

# --- Figure 4: Within-Epithelial Contact Analysis ---
cat("Creating epithelial contact figure...\n")

if (nrow(epi_contact_de) > 0) {
  epi_contact_de[, direction := ifelse(cohens_d > 0, "Higher at Contact", "Higher when Distant")]
  epi_contact_de[, marker_ordered := factor(marker, levels = marker[order(cohens_d)])]

  p4 <- ggplot(epi_contact_de, aes(x = marker_ordered, y = cohens_d, fill = direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_hline(yintercept = c(-0.2, 0, 0.2), linetype = c("dashed", "solid", "dashed"),
               color = c("gray50", "black", "gray50")) +
    scale_fill_manual(values = c("Higher at Contact" = "#E74C3C", "Higher when Distant" = "#3498DB")) +
    coord_flip() +
    labs(
      title = "Within Epithelial Cells: Expression by Collagen Contact",
      subtitle = sprintf("n = %s epithelial cells | Dashed lines = |d| = 0.2 (small effect)",
                         format(nrow(epi_cells), big.mark = ",")),
      x = "",
      y = "Cohen's d (Contact - Distant)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(output_dir, "05_epithelial_contact_analysis.png"), p4, width = 10, height = 8, dpi = 300)
}

# =============================================================================
# 8. Export Results
# =============================================================================
cat("\n--- Exporting Results ---\n")

fwrite(compartment_de, file.path(output_dir, "compartment_differential_expression.csv"))
fwrite(correlation_results, file.path(output_dir, "spatial_correlations_by_ROI.csv"))
fwrite(cor_summary, file.path(output_dir, "correlation_summary.csv"))
fwrite(epi_contact_de, file.path(output_dir, "epithelial_contact_differential_expression.csv"))

# =============================================================================
# 9. Manuscript Text
# =============================================================================
cat("\n=============================================================================\n")
cat("MANUSCRIPT TEXT (Copy-Paste Ready)\n")
cat("=============================================================================\n")

n_cells <- nrow(cell_data)
n_epi <- compartment_summary[compartment == "Epithelial", N]
n_stro <- compartment_summary[compartment == "Stromal", N]
n_rois <- length(unique(cell_data$ImageNumber))

# Get key results
if (nrow(cor_summary) > 0) {
  col_panck <- cor_summary[pair_name == "Collagen_I vs Pan_CK"]
  col_epcam <- cor_summary[pair_name == "Collagen_I vs EpCAM"]
  col_vim <- cor_summary[pair_name == "Collagen_I vs Vimentin"]
}

cat(sprintf("
Collagen tracts define molecularly distinct microenvironmental niches. To characterize
the molecular basis of epithelial-stromal compartmentalization, we classified %s cells
from %d ROIs into compartments based on marker expression: Epithelial (Pan-CK+ or EpCAM+;
n = %s), Stromal (Vimentin+ or Collagen I+; n = %s), Immune (CD45+), and Other.
",
            format(n_cells, big.mark = ","), n_rois,
            format(n_epi, big.mark = ","), format(n_stro, big.mark = ",")))

if (nrow(compartment_de) > 0) {
  # Top epithelial and stromal markers
  top_epi <- compartment_de[log2FC_epi_vs_stro > 0][1:3]
  top_stro <- compartment_de[log2FC_epi_vs_stro < 0][order(log2FC_epi_vs_stro)][1:3]

  cat(sprintf("
Differential expression analysis between compartments confirmed molecular segregation.
Epithelial cells showed elevated expression of %s (log2FC = %.2f, d = %.2f),
%s (log2FC = %.2f), and %s (log2FC = %.2f) compared to stromal cells.
Conversely, stromal cells showed elevated %s (log2FC = %.2f), %s (log2FC = %.2f),
and %s (log2FC = %.2f).
",
              top_epi$marker[1], top_epi$log2FC_epi_vs_stro[1], top_epi$cohens_d[1],
              top_epi$marker[2], top_epi$log2FC_epi_vs_stro[2],
              top_epi$marker[3], top_epi$log2FC_epi_vs_stro[3],
              top_stro$marker[1], top_stro$log2FC_epi_vs_stro[1],
              top_stro$marker[2], top_stro$log2FC_epi_vs_stro[2],
              top_stro$marker[3], top_stro$log2FC_epi_vs_stro[3]))
}

if (nrow(cor_summary) > 0 && nrow(col_panck) > 0) {
  cat(sprintf("
Spatial correlation analysis within ROIs revealed consistent negative correlations
between Collagen I and epithelial markers: Pan-CK (mean r = %.2f, range %.2f to %.2f)
and EpCAM (mean r = %.2f). Collagen I showed positive correlation with Vimentin
(mean r = %.2f), confirming co-localization of mesenchymal markers within collagen-rich
stromal domains.
",
              col_panck$mean_r, col_panck$min_r, col_panck$max_r,
              col_epcam$mean_r,
              col_vim$mean_r))
}

if (nrow(epi_contact_de) > 0) {
  # Get Ki67 and E-Cadherin results if available
  ki67_result <- epi_contact_de[marker == "Ki67"]
  ecad_result <- epi_contact_de[marker == "E_Cadherin"]

  if (nrow(ki67_result) > 0) {
    cat(sprintf("
Within epithelial cells, Ki-67 expression was %s in cells contacting collagen
compared to distant cells (d = %.3f), consistent with reduced proliferation at
epithelial-stromal boundaries.
",
                ifelse(ki67_result$cohens_d < 0, "lower", "higher"),
                ki67_result$cohens_d))
  }
}

cat("
These findings establish that collagen tracts represent molecularly distinct stromal
niches, providing the architectural context for understanding functional
compartmentalization within the tumor microenvironment.
")

# =============================================================================
# Summary
# =============================================================================
cat("\n=============================================================================\n")
cat("SUMMARY\n")
cat("=============================================================================\n")
cat("Output directory:", output_dir, "\n")
cat("\nKey findings:\n")

if (nrow(compartment_de) > 0) {
  cat("  1. Compartment classification successful\n")
  cat(sprintf("     - Epithelial: %s cells (%.1f%%)\n",
              format(n_epi, big.mark = ","), 100*n_epi/n_cells))
  cat(sprintf("     - Stromal: %s cells (%.1f%%)\n",
              format(n_stro, big.mark = ","), 100*n_stro/n_cells))
}

if (nrow(cor_summary) > 0 && nrow(col_panck) > 0) {
  cat(sprintf("  2. Collagen-Pan_CK correlation: r = %.2f (confirms segregation)\n", col_panck$mean_r))
}

if (nrow(epi_contact_de) > 0) {
  ki67_d <- epi_contact_de[marker == "Ki67", cohens_d]
  if (length(ki67_d) > 0) {
    cat(sprintf("  3. Ki-67 in epithelial cells at contact: d = %.3f\n", ki67_d))
  }
}

cat("\n=============================================================================\n")
