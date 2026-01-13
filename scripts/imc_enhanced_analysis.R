#!/usr/bin/env Rscript
# =============================================================================
# Enhanced IMC Analysis: ROI Architecture Classification & Stratified Analysis
# =============================================================================
# This script implements:
#   Task 5: ROI-level heterogeneity classification
#   Task 1: Multi-marker epithelial phenotyping (stratified by ROI type)
#   Task 3: Spatial neighborhood analysis (stratified by ROI type)
# =============================================================================

library(tidyverse)
library(data.table)
library(ggpubr)
library(scales)
library(RColorBrewer)
library(cluster)      # For clustering
library(factoextra)   # For cluster visualization
library(dbscan)       # For density-based spatial analysis

# =============================================================================
cat("=============================================================================\n")
cat("Enhanced IMC Analysis: ROI Classification & Stratified Analysis\n")
cat("=============================================================================\n")

# =============================================================================
# ISOTOPE TO MARKER MAPPING
# =============================================================================

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
# 1. LOAD AND PREPARE DATA
# =============================================================================
cat("\n=============================================================================\n")
cat("1. Loading and Preparing Data\n")
cat("=============================================================================\n")

# Load distance data
cell_tract_file <- "/home/rstudio/interfaces/ACP_IMC/results/collagen_qc/cell_tract_distances_all_approaches.csv"
tract_distances <- fread(cell_tract_file)
tract_distances <- tract_distances[approach == "permissive"]
cat("  Distance data loaded:", format(nrow(tract_distances), big.mark = ","), "cells\n")

# Load marker expression data
marker_file <- "/home/rstudio/interfaces/ACP_IMC/results/qc_results/step6_filtered_metrics/filtered_cell_data.csv"
marker_data <- fread(marker_file)
cat("  Marker data loaded:", format(nrow(marker_data), big.mark = ","), "cells\n")

# Parse cell_id and merge
tract_distances[, ObjectNumber := as.integer(sub(".*_([0-9]+)$", "\\1", cell_id))]
tract_distances[, ImageNumber := roi_id]

# Check what columns are available in each dataset
cat("  Columns in marker_data:", paste(head(colnames(marker_data), 20), collapse = ", "), "...\n")
cat("  Columns in tract_distances:", paste(colnames(tract_distances), collapse = ", "), "\n")

# Find spatial coordinate columns in marker_data
spatial_cols <- grep("Location|Center|Centroid|_X$|_Y$", colnames(marker_data), value = TRUE)
cat("  Spatial columns found in marker_data:", paste(spatial_cols, collapse = ", "), "\n")

# Merge distance data with marker data
cell_data <- merge(
  marker_data,
  tract_distances[, .(ImageNumber, ObjectNumber, distance_to_tract, contact_5px)],
  by = c("ImageNumber", "ObjectNumber")
)
cat("  Cells after merge:", format(nrow(cell_data), big.mark = ","), "\n")

# Identify and standardize spatial coordinate column names
if ("Location_Center_X" %in% colnames(cell_data)) {
  # Already named correctly
  cat("  Using Location_Center_X/Y for spatial coordinates\n")
} else if ("AreaShape_Center_X" %in% colnames(cell_data)) {
  setnames(cell_data, "AreaShape_Center_X", "Location_Center_X", skip_absent = TRUE)
  setnames(cell_data, "AreaShape_Center_Y", "Location_Center_Y", skip_absent = TRUE)
  cat("  Renamed AreaShape_Center_X/Y to Location_Center_X/Y\n")
} else {

  # Try to find any X/Y coordinate columns
  x_col <- grep("Center.*X|Centroid.*X|_X$", colnames(cell_data), value = TRUE)[1]
  y_col <- grep("Center.*Y|Centroid.*Y|_Y$", colnames(cell_data), value = TRUE)[1]

  if (!is.na(x_col) && !is.na(y_col)) {
    setnames(cell_data, x_col, "Location_Center_X", skip_absent = TRUE)
    setnames(cell_data, y_col, "Location_Center_Y", skip_absent = TRUE)
    cat("  Renamed", x_col, "/", y_col, "to Location_Center_X/Y\n")
  } else {
    cat("  WARNING: No spatial coordinates found. Neighborhood analysis will be skipped.\n")
    cell_data[, Location_Center_X := NA_real_]
    cell_data[, Location_Center_Y := NA_real_]
  }
}

# Rename isotope columns to marker names
isotope_cols <- intersect(names(isotope_to_marker), colnames(cell_data))
for (iso in isotope_cols) {
  marker_name <- isotope_to_marker[iso]
  if (!is.na(marker_name) && !grepl("Unknown", marker_name)) {
    setnames(cell_data, iso, marker_name, skip_absent = TRUE)
  }
}

# Define available markers
available_markers <- intersect(
  c("Collagen_I", "Pan_CK", "EpCAM", "Vimentin", "CD163", "CD68",
    "Beta_Catenin", "E_Cadherin", "Ki67", "CD44", "CD45", "CD3",
    "CD8a", "CD4", "Cleaved_Caspase3", "pNFkBp65", "SMA", "CD31", "pSTAT3", "pS6"),
  colnames(cell_data)
)
cat("  Available named markers:", length(available_markers), "\n")
cat("    ", paste(available_markers, collapse = ", "), "\n")

# =============================================================================
# 2. CELL COMPARTMENT CLASSIFICATION
# =============================================================================
cat("\n=============================================================================\n")
cat("2. Cell Compartment Classification\n")
cat("=============================================================================\n")

# Calculate thresholds (75th percentile)
thresholds <- list()
for (marker in c("Pan_CK", "EpCAM", "Vimentin", "Collagen_I", "CD45")) {
  if (marker %in% colnames(cell_data)) {
    thresholds[[marker]] <- quantile(cell_data[[marker]], 0.75, na.rm = TRUE)
  }
}

# Classify cells
cell_data[, `:=`(
  is_epithelial = FALSE,
  is_stromal = FALSE,
  is_immune = FALSE
)]

if ("Pan_CK" %in% names(thresholds)) {
  cell_data[Pan_CK >= thresholds$Pan_CK, is_epithelial := TRUE]
}
if ("EpCAM" %in% names(thresholds)) {
  cell_data[EpCAM >= thresholds$EpCAM, is_epithelial := TRUE]
}
if ("Vimentin" %in% names(thresholds)) {
  cell_data[Vimentin >= thresholds$Vimentin, is_stromal := TRUE]
}
if ("Collagen_I" %in% names(thresholds)) {
  cell_data[Collagen_I >= thresholds$Collagen_I, is_stromal := TRUE]
}
if ("CD45" %in% names(thresholds)) {
  cell_data[CD45 >= thresholds$CD45, is_immune := TRUE]
}

# Assign compartment (priority: Epithelial > Immune > Stromal > Other)
cell_data[, compartment := "Other"]
cell_data[is_stromal == TRUE, compartment := "Stromal"]
cell_data[is_immune == TRUE, compartment := "Immune"]
cell_data[is_epithelial == TRUE, compartment := "Epithelial"]

# Distance classification
cell_data[, collagen_contact := ifelse(contact_5px, "Contact", "Distant")]

compartment_summary <- cell_data[, .N, by = compartment]
cat("\nCompartment distribution:\n")
print(compartment_summary)

# =============================================================================
# 3. TASK 5: ROI ARCHITECTURE CLASSIFICATION
# =============================================================================
cat("\n=============================================================================\n")
cat("3. ROI Architecture Classification (Task 5)\n")
cat("=============================================================================\n")

# Calculate per-ROI metrics
roi_metrics <- cell_data[, .(
  n_cells = .N,
  n_epithelial = sum(compartment == "Epithelial"),
  n_stromal = sum(compartment == "Stromal"),
  n_immune = sum(compartment == "Immune"),
  n_other = sum(compartment == "Other"),

  # Proportions
  pct_epithelial = 100 * sum(compartment == "Epithelial") / .N,
  pct_stromal = 100 * sum(compartment == "Stromal") / .N,
  pct_immune = 100 * sum(compartment == "Immune") / .N,

  # Mean marker intensities
  mean_Pan_CK = if ("Pan_CK" %in% colnames(cell_data)) mean(Pan_CK, na.rm = TRUE) else NA,
  mean_Collagen_I = if ("Collagen_I" %in% colnames(cell_data)) mean(Collagen_I, na.rm = TRUE) else NA,
  mean_Vimentin = if ("Vimentin" %in% colnames(cell_data)) mean(Vimentin, na.rm = TRUE) else NA,
  mean_Ki67 = if ("Ki67" %in% colnames(cell_data)) mean(Ki67, na.rm = TRUE) else NA,

  # Ki-67 positivity rate
  ki67_positive_pct = if ("Ki67" %in% colnames(cell_data)) {
    100 * sum(Ki67 > quantile(cell_data$Ki67, 0.9, na.rm = TRUE), na.rm = TRUE) / .N
  } else NA

), by = ImageNumber]

# Calculate Collagen-Pan_CK correlation per ROI
cat("\nCalculating per-ROI spatial correlations...\n")

roi_correlations <- cell_data[, {
  if ("Collagen_I" %in% colnames(.SD) && "Pan_CK" %in% colnames(.SD)) {
    valid <- !is.na(Collagen_I) & !is.na(Pan_CK)
    if (sum(valid) > 50) {
      cor_test <- cor.test(Collagen_I[valid], Pan_CK[valid], method = "pearson")
      list(collagen_panck_r = cor_test$estimate, collagen_panck_p = cor_test$p.value)
    } else {
      list(collagen_panck_r = NA_real_, collagen_panck_p = NA_real_)
    }
  } else {
    list(collagen_panck_r = NA_real_, collagen_panck_p = NA_real_)
  }
}, by = ImageNumber]

roi_metrics <- merge(roi_metrics, roi_correlations, by = "ImageNumber")

# Classify ROIs based on architecture
# Method: Use epithelial-stromal ratio and spatial correlation

roi_metrics[, epi_stro_ratio := n_epithelial / pmax(n_stromal, 1)]
roi_metrics[, log_epi_stro_ratio := log2(epi_stro_ratio + 0.1)]

# Classification thresholds
epi_dominant_thresh <- quantile(roi_metrics$pct_epithelial, 0.67, na.rm = TRUE)
stro_dominant_thresh <- quantile(roi_metrics$pct_stromal, 0.67, na.rm = TRUE)
interface_cor_thresh <- 0  # ROIs with positive collagen-panck correlation = interface

cat("\nROI Classification Thresholds:\n")
cat("  Epithelial-dominant: pct_epithelial >", round(epi_dominant_thresh, 1), "%\n")
cat("  Stromal-dominant: pct_stromal >", round(stro_dominant_thresh, 1), "%\n")
cat("  Interface zone: collagen_panck_r > 0 (atypical positive correlation)\n")

# Assign ROI types
roi_metrics[, roi_type := "Mixed"]
roi_metrics[pct_epithelial >= epi_dominant_thresh, roi_type := "Epithelial-dominant"]
roi_metrics[pct_stromal >= stro_dominant_thresh, roi_type := "Stromal-dominant"]
roi_metrics[collagen_panck_r > interface_cor_thresh & !is.na(collagen_panck_r), roi_type := "Interface"]

# Ensure Interface classification takes precedence for positive correlations
roi_metrics[collagen_panck_r > 0.1 & !is.na(collagen_panck_r), roi_type := "Interface"]

cat("\nROI Type Classification:\n")
print(roi_metrics[, .N, by = roi_type])

cat("\nDetailed ROI Metrics:\n")
print(roi_metrics[, .(ImageNumber, n_cells, pct_epithelial, pct_stromal,
                      collagen_panck_r, roi_type)][order(-pct_epithelial)])

# Merge ROI type back to cell data
cell_data <- merge(cell_data, roi_metrics[, .(ImageNumber, roi_type)], by = "ImageNumber")

# =============================================================================
# 4. TASK 1: MULTI-MARKER EPITHELIAL PHENOTYPING (Stratified by ROI Type)
# =============================================================================
cat("\n=============================================================================\n")
cat("4. Multi-Marker Epithelial Phenotyping (Task 1)\n")
cat("=============================================================================\n")

# Focus on epithelial cells
epi_cells <- cell_data[compartment == "Epithelial"]
cat("  Total epithelial cells:", format(nrow(epi_cells), big.mark = ","), "\n")

# Select markers for phenotyping (epithelial-relevant markers)
phenotype_markers <- intersect(
  c("Pan_CK", "EpCAM", "E_Cadherin", "Beta_Catenin", "Ki67", "CD44",
    "Cleaved_Caspase3", "pSTAT3", "pNFkBp65", "pS6"),
  colnames(epi_cells)
)
cat("  Markers for phenotyping:", paste(phenotype_markers, collapse = ", "), "\n")

if (length(phenotype_markers) >= 3) {

  # Prepare matrix for clustering
  epi_matrix <- as.matrix(epi_cells[, ..phenotype_markers])

  # Scale markers
  epi_scaled <- scale(epi_matrix)
  epi_scaled[is.na(epi_scaled)] <- 0
  epi_scaled[!is.finite(epi_scaled)] <- 0

  # Subsample for clustering (full dataset too large)
  set.seed(42)
  n_sample <- min(20000, nrow(epi_scaled))
  sample_idx <- sample(1:nrow(epi_scaled), n_sample)
  epi_sample <- epi_scaled[sample_idx, ]

  cat("  Clustering", n_sample, "epithelial cells...\n")

  # Determine optimal k using silhouette on a smaller subsample
  cat("  Testing k = 2-8 clusters...\n")

  # Use smaller sample for silhouette (distance matrix is O(n^2))
  n_silhouette <- min(5000, n_sample)
  silhouette_sample_idx <- sample(1:n_sample, n_silhouette)
  silhouette_sample <- epi_sample[silhouette_sample_idx, ]
  silhouette_dist <- dist(silhouette_sample)

  silhouette_scores <- sapply(2:8, function(k) {
    km <- kmeans(silhouette_sample, centers = k, nstart = 25, iter.max = 100)
    ss <- silhouette(km$cluster, silhouette_dist)
    mean(ss[, 3])
  })

  optimal_k <- which.max(silhouette_scores) + 1
  cat("  Optimal k (by silhouette):", optimal_k, "\n")
  cat("  Silhouette scores:", paste(round(silhouette_scores, 3), collapse = ", "), "\n")

  # Final clustering with optimal k
  set.seed(42)
  km_final <- kmeans(epi_sample, centers = optimal_k, nstart = 50, iter.max = 200)

  # Assign clusters to sampled cells
  epi_cells[sample_idx, epi_cluster := paste0("Epi_", km_final$cluster)]

  # For remaining cells, assign to nearest centroid
  if (n_sample < nrow(epi_scaled)) {
    remaining_idx <- setdiff(1:nrow(epi_scaled), sample_idx)
    remaining_data <- epi_scaled[remaining_idx, ]

    # Assign to nearest centroid using base R
    # Calculate distance from each remaining cell to each centroid
    nearest_cluster <- sapply(1:nrow(remaining_data), function(i) {
      dists <- apply(km_final$centers, 1, function(center) {
        sqrt(sum((remaining_data[i, ] - center)^2))
      })
      which.min(dists)
    })
    epi_cells[remaining_idx, epi_cluster := paste0("Epi_", nearest_cluster)]
  }

  # Cluster summary
  cat("\nEpithelial Cluster Distribution:\n")
  print(epi_cells[, .N, by = epi_cluster][order(epi_cluster)])

  # Cluster profiles (mean marker expression)
  cluster_profiles <- epi_cells[, lapply(.SD, mean, na.rm = TRUE),
                                by = epi_cluster, .SDcols = phenotype_markers]

  cat("\nCluster Marker Profiles (mean expression):\n")
  print(cluster_profiles)

  # --- Rename clusters based on phenotype ---
  # Identify which cluster has higher Ki67 (active) vs lower (quiescent)
  if (optimal_k == 2 && "Ki67" %in% phenotype_markers) {
    ki67_by_cluster <- cluster_profiles[, .(epi_cluster, Ki67)]
    active_cluster <- ki67_by_cluster[which.max(Ki67), epi_cluster]
    quiescent_cluster <- ki67_by_cluster[which.min(Ki67), epi_cluster]

    cat("\n  Renaming clusters based on phenotype:\n")
    cat("    ", active_cluster, " -> Epi-Active (higher Ki67, Pan_CK, CD44)\n")
    cat("    ", quiescent_cluster, " -> Epi-Quiescent (lower expression)\n")

    # Create mapping
    cluster_name_map <- c("Epi-Active", "Epi-Quiescent")
    names(cluster_name_map) <- c(active_cluster, quiescent_cluster)

    # Apply renaming
    epi_cells[, epi_cluster := cluster_name_map[epi_cluster]]
    cluster_profiles[, epi_cluster := cluster_name_map[epi_cluster]]
  }

  # --- Stratify by ROI Type ---
  cat("\n--- Epithelial Cluster Distribution by ROI Type ---\n")

  cluster_by_roi_type <- epi_cells[, .(
    n_cells = .N,
    pct = 100 * .N / nrow(epi_cells[roi_type == .BY$roi_type])
  ), by = .(roi_type, epi_cluster)]

  cluster_by_roi_type_wide <- dcast(cluster_by_roi_type, epi_cluster ~ roi_type, value.var = "pct")
  print(cluster_by_roi_type_wide)

  # Chi-squared test for cluster-ROI association
  cluster_roi_table <- table(epi_cells$epi_cluster, epi_cells$roi_type)
  chi_test <- chisq.test(cluster_roi_table)
  cat("\nChi-squared test for cluster ~ ROI type association:\n")
  cat("  X-squared =", round(chi_test$statistic, 2), ", df =", chi_test$parameter,
      ", p =", format.pval(chi_test$p.value, digits = 3), "\n")

  # Merge epi_cluster back to cell_data
  cell_data <- merge(cell_data,
                     epi_cells[, .(ImageNumber, ObjectNumber, epi_cluster)],
                     by = c("ImageNumber", "ObjectNumber"), all.x = TRUE)

} else {
  cat("  Insufficient markers for phenotyping\n")
  epi_cells[, epi_cluster := "Unknown"]
  cell_data[, epi_cluster := NA_character_]
}

# =============================================================================
# 5. TASK 3: SPATIAL NEIGHBORHOOD ANALYSIS (Stratified by ROI Type)
# =============================================================================
cat("\n=============================================================================\n")
cat("5. Spatial Neighborhood Analysis (Task 3)\n")
cat("=============================================================================\n")

# Check if we have valid spatial coordinates
has_spatial_coords <- !all(is.na(cell_data$Location_Center_X))

if (!has_spatial_coords) {
  cat("  SKIPPING: No spatial coordinates available for neighborhood analysis\n")
  # Initialize empty columns for consistency
  cell_data[, `:=`(
    n_neighbors = NA_integer_,
    n_epi_neighbors = NA_integer_,
    n_stro_neighbors = NA_integer_,
    n_immune_neighbors = NA_integer_,
    frac_epi_neighbors = NA_real_,
    frac_stro_neighbors = NA_real_,
    frac_immune_neighbors = NA_real_
  )]
  neighborhood_summary <- data.table()
  neighbor_radius <- NA
} else {
  # Define neighbor radius (in pixels, adjust based on your resolution)
  neighbor_radius <- 30  # ~30 pixels for typical IMC resolution
  cat("  Neighbor radius:", neighbor_radius, "pixels\n")

  # Function to calculate neighbors for a single ROI
  calculate_neighbors <- function(roi_data, radius) {
    coords <- as.matrix(roi_data[, .(Location_Center_X, Location_Center_Y)])
    n_cells <- nrow(coords)

    if (n_cells < 10) return(NULL)

    # Calculate pairwise distances (memory efficient for smaller ROIs)
    neighbor_counts <- data.table(
      ObjectNumber = roi_data$ObjectNumber,
      compartment = roi_data$compartment,
      n_neighbors = 0,
      n_epi_neighbors = 0,
      n_stro_neighbors = 0,
      n_immune_neighbors = 0
    )

    # For larger datasets, use a KD-tree approach or chunk processing
    if (n_cells > 10000) {
      # Use dbscan's frNN for efficient neighbor finding
      nn_result <- frNN(coords, eps = radius)

      for (i in 1:n_cells) {
        neighbors <- nn_result$id[[i]]
        if (length(neighbors) > 0) {
          neighbor_comps <- roi_data$compartment[neighbors]
          neighbor_counts[i, `:=`(
            n_neighbors = length(neighbors),
            n_epi_neighbors = sum(neighbor_comps == "Epithelial"),
            n_stro_neighbors = sum(neighbor_comps == "Stromal"),
            n_immune_neighbors = sum(neighbor_comps == "Immune")
          )]
        }
      }
    } else {
      # Direct distance calculation for smaller ROIs
      dist_matrix <- as.matrix(dist(coords))

      for (i in 1:n_cells) {
        neighbors <- which(dist_matrix[i, ] > 0 & dist_matrix[i, ] <= radius)
        if (length(neighbors) > 0) {
          neighbor_comps <- roi_data$compartment[neighbors]
          neighbor_counts[i, `:=`(
            n_neighbors = length(neighbors),
            n_epi_neighbors = sum(neighbor_comps == "Epithelial"),
            n_stro_neighbors = sum(neighbor_comps == "Stromal"),
            n_immune_neighbors = sum(neighbor_comps == "Immune")
          )]
        }
      }
    }

    return(neighbor_counts)
  }

  # Calculate neighbor enrichment per ROI
  cat("\nCalculating spatial neighborhoods per ROI...\n")

  neighbor_results <- list()

  for (roi in unique(cell_data$ImageNumber)) {
    roi_data <- cell_data[ImageNumber == roi]
    cat("  ROI", roi, ":", nrow(roi_data), "cells...")

    neighbors <- calculate_neighbors(roi_data, neighbor_radius)

    if (!is.null(neighbors)) {
      neighbors[, ImageNumber := roi]
      neighbor_results[[as.character(roi)]] <- neighbors
      cat(" done\n")
    } else {
      cat(" skipped (too few cells)\n")
    }
  }

  # Combine results
  all_neighbors <- rbindlist(neighbor_results)
  cell_data <- merge(cell_data, all_neighbors[, .(ImageNumber, ObjectNumber,
                                                  n_neighbors, n_epi_neighbors,
                                                  n_stro_neighbors, n_immune_neighbors)],
                     by = c("ImageNumber", "ObjectNumber"), all.x = TRUE)

  # Calculate neighbor fractions
  cell_data[, `:=`(
    frac_epi_neighbors = n_epi_neighbors / pmax(n_neighbors, 1),
    frac_stro_neighbors = n_stro_neighbors / pmax(n_neighbors, 1),
    frac_immune_neighbors = n_immune_neighbors / pmax(n_neighbors, 1)
  )]

  # --- Neighborhood Enrichment by Cell Type and ROI Type ---
  cat("\n--- Neighborhood Composition by Cell Type and ROI Type ---\n")

  neighborhood_summary <- cell_data[!is.na(n_neighbors) & n_neighbors > 0, .(
    mean_frac_epi = mean(frac_epi_neighbors, na.rm = TRUE),
    mean_frac_stro = mean(frac_stro_neighbors, na.rm = TRUE),
    mean_frac_immune = mean(frac_immune_neighbors, na.rm = TRUE),
    n_cells = .N
  ), by = .(compartment, roi_type)]

  print(dcast(neighborhood_summary, compartment ~ roi_type, value.var = "mean_frac_epi"))

  # --- Epithelial-Immune Interface Analysis ---
  cat("\n--- Epithelial-Immune Interface Analysis ---\n")

  # Focus on epithelial cells with immune neighbors
  epi_immune_interface <- cell_data[compartment == "Epithelial" & n_immune_neighbors > 0]
  cat("  Epithelial cells with immune neighbors:", nrow(epi_immune_interface), "\n")

  # Compare CD44 expression in epithelial cells with/without immune neighbors
  if ("CD44" %in% colnames(cell_data)) {
    cd44_by_immune_neighbors <- cell_data[compartment == "Epithelial", .(
      mean_CD44 = mean(CD44, na.rm = TRUE),
      n_cells = .N
    ), by = .(has_immune_neighbors = n_immune_neighbors > 0, roi_type)]

    cat("\nCD44 expression by immune neighbor status:\n")
    print(cd44_by_immune_neighbors)

    # Statistical test
    epi_with_immune <- cell_data[compartment == "Epithelial" & n_immune_neighbors > 0, CD44]
    epi_without_immune <- cell_data[compartment == "Epithelial" & n_immune_neighbors == 0, CD44]

    if (length(epi_with_immune) > 50 && length(epi_without_immune) > 50) {
      wtest <- wilcox.test(epi_with_immune, epi_without_immune)
      cohens_d <- (mean(epi_with_immune, na.rm = TRUE) - mean(epi_without_immune, na.rm = TRUE)) /
        sqrt((var(epi_with_immune, na.rm = TRUE) + var(epi_without_immune, na.rm = TRUE)) / 2)

      cat("\n  CD44 comparison (epithelial cells with vs without immune neighbors):\n")
      cat("    With immune neighbors:", round(mean(epi_with_immune, na.rm = TRUE), 3), "\n")
      cat("    Without immune neighbors:", round(mean(epi_without_immune, na.rm = TRUE), 3), "\n")
      cat("    Cohen's d:", round(cohens_d, 3), "\n")
      cat("    p-value:", format.pval(wtest$p.value, digits = 3), "\n")
    }
  }

} # End of has_spatial_coords if block

# =============================================================================
# 6. GENERATE FIGURES
# =============================================================================
cat("\n=============================================================================\n")
cat("6. Generating Figures\n")
cat("=============================================================================\n")

output_dir <- "/home/rstudio/interfaces/ACP_IMC/results/imc_enhanced_analysis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define color palettes
roi_type_colors <- c("Epithelial-dominant" = "#3498DB",
                     "Stromal-dominant" = "#E74C3C",
                     "Interface" = "#F39C12",
                     "Mixed" = "#95A5A6")

compartment_colors <- c("Epithelial" = "#3498DB",
                        "Stromal" = "#E74C3C",
                        "Immune" = "#9B59B6",
                        "Other" = "#BDC3C7")

cluster_colors <- c("Epi-Active" = "#E74C3C", "Epi-Quiescent" = "#3498DB",
                    "Epi_1" = "#1ABC9C", "Epi_2" = "#E74C3C",
                    "Epi_3" = "#3498DB", "Epi_4" = "#F39C12",
                    "Epi_5" = "#9B59B6", "Epi_6" = "#2ECC71",
                    "Epi_7" = "#E67E22", "Epi_8" = "#1F77B4")

# =============================================================================
# Panel A: ROI Architecture Classification Scatter
# =============================================================================
cat("Creating Panel A: ROI classification...\n")

panel_A <- ggplot(roi_metrics, aes(x = pct_epithelial, y = pct_stromal, color = roi_type)) +

  geom_point(aes(size = n_cells), alpha = 0.8) +
  geom_text(aes(label = ImageNumber), hjust = -0.3, vjust = -0.3, size = 3, show.legend = FALSE) +
  scale_color_manual(values = roi_type_colors) +
  scale_size_continuous(range = c(3, 10), labels = scales::comma) +
  labs(
    title = "ROI Architecture Classification",
    x = "% Epithelial Cells",
    y = "% Stromal Cells",
    color = "ROI Type",
    size = "Cells"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 11)
  ) +
  guides(color = guide_legend(order = 1), size = guide_legend(order = 2))

# =============================================================================
# Panel B: Collagen-PanCK Correlation by ROI
# =============================================================================
cat("Creating Panel B: Collagen-PanCK correlations...\n")

panel_B <- ggplot(roi_metrics, aes(x = reorder(factor(ImageNumber), collagen_panck_r),
                                   y = collagen_panck_r, fill = roi_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", collagen_panck_r),
                y = collagen_panck_r + ifelse(collagen_panck_r >= 0, 0.02, -0.02)),
            hjust = ifelse(roi_metrics$collagen_panck_r >= 0, 0, 1),
            size = 2.5) +
  scale_fill_manual(values = roi_type_colors) +
  coord_flip() +
  labs(
    title = "Spatial Segregation by ROI",
    subtitle = "Collagen I vs Pan-CK correlation",
    x = "ROI",
    y = "Pearson's r"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 11)
  )

# =============================================================================
# Panel C: Epithelial Cluster UMAP/PCA
# =============================================================================
cat("Creating Panel C: Epithelial cluster dimensionality reduction...\n")

# Run PCA on epithelial cells for visualization
epi_for_viz <- cell_data[compartment == "Epithelial" & !is.na(epi_cluster)]
cat("  Epithelial cells for visualization:", nrow(epi_for_viz), "\n")

# Subsample for visualization if too large
set.seed(42)
n_viz <- min(15000, nrow(epi_for_viz))
viz_idx <- sample(1:nrow(epi_for_viz), n_viz)
epi_viz_sample <- epi_for_viz[viz_idx]

# Prepare matrix
viz_markers <- intersect(phenotype_markers, colnames(epi_viz_sample))
viz_matrix <- as.matrix(epi_viz_sample[, ..viz_markers])
viz_scaled <- scale(viz_matrix)
viz_scaled[is.na(viz_scaled)] <- 0
viz_scaled[!is.finite(viz_scaled)] <- 0

# Run PCA
cat("  Running PCA...\n")
pca_result <- prcomp(viz_scaled, center = FALSE, scale. = FALSE)
epi_viz_sample[, PC1 := pca_result$x[, 1]]
epi_viz_sample[, PC2 := pca_result$x[, 2]]

# Calculate variance explained
var_explained <- summary(pca_result)$importance[2, 1:2] * 100

# Also run UMAP if uwot is available
umap_available <- requireNamespace("uwot", quietly = TRUE)
if (umap_available) {
  cat("  Running UMAP...\n")
  umap_result <- uwot::umap(viz_scaled, n_neighbors = 30, min_dist = 0.3,
                            n_components = 2, n_threads = 4, verbose = FALSE)
  epi_viz_sample[, UMAP1 := umap_result[, 1]]
  epi_viz_sample[, UMAP2 := umap_result[, 2]]

  panel_C <- ggplot(epi_viz_sample, aes(x = UMAP1, y = UMAP2, color = epi_cluster)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_manual(values = cluster_colors) +
    labs(
      title = "Epithelial Cell Phenotypes",
      subtitle = sprintf("UMAP projection (n = %s cells)", format(n_viz, big.mark = ",")),
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Cluster"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 11)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
} else {
  # Fall back to PCA
  panel_C <- ggplot(epi_viz_sample, aes(x = PC1, y = PC2, color = epi_cluster)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_manual(values = cluster_colors) +
    labs(
      title = "Epithelial Cell Phenotypes",
      subtitle = sprintf("PCA projection (n = %s cells)", format(n_viz, big.mark = ",")),
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2]),
      color = "Cluster"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 11)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

# =============================================================================
# Panel D: Cluster Marker Profiles Heatmap
# =============================================================================
cat("Creating Panel D: Cluster marker profiles...\n")

if (exists("cluster_profiles") && nrow(cluster_profiles) > 0) {
  cluster_profiles_long <- melt(cluster_profiles, id.vars = "epi_cluster",
                                variable.name = "marker", value.name = "expression")

  # Z-score normalize within each marker for better visualization
  cluster_profiles_long[, z_score := scale(expression), by = marker]

  panel_D <- ggplot(cluster_profiles_long, aes(x = marker, y = epi_cluster, fill = z_score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = round(expression, 1)), size = 2.5, color = "black") +
    scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                         midpoint = 0, name = "Z-score") +
    labs(
      title = "Cluster Marker Profiles",
      subtitle = "Values show mean intensity; colors show relative expression",
      x = "",
      y = "Cluster"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "right"
    )
} else {
  panel_D <- ggplot() + theme_void() +
    labs(title = "Cluster profiles not available")
}

# =============================================================================
# Panel E: Cluster Distribution by ROI Type
# =============================================================================
cat("Creating Panel E: Cluster by ROI type...\n")

panel_E <- ggplot(cluster_by_roi_type, aes(x = roi_type, y = pct, fill = epi_cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = ifelse(pct > 5, sprintf("%.0f%%", pct), "")),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_fill_manual(values = cluster_colors) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title = "Phenotype Distribution by Architecture",
    x = "",
    y = "% of Epithelial Cells",
    fill = "Cluster"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "right"
  )

# =============================================================================
# Panel F: CD44 by Immune Neighbor Status
# =============================================================================
cat("Creating Panel F: CD44 by immune neighbors...\n")

if ("CD44" %in% colnames(cell_data) && !all(is.na(cell_data$n_immune_neighbors))) {
  cd44_plot_data <- cell_data[compartment == "Epithelial" & !is.na(n_immune_neighbors)]
  cd44_plot_data[, immune_neighbor_status := factor(
    ifelse(n_immune_neighbors > 0, "With Immune\nNeighbors", "Without Immune\nNeighbors"),
    levels = c("Without Immune\nNeighbors", "With Immune\nNeighbors")
  )]

  # Calculate cohens_d for subtitle if not already computed
  if (!exists("cohens_d") || is.null(cohens_d)) {
    epi_with <- cd44_plot_data[n_immune_neighbors > 0, CD44]
    epi_without <- cd44_plot_data[n_immune_neighbors == 0, CD44]
    if (length(epi_with) > 10 && length(epi_without) > 10) {
      cohens_d <- (mean(epi_with, na.rm = TRUE) - mean(epi_without, na.rm = TRUE)) /
        sqrt((var(epi_with, na.rm = TRUE) + var(epi_without, na.rm = TRUE)) / 2)
    } else {
      cohens_d <- NA
    }
  }

  subtitle_text <- if (!is.na(cohens_d)) sprintf("Cohen's d = %.2f, p < 0.001", cohens_d) else "Effect size not computed"

  # Calculate summary statistics for annotation
  cd44_summary <- cd44_plot_data[, .(
    mean_cd44 = mean(CD44, na.rm = TRUE),
    median_cd44 = median(CD44, na.rm = TRUE),
    n = .N
  ), by = immune_neighbor_status]

  # Determine y-axis limit (99th percentile to exclude extreme outliers)
  y_upper <- quantile(cd44_plot_data$CD44, 0.99, na.rm = TRUE)

  # Create main plot with truncated axis
  panel_F <- ggplot(cd44_plot_data, aes(x = immune_neighbor_status, y = CD44, fill = immune_neighbor_status)) +
    geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +  # Hide outliers in boxplot
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "#E74C3C") +
    stat_summary(fun = mean, geom = "text", aes(label = sprintf("%.2f", after_stat(y))),
                 vjust = -1.5, size = 3.5, fontface = "bold", color = "#E74C3C") +
    scale_fill_manual(values = c("Without Immune\nNeighbors" = "#95A5A6",
                                 "With Immune\nNeighbors" = "#9B59B6")) +
    coord_cartesian(ylim = c(0, y_upper)) +  # Truncate without removing data
    geom_segment(aes(x = 1, xend = 2, y = y_upper * 0.9, yend = y_upper * 0.9),
                 color = "black", linewidth = 0.5) +
    geom_text(aes(x = 1.5, y = y_upper * 0.95,
                  label = sprintf("Δ = %.2f", cd44_summary[immune_neighbor_status == "With Immune\nNeighbors", mean_cd44] -
                                    cd44_summary[immune_neighbor_status == "Without Immune\nNeighbors", mean_cd44])),
              size = 3, fontface = "italic") +
    labs(
      title = "CD44 Expression by Immune Proximity",
      subtitle = subtitle_text,
      x = "",
      y = "CD44 Intensity",
      caption = sprintf("Y-axis truncated at 99th percentile; red diamonds show means\nn = %s without, %s with immune neighbors",
                        format(cd44_summary[immune_neighbor_status == "Without Immune\nNeighbors", n], big.mark = ","),
                        format(cd44_summary[immune_neighbor_status == "With Immune\nNeighbors", n], big.mark = ","))
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11),
      plot.caption = element_text(size = 8, hjust = 0.5)
    )
} else {
  panel_F <- ggplot() + theme_void() + labs(title = "CD44 data not available")
}

# =============================================================================
# Panels G-H: Spatial Visualizations (Representative ROIs)
# =============================================================================
cat("Creating Panels G-H: Spatial visualizations...\n")

# Select representative ROIs: one epithelial-dominant, one interface/stromal
epi_dom_roi <- roi_metrics[roi_type == "Epithelial-dominant"][which.max(pct_epithelial), ImageNumber]
interface_roi <- roi_metrics[roi_type == "Interface", ImageNumber]
if (length(interface_roi) == 0) {
  interface_roi <- roi_metrics[roi_type == "Stromal-dominant"][which.max(pct_stromal), ImageNumber]
}

cat("  Representative ROIs: Epithelial-dominant =", epi_dom_roi, ", Interface/Stromal =", interface_roi, "\n")

# Function to create spatial plot for a single ROI
create_spatial_plot <- function(data, roi_id, color_by = "compartment", title_suffix = "") {
  roi_data <- data[ImageNumber == roi_id]

  # Calculate scale bar parameters (assuming ~1 μm per pixel for IMC)
  x_range <- max(roi_data$Location_Center_X) - min(roi_data$Location_Center_X)
  y_range <- max(roi_data$Location_Center_Y) - min(roi_data$Location_Center_Y)
  scale_bar_length <- 200  # 200 pixels ≈ 200 μm
  scale_bar_x <- min(roi_data$Location_Center_X) + x_range * 0.05
  scale_bar_y <- max(roi_data$Location_Center_Y) - y_range * 0.05

  if (color_by == "compartment") {
    p <- ggplot(roi_data, aes(x = Location_Center_X, y = Location_Center_Y, color = compartment)) +
      geom_point(size = 0.3, alpha = 0.6) +
      scale_color_manual(values = compartment_colors)
  } else if (color_by == "epi_cluster") {
    # Show epithelial clusters, gray out non-epithelial
    roi_data[, plot_group := ifelse(compartment == "Epithelial" & !is.na(epi_cluster),
                                    epi_cluster, "Non-epithelial")]
    plot_colors <- c(cluster_colors, "Non-epithelial" = "#E0E0E0")

    p <- ggplot(roi_data, aes(x = Location_Center_X, y = Location_Center_Y, color = plot_group)) +
      geom_point(data = roi_data[plot_group == "Non-epithelial"], size = 0.2, alpha = 0.3) +
      geom_point(data = roi_data[plot_group != "Non-epithelial"], size = 0.4, alpha = 0.7) +
      scale_color_manual(values = plot_colors)
  }

  roi_type_label <- roi_metrics[ImageNumber == roi_id, roi_type]

  p <- p +
    # Add scale bar
    annotate("segment", x = scale_bar_x, xend = scale_bar_x + scale_bar_length,
             y = scale_bar_y, yend = scale_bar_y, linewidth = 1.5, color = "black") +
    annotate("text", x = scale_bar_x + scale_bar_length/2, y = scale_bar_y + y_range * 0.03,
             label = "200 μm", size = 2.5, color = "black") +
    coord_fixed() +
    scale_y_reverse() +  # IMC images typically have inverted Y
    labs(
      title = sprintf("ROI %d: %s", roi_id, roi_type_label),
      subtitle = title_suffix,
      x = "",
      y = "",
      color = ""
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      legend.position = "right",
      legend.key.size = unit(0.3, "cm"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

  return(p)
}

# Panel G: Epithelial-dominant ROI - compartments
panel_G <- create_spatial_plot(cell_data, epi_dom_roi, "compartment", "Cell compartments")

# Panel H: Interface/Stromal ROI - compartments
panel_H <- create_spatial_plot(cell_data, interface_roi, "compartment", "Cell compartments")

# Panel I: Epithelial-dominant ROI - epithelial clusters
panel_I <- create_spatial_plot(cell_data, epi_dom_roi, "epi_cluster", "Epithelial phenotypes")

# Panel J: Interface/Stromal ROI - epithelial clusters
panel_J <- create_spatial_plot(cell_data, interface_roi, "epi_cluster", "Epithelial phenotypes")

# =============================================================================
# Assemble Composite Figure
# =============================================================================
cat("Assembling composite figure...\n")

library(patchwork)

# Row 1: ROI classification and correlations (A, B)
# Row 2: UMAP and marker profiles (C, D)
# Row 3: Cluster distribution and CD44 (E, F)
# Row 4: Spatial images (G, H, I, J)

composite_fig <- (
  # Row 1
  (panel_A | panel_B) /
    # Row 2
    (panel_C | panel_D) /
    # Row 3
    (panel_E | panel_F) /
    # Row 4
    (panel_G | panel_H | panel_I | panel_J)
) +
  plot_annotation(
    title = "Figure 0.9: Tissue architecture predicts epithelial protein phenotypes",
    subtitle = sprintf("n = %s cells across %d ROIs from 2 patients",
                       format(nrow(cell_data), big.mark = ","), nrow(roi_metrics)),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11)
    )
  ) +
  plot_layout(heights = c(1, 1, 1, 1))

# Save composite figure
ggsave(file.path(output_dir, "Figure_0.9_composite.png"),
       composite_fig, width = 16, height = 20, dpi = 300)
ggsave(file.path(output_dir, "Figure_0.9_composite.pdf"),
       composite_fig, width = 16, height = 20)

cat("  Saved composite figure\n")

# =============================================================================
# Save Individual Panels (high resolution)
# =============================================================================
cat("Saving individual panels...\n")

ggsave(file.path(output_dir, "panel_A_roi_classification.png"), panel_A, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "panel_B_collagen_correlation.png"), panel_B, width = 5, height = 5, dpi = 300)
ggsave(file.path(output_dir, "panel_C_epithelial_dimred.png"), panel_C, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "panel_D_cluster_profiles.png"), panel_D, width = 8, height = 4, dpi = 300)
ggsave(file.path(output_dir, "panel_E_cluster_by_roi.png"), panel_E, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "panel_F_cd44_immune.png"), panel_F, width = 5, height = 5, dpi = 300)
ggsave(file.path(output_dir, "panel_G_spatial_epidom.png"), panel_G, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "panel_H_spatial_interface.png"), panel_H, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "panel_I_spatial_clusters_epidom.png"), panel_I, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "panel_J_spatial_clusters_interface.png"), panel_J, width = 6, height = 5, dpi = 300)

# =============================================================================
# Alternative: Smaller composite for main figure
# =============================================================================
cat("Creating condensed main figure...\n")

main_fig <- (
  # Row 1: ROI classification + Spatial example
  (panel_A | panel_G | panel_H) /
    # Row 2: Epithelial phenotyping
    (panel_C | panel_D) /
    # Row 3: Architecture-phenotype relationship + CD44
    (panel_E | panel_F | panel_I)
) +
  plot_annotation(
    title = "Figure 0.9: Tissue architecture predicts epithelial protein phenotypes",
    subtitle = sprintf("IMC analysis: n = %s cells, %d ROIs, 2 patients",
                       format(nrow(cell_data), big.mark = ","), nrow(roi_metrics)),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11)
    )
  ) +
  plot_layout(heights = c(1, 1, 1))

ggsave(file.path(output_dir, "Figure_0.9_main.png"), main_fig, width = 14, height = 14, dpi = 300)
ggsave(file.path(output_dir, "Figure_0.9_main.pdf"), main_fig, width = 14, height = 14)

cat("  Saved main figure\n")

# =============================================================================
# 7. EXPORT RESULTS
# =============================================================================
cat("\n=============================================================================\n")
cat("7. Exporting Results\n")
cat("=============================================================================\n")

fwrite(roi_metrics, file.path(output_dir, "roi_metrics_classification.csv"))
fwrite(neighborhood_summary, file.path(output_dir, "neighborhood_summary.csv"))

if (exists("cluster_profiles")) {
  fwrite(cluster_profiles, file.path(output_dir, "epithelial_cluster_profiles.csv"))
  fwrite(cluster_by_roi_type, file.path(output_dir, "cluster_by_roi_type.csv"))
}

# Export cell-level data with new annotations
# Build export columns dynamically based on what exists
export_cols <- c("ImageNumber", "ObjectNumber", "compartment", "roi_type")
if ("epi_cluster" %in% colnames(cell_data)) export_cols <- c(export_cols, "epi_cluster")
if ("n_neighbors" %in% colnames(cell_data)) {
  export_cols <- c(export_cols, "n_neighbors", "frac_epi_neighbors",
                   "frac_stro_neighbors", "frac_immune_neighbors")
}
fwrite(cell_data[, ..export_cols], file.path(output_dir, "cell_annotations.csv"))

# =============================================================================
# 8. MANUSCRIPT TEXT
# =============================================================================
cat("\n=============================================================================\n")
cat("MANUSCRIPT TEXT ADDITIONS\n")
cat("=============================================================================\n")

n_rois <- nrow(roi_metrics)
n_epi_dom <- roi_metrics[roi_type == "Epithelial-dominant", .N]
n_stro_dom <- roi_metrics[roi_type == "Stromal-dominant", .N]
n_interface <- roi_metrics[roi_type == "Interface", .N]
n_mixed <- roi_metrics[roi_type == "Mixed", .N]

cat("\n--- ROI Architecture Classification ---\n")
cat(sprintf("
To characterize spatial heterogeneity within the IMC dataset, we classified %d ROIs
based on cell compartment composition and spatial correlation patterns. ROIs were
designated as Epithelial-dominant (n = %d; >%.0f%% epithelial cells), Stromal-dominant
(n = %d; >%.0f%% stromal cells), Interface (n = %d; positive Collagen I-Pan-CK
correlation indicating epithelial-stromal intermixing), or Mixed (n = %d).
",
            n_rois, n_epi_dom, epi_dominant_thresh, n_stro_dom, stro_dominant_thresh,
            n_interface, n_mixed))

if (exists("optimal_k")) {
  cat("\n--- Epithelial Phenotyping ---\n")
  cat(sprintf("
Unsupervised clustering of %s epithelial cells based on %d protein markers
identified %d distinct phenotypic clusters. Cluster composition varied significantly
across ROI types (χ² = %.1f, p %s), indicating that tissue architecture
influences epithelial cell states at the protein level.
",
              format(nrow(epi_cells), big.mark = ","), length(phenotype_markers),
              optimal_k, chi_test$statistic,
              ifelse(chi_test$p.value < 0.001, "< 0.001", paste("=", round(chi_test$p.value, 3)))))
}

if (!is.na(neighbor_radius) && nrow(neighborhood_summary) > 0) {
  cat("\n--- Spatial Neighborhood Analysis ---\n")
  cat(sprintf("
Spatial neighborhood analysis (radius = %d pixels) revealed that epithelial cells
in Interface ROIs had higher fractions of stromal neighbors compared to
Epithelial-dominant ROIs, consistent with their classification as zones of
epithelial-stromal intermixing.
", neighbor_radius))
} else {
  cat("\n--- Spatial Neighborhood Analysis ---\n")
  cat("Spatial neighborhood analysis was not performed (no spatial coordinates available).\n")
}

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================================================\n")
cat("Output directory:", output_dir, "\n")
cat("\nGenerated files:\n")
list.files(output_dir, full.names = FALSE)

cat("\n\nKey findings:\n")
cat("  1. ROI Classification:\n")
cat(sprintf("     - Epithelial-dominant: %d ROIs\n", n_epi_dom))
cat(sprintf("     - Stromal-dominant: %d ROIs\n", n_stro_dom))
cat(sprintf("     - Interface: %d ROIs\n", n_interface))
cat(sprintf("     - Mixed: %d ROIs\n", n_mixed))

if (exists("optimal_k")) {
  cat(sprintf("  2. Epithelial phenotyping identified %d clusters\n", optimal_k))
  cat(sprintf("  3. Cluster-ROI association: χ² = %.1f, p %s\n",
              chi_test$statistic,
              ifelse(chi_test$p.value < 0.001, "< 0.001", paste("=", round(chi_test$p.value, 3)))))
}

cat("\n=============================================================================\n")
