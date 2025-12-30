# ==============================================================================
# R/trajectory/spatial_trajectory_analysis.R
# ==============================================================================
# Spatial-aware trajectory analysis for 10x Visium data
#
# This module extends standard trajectory analysis by incorporating
# spatial information from Visium spots:
#
#   1. Spatial autocorrelation of pseudotime (Moran's I)
#   2. Spatial smoothing/denoising of pseudotime
#   3. Spatial gradient analysis (direction of differentiation)
#   4. Niche-specific trajectory analysis
#   5. Spatial visualization overlays
#   6. Whorl structure analysis (ACP-specific)
#
# Requires: Seurat (with spatial data), spdep (optional)
# ==============================================================================

#' Run Spatial Trajectory Analysis
#'
#' Comprehensive spatial analysis of pseudotemporal trajectories
#'
#' @param seurat_obj Seurat object with spatial data and pseudotime
#' @param config Configuration list
#' @param pseudotime_col Column containing pseudotime values
#' @param subtype_col Column containing cell type classifications
#' @param k_neighbors Number of neighbors for spatial analyses (default: 6)
#' @param run_smoothing Whether to run spatial smoothing
#' @param run_gradients Whether to compute spatial gradients
#' @param run_niche_analysis Whether to run niche-specific analysis
#' @return List with spatial trajectory analysis results
#' @export
run_spatial_trajectory_analysis <- function(seurat_obj,
                                            config = NULL,
                                            pseudotime_col = "pseudotime_consensus",
                                            subtype_col = "module_score_subtype",
                                            k_neighbors = 6,
                                            run_smoothing = TRUE,
                                            run_gradients = TRUE,
                                            run_niche_analysis = TRUE) {


  message("\n========================================")
  message("Spatial Trajectory Analysis (Visium)")
  message("========================================\n")

  # Validate spatial data
  if (!"Spatial" %in% names(seurat_obj@assays) &&
      length(seurat_obj@images) == 0) {
    stop("No spatial data found in Seurat object")
  }

  if (!pseudotime_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Pseudotime column '%s' not found", pseudotime_col))
  }

  results <- list()

  # Get spatial coordinates
  coords <- get_spatial_coordinates(seurat_obj)
  message(sprintf("Analyzing %d spots with spatial coordinates", nrow(coords)))

  # Get pseudotime
  pseudotime <- seurat_obj@meta.data[[pseudotime_col]]
  names(pseudotime) <- colnames(seurat_obj)

  # ===========================================================================
  # 1. Build Spatial Neighborhood Graph
  # ===========================================================================
  message("\n--- Building Spatial Neighborhood Graph ---")

  spatial_graph <- build_spatial_graph(coords, k = k_neighbors)
  results$spatial_graph <- spatial_graph

  message(sprintf("  K-nearest neighbors: %d", k_neighbors))
  message(sprintf("  Mean neighbor distance: %.2f", mean(spatial_graph$distances)))

  # ===========================================================================
  # 2. Spatial Autocorrelation of Pseudotime
  # ===========================================================================
  message("\n--- Spatial Autocorrelation Analysis ---")

  autocorr <- calculate_spatial_autocorrelation(
    values = pseudotime,
    spatial_graph = spatial_graph,
    coords = coords
  )
  results$autocorrelation <- autocorr

  message(sprintf("  Moran's I: %.4f (p = %.2e)",
                  autocorr$morans_i$observed,
                  autocorr$morans_i$p_value))
  message(sprintf("  Geary's C: %.4f (p = %.2e)",
                  autocorr$gearys_c$observed,
                  autocorr$gearys_c$p_value))
  message(sprintf("  Interpretation: %s", autocorr$interpretation))

  # ===========================================================================
  # 3. Spatial Smoothing of Pseudotime
  # ===========================================================================
  if (run_smoothing) {
    message("\n--- Spatial Smoothing ---")

    smoothed <- spatial_smooth_pseudotime(
      pseudotime = pseudotime,
      spatial_graph = spatial_graph,
      method = "weighted_mean",
      iterations = 1
    )

    results$smoothed_pseudotime <- smoothed$smoothed
    results$smoothing_delta <- smoothed$delta

    # Add to Seurat object
    seurat_obj$pseudotime_spatial_smoothed <- smoothed$smoothed
    seurat_obj$pseudotime_smoothing_delta <- smoothed$delta

    message(sprintf("  Mean change after smoothing: %.4f", mean(abs(smoothed$delta))))
    message(sprintf("  Correlation with original: %.4f",
                    cor(pseudotime, smoothed$smoothed, use = "complete.obs")))
  }

  # ===========================================================================
  # 4. Spatial Gradient Analysis
  # ===========================================================================
  if (run_gradients) {
    message("\n--- Spatial Gradient Analysis ---")

    gradients <- calculate_spatial_gradients(
      pseudotime = pseudotime,
      coords = coords,
      spatial_graph = spatial_graph
    )

    results$gradients <- gradients

    # Add to Seurat object
    seurat_obj$pseudotime_gradient_magnitude <- gradients$magnitude
    seurat_obj$pseudotime_gradient_angle <- gradients$angle

    message(sprintf("  Mean gradient magnitude: %.4f", mean(gradients$magnitude, na.rm = TRUE)))
    message(sprintf("  Dominant direction: %.1f degrees", gradients$dominant_direction))

    # Test for global directionality
    if (gradients$rayleigh_test$p_value < 0.05) {
      message(sprintf("  Significant directional bias detected (p = %.4f)",
                      gradients$rayleigh_test$p_value))
    } else {
      message("  No significant directional bias")
    }
  }

  # ===========================================================================
  # 5. Per-Subtype Spatial Analysis
  # ===========================================================================
  message("\n--- Per-Subtype Spatial Analysis ---")

  subtypes <- seurat_obj@meta.data[[subtype_col]]
  subtype_spatial <- analyze_subtype_spatial_patterns(
    subtypes = subtypes,
    pseudotime = pseudotime,
    coords = coords,
    spatial_graph = spatial_graph
  )

  results$subtype_spatial <- subtype_spatial

  message("  Subtype clustering (Moran's I on indicator):")
  for (st in names(subtype_spatial$clustering)) {
    mi <- subtype_spatial$clustering[[st]]$morans_i
    message(sprintf("    %s: I = %.3f (p = %.3e)", st, mi$observed, mi$p_value))
  }

  # ===========================================================================
  # 6. Niche-Specific Trajectory Analysis
  # ===========================================================================
  if (run_niche_analysis) {
    message("\n--- Niche-Specific Trajectory Analysis ---")

    # Identify spatial niches/domains
    niches <- identify_spatial_niches(
      seurat_obj = seurat_obj,
      coords = coords,
      spatial_graph = spatial_graph,
      n_niches = 4  # Default: 4 spatial domains
    )

    results$niches <- niches
    seurat_obj$spatial_niche <- niches$assignments

    message(sprintf("  Identified %d spatial niches", niches$n_niches))

    # Analyze trajectory within each niche
    niche_trajectories <- analyze_niche_trajectories(
      pseudotime = pseudotime,
      subtypes = subtypes,
      niches = niches$assignments,
      coords = coords
    )

    results$niche_trajectories <- niche_trajectories

    message("  Per-niche trajectory characteristics:")
    for (niche in names(niche_trajectories$per_niche)) {
      nt <- niche_trajectories$per_niche[[niche]]
      message(sprintf("    %s: n=%d, mean_pt=%.3f, dominant=%s",
                      niche, nt$n_spots, nt$mean_pseudotime, nt$dominant_subtype))
    }
  }

  # ===========================================================================
  # 7. Spatial Transition Zones
  # ===========================================================================
  message("\n--- Identifying Spatial Transition Zones ---")

  transitions <- identify_transition_zones(
    pseudotime = pseudotime,
    subtypes = subtypes,
    spatial_graph = spatial_graph,
    coords = coords
  )

  results$transition_zones <- transitions
  seurat_obj$is_transition_zone <- transitions$is_transition
  seurat_obj$transition_score <- transitions$scores

  message(sprintf("  Transition zone spots: %d (%.1f%%)",
                  sum(transitions$is_transition),
                  mean(transitions$is_transition) * 100))

  # ===========================================================================
  # 8. Sample-Specific Spatial Analysis (if multiple samples)
  # ===========================================================================
  if ("sample" %in% colnames(seurat_obj@meta.data)) {
    message("\n--- Per-Sample Spatial Analysis ---")

    samples <- unique(seurat_obj$sample)
    sample_spatial <- list()

    for (samp in samples) {
      cells_in_sample <- which(seurat_obj$sample == samp)

      if (length(cells_in_sample) < 50) {
        message(sprintf("  %s: Skipped (n=%d < 50)", samp, length(cells_in_sample)))
        next
      }

      # Subset coordinates and pseudotime
      samp_coords <- coords[cells_in_sample, ]
      samp_pt <- pseudotime[cells_in_sample]

      # Build sample-specific graph
      samp_graph <- build_spatial_graph(samp_coords, k = k_neighbors)

      # Calculate autocorrelation
      samp_autocorr <- calculate_spatial_autocorrelation(
        values = samp_pt,
        spatial_graph = samp_graph,
        coords = samp_coords
      )

      sample_spatial[[samp]] <- list(
        n_spots = length(cells_in_sample),
        morans_i = samp_autocorr$morans_i$observed,
        morans_p = samp_autocorr$morans_i$p_value,
        interpretation = samp_autocorr$interpretation
      )

      message(sprintf("  %s: n=%d, Moran's I=%.3f (p=%.3e)",
                      samp, length(cells_in_sample),
                      samp_autocorr$morans_i$observed,
                      samp_autocorr$morans_i$p_value))
    }

    results$sample_spatial <- sample_spatial
  }

  # ===========================================================================
  # 9. Summary Statistics
  # ===========================================================================
  message("\n--- Summary ---")

  results$summary <- list(
    n_spots = nrow(coords),
    k_neighbors = k_neighbors,
    pseudotime_spatially_organized = autocorr$morans_i$p_value < 0.05 &&
      autocorr$morans_i$observed > 0,
    spatial_autocorrelation = autocorr$morans_i$observed,
    has_directional_gradient = if (run_gradients) gradients$rayleigh_test$p_value < 0.05 else NA,
    n_transition_spots = sum(transitions$is_transition),
    pct_transition_spots = mean(transitions$is_transition) * 100
  )

  if (results$summary$pseudotime_spatially_organized) {
    message("  ✓ Pseudotime shows significant spatial organization")
  } else {
    message("  ✗ Pseudotime does not show significant spatial organization")
  }

  # Store updated seurat object
  results$seurat_obj <- seurat_obj

  message("\n========================================\n")

  return(results)
}


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Get Spatial Coordinates from Seurat Object
#'
#' @param seurat_obj Seurat object with spatial data
#' @return Data frame with x, y coordinates
#' @keywords internal
get_spatial_coordinates <- function(seurat_obj) {

  # Try different methods to get coordinates
  coords <- NULL

  # Method 1: GetTissueCoordinates (Seurat v5)
  if (is.null(coords)) {
    coords <- tryCatch({
      Seurat::GetTissueCoordinates(seurat_obj)
    }, error = function(e) NULL)
  }

  # Method 2: From images slot
  if (is.null(coords) && length(seurat_obj@images) > 0) {
    coords <- tryCatch({
      img_name <- names(seurat_obj@images)[1]
      Seurat::GetTissueCoordinates(seurat_obj@images[[img_name]])
    }, error = function(e) NULL)
  }

  # Method 3: Direct access to coordinates
  if (is.null(coords) && length(seurat_obj@images) > 0) {
    coords <- tryCatch({
      img_name <- names(seurat_obj@images)[1]
      coord_df <- seurat_obj@images[[img_name]]@coordinates
      data.frame(
        x = coord_df$imagecol,
        y = coord_df$imagerow,
        row.names = rownames(coord_df)
      )
    }, error = function(e) NULL)
  }

  if (is.null(coords)) {
    stop("Could not extract spatial coordinates from Seurat object")
  }

  # Ensure column names are standardized
  if (!"x" %in% colnames(coords)) {
    if ("imagerow" %in% colnames(coords)) {
      coords$x <- coords$imagecol
      coords$y <- coords$imagerow
    } else {
      colnames(coords)[1:2] <- c("x", "y")
    }
  }

  # Filter to cells in the object
  common_cells <- intersect(rownames(coords), colnames(seurat_obj))
  coords <- coords[common_cells, c("x", "y")]

  return(coords)
}


#' Build Spatial Neighborhood Graph
#'
#' @param coords Data frame with x, y coordinates
#' @param k Number of nearest neighbors
#' @param max_dist Maximum distance for neighbors (optional)
#' @return List with adjacency matrix and distances
#' @keywords internal
build_spatial_graph <- function(coords, k = 6, max_dist = NULL) {

  n <- nrow(coords)

  # Calculate pairwise distances
  dist_matrix <- as.matrix(dist(coords[, c("x", "y")]))

  # Build k-NN adjacency
  adj_matrix <- matrix(0, n, n)
  rownames(adj_matrix) <- colnames(adj_matrix) <- rownames(coords)

  neighbor_distances <- c()

  for (i in 1:n) {
    # Get k nearest neighbors (excluding self)
    dists <- dist_matrix[i, ]
    dists[i] <- Inf  # Exclude self

    nn_idx <- order(dists)[1:k]
    nn_dists <- dists[nn_idx]

    # Apply max distance filter if specified
    if (!is.null(max_dist)) {
      nn_idx <- nn_idx[nn_dists <= max_dist]
      nn_dists <- nn_dists[nn_dists <= max_dist]
    }

    adj_matrix[i, nn_idx] <- 1
    neighbor_distances <- c(neighbor_distances, nn_dists)
  }

  # Symmetrize (mutual neighbors)
  adj_matrix <- pmax(adj_matrix, t(adj_matrix))

  # Create neighbor list (for spdep compatibility)
  neighbor_list <- lapply(1:n, function(i) which(adj_matrix[i, ] > 0))
  names(neighbor_list) <- rownames(coords)

  return(list(
    adjacency = adj_matrix,
    distances = neighbor_distances,
    neighbor_list = neighbor_list,
    k = k,
    n = n
  ))
}


#' Calculate Spatial Autocorrelation
#'
#' Computes Moran's I and Geary's C for spatial autocorrelation
#'
#' @param values Numeric vector of values to test
#' @param spatial_graph Spatial graph from build_spatial_graph
#' @param coords Spatial coordinates
#' @param n_permutations Number of permutations for p-value
#' @return List with autocorrelation statistics
#' @keywords internal
calculate_spatial_autocorrelation <- function(values,
                                              spatial_graph,
                                              coords,
                                              n_permutations = 999) {

  n <- length(values)
  W <- spatial_graph$adjacency

  # Normalize weights (row-standardized)
  row_sums <- rowSums(W)
  row_sums[row_sums == 0] <- 1  # Avoid division by zero
  W_norm <- W / row_sums

  # Remove NAs
  valid_idx <- !is.na(values)
  values_valid <- values[valid_idx]
  W_valid <- W_norm[valid_idx, valid_idx]
  n_valid <- sum(valid_idx)

  # ===========================================================================
  # Moran's I
  # ===========================================================================

  # Center values
  y <- values_valid - mean(values_valid)

  # Calculate Moran's I
  numerator <- sum(W_valid * outer(y, y))
  denominator <- sum(y^2)
  S0 <- sum(W_valid)

  morans_i <- (n_valid / S0) * (numerator / denominator)

  # Permutation test for Moran's I
  perm_i <- replicate(n_permutations, {
    y_perm <- sample(y)
    num_perm <- sum(W_valid * outer(y_perm, y_perm))
    (n_valid / S0) * (num_perm / sum(y_perm^2))
  })

  # P-value (two-tailed)
  p_morans <- (sum(abs(perm_i) >= abs(morans_i)) + 1) / (n_permutations + 1)

  # Expected value under null
  E_I <- -1 / (n_valid - 1)

  # ===========================================================================
  # Geary's C
  # ===========================================================================

  # Calculate Geary's C
  diff_sq <- outer(values_valid, values_valid, FUN = function(x, y) (x - y)^2)
  gearys_c <- ((n_valid - 1) * sum(W_valid * diff_sq)) / (2 * S0 * sum(y^2))

  # Permutation test for Geary's C
  perm_c <- replicate(n_permutations, {
    v_perm <- sample(values_valid)
    y_perm <- v_perm - mean(v_perm)
    diff_sq_perm <- outer(v_perm, v_perm, FUN = function(x, y) (x - y)^2)
    ((n_valid - 1) * sum(W_valid * diff_sq_perm)) / (2 * S0 * sum(y_perm^2))
  })

  # P-value for Geary's C (C < 1 indicates positive autocorrelation)
  p_gearys <- (sum(perm_c <= gearys_c) + 1) / (n_permutations + 1)

  # ===========================================================================
  # Interpretation
  # ===========================================================================

  if (p_morans < 0.05) {
    if (morans_i > E_I) {
      interpretation <- "Significant POSITIVE spatial autocorrelation (similar values cluster)"
    } else {
      interpretation <- "Significant NEGATIVE spatial autocorrelation (dissimilar values cluster)"
    }
  } else {
    interpretation <- "No significant spatial autocorrelation (random spatial distribution)"
  }

  return(list(
    morans_i = list(
      observed = morans_i,
      expected = E_I,
      p_value = p_morans,
      permutation_dist = perm_i
    ),
    gearys_c = list(
      observed = gearys_c,
      expected = 1,  # Expected under null
      p_value = p_gearys,
      permutation_dist = perm_c
    ),
    interpretation = interpretation,
    n_valid = n_valid
  ))
}


#' Spatial Smoothing of Pseudotime
#'
#' Denoises pseudotime using spatial neighbors
#'
#' @param pseudotime Pseudotime vector
#' @param spatial_graph Spatial graph
#' @param method Smoothing method: "weighted_mean", "median", "gaussian"
#' @param iterations Number of smoothing iterations
#' @param sigma Bandwidth for gaussian smoothing
#' @return List with smoothed pseudotime and delta
#' @keywords internal
spatial_smooth_pseudotime <- function(pseudotime,
                                      spatial_graph,
                                      method = "weighted_mean",
                                      iterations = 1,
                                      sigma = NULL) {

  n <- length(pseudotime)
  smoothed <- pseudotime

  # Get neighbor list and distances
  neighbors <- spatial_graph$neighbor_list
  adj <- spatial_graph$adjacency

  # Calculate distance-based weights if needed
  if (method == "gaussian" || method == "weighted_mean") {
    # Use adjacency for now; could use actual distances
    weights <- adj
    row_sums <- rowSums(weights)
    row_sums[row_sums == 0] <- 1
    weights <- weights / row_sums
  }

  for (iter in 1:iterations) {
    new_smoothed <- smoothed

    for (i in 1:n) {
      nn_idx <- neighbors[[i]]

      if (length(nn_idx) == 0) {
        next  # No neighbors, keep original
      }

      nn_values <- smoothed[nn_idx]

      # Remove NAs
      valid <- !is.na(nn_values)
      if (sum(valid) == 0) next

      nn_values <- nn_values[valid]

      if (method == "median") {
        # Include self in median
        new_smoothed[i] <- median(c(smoothed[i], nn_values), na.rm = TRUE)

      } else if (method == "weighted_mean") {
        # Weighted by adjacency (could use distance-weighted)
        w <- weights[i, nn_idx[valid]]
        w <- c(1, w)  # Weight 1 for self
        w <- w / sum(w)
        new_smoothed[i] <- sum(w * c(smoothed[i], nn_values))

      } else if (method == "gaussian") {
        # Gaussian kernel smoothing
        if (is.null(sigma)) sigma <- mean(spatial_graph$distances)
        # Simplified: use uniform weights for neighbors
        new_smoothed[i] <- mean(c(smoothed[i], nn_values))
      }
    }

    smoothed <- new_smoothed
  }

  # Calculate change from original
  delta <- smoothed - pseudotime

  return(list(
    smoothed = smoothed,
    delta = delta,
    method = method,
    iterations = iterations
  ))
}


#' Calculate Spatial Gradients
#'
#' Estimates local gradients of pseudotime in spatial coordinates
#'
#' @param pseudotime Pseudotime vector
#' @param coords Spatial coordinates
#' @param spatial_graph Spatial graph
#' @return List with gradient magnitudes, angles, and statistics
#' @keywords internal
calculate_spatial_gradients <- function(pseudotime, coords, spatial_graph) {

  n <- length(pseudotime)
  neighbors <- spatial_graph$neighbor_list

  gradient_x <- rep(NA, n)
  gradient_y <- rep(NA, n)

  for (i in 1:n) {
    nn_idx <- neighbors[[i]]

    if (length(nn_idx) < 2 || is.na(pseudotime[i])) next

    # Get neighbor coordinates and pseudotime
    nn_coords <- coords[nn_idx, ]
    nn_pt <- pseudotime[nn_idx]

    # Remove NAs
    valid <- !is.na(nn_pt)
    if (sum(valid) < 2) next

    nn_coords <- nn_coords[valid, ]
    nn_pt <- nn_pt[valid]

    # Calculate differences from focal spot
    dx <- nn_coords$x - coords$x[i]
    dy <- nn_coords$y - coords$y[i]
    dpt <- nn_pt - pseudotime[i]

    # Fit plane: dpt = gx*dx + gy*dy
    # Using least squares
    X <- cbind(dx, dy)

    tryCatch({
      fit <- lm.fit(X, dpt)
      gradient_x[i] <- fit$coefficients[1]
      gradient_y[i] <- fit$coefficients[2]
    }, error = function(e) {
      # Linear fit failed, use simple difference
      gradient_x[i] <- mean(dpt / dx, na.rm = TRUE)
      gradient_y[i] <- mean(dpt / dy, na.rm = TRUE)
    })
  }

  # Calculate magnitude and angle
  magnitude <- sqrt(gradient_x^2 + gradient_y^2)
  angle <- atan2(gradient_y, gradient_x) * 180 / pi  # In degrees

  # Rayleigh test for directional uniformity
  # Tests whether angles have a preferred direction
  valid_angles <- angle[!is.na(angle)]
  valid_angles_rad <- valid_angles * pi / 180

  # Mean resultant length (R-bar)
  mean_cos <- mean(cos(valid_angles_rad))
  mean_sin <- mean(sin(valid_angles_rad))
  R_bar <- sqrt(mean_cos^2 + mean_sin^2)

  # Rayleigh test statistic
  n_valid <- length(valid_angles)
  z <- n_valid * R_bar^2

  # P-value (approximation)
  p_rayleigh <- exp(-z) * (1 + (2*z - z^2) / (4*n_valid) -
                             (24*z - 132*z^2 + 76*z^3 - 9*z^4) / (288*n_valid^2))

  # Dominant direction
  dominant_direction <- atan2(mean_sin, mean_cos) * 180 / pi
  if (dominant_direction < 0) dominant_direction <- dominant_direction + 360

  return(list(
    gradient_x = gradient_x,
    gradient_y = gradient_y,
    magnitude = magnitude,
    angle = angle,
    dominant_direction = dominant_direction,
    mean_resultant_length = R_bar,
    rayleigh_test = list(
      statistic = z,
      p_value = p_rayleigh,
      significant = p_rayleigh < 0.05
    )
  ))
}


#' Analyze Subtype Spatial Patterns
#'
#' Tests spatial clustering of each cell subtype
#'
#' @param subtypes Subtype labels
#' @param pseudotime Pseudotime vector
#' @param coords Spatial coordinates
#' @param spatial_graph Spatial graph
#' @return List with per-subtype spatial statistics
#' @keywords internal
analyze_subtype_spatial_patterns <- function(subtypes, pseudotime, coords, spatial_graph) {

  unique_subtypes <- unique(subtypes[!is.na(subtypes)])

  clustering <- list()
  spatial_stats <- list()

  for (st in unique_subtypes) {
    # Create binary indicator
    indicator <- as.numeric(subtypes == st)
    indicator[is.na(subtypes)] <- NA

    # Test spatial clustering of this subtype
    autocorr <- calculate_spatial_autocorrelation(
      values = indicator,
      spatial_graph = spatial_graph,
      coords = coords,
      n_permutations = 499
    )

    clustering[[st]] <- list(
      morans_i = autocorr$morans_i,
      n_spots = sum(indicator == 1, na.rm = TRUE),
      clustered = autocorr$morans_i$p_value < 0.05 && autocorr$morans_i$observed > 0
    )

    # Spatial statistics for this subtype
    st_idx <- which(subtypes == st)
    if (length(st_idx) >= 3) {
      st_coords <- coords[st_idx, ]

      # Centroid
      centroid <- colMeans(st_coords)

      # Dispersion (mean distance from centroid)
      dists_from_centroid <- sqrt((st_coords$x - centroid[1])^2 +
                                    (st_coords$y - centroid[2])^2)

      spatial_stats[[st]] <- list(
        centroid = centroid,
        mean_dispersion = mean(dists_from_centroid),
        mean_pseudotime = mean(pseudotime[st_idx], na.rm = TRUE)
      )
    }
  }

  return(list(
    clustering = clustering,
    spatial_stats = spatial_stats
  ))
}


#' Identify Spatial Niches
#'
#' Clusters spots into spatial domains/niches
#'
#' @param seurat_obj Seurat object
#' @param coords Spatial coordinates
#' @param spatial_graph Spatial graph
#' @param n_niches Number of niches to identify
#' @param method Clustering method: "kmeans", "leiden", "expression"
#' @return List with niche assignments
#' @keywords internal
identify_spatial_niches <- function(seurat_obj,
                                    coords,
                                    spatial_graph,
                                    n_niches = 4,
                                    method = "kmeans") {

  if (method == "kmeans") {
    # Simple k-means on spatial coordinates
    km <- kmeans(coords[, c("x", "y")], centers = n_niches, nstart = 25)
    assignments <- paste0("Niche_", km$cluster)
    names(assignments) <- rownames(coords)

  } else if (method == "expression") {
    # Use Seurat clustering on spatially-variable genes
    # This requires FindSpatiallyVariableFeatures to have been run
    if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      assignments <- paste0("Niche_", seurat_obj$seurat_clusters)
      names(assignments) <- colnames(seurat_obj)
    } else {
      # Fallback to kmeans
      km <- kmeans(coords[, c("x", "y")], centers = n_niches, nstart = 25)
      assignments <- paste0("Niche_", km$cluster)
      names(assignments) <- rownames(coords)
    }

  } else if (method == "leiden") {
    # Leiden clustering on spatial graph
    if (requireNamespace("igraph", quietly = TRUE)) {
      g <- igraph::graph_from_adjacency_matrix(spatial_graph$adjacency, mode = "undirected")

      # Use Louvain as fallback if leiden not available
      communities <- tryCatch({
        igraph::cluster_leiden(g, resolution_parameter = 0.5)
      }, error = function(e) {
        igraph::cluster_louvain(g)
      })

      assignments <- paste0("Niche_", igraph::membership(communities))
      names(assignments) <- rownames(coords)
    } else {
      # Fallback
      km <- kmeans(coords[, c("x", "y")], centers = n_niches, nstart = 25)
      assignments <- paste0("Niche_", km$cluster)
      names(assignments) <- rownames(coords)
    }
  }

  return(list(
    assignments = assignments,
    n_niches = length(unique(assignments)),
    method = method
  ))
}


#' Analyze Niche-Specific Trajectories
#'
#' Examines trajectory characteristics within each spatial niche
#'
#' @param pseudotime Pseudotime vector
#' @param subtypes Subtype labels
#' @param niches Niche assignments
#' @param coords Spatial coordinates
#' @return List with per-niche trajectory statistics
#' @keywords internal
analyze_niche_trajectories <- function(pseudotime, subtypes, niches, coords) {

  unique_niches <- unique(niches[!is.na(niches)])

  per_niche <- list()

  for (niche in unique_niches) {
    niche_idx <- which(niches == niche)

    if (length(niche_idx) < 10) next

    niche_pt <- pseudotime[niche_idx]
    niche_subtypes <- subtypes[niche_idx]
    niche_coords <- coords[niche_idx, ]

    # Basic statistics
    per_niche[[niche]] <- list(
      n_spots = length(niche_idx),
      mean_pseudotime = mean(niche_pt, na.rm = TRUE),
      sd_pseudotime = sd(niche_pt, na.rm = TRUE),
      pseudotime_range = diff(range(niche_pt, na.rm = TRUE)),
      subtype_composition = table(niche_subtypes) / length(niche_subtypes),
      dominant_subtype = names(which.max(table(niche_subtypes))),
      centroid = colMeans(niche_coords)
    )
  }

  # Test for differences between niches
  if (length(per_niche) >= 2) {
    # Kruskal-Wallis test
    kw_test <- kruskal.test(pseudotime ~ niches)

    # Pairwise comparisons
    pairwise <- tryCatch({
      pairwise.wilcox.test(pseudotime, niches, p.adjust.method = "BH")
    }, error = function(e) NULL)

  } else {
    kw_test <- NULL
    pairwise <- NULL
  }

  return(list(
    per_niche = per_niche,
    kruskal_wallis = kw_test,
    pairwise_tests = pairwise,
    niches_differ = !is.null(kw_test) && kw_test$p.value < 0.05
  ))
}


#' Identify Spatial Transition Zones
#'
#' Finds spots at boundaries between different subtypes/pseudotime states
#'
#' @param pseudotime Pseudotime vector
#' @param subtypes Subtype labels
#' @param spatial_graph Spatial graph
#' @param coords Spatial coordinates
#' @param threshold Percentile threshold for transition score
#' @return List with transition zone identification
#' @keywords internal
identify_transition_zones <- function(pseudotime,
                                      subtypes,
                                      spatial_graph,
                                      coords,
                                      threshold = 0.75) {

  n <- length(pseudotime)
  neighbors <- spatial_graph$neighbor_list

  transition_scores <- rep(0, n)

  for (i in 1:n) {
    nn_idx <- neighbors[[i]]

    if (length(nn_idx) == 0 || is.na(subtypes[i])) next

    # Score 1: Pseudotime heterogeneity among neighbors
    nn_pt <- pseudotime[nn_idx]
    pt_var <- var(c(pseudotime[i], nn_pt), na.rm = TRUE)

    # Score 2: Subtype heterogeneity among neighbors
    nn_subtypes <- subtypes[nn_idx]
    n_different <- sum(nn_subtypes != subtypes[i], na.rm = TRUE)
    subtype_hetero <- n_different / length(nn_idx)

    # Combined score (weighted)
    transition_scores[i] <- 0.5 * (pt_var / max(var(pseudotime, na.rm = TRUE), 0.001)) +
      0.5 * subtype_hetero
  }

  # Identify transition zones
  score_threshold <- quantile(transition_scores, threshold, na.rm = TRUE)
  is_transition <- transition_scores >= score_threshold

  # Identify specific transition types
  transition_types <- rep(NA, n)

  for (i in which(is_transition)) {
    nn_idx <- neighbors[[i]]
    nn_subtypes <- unique(subtypes[nn_idx])
    nn_subtypes <- nn_subtypes[!is.na(nn_subtypes)]

    if (length(nn_subtypes) > 1) {
      transition_types[i] <- paste(sort(c(subtypes[i], nn_subtypes)), collapse = "-")
    }
  }

  return(list(
    scores = transition_scores,
    is_transition = is_transition,
    transition_types = transition_types,
    threshold = score_threshold,
    n_transition = sum(is_transition),
    transition_type_counts = table(transition_types[is_transition])
  ))
}


# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

#' Plot Spatial Pseudotime
#'
#' Creates spatial plot with pseudotime overlay
#'
#' @param seurat_obj Seurat object with spatial data and pseudotime
#' @param pseudotime_col Pseudotime column name
#' @param pt.size Point size
#' @param alpha Transparency
#' @return ggplot object
#' @export
plot_spatial_pseudotime <- function(seurat_obj,
                                    pseudotime_col = "pseudotime_consensus",
                                    pt.size = 1.5,
                                    alpha = 0.8) {

  require(ggplot2)

  # Get coordinates
  coords <- get_spatial_coordinates(seurat_obj)

  plot_data <- data.frame(
    x = coords$x,
    y = coords$y,
    pseudotime = seurat_obj@meta.data[[pseudotime_col]]
  )

  p <- ggplot(plot_data, aes(x = x, y = -y, color = pseudotime)) +
    geom_point(size = pt.size, alpha = alpha) +
    scale_color_viridis_c(option = "plasma") +
    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    ) +
    labs(title = "Spatial Pseudotime", color = "Pseudotime")

  return(p)
}


#' Plot Spatial Gradients
#'
#' Visualizes pseudotime gradients as arrows on tissue
#'
#' @param seurat_obj Seurat object
#' @param spatial_results Results from run_spatial_trajectory_analysis
#' @param arrow_scale Scale factor for arrows
#' @param subsample Fraction of spots to show arrows for
#' @return ggplot object
#' @export
plot_spatial_gradients <- function(seurat_obj,
                                   spatial_results,
                                   arrow_scale = 50,
                                   subsample = 0.2) {

  require(ggplot2)

  coords <- get_spatial_coordinates(seurat_obj)
  gradients <- spatial_results$gradients

  plot_data <- data.frame(
    x = coords$x,
    y = coords$y,
    gx = gradients$gradient_x * arrow_scale,
    gy = gradients$gradient_y * arrow_scale,
    magnitude = gradients$magnitude,
    pseudotime = seurat_obj$pseudotime_consensus
  )

  # Remove NAs and subsample
  plot_data <- plot_data[!is.na(plot_data$gx), ]

  if (subsample < 1) {
    n_show <- round(nrow(plot_data) * subsample)
    idx <- sample(1:nrow(plot_data), n_show)
    arrow_data <- plot_data[idx, ]
  } else {
    arrow_data <- plot_data
  }

  p <- ggplot() +
    geom_point(data = plot_data, aes(x = x, y = -y, color = pseudotime),
               size = 1, alpha = 0.5) +
    geom_segment(data = arrow_data,
                 aes(x = x, y = -y, xend = x + gx, yend = -(y + gy)),
                 arrow = arrow(length = unit(0.1, "cm")),
                 color = "black", alpha = 0.7) +
    scale_color_viridis_c(option = "plasma") +
    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    ) +
    labs(title = sprintf("Spatial Gradients (Dominant: %.0f°)",
                         spatial_results$gradients$dominant_direction),
         subtitle = sprintf("Rayleigh test p = %.4f",
                            spatial_results$gradients$rayleigh_test$p_value),
         color = "Pseudotime")

  return(p)
}


#' Plot Spatial Transition Zones
#'
#' Highlights spots in transition zones
#'
#' @param seurat_obj Seurat object
#' @param spatial_results Results from run_spatial_trajectory_analysis
#' @return ggplot object
#' @export
plot_transition_zones <- function(seurat_obj, spatial_results) {

  require(ggplot2)

  coords <- get_spatial_coordinates(seurat_obj)
  transitions <- spatial_results$transition_zones

  plot_data <- data.frame(
    x = coords$x,
    y = coords$y,
    score = transitions$scores,
    is_transition = transitions$is_transition,
    subtype = seurat_obj$module_score_subtype
  )

  p <- ggplot(plot_data, aes(x = x, y = -y)) +
    geom_point(aes(color = score), size = 1.5, alpha = 0.7) +
    geom_point(data = plot_data[plot_data$is_transition, ],
               shape = 1, size = 3, color = "red", stroke = 0.5) +
    scale_color_viridis_c(option = "magma") +
    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    ) +
    labs(title = sprintf("Transition Zones (n=%d, %.1f%%)",
                         sum(transitions$is_transition),
                         mean(transitions$is_transition) * 100),
         subtitle = "Red circles = transition zone spots",
         color = "Transition\nScore")

  return(p)
}


#' Generate Spatial Trajectory Report
#'
#' Creates comprehensive spatial analysis plots
#'
#' @param seurat_obj Seurat object
#' @param spatial_results Results from run_spatial_trajectory_analysis
#' @param config Configuration list
#' @param output_dir Output directory for plots
#' @return List of ggplot objects
#' @export
generate_spatial_trajectory_plots <- function(seurat_obj,
                                              spatial_results,
                                              config = NULL,
                                              output_dir = NULL) {

  require(ggplot2)
  require(patchwork)

  plots <- list()

  # Get colors
  if (!is.null(config)) {
    colors <- unlist(config$visualization$colors$epithelial_subtypes)
  } else {
    colors <- NULL
  }

  # 1. Spatial pseudotime
  plots$pseudotime <- plot_spatial_pseudotime(seurat_obj)

  # 2. Spatial pseudotime (smoothed if available)
  if ("pseudotime_spatial_smoothed" %in% colnames(seurat_obj@meta.data)) {
    plots$pseudotime_smoothed <- plot_spatial_pseudotime(
      seurat_obj,
      pseudotime_col = "pseudotime_spatial_smoothed"
    ) + labs(title = "Spatially Smoothed Pseudotime")
  }

  # 3. Spatial gradients
  if (!is.null(spatial_results$gradients)) {
    plots$gradients <- plot_spatial_gradients(seurat_obj, spatial_results)
  }

  # 4. Transition zones
  if (!is.null(spatial_results$transition_zones)) {
    plots$transitions <- plot_transition_zones(seurat_obj, spatial_results)
  }

  # 5. Subtype spatial distribution
  coords <- get_spatial_coordinates(seurat_obj)
  subtype_data <- data.frame(
    x = coords$x,
    y = coords$y,
    subtype = seurat_obj$module_score_subtype
  )

  plots$subtypes <- ggplot(subtype_data, aes(x = x, y = -y, color = subtype)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values = colors) +
    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    ) +
    labs(title = "Spatial Subtype Distribution", color = "Subtype")

  # 6. Spatial niches (if available)
  if ("spatial_niche" %in% colnames(seurat_obj@meta.data)) {
    niche_data <- data.frame(
      x = coords$x,
      y = coords$y,
      niche = seurat_obj$spatial_niche
    )

    plots$niches <- ggplot(niche_data, aes(x = x, y = -y, color = niche)) +
      geom_point(size = 1.5, alpha = 0.7) +
      coord_fixed() +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      ) +
      labs(title = "Spatial Niches", color = "Niche")
  }

  # 7. Autocorrelation summary plot
  if (!is.null(spatial_results$autocorrelation)) {
    ac <- spatial_results$autocorrelation

    perm_data <- data.frame(
      value = ac$morans_i$permutation_dist
    )

    plots$autocorr <- ggplot(perm_data, aes(x = value)) +
      geom_histogram(bins = 30, fill = "gray70", color = "white") +
      geom_vline(xintercept = ac$morans_i$observed,
                 color = "red", linewidth = 1.5) +
      geom_vline(xintercept = ac$morans_i$expected,
                 color = "blue", linetype = "dashed") +
      annotate("text", x = ac$morans_i$observed, y = Inf,
               label = sprintf("Observed\nI = %.3f", ac$morans_i$observed),
               vjust = 2, hjust = -0.1, color = "red") +
      theme_minimal() +
      labs(title = "Spatial Autocorrelation (Moran's I)",
           subtitle = sprintf("p = %.4f | %s",
                              ac$morans_i$p_value,
                              ac$interpretation),
           x = "Moran's I", y = "Permutation Count")
  }

  # Combine into composite figure
  if (length(plots) >= 4) {
    combined <- (plots$pseudotime | plots$subtypes) /
      (plots$gradients | plots$transitions) +
      plot_annotation(title = "Spatial Trajectory Analysis",
                      tag_levels = "A")
  } else {
    combined <- wrap_plots(plots) +
      plot_annotation(title = "Spatial Trajectory Analysis")
  }

  plots$combined <- combined

  # Save if output directory provided
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    ggsave(file.path(output_dir, "spatial_trajectory_combined.pdf"),
           combined, width = 14, height = 12)
    ggsave(file.path(output_dir, "spatial_trajectory_combined.png"),
           combined, width = 14, height = 12, dpi = 300)

    message(sprintf("Saved spatial trajectory plots to: %s", output_dir))
  }

  return(plots)
}
