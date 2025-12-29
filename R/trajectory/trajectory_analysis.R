# ==============================================================================
# R/trajectory/multi_method_trajectory.R
# ==============================================================================
# Multi-method trajectory analysis with consensus determination
#
# Methods implemented:
#   1. Monocle3 - graph-based trajectory inference
#   2. Slingshot - principal curves through clusters
#   3. Diffusion Pseudotime - diffusion map-based ordering
# ==============================================================================

#' Run Multi-Method Trajectory Analysis
#'
#' Performs trajectory inference using multiple methods and determines consensus
#'
#' @param seurat_obj Seurat object with cell type annotations
#' @param config Configuration list
#' @param methods Character vector of methods: "monocle3", "slingshot", "diffusion"
#' @param subtype_column Column name containing cell type classifications
#' @param root_subtype Subtype to use as trajectory root
#' @param reduction Which dimensionality reduction to use (default: "umap")
#' @return List containing results from each method and consensus pseudotime
#' @export
run_multi_method_trajectory <- function(seurat_obj,
                                        config = NULL,
                                        methods = c("monocle3", "slingshot", "diffusion"),
                                        subtype_column = "module_score_subtype",
                                        root_subtype = "Basal-like",
                                        reduction = "umap") {


  message("\n========================================")
  message("Multi-Method Trajectory Analysis")
  message("========================================\n")

  # Get parameters from config
  if (!is.null(config)) {
    if (!is.null(config$trajectory$methods)) {
      methods <- config$trajectory$methods
    }
    if (!is.null(config$trajectory$monocle3$root_subtype)) {
      root_subtype <- config$trajectory$monocle3$root_subtype
    }
  }

  message(sprintf("Methods: %s", paste(methods, collapse = ", ")))
  message(sprintf("Root subtype: %s", root_subtype))
  message(sprintf("Using reduction: %s", reduction))

  # Validate inputs
  if (!subtype_column %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Subtype column '%s' not found in metadata", subtype_column))
  }

  # Check for reduction
  available_reductions <- names(seurat_obj@reductions)
  if (!tolower(reduction) %in% tolower(available_reductions)) {
    message(sprintf("  Available reductions: %s", paste(available_reductions, collapse = ", ")))
    if ("umap" %in% tolower(available_reductions)) {
      reduction <- available_reductions[tolower(available_reductions) == "umap"][1]
    } else if ("pca" %in% tolower(available_reductions)) {
      reduction <- available_reductions[tolower(available_reductions) == "pca"][1]
      message("  Falling back to PCA")
    } else {
      stop("No suitable dimensionality reduction found")
    }
  }

  # Initialize results
  results <- list()
  pseudotime_matrix <- matrix(NA, nrow = ncol(seurat_obj), ncol = length(methods))
  rownames(pseudotime_matrix) <- colnames(seurat_obj)
  colnames(pseudotime_matrix) <- methods
  methods_used <- c()

  # Get subtypes
  subtypes <- seurat_obj@meta.data[[subtype_column]]

  # ===========================================================================
  # Method 1: Monocle3
  # ===========================================================================
  if ("monocle3" %in% methods) {
    message("\n--- Method 1: Monocle3 ---")

    result <- tryCatch({
      monocle_result <- run_monocle3_trajectory(
        seurat_obj,
        config = config,
        subtype_column = subtype_column,
        root_subtype = root_subtype
      )

      results$monocle3 <- monocle_result

      # Extract pseudotime
      pt <- monocle3::pseudotime(monocle_result$cds)

      # Match to seurat cells
      common_cells <- intersect(names(pt), colnames(seurat_obj))
      pseudotime_matrix[common_cells, "monocle3"] <- pt[common_cells]

      methods_used <- c(methods_used, "monocle3")

      message(sprintf("  Pseudotime range: [%.2f, %.2f]",
                      min(pt, na.rm = TRUE), max(pt, na.rm = TRUE)))

      TRUE
    }, error = function(e) {
      warning(paste("  Monocle3 failed:", e$message))
      FALSE
    })
  }

  # ===========================================================================
  # Method 2: Slingshot
  # ===========================================================================
  if ("slingshot" %in% methods) {
    message("\n--- Method 2: Slingshot ---")

    result <- tryCatch({
      slingshot_result <- run_slingshot_trajectory(
        seurat_obj,
        config = config,
        subtype_column = subtype_column,
        root_subtype = root_subtype,
        reduction = reduction
      )

      results$slingshot <- slingshot_result

      # Extract pseudotime (use first lineage or average across lineages)
      pt <- slingshot_result$pseudotime
      pseudotime_matrix[names(pt), "slingshot"] <- pt

      methods_used <- c(methods_used, "slingshot")

      message(sprintf("  Pseudotime range: [%.2f, %.2f]",
                      min(pt, na.rm = TRUE), max(pt, na.rm = TRUE)))
      message(sprintf("  Lineages detected: %d", slingshot_result$n_lineages))

      TRUE
    }, error = function(e) {
      warning(paste("  Slingshot failed:", e$message))
      FALSE
    })
  }

  # ===========================================================================
  # Method 3: Diffusion Pseudotime
  # ===========================================================================
  if ("diffusion" %in% methods) {
    message("\n--- Method 3: Diffusion Pseudotime ---")

    result <- tryCatch({
      diffusion_result <- run_diffusion_pseudotime(
        seurat_obj,
        config = config,
        subtype_column = subtype_column,
        root_subtype = root_subtype,
        reduction = reduction
      )

      results$diffusion <- diffusion_result

      # Extract pseudotime
      pt <- diffusion_result$pseudotime
      pseudotime_matrix[names(pt), "diffusion"] <- pt

      methods_used <- c(methods_used, "diffusion")

      message(sprintf("  Pseudotime range: [%.2f, %.2f]",
                      min(pt, na.rm = TRUE), max(pt, na.rm = TRUE)))

      TRUE
    }, error = function(e) {
      warning(paste("  Diffusion pseudotime failed:", e$message))
      FALSE
    })
  }

  # ===========================================================================
  # Consensus Pseudotime
  # ===========================================================================
  message("\n--- Calculating Consensus Pseudotime ---")
  message(sprintf("Methods used: %s", paste(methods_used, collapse = ", ")))

  if (length(methods_used) == 0) {
    stop("No trajectory methods succeeded")
  }

  consensus <- calculate_consensus_pseudotime(
    pseudotime_matrix[, methods_used, drop = FALSE],
    subtypes = subtypes
  )

  # Store results
  results$consensus <- consensus
  results$methods_used <- methods_used
  results$pseudotime_matrix <- pseudotime_matrix[, methods_used, drop = FALSE]

  # Add pseudotime to seurat object metadata (return updated object)
  results$seurat_obj <- seurat_obj
  results$seurat_obj$pseudotime_consensus <- consensus$pseudotime
  results$seurat_obj$pseudotime_confidence <- consensus$confidence

  for (method in methods_used) {
    col_name <- paste0("pseudotime_", method)
    results$seurat_obj@meta.data[[col_name]] <- pseudotime_matrix[, method]
  }

  message(sprintf("\n  Consensus pseudotime range: [%.2f, %.2f]",
                  min(consensus$pseudotime, na.rm = TRUE),
                  max(consensus$pseudotime, na.rm = TRUE)))
  message(sprintf("  Mean confidence: %.3f", mean(consensus$confidence, na.rm = TRUE)))

  # Calculate method correlations
  if (length(methods_used) > 1) {
    message("\n  Method correlations (Spearman):")
    cor_matrix <- cor(pseudotime_matrix[, methods_used],
                      method = "spearman", use = "pairwise.complete.obs")
    print(round(cor_matrix, 3))
    results$method_correlations <- cor_matrix
  }

  return(results)
}


#' Run Monocle3 Trajectory Analysis
#'
#' @param seurat_obj Seurat object
#' @param config Configuration list
#' @param subtype_column Column with cell types
#' @param root_subtype Root cell type
#' @return List with cds and trajectory info
#' @keywords internal
run_monocle3_trajectory <- function(seurat_obj, config = NULL,
                                    subtype_column = "module_score_subtype",
                                    root_subtype = "Basal-like") {

  if (!requireNamespace("monocle3", quietly = TRUE)) {
    stop("monocle3 package required. Install from: https://cole-trapnell-lab.github.io/monocle3/")
  }
  if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
    stop("SeuratWrappers package required. Install with: remotes::install_github('satijalab/seurat-wrappers')")
  }

  library(monocle3)

  # Convert to cell_data_set
  message("  Converting to cell_data_set...")
  cds <- SeuratWrappers::as.cell_data_set(seurat_obj)

  # Transfer subtype info
  SummarizedExperiment::colData(cds)$subtype <- seurat_obj@meta.data[[subtype_column]]

  # Get parameters
  n_pcs <- 30
  use_partition <- TRUE
  close_loop <- FALSE

  if (!is.null(config)) {
    if (!is.null(config$trajectory$monocle3$n_pcs)) n_pcs <- config$trajectory$monocle3$n_pcs
    if (!is.null(config$trajectory$monocle3$use_partition)) use_partition <- config$trajectory$monocle3$use_partition
    if (!is.null(config$trajectory$monocle3$close_loop)) close_loop <- config$trajectory$monocle3$close_loop
  }

  # Preprocess
  message("  Preprocessing...")
  cds <- monocle3::preprocess_cds(cds, num_dim = n_pcs)

  # Reduce dimensions
  message("  Reducing dimensions...")
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")

  # Cluster
  message("  Clustering...")
  cds <- monocle3::cluster_cells(cds)

  # Learn graph
  message("  Learning trajectory graph...")
  cds <- monocle3::learn_graph(cds,
                               use_partition = use_partition,
                               close_loop = close_loop)

  # Order cells from root
  message("  Ordering cells...")
  subtypes <- SummarizedExperiment::colData(cds)$subtype
  root_cells <- colnames(cds)[subtypes == root_subtype]

  if (length(root_cells) == 0) {
    warning(sprintf("No cells found for root subtype '%s', using automatic root", root_subtype))
    cds <- monocle3::order_cells(cds)
  } else {
    message(sprintf("  Using %d cells as root (%s)", length(root_cells), root_subtype))
    cds <- monocle3::order_cells(cds, root_cells = root_cells)
  }

  # Count branch points
  graph <- monocle3::principal_graph(cds)$UMAP
  degree <- igraph::degree(graph)
  n_branch_points <- sum(degree > 2)

  return(list(
    cds = cds,
    n_branch_points = n_branch_points,
    method = "monocle3"
  ))
}


#' Run Slingshot Trajectory Analysis
#'
#' @param seurat_obj Seurat object
#' @param config Configuration list
#' @param subtype_column Column with cell types
#' @param root_subtype Root cell type
#' @param reduction Dimensionality reduction to use
#' @return List with slingshot results and pseudotime
#' @keywords internal
run_slingshot_trajectory <- function(seurat_obj, config = NULL,
                                     subtype_column = "module_score_subtype",
                                     root_subtype = "Basal-like",
                                     reduction = "umap") {

  if (!requireNamespace("slingshot", quietly = TRUE)) {
    stop("slingshot package required. Install with: BiocManager::install('slingshot')")
  }

  library(slingshot)

  # Get dimensionality reduction coordinates
  rd <- Seurat::Embeddings(seurat_obj, reduction = reduction)

  # Get cluster assignments
  clusters <- seurat_obj@meta.data[[subtype_column]]

  # Determine start cluster
  start_cluster <- root_subtype
  if (!start_cluster %in% unique(clusters)) {
    # Find closest match
    available <- unique(clusters)
    if (any(grepl("Basal", available, ignore.case = TRUE))) {
      start_cluster <- available[grepl("Basal", available, ignore.case = TRUE)][1]
    } else {
      start_cluster <- NULL
      message("  Warning: Could not determine start cluster, using automatic")
    }
  }

  # Run slingshot
  message("  Running slingshot...")

  sds <- slingshot::slingshot(
    data = rd,
    clusterLabels = clusters,
    start.clus = start_cluster,
    stretch = 2,
    extend = "n"
  )

  # Extract pseudotime
  # If multiple lineages, use average or first lineage
  pt_matrix <- slingshot::slingPseudotime(sds)

  if (is.matrix(pt_matrix) && ncol(pt_matrix) > 1) {
    message(sprintf("  Multiple lineages detected: %d", ncol(pt_matrix)))
    # Use average pseudotime across lineages (where defined)
    pseudotime <- rowMeans(pt_matrix, na.rm = TRUE)
    # For cells only in one lineage, use that value
    single_lineage <- rowSums(!is.na(pt_matrix)) == 1
    if (any(single_lineage)) {
      for (i in which(single_lineage)) {
        pseudotime[i] <- pt_matrix[i, !is.na(pt_matrix[i, ])][1]
      }
    }
  } else {
    pseudotime <- as.numeric(pt_matrix)
  }

  names(pseudotime) <- colnames(seurat_obj)

  # Get lineage info
  lineages <- slingshot::slingLineages(sds)
  n_lineages <- length(lineages)

  # Get curves
  curves <- slingshot::slingCurves(sds)

  return(list(
    sds = sds,
    pseudotime = pseudotime,
    n_lineages = n_lineages,
    lineages = lineages,
    curves = curves,
    method = "slingshot"
  ))
}


#' Run Diffusion Pseudotime Analysis
#'
#' Simple diffusion-based pseudotime using destiny package or manual calculation
#'
#' @param seurat_obj Seurat object
#' @param config Configuration list
#' @param subtype_column Column with cell types
#' @param root_subtype Root cell type
#' @param reduction Dimensionality reduction to use
#' @return List with diffusion results and pseudotime
#' @keywords internal
run_diffusion_pseudotime <- function(seurat_obj, config = NULL,
                                     subtype_column = "module_score_subtype",
                                     root_subtype = "Basal-like",
                                     reduction = "pca") {

  # Get PCA or specified reduction
  if (reduction == "pca" && "pca" %in% names(seurat_obj@reductions)) {
    data_matrix <- Seurat::Embeddings(seurat_obj, "pca")[, 1:min(30, ncol(Seurat::Embeddings(seurat_obj, "pca")))]
  } else {
    data_matrix <- Seurat::Embeddings(seurat_obj, reduction)
  }

  # Try using destiny package if available
  if (requireNamespace("destiny", quietly = TRUE)) {
    message("  Using destiny package for diffusion map...")

    dm <- destiny::DiffusionMap(data_matrix, n_pcs = NA)

    # Get diffusion components
    dc <- destiny::eigenvectors(dm)[, 1:2]

    # Find root cells
    subtypes <- seurat_obj@meta.data[[subtype_column]]
    root_cells <- which(subtypes == root_subtype)

    if (length(root_cells) == 0) {
      root_cells <- 1  # Fallback
    }

    # Calculate pseudotime as distance from root in diffusion space
    root_centroid <- colMeans(dc[root_cells, , drop = FALSE])
    pseudotime <- sqrt(rowSums((dc - matrix(root_centroid, nrow = nrow(dc), ncol = 2, byrow = TRUE))^2))

    names(pseudotime) <- colnames(seurat_obj)

    return(list(
      diffusion_map = dm,
      diffusion_components = dc,
      pseudotime = pseudotime,
      method = "diffusion_destiny"
    ))

  } else {
    message("  Using manual diffusion pseudotime calculation...")

    # Simple diffusion-based approach without destiny
    # Calculate transition matrix based on k-NN graph

    k <- 30  # Number of neighbors

    # Compute distances
    dist_matrix <- as.matrix(dist(data_matrix))

    # Build k-NN graph and transition probabilities
    n_cells <- nrow(data_matrix)
    trans_prob <- matrix(0, n_cells, n_cells)

    for (i in 1:n_cells) {
      nn_idx <- order(dist_matrix[i, ])[2:(k+1)]  # Exclude self
      sigma <- dist_matrix[i, nn_idx[k]]  # Adaptive bandwidth
      weights <- exp(-dist_matrix[i, nn_idx]^2 / (2 * sigma^2))
      trans_prob[i, nn_idx] <- weights / sum(weights)
    }

    # Symmetrize
    trans_prob <- (trans_prob + t(trans_prob)) / 2

    # Compute first few eigenvectors of transition matrix
    eigen_result <- eigen(trans_prob, symmetric = TRUE)
    dc <- eigen_result$vectors[, 2:3]  # Skip first (trivial) eigenvector

    # Find root and calculate pseudotime
    subtypes <- seurat_obj@meta.data[[subtype_column]]
    root_cells <- which(subtypes == root_subtype)

    if (length(root_cells) == 0) {
      root_cells <- 1
    }

    root_centroid <- colMeans(dc[root_cells, , drop = FALSE])
    pseudotime <- sqrt(rowSums((dc - matrix(root_centroid, nrow = nrow(dc), ncol = 2, byrow = TRUE))^2))

    names(pseudotime) <- colnames(seurat_obj)

    return(list(
      diffusion_components = dc,
      pseudotime = pseudotime,
      method = "diffusion_manual"
    ))
  }
}


#' Calculate Consensus Pseudotime
#'
#' Combines pseudotime from multiple methods into consensus
#'
#' @param pseudotime_matrix Matrix with cells as rows, methods as columns
#' @param subtypes Cell type annotations (for validation)
#' @return List with consensus pseudotime and confidence
#' @keywords internal
calculate_consensus_pseudotime <- function(pseudotime_matrix, subtypes = NULL) {

  n_methods <- ncol(pseudotime_matrix)
  n_cells <- nrow(pseudotime_matrix)

  if (n_methods == 1) {
    # Single method - just normalize
    pt <- pseudotime_matrix[, 1]
    pt_norm <- (pt - min(pt, na.rm = TRUE)) / (max(pt, na.rm = TRUE) - min(pt, na.rm = TRUE))

    return(list(
      pseudotime = pt_norm,
      confidence = rep(1.0, n_cells),
      n_methods = 1
    ))
  }

  # Normalize each method to [0, 1] for comparison
  pt_normalized <- apply(pseudotime_matrix, 2, function(x) {
    if (all(is.na(x))) return(x)
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  })

  # Check if methods are anti-correlated (reversed direction)
  # If so, flip the reversed ones
  if (n_methods > 1) {
    reference <- pt_normalized[, 1]
    for (i in 2:n_methods) {
      cor_val <- cor(reference, pt_normalized[, i], use = "pairwise.complete.obs")
      if (!is.na(cor_val) && cor_val < -0.5) {
        message(sprintf("  Flipping %s (negative correlation: %.2f)",
                        colnames(pt_normalized)[i], cor_val))
        pt_normalized[, i] <- 1 - pt_normalized[, i]
      }
    }
  }

  # Calculate consensus as mean (robust to missing values)
  consensus_pt <- rowMeans(pt_normalized, na.rm = TRUE)

  # Calculate confidence based on agreement
  # Confidence = 1 - (SD across methods) for each cell
  pt_sd <- apply(pt_normalized, 1, sd, na.rm = TRUE)
  pt_sd[is.na(pt_sd)] <- 0  # Single method = 0 SD

  # Normalize SD to [0, 1] and invert for confidence
  max_possible_sd <- 0.5  # Max SD for uniform [0,1] values
  confidence <- 1 - pmin(pt_sd / max_possible_sd, 1)

  # Count methods with valid values per cell
  n_valid <- rowSums(!is.na(pseudotime_matrix))

  return(list(
    pseudotime = consensus_pt,
    confidence = confidence,
    n_methods = n_methods,
    n_valid_per_cell = n_valid,
    normalized_matrix = pt_normalized
  ))
}


#' Run Per-Sample Multi-Method Trajectory
#'
#' Runs trajectory analysis separately for each sample using multiple methods
#'
#' @param seurat_obj Seurat object
#' @param config Configuration list
#' @param methods Trajectory methods to use
#' @param sample_column Column containing sample IDs
#' @param subtype_column Column containing subtypes
#' @return List with per-sample results and cross-sample concordance
#' @export
run_per_sample_multi_trajectory <- function(seurat_obj, config = NULL,
                                            methods = c("monocle3", "slingshot"),
                                            sample_column = "sample",
                                            subtype_column = "module_score_subtype") {

  message("\n========================================")
  message("Per-Sample Multi-Method Trajectory Analysis")
  message("========================================\n")

  # Get parameters from config
  if (!is.null(config)) {
    if (!is.null(config$trajectory$methods)) {
      methods <- config$trajectory$methods
    }
  }

  # Get unique samples
  samples <- unique(seurat_obj@meta.data[[sample_column]])
  samples <- samples[!is.na(samples)]

  message(sprintf("Samples: %s", paste(samples, collapse = ", ")))
  message(sprintf("Methods: %s", paste(methods, collapse = ", ")))

  # Minimum cells
  min_cells <- 100
  if (!is.null(config$trajectory$per_sample$min_cells_per_sample)) {
    min_cells <- config$trajectory$per_sample$min_cells_per_sample
  }

  # Store results
  all_results <- list()

  for (sample_id in samples) {
    message(sprintf("\n=== Sample: %s ===", sample_id))

    # Subset
    cells_in_sample <- seurat_obj@meta.data[[sample_column]] == sample_id
    cells_in_sample[is.na(cells_in_sample)] <- FALSE
    n_cells <- sum(cells_in_sample)

    message(sprintf("  Cells: %d", n_cells))

    if (n_cells < min_cells) {
      message(sprintf("  Skipping: fewer than %d cells", min_cells))
      all_results[[sample_id]] <- list(
        sample = sample_id,
        n_cells = n_cells,
        skipped = TRUE,
        reason = "insufficient_cells"
      )
      next
    }

    # Subset Seurat
    seurat_sample <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_in_sample])

    # Run multi-method trajectory
    tryCatch({
      sample_result <- run_multi_method_trajectory(
        seurat_sample,
        config = config,
        methods = methods,
        subtype_column = subtype_column
      )

      # Calculate per-subtype mean pseudotime
      subtypes <- seurat_sample@meta.data[[subtype_column]]
      pt_by_subtype <- tapply(sample_result$consensus$pseudotime, subtypes, mean, na.rm = TRUE)
      subtype_rank <- rank(pt_by_subtype)

      all_results[[sample_id]] <- list(
        sample = sample_id,
        n_cells = n_cells,
        skipped = FALSE,
        results = sample_result,
        pseudotime_by_subtype = pt_by_subtype,
        subtype_rank = subtype_rank,
        subtype_order = names(sort(pt_by_subtype)),
        methods_used = sample_result$methods_used
      )

      message(sprintf("  Subtype order: %s",
                      paste(names(sort(pt_by_subtype)), collapse = " → ")))

    }, error = function(e) {
      message(sprintf("  ERROR: %s", e$message))
      all_results[[sample_id]] <- list(
        sample = sample_id,
        n_cells = n_cells,
        skipped = TRUE,
        reason = "trajectory_failed",
        error = e$message
      )
    })
  }

  # Calculate cross-sample concordance
  message("\n=== Cross-Sample Concordance ===")
  concordance <- calculate_cross_sample_concordance(all_results)
  attr(all_results, "concordance") <- concordance

  return(all_results)
}


#' Calculate Cross-Sample Concordance
#'
#' @param results List of per-sample results
#' @return List with concordance metrics
#' @keywords internal
calculate_cross_sample_concordance <- function(results) {

  # Get valid results
  valid_results <- results[sapply(results, function(x) !x$skipped)]

  if (length(valid_results) < 2) {
    message("  Insufficient samples for concordance analysis")
    return(list(mean_correlation = NA, n_samples = length(valid_results)))
  }

  # Build rank matrix
  all_subtypes <- unique(unlist(lapply(valid_results, function(x) names(x$subtype_rank))))
  rank_matrix <- matrix(NA, nrow = length(all_subtypes), ncol = length(valid_results))
  rownames(rank_matrix) <- all_subtypes
  colnames(rank_matrix) <- names(valid_results)

  for (sample_id in names(valid_results)) {
    ranks <- valid_results[[sample_id]]$subtype_rank
    for (subtype in names(ranks)) {
      rank_matrix[subtype, sample_id] <- ranks[subtype]
    }
  }

  # Pairwise Spearman correlations
  n_samples <- ncol(rank_matrix)
  correlations <- c()

  for (i in 1:(n_samples - 1)) {
    for (j in (i + 1):n_samples) {
      shared <- !is.na(rank_matrix[, i]) & !is.na(rank_matrix[, j])
      if (sum(shared) >= 3) {
        rho <- cor(rank_matrix[shared, i], rank_matrix[shared, j],
                   method = "spearman", use = "complete.obs")
        correlations <- c(correlations, rho)
      }
    }
  }

  concordance <- list(
    rank_matrix = rank_matrix,
    mean_correlation = mean(correlations, na.rm = TRUE),
    sd_correlation = sd(correlations, na.rm = TRUE),
    min_correlation = min(correlations, na.rm = TRUE),
    max_correlation = max(correlations, na.rm = TRUE),
    n_comparisons = length(correlations),
    n_samples = n_samples
  )

  message(sprintf("  Samples analyzed: %d", n_samples))
  message(sprintf("  Mean Spearman rho: %.3f ± %.3f",
                  concordance$mean_correlation,
                  ifelse(is.na(concordance$sd_correlation), 0, concordance$sd_correlation)))

  return(concordance)
}


#' Validate Trajectory Against Expected Order
#'
#' @param results Trajectory results (single or per-sample)
#' @param expected_order Expected subtype order (early to late)
#' @return Data frame with validation metrics
#' @export
validate_trajectory_order <- function(results,
                                      expected_order = c("Basal-like", "Transit-Amplifying",
                                                         "Intermediate", "Specialized")) {

  message("\n--- Validating Trajectory Order ---")
  message(sprintf("Expected: %s", paste(expected_order, collapse = " → ")))

  # Handle different result types
  if ("consensus" %in% names(results)) {
    # Single multi-method result
    pt <- results$consensus$pseudotime
    subtypes <- results$seurat_obj@meta.data$module_score_subtype

    mean_pt <- tapply(pt, subtypes, mean, na.rm = TRUE)
    observed_ranks <- rank(mean_pt)
    observed_order <- names(sort(mean_pt))

    # Correlate with expected
    expected_ranks <- seq_along(expected_order)
    names(expected_ranks) <- expected_order

    shared <- intersect(names(observed_ranks), expected_order)

    if (length(shared) >= 3) {
      cor_test <- cor.test(observed_ranks[shared], expected_ranks[shared], method = "spearman")

      validation <- data.frame(
        sample = "pooled",
        spearman_rho = cor_test$estimate,
        p_value = cor_test$p.value,
        n_subtypes = length(shared),
        observed_order = paste(observed_order, collapse = " → "),
        concordant = cor_test$estimate > 0.5 && cor_test$p.value < 0.05
      )
    } else {
      validation <- data.frame(
        sample = "pooled",
        spearman_rho = NA,
        p_value = NA,
        n_subtypes = length(shared),
        observed_order = paste(observed_order, collapse = " → "),
        concordant = NA
      )
    }

  } else {
    # Per-sample results
    validation <- do.call(rbind, lapply(names(results), function(sample_id) {
      res <- results[[sample_id]]

      if (res$skipped) {
        return(data.frame(
          sample = sample_id,
          spearman_rho = NA,
          p_value = NA,
          n_subtypes = NA,
          observed_order = NA,
          concordant = NA
        ))
      }

      observed_ranks <- res$subtype_rank
      observed_order <- res$subtype_order

      expected_ranks <- seq_along(expected_order)
      names(expected_ranks) <- expected_order

      shared <- intersect(names(observed_ranks), expected_order)

      if (length(shared) >= 3) {
        cor_test <- cor.test(observed_ranks[shared], expected_ranks[shared], method = "spearman")

        data.frame(
          sample = sample_id,
          spearman_rho = cor_test$estimate,
          p_value = cor_test$p.value,
          n_subtypes = length(shared),
          observed_order = paste(observed_order, collapse = " → "),
          concordant = cor_test$estimate > 0.5 && cor_test$p.value < 0.05
        )
      } else {
        data.frame(
          sample = sample_id,
          spearman_rho = NA,
          p_value = NA,
          n_subtypes = length(shared),
          observed_order = paste(observed_order, collapse = " → "),
          concordant = NA
        )
      }
    }))
  }

  message(sprintf("\nValidation results:"))
  message(sprintf("  Concordant samples: %d/%d",
                  sum(validation$concordant, na.rm = TRUE),
                  sum(!is.na(validation$concordant))))
  message(sprintf("  Mean Spearman rho: %.3f", mean(validation$spearman_rho, na.rm = TRUE)))

  return(validation)
}


#' Plot Multi-Method Trajectory Comparison
#'
#' Creates diagnostic plots comparing trajectory methods
#'
#' @param results Results from run_multi_method_trajectory
#' @param config Configuration list
#' @return ggplot object
#' @export
plot_trajectory_comparison <- function(results, config = NULL) {

  require(ggplot2)
  require(patchwork)

  plots <- list()

  # Get colors
  if (!is.null(config)) {
    colors <- config$visualization$colors$epithelial_subtypes
    if (!is.null(colors)) colors <- unlist(colors)
  } else {
    colors <- NULL
  }

  seurat_obj <- results$seurat_obj
  methods_used <- results$methods_used

  # Plot 1: Consensus pseudotime by subtype
  pt_data <- data.frame(
    pseudotime = results$consensus$pseudotime,
    confidence = results$consensus$confidence,
    subtype = seurat_obj$module_score_subtype
  )

  plots$consensus_violin <- ggplot(pt_data, aes(x = reorder(subtype, pseudotime),
                                                y = pseudotime, fill = subtype)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Consensus Pseudotime by Subtype",
         x = "", y = "Pseudotime")

  # Plot 2: Method comparison scatter
  if (length(methods_used) >= 2) {
    pt_matrix <- results$pseudotime_matrix

    method1 <- methods_used[1]
    method2 <- methods_used[2]

    compare_data <- data.frame(
      method1 = pt_matrix[, method1],
      method2 = pt_matrix[, method2],
      subtype = seurat_obj$module_score_subtype
    )

    cor_val <- cor(compare_data$method1, compare_data$method2,
                   use = "pairwise.complete.obs", method = "spearman")

    plots$method_scatter <- ggplot(compare_data, aes(x = method1, y = method2, color = subtype)) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = colors) +
      theme_minimal() +
      labs(title = sprintf("%s vs %s (ρ = %.3f)", method1, method2, cor_val),
           x = method1, y = method2, color = "Subtype")
  }

  # Plot 3: Confidence distribution
  plots$confidence <- ggplot(pt_data, aes(x = confidence)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    geom_vline(xintercept = mean(pt_data$confidence), linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = "Pseudotime Confidence Distribution",
         subtitle = sprintf("Mean = %.3f", mean(pt_data$confidence)),
         x = "Confidence", y = "Cells")

  # Plot 4: Confidence by subtype
  plots$confidence_subtype <- ggplot(pt_data, aes(x = subtype, y = confidence, fill = subtype)) +
    geom_boxplot() +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Confidence by Subtype", x = "", y = "Confidence")

  # Combine
  combined <- (plots$consensus_violin | plots$confidence) /
    (plots$method_scatter | plots$confidence_subtype) +
    plot_annotation(title = "Multi-Method Trajectory Comparison",
                    tag_levels = "A")

  return(combined)
}


# ==============================================================================
# TRAJECTORY STATISTICAL EVALUATION
# ==============================================================================

#' Evaluate Trajectory Analysis Statistics
#'
#' Computes comprehensive statistics to evaluate trajectory inference quality,
#' including method agreement, biological validation, and cell cycle assessment
#'
#' @param seurat_obj Seurat object with pseudotime and cell cycle scores
#' @param results Results list from run_multi_method_trajectory
#' @param config Configuration list (optional)
#' @param differentiation_markers Named list of markers expected to change along trajectory
#' @param expected_order Expected order of subtypes along trajectory
#' @return List with statistical evaluation results
#' @export
evaluate_trajectory_statistics <- function(seurat_obj,
                                           results,
                                           config = NULL,
                                           differentiation_markers = NULL,
                                           expected_order = NULL) {

  message("\n========================================")
  message("Trajectory Analysis Statistics")
  message("========================================\n")

  stats <- list()

  # Get expected order from config if not provided
  if (is.null(expected_order) && !is.null(config$trajectory$expected_order)) {
    expected_order <- config$trajectory$expected_order
  }

  # Default differentiation markers if not provided
  if (is.null(differentiation_markers)) {
    differentiation_markers <- list(
      # Markers expected to DECREASE along differentiation (early/progenitor)
      early = c("SOX9", "KRT5", "KRT14", "TP63", "ITGA6", "ITGB1"),
      # Markers expected to INCREASE along differentiation (late/mature)
      late = c("KRT20", "MUC2", "FABP1", "SI", "ALPI", "VIL1")
    )
  }

  # ===========================================================================
  # 1. Method Agreement Statistics
  # ===========================================================================
  message("--- Method Agreement ---")

  pt_matrix <- results$pseudotime_matrix
  methods_used <- colnames(pt_matrix)[colSums(!is.na(pt_matrix)) > 0]
  n_methods <- length(methods_used)

  if (n_methods >= 2) {
    # Pairwise Spearman correlations
    cor_results <- list()
    cor_matrix <- matrix(NA, nrow = n_methods, ncol = n_methods)
    rownames(cor_matrix) <- colnames(cor_matrix) <- methods_used

    for (i in 1:n_methods) {
      for (j in 1:n_methods) {
        if (i != j) {
          rho <- cor(pt_matrix[, methods_used[i]], pt_matrix[, methods_used[j]],
                     use = "pairwise.complete.obs", method = "spearman")
          cor_matrix[i, j] <- rho

          if (i < j) {
            pair_name <- paste(methods_used[i], "vs", methods_used[j])

            # Statistical test
            test <- cor.test(pt_matrix[, methods_used[i]], pt_matrix[, methods_used[j]],
                             method = "spearman")

            cor_results[[pair_name]] <- list(
              rho = rho,
              p_value = test$p.value,
              n = sum(complete.cases(pt_matrix[, c(methods_used[i], methods_used[j])]))
            )

            message(sprintf("  %s: ρ = %.3f (p = %.2e)", pair_name, rho, test$p.value))
          }
        } else {
          cor_matrix[i, j] <- 1
        }
      }
    }

    stats$method_correlations <- cor_results
    stats$correlation_matrix <- cor_matrix

    # Mean pairwise correlation
    upper_tri <- cor_matrix[upper.tri(cor_matrix)]
    stats$mean_correlation <- mean(upper_tri, na.rm = TRUE)
    message(sprintf("  Mean pairwise correlation: ρ = %.3f", stats$mean_correlation))

  } else {
    message("  Only one method - skipping agreement statistics")
    stats$mean_correlation <- NA
  }

  # ===========================================================================
  # 2. Cell Cycle Confounding Assessment
  # ===========================================================================
  message("\n--- Cell Cycle Confounding ---")

  has_cell_cycle <- "cc_consensus" %in% colnames(seurat_obj@meta.data)

  if (has_cell_cycle) {
    pseudotime <- seurat_obj$pseudotime_consensus

    # Convert phase to numeric for correlation
    phase_numeric <- as.numeric(factor(seurat_obj$cc_consensus,
                                       levels = c("G1", "S", "G2M")))

    # Correlation between pseudotime and cell cycle
    cc_cor <- cor(pseudotime, phase_numeric,
                  use = "pairwise.complete.obs", method = "spearman")
    cc_test <- cor.test(pseudotime, phase_numeric, method = "spearman")

    stats$cell_cycle_correlation <- list(
      rho = cc_cor,
      p_value = cc_test$p.value,
      interpretation = ifelse(abs(cc_cor) > 0.3, "POTENTIAL CONFOUNDING", "LOW CONFOUNDING")
    )

    message(sprintf("  Pseudotime ~ Cell Cycle Phase: ρ = %.3f (p = %.2e)",
                    cc_cor, cc_test$p.value))
    message(sprintf("  Interpretation: %s", stats$cell_cycle_correlation$interpretation))

    # Cycling percentage along trajectory
    # Bin pseudotime into quintiles
    pt_bins <- cut(pseudotime, breaks = 5, labels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%"))

    cycling_by_bin <- tapply(seurat_obj$cc_consensus %in% c("S", "G2M"), pt_bins, mean, na.rm = TRUE) * 100

    stats$cycling_along_trajectory <- data.frame(
      pseudotime_bin = names(cycling_by_bin),
      cycling_percent = as.numeric(cycling_by_bin)
    )

    message("  Cycling cells (%) along trajectory:")
    for (i in 1:length(cycling_by_bin)) {
      message(sprintf("    %s: %.1f%%", names(cycling_by_bin)[i], cycling_by_bin[i]))
    }

    # Expected: cycling should decrease along differentiation
    cycling_trend <- cor(1:5, cycling_by_bin, method = "spearman")
    stats$cycling_trend <- cycling_trend

    if (cycling_trend < -0.5) {
      message("  ✓ Cycling decreases along trajectory (expected for differentiation)")
    } else if (cycling_trend > 0.5) {
      message("  ⚠ Cycling increases along trajectory (unexpected)")
    } else {
      message("  ~ No clear cycling trend along trajectory")
    }

    # Use cell cycle confidence
    if ("cc_confidence" %in% colnames(seurat_obj@meta.data)) {
      mean_cc_conf <- mean(seurat_obj$cc_confidence, na.rm = TRUE)
      stats$cell_cycle_confidence <- mean_cc_conf
      message(sprintf("  Mean cell cycle confidence: %.3f", mean_cc_conf))
    }

  } else {
    message("  Cell cycle data not available - skipping")
    stats$cell_cycle_correlation <- NULL
  }

  # ===========================================================================
  # 3. Biological Validation - Marker Correlations
  # ===========================================================================
  message("\n--- Biological Validation (Marker Trends) ---")

  pseudotime <- seurat_obj$pseudotime_consensus

  # Check available markers
  available_early <- intersect(differentiation_markers$early, rownames(seurat_obj))
  available_late <- intersect(differentiation_markers$late, rownames(seurat_obj))

  message(sprintf("  Early markers available: %d/%d",
                  length(available_early), length(differentiation_markers$early)))
  message(sprintf("  Late markers available: %d/%d",
                  length(available_late), length(differentiation_markers$late)))

  marker_results <- list()

  if (length(available_early) > 0 || length(available_late) > 0) {
    # Get expression data
    all_markers <- c(available_early, available_late)
    expr_data <- as.matrix(Seurat::GetAssayData(seurat_obj, layer = "data")[all_markers, , drop = FALSE])

    # Test each marker
    for (marker in all_markers) {
      expected_direction <- ifelse(marker %in% available_early, "negative", "positive")

      rho <- cor(pseudotime, expr_data[marker, ],
                 use = "pairwise.complete.obs", method = "spearman")
      test <- cor.test(pseudotime, expr_data[marker, ], method = "spearman")

      # Check if direction matches expectation
      if (expected_direction == "negative") {
        correct_direction <- rho < 0
      } else {
        correct_direction <- rho > 0
      }

      marker_results[[marker]] <- list(
        rho = rho,
        p_value = test$p.value,
        expected = expected_direction,
        correct_direction = correct_direction,
        significant = test$p.value < 0.05
      )
    }

    # Summary
    marker_df <- do.call(rbind, lapply(names(marker_results), function(m) {
      data.frame(
        marker = m,
        rho = marker_results[[m]]$rho,
        p_value = marker_results[[m]]$p_value,
        expected = marker_results[[m]]$expected,
        correct_direction = marker_results[[m]]$correct_direction,
        significant = marker_results[[m]]$significant
      )
    }))

    stats$marker_validation <- marker_df

    # Print results
    message("\n  Early markers (expected: negative correlation):")
    for (marker in available_early) {
      r <- marker_results[[marker]]
      status <- ifelse(r$correct_direction & r$significant, "✓",
                       ifelse(r$correct_direction, "~", "✗"))
      message(sprintf("    %s %s: ρ = %.3f (p = %.2e)", status, marker, r$rho, r$p_value))
    }

    message("\n  Late markers (expected: positive correlation):")
    for (marker in available_late) {
      r <- marker_results[[marker]]
      status <- ifelse(r$correct_direction & r$significant, "✓",
                       ifelse(r$correct_direction, "~", "✗"))
      message(sprintf("    %s %s: ρ = %.3f (p = %.2e)", status, marker, r$rho, r$p_value))
    }

    # Validation summary
    n_correct <- sum(marker_df$correct_direction)
    n_sig_correct <- sum(marker_df$correct_direction & marker_df$significant)
    n_total <- nrow(marker_df)

    stats$marker_validation_summary <- list(
      n_correct_direction = n_correct,
      n_significant_correct = n_sig_correct,
      n_total = n_total,
      percent_correct = n_correct / n_total * 100,
      validated = n_sig_correct >= n_total / 2
    )

    message(sprintf("\n  Validation: %d/%d correct direction, %d significant",
                    n_correct, n_total, n_sig_correct))
    message(sprintf("  Biological validation: %s",
                    ifelse(stats$marker_validation_summary$validated, "PASSED", "REVIEW NEEDED")))

  } else {
    message("  No differentiation markers found in dataset")
    stats$marker_validation <- NULL
  }

  # ===========================================================================
  # 4. Trajectory Order Validation
  # ===========================================================================
  message("\n--- Trajectory Order Validation ---")

  if (!is.null(expected_order)) {
    subtype_col <- ifelse("module_score_subtype" %in% colnames(seurat_obj@meta.data),
                          "module_score_subtype", "subtype")

    subtypes <- seurat_obj@meta.data[[subtype_col]]

    # Calculate mean pseudotime per subtype
    mean_pt_by_subtype <- tapply(pseudotime, subtypes, mean, na.rm = TRUE)

    # Get observed order
    observed_order <- names(sort(mean_pt_by_subtype))

    # Filter to only subtypes in expected order
    observed_in_expected <- observed_order[observed_order %in% expected_order]
    expected_filtered <- expected_order[expected_order %in% observed_order]

    message(sprintf("  Expected order: %s", paste(expected_filtered, collapse = " → ")))
    message(sprintf("  Observed order: %s", paste(observed_in_expected, collapse = " → ")))

    # Spearman correlation of ranks
    if (length(observed_in_expected) >= 3) {
      expected_ranks <- match(observed_in_expected, expected_filtered)
      observed_ranks <- 1:length(observed_in_expected)

      rank_cor <- cor(expected_ranks, observed_ranks, method = "spearman")

      stats$order_validation <- list(
        expected = expected_filtered,
        observed = observed_in_expected,
        rank_correlation = rank_cor,
        order_preserved = rank_cor > 0.7
      )

      message(sprintf("  Rank correlation: ρ = %.3f", rank_cor))
      message(sprintf("  Order preserved: %s", ifelse(rank_cor > 0.7, "YES", "NO")))
    } else {
      message("  Too few shared subtypes for order validation")
      stats$order_validation <- NULL
    }

    # Mean pseudotime table
    stats$mean_pseudotime_by_subtype <- data.frame(
      subtype = names(mean_pt_by_subtype),
      mean_pseudotime = as.numeric(mean_pt_by_subtype),
      sd_pseudotime = as.numeric(tapply(pseudotime, subtypes, sd, na.rm = TRUE)),
      n_cells = as.numeric(table(subtypes)[names(mean_pt_by_subtype)])
    )
    stats$mean_pseudotime_by_subtype <- stats$mean_pseudotime_by_subtype[
      order(stats$mean_pseudotime_by_subtype$mean_pseudotime), ]

  } else {
    message("  No expected order provided - skipping")
    stats$order_validation <- NULL
  }

  # ===========================================================================
  # 5. Confidence and Uncertainty Statistics
  # ===========================================================================
  message("\n--- Confidence Statistics ---")

  if ("pseudotime_confidence" %in% colnames(seurat_obj@meta.data)) {
    conf <- seurat_obj$pseudotime_confidence

    stats$confidence_stats <- list(
      mean = mean(conf, na.rm = TRUE),
      median = median(conf, na.rm = TRUE),
      sd = sd(conf, na.rm = TRUE),
      high_confidence_pct = mean(conf >= 0.7, na.rm = TRUE) * 100,
      low_confidence_pct = mean(conf < 0.4, na.rm = TRUE) * 100
    )

    message(sprintf("  Mean confidence: %.3f ± %.3f",
                    stats$confidence_stats$mean, stats$confidence_stats$sd))
    message(sprintf("  High confidence (≥0.7): %.1f%%", stats$confidence_stats$high_confidence_pct))
    message(sprintf("  Low confidence (<0.4): %.1f%%", stats$confidence_stats$low_confidence_pct))
  }

  # ===========================================================================
  # 6. Per-Sample Concordance (if available)
  # ===========================================================================
  if (!is.null(results$per_sample) && !is.null(results$per_sample$concordance)) {
    message("\n--- Per-Sample Concordance ---")

    conc <- results$per_sample$concordance

    # Mean off-diagonal correlation
    if (is.matrix(conc) && nrow(conc) > 1) {
      off_diag <- conc[upper.tri(conc)]
      stats$sample_concordance <- list(
        mean_correlation = mean(off_diag, na.rm = TRUE),
        min_correlation = min(off_diag, na.rm = TRUE),
        max_correlation = max(off_diag, na.rm = TRUE)
      )

      message(sprintf("  Mean cross-sample correlation: %.3f (range: %.3f - %.3f)",
                      stats$sample_concordance$mean_correlation,
                      stats$sample_concordance$min_correlation,
                      stats$sample_concordance$max_correlation))
    }
  }

  # ===========================================================================
  # 7. Overall Quality Score
  # ===========================================================================
  message("\n--- Overall Quality Score ---")

  quality_components <- c()

  # Component 1: Method agreement (weight: 25%)
  if (!is.na(stats$mean_correlation)) {
    # Correlation of 0.7+ = max score
    agreement_score <- min(stats$mean_correlation / 0.7, 1)
    quality_components["method_agreement"] <- agreement_score * 0.25
  }

  # Component 2: Biological validation (weight: 30%)
  if (!is.null(stats$marker_validation_summary)) {
    marker_score <- stats$marker_validation_summary$percent_correct / 100
    quality_components["marker_validation"] <- marker_score * 0.30
  }

  # Component 3: Order preservation (weight: 25%)
  if (!is.null(stats$order_validation)) {
    # Correlation > 0.7 = full score
    order_score <- max(0, (stats$order_validation$rank_correlation + 1) / 2)
    quality_components["order_preservation"] <- order_score * 0.25
  }

  # Component 4: Low cell cycle confounding (weight: 10%)
  if (!is.null(stats$cell_cycle_correlation)) {
    # Lower correlation = better (less confounding)
    cc_score <- 1 - min(abs(stats$cell_cycle_correlation$rho), 1)
    quality_components["cell_cycle"] <- cc_score * 0.10
  }

  # Component 5: Confidence (weight: 10%)
  if (!is.null(stats$confidence_stats)) {
    conf_score <- stats$confidence_stats$mean
    quality_components["confidence"] <- conf_score * 0.10
  }

  overall_score <- sum(quality_components, na.rm = TRUE)

  # Normalize to available components
  max_possible <- sum(c(0.25, 0.30, 0.25, 0.10, 0.10)[!is.na(match(
    c("method_agreement", "marker_validation", "order_preservation", "cell_cycle", "confidence"),
    names(quality_components)
  ))])

  if (max_possible > 0) {
    overall_score <- overall_score / max_possible
  }

  stats$quality_score <- list(
    components = quality_components,
    overall = overall_score,
    grade = ifelse(overall_score >= 0.8, "Excellent",
                   ifelse(overall_score >= 0.6, "Good",
                          ifelse(overall_score >= 0.4, "Acceptable", "Review Needed")))
  )

  message(sprintf("\n  Quality Score: %.2f / 1.00 (%s)", overall_score, stats$quality_score$grade))
  message("  Components:")
  for (comp in names(quality_components)) {
    message(sprintf("    %s: %.3f", comp, quality_components[comp]))
  }

  message("\n========================================\n")

  return(stats)
}


#' Generate Trajectory Statistics Report
#'
#' Creates a formatted markdown report of trajectory analysis statistics
#'
#' @param stats_results Results from evaluate_trajectory_statistics
#' @param output_path Path to save the report (optional)
#' @return Character string with formatted report
#' @export
generate_trajectory_report <- function(stats_results, output_path = NULL) {

  report <- c(
    "# Trajectory Analysis Evaluation Report",
    paste0("Generated: ", Sys.time()),
    "",
    "## Quality Score",
    sprintf("**Overall: %.2f / 1.00 (%s)**",
            stats_results$quality_score$overall,
            stats_results$quality_score$grade),
    ""
  )

  # Method agreement section
  if (!is.null(stats_results$method_correlations)) {
    report <- c(report, "## Method Agreement", "")
    report <- c(report, "| Comparison | Spearman ρ | p-value |")
    report <- c(report, "|------------|------------|---------|")

    for (name in names(stats_results$method_correlations)) {
      r <- stats_results$method_correlations[[name]]
      report <- c(report, sprintf("| %s | %.3f | %.2e |", name, r$rho, r$p_value))
    }
    report <- c(report, "", sprintf("**Mean correlation: %.3f**", stats_results$mean_correlation), "")
  }

  # Cell cycle section
  if (!is.null(stats_results$cell_cycle_correlation)) {
    cc <- stats_results$cell_cycle_correlation
    report <- c(report, "## Cell Cycle Confounding", "",
                sprintf("- Pseudotime ~ Phase correlation: ρ = %.3f (p = %.2e)", cc$rho, cc$p_value),
                sprintf("- Interpretation: **%s**", cc$interpretation), "")

    if (!is.null(stats_results$cycling_along_trajectory)) {
      report <- c(report, "### Cycling Cells Along Trajectory", "")
      report <- c(report, "| Pseudotime Bin | % Cycling |")
      report <- c(report, "|----------------|-----------|")
      ct <- stats_results$cycling_along_trajectory
      for (i in 1:nrow(ct)) {
        report <- c(report, sprintf("| %s | %.1f%% |", ct$pseudotime_bin[i], ct$cycling_percent[i]))
      }
      report <- c(report, "")
    }
  }

  # Marker validation
  if (!is.null(stats_results$marker_validation)) {
    report <- c(report, "## Biological Validation (Marker Trends)", "")
    report <- c(report, "| Marker | Spearman ρ | p-value | Expected | Correct |")
    report <- c(report, "|--------|------------|---------|----------|---------|")

    mv <- stats_results$marker_validation
    for (i in 1:nrow(mv)) {
      report <- c(report, sprintf("| %s | %.3f | %.2e | %s | %s |",
                                  mv$marker[i], mv$rho[i], mv$p_value[i],
                                  mv$expected[i],
                                  ifelse(mv$correct_direction[i], "✓", "✗")))
    }

    if (!is.null(stats_results$marker_validation_summary)) {
      mvs <- stats_results$marker_validation_summary
      report <- c(report, "",
                  sprintf("**Summary:** %d/%d correct direction (%.1f%%), Validation: %s",
                          mvs$n_correct_direction, mvs$n_total, mvs$percent_correct,
                          ifelse(mvs$validated, "PASSED", "REVIEW NEEDED")), "")
    }
  }

  # Order validation
  if (!is.null(stats_results$order_validation)) {
    ov <- stats_results$order_validation
    report <- c(report, "## Trajectory Order Validation", "",
                sprintf("- Expected: %s", paste(ov$expected, collapse = " → ")),
                sprintf("- Observed: %s", paste(ov$observed, collapse = " → ")),
                sprintf("- Rank correlation: ρ = %.3f", ov$rank_correlation),
                sprintf("- Order preserved: **%s**", ifelse(ov$order_preserved, "YES", "NO")), "")
  }

  # Pseudotime by subtype
  if (!is.null(stats_results$mean_pseudotime_by_subtype)) {
    report <- c(report, "## Mean Pseudotime by Subtype", "")
    report <- c(report, "| Subtype | Mean PT | SD | N Cells |")
    report <- c(report, "|---------|---------|-------|---------|")

    pt <- stats_results$mean_pseudotime_by_subtype
    for (i in 1:nrow(pt)) {
      report <- c(report, sprintf("| %s | %.3f | %.3f | %d |",
                                  pt$subtype[i], pt$mean_pseudotime[i],
                                  pt$sd_pseudotime[i], pt$n_cells[i]))
    }
    report <- c(report, "")
  }

  report_text <- paste(report, collapse = "\n")

  if (!is.null(output_path)) {
    writeLines(report_text, output_path)
    message(sprintf("Report saved to: %s", output_path))
  }

  return(report_text)
}


#' Plot Cell Cycle Along Trajectory
#'
#' Visualizes the relationship between pseudotime and cell cycle
#'
#' @param seurat_obj Seurat object with pseudotime and cell cycle
#' @param config Configuration list (optional)
#' @return ggplot object
#' @export
plot_cell_cycle_trajectory <- function(seurat_obj, config = NULL) {

  require(ggplot2)
  require(patchwork)

  if (!"pseudotime_consensus" %in% colnames(seurat_obj@meta.data) ||
      !"cc_consensus" %in% colnames(seurat_obj@meta.data)) {
    stop("Requires both pseudotime_consensus and cc_consensus in metadata")
  }

  # Get colors
  if (!is.null(config) && !is.null(config$visualization$colors$cell_cycle)) {
    cc_colors <- unlist(config$visualization$colors$cell_cycle)
  } else {
    cc_colors <- c(G1 = "#F8766D", S = "#00BA38", G2M = "#619CFF")
  }

  plot_data <- data.frame(
    pseudotime = seurat_obj$pseudotime_consensus,
    phase = seurat_obj$cc_consensus,
    confidence = seurat_obj$cc_confidence
  )

  # Plot 1: Pseudotime distribution by phase
  p1 <- ggplot(plot_data, aes(x = phase, y = pseudotime, fill = phase)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = cc_colors) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Pseudotime by Cell Cycle Phase", x = "", y = "Pseudotime")

  # Plot 2: Density along pseudotime by phase
  p2 <- ggplot(plot_data, aes(x = pseudotime, fill = phase)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = cc_colors) +
    theme_minimal() +
    labs(title = "Cell Cycle Phase Density Along Trajectory",
         x = "Pseudotime", y = "Density", fill = "Phase")

  # Plot 3: Proportion of cycling cells along trajectory
  plot_data$pt_bin <- cut(plot_data$pseudotime, breaks = 10, labels = FALSE)

  cycling_prop <- aggregate(
    cycling ~ pt_bin,
    data = transform(plot_data, cycling = phase %in% c("S", "G2M")),
    FUN = mean
  )
  cycling_prop$pt_mid <- (cycling_prop$pt_bin - 0.5) / 10

  p3 <- ggplot(cycling_prop, aes(x = pt_mid, y = cycling * 100)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(color = "steelblue", size = 2) +
    geom_smooth(method = "loess", se = TRUE, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = "Proliferation Along Trajectory",
         x = "Pseudotime (normalized)", y = "% Cycling (S + G2M)")

  # Plot 4: UMAP colored by cell cycle (if UMAP available)
  if ("umap" %in% tolower(names(seurat_obj@reductions))) {
    umap_coords <- Seurat::Embeddings(seurat_obj, "umap")
    plot_data$UMAP_1 <- umap_coords[, 1]
    plot_data$UMAP_2 <- umap_coords[, 2]

    p4 <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = phase)) +
      geom_point(size = 0.3, alpha = 0.5) +
      scale_color_manual(values = cc_colors) +
      theme_minimal() +
      theme(legend.position = "bottom") +
      labs(title = "UMAP: Cell Cycle Phase", color = "Phase")
  } else {
    p4 <- ggplot() + theme_void() +
      labs(title = "UMAP not available")
  }

  # Combine
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(title = "Cell Cycle and Trajectory Relationship",
                    tag_levels = "A")

  return(combined)
}
