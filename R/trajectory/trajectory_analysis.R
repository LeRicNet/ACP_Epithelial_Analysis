# ==============================================================================
# R/trajectory/trajectory_analysis.R
# ==============================================================================
# Multi-method trajectory analysis with consensus determination
# Updated with HPC parallelization and comprehensive statistical evaluation
#
# Methods implemented:
#   1. Monocle3 - graph-based trajectory inference
#   2. Slingshot - principal curves through clusters
#   3. Diffusion Pseudotime - diffusion map-based ordering
#
# Statistical Framework:
#   - Per-method quality metrics
#   - Per-sample heterogeneity testing
#   - Bootstrap confidence intervals
#   - Permutation significance testing
# ==============================================================================

# ==============================================================================
# PARALLELIZATION SETUP
# ==============================================================================

#' Setup Parallel Backend
#'
#' Configures BiocParallel and future backends for HPC execution
#'
#' @param n_cores Number of cores to use (default: auto-detect)
#' @param type Backend type: "multicore" (Unix), "snow" (Windows), "serial"
#' @param memory_per_core Memory limit per core in GB (for SLURM)
#' @return List with configured backends
#' @export
setup_parallel_backend <- function(n_cores = NULL,
                                   type = "auto",
                                   memory_per_core = NULL) {


  message("\n--- Setting up parallel backend ---")

  # Auto-detect cores if not specified
  if (is.null(n_cores)) {
    # Check for SLURM environment
    slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)
    if (!is.na(slurm_cpus)) {
      n_cores <- as.integer(slurm_cpus)
      message(sprintf("  Detected SLURM environment: %d CPUs", n_cores))
    } else {
      n_cores <- max(1, parallel::detectCores() - 1)
      message(sprintf("  Auto-detected %d cores (reserving 1)", n_cores))
    }
  }

  # Ensure at least 1 core
  n_cores <- max(1, n_cores)

  # Determine backend type
  if (type == "auto") {
    type <- if (.Platform$OS.type == "unix") "multicore" else "snow"
  }

  backends <- list(
    n_cores = n_cores,
    type = type
  )

  # Setup BiocParallel
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    if (type == "multicore" && n_cores > 1) {
      backends$bioc <- BiocParallel::MulticoreParam(
        workers = n_cores,
        progressbar = TRUE
      )
    } else if (type == "snow" && n_cores > 1) {
      backends$bioc <- BiocParallel::SnowParam(
        workers = n_cores,
        progressbar = TRUE
      )
    } else {
      backends$bioc <- BiocParallel::SerialParam()
    }
    BiocParallel::register(backends$bioc)
    message(sprintf("  BiocParallel: %s with %d workers", type, n_cores))
  }

  # Setup future for future.apply
  if (requireNamespace("future", quietly = TRUE)) {
    if (type == "multicore" && n_cores > 1) {
      future::plan(future::multicore, workers = n_cores)
    } else if (type == "snow" && n_cores > 1) {
      future::plan(future::multisession, workers = n_cores)
    } else {
      future::plan(future::sequential)
    }

    # Set memory limits if specified
    if (!is.null(memory_per_core)) {
      options(future.globals.maxSize = memory_per_core * 1024^3)
    }

    backends$future_plan <- future::plan()
    message(sprintf("  future: %s", class(backends$future_plan)[1]))
  }

  # Store in options for later retrieval
  options(trajectory.parallel.backend = backends)

  return(backends)
}

#' Get Current Parallel Backend
#'
#' @return List with backend configuration
#' @export
get_parallel_backend <- function() {
  backend <- getOption("trajectory.parallel.backend")
  if (is.null(backend)) {
    backend <- list(n_cores = 1, type = "serial")
  }
  return(backend)
}

# ==============================================================================
# MAIN TRAJECTORY FUNCTIONS
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
#' @param n_cores Number of cores for parallel processing
#' @param bootstrap_ci Whether to compute bootstrap confidence intervals
#' @param n_bootstrap Number of bootstrap iterations
#' @return List containing results from each method and consensus pseudotime
#' @export
run_multi_method_trajectory <- function(seurat_obj,
                                        config = NULL,
                                        methods = c("monocle3", "slingshot", "diffusion"),
                                        subtype_column = "module_score_subtype",
                                        root_subtype = "Basal-like",
                                        reduction = "umap",
                                        n_cores = NULL,
                                        bootstrap_ci = FALSE,
                                        n_bootstrap = 100) {

  message("\n========================================")
  message("Multi-Method Trajectory Analysis")
  message("========================================\n")

  # Setup parallelization
  if (is.null(n_cores) && !is.null(config$reproducibility$n_cores)) {
    n_cores <- config$reproducibility$n_cores
  }
  if (!is.null(n_cores) && n_cores > 1) {
    setup_parallel_backend(n_cores)
  }

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
  method_quality <- list()

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

      # Calculate method quality metrics
      method_quality$monocle3 <- evaluate_method_quality(
        pseudotime = pt[common_cells],
        subtypes = subtypes[match(common_cells, colnames(seurat_obj))],
        method = "monocle3",
        cds = monocle_result$cds
      )

      message(sprintf("  Pseudotime range: [%.2f, %.2f]",
                      min(pt, na.rm = TRUE), max(pt, na.rm = TRUE)))
      message(sprintf("  Quality score: %.3f", method_quality$monocle3$overall_quality))

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

      # Calculate method quality metrics
      method_quality$slingshot <- evaluate_method_quality(
        pseudotime = pt,
        subtypes = subtypes[match(names(pt), colnames(seurat_obj))],
        method = "slingshot",
        sds = slingshot_result$sds
      )

      message(sprintf("  Pseudotime range: [%.2f, %.2f]",
                      min(pt, na.rm = TRUE), max(pt, na.rm = TRUE)))
      message(sprintf("  Lineages detected: %d", slingshot_result$n_lineages))
      message(sprintf("  Quality score: %.3f", method_quality$slingshot$overall_quality))

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
        reduction = reduction,
        parallel = n_cores > 1
      )

      results$diffusion <- diffusion_result

      # Extract pseudotime
      pt <- diffusion_result$pseudotime
      pseudotime_matrix[names(pt), "diffusion"] <- pt

      methods_used <- c(methods_used, "diffusion")

      # Calculate method quality metrics
      method_quality$diffusion <- evaluate_method_quality(
        pseudotime = pt,
        subtypes = subtypes[match(names(pt), colnames(seurat_obj))],
        method = "diffusion"
      )

      message(sprintf("  Pseudotime range: [%.2f, %.2f]",
                      min(pt, na.rm = TRUE), max(pt, na.rm = TRUE)))
      message(sprintf("  Quality score: %.3f", method_quality$diffusion$overall_quality))

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
  results$method_quality <- method_quality
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

  # ===========================================================================
  # Bootstrap Confidence Intervals (if requested)
  # ===========================================================================
  if (bootstrap_ci && length(methods_used) > 1) {
    message("\n--- Computing Bootstrap Confidence Intervals ---")

    bootstrap_results <- bootstrap_pseudotime_ci(
      seurat_obj = seurat_obj,
      methods = methods_used,
      config = config,
      subtype_column = subtype_column,
      root_subtype = root_subtype,
      reduction = reduction,
      n_bootstrap = n_bootstrap,
      parallel = n_cores > 1
    )

    results$bootstrap <- bootstrap_results
    results$seurat_obj$pseudotime_ci_lower <- bootstrap_results$ci_lower
    results$seurat_obj$pseudotime_ci_upper <- bootstrap_results$ci_upper
    results$seurat_obj$pseudotime_ci_width <- bootstrap_results$ci_width

    message(sprintf("  Mean CI width: %.3f", mean(bootstrap_results$ci_width, na.rm = TRUE)))
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

  # Calculate graph metrics
  graph_metrics <- list(
    n_nodes = igraph::vcount(graph),
    n_edges = igraph::ecount(graph),
    n_branch_points = n_branch_points,
    n_tips = sum(degree == 1),
    graph_diameter = igraph::diameter(graph)
  )

  return(list(
    cds = cds,
    n_branch_points = n_branch_points,
    graph_metrics = graph_metrics,
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

  # Calculate curve metrics
  curve_metrics <- list(
    n_lineages = n_lineages,
    lineage_clusters = lineages,
    curve_lengths = sapply(curves, function(c) sum(c$lambda))
  )

  return(list(
    sds = sds,
    pseudotime = pseudotime,
    pseudotime_matrix = pt_matrix,
    n_lineages = n_lineages,
    lineages = lineages,
    curves = curves,
    curve_metrics = curve_metrics,
    method = "slingshot"
  ))
}


#' Run Diffusion Pseudotime Analysis (Parallelized)
#'
#' Simple diffusion-based pseudotime using destiny package or manual calculation
#'
#' @param seurat_obj Seurat object
#' @param config Configuration list
#' @param subtype_column Column with cell types
#' @param root_subtype Root cell type
#' @param reduction Dimensionality reduction to use
#' @param parallel Whether to use parallel processing
#' @return List with diffusion results and pseudotime
#' @keywords internal
run_diffusion_pseudotime <- function(seurat_obj, config = NULL,
                                     subtype_column = "module_score_subtype",
                                     root_subtype = "Basal-like",
                                     reduction = "pca",
                                     parallel = FALSE) {

  # Get PCA or specified reduction
  if (reduction == "pca" && "pca" %in% names(seurat_obj@reductions)) {
    data_matrix <- Seurat::Embeddings(seurat_obj, "pca")[, 1:min(30, ncol(Seurat::Embeddings(seurat_obj, "pca")))]
  } else {
    data_matrix <- Seurat::Embeddings(seurat_obj, reduction)
  }

  n_cells <- nrow(data_matrix)
  message(sprintf("  Processing %d cells...", n_cells))

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

    # Compute distances - PARALLELIZED for large datasets
    if (parallel && n_cells > 1000 && requireNamespace("future.apply", quietly = TRUE)) {
      message("  Computing distances (parallel)...")

      # Split into chunks for parallel processing
      chunk_size <- ceiling(n_cells / get_parallel_backend()$n_cores)
      chunks <- split(1:n_cells, ceiling(1:n_cells / chunk_size))

      dist_list <- future.apply::future_lapply(chunks, function(idx) {
        as.matrix(dist(data_matrix))[idx, , drop = FALSE]
      }, future.seed = TRUE)

      dist_matrix <- do.call(rbind, dist_list)

    } else {
      message("  Computing distances (serial)...")
      dist_matrix <- as.matrix(dist(data_matrix))
    }

    # Build k-NN graph and transition probabilities
    message("  Building transition matrix...")
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
    message("  Computing eigenvectors...")
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
      transition_matrix = trans_prob,
      eigenvalues = eigen_result$values[1:10],
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


# ==============================================================================
# PER-SAMPLE ANALYSIS (PARALLELIZED)
# ==============================================================================

#' Run Per-Sample Multi-Method Trajectory (Parallelized)
#'
#' Runs trajectory analysis separately for each sample using multiple methods
#'
#' @param seurat_obj Seurat object
#' @param config Configuration list
#' @param methods Trajectory methods to use
#' @param sample_column Column containing sample IDs
#' @param subtype_column Column containing subtypes
#' @param parallel Whether to run samples in parallel
#' @return List with per-sample results and cross-sample concordance
#' @export
run_per_sample_multi_trajectory <- function(seurat_obj, config = NULL,
                                            methods = c("monocle3", "slingshot"),
                                            sample_column = "sample",
                                            subtype_column = "module_score_subtype",
                                            parallel = TRUE) {

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

  # Determine if we should run in parallel
  backend <- get_parallel_backend()
  use_parallel <- parallel && backend$n_cores > 1 && length(samples) > 1

  if (use_parallel) {
    message(sprintf("Running %d samples in parallel across %d cores",
                    length(samples), backend$n_cores))
  }

  # Function to process single sample
  process_sample <- function(sample_id) {
    # Subset
    cells_in_sample <- seurat_obj@meta.data[[sample_column]] == sample_id
    cells_in_sample[is.na(cells_in_sample)] <- FALSE
    n_cells <- sum(cells_in_sample)

    if (n_cells < min_cells) {
      return(list(
        sample = sample_id,
        n_cells = n_cells,
        skipped = TRUE,
        reason = "insufficient_cells"
      ))
    }

    # Subset Seurat
    seurat_sample <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_in_sample])

    # Run multi-method trajectory (without further parallelization to avoid nesting)
    tryCatch({
      sample_result <- run_multi_method_trajectory(
        seurat_sample,
        config = config,
        methods = methods,
        subtype_column = subtype_column,
        n_cores = 1  # No nested parallelization
      )

      # Calculate per-subtype mean pseudotime
      subtypes <- seurat_sample@meta.data[[subtype_column]]
      pt_by_subtype <- tapply(sample_result$consensus$pseudotime, subtypes, mean, na.rm = TRUE)
      subtype_rank <- rank(pt_by_subtype)

      list(
        sample = sample_id,
        n_cells = n_cells,
        skipped = FALSE,
        results = sample_result,
        pseudotime_by_subtype = pt_by_subtype,
        subtype_rank = subtype_rank,
        subtype_order = names(sort(pt_by_subtype)),
        methods_used = sample_result$methods_used,
        method_quality = sample_result$method_quality
      )

    }, error = function(e) {
      list(
        sample = sample_id,
        n_cells = n_cells,
        skipped = TRUE,
        reason = "trajectory_failed",
        error = e$message
      )
    })
  }

  # Run analysis
  if (use_parallel && requireNamespace("future.apply", quietly = TRUE)) {
    all_results <- future.apply::future_lapply(
      samples,
      process_sample,
      future.seed = TRUE
    )
    names(all_results) <- samples
  } else {
    # Sequential processing with progress
    all_results <- list()
    for (i in seq_along(samples)) {
      sample_id <- samples[i]
      message(sprintf("\n=== Sample %d/%d: %s ===", i, length(samples), sample_id))
      all_results[[sample_id]] <- process_sample(sample_id)

      if (!all_results[[sample_id]]$skipped) {
        message(sprintf("  Subtype order: %s",
                        paste(all_results[[sample_id]]$subtype_order, collapse = " → ")))
      }
    }
  }

  # Calculate cross-sample concordance
  message("\n=== Cross-Sample Concordance ===")
  concordance <- calculate_cross_sample_concordance(all_results)
  attr(all_results, "concordance") <- concordance

  # Calculate sample heterogeneity statistics
  message("\n=== Sample Heterogeneity Statistics ===")
  heterogeneity <- evaluate_sample_heterogeneity(all_results)
  attr(all_results, "heterogeneity") <- heterogeneity

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
  pairwise_results <- list()

  for (i in 1:(n_samples - 1)) {
    for (j in (i + 1):n_samples) {
      shared <- !is.na(rank_matrix[, i]) & !is.na(rank_matrix[, j])
      if (sum(shared) >= 3) {
        rho <- cor(rank_matrix[shared, i], rank_matrix[shared, j],
                   method = "spearman", use = "complete.obs")
        correlations <- c(correlations, rho)

        pair_name <- paste(colnames(rank_matrix)[i], colnames(rank_matrix)[j], sep = "_vs_")
        pairwise_results[[pair_name]] <- list(
          sample1 = colnames(rank_matrix)[i],
          sample2 = colnames(rank_matrix)[j],
          correlation = rho,
          n_shared = sum(shared)
        )
      }
    }
  }

  concordance <- list(
    rank_matrix = rank_matrix,
    mean_correlation = mean(correlations, na.rm = TRUE),
    sd_correlation = sd(correlations, na.rm = TRUE),
    min_correlation = min(correlations, na.rm = TRUE),
    max_correlation = max(correlations, na.rm = TRUE),
    all_correlations = correlations,
    pairwise_results = pairwise_results,
    n_comparisons = length(correlations),
    n_samples = n_samples
  )

  message(sprintf("  Samples analyzed: %d", n_samples))
  message(sprintf("  Mean Spearman rho: %.3f ± %.3f",
                  concordance$mean_correlation,
                  ifelse(is.na(concordance$sd_correlation), 0, concordance$sd_correlation)))

  return(concordance)
}


# ==============================================================================
# METHOD QUALITY EVALUATION
# ==============================================================================

#' Evaluate Method Quality
#'
#' Computes quality metrics for a single trajectory method
#'
#' @param pseudotime Pseudotime vector
#' @param subtypes Cell type labels
#' @param method Method name
#' @param cds Monocle3 cds object (optional)
#' @param sds Slingshot object (optional)
#' @return List with quality metrics
#' @export
evaluate_method_quality <- function(pseudotime, subtypes, method,
                                    cds = NULL, sds = NULL) {

  quality <- list(method = method)

  # 1. Pseudotime spread (should be well distributed)
  pt_range <- diff(range(pseudotime, na.rm = TRUE))
  pt_iqr <- IQR(pseudotime, na.rm = TRUE)
  quality$spread <- list(
    range = pt_range,
    iqr = pt_iqr,
    cv = sd(pseudotime, na.rm = TRUE) / mean(pseudotime, na.rm = TRUE)
  )

  # 2. Subtype separation (should show clear ordering)
  if (!is.null(subtypes) && length(unique(subtypes)) > 1) {
    pt_by_subtype <- tapply(pseudotime, subtypes, mean, na.rm = TRUE)

    # Between-group variance / total variance
    total_var <- var(pseudotime, na.rm = TRUE)
    between_var <- var(pt_by_subtype)
    quality$separation <- between_var / total_var

    # Kruskal-Wallis test
    kw_test <- kruskal.test(pseudotime ~ subtypes)
    quality$kruskal_wallis <- list(
      statistic = kw_test$statistic,
      p_value = kw_test$p.value
    )
  }

  # 3. Trajectory continuity (cells should have similar pseudotime to neighbors)
  # This requires additional computation, so we'll estimate from variance
  quality$continuity <- 1 - (sd(diff(sort(pseudotime)), na.rm = TRUE) / pt_range)

  # 4. Method-specific metrics
  if (!is.null(cds) && method == "monocle3") {
    # Graph coherence
    graph <- monocle3::principal_graph(cds)$UMAP
    degree_dist <- table(igraph::degree(graph))
    quality$graph_metrics <- list(
      mean_degree = mean(igraph::degree(graph)),
      n_branch_points = sum(igraph::degree(graph) > 2),
      n_tips = sum(igraph::degree(graph) == 1)
    )
  }

  if (!is.null(sds) && method == "slingshot") {
    quality$lineage_metrics <- list(
      n_lineages = length(slingshot::slingLineages(sds)),
      curve_coverage = mean(!is.na(slingshot::slingPseudotime(sds)))
    )
  }

  # 5. Overall quality score (weighted combination)
  score_components <- c()

  # Separation (higher is better)
  if (!is.null(quality$separation)) {
    score_components["separation"] <- min(quality$separation, 1)
  }

  # Continuity (higher is better)
  if (!is.null(quality$continuity) && !is.na(quality$continuity)) {
    score_components["continuity"] <- max(0, quality$continuity)
  }

  # Statistical significance (p < 0.05 = good)
  if (!is.null(quality$kruskal_wallis)) {
    score_components["significance"] <- ifelse(quality$kruskal_wallis$p_value < 0.05, 1, 0.5)
  }

  quality$overall_quality <- mean(score_components, na.rm = TRUE)
  quality$quality_components <- score_components

  return(quality)
}


# ==============================================================================
# SAMPLE HETEROGENEITY EVALUATION
# ==============================================================================

#' Evaluate Sample Heterogeneity
#'
#' Statistical tests for heterogeneity across samples
#'
#' @param results Per-sample trajectory results
#' @return List with heterogeneity statistics
#' @export
evaluate_sample_heterogeneity <- function(results) {

  valid_results <- results[sapply(results, function(x) !x$skipped)]

  if (length(valid_results) < 2) {
    return(list(
      n_samples = length(valid_results),
      heterogeneity_test = NA,
      message = "Insufficient samples for heterogeneity testing"
    ))
  }

  heterogeneity <- list(n_samples = length(valid_results))

  # 1. Concordance-based heterogeneity
  concordance <- attr(results, "concordance")
  if (!is.null(concordance) && !is.na(concordance$mean_correlation)) {
    # High correlation = low heterogeneity
    heterogeneity$concordance_based <- list(
      mean_correlation = concordance$mean_correlation,
      heterogeneity_index = 1 - concordance$mean_correlation,
      interpretation = ifelse(concordance$mean_correlation > 0.7, "Low heterogeneity",
                              ifelse(concordance$mean_correlation > 0.4, "Moderate heterogeneity",
                                     "High heterogeneity"))
    )
  }

  # 2. Ranking consistency test (Friedman test if >= 3 samples)
  if (length(valid_results) >= 3) {
    # Build data for Friedman test
    rank_matrix <- concordance$rank_matrix

    # Remove subtypes with missing values
    complete_subtypes <- rowSums(is.na(rank_matrix)) == 0
    if (sum(complete_subtypes) >= 3) {
      rank_complete <- rank_matrix[complete_subtypes, ]

      tryCatch({
        # Reshape for Friedman test
        friedman_data <- as.data.frame(t(rank_complete))
        friedman_result <- friedman.test(as.matrix(friedman_data))

        heterogeneity$friedman_test <- list(
          statistic = friedman_result$statistic,
          p_value = friedman_result$p.value,
          significant = friedman_result$p.value < 0.05,
          interpretation = ifelse(friedman_result$p.value < 0.05,
                                  "Significant ordering difference across samples",
                                  "Consistent ordering across samples")
        )
      }, error = function(e) {
        heterogeneity$friedman_test <- list(error = e$message)
      })
    }
  }

  # 3. Per-sample quality comparison
  sample_qualities <- sapply(valid_results, function(x) {
    if (!is.null(x$method_quality)) {
      mean(sapply(x$method_quality, function(m) m$overall_quality), na.rm = TRUE)
    } else {
      NA
    }
  })

  heterogeneity$quality_comparison <- list(
    mean_quality = mean(sample_qualities, na.rm = TRUE),
    sd_quality = sd(sample_qualities, na.rm = TRUE),
    min_quality = min(sample_qualities, na.rm = TRUE),
    max_quality = max(sample_qualities, na.rm = TRUE),
    per_sample = sample_qualities
  )

  # 4. Outlier sample detection
  if (length(sample_qualities) >= 3) {
    # Use IQR-based outlier detection
    q1 <- quantile(sample_qualities, 0.25, na.rm = TRUE)
    q3 <- quantile(sample_qualities, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower_bound <- q1 - 1.5 * iqr
    upper_bound <- q3 + 1.5 * iqr

    outliers <- names(sample_qualities)[sample_qualities < lower_bound |
                                          sample_qualities > upper_bound]

    heterogeneity$outlier_samples <- list(
      outliers = outliers,
      n_outliers = length(outliers),
      bounds = c(lower = lower_bound, upper = upper_bound)
    )
  }

  # 5. Overall heterogeneity summary
  heterogeneity$summary <- list(
    is_heterogeneous = FALSE,
    reasons = c()
  )

  if (!is.null(heterogeneity$concordance_based) &&
      heterogeneity$concordance_based$mean_correlation < 0.5) {
    heterogeneity$summary$is_heterogeneous <- TRUE
    heterogeneity$summary$reasons <- c(heterogeneity$summary$reasons,
                                       "Low cross-sample concordance")
  }

  if (!is.null(heterogeneity$friedman_test) &&
      !is.null(heterogeneity$friedman_test$significant) &&
      heterogeneity$friedman_test$significant) {
    heterogeneity$summary$is_heterogeneous <- TRUE
    heterogeneity$summary$reasons <- c(heterogeneity$summary$reasons,
                                       "Significant Friedman test")
  }

  if (!is.null(heterogeneity$outlier_samples) &&
      heterogeneity$outlier_samples$n_outliers > 0) {
    heterogeneity$summary$is_heterogeneous <- TRUE
    heterogeneity$summary$reasons <- c(heterogeneity$summary$reasons,
                                       paste("Outlier samples:",
                                             paste(heterogeneity$outlier_samples$outliers,
                                                   collapse = ", ")))
  }

  message(sprintf("  Heterogeneity: %s",
                  ifelse(heterogeneity$summary$is_heterogeneous, "DETECTED", "NOT DETECTED")))
  if (length(heterogeneity$summary$reasons) > 0) {
    for (reason in heterogeneity$summary$reasons) {
      message(sprintf("    - %s", reason))
    }
  }

  return(heterogeneity)
}


# ==============================================================================
# BOOTSTRAP CONFIDENCE INTERVALS
# ==============================================================================

#' Bootstrap Pseudotime Confidence Intervals
#'
#' Computes bootstrap CIs for consensus pseudotime
#'
#' @param seurat_obj Seurat object
#' @param methods Methods to use
#' @param config Configuration list
#' @param subtype_column Subtype column name
#' @param root_subtype Root subtype
#' @param reduction Reduction to use
#' @param n_bootstrap Number of bootstrap iterations
#' @param ci_level Confidence interval level (default 0.95)
#' @param parallel Whether to run in parallel
#' @return List with confidence intervals
#' @export
bootstrap_pseudotime_ci <- function(seurat_obj,
                                    methods,
                                    config = NULL,
                                    subtype_column = "module_score_subtype",
                                    root_subtype = "Basal-like",
                                    reduction = "umap",
                                    n_bootstrap = 100,
                                    ci_level = 0.95,
                                    parallel = TRUE) {

  n_cells <- ncol(seurat_obj)
  alpha <- 1 - ci_level

  message(sprintf("  Running %d bootstrap iterations...", n_bootstrap))

  # Function for single bootstrap iteration
  run_bootstrap_iter <- function(iter) {
    # Sample cells with replacement
    boot_idx <- sample(1:n_cells, n_cells, replace = TRUE)
    boot_cells <- colnames(seurat_obj)[boot_idx]

    # Create bootstrap Seurat object
    boot_seurat <- subset(seurat_obj, cells = unique(boot_cells))

    # Run trajectory (lightweight - single method for speed)
    tryCatch({
      # Use fastest method for bootstrap
      boot_result <- run_multi_method_trajectory(
        boot_seurat,
        config = config,
        methods = methods[1],  # Use first method only for speed
        subtype_column = subtype_column,
        root_subtype = root_subtype,
        n_cores = 1,
        bootstrap_ci = FALSE
      )

      # Return pseudotime for original cells
      pt <- boot_result$consensus$pseudotime
      names(pt) <- colnames(boot_seurat)

      # Map back to original cell indices
      pt_full <- rep(NA, n_cells)
      names(pt_full) <- colnames(seurat_obj)

      # Average for cells that appear multiple times
      for (cell in unique(boot_cells)) {
        if (cell %in% names(pt)) {
          pt_full[cell] <- pt[cell]
        }
      }

      pt_full

    }, error = function(e) {
      rep(NA, n_cells)
    })
  }

  # Run bootstrap iterations
  backend <- get_parallel_backend()

  if (parallel && backend$n_cores > 1 && requireNamespace("future.apply", quietly = TRUE)) {
    message(sprintf("  Using %d cores for bootstrap...", backend$n_cores))

    boot_results <- future.apply::future_lapply(
      1:n_bootstrap,
      run_bootstrap_iter,
      future.seed = TRUE
    )
  } else {
    boot_results <- lapply(1:n_bootstrap, function(i) {
      if (i %% 10 == 0) message(sprintf("    Iteration %d/%d", i, n_bootstrap))
      run_bootstrap_iter(i)
    })
  }

  # Combine results into matrix
  boot_matrix <- do.call(cbind, boot_results)
  rownames(boot_matrix) <- colnames(seurat_obj)

  # Calculate confidence intervals
  ci_lower <- apply(boot_matrix, 1, quantile, probs = alpha/2, na.rm = TRUE)
  ci_upper <- apply(boot_matrix, 1, quantile, probs = 1 - alpha/2, na.rm = TRUE)
  ci_width <- ci_upper - ci_lower

  # Calculate bootstrap standard error
  boot_se <- apply(boot_matrix, 1, sd, na.rm = TRUE)

  # Calculate number of successful iterations per cell
  n_valid <- rowSums(!is.na(boot_matrix))

  message(sprintf("  Mean CI width: %.3f", mean(ci_width, na.rm = TRUE)))
  message(sprintf("  Mean bootstrap SE: %.3f", mean(boot_se, na.rm = TRUE)))

  return(list(
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    ci_width = ci_width,
    boot_se = boot_se,
    n_valid = n_valid,
    boot_matrix = boot_matrix,
    ci_level = ci_level,
    n_bootstrap = n_bootstrap
  ))
}


# ==============================================================================
# PERMUTATION SIGNIFICANCE TESTING
# ==============================================================================

#' Permutation Test for Trajectory Significance
#'
#' Tests whether the observed trajectory ordering is significantly different from random
#'
#' @param seurat_obj Seurat object with pseudotime
#' @param subtype_column Column with subtypes
#' @param expected_order Expected subtype ordering
#' @param n_permutations Number of permutations
#' @param parallel Whether to run in parallel
#' @return List with permutation test results
#' @export
permutation_test_trajectory <- function(seurat_obj,
                                        subtype_column = "module_score_subtype",
                                        expected_order = NULL,
                                        n_permutations = 1000,
                                        parallel = TRUE) {

  message("\n--- Permutation Test for Trajectory Significance ---")

  pseudotime <- seurat_obj$pseudotime_consensus
  subtypes <- seurat_obj@meta.data[[subtype_column]]

  if (is.null(expected_order)) {
    expected_order <- c("Basal-like", "Transit-Amplifying", "Intermediate", "Specialized")
  }

  # Calculate observed statistic
  # Use Spearman correlation between mean pseudotime rank and expected rank
  observed_mean_pt <- tapply(pseudotime, subtypes, mean, na.rm = TRUE)

  # Filter to subtypes in expected order
  common_subtypes <- intersect(names(observed_mean_pt), expected_order)

  if (length(common_subtypes) < 3) {
    warning("Too few common subtypes for permutation test")
    return(list(
      p_value = NA,
      message = "Insufficient common subtypes"
    ))
  }

  # Expected ranks
  expected_ranks <- match(common_subtypes, expected_order)

  # Observed correlation
  observed_cor <- cor(rank(observed_mean_pt[common_subtypes]), expected_ranks,
                      method = "spearman")

  message(sprintf("  Observed correlation: %.3f", observed_cor))
  message(sprintf("  Running %d permutations...", n_permutations))

  # Permutation function
  run_permutation <- function(i) {
    # Permute pseudotime values
    perm_pt <- sample(pseudotime)

    # Calculate mean pseudotime by subtype
    perm_mean_pt <- tapply(perm_pt, subtypes, mean, na.rm = TRUE)

    # Calculate correlation with expected order
    cor(rank(perm_mean_pt[common_subtypes]), expected_ranks, method = "spearman")
  }

  # Run permutations
  backend <- get_parallel_backend()

  if (parallel && backend$n_cores > 1 && requireNamespace("future.apply", quietly = TRUE)) {
    perm_cors <- future.apply::future_sapply(
      1:n_permutations,
      run_permutation,
      future.seed = TRUE
    )
  } else {
    perm_cors <- sapply(1:n_permutations, run_permutation)
  }

  # Calculate p-value (one-tailed: observed should be higher than random)
  p_value <- mean(perm_cors >= observed_cor)

  # Also calculate two-tailed p-value
  p_value_two_tailed <- mean(abs(perm_cors) >= abs(observed_cor))

  message(sprintf("  P-value (one-tailed): %.4f", p_value))
  message(sprintf("  P-value (two-tailed): %.4f", p_value_two_tailed))

  # Effect size (how many SDs above random)
  effect_size <- (observed_cor - mean(perm_cors)) / sd(perm_cors)
  message(sprintf("  Effect size (Cohen's d): %.2f", effect_size))

  return(list(
    observed_correlation = observed_cor,
    permutation_correlations = perm_cors,
    p_value = p_value,
    p_value_two_tailed = p_value_two_tailed,
    effect_size = effect_size,
    n_permutations = n_permutations,
    significant = p_value < 0.05,
    interpretation = ifelse(p_value < 0.001, "Highly significant trajectory",
                            ifelse(p_value < 0.01, "Very significant trajectory",
                                   ifelse(p_value < 0.05, "Significant trajectory",
                                          "Non-significant trajectory")))
  ))
}


#' Permutation Test for Method Concordance
#'
#' Tests whether method agreement is significantly better than chance
#'
#' @param pseudotime_matrix Matrix of pseudotime values (cells x methods)
#' @param n_permutations Number of permutations
#' @param parallel Whether to run in parallel
#' @return List with permutation test results
#' @export
permutation_test_method_concordance <- function(pseudotime_matrix,
                                                n_permutations = 1000,
                                                parallel = TRUE) {

  message("\n--- Permutation Test for Method Concordance ---")

  n_methods <- ncol(pseudotime_matrix)

  if (n_methods < 2) {
    return(list(p_value = NA, message = "Need at least 2 methods"))
  }

  # Observed mean pairwise correlation
  cor_matrix <- cor(pseudotime_matrix, method = "spearman", use = "pairwise.complete.obs")
  observed_mean_cor <- mean(cor_matrix[upper.tri(cor_matrix)])

  message(sprintf("  Observed mean correlation: %.3f", observed_mean_cor))

  # Permutation function
  run_permutation <- function(i) {
    # Permute each method's pseudotime independently
    perm_matrix <- apply(pseudotime_matrix, 2, sample)

    # Calculate correlation
    perm_cor <- cor(perm_matrix, method = "spearman", use = "pairwise.complete.obs")
    mean(perm_cor[upper.tri(perm_cor)])
  }

  # Run permutations
  backend <- get_parallel_backend()

  if (parallel && backend$n_cores > 1 && requireNamespace("future.apply", quietly = TRUE)) {
    perm_cors <- future.apply::future_sapply(
      1:n_permutations,
      run_permutation,
      future.seed = TRUE
    )
  } else {
    perm_cors <- sapply(1:n_permutations, run_permutation)
  }

  # Calculate p-value
  p_value <- mean(perm_cors >= observed_mean_cor)

  message(sprintf("  P-value: %.4f", p_value))

  return(list(
    observed_correlation = observed_mean_cor,
    permutation_correlations = perm_cors,
    p_value = p_value,
    significant = p_value < 0.05,
    n_permutations = n_permutations
  ))
}


# ==============================================================================
# TRAJECTORY VALIDATION
# ==============================================================================

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


# ==============================================================================
# COMPREHENSIVE STATISTICAL EVALUATION
# ==============================================================================

#' Evaluate Trajectory Analysis Statistics
#'
#' Computes comprehensive statistics to evaluate trajectory inference quality
#'
#' @param seurat_obj Seurat object with pseudotime
#' @param results Results list from run_multi_method_trajectory
#' @param config Configuration list
#' @param differentiation_markers Named list of markers
#' @param expected_order Expected order of subtypes
#' @param run_permutation Whether to run permutation tests
#' @param n_permutations Number of permutations
#' @return List with statistical evaluation results
#' @export
evaluate_trajectory_statistics <- function(seurat_obj,
                                           results,
                                           config = NULL,
                                           differentiation_markers = NULL,
                                           expected_order = NULL,
                                           run_permutation = TRUE,
                                           n_permutations = 1000) {

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
      early = c("SOX9", "KRT5", "KRT14", "TP63", "ITGA6", "ITGB1"),
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
    upper_tri <- cor_matrix[upper.tri(cor_matrix)]
    stats$mean_correlation <- mean(upper_tri, na.rm = TRUE)
    message(sprintf("  Mean pairwise correlation: ρ = %.3f", stats$mean_correlation))

    # Permutation test for method concordance
    if (run_permutation) {
      concordance_perm <- permutation_test_method_concordance(
        pt_matrix[, methods_used],
        n_permutations = n_permutations,
        parallel = get_parallel_backend()$n_cores > 1
      )
      stats$method_concordance_permutation <- concordance_perm
    }

  } else {
    message("  Only one method - skipping agreement statistics")
    stats$mean_correlation <- NA
  }

  # ===========================================================================
  # 2. Per-Method Quality Scores
  # ===========================================================================
  message("\n--- Per-Method Quality ---")

  if (!is.null(results$method_quality)) {
    stats$method_quality <- results$method_quality

    for (method in names(results$method_quality)) {
      mq <- results$method_quality[[method]]
      message(sprintf("  %s: overall=%.3f, separation=%.3f",
                      method,
                      mq$overall_quality,
                      ifelse(is.null(mq$separation), NA, mq$separation)))
    }
  }

  # ===========================================================================
  # 3. Cell Cycle Confounding Assessment
  # ===========================================================================
  message("\n--- Cell Cycle Confounding ---")

  has_cell_cycle <- "cc_consensus" %in% colnames(seurat_obj@meta.data) ||
    "Phase" %in% colnames(seurat_obj@meta.data)

  if (has_cell_cycle) {
    pseudotime <- seurat_obj$pseudotime_consensus

    phase_col <- ifelse("cc_consensus" %in% colnames(seurat_obj@meta.data),
                        "cc_consensus", "Phase")

    phase_numeric <- as.numeric(factor(seurat_obj@meta.data[[phase_col]],
                                       levels = c("G1", "S", "G2M")))

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
    pt_bins <- cut(pseudotime, breaks = 5, labels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%"))
    cycling_by_bin <- tapply(seurat_obj@meta.data[[phase_col]] %in% c("S", "G2M"),
                             pt_bins, mean, na.rm = TRUE) * 100

    stats$cycling_along_trajectory <- data.frame(
      pseudotime_bin = names(cycling_by_bin),
      cycling_percent = as.numeric(cycling_by_bin)
    )

  } else {
    message("  Cell cycle data not available - skipping")
    stats$cell_cycle_correlation <- NULL
  }

  # ===========================================================================
  # 4. Biological Validation - Marker Correlations
  # ===========================================================================
  message("\n--- Biological Validation (Marker Trends) ---")

  pseudotime <- seurat_obj$pseudotime_consensus
  available_early <- intersect(differentiation_markers$early, rownames(seurat_obj))
  available_late <- intersect(differentiation_markers$late, rownames(seurat_obj))

  message(sprintf("  Early markers available: %d/%d",
                  length(available_early), length(differentiation_markers$early)))
  message(sprintf("  Late markers available: %d/%d",
                  length(available_late), length(differentiation_markers$late)))

  if (length(available_early) > 0 || length(available_late) > 0) {
    all_markers <- c(available_early, available_late)
    expr_data <- as.matrix(Seurat::GetAssayData(seurat_obj, layer = "data")[all_markers, , drop = FALSE])

    marker_results <- lapply(all_markers, function(marker) {
      expected_direction <- ifelse(marker %in% available_early, "negative", "positive")

      rho <- cor(pseudotime, expr_data[marker, ],
                 use = "pairwise.complete.obs", method = "spearman")
      test <- cor.test(pseudotime, expr_data[marker, ], method = "spearman")

      correct_direction <- ifelse(expected_direction == "negative", rho < 0, rho > 0)

      list(
        marker = marker,
        rho = rho,
        p_value = test$p.value,
        expected = expected_direction,
        correct_direction = correct_direction,
        significant = test$p.value < 0.05
      )
    })

    marker_df <- do.call(rbind, lapply(marker_results, as.data.frame))
    stats$marker_validation <- marker_df

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
  }

  # ===========================================================================
  # 5. Trajectory Order Validation
  # ===========================================================================
  message("\n--- Trajectory Order Validation ---")

  if (!is.null(expected_order)) {
    subtype_col <- ifelse("module_score_subtype" %in% colnames(seurat_obj@meta.data),
                          "module_score_subtype", "subtype")
    subtypes <- seurat_obj@meta.data[[subtype_col]]
    mean_pt_by_subtype <- tapply(pseudotime, subtypes, mean, na.rm = TRUE)
    observed_order <- names(sort(mean_pt_by_subtype))
    observed_in_expected <- observed_order[observed_order %in% expected_order]
    expected_filtered <- expected_order[expected_order %in% observed_order]

    message(sprintf("  Expected order: %s", paste(expected_filtered, collapse = " → ")))
    message(sprintf("  Observed order: %s", paste(observed_in_expected, collapse = " → ")))

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
    }

    stats$mean_pseudotime_by_subtype <- data.frame(
      subtype = names(mean_pt_by_subtype),
      mean_pseudotime = as.numeric(mean_pt_by_subtype),
      sd_pseudotime = as.numeric(tapply(pseudotime, subtypes, sd, na.rm = TRUE)),
      n_cells = as.numeric(table(subtypes)[names(mean_pt_by_subtype)])
    )
  }

  # ===========================================================================
  # 6. Permutation Test for Trajectory Significance
  # ===========================================================================
  if (run_permutation && !is.null(expected_order)) {
    perm_test <- permutation_test_trajectory(
      seurat_obj,
      subtype_column = ifelse("module_score_subtype" %in% colnames(seurat_obj@meta.data),
                              "module_score_subtype", "subtype"),
      expected_order = expected_order,
      n_permutations = n_permutations,
      parallel = get_parallel_backend()$n_cores > 1
    )
    stats$permutation_test <- perm_test
  }

  # ===========================================================================
  # 7. Bootstrap Results Summary (if available)
  # ===========================================================================
  if (!is.null(results$bootstrap)) {
    message("\n--- Bootstrap Confidence Intervals ---")

    stats$bootstrap_summary <- list(
      mean_ci_width = mean(results$bootstrap$ci_width, na.rm = TRUE),
      median_ci_width = median(results$bootstrap$ci_width, na.rm = TRUE),
      mean_se = mean(results$bootstrap$boot_se, na.rm = TRUE),
      coverage = mean(results$bootstrap$n_valid >= results$bootstrap$n_bootstrap * 0.9)
    )

    message(sprintf("  Mean CI width: %.3f", stats$bootstrap_summary$mean_ci_width))
    message(sprintf("  Mean bootstrap SE: %.3f", stats$bootstrap_summary$mean_se))
  }

  # ===========================================================================
  # 8. Overall Quality Score
  # ===========================================================================
  message("\n--- Overall Quality Score ---")

  quality_components <- c()

  if (!is.na(stats$mean_correlation)) {
    agreement_score <- min(stats$mean_correlation / 0.7, 1)
    quality_components["method_agreement"] <- agreement_score * 0.20
  }

  if (!is.null(stats$marker_validation_summary)) {
    marker_score <- stats$marker_validation_summary$percent_correct / 100
    quality_components["marker_validation"] <- marker_score * 0.25
  }

  if (!is.null(stats$order_validation)) {
    order_score <- max(0, (stats$order_validation$rank_correlation + 1) / 2)
    quality_components["order_preservation"] <- order_score * 0.20
  }

  if (!is.null(stats$cell_cycle_correlation)) {
    cc_score <- 1 - min(abs(stats$cell_cycle_correlation$rho), 1)
    quality_components["cell_cycle"] <- cc_score * 0.10
  }

  if (!is.null(stats$permutation_test)) {
    perm_score <- ifelse(stats$permutation_test$p_value < 0.05, 1, 0.5)
    quality_components["significance"] <- perm_score * 0.15
  }

  if (!is.null(results$consensus$confidence)) {
    conf_score <- mean(results$consensus$confidence, na.rm = TRUE)
    quality_components["confidence"] <- conf_score * 0.10
  }

  overall_score <- sum(quality_components, na.rm = TRUE)
  max_possible <- sum(c(0.20, 0.25, 0.20, 0.10, 0.15, 0.10)[
    c("method_agreement", "marker_validation", "order_preservation",
      "cell_cycle", "significance", "confidence") %in% names(quality_components)
  ])

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


# ==============================================================================
# REPORTING FUNCTIONS
# ==============================================================================

#' Generate Trajectory Statistics Report
#'
#' Creates a formatted markdown report
#'
#' @param stats_results Results from evaluate_trajectory_statistics
#' @param output_path Path to save the report
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

    # Add permutation test result if available
    if (!is.null(stats_results$method_concordance_permutation)) {
      perm <- stats_results$method_concordance_permutation
      report <- c(report, sprintf("**Concordance permutation test: p = %.4f (%s)**",
                                  perm$p_value,
                                  ifelse(perm$significant, "significant", "not significant")), "")
    }
  }

  # Per-method quality
  if (!is.null(stats_results$method_quality)) {
    report <- c(report, "## Per-Method Quality", "")
    report <- c(report, "| Method | Overall Quality | Separation | Continuity |")
    report <- c(report, "|--------|-----------------|------------|------------|")

    for (method in names(stats_results$method_quality)) {
      mq <- stats_results$method_quality[[method]]
      report <- c(report, sprintf("| %s | %.3f | %.3f | %.3f |",
                                  method,
                                  mq$overall_quality,
                                  ifelse(is.null(mq$separation), NA, mq$separation),
                                  ifelse(is.null(mq$continuity), NA, mq$continuity)))
    }
    report <- c(report, "")
  }

  # Permutation test results
  if (!is.null(stats_results$permutation_test)) {
    perm <- stats_results$permutation_test
    report <- c(report, "## Trajectory Significance (Permutation Test)", "",
                sprintf("- Observed correlation: %.3f", perm$observed_correlation),
                sprintf("- P-value: %.4f", perm$p_value),
                sprintf("- Effect size (Cohen's d): %.2f", perm$effect_size),
                sprintf("- Interpretation: **%s**", perm$interpretation), "")
  }

  # Cell cycle section
  if (!is.null(stats_results$cell_cycle_correlation)) {
    cc <- stats_results$cell_cycle_correlation
    report <- c(report, "## Cell Cycle Confounding", "",
                sprintf("- Pseudotime ~ Phase correlation: ρ = %.3f (p = %.2e)", cc$rho, cc$p_value),
                sprintf("- Interpretation: **%s**", cc$interpretation), "")
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

  # Bootstrap summary
  if (!is.null(stats_results$bootstrap_summary)) {
    bs <- stats_results$bootstrap_summary
    report <- c(report, "## Bootstrap Confidence Intervals", "",
                sprintf("- Mean CI width: %.3f", bs$mean_ci_width),
                sprintf("- Mean bootstrap SE: %.3f", bs$mean_se),
                sprintf("- Coverage: %.1f%%", bs$coverage * 100), "")
  }

  # Pseudotime by subtype
  if (!is.null(stats_results$mean_pseudotime_by_subtype)) {
    report <- c(report, "## Mean Pseudotime by Subtype", "")
    report <- c(report, "| Subtype | Mean PT | SD | N Cells |")
    report <- c(report, "|---------|---------|-------|---------|")

    pt <- stats_results$mean_pseudotime_by_subtype
    pt <- pt[order(pt$mean_pseudotime), ]
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


# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

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


#' Plot Cell Cycle Along Trajectory
#'
#' Visualizes the relationship between pseudotime and cell cycle
#'
#' @param seurat_obj Seurat object with pseudotime and cell cycle
#' @param config Configuration list
#' @return ggplot object
#' @export
plot_cell_cycle_trajectory <- function(seurat_obj, config = NULL) {

  require(ggplot2)
  require(patchwork)

  # Determine cell cycle column
  cc_col <- ifelse("cc_consensus" %in% colnames(seurat_obj@meta.data),
                   "cc_consensus", "Phase")

  if (!cc_col %in% colnames(seurat_obj@meta.data) ||
      !"pseudotime_consensus" %in% colnames(seurat_obj@meta.data)) {
    stop("Requires both pseudotime_consensus and cell cycle in metadata")
  }

  # Get colors
  if (!is.null(config) && !is.null(config$visualization$colors$cell_cycle)) {
    cc_colors <- unlist(config$visualization$colors$cell_cycle)
  } else {
    cc_colors <- c(G1 = "#F8766D", S = "#00BA38", G2M = "#619CFF")
  }

  plot_data <- data.frame(
    pseudotime = seurat_obj$pseudotime_consensus,
    phase = seurat_obj@meta.data[[cc_col]]
  )

  if ("pseudotime_confidence" %in% colnames(seurat_obj@meta.data)) {
    plot_data$confidence <- seurat_obj$pseudotime_confidence
  }

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

  # Plot 4: UMAP colored by cell cycle
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
