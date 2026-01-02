# ==============================================================================
# R/trajectory/trajectory_analysis.R
# ==============================================================================
# Multi-method trajectory analysis with consensus determination
# ==============================================================================

# ==============================================================================
# PARALLELIZATION SETUP
# ==============================================================================

#' Setup Parallel Backend
#' @export
setup_parallel_backend <- function(n_cores = NULL, type = "auto", memory_per_core = NULL) {

  message("\n--- Setting up parallel backend ---")

  if (is.null(n_cores)) {
    slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)
    if (!is.na(slurm_cpus)) {
      n_cores <- as.integer(slurm_cpus)
      message(sprintf("  Detected SLURM environment: %d CPUs", n_cores))
    } else {
      n_cores <- max(1, parallel::detectCores() - 1)
      message(sprintf("  Auto-detected %d cores (reserving 1)", n_cores))
    }
  }

  n_cores <- max(1, n_cores)

  if (type == "auto") {
    type <- if (.Platform$OS.type == "unix") "multicore" else "snow"
  }

  backends <- list(n_cores = n_cores, type = type)

  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    if (type == "multicore" && n_cores > 1) {
      backends$bioc <- BiocParallel::MulticoreParam(workers = n_cores, progressbar = TRUE)
    } else if (type == "snow" && n_cores > 1) {
      backends$bioc <- BiocParallel::SnowParam(workers = n_cores, progressbar = TRUE)
    } else {
      backends$bioc <- BiocParallel::SerialParam()
    }
    BiocParallel::register(backends$bioc)
    message(sprintf("  BiocParallel: %s with %d workers", type, n_cores))
  }

  if (requireNamespace("future", quietly = TRUE)) {
    if (type == "multicore" && n_cores > 1) {
      future::plan(future::multicore, workers = n_cores)
    } else if (type == "snow" && n_cores > 1) {
      future::plan(future::multisession, workers = n_cores)
    } else {
      future::plan(future::sequential)
    }

    if (!is.null(memory_per_core)) {
      options(future.globals.maxSize = memory_per_core * 1024^3)
    }

    backends$future_plan <- future::plan()
    message(sprintf("  future: %s", class(backends$future_plan)[1]))
  }

  options(trajectory.parallel.backend = backends)
  return(backends)
}

#' Get Current Parallel Backend
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
#' @export
run_multi_method_trajectory <- function(seurat_obj,
                                        config = NULL,
                                        methods = c("monocle3", "slingshot"),
                                        subtype_column = "module_score_subtype",
                                        root_subtype = "Basal-like",
                                        reduction = "umap",
                                        n_cores = NULL,
                                        bootstrap_ci = FALSE,
                                        n_bootstrap = 100) {

  message("\n========================================")
  message("Multi-Method Trajectory Analysis")
  message("========================================\n")

  if (is.null(n_cores) && !is.null(config$reproducibility$n_cores)) {
    n_cores <- config$reproducibility$n_cores
  }
  if (!is.null(n_cores) && n_cores > 1) {
    setup_parallel_backend(n_cores)
  }

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

  if (!subtype_column %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Subtype column '%s' not found in metadata", subtype_column))
  }

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

  results <- list()
  pseudotime_matrix <- matrix(NA, nrow = ncol(seurat_obj), ncol = length(methods))
  rownames(pseudotime_matrix) <- colnames(seurat_obj)
  colnames(pseudotime_matrix) <- methods
  methods_used <- c()
  method_quality <- list()

  subtypes <- seurat_obj@meta.data[[subtype_column]]

  # Method 1: Monocle3
  if ("monocle3" %in% methods) {
    message("\n--- Method 1: Monocle3 ---")

    result <- tryCatch({
      monocle_result <- run_monocle3_trajectory(
        seurat_obj, config = config, subtype_column = subtype_column, root_subtype = root_subtype
      )

      results$monocle3 <- monocle_result
      pt <- monocle3::pseudotime(monocle_result$cds)
      common_cells <- intersect(names(pt), colnames(seurat_obj))
      pseudotime_matrix[common_cells, "monocle3"] <- pt[common_cells]
      methods_used <- c(methods_used, "monocle3")

      method_quality$monocle3 <- evaluate_method_quality(
        pseudotime = pt[common_cells],
        subtypes = subtypes[match(common_cells, colnames(seurat_obj))],
        method = "monocle3",
        cds = monocle_result$cds
      )

      message(sprintf("  Pseudotime range: [%.2f, %.2f]", min(pt, na.rm = TRUE), max(pt, na.rm = TRUE)))
      message(sprintf("  Quality score: %.3f", method_quality$monocle3$overall_quality))
      TRUE
    }, error = function(e) {
      warning(paste("  Monocle3 failed:", e$message))
      FALSE
    })
  }

  # Method 2: Slingshot
  if ("slingshot" %in% methods) {
    message("\n--- Method 2: Slingshot ---")

    result <- tryCatch({
      slingshot_result <- run_slingshot_trajectory(
        seurat_obj, config = config, subtype_column = subtype_column,
        root_subtype = root_subtype, reduction = reduction
      )

      results$slingshot <- slingshot_result
      pt <- slingshot_result$pseudotime
      pseudotime_matrix[names(pt), "slingshot"] <- pt
      methods_used <- c(methods_used, "slingshot")

      method_quality$slingshot <- evaluate_method_quality(
        pseudotime = pt,
        subtypes = subtypes[match(names(pt), colnames(seurat_obj))],
        method = "slingshot",
        sds = slingshot_result$sds
      )

      message(sprintf("  Pseudotime range: [%.2f, %.2f]", min(pt, na.rm = TRUE), max(pt, na.rm = TRUE)))
      message(sprintf("  Lineages detected: %d", slingshot_result$n_lineages))
      message(sprintf("  Quality score: %.3f", method_quality$slingshot$overall_quality))
      TRUE
    }, error = function(e) {
      warning(paste("  Slingshot failed:", e$message))
      FALSE
    })
  }

  if (length(methods_used) == 0) {
    stop("All trajectory methods failed!")
  }

  results$methods_used <- methods_used
  results$pseudotime_matrix <- pseudotime_matrix[, methods_used, drop = FALSE]
  results$method_quality <- method_quality

  # Calculate consensus
  message("\n--- Calculating Consensus Pseudotime ---")
  consensus <- calculate_consensus_pseudotime(results$pseudotime_matrix)
  results$consensus <- consensus

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
#' @keywords internal
run_monocle3_trajectory <- function(seurat_obj, config = NULL,
                                    subtype_column = "module_score_subtype",
                                    root_subtype = "Basal-like") {

  if (!requireNamespace("monocle3", quietly = TRUE)) {
    stop("monocle3 package required")
  }
  if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
    stop("SeuratWrappers package required")
  }

  library(monocle3)

  message("  Converting to cell_data_set...")
  cds <- SeuratWrappers::as.cell_data_set(seurat_obj)
  SummarizedExperiment::colData(cds)$subtype <- seurat_obj@meta.data[[subtype_column]]

  n_pcs <- 30
  use_partition <- TRUE
  close_loop <- FALSE

  if (!is.null(config)) {
    if (!is.null(config$trajectory$monocle3$n_pcs)) n_pcs <- config$trajectory$monocle3$n_pcs
    if (!is.null(config$trajectory$monocle3$use_partition)) use_partition <- config$trajectory$monocle3$use_partition
    if (!is.null(config$trajectory$monocle3$close_loop)) close_loop <- config$trajectory$monocle3$close_loop
  }

  message("  Preprocessing...")
  cds <- monocle3::preprocess_cds(cds, num_dim = n_pcs)

  message("  Reducing dimensions...")
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")

  message("  Clustering...")
  cds <- monocle3::cluster_cells(cds)

  message("  Learning trajectory graph...")
  cds <- monocle3::learn_graph(cds, use_partition = use_partition, close_loop = close_loop)

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

  graph <- monocle3::principal_graph(cds)$UMAP
  degree <- igraph::degree(graph)
  n_branch_points <- sum(degree > 2)

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
#' @keywords internal
run_slingshot_trajectory <- function(seurat_obj, config = NULL,
                                     subtype_column = "module_score_subtype",
                                     root_subtype = "Basal-like",
                                     reduction = "umap") {

  if (!requireNamespace("slingshot", quietly = TRUE)) {
    stop("slingshot package required")
  }

  library(slingshot)

  rd <- Seurat::Embeddings(seurat_obj, reduction = reduction)
  cl <- seurat_obj@meta.data[[subtype_column]]

  message("  Running slingshot...")
  sds <- slingshot::slingshot(rd, clusterLabels = cl, start.clus = root_subtype)

  pt <- slingshot::slingPseudotime(sds)

  if (ncol(pt) > 1) {
    message(sprintf("  Multiple lineages detected (%d), using average pseudotime", ncol(pt)))
    pt_avg <- rowMeans(pt, na.rm = TRUE)
  } else {
    pt_avg <- pt[, 1]
  }

  names(pt_avg) <- colnames(seurat_obj)

  return(list(
    sds = sds,
    pseudotime = pt_avg,
    pseudotime_matrix = pt,
    n_lineages = ncol(pt),
    method = "slingshot"
  ))
}

#' Calculate Consensus Pseudotime
#' @export
calculate_consensus_pseudotime <- function(pseudotime_matrix) {

  n_methods <- ncol(pseudotime_matrix)
  n_cells <- nrow(pseudotime_matrix)

  if (n_methods == 1) {
    pt <- pseudotime_matrix[, 1]
    pt_norm <- (pt - min(pt, na.rm = TRUE)) / (max(pt, na.rm = TRUE) - min(pt, na.rm = TRUE))
    return(list(
      pseudotime = pt_norm,
      confidence = rep(1.0, n_cells),
      n_methods = 1
    ))
  }

  pt_normalized <- apply(pseudotime_matrix, 2, function(x) {
    if (all(is.na(x))) return(x)
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  })

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

  consensus_pt <- rowMeans(pt_normalized, na.rm = TRUE)
  pt_sd <- apply(pt_normalized, 1, sd, na.rm = TRUE)
  pt_sd[is.na(pt_sd)] <- 0

  max_possible_sd <- 0.5
  confidence <- 1 - pmin(pt_sd / max_possible_sd, 1)
  n_valid <- rowSums(!is.na(pseudotime_matrix))

  return(list(
    pseudotime = consensus_pt,
    confidence = confidence,
    n_methods = n_methods,
    n_valid_per_cell = n_valid,
    normalized_matrix = pt_normalized
  ))
}

#' Evaluate Method Quality
#' @export
evaluate_method_quality <- function(pseudotime, subtypes, method, cds = NULL, sds = NULL) {

  quality <- list(method = method)

  pt_range <- diff(range(pseudotime, na.rm = TRUE))
  pt_iqr <- IQR(pseudotime, na.rm = TRUE)
  quality$spread <- list(
    range = pt_range,
    iqr = pt_iqr,
    cv = sd(pseudotime, na.rm = TRUE) / mean(pseudotime, na.rm = TRUE)
  )

  if (!is.null(subtypes) && length(unique(subtypes)) > 1) {
    pt_by_subtype <- tapply(pseudotime, subtypes, mean, na.rm = TRUE)
    total_var <- var(pseudotime, na.rm = TRUE)
    between_var <- var(pt_by_subtype)
    quality$separation <- between_var / total_var

    kw_test <- kruskal.test(pseudotime ~ subtypes)
    quality$kruskal_wallis <- list(
      statistic = kw_test$statistic,
      p_value = kw_test$p.value
    )
  } else {
    quality$separation <- NA
  }

  quality$continuity <- 1 - (sd(diff(sort(pseudotime)), na.rm = TRUE) / pt_range)

  # Calculate overall quality score
  quality$quality_components <- c(
    separation = ifelse(is.na(quality$separation), 0.5, quality$separation),
    continuity = quality$continuity,
    spread = min(1, quality$spread$cv)
  )

  quality$overall_quality <- mean(quality$quality_components, na.rm = TRUE)

  return(quality)
}

#' Validate Trajectory Against Expected Order
#' @export
validate_trajectory_order <- function(results, expected_order = c("Basal-like", "Transit-Amplifying",
                                                                  "Intermediate", "Specialized")) {

  message("\n--- Validating Trajectory Order ---")
  message(sprintf("  Expected: %s", paste(expected_order, collapse = " → ")))

  if (is.list(results) && "seurat_obj" %in% names(results)) {
    seurat_obj <- results$seurat_obj
    pt_by_subtype <- tapply(seurat_obj$pseudotime_consensus,
                            seurat_obj$module_score_subtype, mean, na.rm = TRUE)
  } else {
    stop("Expected results from run_multi_method_trajectory")
  }

  observed_order <- names(sort(pt_by_subtype))
  message(sprintf("  Observed: %s", paste(observed_order, collapse = " → ")))

  shared <- intersect(expected_order, observed_order)
  if (length(shared) < 3) {
    warning("Fewer than 3 shared subtypes for validation")
    return(data.frame(
      expected_order = paste(expected_order, collapse = " → "),
      observed_order = paste(observed_order, collapse = " → "),
      spearman_rho = NA,
      p_value = NA,
      n_subtypes = length(shared),
      concordant = NA
    ))
  }

  exp_ranks <- match(shared, expected_order)
  obs_ranks <- match(shared, observed_order)

  cor_test <- cor.test(exp_ranks, obs_ranks, method = "spearman")

  concordant <- cor_test$estimate > 0.5 && cor_test$p.value < 0.05

  message(sprintf("  Spearman rho: %.3f (p = %.4f)", cor_test$estimate, cor_test$p.value))
  message(sprintf("  Concordant: %s", ifelse(concordant, "YES", "NO")))

  return(data.frame(
    expected_order = paste(expected_order, collapse = " → "),
    observed_order = paste(observed_order, collapse = " → "),
    spearman_rho = as.numeric(cor_test$estimate),
    p_value = cor_test$p.value,
    n_subtypes = length(shared),
    concordant = concordant
  ))
}

#' Permutation Test for Trajectory Significance
#' @export
permutation_test_trajectory <- function(seurat_obj,
                                        subtype_column = "module_score_subtype",
                                        expected_order = c("Basal-like", "Transit-Amplifying",
                                                           "Intermediate", "Specialized"),
                                        n_permutations = 1000,
                                        parallel = TRUE) {

  message("\n--- Permutation Test for Trajectory Significance ---")

  subtypes <- seurat_obj@meta.data[[subtype_column]]
  pseudotime <- seurat_obj$pseudotime_consensus

  pt_by_subtype <- tapply(pseudotime, subtypes, mean, na.rm = TRUE)
  observed_order <- names(sort(pt_by_subtype))

  shared <- intersect(expected_order, observed_order)
  if (length(shared) < 3) {
    return(list(p_value = NA, message = "Insufficient shared subtypes"))
  }

  exp_ranks <- match(shared, expected_order)
  obs_ranks <- match(shared, observed_order)
  observed_cor <- cor(exp_ranks, obs_ranks, method = "spearman")

  message(sprintf("  Observed Spearman rho: %.3f", observed_cor))
  message(sprintf("  Running %d permutations...", n_permutations))

  run_permutation <- function(i) {
    perm_pt <- sample(pseudotime)
    perm_by_subtype <- tapply(perm_pt, subtypes, mean, na.rm = TRUE)
    perm_order <- names(sort(perm_by_subtype))
    perm_ranks <- match(shared, perm_order)
    cor(exp_ranks, perm_ranks, method = "spearman")
  }

  perm_cors <- sapply(1:n_permutations, run_permutation)

  p_value <- mean(perm_cors >= observed_cor)
  p_value_two_tailed <- mean(abs(perm_cors) >= abs(observed_cor))

  effect_size <- (observed_cor - mean(perm_cors)) / sd(perm_cors)

  if (p_value < 0.001) {
    interpretation <- "Highly significant trajectory order"
  } else if (p_value < 0.01) {
    interpretation <- "Significant trajectory order"
  } else if (p_value < 0.05) {
    interpretation <- "Marginally significant trajectory order"
  } else {
    interpretation <- "Non-significant trajectory order"
  }

  message(sprintf("  P-value (one-tailed): %.4f", p_value))
  message(sprintf("  Effect size (d): %.2f", effect_size))
  message(sprintf("  Interpretation: %s", interpretation))

  return(list(
    observed_correlation = observed_cor,
    permutation_correlations = perm_cors,
    p_value = p_value,
    p_value_two_tailed = p_value_two_tailed,
    effect_size = effect_size,
    n_permutations = n_permutations,
    significant = p_value < 0.05,
    interpretation = interpretation
  ))
}

#' Permutation Test for Method Concordance
#' @export
permutation_test_method_concordance <- function(pseudotime_matrix,
                                                n_permutations = 1000,
                                                parallel = TRUE) {

  message("\n--- Permutation Test for Method Concordance ---")

  n_methods <- ncol(pseudotime_matrix)

  if (n_methods < 2) {
    return(list(p_value = NA, message = "Need at least 2 methods"))
  }

  cor_matrix <- cor(pseudotime_matrix, method = "spearman", use = "pairwise.complete.obs")
  observed_mean_cor <- mean(cor_matrix[upper.tri(cor_matrix)])

  message(sprintf("  Observed mean correlation: %.3f", observed_mean_cor))

  run_permutation <- function(i) {
    perm_matrix <- apply(pseudotime_matrix, 2, sample)
    perm_cor <- cor(perm_matrix, method = "spearman", use = "pairwise.complete.obs")
    mean(perm_cor[upper.tri(perm_cor)])
  }

  perm_cors <- sapply(1:n_permutations, run_permutation)
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

#' Evaluate Trajectory Analysis Statistics
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

  if (is.null(expected_order) && !is.null(config$trajectory$expected_order)) {
    expected_order <- config$trajectory$expected_order
  }

  if (is.null(differentiation_markers)) {
    differentiation_markers <- list(
      early = c("SOX9", "KRT5", "KRT14", "TP63", "ITGA6", "ITGB1"),
      late = c("KRT20", "MUC2", "FABP1", "SI", "ALPI", "VIL1")
    )
  }

  # Mean pseudotime by subtype
  pt_by_subtype <- seurat_obj@meta.data %>%
    dplyr::group_by(module_score_subtype) %>%
    dplyr::summarise(
      mean_pseudotime = mean(pseudotime_consensus, na.rm = TRUE),
      sd_pseudotime = sd(pseudotime_consensus, na.rm = TRUE),
      n_cells = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(mean_pseudotime)

  stats$mean_pseudotime_by_subtype <- pt_by_subtype

  message("Mean pseudotime by subtype:")
  print(pt_by_subtype)

  # Method correlations
  if (length(results$methods_used) > 1) {
    stats$mean_correlation <- mean(results$method_correlations[upper.tri(results$method_correlations)])
    message(sprintf("\nMean method correlation: %.3f", stats$mean_correlation))
  } else {
    stats$mean_correlation <- NA
  }

  # Order validation
  observed_order <- pt_by_subtype$module_score_subtype

  if (!is.null(expected_order)) {
    shared <- intersect(expected_order, observed_order)
    if (length(shared) >= 3) {
      exp_ranks <- match(shared, expected_order)
      obs_ranks <- match(shared, observed_order)
      rank_cor <- cor(exp_ranks, obs_ranks, method = "spearman")

      stats$order_validation <- list(
        expected = expected_order,
        observed = observed_order,
        rank_correlation = rank_cor,
        order_preserved = rank_cor > 0.5
      )

      message(sprintf("\nOrder validation: ρ = %.3f (%s)",
                      rank_cor,
                      ifelse(rank_cor > 0.5, "CONCORDANT", "NOT CONCORDANT")))
    }
  }

  # Marker validation
  message("\n--- Marker-Pseudotime Correlations ---")

  available_genes <- rownames(seurat_obj)
  early_markers <- intersect(differentiation_markers$early, available_genes)
  late_markers <- intersect(differentiation_markers$late, available_genes)

  marker_cors <- c()
  marker_validation <- data.frame(
    marker = character(),
    type = character(),
    correlation = numeric(),
    expected_direction = character(),
    correct_direction = logical(),
    stringsAsFactors = FALSE
  )

  if (length(early_markers) > 0) {
    for (marker in early_markers) {
      expr <- Seurat::GetAssayData(seurat_obj, layer = "data")[marker, ]
      cor_val <- cor(expr, seurat_obj$pseudotime_consensus, method = "spearman", use = "complete.obs")
      marker_cors <- c(marker_cors, cor_val)
      marker_validation <- rbind(marker_validation, data.frame(
        marker = marker,
        type = "early",
        correlation = cor_val,
        expected_direction = "negative",
        correct_direction = cor_val < 0
      ))
    }
  }

  if (length(late_markers) > 0) {
    for (marker in late_markers) {
      expr <- Seurat::GetAssayData(seurat_obj, layer = "data")[marker, ]
      cor_val <- cor(expr, seurat_obj$pseudotime_consensus, method = "spearman", use = "complete.obs")
      marker_cors <- c(marker_cors, cor_val)
      marker_validation <- rbind(marker_validation, data.frame(
        marker = marker,
        type = "late",
        correlation = cor_val,
        expected_direction = "positive",
        correct_direction = cor_val > 0
      ))
    }
  }

  stats$marker_validation <- marker_validation

  n_correct <- sum(marker_validation$correct_direction, na.rm = TRUE)
  n_total <- nrow(marker_validation)

  stats$marker_validation_summary <- list(
    n_correct_direction = n_correct,
    n_total = n_total,
    percent_correct = n_correct / n_total * 100,
    validated = n_correct / n_total >= 0.5
  )

  message(sprintf("  %d/%d markers show expected direction (%.0f%%)",
                  n_correct, n_total, n_correct / n_total * 100))

  # Cell cycle correlation
  if ("Phase" %in% colnames(seurat_obj@meta.data)) {
    cycling <- seurat_obj$Phase %in% c("S", "G2M")
    cor_cc <- cor(as.numeric(cycling), seurat_obj$pseudotime_consensus,
                  method = "spearman", use = "complete.obs")

    stats$cell_cycle_correlation <- list(
      rho = cor_cc,
      interpretation = ifelse(abs(cor_cc) > 0.3, "Potentially confounded", "Low confounding")
    )

    message(sprintf("\nCell cycle correlation: ρ = %.3f (%s)",
                    cor_cc, stats$cell_cycle_correlation$interpretation))
  }

  # Quality score
  quality_components <- c(
    method_agreement = ifelse(is.na(stats$mean_correlation), 0.5, stats$mean_correlation),
    order_preservation = ifelse(is.null(stats$order_validation), 0.5,
                                ifelse(stats$order_validation$rank_correlation > 0,
                                       stats$order_validation$rank_correlation, 0)),
    marker_validation = stats$marker_validation_summary$percent_correct / 100
  )

  overall_quality <- mean(quality_components, na.rm = TRUE)

  grade <- if (overall_quality >= 0.8) "Excellent"
  else if (overall_quality >= 0.6) "Good"
  else if (overall_quality >= 0.4) "Moderate"
  else "Review needed"

  stats$quality_score <- list(
    components = quality_components,
    overall = overall_quality,
    grade = grade
  )

  message(sprintf("\n=== Overall Quality Score: %.2f (%s) ===", overall_quality, grade))

  return(stats)
}

#' Generate Trajectory Report
#' @export
generate_trajectory_report <- function(stats_results, output_path = NULL) {

  report <- c(
    "# Trajectory Analysis Evaluation Report",
    "",
    sprintf("Generated: %s", Sys.time()),
    ""
  )

  # Quality score
  if (!is.null(stats_results$quality_score)) {
    qs <- stats_results$quality_score
    report <- c(report, "## Overall Quality", "",
                sprintf("**Score: %.2f (%s)**", qs$overall, qs$grade), "",
                "Components:", "")
    for (name in names(qs$components)) {
      report <- c(report, sprintf("- %s: %.3f", name, qs$components[name]))
    }
    report <- c(report, "")
  }

  # Marker validation
  if (!is.null(stats_results$marker_validation)) {
    report <- c(report, "## Marker-Pseudotime Correlations", "")
    report <- c(report, "| Marker | Type | Correlation | Expected | Correct |")
    report <- c(report, "|--------|------|-------------|----------|---------|")

    mv <- stats_results$marker_validation
    for (i in 1:nrow(mv)) {
      report <- c(report, sprintf("| %s | %s | %.3f | %s | %s |",
                                  mv$marker[i], mv$type[i], mv$correlation[i],
                                  mv$expected_direction[i],
                                  ifelse(mv$correct_direction[i], "✓", "✗")))
    }
    report <- c(report, "")
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

  report_text <- paste(report, collapse = "\n")

  if (!is.null(output_path)) {
    writeLines(report_text, output_path)
    message(sprintf("Report saved to: %s", output_path))
  }

  return(report_text)
}

#' Plot Cell Cycle Trajectory
#' @export
plot_cell_cycle_trajectory <- function(seurat_obj, config = NULL) {

  require(ggplot2)
  require(patchwork)

  cc_col <- ifelse("cc_consensus" %in% colnames(seurat_obj@meta.data), "cc_consensus", "Phase")

  if (!cc_col %in% colnames(seurat_obj@meta.data)) {
    return(NULL)
  }

  cc_colors <- c(G1 = "#F8766D", S = "#00BA38", G2M = "#619CFF")
  if (!is.null(config$visualization$colors$cell_cycle)) {
    cc_colors <- unlist(config$visualization$colors$cell_cycle)
  }

  pt_data <- data.frame(
    pseudotime = seurat_obj$pseudotime_consensus,
    phase = seurat_obj@meta.data[[cc_col]],
    subtype = seurat_obj$module_score_subtype
  )

  p1 <- ggplot(pt_data, aes(x = phase, y = pseudotime, fill = phase)) +
    geom_violin() +
    geom_boxplot(width = 0.2, fill = "white") +
    scale_fill_manual(values = cc_colors) +
    theme_minimal() +
    labs(title = "Pseudotime by Cell Cycle Phase", x = "Phase", y = "Pseudotime")

  p2 <- ggplot(pt_data, aes(x = pseudotime, fill = phase)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = cc_colors) +
    theme_minimal() +
    labs(title = "Cell Cycle Phase Distribution Along Trajectory",
         x = "Pseudotime", y = "Density")

  return(p1 / p2)
}
