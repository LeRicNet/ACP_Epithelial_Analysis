# ==============================================================================
# R/trajectory/trajectory_analysis.R
# ==============================================================================
# Monocle3-based trajectory analysis with per-sample validation
# ==============================================================================

#' Run Monocle3 Trajectory Analysis
#'
#' Performs trajectory inference using Monocle3 with configurable root selection
#'
#' @param seurat_obj Seurat object with classifications
#' @param config Configuration list
#' @param subtype_column Column name containing cell type classifications
#' @return cell_data_set object with trajectory
#' @export
run_trajectory_analysis <- function(seurat_obj, config, 
                                    subtype_column = "module_score_subtype") {
  
  message("\n========================================")
  message("Running Trajectory Analysis")
  message("========================================\n")
  
  # Load Monocle3
  if (!requireNamespace("monocle3", quietly = TRUE)) {
    stop("Package 'monocle3' required. Install from: https://cole-trapnell-lab.github.io/monocle3/")
  }
  
  library(monocle3)
  library(SeuratWrappers)
  
  # Convert Seurat to cell_data_set
  message("Converting Seurat object to cell_data_set...")
  cds <- SeuratWrappers::as.cell_data_set(seurat_obj)
  
  # Transfer metadata
  SummarizedExperiment::colData(cds)$subtype <- seurat_obj@meta.data[[subtype_column]]
  
  # Preprocess
  message("Preprocessing...")
  n_pcs <- config$trajectory$monocle3$n_pcs
  if (is.null(n_pcs)) n_pcs <- 30
  
  cds <- monocle3::preprocess_cds(cds, num_dim = n_pcs)
  
  # Reduce dimensions
  message("Reducing dimensions...")
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
  
  # Cluster cells
  message("Clustering cells...")
  cds <- monocle3::cluster_cells(cds)
  
  # Learn graph
  message("Learning trajectory graph...")
  use_partition <- config$trajectory$monocle3$use_partition
  if (is.null(use_partition)) use_partition <- TRUE
  
  close_loop <- config$trajectory$monocle3$close_loop
  if (is.null(close_loop)) close_loop <- FALSE
  
  cds <- monocle3::learn_graph(cds, 
                                use_partition = use_partition,
                                close_loop = close_loop)
  
  # Order cells (select root)
  message("Ordering cells along trajectory...")
  root_subtype <- config$trajectory$monocle3$root_subtype
  if (is.null(root_subtype)) root_subtype <- "Basal-like"
  
  root_method <- config$trajectory$monocle3$root_selection_method
  if (is.null(root_method)) root_method <- "all_cells"
  
  cds <- order_cells_by_subtype(cds, root_subtype, root_method, config)
  
  message("\nTrajectory analysis complete!")
  
  return(cds)
}

#' Order Cells by Subtype
#'
#' Orders cells along trajectory using specified root subtype
#'
#' @param cds cell_data_set object
#' @param root_subtype Subtype to use as root
#' @param method Root selection method
#' @param config Configuration list
#' @return cell_data_set with pseudotime
#' @export
order_cells_by_subtype <- function(cds, root_subtype, method = "all_cells", config = NULL) {
  
  library(monocle3)
  
  # Get cells of root subtype
  subtypes <- SummarizedExperiment::colData(cds)$subtype
  root_cells <- colnames(cds)[subtypes == root_subtype]
  
  if (length(root_cells) == 0) {
    warning(paste("No cells found for root subtype:", root_subtype))
    warning("Using first partition node as root")
    cds <- monocle3::order_cells(cds)
    return(cds)
  }
  
  message(sprintf("  Found %d cells of root subtype: %s", length(root_cells), root_subtype))
  
  if (method == "all_cells") {
    # Use all cells of the root subtype
    cds <- monocle3::order_cells(cds, root_cells = root_cells)
    
  } else if (method == "random_sample") {
    # Random sample of root cells
    n_sample <- config$trajectory$monocle3$n_root_cells
    if (is.null(n_sample)) n_sample <- 50
    n_sample <- min(n_sample, length(root_cells))
    
    set.seed(config$reproducibility$seed)
    sampled_roots <- sample(root_cells, n_sample)
    cds <- monocle3::order_cells(cds, root_cells = sampled_roots)
    
  } else if (method == "graph_node") {
    # Find graph node closest to root subtype centroid
    # Get UMAP coordinates
    umap_coords <- SingleCellExperiment::reducedDims(cds)$UMAP
    
    # Calculate centroid of root cells
    root_idx <- which(colnames(cds) %in% root_cells)
    centroid <- colMeans(umap_coords[root_idx, , drop = FALSE])
    
    # Get principal graph nodes
    graph <- monocle3::principal_graph(cds)$UMAP
    node_coords <- t(igraph::V(graph)$Y)  # Node coordinates
    
    # Find nearest node to centroid
    distances <- sqrt(rowSums((t(node_coords) - centroid)^2))
    nearest_node <- names(distances)[which.min(distances)]
    
    # Order from this node
    cds <- monocle3::order_cells(cds, root_pr_nodes = nearest_node)
  }
  
  return(cds)
}

#' Run Per-Sample Trajectory Analysis
#'
#' Runs trajectory analysis separately for each sample
#'
#' @param seurat_obj Seurat object
#' @param config Configuration list
#' @param sample_column Column containing sample IDs
#' @param subtype_column Column containing subtypes
#' @return List with per-sample results
#' @export
run_per_sample_trajectory <- function(seurat_obj, config,
                                      sample_column = "sample",
                                      subtype_column = "module_score_subtype") {
  
  message("\n========================================")
  message("Running Per-Sample Trajectory Analysis")
  message("========================================\n")
  
  # Get unique samples
  samples <- unique(seurat_obj@meta.data[[sample_column]])
  samples <- samples[!is.na(samples)]
  
  message(sprintf("Found %d samples: %s", length(samples), paste(samples, collapse = ", ")))
  
  # Minimum cells per sample

  min_cells <- config$trajectory$per_sample$min_cells_per_sample
  if (is.null(min_cells)) min_cells <- 50
  
  results <- list()
  
  for (sample_id in samples) {
    message(sprintf("\n--- Processing sample: %s ---", sample_id))
    
    # Subset to this sample
    cells_in_sample <- seurat_obj@meta.data[[sample_column]] == sample_id
    cells_in_sample[is.na(cells_in_sample)] <- FALSE
    
    n_cells <- sum(cells_in_sample)
    message(sprintf("  Cells in sample: %d", n_cells))
    
    if (n_cells < min_cells) {
      message(sprintf("  Skipping: fewer than %d cells", min_cells))
      results[[sample_id]] <- list(
        sample = sample_id,
        n_cells = n_cells,
        skipped = TRUE,
        reason = "insufficient_cells"
      )
      next
    }
    
    # Subset Seurat object
    seurat_sample <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_in_sample])
    
    # Check subtype distribution
    subtype_counts <- table(seurat_sample@meta.data[[subtype_column]])
    message("  Subtype distribution:")
    print(subtype_counts)
    
    # Run trajectory
    tryCatch({
      cds <- run_trajectory_analysis(seurat_sample, config, subtype_column)
      
      # Extract results
      pseudotime <- monocle3::pseudotime(cds)
      subtypes <- SummarizedExperiment::colData(cds)$subtype
      
      # Calculate mean pseudotime by subtype
      pseudotime_by_subtype <- tapply(pseudotime, subtypes, mean, na.rm = TRUE)
      
      # Determine subtype ordering
      subtype_order <- names(sort(pseudotime_by_subtype))
      
      # Count branch points
      graph <- monocle3::principal_graph(cds)$UMAP
      degree <- igraph::degree(graph)
      n_branch_points <- sum(degree > 2)
      
      results[[sample_id]] <- list(
        sample = sample_id,
        n_cells = n_cells,
        skipped = FALSE,
        cds = cds,
        pseudotime = pseudotime,
        pseudotime_by_subtype = pseudotime_by_subtype,
        subtype_order = subtype_order,
        subtype_rank = rank(pseudotime_by_subtype),
        n_branch_points = n_branch_points,
        subtype_counts = subtype_counts
      )
      
      message(sprintf("  Branch points detected: %d", n_branch_points))
      message(sprintf("  Subtype ordering: %s", paste(subtype_order, collapse = " → ")))
      
    }, error = function(e) {
      message(sprintf("  ERROR: %s", e$message))
      results[[sample_id]] <- list(
        sample = sample_id,
        n_cells = n_cells,
        skipped = TRUE,
        reason = "trajectory_failed",
        error = e$message
      )
    })
  }
  
  # Calculate cross-sample concordance
  message("\n--- Cross-Sample Concordance ---")
  results <- calculate_concordance(results, config)
  
  return(results)
}

#' Calculate Cross-Sample Concordance
#'
#' Calculates Spearman correlation of subtype rankings across samples
#'
#' @param results List of per-sample results
#' @param config Configuration list
#' @return Updated results list with concordance metrics
#' @export
calculate_concordance <- function(results, config) {
  
  # Get valid results
  valid_results <- results[sapply(results, function(x) !x$skipped)]
  
  if (length(valid_results) < 2) {
    message("  Insufficient valid samples for concordance analysis")
    return(results)
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
  
  # Calculate pairwise Spearman correlations
  n_samples <- ncol(rank_matrix)
  correlations <- c()
  
  for (i in 1:(n_samples - 1)) {
    for (j in (i + 1):n_samples) {
      # Get shared subtypes
      shared <- !is.na(rank_matrix[, i]) & !is.na(rank_matrix[, j])
      if (sum(shared) >= 3) {
        rho <- cor(rank_matrix[shared, i], rank_matrix[shared, j], 
                   method = "spearman", use = "complete.obs")
        correlations <- c(correlations, rho)
      }
    }
  }
  
  # Summary statistics
  concordance_summary <- list(
    rank_matrix = rank_matrix,
    mean_correlation = mean(correlations, na.rm = TRUE),
    sd_correlation = sd(correlations, na.rm = TRUE),
    min_correlation = min(correlations, na.rm = TRUE),
    max_correlation = max(correlations, na.rm = TRUE),
    n_comparisons = length(correlations),
    pairwise_correlations = correlations
  )
  
  message(sprintf("  Mean Spearman rho: %.3f ± %.3f", 
                  concordance_summary$mean_correlation,
                  concordance_summary$sd_correlation))
  
  attr(results, "concordance") <- concordance_summary
  
  return(results)
}

#' Extract Pseudotime Summary Statistics
#'
#' Creates summary data frame from trajectory results
#'
#' @param cds cell_data_set object OR list of per-sample results
#' @param subtype_column Column name for subtypes
#' @return Data frame with pseudotime statistics
#' @export
summarize_pseudotime <- function(cds, subtype_column = "subtype") {
  
  # Handle list of per-sample results
  if (is.list(cds) && !inherits(cds, "cell_data_set")) {
    # This is per-sample results
    results_list <- cds
    
    summary_df <- do.call(rbind, lapply(names(results_list), function(sample_id) {
      res <- results_list[[sample_id]]
      if (res$skipped) return(NULL)
      
      pseudotime <- res$pseudotime
      subtypes <- SummarizedExperiment::colData(res$cds)$subtype
      
      df <- data.frame(
        sample = sample_id,
        subtype = subtypes,
        pseudotime = pseudotime,
        stringsAsFactors = FALSE
      )
      
      return(df)
    }))
    
  } else {
    # Single cds object
    pseudotime <- monocle3::pseudotime(cds)
    subtypes <- SummarizedExperiment::colData(cds)[[subtype_column]]
    
    summary_df <- data.frame(
      subtype = subtypes,
      pseudotime = pseudotime,
      stringsAsFactors = FALSE
    )
  }
  
  return(summary_df)
}

#' Validate Trajectory with Expected Ordering
#'
#' Tests whether observed ordering matches expected differentiation sequence
#'
#' @param results Trajectory results (cds or list)
#' @param expected_order Expected subtype ordering (early to late)
#' @param config Configuration list
#' @return Data frame with validation metrics
#' @export
validate_trajectory_ordering <- function(results, 
                                         expected_order = c("Basal-like", "Transit-Amplifying", 
                                                            "Intermediate", "Specialized"),
                                         config = NULL) {
  
  message("\n--- Validating Trajectory Ordering ---")
  message(sprintf("Expected order: %s", paste(expected_order, collapse = " → ")))
  
  # Handle per-sample results
  if (is.list(results) && !inherits(results, "cell_data_set")) {
    
    validation_df <- do.call(rbind, lapply(names(results), function(sample_id) {
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
      
      observed_order <- res$subtype_order
      observed_ranks <- res$subtype_rank
      
      # Get shared subtypes
      shared <- intersect(names(observed_ranks), expected_order)
      
      if (length(shared) < 3) {
        return(data.frame(
          sample = sample_id,
          spearman_rho = NA,
          p_value = NA,
          n_subtypes = length(shared),
          observed_order = paste(observed_order, collapse = " → "),
          concordant = NA
        ))
      }
      
      # Calculate Spearman correlation
      expected_ranks <- seq_along(expected_order)
      names(expected_ranks) <- expected_order
      
      obs <- observed_ranks[shared]
      exp <- expected_ranks[shared]
      
      cor_test <- cor.test(obs, exp, method = "spearman")
      
      data.frame(
        sample = sample_id,
        spearman_rho = cor_test$estimate,
        p_value = cor_test$p.value,
        n_subtypes = length(shared),
        observed_order = paste(observed_order, collapse = " → "),
        concordant = cor_test$estimate > 0.5 && cor_test$p.value < 0.05
      )
    }))
    
  } else {
    # Single cds
    pseudotime <- monocle3::pseudotime(results)
    subtypes <- SummarizedExperiment::colData(results)$subtype
    
    mean_pt <- tapply(pseudotime, subtypes, mean, na.rm = TRUE)
    observed_order <- names(sort(mean_pt))
    observed_ranks <- rank(mean_pt)
    
    shared <- intersect(names(observed_ranks), expected_order)
    
    expected_ranks <- seq_along(expected_order)
    names(expected_ranks) <- expected_order
    
    obs <- observed_ranks[shared]
    exp <- expected_ranks[shared]
    
    cor_test <- cor.test(obs, exp, method = "spearman")
    
    validation_df <- data.frame(
      sample = "pooled",
      spearman_rho = cor_test$estimate,
      p_value = cor_test$p.value,
      n_subtypes = length(shared),
      observed_order = paste(observed_order, collapse = " → "),
      concordant = cor_test$estimate > 0.5 && cor_test$p.value < 0.05
    )
  }
  
  # Print summary
  message(sprintf("\nValidation Summary:"))
  message(sprintf("  Samples concordant with expected: %d/%d", 
                  sum(validation_df$concordant, na.rm = TRUE),
                  sum(!is.na(validation_df$concordant))))
  message(sprintf("  Mean Spearman rho: %.3f", 
                  mean(validation_df$spearman_rho, na.rm = TRUE)))
  
  return(validation_df)
}
