# ==============================================================================
# R/classification/module_score_classification.R
# ==============================================================================
# Module score-based epithelial subtype classification
# Replaces binary threshold approach with continuous scoring
# ==============================================================================

#' Calculate Module Scores for All Epithelial Signatures
#'
#' Calculates UCell or Seurat module scores for all epithelial subtype signatures
#'
#' @param seurat_obj Seurat object
#' @param config Configuration list
#' @param method Scoring method: "UCell" (default) or "Seurat"
#' @return Seurat object with module scores added to metadata
#' @export
calculate_epithelial_scores <- function(seurat_obj, config, method = "UCell") {
  
  message("Calculating epithelial subtype module scores...")
  
  # Get available genes
  available_genes <- rownames(seurat_obj)
  
  # Get epithelial subtype signatures
  sigs <- list()
  for (name in names(config$signatures$epithelial_subtypes)) {
    genes <- config$signatures$epithelial_subtypes[[name]]
    present <- intersect(genes, available_genes)
    
    if (length(present) >= 3) {
      sigs[[name]] <- present
      message(sprintf("  %s: %d/%d genes present", name, length(present), length(genes)))
    } else {
      warning(sprintf("  %s: Only %d genes present, skipping", name, length(present)))
    }
  }
  
  if (length(sigs) == 0) {
    stop("No valid signatures found!")
  }
  
  # Calculate scores
  if (method == "UCell") {
    if (!requireNamespace("UCell", quietly = TRUE)) {
      message("UCell not available, falling back to Seurat AddModuleScore")
      method <- "Seurat"
    } else {
      seurat_obj <- UCell::AddModuleScore_UCell(
        seurat_obj, 
        features = sigs, 
        name = ""
      )
      # Rename columns for consistency
      for (name in names(sigs)) {
        old_name <- paste0(name, "_UCell")
        new_name <- paste0("score_", name)
        if (old_name %in% colnames(seurat_obj@meta.data)) {
          seurat_obj@meta.data[[new_name]] <- seurat_obj@meta.data[[old_name]]
        }
      }
    }
  }
  
  if (method == "Seurat") {
    for (i in seq_along(sigs)) {
      name <- names(sigs)[i]
      seurat_obj <- Seurat::AddModuleScore(
        seurat_obj,
        features = list(sigs[[name]]),
        name = paste0("score_", name),
        seed = config$reproducibility$seed
      )
      # Seurat appends "1" to the name
      old_name <- paste0("score_", name, "1")
      new_name <- paste0("score_", name)
      if (old_name %in% colnames(seurat_obj@meta.data)) {
        seurat_obj@meta.data[[new_name]] <- seurat_obj@meta.data[[old_name]]
        seurat_obj@meta.data[[old_name]] <- NULL
      }
    }
  }
  
  # Also calculate proliferation score if not already done
  prolif_sig <- config$signatures$epithelial_subtypes$Transit_Amplifying
  if (!is.null(prolif_sig)) {
    prolif_present <- intersect(prolif_sig, available_genes)
    if (length(prolif_present) >= 3) {
      if (method == "UCell" && requireNamespace("UCell", quietly = TRUE)) {
        seurat_obj <- UCell::AddModuleScore_UCell(
          seurat_obj,
          features = list(Proliferation = prolif_present),
          name = ""
        )
        seurat_obj@meta.data$score_Proliferation <- seurat_obj@meta.data$Proliferation_UCell
      } else {
        seurat_obj <- Seurat::AddModuleScore(
          seurat_obj,
          features = list(prolif_present),
          name = "score_Proliferation",
          seed = config$reproducibility$seed
        )
        seurat_obj@meta.data$score_Proliferation <- seurat_obj@meta.data$score_Proliferation1
        seurat_obj@meta.data$score_Proliferation1 <- NULL
      }
    }
  }
  
  message("Module scores calculated successfully")
  return(seurat_obj)
}

#' Classify Cells by Highest Module Score
#'
#' Assigns epithelial subtypes based on highest module score above threshold
#'
#' @param seurat_obj Seurat object with module scores calculated
#' @param config Configuration list
#' @param include_proliferation Whether to use proliferation score for TA identification
#' @return Seurat object with classification added to metadata
#' @export
classify_by_module_score <- function(seurat_obj, config, include_proliferation = TRUE) {
  
  message("Classifying cells by module score...")
  
  # Get score columns
  score_cols <- grep("^score_", colnames(seurat_obj@meta.data), value = TRUE)
  
  # Identify core subtype scores (exclude Proliferation for initial classification)
  core_scores <- score_cols[!grepl("Proliferation|Transit", score_cols)]
  
  if (length(core_scores) == 0) {
    stop("No module score columns found. Run calculate_epithelial_scores first.")
  }
  
  # Get parameters from config
  min_threshold <- config$classification$min_score_threshold
  if (is.null(min_threshold)) min_threshold <- 0.1
  
  # Extract score matrix
  score_matrix <- seurat_obj@meta.data[, core_scores, drop = FALSE]
  colnames(score_matrix) <- gsub("^score_", "", colnames(score_matrix))
  
  # Find max score and corresponding subtype for each cell
  max_scores <- apply(score_matrix, 1, max, na.rm = TRUE)
  max_subtype <- apply(score_matrix, 1, function(x) {
    if (all(is.na(x))) return(NA)
    colnames(score_matrix)[which.max(x)]
  })
  
  # Assign classifications
  classification <- ifelse(max_scores >= min_threshold, max_subtype, "Unassigned")
  
  # Handle proliferating/Transit-Amplifying cells
  if (include_proliferation && "score_Proliferation" %in% colnames(seurat_obj@meta.data)) {
    
    prolif_threshold <- config$classification$transit_amplifying$proliferation_threshold
    if (is.null(prolif_threshold)) prolif_threshold <- 0.15
    
    override_threshold <- config$classification$transit_amplifying$override_threshold
    if (is.null(override_threshold)) override_threshold <- 0.3
    
    require_cycling <- config$classification$transit_amplifying$require_cycling
    if (is.null(require_cycling)) require_cycling <- TRUE
    
    prolif_scores <- seurat_obj@meta.data$score_Proliferation
    
    # Identify cycling cells if cell cycle info available
    is_cycling <- rep(TRUE, ncol(seurat_obj))
    if (require_cycling && "Phase" %in% colnames(seurat_obj@meta.data)) {
      is_cycling <- seurat_obj@meta.data$Phase %in% c("S", "G2M")
    }
    
    # Override classification for highly proliferating cells
    override_as_ta <- config$classification$transit_amplifying$override_on_high_proliferation
    if (is.null(override_as_ta)) override_as_ta <- TRUE
    
    if (override_as_ta) {
      # High proliferation score overrides other classification
      high_prolif <- prolif_scores >= override_threshold & is_cycling
      classification[high_prolif] <- "Transit-Amplifying"
    }
    
    # Also label Unassigned cycling cells as TA if above threshold
    unassigned_cycling <- classification == "Unassigned" & 
                          prolif_scores >= prolif_threshold & 
                          is_cycling
    classification[unassigned_cycling] <- "Transit-Amplifying"
  }
  
  # Handle hybrid states (cells with high scores in multiple categories)
  basal_scores <- seurat_obj@meta.data$score_Basal_like
  int_scores <- seurat_obj@meta.data$score_Intermediate
  spec_scores <- seurat_obj@meta.data$score_Specialized
  
  if (!is.null(basal_scores) && !is.null(int_scores)) {
    # Basal-Intermediate: high in both Basal and Intermediate
    hybrid_threshold <- min_threshold * 0.8  # Slightly lower threshold for hybrids
    basal_int <- basal_scores >= hybrid_threshold & 
                 int_scores >= hybrid_threshold &
                 classification %in% c("Basal_like", "Intermediate")
    # Only assign if both scores are relatively close
    score_ratio <- abs(basal_scores - int_scores) / pmax(basal_scores, int_scores)
    close_scores <- score_ratio < 0.5
    classification[basal_int & close_scores] <- "Basal-Intermediate"
  }
  
  if (!is.null(int_scores) && !is.null(spec_scores)) {
    # Intermediate-Specialized: high in both
    int_spec <- int_scores >= hybrid_threshold & 
                spec_scores >= hybrid_threshold &
                classification %in% c("Intermediate", "Specialized")
    score_ratio <- abs(int_scores - spec_scores) / pmax(int_scores, spec_scores)
    close_scores <- score_ratio < 0.5
    classification[int_spec & close_scores] <- "Intermediate-Specialized"
  }
  
  # Standardize subtype names
  classification <- gsub("_", "-", classification)
  classification <- gsub("Basal-like", "Basal-like", classification)
  
  # Add to metadata
  seurat_obj$module_score_subtype <- classification
  
  # Report results
  message("\nClassification results:")
  print(table(classification))
  
  return(seurat_obj)
}

#' Compare Classification Methods
#'
#' Compares module score classification with binary threshold classification
#'
#' @param seurat_obj Seurat object with both classifications
#' @param config Configuration list
#' @return Data frame with comparison metrics
#' @export
compare_classifications <- function(seurat_obj, config) {
  
  # Check if both classifications exist
  has_module <- "module_score_subtype" %in% colnames(seurat_obj@meta.data)
  has_binary <- "epithelial_subtype" %in% colnames(seurat_obj@meta.data)
  
  if (!has_module || !has_binary) {
    stop("Both module_score_subtype and epithelial_subtype required for comparison")
  }
  
  module_class <- seurat_obj$module_score_subtype
  binary_class <- seurat_obj$epithelial_subtype
  
  # Create confusion matrix
  confusion <- table(Binary = binary_class, ModuleScore = module_class)
  
  # Calculate agreement
  agreement <- sum(diag(as.matrix(confusion))) / sum(confusion)
  
  # Calculate Cohen's Kappa
  # (agreement beyond chance)
  n <- sum(confusion)
  p_observed <- sum(diag(as.matrix(confusion))) / n
  
  row_totals <- rowSums(confusion)
  col_totals <- colSums(confusion)
  p_expected <- sum(row_totals * col_totals) / (n^2)
  
  kappa <- (p_observed - p_expected) / (1 - p_expected)
  
  # Per-class metrics
  classes <- unique(c(binary_class, module_class))
  per_class_metrics <- data.frame(
    class = classes,
    n_binary = sapply(classes, function(x) sum(binary_class == x)),
    n_module = sapply(classes, function(x) sum(module_class == x)),
    pct_change = NA
  )
  per_class_metrics$pct_change <- 
    (per_class_metrics$n_module - per_class_metrics$n_binary) / 
    pmax(per_class_metrics$n_binary, 1) * 100
  
  # Summary
  results <- list(
    confusion_matrix = confusion,
    overall_agreement = agreement,
    cohens_kappa = kappa,
    per_class_metrics = per_class_metrics,
    binary_counts = table(binary_class),
    module_counts = table(module_class)
  )
  
  message(sprintf("Overall agreement: %.1f%%", agreement * 100))
  message(sprintf("Cohen's Kappa: %.3f", kappa))
  
  return(results)
}

#' Visualize Module Scores
#'
#' Creates diagnostic plots of module scores
#'
#' @param seurat_obj Seurat object with module scores
#' @param config Configuration list
#' @return List of ggplot objects
#' @export
visualize_module_scores <- function(seurat_obj, config) {
  
  require(ggplot2)
  require(tidyr)
  require(dplyr)
  
  # Get score columns
  score_cols <- grep("^score_", colnames(seurat_obj@meta.data), value = TRUE)
  
  # Prepare data
  score_data <- seurat_obj@meta.data[, c(score_cols, "module_score_subtype")]
  colnames(score_data) <- gsub("^score_", "", colnames(score_data))
  
  score_long <- score_data %>%
    tidyr::pivot_longer(
      cols = -module_score_subtype,
      names_to = "Signature",
      values_to = "Score"
    )
  
  # Get colors
  colors <- get_colors(config, "epithelial_subtypes")
  
  plots <- list()
  
  # 1. Score distributions by signature
  plots$score_distributions <- ggplot(score_long, aes(x = Score, fill = Signature)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Signature, scales = "free_y") +
    geom_vline(xintercept = config$classification$min_score_threshold, 
               linetype = "dashed", color = "red") +
    theme_bw() +
    labs(title = "Module Score Distributions",
         subtitle = paste("Red line = classification threshold:", 
                          config$classification$min_score_threshold))
  
  # 2. Scores by assigned subtype
  plots$scores_by_subtype <- ggplot(score_long, 
                                     aes(x = module_score_subtype, y = Score, 
                                         fill = module_score_subtype)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~Signature, scales = "free_y") +
    scale_fill_manual(values = colors, na.value = "gray50") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Module Scores by Assigned Subtype",
         x = "Assigned Subtype", y = "Module Score")
  
  # 3. Score heatmap
  score_means <- score_data %>%
    group_by(module_score_subtype) %>%
    summarise(across(-module_score_subtype, mean, na.rm = TRUE)) %>%
    tidyr::pivot_longer(-module_score_subtype, names_to = "Signature", values_to = "Score")
  
  plots$score_heatmap <- ggplot(score_means, 
                                 aes(x = Signature, y = module_score_subtype, fill = Score)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", Score)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Mean Module Scores by Subtype",
         x = "Signature", y = "Assigned Subtype")
  
  # 4. UMAP with scores (if available)
  if ("umap" %in% names(seurat_obj@reductions) || "UMAP" %in% names(seurat_obj@reductions)) {
    for (sig in score_cols) {
      sig_name <- gsub("^score_", "", sig)
      p <- Seurat::FeaturePlot(seurat_obj, features = sig) +
        scale_color_viridis_c() +
        labs(title = paste(sig_name, "Score"))
      plots[[paste0("umap_", sig_name)]] <- p
    }
  }
  
  return(plots)
}

#' Run Complete Classification Pipeline
#'
#' Runs the full classification pipeline from raw data to classified object
#'
#' @param seurat_obj Seurat object
#' @param config Configuration list
#' @param compare_with_binary Also run binary classification for comparison
#' @return Seurat object with classifications
#' @export
run_classification_pipeline <- function(seurat_obj, config, compare_with_binary = TRUE) {
  
  message("\n========================================")
  message("Running Classification Pipeline")
  message("========================================\n")
  
  # Step 1: Calculate module scores
  seurat_obj <- calculate_epithelial_scores(seurat_obj, config)
  
  # Step 2: Classify by module score
  seurat_obj <- classify_by_module_score(seurat_obj, config)
  
  # Step 3: Compare with binary classification if exists
  if (compare_with_binary && "epithelial_subtype" %in% colnames(seurat_obj@meta.data)) {
    comparison <- compare_classifications(seurat_obj, config)
    attr(seurat_obj, "classification_comparison") <- comparison
  }
  
  # Step 4: Create diagnostic visualizations
  plots <- visualize_module_scores(seurat_obj, config)
  attr(seurat_obj, "classification_plots") <- plots
  
  message("\n========================================")
  message("Classification Complete")
  message("========================================\n")
  
  return(seurat_obj)
}
