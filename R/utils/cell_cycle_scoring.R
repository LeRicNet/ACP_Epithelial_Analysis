# ==============================================================================
# R/utils/cell_cycle_scoring.R
# ==============================================================================
# Multi-method cell cycle scoring with consensus determination
#
# Methods implemented:
#   1. Seurat CellCycleScoring - regression-based approach
#   2. Module score threshold - AddModuleScore with threshold classification
#   3. Cyclone (scran) - paired marker gene approach (if available)
# ==============================================================================

#' Score Cell Cycle Using Multiple Methods
#'
#' Applies multiple cell cycle scoring methods and determines consensus
#'
#' @param seurat_obj Seurat object
#' @param config Configuration list (optional)
#' @param methods Character vector of methods to use. Options: "seurat", "module_score", "cyclone"
#' @param s_genes S phase genes (uses Seurat defaults if NULL)
#' @param g2m_genes G2/M phase genes (uses Seurat defaults if NULL)
#' @param min_genes Minimum genes required per phase (default: 5)
#' @param seed Random seed for reproducibility
#' @return Seurat object with cell cycle scores and consensus phase
#' @export
score_cell_cycle_consensus <- function(seurat_obj,
                                       config = NULL,
                                       methods = c("seurat", "module_score", "cyclone"),
                                       s_genes = NULL,
                                       g2m_genes = NULL,
                                       min_genes = 5,
                                       seed = 42) {

  message("\n========================================")
  message("Multi-Method Cell Cycle Scoring")
  message("========================================\n")

  # Set seed
  set.seed(seed)

  # Get parameters from config if provided
  if (!is.null(config)) {
    if (!is.null(config$cell_cycle$methods)) {
      methods <- config$cell_cycle$methods
    }
    if (!is.null(config$cell_cycle$min_genes)) {
      min_genes <- config$cell_cycle$min_genes
    }
    if (!is.null(config$reproducibility$seed)) {
      seed <- config$reproducibility$seed
      set.seed(seed)
    }
  }

  # Get gene lists
  if (is.null(s_genes)) {
    s_genes <- Seurat::cc.genes.updated.2019$s.genes
  }
  if (is.null(g2m_genes)) {
    g2m_genes <- Seurat::cc.genes.updated.2019$g2m.genes
  }

  # Check gene availability
  available_genes <- rownames(seurat_obj)
  s_present <- intersect(s_genes, available_genes)
  g2m_present <- intersect(g2m_genes, available_genes)

  message(sprintf("S phase genes: %d/%d present", length(s_present), length(s_genes)))
  message(sprintf("G2M phase genes: %d/%d present", length(g2m_present), length(g2m_genes)))

  if (length(s_present) < min_genes || length(g2m_present) < min_genes) {
    warning(sprintf("Insufficient cell cycle genes (need >= %d per phase)", min_genes))
    seurat_obj$Phase <- "Unknown"
    seurat_obj$cc_consensus <- "Unknown"
    seurat_obj$cc_confidence <- 0
    return(seurat_obj)
  }

  # Initialize results storage
  n_cells <- ncol(seurat_obj)
  phase_calls <- list()
  methods_used <- c()

  # ===========================================================================
  # Method 1: Seurat CellCycleScoring
  # ===========================================================================
  if ("seurat" %in% methods) {
    message("\n--- Method 1: Seurat CellCycleScoring ---")

    result <- tryCatch({
      # CellCycleScoring uses AddModuleScore internally, which can fail with
      # default nbin=24 on smaller datasets. We'll try with adjusted nbin.
      # Unfortunately CellCycleScoring doesn't expose nbin directly, so we
      # pre-calculate the scores ourselves if needed.

      n_genes <- nrow(seurat_obj)

      # Try standard CellCycleScoring first
      seurat_obj <- tryCatch({
        Seurat::CellCycleScoring(
          seurat_obj,
          s.features = s_present,
          g2m.features = g2m_present,
          set.ident = FALSE
        )
      }, error = function(e) {
        # If it fails (likely nbin issue), calculate scores manually
        message("  Standard scoring failed, using manual calculation...")

        # Calculate S and G2M scores with safe nbin
        obj <- score_module_safe(seurat_obj, s_present, "S.Score", seed)
        obj <- score_module_safe(obj, g2m_present, "G2M.Score", seed)

        # Assign phase based on scores (matching Seurat's logic)
        s_score <- obj@meta.data$S.Score
        g2m_score <- obj@meta.data$G2M.Score

        phase <- rep("G1", ncol(obj))
        phase[s_score > g2m_score & s_score > 0] <- "S"
        phase[g2m_score > s_score & g2m_score > 0] <- "G2M"

        obj$Phase <- phase
        obj
      })

      phase_calls$seurat <- seurat_obj$Phase
      methods_used <- c(methods_used, "seurat")

      message("  Seurat phase distribution:")
      print(table(seurat_obj$Phase))

      TRUE
    }, error = function(e) {
      warning(paste("  Seurat scoring failed:", e$message))
      FALSE
    })
  }

  # ===========================================================================
  # Method 2: Module Score Threshold
  # ===========================================================================
  if ("module_score" %in% methods) {
    message("\n--- Method 2: Module Score Threshold ---")

    result <- tryCatch({
      # Calculate S phase score
      seurat_obj <- score_module_safe(
        seurat_obj,
        features = s_present,
        name = "cc_S_score",
        seed = seed
      )

      # Calculate G2M phase score
      seurat_obj <- score_module_safe(
        seurat_obj,
        features = g2m_present,
        name = "cc_G2M_score",
        seed = seed
      )

      # Get scores
      s_scores <- seurat_obj@meta.data$cc_S_score
      g2m_scores <- seurat_obj@meta.data$cc_G2M_score

      # Use adaptive thresholding based on score distribution
      # Default threshold from config, but adapt if scores are centered below it
      base_threshold <- 0.1
      if (!is.null(config) && !is.null(config$cell_cycle$score_threshold)) {
        base_threshold <- config$cell_cycle$score_threshold
      }

      # Adaptive: use mean + 0.5*SD as threshold if scores are low
      # This captures cells that are notably above average
      s_adaptive <- mean(s_scores, na.rm = TRUE) + 0.5 * sd(s_scores, na.rm = TRUE)
      g2m_adaptive <- mean(g2m_scores, na.rm = TRUE) + 0.5 * sd(g2m_scores, na.rm = TRUE)

      # Use the more permissive of base or adaptive threshold
      s_threshold <- min(base_threshold, s_adaptive)
      g2m_threshold <- min(base_threshold, g2m_adaptive)

      message(sprintf("  S threshold: %.3f (adaptive: %.3f)", s_threshold, s_adaptive))
      message(sprintf("  G2M threshold: %.3f (adaptive: %.3f)", g2m_threshold, g2m_adaptive))

      # Classify based on scores
      phase_module <- rep("G1", n_cells)

      # S phase: S score > threshold AND S > G2M
      s_phase <- s_scores > s_threshold & s_scores > g2m_scores
      phase_module[s_phase] <- "S"

      # G2M phase: G2M score > threshold AND G2M > S
      g2m_phase <- g2m_scores > g2m_threshold & g2m_scores > s_scores
      phase_module[g2m_phase] <- "G2M"

      phase_calls$module_score <- phase_module
      methods_used <- c(methods_used, "module_score")

      message("  Module score phase distribution:")
      print(table(phase_module))

      TRUE
    }, error = function(e) {
      warning(paste("  Module score method failed:", e$message))
      FALSE
    })
  }

  # ===========================================================================
  # Method 3: Cyclone (scran) - optional
  # ===========================================================================
  if ("cyclone" %in% methods) {
    message("\n--- Method 3: Cyclone (scran) ---")

    if (!requireNamespace("scran", quietly = TRUE)) {
      message("  scran package not available - skipping Cyclone method")
      message("  Install with: BiocManager::install('scran')")
    } else {
      result <- tryCatch({
        # Get expression matrix (SeuratObject 5.0+ compatible)
        expr_matrix <- tryCatch({
          Seurat::GetAssayData(seurat_obj, layer = "counts")
        }, error = function(e) {
          # Fallback for older Seurat versions
          Seurat::GetAssayData(seurat_obj, slot = "counts")
        })

        # Get human cell cycle markers from scran
        hs_pairs <- tryCatch({
          readRDS(system.file("exdata", "human_cycle_markers.rds",
                              package = "scran"))
        }, error = function(e) {
          stop("Could not load cell cycle markers from scran")
        })

        # Check if we need to convert gene symbols to Ensembl IDs
        # scran markers use Ensembl IDs (ENSG...)
        current_genes <- rownames(expr_matrix)
        uses_symbols <- !any(grepl("^ENSG", current_genes[1:min(100, length(current_genes))]))

        if (uses_symbols) {
          message("  Converting gene symbols to Ensembl IDs...")

          # Try to load annotation package
          if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
            stop("org.Hs.eg.db package required for gene ID conversion. ",
                 "Install with: BiocManager::install('org.Hs.eg.db')")
          }

          if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
            stop("AnnotationDbi package required for gene ID conversion. ",
                 "Install with: BiocManager::install('AnnotationDbi')")
          }

          # Map symbols to Ensembl IDs
          symbol_to_ensembl <- tryCatch({
            AnnotationDbi::mapIds(
              org.Hs.eg.db::org.Hs.eg.db,
              keys = current_genes,
              keytype = "SYMBOL",
              column = "ENSEMBL",
              multiVals = "first"
            )
          }, error = function(e) {
            stop("Failed to map gene symbols to Ensembl IDs: ", e$message)
          })

          # Count successful mappings
          n_mapped <- sum(!is.na(symbol_to_ensembl))
          message(sprintf("  Mapped %d/%d genes (%.1f%%)",
                          n_mapped, length(current_genes),
                          n_mapped/length(current_genes) * 100))

          if (n_mapped < 1000) {
            stop("Too few genes mapped to Ensembl IDs (", n_mapped, ")")
          }

          # Create new expression matrix with Ensembl IDs
          # Keep only genes that mapped successfully
          mapped_idx <- !is.na(symbol_to_ensembl)
          expr_matrix_mapped <- expr_matrix[mapped_idx, , drop = FALSE]
          rownames(expr_matrix_mapped) <- symbol_to_ensembl[mapped_idx]

          # Remove any duplicate Ensembl IDs (keep first occurrence)
          dup_ids <- duplicated(rownames(expr_matrix_mapped))
          if (any(dup_ids)) {
            message(sprintf("  Removing %d duplicate Ensembl IDs", sum(dup_ids)))
            expr_matrix_mapped <- expr_matrix_mapped[!dup_ids, , drop = FALSE]
          }

          expr_matrix <- expr_matrix_mapped
          message(sprintf("  Final matrix: %d genes x %d cells",
                          nrow(expr_matrix), ncol(expr_matrix)))
        }

        # Check overlap with marker pairs
        # Handle different scran versions - pairs structure may vary
        all_marker_genes <- tryCatch({
          # Try newer format (list of data frames)
          if (is.data.frame(hs_pairs$G1) || is.matrix(hs_pairs$G1)) {
            unique(c(
              as.character(hs_pairs$G1[,1]), as.character(hs_pairs$G1[,2]),
              as.character(hs_pairs$S[,1]), as.character(hs_pairs$S[,2]),
              as.character(hs_pairs$G2M[,1]), as.character(hs_pairs$G2M[,2])
            ))
          } else if (is.list(hs_pairs$G1)) {
            # Older format (list of lists with $first/$second)
            unique(c(
              unlist(lapply(hs_pairs$G1, function(x) c(x$first, x$second))),
              unlist(lapply(hs_pairs$S, function(x) c(x$first, x$second))),
              unlist(lapply(hs_pairs$G2M, function(x) c(x$first, x$second)))
            ))
          } else {
            # Fallback: try to extract all character elements
            unique(unlist(hs_pairs))
          }
        }, error = function(e) {
          message("  Warning: Could not parse marker pair structure")
          character(0)
        })

        if (length(all_marker_genes) > 0) {
          marker_overlap <- sum(all_marker_genes %in% rownames(expr_matrix))
          message(sprintf("  Marker gene overlap: %d/%d (%.1f%%)",
                          marker_overlap, length(all_marker_genes),
                          marker_overlap/length(all_marker_genes) * 100))

          if (marker_overlap < 50) {
            stop("Insufficient marker gene overlap (", marker_overlap, " genes)")
          }
        }

        # Run cyclone - convert to SingleCellExperiment if needed
        message("  Running cyclone classifier...")

        # Set up parallel processing for HPC
        n_cores <- 1
        if (!is.null(config) && !is.null(config$reproducibility$n_cores)) {
          n_cores <- config$reproducibility$n_cores
        } else {
          # Auto-detect available cores (use half to be safe)
          n_cores <- max(1, floor(parallel::detectCores() / 2))
        }

        # Set up BiocParallel backend
        BPPARAM <- BiocParallel::SerialParam()  # Default: serial

        if (n_cores > 1 && requireNamespace("BiocParallel", quietly = TRUE)) {
          message(sprintf("  Using %d cores for parallel processing", n_cores))

          # Choose appropriate backend
          if (.Platform$OS.type == "unix") {
            # Unix/Linux/Mac - use MulticoreParam (fork-based, faster)
            BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores)
          } else {
            # Windows - use SnowParam (socket-based)
            BPPARAM <- BiocParallel::SnowParam(workers = n_cores)
          }
        } else {
          message("  Running in serial mode (set n_cores in config for parallel)")
        }

        # scran::cyclone expects a SingleCellExperiment or matrix
        if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
          sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = expr_matrix)
          )
          cyclone_result <- scran::cyclone(sce, pairs = hs_pairs, BPPARAM = BPPARAM)
        } else {
          # Pass matrix directly
          cyclone_result <- scran::cyclone(expr_matrix, pairs = hs_pairs, BPPARAM = BPPARAM)
        }

        # Get phase assignments - handle different result formats
        phase_cyclone <- if ("phases" %in% names(cyclone_result)) {
          cyclone_result$phases
        } else if ("phase" %in% names(cyclone_result)) {
          cyclone_result$phase
        } else {
          # Try to determine phase from scores
          scores_df <- cyclone_result$scores
          if (is.null(scores_df)) {
            scores_df <- cyclone_result$normalized.scores
          }
          if (!is.null(scores_df)) {
            apply(scores_df, 1, function(x) names(x)[which.max(x)])
          } else {
            stop("Could not extract phases from cyclone result")
          }
        }

        # Check if cyclone returned valid results
        if (is.null(phase_cyclone) || length(phase_cyclone) == 0 ||
            all(is.na(phase_cyclone)) || length(unique(na.omit(phase_cyclone))) == 0) {
          stop("Cyclone returned empty or invalid phase assignments")
        }

        phase_calls$cyclone <- phase_cyclone
        methods_used <- c(methods_used, "cyclone")

        # Store cyclone scores - handle different result formats
        scores_df <- cyclone_result$scores
        if (is.null(scores_df)) {
          scores_df <- cyclone_result$normalized.scores
        }

        if (!is.null(scores_df) && is.data.frame(scores_df)) {
          if ("G1" %in% colnames(scores_df)) seurat_obj$cc_cyclone_G1 <- scores_df$G1
          if ("S" %in% colnames(scores_df)) seurat_obj$cc_cyclone_S <- scores_df$S
          if ("G2M" %in% colnames(scores_df)) seurat_obj$cc_cyclone_G2M <- scores_df$G2M
        }

        message("  Cyclone phase distribution:")
        print(table(phase_cyclone))

        TRUE
      }, error = function(e) {
        warning(paste("  Cyclone method failed:", e$message))
        FALSE
      })
    }
  }

  # ===========================================================================
  # Consensus Determination
  # ===========================================================================
  message("\n--- Determining Consensus ---")
  message(sprintf("Methods used: %s", paste(methods_used, collapse = ", ")))

  if (length(methods_used) == 0) {
    warning("No cell cycle scoring methods succeeded")
    seurat_obj$Phase <- "Unknown"
    seurat_obj$cc_consensus <- "Unknown"
    seurat_obj$cc_confidence <- 0
    return(seurat_obj)
  }

  if (length(methods_used) == 1) {
    # Only one method - use it directly
    seurat_obj$cc_consensus <- phase_calls[[methods_used[1]]]
    seurat_obj$cc_confidence <- 1.0
    seurat_obj$cc_n_methods <- 1
    message("  Single method used - no consensus needed")

  } else {
    # Multiple methods - determine consensus
    consensus_result <- determine_consensus(phase_calls, methods_used)

    seurat_obj$cc_consensus <- consensus_result$consensus
    seurat_obj$cc_confidence <- consensus_result$confidence
    seurat_obj$cc_n_methods <- length(methods_used)

    message(sprintf("  Consensus phase distribution:"))
    print(table(consensus_result$consensus))

    message(sprintf("\n  Mean confidence: %.2f", mean(consensus_result$confidence)))
    message(sprintf("  Cells with full agreement: %d (%.1f%%)",
                    sum(consensus_result$confidence == 1.0),
                    mean(consensus_result$confidence == 1.0) * 100))
  }

  # Store individual method calls in metadata
  for (method in methods_used) {
    col_name <- paste0("cc_", method)
    seurat_obj@meta.data[[col_name]] <- phase_calls[[method]]
  }

  # Set Phase to consensus (for compatibility with downstream code)
  seurat_obj$Phase <- seurat_obj$cc_consensus

  message("\n  Final Phase assignment (consensus):")
  print(table(seurat_obj$Phase))

  return(seurat_obj)
}


#' Determine Consensus Phase from Multiple Methods
#'
#' @param phase_calls Named list of phase assignments from each method
#' @param methods_used Character vector of method names used
#' @return List with consensus phase and confidence score
#' @keywords internal
determine_consensus <- function(phase_calls, methods_used) {

  n_cells <- length(phase_calls[[1]])
  n_methods <- length(methods_used)

  consensus <- character(n_cells)
  confidence <- numeric(n_cells)

  for (i in seq_len(n_cells)) {
    # Get calls for this cell
    calls <- sapply(methods_used, function(m) phase_calls[[m]][i])

    # Count votes for each phase
    vote_counts <- table(calls)
    max_votes <- max(vote_counts)
    winner <- names(vote_counts)[which.max(vote_counts)]

    consensus[i] <- winner
    confidence[i] <- max_votes / n_methods
  }

  return(list(
    consensus = consensus,
    confidence = confidence
  ))
}


#' Safe Module Scoring with Automatic nbin Adjustment
#'
#' @param seurat_obj Seurat object
#' @param features Gene list
#' @param name Score name
#' @param seed Random seed
#' @return Seurat object with score added
#' @keywords internal
score_module_safe <- function(seurat_obj, features, name, seed = 42) {

  n_genes <- nrow(seurat_obj)
  nbin_values <- c(24, 12, 8, 5)
  nbin_values <- nbin_values[nbin_values < n_genes / 10]
  if (length(nbin_values) == 0) nbin_values <- 5

  for (nbin in nbin_values) {
    result <- tryCatch({
      obj <- Seurat::AddModuleScore(
        seurat_obj,
        features = list(features),
        name = name,
        nbin = nbin,
        seed = seed
      )
      list(success = TRUE, obj = obj)
    }, error = function(e) {
      list(success = FALSE)
    })

    if (result$success) {
      # Rename column (Seurat appends "1")
      old_name <- paste0(name, "1")
      if (old_name %in% colnames(result$obj@meta.data)) {
        result$obj@meta.data[[name]] <- result$obj@meta.data[[old_name]]
        result$obj@meta.data[[old_name]] <- NULL
      }
      return(result$obj)
    }
  }

  stop(paste("Failed to calculate module score for", name))
}


#' Generate Cell Cycle Diagnostic Plot
#'
#' Creates diagnostic visualization of cell cycle scoring results
#'
#' @param seurat_obj Seurat object with cell cycle scores
#' @param config Configuration list (optional, for colors)
#' @return ggplot object
#' @export
plot_cell_cycle_diagnostics <- function(seurat_obj, config = NULL) {

  require(ggplot2)
  require(patchwork)

  # Get colors
  if (!is.null(config) && !is.null(config$visualization$colors$cell_cycle)) {
    cycle_colors <- unlist(config$visualization$colors$cell_cycle)
  } else {
    cycle_colors <- c(G1 = "#F8766D", S = "#00BA38", G2M = "#619CFF")
  }

  plots <- list()

  # Plot 1: Consensus phase distribution
  phase_counts <- as.data.frame(table(seurat_obj$cc_consensus))
  colnames(phase_counts) <- c("Phase", "Count")
  phase_counts$Pct <- phase_counts$Count / sum(phase_counts$Count) * 100

  plots$consensus <- ggplot(phase_counts, aes(x = Phase, y = Count, fill = Phase)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f%%", Pct)), vjust = -0.3) +
    scale_fill_manual(values = cycle_colors) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Consensus Phase Distribution", x = "", y = "Cells")

  # Plot 2: Confidence distribution
  if ("cc_confidence" %in% colnames(seurat_obj@meta.data)) {
    plots$confidence <- ggplot(seurat_obj@meta.data, aes(x = cc_confidence)) +
      geom_histogram(bins = 20, fill = "steelblue", color = "white") +
      theme_minimal() +
      labs(title = "Consensus Confidence", x = "Confidence", y = "Cells")
  }

  # Plot 3: S vs G2M scores (if available)
  if (all(c("S.Score", "G2M.Score") %in% colnames(seurat_obj@meta.data))) {
    plots$scores <- ggplot(seurat_obj@meta.data,
                           aes(x = S.Score, y = G2M.Score, color = cc_consensus)) +
      geom_point(size = 0.5, alpha = 0.5) +
      scale_color_manual(values = cycle_colors) +
      theme_minimal() +
      labs(title = "Seurat S vs G2M Scores", color = "Phase")
  }

  # Plot 4: Method agreement (if multiple methods)
  method_cols <- grep("^cc_seurat$|^cc_module_score$|^cc_cyclone$",
                      colnames(seurat_obj@meta.data), value = TRUE)

  if (length(method_cols) >= 2) {
    # Calculate agreement matrix
    agreement_data <- data.frame(
      Cell = colnames(seurat_obj)
    )
    for (col in method_cols) {
      agreement_data[[col]] <- seurat_obj@meta.data[[col]]
    }

    # Simple bar plot of agreement
    n_methods <- length(method_cols)
    agreement_counts <- sapply(seq_len(ncol(seurat_obj)), function(i) {
      calls <- sapply(method_cols, function(m) seurat_obj@meta.data[[m]][i])
      max(table(calls))
    })

    agree_df <- data.frame(
      Agreement = factor(agreement_counts, levels = 1:n_methods),
      Count = 1
    ) %>%
      dplyr::group_by(Agreement) %>%
      dplyr::summarise(Count = dplyr::n(), .groups = "drop")

    plots$agreement <- ggplot(agree_df, aes(x = Agreement, y = Count)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      labs(title = "Method Agreement",
           x = paste("Methods Agreeing (of", n_methods, ")"),
           y = "Cells")
  }

  # Combine plots
  if (length(plots) >= 4) {
    combined <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]]) +
      plot_annotation(title = "Cell Cycle Scoring Diagnostics")
  } else if (length(plots) >= 2) {
    combined <- wrap_plots(plots) +
      plot_annotation(title = "Cell Cycle Scoring Diagnostics")
  } else {
    combined <- plots[[1]]
  }

  return(combined)
}


#' Evaluate Cell Cycle Scoring Statistics
#'
#' Computes comprehensive statistics to evaluate cell cycle scoring quality
#'
#' @param seurat_obj Seurat object with cell cycle scores
#' @param config Configuration list (optional)
#' @param proliferation_markers Genes to validate against (default: MKI67, TOP2A, PCNA)
#' @return List with statistical evaluation results
#' @export
evaluate_cell_cycle_statistics <- function(seurat_obj,
                                           config = NULL,
                                           proliferation_markers = c("MKI67", "TOP2A", "PCNA", "STMN1")) {

  message("\n========================================")
  message("Cell Cycle Scoring Statistics")
  message("========================================\n")

  results <- list()

  # Get method columns

  method_cols <- grep("^cc_seurat$|^cc_module_score$|^cc_cyclone$",
                      colnames(seurat_obj@meta.data), value = TRUE)
  n_methods <- length(method_cols)

  message(sprintf("Methods evaluated: %s", paste(gsub("cc_", "", method_cols), collapse = ", ")))

  # ===========================================================================
  # 1. Inter-rater Agreement Statistics
  # ===========================================================================
  message("\n--- Inter-Method Agreement ---")

  if (n_methods >= 2) {
    # Build classification matrix
    class_matrix <- as.data.frame(seurat_obj@meta.data[, method_cols])

    # Pairwise Cohen's Kappa
    kappa_results <- list()

    for (i in 1:(n_methods - 1)) {
      for (j in (i + 1):n_methods) {
        method1 <- method_cols[i]
        method2 <- method_cols[j]

        # Confusion matrix
        tbl <- table(class_matrix[[method1]], class_matrix[[method2]])

        # Cohen's Kappa calculation
        n <- sum(tbl)
        p_observed <- sum(diag(tbl)) / n

        row_marginals <- rowSums(tbl) / n
        col_marginals <- colSums(tbl) / n
        p_expected <- sum(row_marginals * col_marginals)

        kappa <- (p_observed - p_expected) / (1 - p_expected)

        # Percent agreement
        pct_agree <- p_observed * 100

        pair_name <- paste(gsub("cc_", "", method1), "vs", gsub("cc_", "", method2))
        kappa_results[[pair_name]] <- list(
          kappa = kappa,
          percent_agreement = pct_agree,
          confusion_matrix = tbl
        )

        message(sprintf("  %s: κ = %.3f, Agreement = %.1f%%",
                        pair_name, kappa, pct_agree))
      }
    }

    results$pairwise_kappa <- kappa_results

    # Fleiss' Kappa for all methods (if 3+)
    if (n_methods >= 3) {
      fleiss_kappa <- tryCatch({
        calculate_fleiss_kappa(class_matrix)
      }, error = function(e) {
        message(sprintf("  Could not calculate Fleiss' Kappa: %s", e$message))
        NA
      })

      if (!is.na(fleiss_kappa)) {
        results$fleiss_kappa <- fleiss_kappa
        message(sprintf("  Fleiss' Kappa (all methods): κ = %.3f", fleiss_kappa))
      }
    }

    # Overall agreement rate
    all_agree <- apply(class_matrix, 1, function(x) length(unique(x)) == 1)
    results$full_agreement_rate <- mean(all_agree)
    message(sprintf("  Full agreement (all methods): %.1f%%", mean(all_agree) * 100))

    # Check for method outliers (one method disagreeing with others)
    if (n_methods >= 3) {
      # Calculate mean kappa for each method with others
      method_agreement <- sapply(method_cols, function(m) {
        other_kappas <- sapply(names(kappa_results), function(pair) {
          if (grepl(gsub("cc_", "", m), pair)) {
            return(kappa_results[[pair]]$kappa)
          }
          return(NA)
        })
        mean(other_kappas, na.rm = TRUE)
      })
      names(method_agreement) <- gsub("cc_", "", method_cols)

      # Identify outlier method (mean kappa < 0.2)
      outlier_methods <- names(method_agreement)[method_agreement < 0.2]
      if (length(outlier_methods) > 0) {
        results$outlier_methods <- outlier_methods
        message(sprintf("\n  ⚠ Method disagreement detected: %s has low agreement with others",
                        paste(outlier_methods, collapse = ", ")))
        message("    Consider running without this method or investigating differences")
      }
    }

  } else {
    message("  Only one method - skipping agreement statistics")
  }

  # ===========================================================================
  # 2. Biological Validation
  # ===========================================================================
  message("\n--- Biological Validation ---")

  # Check proliferation marker expression by phase
  available_markers <- intersect(proliferation_markers, rownames(seurat_obj))

  if (length(available_markers) > 0) {
    message(sprintf("  Validating against: %s", paste(available_markers, collapse = ", ")))

    # Get expression data (handle both Seurat v4 and v5)
    expr_data <- tryCatch({
      as.matrix(Seurat::GetAssayData(seurat_obj, layer = "data")[available_markers, , drop = FALSE])
    }, error = function(e) {
      tryCatch({
        as.matrix(Seurat::GetAssayData(seurat_obj, slot = "data")[available_markers, , drop = FALSE])
      }, error = function(e2) {
        NULL
      })
    })

    if (!is.null(expr_data)) {
      # Mean expression by phase
      phases <- as.character(seurat_obj$cc_consensus)

      # Check we have valid phases
      valid_phases <- phases %in% c("G1", "S", "G2M")
      if (sum(valid_phases) < 100) {
        message("  Not enough cells with valid phase assignments")
      } else {
        marker_by_phase <- lapply(available_markers, function(marker) {
          tapply(expr_data[marker, valid_phases], phases[valid_phases], mean, na.rm = TRUE)
        })
        names(marker_by_phase) <- available_markers

        results$marker_expression_by_phase <- marker_by_phase

        # Calculate fold change and test for each marker
        validation_list <- lapply(available_markers, function(marker) {
          expr_by_phase <- marker_by_phase[[marker]]

          # Get values safely
          s_val <- if ("S" %in% names(expr_by_phase)) expr_by_phase["S"] else NA
          g2m_val <- if ("G2M" %in% names(expr_by_phase)) expr_by_phase["G2M"] else NA
          g1_val <- if ("G1" %in% names(expr_by_phase)) expr_by_phase["G1"] else NA

          cycling_mean <- mean(c(s_val, g2m_val), na.rm = TRUE)

          # Calculate fold change
          fold_change <- NA
          if (!is.na(g1_val) && !is.na(cycling_mean) && g1_val > 0) {
            fold_change <- cycling_mean / g1_val
          }

          # Statistical test: Wilcoxon
          p_value <- NA
          cycling_cells <- which(phases %in% c("S", "G2M"))
          g1_cells <- which(phases == "G1")

          if (length(cycling_cells) >= 10 && length(g1_cells) >= 10) {
            test <- tryCatch({
              wilcox.test(expr_data[marker, cycling_cells],
                          expr_data[marker, g1_cells])
            }, error = function(e) NULL)
            if (!is.null(test)) p_value <- test$p.value
          }

          data.frame(
            marker = marker,
            fold_change = fold_change,
            p_value = p_value,
            stringsAsFactors = FALSE
          )
        })

        marker_df <- do.call(rbind, validation_list)
        marker_df$significant <- !is.na(marker_df$p_value) & marker_df$p_value < 0.05
        rownames(marker_df) <- marker_df$marker

        results$marker_validation <- marker_df

        message("\n  Proliferation marker enrichment in cycling cells:")
        for (i in 1:nrow(marker_df)) {
          fc <- marker_df$fold_change[i]
          pval <- marker_df$p_value[i]
          marker <- marker_df$marker[i]
          if (!is.na(fc) && !is.na(pval)) {
            sig <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "")))
            message(sprintf("    %s: FC = %.2f, p = %.2e %s", marker, fc, pval, sig))
          } else {
            message(sprintf("    %s: Could not calculate", marker))
          }
        }

        # Overall validation score
        mean_fc <- mean(marker_df$fold_change, na.rm = TRUE)
        n_significant <- sum(marker_df$significant, na.rm = TRUE)

        is_validated <- FALSE
        if (!is.na(mean_fc) && mean_fc > 1.5) {
          if (n_significant >= length(available_markers) / 2) {
            is_validated <- TRUE
          }
        }

        results$validation_summary <- list(
          mean_fold_change = mean_fc,
          n_markers_significant = n_significant,
          n_markers_tested = length(available_markers),
          validated = is_validated
        )

        message(sprintf("\n  Validation summary: Mean FC = %.2f, %d/%d markers significant",
                        mean_fc, n_significant, length(available_markers)))
        message(sprintf("  Biological validation: %s",
                        ifelse(is_validated, "PASSED", "REVIEW NEEDED")))
      }

    } else {
      message("  Could not retrieve expression data")
    }

  } else {
    message("  No proliferation markers found in dataset")
  }

  # ===========================================================================
  # 3. Phase Distribution Statistics
  # ===========================================================================
  message("\n--- Phase Distribution ---")

  phase_table <- table(seurat_obj$cc_consensus)
  phase_pct <- prop.table(phase_table) * 100

  results$phase_distribution <- data.frame(
    Phase = names(phase_table),
    Count = as.numeric(phase_table),
    Percent = as.numeric(phase_pct)
  )

  # Expected ranges for typical tissues (rough guidelines)
  expected_ranges <- list(
    G1 = c(50, 80),
    S = c(10, 30),
    G2M = c(10, 25)
  )

  message("  Phase distribution vs expected ranges:")
  for (phase in c("G1", "S", "G2M")) {
    pct <- phase_pct[phase]
    expected <- expected_ranges[[phase]]
    in_range <- pct >= expected[1] && pct <= expected[2]
    status <- ifelse(in_range, "OK", ifelse(pct < expected[1], "LOW", "HIGH"))
    message(sprintf("    %s: %.1f%% (expected: %d-%d%%) [%s]",
                    phase, pct, expected[1], expected[2], status))
  }

  # ===========================================================================
  # 4. Method-Specific Statistics
  # ===========================================================================
  message("\n--- Method-Specific Statistics ---")

  method_stats <- list()

  for (method in method_cols) {
    method_name <- gsub("cc_", "", method)
    phases <- seurat_obj@meta.data[[method]]

    # Phase distribution
    dist <- table(phases)
    pct <- prop.table(dist) * 100

    # Cycling percentage
    cycling_pct <- sum(pct[names(pct) %in% c("S", "G2M")])

    method_stats[[method_name]] <- list(
      distribution = dist,
      percent = pct,
      cycling_percent = cycling_pct
    )

    message(sprintf("  %s: G1=%.1f%%, S=%.1f%%, G2M=%.1f%% (Cycling=%.1f%%)",
                    method_name,
                    pct["G1"], pct["S"], pct["G2M"], cycling_pct))
  }

  results$method_statistics <- method_stats

  # ===========================================================================
  # 5. Confidence Statistics
  # ===========================================================================
  if ("cc_confidence" %in% colnames(seurat_obj@meta.data)) {
    message("\n--- Confidence Statistics ---")

    conf <- seurat_obj$cc_confidence

    results$confidence_stats <- list(
      mean = mean(conf, na.rm = TRUE),
      median = median(conf, na.rm = TRUE),
      sd = sd(conf, na.rm = TRUE),
      q25 = quantile(conf, 0.25, na.rm = TRUE),
      q75 = quantile(conf, 0.75, na.rm = TRUE),
      high_confidence_pct = mean(conf >= 0.8, na.rm = TRUE) * 100,
      low_confidence_pct = mean(conf < 0.5, na.rm = TRUE) * 100
    )

    message(sprintf("  Mean confidence: %.3f ± %.3f",
                    results$confidence_stats$mean, results$confidence_stats$sd))
    message(sprintf("  High confidence (≥0.8): %.1f%%", results$confidence_stats$high_confidence_pct))
    message(sprintf("  Low confidence (<0.5): %.1f%%", results$confidence_stats$low_confidence_pct))

    # Confidence by phase
    conf_by_phase <- tapply(conf, seurat_obj$cc_consensus, mean, na.rm = TRUE)
    results$confidence_by_phase <- conf_by_phase
    message("  Confidence by phase:")
    for (phase in names(conf_by_phase)) {
      message(sprintf("    %s: %.3f", phase, conf_by_phase[phase]))
    }
  }

  # ===========================================================================
  # 6. Summary Score
  # ===========================================================================
  message("\n--- Overall Quality Score ---")

  quality_components <- c()

  # Component 1: Method agreement (weight: 30%)
  if (!is.null(results$full_agreement_rate)) {
    agreement_score <- min(results$full_agreement_rate / 0.5, 1)  # 50% agreement = max score
    quality_components["agreement"] <- agreement_score * 0.3
  }

  # Component 2: Biological validation (weight: 40%)
  if (!is.null(results$validation_summary)) {
    fc_score <- min((results$validation_summary$mean_fold_change - 1) / 1, 1)  # FC of 2 = max score
    fc_score <- max(fc_score, 0)
    sig_score <- results$validation_summary$n_markers_significant / results$validation_summary$n_markers_tested
    validation_score <- (fc_score + sig_score) / 2
    quality_components["validation"] <- validation_score * 0.4
  }

  # Component 3: Confidence (weight: 30%)
  if (!is.null(results$confidence_stats)) {
    conf_score <- results$confidence_stats$mean
    quality_components["confidence"] <- conf_score * 0.3
  }

  overall_score <- sum(quality_components, na.rm = TRUE)
  results$quality_score <- list(
    components = quality_components,
    overall = overall_score,
    grade = ifelse(overall_score >= 0.8, "Excellent",
                   ifelse(overall_score >= 0.6, "Good",
                          ifelse(overall_score >= 0.4, "Acceptable", "Review Needed")))
  )

  message(sprintf("\n  Quality Score: %.2f / 1.00 (%s)", overall_score, results$quality_score$grade))
  message("  Components:")
  for (comp in names(quality_components)) {
    message(sprintf("    %s: %.3f", comp, quality_components[comp]))
  }

  message("\n========================================\n")

  return(results)
}


#' Calculate Fleiss' Kappa for Multiple Raters
#'
#' @param ratings Data frame with raters as columns, subjects as rows
#' @return Fleiss' kappa value
#' @keywords internal
calculate_fleiss_kappa <- function(ratings) {

  # Convert to character matrix and handle NAs
  ratings <- as.data.frame(lapply(ratings, as.character))

  # Remove rows with any NA values
  complete_rows <- complete.cases(ratings)
  if (sum(complete_rows) < nrow(ratings) * 0.5) {
    warning("More than 50% of rows have NA values - Fleiss' kappa may be unreliable")
  }
  ratings <- ratings[complete_rows, , drop = FALSE]

  if (nrow(ratings) < 10) {
    warning("Too few complete cases for Fleiss' kappa")
    return(NA)
  }

  # Get unique categories (excluding NA)
  categories <- unique(unlist(ratings))
  categories <- categories[!is.na(categories)]

  n_subjects <- nrow(ratings)
  n_raters <- ncol(ratings)

  # Build count matrix (subjects x categories)
  count_matrix <- matrix(0, nrow = n_subjects, ncol = length(categories))
  colnames(count_matrix) <- categories

  for (i in 1:n_subjects) {
    for (j in 1:n_raters) {
      cat <- as.character(ratings[i, j])
      if (!is.na(cat) && cat %in% categories) {
        count_matrix[i, cat] <- count_matrix[i, cat] + 1
      }
    }
  }

  # Calculate P_i (agreement for each subject)
  P_i <- (rowSums(count_matrix^2) - n_raters) / (n_raters * (n_raters - 1))
  P_bar <- mean(P_i)

  # Calculate P_j (proportion in each category)
  p_j <- colSums(count_matrix) / (n_subjects * n_raters)
  P_e <- sum(p_j^2)

  # Fleiss' kappa
  if (P_e >= 1) {
    return(NA)  # Cannot calculate if expected agreement is 100%
  }

  kappa <- (P_bar - P_e) / (1 - P_e)

  return(kappa)
}


#' Generate Cell Cycle Statistics Report
#'
#' Creates a formatted report of cell cycle scoring statistics
#'
#' @param stats_results Results from evaluate_cell_cycle_statistics
#' @param output_path Path to save the report (optional)
#' @return Character string with formatted report
#' @export
generate_cell_cycle_report <- function(stats_results, output_path = NULL) {

  report <- c(
    "# Cell Cycle Scoring Evaluation Report",
    paste0("Generated: ", Sys.time()),
    "",
    "## Quality Score",
    sprintf("**Overall: %.2f / 1.00 (%s)**",
            stats_results$quality_score$overall,
            stats_results$quality_score$grade),
    ""
  )

  # Method agreement section
  if (!is.null(stats_results$pairwise_kappa)) {
    report <- c(report, "## Inter-Method Agreement", "")
    report <- c(report, "| Comparison | Cohen's κ | Agreement |")
    report <- c(report, "|------------|-----------|-----------|")

    for (name in names(stats_results$pairwise_kappa)) {
      k <- stats_results$pairwise_kappa[[name]]
      report <- c(report, sprintf("| %s | %.3f | %.1f%% |", name, k$kappa, k$percent_agreement))
    }

    if (!is.null(stats_results$fleiss_kappa)) {
      report <- c(report, "", sprintf("**Fleiss' Kappa (all methods): %.3f**", stats_results$fleiss_kappa))
    }
    report <- c(report, "")
  }

  # Biological validation
  if (!is.null(stats_results$marker_validation)) {
    report <- c(report, "## Biological Validation", "")
    report <- c(report, "| Marker | Fold Change | p-value | Significant |")
    report <- c(report, "|--------|-------------|---------|-------------|")

    mv <- stats_results$marker_validation
    for (i in 1:nrow(mv)) {
      report <- c(report, sprintf("| %s | %.2f | %.2e | %s |",
                                  mv$marker[i], mv$fold_change[i], mv$p_value[i],
                                  ifelse(mv$significant[i], "Yes", "No")))
    }
    report <- c(report, "")
  }

  # Phase distribution
  report <- c(report, "## Phase Distribution", "")
  report <- c(report, "| Phase | Count | Percent |")
  report <- c(report, "|-------|-------|---------|")

  pd <- stats_results$phase_distribution
  for (i in 1:nrow(pd)) {
    report <- c(report, sprintf("| %s | %d | %.1f%% |", pd$Phase[i], pd$Count[i], pd$Percent[i]))
  }
  report <- c(report, "")

  # Confidence
  if (!is.null(stats_results$confidence_stats)) {
    cs <- stats_results$confidence_stats
    report <- c(report, "## Confidence Statistics", "",
                sprintf("- Mean: %.3f ± %.3f", cs$mean, cs$sd),
                sprintf("- Median: %.3f", cs$median),
                sprintf("- High confidence (≥0.8): %.1f%%", cs$high_confidence_pct),
                sprintf("- Low confidence (<0.5): %.1f%%", cs$low_confidence_pct),
                "")
  }

  report_text <- paste(report, collapse = "\n")

  # Save if path provided
  if (!is.null(output_path)) {
    writeLines(report_text, output_path)
    message(sprintf("Report saved to: %s", output_path))
  }

  return(report_text)
}
