#!/usr/bin/env Rscript
# ==============================================================================
# scripts/01_cell_type_annotation.R
# ==============================================================================
# Cell Type Annotation for ACP Spatial Transcriptomics Data
#
# This script performs module score-based classification of epithelial subtypes,
# replacing the original binary threshold approach.
#
# Key improvements over binary classification:
#   1. Continuous module scores instead of arbitrary thresholds
#   2. Explicit Transit-Amplifying identification via proliferation markers
#   3. Hybrid state detection for cells with high scores in multiple categories
#   4. Reduced "Unknown" category through more sensitive classification
#
# Input:
#   - Raw spatial Seurat object (config$paths$spatial_object)
#
# Output:
#   - Annotated Seurat object with module scores and classifications
#   - Classification summary tables
#   - Diagnostic plots
#
# Usage:
#   Rscript scripts/01_cell_type_annotation.R
#   Rscript scripts/01_cell_type_annotation.R --config path/to/config.yaml
#
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL

if (length(args) > 0) {
  if (args[1] == "--config" && length(args) > 1) {
    config_path <- args[2]
  } else if (!startsWith(args[1], "--")) {
    config_path <- args[1]
  }
}

# Set working directory to project root
# Try multiple methods to find project root
find_project_root <- function() {
  # Method 1: Already in project root
  if (file.exists("config/config.yaml")) {
    return(getwd())
  }

  # Method 2: In scripts/ subdirectory
  if (file.exists("../config/config.yaml")) {
    return(normalizePath(".."))
  }

  # Method 3: Use 'here' package if available
  if (requireNamespace("here", quietly = TRUE)) {
    root <- here::here()
    if (file.exists(file.path(root, "config/config.yaml"))) {
      return(root)
    }
  }

  # Method 4: Search up directory tree
  current <- getwd()
  for (i in 1:5) {
    if (file.exists(file.path(current, "config/config.yaml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (parent == current) break  # Reached filesystem root
    current <- parent
  }

  return(NULL)
}

project_root <- find_project_root()
if (is.null(project_root)) {
  stop("Could not find project root directory.\n",
       "Please run this script from the project root or scripts/ directory.\n",
       "Expected to find: config/config.yaml")
}
setwd(project_root)

message("\n")
message("================================================================")
message("  Step 1: Cell Type Annotation")
message("================================================================")
message(paste("  Started:", Sys.time()))
message(paste("  Working directory:", getwd()))
message("\n")

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# Source utility functions
source("R/utils/config.R")

# Load configuration
config <- load_config(config_path)
print_config_summary(config)

# Set seed for reproducibility
set.seed(config$reproducibility$seed)

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("\n========================================")
message("Loading Data")
message("========================================\n")

# Determine input file
input_path <- get_path(config, config$paths$spatial_object)

# Check if file with existing subtypes exists (for comparison)
alt_path <- get_path(config, config$paths$spatial_object_with_subtypes)
if (file.exists(alt_path)) {
  message("Found object with existing classifications - will compare methods")
  input_path <- alt_path
}

if (!file.exists(input_path)) {
  stop(paste("Input file not found:", input_path,
             "\nPlease check paths in config/config.yaml"))
}

message(paste("Loading:", input_path))
seurat_obj <- readRDS(input_path)

message(sprintf("  Cells: %d", ncol(seurat_obj)))
message(sprintf("  Genes: %d", nrow(seurat_obj)))

# Report existing metadata
if ("epithelial_subtype" %in% colnames(seurat_obj@meta.data)) {
  message("\nExisting classification (epithelial_subtype):")
  print(table(seurat_obj$epithelial_subtype))
}

if ("sample" %in% colnames(seurat_obj@meta.data)) {
  message("\nSamples:")
  print(table(seurat_obj$sample))
}

# ==============================================================================
# CALCULATE MODULE SCORES
# ==============================================================================

message("\n========================================")
message("Calculating Module Scores")
message("========================================\n")

# Get available genes
available_genes <- rownames(seurat_obj)
n_genes <- length(available_genes)
message(sprintf("Genes available in dataset: %d", n_genes))

# Calculate appropriate nbin for AddModuleScore
# Default is 24, but smaller datasets need fewer bins
# Rule of thumb: nbin should be < n_genes / 30
optimal_nbin <- min(24, max(5, floor(n_genes / 30)))
message(sprintf("Using nbin = %d for module scoring (based on %d genes)", optimal_nbin, n_genes))

# Prepare signatures from config
signatures <- list()

# Core epithelial subtype signatures
message("\nPreparing epithelial subtype signatures:")
for (sig_name in names(config$signatures$epithelial_subtypes)) {
  sig_genes <- config$signatures$epithelial_subtypes[[sig_name]]
  present_genes <- intersect(sig_genes, available_genes)

  message(sprintf("  %s: %d/%d genes present",
                  sig_name, length(present_genes), length(sig_genes)))

  if (length(present_genes) >= 3) {
    signatures[[sig_name]] <- present_genes
  } else {
    warning(sprintf("  Skipping %s: insufficient genes", sig_name))
  }
}

# Additional signatures for context
message("\nPreparing additional signatures:")
additional_sigs <- c("Whorl_Wnt", "Senescence", "SASP")
for (sig_name in additional_sigs) {
  if (sig_name %in% names(config$signatures)) {
    sig_genes <- config$signatures[[sig_name]]
    present_genes <- intersect(sig_genes, available_genes)

    message(sprintf("  %s: %d/%d genes present",
                    sig_name, length(present_genes), length(sig_genes)))

    if (length(present_genes) >= 3) {
      signatures[[sig_name]] <- present_genes
    }
  }
}

# Helper function for robust module scoring with fallback
score_module_robust <- function(seurat_obj, features, name, seed, nbin_start = 24) {
  # Try with decreasing nbin values until successful
  nbin_values <- c(nbin_start, 12, 8, 5)
  nbin_values <- nbin_values[nbin_values <= nbin_start]

  for (nbin in nbin_values) {
    result <- tryCatch({
      obj <- AddModuleScore(
        seurat_obj,
        features = list(features),
        name = name,
        nbin = nbin,
        seed = seed
      )
      list(success = TRUE, obj = obj, nbin = nbin)
    }, error = function(e) {
      list(success = FALSE, error = e$message)
    })

    if (result$success) {
      return(result)
    }
  }

  # If all nbin values fail, return failure
  return(list(success = FALSE, error = "All nbin values failed"))
}

# Calculate module scores using Seurat's AddModuleScore
message("\nCalculating module scores...")

for (sig_name in names(signatures)) {
  message(sprintf("  Scoring: %s", sig_name))

  result <- score_module_robust(
    seurat_obj,
    features = signatures[[sig_name]],
    name = paste0("score_", sig_name),
    seed = config$reproducibility$seed,
    nbin_start = optimal_nbin
  )

  if (result$success) {
    seurat_obj <- result$obj
    if (result$nbin != optimal_nbin) {
      message(sprintf("    (used nbin=%d)", result$nbin))
    }

    # Seurat appends "1" to the column name - rename for clarity
    old_name <- paste0("score_", sig_name, "1")
    new_name <- paste0("score_", sig_name)

    if (old_name %in% colnames(seurat_obj@meta.data)) {
      seurat_obj@meta.data[[new_name]] <- seurat_obj@meta.data[[old_name]]
      seurat_obj@meta.data[[old_name]] <- NULL
    }
  } else {
    warning(sprintf("    Failed to score %s: %s", sig_name, result$error))
  }
}

# Report score distributions
message("\nModule score summary statistics:")
score_cols <- grep("^score_", colnames(seurat_obj@meta.data), value = TRUE)

if (length(score_cols) == 0) {
  stop("No module scores were successfully calculated. Check gene availability and data structure.")
}

for (col in score_cols) {
  scores <- seurat_obj@meta.data[[col]]
  message(sprintf("  %s: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]",
                  gsub("^score_", "", col),
                  mean(scores, na.rm = TRUE),
                  sd(scores, na.rm = TRUE),
                  min(scores, na.rm = TRUE),
                  max(scores, na.rm = TRUE)))
}

# ==============================================================================
# CELL CYCLE SCORING (Multi-Method Consensus)
# ==============================================================================

message("\n========================================")
message("Cell Cycle Scoring")
message("========================================\n")

# Source the cell cycle scoring utility
source("R/utils/cell_cycle_scoring.R")

# Check if cell cycle already scored with consensus
if ("cc_consensus" %in% colnames(seurat_obj@meta.data)) {
  message("Cell cycle consensus already present in metadata")
  message("\nConsensus phase distribution:")
  print(table(seurat_obj$cc_consensus))
  if ("cc_confidence" %in% colnames(seurat_obj@meta.data)) {
    message(sprintf("Mean confidence: %.2f", mean(seurat_obj$cc_confidence)))
  }
} else if ("Phase" %in% colnames(seurat_obj@meta.data) &&
           !"cc_seurat" %in% colnames(seurat_obj@meta.data)) {
  # Legacy single-method scoring exists - keep it but note it
  message("Single-method cell cycle phase already present in metadata")
  print(table(seurat_obj$Phase))
  message("\nTo upgrade to consensus scoring, remove 'Phase' column and re-run")

} else {
  message("Scoring cell cycle phases using multiple methods...")

  # Get methods from config or use defaults
  cc_methods <- c("seurat", "module_score")
  if (!is.null(config$cell_cycle$methods)) {
    cc_methods <- config$cell_cycle$methods
  }
  message(sprintf("Methods: %s", paste(cc_methods, collapse = ", ")))

  # Run multi-method consensus scoring
  seurat_obj <- score_cell_cycle_consensus(
    seurat_obj,
    config = config,
    methods = cc_methods,
    seed = config$reproducibility$seed
  )
}

# ==============================================================================
# CLASSIFY CELLS BY MODULE SCORE
# ==============================================================================

message("\n========================================")
message("Classifying Cells")
message("========================================\n")

# Get classification parameters from config
min_threshold <- config$classification$min_score_threshold
if (is.null(min_threshold)) min_threshold <- 0.1
message(sprintf("Minimum score threshold: %.2f", min_threshold))

# Identify core epithelial signature scores (exclude auxiliary signatures)
core_signatures <- c("Basal_like", "Intermediate", "Specialized")
core_score_cols <- paste0("score_", core_signatures)
core_score_cols <- core_score_cols[core_score_cols %in% colnames(seurat_obj@meta.data)]

if (length(core_score_cols) == 0) {
  stop("No core epithelial signatures were scored. Cannot proceed with classification.")
}

message(sprintf("Core signatures for classification: %s",
                paste(gsub("score_", "", core_score_cols), collapse = ", ")))

# Extract score matrix for classification
score_matrix <- seurat_obj@meta.data[, core_score_cols, drop = FALSE]
colnames(score_matrix) <- gsub("^score_", "", colnames(score_matrix))

# Find maximum score and corresponding signature for each cell
max_scores <- apply(score_matrix, 1, max, na.rm = TRUE)
max_signature <- apply(score_matrix, 1, function(x) {
  if (all(is.na(x))) return(NA)
  names(x)[which.max(x)]
})

# Initial classification based on highest score
classification <- ifelse(max_scores >= min_threshold, max_signature, "Unassigned")
message(sprintf("\nInitial classification (threshold = %.2f):", min_threshold))
print(table(classification))

# ==============================================================================
# IDENTIFY TRANSIT-AMPLIFYING CELLS
# ==============================================================================

message("\n--- Identifying Transit-Amplifying Cells ---")

# Get Transit-Amplifying parameters
ta_config <- config$classification$transit_amplifying
prolif_threshold <- ta_config$proliferation_threshold
if (is.null(prolif_threshold)) prolif_threshold <- 0.15

override_threshold <- ta_config$override_threshold
if (is.null(override_threshold)) override_threshold <- 0.3

require_cycling <- ta_config$require_cycling
if (is.null(require_cycling)) require_cycling <- TRUE

message(sprintf("  Proliferation threshold: %.2f", prolif_threshold))
message(sprintf("  Override threshold: %.2f", override_threshold))
message(sprintf("  Require cycling (S/G2M): %s", require_cycling))

# Get proliferation score
if ("score_Transit_Amplifying" %in% colnames(seurat_obj@meta.data)) {
  prolif_scores <- seurat_obj@meta.data$score_Transit_Amplifying
} else {
  warning("  Transit_Amplifying score not found - creating from proliferation markers")
  # Fall back to using proliferation markers directly
  prolif_markers <- config$markers$proliferation
  prolif_present <- intersect(prolif_markers, available_genes)

  if (length(prolif_present) >= 2) {
    result <- score_module_robust(
      seurat_obj,
      features = prolif_present,
      name = "score_Proliferation",
      seed = config$reproducibility$seed,
      nbin_start = optimal_nbin
    )

    if (result$success) {
      seurat_obj <- result$obj
      seurat_obj@meta.data$score_Proliferation <- seurat_obj@meta.data$score_Proliferation1
      seurat_obj@meta.data$score_Proliferation1 <- NULL
      prolif_scores <- seurat_obj@meta.data$score_Proliferation
    } else {
      prolif_scores <- rep(0, ncol(seurat_obj))
      warning("  Failed to score proliferation markers")
    }
  } else {
    prolif_scores <- rep(0, ncol(seurat_obj))
    warning("  Insufficient proliferation markers - cannot identify TA cells")
  }
}

# Identify cycling cells
if (require_cycling && "Phase" %in% colnames(seurat_obj@meta.data)) {
  is_cycling <- seurat_obj$Phase %in% c("S", "G2M")
  message(sprintf("  Cycling cells (S/G2M): %d (%.1f%%)",
                  sum(is_cycling), mean(is_cycling) * 100))
} else {
  is_cycling <- rep(TRUE, ncol(seurat_obj))
  message("  Cell cycle not required or not available")
}

# Classify Transit-Amplifying cells
# Strategy 1: High proliferation overrides other classification
if (!is.null(ta_config$override_on_high_proliferation) && ta_config$override_on_high_proliferation) {
  high_prolif <- prolif_scores >= override_threshold & is_cycling
  n_override <- sum(high_prolif & classification != "Transit-Amplifying")
  classification[high_prolif] <- "Transit-Amplifying"
  message(sprintf("  High proliferation override: %d cells reclassified", n_override))
}

# Strategy 2: Unassigned cycling cells with moderate proliferation
moderate_prolif <- prolif_scores >= prolif_threshold & is_cycling
unassigned_ta <- classification == "Unassigned" & moderate_prolif
classification[unassigned_ta] <- "Transit-Amplifying"
message(sprintf("  Unassigned → Transit-Amplifying: %d cells", sum(unassigned_ta)))

# ==============================================================================
# DETECT HYBRID STATES
# ==============================================================================

message("\n--- Detecting Hybrid States ---")

# Hybrid detection parameters
hybrid_threshold <- min_threshold * 0.8
score_ratio_max <- 0.5  # Scores must be within 50% of each other

# Get individual scores
basal_scores <- seurat_obj@meta.data$score_Basal_like
int_scores <- seurat_obj@meta.data$score_Intermediate
spec_scores <- seurat_obj@meta.data$score_Specialized

# Basal-Intermediate hybrids
if (!is.null(basal_scores) && !is.null(int_scores)) {
  # Both scores above threshold
  both_high <- basal_scores >= hybrid_threshold & int_scores >= hybrid_threshold

  # Scores are close to each other (neither dominates)
  score_diff <- abs(basal_scores - int_scores)
  max_score <- pmax(basal_scores, int_scores)
  close_scores <- (score_diff / max_score) < score_ratio_max

  # Currently classified as either Basal or Intermediate
  eligible <- classification %in% c("Basal_like", "Intermediate")

  basal_int_hybrid <- both_high & close_scores & eligible
  classification[basal_int_hybrid] <- "Basal-Intermediate"
  message(sprintf("  Basal-Intermediate hybrids detected: %d", sum(basal_int_hybrid)))
}

# Intermediate-Specialized hybrids
if (!is.null(int_scores) && !is.null(spec_scores)) {
  both_high <- int_scores >= hybrid_threshold & spec_scores >= hybrid_threshold

  score_diff <- abs(int_scores - spec_scores)
  max_score <- pmax(int_scores, spec_scores)
  close_scores <- (score_diff / max_score) < score_ratio_max

  eligible <- classification %in% c("Intermediate", "Specialized")

  int_spec_hybrid <- both_high & close_scores & eligible
  classification[int_spec_hybrid] <- "Intermediate-Specialized"
  message(sprintf("  Intermediate-Specialized hybrids detected: %d", sum(int_spec_hybrid)))
}

# ==============================================================================
# FINALIZE CLASSIFICATION
# ==============================================================================

message("\n--- Finalizing Classification ---")

# Standardize names (replace underscores with hyphens for consistency)
classification <- gsub("Basal_like", "Basal-like", classification)

# Add to metadata
seurat_obj$module_score_subtype <- classification

# Create a factor with ordered levels
subtype_levels <- c(
  "Basal-like",
  "Transit-Amplifying",
  "Intermediate",
  "Specialized",
  "Basal-Intermediate",
  "Intermediate-Specialized",
  "Unassigned"
)
subtype_levels <- subtype_levels[subtype_levels %in% unique(classification)]
seurat_obj$module_score_subtype <- factor(
  seurat_obj$module_score_subtype,
  levels = subtype_levels
)

# Final classification summary
message("\nFinal classification:")
final_counts <- table(seurat_obj$module_score_subtype)
print(final_counts)

# Calculate percentages
final_pct <- prop.table(final_counts) * 100
message("\nPercentages:")
for (i in seq_along(final_counts)) {
  message(sprintf("  %s: %.1f%%", names(final_counts)[i], final_pct[i]))
}

# ==============================================================================
# COMPARE WITH ORIGINAL CLASSIFICATION
# ==============================================================================

if ("epithelial_subtype" %in% colnames(seurat_obj@meta.data)) {

  message("\n========================================")
  message("Comparison with Original Classification")
  message("========================================\n")

  original <- seurat_obj$epithelial_subtype
  new_class <- seurat_obj$module_score_subtype

  # Confusion matrix
  confusion <- table(Original = original, ModuleScore = new_class)

  # Calculate agreement
  # Need to handle different category names
  common_categories <- intersect(rownames(confusion), colnames(confusion))
  if (length(common_categories) > 0) {
    agreement_count <- sum(sapply(common_categories, function(cat) {
      if (cat %in% rownames(confusion) && cat %in% colnames(confusion)) {
        confusion[cat, cat]
      } else {
        0
      }
    }))
    total <- sum(confusion)
    agreement <- agreement_count / total * 100
    message(sprintf("Overall agreement: %.1f%%", agreement))
  }

  # Key changes
  message("\nKey reclassifications:")

  # Unknown → Transit-Amplifying
  if ("Unknown" %in% rownames(confusion) && "Transit-Amplifying" %in% colnames(confusion)) {
    n_unknown_to_ta <- confusion["Unknown", "Transit-Amplifying"]
    message(sprintf("  Unknown → Transit-Amplifying: %d cells", n_unknown_to_ta))
  }

  # Unknown → Other categories
  if ("Unknown" %in% rownames(confusion)) {
    unknown_row <- confusion["Unknown", ]
    for (cat in names(unknown_row)) {
      if (cat != "Unknown" && cat != "Unassigned" && unknown_row[cat] > 0) {
        message(sprintf("  Unknown → %s: %d cells", cat, unknown_row[cat]))
      }
    }
  }

  # Store comparison for later reference
  comparison_results <- list(
    confusion_matrix = confusion,
    original_counts = table(original),
    new_counts = table(new_class)
  )
  attr(seurat_obj, "classification_comparison") <- comparison_results
}

# ==============================================================================
# GENERATE DIAGNOSTIC PLOTS
# ==============================================================================

message("\n========================================")
message("Generating Diagnostic Plots")
message("========================================\n")

# Create output directory
fig_dir <- get_path(config, config$paths$figures_supp_dir)
ensure_dir(fig_dir)

# Get colors
colors <- get_colors(config, "epithelial_subtypes")

# Plot 1: Module score distributions
message("Creating score distribution plot...")
score_data <- seurat_obj@meta.data[, c(score_cols, "module_score_subtype")]
score_long <- tidyr::pivot_longer(
  score_data,
  cols = -module_score_subtype,
  names_to = "Signature",
  values_to = "Score"
)
score_long$Signature <- gsub("^score_", "", score_long$Signature)

p1 <- ggplot(score_long, aes(x = Score, fill = Signature)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = min_threshold, linetype = "dashed", color = "red") +
  facet_wrap(~Signature, scales = "free_y", ncol = 2) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Module Score Distributions",
       subtitle = sprintf("Red line = classification threshold (%.2f)", min_threshold),
       x = "Module Score", y = "Density")

# Plot 2: Scores by assigned subtype
message("Creating scores by subtype plot...")
p2 <- ggplot(score_long, aes(x = module_score_subtype, y = Score, fill = module_score_subtype)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~Signature, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = colors, na.value = "gray50") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Module Scores by Assigned Subtype",
       x = "", y = "Score")

# Plot 3: Classification comparison (if applicable)
if ("epithelial_subtype" %in% colnames(seurat_obj@meta.data)) {
  message("Creating classification comparison plot...")

  comparison_data <- seurat_obj@meta.data %>%
    count(epithelial_subtype, module_score_subtype) %>%
    group_by(epithelial_subtype) %>%
    mutate(pct = n / sum(n) * 100)

  p3 <- ggplot(comparison_data, aes(x = epithelial_subtype, y = pct, fill = module_score_subtype)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors, na.value = "gray50") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Reclassification: Original → Module Score",
         x = "Original Classification", y = "Percentage",
         fill = "New Classification")
} else {
  p3 <- NULL
}

# Plot 4: Cell counts by subtype
message("Creating cell counts plot...")
count_data <- as.data.frame(table(seurat_obj$module_score_subtype))
colnames(count_data) <- c("Subtype", "Count")
count_data$Percentage <- count_data$Count / sum(count_data$Count) * 100

p4 <- ggplot(count_data, aes(x = reorder(Subtype, -Count), y = Count, fill = Subtype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = colors, na.value = "gray50") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Cell Counts by Subtype",
       x = "", y = "Number of Cells")

# Plot 5: Proliferation vs subtype
message("Creating proliferation analysis plot...")
if ("score_Transit_Amplifying" %in% colnames(seurat_obj@meta.data) ||
    "score_Proliferation" %in% colnames(seurat_obj@meta.data)) {

  prolif_col <- ifelse("score_Transit_Amplifying" %in% colnames(seurat_obj@meta.data),
                       "score_Transit_Amplifying", "score_Proliferation")

  p5 <- ggplot(seurat_obj@meta.data,
               aes(x = module_score_subtype, y = .data[[prolif_col]],
                   fill = module_score_subtype)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
    geom_hline(yintercept = prolif_threshold, linetype = "dashed", color = "red") +
    scale_fill_manual(values = colors, na.value = "gray50") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Proliferation Score by Subtype",
         subtitle = sprintf("Red line = TA threshold (%.2f)", prolif_threshold),
         x = "", y = "Proliferation Score")
} else {
  p5 <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Proliferation score not available")
}

# Plot 6: Cell cycle by subtype
message("Creating cell cycle plot...")
if ("Phase" %in% colnames(seurat_obj@meta.data) &&
    !all(seurat_obj$Phase == "Unknown")) {

  cycle_data <- seurat_obj@meta.data %>%
    count(module_score_subtype, Phase) %>%
    group_by(module_score_subtype) %>%
    mutate(Proportion = n / sum(n))

  cycle_colors <- get_colors(config, "cell_cycle")

  p6 <- ggplot(cycle_data, aes(x = module_score_subtype, y = Proportion, fill = Phase)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cycle_colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Cell Cycle Distribution by Subtype",
         x = "", y = "Proportion")
} else {
  p6 <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Cell cycle not scored")
}

# Plot 7: Cell cycle consensus diagnostics (if multi-method scoring was used)
if ("cc_consensus" %in% colnames(seurat_obj@meta.data)) {
  message("Creating cell cycle consensus diagnostic plot...")
  cc_diagnostic_plot <- tryCatch({
    plot_cell_cycle_diagnostics(seurat_obj, config)
  }, error = function(e) {
    message(sprintf("  Could not generate diagnostic plot: %s", e$message))
    NULL
  })

  if (!is.null(cc_diagnostic_plot)) {
    cc_plot_path <- file.path(fig_dir, "01_cell_cycle_consensus_diagnostics")
    ggsave(paste0(cc_plot_path, ".pdf"), cc_diagnostic_plot, width = 12, height = 10)
    ggsave(paste0(cc_plot_path, ".png"), cc_diagnostic_plot, width = 12, height = 10, dpi = 300)
    message(sprintf("Saved: %s.pdf/.png", cc_plot_path))
  }
}


# Combine plots
message("Combining diagnostic plots...")
if (!is.null(p3)) {
  diagnostic_plot <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
    plot_annotation(
      title = "Cell Type Annotation Diagnostics",
      subtitle = sprintf("Module score classification with TA identification | n = %d cells",
                         ncol(seurat_obj)),
      tag_levels = "A"
    )
} else {
  diagnostic_plot <- (p1 | p2) / (p4 | p5) / (p6 | plot_spacer()) +
    plot_annotation(
      title = "Cell Type Annotation Diagnostics",
      subtitle = sprintf("Module score classification with TA identification | n = %d cells",
                         ncol(seurat_obj)),
      tag_levels = "A"
    )
}

# Save diagnostic plot
diagnostic_path <- file.path(fig_dir, "01_cell_type_annotation_diagnostics")
ggsave(paste0(diagnostic_path, ".pdf"), diagnostic_plot, width = 14, height = 16)
ggsave(paste0(diagnostic_path, ".png"), diagnostic_plot, width = 14, height = 16, dpi = 300)
message(sprintf("Saved: %s.pdf/.png", diagnostic_path))

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

message("\n========================================")
message("Saving Outputs")
message("========================================\n")

# Create output directories
objects_dir <- get_path(config, config$paths$objects_dir)
tables_dir <- get_path(config, config$paths$tables_dir)
ensure_dir(objects_dir)
ensure_dir(tables_dir)

# Save annotated Seurat object
output_path <- file.path(objects_dir, "01_seurat_annotated.rds")
saveRDS(seurat_obj, output_path)
message(sprintf("Saved annotated object: %s", output_path))

# Save classification summary table
summary_table <- seurat_obj@meta.data %>%
  group_by(module_score_subtype) %>%
  summarise(
    n_cells = n(),
    pct_total = n() / nrow(seurat_obj@meta.data) * 100,
    mean_prolif_score = if ("score_Transit_Amplifying" %in% names(.)) {
      mean(score_Transit_Amplifying, na.rm = TRUE)
    } else if ("score_Proliferation" %in% names(.)) {
      mean(score_Proliferation, na.rm = TRUE)
    } else NA,
    pct_cycling = if ("Phase" %in% names(.)) {
      mean(Phase %in% c("S", "G2M"), na.rm = TRUE) * 100
    } else NA,
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

write.csv(summary_table,
          file.path(tables_dir, "01_classification_summary.csv"),
          row.names = FALSE)
message(sprintf("Saved classification summary: %s",
                file.path(tables_dir, "01_classification_summary.csv")))

# Save per-sample breakdown if applicable
if ("sample" %in% colnames(seurat_obj@meta.data)) {
  sample_summary <- seurat_obj@meta.data %>%
    count(sample, module_score_subtype) %>%
    group_by(sample) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()

  write.csv(sample_summary,
            file.path(tables_dir, "01_classification_by_sample.csv"),
            row.names = FALSE)
  message(sprintf("Saved per-sample summary: %s",
                  file.path(tables_dir, "01_classification_by_sample.csv")))
}

# Save gene signatures used
signatures_df <- do.call(rbind, lapply(names(signatures), function(sig_name) {
  data.frame(
    signature = sig_name,
    gene = signatures[[sig_name]],
    stringsAsFactors = FALSE
  )
}))
write.csv(signatures_df,
          file.path(tables_dir, "01_gene_signatures_used.csv"),
          row.names = FALSE)
message(sprintf("Saved gene signatures: %s",
                file.path(tables_dir, "01_gene_signatures_used.csv")))

# Evaluate cell cycle scoring statistics
if ("cc_consensus" %in% colnames(seurat_obj@meta.data)) {
  message("\nEvaluating cell cycle scoring statistics...")

  cc_stats <- tryCatch({
    evaluate_cell_cycle_statistics(
      seurat_obj,
      config = config,
      proliferation_markers = c("MKI67", "TOP2A", "PCNA", "STMN1", "CDK1", "CCNB1")
    )
  }, error = function(e) {
    message(sprintf("  Could not evaluate statistics: %s", e$message))
    NULL
  })

  if (!is.null(cc_stats)) {
    # Determine table directory (handle both table_dir and tables_dir variable names)
    cc_table_dir <- if (exists("table_dir")) {
      table_dir
    } else if (exists("tables_dir")) {
      tables_dir
    } else {
      get_path(config, config$paths$tables_dir)
    }

    # Ensure directory exists
    if (!dir.exists(cc_table_dir)) {
      dir.create(cc_table_dir, recursive = TRUE)
    }

    # Save statistics as RDS
    stats_path <- file.path(cc_table_dir, "01_cell_cycle_statistics.rds")
    saveRDS(cc_stats, stats_path)
    message(sprintf("Saved cell cycle statistics: %s", stats_path))

    # Generate and save markdown report
    report_path <- file.path(cc_table_dir, "01_cell_cycle_evaluation_report.md")
    tryCatch({
      generate_cell_cycle_report(cc_stats, report_path)
    }, error = function(e) {
      message(sprintf("  Could not generate report: %s", e$message))
    })

    # Save summary table as CSV
    if (!is.null(cc_stats$marker_validation)) {
      marker_path <- file.path(cc_table_dir, "01_cell_cycle_marker_validation.csv")
      write.csv(cc_stats$marker_validation, marker_path, row.names = FALSE)
      message(sprintf("Saved marker validation: %s", marker_path))
    }

    # Save pairwise agreement as CSV
    if (!is.null(cc_stats$pairwise_kappa)) {
      kappa_df <- do.call(rbind, lapply(names(cc_stats$pairwise_kappa), function(pair) {
        data.frame(
          comparison = pair,
          cohens_kappa = cc_stats$pairwise_kappa[[pair]]$kappa,
          percent_agreement = cc_stats$pairwise_kappa[[pair]]$percent_agreement
        )
      }))
      kappa_path <- file.path(cc_table_dir, "01_cell_cycle_method_agreement.csv")
      write.csv(kappa_df, kappa_path, row.names = FALSE)
      message(sprintf("Saved method agreement: %s", kappa_path))
    }
  }
}

# ==============================================================================
# COMPLETION
# ==============================================================================

message("\n")
message("================================================================")
message("  Cell Type Annotation Complete!")
message("================================================================")
message(sprintf("  Finished: %s", Sys.time()))
message(sprintf("  Total cells: %d", ncol(seurat_obj)))
message(sprintf("  Subtypes identified: %d", length(unique(seurat_obj$module_score_subtype))))
message("\nOutputs:")
message(sprintf("  Object: %s", output_path))
message(sprintf("  Tables: %s", tables_dir))
message(sprintf("  Figures: %s", fig_dir))
message("\nNext step: Run 02_trajectory_analysis.R")
message("================================================================\n")
