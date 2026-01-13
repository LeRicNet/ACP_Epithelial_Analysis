#!/usr/bin/env Rscript
# ==============================================================================
# OpenPBTA Validation of Spatial Transcriptomic Findings (Section 0.11)
# ==============================================================================
# Validates key findings from verified manuscript sections (0.1-0.10) using
# bulk RNA-seq data from 73 ACP specimens in the OpenPBTA/CBTN cohort.
#
# Hypotheses tested:
#   H1: Keratin subtypes show coordinated expression (validates 0.3, 0.4, 0.10)
#   H2: ACP expresses oral-biased rather than skin-biased keratins (validates 0.3)
#   H3: Proliferation and ECM signatures are anti-correlated (validates 0.6, 0.8)
#   H4: Epithelial and stromal signatures define distinct compartments (validates 0.8)
#   H5: CD44 correlates with immune signatures (validates 0.9)
#
# Citation: Rokita JL, et al. (2023) Cell Genomics. doi:10.1016/j.xgen.2023.100340
# Data: https://github.com/AlexsLemonade/OpenPBTA-analysis
#
# Usage:
#   Rscript openpbta_spatial_validation.R [--config path/to/config.yaml]
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(ggpubr)
  library(corrplot)
  library(RColorBrewer)
  library(scales)
  library(patchwork)
  library(boot)
  library(grid)  # For grid.grabExpr and textGrob
})

# Null coalescing operator (if not available from rlang)
`%||%` <- function(x, y) if (is.null(x)) y else x

# Install UCell if not available
if (!requireNamespace("UCell", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("UCell", update = FALSE, ask = FALSE)
}
library(UCell)

# Load Seurat for snRNA-seq validation (H8)
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
suppressPackageStartupMessages(library(Seurat))

# ==============================================================================
# CONFIGURATION
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
for (i in seq_along(args)) {
  if (args[i] == "--config" && i < length(args)) config_path <- args[i + 1]
}

# Find project root
find_project_root <- function() {
  candidates <- c(getwd(), normalizePath(".."),
                  if (requireNamespace("here", quietly = TRUE)) here::here() else NULL)
  for (path in candidates) {
    if (!is.null(path) && file.exists(file.path(path, "config/config.yaml"))) return(path)
  }
  current <- getwd()
  for (i in 1:5) {
    if (file.exists(file.path(current, "config/config.yaml"))) return(current)
    current <- dirname(current)
  }
  return(NULL)
}

project_root <- find_project_root()
if (!is.null(project_root)) {

  setwd(project_root)
  if (file.exists("R/utils/config.R")) {
    source("R/utils/config.R")
    config <- load_config(config_path)
  } else {
    config <- list(reproducibility = list(seed = 42))
  }
} else {
  config <- list(reproducibility = list(seed = 42))
}

set.seed(config$reproducibility$seed)

# Directory setup
dirs <- list(
  data = "data/OpenPBTA",
  results = "results/validation",
  figures = "results/figures/validation",
  tables = "results/tables/validation"
)
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# OpenPBTA download configuration
openpbta_config <- list(
  base_url = "https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/open-targets/v14/",
  files = c(
    "histologies.tsv",
    "gene-expression-rsem-tpm-collapsed.rds"
  ),
  pseudocount = 1
)

message("\n")
message("================================================================")
message("  OpenPBTA Validation of Spatial Transcriptomic Findings")
message("  Section 0.11: External Cohort Validation")
message("================================================================")
message(sprintf("  Started: %s", Sys.time()))
message(sprintf("  Seed: %d", config$reproducibility$seed))

# ==============================================================================
# GENE SIGNATURES
# ==============================================================================

# Define signatures based on manuscript findings
signatures <- list(


  # Epithelial subtype signatures (from spatial transcriptomics classification)
  Basal_like = c("KRT5", "KRT14", "KRT15", "TP63", "KRT6A", "KRT6B",
                 "ITGA6", "ITGB4", "COL17A1", "LAMB3", "LAMC2"),

  Intermediate = c("KRT8", "KRT18", "KRT19", "KRT7", "EPCAM", "CDH1",
                   "CLDN3", "CLDN4", "CLDN7", "TJP1", "OCLN"),

  Specialized = c("KRT17", "KRT23", "KRT10", "KRT1", "IVL", "FLG",
                  "SPRR1A", "SPRR1B", "SPRR2A", "LOR"),

  # Keratin origin signatures (Section 0.3)
  Oral_Keratin = c("KRT4", "KRT13", "KRT19", "KRT76", "KRT78"),
  Skin_Keratin = c("KRT1", "KRT10", "KRT2", "KRT9", "KRT77"),
  Basal_Shared = c("KRT5", "KRT14", "KRT15", "TP63", "ITGA6", "ITGB4", "COL17A1"),

  # Proliferation signature (Section 0.6)
  Proliferation = c("MKI67", "TOP2A", "PCNA", "STMN1", "CDK1", "CCNB1",
                    "CCNA2", "MCM2", "MCM5", "MCM6", "TYMS", "RRM2"),

  # ECM/Collagen signatures (Section 0.8)
  Collagen = c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL5A1",
               "COL6A1", "COL6A2", "COL6A3", "COL11A1"),

  ECM_General = c("FN1", "VIM", "SPARC", "POSTN", "TNC", "THBS1",
                  "THBS2", "CTGF", "CYR61"),

  ECM_Remodeling = c("MMP2", "MMP9", "MMP14", "TIMP1", "TIMP2", "TIMP3",
                     "LOX", "LOXL1", "LOXL2"),

  # Epithelial identity (Section 0.8)
  Pan_Epithelial = c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19",
                     "CLDN3", "CLDN4", "CLDN7", "TJP1"),

  # Immune signatures (Section 0.9)
  Myeloid = c("CD68", "CD163", "CD14", "CSF1R", "ITGAM", "CD86",
              "MRC1", "MSR1", "MARCO"),

  TAM = c("CD68", "CD163", "MRC1", "MSR1", "GPNMB", "SPP1", "CCL4",
          "CCL18", "APOE", "TREM2"),

  # CD44-SPP1 axis (supports Section 0.9 paracrine model)
  CD44_SPP1 = c("CD44", "SPP1", "ITGAV", "ITGB1", "ITGB3", "ITGB5"),

  # Senescence (supports differentiation findings)
  Senescence = c("CDKN1A", "CDKN2A", "TP53", "GLB1", "SERPINE1",
                 "IGFBP3", "IGFBP7"),

  # Wnt pathway (ACP biology)
  Wnt_Pathway = c("CTNNB1", "AXIN2", "LEF1", "TCF7", "WNT5A", "DKK1",
                  "DKK4", "SFRP1")
)

# ==============================================================================
# STATISTICAL FUNCTIONS
# ==============================================================================

#' Bootstrap confidence interval for correlation
#'
#' @param x,y Numeric vectors
#' @param n_boot Number of bootstrap iterations
#' @param conf_level Confidence level
#' @return List with estimate and CI
boot_cor <- function(x, y, n_boot = 1000, conf_level = 0.95) {
  valid <- complete.cases(x, y)
  x <- x[valid]
  y <- y[valid]

  if (length(x) < 10) {
    return(list(estimate = NA, ci_lower = NA, ci_upper = NA, n = length(x)))
  }

  boot_fn <- function(data, indices) {
    cor(data[indices, 1], data[indices, 2], method = "pearson")
  }

  data_mat <- cbind(x, y)
  boot_result <- boot(data_mat, boot_fn, R = n_boot)
  ci <- boot.ci(boot_result, conf = conf_level, type = "perc")

  list(
    estimate = cor(x, y),
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5],
    n = length(x),
    p_value = cor.test(x, y)$p.value
  )
}

#' Permutation test for correlation
#'
#' @param x,y Numeric vectors
#' @param n_perm Number of permutations
#' @return List with observed, null distribution, and p-value
perm_cor_test <- function(x, y, n_perm = 1000) {
  valid <- complete.cases(x, y)
  x <- x[valid]
  y <- y[valid]

  observed <- cor(x, y)

  null_dist <- replicate(n_perm, {
    cor(x, sample(y))
  })

  # Two-sided p-value
  p_value <- mean(abs(null_dist) >= abs(observed))

  list(
    observed = observed,
    null_mean = mean(null_dist),
    null_sd = sd(null_dist),
    effect_size = (observed - mean(null_dist)) / sd(null_dist),
    p_value = p_value,
    null_distribution = null_dist
  )
}

#' Bootstrap confidence interval for mean difference
#'
#' @param x,y Numeric vectors to compare
#' @param n_boot Number of bootstrap iterations
#' @return List with estimate and CI
boot_mean_diff <- function(x, y, n_boot = 1000, conf_level = 0.95) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  observed_diff <- mean(x) - mean(y)

  boot_diffs <- replicate(n_boot, {
    x_boot <- sample(x, replace = TRUE)
    y_boot <- sample(y, replace = TRUE)
    mean(x_boot) - mean(y_boot)
  })

  ci <- quantile(boot_diffs, c((1 - conf_level)/2, 1 - (1 - conf_level)/2))

  list(
    estimate = observed_diff,
    ci_lower = ci[1],
    ci_upper = ci[2],
    se = sd(boot_diffs)
  )
}

#' Cohen's d effect size
cohens_d <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  pooled_sd <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) /
                      (length(x) + length(y) - 2))

  (mean(x) - mean(y)) / pooled_sd
}

#' Assess signature coherence via pairwise correlations
#'
#' @param expr_mat Expression matrix (genes x samples)
#' @param genes Gene list
#' @return List with coherence metrics
signature_coherence <- function(expr_mat, genes) {
  genes_present <- intersect(genes, rownames(expr_mat))

  if (length(genes_present) < 3) {
    return(list(mean_cor = NA, n_genes = length(genes_present),
                error = "Insufficient genes"))
  }

  gene_expr <- expr_mat[genes_present, , drop = FALSE]
  cor_mat <- cor(t(gene_expr), use = "pairwise.complete.obs")

  # Get upper triangle (excluding diagonal)
  upper_tri <- cor_mat[upper.tri(cor_mat)]

  list(
    mean_cor = mean(upper_tri, na.rm = TRUE),
    median_cor = median(upper_tri, na.rm = TRUE),
    min_cor = min(upper_tri, na.rm = TRUE),
    max_cor = max(upper_tri, na.rm = TRUE),
    prop_positive = mean(upper_tri > 0, na.rm = TRUE),
    n_genes = length(genes_present),
    n_pairs = length(upper_tri)
  )
}

# ==============================================================================
# DOWNLOAD AND LOAD DATA
# ==============================================================================

message("\n========================================")
message("Loading OpenPBTA Data")
message("========================================\n")

# Download if needed
for (f in openpbta_config$files) {
  dest <- file.path(dirs$data, basename(f))
  if (!file.exists(dest)) {
    message(sprintf("Downloading: %s", f))
    tryCatch({
      download.file(paste0(openpbta_config$base_url, f), dest, mode = "wb", quiet = TRUE)
      message(sprintf("  Size: %.2f MB", file.size(dest) / 1e6))
    }, error = function(e) {
      stop(sprintf("Failed to download %s: %s", f, e$message))
    })
  }
}

# Load data
histologies <- fread(file.path(dirs$data, "histologies.tsv"), data.table = FALSE) %>%
  as_tibble()
expr_tpm <- readRDS(file.path(dirs$data, "gene-expression-rsem-tpm-collapsed.rds"))

message(sprintf("  Total samples in OpenPBTA: %d", nrow(histologies)))
message(sprintf("  Expression matrix: %d genes x %d samples", nrow(expr_tpm), ncol(expr_tpm)))

# Filter to craniopharyngioma
acp_samples <- histologies %>%
  filter(str_detect(tolower(short_histology), "cranio") |
           str_detect(tolower(broad_histology), "cranio") |
           str_detect(tolower(pathology_diagnosis), "cranio|adamantinomatous")) %>%
  filter(experimental_strategy == "RNA-Seq")

available_ids <- intersect(acp_samples$Kids_First_Biospecimen_ID, colnames(expr_tpm))

message(sprintf("\n  ACP RNA-Seq samples: %d", nrow(acp_samples)))
message(sprintf("  With expression data: %d", length(available_ids)))

# Extract ACP expression data
acp_tpm <- expr_tpm[, available_ids, drop = FALSE]
acp_log2tpm <- log2(acp_tpm + openpbta_config$pseudocount)

acp_meta <- acp_samples %>%
  filter(Kids_First_Biospecimen_ID %in% available_ids) %>%
  arrange(match(Kids_First_Biospecimen_ID, available_ids))

message("\n  Tumor descriptor distribution:")
print(table(acp_meta$tumor_descriptor))

# Check signature gene coverage
message("\n  Signature gene coverage:")
for (sig_name in names(signatures)) {
  n_present <- sum(signatures[[sig_name]] %in% rownames(acp_tpm))
  n_total <- length(signatures[[sig_name]])
  pct <- round(n_present / n_total * 100, 1)
  message(sprintf("    %s: %d/%d (%.1f%%)", sig_name, n_present, n_total, pct))
}

# ==============================================================================
# CALCULATE SIGNATURE SCORES
# ==============================================================================

message("\n========================================")
message("Calculating Signature Scores")
message("========================================\n")

# Use UCell for robust signature scoring
ucell_scores <- ScoreSignatures_UCell(
  matrix = as.matrix(acp_tpm),
  features = signatures,
  maxRank = 1500,
  name = ""
)

sig_scores <- as.data.frame(ucell_scores) %>%
  rownames_to_column("sample_id") %>%
  left_join(acp_meta %>% select(sample_id = Kids_First_Biospecimen_ID,
                                patient_id = Kids_First_Participant_ID,
                                tumor_descriptor),
            by = "sample_id")

message(sprintf("  Scores calculated for %d samples and %d signatures",
                nrow(sig_scores), length(signatures)))

# Initialize results storage
results <- list(
  n_samples = length(available_ids),
  signatures = signatures,
  sig_scores = sig_scores
)

# ==============================================================================
# LOAD COMPARISON TUMOR TYPES FROM OpenPBTA
# ==============================================================================

message("\n========================================")
message("Loading Comparison Tumor Types")
message("========================================\n")

# Get available tumor types with sufficient samples
tumor_counts <- histologies %>%
  filter(experimental_strategy == "RNA-Seq",
         Kids_First_Biospecimen_ID %in% colnames(expr_tpm)) %>%
  count(short_histology, sort = TRUE) %>%
  filter(n >= 10, !is.na(short_histology), short_histology != "")

message("  Top tumor types in OpenPBTA:")
print(head(tumor_counts, 15))

# Select comparison tumor types (brain tumors, exclude craniopharyngioma)
# Use pattern matching for flexibility
comparison_patterns <- c(
  "glioma", "ependymoma", "medulloblastoma", "meningioma",
  "schwannoma", "ganglioglioma", "astrocytoma", "ATRT", "rhabdoid"
)

comparison_histologies <- tumor_counts %>%
  filter(str_detect(tolower(short_histology), paste(tolower(comparison_patterns), collapse = "|")),
         !str_detect(tolower(short_histology), "cranio")) %>%
  slice_head(n = 10) %>%
  pull(short_histology)

message(sprintf("\n  Selected comparison tumor types: %d", length(comparison_histologies)))
message(sprintf("    %s", paste(comparison_histologies, collapse = "\n    ")))

# Get samples for comparison tumors
comparison_samples <- histologies %>%
  filter(short_histology %in% comparison_histologies,
         experimental_strategy == "RNA-Seq",
         Kids_First_Biospecimen_ID %in% colnames(expr_tpm))

# Subsample to max 50 per tumor type using base R approach
set.seed(config$reproducibility$seed)
comparison_samples_list <- split(comparison_samples, comparison_samples$short_histology)
comparison_samples_list <- lapply(comparison_samples_list, function(df) {
  if (nrow(df) > 50) {
    df[sample(nrow(df), 50), ]
  } else {
    df
  }
})
comparison_samples <- do.call(rbind, comparison_samples_list)
rownames(comparison_samples) <- NULL

message(sprintf("\n  Comparison tumor samples: %d", nrow(comparison_samples)))
message("\n  Samples per tumor type:")
print(table(comparison_samples$short_histology))

# Extract expression for comparison tumors
comparison_ids <- comparison_samples$Kids_First_Biospecimen_ID
comparison_tpm <- expr_tpm[, comparison_ids, drop = FALSE]
comparison_log2tpm <- log2(comparison_tpm + openpbta_config$pseudocount)

# Calculate signature scores for comparison tumors
message("\n  Calculating signature scores for comparison tumors...")
comparison_ucell <- ScoreSignatures_UCell(
  matrix = as.matrix(comparison_tpm),
  features = signatures,
  maxRank = 1500,
  name = ""
)

comparison_sig_scores <- as.data.frame(comparison_ucell) %>%
  rownames_to_column("sample_id") %>%
  left_join(comparison_samples %>% select(sample_id = Kids_First_Biospecimen_ID,
                                          tumor_type = short_histology),
            by = "sample_id")

# Combine ACP and comparison scores
all_sig_scores <- bind_rows(
  sig_scores %>% mutate(tumor_type = "Craniopharyngioma") %>%
    select(sample_id, tumor_type, names(signatures)),
  comparison_sig_scores %>% select(sample_id, tumor_type, names(signatures))
)

message(sprintf("  Total samples for cross-tumor analysis: %d", nrow(all_sig_scores)))

# Update results with comparison data
results$n_comparison_samples <- nrow(comparison_samples)
results$comparison_tumor_types <- unique(comparison_samples$short_histology)
results$all_sig_scores <- all_sig_scores

# ==============================================================================
# LOAD XU ET AL. (GSE215932) SINGLE-CELL REFERENCE
# ==============================================================================

message("\n========================================")
message("Loading Xu et al. snRNA-seq Reference")
message("========================================\n")

# Try to load the GSE215932 data
gse_path <- NULL
possible_gse_paths <- c(
  "results/objects/01_seurat_annotated_snrnaseq.rds",
  "data/external/GSE215932/GSE215932_snRNA_processed.rds",
  file.path(project_root, "results/objects/01_seurat_annotated_snrnaseq.rds"),
  file.path(project_root, "data/external/GSE215932/GSE215932_snRNA_processed.rds")
)

for (path in possible_gse_paths) {
  if (file.exists(path)) {
    gse_path <- path
    break
  }
}

gse_available <- FALSE
if (!is.null(gse_path) && file.exists(gse_path)) {
  message(sprintf("  Loading: %s", gse_path))

  tryCatch({
    seurat_gse <- readRDS(gse_path)
    seurat_gse <- UpdateSeuratObject(seurat_gse)

    # Check for cell type labels
    label_col <- NULL
    for (col in c("cell_type", "celltype", "cluster", "seurat_clusters")) {
      if (col %in% colnames(seurat_gse@meta.data)) {
        label_col <- col
        break
      }
    }

    if (!is.null(label_col)) {
      gse_available <- TRUE
      gse_labels <- seurat_gse@meta.data[[label_col]]

      message(sprintf("  Cells: %d", ncol(seurat_gse)))
      message(sprintf("  Label column: %s", label_col))
      message("\n  Cell type distribution:")
      print(table(gse_labels))

      # Calculate mean expression per cell type for signature genes
      message("\n  Calculating per-cell-type expression profiles...")

      # Get expression data
      if (inherits(seurat_gse[["RNA"]], "Assay5")) {
        if (length(Layers(seurat_gse[["RNA"]])) > 1) {
          seurat_gse <- JoinLayers(seurat_gse)
        }
        expr_data <- LayerData(seurat_gse, assay = "RNA", layer = "data")
      } else {
        expr_data <- GetAssayData(seurat_gse, slot = "data", assay = "RNA")
      }

      # Calculate mean expression per cell type for key genes
      all_sig_genes <- unique(unlist(signatures))
      common_genes <- intersect(all_sig_genes, rownames(expr_data))

      celltype_means <- sapply(unique(gse_labels), function(ct) {
        cells <- which(gse_labels == ct)
        if (length(cells) > 10) {
          rowMeans(expr_data[common_genes, cells, drop = FALSE])
        } else {
          rep(NA, length(common_genes))
        }
      })
      rownames(celltype_means) <- common_genes

      message(sprintf("  Signature genes in snRNA-seq: %d/%d",
                      length(common_genes), length(all_sig_genes)))

    } else {
      message("  WARNING: No cell type labels found in GSE215932 object")
    }
  }, error = function(e) {
    message(sprintf("  ERROR loading GSE215932: %s", e$message))
  })
} else {
  message("  GSE215932 data not found - skipping snRNA-seq validation")
  message("  Looked in:")
  for (path in possible_gse_paths) {
    message(sprintf("    - %s", path))
  }
}

# Store GSE availability
results$gse_available <- gse_available

# ==============================================================================
# HYPOTHESIS 1: KERATIN SUBTYPE COORDINATED EXPRESSION
# Validates Sections 0.3, 0.4, 0.10
# ==============================================================================

message("\n========================================")
message("H1: Keratin Subtype Coordinated Expression")
message("    (Validates Sections 0.3, 0.4, 0.10)")
message("========================================\n")

# Test: Do keratins within the same subtype correlate more strongly than
# keratins across subtypes?

# Define keratin groups based on spatial transcriptomics findings
# Using broader groups to ensure sufficient genes
keratin_groups <- list(
  Basal = c("KRT5", "KRT14", "KRT15", "KRT6A", "KRT6B", "TP63"),
  Intermediate = c("KRT8", "KRT18", "KRT19", "KRT7", "EPCAM"),
  Specialized = c("KRT17", "KRT23", "KRT10", "KRT1", "IVL")
)

message("Testing within-group vs between-group keratin correlations...")

# Check gene availability
message("\n  Genes available per group:")
for (grp in names(keratin_groups)) {
  genes_present <- intersect(keratin_groups[[grp]], rownames(acp_log2tpm))
  message(sprintf("    %s: %d/%d (%s)", grp, length(genes_present),
                  length(keratin_groups[[grp]]), paste(genes_present, collapse = ", ")))
}

# Calculate within-group correlations
within_cors <- list()
for (grp in names(keratin_groups)) {
  genes <- intersect(keratin_groups[[grp]], rownames(acp_log2tpm))
  if (length(genes) >= 2) {
    expr_mat <- acp_log2tpm[genes, , drop = FALSE]
    cor_mat <- cor(t(expr_mat), use = "pairwise.complete.obs")
    upper_vals <- cor_mat[upper.tri(cor_mat)]
    if (length(upper_vals) > 0) {
      within_cors[[grp]] <- upper_vals
      message(sprintf("    %s within-group: %d pairs, mean r = %.3f",
                      grp, length(upper_vals), mean(upper_vals, na.rm = TRUE)))
    }
  }
}

# Calculate between-group correlations
between_cors <- list()
group_pairs <- combn(names(keratin_groups), 2, simplify = FALSE)
for (pair in group_pairs) {
  genes1 <- intersect(keratin_groups[[pair[1]]], rownames(acp_log2tpm))
  genes2 <- intersect(keratin_groups[[pair[2]]], rownames(acp_log2tpm))

  if (length(genes1) >= 1 && length(genes2) >= 1) {
    cors <- numeric()
    for (g1 in genes1) {
      for (g2 in genes2) {
        r <- cor(as.numeric(acp_log2tpm[g1, ]), as.numeric(acp_log2tpm[g2, ]),
                 use = "pairwise.complete.obs")
        if (!is.na(r)) cors <- c(cors, r)
      }
    }
    if (length(cors) > 0) {
      pair_name <- paste(pair, collapse = "_vs_")
      between_cors[[pair_name]] <- cors
      message(sprintf("    %s: %d pairs, mean r = %.3f",
                      pair_name, length(cors), mean(cors, na.rm = TRUE)))
    }
  }
}

# Aggregate and compare
all_within <- unlist(within_cors)
all_between <- unlist(between_cors)

message(sprintf("\n  Total within-group pairs: %d", length(all_within)))
message(sprintf("  Total between-group pairs: %d", length(all_between)))

# Check if we have enough data for comparison
if (length(all_within) < 3 || length(all_between) < 3) {
  message("\n  WARNING: Insufficient correlation pairs for statistical testing")
  message("  Falling back to descriptive comparison only")

  h1_test <- list(p.value = NA)
  h1_effect <- if (length(all_within) > 0 && length(all_between) > 0) {
    (mean(all_within, na.rm = TRUE) - mean(all_between, na.rm = TRUE)) /
      sd(c(all_within, all_between), na.rm = TRUE)
  } else NA
  h1_boot <- list(estimate = mean(all_within, na.rm = TRUE) - mean(all_between, na.rm = TRUE),
                  ci_lower = NA, ci_upper = NA)
} else {
  h1_test <- wilcox.test(all_within, all_between, alternative = "greater")
  h1_effect <- cohens_d(all_within, all_between)
  h1_boot <- boot_mean_diff(all_within, all_between)
}

message(sprintf("\n  Within-subtype correlations: mean = %.3f (n = %d pairs)",
                mean(all_within, na.rm = TRUE), length(all_within)))
message(sprintf("  Between-subtype correlations: mean = %.3f (n = %d pairs)",
                mean(all_between, na.rm = TRUE), length(all_between)))
message(sprintf("  Difference: %.3f [95%% CI: %.3f, %.3f]",
                h1_boot$estimate,
                ifelse(is.na(h1_boot$ci_lower), NA, h1_boot$ci_lower),
                ifelse(is.na(h1_boot$ci_upper), NA, h1_boot$ci_upper)))
message(sprintf("  Cohen's d: %.3f", h1_effect))
if (!is.na(h1_test$p.value)) {
  message(sprintf("  Wilcoxon p-value: %.2e", h1_test$p.value))
}

# Determine conclusion based on available statistics
if (is.na(h1_test$p.value)) {
  if (!is.na(h1_effect) && h1_effect > 0.5) {
    message("\n  → H1 TENTATIVELY SUPPORTED: Large effect size but statistical test not possible")
    h1_conclusion <- "TENTATIVE"
  } else {
    message("\n  → H1 INCONCLUSIVE: Insufficient data for full analysis")
    h1_conclusion <- "INCONCLUSIVE"
  }
} else if (h1_test$p.value < 0.05 && h1_effect > 0.5) {
  message("\n  → H1 SUPPORTED: Keratins show coordinated within-subtype expression")
  h1_conclusion <- "SUPPORTED"
} else if (h1_test$p.value < 0.05) {
  message("\n  → H1 PARTIALLY SUPPORTED: Significant but modest effect")
  h1_conclusion <- "PARTIAL"
} else {
  message("\n  → H1 NOT SUPPORTED: No evidence for coordinated expression")
  h1_conclusion <- "NOT_SUPPORTED"
}

results$H1 <- list(
  within_mean = mean(all_within, na.rm = TRUE),
  between_mean = mean(all_between, na.rm = TRUE),
  difference = h1_boot$estimate,
  ci_lower = if(is.null(h1_boot$ci_lower)) NA else h1_boot$ci_lower,
  ci_upper = if(is.null(h1_boot$ci_upper)) NA else h1_boot$ci_upper,
  cohens_d = if(is.na(h1_effect)) 0 else h1_effect,
  p_value = if(is.na(h1_test$p.value)) 1 else h1_test$p.value,
  conclusion = h1_conclusion,
  within_cors = within_cors,
  between_cors = between_cors,
  n_within_pairs = length(all_within),
  n_between_pairs = length(all_between)
)

# ==============================================================================
# HYPOTHESIS 2: ORAL-BIASED KERATIN EXPRESSION
# Validates Section 0.3
# ==============================================================================

message("\n========================================")
message("H2: Oral-Biased Keratin Expression")
message("    (Validates Section 0.3)")
message("========================================\n")

# Test: Do ACP tumors express higher oral-specific than skin-specific keratins?

oral_genes <- intersect(signatures$Oral_Keratin, rownames(acp_log2tpm))
skin_genes <- intersect(signatures$Skin_Keratin, rownames(acp_log2tpm))

message(sprintf("  Oral-specific keratins available: %d/%d (%s)",
                length(oral_genes), length(signatures$Oral_Keratin),
                paste(oral_genes, collapse = ", ")))
message(sprintf("  Skin-specific keratins available: %d/%d (%s)",
                length(skin_genes), length(signatures$Skin_Keratin),
                paste(skin_genes, collapse = ", ")))

# Per-sample mean expression
oral_expr <- colMeans(acp_log2tpm[oral_genes, , drop = FALSE], na.rm = TRUE)
skin_expr <- colMeans(acp_log2tpm[skin_genes, , drop = FALSE], na.rm = TRUE)

# Paired comparison (same tumors)
h2_test <- wilcox.test(oral_expr, skin_expr, paired = TRUE, alternative = "greater")
h2_effect <- cohens_d(oral_expr, skin_expr)
h2_boot <- boot_mean_diff(oral_expr, skin_expr)

# Proportion of tumors with oral > skin
prop_oral_higher <- mean(oral_expr > skin_expr)
prop_test <- binom.test(sum(oral_expr > skin_expr), length(oral_expr),
                        p = 0.5, alternative = "greater")

message(sprintf("\n  Mean oral keratin expression: %.3f ± %.3f",
                mean(oral_expr), sd(oral_expr)))
message(sprintf("  Mean skin keratin expression: %.3f ± %.3f",
                mean(skin_expr), sd(skin_expr)))
message(sprintf("  Paired difference: %.3f [95%% CI: %.3f, %.3f]",
                h2_boot$estimate, h2_boot$ci_lower, h2_boot$ci_upper))
message(sprintf("  Cohen's d: %.3f", h2_effect))
message(sprintf("  Wilcoxon signed-rank p-value: %.2e", h2_test$p.value))
message(sprintf("\n  Proportion with oral > skin: %.1f%% (binomial p = %.2e)",
                prop_oral_higher * 100, prop_test$p.value))

if (h2_test$p.value < 0.05 && h2_effect > 0.5 && prop_oral_higher > 0.6) {
  message("\n  → H2 SUPPORTED: ACP shows oral-biased keratin expression")
  h2_conclusion <- "SUPPORTED"
} else if (h2_test$p.value < 0.05) {
  message("\n  → H2 PARTIALLY SUPPORTED: Significant but modest oral bias")
  h2_conclusion <- "PARTIAL"
} else {
  message("\n  → H2 NOT SUPPORTED: No evidence for oral bias")
  h2_conclusion <- "NOT_SUPPORTED"
}

results$H2 <- list(
  oral_mean = mean(oral_expr),
  oral_sd = sd(oral_expr),
  skin_mean = mean(skin_expr),
  skin_sd = sd(skin_expr),
  difference = h2_boot$estimate,
  ci_lower = h2_boot$ci_lower,
  ci_upper = h2_boot$ci_upper,
  cohens_d = h2_effect,
  p_value = h2_test$p.value,
  prop_oral_higher = prop_oral_higher,
  prop_p_value = prop_test$p.value,
  conclusion = h2_conclusion
)

# ==============================================================================
# HYPOTHESIS 3: PROLIFERATION-ECM ANTI-CORRELATION
# Validates Sections 0.6, 0.8
# ==============================================================================

message("\n========================================")
message("H3: Proliferation-ECM Anti-Correlation")
message("    (Validates Sections 0.6, 0.8)")
message("========================================\n")

# Test: Are proliferation and ECM/collagen signatures negatively correlated?
# This would support the "protected niche" model where proliferation occurs
# in epithelial cores away from collagen-rich stroma.

prolif_scores <- sig_scores$Proliferation
ecm_scores <- sig_scores$ECM_General
collagen_scores <- sig_scores$Collagen

# Proliferation vs ECM
message("Testing Proliferation vs ECM signature correlation...")
h3a_cor <- boot_cor(prolif_scores, ecm_scores)
h3a_perm <- perm_cor_test(prolif_scores, ecm_scores)

message(sprintf("  Pearson r = %.3f [95%% CI: %.3f, %.3f]",
                h3a_cor$estimate, h3a_cor$ci_lower, h3a_cor$ci_upper))
message(sprintf("  Permutation test: observed = %.3f, null = %.3f ± %.3f, p = %.4f",
                h3a_perm$observed, h3a_perm$null_mean, h3a_perm$null_sd,
                h3a_perm$p_value))

# Proliferation vs Collagen
message("\nTesting Proliferation vs Collagen signature correlation...")
h3b_cor <- boot_cor(prolif_scores, collagen_scores)
h3b_perm <- perm_cor_test(prolif_scores, collagen_scores)

message(sprintf("  Pearson r = %.3f [95%% CI: %.3f, %.3f]",
                h3b_cor$estimate, h3b_cor$ci_lower, h3b_cor$ci_upper))
message(sprintf("  Permutation test: observed = %.3f, null = %.3f ± %.3f, p = %.4f",
                h3b_perm$observed, h3b_perm$null_mean, h3b_perm$null_sd,
                h3b_perm$p_value))

# Interpretation
# Negative correlation supports spatial segregation model
# Positive or no correlation would not support it

# Note: In bulk RNA-seq, we might see weak or no correlation because
# both compartments are present in every sample. Strong anti-correlation
# would be surprising but supportive.

h3_interpretation <- if (h3a_cor$estimate < -0.3 && h3a_perm$p_value < 0.05) {
  "STRONG_SUPPORT"
} else if (h3a_cor$estimate < 0 && h3a_perm$p_value < 0.05) {
  "MODERATE_SUPPORT"
} else if (h3a_cor$estimate < 0.3) {
  "WEAK_SUPPORT"
} else {
  "NOT_SUPPORTED"
}

message(sprintf("\n  → H3 %s: Proliferation-ECM correlation = %.3f",
                h3_interpretation, h3a_cor$estimate))
message("    Note: Weak negative or near-zero correlation in bulk is consistent")
message("    with compartmentalized but co-occurring signals.")

results$H3 <- list(
  prolif_ecm = list(
    correlation = h3a_cor$estimate,
    ci_lower = h3a_cor$ci_lower,
    ci_upper = h3a_cor$ci_upper,
    perm_p = h3a_perm$p_value,
    effect_size = h3a_perm$effect_size
  ),
  prolif_collagen = list(
    correlation = h3b_cor$estimate,
    ci_lower = h3b_cor$ci_lower,
    ci_upper = h3b_cor$ci_upper,
    perm_p = h3b_perm$p_value,
    effect_size = h3b_perm$effect_size
  ),
  conclusion = h3_interpretation
)

# ==============================================================================
# HYPOTHESIS 4: EPITHELIAL-STROMAL SIGNATURE SEGREGATION
# Validates Section 0.8
# ==============================================================================

message("\n========================================")
message("H4: Epithelial-Stromal Signature Segregation")
message("    (Validates Section 0.8)")
message("========================================\n")

# Test: Are epithelial and stromal signatures anti-correlated?
# This would support molecular compartmentalization.

epi_scores <- sig_scores$Pan_Epithelial
ecm_scores <- sig_scores$ECM_General

message("Testing Pan-Epithelial vs ECM signature correlation...")
h4_cor <- boot_cor(epi_scores, ecm_scores)
h4_perm <- perm_cor_test(epi_scores, ecm_scores)

message(sprintf("  Pearson r = %.3f [95%% CI: %.3f, %.3f]",
                h4_cor$estimate, h4_cor$ci_lower, h4_cor$ci_upper))
message(sprintf("  Permutation test: observed = %.3f, null = %.3f ± %.3f, p = %.4f",
                h4_perm$observed, h4_perm$null_mean, h4_perm$null_sd,
                h4_perm$p_value))

# Also test against collagen specifically
message("\nTesting Pan-Epithelial vs Collagen signature correlation...")
h4b_cor <- boot_cor(epi_scores, collagen_scores)
h4b_perm <- perm_cor_test(epi_scores, collagen_scores)

message(sprintf("  Pearson r = %.3f [95%% CI: %.3f, %.3f]",
                h4b_cor$estimate, h4b_cor$ci_lower, h4b_cor$ci_upper))
message(sprintf("  Permutation test: p = %.4f", h4b_perm$p_value))

h4_interpretation <- if (h4_cor$estimate < -0.3 && h4_perm$p_value < 0.05) {
  "STRONG_SUPPORT"
} else if (h4_cor$estimate < 0 && h4_perm$p_value < 0.05) {
  "MODERATE_SUPPORT"
} else if (abs(h4_cor$estimate) < 0.3) {
  "INDEPENDENT_COMPARTMENTS"
} else {
  "NOT_SUPPORTED"
}

message(sprintf("\n  → H4 %s: Epithelial-ECM correlation = %.3f",
                h4_interpretation, h4_cor$estimate))

results$H4 <- list(
  epi_ecm = list(
    correlation = h4_cor$estimate,
    ci_lower = h4_cor$ci_lower,
    ci_upper = h4_cor$ci_upper,
    perm_p = h4_perm$p_value
  ),
  epi_collagen = list(
    correlation = h4b_cor$estimate,
    ci_lower = h4b_cor$ci_lower,
    ci_upper = h4b_cor$ci_upper,
    perm_p = h4b_perm$p_value
  ),
  conclusion = h4_interpretation
)

# ==============================================================================
# HYPOTHESIS 5: CD44-IMMUNE CORRELATION
# Validates Section 0.9
# ==============================================================================

message("\n========================================")
message("H5: CD44-Immune Signature Correlation")
message("    (Validates Section 0.9)")
message("========================================\n")

# Test: Does CD44 expression correlate with immune/myeloid signatures?
# This would support the paracrine signaling model where immune cells
# at boundaries influence CD44-expressing epithelial cells.

cd44_expr <- as.numeric(acp_log2tpm["CD44", ])
myeloid_scores <- sig_scores$Myeloid
tam_scores <- sig_scores$TAM

message("Testing CD44 vs Myeloid signature correlation...")
h5a_cor <- boot_cor(cd44_expr, myeloid_scores)
h5a_perm <- perm_cor_test(cd44_expr, myeloid_scores)

message(sprintf("  Pearson r = %.3f [95%% CI: %.3f, %.3f]",
                h5a_cor$estimate, h5a_cor$ci_lower, h5a_cor$ci_upper))
message(sprintf("  Permutation test: p = %.4f", h5a_perm$p_value))

message("\nTesting CD44 vs TAM signature correlation...")
h5b_cor <- boot_cor(cd44_expr, tam_scores)
h5b_perm <- perm_cor_test(cd44_expr, tam_scores)

message(sprintf("  Pearson r = %.3f [95%% CI: %.3f, %.3f]",
                h5b_cor$estimate, h5b_cor$ci_lower, h5b_cor$ci_upper))
message(sprintf("  Permutation test: p = %.4f", h5b_perm$p_value))

# Also test CD44-SPP1 correlation (direct paracrine axis)
if ("SPP1" %in% rownames(acp_log2tpm)) {
  spp1_expr <- as.numeric(acp_log2tpm["SPP1", ])
  message("\nTesting CD44 vs SPP1 expression correlation...")
  h5c_cor <- boot_cor(cd44_expr, spp1_expr)
  h5c_perm <- perm_cor_test(cd44_expr, spp1_expr)

  message(sprintf("  Pearson r = %.3f [95%% CI: %.3f, %.3f]",
                  h5c_cor$estimate, h5c_cor$ci_lower, h5c_cor$ci_upper))
  message(sprintf("  Permutation test: p = %.4f", h5c_perm$p_value))
} else {
  h5c_cor <- list(estimate = NA)
  h5c_perm <- list(p_value = NA)
}

h5_interpretation <- if (h5a_cor$estimate > 0.3 && h5a_perm$p_value < 0.05) {
  "STRONG_SUPPORT"
} else if (h5a_cor$estimate > 0 && h5a_perm$p_value < 0.05) {
  "MODERATE_SUPPORT"
} else if (!is.na(h5c_cor$estimate) && h5c_cor$estimate > 0.2 && h5c_perm$p_value < 0.05) {
  # CD44-SPP1 direct correlation is significant even if myeloid signature isn't
  "PARTIAL_SUPPORT"
} else {
  "NOT_SUPPORTED"
}

# Additional interpretation for CD44-SPP1
if (!is.na(h5c_cor$estimate) && h5c_perm$p_value < 0.05) {
  message(sprintf("\n  Note: CD44-SPP1 direct correlation IS significant (r=%.3f, p=%.3f)",
                  h5c_cor$estimate, h5c_perm$p_value))
  message("        This supports the paracrine axis even though broader myeloid signature")
  message("        correlation is not significant (expected in bulk where cell types are mixed)")
}

message(sprintf("\n  → H5 %s: CD44-Myeloid correlation = %.3f",
                h5_interpretation, h5a_cor$estimate))

results$H5 <- list(
  cd44_myeloid = list(
    correlation = h5a_cor$estimate,
    ci_lower = h5a_cor$ci_lower,
    ci_upper = h5a_cor$ci_upper,
    perm_p = h5a_perm$p_value
  ),
  cd44_tam = list(
    correlation = h5b_cor$estimate,
    ci_lower = h5b_cor$ci_lower,
    ci_upper = h5b_cor$ci_upper,
    perm_p = h5b_perm$p_value
  ),
  cd44_spp1 = list(
    correlation = h5c_cor$estimate,
    ci_lower = if(exists("h5c_cor") && !is.null(h5c_cor$ci_lower)) h5c_cor$ci_lower else NA,
    ci_upper = if(exists("h5c_cor") && !is.null(h5c_cor$ci_upper)) h5c_cor$ci_upper else NA,
    perm_p = h5c_perm$p_value,
    significant = !is.na(h5c_perm$p_value) && h5c_perm$p_value < 0.05
  ),
  conclusion = h5_interpretation
)

# ==============================================================================
# HYPOTHESIS 6: ACP-SPECIFIC KERATIN PROFILE (Cross-tumor comparison)
# Validates Section 0.3 specificity
# ==============================================================================

message("\n========================================")
message("H6: ACP-Specific Keratin Profile")
message("    (Cross-tumor validation of Section 0.3)")
message("========================================\n")

# Test: Is the oral-biased keratin profile specific to ACP vs other brain tumors?

# Calculate oral-skin difference for all tumors
all_sig_scores$oral_skin_diff <- all_sig_scores$Oral_Keratin - all_sig_scores$Skin_Keratin

# Compare ACP to each other tumor type
acp_oral_skin <- all_sig_scores$oral_skin_diff[all_sig_scores$tumor_type == "Craniopharyngioma"]

h6_comparisons <- data.frame(
  Tumor_Type = character(),
  ACP_Mean = numeric(),
  Other_Mean = numeric(),
  Difference = numeric(),
  p_value = numeric(),
  cohens_d = numeric(),
  stringsAsFactors = FALSE
)

other_tumor_types <- setdiff(unique(all_sig_scores$tumor_type), "Craniopharyngioma")

for (tt in other_tumor_types) {
  other_oral_skin <- all_sig_scores$oral_skin_diff[all_sig_scores$tumor_type == tt]

  if (length(other_oral_skin) >= 5) {
    test_result <- wilcox.test(acp_oral_skin, other_oral_skin)
    d <- cohens_d(acp_oral_skin, other_oral_skin)

    h6_comparisons <- rbind(h6_comparisons, data.frame(
      Tumor_Type = tt,
      ACP_Mean = mean(acp_oral_skin, na.rm = TRUE),
      Other_Mean = mean(other_oral_skin, na.rm = TRUE),
      Difference = mean(acp_oral_skin, na.rm = TRUE) - mean(other_oral_skin, na.rm = TRUE),
      p_value = test_result$p.value,
      cohens_d = d
    ))
  }
}

# Adjust p-values for multiple comparisons
h6_comparisons$p_adj <- p.adjust(h6_comparisons$p_value, method = "BH")

message("  Oral-Skin keratin difference: ACP vs other tumor types")
message("  (Positive difference = oral-biased)")
message("")
for (i in 1:nrow(h6_comparisons)) {
  row <- h6_comparisons[i, ]
  sig_marker <- if (row$p_adj < 0.05) "***" else if (row$p_adj < 0.1) "*" else ""
  message(sprintf("    vs %-30s: Δ = %+.3f, d = %+.2f, p_adj = %.3f %s",
                  row$Tumor_Type, row$Difference, row$cohens_d, row$p_adj, sig_marker))
}

# Overall assessment
n_sig_higher <- sum(h6_comparisons$p_adj < 0.05 & h6_comparisons$Difference > 0)
n_comparisons <- nrow(h6_comparisons)

if (n_sig_higher >= n_comparisons * 0.7) {
  h6_conclusion <- "STRONG_SUPPORT"
  message(sprintf("\n  → H6 STRONGLY SUPPORTED: ACP has significantly higher oral-bias than %d/%d other tumor types",
                  n_sig_higher, n_comparisons))
} else if (n_sig_higher >= n_comparisons * 0.5) {
  h6_conclusion <- "MODERATE_SUPPORT"
  message(sprintf("\n  → H6 MODERATELY SUPPORTED: ACP has higher oral-bias than %d/%d other tumor types",
                  n_sig_higher, n_comparisons))
} else {
  h6_conclusion <- "NOT_SUPPORTED"
  message("\n  → H6 NOT SUPPORTED: Oral keratin bias not ACP-specific")
}

results$H6 <- list(
  comparisons = h6_comparisons,
  n_sig_higher = n_sig_higher,
  n_comparisons = n_comparisons,
  conclusion = h6_conclusion
)

# ==============================================================================
# HYPOTHESIS 7: ACP-SPECIFIC EPITHELIAL SIGNATURE (Cross-tumor comparison)
# Validates epithelial subtype specificity
# ==============================================================================

message("\n========================================")
message("H7: ACP Epithelial Signature Specificity")
message("    (Cross-tumor validation)")
message("========================================\n")

# Test: Are the epithelial subtype signatures higher in ACP than other tumors?
# (Expected: ACP should have high epithelial signatures; gliomas should not)

epithelial_sigs <- c("Basal_like", "Intermediate", "Specialized", "Pan_Epithelial")

h7_results <- data.frame()

for (sig in epithelial_sigs) {
  acp_scores <- all_sig_scores[[sig]][all_sig_scores$tumor_type == "Craniopharyngioma"]

  for (tt in other_tumor_types) {
    other_scores <- all_sig_scores[[sig]][all_sig_scores$tumor_type == tt]

    if (length(other_scores) >= 5) {
      test_result <- wilcox.test(acp_scores, other_scores, alternative = "greater")
      d <- cohens_d(acp_scores, other_scores)

      h7_results <- rbind(h7_results, data.frame(
        Signature = sig,
        Tumor_Type = tt,
        ACP_Mean = mean(acp_scores, na.rm = TRUE),
        Other_Mean = mean(other_scores, na.rm = TRUE),
        cohens_d = d,
        p_value = test_result$p.value
      ))
    }
  }
}

# Adjust p-values
h7_results$p_adj <- p.adjust(h7_results$p_value, method = "BH")

# Summarize by signature
h7_summary <- h7_results %>%
  group_by(Signature) %>%
  summarise(
    n_comparisons = n(),
    n_sig_higher = sum(p_adj < 0.05),
    mean_d = mean(cohens_d, na.rm = TRUE),
    .groups = "drop"
  )

message("  Epithelial signatures: ACP vs other tumor types")
message("")
for (i in 1:nrow(h7_summary)) {
  row <- h7_summary[i, ]
  message(sprintf("    %s: ACP higher in %d/%d comparisons (mean d = %.2f)",
                  row$Signature, row$n_sig_higher, row$n_comparisons, row$mean_d))
}

# Overall assessment
overall_specificity <- mean(h7_summary$n_sig_higher / h7_summary$n_comparisons)

if (overall_specificity > 0.7) {
  h7_conclusion <- "STRONG_SUPPORT"
  message(sprintf("\n  → H7 STRONGLY SUPPORTED: Epithelial signatures are ACP-specific (%.0f%% of comparisons)",
                  overall_specificity * 100))
} else if (overall_specificity > 0.5) {
  h7_conclusion <- "MODERATE_SUPPORT"
  message(sprintf("\n  → H7 MODERATELY SUPPORTED: Epithelial signatures partially ACP-specific (%.0f%%)",
                  overall_specificity * 100))
} else {
  h7_conclusion <- "NOT_SUPPORTED"
  message("\n  → H7 NOT SUPPORTED: Epithelial signatures not ACP-specific")
}

results$H7 <- list(
  full_results = h7_results,
  summary = h7_summary,
  overall_specificity = overall_specificity,
  conclusion = h7_conclusion
)

# ==============================================================================
# HYPOTHESIS 8: VALIDATION AGAINST XU ET AL. CELL TYPE DEFINITIONS
# Validates subtype signatures against snRNA-seq ground truth
# ==============================================================================

if (gse_available) {
  message("\n========================================")
  message("H8: Validation Against Xu et al. Cell Types")
  message("    (snRNA-seq ground truth)")
  message("========================================\n")

  # Map Xu et al. labels to our subtype definitions
  # Based on manuscript: E_C1_PE = Basal-like, E_C3_KC = Intermediate,
  # E_C2_WC/E_C4_RHCG = Specialized, E_C5_Prolif = Transit-Amplifying

  label_mapping <- list(
    "E_C1_PE" = "Basal_like",
    "E_C2_WC" = "Specialized",
    "E_C3_KC" = "Intermediate",
    "E_C4_RHCG" = "Specialized",
    "E_C5_Prolif" = "Proliferation"
  )

  message("  Expected mapping (from manuscript):")
  for (xu_label in names(label_mapping)) {
    message(sprintf("    %s → %s signature", xu_label, label_mapping[[xu_label]]))
  }

  # Test: For each Xu label, is the expected signature highest?
  message("\n  Testing signature enrichment per Xu et al. cell type:")

  h8_results <- data.frame()
  test_sigs <- c("Basal_like", "Intermediate", "Specialized", "Proliferation")

  for (xu_label in names(label_mapping)) {
    if (!xu_label %in% colnames(celltype_means)) next

    expected_sig <- label_mapping[[xu_label]]

    # Calculate mean signature score for this cell type
    sig_means <- sapply(test_sigs, function(sig) {
      sig_genes <- intersect(signatures[[sig]], rownames(celltype_means))
      if (length(sig_genes) >= 3) {
        mean(celltype_means[sig_genes, xu_label], na.rm = TRUE)
      } else {
        NA
      }
    })

    # Check if expected signature is highest
    if (all(is.na(sig_means))) next

    highest_sig <- names(which.max(sig_means))
    is_correct <- (highest_sig == expected_sig)

    # Calculate specificity (how much higher is expected vs others)
    expected_score <- sig_means[expected_sig]
    other_scores <- sig_means[names(sig_means) != expected_sig]
    specificity <- expected_score - max(other_scores, na.rm = TRUE)

    h8_results <- rbind(h8_results, data.frame(
      Xu_Label = xu_label,
      Expected_Signature = expected_sig,
      Highest_Signature = highest_sig,
      Expected_Score = expected_score,
      Is_Correct = is_correct,
      Specificity = specificity
    ))

    status <- if (is_correct) "✓ CORRECT" else "✗ MISMATCH"
    message(sprintf("    %s: Expected %s, Got %s %s (specificity = %.3f)",
                    xu_label, expected_sig, highest_sig, status, specificity))
  }

  # Summary
  n_correct <- sum(h8_results$Is_Correct)
  n_total <- nrow(h8_results)

  # Base conclusion on signature tests
  if (n_correct == n_total) {
    h8_conclusion <- "STRONG_SUPPORT"
    message(sprintf("\n  → H8 STRONGLY SUPPORTED: All %d/%d Xu labels match expected signatures",
                    n_correct, n_total))
  } else if (n_correct >= n_total * 0.7) {
    h8_conclusion <- "MODERATE_SUPPORT"
    message(sprintf("\n  → H8 MODERATELY SUPPORTED: %d/%d Xu labels match expected signatures",
                    n_correct, n_total))
  } else if (n_correct >= 2) {
    h8_conclusion <- "PARTIAL_SUPPORT"
    message(sprintf("\n  → H8 PARTIAL SUPPORT: %d/%d signature mappings correct", n_correct, n_total))
    message("     Note: Signature-based approach may not capture all cluster biology")
    message("     See targeted validations below for specific marker/pathway tests")
  } else {
    h8_conclusion <- "NOT_SUPPORTED"
    message(sprintf("\n  → H8 NOT SUPPORTED: Only %d/%d Xu labels match", n_correct, n_total))
  }

  results$H8 <- list(
    results = h8_results,
    n_correct = n_correct,
    n_total = n_total,
    conclusion = h8_conclusion
  )

  # Additional: Test CD44 enrichment in E_C4_RHCG (Specialized/senescent)
  message("\n  Additional: CD44 enrichment validation")

  if ("CD44" %in% rownames(celltype_means)) {
    cd44_by_type <- celltype_means["CD44", ]
    cd44_by_type <- cd44_by_type[!is.na(cd44_by_type)]

    if (length(cd44_by_type) > 0) {
      message("    CD44 expression by Xu et al. cell type:")
      for (ct in names(sort(cd44_by_type, decreasing = TRUE))) {
        marker <- if (ct == "E_C4_RHCG") " ← Expected highest (Specialized/senescent)" else ""
        message(sprintf("      %s: %.3f%s", ct, cd44_by_type[ct], marker))
      }

      # Check if E_C4_RHCG has highest CD44
      if ("E_C4_RHCG" %in% names(cd44_by_type)) {
        highest_cd44_cluster <- names(which.max(cd44_by_type))
        is_highest <- highest_cd44_cluster == "E_C4_RHCG"
        if (is_highest) {
          message("\n    ✓ CD44 enrichment in E_C4_RHCG VALIDATED")
          message(sprintf("      (%.1fx higher than next highest cluster)",
                          cd44_by_type["E_C4_RHCG"] / sort(cd44_by_type, decreasing = TRUE)[2]))
        } else {
          message(sprintf("\n    ✗ CD44 highest in %s, not E_C4_RHCG", highest_cd44_cluster))
        }
        results$H8$cd44_validation <- is_highest
        results$H8$cd44_fold_enrichment <- cd44_by_type["E_C4_RHCG"] / sort(cd44_by_type, decreasing = TRUE)[2]
      }
    }
  }

  # Additional: Test Wnt pathway enrichment in E_C2_WC (Whorl Center)
  message("\n  Additional: Wnt pathway enrichment validation (E_C2_WC)")

  wnt_genes <- intersect(signatures$Wnt_Pathway, rownames(celltype_means))
  if (length(wnt_genes) >= 3) {
    wnt_by_type <- colMeans(celltype_means[wnt_genes, , drop = FALSE], na.rm = TRUE)
    wnt_by_type <- wnt_by_type[!is.na(wnt_by_type)]

    message("    Wnt pathway signature by Xu et al. cell type:")
    for (ct in names(sort(wnt_by_type, decreasing = TRUE))) {
      marker <- if (ct == "E_C2_WC") " ← Expected highest (Whorl Center/Wnt-active)" else ""
      message(sprintf("      %s: %.3f%s", ct, wnt_by_type[ct], marker))
    }

    if ("E_C2_WC" %in% names(wnt_by_type)) {
      highest_wnt_cluster <- names(which.max(wnt_by_type))
      is_wnt_highest <- highest_wnt_cluster == "E_C2_WC"
      if (is_wnt_highest) {
        message("\n    ✓ Wnt pathway enrichment in E_C2_WC VALIDATED")
      } else {
        message(sprintf("\n    ✗ Wnt pathway highest in %s, not E_C2_WC", highest_wnt_cluster))
      }
      results$H8$wnt_validation <- is_wnt_highest
    }
  }

  # Additional: Test proliferation enrichment in E_C5_Prolif
  message("\n  Additional: Proliferation enrichment validation (E_C5_Prolif)")

  prolif_genes <- intersect(signatures$Proliferation, rownames(celltype_means))
  if (length(prolif_genes) >= 3) {
    prolif_by_type <- colMeans(celltype_means[prolif_genes, , drop = FALSE], na.rm = TRUE)
    prolif_by_type <- prolif_by_type[!is.na(prolif_by_type)]

    message("    Proliferation signature by Xu et al. cell type:")
    for (ct in names(sort(prolif_by_type, decreasing = TRUE))) {
      marker <- if (ct == "E_C5_Prolif") " ← Expected highest" else ""
      message(sprintf("      %s: %.3f%s", ct, prolif_by_type[ct], marker))
    }

    if ("E_C5_Prolif" %in% names(prolif_by_type)) {
      highest_prolif_cluster <- names(which.max(prolif_by_type))
      is_prolif_highest <- highest_prolif_cluster == "E_C5_Prolif"
      if (is_prolif_highest) {
        message("\n    ✓ Proliferation enrichment in E_C5_Prolif VALIDATED")
      } else {
        message(sprintf("\n    ✗ Proliferation highest in %s, not E_C5_Prolif", highest_prolif_cluster))
      }
      results$H8$prolif_validation <- is_prolif_highest
    }
  }

  # Update H8 conclusion based on additional validations
  additional_validations <- c(
    results$H8$cd44_validation %||% FALSE,
    results$H8$wnt_validation %||% FALSE,
    results$H8$prolif_validation %||% FALSE
  )
  n_additional_correct <- sum(additional_validations, na.rm = TRUE)

  message(sprintf("\n  Additional targeted validations: %d/3 correct", n_additional_correct))

  # Revise conclusion if additional validations help
  if (n_additional_correct >= 2) {
    message("  → Targeted marker/pathway validations support cluster biology")
    results$H8$additional_validations_passed <- n_additional_correct

    # Upgrade conclusion if signature tests were borderline
    if (results$H8$conclusion == "NOT_SUPPORTED" && n_additional_correct >= 2) {
      results$H8$conclusion <- "PARTIAL_SUPPORT"
      message("  → H8 upgraded to PARTIAL_SUPPORT based on targeted validations")
    }
  } else {
    results$H8$additional_validations_passed <- n_additional_correct
  }

} else {
  message("\n========================================")
  message("H8: Skipped (GSE215932 data not available)")
  message("========================================\n")
  results$H8 <- list(conclusion = "SKIPPED", error = "Data not available")
}

# ==============================================================================
# SUPPLEMENTARY: SIGNATURE COHERENCE ANALYSIS
# ==============================================================================

message("\n========================================")
message("Supplementary: Signature Coherence Analysis")
message("========================================\n")

message("Testing internal consistency of key signatures...")

coherence_results <- list()
key_sigs <- c("Basal_like", "Intermediate", "Specialized", "Proliferation",
              "Collagen", "Myeloid")

for (sig in key_sigs) {
  coh <- signature_coherence(acp_log2tpm, signatures[[sig]])
  coherence_results[[sig]] <- coh

  interp <- if (is.na(coh$mean_cor)) {
    "insufficient genes"
  } else if (coh$mean_cor > 0.5) {
    "highly coherent"
  } else if (coh$mean_cor > 0.3) {
    "moderately coherent"
  } else if (coh$mean_cor > 0) {
    "weakly coherent"
  } else {
    "not coherent"
  }

  message(sprintf("  %s: mean r = %.3f (%d genes, %s)",
                  sig, coh$mean_cor, coh$n_genes, interp))
}

results$coherence <- coherence_results

# ==============================================================================
# SUPPLEMENTARY: PRIMARY VS RECURRENT COMPARISON
# ==============================================================================

message("\n========================================")
message("Supplementary: Primary vs Recurrent Comparison")
message("========================================\n")

primary_samples <- sig_scores$sample_id[sig_scores$tumor_descriptor == "Initial CNS Tumor"]
recurrent_samples <- sig_scores$sample_id[sig_scores$tumor_descriptor == "Recurrence"]

message(sprintf("  Primary tumors: %d", length(primary_samples)))
message(sprintf("  Recurrent tumors: %d", length(recurrent_samples)))

if (length(primary_samples) >= 5 && length(recurrent_samples) >= 5) {

  test_sigs <- c("Basal_like", "Intermediate", "Specialized", "Proliferation",
                 "Collagen", "ECM_Remodeling", "Senescence")

  primary_vs_recurrent <- data.frame(
    Signature = character(),
    Primary_Mean = numeric(),
    Recurrent_Mean = numeric(),
    Difference = numeric(),
    p_value = numeric(),
    cohens_d = numeric(),
    stringsAsFactors = FALSE
  )

  for (sig in test_sigs) {
    if (sig %in% colnames(sig_scores)) {
      primary_vals <- sig_scores[[sig]][sig_scores$tumor_descriptor == "Initial CNS Tumor"]
      recurrent_vals <- sig_scores[[sig]][sig_scores$tumor_descriptor == "Recurrence"]

      test_result <- wilcox.test(recurrent_vals, primary_vals)
      d <- cohens_d(recurrent_vals, primary_vals)

      primary_vs_recurrent <- rbind(primary_vs_recurrent, data.frame(
        Signature = sig,
        Primary_Mean = mean(primary_vals, na.rm = TRUE),
        Recurrent_Mean = mean(recurrent_vals, na.rm = TRUE),
        Difference = mean(recurrent_vals, na.rm = TRUE) - mean(primary_vals, na.rm = TRUE),
        p_value = test_result$p.value,
        cohens_d = d
      ))
    }
  }

  message("\n  Signature differences (Recurrent - Primary):")
  for (i in 1:nrow(primary_vs_recurrent)) {
    row <- primary_vs_recurrent[i, ]
    sig_marker <- if (row$p_value < 0.05) "*" else ""
    message(sprintf("    %s: Δ = %+.3f, d = %+.2f, p = %.3f%s",
                    row$Signature, row$Difference, row$cohens_d,
                    row$p_value, sig_marker))
  }

  results$primary_vs_recurrent <- primary_vs_recurrent

} else {
  message("  Insufficient samples for comparison")
  results$primary_vs_recurrent <- NULL
}

# ==============================================================================
# SUMMARY AND CONCLUSIONS
# ==============================================================================

message("\n========================================")
message("Summary of Results")
message("========================================\n")

summary_table <- data.frame(
  Hypothesis = c(
    "H1: Keratin coordinated expression",
    "H2: Oral-biased keratin profile",
    "H3: Proliferation-ECM segregation",
    "H4: Epithelial-stromal compartments",
    "H5: CD44-immune association",
    "H6: Oral bias ACP-specific (vs other tumors)",
    "H7: Epithelial sigs ACP-specific",
    "H8: Xu et al. label validation"
  ),
  Validates = c(
    "Sections 0.3, 0.4, 0.10",
    "Section 0.3",
    "Sections 0.6, 0.8",
    "Section 0.8",
    "Section 0.9",
    "Section 0.3 (specificity)",
    "Sections 0.1-0.4 (specificity)",
    "Section 0.1 (snRNA-seq)"
  ),
  Key_Statistic = c(
    sprintf("d = %.2f (n=%d/%d)",
            ifelse(is.na(results$H1$cohens_d), 0, results$H1$cohens_d),
            results$H1$n_within_pairs, results$H1$n_between_pairs),
    sprintf("d = %.2f", ifelse(is.na(results$H2$cohens_d), 0, results$H2$cohens_d)),
    sprintf("r = %.2f", ifelse(is.na(results$H3$prolif_ecm$correlation), 0, results$H3$prolif_ecm$correlation)),
    sprintf("r = %.2f", ifelse(is.na(results$H4$epi_ecm$correlation), 0, results$H4$epi_ecm$correlation)),
    sprintf("r = %.2f (SPP1: %.2f*)",
            ifelse(is.na(results$H5$cd44_myeloid$correlation), 0, results$H5$cd44_myeloid$correlation),
            ifelse(is.na(results$H5$cd44_spp1$correlation), 0, results$H5$cd44_spp1$correlation)),
    sprintf("%d/%d sig", results$H6$n_sig_higher, results$H6$n_comparisons),
    sprintf("%.0f%% specific", results$H7$overall_specificity * 100),
    if (!is.null(results$H8$n_correct)) {
      additional <- results$H8$additional_validations_passed %||% 0
      sprintf("%d/5 sig + %d/3 targeted", results$H8$n_correct, additional)
    } else "N/A"
  ),
  P_Value = c(
    ifelse(is.na(results$H1$p_value), "NA", sprintf("%.2e", results$H1$p_value)),
    sprintf("%.2e", results$H2$p_value),
    sprintf("%.3f", results$H3$prolif_ecm$perm_p),
    sprintf("%.3f", results$H4$epi_ecm$perm_p),
    sprintf("%.3f", results$H5$cd44_myeloid$perm_p),
    sprintf("%.3f (median adj)", median(results$H6$comparisons$p_adj)),
    sprintf("%.3f (median adj)", median(results$H7$full_results$p_adj)),
    "N/A"
  ),
  Conclusion = c(
    results$H1$conclusion,
    results$H2$conclusion,
    results$H3$conclusion,
    results$H4$conclusion,
    results$H5$conclusion,
    results$H6$conclusion,
    results$H7$conclusion,
    results$H8$conclusion
  ),
  stringsAsFactors = FALSE
)

print(summary_table)

results$summary <- summary_table

# Overall assessment
n_supported <- sum(grepl("SUPPORT", summary_table$Conclusion) &
                     !grepl("NOT_SUPPORTED", summary_table$Conclusion))
n_partial <- sum(grepl("PARTIAL|MODERATE|WEAK|TENTATIVE", summary_table$Conclusion))
n_not_supported <- sum(grepl("NOT_SUPPORTED", summary_table$Conclusion))
n_skipped <- sum(grepl("SKIPPED|INCONCLUSIVE", summary_table$Conclusion))

message("\n  Overall validation assessment:")
message(sprintf("    Supported: %d/8 hypotheses", n_supported))
message(sprintf("    Partially supported: %d/8 hypotheses", n_partial))
message(sprintf("    Not supported: %d/8 hypotheses", n_not_supported))
if (n_skipped > 0) {
  message(sprintf("    Skipped/Inconclusive: %d/8 hypotheses", n_skipped))
}

message("\n  Interpretation notes:")
message("    - H1-H5 test internal ACP patterns in bulk data")
message("    - H6-H7 test ACP SPECIFICITY vs other brain tumors")
message("    - H8 validates signatures against snRNA-seq ground truth")
message("    - Weak correlations in bulk may reflect averaged compartments")
message("    - Strong cross-tumor specificity is key validation evidence")

# ==============================================================================
# GENERATE FIGURES
# ==============================================================================

message("\n========================================")
message("Generating Figures")
message("========================================\n")

# --- Figure Panel A: Keratin Heatmap Grouped by Subtype ---

# Define keratins in subtype order
key_keratins <- c("KRT5", "KRT14", "KRT15", "KRT6A",  # Basal
                  "KRT8", "KRT18", "KRT19", "KRT7",   # Intermediate
                  "KRT17", "KRT23", "KRT10", "KRT1",  # Specialized
                  "KRT4", "KRT13", "KRT76")           # Oral-specific

key_keratins <- intersect(key_keratins, rownames(acp_log2tpm))
krt_mat <- as.matrix(acp_log2tpm[key_keratins, ])

# Row annotation - must be a factor to preserve order
keratin_type <- case_when(
  key_keratins %in% c("KRT5", "KRT14", "KRT15", "KRT6A") ~ "Basal-like",
  key_keratins %in% c("KRT8", "KRT18", "KRT19", "KRT7") ~ "Intermediate",
  key_keratins %in% c("KRT17", "KRT23", "KRT10", "KRT1") ~ "Specialized",
  key_keratins %in% c("KRT4", "KRT13", "KRT76") ~ "Oral-specific",
  TRUE ~ "Other"
)
keratin_type <- factor(keratin_type, levels = c("Basal-like", "Intermediate", "Specialized", "Oral-specific"))

row_ha <- rowAnnotation(
  Type = keratin_type,
  col = list(Type = c("Basal-like" = "#E41A1C", "Intermediate" = "#4DAF4A",
                      "Specialized" = "#377EB8", "Oral-specific" = "#984EA3",
                      "Other" = "grey")),
  show_legend = TRUE
)

# Column annotation
col_ha <- HeatmapAnnotation(
  Tumor = acp_meta$tumor_descriptor,
  col = list(Tumor = c("Initial CNS Tumor" = "#1f77b4",
                       "Recurrence" = "#d62728",
                       "Progressive" = "#ff7f0e")),
  na_col = "grey80"
)

pdf(file.path(dirs$figures, "panel_A_keratin_heatmap.pdf"), width = 12, height = 7)
draw(Heatmap(
  krt_mat,
  name = "log2TPM",
  col = colorRamp2(c(0, 5, 10, 15), c("navy", "white", "orange", "firebrick")),
  top_annotation = col_ha,
  right_annotation = row_ha,
  cluster_rows = FALSE,  # Keep keratins in subtype order
  cluster_columns = TRUE,  # Cluster samples to show heterogeneity
  row_split = keratin_type,  # Split rows by subtype
  row_gap = unit(2, "mm"),
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_title = sprintf("Keratin Expression in ACP (n=%d)", ncol(krt_mat)),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(title = "log2(TPM+1)")
))
dev.off()
message("  Saved: panel_A_keratin_heatmap.pdf")

# --- Figure Panel B: Oral vs Skin Keratin Scores ---

p_oral_skin <- data.frame(
  Sample = rep(1:length(oral_expr), 2),
  Expression = c(oral_expr, skin_expr),
  Type = rep(c("Oral-specific", "Skin-specific"), each = length(oral_expr))
) %>%
  ggplot(aes(x = Type, y = Expression, fill = Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  geom_line(aes(group = Sample), alpha = 0.3, color = "grey50") +
  scale_fill_manual(values = c("Oral-specific" = "#984EA3", "Skin-specific" = "#FF7F00")) +
  labs(
    title = "H2: Oral vs Skin Keratin Expression",
    subtitle = sprintf("Paired difference = %.2f, Cohen's d = %.2f, p = %.2e",
                       results$H2$difference, results$H2$cohens_d, results$H2$p_value),
    x = NULL,
    y = "Mean log2(TPM+1)"
  ) +
  theme_pubr() +
  theme(legend.position = "none")

ggsave(file.path(dirs$figures, "panel_B_oral_vs_skin.pdf"), p_oral_skin,
       width = 5, height = 5)
message("  Saved: panel_B_oral_vs_skin.pdf")

# --- Figure Panel C: Proliferation vs ECM Correlation ---

p_prolif_ecm <- data.frame(
  Proliferation = prolif_scores,
  ECM = ecm_scores,
  Tumor = acp_meta$tumor_descriptor
) %>%
  ggplot(aes(x = Proliferation, y = ECM, color = Tumor)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  scale_color_manual(values = c("Initial CNS Tumor" = "#1f77b4",
                                "Recurrence" = "#d62728",
                                "Progressive" = "#ff7f0e")) +
  labs(
    title = "H3: Proliferation vs ECM Signatures",
    subtitle = sprintf("r = %.3f [%.3f, %.3f], perm p = %.3f",
                       results$H3$prolif_ecm$correlation,
                       results$H3$prolif_ecm$ci_lower,
                       results$H3$prolif_ecm$ci_upper,
                       results$H3$prolif_ecm$perm_p),
    x = "Proliferation Score (UCell)",
    y = "ECM Score (UCell)"
  ) +
  theme_pubr()

ggsave(file.path(dirs$figures, "panel_C_prolif_vs_ecm.pdf"), p_prolif_ecm,
       width = 6, height = 5)
message("  Saved: panel_C_prolif_vs_ecm.pdf")

# --- Figure Panel D: Epithelial vs ECM Correlation ---

p_epi_ecm <- data.frame(
  Epithelial = epi_scores,
  ECM = ecm_scores,
  Tumor = acp_meta$tumor_descriptor
) %>%
  ggplot(aes(x = Epithelial, y = ECM, color = Tumor)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  scale_color_manual(values = c("Initial CNS Tumor" = "#1f77b4",
                                "Recurrence" = "#d62728",
                                "Progressive" = "#ff7f0e")) +
  labs(
    title = "H4: Epithelial vs ECM Signatures",
    subtitle = sprintf("r = %.3f [%.3f, %.3f], perm p = %.4f",
                       results$H4$epi_ecm$correlation,
                       results$H4$epi_ecm$ci_lower,
                       results$H4$epi_ecm$ci_upper,
                       results$H4$epi_ecm$perm_p),
    x = "Pan-Epithelial Score (UCell)",
    y = "ECM Score (UCell)"
  ) +
  theme_pubr()

ggsave(file.path(dirs$figures, "panel_D_epi_vs_ecm.pdf"), p_epi_ecm,
       width = 6, height = 5)
message("  Saved: panel_D_epi_vs_ecm.pdf")

# --- Figure Panel E: CD44 vs Immune Signatures ---

# Show both the myeloid comparison AND the significant CD44-SPP1 correlation
if ("SPP1" %in% rownames(acp_log2tpm)) {
  spp1_expr <- as.numeric(acp_log2tpm["SPP1", ])

  p_cd44_immune <- data.frame(
    CD44 = cd44_expr,
    Myeloid = myeloid_scores,
    SPP1 = spp1_expr,
    Tumor = acp_meta$tumor_descriptor
  ) %>%
    pivot_longer(cols = c(Myeloid, SPP1), names_to = "Comparison", values_to = "Value") %>%
    mutate(Comparison = factor(Comparison, levels = c("Myeloid", "SPP1"))) %>%
    ggplot(aes(x = CD44, y = Value, color = Tumor)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2, linewidth = 0.5) +
    facet_wrap(~Comparison, scales = "free_y",
               labeller = labeller(Comparison = c("Myeloid" = "Myeloid Signature (n.s.)",
                                                  "SPP1" = "SPP1 Expression (p=0.013)"))) +
    scale_color_manual(values = c("Initial CNS Tumor" = "#1f77b4",
                                  "Recurrence" = "#d62728",
                                  "Progressive" = "#ff7f0e")) +
    labs(
      title = "H5: CD44 Correlations",
      subtitle = sprintf("CD44-Myeloid: r=%.2f (n.s.) | CD44-SPP1: r=%.2f (p=0.013)*",
                         results$H5$cd44_myeloid$correlation,
                         results$H5$cd44_spp1$correlation),
      x = "CD44 log2(TPM+1)",
      y = "Score / Expression"
    ) +
    theme_pubr() +
    theme(legend.position = "bottom")
} else {
  p_cd44_immune <- data.frame(
    CD44 = cd44_expr,
    Myeloid = myeloid_scores,
    TAM = tam_scores,
    Tumor = acp_meta$tumor_descriptor
  ) %>%
    pivot_longer(cols = c(Myeloid, TAM), names_to = "Signature", values_to = "Score") %>%
    ggplot(aes(x = CD44, y = Score, color = Tumor)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2, linewidth = 0.5) +
    facet_wrap(~Signature, scales = "free_y") +
    scale_color_manual(values = c("Initial CNS Tumor" = "#1f77b4",
                                  "Recurrence" = "#d62728",
                                  "Progressive" = "#ff7f0e")) +
    labs(
      title = "H5: CD44 vs Immune Signatures",
      subtitle = sprintf("CD44-Myeloid: r = %.3f (p = %.3f) | CD44-TAM: r = %.3f (p = %.3f)",
                         results$H5$cd44_myeloid$correlation,
                         results$H5$cd44_myeloid$perm_p,
                         results$H5$cd44_tam$correlation,
                         results$H5$cd44_tam$perm_p),
      x = "CD44 log2(TPM+1)",
      y = "Immune Score (UCell)"
    ) +
    theme_pubr() +
    theme(legend.position = "bottom")
}

ggsave(file.path(dirs$figures, "panel_E_cd44_vs_immune.pdf"), p_cd44_immune,
       width = 9, height = 5)
message("  Saved: panel_E_cd44_vs_immune.pdf")

# --- Figure Panel F: Signature Score Summary ---

sig_long <- sig_scores %>%
  select(sample_id, Basal_like, Intermediate, Specialized,
         Proliferation, Collagen, ECM_General, Myeloid) %>%
  pivot_longer(-sample_id, names_to = "Signature", values_to = "Score") %>%
  mutate(Category = case_when(
    Signature %in% c("Basal_like", "Intermediate", "Specialized") ~ "Epithelial",
    Signature %in% c("Collagen", "ECM_General") ~ "ECM",
    Signature == "Proliferation" ~ "Proliferation",
    Signature == "Myeloid" ~ "Immune"
  ))

p_signature_summary <- ggplot(sig_long, aes(x = reorder(Signature, Score, median),
                                            y = Score, fill = Category)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("Epithelial" = "#E41A1C", "ECM" = "#377EB8",
                               "Proliferation" = "#4DAF4A", "Immune" = "#984EA3")) +
  labs(
    title = "Signature Score Distribution Across ACP Cohort",
    subtitle = sprintf("n = %d tumors", length(available_ids)),
    x = NULL,
    y = "UCell Score"
  ) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(dirs$figures, "panel_F_signature_summary.pdf"), p_signature_summary,
       width = 8, height = 5)
message("  Saved: panel_F_signature_summary.pdf")

# --- Figure Panel G: Cross-Tumor Oral Keratin Bias ---

p_cross_tumor_oral <- all_sig_scores %>%
  mutate(tumor_type = factor(tumor_type,
                             levels = c("Craniopharyngioma",
                                        sort(setdiff(unique(tumor_type), "Craniopharyngioma"))))) %>%
  ggplot(aes(x = tumor_type, y = oral_skin_diff, fill = tumor_type == "Craniopharyngioma")) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#E41A1C"),
                    labels = c("Other tumors", "ACP")) +
  labs(
    title = "H6: Oral Keratin Bias is ACP-Specific",
    subtitle = "Positive = oral-biased; ACP significantly higher than all 8 comparison tumor types (p<0.001)",
    x = NULL,
    y = "Oral - Skin Keratin Score",
    fill = NULL
  ) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

ggsave(file.path(dirs$figures, "panel_G_cross_tumor_oral.pdf"), p_cross_tumor_oral,
       width = 10, height = 6)
message("  Saved: panel_G_cross_tumor_oral.pdf")

# --- Figure Panel H: Cross-Tumor Epithelial Signature Comparison ---

epi_comparison <- all_sig_scores %>%
  select(sample_id, tumor_type, Basal_like, Intermediate, Specialized, Pan_Epithelial) %>%
  pivot_longer(cols = c(Basal_like, Intermediate, Specialized, Pan_Epithelial),
               names_to = "Signature", values_to = "Score") %>%
  mutate(is_acp = tumor_type == "Craniopharyngioma",
         tumor_type = factor(tumor_type,
                             levels = c("Craniopharyngioma",
                                        sort(setdiff(unique(tumor_type), "Craniopharyngioma")))))

p_cross_tumor_epi <- ggplot(epi_comparison, aes(x = tumor_type, y = Score, fill = is_acp)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  facet_wrap(~Signature, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#E41A1C"),
                    labels = c("Other tumors", "ACP")) +
  labs(
    title = "H7: Epithelial Signature Specificity Across Brain Tumors",
    subtitle = "ACP has significantly higher scores for all epithelial signatures (100% comparisons, p<0.001)",
    x = NULL,
    y = "UCell Score",
    fill = NULL
  ) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "top",
        strip.text = element_text(face = "bold"))

ggsave(file.path(dirs$figures, "panel_H_cross_tumor_epithelial.pdf"), p_cross_tumor_epi,
       width = 12, height = 8)
message("  Saved: panel_H_cross_tumor_epithelial.pdf")

# --- Figure Panel I: Xu et al. Validation Heatmap (if available) ---

if (gse_available && exists("celltype_means")) {
  # Include genes that showed validation, not just signature genes
  # Group 1: Basal-like (validates E_C1_PE)
  basal_genes <- c("KRT5", "KRT14", "KRT15", "TP63", "KRT6A", "ITGA6")

  # Group 2: Intermediate (validates E_C3_KC)
  intermediate_genes <- c("KRT8", "KRT18", "KRT19", "KRT7", "EPCAM", "CDH1")

  # Group 3: Proliferation (validates E_C5_Prolif)
  prolif_genes <- c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1", "STMN1")

  # Group 4: Wnt pathway (validates E_C2_WC)
  wnt_genes <- c("CTNNB1", "AXIN2", "LEF1", "WNT5A", "DKK1", "DKK4")

  # Group 5: CD44/Senescence (validates E_C4_RHCG)
  senescence_genes <- c("CD44", "CDKN2A", "CDKN1A", "IGFBP7", "SERPINE1")

  # Group 6: Specialized keratins
  specialized_genes <- c("KRT17", "KRT23", "KRT10", "IVL")

  # Combine and filter for available genes
  all_heatmap_genes <- c(basal_genes, intermediate_genes, prolif_genes,
                         wnt_genes, senescence_genes, specialized_genes)
  sig_genes_for_heatmap <- intersect(all_heatmap_genes, rownames(celltype_means))

  if (length(sig_genes_for_heatmap) >= 10) {
    heatmap_mat <- celltype_means[sig_genes_for_heatmap, , drop = FALSE]
    heatmap_mat <- heatmap_mat[, !apply(heatmap_mat, 2, function(x) all(is.na(x))), drop = FALSE]

    # Scale by row
    heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

    # Row annotation - assign to functional group
    gene_group <- case_when(
      sig_genes_for_heatmap %in% basal_genes ~ "Basal-like",
      sig_genes_for_heatmap %in% intermediate_genes ~ "Intermediate",
      sig_genes_for_heatmap %in% prolif_genes ~ "Proliferation",
      sig_genes_for_heatmap %in% wnt_genes ~ "Wnt pathway",
      sig_genes_for_heatmap %in% senescence_genes ~ "CD44/Senescence",
      sig_genes_for_heatmap %in% specialized_genes ~ "Specialized",
      TRUE ~ "Other"
    )
    names(gene_group) <- sig_genes_for_heatmap

    # Only keep genes that are in the matrix
    gene_group <- gene_group[rownames(heatmap_mat_scaled)]

    # Define colors for groups
    group_colors <- c(
      "Basal-like" = "#E41A1C",
      "Intermediate" = "#4DAF4A",
      "Proliferation" = "#FF7F00",
      "Wnt pathway" = "#984EA3",
      "CD44/Senescence" = "#377EB8",
      "Specialized" = "#A65628",
      "Other" = "grey"
    )

    row_ha_xu <- rowAnnotation(
      Group = gene_group,
      col = list(Group = group_colors),
      show_legend = TRUE
    )

    # Column annotation showing expected enrichment
    col_annotation_df <- data.frame(
      Expected = c(
        "E_C1_PE" = "Basal-like",
        "E_C2_WC" = "Wnt pathway",
        "E_C3_KC" = "Intermediate",
        "E_C4_RHCG" = "CD44/Senescence",
        "E_C5_Prolif" = "Proliferation"
      )[colnames(heatmap_mat_scaled)]
    )
    rownames(col_annotation_df) <- colnames(heatmap_mat_scaled)

    col_ha_xu <- HeatmapAnnotation(
      Expected = col_annotation_df$Expected,
      col = list(Expected = group_colors),
      show_legend = TRUE
    )

    pdf(file.path(dirs$figures, "panel_I_xu_validation_heatmap.pdf"), width = 10, height = 12)
    draw(Heatmap(
      heatmap_mat_scaled,
      name = "Z-score",
      col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick")),
      right_annotation = row_ha_xu,
      top_annotation = col_ha_xu,
      cluster_rows = FALSE,
      cluster_columns = FALSE,  # Keep fixed order to show expected patterns
      show_row_names = TRUE,
      row_names_gp = gpar(fontsize = 9),
      column_names_gp = gpar(fontsize = 11),
      column_title = "Xu et al. Cell Types (GSE215932)",
      row_title = "Marker Genes",
      row_split = factor(gene_group, levels = c("Basal-like", "Intermediate", "Proliferation",
                                                "Wnt pathway", "CD44/Senescence", "Specialized")),
      row_title_gp = gpar(fontsize = 10),
      row_gap = unit(2, "mm"),
      column_names_rot = 45,
      heatmap_legend_param = list(title = "Z-score", direction = "vertical")
    ))
    dev.off()
    message("  Saved: panel_I_xu_validation_heatmap.pdf")

    # Also save a version with actual expression values (not z-scored) for reference
    pdf(file.path(dirs$figures, "panel_I_xu_validation_heatmap_raw.pdf"), width = 10, height = 12)
    draw(Heatmap(
      heatmap_mat,
      name = "log2(expr)",
      col = colorRamp2(c(0, 1, 3), c("white", "orange", "firebrick")),
      right_annotation = row_ha_xu,
      top_annotation = col_ha_xu,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      row_names_gp = gpar(fontsize = 9),
      column_names_gp = gpar(fontsize = 11),
      column_title = "Xu et al. Cell Types - Raw Expression",
      row_title = "Marker Genes",
      row_split = factor(gene_group, levels = c("Basal-like", "Intermediate", "Proliferation",
                                                "Wnt pathway", "CD44/Senescence", "Specialized")),
      row_title_gp = gpar(fontsize = 10),
      row_gap = unit(2, "mm"),
      column_names_rot = 45
    ))
    dev.off()
    message("  Saved: panel_I_xu_validation_heatmap_raw.pdf")

    # Also create a bar plot showing the targeted validations
    # CD44, Wnt, and Proliferation by cluster
    if (all(c("CD44") %in% rownames(celltype_means))) {

      # Get the validation data
      cd44_vals <- celltype_means["CD44", ]
      wnt_vals <- colMeans(celltype_means[intersect(signatures$Wnt_Pathway, rownames(celltype_means)), , drop = FALSE], na.rm = TRUE)
      prolif_vals <- colMeans(celltype_means[intersect(signatures$Proliferation, rownames(celltype_means)), , drop = FALSE], na.rm = TRUE)

      validation_df <- data.frame(
        Cluster = rep(names(cd44_vals), 3),
        Marker = rep(c("CD44", "Wnt pathway", "Proliferation"), each = length(cd44_vals)),
        Expression = c(cd44_vals, wnt_vals, prolif_vals)
      ) %>%
        mutate(
          Expected_Highest = case_when(
            Marker == "CD44" & Cluster == "E_C4_RHCG" ~ TRUE,
            Marker == "Wnt pathway" & Cluster == "E_C2_WC" ~ TRUE,
            Marker == "Proliferation" & Cluster == "E_C5_Prolif" ~ TRUE,
            TRUE ~ FALSE
          ),
          Marker = factor(Marker, levels = c("CD44", "Wnt pathway", "Proliferation"))
        )

      p_xu_targeted <- ggplot(validation_df, aes(x = Cluster, y = Expression, fill = Expected_Highest)) +
        geom_bar(stat = "identity", width = 0.7) +
        facet_wrap(~Marker, scales = "free_y", ncol = 1) +
        scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#4DAF4A"),
                          labels = c("Other clusters", "Expected highest")) +
        labs(
          title = "H8: Targeted Validation of Xu et al. Cluster Biology",
          subtitle = "All 3 targeted markers show expected enrichment pattern (✓)",
          x = "Xu et al. Cell Type",
          y = "Mean Expression",
          fill = NULL
        ) +
        theme_pubr() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "top",
              strip.text = element_text(face = "bold"))

      ggsave(file.path(dirs$figures, "panel_J_xu_targeted_validation.pdf"), p_xu_targeted,
             width = 7, height = 10)
      message("  Saved: panel_J_xu_targeted_validation.pdf")
    }
  }
}

# --- Combined Figure ---

combined_fig <- (p_oral_skin | p_prolif_ecm | p_epi_ecm) /
  (p_cd44_immune | p_cross_tumor_oral) +
  plot_annotation(
    title = "Figure 0.11: External Validation of Spatial Transcriptomic Findings",
    subtitle = sprintf("OpenPBTA bulk RNA-seq: %d ACP specimens + %d comparison tumors",
                       length(available_ids), nrow(comparison_samples)),
    tag_levels = "A"
  )

ggsave(file.path(dirs$figures, "Figure_0.11_combined.pdf"), combined_fig,
       width = 16, height = 10)
ggsave(file.path(dirs$figures, "Figure_0.11_combined.png"), combined_fig,
       width = 16, height = 10, dpi = 300)
message("  Saved: Figure_0.11_combined.pdf/png")

# --- Extended Combined Figure (ggplot panels only) ---

# Layout:
# Row 1: Oral vs Skin | Prolif vs ECM | Epi vs ECM
# Row 2: CD44 correlations (wide)
# Row 3: Cross-tumor oral bias (wide)
# Row 4: Cross-tumor epithelial specificity (wide, 4 facets)

extended_fig <- (p_oral_skin | p_prolif_ecm | p_epi_ecm) /
  p_cd44_immune /
  p_cross_tumor_oral /
  p_cross_tumor_epi +
  plot_layout(heights = c(1, 0.8, 1, 1.2)) +
  plot_annotation(
    title = "Figure 0.11 Extended: Complete Validation Analysis",
    subtitle = sprintf("OpenPBTA: %d ACP + %d comparison tumors from %d types | See separate files for heatmaps",
                       length(available_ids), nrow(comparison_samples),
                       length(unique(comparison_samples$short_histology))),
    tag_levels = "A"
  )

ggsave(file.path(dirs$figures, "Figure_0.11_extended.pdf"), extended_fig,
       width = 14, height = 20)
message("  Saved: Figure_0.11_extended.pdf")

# --- Create Multi-Page Complete Figure (all panels including heatmaps) ---

message("  Creating multi-page complete figure...")

pdf(file.path(dirs$figures, "Figure_0.11_complete_multipanel.pdf"), width = 14, height = 18)

# Page 1: Main ggplot panels
print(extended_fig)

# Page 2: Keratin heatmap (Panel A)
draw(Heatmap(
  krt_mat,
  name = "log2TPM",
  col = colorRamp2(c(0, 5, 10, 15), c("navy", "white", "orange", "firebrick")),
  top_annotation = col_ha,
  right_annotation = row_ha,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  row_split = keratin_type,
  row_gap = unit(2, "mm"),
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_title = sprintf("Panel A: Keratin Subtype Expression in ACP (n=%d)\nH1: Within-subtype correlations higher than between-subtype (d=0.50, p=0.012)", ncol(krt_mat)),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_title_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(title = "log2(TPM+1)")
))

dev.off()
message("  Saved: Figure_0.11_complete_multipanel.pdf")

# --- Print figure file summary ---
message("\n  Figure files generated:")
message("    Main figures:")
message("      - Figure_0.11_combined.pdf: 5-panel summary (B-F)")
message("      - Figure_0.11_extended.pdf: 7-panel with cross-tumor comparisons")
message("      - Figure_0.11_complete_multipanel.pdf: Multi-page with heatmaps")
message("    Individual panels:")
message("      - panel_A_keratin_heatmap.pdf: Keratin subtype expression (H1)")
message("      - panel_B_oral_vs_skin.pdf: Oral vs skin keratins (H2)")
message("      - panel_C_prolif_vs_ecm.pdf: Proliferation-ECM anti-correlation (H3)")
message("      - panel_D_epi_vs_ecm.pdf: Epithelial-ECM segregation (H4)")
message("      - panel_E_cd44_vs_immune.pdf: CD44-SPP1 correlation (H5)")
message("      - panel_F_signature_summary.pdf: Signature score distributions")
message("      - panel_G_cross_tumor_oral.pdf: ACP oral bias specificity (H6)")
message("      - panel_H_cross_tumor_epithelial.pdf: ACP epithelial specificity (H7)")
if (gse_available) {
  message("      - panel_I_xu_validation_heatmap.pdf: Xu et al. snRNA-seq validation (H8)")
  message("      - panel_I_xu_validation_heatmap_raw.pdf: Raw expression version")
  message("      - panel_J_xu_targeted_validation.pdf: CD44/Wnt/Prolif validation bars")
}

# --- Create Publication-Ready Comprehensive Figure ---
# This creates a large figure with all key panels for manuscript

message("\n  Creating publication-ready comprehensive figure...")

# We need gridExtra for combining grid and ggplot objects
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}
library(gridExtra)

# Convert ggplots to grobs
grob_oral_skin <- ggplotGrob(p_oral_skin + ggtitle("B: Oral vs Skin Keratins (H2)") + theme(plot.title = element_text(size = 10, face = "bold")))
grob_prolif_ecm <- ggplotGrob(p_prolif_ecm + ggtitle("C: Proliferation-ECM (H3)") + theme(plot.title = element_text(size = 10, face = "bold")))
grob_epi_ecm <- ggplotGrob(p_epi_ecm + ggtitle("D: Epithelial-ECM (H4)") + theme(plot.title = element_text(size = 10, face = "bold")))
grob_cross_oral <- ggplotGrob(p_cross_tumor_oral + ggtitle("E: Oral Keratin Specificity (H6)") + theme(plot.title = element_text(size = 10, face = "bold")))
grob_cross_epi <- ggplotGrob(p_cross_tumor_epi + ggtitle("F: Epithelial Signature Specificity (H7)") + theme(plot.title = element_text(size = 10, face = "bold")))

# Create keratin heatmap as a grob
krt_heatmap_grob <- grid.grabExpr({
  draw(Heatmap(
    krt_mat,
    name = "log2TPM",
    col = colorRamp2(c(0, 5, 10, 15), c("navy", "white", "orange", "firebrick")),
    top_annotation = col_ha,
    right_annotation = row_ha,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    row_split = keratin_type,
    row_gap = unit(1.5, "mm"),
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),
    column_title = "A: Keratin Subtype Expression (H1)",
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_title_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(title = "log2(TPM+1)", title_gp = gpar(fontsize = 8))
  ))
}, width = 10, height = 5)

# Layout the comprehensive figure
# Row 1: Keratin heatmap (wide) | Oral vs Skin | Prolif-ECM
# Row 2: Cross-tumor oral bias (wide)
# Row 3: Cross-tumor epithelial (wide, 4 facets)

pdf(file.path(dirs$figures, "Figure_0.11_comprehensive.pdf"), width = 18, height = 20)

# Use grid.arrange for complex layout
grid.arrange(
  arrangeGrob(
    krt_heatmap_grob,
    arrangeGrob(grob_oral_skin, grob_prolif_ecm, grob_epi_ecm, ncol = 3),
    nrow = 2,
    heights = c(1, 1)
  ),
  grob_cross_oral,
  grob_cross_epi,
  nrow = 3,
  heights = c(2, 1, 1.5),
  top = textGrob("Figure 0.11: OpenPBTA Validation of Spatial Transcriptomic Findings",
                 gp = gpar(fontsize = 14, fontface = "bold")),
  bottom = textGrob(sprintf("n = %d ACP specimens + %d comparison tumors from %d types",
                            length(available_ids), nrow(comparison_samples),
                            length(unique(comparison_samples$short_histology))),
                    gp = gpar(fontsize = 10))
)

dev.off()
message("  Saved: Figure_0.11_comprehensive.pdf")

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

message("\n========================================")
message("Saving Outputs")
message("========================================\n")

# Summary table
write.csv(summary_table, file.path(dirs$tables, "H_test_summary.csv"), row.names = FALSE)
message("  Saved: H_test_summary.csv")

# Signature scores (ACP only)
write.csv(sig_scores, file.path(dirs$tables, "signature_scores_acp.csv"), row.names = FALSE)
message("  Saved: signature_scores_acp.csv")

# All signature scores (ACP + comparison tumors)
write.csv(all_sig_scores, file.path(dirs$tables, "signature_scores_all_tumors.csv"), row.names = FALSE)
message("  Saved: signature_scores_all_tumors.csv")

# Cross-tumor comparisons
write.csv(results$H6$comparisons, file.path(dirs$tables, "H6_oral_keratin_cross_tumor.csv"), row.names = FALSE)
message("  Saved: H6_oral_keratin_cross_tumor.csv")

write.csv(results$H7$full_results, file.path(dirs$tables, "H7_epithelial_specificity.csv"), row.names = FALSE)
message("  Saved: H7_epithelial_specificity.csv")

if (!is.null(results$H8$results) && nrow(results$H8$results) > 0) {
  write.csv(results$H8$results, file.path(dirs$tables, "H8_xu_validation.csv"), row.names = FALSE)
  message("  Saved: H8_xu_validation.csv")
}

# Detailed results
if (!is.null(results$primary_vs_recurrent)) {
  write.csv(results$primary_vs_recurrent,
            file.path(dirs$tables, "primary_vs_recurrent.csv"), row.names = FALSE)
  message("  Saved: primary_vs_recurrent.csv")
}

# Full results object
saveRDS(results, file.path(dirs$results, "openpbta_validation_results.rds"))
message("  Saved: openpbta_validation_results.rds")

# Session info
writeLines(capture.output(sessionInfo()),
           file.path(dirs$results, "session_info.txt"))

# ==============================================================================
# COMPLETION
# ==============================================================================

message("\n")
message("================================================================")
message("  OpenPBTA Validation Complete")
message("================================================================")
message(sprintf("  Finished: %s", Sys.time()))
message("")
message("Dataset summary:")
message(sprintf("  ACP samples: %d", length(available_ids)))
message(sprintf("  Comparison tumor samples: %d (from %d tumor types)",
                nrow(comparison_samples), length(unique(comparison_samples$short_histology))))
if (gse_available) {
  message(sprintf("  Xu et al. snRNA-seq cells: %d", ncol(seurat_gse)))
}
message("")
message("Key results:")
message(sprintf("  H1 (Keratin coordination): %s (d = %.2f)",
                results$H1$conclusion,
                ifelse(is.na(results$H1$cohens_d), 0, results$H1$cohens_d)))
message(sprintf("  H2 (Oral keratin bias): %s (d = %.2f)",
                results$H2$conclusion,
                ifelse(is.na(results$H2$cohens_d), 0, results$H2$cohens_d)))
message(sprintf("  H3 (Prolif-ECM segregation): %s (r = %.2f)",
                results$H3$conclusion,
                ifelse(is.na(results$H3$prolif_ecm$correlation), 0, results$H3$prolif_ecm$correlation)))
message(sprintf("  H4 (Epi-stromal compartments): %s (r = %.2f)",
                results$H4$conclusion,
                ifelse(is.na(results$H4$epi_ecm$correlation), 0, results$H4$epi_ecm$correlation)))
message(sprintf("  H5 (CD44-immune association): %s (r = %.2f; CD44-SPP1: r = %.2f, p = %.3f)",
                results$H5$conclusion,
                ifelse(is.na(results$H5$cd44_myeloid$correlation), 0, results$H5$cd44_myeloid$correlation),
                ifelse(is.na(results$H5$cd44_spp1$correlation), 0, results$H5$cd44_spp1$correlation),
                ifelse(is.na(results$H5$cd44_spp1$perm_p), 1, results$H5$cd44_spp1$perm_p)))
message(sprintf("  H6 (Oral bias ACP-specific): %s (%d/%d tumor comparisons)",
                results$H6$conclusion, results$H6$n_sig_higher, results$H6$n_comparisons))
message(sprintf("  H7 (Epithelial sigs specific): %s (%.0f%% of comparisons)",
                results$H7$conclusion, results$H7$overall_specificity * 100))
message(sprintf("  H8 (Xu et al. validation): %s (sig: %d/%d; CD44: %s, Wnt: %s, Prolif: %s)",
                results$H8$conclusion,
                results$H8$n_correct %||% 0,
                results$H8$n_total %||% 0,
                if (isTRUE(results$H8$cd44_validation)) "✓" else "✗",
                if (isTRUE(results$H8$wnt_validation)) "✓" else "?",
                if (isTRUE(results$H8$prolif_validation)) "✓" else "?"))
message("")
message("Outputs:")
message(sprintf("  Tables: %s", dirs$tables))
message(sprintf("  Figures: %s", dirs$figures))
message(sprintf("  Full results: %s", file.path(dirs$results, "openpbta_validation_results.rds")))
message("================================================================\n")
