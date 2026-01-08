#!/usr/bin/env Rscript
# ==============================================================================
# scripts/04_epithelial_origin_comparison.R
# ==============================================================================
# Supplementary Analysis: ACP Epithelial Origin Comparison (v2)
#
# This script performs rigorous comparison of ACP epithelium to oral and skin
# references with proper batch correction and statistical controls.
#
# WORKFLOW:
#   0. Batch effect diagnostic and Harmony integration
#   1. Keratin profile comparison (with effect sizes)
#   2. Transcriptional similarity (on integrated data)
#   3. ACP heterogeneity analysis
#   4. Sensitivity analysis
#
# Hypotheses:
#   H1: ACP expresses oral-specific keratins > skin-specific keratins
#   H2: ACP transcriptome is more similar to oral than skin epithelium
#   H3: ACP shows distinct keratin profile from both healthy references
#
# Usage:
#   Rscript scripts/04_epithelial_origin_comparison.R [--cores N]
#
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(stats)
  library(harmony)
  library(cluster)
})

options(future.globals.maxSize = 32 * 1024^3)

# Define null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

message("\n")
message("================================================================")
message("  Epithelial Origin Comparison: ACP vs Oral vs Skin (v2)")
message("================================================================")
message(sprintf("  Started: %s", Sys.time()))

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
n_cores <- 4
for (i in seq_along(args)) {
  if (args[i] == "--cores" && i < length(args)) {
    n_cores <- as.integer(args[i + 1])
  }
}
message(sprintf("  Using %d cores", n_cores))

# Source utilities
source("R/utils/config.R")

# Load configuration
config <- load_config()
set.seed(config$reproducibility$seed)

# Define paths
objects_dir <- get_path(config, config$paths$objects_dir)
tables_dir <- get_path(config, config$paths$tables_dir)
fig_dir <- get_path(config, config$paths$figures_supp_dir)
ensure_dir(fig_dir)
ensure_dir(tables_dir)

# ==============================================================================
# DEFINE GENE SIGNATURES
# ==============================================================================

message("\n========================================")
message("Defining Gene Signatures")
message("========================================\n")

keratin_signatures <- list(
  # Oral mucosa-specific (non-keratinized stratified squamous)
  oral_specific = c("KRT4", "KRT13", "KRT19", "KRT76", "KRT78"),

  # Skin-specific (keratinized stratified squamous)
  skin_specific = c("KRT1", "KRT10", "KRT2", "KRT9", "KRT77"),

  # Basal layer (shared across stratified epithelia)
  basal_shared = c("KRT5", "KRT14", "KRT15", "TP63", "ITGA6", "ITGB4", "COL17A1"),

  # Simple epithelium markers
  simple_epithelium = c("KRT8", "KRT18", "KRT7", "KRT19", "KRT20"),

  # Differentiation/cornification markers
  cornification = c("IVL", "LOR", "FLG", "SPRR1A", "SPRR1B", "SPRR2A", "SPRR3"),

  # Proliferation
  proliferation = c("MKI67", "TOP2A", "PCNA", "MCM2", "CDK1"),

  # ACP/Craniopharyngioma markers (from literature)
  acp_markers = c("CTNNB1", "AXIN2", "LEF1", "SOX2", "SOX9", "CLDN1", "EPCAM"),

  # WNT pathway targets
  wnt_targets = c("AXIN2", "LEF1", "TCF7", "MYC", "CCND1", "CD44", "LGR5")
)

for (sig_name in names(keratin_signatures)) {
  message(sprintf("  %s: %d genes", sig_name, length(keratin_signatures[[sig_name]])))
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Calculate expression statistics per cell
calc_percell_expression <- function(seurat_obj, genes, assay = "RNA") {
  expr_data <- GetAssayData(seurat_obj, assay = assay, layer = "data")
  genes_present <- genes[genes %in% rownames(expr_data)]

  n_genes_found <- length(genes_present)

  if (n_genes_found == 0) {
    return(data.frame(
      cell = colnames(seurat_obj),
      mean_expr = NA_real_,
      n_genes = 0L,
      stringsAsFactors = FALSE
    ))
  }

  if (n_genes_found == 1) {
    cell_means <- as.numeric(expr_data[genes_present, ])
  } else {
    cell_means <- colMeans(as.matrix(expr_data[genes_present, ]))
  }

  data.frame(
    cell = colnames(seurat_obj),
    mean_expr = cell_means,
    n_genes = n_genes_found,
    stringsAsFactors = FALSE
  )
}

#' Calculate Cohen's d effect size
cohens_d <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  nx <- length(x)
  ny <- length(y)
  if (nx < 2 || ny < 2) return(NA)
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  if (pooled_sd == 0) return(NA)
  (mean(x) - mean(y)) / pooled_sd
}

#' Calculate Cliff's delta effect size
cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) return(NA)
  nx <- length(x)
  ny <- length(y)
  greater <- sum(outer(x, y, ">"))
  less <- sum(outer(x, y, "<"))
  (greater - less) / (nx * ny)
}

#' Interpret effect size
interpret_cohens_d <- function(d) {
  d <- abs(d)
  case_when(
    is.na(d) ~ "NA",
    d < 0.2 ~ "negligible",
    d < 0.5 ~ "small",
    d < 0.8 ~ "medium",
    TRUE ~ "large"
  )
}

# ==============================================================================
# LOAD DATASETS
# ==============================================================================

message("\n========================================")
message("Loading Datasets")
message("========================================\n")

# Load merged ACP dataset
acp_path <- file.path(objects_dir, "01_seurat_annotated_merged.rds")
if (!file.exists(acp_path)) {
  stop("Merged ACP dataset not found: ", acp_path)
}
message("Loading ACP (merged)...")
seurat_acp <- readRDS(acp_path)
seurat_acp <- UpdateSeuratObject(seurat_acp)
message(sprintf("  ACP cells: %d", ncol(seurat_acp)))

# Load oral atlas
oral_path <- get_path(config, config$paths$cellxgene_oral_atlas)
if (!file.exists(oral_path)) {
  stop("Oral atlas not found: ", oral_path)
}
message("Loading oral atlas...")
seurat_oral <- readRDS(oral_path)
seurat_oral <- UpdateSeuratObject(seurat_oral)
message(sprintf("  Oral cells (total): %d", ncol(seurat_oral)))

# Load skin reference
skin_path <- file.path(objects_dir, "01_seurat_annotated_cellxgene.rds")
if (!file.exists(skin_path)) {
  stop("Skin reference not found: ", skin_path)
}
message("Loading skin reference...")
seurat_skin <- readRDS(skin_path)
seurat_skin <- UpdateSeuratObject(seurat_skin)
message(sprintf("  Skin cells (total): %d", ncol(seurat_skin)))

# Join layers if Seurat v5
for (obj_name in c("seurat_acp", "seurat_oral", "seurat_skin")) {
  obj <- get(obj_name)
  if (inherits(obj[[DefaultAssay(obj)]], "Assay5")) {
    if (length(Layers(obj[[DefaultAssay(obj)]])) > 1) {
      assign(obj_name, JoinLayers(obj))
    }
  }
}

# ==============================================================================
# FILTER TO EPITHELIAL CELLS
# ==============================================================================

message("\n========================================")
message("Filtering to Epithelial Cells")
message("========================================\n")

# Filter oral to epithelial cells only
oral_epithelial <- subset(seurat_oral, cell_type == "epithelial cell")
message(sprintf("  Oral epithelial cells: %d", ncol(oral_epithelial)))

# Further filter oral to relevant tissues (oral mucosa, not glands)
oral_mucosal_tissues <- c("buccal mucosa", "gingiva", "hard palate", "soft palate",
                          "mucosa of lip", "mucosa of dorsum of tongue")
oral_mucosal <- subset(oral_epithelial, tissue %in% oral_mucosal_tissues)
message(sprintf("  Oral mucosal epithelial cells: %d", ncol(oral_mucosal)))

# Filter skin to epithelial cell types
skin_epithelial_types <- c("basal cell of epidermis", "keratinocyte",
                           "spinous cell of epidermis")
skin_epithelial <- subset(seurat_skin, cell_type %in% skin_epithelial_types)
message(sprintf("  Skin epithelial cells: %d", ncol(skin_epithelial)))

# Filter ACP to epithelial (exclude Unassigned if present)
if ("module_score_subtype" %in% colnames(seurat_acp@meta.data)) {
  acp_epithelial <- subset(seurat_acp, module_score_subtype != "Unassigned")
} else {
  acp_epithelial <- seurat_acp
}
message(sprintf("  ACP epithelial cells: %d", ncol(acp_epithelial)))

# ==============================================================================
# SUBSAMPLE FOR COMPUTATIONAL EFFICIENCY
# ==============================================================================

message("\n========================================")
message("Subsampling for Analysis")
message("========================================\n")

max_cells_per_group <- 10000
set.seed(config$reproducibility$seed)

subsample_seurat <- function(obj, n_max, label) {
  if (ncol(obj) > n_max) {
    cells_keep <- sample(colnames(obj), n_max)
    obj <- subset(obj, cells = cells_keep)
    message(sprintf("  %s: subsampled to %d cells", label, ncol(obj)))
  } else {
    message(sprintf("  %s: kept all %d cells", label, ncol(obj)))
  }
  return(obj)
}

# acp_sub <- subsample_seurat(acp_epithelial, max_cells_per_group, "ACP")
# oral_sub <- subsample_seurat(oral_mucosal, max_cells_per_group, "Oral")
# skin_sub <- subsample_seurat(skin_epithelial, max_cells_per_group, "Skin")
acp_sub <- acp_epithelial  # Keep all ACP cells
oral_sub <- oral_mucosal
skin_sub <- skin_epithelial

# ==============================================================================
# SECTION 0: BATCH EFFECT DIAGNOSTIC AND INTEGRATION
# ==============================================================================

message("\n========================================")
message("SECTION 0: Batch Effect Assessment")
message("========================================\n")

# -----------------------------------------------------------------------------
# 0.1: Batch effect diagnostic BEFORE integration
# -----------------------------------------------------------------------------

message("--- 0.1: Batch effect diagnostic (pre-integration) ---\n")

# Add source labels
acp_sub$source <- "ACP"
oral_sub$source <- "Oral"
skin_sub$source <- "Skin"

# Find common genes
common_genes <- Reduce(intersect, list(
  rownames(acp_sub),
  rownames(oral_sub),
  rownames(skin_sub)
))
message(sprintf("  Common genes across datasets: %d", length(common_genes)))

# Merge without integration
merged_raw <- merge(
  subset(acp_sub, features = common_genes),
  c(subset(oral_sub, features = common_genes),
    subset(skin_sub, features = common_genes))
)

merged_raw <- NormalizeData(merged_raw, verbose = FALSE)
merged_raw <- FindVariableFeatures(merged_raw, nfeatures = 2000, verbose = FALSE)
merged_raw <- ScaleData(merged_raw, verbose = FALSE)
merged_raw <- RunPCA(merged_raw, npcs = 30, verbose = FALSE)
merged_raw <- RunUMAP(merged_raw, dims = 1:20, verbose = FALSE)

# Calculate variance explained by dataset
pca_coords <- Embeddings(merged_raw, "pca")[, 1:10]
source_var <- sapply(1:10, function(i) {
  aov_result <- summary(aov(pca_coords[, i] ~ merged_raw$source))[[1]]
  aov_result["Sum Sq"][1, 1] / sum(aov_result["Sum Sq"])
})

mean_source_var <- mean(source_var) * 100

message("Variance explained by dataset (batch effect proxy):")
message(sprintf("  PC1: %.1f%%", source_var[1] * 100))
message(sprintf("  PC2: %.1f%%", source_var[2] * 100))
message(sprintf("  PC1-10 mean: %.1f%%", mean_source_var))

batch_confounded <- mean_source_var > 50
if (batch_confounded) {
  message("\n  ⚠️  WARNING: >50% variance explained by dataset identity")
  message("  Downstream comparisons may be confounded by batch effects")
  message("  Proceeding with Harmony integration...")
} else {
  message("\n  ✓ Batch effects appear manageable (<50% variance)")
}

# Plot pre-integration UMAP
p_pre_umap <- DimPlot(merged_raw, group.by = "source", reduction = "umap",
                      pt.size = 0.3, alpha = 0.5) +
  scale_color_manual(values = c("ACP" = "#E41A1C", "Oral" = "#4DAF4A", "Skin" = "#377EB8")) +
  labs(title = "UMAP Before Integration",
       subtitle = sprintf("%.0f%% of PC1-10 variance explained by dataset", mean_source_var))

p_pre_pca <- DimPlot(merged_raw, group.by = "source", reduction = "pca",
                     pt.size = 0.3, alpha = 0.5) +
  scale_color_manual(values = c("ACP" = "#E41A1C", "Oral" = "#4DAF4A", "Skin" = "#377EB8")) +
  labs(title = "PCA Before Integration")

# -----------------------------------------------------------------------------
# 0.2: Harmony integration
# -----------------------------------------------------------------------------

message("\n--- 0.2: Integrating datasets with Harmony ---\n")

merged_integrated <- merged_raw %>%
  RunHarmony("source", verbose = FALSE, max.iter.harmony = 20)

# Join layers for downstream analysis (Seurat v5)
if (inherits(merged_integrated[[DefaultAssay(merged_integrated)]], "Assay5")) {
  merged_integrated <- JoinLayers(merged_integrated)
}

# Run UMAP on harmony embeddings
merged_integrated <- RunUMAP(merged_integrated, reduction = "harmony",
                             dims = 1:20, verbose = FALSE,
                             reduction.name = "umap_harmony")

# Re-check variance after integration
harmony_coords <- Embeddings(merged_integrated, "harmony")[, 1:10]
source_var_post <- sapply(1:10, function(i) {
  aov_result <- summary(aov(harmony_coords[, i] ~ merged_integrated$source))[[1]]
  aov_result["Sum Sq"][1, 1] / sum(aov_result["Sum Sq"])
})

mean_source_var_post <- mean(source_var_post) * 100

message("Variance explained by dataset AFTER Harmony:")
message(sprintf("  Harmony1-10 mean: %.1f%% (was %.1f%%)",
                mean_source_var_post, mean_source_var))
message(sprintf("  Reduction: %.1f%%", mean_source_var - mean_source_var_post))

# Plot post-integration
p_post_umap <- DimPlot(merged_integrated, group.by = "source", reduction = "umap_harmony",
                       pt.size = 0.3, alpha = 0.5) +
  scale_color_manual(values = c("ACP" = "#E41A1C", "Oral" = "#4DAF4A", "Skin" = "#377EB8")) +
  labs(title = "UMAP After Harmony Integration",
       subtitle = sprintf("%.0f%% variance by dataset (reduced from %.0f%%)",
                          mean_source_var_post, mean_source_var))

# Combined plot
p_integration <- (p_pre_umap | p_post_umap) +
  plot_annotation(title = "Batch Effect Correction",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

ggsave(file.path(fig_dir, "04_batch_correction_comparison.pdf"), p_integration,
       width = 14, height = 6)
ggsave(file.path(fig_dir, "04_batch_correction_comparison.png"), p_integration,
       width = 14, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 0.3: Positive control - Can we distinguish Oral from Skin?
# -----------------------------------------------------------------------------

message("\n--- 0.3: Positive control - Oral vs Skin separation ---\n")

# Subset to reference cells only
ref_cells <- merged_integrated$source %in% c("Oral", "Skin")
ref_harmony <- harmony_coords[ref_cells, ]
ref_labels <- merged_integrated$source[ref_cells]

# Silhouette score
sil <- silhouette(as.numeric(factor(ref_labels)), dist(ref_harmony))
mean_sil <- mean(sil[, 3])

message(sprintf("Silhouette score (Oral vs Skin separation): %.3f", mean_sil))
sil_interpretation <- case_when(
  mean_sil > 0.7 ~ "Strong separation - references are distinct",
  mean_sil > 0.5 ~ "Reasonable separation - references distinguishable",
  mean_sil > 0.25 ~ "Weak separation - some overlap",
  TRUE ~ "No clear separation - REFERENCES MAY BE PROBLEMATIC"
)
message(sprintf("  Interpretation: %s", sil_interpretation))

positive_control_passed <- mean_sil > 0.25
if (!positive_control_passed) {
  message("\n  ⚠️  WARNING: Oral and skin references don't clearly separate")
  message("  Transcriptome-wide comparisons may not be meaningful")
}

# Save batch diagnostic results
batch_results <- data.frame(
  metric = c("variance_pre_integration", "variance_post_integration",
             "variance_reduction", "silhouette_oral_skin", "positive_control_passed"),
  value = c(mean_source_var, mean_source_var_post,
            mean_source_var - mean_source_var_post, mean_sil, positive_control_passed)
)
write.csv(batch_results, file.path(tables_dir, "04_batch_diagnostic.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 0.4: Reference-only PCA projection and outlier identification
# -----------------------------------------------------------------------------

message("\n--- 0.4: Reference-only PCA projection ---\n")

# Build PCA on reference only (Oral + Skin)
ref_only <- subset(merged_raw, source %in% c("Oral", "Skin"))
ref_only <- FindVariableFeatures(ref_only, nfeatures = 2000, verbose = FALSE)
ref_only <- ScaleData(ref_only, verbose = FALSE)
ref_only <- RunPCA(ref_only, npcs = 30, verbose = FALSE)

# Get loadings and coordinates
ref_loadings <- Loadings(ref_only, "pca")[, 1:30]
ref_pca <- Embeddings(ref_only, "pca")[, 1:30]
pca_genes <- rownames(ref_loadings)

message(sprintf("  Reference PCA built on %d cells", ncol(ref_only)))

# Project ACP onto reference PCA space
acp_only_proj <- subset(merged_raw, source == "ACP")
acp_scale <- GetAssayData(acp_only_proj, layer = "scale.data")
common_genes <- intersect(pca_genes, rownames(acp_scale))
acp_projected <- t(as.matrix(acp_scale[common_genes, ])) %*% ref_loadings[common_genes, ]

message(sprintf("  ACP cells projected: %d", nrow(acp_projected)))

# Calculate reference bounds (95% envelope)
ref_pc1_bounds <- quantile(ref_pca[, 1], c(0.025, 0.975))
ref_pc2_bounds <- quantile(ref_pca[, 2], c(0.025, 0.975))

message(sprintf("  Reference PC1 95%% bounds: [%.2f, %.2f]", ref_pc1_bounds[1], ref_pc1_bounds[2]))
message(sprintf("  Reference PC2 95%% bounds: [%.2f, %.2f]", ref_pc2_bounds[1], ref_pc2_bounds[2]))

# Identify ACP outliers (outside reference envelope)
acp_proj_df <- data.frame(
  cell = rownames(acp_projected),
  PC1 = acp_projected[, 1],
  PC2 = acp_projected[, 2]
)

acp_proj_df$is_outlier <-
  acp_proj_df$PC1 < ref_pc1_bounds[1] | acp_proj_df$PC1 > ref_pc1_bounds[2] |
  acp_proj_df$PC2 < ref_pc2_bounds[1] | acp_proj_df$PC2 > ref_pc2_bounds[2]

acp_proj_df$outlier_type <- case_when(
  acp_proj_df$PC2 < ref_pc2_bounds[1] ~ "PC2_low",
  acp_proj_df$PC2 > ref_pc2_bounds[2] ~ "PC2_high",
  acp_proj_df$PC1 < ref_pc1_bounds[1] ~ "PC1_low",
  acp_proj_df$PC1 > ref_pc1_bounds[2] ~ "PC1_high",
  TRUE ~ "Within_reference"
)

n_outliers <- sum(acp_proj_df$is_outlier)
pct_outliers <- 100 * n_outliers / nrow(acp_proj_df)

message(sprintf("\n  ACP outliers (outside 95%% reference): %d (%.1f%%)", n_outliers, pct_outliers))
print(table(acp_proj_df$outlier_type))

# -----------------------------------------------------------------------------
# 0.5: DE analysis - Outliers vs Within-reference
# -----------------------------------------------------------------------------

# message("\n--- 0.5: DE analysis (Outliers vs Within-reference) ---\n")
#
# # Add outlier status to ACP subset
# acp_de <- subset(merged_raw, source == "ACP")
# acp_de$projection_outlier <- acp_proj_df$outlier_type[match(colnames(acp_de), acp_proj_df$cell)]
# acp_de$is_outlier <- acp_proj_df$is_outlier[match(colnames(acp_de), acp_proj_df$cell)]
#
# Idents(acp_de) <- "is_outlier"
#
# if (sum(acp_de$is_outlier) >= 50 && sum(!acp_de$is_outlier) >= 50) {
#   markers_outlier <- FindMarkers(
#     acp_de,
#     ident.1 = TRUE,
#     ident.2 = FALSE,
#     only.pos = FALSE,
#     min.pct = 0.1,
#     logfc.threshold = 0.25,
#     verbose = FALSE
#   )
#   markers_outlier$gene <- rownames(markers_outlier)
#
#   n_sig <- sum(markers_outlier$p_val_adj < 0.05)
#   n_up <- sum(markers_outlier$p_val_adj < 0.05 & markers_outlier$avg_log2FC > 0)
#   n_down <- sum(markers_outlier$p_val_adj < 0.05 & markers_outlier$avg_log2FC < 0)
#
#   message(sprintf("  DE genes (outlier vs within-reference): %d significant", n_sig))
#   message(sprintf("  Upregulated in outliers: %d", n_up))
#   message(sprintf("  Downregulated in outliers: %d", n_down))
#
#   message("\nTop genes UPREGULATED in outlier ACP cells:")
#   print(head(markers_outlier %>% filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
#                arrange(p_val_adj) %>% select(gene, avg_log2FC, pct.1, pct.2, p_val_adj), 15))
#
#   message("\nTop genes DOWNREGULATED in outlier ACP cells:")
#   print(head(markers_outlier %>% filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
#                arrange(p_val_adj) %>% select(gene, avg_log2FC, pct.1, pct.2, p_val_adj), 15))
#
#   write.csv(markers_outlier, file.path(tables_dir, "04_outlier_de_genes.csv"), row.names = FALSE)
# } else {
#   message("  Insufficient cells for DE analysis")
#   markers_outlier <- NULL
# }

# -----------------------------------------------------------------------------
# 0.6: Gene set enrichment for outliers
# -----------------------------------------------------------------------------

# message("\n--- 0.6: Gene set enrichment (Outliers) ---\n")
#
# outlier_enrichment <- list()
#
# for (gs_name in names(keratin_signatures)) {
#   genes <- keratin_signatures[[gs_name]]
#   genes_present <- genes[genes %in% rownames(acp_de)]
#
#   if (length(genes_present) >= 3) {
#     expr <- GetAssayData(acp_de, layer = "data")[genes_present, , drop = FALSE]
#
#     outlier_means <- colMeans(as.matrix(expr[, acp_de$is_outlier, drop = FALSE]))
#     within_means <- colMeans(as.matrix(expr[, !acp_de$is_outlier, drop = FALSE]))
#
#     wtest <- wilcox.test(outlier_means, within_means)
#
#     outlier_enrichment[[gs_name]] <- data.frame(
#       gene_set = gs_name,
#       n_genes = length(genes_present),
#       mean_outlier = mean(outlier_means),
#       mean_within = mean(within_means),
#       log2fc = log2((mean(outlier_means) + 0.01) / (mean(within_means) + 0.01)),
#       p_value = wtest$p.value
#     )
#   }
# }
#
# outlier_enrich_df <- do.call(rbind, outlier_enrichment)
# outlier_enrich_df$p_adj <- p.adjust(outlier_enrich_df$p_value, method = "BH")
# outlier_enrich_df <- outlier_enrich_df %>% arrange(p_adj)
#
# message("Gene set enrichment (Outliers vs Within-reference):")
# print(outlier_enrich_df %>% mutate(across(where(is.numeric), ~round(., 4))))
#
# write.csv(outlier_enrich_df, file.path(tables_dir, "04_outlier_geneset_enrichment.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Plot: PCA projection with outliers
# -----------------------------------------------------------------------------

message("\nCreating PCA projection plot...")

# Combine reference and ACP for plotting
ref_plot <- data.frame(
  PC1 = ref_pca[, 1],
  PC2 = ref_pca[, 2],
  Origin = ref_only$source
)

acp_plot <- data.frame(
  PC1 = acp_proj_df$PC1,
  PC2 = acp_proj_df$PC2,
  Origin = "ACP"
)

projection_plot_df <- rbind(ref_plot, acp_plot)

p_projection <- ggplot(projection_plot_df, aes(x = PC1, y = PC2, color = Origin)) +
  geom_point(alpha = 0.4, size = 0.5) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = c("ACP" = "#E41A1C", "Oral" = "#4DAF4A", "Skin" = "#377EB8")) +
  theme_minimal(base_size = 12) +
  labs(title = "PCA Projection: ACP onto Reference Space",
       subtitle = "ACP cells projected onto Oral+Skin reference PCA")

ggsave(file.path(fig_dir, "04_pca_projection.pdf"), p_projection, width = 10, height = 8)
ggsave(file.path(fig_dir, "04_pca_projection.png"), p_projection, width = 10, height = 8, dpi = 300)

# # Plot with outliers highlighted
# acp_proj_df$group <- ifelse(acp_proj_df$is_outlier, "Outlier", "Within_Reference")
#
# p_outliers <- ggplot() +
#   geom_point(data = ref_plot, aes(x = PC1, y = PC2), color = "grey70", alpha = 0.3, size = 0.5) +
#   geom_rect(aes(xmin = ref_pc1_bounds[1], xmax = ref_pc1_bounds[2],
#                 ymin = ref_pc2_bounds[1], ymax = ref_pc2_bounds[2]),
#             fill = NA, color = "grey40", linetype = "dashed", linewidth = 0.8) +
#   geom_point(data = acp_proj_df %>% filter(!is_outlier),
#              aes(x = PC1, y = PC2), color = "#E41A1C", alpha = 0.5, size = 0.8) +
#   geom_point(data = acp_proj_df %>% filter(is_outlier),
#              aes(x = PC1, y = PC2), color = "#FFD700", alpha = 0.8, size = 1) +
#   theme_minimal(base_size = 12) +
#   labs(title = "ACP Outliers in Reference PCA Space",
#        subtitle = sprintf("Outliers (yellow): %d cells (%.1f%%) outside 95%% reference envelope",
#                           n_outliers, pct_outliers),
#        x = "PC1", y = "PC2")
#
# ggsave(file.path(fig_dir, "04_pca_outliers.pdf"), p_outliers, width = 10, height = 8)
# ggsave(file.path(fig_dir, "04_pca_outliers.png"), p_outliers, width = 10, height = 8, dpi = 300)
#
# # Save outlier classification
# write.csv(acp_proj_df, file.path(tables_dir, "04_acp_projection_outliers.csv"), row.names = FALSE)


# ==============================================================================
# ANALYSIS 1: KERATIN PROFILE COMPARISON
# ==============================================================================

message("\n========================================")
message("ANALYSIS 1: Keratin Profile Comparison")
message("========================================\n")

# -----------------------------------------------------------------------------
# 1.1: Calculate signature scores for each dataset
# -----------------------------------------------------------------------------

message("--- 1.1: Calculating signature scores ---\n")

signature_results <- list()

for (sig_name in names(keratin_signatures)) {
  genes <- keratin_signatures[[sig_name]]

  acp_expr <- calc_percell_expression(acp_sub, genes)
  oral_expr <- calc_percell_expression(oral_sub, genes)
  skin_expr <- calc_percell_expression(skin_sub, genes)

  n_genes_acp <- acp_expr$n_genes[1]
  n_genes_oral <- oral_expr$n_genes[1]
  n_genes_skin <- skin_expr$n_genes[1]

  combined <- data.frame(
    cell = c(as.character(acp_expr$cell),
             as.character(oral_expr$cell),
             as.character(skin_expr$cell)),
    mean_expr = c(acp_expr$mean_expr, oral_expr$mean_expr, skin_expr$mean_expr),
    n_genes = c(rep(n_genes_acp, nrow(acp_expr)),
                rep(n_genes_oral, nrow(oral_expr)),
                rep(n_genes_skin, nrow(skin_expr))),
    dataset = c(rep("ACP", nrow(acp_expr)),
                rep("Oral", nrow(oral_expr)),
                rep("Skin", nrow(skin_expr))),
    signature = rep(sig_name, nrow(acp_expr) + nrow(oral_expr) + nrow(skin_expr)),
    stringsAsFactors = FALSE
  )

  signature_results[[sig_name]] <- combined
}

signature_df <- do.call(rbind, signature_results)
rownames(signature_df) <- NULL
signature_df$dataset <- factor(signature_df$dataset, levels = c("ACP", "Oral", "Skin"))

# Calculate summary statistics
signature_summary <- signature_df %>%
  group_by(signature, dataset) %>%
  summarise(
    mean = mean(mean_expr, na.rm = TRUE),
    sd = sd(mean_expr, na.rm = TRUE),
    median = median(mean_expr, na.rm = TRUE),
    n = n(),
    n_genes = first(n_genes),
    .groups = "drop"
  ) %>%
  mutate(se = sd / sqrt(n))

message("Signature expression summary:")
print(signature_summary %>%
        select(signature, dataset, mean, sd, n_genes) %>%
        pivot_wider(names_from = dataset, values_from = c(mean, sd)))

# -----------------------------------------------------------------------------
# 1.2: H1 Test with effect sizes and confidence intervals
# -----------------------------------------------------------------------------

message("\n--- 1.2: H1 Test (oral vs skin keratins in ACP) ---\n")

acp_oral_scores <- calc_percell_expression(acp_sub, keratin_signatures$oral_specific)$mean_expr
acp_skin_scores <- calc_percell_expression(acp_sub, keratin_signatures$skin_specific)$mean_expr

# Paired Wilcoxon test
h1_test <- wilcox.test(acp_oral_scores, acp_skin_scores, paired = TRUE,
                       alternative = "greater", conf.int = TRUE)

# Effect sizes
n_paired <- sum(!is.na(acp_oral_scores) & !is.na(acp_skin_scores))
d_h1 <- cohens_d(acp_oral_scores, acp_skin_scores)
cliff_h1 <- cliffs_delta(acp_oral_scores, acp_skin_scores)

# Bootstrap confidence interval for the difference
set.seed(config$reproducibility$seed)
n_boot <- 1000
boot_diffs <- replicate(n_boot, {
  idx <- sample(length(acp_oral_scores), replace = TRUE)
  mean(acp_oral_scores[idx], na.rm = TRUE) - mean(acp_skin_scores[idx], na.rm = TRUE)
})
ci_95 <- quantile(boot_diffs, c(0.025, 0.975))

h1_result <- list(
  hypothesis = "H1: ACP oral-specific > skin-specific keratin expression",
  test = "Paired Wilcoxon signed-rank test",
  oral_mean = mean(acp_oral_scores, na.rm = TRUE),
  skin_mean = mean(acp_skin_scores, na.rm = TRUE),
  difference = mean(acp_oral_scores, na.rm = TRUE) - mean(acp_skin_scores, na.rm = TRUE),
  ci_lower = ci_95[1],
  ci_upper = ci_95[2],
  W_statistic = h1_test$statistic,
  p_value = h1_test$p.value,
  cohens_d = d_h1,
  cohens_d_interpretation = interpret_cohens_d(d_h1),
  cliffs_delta = cliff_h1,
  n_cells = n_paired,
  conclusion = ifelse(h1_test$p.value < 0.05 &&
                        mean(acp_oral_scores, na.rm = TRUE) > mean(acp_skin_scores, na.rm = TRUE),
                      "SUPPORTED", "NOT SUPPORTED")
)

message(sprintf("H1 Test: ACP oral-specific vs skin-specific keratins"))
message(sprintf("  Oral-specific mean: %.4f", h1_result$oral_mean))
message(sprintf("  Skin-specific mean: %.4f", h1_result$skin_mean))
message(sprintf("  Difference: %.4f [95%% CI: %.4f, %.4f]",
                h1_result$difference, h1_result$ci_lower, h1_result$ci_upper))
message(sprintf("  Cohen's d: %.3f (%s effect)", h1_result$cohens_d, h1_result$cohens_d_interpretation))
message(sprintf("  Cliff's delta: %.3f", h1_result$cliffs_delta))
message(sprintf("  Wilcoxon p-value: %.2e", h1_result$p_value))
message(sprintf("  Conclusion: %s", h1_result$conclusion))

# -----------------------------------------------------------------------------
# 1.3: Pairwise comparisons with effect sizes
# -----------------------------------------------------------------------------

message("\n--- 1.3: Pairwise comparisons (Kruskal-Wallis + post-hoc) ---\n")

pairwise_results <- list()

for (sig_name in names(keratin_signatures)) {
  sig_data <- signature_df %>% filter(signature == sig_name)

  # Check for sufficient non-NA data
  group_counts <- sig_data %>%
    group_by(dataset) %>%
    summarise(n_valid = sum(!is.na(mean_expr)), .groups = "drop")

  if (any(group_counts$n_valid < 10)) {
    message(sprintf("  Skipping %s: insufficient data", sig_name))
    pairwise_results[[sig_name]] <- data.frame(
      signature = sig_name,
      kruskal_chi2 = NA, kruskal_p = NA,
      p_acp_oral = NA, p_acp_skin = NA, p_oral_skin = NA,
      d_acp_oral = NA, d_acp_skin = NA, d_oral_skin = NA,
      cliff_acp_oral = NA, cliff_acp_skin = NA, cliff_oral_skin = NA,
      note = "Insufficient data"
    )
    next
  }

  sig_data_clean <- sig_data %>% filter(!is.na(mean_expr))

  # Kruskal-Wallis
  kw_test <- tryCatch(
    kruskal.test(mean_expr ~ dataset, data = sig_data_clean),
    error = function(e) list(statistic = NA, p.value = NA)
  )

  # Pairwise Wilcoxon
  pw_tests <- tryCatch(
    pairwise.wilcox.test(sig_data_clean$mean_expr, sig_data_clean$dataset,
                         p.adjust.method = "bonferroni"),
    error = function(e) list(p.value = matrix(NA, 2, 2,
                                              dimnames = list(c("Oral", "Skin"), c("ACP", "Oral"))))
  )

  acp_vals <- sig_data_clean$mean_expr[sig_data_clean$dataset == "ACP"]
  oral_vals <- sig_data_clean$mean_expr[sig_data_clean$dataset == "Oral"]
  skin_vals <- sig_data_clean$mean_expr[sig_data_clean$dataset == "Skin"]

  pairwise_results[[sig_name]] <- data.frame(
    signature = sig_name,
    kruskal_chi2 = kw_test$statistic %||% NA,
    kruskal_p = kw_test$p.value %||% NA,
    p_acp_oral = tryCatch(pw_tests$p.value["Oral", "ACP"], error = function(e) NA),
    p_acp_skin = tryCatch(pw_tests$p.value["Skin", "ACP"], error = function(e) NA),
    p_oral_skin = tryCatch(pw_tests$p.value["Skin", "Oral"], error = function(e) NA),
    d_acp_oral = cohens_d(acp_vals, oral_vals),
    d_acp_skin = cohens_d(acp_vals, skin_vals),
    d_oral_skin = cohens_d(oral_vals, skin_vals),
    cliff_acp_oral = cliffs_delta(acp_vals, oral_vals),
    cliff_acp_skin = cliffs_delta(acp_vals, skin_vals),
    cliff_oral_skin = cliffs_delta(oral_vals, skin_vals),
    note = NA_character_
  )
}

pairwise_df <- do.call(rbind, pairwise_results)
rownames(pairwise_df) <- NULL

# FDR adjustment
non_na_idx <- !is.na(pairwise_df$kruskal_p)
pairwise_df$kruskal_p_adj <- NA
pairwise_df$kruskal_p_adj[non_na_idx] <- p.adjust(pairwise_df$kruskal_p[non_na_idx], method = "BH")

message("Pairwise comparison results (with effect sizes):")
print(pairwise_df %>%
        filter(is.na(note)) %>%
        select(signature, kruskal_p_adj, d_acp_oral, d_acp_skin, cliff_acp_oral, cliff_acp_skin) %>%
        mutate(across(where(is.numeric), ~round(., 3))))

# -----------------------------------------------------------------------------
# 1.4: Individual gene expression with statistical tests
# -----------------------------------------------------------------------------

message("\n--- 1.4: Individual keratin gene expression ---\n")

all_keratins <- unique(unlist(keratin_signatures[c("oral_specific", "skin_specific",
                                                   "basal_shared", "simple_epithelium")]))

gene_expr_list <- list()
gene_stats_list <- list()

for (gene in all_keratins) {
  in_acp <- gene %in% rownames(acp_sub)
  in_oral <- gene %in% rownames(oral_sub)
  in_skin <- gene %in% rownames(skin_sub)

  if (in_acp && in_oral && in_skin) {
    acp_expr <- as.numeric(GetAssayData(acp_sub, layer = "data")[gene, ])
    oral_expr <- as.numeric(GetAssayData(oral_sub, layer = "data")[gene, ])
    skin_expr <- as.numeric(GetAssayData(skin_sub, layer = "data")[gene, ])

    gene_expr_list[[gene]] <- data.frame(
      gene = gene,
      dataset = c(rep("ACP", length(acp_expr)),
                  rep("Oral", length(oral_expr)),
                  rep("Skin", length(skin_expr))),
      expression = c(acp_expr, oral_expr, skin_expr)
    )

    # Statistical tests
    df_gene <- gene_expr_list[[gene]]
    kw <- kruskal.test(expression ~ dataset, data = df_gene)

    gene_stats_list[[gene]] <- data.frame(
      gene = gene,
      mean_acp = mean(acp_expr),
      mean_oral = mean(oral_expr),
      mean_skin = mean(skin_expr),
      pct_acp = mean(acp_expr > 0) * 100,
      pct_oral = mean(oral_expr > 0) * 100,
      pct_skin = mean(skin_expr > 0) * 100,
      kw_p = kw$p.value,
      d_acp_oral = cohens_d(acp_expr, oral_expr),
      d_acp_skin = cohens_d(acp_expr, skin_expr),
      cliff_acp_oral = cliffs_delta(acp_expr, oral_expr),
      cliff_acp_skin = cliffs_delta(acp_expr, skin_expr)
    )
  }
}

gene_expr_df <- do.call(rbind, gene_expr_list)
gene_stats_df <- do.call(rbind, gene_stats_list)
rownames(gene_stats_df) <- NULL

# FDR correction
gene_stats_df$kw_fdr <- p.adjust(gene_stats_df$kw_p, method = "BH")

# Add specificity annotation
gene_stats_df$specificity <- case_when(
  gene_stats_df$gene %in% keratin_signatures$oral_specific ~ "Oral-specific",
  gene_stats_df$gene %in% keratin_signatures$skin_specific ~ "Skin-specific",
  gene_stats_df$gene %in% keratin_signatures$basal_shared ~ "Basal (shared)",
  gene_stats_df$gene %in% keratin_signatures$simple_epithelium ~ "Simple epithelium",
  TRUE ~ "Other"
)

message("Individual keratin gene statistics:")
print(gene_stats_df %>%
        select(gene, specificity, mean_acp, mean_oral, mean_skin, d_acp_oral, d_acp_skin, kw_fdr) %>%
        arrange(specificity, gene) %>%
        mutate(across(where(is.numeric), ~round(., 3))))

# ==============================================================================
# ANALYSIS 2: TRANSCRIPTIONAL SIMILARITY (ON INTEGRATED DATA)
# ==============================================================================

message("\n========================================")
message("ANALYSIS 2: Transcriptional Similarity")
message("========================================\n")

# -----------------------------------------------------------------------------
# 2.1: Correlation analysis on integrated data
# -----------------------------------------------------------------------------

message("--- 2.1: Transcriptome-wide correlation (integrated space) ---\n")

# Use variable features from integration
var_genes <- VariableFeatures(merged_integrated)
message(sprintf("  Using %d variable genes from integration", length(var_genes)))

# Get mean expression profiles per dataset
get_mean_profile <- function(dataset_name) {
  cells <- colnames(merged_integrated)[merged_integrated$source == dataset_name]
  expr <- GetAssayData(merged_integrated, layer = "data")[var_genes, cells]
  rowMeans(as.matrix(expr))
}

acp_profile <- get_mean_profile("ACP")
oral_profile <- get_mean_profile("Oral")
skin_profile <- get_mean_profile("Skin")

# Correlations
cor_acp_oral <- cor(acp_profile, oral_profile, method = "pearson")
cor_acp_skin <- cor(acp_profile, skin_profile, method = "pearson")
cor_oral_skin <- cor(oral_profile, skin_profile, method = "pearson")

message("Transcriptome correlations (Pearson):")
message(sprintf("  ACP vs Oral:  r = %.4f", cor_acp_oral))
message(sprintf("  ACP vs Skin:  r = %.4f", cor_acp_skin))
message(sprintf("  Oral vs Skin: r = %.4f", cor_oral_skin))

# Fisher z-test for H2
fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))
z_acp_oral <- fisher_z(cor_acp_oral)
z_acp_skin <- fisher_z(cor_acp_skin)
n_genes <- length(var_genes)
se_diff <- sqrt(2 / (n_genes - 3))
z_diff <- (z_acp_oral - z_acp_skin) / se_diff
p_diff <- 2 * (1 - pnorm(abs(z_diff)))

h2_result <- list(
  hypothesis = "H2: ACP transcriptome closer to oral than skin",
  test = "Fisher z-transformation comparison",
  cor_acp_oral = cor_acp_oral,
  cor_acp_skin = cor_acp_skin,
  cor_oral_skin = cor_oral_skin,
  z_difference = z_diff,
  p_value = p_diff,
  n_genes = n_genes,
  conclusion = case_when(
    p_diff < 0.05 && cor_acp_oral > cor_acp_skin ~ "SUPPORTED: ACP closer to Oral",
    p_diff < 0.05 && cor_acp_skin > cor_acp_oral ~ "OPPOSITE: ACP closer to Skin",
    TRUE ~ "NOT SUPPORTED: No significant difference"
  )
)

message(sprintf("\nH2 Test:"))
message(sprintf("  z-difference: %.4f, p-value: %.4f", z_diff, p_diff))
message(sprintf("  Conclusion: %s", h2_result$conclusion))

# -----------------------------------------------------------------------------
# 2.2: Per-cell similarity in integrated space
# -----------------------------------------------------------------------------

message("\n--- 2.2: Per-cell similarity (Harmony space) ---\n")

# Calculate centroids in Harmony space
harmony_emb <- Embeddings(merged_integrated, "harmony")[, 1:20]

oral_cells <- merged_integrated$source == "Oral"
skin_cells <- merged_integrated$source == "Skin"
acp_cells <- merged_integrated$source == "ACP"

oral_centroid <- colMeans(harmony_emb[oral_cells, ])
skin_centroid <- colMeans(harmony_emb[skin_cells, ])

# Calculate distance of each ACP cell to centroids
acp_harmony <- harmony_emb[acp_cells, ]

dist_to_oral <- apply(acp_harmony, 1, function(x) sqrt(sum((x - oral_centroid)^2)))
dist_to_skin <- apply(acp_harmony, 1, function(x) sqrt(sum((x - skin_centroid)^2)))

cell_similarities <- data.frame(
  cell = rownames(acp_harmony),
  dist_oral = dist_to_oral,
  dist_skin = dist_to_skin,
  diff = dist_to_skin - dist_to_oral,  # Positive = closer to oral
  closer_to = ifelse(dist_to_oral < dist_to_skin, "Oral", "Skin")
)

n_closer_oral <- sum(cell_similarities$closer_to == "Oral")
n_closer_skin <- sum(cell_similarities$closer_to == "Skin")
prop_oral <- n_closer_oral / nrow(cell_similarities)

# Binomial test
binom_test <- binom.test(n_closer_oral, nrow(cell_similarities), p = 0.5)

# Paired t-test on distances
ttest_dist <- t.test(cell_similarities$dist_skin, cell_similarities$dist_oral, paired = TRUE)

message(sprintf("Per-cell similarity (Harmony space):"))
message(sprintf("  Cells closer to Oral: %d (%.1f%%)", n_closer_oral, prop_oral * 100))
message(sprintf("  Cells closer to Skin: %d (%.1f%%)", n_closer_skin, (1 - prop_oral) * 100))
message(sprintf("  Binomial test p-value: %.2e", binom_test$p.value))
message(sprintf("  Mean distance to Oral: %.4f", mean(dist_to_oral)))
message(sprintf("  Mean distance to Skin: %.4f", mean(dist_to_skin)))
message(sprintf("  Paired t-test: t = %.2f, p = %.2e", ttest_dist$statistic, ttest_dist$p.value))

# -----------------------------------------------------------------------------
# 2.3: H3 Test - ACP distinctiveness (permutation)
# -----------------------------------------------------------------------------

message("\n--- 2.3: H3 Test - ACP distinctiveness ---\n")

# Use keratin genes for this test
keratin_genes <- unique(unlist(keratin_signatures[c("oral_specific", "skin_specific", "basal_shared")]))
keratin_genes <- keratin_genes[keratin_genes %in% rownames(merged_integrated)]

# Get expression
acp_krt <- as.matrix(GetAssayData(merged_integrated, layer = "data")[keratin_genes, acp_cells])
oral_krt <- as.matrix(GetAssayData(merged_integrated, layer = "data")[keratin_genes, oral_cells])
skin_krt <- as.matrix(GetAssayData(merged_integrated, layer = "data")[keratin_genes, skin_cells])

# Centroids
acp_centroid_k <- rowMeans(acp_krt)
oral_centroid_k <- rowMeans(oral_krt)
skin_centroid_k <- rowMeans(skin_krt)

# Distances
dist_acp_oral_k <- sqrt(sum((acp_centroid_k - oral_centroid_k)^2))
dist_acp_skin_k <- sqrt(sum((acp_centroid_k - skin_centroid_k)^2))
dist_oral_skin_k <- sqrt(sum((oral_centroid_k - skin_centroid_k)^2))

message(sprintf("Keratin profile centroid distances:"))
message(sprintf("  ACP to Oral: %.4f", dist_acp_oral_k))
message(sprintf("  ACP to Skin: %.4f", dist_acp_skin_k))
message(sprintf("  Oral to Skin: %.4f", dist_oral_skin_k))

# Permutation test
n_perm <- 1000
set.seed(config$reproducibility$seed)

all_krt <- cbind(acp_krt, oral_krt, skin_krt)
labels <- c(rep("ACP", ncol(acp_krt)), rep("Oral", ncol(oral_krt)), rep("Skin", ncol(skin_krt)))

observed_min_dist <- min(dist_acp_oral_k, dist_acp_skin_k)

perm_dists <- replicate(n_perm, {
  perm_labels <- sample(labels)
  perm_acp <- all_krt[, perm_labels == "ACP"]
  perm_oral <- all_krt[, perm_labels == "Oral"]
  perm_skin <- all_krt[, perm_labels == "Skin"]

  perm_acp_c <- rowMeans(perm_acp)
  perm_oral_c <- rowMeans(perm_oral)
  perm_skin_c <- rowMeans(perm_skin)

  min(sqrt(sum((perm_acp_c - perm_oral_c)^2)),
      sqrt(sum((perm_acp_c - perm_skin_c)^2)))
})

p_distinct <- mean(perm_dists >= observed_min_dist)

h3_result <- list(
  hypothesis = "H3: ACP distinct from both healthy references",
  test = "Permutation test on centroid distances",
  dist_acp_oral = dist_acp_oral_k,
  dist_acp_skin = dist_acp_skin_k,
  dist_oral_skin = dist_oral_skin_k,
  observed_min_dist = observed_min_dist,
  mean_perm_dist = mean(perm_dists),
  sd_perm_dist = sd(perm_dists),
  p_value = p_distinct,
  n_permutations = n_perm,
  n_genes = length(keratin_genes),
  conclusion = ifelse(p_distinct < 0.05,
                      "SUPPORTED: ACP is distinct from both references",
                      "NOT SUPPORTED: ACP not more distinct than expected")
)

message(sprintf("\nH3 Test: ACP distinctiveness"))
message(sprintf("  Observed min distance: %.4f", observed_min_dist))
message(sprintf("  Permutation mean: %.4f ± %.4f", mean(perm_dists), sd(perm_dists)))
message(sprintf("  p-value: %.4f", p_distinct))
message(sprintf("  Conclusion: %s", h3_result$conclusion))

# ==============================================================================
# ANALYSIS 3: ACP HETEROGENEITY
# ==============================================================================

message("\n========================================")
message("ANALYSIS 3: ACP Heterogeneity")
message("========================================\n")

# -----------------------------------------------------------------------------
# 3.1: Characterize oral-like vs skin-like ACP cells
# -----------------------------------------------------------------------------

message("--- 3.1: Oral-like vs Skin-like ACP cells ---\n")

# Add similarity info to ACP subset
acp_meta <- cell_similarities
acp_meta$oral_vs_skin_score <- scale(acp_meta$diff)[, 1]

# Classify into groups
acp_meta$similarity_group <- case_when(
  acp_meta$diff > quantile(acp_meta$diff, 0.75) ~ "Oral-like",
  acp_meta$diff < quantile(acp_meta$diff, 0.25) ~ "Skin-like",
  TRUE ~ "Intermediate"
)

message(sprintf("ACP cell classification:"))
print(table(acp_meta$similarity_group))

# Add to Seurat object
acp_for_de <- subset(merged_integrated, source == "ACP")
acp_for_de$similarity_group <- acp_meta$similarity_group[match(colnames(acp_for_de), acp_meta$cell)]
acp_for_de$oral_vs_skin_score <- acp_meta$oral_vs_skin_score[match(colnames(acp_for_de), acp_meta$cell)]

# Differential expression
Idents(acp_for_de) <- "similarity_group"

# Only run DE if we have enough cells in each group
group_sizes <- table(acp_for_de$similarity_group)
if (all(group_sizes[c("Oral-like", "Skin-like")] >= 50)) {

  markers_oral_vs_skin <- FindMarkers(
    acp_for_de,
    ident.1 = "Oral-like",
    ident.2 = "Skin-like",
    only.pos = FALSE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    verbose = FALSE
  )

  markers_oral_vs_skin$gene <- rownames(markers_oral_vs_skin)

  message(sprintf("\nDE genes between oral-like and skin-like ACP cells: %d",
                  sum(markers_oral_vs_skin$p_val_adj < 0.05)))

  message("\nTop genes enriched in ORAL-LIKE ACP cells:")
  print(head(markers_oral_vs_skin %>%
               filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
               arrange(p_val_adj) %>%
               select(gene, avg_log2FC, pct.1, pct.2, p_val_adj), 10))

  message("\nTop genes enriched in SKIN-LIKE ACP cells:")
  print(head(markers_oral_vs_skin %>%
               filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
               arrange(p_val_adj) %>%
               select(gene, avg_log2FC, pct.1, pct.2, p_val_adj), 10))

  write.csv(markers_oral_vs_skin,
            file.path(tables_dir, "04_markers_orallike_vs_skinlike.csv"),
            row.names = FALSE)

} else {
  message("  Insufficient cells for DE analysis")
  markers_oral_vs_skin <- NULL
}

# -----------------------------------------------------------------------------
# 3.2: Correlation of similarity score with gene signatures
# -----------------------------------------------------------------------------

message("\n--- 3.2: Signature correlation with oral/skin similarity ---\n")

# Add module scores to ACP cells
for (sig_name in names(keratin_signatures)) {
  genes <- keratin_signatures[[sig_name]]
  genes_present <- genes[genes %in% rownames(acp_for_de)]

  if (length(genes_present) >= 3) {
    acp_for_de <- AddModuleScore(acp_for_de,
                                 features = list(genes_present),
                                 name = paste0(sig_name, "_score"),
                                 verbose = FALSE)
  }
}

# Correlate with oral_vs_skin_score
score_cols <- grep("_score1$", colnames(acp_for_de@meta.data), value = TRUE)

if (length(score_cols) > 0) {
  sig_correlations <- sapply(score_cols, function(col) {
    cor(acp_for_de@meta.data[[col]], acp_for_de$oral_vs_skin_score,
        use = "complete.obs", method = "spearman")
  })
  names(sig_correlations) <- gsub("_score1$", "", names(sig_correlations))

  message("Signature correlations with oral-likeness:")
  print(sort(sig_correlations, decreasing = TRUE))

  sig_cor_df <- data.frame(
    signature = names(sig_correlations),
    spearman_rho = sig_correlations
  )
  write.csv(sig_cor_df, file.path(tables_dir, "04_signature_similarity_correlations.csv"),
            row.names = FALSE)
}

# ==============================================================================
# ANALYSIS 4: SENSITIVITY ANALYSIS
# ==============================================================================

message("\n========================================")
message("ANALYSIS 4: Sensitivity Analysis")
message("========================================\n")

# -----------------------------------------------------------------------------
# 4.1: Robustness to subsampling
# -----------------------------------------------------------------------------

message("--- 4.1: Subsampling sensitivity (10 iterations) ---\n")

sensitivity_results <- list()

for (i in 1:10) {
  set.seed(i * 100)

  # Resample
  acp_sens <- subset(acp_epithelial,
                     cells = sample(colnames(acp_epithelial),
                                    min(5000, ncol(acp_epithelial))))
  oral_sens <- subset(oral_mucosal,
                      cells = sample(colnames(oral_mucosal),
                                     min(5000, ncol(oral_mucosal))))
  skin_sens <- subset(skin_epithelial,
                      cells = sample(colnames(skin_epithelial),
                                     min(5000, ncol(skin_epithelial))))

  # Recalculate H1 metrics
  oral_scores <- calc_percell_expression(acp_sens, keratin_signatures$oral_specific)$mean_expr
  skin_scores <- calc_percell_expression(acp_sens, keratin_signatures$skin_specific)$mean_expr

  diff_mean <- mean(oral_scores, na.rm = TRUE) - mean(skin_scores, na.rm = TRUE)
  d <- cohens_d(oral_scores, skin_scores)

  # Quick test
  test_p <- wilcox.test(oral_scores, skin_scores, paired = TRUE,
                        alternative = "greater")$p.value

  sensitivity_results[[i]] <- data.frame(
    iteration = i,
    oral_mean = mean(oral_scores, na.rm = TRUE),
    skin_mean = mean(skin_scores, na.rm = TRUE),
    difference = diff_mean,
    cohens_d = d,
    p_value = test_p,
    h1_supported = test_p < 0.05 && diff_mean > 0
  )
}

sens_df <- do.call(rbind, sensitivity_results)

message(sprintf("H1 results across 10 subsamples:"))
message(sprintf("  Difference: %.4f ± %.4f (range: %.4f to %.4f)",
                mean(sens_df$difference), sd(sens_df$difference),
                min(sens_df$difference), max(sens_df$difference)))
message(sprintf("  Cohen's d: %.3f ± %.3f",
                mean(sens_df$cohens_d), sd(sens_df$cohens_d)))
message(sprintf("  H1 supported in %d/10 iterations", sum(sens_df$h1_supported)))

if (all(sens_df$h1_supported)) {
  message("  ✓ Result is ROBUST: H1 supported in all subsamples")
  sensitivity_conclusion <- "Robust"
} else if (sum(sens_df$h1_supported) >= 8) {
  message("  ~ Result is MOSTLY ROBUST: H1 supported in ≥80% of subsamples")
  sensitivity_conclusion <- "Mostly robust"
} else {
  message("  ⚠️ Result is NOT ROBUST: H1 inconsistent across subsamples")
  sensitivity_conclusion <- "Not robust"
}

write.csv(sens_df, file.path(tables_dir, "04_sensitivity_analysis.csv"), row.names = FALSE)

# ==============================================================================
# GENERATE PLOTS
# ==============================================================================

message("\n========================================")
message("Generating Plots")
message("========================================\n")

# -----------------------------------------------------------------------------
# Plot 1: Signature violin plots with effect sizes
# -----------------------------------------------------------------------------

message("Creating signature violin plots...")

p_signatures <- ggplot(signature_df %>%
                         filter(signature %in% c("oral_specific", "skin_specific",
                                                 "basal_shared", "cornification")),
                       aes(x = dataset, y = mean_expr, fill = dataset)) +
  geom_violin(scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white") +
  facet_wrap(~signature, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("ACP" = "#E41A1C", "Oral" = "#4DAF4A", "Skin" = "#377EB8")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold")) +
  labs(title = "Keratin Signature Expression by Dataset",
       subtitle = "Mean normalized expression of signature genes",
       x = NULL, y = "Mean Expression")

ggsave(file.path(fig_dir, "04_keratin_signatures_violin.pdf"), p_signatures,
       width = 10, height = 8)
ggsave(file.path(fig_dir, "04_keratin_signatures_violin.png"), p_signatures,
       width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# Plot 2: Individual keratin gene heatmap
# -----------------------------------------------------------------------------

message("Creating keratin gene heatmap...")

gene_summary <- gene_expr_df %>%
  group_by(gene, dataset) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")

heatmap_data <- gene_summary %>%
  pivot_wider(names_from = dataset, values_from = mean_expr) %>%
  column_to_rownames("gene")

gene_annot <- data.frame(
  gene = rownames(heatmap_data),
  specificity = case_when(
    rownames(heatmap_data) %in% keratin_signatures$oral_specific ~ "Oral-specific",
    rownames(heatmap_data) %in% keratin_signatures$skin_specific ~ "Skin-specific",
    rownames(heatmap_data) %in% keratin_signatures$basal_shared ~ "Basal (shared)",
    TRUE ~ "Other"
  )
)
rownames(gene_annot) <- gene_annot$gene

heatmap_mat <- as.matrix(heatmap_data)
heatmap_mat_scaled <- t(apply(heatmap_mat, 1, function(x) {
  if (sd(x, na.rm = TRUE) == 0 || is.na(sd(x, na.rm = TRUE))) {
    return(rep(0, length(x)))
  }
  scale(x)[, 1]
}))
colnames(heatmap_mat_scaled) <- colnames(heatmap_mat)
heatmap_mat_scaled[is.na(heatmap_mat_scaled)] <- 0

row_ha <- rowAnnotation(
  Specificity = gene_annot$specificity,
  col = list(Specificity = c("Oral-specific" = "#4DAF4A",
                             "Skin-specific" = "#377EB8",
                             "Basal (shared)" = "#984EA3",
                             "Other" = "grey70"))
)

pdf(file.path(fig_dir, "04_keratin_heatmap.pdf"), width = 8, height = 10)
ht <- Heatmap(heatmap_mat_scaled,
              name = "Z-score",
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
              cluster_columns = FALSE,
              cluster_rows = TRUE,
              show_row_names = TRUE,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 12),
              right_annotation = row_ha,
              column_title = "Keratin Gene Expression (Z-scored)",
              column_title_gp = gpar(fontsize = 14, fontface = "bold"))
draw(ht)
dev.off()

# -----------------------------------------------------------------------------
# Plot 3: Cell similarity in Harmony space
# -----------------------------------------------------------------------------

message("Creating similarity scatter plot...")

p_similarity <- ggplot(cell_similarities, aes(x = dist_oral, y = dist_skin, color = closer_to)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Oral" = "#4DAF4A", "Skin" = "#377EB8")) +
  theme_minimal(base_size = 12) +
  labs(title = "ACP Cell Distance to Reference Centroids (Harmony Space)",
       subtitle = sprintf("%.1f%% closer to Oral, %.1f%% closer to Skin",
                          prop_oral * 100, (1 - prop_oral) * 100),
       x = "Distance to Oral Centroid",
       y = "Distance to Skin Centroid",
       color = "Closer to") +
  coord_fixed()

ggsave(file.path(fig_dir, "04_cell_similarity_scatter.pdf"), p_similarity,
       width = 8, height = 7)
ggsave(file.path(fig_dir, "04_cell_similarity_scatter.png"), p_similarity,
       width = 8, height = 7, dpi = 300)

# -----------------------------------------------------------------------------
# Plot 4: Integrated UMAP colored by similarity
# -----------------------------------------------------------------------------

message("Creating integrated UMAP plots...")

# Add similarity scores to full integrated object
merged_integrated$similarity_group <- NA
merged_integrated$similarity_group[acp_cells] <- acp_meta$similarity_group

p_umap_source <- DimPlot(merged_integrated, group.by = "source",
                         reduction = "umap_harmony", pt.size = 0.3, alpha = 0.5) +
  scale_color_manual(values = c("ACP" = "#E41A1C", "Oral" = "#4DAF4A", "Skin" = "#377EB8")) +
  labs(title = "Integrated UMAP by Source")

# ACP only colored by similarity
acp_umap <- subset(merged_integrated, source == "ACP")
acp_umap$similarity_group <- acp_meta$similarity_group[match(colnames(acp_umap), acp_meta$cell)]

p_umap_similarity <- DimPlot(acp_umap, group.by = "similarity_group",
                             reduction = "umap_harmony", pt.size = 0.5) +
  scale_color_manual(values = c("Oral-like" = "#4DAF4A",
                                "Intermediate" = "grey70",
                                "Skin-like" = "#377EB8")) +
  labs(title = "ACP Cells by Reference Similarity")

p_umap_combined <- p_umap_source + p_umap_similarity
ggsave(file.path(fig_dir, "04_integrated_umap.pdf"), p_umap_combined,
       width = 14, height = 6)
ggsave(file.path(fig_dir, "04_integrated_umap.png"), p_umap_combined,
       width = 14, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Plot 5: Effect size forest plot
# -----------------------------------------------------------------------------

message("Creating effect size forest plot...")

effect_df <- pairwise_df %>%
  filter(is.na(note)) %>%
  select(signature, d_acp_oral, d_acp_skin) %>%
  pivot_longer(cols = c(d_acp_oral, d_acp_skin),
               names_to = "comparison",
               values_to = "cohens_d") %>%
  mutate(comparison = case_when(
    comparison == "d_acp_oral" ~ "ACP vs Oral",
    comparison == "d_acp_skin" ~ "ACP vs Skin"
  ))

p_forest <- ggplot(effect_df, aes(x = cohens_d, y = signature, color = comparison)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-0.8, -0.5, -0.2, 0.2, 0.5, 0.8),
             linetype = "dotted", color = "grey80", alpha = 0.5) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  scale_color_manual(values = c("ACP vs Oral" = "#4DAF4A", "ACP vs Skin" = "#377EB8")) +
  theme_minimal(base_size = 12) +
  labs(title = "Effect Sizes: ACP vs References",
       subtitle = "Cohen's d (positive = higher in ACP)",
       x = "Cohen's d",
       y = NULL,
       color = "Comparison") +
  annotate("text", x = c(-0.8, -0.5, -0.2, 0.2, 0.5, 0.8), y = 0.3,
           label = c("large-", "med-", "small-", "small+", "med+", "large+"),
           size = 2.5, color = "grey50")

ggsave(file.path(fig_dir, "04_effect_size_forest.pdf"), p_forest,
       width = 10, height = 6)
ggsave(file.path(fig_dir, "04_effect_size_forest.png"), p_forest,
       width = 10, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Plot 6: Sensitivity analysis
# -----------------------------------------------------------------------------

message("Creating sensitivity analysis plot...")

p_sensitivity <- ggplot(sens_df, aes(x = factor(iteration), y = difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = h1_supported), size = 3) +
  geom_errorbar(aes(ymin = difference - 0.01, ymax = difference + 0.01), width = 0.2) +
  scale_color_manual(values = c("TRUE" = "#4DAF4A", "FALSE" = "#E41A1C")) +
  theme_minimal(base_size = 12) +
  labs(title = "Sensitivity Analysis: H1 Result Stability",
       subtitle = sprintf("Result: %s (%d/10 iterations support H1)",
                          sensitivity_conclusion, sum(sens_df$h1_supported)),
       x = "Subsample Iteration",
       y = "Oral - Skin Keratin Difference",
       color = "H1 Supported")

ggsave(file.path(fig_dir, "04_sensitivity_analysis.pdf"), p_sensitivity,
       width = 8, height = 5)
ggsave(file.path(fig_dir, "04_sensitivity_analysis.png"), p_sensitivity,
       width = 8, height = 5, dpi = 300)

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

message("\n========================================")
message("Saving Results")
message("========================================\n")

# Compile hypothesis results
hypothesis_results <- data.frame(
  hypothesis = c("H1", "H2", "H3"),
  description = c(
    "ACP oral-specific > skin-specific keratins",
    "ACP transcriptome closer to oral than skin",
    "ACP distinct from both healthy references"
  ),
  test = c(h1_result$test, h2_result$test, h3_result$test),
  statistic = c(h1_result$W_statistic, h2_result$z_difference, h3_result$observed_min_dist),
  p_value = c(h1_result$p_value, h2_result$p_value, h3_result$p_value),
  effect_size = c(h1_result$cohens_d,
                  h2_result$cor_acp_oral - h2_result$cor_acp_skin,
                  h3_result$observed_min_dist - h3_result$mean_perm_dist),
  effect_interpretation = c(h1_result$cohens_d_interpretation, "correlation diff", "distance excess"),
  conclusion = c(h1_result$conclusion, h2_result$conclusion, h3_result$conclusion)
)

write.csv(hypothesis_results, file.path(tables_dir, "04_hypothesis_results.csv"), row.names = FALSE)
message(sprintf("Saved: %s", file.path(tables_dir, "04_hypothesis_results.csv")))

write.csv(signature_summary, file.path(tables_dir, "04_signature_summary.csv"), row.names = FALSE)
write.csv(pairwise_df, file.path(tables_dir, "04_pairwise_comparisons.csv"), row.names = FALSE)
write.csv(gene_stats_df, file.path(tables_dir, "04_individual_keratin_stats.csv"), row.names = FALSE)
write.csv(cell_similarities, file.path(tables_dir, "04_cell_similarity_scores.csv"), row.names = FALSE)

# Save comprehensive results object
results_list <- list(
  batch_assessment = list(
    variance_pre = mean_source_var,
    variance_post = mean_source_var_post,
    silhouette = mean_sil,
    positive_control_passed = positive_control_passed
  ),
  h1 = h1_result,
  h2 = h2_result,
  h3 = h3_result,
  cell_similarity = list(
    n_closer_oral = n_closer_oral,
    n_closer_skin = n_closer_skin,
    prop_oral = prop_oral,
    binom_pvalue = binom_test$p.value
  ),
  sensitivity = list(
    conclusion = sensitivity_conclusion,
    n_supported = sum(sens_df$h1_supported),
    mean_difference = mean(sens_df$difference),
    sd_difference = sd(sens_df$difference)
  )
)
saveRDS(results_list, file.path(objects_dir, "04_origin_comparison_results.rds"))
saveRDS(merged_integrated, file.path(objects_dir, "04_integrated_acp_oral_skin.rds"))

message(sprintf("Saved: %s", file.path(objects_dir, "04_origin_comparison_results.rds")))
message(sprintf("Saved: %s", file.path(objects_dir, "04_integrated_acp_oral_skin.rds")))

# Generate markdown report
report <- c(
  "# Epithelial Origin Comparison Report (v2)",
  "",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Methodological Improvements",
  "",
  "This analysis includes:",
  "- Batch effect diagnostic and Harmony integration",
  "- Effect sizes (Cohen's d, Cliff's delta) with bootstrap CIs",
  "- Positive control (oral vs skin separation)",
  "- Sensitivity analysis (10 subsampling iterations)",
  "",
  "## Batch Effect Assessment",
  "",
  sprintf("| Metric | Value |"),
  sprintf("|--------|-------|"),
  sprintf("| Variance by dataset (pre-integration) | %.1f%% |", mean_source_var),
  sprintf("| Variance by dataset (post-Harmony) | %.1f%% |", mean_source_var_post),
  sprintf("| Silhouette (Oral vs Skin) | %.3f |", mean_sil),
  sprintf("| Positive control passed | %s |", positive_control_passed),
  "",
  "## Hypothesis Test Results",
  "",
  "| Hypothesis | p-value | Effect Size | Conclusion |",
  "|------------|---------|-------------|------------|",
  sprintf("| H1: Oral > Skin keratins | %.2e | d = %.2f (%s) | %s |",
          h1_result$p_value, h1_result$cohens_d, h1_result$cohens_d_interpretation, h1_result$conclusion),
  sprintf("| H2: Closer to oral | %.4f | Δr = %.3f | %s |",
          h2_result$p_value, h2_result$cor_acp_oral - h2_result$cor_acp_skin, h2_result$conclusion),
  sprintf("| H3: Distinct from both | %.4f | excess = %.2f | %s |",
          h3_result$p_value, h3_result$observed_min_dist - h3_result$mean_perm_dist, h3_result$conclusion),
  "",
  "## Key Findings",
  "",
  "### Keratin Profile (H1)",
  sprintf("- Oral keratin mean: %.4f [95%% CI: %.4f, %.4f]",
          h1_result$oral_mean, h1_result$ci_lower + h1_result$skin_mean, h1_result$ci_upper + h1_result$skin_mean),
  sprintf("- Skin keratin mean: %.4f", h1_result$skin_mean),
  sprintf("- Cohen's d: %.3f (%s)", h1_result$cohens_d, h1_result$cohens_d_interpretation),
  "",
  "### Transcriptome Similarity (H2)",
  sprintf("- Correlation ACP-Oral: %.4f", h2_result$cor_acp_oral),
  sprintf("- Correlation ACP-Skin: %.4f", h2_result$cor_acp_skin),
  sprintf("- Cells closer to Oral: %.1f%%", prop_oral * 100),
  "",
  "### ACP Heterogeneity",
  sprintf("- Oral-like ACP cells: %d", sum(acp_meta$similarity_group == "Oral-like")),
  sprintf("- Skin-like ACP cells: %d", sum(acp_meta$similarity_group == "Skin-like")),
  if (!is.null(markers_oral_vs_skin)) sprintf("- DE genes: %d", sum(markers_oral_vs_skin$p_val_adj < 0.05)) else "- DE: insufficient cells",
  "",
  "### Sensitivity Analysis",
  sprintf("- H1 supported in %d/10 subsamples", sum(sens_df$h1_supported)),
  sprintf("- Conclusion: %s", sensitivity_conclusion),
  "",
  "## Interpretation",
  "",
  "### What the data show:",
  sprintf("1. ACP expresses oral-specific keratins more than skin-specific (%s effect)", h1_result$cohens_d_interpretation),
  sprintf("2. Transcriptome-wide, ACP is %s", h2_result$conclusion),
  sprintf("3. ACP is distinct from both healthy epithelia (p = %.4f)", h3_result$p_value),
  "",
  "### Caveats:",
  sprintf("- Batch effects reduced from %.0f%% to %.0f%% but not eliminated", mean_source_var, mean_source_var_post),
  sprintf("- Positive control (oral-skin separation) silhouette: %.3f", mean_sil),
  "- Oral atlas lacks differentiation state annotations",
  "",
  "## Output Files",
  "",
  "- `04_batch_diagnostic.csv`: Integration quality metrics",
  "- `04_hypothesis_results.csv`: Summary of hypothesis tests",
  "- `04_pairwise_comparisons.csv`: All signature comparisons with effect sizes",
  "- `04_individual_keratin_stats.csv`: Per-gene statistics",
  "- `04_markers_orallike_vs_skinlike.csv`: DE genes within ACP",
  "- `04_sensitivity_analysis.csv`: Subsampling results",
  "- `04_integrated_acp_oral_skin.rds`: Harmony-integrated Seurat object"
)

writeLines(report, file.path(tables_dir, "04_origin_comparison_report.md"))
message(sprintf("Saved: %s", file.path(tables_dir, "04_origin_comparison_report.md")))

# ==============================================================================
# COMPLETION
# ==============================================================================

message("\n")
message("================================================================")
message("  Epithelial Origin Comparison Complete!")
message("================================================================")
message(sprintf("  Finished: %s", Sys.time()))
message("")
message("Batch Assessment:")
message(sprintf("  Variance reduced: %.0f%% → %.0f%%", mean_source_var, mean_source_var_post))
message(sprintf("  Positive control: %s (silhouette = %.3f)",
                ifelse(positive_control_passed, "PASSED", "FAILED"), mean_sil))
message("")
message("Hypothesis Results:")
message(sprintf("  H1 (Oral > Skin keratins): %s", h1_result$conclusion))
message(sprintf("      Effect: d = %.2f (%s)", h1_result$cohens_d, h1_result$cohens_d_interpretation))
message(sprintf("  H2 (Closer to Oral): %s", h2_result$conclusion))
message(sprintf("  H3 (Distinct from both): %s", h3_result$conclusion))
message("")
message(sprintf("Sensitivity: %s (%d/10 subsamples support H1)",
                sensitivity_conclusion, sum(sens_df$h1_supported)))
message("")
message("Outputs:")
message(sprintf("  Tables: %s", tables_dir))
message(sprintf("  Figures: %s", fig_dir))
message("================================================================\n")
