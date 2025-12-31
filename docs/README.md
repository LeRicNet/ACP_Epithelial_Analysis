# Documentation

Detailed documentation for the ACP Epithelial Differentiation Analysis pipeline.

## Scripts

| Script | Description |
|--------|-------------|
| [download_cellxgene_reference.md](download_cellxgene_reference.md) | Download skin keratinocyte reference from CELLxGENE Census |
| [00_merge_datasets.md](00_merge_datasets.md) | Merge ACP and GSE datasets, extract raw counts, harmonize metadata |
| [01_cell_type_annotation.md](01_cell_type_annotation.md) | Epithelial identification, module scoring, subtype classification, per-sample statistics |
| [01b_reference_validation.md](01b_reference_validation.md) | Reference-based validation, enhanced statistics, bootstrap CIs, effect sizes |
| [02_trajectory_analysis.md](02_trajectory_analysis.md) | Multi-method trajectory inference, consensus pseudotime, validation |

## Guides

| Guide | Description |
|-------|-------------|
| [configuration.md](configuration.md) | All `config.yaml` parameters explained |
| [troubleshooting.md](troubleshooting.md) | Common issues and solutions |
| [adding_datasets.md](adding_datasets.md) | How to add support for new datasets |

## Quick Reference

### Running the Pipeline

```bash
# Download reference data (optional - for validation)
Rscript scripts/download_cellxgene_reference.R           # Downloads + converts to Seurat
Rscript scripts/download_cellxgene_reference.R --force   # Force re-download

# Merge datasets (optional - for combined analysis)
Rscript scripts/00_merge_datasets.R              # Basic merge
Rscript scripts/00_merge_datasets.R --integrate  # With batch correction

# Individual dataset analysis
Rscript scripts/01_cell_type_annotation.R --dataset spatial
Rscript scripts/01_cell_type_annotation.R --dataset snrnaseq
Rscript scripts/01_cell_type_annotation.R --dataset acp_scn

# Combined dataset analysis
Rscript scripts/01_cell_type_annotation.R --dataset merged

# Validation and enhanced statistics
Rscript scripts/01b_reference_validation_enhanced_stats.R --dataset merged
Rscript scripts/01b_reference_validation_enhanced_stats.R --dataset merged --reference custom --ref_path data/reference/cellxgene_skin_keratinocytes.rds

# Trajectory analysis
Rscript scripts/02_trajectory_analysis.R --dataset merged --cores 16
```

### Supported Datasets

| Dataset | Flag | Description |
|---------|------|-------------|
| `spatial` | `--dataset spatial` | 10x Visium spatial transcriptomics |
| `snrnaseq` | `--dataset snrnaseq` | GSE215932 snRNA-seq |
| `acp_scn` | `--dataset acp_scn` | ACP scRNA-seq annotated |
| `merged` | `--dataset merged` | Combined ACP + GSE dataset |

### Key Outputs

| File | Location |
|------|----------|
| Reference data (h5ad) | `data/reference/cellxgene_skin_keratinocytes.h5ad` |
| Reference data (Seurat) | `data/reference/cellxgene_skin_keratinocytes.rds` |
| Merged dataset | `data/processed/merged_acp_gse.rds` |
| Classified cells | `results/objects/01_seurat_annotated_{dataset}.rds` |
| Validated cells | `results/objects/01b_seurat_validated_{dataset}.rds` |
| Pseudotime | `results/objects/02_seurat_with_pseudotime_{dataset}.rds` |
| Summary tables | `results/tables/` |
| Figures | `results/figures/main/` and `results/figures/supplementary/` |

### Per-Sample Analysis Outputs

| File | Description |
|------|-------------|
| `01_classification_by_sample_{dataset}.csv` | Cell counts by sample and subtype |
| `01_sample_proportions_wide_{dataset}.csv` | Subtype proportions per sample |
| `01_statistical_tests_{dataset}.csv` | Chi-square, Kruskal-Wallis, JS divergence |
| `01_sample_similarity_matrix_{dataset}.csv` | Pairwise sample similarity |
| `01_classification_by_dataset_merged.csv` | ACP vs GSE comparison (merged only) |

### Validation Outputs (01b script)

| File | Description |
|------|-------------|
| `01b_bootstrap_ci_{dataset}.csv` | Subtype proportions with 95% bootstrap CIs |
| `01b_effect_sizes_{dataset}.csv` | Cohen's d and Cliff's delta for each subtype |
| `01b_confidence_metrics_{dataset}.csv` | Classification confidence (score margins) |
| `01b_silhouette_scores_{dataset}.csv` | Cluster quality metrics |
| `01b_reference_concordance_{dataset}.csv` | Concordance with reference signatures |
| `01b_enhanced_validation_{dataset}.pdf` | Publication-quality validation figure |

### Configuration Quick Reference

```yaml
# Key paths
paths:
  acp_scn_annotated: "data/raw/acp_scn_annotated.rds"
  snrnaseq_processed: "data/external/GSE215932/GSE215932_snRNA_processed.rds"
  merged_dataset: "data/processed/merged_acp_gse.rds"
  
  # Reference data (created by download_cellxgene_reference.R)
  reference_dir: "data/reference"
  cellxgene_keratinocytes: "data/reference/cellxgene_skin_keratinocytes.h5ad"
  cellxgene_keratinocytes_rds: "data/reference/cellxgene_skin_keratinocytes.rds"

# Key classification parameters
classification:
  min_score_threshold: 0.10      # Minimum score for subtype assignment
  transit_amplifying:
    proliferation_threshold: 0.15  # TA identification threshold
    require_cycling: true          # Require S/G2M phase

# Key validation parameters
validation:
  n_bootstrap: 1000               # Bootstrap iterations for CIs
  reference_type: "signatures"    # signatures, custom, or none

# Key trajectory parameters  
trajectory:
  methods: [monocle3, slingshot]
  monocle3:
    root_subtype: "Basal-like"
```

## Typical Workflow

### Single Dataset Analysis
```bash
# Classify a single dataset
Rscript scripts/01_cell_type_annotation.R --dataset snrnaseq

# Run trajectory analysis
Rscript scripts/02_trajectory_analysis.R --dataset snrnaseq
```

### Combined Dataset Analysis
```bash
# Step 1: Merge ACP and GSE datasets
Rscript scripts/00_merge_datasets.R

# Step 2: Classify merged dataset
Rscript scripts/01_cell_type_annotation.R --dataset merged

# Step 3: Trajectory analysis on merged data
Rscript scripts/02_trajectory_analysis.R --dataset merged
```

### With Reference-Based Validation
```bash
# Step 1: Download reference data (one-time)
Rscript scripts/download_cellxgene_reference.R

# Step 2: Classify cells
Rscript scripts/01_cell_type_annotation.R --dataset merged

# Step 3: Validate against reference + enhanced statistics
Rscript scripts/01b_reference_validation_enhanced_stats.R \
  --dataset merged \
  --reference custom \
  --ref_path data/reference/cellxgene_skin_keratinocytes.rds
```

### With Batch Correction
```bash
# Merge with Seurat integration
Rscript scripts/00_merge_datasets.R --integrate

# Continue with classification
Rscript scripts/01_cell_type_annotation.R --dataset merged
```

## Statistical Analysis Features

The pipeline includes comprehensive per-sample statistical analysis:

| Analysis | Description |
|----------|-------------|
| **Chi-square test** | Tests subtype proportion independence across samples |
| **Kruskal-Wallis test** | Tests module score differences across samples |
| **Jensen-Shannon divergence** | Measures pairwise sample similarity |
| **Coefficient of Variation** | Identifies high-variability subtypes |
| **Dataset comparison** | ACP vs GSE breakdown (merged only) |

### Enhanced Validation (01b script)

| Metric | Description |
|--------|-------------|
| **Bootstrap CIs** | 95% confidence intervals for proportions |
| **Cohen's d** | Effect size for module score separation |
| **Cliff's delta** | Non-parametric effect size |
| **Silhouette scores** | Cluster quality assessment |
| **Score margin** | Classification confidence (top vs 2nd score) |
| **Reference concordance** | Agreement with published signatures |

## Troubleshooting Quick Tips

| Issue | Solution |
|-------|----------|
| High unassigned rate (>50%) | Lower `min_score_threshold` to 0.05 |
| Signature genes skipped | Use `00_merge_datasets.R` to extract full counts |
| Memory issues | Increase `future.globals.maxSize` |
| Cyclone warning | Expected with Seurat v5; uses 2-method consensus |
| Python environment issues | Run `setup_python_env()` or `setup_python_env(force_install = TRUE)` |
| Low effect sizes | Subtypes may overlap; check silhouette scores |
| Reference download fails | Check internet; run with `--force` to retry |
