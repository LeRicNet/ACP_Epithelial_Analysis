# !! UNDER CONSTRUCTION !!

# ACP Epithelial Differentiation Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Analysis code for the manuscript: **"Epithelial Differentiation Trajectories in Adamantinomatous Craniopharyngioma Revealed by Spatial Transcriptomics"**

## Overview

This repository contains all code required to reproduce the epithelial subtype classification and trajectory analysis presented in our manuscript. The analysis pipeline:

1. **Downloads reference data** from CELLxGENE Census (Reynolds et al. 2021 Human Skin Atlas)
2. **Merges datasets** from multiple sources (ACP scRNA-seq, GSE215932 snRNA-seq)
3. **Identifies and filters epithelial cells** from mixed-population datasets
4. **Classifies epithelial cells** using module score-based approach (replacing binary thresholds)
5. **Validates classifications** against reference signatures with enhanced statistics (bootstrap CIs, effect sizes)
6. **Identifies Transit-Amplifying populations** via proliferation markers and cell cycle state
7. **Reconstructs differentiation trajectories** using multi-method consensus (Monocle3, Slingshot)
8. **Validates findings** across individual samples and independent datasets

## Quick Start

```bash
# Clone the repository
git clone https://github.com/lericnet/ACP_Epithelial_Analysis.git
cd ACP_Epithelial_Analysis

# Install dependencies using renv
Rscript -e "install.packages('renv'); renv::restore()"

# Run the full analysis pipeline
Rscript scripts/00_run_full_pipeline.R

# Or run with options
Rscript scripts/00_run_full_pipeline.R --dataset merged --skip-reference
Rscript scripts/00_run_full_pipeline.R --dry-run  # Preview steps without running

# Or run individual steps
Rscript scripts/00a_download_skin_reference.R --cell_type keratinocyte
Rscript scripts/00_merge_datasets.R
Rscript scripts/01_cell_type_annotation.R --dataset merged
Rscript scripts/01b_reference_validation_enhanced_stats.R --dataset merged
Rscript scripts/02_trajectory_analysis.R --dataset merged
```

## Project Structure

```
ACP_Epithelial_Analysis/
├── config/
│   └── config.yaml          # Central configuration (paths, signatures, parameters)
├── data/
│   ├── raw/                  # Original data files (not tracked in git)
│   ├── processed/            # Processed intermediates (merged datasets)
│   ├── external/             # External validation datasets (GSE215932)
│   └── reference/            # Reference atlas data (CELLxGENE Census)
├── R/
│   ├── utils/
│   │   ├── config.R          # Configuration loading utilities
│   │   └── cell_cycle_scoring.R  # Multi-method cell cycle scoring
│   ├── classification/
│   │   └── module_score_classification.R
│   ├── trajectory/
│   │   └── trajectory_analysis.R
│   └── figures/
│       └── figure_generators.R
├── scripts/
│   ├── 00_run_full_pipeline.R              # Master orchestration script
│   ├── 00a_download_skin_reference.R       # Download reference atlas
│   ├── 00_merge_datasets.R                 # Merge ACP + GSE datasets
│   ├── 01_cell_type_annotation.R           # Cell type classification
│   ├── 01b_reference_validation_enhanced_stats.R  # Validation & statistics
│   ├── 02_trajectory_analysis.R            # Trajectory inference
│   └── generate_figures.R                  # Figure regeneration script
├── results/
│   ├── tables/               # Summary tables (CSV)
│   ├── figures/
│   │   ├── main/            # Main manuscript figures
│   │   └── supplementary/   # Supplementary figures
│   └── objects/             # Intermediate R objects
├── docs/                     # Additional documentation
├── renv.lock                 # Package versions for reproducibility
└── README.md
```

## Configuration

All analysis parameters are centralized in `config/config.yaml`:

- **File paths**: Input data locations and output directories
- **Gene signatures**: Epithelial subtype markers, proliferation genes, etc.
- **Classification parameters**: Score thresholds, hybrid detection settings
- **Trajectory parameters**: Monocle3 settings, root selection method
- **Visualization**: Color palettes, figure dimensions

To modify the analysis, edit `config/config.yaml` rather than the R code.

---

## Script Documentation

### 00_run_full_pipeline.R

Master orchestration script that runs the complete analysis pipeline.

#### Usage

```bash
# Run full pipeline with defaults
Rscript scripts/00_run_full_pipeline.R

# Run on merged dataset, skip reference download
Rscript scripts/00_run_full_pipeline.R --dataset merged --skip-reference

# Preview steps without executing
Rscript scripts/00_run_full_pipeline.R --dry-run

# Run with batch correction
Rscript scripts/00_run_full_pipeline.R --integrate
```

#### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--dataset` | `merged` | Dataset to analyze: spatial, snrnaseq, acp_scn, merged |
| `--skip-reference` | FALSE | Skip reference data download |
| `--skip-merge` | FALSE | Skip dataset merging |
| `--skip-validation` | FALSE | Skip reference validation |
| `--skip-trajectory` | FALSE | Skip trajectory analysis |
| `--integrate` | FALSE | Use Seurat integration for batch correction |
| `--dry-run` | FALSE | Print steps without executing |

---

### 00a_download_skin_reference.R

Downloads human skin reference atlas from CELLxGENE Census for validation.

#### Usage

```bash
# Download keratinocytes only (recommended)
Rscript scripts/00a_download_skin_reference.R --cell_type keratinocyte

# Limit download size
Rscript scripts/00a_download_skin_reference.R --max_cells 50000
```

#### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--cell_type` | `all` | Filter: all, keratinocyte, epithelial |
| `--max_cells` | `100000` | Maximum cells to download |
| `--assay` | `10x` | Assay type: 10x, smartseq, all |
| `--healthy_only` | TRUE | Only normal/healthy tissue |

#### Outputs

| File | Description |
|------|-------------|
| `data/reference/skin_reference_census.rds` | Seurat object with reference cells |
| `data/reference/skin_reference_signatures.rds` | Gene signature lists |

---

### 00_merge_datasets.R

Merges ACP scRNA-seq and GSE215932 snRNA-seq datasets with full gene counts.

#### Usage

```bash
# Basic merge
Rscript scripts/00_merge_datasets.R

# With batch correction
Rscript scripts/00_merge_datasets.R --integrate
```

#### Outputs

| File | Description |
|------|-------------|
| `data/processed/merged_acp_gse.rds` | Merged Seurat object |
| `data/processed/merged_acp_gse_integrated.rds` | With batch correction (if --integrate) |

---

### 01_cell_type_annotation.R

Epithelial cell type annotation pipeline supporting multiple single-cell and spatial transcriptomics datasets.

#### Supported Datasets

| Dataset | Flag | Input Source | Epithelial Column | Match Type |
|---------|------|--------------|-------------------|------------|
| 10x Visium Spatial | `--dataset spatial` | `config$paths$spatial_object` | `cellstates` | `"Epithelial"` |
| snRNA-seq (GSE215932) | `--dataset snrnaseq` | `config$paths$snrnaseq_processed` | `is_epithelial` | `TRUE` (boolean) |
| ACP scRNA-seq | `--dataset acp_scn` | `config$paths$acp_scn_annotated` | `celltypes` | `"Epithelial"` |

#### Usage

```bash
# Default: spatial dataset
Rscript scripts/01_cell_type_annotation.R

# Explicit dataset selection
Rscript scripts/01_cell_type_annotation.R --dataset spatial
Rscript scripts/01_cell_type_annotation.R --dataset snrnaseq
Rscript scripts/01_cell_type_annotation.R --dataset acp_scn

# With custom config
Rscript scripts/01_cell_type_annotation.R --dataset snrnaseq --config path/to/config.yaml
```

#### Pipeline Steps

1. **Epithelial Cell Identification**
   - Filters to epithelial cells using dataset-specific metadata columns
   - Creates standardized `is_epithelial_cell` flag
   - Reports cell counts and original annotation distributions

2. **Module Score Calculation**
   - Scores cells for epithelial subtype signatures (Basal-like, Intermediate, Specialized, Transit-Amplifying)
   - Additional signatures: Whorl_Wnt, Senescence, SASP

3. **Cell Cycle Scoring**
   - Multi-method consensus approach (Seurat, module score, Cyclone)
   - Reports phase distribution and method agreement

4. **Subtype Classification**
   - Assigns cells to epithelial subtypes based on highest module score above threshold
   - Hybrid state detection (Basal-Intermediate, Intermediate-Specialized)

5. **Transit-Amplifying Identification**
   - Requires cycling phase (S/G2M) if `require_cycling: true`
   - Proliferation thresholds: 0.15 (moderate), 0.30 (high/override)

6. **Diagnostic Visualization**
   - Score distributions, classification comparisons, cell cycle plots

#### Outputs

All outputs include dataset type suffix (e.g., `_spatial`, `_snrnaseq`, `_acp_scn`):

| Output | Location | Description |
|--------|----------|-------------|
| Annotated Seurat object | `results/objects/01_seurat_annotated_{dataset}.rds` | Epithelial cells with classifications |
| Classification summary | `results/tables/01_classification_summary_{dataset}.csv` | Cell counts per subtype |
| Per-sample summary | `results/tables/01_classification_by_sample_{dataset}.csv` | Breakdown by sample |
| Gene signatures used | `results/tables/01_gene_signatures_used_{dataset}.csv` | Genes present per signature |
| Diagnostic plots | `results/figures/supplementary/01_cell_type_annotation_diagnostics_{dataset}.pdf` | QC visualizations |

#### Configuration Parameters

Key parameters in `config/config.yaml`:

```yaml
paths:
  spatial_object: "data/processed/spatial_object.rds"
  snrnaseq_processed: "data/external/GSE215932/GSE215932_snRNA_processed.rds"
  acp_scn_annotated: "data/raw/acp_scn_annotated.rds"

classification:
  min_score_threshold: 0.10
  transit_amplifying:
    proliferation_threshold: 0.15
    override_threshold: 0.30
    require_cycling: true
    override_on_high_proliferation: true

signatures:
  epithelial_subtypes:
    Basal_like: [KRT5, KRT14, ...]
    Intermediate: [KRT8, KRT18, ...]
    Specialized: [KRT17, KRT23, ...]
    Transit_Amplifying: [MKI67, TOP2A, ...]
```

#### Notes

**Gene Coverage**: Datasets with fewer genes (e.g., integrated objects with ~3,000 genes) may have incomplete signature coverage. The script reports gene overlap for each signature and skips signatures with insufficient genes.

**Threshold Sensitivity**: The default threshold (0.10) may be too stringent for some datasets. snRNA-seq data often has lower module scores than spatial data. Examine score distributions in diagnostic plots and adjust thresholds if unassigned rate is high.

**Preserved Metadata**: Original annotations are preserved (not overwritten). The full object (including non-epithelial cells) is available in memory as `seurat_obj_full` during execution.

---

### 01b_reference_validation_enhanced_stats.R

Reference-based validation with enhanced statistical measures.

#### Usage

```bash
# Signature-based validation (no reference object needed)
Rscript scripts/01b_reference_validation_enhanced_stats.R --dataset merged

# With downloaded reference
Rscript scripts/01b_reference_validation_enhanced_stats.R \
  --dataset merged \
  --reference custom \
  --ref_path data/reference/skin_reference_census.rds
```

#### Features

- **Bootstrap confidence intervals** (95% CIs for subtype proportions)
- **Effect sizes** (Cohen's d, Cliff's delta)
- **Silhouette scores** (cluster quality assessment)
- **Classification confidence** (score margins)
- **Reference concordance** (agreement with published signatures)

#### Outputs

| File | Description |
|------|-------------|
| `01b_bootstrap_ci_{dataset}.csv` | Proportions with 95% CIs |
| `01b_effect_sizes_{dataset}.csv` | Cohen's d and Cliff's delta |
| `01b_silhouette_scores_{dataset}.csv` | Cluster quality metrics |
| `01b_enhanced_validation_{dataset}.pdf` | 4-panel validation figure |

#### Interpretation Thresholds

| Metric | Good | Concerning |
|--------|------|------------|
| Effect size (Cohen's d) | > 0.8 | < 0.5 |
| Silhouette score | > 0.5 | < 0.25 |
| High-confidence cells | > 70% | < 50% |

---

### 02_trajectory_analysis.R

Multi-method trajectory inference with consensus pseudotime and comprehensive statistical evaluation.

#### Usage

```bash
Rscript scripts/02_trajectory_analysis.R
Rscript scripts/02_trajectory_analysis.R --cores 16
Rscript scripts/02_trajectory_analysis.R --no-bootstrap  # Skip bootstrap (faster)
```

#### Features

- **Multi-method consensus**: Monocle3, Slingshot, Diffusion pseudotime
- **Per-sample analysis**: Tests trajectory reproducibility across samples
- **Bootstrap confidence intervals**: Uncertainty quantification
- **Permutation testing**: Statistical significance of trajectory ordering
- **HPC parallelization**: Automatic SLURM detection

---

## Requirements

### R Version
- R >= 4.2.0

### Key Packages
- Seurat (>= 5.0.0)
- Monocle3
- SeuratWrappers
- slingshot
- cellxgene.census (for reference download)
- tiledbsoma (for reference download)
- cluster (for silhouette scores)
- effsize (for effect size calculations)
- UCell (optional, for module scoring)
- ggplot2
- patchwork
- dplyr/tidyr
- yaml

All package versions are locked via `renv.lock` for exact reproducibility.

### CELLxGENE Census Installation

The reference download script requires special packages from R-universe:

```r
install.packages('tiledbsoma', 
  repos = c('https://tiledb-inc.r-universe.dev', 'https://cloud.r-project.org'))
install.packages("cellxgene.census", 
  repos = c('https://chanzuckerberg.r-universe.dev', 'https://cloud.r-project.org'))
```

## Data Availability

### Input Data
- **Spatial transcriptomics**: 10x Visium data from 4 ACP tumor samples
- **snRNA-seq validation**: GSE215932 (Xu et al., 2024)
- **ACP scRNA-seq**: Annotated single-cell data

Due to size constraints, raw data is not included in this repository. Data can be obtained from:
- GEO: [GSE215932](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215932)
- [Link to spatial data repository]

After downloading, place files according to paths in `config/config.yaml`.

### Reference Data

Reference data for validation is downloaded automatically from CELLxGENE Census:

- **Source**: Reynolds et al. (2021) Human Skin Cell Atlas, *Science*
- **Download**: `Rscript scripts/00a_download_skin_reference.R --cell_type keratinocyte`
- **Output**: `data/reference/skin_reference_census.rds`

The script also saves gene signatures that can be used without the full reference download.

## Key Analysis Components

### 1. Module Score Classification
Replaces arbitrary binary thresholds with continuous module scores:
```r
# Calculate module scores for epithelial signatures
seurat_obj <- calculate_epithelial_scores(seurat_obj, config)

# Classify by highest score above threshold
seurat_obj <- classify_by_module_score(seurat_obj, config)
```

### 2. Transit-Amplifying Identification
Identifies proliferating progenitors using proliferation markers:
- Cells with high proliferation score AND in S/G2M phase
- Corresponds to "Unknown" population in original classification

### 3. Trajectory Analysis
Multi-method pseudotemporal ordering with consensus:
```r
# Run multi-method trajectory analysis
results <- run_multi_method_trajectory(seurat_obj, config)

# Per-sample validation with heterogeneity testing
per_sample_results <- run_per_sample_multi_trajectory(seurat_obj, config)
```

## Reproducing Manuscript Figures

| Figure | Description | Command |
|--------|-------------|---------|
| Fig 1 | Classification overview | `Rscript scripts/generate_figures.R 1` |
| Fig 2 | Transit-Amplifying characterization | `Rscript scripts/generate_figures.R 2` |
| Fig 3 | Trajectory analysis | `Rscript scripts/generate_figures.R 3` |
| Fig 4 | Gene expression dynamics | `Rscript scripts/generate_figures.R 4` |
| Supp | All supplementary figures | `Rscript scripts/generate_figures.R supp` |

## Citation

If you use this code, please cite:

```
[Author et al.], "Epithelial Differentiation Trajectories in Adamantinomatous 
Craniopharyngioma Revealed by Spatial Transcriptomics", [Journal], [Year].
```

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file.

## Contact

- **Corresponding Author**: [email]
- **Issues**: Please use GitHub Issues for bug reports and questions
