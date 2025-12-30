# 01_cell_type_annotation.R

Epithelial cell type annotation pipeline supporting multiple single-cell and spatial transcriptomics datasets.

## Overview

This script performs:
1. **Epithelial cell identification** - Filters to epithelial cells using dataset-specific metadata columns
2. **Module score calculation** - Scores cells for epithelial subtype signatures
3. **Cell cycle scoring** - Multi-method consensus approach (Seurat, module score, Cyclone)
4. **Subtype classification** - Assigns cells to epithelial subtypes including hybrid states
5. **Transit-Amplifying identification** - Identifies proliferating progenitors
6. **Per-sample statistical analysis** - Chi-square tests, Kruskal-Wallis tests, Jensen-Shannon divergence
7. **Diagnostic visualization** - Generates QC plots for classification assessment

## Supported Datasets

| Dataset | Flag | Input Source | Epithelial Column | Match Type |
|---------|------|--------------|-------------------|------------|
| 10x Visium Spatial | `--dataset spatial` | `config$paths$spatial_object` | `cellstates` | `"Epithelial"` |
| snRNA-seq (GSE215932) | `--dataset snrnaseq` | `config$paths$snrnaseq_processed` | `is_epithelial` | `TRUE` (boolean) |
| ACP scRNA-seq | `--dataset acp_scn` | `config$paths$acp_scn_annotated` | `celltypes` | `"Epithelial"` |
| **Merged ACP+GSE** | `--dataset merged` | `data/processed/merged_acp_gse.rds` | `is_epithelial` | `TRUE` (boolean) |

## Usage

```bash
# Default: spatial dataset
Rscript scripts/01_cell_type_annotation.R

# Explicit dataset selection
Rscript scripts/01_cell_type_annotation.R --dataset spatial
Rscript scripts/01_cell_type_annotation.R --dataset snrnaseq
Rscript scripts/01_cell_type_annotation.R --dataset acp_scn
Rscript scripts/01_cell_type_annotation.R --dataset merged

# With custom config
Rscript scripts/01_cell_type_annotation.R --dataset snrnaseq --config path/to/config.yaml
```

## Pipeline Details

### Step 1: Epithelial Cell Identification

Each dataset uses a different metadata column to identify epithelial cells:

| Dataset | Column | Filter |
|---------|--------|--------|
| Spatial (10x Visium) | `cellstates` | `== "Epithelial"` |
| snRNA-seq (GSE215932) | `is_epithelial` | `== TRUE` |
| ACP scRNA-seq | `celltypes` | `== "Epithelial"` |
| Merged | `is_epithelial` | `== TRUE` |

The script reports per-sample epithelial cell counts before filtering.

### Step 2: Module Score Calculation

Cells are scored against gene signatures defined in `config$signatures$epithelial_subtypes`:

| Signature | Description | Key Genes |
|-----------|-------------|-----------|
| Basal_like | Basal/stem-like markers | KRT5, KRT14, TP63, ITGA6 |
| Intermediate | Transitional state | KRT8, KRT18, KRT19, EPCAM |
| Specialized | Terminally differentiated | KRT17, KRT23, IVL, FLG |
| Transit_Amplifying | Proliferation + early differentiation | MKI67, TOP2A, PCNA, CDK1 |

Additional signatures: Whorl_Wnt, Senescence, SASP

**Gene Coverage Report**: The script reports signature gene coverage and warns if coverage is low (<80%) or critical (<50%). Signatures with <3 genes present are skipped.

### Step 3: Cell Cycle Scoring

Multi-method consensus using up to three approaches:
- **Seurat** `CellCycleScoring()` - Standard approach
- **Module score thresholding** - S.Score/G2M.Score with adaptive thresholds
- **Cyclone classifier** (scran) - When Ensembl ID mapping succeeds

Consensus phase is determined by majority vote with confidence scoring.

### Step 4: Subtype Classification

1. **Initial classification**: Highest-scoring subtype above threshold (default: 0.10)
2. **Hybrid state detection**: Cells scoring highly on adjacent subtypes
   - Basal-Intermediate: High in both Basal_like and Intermediate
   - Intermediate-Specialized: High in both Intermediate and Specialized

### Step 5: Transit-Amplifying Identification

Criteria for TA classification:
- Proliferation score ≥ 0.15 (moderate) or ≥ 0.30 (high/override)
- In S or G2M phase (if `require_cycling: true`)
- High proliferation can override other classifications
- Unassigned cycling cells can be rescued to TA

### Step 6: Per-Sample Statistical Analysis

Comprehensive statistical analysis across samples:

#### Sample Detection
Auto-detects sample column from: `sample`, `Sample`, `sample_id`, `Sample_ID`, `orig.ident`, `patient`, `Patient`, `donor`, `Donor`

#### Statistical Tests

| Test | Purpose | Interpretation |
|------|---------|----------------|
| **Chi-square** | Tests if subtype proportions are independent of sample | p < 0.05 = significant heterogeneity |
| **Kruskal-Wallis** | Tests if module scores differ across samples | p < 0.05 = significant score differences |
| **Jensen-Shannon divergence** | Measures similarity between sample distributions | < 0.1 = highly consistent, 0.1-0.3 = moderate, > 0.3 = heterogeneous |

#### Subtype Variability Metrics

For each subtype across samples:
- Mean percentage
- Standard deviation
- Min/Max range
- **Coefficient of Variation (CV)**: CV > 50% flagged as high variability

#### Dataset Comparison (Merged Only)

For merged datasets, additional analysis compares ACP vs GSE:
- Chi-square test for dataset × subtype independence
- Side-by-side proportion comparison
- Saved to dedicated output files

## Outputs

All outputs include dataset type suffix (e.g., `_spatial`, `_snrnaseq`, `_acp_scn`, `_merged`):

### Core Outputs

| Output | Location | Description |
|--------|----------|-------------|
| Annotated Seurat object | `results/objects/01_seurat_annotated_{dataset}.rds` | Epithelial cells with classifications |
| Classification summary | `results/tables/01_classification_summary_{dataset}.csv` | Cell counts per subtype |
| Gene signatures used | `results/tables/01_gene_signatures_used_{dataset}.csv` | Genes present per signature |
| Signature coverage | `results/tables/01_signature_coverage_{dataset}.csv` | Coverage % per signature |
| Diagnostic plots | `results/figures/supplementary/01_cell_type_annotation_diagnostics_{dataset}.pdf` | Score distributions, cell counts |

### Per-Sample Analysis Outputs

| Output | Description |
|--------|-------------|
| `01_classification_by_sample_{dataset}.csv` | Long format: sample, subtype, n_cells, pct_of_sample |
| `01_sample_proportions_wide_{dataset}.csv` | Wide format: samples as rows, subtypes as columns |
| `01_sample_score_statistics_{dataset}.csv` | Mean/SD of module scores per sample |
| `01_subtype_variability_{dataset}.csv` | CV, range, mean/SD for each subtype |
| `01_statistical_tests_{dataset}.csv` | Chi-square, Kruskal-Wallis results, JS divergence |
| `01_sample_similarity_matrix_{dataset}.csv` | Pairwise JS divergence matrix |
| `01_per_sample_diagnostics_{dataset}.pdf` | Per-sample visualization panels |

### Merged Dataset Additional Outputs

| Output | Description |
|--------|-------------|
| `01_classification_by_dataset_{dataset}.csv` | ACP vs GSE breakdown |
| `01_dataset_comparison_{dataset}.csv` | Side-by-side percentages |

## Diagnostic Plots

### Main Diagnostic Plot
5-panel figure showing:
- A: Score distributions by subtype (density plots)
- B: Mean scores per subtype (bar plot)
- C: Cell counts per subtype (bar plot)
- D: Proliferation score vs classification
- E: Cell cycle phase distribution

### Per-Sample Diagnostic Plot
5-panel figure showing:
- A: Stacked bar plot of subtype composition by sample
- B: Grouped bar plot comparing subtype proportions
- C: Violin plots of module scores by sample
- D: Heatmap of subtype proportions with percentages
- E: Cell cycle distribution by sample

## Configuration

Key parameters in `config/config.yaml`:

```yaml
paths:
  spatial_object: "data/processed/spatial_object.rds"
  snrnaseq_processed: "data/external/GSE215932/GSE215932_snRNA_processed.rds"
  acp_scn_annotated: "data/raw/acp_scn_annotated.rds"
  merged_dataset: "data/processed/merged_acp_gse.rds"  # Optional

classification:
  min_score_threshold: 0.10
  transit_amplifying:
    proliferation_threshold: 0.15
    override_threshold: 0.30
    require_cycling: true
    override_on_high_proliferation: true

signatures:
  epithelial_subtypes:
    Basal_like: [KRT5, KRT14, KRT15, TP63, ...]
    Intermediate: [KRT8, KRT18, KRT19, KRT7, ...]
    Specialized: [KRT17, KRT23, KRT10, IVL, ...]
    Transit_Amplifying: [MKI67, TOP2A, STMN1, PCNA, ...]
```

## Example Output

```
================================================================
  Step 1: Cell Type Annotation
  Dataset type: MERGED
================================================================

=== Gene Signature Coverage Report ===
Total unique genes in signatures: 42
Genes present in data: 40 (95.2%)

Final classification:
        Basal-like: 5199 (22.3%)
  Transit-Amplifying: 550 (2.4%)
        Intermediate: 1995 (8.5%)
         Specialized: 1445 (6.2%)
  Basal-Intermediate: 1699 (7.3%)
Intermediate-Specialized: 727 (3.1%)
          Unassigned: 11741 (50.3%)

--- Dataset Source Breakdown (ACP vs GSE) ---
           Subtype
Dataset     Basal-like Specialized Unassigned
  ACP_SCN        43.3%       14.6%      15.6%
  GSE215932      19.3%        5.0%      55.1%

Chi-square test (dataset × subtype): X² = 1956.57, p = 0.00
  → Significant differences between ACP and GSE classification patterns

--- Sample Consistency Metrics ---
Mean pairwise JS divergence: 0.1952
  → Samples show moderate consistency

================================================================
  Cell Type Annotation Complete!
================================================================
  Total cells: 23356
  Subtypes identified: 7
  Samples analyzed: 19

  Dataset breakdown:
    ACP_SCN: 2879 cells
    GSE215932: 20477 cells
```

## Interpreting Results

### High Unassigned Rate

| Unassigned % | Interpretation | Action |
|--------------|----------------|--------|
| < 20% | Good classification | None needed |
| 20-40% | Moderate | Check score distributions |
| 40-60% | High | Consider lowering threshold |
| > 60% | Very high | Lower threshold to 0.05 or investigate data quality |

**Dataset-specific considerations:**
- **snRNA-seq/GSE**: Often has lower scores due to normalization; may need threshold 0.05
- **ACP scRNA-seq**: Generally works well with 0.10 threshold
- **Merged**: Trade-off between datasets; ACP may be over-classified if threshold lowered

### Jensen-Shannon Divergence Interpretation

| Mean JS | Interpretation |
|---------|----------------|
| < 0.10 | Highly consistent across samples |
| 0.10 - 0.20 | Good consistency |
| 0.20 - 0.30 | Moderate heterogeneity |
| > 0.30 | Substantial sample-to-sample variation |

### High CV Subtypes

Coefficient of Variation > 50% indicates:
- Biological heterogeneity (real differences between samples)
- Technical artifacts (batch effects)
- Small sample sizes in some groups

## Troubleshooting

### High Unassigned Rate

1. Check score distributions in diagnostic plots
2. Lower `min_score_threshold` (try 0.05 or 0.02)
3. Verify gene coverage - integrated/filtered objects may have too few genes

### Gene Coverage Issues

The script reports gene overlap for each signature. If signatures are skipped:
- Check if your object has been heavily filtered
- Consider using `00_merge_datasets.R` to recover full gene counts
- Signatures with <3 genes present are automatically skipped

### Dataset-Specific Issues

**snRNA-seq**: Often has lower module scores than spatial data. Consider dataset-specific thresholds.

**ACP scRNA-seq**: Integrated objects may have reduced gene counts (~3,000 genes). Use `00_merge_datasets.R` to extract full RNA counts.

**Merged**: Significant differences between ACP and GSE patterns are expected. The chi-square test quantifies this.

### Cyclone Warning

```
Cyclone method failed: Please use the `layer` argument instead.
```

This is a Seurat v5 compatibility issue with scran. The script falls back to Seurat + module score consensus (2 methods instead of 3). Results are still valid.

## Adding New Datasets

To add support for a new dataset:

1. Add the input path to `config/config.yaml`:
   ```yaml
   paths:
     new_dataset: "data/path/to/data.rds"
   ```

2. Add to valid datasets list (line ~49):
   ```r
   valid_datasets <- c("spatial", "snrnaseq", "acp_scn", "merged", "new_dataset")
   ```

3. Add configuration block (after line ~143):
   ```r
   } else if (dataset_type == "new_dataset") {
     input_path <- get_path(config, config$paths$new_dataset)
     alt_path <- NULL
     count_col <- "nCount_RNA"
     feature_col <- "nFeature_RNA"
     default_assay <- "RNA"
     epi_column <- "your_epithelial_column"
     epi_values <- "Epithelial"  # or TRUE for boolean
     epi_is_boolean <- FALSE     # or TRUE
     check_full_data <- FALSE    # TRUE if object may be filtered
     message("Configured for new dataset")
   }
   ```

## Workflow Examples

### Standard Single Dataset
```bash
Rscript scripts/01_cell_type_annotation.R --dataset snrnaseq
```

### Merged Dataset Analysis
```bash
# First merge the datasets
Rscript scripts/00_merge_datasets.R

# Then classify
Rscript scripts/01_cell_type_annotation.R --dataset merged
```

### With Batch Correction
```bash
# Merge with integration
Rscript scripts/00_merge_datasets.R --integrate

# Classify integrated data
Rscript scripts/01_cell_type_annotation.R --dataset merged
```
