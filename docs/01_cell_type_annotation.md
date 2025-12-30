# 01_cell_type_annotation.R

Epithelial cell type annotation pipeline supporting multiple single-cell and spatial transcriptomics datasets.

## Overview

This script performs:
1. **Epithelial cell identification** - Filters to epithelial cells using dataset-specific metadata columns
2. **Module score calculation** - Scores cells for epithelial subtype signatures
3. **Cell cycle scoring** - Multi-method consensus approach (Seurat, module score, Cyclone)
4. **Subtype classification** - Assigns cells to epithelial subtypes including hybrid states
5. **Transit-Amplifying identification** - Identifies proliferating progenitors
6. **Diagnostic visualization** - Generates QC plots for classification assessment

## Supported Datasets

| Dataset | Flag | Input Source | Epithelial Column | Match Type |
|---------|------|--------------|-------------------|------------|
| 10x Visium Spatial | `--dataset spatial` | `config$paths$spatial_object` | `cellstates` | `"Epithelial"` |
| snRNA-seq (GSE215932) | `--dataset snrnaseq` | `config$paths$snrnaseq_processed` | `is_epithelial` | `TRUE` (boolean) |
| ACP scRNA-seq | `--dataset acp_scn` | `config$paths$acp_scn_annotated` | `celltypes` | `"Epithelial"` |

## Usage

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

## Pipeline Details

### Step 1: Epithelial Cell Identification

Each dataset uses a different metadata column to identify epithelial cells:

- **Spatial (10x Visium)**: Column `cellstates`, filter where value equals `"Epithelial"`
- **snRNA-seq (GSE215932)**: Column `is_epithelial`, filter where value is `TRUE`
- **ACP scRNA-seq**: Column `celltypes`, filter where value equals `"Epithelial"`

The script creates a standardized `is_epithelial_cell` flag and subsets to epithelial cells only. The full object (including non-epithelial cells) is preserved in memory as `seurat_obj_full`.

### Step 2: Module Score Calculation

Cells are scored against gene signatures defined in `config$signatures$epithelial_subtypes`:

| Signature | Description | Key Genes |
|-----------|-------------|-----------|
| Basal_like | Basal/stem-like markers | KRT5, KRT14, TP63, ITGA6 |
| Intermediate | Transitional state | KRT8, KRT18, KRT19, EPCAM |
| Specialized | Terminally differentiated | KRT17, KRT23, IVL, FLG |
| Transit_Amplifying | Proliferation + early differentiation | MKI67, TOP2A, PCNA, CDK1 |

Additional signatures scored: Whorl_Wnt, Senescence, SASP

### Step 3: Cell Cycle Scoring

Multi-method consensus using up to three approaches:
- Seurat `CellCycleScoring()`
- Module score thresholding (S.Score, G2M.Score)
- Cyclone classifier (scran)

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

## Outputs

All outputs include dataset type suffix (e.g., `_spatial`, `_snrnaseq`, `_acp_scn`):

| Output | Location | Description |
|--------|----------|-------------|
| Annotated Seurat object | `results/objects/01_seurat_annotated_{dataset}.rds` | Epithelial cells with classifications |
| Classification summary | `results/tables/01_classification_summary_{dataset}.csv` | Cell counts per subtype |
| Per-sample summary | `results/tables/01_classification_by_sample_{dataset}.csv` | Breakdown by sample |
| Gene signatures used | `results/tables/01_gene_signatures_used_{dataset}.csv` | Genes present per signature |
| Diagnostic plots | `results/figures/supplementary/01_cell_type_annotation_diagnostics_{dataset}.pdf` | QC visualizations |

## Configuration

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
    Basal_like: [KRT5, KRT14, KRT15, TP63, ...]
    Intermediate: [KRT8, KRT18, KRT19, KRT7, ...]
    Specialized: [KRT17, KRT23, KRT10, IVL, ...]
    Transit_Amplifying: [MKI67, TOP2A, STMN1, PCNA, ...]
```

## Example Output

```
=== Epithelial Cell Identification ===
  Using boolean column 'is_epithelial'
  Total cells: 59432
  Epithelial cells: 20477 (34.5%)
  Non-epithelial cells: 38955 (65.5%)

Module score summary statistics:
  Basal_like: mean=0.12, sd=0.45, range=[-0.89, 2.87]
  Intermediate: mean=0.08, sd=0.31, range=[-0.67, 1.36]
  Specialized: mean=-0.02, sd=0.28, range=[-0.54, 1.89]

Final classification:
        Basal-like: 4032 (19.7%)
  Transit-Amplifying: 539 (2.6%)
        Intermediate: 1352 (6.6%)
         Specialized: 1147 (5.6%)
  Basal-Intermediate: 1044 (5.1%)
          Unassigned: 11863 (57.9%)
```

## Troubleshooting

### High Unassigned Rate

If >50% of cells are unassigned:
1. Check score distributions in diagnostic plots
2. Lower `min_score_threshold` (try 0.05 or 0.02)
3. Verify gene coverage - integrated/filtered objects may have too few genes

### Gene Coverage Issues

The script reports gene overlap for each signature. If signatures are skipped:
- Check if your object has been heavily filtered
- Consider using a less-filtered version of the data
- Signatures with <3 genes present are automatically skipped

### Dataset-Specific Issues

**snRNA-seq**: Often has lower module scores than spatial data due to normalization differences. Consider dataset-specific thresholds.

**ACP scRNA-seq**: Integrated objects may have reduced gene counts (~3,000 genes). Check for gene coverage warnings.

## Adding New Datasets

To add support for a new dataset:

1. Add the input path to `config/config.yaml`:
   ```yaml
   paths:
     new_dataset: "data/path/to/data.rds"
   ```

2. Add configuration block in the script (around line 100):
   ```r
   } else if (dataset_type == "new_dataset") {
     input_path <- get_path(config, config$paths$new_dataset)
     count_col <- "nCount_RNA"
     feature_col <- "nFeature_RNA"
     default_assay <- "RNA"
     epi_column <- "your_epithelial_column"
     epi_values <- "Epithelial"  # or TRUE for boolean
     epi_is_boolean <- FALSE     # or TRUE
   }
   ```

3. Add to valid datasets list:
   ```r
   valid_datasets <- c("spatial", "snrnaseq", "acp_scn", "new_dataset")
   ```
