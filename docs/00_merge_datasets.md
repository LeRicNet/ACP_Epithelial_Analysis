# 00_merge_datasets.R

Preprocessing script to merge multiple single-cell datasets for combined analysis.

## Overview

This script merges the ACP scRNA-seq and GSE215932 snRNA-seq datasets into a single Seurat object for unified classification and trajectory analysis. Key features:

1. **Raw counts extraction** - Recovers full gene counts from ACP object (not variable features)
2. **Metadata harmonization** - Standardizes sample, epithelial, and dataset columns
3. **Common gene identification** - Finds intersection of genes across datasets
4. **Optional integration** - Batch correction using Seurat's integration workflow

## Usage

```bash
# Basic merge (no batch correction)
Rscript scripts/00_merge_datasets.R

# With batch correction (Seurat integration)
Rscript scripts/00_merge_datasets.R --integrate

# With custom config
Rscript scripts/00_merge_datasets.R --config path/to/config.yaml
```

## Input Requirements

| Dataset | Default Path | Required Columns |
|---------|--------------|------------------|
| ACP_SCN | `data/raw/acp_scn_annotated.rds` | `sample_id`, `celltypes` |
| GSE215932 | `data/external/GSE215932/GSE215932_snRNA_processed.rds` | `sample`, `is_epithelial` |

## Pipeline Details

### Step 1: Load Datasets

Both Seurat objects are loaded and updated to the current Seurat version (v5 compatible).

### Step 2: Extract Raw Counts

The script intelligently extracts the **full gene count matrix** from each dataset:

**ACP_SCN extraction priority:**
1. RNA assay → counts layer (preferred)
2. RNA assay → data layer (log-normalized fallback)
3. SCT assay → counts layer
4. SCT assay → data layer
5. Any other assay with more genes

This is critical because integrated ACP objects often have their default assay set to `integrated` with only ~3,000 variable features, while the RNA assay contains the full ~24,000 genes.

**GSE215932 extraction:**
- RNA assay → counts layer (Seurat v5 Assay5 format)

### Step 3: Harmonize Metadata

Metadata columns are standardized across datasets:

| Standardized Column | ACP Source | GSE Source |
|---------------------|------------|------------|
| `sample` | `sample_id` | `sample` |
| `dataset` | Added: `"ACP_SCN"` | Added: `"GSE215932"` |
| `dataset_source` | Added: `"ACP scRNA-seq"` | Added: `"GSE215932 snRNA-seq"` |
| `is_epithelial` | `celltypes == "Epithelial"` | `is_epithelial` (existing) |
| `original_celltype` | `celltypes` | `cell_type` |

### Step 4: Find Common Genes

The script identifies genes present in both datasets:

```
ACP genes: 23,841
GSE genes: 30,828
Common genes: 23,459
Signature coverage: 40/42 (95.2%)
```

Only common genes are retained in the merged object to ensure consistent scoring.

### Step 5: Create Merged Object

New Seurat objects are created from the common-gene count matrices and merged:

```r
merged_obj <- merge(acp_new, gse_new, add.cell.ids = c("ACP", "GSE"))
```

Cell barcodes are prefixed with dataset identifiers to prevent collisions.

### Step 6: Normalize Data

Standard Seurat normalization pipeline:
- `NormalizeData()` - Log-normalization
- `FindVariableFeatures()` - 3,000 variable features
- `ScaleData()` - Z-score scaling

### Step 7: Dimensionality Reduction

- PCA (50 components)
- UMAP (30 dimensions)

### Optional: Integration (Batch Correction)

With `--integrate` flag:

```r
anchors <- FindIntegrationAnchors(object.list, dims = 1:30)
merged_obj <- IntegrateData(anchorset = anchors, dims = 1:30)
```

This creates an `integrated` assay with batch-corrected values.

## Outputs

| Output | Location | Description |
|--------|----------|-------------|
| Merged Seurat object | `data/processed/merged_acp_gse.rds` | Combined dataset |
| Integrated object | `data/processed/merged_acp_gse_integrated.rds` | With batch correction (if `--integrate`) |
| Summary CSV | `data/processed/merged_dataset_summary.csv` | Merge statistics |

### Summary CSV Contents

```csv
metric,value
Total cells,78332
Total genes,23459
ACP cells,18900
GSE cells,59432
ACP epithelial,2879
GSE epithelial,20477
Common genes,23459
Signature gene coverage,40/42 (95.2%)
```

## Configuration

Add to `config/config.yaml`:

```yaml
paths:
  # Input datasets
  acp_scn_annotated: "data/raw/acp_scn_annotated.rds"
  snrnaseq_processed: "data/external/GSE215932/GSE215932_snRNA_processed.rds"
  
  # Merged output (optional - script will find automatically)
  merged_dataset: "data/processed/merged_acp_gse.rds"
```

## Workflow

```bash
# Step 1: Merge datasets
Rscript scripts/00_merge_datasets.R

# Step 2: Run classification on merged data
Rscript scripts/01_cell_type_annotation.R --dataset merged

# Step 3: Continue with trajectory analysis
Rscript scripts/02_trajectory_analysis.R --dataset merged
```

## Technical Notes

### Seurat v5 API Compatibility

The script uses `LayerData()` instead of the deprecated `GetAssayData(slot=...)` for Seurat v5 compatibility:

```r
# Seurat v5 way (used in this script)
counts <- LayerData(obj, assay = "RNA", layer = "counts")

# Deprecated Seurat v4 way (not used)
counts <- GetAssayData(obj, slot = "counts", assay = "RNA")
```

### Memory Considerations

The merged object is large (~78K cells × 23K genes). Expect:
- ~5-10 GB RAM during processing
- ~2-3 GB saved object size
- ~7 minutes runtime on typical hardware

### Why Extract Raw Counts?

Integrated Seurat objects often have:
- Default assay = `integrated` (3,000 variable features only)
- SCT assay with SCTransform-corrected counts
- RNA assay with **full original counts** (what we want)

Using variable features limits signature gene coverage to ~50%, while extracting full RNA counts provides ~95% coverage.

## Troubleshooting

### "Could not extract raw counts" Error

The script will show diagnostic output listing all assays and their layers. Common causes:
- Object was saved without counts (data only)
- Non-standard assay structure

**Solution**: Check your input object structure:
```r
obj <- readRDS("your_object.rds")
Assays(obj)  # List assays
Layers(obj[["RNA"]])  # List layers in RNA assay
```

### Low Common Gene Count

If common genes < 10,000:
- Check if one dataset was heavily filtered
- Verify both use the same gene annotation (symbols vs Ensembl)

### Memory Issues

For very large datasets:
```r
# Increase memory limit (if needed)
options(future.globals.maxSize = 8000 * 1024^2)  # 8 GB
```

## Example Output

```
================================================================
  Dataset Merge: ACP_SCN + GSE215932
================================================================
  Started: 2025-12-30 19:31:18
  Seurat version: 5.4.0

========================================
Extracting Raw Counts from ACP_SCN
========================================

Available assays in ACP_SCN:
  RNA (Assay): 23841 features
  SCT (Assay): 21973 features
  integrated (Assay): 3000 features

  RNA assay class: Assay
    Available layers: counts, data
    Counts layer dimensions: 23841 x 18900
    Successfully extracted 23841 genes from counts layer

========================================
Finding Common Genes
========================================

ACP genes: 23841
GSE genes: 30828
Common genes: 23459
Signature genes in common set: 40/42 (95.2%)

========================================
Creating Merged Seurat Object
========================================

Merged object:
  Total cells: 78332
  Genes: 23459
  ACP cells: 18900
  GSE cells: 59432

================================================================
  Dataset Merge Complete!
================================================================
  Total cells: 78332
  Total genes: 23459

To run classification on merged data:
  Rscript scripts/01_cell_type_annotation.R --dataset merged
================================================================
```
