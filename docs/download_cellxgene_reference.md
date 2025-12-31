# Download CELLxGENE Reference (download_cellxgene_reference.R)

Downloads human skin keratinocyte reference data from CELLxGENE Census and converts to Seurat format.

## Overview

This script provides an end-to-end solution for obtaining skin keratinocyte reference data:
1. Sets up a Python environment with required packages
2. Downloads data from CELLxGENE Census via Python
3. Converts h5ad to Seurat object
4. Caches results for fast subsequent access

## Reference Data

**Source**: [CELLxGENE Census](https://cellxgene.cziscience.com/) - a standardized collection of single-cell datasets.

**Cell types included**:
| Cell Type | Typical Count |
|-----------|---------------|
| keratinocyte | ~72,000 |
| basal cell of epidermis | ~29,000 |
| spinous cell of epidermis | ~26,000 |

**Total**: ~128,000 cells × 61,497 genes from healthy human skin.

## Requirements

### R Packages

```r
library(Seurat)
library(reticulate)
library(Matrix)
library(yaml)
```

### Python (managed automatically)

The script creates a Python virtualenv `acp_epi` with:
- `cellxgene-census` - CELLxGENE API
- `scanpy`, `anndata` - Data handling
- `numpy<2.4` - Numba compatibility
- `scipy`, `pandas`, `h5py` - Dependencies

### System Requirements

- Internet connection
- ~8 GB RAM for conversion
- ~4 GB disk space (2.1 GB h5ad + 1.5 GB RDS)

## File Structure

```
project/
├── scripts/
│   └── download_cellxgene_reference.R   # This script
├── python/
│   └── download_cellxgene_reference.py  # Python download script
├── config/
│   └── config.yaml                      # Path configuration
└── data/
    └── reference/
        ├── cellxgene_skin_keratinocytes.h5ad  # Downloaded data
        └── cellxgene_skin_keratinocytes.rds   # Seurat object
```

## Usage

### Interactive R Session (Recommended)

```r
source("scripts/download_cellxgene_reference.R")

# First time: sets up Python, downloads, converts, caches
skin_ref <- get_cellxgene_reference()

# Subsequent runs: loads from RDS cache (fast)
skin_ref <- get_cellxgene_reference()

# Force fresh download
skin_ref <- get_cellxgene_reference(force_download = TRUE)
```

### Command Line

```bash
# Download and convert (if not already done)
Rscript scripts/download_cellxgene_reference.R

# Force re-download
Rscript scripts/download_cellxgene_reference.R --force

# Show help
Rscript scripts/download_cellxgene_reference.R --help
```

### One-Time Python Setup

If you need to set up the Python environment manually:

```r
source("scripts/download_cellxgene_reference.R")
setup_python_env()  # Creates 'acp_epi' virtualenv

# Force reinstall packages
setup_python_env(force_install = TRUE)
```

## Functions

### `get_cellxgene_reference()`

**Main function** - handles everything automatically.

```r
get_cellxgene_reference(
  force_download = FALSE,  # Re-download h5ad
force_convert = FALSE,   # Re-convert to Seurat
  max_cells = 150000,      # Max cells to download
  save_rds = TRUE,         # Cache as RDS
  h5ad_path = NULL,        # Override h5ad path
  rds_path = NULL          # Override RDS path
)
```

**Returns**: Seurat object

**Behavior**:
1. If RDS exists → loads it (fastest)
2. If h5ad exists but no RDS → converts and saves
3. If nothing exists → downloads, converts, saves

### `download_cellxgene_reference()`

Downloads h5ad file only.

```r
download_cellxgene_reference(
  h5ad_path = NULL,   # Output path (default: from config)
  max_cells = 150000, # Max cells
  force = FALSE,      # Re-download if exists
  script_path = "../python/download_cellxgene_reference.py"
)
```

### `convert_h5ad_to_seurat()`
 
Converts any h5ad file to Seurat.

```r
convert_h5ad_to_seurat(
  h5ad_path,            # Input h5ad file
  assay_name = "RNA",   # Seurat assay name
  save_rds = FALSE,     # Save as RDS
  output_path = NULL,   # RDS output path
  envname = "acp_epi"   # Python env name
)
```

### `setup_python_env()`

Sets up Python environment.

```r
setup_python_env(
  envname = "acp_epi",    # Virtualenv name
  force_install = FALSE   # Reinstall packages
)
```

## Configuration

Add to `config/config.yaml`:

```yaml
paths:
  cellxgene_keratinocytes: "data/reference/cellxgene_skin_keratinocytes.h5ad"
  cellxgene_keratinocytes_rds: "data/reference/cellxgene_skin_keratinocytes.rds"
```

If not specified, defaults are used.

## Output

### Seurat Object Structure

```r
skin_ref
# An object of class Seurat
# 61497 features across 128072 samples within 1 assay
# Active assay: RNA (61497 features, 0 variable features)

# Metadata columns:
colnames(skin_ref@meta.data)
# [1] "cell_type"                  "cell_type_ontology_term_id"
# [3] "tissue"                     "tissue_general"
# [5] "disease"                    "assay"
# [7] "donor_id"                   "dataset_id"
# [9] "sex"                        "development_stage"

# Cell type distribution:
table(skin_ref$cell_type)
# basal cell of epidermis       keratinocyte spinous cell of epidermis
#                   29185              72566                     26321
```

## Using the Reference

### Label Transfer

```r
# Load reference and query
skin_ref <- get_cellxgene_reference()
query <- readRDS("results/objects/01_seurat_annotated_merged.rds")

# Normalize reference if needed
skin_ref <- NormalizeData(skin_ref)
skin_ref <- FindVariableFeatures(skin_ref)
skin_ref <- ScaleData(skin_ref)
skin_ref <- RunPCA(skin_ref)

# Find transfer anchors
anchors <- FindTransferAnchors(
  reference = skin_ref,
  query = query,
  dims = 1:30,
  reference.reduction = "pca"
)

# Transfer labels
predictions <- TransferData(
  anchorset = anchors,
  refdata = skin_ref$cell_type,
  dims = 1:30
)

query$predicted_celltype <- predictions$predicted.id
query$prediction_score <- predictions$prediction.score.max
```

### Signature Scoring

```r
# Define signatures from reference
basal_genes <- c("KRT5", "KRT14", "KRT15", "COL17A1", "ITGA6", "TP63")
spinous_genes <- c("KRT1", "KRT10", "DSC1", "DSG1")
differentiated_genes <- c("IVL", "LOR", "FLG", "SBSN")

signatures <- list(
  Basal = basal_genes,
  Spinous = spinous_genes,
  Differentiated = differentiated_genes
)

# Score cells
query <- AddModuleScore(query, features = signatures, name = "KC_")
```

## Troubleshooting

### Python Environment Issues

**Error**: `Python environment 'acp_epi' not found`
```r
setup_python_env()  # Creates the environment
```

**Error**: `scanpy not available`
```r
setup_python_env(force_install = TRUE)
```

**Error**: `RETICULATE_PYTHON conflicts`

The script automatically clears this, but if issues persist:
```r
Sys.unsetenv("RETICULATE_PYTHON")
```

Or remove from `.Renviron`:
```r
usethis::edit_r_environ("project")
# Delete/comment the RETICULATE_PYTHON line
```

### NumPy Compatibility

**Error**: `Numba needs NumPy 2.3 or less`

The script installs `numpy<2.4` automatically. If needed manually:
```r
reticulate::py_install("numpy<2.4", envname = "acp_epi", pip = TRUE)
```

### Download Issues

**Error**: `Download script not found`

Ensure the Python script is at `python/download_cellxgene_reference.py` relative to project root.

**Error**: `Download failed with exit code`

Check internet connection and try again. The CELLxGENE servers may be temporarily unavailable.

### Memory Issues

**Error**: `cannot allocate vector of size...`

The full dataset is large (~128K cells). Options:
1. Increase R memory limit
2. Download fewer cells: modify `max_cells` in Python script
3. Use a machine with more RAM (8+ GB recommended)

### Conversion Issues

**Error**: `File not found` after download

The h5ad file should be at `data/reference/cellxgene_skin_keratinocytes.h5ad`. Check:
```r
file.exists("data/reference/cellxgene_skin_keratinocytes.h5ad")
```

## Notes

- First-time setup takes 15-30 minutes (Python packages + download + conversion)
- Subsequent runs load from RDS cache in seconds
- The reference includes cells from multiple body sites (trunk, scalp, leg, etc.)
- All cells are from healthy/normal tissue (disease == 'normal')
- Data is from CELLxGENE Census stable release (updated periodically)
