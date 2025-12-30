# Configuration Guide

All analysis parameters are centralized in `config/config.yaml`. This document explains each section.

## Table of Contents

- [Project Metadata](#project-metadata)
- [File Paths](#file-paths)
- [Gene Signatures](#gene-signatures)
- [Classification Parameters](#classification-parameters)
- [Trajectory Analysis](#trajectory-analysis)
- [Cell Cycle Scoring](#cell-cycle-scoring)
- [Validation](#validation)
- [Visualization](#visualization)
- [Reproducibility](#reproducibility)
- [Output Formats](#output-formats)

---

## Project Metadata

```yaml
project:
  name: "ACP Epithelial Differentiation Analysis"
  version: "1.0.0"
  author: "Your Name"
  date_created: "2025-12-29"
  description: >
    Analysis of epithelial differentiation trajectories...
```

Metadata for documentation purposes. Not used by analysis code.

---

## File Paths

```yaml
paths:
  # Input data
  spatial_object: "data/raw/spatial-object.rds"
  spatial_object_with_subtypes: "data/raw/spatial-object-with-epithelial-subtypes.rds"
  snrnaseq_processed: "data/external/GSE215932/GSE215932_snRNA_processed.rds"
  acp_scn_annotated: "data/raw/acp_scn_annotated.rds"
  
  # Merged dataset (created by 00_merge_datasets.R)
  merged_dataset: "data/processed/merged_acp_gse.rds"

  # Output directories
  tables_dir: "results/tables"
  figures_main_dir: "results/figures/main"
  figures_supp_dir: "results/figures/supplementary"
  objects_dir: "results/objects"
```

| Path | Description |
|------|-------------|
| `spatial_object` | Primary 10x Visium Seurat object |
| `spatial_object_with_subtypes` | Visium object with existing classifications (for comparison) |
| `snrnaseq_processed` | Processed snRNA-seq validation data (GSE215932) |
| `acp_scn_annotated` | ACP scRNA-seq annotated object |
| `merged_dataset` | Combined ACP + GSE dataset (created by `00_merge_datasets.R`) |
| `*_dir` | Output directories (created automatically) |

**Note**: Paths are relative to project root. Use `get_path(config, config$paths$...)` in R code.

### Merged Dataset Notes

The merged dataset is created by `00_merge_datasets.R` which:
- Extracts **raw counts** from ACP (not variable features)
- Finds common genes between ACP and GSE
- Harmonizes metadata (sample, dataset, is_epithelial columns)
- Optionally integrates (batch correction) with `--integrate` flag

Two versions may exist:
- `merged_acp_gse.rds` - Basic merge (no batch correction)
- `merged_acp_gse_integrated.rds` - With Seurat integration

The script auto-detects which version exists.

---

## Gene Signatures

### Epithelial Subtypes

```yaml
signatures:
  epithelial_subtypes:
    Basal_like:
      - KRT5
      - KRT14
      - KRT15
      - TP63
      # ... more genes
    
    Intermediate:
      - KRT8
      - KRT18
      - KRT19
      # ...
    
    Specialized:
      - KRT17
      - KRT23
      - IVL
      # ...
    
    Transit_Amplifying:
      - MKI67
      - TOP2A
      - PCNA
      # ...
```

These signatures are used for module scoring. Each signature should have 8-15 genes for robust scoring.

**Minimum genes**: Signatures with <3 genes present in the dataset are skipped.

### Additional Signatures

```yaml
signatures:
  Whorl_Wnt:        # Wnt pathway activity
  Senescence:       # Senescence markers
  SASP:             # Senescence-associated secretory phenotype
  CD44_SPP1:        # CD44-SPP1 axis
  Collagen:         # ECM collagen genes
  Matricellular:    # Matricellular proteins
```

### Individual Markers

```yaml
markers:
  keratins: [KRT5, KRT8, KRT14, KRT17, KRT18, KRT19, KRT23]
  proliferation: [MKI67, TOP2A, PCNA, STMN1]
  whorl: [SHH, DKK4, CTNNB1, LEF1]
  senescence: [CDKN1A, CDKN2A, GLB1]
```

Used for visualization (violin plots, feature plots).

---

## Classification Parameters

```yaml
classification:
  method: "module_score"  # Options: "module_score", "binary_threshold", "hybrid"
  min_score_threshold: 0.1
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `"module_score"` | Classification approach |
| `min_score_threshold` | `0.1` | Minimum score to assign a subtype (cells below = "Unassigned") |

### Dataset-Specific Threshold Recommendations

| Dataset | Recommended Threshold | Notes |
|---------|----------------------|-------|
| `spatial` | 0.10 | Standard threshold works well |
| `acp_scn` | 0.10 | Works well with full gene counts |
| `snrnaseq` | 0.05 | Lower scores typical; reduce threshold |
| `merged` | 0.10 | Trade-off; GSE will have higher unassigned |

### Transit-Amplifying Identification

```yaml
classification:
  transit_amplifying:
    proliferation_threshold: 0.15
    override_threshold: 0.30
    require_cycling: true
    override_on_high_proliferation: true
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `proliferation_threshold` | `0.15` | Minimum proliferation score for TA (with cycling) |
| `override_threshold` | `0.30` | High proliferation overrides other classifications |
| `require_cycling` | `true` | Require S/G2M phase for TA assignment |
| `override_on_high_proliferation` | `true` | Allow high proliferation to override other subtypes |

### Per-Sample Analysis Settings

The script automatically performs per-sample statistical analysis:

| Analysis | Description |
|----------|-------------|
| Chi-square test | Tests subtype × sample independence |
| Kruskal-Wallis test | Tests module score differences across samples |
| Jensen-Shannon divergence | Pairwise sample similarity matrix |
| Coefficient of Variation | Flags subtypes with CV > 50% |

Sample column auto-detection order: `sample`, `Sample`, `sample_id`, `Sample_ID`, `orig.ident`, `patient`, `Patient`, `donor`, `Donor`

**Tuning tips**:
- High unassigned rate → lower `min_score_threshold` (try 0.05)
- Too many TA cells → increase `proliferation_threshold`
- snRNA-seq often needs lower thresholds than spatial data
- For merged datasets: check `01_dataset_comparison_merged.csv` to see ACP vs GSE differences

---

## Trajectory Analysis

```yaml
trajectory:
  methods:
    - monocle3
    - slingshot
    - diffusion
```

Available methods:
- `monocle3`: Graph-based trajectory (recommended)
- `slingshot`: Principal curves through clusters
- `diffusion`: Diffusion pseudotime

### Monocle3 Settings

```yaml
trajectory:
  monocle3:
    root_subtype: "Basal-like"
    root_selection_method: "all_cells"  # Options: "all_cells", "random_sample", "graph_node"
    n_root_cells: 50
    use_partition: true
    close_loop: false
    n_pcs: 30
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `root_subtype` | `"Basal-like"` | Subtype to use as trajectory start |
| `use_partition` | `true` | Learn separate trajectories per partition |
| `close_loop` | `false` | Allow circular trajectories |
| `n_pcs` | `30` | PCs for dimensionality reduction |

### Per-Sample Analysis

```yaml
trajectory:
  per_sample:
    enabled: true
    min_cells_per_sample: 100
    min_cells_per_subtype: 10
```

### Expected Ordering

```yaml
trajectory:
  expected_order:
    - "Basal-like"
    - "Transit-Amplifying"
    - "Intermediate"
    - "Specialized"
```

Used for validation - checks if observed pseudotime ordering matches expected biological progression.

### Statistical Settings

```yaml
trajectory:
  statistics:
    n_permutations: 1000
    alpha: 0.05
```

---

## Cell Cycle Scoring

```yaml
cell_cycle:
  methods:
    - seurat
    - module_score
    - cyclone
  min_genes: 5
  score_threshold: 0.1
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `methods` | `[seurat, module_score, cyclone]` | Scoring methods to use |
| `min_genes` | `5` | Minimum genes per phase required |
| `score_threshold` | `0.1` | Threshold for module_score method |

**Consensus determination**: Phase assigned by majority vote across methods. Confidence = proportion of methods agreeing.

---

## Validation

### snRNA-seq Cluster Mapping

```yaml
validation:
  snrnaseq:
    cluster_mapping:
      E_C1_PE: "Basal-like"
      E_C2_WC: "Specialized-Wnt"
      E_C3_KC: "Intermediate"
      E_C4_RHCG: "Specialized-Senescent"
      E_C5_Prolif: "Transit-Amplifying"
```

Maps GSE215932 cluster names to our nomenclature for validation comparisons.

---

## Visualization

### Color Palettes

```yaml
visualization:
  colors:
    epithelial_subtypes:
      Basal-like: "#E41A1C"
      Transit-Amplifying: "#FF7F00"
      Intermediate: "#4DAF4A"
      Specialized: "#377EB8"
      Basal-Intermediate: "#984EA3"
      Intermediate-Specialized: "#A65628"
      Unknown: "#999999"
    
    cell_cycle:
      G1: "#F8766D"
      S: "#00BA38"
      G2M: "#619CFF"
    
    samples:
      sampleA: "#1B9E77"
      sampleB: "#D95F02"
      # ...
```

### Figure Sizes

```yaml
visualization:
  figure_sizes:
    single_panel:
      width: 6
      height: 5
    double_panel:
      width: 12
      height: 5
    multi_panel:
      width: 14
      height: 10
```

### Plot Settings

```yaml
visualization:
  plot_settings:
    point_size: 1.5
    alpha: 0.7
    dpi: 300
    font_family: "Arial"
    base_size: 12
```

---

## Reproducibility

```yaml
reproducibility:
  seed: 42
  n_cores: 28
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `seed` | `42` | Random seed for all analyses |
| `n_cores` | `28` | Cores for parallel processing (auto-detected from SLURM if available) |

---

## Output Formats

```yaml
output:
  tables:
    format: "csv"  # Options: "csv", "tsv", "xlsx"
    include_excel_workbook: true

  figures:
    formats:
      - "pdf"
      - "png"
    png_dpi: 300

  objects:
    compress: true
```

---

## Accessing Config in R

```r
# Load config
config <- load_config()
config <- load_config("path/to/custom/config.yaml")

# Get full path
input_file <- get_path(config, config$paths$spatial_object)

# Get signature genes (filtered for available genes)
basal_genes <- get_signature(config, "Basal_like", rownames(seurat_obj))

# Get colors
colors <- get_colors(config, "epithelial_subtypes")

# Save with config settings
save_figure(plot, "my_figure", config, subdir = "main")
save_table(df, "my_table", config)
```

---

## Environment Variables

These environment variables override config settings:

| Variable | Overrides |
|----------|-----------|
| `RENV_PROFILE` | renv profile |
| `SLURM_CPUS_PER_TASK` | `reproducibility.n_cores` |

---

## Example: Adjusting for Different Datasets

### snRNA-seq (lower scores typical)

```yaml
classification:
  min_score_threshold: 0.05  # Lower threshold
  transit_amplifying:
    proliferation_threshold: 0.10
```

### Merged Dataset Analysis

```yaml
paths:
  merged_dataset: "data/processed/merged_acp_gse.rds"

classification:
  min_score_threshold: 0.10  # Balance between ACP and GSE
```

**Note**: For merged datasets, ACP samples typically classify better than GSE samples at the same threshold. Review `01_classification_by_dataset_merged.csv` to check for systematic differences.

### Highly filtered data (few genes)

```yaml
# Check gene coverage first, then adjust if needed
signatures:
  epithelial_subtypes:
    Basal_like:
      # Use only most robust markers
      - KRT5
      - KRT14
      - TP63
```

### Using Integrated Merged Data

If batch effects are a concern:

```bash
# Create integrated merged dataset
Rscript scripts/00_merge_datasets.R --integrate
```

This creates `merged_acp_gse_integrated.rds` with batch-corrected values. The script auto-detects and uses the integrated version if present.
