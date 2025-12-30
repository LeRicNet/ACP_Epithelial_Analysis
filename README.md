# !! UNDER CONSTRUCTION !!

# ACP Epithelial Differentiation Analysis

Analysis code for epithelial differentiation trajectories in adamantinomatous craniopharyngioma using spatial transcriptomics and snRNA-seq validation.

## Quick Start

```bash
# Clone and setup
git clone https://github.com/lericnet/ACP_Epithelial_Analysis.git
cd ACP_Epithelial_Analysis
Rscript scripts/setup_environment.R

# Run full pipeline
source("scripts/00_run_full_pipeline.R")

# Or run individual steps
Rscript scripts/01_cell_type_annotation.R --dataset spatial
Rscript scripts/02_trajectory_analysis.R
```

## Project Structure

```
├── config/config.yaml       # All parameters and paths
├── scripts/
│   ├── 01_cell_type_annotation.R   # Epithelial classification
│   ├── 02_trajectory_analysis.R    # Trajectory inference
│   └── generate_figures.R          # Figure regeneration
├── R/                       # Core functions
├── results/                 # Output tables, figures, objects
└── docs/                    # Detailed documentation
```

## Pipeline Overview

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_cell_type_annotation.R` | Filter epithelial cells, classify subtypes (Basal-like, Intermediate, Specialized, Transit-Amplifying) |
| 2 | `02_trajectory_analysis.R` | Multi-method trajectory inference (Monocle3, Slingshot), consensus pseudotime |
| 3 | `generate_figures.R` | Publication-ready figures |

## Supported Datasets

```bash
Rscript scripts/01_cell_type_annotation.R --dataset spatial   # 10x Visium
Rscript scripts/01_cell_type_annotation.R --dataset snrnaseq  # GSE215932
Rscript scripts/01_cell_type_annotation.R --dataset acp_scn   # ACP scRNA-seq
```

## Documentation

- [01_cell_type_annotation.R](docs/01_cell_type_annotation.md) - Epithelial filtering and subtype classification
- [02_trajectory_analysis.R](docs/02_trajectory_analysis.md) - Trajectory inference and validation
- [Configuration Guide](docs/configuration.md) - All parameters explained

## Requirements

- R ≥ 4.2.0
- Seurat ≥ 5.0.0
- Monocle3, slingshot

Run `renv::restore()` for exact package versions.

## Data

- **Spatial**: 10x Visium from 4 ACP samples
- **Validation**: [GSE215932](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215932)

## Citation

```
[Author et al.], "Epithelial Differentiation Trajectories in Adamantinomatous 
Craniopharyngioma Revealed by Spatial Transcriptomics", [Journal], [Year].
```

## License

MIT
