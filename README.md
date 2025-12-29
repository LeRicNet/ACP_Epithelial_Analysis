# !! UNDER CONSTRUCTION !!

# ACP Epithelial Differentiation Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Analysis code for the manuscript: **"Epithelial Differentiation Trajectories in Adamantinomatous Craniopharyngioma Revealed by Spatial Transcriptomics"**

## Overview

This repository contains all code required to reproduce the epithelial subtype classification and trajectory analysis presented in our manuscript. The analysis pipeline:

1. **Classifies epithelial cells** using module score-based approach (replacing binary thresholds)
2. **Identifies Transit-Amplifying populations** via proliferation markers
3. **Reconstructs differentiation trajectories** using Monocle3
4. **Validates findings** across individual samples and independent datasets

## Quick Start

```r
# Clone the repository
git clone https://github.com/your-username/ACP_Epithelial_Analysis.git
cd ACP_Epithelial_Analysis

# Install dependencies using renv
install.packages("renv")
renv::restore()

# Run the full analysis pipeline
source("scripts/00_run_full_pipeline.R")

# Or regenerate specific figures
source("scripts/generate_figures.R")  # All figures
Rscript scripts/generate_figures.R 1  # Figure 1 only
```

## Project Structure

```
ACP_Epithelial_Analysis/
├── config/
│   └── config.yaml          # Central configuration (paths, signatures, parameters)
├── data/
│   ├── raw/                  # Original data files (not tracked in git)
│   ├── processed/            # Processed intermediates
│   └── external/             # External validation datasets
├── R/
│   ├── utils/
│   │   └── config.R          # Configuration loading utilities
│   ├── classification/
│   │   └── module_score_classification.R
│   ├── trajectory/
│   │   └── trajectory_analysis.R
│   └── figures/
│       └── figure_generators.R
├── scripts/
│   ├── 00_run_full_pipeline.R    # Master analysis script
│   └── generate_figures.R         # Figure regeneration script
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

## Requirements

### R Version
- R >= 4.2.0

### Key Packages
- Seurat (>= 5.0.0)
- Monocle3
- SeuratWrappers
- UCell (optional, for module scoring)
- ggplot2
- patchwork
- dplyr/tidyr
- yaml

All package versions are locked via `renv.lock` for exact reproducibility.

## Data Availability

### Input Data
- **Spatial transcriptomics**: 10x Visium data from 4 ACP tumor samples
- **snRNA-seq validation**: GSE215932 (Xu et al., 2024)

Due to size constraints, raw data is not included in this repository. Data can be obtained from:
- GEO: [GSE215932](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215932)
- [Link to spatial data repository]

After downloading, place files according to paths in `config/config.yaml`.

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
Monocle3-based pseudotemporal ordering:
```r
# Run trajectory analysis
cds <- run_trajectory_analysis(seurat_obj, config)

# Per-sample validation
per_sample_results <- run_per_sample_trajectory(seurat_obj, config)
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
