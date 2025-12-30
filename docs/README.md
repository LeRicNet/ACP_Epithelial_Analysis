# Documentation

Detailed documentation for the ACP Epithelial Differentiation Analysis pipeline.

## Scripts

| Script | Description |
|--------|-------------|
| [01_cell_type_annotation.md](01_cell_type_annotation.md) | Epithelial cell identification, module scoring, subtype classification |
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
# Full pipeline
source("scripts/00_run_full_pipeline.R")

# Individual steps
Rscript scripts/01_cell_type_annotation.R --dataset spatial
Rscript scripts/02_trajectory_analysis.R --cores 16
Rscript scripts/generate_figures.R 1 2 3
```

### Key Outputs

| File | Location |
|------|----------|
| Classified cells | `results/objects/01_seurat_annotated_{dataset}.rds` |
| Pseudotime | `results/objects/02_seurat_with_pseudotime.rds` |
| Summary tables | `results/tables/` |
| Figures | `results/figures/main/` and `results/figures/supplementary/` |

### Configuration Quick Reference

```yaml
# Key classification parameters
classification:
  min_score_threshold: 0.10      # Minimum score for subtype assignment
  transit_amplifying:
    proliferation_threshold: 0.15  # TA identification threshold
    require_cycling: true          # Require S/G2M phase

# Key trajectory parameters  
trajectory:
  methods: [monocle3, slingshot, diffusion]
  monocle3:
    root_subtype: "Basal-like"
```
