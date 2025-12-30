# Troubleshooting

Common issues and solutions for the ACP Epithelial Differentiation Analysis pipeline.

## Table of Contents

- [Classification Issues](#classification-issues)
- [Gene Coverage Problems](#gene-coverage-problems)
- [Trajectory Analysis Issues](#trajectory-analysis-issues)
- [Memory and Performance](#memory-and-performance)
- [Data Format Errors](#data-format-errors)
- [Visualization Problems](#visualization-problems)

---

## Classification Issues

### High Unassigned Rate (>30%)

**Symptoms**: Large proportion of cells classified as "Unassigned" despite being epithelial.

**Causes & Solutions**:

1. **Threshold too stringent**
   ```yaml
   # config.yaml - lower the threshold
   classification:
     min_score_threshold: 0.05  # default is 0.10
   ```

2. **Score distributions differ by dataset**
   
   Check diagnostic plots in `results/figures/supplementary/01_cell_type_annotation_diagnostics_{dataset}.pdf`. If scores are clustered near zero:
   - snRNA-seq typically has lower module scores than spatial data
   - Try threshold of 0.02-0.05 for snRNA-seq

3. **Normalization differences**
   
   Module scores depend on normalization. Verify your object uses log-normalized counts:
   ```r
   # Check
   seurat_obj@assays$RNA@data[1:5, 1:5]
   
   # Re-normalize if needed
   seurat_obj <- NormalizeData(seurat_obj)
   ```

### Too Many Transit-Amplifying Cells

**Symptoms**: >20% of cells classified as TA when biological expectation is <10%.

**Solutions**:

1. **Increase proliferation threshold**
   ```yaml
   classification:
     transit_amplifying:
       proliferation_threshold: 0.20  # default is 0.15
       override_threshold: 0.40       # default is 0.30
   ```

2. **Require cell cycle phase**
   ```yaml
   classification:
     transit_amplifying:
       require_cycling: true  # only S/G2M cells can be TA
   ```

3. **Check for batch effects** - proliferation signatures can be inflated in certain batches

### No Hybrid States Detected

**Symptoms**: No Basal-Intermediate or Intermediate-Specialized cells.

**Causes**:
- Hybrid detection requires cells to score above threshold on both adjacent subtypes
- Check if `hybrid_detection.enabled: true` in config
- May indicate clean separation between subtypes (not necessarily a problem)

---

## Gene Coverage Problems

### "Skipping signature X: insufficient genes"

**Symptoms**: Warning that a signature has <3 genes present.

**Diagnosis**:
```r
# Check how many signature genes are in your data
signature_genes <- config$signatures$epithelial_subtypes$Basal_like
present <- signature_genes %in% rownames(seurat_obj)
cat(sum(present), "of", length(signature_genes), "genes present\n")
print(signature_genes[!present])  # Missing genes
```

**Solutions**:

1. **Use less-filtered data**
   - Integrated objects often have reduced gene counts
   - Try using the object before integration or heavy filtering

2. **Check gene naming**
   - Ensure gene names match (e.g., HGNC symbols vs Ensembl IDs)
   ```r
   # Convert if needed
   library(biomaRt)
   # ... conversion code
   ```

3. **Update signatures for available genes**
   - Add alternative markers that are present
   - Minimum 5-8 genes recommended per signature

### Low Gene Count Warning

**Symptoms**: "Object has only X genes (expected >15,000)"

**Common with**: Integrated/batch-corrected objects, heavily filtered data

**Solutions**:
- Use raw or less-filtered object for classification
- Transfer labels back to integrated object after classification:
  ```r
  # Classify on full object
  full_obj <- readRDS("full_object.rds")
  # ... run classification ...
  
  # Transfer to integrated
  integrated_obj$epithelial_subtype <- full_obj$epithelial_subtype[Cells(integrated_obj)]
  ```

---

## Trajectory Analysis Issues

### "No root cells found in Basal-like"

**Symptoms**: Monocle3 fails because no cells match root subtype.

**Solutions**:

1. **Check subtype names match exactly**
   ```r
   table(seurat_obj$epithelial_subtype)
   # Ensure "Basal-like" exists (check hyphen vs underscore)
   ```

2. **Change root subtype**
   ```yaml
   trajectory:
     monocle3:
       root_subtype: "Basal_like"  # match your naming
   ```

3. **Lower classification threshold** to get more Basal-like cells

### Disconnected Trajectory Graph

**Symptoms**: Multiple disconnected components, trajectory doesn't span all subtypes.

**Solutions**:

1. **Adjust clustering resolution**
   ```yaml
   trajectory:
     monocle3:
       cluster_resolution: 0.001  # lower = fewer clusters
   ```

2. **Disable partitioning**
   ```yaml
   trajectory:
     monocle3:
       use_partition: false
   ```

3. **Increase k for nearest neighbors**
   ```r
   cds <- learn_graph(cds, learn_graph_control = list(nn.k = 25))
   ```

### Pseudotime Values All NA

**Symptoms**: `pseudotime` column is all NA after trajectory inference.

**Causes**:
- Root not set properly
- Cells not connected to root in graph

**Solutions**:
```r
# Check root selection
root_cells <- WhichCells(seurat_obj, expression = epithelial_subtype == "Basal-like")
length(root_cells)  # Should be >0

# Manually set root if needed
cds <- order_cells(cds, root_cells = root_cells[1:10])
```

---

## Memory and Performance

### R Session Crashes / Out of Memory

**Symptoms**: R crashes when loading data or during analysis.

**Solutions**:

1. **Increase memory allocation**
   ```bash
   # SLURM
   #SBATCH --mem=64G
   
   # Local R
   R --max-mem-size=32G
   ```

2. **Process in chunks**
   ```r
   # Split by sample
   samples <- unique(seurat_obj$sample_id)
   for (s in samples) {
     subset_obj <- subset(seurat_obj, sample_id == s)
     # ... process ...
   }
   ```

3. **Reduce object size before analysis**
   ```r
   # Remove unnecessary assays
   seurat_obj[["SCT"]] <- NULL
   seurat_obj[["integrated"]] <- NULL
   
   # Slim down
   seurat_obj <- DietSeurat(seurat_obj, assays = "RNA")
   ```

### Slow Module Score Calculation

**Symptoms**: `AddModuleScore()` takes >10 minutes.

**Solutions**:

1. **Reduce control genes**
   ```r
   seurat_obj <- AddModuleScore(seurat_obj, features = gene_list, ctrl = 50)  # default is 100
   ```

2. **Parallelize** (if running multiple signatures)
   ```r
   library(future)
   plan(multisession, workers = 4)
   ```

---

## Data Format Errors

### "Column X not found in metadata"

**Symptoms**: Script fails looking for epithelial identification column.

**Solutions**:

1. **Check available columns**
   ```r
   colnames(seurat_obj@meta.data)
   ```

2. **Update config for your column name**
   - Modify the dataset configuration block in `01_cell_type_annotation.R` (~line 100)

3. **Create the expected column**
   ```r
   # Example: create is_epithelial from cell_type
   seurat_obj$is_epithelial <- seurat_obj$cell_type == "Epithelial"
   ```

### "Assay RNA not found"

**Symptoms**: Default assay doesn't exist.

**Solutions**:
```r
# Check available assays
Assays(seurat_obj)

# Set correct default
DefaultAssay(seurat_obj) <- "RNA"  # or "Spatial", "SCT", etc.
```

### Seurat v4 vs v5 Compatibility

**Symptoms**: Errors about layers, or `@data` slot access failing.

**Solutions**:

```r
# Check version
packageVersion("Seurat")

# Seurat v5: access normalized data
LayerData(seurat_obj, assay = "RNA", layer = "data")

# Seurat v4: access normalized data
seurat_obj@assays$RNA@data

# Convert v5 to v4-style if needed
seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")
```

---

## Visualization Problems

### Plots Not Rendering / Blank PDFs

**Solutions**:

1. **Check output directory exists**
   ```r
   dir.create("results/figures/main", recursive = TRUE, showWarnings = FALSE)
   ```

2. **Close open devices**
   ```r
   dev.off()  # Run multiple times if needed
   graphics.off()
   ```

3. **Specify device explicitly**
   ```r
   ggsave("plot.pdf", plot, device = cairo_pdf)
   ```

### Colors Don't Match Between Plots

**Solutions**:

1. **Use consistent color mapping**
   ```r
   colors <- get_colors(config, "epithelial_subtypes")
   
   # Force factor levels
   seurat_obj$epithelial_subtype <- factor(
     seurat_obj$epithelial_subtype,
     levels = names(colors)
   )
   ```

2. **Check for NA values creating extra categories**
   ```r
   table(seurat_obj$epithelial_subtype, useNA = "ifany")
   ```

### UMAP/Spatial Coordinates Missing

**Symptoms**: "DimReduc object umap not found"

**Solutions**:
```r
# Check available reductions
Reductions(seurat_obj)

# Run UMAP if missing
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# For spatial data, check image
Images(seurat_obj)
SpatialDimPlot(seurat_obj)  # Test
```

---

## Getting Help

### Generating a Debug Report

Run this to create a report for troubleshooting:

```r
# Save session info
sink("debug_report.txt")
cat("=== Session Info ===\n")
sessionInfo()
cat("\n=== Object Summary ===\n")
print(seurat_obj)
cat("\n=== Metadata Columns ===\n")
print(colnames(seurat_obj@meta.data))
cat("\n=== Assays ===\n")
print(Assays(seurat_obj))
cat("\n=== Gene Count ===\n")
print(nrow(seurat_obj))
sink()
```

### Common Log Messages Explained

| Message | Meaning | Action |
|---------|---------|--------|
| "Using column X for epithelial identification" | Info | None needed |
| "X cells identified as epithelial (Y%)" | Info | Check if % is reasonable |
| "Signature X: using N of M genes" | Some genes missing | OK if N ≥ 5 |
| "Skipping signature X" | <3 genes present | Check gene coverage |
| "High unassigned rate" | >30% unassigned | Lower threshold |
| "No cells in cycling phase" | Cell cycle scoring issue | Check phase assignment |

---

## Still Stuck?

1. Check the diagnostic PDF outputs for visual clues
2. Verify input data format matches expectations
3. Run scripts interactively to identify exact failure point
4. Open an issue on GitHub with debug report
