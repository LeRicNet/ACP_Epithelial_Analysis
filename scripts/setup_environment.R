#!/usr/bin/env Rscript
# ==============================================================================
# scripts/setup_environment.R
# ==============================================================================
# Initialize renv and install all required packages
#
# Run this script ONCE after cloning the repository:
#   Rscript scripts/setup_environment.R
# ==============================================================================

message("================================================================")
message("  Setting up ACP Epithelial Analysis Environment")
message("================================================================\n")

# Check R version
r_version <- paste(R.version$major, R.version$minor, sep = ".")
message(paste("R version:", r_version))

if (as.numeric(R.version$major) < 4 || 
    (as.numeric(R.version$major) == 4 && as.numeric(R.version$minor) < 2)) {
  warning("R version 4.2.0 or higher is recommended")
}

# Initialize renv if not already done
if (!file.exists("renv.lock")) {
  message("\nInitializing renv...")
  
  # Install renv if needed
  if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
  }
  
  # Initialize renv
  renv::init(bare = TRUE)
  
  # Install core packages
  message("\nInstalling core packages...")
  
  packages <- c(
    # CRAN packages
    "Seurat",
    "ggplot2",
    "patchwork",
    "dplyr",
    "tidyr",
    "yaml",
    "here",
    "pheatmap",
    "viridis",
    "RColorBrewer",
    "writexl",
    "reshape2",
    
    # Bioconductor packages
    "BiocManager"
  )
  
  # Install CRAN packages
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing:", pkg))
      tryCatch({
        install.packages(pkg)
      }, error = function(e) {
        message(paste("  Failed to install", pkg, ":", e$message))
      })
    }
  }
  
  # Install Bioconductor packages
  message("\nInstalling Bioconductor packages...")
  
  bioc_packages <- c(
    "SingleCellExperiment",
    "SummarizedExperiment",
    "UCell"
  )
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing:", pkg))
      tryCatch({
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      }, error = function(e) {
        message(paste("  Failed to install", pkg, ":", e$message))
      })
    }
  }
  
  # Install Monocle3 (requires special installation)
  message("\nInstalling Monocle3...")
  if (!requireNamespace("monocle3", quietly = TRUE)) {
    tryCatch({
      BiocManager::install(c("BiocGenerics", "DelayedArray", "DelayedMatrixStats",
                             "limma", "lme4", "S4Vectors", "SingleCellExperiment",
                             "SummarizedExperiment", "batchelor", "HDF5Array",
                             "terra", "ggrastr"),
                           ask = FALSE, update = FALSE)
      
      if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools")
      }
      devtools::install_github("cole-trapnell-lab/monocle3")
    }, error = function(e) {
      message(paste("Failed to install Monocle3:", e$message))
      message("Please install manually from: https://cole-trapnell-lab.github.io/monocle3/")
    })
  }
  
  # Install SeuratWrappers
  message("\nInstalling SeuratWrappers...")
  if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
    tryCatch({
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes")
      }
      remotes::install_github("satijalab/seurat-wrappers")
    }, error = function(e) {
      message(paste("Failed to install SeuratWrappers:", e$message))
    })
  }
  
  # Snapshot the environment
  message("\nCreating renv snapshot...")
  renv::snapshot()
  
} else {
  message("renv.lock found - restoring environment...")
  renv::restore()
}

# Verify installation
message("\n================================================================")
message("  Verifying Installation")
message("================================================================\n")

required_packages <- c(
  "Seurat", "ggplot2", "patchwork", "dplyr", "tidyr", "yaml",
  "monocle3", "SeuratWrappers"
)

all_installed <- TRUE
for (pkg in required_packages) {
  installed <- requireNamespace(pkg, quietly = TRUE)
  status <- if (installed) "OK" else "MISSING"
  message(sprintf("  %-20s %s", pkg, status))
  if (!installed) all_installed <- FALSE
}

if (all_installed) {
  message("\n✓ All required packages installed successfully!")
} else {
  message("\n✗ Some packages are missing. Please install manually.")
}

message("\n================================================================")
message("  Setup Complete")
message("================================================================")
message("\nTo run the analysis:")
message("  source('scripts/00_run_full_pipeline.R')")
message("\nTo regenerate figures:")
message("  source('scripts/generate_figures.R')")
