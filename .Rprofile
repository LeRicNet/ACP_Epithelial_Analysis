# ==============================================================================
# .Rprofile - Project-level R startup configuration
# ==============================================================================
# This file is automatically sourced when R starts in this project directory.
# It activates renv to ensure reproducible package management.
# ==============================================================================

# Activate renv
source("renv/activate.R")

# Set options
options(
  # Don't ask about saving workspace
  save.defaults = list(ask = FALSE),
  
  # CRAN mirror
  repos = c(CRAN = "https://cloud.r-project.org"),
  
  # Warn on partial matching
  warnPartialMatchArgs = TRUE,
  warnPartialMatchAttr = TRUE,
  warnPartialMatchDollar = TRUE,
  
  # Show more warnings
  warn = 1
)

# Welcome message
if (interactive()) {
  message("\n")
  message("==============================================")
  message("  ACP Epithelial Differentiation Analysis")
  message("==============================================")
  message("  Type: source('scripts/00_run_full_pipeline.R') to run analysis")
  message("  Type: source('scripts/generate_figures.R') to regenerate figures")
  message("\n")
}

# Python environment for reticulate
Sys.setenv(RETICULATE_PYTHON = "/home/rstudio/interfaces/ACP_Epithelial_Analysis/.venv/bin/python")
