# renv activation script
# This will be replaced when renv::init() is run

local({
  # Check if renv is installed
  if (!requireNamespace("renv", quietly = TRUE)) {
    message("renv not found. Please run: install.packages('renv')")
    message("Then run: renv::restore()")
    return(invisible())
  }
  
  # Activate renv
  renv::activate()
})
