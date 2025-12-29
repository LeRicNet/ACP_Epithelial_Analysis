# ==============================================================================
# R/utils/config.R - Configuration Loading and Management
# ==============================================================================
# Functions for loading and accessing the YAML configuration file
# ==============================================================================

#' Load Project Configuration
#'
#' Loads the YAML configuration file and returns a list object.
#' Also sets up the project paths and random seed.
#'
#' @param config_path Path to config.yaml (default: "config/config.yaml")
#' @return List containing all configuration parameters
#' @export
load_config <- function(config_path = NULL) {
  
  # Find config file
 if (is.null(config_path)) {
    # Try multiple locations
    possible_paths <- c(
      "config/config.yaml",
      "../config/config.yaml",
      file.path(here::here(), "config/config.yaml")
    )
    
    for (path in possible_paths) {
      if (file.exists(path)) {
        config_path <- path
        break
      }
    }
    
    if (is.null(config_path)) {
      stop("Could not find config.yaml. Please specify path explicitly.")
    }
  }
  
  if (!file.exists(config_path)) {
    stop(paste("Configuration file not found:", config_path))
  }
  
  # Load YAML
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install with: install.packages('yaml')")
  }
  
  config <- yaml::read_yaml(config_path)
  
  # Set random seed for reproducibility
  if (!is.null(config$reproducibility$seed)) {
    set.seed(config$reproducibility$seed)
    message(paste("Random seed set to:", config$reproducibility$seed))
  }
  
  # Store config path for reference
  config$.config_path <- normalizePath(config_path)
  config$.project_root <- dirname(dirname(normalizePath(config_path)))
  
  # Validate required sections exist
  required_sections <- c("paths", "signatures", "classification", "visualization")
  missing <- setdiff(required_sections, names(config))
  if (length(missing) > 0) {
    warning(paste("Missing config sections:", paste(missing, collapse = ", ")))
  }
  
  message(paste("Configuration loaded from:", config_path))
  
  return(config)
}

#' Get Full Path
#'
#' Converts a relative path from config to an absolute path
#'
#' @param config Configuration list
#' @param relative_path Relative path string
#' @return Absolute path
#' @export
get_path <- function(config, relative_path) {
  file.path(config$.project_root, relative_path)
}

#' Ensure Directory Exists
#'
#' Creates directory if it doesn't exist
#'
#' @param path Directory path
#' @return Path (invisibly)
#' @export
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    message(paste("Created directory:", path))
  }
  invisible(path)
}

#' Get Signature Genes
#'
#' Retrieves a gene signature from config, filtering for genes present in data
#'
#' @param config Configuration list
#' @param signature_name Name of the signature
#' @param available_genes Character vector of available genes (e.g., rownames(seurat))
#' @param warn_missing Warn about missing genes (default TRUE)
#' @return Character vector of genes
#' @export
get_signature <- function(config, signature_name, available_genes = NULL, warn_missing = TRUE) {
  
  # Try to find signature in different locations
  sig <- NULL
  
  # Check epithelial_subtypes first
  if (signature_name %in% names(config$signatures$epithelial_subtypes)) {
    sig <- config$signatures$epithelial_subtypes[[signature_name]]
  } else if (signature_name %in% names(config$signatures)) {
    sig <- config$signatures[[signature_name]]
  } else if (signature_name %in% names(config$markers)) {
    sig <- config$markers[[signature_name]]
  }
  
  if (is.null(sig)) {
    stop(paste("Signature not found:", signature_name))
  }
  
  # Filter for available genes if provided
  if (!is.null(available_genes)) {
    present <- intersect(sig, available_genes)
    missing <- setdiff(sig, available_genes)
    
    if (warn_missing && length(missing) > 0) {
      warning(paste0(
        signature_name, ": ", length(present), "/", length(sig), " genes present. ",
        "Missing: ", paste(missing, collapse = ", ")
      ))
    }
    
    return(present)
  }
  
  return(sig)
}

#' Get All Signatures
#'
#' Returns all signatures as a named list, optionally filtered for available genes
#'
#' @param config Configuration list
#' @param available_genes Character vector of available genes
#' @param min_genes Minimum genes required to include signature
#' @return Named list of gene signatures
#' @export
get_all_signatures <- function(config, available_genes = NULL, min_genes = 3) {
  
  all_sigs <- list()
  
  # Get epithelial subtype signatures
  for (name in names(config$signatures$epithelial_subtypes)) {
    all_sigs[[name]] <- config$signatures$epithelial_subtypes[[name]]
  }
  
  # Get other signatures
  other_sigs <- setdiff(names(config$signatures), "epithelial_subtypes")
  for (name in other_sigs) {
    all_sigs[[name]] <- config$signatures[[name]]
  }
  
  # Filter for available genes
  if (!is.null(available_genes)) {
    all_sigs <- lapply(all_sigs, function(sig) {
      intersect(sig, available_genes)
    })
    
    # Remove signatures with too few genes
    all_sigs <- all_sigs[sapply(all_sigs, length) >= min_genes]
  }
  
  return(all_sigs)
}

#' Get Color Palette
#'
#' Retrieves color palette from config
#'
#' @param config Configuration list
#' @param palette_name Name of palette (e.g., "epithelial_subtypes", "samples")
#' @return Named character vector of colors
#' @export
get_colors <- function(config, palette_name) {
  colors <- config$visualization$colors[[palette_name]]
  if (is.null(colors)) {
    stop(paste("Color palette not found:", palette_name))
  }
  return(unlist(colors))
}

#' Get Figure Size
#'
#' Retrieves figure dimensions from config
#'
#' @param config Configuration list
#' @param size_name Name of size preset
#' @return Named list with width and height
#' @export
get_figure_size <- function(config, size_name) {
  size <- config$visualization$figure_sizes[[size_name]]
  if (is.null(size)) {
    warning(paste("Figure size not found:", size_name, "- using defaults"))
    return(list(width = 8, height = 6))
  }
  return(size)
}

#' Save Figure
#'
#' Saves a ggplot figure in multiple formats based on config
#'
#' @param plot ggplot object
#' @param filename Base filename (without extension)
#' @param config Configuration list
#' @param subdir Subdirectory under figures (e.g., "main", "supplementary")
#' @param width Width in inches (or NULL to use config default)
#' @param height Height in inches (or NULL to use config default)
#' @param size_preset Size preset name from config (overrides width/height)
#' @return Invisibly returns the file paths
#' @export
save_figure <- function(plot, filename, config, subdir = "main", 
                        width = NULL, height = NULL, size_preset = NULL) {
  
  # Determine output directory
  if (subdir == "main") {
    out_dir <- get_path(config, config$paths$figures_main_dir)
  } else if (subdir == "supplementary") {
    out_dir <- get_path(config, config$paths$figures_supp_dir)
  } else {
    out_dir <- get_path(config, file.path("results/figures", subdir))
  }
  ensure_dir(out_dir)
  
  # Determine dimensions
  if (!is.null(size_preset)) {
    dims <- get_figure_size(config, size_preset)
    width <- dims$width
    height <- dims$height
  }
  
  if (is.null(width)) width <- 8
  if (is.null(height)) height <- 6
  
  # Get output formats
  formats <- config$output$figures$formats
  if (is.null(formats)) formats <- c("pdf", "png")
  
  dpi <- config$output$figures$png_dpi
  if (is.null(dpi)) dpi <- 300
  
  # Save in each format
  saved_paths <- c()
  for (fmt in formats) {
    filepath <- file.path(out_dir, paste0(filename, ".", fmt))
    
    if (fmt == "pdf") {
      ggplot2::ggsave(filepath, plot, width = width, height = height, device = "pdf")
    } else if (fmt == "png") {
      ggplot2::ggsave(filepath, plot, width = width, height = height, dpi = dpi, device = "png")
    }
    
    saved_paths <- c(saved_paths, filepath)
    message(paste("Saved:", filepath))
  }
  
  invisible(saved_paths)
}

#' Save Table
#'
#' Saves a data frame in the format specified by config
#'
#' @param df Data frame to save
#' @param filename Base filename (without extension)
#' @param config Configuration list
#' @return Invisibly returns the file path
#' @export
save_table <- function(df, filename, config) {
  
  out_dir <- get_path(config, config$paths$tables_dir)
  ensure_dir(out_dir)
  
  format <- config$output$tables$format
  if (is.null(format)) format <- "csv"
  
  filepath <- file.path(out_dir, paste0(filename, ".", format))
  
  if (format == "csv") {
    write.csv(df, filepath, row.names = FALSE)
  } else if (format == "tsv") {
    write.table(df, filepath, sep = "\t", row.names = FALSE, quote = FALSE)
  } else if (format == "xlsx") {
    if (requireNamespace("writexl", quietly = TRUE)) {
      writexl::write_xlsx(df, filepath)
    } else {
      warning("writexl not installed, saving as CSV instead")
      filepath <- file.path(out_dir, paste0(filename, ".csv"))
      write.csv(df, filepath, row.names = FALSE)
    }
  }
  
  message(paste("Saved:", filepath))
  invisible(filepath)
}

#' Print Configuration Summary
#'
#' Prints a summary of the loaded configuration
#'
#' @param config Configuration list
#' @export
print_config_summary <- function(config) {
  cat("\n")
  cat("=======================================================\n")
  cat("  ACP Epithelial Differentiation Analysis\n")
  cat("=======================================================\n")
  cat("Project:", config$project$name, "\n")
  cat("Version:", config$project$version, "\n")
  cat("Config:", config$.config_path, "\n")
  cat("Root:", config$.project_root, "\n")
  cat("\nSignatures loaded:", length(get_all_signatures(config)), "\n")
  cat("Classification method:", config$classification$method, "\n")
  cat("Root subtype for trajectory:", config$trajectory$monocle3$root_subtype, "\n")
  cat("Random seed:", config$reproducibility$seed, "\n")
  cat("=======================================================\n\n")
}
