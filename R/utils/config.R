# ==============================================================================
# R/utils/config.R - Configuration Loading and Management
# ==============================================================================

#' Load Project Configuration
#'
#' @param config_path Path to config.yaml (default: "config/config.yaml")
#' @return List containing all configuration parameters
#' @export
load_config <- function(config_path = NULL) {

  if (is.null(config_path)) {
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

  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install with: install.packages('yaml')")
  }

  config <- yaml::read_yaml(config_path)

  if (!is.null(config$reproducibility$seed)) {
    set.seed(config$reproducibility$seed)
    message(paste("Random seed set to:", config$reproducibility$seed))
  }

  config$.config_path <- normalizePath(config_path)
  config$.project_root <- dirname(dirname(normalizePath(config_path)))

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
#' @param config Configuration list
#' @param relative_path Relative path string
#' @return Absolute path
#' @export
get_path <- function(config, relative_path) {
  file.path(config$.project_root, relative_path)
}

#' Ensure Directory Exists
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
#' @param config Configuration list
#' @param signature_name Name of the signature
#' @param available_genes Character vector of available genes
#' @param warn_missing Warn about missing genes (default TRUE)
#' @return Character vector of genes
#' @export
get_signature <- function(config, signature_name, available_genes = NULL, warn_missing = TRUE) {

  sig <- NULL

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
#' @param config Configuration list
#' @param available_genes Character vector of available genes
#' @param min_genes Minimum genes required to include signature
#' @return Named list of gene signatures
#' @export
get_all_signatures <- function(config, available_genes = NULL, min_genes = 3) {

  all_sigs <- list()

  for (name in names(config$signatures$epithelial_subtypes)) {
    all_sigs[[name]] <- config$signatures$epithelial_subtypes[[name]]
  }

  other_sigs <- setdiff(names(config$signatures), "epithelial_subtypes")
  for (name in other_sigs) {
    all_sigs[[name]] <- config$signatures[[name]]
  }

  if (!is.null(available_genes)) {
    all_sigs <- lapply(all_sigs, function(sig) {
      intersect(sig, available_genes)
    })

    all_sigs <- all_sigs[sapply(all_sigs, length) >= min_genes]
  }

  return(all_sigs)
}

#' Get Color Palette
#'
#' @param config Configuration list
#' @param palette_name Name of palette
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

#' Print Configuration Summary
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
