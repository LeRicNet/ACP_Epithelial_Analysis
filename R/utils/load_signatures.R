#!/usr/bin/env Rscript
# ==============================================================================
# R/utils/load_signatures.R
# ==============================================================================
# Utility functions for loading gene signatures from external files
# Supports: GMT, GRP, TXT, CSV formats and MSigDB download
#
# Usage:
#   source("R/utils/load_signatures.R")
#   sigs <- load_signatures_from_config(config)
#   senescence_genes <- load_gmt("config/signatures/SAUL_SEN_MAYO.gmt")[[1]]
#
# Author: [Your Name]
# Date: 2025-12-30
# ==============================================================================

# --- Load GMT File ------------------------------------------------------------
#' Load gene sets from GMT format (MSigDB standard)
#'
#' @param gmt_file Path to .gmt file
#' @return Named list of gene vectors
load_gmt <- function(gmt_file) {
  if (!file.exists(gmt_file)) {
    stop("GMT file not found: ", gmt_file)
  }

  lines <- readLines(gmt_file)
  gene_sets <- list()

  for (line in lines) {
    if (nchar(trimws(line)) == 0) next

    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 3) next

    set_name <- fields[1]
    # fields[2] is description, skip it
    genes <- fields[3:length(fields)]
    genes <- genes[genes != "" & !is.na(genes)]

    gene_sets[[set_name]] <- genes
  }

  message(sprintf("Loaded %d gene sets from %s", length(gene_sets), basename(gmt_file)))
  return(gene_sets)
}

# --- Load GRP File ------------------------------------------------------------
#' Load gene set from GRP format (one gene per line, first line is set name)
#'
#' @param grp_file Path to .grp file
#' @return Character vector of genes
load_grp <- function(grp_file) {
  if (!file.exists(grp_file)) {
    stop("GRP file not found: ", grp_file)
  }

  lines <- readLines(grp_file)
  lines <- trimws(lines)
  lines <- lines[lines != "" & !grepl("^#", lines)]

  message(sprintf("Loaded %d genes from %s", length(lines), basename(grp_file)))
  return(lines)
}

# --- Load TXT File ------------------------------------------------------------
#' Load gene set from simple text file (one gene per line)
#'
#' @param txt_file Path to .txt file
#' @param skip_header Logical, skip first line if TRUE
#' @return Character vector of genes
load_txt <- function(txt_file, skip_header = FALSE) {
  if (!file.exists(txt_file)) {
    stop("TXT file not found: ", txt_file)
  }

  lines <- readLines(txt_file)
  if (skip_header && length(lines) > 1) {
    lines <- lines[-1]
  }

  lines <- trimws(lines)
  lines <- lines[lines != "" & !grepl("^#", lines)]

  message(sprintf("Loaded %d genes from %s", length(lines), basename(txt_file)))
  return(lines)
}

# --- Load CSV File ------------------------------------------------------------
#' Load gene set from CSV file
#'
#' @param csv_file Path to .csv file
#' @param gene_column Column name or index containing gene symbols
#' @return Character vector of genes
load_csv <- function(csv_file, gene_column = 1) {
  if (!file.exists(csv_file)) {
    stop("CSV file not found: ", csv_file)
  }

  df <- read.csv(csv_file, stringsAsFactors = FALSE)

  if (is.character(gene_column)) {
    if (!gene_column %in% colnames(df)) {
      stop("Column '", gene_column, "' not found in CSV")
    }
    genes <- df[[gene_column]]
  } else {
    genes <- df[[gene_column]]
  }

  genes <- trimws(genes)
  genes <- genes[genes != "" & !is.na(genes)]
  genes <- unique(genes)

  message(sprintf("Loaded %d genes from %s", length(genes), basename(csv_file)))
  return(genes)
}

# --- Auto-detect Format -------------------------------------------------------
#' Load gene set from file, auto-detecting format
#'
#' @param file_path Path to gene set file
#' @return Gene set(s) - list for GMT, vector for others
load_signature_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Signature file not found: ", file_path)
  }

  ext <- tolower(tools::file_ext(file_path))

  result <- switch(ext,
                   "gmt" = load_gmt(file_path),
                   "grp" = load_grp(file_path),
                   "txt" = load_txt(file_path),
                   "csv" = load_csv(file_path),
                   "tsv" = {
                     df <- read.delim(file_path, stringsAsFactors = FALSE)
                     unique(trimws(df[[1]]))
                   },
                   stop("Unsupported file format: ", ext)
  )

  return(result)
}

# --- Download from MSigDB -----------------------------------------------------
#' Download gene set from MSigDB
#'
#' @param set_name Gene set name (e.g., "SAUL_SEN_MAYO")
#' @param collection Collection code (e.g., "C2", "H")
#' @param output_dir Directory to save downloaded file
#' @param format Download format: "gmt", "grp", "json"
#' @return Path to downloaded file
download_msigdb_geneset <- function(set_name,
                                    collection = NULL,
                                    output_dir = "config/signatures",
                                    format = "grp") {

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # MSigDB download URL pattern
  # Example: https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=SAUL_SEN_MAYO&fileType=grp
  base_url <- "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp"
  url <- sprintf("%s?geneSetName=%s&fileType=%s", base_url, set_name, format)

  output_file <- file.path(output_dir, paste0(set_name, ".", format))

  message("Downloading ", set_name, " from MSigDB...")

  tryCatch({
    download.file(url, output_file, mode = "w", quiet = TRUE)
    message("Saved to: ", output_file)
    return(output_file)
  }, error = function(e) {
    warning("Download failed: ", e$message)
    warning("Please download manually from: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/", set_name, ".html")
    return(NULL)
  })
}

# --- Load Signatures from Config ----------------------------------------------
#' Load all signatures from config, including external files
#'
#' @param config Config list loaded from config.yaml
#' @param project_root Project root directory for resolving relative paths
#' @return Named list of gene signatures
load_signatures_from_config <- function(config, project_root = ".") {
  signatures <- list()

  # Load inline signatures from config
  if (!is.null(config$signatures$epithelial_subtypes)) {
    for (name in names(config$signatures$epithelial_subtypes)) {
      signatures[[name]] <- config$signatures$epithelial_subtypes[[name]]
    }
  }

  # Load additional inline signatures
  inline_sigs <- c("Whorl_Wnt", "Senescence", "SASP", "CD44_SPP1",
                   "Collagen", "Matricellular", "Proliferation")
  for (sig_name in inline_sigs) {
    if (!is.null(config$signatures[[sig_name]])) {
      signatures[[sig_name]] <- config$signatures[[sig_name]]
    }
  }

  # Load external signature files
  if (!is.null(config$signatures$external_files)) {
    for (name in names(config$signatures$external_files)) {
      file_path <- config$signatures$external_files[[name]]

      # Resolve relative paths
      if (!startsWith(file_path, "/")) {
        file_path <- file.path(project_root, file_path)
      }

      if (file.exists(file_path)) {
        loaded <- load_signature_file(file_path)

        # If GMT returned multiple sets, prefix with file name
        if (is.list(loaded) && !is.null(names(loaded))) {
          if (length(loaded) == 1) {
            signatures[[name]] <- loaded[[1]]
          } else {
            for (set_name in names(loaded)) {
              signatures[[paste0(name, "_", set_name)]] <- loaded[[set_name]]
            }
          }
        } else {
          signatures[[name]] <- loaded
        }
      } else {
        warning("External signature file not found: ", file_path)
      }
    }
  }

  message(sprintf("Loaded %d total signatures", length(signatures)))
  return(signatures)
}

# --- Convenience: SenMayo Signature -------------------------------------------
#' Get SenMayo senescence signature
#'
#' Returns the SAUL_SEN_MAYO signature (124 genes)
#' Either loads from file or returns built-in version
#'
#' @param file_path Optional path to downloaded gene set file
#' @return Character vector of senescence genes
get_senmayo_signature <- function(file_path = NULL) {

  if (!is.null(file_path) && file.exists(file_path)) {
    return(load_signature_file(file_path))
  }

  # Built-in SenMayo signature (SAUL_SEN_MAYO from MSigDB)
  # Source: Saul et al. 2022, PMID 35974106
  senmayo <- c(
    "ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", "BMP2",
    "BMP6", "C3", "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", "CCL24",
    "CCL26", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", "CD55",
    "CD9", "CSF1", "CSF2", "CSF2RB", "CST4", "CTNNB1", "CTSB", "CXCL1",
    "CXCL10", "CXCL12", "CXCL16", "CXCL2", "CXCL3", "CXCL8", "CXCR2",
    "DKK1", "EDN1", "EGFR", "EGR1", "EGR3", "EREG", "ESM1", "ETS2",
    "FAS", "FGF1", "FGF2", "FGF7", "GDF15", "GEM", "GMFG", "HGF",
    "HMGB1", "ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2", "IGFBP3",
    "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15",
    "IL18", "IL1A", "IL1B", "IL2", "IL32", "IL6", "IL6ST", "IL7",
    "INHA", "INHBA", "IQGAP2", "ITGA2", "ITPKA", "JUN", "KITLG",
    "LCP1", "MIF", "MMP1", "MMP10", "MMP12", "MMP13", "MMP14", "MMP2",
    "MMP3", "MMP9", "NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF",
    "PIGF", "PLAT", "PLAU", "PLAUR", "PTBP1", "PTGER2", "PTGES",
    "RPS6KA5", "SCAMP4", "SELPLG", "SEMA3F", "SERPINB4", "SERPINE1",
    "SERPINE2", "SPP1", "SPX", "TIMP2", "TNFRSF10C", "TNFRSF11B",
    "TNFRSF1A", "TNFRSF1B", "TUBGCP2", "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2"
  )

  message("Using built-in SenMayo signature (124 genes)")
  return(senmayo)
}

# --- Example Usage ------------------------------------------------------------
if (FALSE) {
  # Load from downloaded file
  senescence <- load_grp("config/signatures/SAUL_SEN_MAYO.grp")

  # Or use built-in
  senescence <- get_senmayo_signature()

  # Download from MSigDB
  download_msigdb_geneset("SAUL_SEN_MAYO", output_dir = "config/signatures")

  # Load all signatures from config
  config <- yaml::read_yaml("config/config.yaml")
  all_sigs <- load_signatures_from_config(config)
}
