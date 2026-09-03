# Shared path configuration for the Cell Metabolism transfer workspace.

CONFIG_FILE <- if (exists(".local_config", inherits = TRUE)) {
  normalizePath(get(".local_config", inherits = TRUE), winslash = "/", mustWork = TRUE)
} else {
  config_candidates <- c(
    file.path(getwd(), "code", "Fig00_Config.R"),
    file.path(getwd(), "Fig00_Config.R")
  )
  config_candidate <- config_candidates[file.exists(config_candidates)][1]
  if (is.na(config_candidate)) {
    stop("Cannot locate code/Fig00_Config.R from the current working directory.", call. = FALSE)
  }
  normalizePath(config_candidate, winslash = "/", mustWork = TRUE)
}

PROJECT_ROOT_ENV <- Sys.getenv("ATM_PROJECT_ROOT", unset = NA_character_)
if (is.na(PROJECT_ROOT_ENV) || !nzchar(PROJECT_ROOT_ENV)) {
  stop("Set ATM_PROJECT_ROOT to the analysis project root before running release scripts.", call. = FALSE)
}
PROJECT_ROOT <- normalizePath(PROJECT_ROOT_ENV, winslash = "/", mustWork = TRUE)
WORKSPACE_ROOT <- normalizePath(file.path(PROJECT_ROOT, ".."), winslash = "/", mustWork = TRUE)

path_project <- function(...) file.path(PROJECT_ROOT, ...)
path_workspace <- function(...) file.path(WORKSPACE_ROOT, ...)
path_raw <- function(...) file.path(DIR_RAW, ...)
path_metadata <- function(...) file.path(DIR_METADATA, ...)
path_processed <- function(...) file.path(DIR_PROCESSED, ...)
path_intermediate <- function(...) file.path(DIR_INTERMEDIATE, ...)
path_analysis_ready <- function(...) file.path(DIR_ANALYSIS_READY, ...)
path_analysis_ready_source <- function(...) file.path(DIR_ANALYSIS_READY_SOURCE, ...)
path_ancillary <- function(...) file.path(DIR_ANCILLARY, ...)
path_qc <- function(...) file.path(DIR_QC, ...)
path_qc_source <- function(...) file.path(DIR_QC_SOURCE, ...)
path_results <- function(...) file.path(DIR_RESULTS, ...)

DIR_DATA <- path_workspace("Data")
DIR_RAW <- file.path(DIR_DATA, "Raw")
DIR_METADATA <- file.path(DIR_DATA, "Metadata")
DIR_PROCESSED <- file.path(DIR_DATA, "Processed")
DIR_INTERMEDIATE <- file.path(DIR_PROCESSED, "Intermediate")
DIR_ANALYSIS_READY <- file.path(DIR_PROCESSED, "AnalysisReady")
DIR_ANALYSIS_READY_SOURCE <- file.path(DIR_INTERMEDIATE, "AnalysisReady_Source")
DIR_ANCILLARY <- file.path(DIR_PROCESSED, "Ancillary")
DIR_QC <- file.path(DIR_PROCESSED, "QC")
DIR_QC_SOURCE <- file.path(DIR_INTERMEDIATE, "QC_Source")
DIR_RESULTS <- path_project("results")
DIR_R_LIBRARY <- path_project("R_Library")
DIR_R_CACHE <- path_project("R_Cache")
HISTORICAL_R_LIBS <- c(
  Sys.getenv("ATM_HISTORICAL_R_LIBRARY", unset = NA_character_),
  Sys.getenv("ATM_LEGACY_R_LIBRARY", unset = NA_character_)
)

if (!dir.exists(DIR_RAW)) {
  stop("Canonical immutable input directory not found: ", DIR_RAW, call. = FALSE)
}
if (!dir.exists(DIR_METADATA)) {
  stop("Canonical metadata directory not found: ", DIR_METADATA, call. = FALSE)
}
dir.create(DIR_PROCESSED, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_INTERMEDIATE, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_ANALYSIS_READY, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_ANALYSIS_READY_SOURCE, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_ANCILLARY, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_QC, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_QC_SOURCE, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_RESULTS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_R_LIBRARY, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_R_CACHE, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_USER_CACHE_DIR = DIR_R_CACHE)

is_readable_r_library <- function(path) {
  if (is.na(path) || !dir.exists(path)) return(FALSE)
  package_dirs <- list.dirs(path, full.names = TRUE, recursive = FALSE)
  if (length(package_dirs) == 0) return(FALSE)
  any(file.exists(file.path(package_dirs, "DESCRIPTION")))
}

HISTORICAL_R_LIBS <- unique(HISTORICAL_R_LIBS[
  !is.na(HISTORICAL_R_LIBS) &
    vapply(HISTORICAL_R_LIBS, is_readable_r_library, logical(1))
])
.libPaths(unique(c(DIR_R_LIBRARY, HISTORICAL_R_LIBS, .libPaths())))

BODY_FAT_OBESE_THRESHOLD <- 30

classify_body_fat <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x > BODY_FAT_OBESE_THRESHOLD ~ "Obese",
    TRUE ~ "Lean"
  )
}

p_load <- function(...) {
  pkgs <- as.character(substitute(list(...)))[-1]
  pkgs <- pkgs[nzchar(pkgs)]
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(lapply(pkgs, function(pkg) suppressPackageStartupMessages(library(pkg, character.only = TRUE))))
}

load_hallmark_gene_sets <- function(cache_file = Sys.getenv(
  "ATM_MSIGDB_HALLMARK_RDS", unset = ""
)) {
  cache_candidates <- character()
  if (nzchar(cache_file)) cache_candidates <- cache_file

  project_cache_dir <- file.path(DIR_R_CACHE, "R", "msigdbr")
  if (dir.exists(project_cache_dir)) {
    cache_candidates <- c(
      cache_candidates,
      list.files(
        project_cache_dir,
        pattern = "^msigdb[.].*[.]Hs[.]H[.]rds$",
        full.names = TRUE
      )
    )
  }
  cache_candidates <- unique(cache_candidates[file.exists(cache_candidates)])

  if (length(cache_candidates) > 0L) {
    cache_path <- normalizePath(cache_candidates[[1]], winslash = "/", mustWork = TRUE)
    hallmark <- readRDS(cache_path)
    message("Loaded Hallmark gene sets from local cache: ", cache_path)
  } else {
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      stop(
        "Hallmark gene sets are unavailable. Install msigdbr or set ",
        "ATM_MSIGDB_HALLMARK_RDS to an authorized local cache file.",
        call. = FALSE
      )
    }
    hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  }

  gene_column <- intersect(c("gene_symbol", "db_gene_symbol"), names(hallmark))[1]
  if (is.na(gene_column) || !"gs_name" %in% names(hallmark)) {
    stop("Hallmark gene-set data lack gs_name or gene-symbol columns.", call. = FALSE)
  }

  hallmark <- hallmark[
    !is.na(hallmark$gs_name) & !is.na(hallmark[[gene_column]]) &
      nzchar(hallmark$gs_name) & nzchar(hallmark[[gene_column]]),
    c("gs_name", gene_column),
    drop = FALSE
  ]
  names(hallmark)[2] <- "gene_symbol"
  hallmark <- unique(hallmark)
  split(hallmark$gene_symbol, hallmark$gs_name)
}

geom_point_rast <- function(..., raster.dpi = 300) {
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::geom_point_rast(..., raster.dpi = raster.dpi)
  } else {
    ggplot2::geom_point(...)
  }
}

geom_jitter_rast <- function(..., raster.dpi = 300) {
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::geom_jitter_rast(..., raster.dpi = raster.dpi)
  } else {
    ggplot2::geom_jitter(...)
  }
}

rasterise <- function(plot, dpi = 300, ...) {
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::rasterise(plot, dpi = dpi, ...)
  } else {
    plot
  }
}

partial_spearman_test <- function(x, y, covariates = NULL) {
  df <- data.frame(.x = x, .y = y, covariates, check.names = FALSE)
  df <- df[stats::complete.cases(df), , drop = FALSE]
  if (nrow(df) < 4 || stats::sd(df$.x) == 0 || stats::sd(df$.y) == 0) {
    return(list(estimate = NA_real_, p.value = NA_real_))
  }

  ranked <- as.data.frame(lapply(df, function(v) {
    if (is.numeric(v)) {
      rank(v, ties.method = "average")
    } else {
      rank(as.numeric(as.factor(v)), ties.method = "average")
    }
  }))

  covar_names <- setdiff(names(ranked), c(".x", ".y"))
  if (length(covar_names) > 0) {
    form <- stats::as.formula(paste(".value ~", paste(covar_names, collapse = " + ")))
    rx <- stats::resid(stats::lm(form, data = transform(ranked, .value = .x)))
    ry <- stats::resid(stats::lm(form, data = transform(ranked, .value = .y)))
  } else {
    rx <- ranked$.x
    ry <- ranked$.y
  }

  test <- stats::cor.test(rx, ry, method = "pearson")
  list(estimate = unname(test$estimate), p.value = test$p.value)
}

analysis_status_is_released <- function(status) {
  grepl("^Released", as.character(status))
}

read_analysis_matrix <- function(input_path) {
  if (!file.exists(input_path)) {
    stop("Analysis matrix not found: ", input_path, call. = FALSE)
  }

  extension <- tolower(tools::file_ext(input_path))
  if (extension == "rds") {
    output <- readRDS(input_path)
    if (!is.matrix(output) && !is.data.frame(output)) {
      stop("RDS analysis input must contain a matrix or data frame: ", input_path,
           call. = FALSE)
    }
    output <- as.matrix(output)
    suppressWarnings(storage.mode(output) <- "double")
    if (is.null(rownames(output)) || is.null(colnames(output))) {
      stop("RDS analysis input lacks row or column identifiers: ", input_path,
           call. = FALSE)
    }
    rownames(output) <- make.unique(as.character(rownames(output)))
    colnames(output) <- tolower(as.character(colnames(output)))
    return(output)
  }

  if (extension %in% c("csv", "gz")) {
    frame <- readr::read_csv(
      input_path,
      show_col_types = FALSE,
      progress = FALSE
    )
    if (ncol(frame) < 2L) {
      stop("CSV analysis input has no sample columns: ", input_path,
           call. = FALSE)
    }
    output <- data.matrix(frame[, -1, drop = FALSE])
    suppressWarnings(storage.mode(output) <- "double")
    rownames(output) <- make.unique(as.character(frame[[1]]))
    colnames(output) <- tolower(colnames(frame)[-1])
    return(output)
  }

  stop("Unsupported analysis-input extension: ", extension, call. = FALSE)
}

prepare_pca_matrix <- function(mat, dataset_name = "PCA matrix") {
  mat <- as.matrix(mat)
  suppressWarnings(storage.mode(mat) <- "numeric")

  finite_features <- apply(mat, 2, function(x) all(is.finite(x)))
  if (any(!finite_features)) {
    message("  [PCA] ", dataset_name, ": removed ", sum(!finite_features),
            " feature(s) with NA/Inf values.")
  }
  mat <- mat[, finite_features, drop = FALSE]

  feature_sd <- apply(mat, 2, stats::sd)
  variable_features <- is.finite(feature_sd) & feature_sd > 0
  if (any(!variable_features)) {
    message("  [PCA] ", dataset_name, ": removed ", sum(!variable_features),
            " zero-variance feature(s).")
  }
  mat <- mat[, variable_features, drop = FALSE]

  if (nrow(mat) < 2 || ncol(mat) < 2) {
    stop("Not enough finite, variable data for PCA in ", dataset_name, call. = FALSE)
  }

  mat
}
