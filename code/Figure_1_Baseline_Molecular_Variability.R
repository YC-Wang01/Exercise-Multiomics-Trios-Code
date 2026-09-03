# Figure 1D: Cell-submission-style baseline molecular variability.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.local_config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  file.path(getwd(), "Fig00_Config.R"),
  if (!is.na(.local_script)) {
    file.path(dirname(normalizePath(.local_script, winslash = "/", mustWork = FALSE)), "Fig00_Config.R")
  } else {
    NA_character_
  }
))
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (is.na(.local_config)) stop("Cannot find CellMetabolism_Transfer code/Fig00_Config.R.")
source(.local_config)
rm(.local_script, .local_config_candidates, .local_config)

p_load(readr, dplyr, tidyr, ggplot2, matrixStats, openxlsx, cowplot, digest)
options(cli.unicode = FALSE, pillar.unicode = FALSE, width = 180)

SEED <- 20260804L
OBSERVED_FRACTION <- 0.70
MIN_SAMPLES <- 3L
MAX_DISPLAY_POINTS <- 5000L
Y_LIMIT <- 150
GROUP_LEVELS <- c("Lean", "All participants", "Obese")

OUT_DIR <- path_project("fig01", "fig01D")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_COMBINED <- file.path(OUT_DIR, "Fig1D_Baseline_Molecular_Variability_CellStyle.pdf")
OUT_SOURCE <- file.path(OUT_DIR, "Fig1D_Baseline_Molecular_Variability_SourceData.xlsx")
OUT_FEATURES <- file.path(OUT_DIR, "Fig1D_Baseline_Molecular_Variability_FeatureData.rds")
OUT_LOG <- file.path(OUT_DIR, "Fig1D_Baseline_Molecular_Variability_Run.log")

DATASETS <- tibble::tribble(
  ~Dataset, ~Tissue, ~Layer, ~Label, ~Order, ~InputFile, ~InputScale, ~CVScale,
  "Serum_Proteomics", "Serum", "Proteomics", "Proteomics", 1L, "Analysis_Serum_Proteomics.csv", "log2 relative abundance", "2^(log2 relative abundance)",
  "Serum_Metabonomics", "Serum", "Metabonomics", "Metabonomics", 2L, "Analysis_Serum_Metabonomics.csv", "positive supplied scale", "supplied positive scale",
  "Adipose_Microarray", "Adipose", "Microarray", "Microarray", 3L, "Analysis_Adipose_Microarray.csv", "RMA log2", "2^(RMA log2)",
  "Adipose_Methylation", "Adipose", "Methylation", "Methylation", 4L, "Reporting_Adipose_Methylation_Beta.rds", "ssNoob Beta value", "Beta value",
  "Adipose_Proteomics", "Adipose", "Proteomics", "Proteomics", 5L, "Analysis_Adipose_Proteomics.csv", "log2 relative abundance", "2^(log2 relative abundance)",
  "Muscle_Microarray", "Muscle", "Microarray", "Microarray", 6L, "Analysis_Muscle_Microarray.csv", "RMA log2", "2^(RMA log2)",
  "Muscle_Methylation", "Muscle", "Methylation", "Methylation", 7L, "Reporting_Muscle_Methylation_Beta.rds", "ssNoob Beta value", "Beta value",
  "Muscle_Proteomics", "Muscle", "Proteomics", "Proteomics", 8L, "Analysis_Muscle_Proteomics.csv", "log2 relative abundance", "2^(log2 relative abundance)"
)

LAYER_COLORS <- c(
  Proteomics = "#8D6E9F",
  Metabonomics = "#D66B73",
  Methylation = "#52749A",
  Microarray = "#7F936C"
)

sample_sheet <- readr::read_csv(
  path_analysis_ready("Metadata", "Project_Sample_Sheet.csv"),
  show_col_types = FALSE,
  progress = FALSE
) %>%
  dplyr::mutate(
    SampleID = tolower(as.character(SampleID)),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort)
  )

clinical <- readr::read_csv(
  path_analysis_ready("Metadata", "Clinical_Observed_Primary.csv"),
  show_col_types = FALSE,
  progress = FALSE
) %>%
  dplyr::transmute(
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    BodyFat_percent = as.numeric(Fat_p_G),
    Adiposity = classify_body_fat(BodyFat_percent)
  ) %>%
  dplyr::distinct(Clinical_Subject_ID, .keep_all = TRUE)

to_cv_scale <- function(mat, input_scale) {
  output <- if (input_scale %in% c(
    "RMA log2",
    "log2 relative abundance"
  )) {
    2^mat
  } else {
    mat
  }
  output[!is.finite(output)] <- NA_real_
  output
}

calculate_cv <- function(mat, group_name, dataset_row, participant_n) {
  required_n <- max(MIN_SAMPLES, ceiling(participant_n * OBSERVED_FRACTION))
  observed_n <- matrixStats::rowCounts(is.finite(mat), value = TRUE)
  means <- matrixStats::rowMeans2(mat, na.rm = TRUE)
  sds <- matrixStats::rowSds(mat, na.rm = TRUE)
  keep <- observed_n >= required_n & is.finite(means) & means > 1e-6 & is.finite(sds)
  tibble::tibble(
    Group = group_name,
    Dataset = dataset_row$Dataset,
    Tissue = dataset_row$Tissue,
    Layer = dataset_row$Layer,
    Label = dataset_row$Label,
    Order = dataset_row$Order,
    FeatureID = rownames(mat)[keep],
    ObservedN = observed_n[keep],
    ParticipantN = participant_n,
    CVPercent = 100 * sds[keep] / means[keep]
  ) %>%
    dplyr::filter(is.finite(CVPercent), CVPercent >= 0)
}

analyse_dataset <- function(dataset_row) {
  message("Processing ", dataset_row$Dataset)
  mat <- read_analysis_matrix(path_analysis_ready(dataset_row$InputFile))
  mat <- to_cv_scale(mat, dataset_row$InputScale)

  metadata <- sample_sheet %>%
    dplyr::filter(
      Dataset == dataset_row$Dataset,
      AnalysisSet == "Primary",
      Is_Primary_Cohort,
      Timepoint == "Pre",
      SampleID %in% colnames(mat)
    ) %>%
    dplyr::left_join(clinical, by = "Clinical_Subject_ID") %>%
    dplyr::filter(!is.na(Adiposity)) %>%
    dplyr::distinct(SampleID, .keep_all = TRUE)

  output <- list()
  counter <- 0L
  for (group_name in GROUP_LEVELS) {
    selected <- if (group_name == "All participants") {
      metadata
    } else {
      metadata %>% dplyr::filter(Adiposity == group_name)
    }
    if (nrow(selected) < MIN_SAMPLES) next
    counter <- counter + 1L
    output[[counter]] <- calculate_cv(
      mat[, selected$SampleID, drop = FALSE],
      group_name,
      dataset_row,
      nrow(selected)
    )
  }
  rm(mat)
  invisible(gc())
  dplyr::bind_rows(output)
}

feature_data <- dplyr::bind_rows(lapply(seq_len(nrow(DATASETS)), function(i) {
  analyse_dataset(DATASETS[i, ])
})) %>%
  dplyr::mutate(
    Group = factor(Group, levels = GROUP_LEVELS),
    Dataset = factor(Dataset, levels = DATASETS$Dataset),
    Layer = factor(Layer, levels = names(LAYER_COLORS))
  ) %>%
  dplyr::arrange(Group, Dataset, FeatureID)

if (nrow(feature_data) == 0L) stop("No eligible feature-level CV values were generated.")

box_stats <- feature_data %>%
  dplyr::group_by(Group, Dataset, Tissue, Layer, Label, Order, ParticipantN) %>%
  dplyr::summarise(
    FeatureN = dplyr::n(),
    MeanCV = mean(CVPercent),
    MedianCV = median(CVPercent),
    SDCV = stats::sd(CVPercent),
    Q1 = stats::quantile(CVPercent, 0.25),
    Q3 = stats::quantile(CVPercent, 0.75),
    P95 = stats::quantile(CVPercent, 0.95),
    P99 = stats::quantile(CVPercent, 0.99),
    Maximum = max(CVPercent),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Group, Order)

set.seed(SEED)
plot_points <- feature_data %>%
  dplyr::group_by(Group, Dataset) %>%
  dplyr::group_modify(function(.x, .y) {
    dplyr::slice_sample(.x, n = min(nrow(.x), MAX_DISPLAY_POINTS))
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(DisplayCV = pmin(CVPercent, Y_LIMIT))

make_plot <- function(groups, title = NULL) {
  points <- plot_points %>% dplyr::filter(as.character(Group) %in% groups)
  plot <- ggplot2::ggplot(
    points,
    ggplot2::aes(x = Dataset, y = DisplayCV, color = Layer)
  ) +
    geom_jitter_rast(width = 0.24, height = 0, size = 0.32, alpha = 0.26, raster.dpi = 450) +
    ggplot2::geom_boxplot(
      fill = NA,
      width = 0.56,
      linewidth = 0.62,
      outlier.shape = NA
    ) +
    ggplot2::scale_color_manual(values = LAYER_COLORS, drop = FALSE) +
    ggplot2::scale_x_discrete(
      labels = setNames(
        c("", "Serum", "", "Adipose", "", "", "Muscle", ""),
        DATASETS$Dataset
      ),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 25, 50, 75, 100, 125, 150),
      limits = c(0, Y_LIMIT),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = "Coefficient of variation (CV, %)",
      color = NULL
    ) +
    ggplot2::theme_classic(base_family = "Arial", base_size = 8.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 9.5, hjust = 0),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 9.5, hjust = 0),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.65),
      axis.text.x = ggplot2::element_text(color = "black", size = 8.2, face = "bold"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(color = "black", size = 7.5),
      axis.title.y = ggplot2::element_text(size = 8.4, margin = ggplot2::margin(r = 4)),
      strip.text.y.right = ggplot2::element_text(angle = 0, face = "bold", size = 9.5),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 7.5),
      plot.margin = ggplot2::margin(4, 5, 4, 5)
    )
  if (length(groups) > 1L) {
    plot <- plot + ggplot2::facet_grid(rows = ggplot2::vars(Group))
  }
  plot
}

for (group_name in GROUP_LEVELS) {
  slug <- gsub(" ", "_", group_name)
  plot <- make_plot(group_name, group_name)
  pdf_path <- file.path(OUT_DIR, paste0("Fig1D_Baseline_Molecular_Variability_", slug, ".pdf"))
  ggplot2::ggsave(pdf_path, plot, width = 160, height = 72, units = "mm", device = grDevices::cairo_pdf)
  ggplot2::ggsave(sub("\\.pdf$", ".png", pdf_path), plot, width = 160, height = 72, units = "mm", dpi = 450, bg = "white")
}

combined_plot <- make_plot(GROUP_LEVELS)
ggplot2::ggsave(OUT_COMBINED, combined_plot, width = 160, height = 178, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(sub("\\.pdf$", ".png", OUT_COMBINED), combined_plot, width = 160, height = 178, units = "mm", dpi = 450, bg = "white")

sample_counts <- feature_data %>%
  dplyr::distinct(Group, Dataset, Tissue, Layer, ParticipantN) %>%
  dplyr::arrange(Group, Dataset)

parameters <- tibble::tibble(
  Parameter = c(
    "Figure definition", "Baseline", "Groups", "Adiposity threshold",
    "CV formula", "Feature observed fraction", "Outlier policy",
    "Methylation scale", "Serum proteomics caveat", "Point display",
    "Boxplot data", "Displayed range", "Inference"
  ),
  Value = c(
    "Cell-submission Figure 1D left-column layout",
    "Pre only",
    "All participants, Lean, Obese (individual adiposity; All is not Mixed families)",
    "Obese: body fat >30%; Lean: body fat <=30%",
    "100 x sample SD / sample mean, feature-wise across participants",
    as.character(OBSERVED_FRACTION),
    "No participant excluded; PCA-extreme flags are review-only",
    "ssNoob Beta values",
    "Log2 expression and relative-abundance values were back-transformed before CV calculation",
    paste0("Deterministic random sample of up to ", MAX_DISPLAY_POINTS, " features per dataset and group"),
    "All eligible features",
    paste0("0-", Y_LIMIT, "%; larger values retained in FeatureData.rds"),
    "Descriptive only; no Lean-Obese hypothesis test"
  )
)

input_manifest <- DATASETS %>%
  dplyr::mutate(
    Path = vapply(InputFile, function(x) normalizePath(path_analysis_ready(x), winslash = "/", mustWork = TRUE), character(1)),
    SHA256 = vapply(Path, function(x) digest::digest(x, algo = "sha256", file = TRUE), character(1))
  )

workbook <- openxlsx::createWorkbook()
write_sheet <- function(name, data) {
  openxlsx::addWorksheet(workbook, name)
  openxlsx::writeData(workbook, name, data)
  openxlsx::freezePane(workbook, name, firstRow = TRUE)
  openxlsx::setColWidths(workbook, name, cols = seq_len(ncol(data)), widths = "auto")
}
write_sheet("Parameters", parameters)
write_sheet("Inputs", input_manifest)
write_sheet("SampleCounts", sample_counts)
write_sheet("Summary", box_stats)
openxlsx::saveWorkbook(workbook, OUT_SOURCE, overwrite = TRUE)
saveRDS(feature_data, OUT_FEATURES, compress = "gzip")

log_lines <- c(
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R version: ", R.version.string),
  "Script: code/Figure_1_Baseline_Molecular_Variability.R",
  paste0("Output: ", normalizePath(OUT_COMBINED, winslash = "/", mustWork = FALSE)),
  paste0("Feature rows: ", nrow(feature_data)),
  paste0("Groups: ", paste(GROUP_LEVELS, collapse = ", ")),
  paste0("Observed fraction: ", OBSERVED_FRACTION),
  "No participant-level automatic outlier exclusion was applied.",
  "Serum proteomics CV uses back-transformed log2 sample-median-centred relative abundance."
)
writeLines(log_lines, OUT_LOG, useBytes = TRUE)

message("Figure 1D baseline CV outputs written to: ", OUT_DIR)
