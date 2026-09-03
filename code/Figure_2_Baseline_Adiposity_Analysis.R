# Figure 2A/C/D/F baseline adiposity analysis across eight omics platforms.
#
# Primary model: Feature ~ Adiposity + Role, with FamilyID correlation blocking.
# Figure 2A and Figure 2B preserve the submitted nominal-P visual language;
# BH-FDR values are retained in all machine-readable outputs.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.local_config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  file.path(getwd(), "Fig00_Config.R"),
  if (!is.na(.local_script)) {
    file.path(
      dirname(normalizePath(.local_script, winslash = "/", mustWork = FALSE)),
      "Fig00_Config.R"
    )
  } else {
    NA_character_
  }
))
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (is.na(.local_config)) stop("Cannot find CellMetabolism_Transfer code/Fig00_Config.R.")
source(.local_config)
rm(.local_script, .local_config_candidates, .local_config)

p_load(
  readr,
  dplyr,
  tidyr,
  tibble,
  ggplot2,
  ggrepel,
  limma,
  digest,
  png,
  ragg,
  patchwork,
  writexl
)

set.seed(20260804)
select <- dplyr::select

analysis_id <- "Fig02AB_Baseline_Adiposity_2026-08-09"
release_id <- "2026-08-09_v6"
nominal_p_threshold <- 0.05
fdr_threshold <- 0.05
minimum_feature_observations <- 10L
volcano_y_cap <- 8
labels_per_side <- 3L
role_levels <- c("Daughter", "Mother", "Father")
model_description <- "Feature ~ Adiposity + Role; FamilyID duplicateCorrelation block"

direction_colors <- c(
  "Not significant" = "#D7DEE2",
  "Higher in Lean" = "#3FA0A8",
  "Higher in Obese" = "#F59646"
)

fig02_root <- path_project("fig02")
fig02a_dir <- file.path(fig02_root, "fig02A")
fig02c_dir <- file.path(fig02_root, "fig02C")
fig02d_dir <- file.path(fig02_root, "fig02D")
fig02f_dir <- file.path(fig02_root, "fig02F")
dir.create(fig02a_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig02c_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig02d_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig02f_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(DIR_R_CACHE, "Fig02AB")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

metadata_dir <- path_analysis_ready("Metadata")
clinical_file <- file.path(metadata_dir, "Clinical_Observed_Primary.csv")
sample_sheet_file <- file.path(metadata_dir, "Project_Sample_Sheet.csv")
manifest_file <- file.path(metadata_dir, "Matrix_Manifest.csv")

dataset_registry <- tibble::tribble(
  ~Dataset, ~InputFile, ~DisplayLabel, ~PanelTitle, ~DisplayOrder,
  "Serum_Proteomics", "Analysis_Serum_Proteomics.csv", "Serum\n(Proteomics)", "Serum proteomics", 1L,
  "Serum_Metabonomics", "Analysis_Serum_Metabonomics.csv", "Serum\n(Metabolomics)", "Serum metabolomics", 2L,
  "Muscle_Proteomics", "Analysis_Muscle_Proteomics.csv", "Muscle\n(Proteomics)", "Muscle proteomics", 3L,
  "Muscle_Microarray", "Analysis_Muscle_Microarray.csv", "Muscle\n(Transcriptomics)", "Muscle transcriptomics", 4L,
  "Muscle_Methylation", "Analysis_Muscle_Methylation_MValue.rds", "Muscle\n(DNA methylation)", "Muscle DNA methylation", 5L,
  "Adipose_Proteomics", "Analysis_Adipose_Proteomics.csv", "Adipose\n(Proteomics)", "Adipose proteomics", 6L,
  "Adipose_Microarray", "Analysis_Adipose_Microarray.csv", "Adipose\n(Transcriptomics)", "Adipose transcriptomics", 7L,
  "Adipose_Methylation", "Analysis_Adipose_Methylation_MValue.rds", "Adipose\n(DNA methylation)", "Adipose DNA methylation", 8L
)

manifest <- readr::read_csv(manifest_file, show_col_types = FALSE) %>%
  filter(MatrixRole == "Analysis")
dataset_registry <- dataset_registry %>%
  left_join(
    manifest %>% select(Dataset, UseStatus, SHA256, Features, Samples, MissingRate),
    by = "Dataset"
  )
if (any(is.na(dataset_registry$UseStatus))) {
  stop("Dataset registry does not match Matrix_Manifest.csv.", call. = FALSE)
}
if (any(!analysis_status_is_released(dataset_registry$UseStatus))) {
  stop("At least one requested matrix is not released.", call. = FALSE)
}

clinical <- readr::read_csv(clinical_file, show_col_types = FALSE) %>%
  transmute(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = factor(Role, levels = role_levels),
    FatPercent = suppressWarnings(as.numeric(Fat_p_G)),
    Adiposity = factor(classify_body_fat(FatPercent), levels = c("Lean", "Obese"))
  ) %>%
  filter(!is.na(FamilyID), !is.na(Role), is.finite(FatPercent), !is.na(Adiposity)) %>%
  distinct(FamilyID, Role, .keep_all = TRUE)

if (
  nrow(clinical) != 81L ||
    dplyr::n_distinct(clinical$FamilyID) != 27L ||
    any(table(clinical$FamilyID) != 3L)
) {
  stop("Clinical metadata is not the locked 27-family cohort.", call. = FALSE)
}

sample_sheet <- readr::read_csv(sample_sheet_file, show_col_types = FALSE) %>%
  mutate(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    SampleID = tolower(as.character(SampleID)),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort)
  ) %>%
  filter(AnalysisSet == "Primary", Is_Primary_Cohort)

prepare_baseline <- function(dataset_name, expression_matrix) {
  metadata <- sample_sheet %>%
    filter(
      Dataset == dataset_name,
      Timepoint == "Pre",
      SampleID %in% colnames(expression_matrix)
    ) %>%
    inner_join(clinical, by = c("FamilyID", "Clinical_Subject_ID", "Role")) %>%
    filter(!is.na(Adiposity), !is.na(Role)) %>%
    distinct(Clinical_Subject_ID, .keep_all = TRUE) %>%
    arrange(FamilyID, Role)

  matrix_index <- match(metadata$SampleID, tolower(colnames(expression_matrix)))
  if (anyNA(matrix_index)) stop("Sample alignment failed for ", dataset_name)
  expression_matrix <- expression_matrix[, matrix_index, drop = FALSE]
  colnames(expression_matrix) <- metadata$Clinical_Subject_ID

  keep <- rowSums(is.finite(expression_matrix)) >= minimum_feature_observations
  variable <- apply(expression_matrix, 1L, function(values) {
    observed <- values[is.finite(values)]
    length(observed) >= minimum_feature_observations &&
      is.finite(stats::sd(observed)) && stats::sd(observed) > 0
  })

  list(
    matrix = expression_matrix[keep & variable, , drop = FALSE],
    metadata = metadata
  )
}

fit_baseline <- function(expression_matrix, metadata, dataset_name) {
  design <- stats::model.matrix(~ Adiposity + Role, data = metadata)
  if (qr(design)$rank != ncol(design)) {
    stop("Rank-deficient design for ", dataset_name)
  }
  if (nrow(design) - ncol(design) < 5L) {
    stop("Insufficient residual degrees of freedom for ", dataset_name)
  }

  correlation_fit <- limma::duplicateCorrelation(
    expression_matrix,
    design,
    block = metadata$FamilyID
  )
  family_correlation <- correlation_fit$consensus.correlation
  if (!is.finite(family_correlation)) {
    stop("Non-finite family correlation for ", dataset_name)
  }

  fitted <- limma::lmFit(
    expression_matrix,
    design,
    block = metadata$FamilyID,
    correlation = family_correlation
  )
  fitted <- limma::eBayes(fitted, trend = TRUE, robust = TRUE)
  coefficient_index <- match("AdiposityObese", colnames(design))
  result <- limma::topTable(
    fitted,
    coef = coefficient_index,
    number = Inf,
    sort.by = "none"
  )

  tibble(
    Dataset = dataset_name,
    FeatureID = rownames(result),
    Effect = result$logFC,
    ModeratedT = result$t,
    P_Value = result$P.Value,
    BH_FDR = result$adj.P.Val,
    AverageExpression = result$AveExpr,
    ObservedN = rowSums(is.finite(expression_matrix)),
    Samples = nrow(metadata),
    Families = dplyr::n_distinct(metadata$FamilyID),
    LeanN = sum(metadata$Adiposity == "Lean"),
    ObeseN = sum(metadata$Adiposity == "Obese"),
    FamilyCorrelation = family_correlation
  )
}

read_annotation <- function(dataset_name) {
  if (dataset_name %in% c("Adipose_Microarray", "Muscle_Microarray")) {
    readr::read_csv(
      file.path(metadata_dir, paste0("Feature_Annotation_", dataset_name, ".csv")),
      show_col_types = FALSE
    ) %>%
      transmute(
        FeatureID = as.character(ProbeSetID),
        FeatureLabel = if_else(!is.na(SYMBOL) & nzchar(SYMBOL), SYMBOL, ProbeSetID)
      ) %>%
      distinct(FeatureID, .keep_all = TRUE)
  } else if (dataset_name %in% c("Adipose_Proteomics", "Muscle_Proteomics")) {
    readr::read_csv(
      file.path(metadata_dir, paste0("Feature_Annotation_", dataset_name, ".csv")),
      show_col_types = FALSE
    ) %>%
      transmute(
        FeatureID = as.character(FeatureID),
        FeatureLabel = if_else(
          !is.na(IntegrationGene) & nzchar(IntegrationGene),
          IntegrationGene,
          if_else(!is.na(PG.Genes) & nzchar(PG.Genes), PG.Genes, FeatureID)
        )
      ) %>%
      distinct(FeatureID, .keep_all = TRUE)
  } else if (dataset_name == "Serum_Proteomics") {
    readr::read_csv(
      file.path(metadata_dir, "Feature_Dictionary_Serum_Proteomics.csv"),
      show_col_types = FALSE
    ) %>%
      transmute(
        FeatureID = as.character(FeatureID),
        FeatureLabel = if_else(!is.na(GeneSymbol) & nzchar(GeneSymbol), GeneSymbol, FeatureID)
      ) %>%
      distinct(FeatureID, .keep_all = TRUE)
  } else if (dataset_name == "Serum_Metabonomics") {
    readr::read_csv(
      file.path(metadata_dir, "Feature_Dictionary_Serum_Metabonomics.csv"),
      show_col_types = FALSE
    ) %>%
      transmute(
        FeatureID = as.character(FeatureID),
        FeatureLabel = if_else(!is.na(LongName) & nzchar(LongName), LongName, FeatureID)
      ) %>%
      distinct(FeatureID, .keep_all = TRUE)
  } else {
    tibble(FeatureID = character(), FeatureLabel = character())
  }
}

result_list <- list()
summary_list <- list()
input_manifest <- list()
cache_manifest <- list()

for (dataset_index in seq_len(nrow(dataset_registry))) {
  dataset_row <- dataset_registry[dataset_index, ]
  input_path <- path_analysis_ready(dataset_row$InputFile)
  input_hash <- digest::digest(file = input_path, algo = "sha256", serialize = FALSE)
  if (!identical(tolower(input_hash), tolower(dataset_row$SHA256))) {
    stop("Input hash differs from Matrix_Manifest.csv for ", dataset_row$Dataset)
  }

  cache_file <- file.path(
    cache_dir,
    paste0("FullResults_", dataset_row$Dataset, "_CategoricalRole.rds")
  )
  cache_reused <- FALSE
  if (file.exists(cache_file)) {
    cache <- tryCatch(readRDS(cache_file), error = function(error) NULL)
    cache_reused <- !is.null(cache) &&
      cache$ReleaseID %in% c(release_id, "2026-08-01_v5") &&
      identical(tolower(cache$InputSHA256), tolower(input_hash)) &&
      identical(cache$Model, model_description) &&
      identical(cache$SampleWeights, "None") &&
      is.data.frame(cache$Results)
  }

  if (cache_reused) {
    message("Loading validated cache for ", dataset_row$Dataset)
    result <- tibble::as_tibble(cache$Results)
  } else {
    message("Reading and fitting ", dataset_row$Dataset)
    expression_matrix <- read_analysis_matrix(input_path)
    baseline <- prepare_baseline(dataset_row$Dataset, expression_matrix)
    result <- fit_baseline(baseline$matrix, baseline$metadata, dataset_row$Dataset)
    annotation <- read_annotation(dataset_row$Dataset)
    if (nrow(annotation) > 0L) {
      result <- result %>% left_join(annotation, by = "FeatureID")
    } else {
      result$FeatureLabel <- result$FeatureID
    }
    saveRDS(
      list(
        AnalysisID = analysis_id,
        ReleaseID = release_id,
        InputSHA256 = input_hash,
        Model = model_description,
        SampleWeights = "None",
        Results = result
      ),
      cache_file,
      compress = "gzip"
    )
    rm(expression_matrix, baseline, annotation)
  }

  if (!"FeatureLabel" %in% names(result)) result$FeatureLabel <- result$FeatureID
  if (dataset_row$Dataset == "Serum_Metabonomics") {
    result$FeatureLabel <- result$FeatureID
  }

  result <- result %>%
    mutate(
      FeatureLabel = if_else(
        is.na(FeatureLabel) | !nzchar(FeatureLabel),
        FeatureID,
        FeatureLabel
      ),
      NominalDirection = case_when(
        P_Value < nominal_p_threshold & Effect < 0 ~ "Higher in Lean",
        P_Value < nominal_p_threshold & Effect > 0 ~ "Higher in Obese",
        TRUE ~ "Not significant"
      ),
      FDRSignificant = is.finite(BH_FDR) & BH_FDR < fdr_threshold
    )

  result_list[[dataset_row$Dataset]] <- result
  summary_list[[dataset_row$Dataset]] <- tibble(
    Dataset = dataset_row$Dataset,
    DisplayLabel = dataset_row$DisplayLabel,
    PanelTitle = dataset_row$PanelTitle,
    DisplayOrder = dataset_row$DisplayOrder,
    Samples = unique(result$Samples),
    Families = unique(result$Families),
    LeanN = unique(result$LeanN),
    ObeseN = unique(result$ObeseN),
    FeaturesTested = nrow(result),
    HigherLeanNominal = sum(result$NominalDirection == "Higher in Lean"),
    HigherObeseNominal = sum(result$NominalDirection == "Higher in Obese"),
    SignificantFDR05 = sum(result$FDRSignificant, na.rm = TRUE),
    MinimumP = min(result$P_Value, na.rm = TRUE),
    MinimumBH_FDR = min(result$BH_FDR, na.rm = TRUE),
    FamilyCorrelation = unique(result$FamilyCorrelation)
  )
  input_manifest[[dataset_row$Dataset]] <- tibble(
    Dataset = dataset_row$Dataset,
    InputFile = dataset_row$InputFile,
    InputSHA256 = input_hash,
    FeaturesInRelease = dataset_row$Features,
    SamplesInRelease = dataset_row$Samples
  )
  cache_manifest[[dataset_row$Dataset]] <- tibble(
    Dataset = dataset_row$Dataset,
    CacheReused = cache_reused,
    CacheFile = basename(cache_file),
    CacheSHA256 = digest::digest(file = cache_file, algo = "sha256", serialize = FALSE)
  )

  rm(result)
  invisible(gc())
}

summary_table <- bind_rows(summary_list) %>%
  mutate(
    PercentHigherLeanNominal = 100 * HigherLeanNominal / FeaturesTested,
    PercentHigherObeseNominal = 100 * HigherObeseNominal / FeaturesTested
  ) %>%
  arrange(DisplayOrder)
input_manifest <- bind_rows(input_manifest) %>% arrange(match(Dataset, dataset_registry$Dataset))
cache_manifest <- bind_rows(cache_manifest) %>% arrange(match(Dataset, dataset_registry$Dataset))

# Figure 2A: submitted mirrored landscape, recalculated from the v5 release.
landscape <- summary_table %>%
  select(
    Dataset,
    DisplayLabel,
    DisplayOrder,
    FeaturesTested,
    HigherLeanNominal,
    HigherObeseNominal,
    PercentHigherLeanNominal,
    PercentHigherObeseNominal
  ) %>%
  pivot_longer(
    cols = c(PercentHigherLeanNominal, PercentHigherObeseNominal),
    names_to = "DirectionCode",
    values_to = "Percentage"
  ) %>%
  mutate(
    Direction = recode(
      DirectionCode,
      PercentHigherLeanNominal = "Higher in Lean",
      PercentHigherObeseNominal = "Higher in Obese"
    ),
    NominalCount = if_else(
      Direction == "Higher in Lean",
      HigherLeanNominal,
      HigherObeseNominal
    ),
    PlotValue = if_else(Direction == "Higher in Lean", -Percentage, Percentage),
    DisplayLabel = factor(DisplayLabel, levels = rev(dataset_registry$DisplayLabel))
  )

landscape_limit <- max(abs(landscape$PlotValue), na.rm = TRUE)
if (!is.finite(landscape_limit) || landscape_limit <= 0) landscape_limit <- 1
landscape_limit <- landscape_limit * 1.18

fig2a <- ggplot(landscape, aes(x = PlotValue, y = DisplayLabel, fill = Direction)) +
  geom_col(width = 0.58, color = "#26343C", linewidth = 0.28) +
  geom_vline(xintercept = 0, color = "#26343C", linewidth = 0.55) +
  scale_fill_manual(values = direction_colors[c("Higher in Lean", "Higher in Obese")]) +
  scale_x_continuous(
    limits = c(-landscape_limit, landscape_limit),
    labels = function(x) abs(x),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Global landscape of baseline differences",
    x = "Nominally associated features (%)",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 9.5, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 9.5),
    axis.text.y = element_text(size = 8.5, color = "#26343C"),
    axis.text.x = element_text(size = 8, color = "#26343C"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = grid::unit(0.42, "cm"),
    plot.margin = margin(7, 10, 7, 7)
  )

ggsave(
  file.path(fig02a_dir, "Fig02A_Baseline_EightOmics_Landscape.pdf"),
  fig2a,
  width = 7.3,
  height = 5.4,
  device = grDevices::cairo_pdf
)
ggsave(
  file.path(fig02a_dir, "Fig02A_Baseline_EightOmics_Landscape.png"),
  fig2a,
  width = 7.3,
  height = 5.4,
  dpi = 450,
  bg = "white"
)

parameters <- tibble(
  Parameter = c(
    "AnalysisID",
    "ReleaseID",
    "Baseline",
    "AdiposityDefinition",
    "Model",
    "FamilyDependence",
    "SampleWeights",
    "VarianceModeration",
    "Figure2AThreshold",
    "Figure2BThreshold",
    "MultipleTesting",
    "VolcanoYCap"
  ),
  Value = c(
    analysis_id,
    release_id,
    "Pre",
    "Obese if body fat percentage >30%; Lean otherwise",
    "Feature ~ Adiposity + Role",
    "limma duplicateCorrelation block = FamilyID",
    "None",
    "limma eBayes trend=TRUE, robust=TRUE",
    "Nominal P <0.05; direction from the Obese-Lean effect",
    "Y=-log10(nominal P); color denotes nominal P <0.05 and effect direction",
    "Benjamini-Hochberg within each omics platform; values retained in source data",
    as.character(volcano_y_cap)
  )
)

writexl::write_xlsx(
  list(
    Omics_Summary = summary_table,
    Plot_Data = landscape,
    Input_Manifest = input_manifest,
    Cache_Manifest = cache_manifest,
    Parameters = parameters
  ),
  file.path(fig02a_dir, "Fig02A_SourceData.xlsx")
)

select_volcano_labels <- function(data, label_map = NULL) {
  if (!is.null(label_map)) {
    data <- data %>%
      select(-FeatureLabel) %>%
      inner_join(
        label_map %>% select(
          FeatureID,
          GeneLabel,
          UCSC_RefGene_Name,
          UCSC_RefGene_Group,
          Chromosome,
          Position_hg19,
          AnnotationPackage,
          AnnotationPackageVersion
        ),
        by = "FeatureID"
      ) %>%
      mutate(FeatureLabel = GeneLabel)
  }
  candidates <- data %>%
    filter(
      is.finite(Effect),
      is.finite(P_Value),
      P_Value < nominal_p_threshold,
      NominalDirection != "Not significant",
      !is.na(FeatureLabel),
      nzchar(FeatureLabel)
    ) %>%
    mutate(EffectSide = if_else(Effect < 0, "Higher in Lean", "Higher in Obese")) %>%
    group_by(EffectSide, FeatureLabel) %>%
    arrange(desc(FDRSignificant), P_Value, desc(abs(Effect)), .by_group = TRUE) %>%
    slice_head(n = 1L) %>%
    ungroup() %>%
    group_by(EffectSide) %>%
    arrange(desc(FDRSignificant), P_Value, desc(abs(Effect)), .by_group = TRUE) %>%
    slice_head(n = labels_per_side) %>%
    ungroup()
  candidates
}

make_raster_volcano <- function(data, label_data, title_text) {
  data <- data %>%
    filter(is.finite(Effect), is.finite(P_Value)) %>%
    mutate(
      PlotY = pmin(-log10(pmax(P_Value, 1e-300)), volcano_y_cap),
      PlotDirection = factor(
        NominalDirection,
        levels = c("Not significant", "Higher in Lean", "Higher in Obese")
      )
    )
  x_limit <- max(abs(data$Effect), na.rm = TRUE)
  if (!is.finite(x_limit) || x_limit <= 0) x_limit <- 1
  x_limit <- x_limit * 1.05

  raster_file <- tempfile(fileext = ".png")
  ragg::agg_png(
    raster_file,
    width = 1300,
    height = 1000,
    res = 220,
    background = "transparent"
  )
  graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  graphics::plot.new()
  graphics::plot.window(xlim = c(-x_limit, x_limit), ylim = c(0, volcano_y_cap))
  for (direction in c("Not significant", "Higher in Lean", "Higher in Obese")) {
    subset_data <- data[data$NominalDirection == direction, , drop = FALSE]
    if (nrow(subset_data) == 0L) next
    point_color <- grDevices::adjustcolor(
      direction_colors[[direction]],
      alpha.f = if (direction == "Not significant") 0.40 else 0.78
    )
    graphics::points(
      subset_data$Effect,
      subset_data$PlotY,
      pch = 16,
      cex = if (direction == "Not significant") 0.24 else 0.40,
      col = point_color
    )
  }
  grDevices::dev.off()
  raster_image <- png::readPNG(raster_file)
  unlink(raster_file)

  label_data <- label_data %>%
    mutate(PlotY = pmin(-log10(pmax(P_Value, 1e-300)), volcano_y_cap))
  legend_data <- tibble(
    Effect = 0,
    PlotY = 0,
    PlotDirection = factor(
      c("Not significant", "Higher in Lean", "Higher in Obese"),
      levels = c("Not significant", "Higher in Lean", "Higher in Obese")
    )
  )

  ggplot() +
    annotation_raster(
      raster_image,
      xmin = -x_limit,
      xmax = x_limit,
      ymin = 0,
      ymax = volcano_y_cap,
      interpolate = FALSE
    ) +
    geom_point(
      data = legend_data,
      aes(x = Effect, y = PlotY, color = PlotDirection),
      alpha = 0,
      size = 0.01,
      show.legend = TRUE
    ) +
    geom_hline(
      yintercept = -log10(nominal_p_threshold),
      linetype = "22",
      linewidth = 0.42,
      color = "#62717A"
    ) +
    geom_vline(xintercept = 0, linewidth = 0.32, color = "#62717A") +
    ggrepel::geom_text_repel(
      data = label_data,
      aes(x = Effect, y = PlotY, label = FeatureLabel),
      size = 2.15,
      family = "Arial",
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.color = "#68747C",
      segment.size = 0.24,
      box.padding = 0.20,
      point.padding = 0.10,
      seed = 20260804
    ) +
    scale_color_manual(
      values = direction_colors,
      breaks = c("Not significant", "Higher in Lean", "Higher in Obese"),
      name = NULL,
      guide = guide_legend(override.aes = list(alpha = 1, size = 2.4))
    ) +
    scale_y_continuous(
      breaks = seq(0, volcano_y_cap, by = 2),
      limits = c(0, volcano_y_cap),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(
      xlim = c(-x_limit, x_limit),
      ylim = c(0, volcano_y_cap),
      expand = FALSE
    ) +
    labs(
      title = title_text,
      x = "Effect estimate (Obese - Lean)",
      y = expression(-log[10]~"P")
    ) +
    theme_classic(base_size = 9.5, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = 10.5, hjust = 0.5),
      axis.title = element_text(size = 9.2),
      axis.text = element_text(size = 8.5, color = "#26343C"),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.margin = margin(6, 6, 6, 6)
    )
}

methylation_annotation_file <- file.path(
  fig02a_dir,
  "Fig02A_Methylation_Label_Annotation.csv.gz"
)
statistics_file <- file.path(fig02a_dir, "Fig02A_Full_Statistics.rds")
saveRDS(
  list(
    AnalysisID = analysis_id,
    ReleaseID = release_id,
    Model = "Feature ~ Adiposity + Role; FamilyID duplicateCorrelation block",
    Results = result_list,
    Summary = summary_table,
    InputManifest = input_manifest,
    CacheManifest = cache_manifest,
    Parameters = parameters
  ),
  statistics_file,
  compress = "gzip"
)
if (!file.exists(methylation_annotation_file)) {
  stop(
    "Missing controlled EPIC label annotation: ", annotation_file,
    ". Provide the documented probe-to-gene annotation before running this script."
  )
}
methylation_annotation <- readr::read_csv(
  methylation_annotation_file,
  show_col_types = FALSE
)
if (!all(c("Dataset", "FeatureID", "GeneLabel") %in% names(methylation_annotation))) {
  stop("EPIC label annotation lacks required columns.")
}

label_list <- lapply(names(result_list), function(dataset_name) {
  if (dataset_name %in% c("Muscle_Methylation", "Adipose_Methylation")) {
    select_volcano_labels(
      result_list[[dataset_name]],
      methylation_annotation %>% filter(Dataset == dataset_name)
    )
  } else {
    select_volcano_labels(result_list[[dataset_name]])
  }
})
names(label_list) <- names(result_list)

volcano_plots <- lapply(seq_len(nrow(dataset_registry)), function(index) {
  dataset_name <- dataset_registry$Dataset[index]
  make_raster_volcano(
    result_list[[dataset_name]],
    label_list[[dataset_name]],
    dataset_registry$PanelTitle[index]
  )
})

names(volcano_plots) <- dataset_registry$Dataset

assemble_tissue_volcanoes <- function(dataset_names, title_text) {
  patchwork::wrap_plots(volcano_plots[dataset_names], nrow = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = title_text,
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, family = "Arial")
      )
    ) &
    theme(legend.position = "right")
}

fig2b_serum <- assemble_tissue_volcanoes(
  c("Serum_Proteomics", "Serum_Metabonomics"),
  "Baseline serum differences between lean and obese participants"
)
fig2b_adipose <- assemble_tissue_volcanoes(
  c("Adipose_Proteomics", "Adipose_Microarray", "Adipose_Methylation"),
  "Baseline adipose differences between lean and obese participants"
)
fig2b_muscle <- assemble_tissue_volcanoes(
  c("Muscle_Proteomics", "Muscle_Microarray", "Muscle_Methylation"),
  "Baseline muscle differences between lean and obese participants"
)

tissue_figures <- list(
  Serum = list(
    plot = fig2b_serum,
    width = 7.4,
    output_dir = fig02c_dir,
    output_stem = "Fig02C_Serum_Omics_Volcanoes"
  ),
  Adipose = list(
    plot = fig2b_adipose,
    width = 10.4,
    output_dir = fig02d_dir,
    output_stem = "Fig02D_Adipose_Omics_Volcanoes"
  ),
  Muscle = list(
    plot = fig2b_muscle,
    width = 10.4,
    output_dir = fig02f_dir,
    output_stem = "Fig02F_Muscle_Omics_Volcanoes"
  )
)
for (tissue_name in names(tissue_figures)) {
  figure <- tissue_figures[[tissue_name]]
  ggsave(
    file.path(figure$output_dir, paste0(figure$output_stem, ".pdf")),
    figure$plot,
    width = figure$width,
    height = 3.9,
    device = grDevices::cairo_pdf
  )
  ggsave(
    file.path(figure$output_dir, paste0(figure$output_stem, ".png")),
    figure$plot,
    width = figure$width,
    height = 3.9,
    dpi = 450,
    bg = "white"
  )
}

all_results <- bind_rows(result_list)
top_results <- all_results %>%
  group_by(Dataset) %>%
  arrange(P_Value, BH_FDR, desc(abs(Effect)), .by_group = TRUE) %>%
  slice_head(n = 100L) %>%
  ungroup()
label_table <- bind_rows(label_list)

saveRDS(
  list(
    AnalysisID = analysis_id,
    ReleaseID = release_id,
    Model = "Feature ~ Adiposity + Role; FamilyID duplicateCorrelation block",
    Results = result_list,
    Summary = summary_table,
    InputManifest = input_manifest,
    CacheManifest = cache_manifest,
    Parameters = parameters
  ),
  statistics_file,
  compress = "gzip"
)

writexl::write_xlsx(
  list(
    Omics_Summary = summary_table,
    Volcano_Labels = label_table,
    Top100_Per_Platform = top_results,
    Input_Manifest = input_manifest,
    Cache_Manifest = cache_manifest,
    Parameters = parameters
  ),
  file.path(fig02a_dir, "Fig02A_Full_SourceData.xlsx")
)

notes <- c(
  "# Figure 2 baseline adiposity analysis",
  "",
  paste0("Analysis ID: `", analysis_id, "`"),
  "",
  "## Analysis definition",
  "",
  "- Baseline is strictly `Pre`.",
  "- Obese is body fat percentage >30%; Lean is <=30%.",
  "- The model is `Feature ~ Adiposity + Role`.",
  "- `FamilyID` is supplied to limma `duplicateCorrelation`.",
  "- Age, sex and sample-quality weights are not included.",
  "- Robust, trend-aware empirical Bayes moderation is used.",
  "",
  "## Figure definitions",
  "",
  "- Figure 2A preserves the submitted mirrored layout but uses v5 data and exact `Pre` samples.",
  "- Figure 2A bars show the percentage of tested features with nominal P <0.05 in each effect direction.",
  "- Figure 2C presents serum proteomics and metabolomics.",
  "- Figure 2D presents the three adipose omics platforms.",
  "- Figure 2F presents the three muscle omics platforms.",
  "- Volcano-plot y axes show -log10(nominal P).",
  "- Color denotes nominal P <0.05 and effect direction. BH-FDR is retained in the source tables.",
  "- Methylation labels use UCSC RefGene names from the official Illumina EPIC ilm10b4 hg19 annotation package; CpG identifiers remain in source data.",
  paste0("- Volcano values above -log10(P) = ", volcano_y_cap, " are displayed at the axis ceiling; exact values remain in the source data."),
  "",
  "## Interpretation boundary",
  "",
  "- Nominal P-based panels describe the signal landscape and do not define multiplicity-controlled discoveries.",
  "- Confirmatory feature claims must use the platform-specific BH-FDR values.",
  "- Effect estimates are on platform-specific processed scales and must not be compared numerically across platforms."
)
writeLines(notes, file.path(fig02a_dir, "Fig02_AB_NOTES.md"), useBytes = TRUE)

output_files <- c(
  file.path(fig02a_dir, "Fig02A_Baseline_EightOmics_Landscape.pdf"),
  file.path(fig02a_dir, "Fig02A_Baseline_EightOmics_Landscape.png"),
  file.path(fig02a_dir, "Fig02A_SourceData.xlsx"),
  file.path(fig02c_dir, "Fig02C_Serum_Omics_Volcanoes.pdf"),
  file.path(fig02c_dir, "Fig02C_Serum_Omics_Volcanoes.png"),
  file.path(fig02d_dir, "Fig02D_Adipose_Omics_Volcanoes.pdf"),
  file.path(fig02d_dir, "Fig02D_Adipose_Omics_Volcanoes.png"),
  file.path(fig02f_dir, "Fig02F_Muscle_Omics_Volcanoes.pdf"),
  file.path(fig02f_dir, "Fig02F_Muscle_Omics_Volcanoes.png"),
  methylation_annotation_file,
  file.path(fig02a_dir, "Fig02A_Full_Statistics.rds"),
  file.path(fig02a_dir, "Fig02A_Full_SourceData.xlsx"),
  file.path(fig02a_dir, "Fig02_AB_NOTES.md")
)

run_log <- c(
  paste0("AnalysisID=", analysis_id),
  paste0("ReleaseID=", release_id),
  paste0("RunTime=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("RVersion=", R.version.string),
  "Model=Feature ~ Adiposity + Role",
  "FamilyDependence=duplicateCorrelation(block=FamilyID)",
  "Baseline=Pre",
  "AgeIncluded=FALSE",
  "SexIncluded=FALSE",
  "SampleWeights=FALSE",
  "FigureThreshold=Nominal P <0.05",
  "MultipleTesting=BH within platform retained in source data",
  "",
  "Inputs:",
  paste(input_manifest$Dataset, input_manifest$InputFile, input_manifest$InputSHA256, sep = "\t"),
  "",
  "Validated caches:",
  paste(cache_manifest$Dataset, cache_manifest$CacheReused, cache_manifest$CacheFile, cache_manifest$CacheSHA256, sep = "\t"),
  "",
  "Outputs:",
  vapply(output_files[file.exists(output_files)], function(path) {
    paste0(
      basename(path),
      "\t",
      digest::digest(file = path, algo = "sha256", serialize = FALSE),
      "\t",
      file.info(path)$size
    )
  }, character(1)),
  "",
  capture.output(sessionInfo())
)
writeLines(
  run_log,
  file.path(fig02a_dir, "Fig02_AB_RunLog.txt"),
  useBytes = TRUE
)

message("Completed Figure 2A/C/D/F baseline adiposity analysis.")
