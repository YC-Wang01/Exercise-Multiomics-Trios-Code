# Figure S1: platform-specific PCA and adiposity-stratified cohort distributions.

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

p_load(readr, dplyr, tidyr, tibble, ggplot2, cowplot, openxlsx, matrixStats, digest)
options(cli.unicode = FALSE, pillar.unicode = FALSE, width = 220)
set.seed(20260802)

OUT_DIR <- path_project("fig01", "figS1")
FIG1_DIR <- path_project("fig01")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG1_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_COMBINED <- file.path(OUT_DIR, "FigS1_Cell_Layout_BH_FDR.pdf")
OUT_COMBINED_PNG <- file.path(OUT_DIR, "FigS1_Cell_Layout_BH_FDR.png")
OUT_SERUM <- file.path(OUT_DIR, "FigS01A_Serum_PCA.pdf")
OUT_ADIPOSE <- file.path(OUT_DIR, "FigS01B_Adipose_PCA.pdf")
OUT_MUSCLE <- file.path(OUT_DIR, "FigS01C_Muscle_PCA.pdf")
OUT_CLINICAL <- file.path(OUT_DIR, "FigS01D-F_Cell_Boxplots_BH_FDR.pdf")
OUT_FIG1_BMI <- file.path(FIG1_DIR, "Fig1_Original_BMI_BH_FDR.pdf")
OUT_FIG1_FITNESS <- file.path(FIG1_DIR, "Fig1_Original_Aerobic_Fitness_BH_FDR.pdf")
OUT_FIG1_VISCERAL <- file.path(FIG1_DIR, "Fig1_Original_Visceral_Fat_BH_FDR.pdf")
OUT_FIG1_LIVER <- file.path(FIG1_DIR, "Fig1_Original_Liver_Fat_BH_FDR.pdf")
OUT_XLSX <- file.path(OUT_DIR, "FigS1_and_Fig1_Original4_SourceData.xlsx")
OUT_RDS <- file.path(OUT_DIR, "FigS01_PCA_SelectedFeatures.rds")
OUT_LOG <- file.path(OUT_DIR, "FigS01_Run.log")

MIN_OBSERVED_FRACTION <- 0.70
ROLE_LEVELS <- c("Daughter", "Mother", "Father")
ADIPOSITY_LEVELS <- c("Obese", "Lean")
GROUP_LEVELS <- c("DO", "DL", "MO", "ML", "FO", "FL")

# Historical Figure S1 visual encoding retained.
ROLE_COLORS <- c(
  Daughter = "#1B9E77",
  Mother = "#E64B35",
  Father = "#4DBBD5"
)
PHASE_SHAPES <- c(Baseline = 16, Exercise = 17)
PHASE_ELLIPSE_COLORS <- c(Baseline = "#4E79A7", Exercise = "#F28E6B")
ADIPOSITY_COLORS <- c(Obese = "#F79647", Lean = "#3DA6AE")

DATASETS <- tibble::tribble(
  ~Dataset, ~Compartment, ~Display, ~InputFile, ~InputScale, ~MaxFeatures, ~Timepoints, ~Order,
  "Serum_Proteomics", "Serum", "Serum proteomics", "Analysis_Serum_Proteomics.csv", "Log2 sample-median-centred abundance", 227L, "Pre|Post1h|Post3h", 1L,
  "Serum_Metabonomics", "Serum", "Serum metabonomics", "IntegrationZ_Serum_Metabonomics.csv", "Feature-wise z score for mixed-unit PCA", 137L, "Pre|Post1h|Post3h", 2L,
  "Adipose_Microarray", "Adipose", "Adipose microarray", "Analysis_Adipose_Microarray.csv", "RMA log2 expression", 5000L, "Pre", 3L,
  "Adipose_Methylation", "Adipose", "Adipose methylation", "Analysis_Adipose_Methylation_MValue.rds", "ssNoob M value", 5000L, "Pre|Post3h", 4L,
  "Adipose_Proteomics", "Adipose", "Adipose proteomics", "Analysis_Adipose_Proteomics.csv", "Log2 sample-median-centred abundance", 3000L, "Pre|Post3h", 5L,
  "Muscle_Microarray", "Muscle", "Muscle microarray", "Analysis_Muscle_Microarray.csv", "RMA log2 expression", 5000L, "Pre", 6L,
  "Muscle_Methylation", "Muscle", "Muscle methylation", "Analysis_Muscle_Methylation_MValue.rds", "ssNoob M value", 5000L, "Pre|Post3h", 7L,
  "Muscle_Proteomics", "Muscle", "Muscle proteomics", "Analysis_Muscle_Proteomics.csv", "Log2 sample-median-centred abundance", 3000L, "Pre|Post3h", 8L
)

sample_sheet <- readr::read_csv(
  path_analysis_ready("Metadata", "Project_Sample_Sheet.csv"),
  show_col_types = FALSE,
  progress = FALSE
) %>%
  dplyr::mutate(
    SampleID = tolower(as.character(SampleID)),
    Dataset = as.character(Dataset),
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = factor(as.character(Role), levels = ROLE_LEVELS),
    Timepoint = as.character(Timepoint)
  )

hash_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

safe_row_vars <- function(matrix) {
  suppressWarnings(matrixStats::rowVars(matrix, na.rm = TRUE))
}

fit_dataset_pca <- function(spec) {
  dataset <- spec$Dataset[[1]]
  input_path <- path_analysis_ready(spec$InputFile[[1]])
  allowed_timepoints <- strsplit(spec$Timepoints[[1]], "\\|", fixed = FALSE)[[1]]
  message("PCA: ", dataset, " [", paste(allowed_timepoints, collapse = ", "), "]")
  matrix <- read_analysis_matrix(input_path)

  metadata <- sample_sheet %>%
    dplyr::filter(
      Dataset == dataset,
      Timepoint %in% allowed_timepoints,
      SampleID %in% colnames(matrix)
    ) %>%
    dplyr::distinct(SampleID, .keep_all = TRUE)
  common_samples <- colnames(matrix)[colnames(matrix) %in% metadata$SampleID]
  if (length(common_samples) < 3L) stop(dataset, ": fewer than three mapped display samples.")
  matrix <- matrix[, common_samples, drop = FALSE]
  metadata <- metadata[match(common_samples, metadata$SampleID), , drop = FALSE]
  if (any(is.na(metadata$SampleID))) stop(dataset, ": sample metadata alignment failed.")

  observed_fraction <- rowMeans(is.finite(matrix))
  row_variance <- safe_row_vars(matrix)
  eligible <- is.finite(row_variance) & row_variance > 0 &
    observed_fraction >= MIN_OBSERVED_FRACTION
  eligible_index <- which(eligible)
  eligible_index <- eligible_index[order(row_variance[eligible_index], decreasing = TRUE)]
  selected_index <- head(eligible_index, spec$MaxFeatures[[1]])
  if (length(selected_index) < 2L) stop(dataset, ": insufficient eligible features for PCA.")

  selected <- matrix[selected_index, , drop = FALSE]
  row_medians <- matrixStats::rowMedians(selected, na.rm = TRUE)
  missing_index <- which(!is.finite(selected), arr.ind = TRUE)
  if (nrow(missing_index) > 0L) selected[missing_index] <- row_medians[missing_index[, 1]]
  if (any(!is.finite(selected))) stop(dataset, ": non-finite value remained after PCA-only filling.")

  pca <- stats::prcomp(t(selected), center = TRUE, scale. = FALSE, rank. = 5)
  variance_fraction <- pca$sdev^2 / sum(pca$sdev^2)
  score_count <- min(5L, ncol(pca$x))
  score_frame <- as.data.frame(pca$x[, seq_len(score_count), drop = FALSE])
  score_frame$SampleID <- rownames(score_frame)
  scores <- metadata %>%
    dplyr::select(
      Dataset, SampleID, FamilyID, Role, Timepoint, Clinical_Subject_ID,
      QCReviewFlag, PairContainsQCReviewSample, RecommendedSensitivityExclusion
    ) %>%
    dplyr::left_join(score_frame, by = "SampleID") %>%
    dplyr::mutate(
      Display = spec$Display[[1]],
      Compartment = spec$Compartment[[1]],
      Phase = factor(ifelse(Timepoint == "Pre", "Baseline", "Exercise"), levels = names(PHASE_SHAPES)),
      TimepointOrder = match(Timepoint, c("Pre", "Post1h", "Post3h")),
      PC1VariancePercent = 100 * variance_fraction[[1]],
      PC2VariancePercent = 100 * variance_fraction[[2]],
      .before = 1
    ) %>%
    dplyr::arrange(Clinical_Subject_ID, TimepointOrder)

  variance <- tibble::tibble(
    Dataset = dataset,
    Display = spec$Display[[1]],
    PC = paste0("PC", seq_along(variance_fraction)),
    VariancePercent = 100 * variance_fraction
  )
  features <- tibble::tibble(
    Dataset = dataset,
    FeatureID = rownames(matrix)[selected_index],
    Variance = row_variance[selected_index],
    ObservedFraction = observed_fraction[selected_index],
    VarianceRank = seq_along(selected_index)
  )
  audit <- tibble::tibble(
    Dataset = dataset,
    Compartment = spec$Compartment[[1]],
    Display = spec$Display[[1]],
    InputFile = spec$InputFile[[1]],
    InputScale = spec$InputScale[[1]],
    DisplayTimepoints = paste(allowed_timepoints, collapse = ", "),
    InputSHA256 = hash_file(input_path),
    SamplesUsed = ncol(matrix),
    FeaturesInput = nrow(matrix),
    FeaturesEligible = length(eligible_index),
    FeaturesUsed = length(selected_index),
    MissingValuesFilledForPCA = nrow(missing_index),
    MinimumObservedFraction = MIN_OBSERVED_FRACTION,
    PCAFeatureCentering = TRUE,
    PCAFeatureScaling = FALSE,
    MissingValuePolicy = "Feature median filling for PCA projection only; source matrices unchanged"
  )

  rm(matrix, selected, pca)
  invisible(gc())
  list(scores = scores, variance = variance, features = features, audit = audit)
}

pca_runs <- lapply(seq_len(nrow(DATASETS)), function(index) {
  fit_dataset_pca(DATASETS[index, , drop = FALSE])
})
pca_scores <- dplyr::bind_rows(lapply(pca_runs, `[[`, "scores"))
pca_variance <- dplyr::bind_rows(lapply(pca_runs, `[[`, "variance"))
pca_features <- dplyr::bind_rows(lapply(pca_runs, `[[`, "features"))
pca_audit <- dplyr::bind_rows(lapply(pca_runs, `[[`, "audit"))
rm(pca_runs)
invisible(gc())

pca_theme <- ggplot2::theme_bw(base_family = "sans", base_size = 7.2) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "plain", size = 8.2, hjust = 0.5),
    axis.title = ggplot2::element_text(size = 7.0),
    axis.text = ggplot2::element_text(size = 7.0, color = "#222222"),
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = "#222222", linewidth = 0.35),
    axis.ticks = ggplot2::element_line(color = "#222222", linewidth = 0.30),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 7),
    legend.key.height = grid::unit(3.4, "mm"),
    legend.key.width = grid::unit(4.0, "mm"),
    plot.margin = ggplot2::margin(2, 3, 2, 3)
  )

plot_pca <- function(dataset, show_legend = FALSE) {
  data <- pca_scores %>% dplyr::filter(Dataset == dataset)
  display <- unique(data$Display)
  pc1 <- unique(data$PC1VariancePercent)
  pc2 <- unique(data$PC2VariancePercent)

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_path(
      data = data %>% dplyr::filter(dplyr::n_distinct(Timepoint) > 1),
      ggplot2::aes(group = Clinical_Subject_ID, color = Role),
      linewidth = 0.30,
      alpha = 0.24,
      arrow = grid::arrow(length = grid::unit(1.1, "mm"), type = "closed"),
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = Role, shape = Phase),
      size = 1.55,
      alpha = 0.90,
      stroke = 0.35
    ) +
    ggplot2::scale_color_manual(values = ROLE_COLORS, drop = FALSE) +
    ggplot2::scale_shape_manual(values = PHASE_SHAPES, drop = FALSE) +
    ggplot2::labs(
      title = display,
      x = sprintf("PC1 (%.1f%%)", pc1[[1]]),
      y = sprintf("PC2 (%.1f%%)", pc2[[1]]),
      color = NULL,
      shape = NULL
    ) +
    pca_theme +
    ggplot2::theme(legend.position = if (show_legend) "right" else "none")

  for (phase_name in names(PHASE_ELLIPSE_COLORS)) {
    phase_data <- data %>% dplyr::filter(Phase == phase_name)
    if (nrow(phase_data) >= 4L) {
      plot <- plot + ggplot2::stat_ellipse(
        data = phase_data,
        ggplot2::aes(x = PC1, y = PC2),
        inherit.aes = FALSE,
        type = "norm",
        level = 0.90,
        linewidth = 0.42,
        linetype = "dashed",
        color = PHASE_ELLIPSE_COLORS[[phase_name]],
        na.rm = TRUE
      )
    }
  }
  plot
}

pca_legend <- cowplot::get_legend(
  plot_pca("Serum_Proteomics", show_legend = TRUE) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(order = 1, override.aes = list(color = "#222222")),
      color = ggplot2::guide_legend(
        order = 2,
        override.aes = list(shape = 16, linewidth = 0, linetype = 0, alpha = 1)
      )
    ) +
    ggplot2::theme(legend.position = "right")
)

make_pca_section <- function(include_headings = TRUE) {
  row1 <- cowplot::plot_grid(
    plot_pca("Serum_Proteomics"),
    plot_pca("Adipose_Microarray"),
    plot_pca("Adipose_Methylation"),
    plot_pca("Adipose_Proteomics"),
    pca_legend,
    nrow = 1,
    rel_widths = c(1, 1, 1, 1, 0.53)
  )
  row2 <- cowplot::plot_grid(
    plot_pca("Serum_Metabonomics"),
    plot_pca("Muscle_Microarray"),
    plot_pca("Muscle_Methylation"),
    plot_pca("Muscle_Proteomics"),
    pca_legend,
    nrow = 1,
    rel_widths = c(1, 1, 1, 1, 0.53)
  )
  if (!include_headings) {
    return(cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1)))
  }
  serum_title <- cowplot::ggdraw() +
    cowplot::draw_label("A  Sample distribution in serum", x = 0.02, hjust = 0, fontface = "bold", size = 10)
  adipose_title <- cowplot::ggdraw() +
    cowplot::draw_label("B  Sample distribution in adipose", x = 0.5, hjust = 0.5, fontface = "bold", size = 10)
  muscle_title <- cowplot::ggdraw() +
    cowplot::draw_label("C  Sample distribution in muscle", x = 0.5, hjust = 0.5, fontface = "bold", size = 10)
  title_row1 <- cowplot::plot_grid(
    serum_title, adipose_title, cowplot::ggdraw(),
    nrow = 1,
    rel_widths = c(1, 3, 0.53)
  )
  title_row2 <- cowplot::plot_grid(
    cowplot::ggdraw(), muscle_title, cowplot::ggdraw(),
    nrow = 1,
    rel_widths = c(1, 3, 0.53)
  )
  cowplot::plot_grid(
    title_row1, row1, title_row2, row2,
    ncol = 1,
    rel_heights = c(0.13, 1, 0.13, 1)
  )
}

clinical_core_path <- path_analysis_ready("Metadata", "Clinical_Core.csv")
clinical_phenotypes_path <- path_analysis_ready("Metadata", "Clinical_Phenotypes.csv")
clinical_diet_path <- path_analysis_ready("Metadata", "Clinical_Diet.csv")

clinical_core <- readr::read_csv(clinical_core_path, show_col_types = FALSE, progress = FALSE)
clinical_phenotypes <- readr::read_csv(clinical_phenotypes_path, show_col_types = FALSE, progress = FALSE)
clinical_diet <- readr::read_csv(clinical_diet_path, show_col_types = FALSE, progress = FALSE)

phenotype_variables <- c(
  "SubcutaneousFatKg_G", "Total_FM_2010_G", "Total_LM_2010_G", "Total_BM_2010_G",
  "DLW_TEE_G", "BMR_REE30min_G", "DLW_AEE_G", "WAISTLINEcm_84M_G",
  "BMI_G", "MOU_VO2maxkg_G", "IntraperitonealFatKg_G", "Liverfat_G"
)
diet_variables <- c("Energy_kcal", "Ch_E_pros", "Prot_E_pros", "Fat_E_pros")
required_phenotype_variables <- c("Subject_UID", phenotype_variables)
required_diet_variables <- c("Subject_UID", diet_variables)
if (!all(required_phenotype_variables %in% names(clinical_phenotypes))) {
  stop("Clinical_Phenotypes.csv is missing Figure S1 variables: ",
       paste(setdiff(required_phenotype_variables, names(clinical_phenotypes)), collapse = ", "))
}
if (!all(required_diet_variables %in% names(clinical_diet))) {
  stop("Clinical_Diet.csv is missing Figure S1 variables: ",
       paste(setdiff(required_diet_variables, names(clinical_diet)), collapse = ", "))
}

clinical_wide <- clinical_core %>%
  dplyr::select(Family_UID, Subject_UID, Role, Adiposity_Group, BodyFat_percent) %>%
  dplyr::left_join(
    clinical_phenotypes %>% dplyr::select(dplyr::all_of(required_phenotype_variables)),
    by = "Subject_UID"
  ) %>%
  dplyr::left_join(
    clinical_diet %>% dplyr::select(dplyr::all_of(required_diet_variables)),
    by = "Subject_UID"
  ) %>%
  dplyr::mutate(
    Role = factor(Role, levels = ROLE_LEVELS),
    Adiposity_Group = factor(Adiposity_Group, levels = ADIPOSITY_LEVELS),
    Group = factor(
      paste0(substr(as.character(Role), 1, 1), substr(as.character(Adiposity_Group), 1, 1)),
      levels = GROUP_LEVELS
    ),
    DietAllZero = dplyr::if_all(dplyr::all_of(diet_variables), ~ is.na(.x) | .x == 0),
    TotalFatMass_kg = Total_FM_2010_G / 1000,
    TotalLeanMass_kg = Total_LM_2010_G / 1000,
    TotalBoneMass_kg = Total_BM_2010_G / 1000,
    EnergyIntake_kcal_day = dplyr::if_else(DietAllZero, NA_real_, Energy_kcal),
    Carbohydrate_E_percent = dplyr::if_else(DietAllZero, NA_real_, Ch_E_pros),
    Protein_E_percent = dplyr::if_else(DietAllZero, NA_real_, Prot_E_pros),
    Fat_E_percent = dplyr::if_else(DietAllZero, NA_real_, Fat_E_pros)
  )

CLINICAL_SPECS <- tibble::tribble(
  ~Variable, ~SourceVariable, ~Panel, ~RowOrder, ~ColumnOrder, ~Title, ~YAxis, ~DisplayTransform,
  "SubcutaneousFatKg_G", "SubcutaneousFatKg_G", "D", 1L, 1L, "Subcutaneous fat (kg)", "Subcutaneous fat (kg)", "None",
  "TotalFatMass_kg", "Total_FM_2010_G", "D", 1L, 2L, "Fat mass (kg)", "Fat mass (kg)", "g to kg (/1000)",
  "TotalLeanMass_kg", "Total_LM_2010_G", "D", 1L, 3L, "Lean mass (kg)", "Lean mass (kg)", "g to kg (/1000)",
  "TotalBoneMass_kg", "Total_BM_2010_G", "D", 1L, 4L, "Bone mass (kg)", "Bone mass (kg)", "g to kg (/1000)",
  "DLW_TEE_G", "DLW_TEE_G", "E", 2L, 1L, "TEE (kcal/day)", "TEE (kcal/day)", "None",
  "BMR_REE30min_G", "BMR_REE30min_G", "E", 2L, 2L, "REE (kcal/day)", "REE (kcal/day)", "None",
  "DLW_AEE_G", "DLW_AEE_G", "E", 2L, 3L, "AEE (kcal/day)", "AEE (kcal/day)", "None",
  "WAISTLINEcm_84M_G", "WAISTLINEcm_84M_G", "E", 2L, 4L, "Waist circumference (cm)", "Waist circumference (cm)", "None",
  "EnergyIntake_kcal_day", "Energy_kcal", "F", 3L, 1L, "Energy intake (kcal/day)", "Energy intake (kcal/day)", "All-zero diet record to missing",
  "Carbohydrate_E_percent", "Ch_E_pros", "F", 3L, 2L, "Carbohydrate (E%)", "Carbohydrate (E%)", "All-zero diet record to missing",
  "Protein_E_percent", "Prot_E_pros", "F", 3L, 3L, "Protein (E%)", "Protein (E%)", "All-zero diet record to missing",
  "Fat_E_percent", "Fat_E_pros", "F", 3L, 4L, "Fat (E%)", "Fat (E%)", "All-zero diet record to missing",
  "BMI_G", "BMI_G", "G", 4L, 1L, "BMI", "BMI (kg/m2)", "None",
  "MOU_VO2maxkg_G", "MOU_VO2maxkg_G", "G", 4L, 2L, "Aerobic fitness", "VO2max (ml/kg/min)", "None",
  "IntraperitonealFatKg_G", "IntraperitonealFatKg_G", "G", 4L, 3L, "Visceral fat", "Visceral fat (kg)", "None",
  "Liverfat_G", "Liverfat_G", "G", 4L, 4L, "Liver fat", "Liver fat (%)", "None"
)

clinical_long <- clinical_wide %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(CLINICAL_SPECS$Variable),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  dplyr::left_join(CLINICAL_SPECS, by = "Variable") %>%
  dplyr::mutate(
    Variable = factor(Variable, levels = CLINICAL_SPECS$Variable),
    Group = factor(Group, levels = GROUP_LEVELS)
  )

clinical_counts <- clinical_long %>%
  dplyr::group_by(Panel, RowOrder, ColumnOrder, Variable, Title, Role, Adiposity_Group, Group) %>%
  dplyr::summarise(
    NObserved = sum(is.finite(Value)),
    NMissing = sum(!is.finite(Value)),
    .groups = "drop"
  )

run_role_adiposity_test <- function(variable, role) {
  spec <- CLINICAL_SPECS %>% dplyr::filter(Variable == variable)
  data <- clinical_long %>%
    dplyr::filter(Variable == variable, Role == role, is.finite(Value))
  obese <- data$Value[data$Adiposity_Group == "Obese"]
  lean <- data$Value[data$Adiposity_Group == "Lean"]
  p_value <- if (length(obese) >= 3L && length(lean) >= 3L) {
    stats::wilcox.test(
      obese,
      lean,
      alternative = "two.sided",
      exact = FALSE,
      correct = FALSE
    )$p.value
  } else {
    NA_real_
  }
  tibble::tibble(
    Panel = spec$Panel[[1]],
    RowOrder = spec$RowOrder[[1]],
    ColumnOrder = spec$ColumnOrder[[1]],
    Variable = variable,
    Title = spec$Title[[1]],
    Role = role,
    NObese = length(obese),
    NLean = length(lean),
    MedianObese = if (length(obese) > 0L) stats::median(obese) else NA_real_,
    MedianLean = if (length(lean) > 0L) stats::median(lean) else NA_real_,
    MedianDifference_ObeseMinusLean = if (length(obese) > 0L && length(lean) > 0L) {
      stats::median(obese) - stats::median(lean)
    } else {
      NA_real_
    },
    PValue = p_value
  )
}

clinical_group_tests <- dplyr::bind_rows(lapply(CLINICAL_SPECS$Variable, function(variable) {
  dplyr::bind_rows(lapply(ROLE_LEVELS, function(role) run_role_adiposity_test(variable, role)))
})) %>%
  dplyr::mutate(
    BH_FDR = stats::p.adjust(PValue, method = "BH"),
    P_lt_0_05 = !is.na(PValue) & PValue < 0.05,
    BH_FDR_lt_0_05 = !is.na(BH_FDR) & BH_FDR < 0.05,
    BH_FDR_lt_0_10 = !is.na(BH_FDR) & BH_FDR < 0.10,
    BH_Significance = dplyr::case_when(
      is.na(BH_FDR) ~ NA_character_,
      BH_FDR < 0.001 ~ "***",
      BH_FDR < 0.01 ~ "**",
      BH_FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::arrange(BH_FDR, PValue, RowOrder, ColumnOrder, Role)

clinical_test_summary <- tibble::tibble(
  TestFamily = "Sixteen endpoints x three Roles; Obese versus Lean within Role",
  Test = "Two-sided Wilcoxon rank-sum test without continuity correction",
  MultipleTesting = "Benjamini-Hochberg correction across all available tests",
  TestsAttempted = nrow(clinical_group_tests),
  TestsWithPValue = sum(is.finite(clinical_group_tests$PValue)),
  RawP_lt_0_05 = sum(clinical_group_tests$P_lt_0_05, na.rm = TRUE),
  BH_FDR_lt_0_05 = sum(clinical_group_tests$BH_FDR_lt_0_05, na.rm = TRUE),
  BH_FDR_lt_0_10 = sum(clinical_group_tests$BH_FDR_lt_0_10, na.rm = TRUE)
)

box_theme <- ggplot2::theme_bw(base_family = "sans", base_size = 7.3) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 7.5, hjust = 0.5),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(face = "bold", size = 7.0),
    axis.text.x = ggplot2::element_text(face = "bold", size = 7.0, color = "#222222"),
    axis.text.y = ggplot2::element_text(size = 7.0, color = "#222222"),
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = "#222222", linewidth = 0.35),
    axis.ticks = ggplot2::element_line(color = "#222222", linewidth = 0.30),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 7),
    legend.key.height = grid::unit(4.0, "mm"),
    legend.key.width = grid::unit(4.0, "mm"),
    plot.margin = ggplot2::margin(2, 3, 2, 3)
  )

natural_breaks_5 <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) return(list(breaks = 0:4, accuracy = 1))
  value_min <- min(values)
  value_max <- max(values)
  value_range <- value_max - value_min
  if (!is.finite(value_range) || value_range <= 0) value_range <- max(abs(values), 1)
  base_power <- floor(log10(value_range / 4))
  multipliers <- tibble::tibble(
    Multiplier = c(1, 2, 2.5, 5, 10, 1.5, 4, 6, 8),
    SimplicityPenalty = c(rep(0, 5), rep(0.08, 4))
  )
  candidates <- dplyr::bind_rows(lapply(seq(base_power - 2, base_power + 2), function(power) {
    multipliers %>% dplyr::mutate(Step = Multiplier * 10^power)
  })) %>%
    dplyr::filter(is.finite(Step), Step > 0) %>%
    dplyr::distinct(Step, .keep_all = TRUE)
  options <- dplyr::bind_rows(lapply(seq_len(nrow(candidates)), function(index) {
    step <- candidates$Step[[index]]
    lower <- floor((value_min + 1e-12 * value_range) / step) * step
    dplyr::bind_rows(lapply(5L, function(count) {
      upper <- lower + (count - 1L) * step
      if (upper < value_max - 1e-10 * value_range) return(NULL)
      tibble::tibble(
        Step = step,
        Count = count,
        Lower = lower,
        Upper = upper,
        Score = (upper - lower) / value_range - 1 +
          candidates$SimplicityPenalty[[index]]
      )
    }))
  }))
  if (nrow(options) == 0L) stop("Unable to construct a five-tick natural y-axis.")
  selected <- options %>% dplyr::arrange(Score, Count, Step) %>% dplyr::slice(1)
  step <- selected$Step[[1]]
  breaks <- selected$Lower[[1]] + seq_len(selected$Count[[1]]) * step - step
  accuracy <- if (step >= 1 && abs(step - round(step)) < 1e-9) {
    1
  } else if (step >= 0.1) {
    0.1
  } else if (step >= 0.01) {
    0.01
  } else {
    0.001
  }
  list(breaks = breaks, accuracy = accuracy)
}

plot_clinical <- function(variable, show_legend = FALSE) {
  spec <- CLINICAL_SPECS %>% dplyr::filter(Variable == variable)
  data <- clinical_long %>% dplyr::filter(Variable == variable, is.finite(Value))
  annotation <- clinical_group_tests %>%
    dplyr::filter(Variable == variable, BH_FDR_lt_0_05) %>%
    dplyr::mutate(
      XStart = dplyr::recode(Role, Daughter = 1, Mother = 3, Father = 5),
      XEnd = XStart + 1
    )
  value_range <- diff(range(data$Value, na.rm = TRUE))
  if (!is.finite(value_range) || value_range <= 0) value_range <- max(abs(data$Value), na.rm = TRUE)
  if (!is.finite(value_range) || value_range <= 0) value_range <- 1
  if (nrow(annotation) > 0L) {
    annotation <- annotation %>%
      dplyr::mutate(
        Y = max(data$Value, na.rm = TRUE) + 0.075 * value_range,
        YTick = Y - 0.025 * value_range,
        YLabel = Y + 0.020 * value_range
      )
  }
  annotation_values <- if (nrow(annotation) > 0L) annotation$YLabel else numeric()
  axis_spec <- natural_breaks_5(c(data$Value, annotation_values))
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = Group, y = Value, fill = Adiposity_Group)) +
    ggplot2::stat_boxplot(
      geom = "errorbar",
      width = 0.28,
      color = "#333333",
      linewidth = 0.36
    ) +
    ggplot2::geom_boxplot(
      width = 0.62,
      alpha = 0.55,
      outlier.shape = NA,
      color = "#333333",
      linewidth = 0.38
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = Adiposity_Group),
      size = 1.15,
      alpha = 0.68,
      stroke = 0,
      position = ggplot2::position_jitter(width = 0.13, height = 0, seed = 20260802)
    ) +
    ggplot2::scale_fill_manual(values = ADIPOSITY_COLORS, drop = FALSE) +
    ggplot2::scale_color_manual(values = ADIPOSITY_COLORS, drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_continuous(
      breaks = axis_spec$breaks,
      labels = scales::label_number(accuracy = axis_spec$accuracy, big.mark = ""),
      limits = range(axis_spec$breaks),
      expand = ggplot2::expansion(mult = c(0.05, 0.08))
    ) +
    ggplot2::labs(
      title = spec$Title[[1]],
      y = spec$YAxis[[1]],
      fill = NULL,
      color = NULL
    ) +
    box_theme +
    ggplot2::theme(legend.position = if (show_legend) "right" else "none")
  if (nrow(annotation) > 0L) {
    plot <- plot +
      ggplot2::geom_segment(
        data = annotation,
        ggplot2::aes(x = XStart, xend = XEnd, y = Y, yend = Y),
        inherit.aes = FALSE,
        color = "#37474F",
        linewidth = 0.42
      ) +
      ggplot2::geom_segment(
        data = annotation,
        ggplot2::aes(x = XStart, xend = XStart, y = YTick, yend = Y),
        inherit.aes = FALSE,
        color = "#37474F",
        linewidth = 0.42
      ) +
      ggplot2::geom_segment(
        data = annotation,
        ggplot2::aes(x = XEnd, xend = XEnd, y = YTick, yend = Y),
        inherit.aes = FALSE,
        color = "#37474F",
        linewidth = 0.42
      ) +
      ggplot2::geom_text(
        data = annotation,
        ggplot2::aes(x = (XStart + XEnd) / 2, y = YLabel, label = BH_Significance),
        inherit.aes = FALSE,
        color = "#263238",
        size = 2.45,
        fontface = "bold",
        vjust = 0
      )
  }
  plot
}

clinical_legend <- cowplot::get_legend(
  plot_clinical("SubcutaneousFatKg_G", show_legend = TRUE) +
    ggplot2::guides(color = "none", fill = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
    ggplot2::theme(legend.position = "right")
)

ROW_TITLES <- c(
  D = "Body composition by adiposity group",
  E = "Energy expenditure and waist circumference by adiposity group",
  F = "Dietary intake by adiposity group",
  G = "Adiposity, aerobic fitness and ectopic fat"
)

make_clinical_row <- function(panel, include_heading = TRUE) {
  specs <- CLINICAL_SPECS %>% dplyr::filter(Panel == panel) %>% dplyr::arrange(ColumnOrder)
  plots <- lapply(specs$Variable, plot_clinical)
  row <- cowplot::plot_grid(
    plotlist = c(plots, list(clinical_legend)),
    nrow = 1,
    rel_widths = c(rep(1, 4), 0.53)
  )
  if (!include_heading) return(row)
  title_strip <- cowplot::ggdraw() +
    cowplot::draw_label(panel, x = 0.008, hjust = 0, fontface = "bold", size = 11) +
    cowplot::draw_label(ROW_TITLES[[panel]], x = 0.445, hjust = 0.5, fontface = "bold", size = 9.5)
  cowplot::plot_grid(title_strip, row, ncol = 1, rel_heights = c(0.13, 1))
}

serum_export <- cowplot::plot_grid(
  plot_pca("Serum_Proteomics"),
  plot_pca("Serum_Metabonomics"),
  pca_legend,
  nrow = 1,
  rel_widths = c(1, 1, 0.55)
)
adipose_export <- cowplot::plot_grid(
  plot_pca("Adipose_Microarray"),
  plot_pca("Adipose_Methylation"),
  plot_pca("Adipose_Proteomics"),
  pca_legend,
  nrow = 1,
  rel_widths = c(1, 1, 1, 0.55)
)
muscle_export <- cowplot::plot_grid(
  plot_pca("Muscle_Microarray"),
  plot_pca("Muscle_Methylation"),
  plot_pca("Muscle_Proteomics"),
  pca_legend,
  nrow = 1,
  rel_widths = c(1, 1, 1, 0.55)
)

pca_section <- make_pca_section(include_headings = TRUE)
clinical_rows <- lapply(c("D", "E", "F"), make_clinical_row)
axis_break_audit <- dplyr::bind_rows(lapply(CLINICAL_SPECS$Variable, function(variable) {
  built <- ggplot2::ggplot_build(plot_clinical(variable))
  breaks <- built$layout$panel_params[[1]]$y$breaks
  breaks <- breaks[is.finite(breaks)]
  spec <- CLINICAL_SPECS %>% dplyr::filter(Variable == variable)
  tibble::tibble(
    Panel = spec$Panel[[1]],
    Variable = variable,
    Title = spec$Title[[1]],
    BreakCount = length(breaks),
    Breaks = paste(format(breaks, trim = TRUE, scientific = FALSE), collapse = ", ")
  )
}))
if (any(axis_break_audit$BreakCount != 5L)) {
  stop("Figure S1 clinical y-axis audit found a panel without exactly five major tick labels.")
}
clinical_export <- cowplot::plot_grid(plotlist = clinical_rows, ncol = 1, rel_heights = rep(1, 3))
make_original_fig1_plot <- function(variable) {
  plot_clinical(variable, show_legend = TRUE) +
    ggplot2::guides(
      color = "none",
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    ggplot2::theme(legend.position = "right")
}
original_fig1_plots <- list(
  BMI_G = make_original_fig1_plot("BMI_G"),
  MOU_VO2maxkg_G = make_original_fig1_plot("MOU_VO2maxkg_G"),
  IntraperitonealFatKg_G = make_original_fig1_plot("IntraperitonealFatKg_G"),
  Liverfat_G = make_original_fig1_plot("Liverfat_G")
)
combined <- cowplot::plot_grid(
  pca_section,
  clinical_rows[[1]],
  clinical_rows[[2]],
  clinical_rows[[3]],
  ncol = 1,
  rel_heights = c(2.0, 1, 1, 1)
)

ggplot2::ggsave(OUT_SERUM, serum_export, width = 180, height = 69, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(OUT_ADIPOSE, adipose_export, width = 180, height = 69, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(OUT_MUSCLE, muscle_export, width = 180, height = 69, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(OUT_CLINICAL, clinical_export, width = 180, height = 126, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(OUT_FIG1_BMI, original_fig1_plots$BMI_G, width = 72, height = 56, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(OUT_FIG1_FITNESS, original_fig1_plots$MOU_VO2maxkg_G, width = 72, height = 56, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(OUT_FIG1_VISCERAL, original_fig1_plots$IntraperitonealFatKg_G, width = 72, height = 56, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(OUT_FIG1_LIVER, original_fig1_plots$Liverfat_G, width = 72, height = 56, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(OUT_COMBINED, combined, width = 180, height = 210, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(OUT_COMBINED_PNG, combined, width = 180, height = 210, units = "mm", dpi = 450, bg = "white")

saveRDS(split(pca_features, pca_features$Dataset), OUT_RDS, compress = "gzip")

palette_table <- dplyr::bind_rows(
  tibble::tibble(Encoding = "Role", Level = names(ROLE_COLORS), HexColor = unname(ROLE_COLORS), Shape = NA_integer_),
  tibble::tibble(Encoding = "Phase", Level = names(PHASE_SHAPES), HexColor = unname(PHASE_ELLIPSE_COLORS), Shape = unname(PHASE_SHAPES)),
  tibble::tibble(Encoding = "Adiposity", Level = names(ADIPOSITY_COLORS), HexColor = unname(ADIPOSITY_COLORS), Shape = NA_integer_)
)

input_manifest <- dplyr::bind_rows(
  pca_audit %>% dplyr::select(Dataset, InputFile, InputSHA256, InputScale),
  tibble::tibble(
    Dataset = c("Clinical_Core", "Clinical_Phenotypes", "Clinical_Diet"),
    InputFile = c("Metadata/Clinical_Core.csv", "Metadata/Clinical_Phenotypes.csv", "Metadata/Clinical_Diet.csv"),
    InputSHA256 = c(hash_file(clinical_core_path), hash_file(clinical_phenotypes_path), hash_file(clinical_diet_path)),
    InputScale = "Observed clinical or dietary values; no imputation"
  )
)

parameters <- tibble::tibble(
  Parameter = c(
    "Analysis-ready release", "PCA observed-feature threshold", "PCA feature selection",
    "PCA centering", "PCA scaling", "PCA missing-value handling", "PCA displayed baseline",
    "PCA displayed exercise samples", "Serum fasting sample handling", "Methylation PCA scale",
    "Serum metabolomics PCA scale", "Clinical missing-value handling", "Diet all-zero record handling",
    "DXA display conversion", "Adiposity definition", "Clinical inference", "Figure S1 scope"
  ),
  Value = c(
    "2026-08-01_v5", as.character(MIN_OBSERVED_FRACTION),
    "Highest-variance eligible features up to the dataset-specific maximum",
    "Feature centering enabled", "No additional feature scaling",
    "Feature median filling for PCA projection only; no matrix is overwritten", "Pre",
    "Post1h and Post3h for serum; Post3h for tissue platforms when available",
    "Fast is excluded from Figure S1 PCA and is not merged with Pre", "ssNoob M values",
    "Released IntegrationZ feature-wise z scores", "Observed values only; no imputation",
    "Rows with energy and all three macronutrient E% values equal to zero are displayed as missing",
    "Total fat, lean and bone mass source values divided by 1000 for kg display",
    "Lean <=30% body fat; Obese >30% body fat", "Role-stratified Wilcoxon tests; significance marks use BH FDR across 48 tests",
    "Figure S1 panels A-C: eight platform-specific PCA plots; panels D-F: the 12 Cell-format cohort boxplots. Four former main-Figure-1 endpoints are exported separately."
  )
)

workbook <- openxlsx::createWorkbook()
write_sheet <- function(name, data) {
  openxlsx::addWorksheet(workbook, name)
  openxlsx::writeData(workbook, name, data, withFilter = TRUE)
  openxlsx::freezePane(workbook, name, firstRow = TRUE)
  if (ncol(data) > 0L) openxlsx::setColWidths(workbook, name, cols = seq_len(ncol(data)), widths = "auto")
}
write_sheet("PCA_Scores", pca_scores)
write_sheet("PCA_Variance", pca_variance)
write_sheet("PCA_Audit", pca_audit)
write_sheet("PCA_Selected_Features", pca_features)
write_sheet("Clinical_Display_Values", clinical_long)
write_sheet("Clinical_Wide_Audit", clinical_wide)
write_sheet("Clinical_Counts", clinical_counts)
write_sheet("Clinical_Group_Tests", clinical_group_tests)
write_sheet("Clinical_Test_Summary", clinical_test_summary)
write_sheet("Clinical_Axis_Breaks", axis_break_audit)
write_sheet("Clinical_Specification", CLINICAL_SPECS)
write_sheet("Parameters", parameters)
write_sheet("Palette", palette_table)
write_sheet("InputManifest", input_manifest)
openxlsx::saveWorkbook(workbook, OUT_XLSX, overwrite = TRUE)

log_lines <- c(
  paste("Run timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("R:", R.version.string),
  paste("Combined figure:", OUT_COMBINED),
  paste("PCA score rows:", nrow(pca_scores)),
  paste("Selected feature rows:", nrow(pca_features)),
  paste("Clinical display rows:", nrow(clinical_long)),
  paste("Diet all-zero records treated as missing:", sum(clinical_wide$DietAllZero, na.rm = TRUE)),
  paste("Role-specific Obese-versus-Lean tests:", sum(is.finite(clinical_group_tests$PValue))),
  paste("Raw P < 0.05:", sum(clinical_group_tests$P_lt_0_05, na.rm = TRUE)),
  paste("BH FDR < 0.05:", sum(clinical_group_tests$BH_FDR_lt_0_05, na.rm = TRUE)),
  paste("BH FDR < 0.10:", sum(clinical_group_tests$BH_FDR_lt_0_10, na.rm = TRUE)),
  paste("Clinical y-axis major ticks per panel:", paste(sort(unique(axis_break_audit$BreakCount)), collapse = ", ")),
  "",
  capture.output(print(pca_audit)),
  "",
  capture.output(print(clinical_counts)),
  "",
  capture.output(print(clinical_test_summary)),
  "",
  capture.output(print(clinical_group_tests %>% dplyr::filter(P_lt_0_05 | BH_FDR_lt_0_10))),
  "",
  capture.output(print(axis_break_audit)),
  "",
  capture.output(utils::sessionInfo())
)
writeLines(log_lines, OUT_LOG, useBytes = TRUE)

message("Saved historical-format Figure S1 and traceable source data under ", OUT_DIR)
