# Figure 3A-B: acute serum multi-omics response and temporal trajectories.
#
# Pre is the formal exercise baseline. Subject-level changes are modelled with
# centered adiposity and role covariates; FamilyID is the limma correlation block.

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

p_load(readr, dplyr, tidyr, tibble, ggplot2, limma, patchwork, writexl, digest)
select <- dplyr::select
set.seed(20260805)

analysis_id <- "Fig03AB_Acute_Serum_Response_2026-08-05"
release_id <- "2026-08-01_v5"
fdr_threshold <- 0.05
minimum_pairs <- 10L
display_features <- 22L
time_levels <- c("Fast", "Pre", "Post1h", "Post3h")
contrast_levels <- c("Fast_vs_Pre", "Post1h_vs_Pre", "Post3h_vs_Pre")
acute_contrasts <- c("Post1h_vs_Pre", "Post3h_vs_Pre")
role_levels <- c("Daughter", "Mother", "Father")

lean_color <- "#3FA0A8"
obese_color <- "#F59646"
decrease_color <- "#2F6FA3"
neutral_color <- "#F7F7F4"
increase_color <- "#C94B45"

fig03_root <- path_project("fig03")
figure3_dir <- file.path(fig03_root, "Figure3")
fig03a_dir <- figure3_dir
fig03b_dir <- figure3_dir
dir.create(fig03a_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig03b_dir, recursive = TRUE, showWarnings = FALSE)

metadata_dir <- path_analysis_ready("Metadata")
sample_sheet_file <- file.path(metadata_dir, "Project_Sample_Sheet.csv")
clinical_file <- file.path(metadata_dir, "Clinical_Observed_Primary.csv")
manifest_file <- file.path(metadata_dir, "Matrix_Manifest.csv")

dataset_registry <- tibble::tribble(
  ~Dataset, ~InputFile, ~DictionaryFile, ~FeatureLabelColumn, ~PanelTitle,
  "Serum_Metabonomics", "Analysis_Serum_Metabonomics.csv",
  "Feature_Dictionary_Serum_Metabonomics.csv", "FeatureID", "Serum metabolomics",
  "Serum_Proteomics", "Analysis_Serum_Proteomics.csv",
  "Feature_Dictionary_Serum_Proteomics.csv", "GeneSymbol", "Serum proteomics"
)

manifest <- readr::read_csv(manifest_file, show_col_types = FALSE)
dataset_registry <- dataset_registry %>%
  left_join(
    manifest %>% select(Dataset, UseStatus, SHA256, Features, Samples, MissingRate),
    by = "Dataset"
  )
if (any(!analysis_status_is_released(dataset_registry$UseStatus))) {
  stop("One or more serum analysis matrices are not released.")
}

clinical <- readr::read_csv(clinical_file, show_col_types = FALSE) %>%
  transmute(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = factor(Role, levels = role_levels),
    FatPercent = suppressWarnings(as.numeric(Fat_p_G)),
    Adiposity = factor(classify_body_fat(FatPercent), levels = c("Lean", "Obese"))
  ) %>%
  filter(!is.na(FamilyID), !is.na(Clinical_Subject_ID), !is.na(Role), !is.na(Adiposity)) %>%
  distinct(FamilyID, Clinical_Subject_ID, Role, .keep_all = TRUE)

sample_sheet <- readr::read_csv(sample_sheet_file, show_col_types = FALSE) %>%
  mutate(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    SampleID = tolower(as.character(SampleID)),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort),
    RetainInPrimaryInput = as.logical(RetainInPrimaryInput),
    Timepoint = factor(Timepoint, levels = time_levels)
  ) %>%
  filter(
    AnalysisSet == "Primary",
    Is_Primary_Cohort,
    is.na(RetainInPrimaryInput) | RetainInPrimaryInput
  ) %>%
  inner_join(clinical, by = c("FamilyID", "Clinical_Subject_ID", "Role"))

read_dictionary <- function(dataset_row) {
  dictionary <- readr::read_csv(
    file.path(metadata_dir, dataset_row$DictionaryFile),
    show_col_types = FALSE
  )
  label_column <- dataset_row$FeatureLabelColumn
  if (!label_column %in% names(dictionary)) stop("Dictionary label column is missing.")
  labels <- as.character(dictionary[[label_column]])
  labels[is.na(labels) | !nzchar(trimws(labels))] <- as.character(dictionary$FeatureID)[
    is.na(labels) | !nzchar(trimws(labels))
  ]
  tibble(
    FeatureID = as.character(dictionary$FeatureID),
    FeatureLabel = make.unique(labels)
  )
}

prepare_dataset <- function(dataset_row) {
  input_path <- path_analysis_ready(dataset_row$InputFile)
  matrix <- read_analysis_matrix(input_path)
  metadata <- sample_sheet %>%
    filter(Dataset == dataset_row$Dataset, SampleID %in% colnames(matrix)) %>%
    distinct(Clinical_Subject_ID, Timepoint, .keep_all = TRUE) %>%
    arrange(Clinical_Subject_ID, Timepoint)
  index <- match(metadata$SampleID, colnames(matrix))
  if (anyNA(index)) stop("Sample alignment failed for ", dataset_row$Dataset)
  matrix <- matrix[, index, drop = FALSE]
  colnames(matrix) <- metadata$SampleID
  list(matrix = matrix, metadata = metadata, dictionary = read_dictionary(dataset_row))
}

fit_change <- function(dataset_name, prepared, post_timepoint) {
  metadata <- prepared$metadata
  pre <- metadata %>%
    filter(Timepoint == "Pre") %>%
    select(FamilyID, Clinical_Subject_ID, Role, Adiposity, PreSample = SampleID)
  post <- metadata %>%
    filter(Timepoint == post_timepoint) %>%
    select(FamilyID, Clinical_Subject_ID, Role, Adiposity, PostSample = SampleID)
  pairs <- inner_join(
    pre,
    post,
    by = c("FamilyID", "Clinical_Subject_ID", "Role", "Adiposity")
  ) %>%
    arrange(FamilyID, Role, Clinical_Subject_ID)
  if (nrow(pairs) < minimum_pairs) {
    stop("Too few paired samples for ", dataset_name, " ", post_timepoint)
  }
  pre_index <- match(pairs$PreSample, colnames(prepared$matrix))
  post_index <- match(pairs$PostSample, colnames(prepared$matrix))
  delta <- prepared$matrix[, post_index, drop = FALSE] - prepared$matrix[, pre_index, drop = FALSE]
  colnames(delta) <- pairs$Clinical_Subject_ID

  adiposity <- as.numeric(pairs$Adiposity == "Obese")
  role_mother <- as.numeric(pairs$Role == "Mother")
  role_father <- as.numeric(pairs$Role == "Father")
  design <- cbind(
    Overall = 1,
    Adiposity = adiposity - mean(adiposity),
    RoleMother = role_mother - mean(role_mother),
    RoleFather = role_father - mean(role_father)
  )
  keep <- rowSums(is.finite(delta)) >= minimum_pairs
  variable <- apply(delta, 1L, function(values) {
    values <- values[is.finite(values)]
    length(values) >= minimum_pairs && is.finite(stats::sd(values)) && stats::sd(values) > 0
  })
  delta <- delta[keep & variable, , drop = FALSE]
  correlation <- tryCatch(
    limma::duplicateCorrelation(delta, design, block = pairs$FamilyID)$consensus.correlation,
    error = function(error) NA_real_
  )
  if (!is.finite(correlation)) correlation <- 0
  fit <- limma::lmFit(delta, design, block = pairs$FamilyID, correlation = correlation)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

  extract_coef <- function(coefficient, effect_name) {
    result <- limma::topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
    tibble(
      Dataset = dataset_name,
      Contrast = paste0(post_timepoint, "_vs_Pre"),
      EffectType = effect_name,
      FeatureID = rownames(result),
      Effect = result$logFC,
      ModeratedT = result$t,
      PValue = result$P.Value,
      BH_FDR = result$adj.P.Val,
      NPairs = rowSums(is.finite(delta)),
      FamilyCorrelation = correlation
    )
  }
  list(
    results = bind_rows(
      extract_coef("Overall", "Overall response"),
      extract_coef("Adiposity", "Time by adiposity interaction")
    ),
    pairs = pairs
  )
}

prepared_list <- list()
result_list <- list()
pair_list <- list()
for (index in seq_len(nrow(dataset_registry))) {
  dataset_row <- dataset_registry[index, ]
  message("Preparing ", dataset_row$Dataset)
  prepared <- prepare_dataset(dataset_row)
  prepared_list[[dataset_row$Dataset]] <- prepared
  for (timepoint in c("Fast", "Post1h", "Post3h")) {
    fitted <- fit_change(dataset_row$Dataset, prepared, timepoint)
    result_list[[paste(dataset_row$Dataset, timepoint, sep = "__")]] <- fitted$results
    pair_list[[paste(dataset_row$Dataset, timepoint, sep = "__")]] <- fitted$pairs %>%
      mutate(Dataset = dataset_row$Dataset, Contrast = paste0(timepoint, "_vs_Pre"))
  }
}

all_results <- bind_rows(result_list)
all_pairs <- bind_rows(pair_list) %>%
  distinct(Dataset, Contrast, FamilyID, Clinical_Subject_ID, .keep_all = TRUE)

label_table <- bind_rows(lapply(names(prepared_list), function(dataset_name) {
  prepared_list[[dataset_name]]$dictionary %>% mutate(Dataset = dataset_name)
}))
all_results <- all_results %>% left_join(label_table, by = c("Dataset", "FeatureID"))

selection <- all_results %>%
  filter(EffectType == "Overall response", Contrast %in% acute_contrasts) %>%
  group_by(Dataset, FeatureID, FeatureLabel) %>%
  summarise(
    AnyFDR05 = any(BH_FDR < fdr_threshold, na.rm = TRUE),
    MinimumFDR = min(BH_FDR, na.rm = TRUE),
    MinimumP = min(PValue, na.rm = TRUE),
    MaximumAbsT = max(abs(ModeratedT), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Dataset) %>%
  arrange(desc(AnyFDR05), MinimumFDR, MinimumP, desc(MaximumAbsT), .by_group = TRUE) %>%
  slice_head(n = display_features) %>%
  mutate(DisplayRank = row_number()) %>%
  ungroup()

circle_data <- all_results %>%
  filter(
    EffectType == "Overall response",
    Contrast %in% contrast_levels
  ) %>%
  inner_join(selection %>% select(Dataset, FeatureID, FeatureLabel, DisplayRank, AnyFDR05),
             by = c("Dataset", "FeatureID", "FeatureLabel")) %>%
  mutate(
    Contrast = factor(Contrast, levels = contrast_levels),
    Ring = as.integer(Contrast),
    DisplayT = pmax(-4, pmin(4, ModeratedT))
  )

make_circle <- function(dataset_name, title_text) {
  data <- circle_data %>%
    filter(Dataset == dataset_name) %>%
    arrange(DisplayRank, Contrast) %>%
    mutate(
      FeatureLabel = factor(FeatureLabel, levels = rev(unique(FeatureLabel[order(DisplayRank, decreasing = TRUE)]))),
      FeatureIndex = as.numeric(FeatureLabel)
    )
  feature_key <- data %>%
    distinct(FeatureLabel, FeatureIndex, AnyFDR05) %>%
    arrange(FeatureIndex) %>%
    mutate(
      Angle = 90 - 360 * (FeatureIndex - 0.5) / n(),
      HJust = ifelse(Angle < -90, 1, 0),
      LabelAngle = ifelse(Angle < -90, Angle + 180, Angle)
    )
  ggplot(data, aes(x = FeatureIndex, y = Ring, fill = DisplayT)) +
    geom_tile(width = 0.94, height = 0.9, color = "white", linewidth = 0.18) +
    geom_point(
      data = feature_key %>% filter(AnyFDR05),
      aes(x = FeatureIndex, y = 3.58, shape = "BH-FDR < 0.05"),
      inherit.aes = FALSE,
      size = 1.5, stroke = 0.35, fill = "black", color = "black"
    ) +
    geom_text(
      data = feature_key,
      aes(
        x = FeatureIndex,
        y = 3.82,
        label = FeatureLabel,
        angle = LabelAngle,
        hjust = HJust
      ),
      inherit.aes = FALSE,
      family = "Arial", size = 6 / ggplot2::.pt, color = "black"
    ) +
    coord_polar(start = 0, clip = "off") +
    scale_fill_gradient2(
      low = decrease_color, mid = neutral_color, high = increase_color,
      midpoint = 0, limits = c(-4, 4), name = "Moderated t"
    ) +
    scale_shape_manual(
      values = c("BH-FDR < 0.05" = 21),
      limits = "BH-FDR < 0.05", drop = FALSE, name = NULL
    ) +
    scale_y_continuous(
      limits = c(0.45, 4.18), breaks = 1:3,
      labels = NULL
    ) +
    labs(
      title = title_text,
      subtitle = "Inner: Fast-Pre   Middle: Post1h-Pre   Outer: Post3h-Pre",
      x = NULL, y = NULL
    ) +
    theme_void(base_family = "Arial", base_size = 7) +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 6.5, hjust = 0.5, color = "black"),
      axis.text.y = element_blank(),
      legend.title = element_text(size = 8, face = "bold", color = "black"),
      legend.text = element_text(size = 7, color = "black"),
      plot.margin = margin(16, 32, 16, 32)
    )
}

fig3a <- make_circle("Serum_Metabonomics", "Serum metabolomics") +
  make_circle("Serum_Proteomics", "Serum proteomics") +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(
  file.path(fig03a_dir, "Fig03A_Temporal_Serum_Response_Rings.pdf"),
  fig3a, width = 11.8, height = 5.9, device = grDevices::cairo_pdf
)
ggsave(
  file.path(fig03a_dir, "Fig03A_Temporal_Serum_Response_Rings.png"),
  fig3a, width = 11.8, height = 5.9, dpi = 500, bg = "white"
)

build_trajectory_data <- function(dataset_name) {
  prepared <- prepared_list[[dataset_name]]
  selected <- selection %>% filter(Dataset == dataset_name)
  metadata <- prepared$metadata %>%
    filter(Timepoint %in% time_levels, SampleID %in% colnames(prepared$matrix)) %>%
    arrange(Clinical_Subject_ID, Timepoint)
  matrix <- prepared$matrix[selected$FeatureID, metadata$SampleID, drop = FALSE]
  row_means <- rowMeans(matrix, na.rm = TRUE)
  row_sds <- apply(matrix, 1L, stats::sd, na.rm = TRUE)
  row_sds[!is.finite(row_sds) | row_sds == 0] <- 1
  z_matrix <- sweep(sweep(matrix, 1L, row_means, "-"), 1L, row_sds, "/")
  long <- as.data.frame(z_matrix) %>%
    rownames_to_column("FeatureID") %>%
    pivot_longer(-FeatureID, names_to = "SampleID", values_to = "Z") %>%
    left_join(
      metadata %>% select(SampleID, Clinical_Subject_ID, Timepoint, Adiposity),
      by = "SampleID"
    ) %>%
    left_join(selected %>% select(FeatureID, FeatureLabel), by = "FeatureID")

  profile <- long %>%
    group_by(FeatureID, Timepoint) %>%
    summarise(MeanZ = mean(Z, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Timepoint, values_from = MeanZ) %>%
    select(FeatureID, all_of(time_levels))
  profile_matrix <- as.matrix(profile[, time_levels, drop = FALSE])
  rownames(profile_matrix) <- profile$FeatureID
  k <- min(4L, nrow(profile_matrix))
  clusters <- stats::cutree(stats::hclust(stats::dist(profile_matrix), method = "ward.D2"), k = k)
  cluster_key <- tibble(FeatureID = names(clusters), Cluster = paste("Cluster", clusters))

  long %>%
    left_join(cluster_key, by = "FeatureID") %>%
    group_by(FeatureID, FeatureLabel, Cluster, Timepoint) %>%
    summarise(MeanZ = mean(Z, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Dataset = dataset_name,
      Timepoint = factor(Timepoint, levels = time_levels),
      Cluster = factor(Cluster, levels = paste("Cluster", seq_len(k)))
    )
}

trajectory_data <- bind_rows(
  build_trajectory_data("Serum_Metabonomics"),
  build_trajectory_data("Serum_Proteomics")
)

make_trajectory_plot <- function(dataset_name, title_text) {
  data <- trajectory_data %>% filter(Dataset == dataset_name)
  centroids <- data %>%
    group_by(Cluster, Timepoint) %>%
    summarise(MeanZ = mean(MeanZ, na.rm = TRUE), .groups = "drop")
  ggplot(data, aes(x = Timepoint, y = MeanZ)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    geom_line(
      aes(group = FeatureID),
      linewidth = 0.35, alpha = 0.32, color = "#C8AEC3"
    ) +
    geom_line(
      data = centroids,
      aes(group = Cluster),
      linewidth = 1.05, color = "#7A4E74"
    ) +
    geom_point(data = centroids, size = 1.7, stroke = 0.25, color = "#7A4E74") +
    facet_wrap(~Cluster, nrow = 1) +
    coord_cartesian(ylim = c(-1.25, 1.25)) +
    labs(
      title = title_text,
      x = if (dataset_name == "Serum_Proteomics") "Time point" else NULL,
      y = "Standardized abundance"
    ) +
    theme_classic(base_family = "Arial", base_size = 7) +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0.5, color = "black"),
      axis.title = element_text(size = 8, color = "black"),
      axis.text = element_text(size = 7, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 8, color = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(5, 6, 5, 5)
    )
}

fig3b <- make_trajectory_plot("Serum_Metabonomics", "Serum metabolomics") /
  make_trajectory_plot("Serum_Proteomics", "Serum proteomics") +
  patchwork::plot_layout(heights = c(1, 1))

ggsave(
  file.path(fig03b_dir, "Fig03B_Temporal_Serum_Response_Clusters.pdf"),
  fig3b, width = 11.8, height = 6.8, device = grDevices::cairo_pdf
)
ggsave(
  file.path(fig03b_dir, "Fig03B_Temporal_Serum_Response_Clusters.png"),
  fig3b, width = 11.8, height = 6.8, dpi = 500, bg = "white"
)

summary_table <- all_results %>%
  filter(EffectType == "Overall response", Contrast %in% acute_contrasts) %>%
  group_by(Dataset, Contrast) %>%
  summarise(
    TestedFeatures = n(),
    NominalP05 = sum(PValue < 0.05, na.rm = TRUE),
    BH_FDR05 = sum(BH_FDR < fdr_threshold, na.rm = TRUE),
    MinimumBH_FDR = min(BH_FDR, na.rm = TRUE),
    .groups = "drop"
  )
interaction_summary <- all_results %>%
  filter(EffectType == "Time by adiposity interaction", Contrast %in% acute_contrasts) %>%
  group_by(Dataset, Contrast) %>%
  summarise(
    TestedFeatures = n(),
    NominalP05 = sum(PValue < 0.05, na.rm = TRUE),
    BH_FDR05 = sum(BH_FDR < fdr_threshold, na.rm = TRUE),
    MinimumBH_FDR = min(BH_FDR, na.rm = TRUE),
    .groups = "drop"
  )
pair_summary <- all_pairs %>%
  group_by(Dataset, Contrast) %>%
  summarise(
    PairedSubjects = n_distinct(Clinical_Subject_ID),
    Families = n_distinct(FamilyID),
    Lean = sum(Adiposity == "Lean"),
    Obese = sum(Adiposity == "Obese"),
    .groups = "drop"
  )
parameters <- tibble(
  Parameter = c(
    "AnalysisID", "ReleaseID", "FormalBaseline", "Model",
    "FamilyCorrelation", "FeatureFDR", "DisplayFeatureRule", "DisplayTCap"
  ),
  Value = c(
    analysis_id, release_id, "Pre",
    "Subject-level delta ~ centered Adiposity + centered Role; FamilyID duplicateCorrelation block",
    "Platform- and contrast-specific limma consensus correlation",
    "Benjamini-Hochberg within platform x contrast x effect type",
    paste0("Top ", display_features, " per serum platform; prioritize any acute BH-FDR <0.05, then minimum BH-FDR and P"),
    "Moderated t clipped to [-4, 4] for ring color only"
  )
)

writexl::write_xlsx(
  list(
    Overall_Summary = summary_table,
    Interaction_Summary = interaction_summary,
    Pair_Summary = pair_summary,
    Display_Selection = selection,
    Ring_Data = circle_data,
    Trajectory_Data = trajectory_data,
    Complete_Statistics = all_results,
    Input_Manifest = dataset_registry,
    Parameters = parameters
  ),
  file.path(fig03a_dir, "Fig03AB_SourceData.xlsx")
)

notes <- c(
  "# Figure 3A-B acute serum response",
  "",
  paste0("Analysis ID: `", analysis_id, "`"),
  "",
  "## Statistical definition",
  "",
  "- `Pre` is the formal exercise baseline; `Fast-Pre` is displayed only as a pre-exercise control contrast.",
  "- Changes are calculated within subject before modelling.",
  "- The model adjusts for adiposity and family role, with FamilyID as the limma correlation block.",
  "- BH-FDR is controlled separately within each platform, contrast and effect type.",
  "",
  "## Figure definition",
  "",
  "- Figure 3A color is the moderated t statistic, capped at +/-4 for display.",
  "- A black outer marker identifies a feature with BH-FDR <0.05 in at least one acute contrast.",
  "- Figure 3B shows all-participant standardized descriptive trajectories for the same ranked features; thin lines are features and thick lines are cluster means.",
  "- Figure 3B is descriptive and does not replace the formal time-by-adiposity tests retained in the source workbook.",
  "",
  "## Interpretation boundary",
  "",
  "- If a selected feature does not have BH-FDR <0.05, it is a ranked exploratory display feature and not a confirmed molecular response.",
  "- Serum proteomics effects are on the supplied SD scale; serum metabolomics features retain heterogeneous supplied scales.",
  "- No fold-change comparison is made across platforms."
)
writeLines(notes, file.path(figure3_dir, "METHODS_AND_INTERPRETATION.md"), useBytes = TRUE)

run_log <- c(
  paste0("AnalysisID: ", analysis_id),
  paste0("RunTime: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("RVersion: ", R.version.string),
  paste0("limma: ", packageVersion("limma")),
  paste0("ReleaseID: ", release_id),
  "",
  "Overall response summary:",
  capture.output(print(as.data.frame(summary_table), row.names = FALSE)),
  "",
  "Interaction summary:",
  capture.output(print(as.data.frame(interaction_summary), row.names = FALSE)),
  "",
  "Pair summary:",
  capture.output(print(as.data.frame(pair_summary), row.names = FALSE)),
  "",
  "Input SHA256:",
  capture.output(print(as.data.frame(dataset_registry %>% select(Dataset, SHA256)), row.names = FALSE)),
  "",
  capture.output(sessionInfo())
)
writeLines(run_log, file.path(figure3_dir, "Fig03AB_RunLog.txt"), useBytes = TRUE)

message("Completed Figure 3A-B acute serum response analysis.")
