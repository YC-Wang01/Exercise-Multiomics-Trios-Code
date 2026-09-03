# Figure 4 candidate scatterplots for Post3h molecular profiles associated
# with baseline Matsuda ISI. Raw/source inputs are read only.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  file.path(getwd(), "Fig00_Config.R"),
  if (!is.na(.local_script)) file.path(dirname(normalizePath(
    .local_script, winslash = "/", mustWork = FALSE
  )), "Fig00_Config.R") else NA_character_
))
.config_file <- .config_candidates[file.exists(.config_candidates)][1]
if (is.na(.config_file)) stop("Cannot locate Fig00_Config.R.")
source(.config_file)
rm(.local_script, .config_candidates, .config_file)

p_load(dplyr, ggplot2, patchwork, readr, scales, stringr, tibble, tidyr, writexl)

analysis_id <- "Fig04D_Post3h_Candidate_Scatterplots_2026-08-28"
output_dir <- path_project("fig04", "fig04D", "Post3h_Candidate_Scatterplots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

role_levels <- c("Daughter", "Mother", "Father")
role_shapes <- c(Daughter = 24, Mother = 21, Father = 22)
body_fat_colors <- c("#6D9FC4", "#F4F1EC", "#D77A72")

rank_normalize <- function(x) {
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  n <- sum(keep)
  if (n >= 3L) {
    out[keep] <- qnorm((rank(x[keep], ties.method = "average") - 0.5) / n)
  }
  out
}

clinical_file <- path_analysis_ready("Metadata", "Clinical_Observed_Primary.csv")
sample_sheet_file <- path_analysis_ready("Metadata", "Project_Sample_Sheet.csv")
formal_results_file <- path_project(
  "Pre_Post_Response_Analysis", "Results", "02_NonMethylation",
  "NonMethylation_Feature_Selected.csv"
)

required_files <- c(clinical_file, sample_sheet_file, formal_results_file)
if (any(!file.exists(required_files))) {
  stop("Missing required input(s): ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}

clinical <- read_csv(clinical_file, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = factor(Role, levels = role_levels),
    BodyFat = suppressWarnings(as.numeric(Fat_p_G)),
    Matsuda_ISI = suppressWarnings(as.numeric(ISI_G))
  ) %>%
  distinct(Clinical_Subject_ID, .keep_all = TRUE) %>%
  mutate(
    BodyFat_Z = rank_normalize(BodyFat),
    Matsuda_ISI_Z = rank_normalize(if_else(Matsuda_ISI > 0, Matsuda_ISI, NA_real_))
  )

if (nrow(clinical) != 81L || n_distinct(clinical$FamilyID) != 27L) {
  stop("Clinical metadata does not match the locked 81-participant, 27-family cohort.")
}

sample_sheet <- read_csv(sample_sheet_file, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    SampleID = tolower(as.character(SampleID)),
    Role = factor(Role, levels = role_levels),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort)
  ) %>%
  filter(AnalysisSet == "Primary", Is_Primary_Cohort)

candidates <- tribble(
  ~Panel, ~Dataset, ~InputFile, ~FeatureID, ~DisplayLabel, ~LayerLabel, ~YAxisLabel,
  "A", "Serum_Metabonomics", "Analysis_Serum_Metabonomics.csv", "TGtoPG", "TG/PG ratio", "Serum metabolomics", "Role-adjusted Post3h TG/PG ratio",
  "B", "Serum_Proteomics", "Analysis_Serum_Proteomics.csv", "P08603;Q9NSY1", "CFH", "Serum proteomics", "Role-adjusted Post3h log2 abundance",
  "C", "Adipose_Proteomics", "Analysis_Adipose_Proteomics.csv", "Q9H4G4", "GLIPR2", "Adipose proteomics", "Role-adjusted Post3h log2 abundance",
  "D", "Adipose_Proteomics", "Analysis_Adipose_Proteomics.csv", "Q5ZPR3", "CD276", "Adipose proteomics", "Role-adjusted Post3h log2 abundance"
)

formal_results <- read_csv(formal_results_file, show_col_types = FALSE, progress = FALSE) %>%
  filter(
    State == "Post3hAbsolute",
    Predictor == "Matsuda_ISI"
  ) %>%
  mutate(
    Effect = as.numeric(Effect),
    P_Value = as.numeric(P_Value),
    BH_FDR = as.numeric(BH_FDR),
    ObservedN = as.integer(ObservedN),
    Subjects = as.integer(Subjects),
    Families = as.integer(Families)
  ) %>%
  select(Dataset, FeatureID, Effect, P_Value, BH_FDR, ObservedN, Subjects, Families)

format_probability <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < 0.001) return(format(x, scientific = TRUE, digits = 2))
  sprintf("%.3f", x)
}

nice_breaks_exact <- function(values, n = 5L) {
  values <- values[is.finite(values)]
  if (length(values) < 2L || diff(range(values)) == 0) return(pretty(values, n = n))
  value_range <- range(values)
  exponent <- floor(log10(diff(value_range) / (n - 1L)))
  candidate_steps <- sort(unique(as.vector(outer(
    c(1, 1.5, 2, 2.5, 5),
    10 ^ seq(exponent - 1L, exponent + 2L),
    `*`
  ))))
  for (step in candidate_steps) {
    start <- floor(value_range[1] / step) * step
    breaks <- start + seq.int(0L, n - 1L) * step
    if (breaks[n] >= value_range[2]) return(breaks)
  }
  pretty(values, n = n)
}

extract_candidate <- function(spec) {
  matrix_file <- path_analysis_ready(spec$InputFile)
  if (!file.exists(matrix_file)) stop("Missing matrix: ", matrix_file)
  matrix <- read_analysis_matrix(matrix_file)
  if (!spec$FeatureID %in% rownames(matrix)) {
    stop("Feature not present in ", spec$Dataset, ": ", spec$FeatureID)
  }

  available <- sample_sheet %>%
    filter(
      Dataset == spec$Dataset,
      Timepoint %in% c("Pre", "Post3h"),
      SampleID %in% colnames(matrix)
    ) %>%
    inner_join(clinical, by = c("FamilyID", "Clinical_Subject_ID", "Role")) %>%
    distinct(Clinical_Subject_ID, Timepoint, .keep_all = TRUE) %>%
    select(
      FamilyID, Clinical_Subject_ID, Role, BodyFat, BodyFat_Z,
      Matsuda_ISI, Matsuda_ISI_Z, Timepoint, SampleID
    ) %>%
    pivot_wider(names_from = Timepoint, values_from = SampleID) %>%
    filter(!is.na(Pre), !is.na(Post3h)) %>%
    arrange(FamilyID, Role, Clinical_Subject_ID)

  values <- as.numeric(matrix[spec$FeatureID, available$Post3h])
  data <- available %>%
    mutate(Post3hValue = values) %>%
    filter(is.finite(Post3hValue), is.finite(Matsuda_ISI_Z), !is.na(Role), !is.na(FamilyID))

  if (nrow(data) < 10L) stop("Insufficient observations for ", spec$DisplayLabel)

  role_fit <- lm(Post3hValue ~ Role, data = data)
  data <- data %>%
    mutate(
      RoleAdjustedPost3h = residuals(role_fit) + mean(Post3hValue, na.rm = TRUE),
      Panel = spec$Panel,
      Dataset = spec$Dataset,
      FeatureID = spec$FeatureID,
      DisplayLabel = spec$DisplayLabel,
      LayerLabel = spec$LayerLabel,
      YAxisLabel = spec$YAxisLabel
    )

  stats <- formal_results %>%
    filter(Dataset == spec$Dataset, FeatureID == spec$FeatureID) %>%
    slice_head(n = 1L)
  if (nrow(stats) != 1L) stop("Formal model result not found for ", spec$DisplayLabel)
  if (nrow(data) != stats$ObservedN) {
    stop(
      "Observed N mismatch for ", spec$DisplayLabel,
      ": plot data n=", nrow(data), ", formal model n=", stats$ObservedN
    )
  }

  stats <- stats %>%
    mutate(
      Panel = spec$Panel,
      DisplayLabel = spec$DisplayLabel,
      LayerLabel = spec$LayerLabel,
      ObservedFamilies = n_distinct(data$FamilyID),
      Intercept = mean(data$RoleAdjustedPost3h) - Effect * mean(data$Matsuda_ISI_Z),
      Annotation = paste0(
        "\u03b2 = ", sprintf("%+.3f", Effect),
        "\nBH-FDR = ", format_probability(BH_FDR)
      )
    )
  list(data = data, stats = stats)
}

extracted <- lapply(seq_len(nrow(candidates)), function(i) extract_candidate(candidates[i, ]))
plot_data <- bind_rows(lapply(extracted, `[[`, "data"))
model_stats <- bind_rows(lapply(extracted, `[[`, "stats"))

make_panel <- function(panel_id) {
  dat <- plot_data %>% filter(Panel == panel_id)
  stat <- model_stats %>% filter(Panel == panel_id)
  y_breaks <- nice_breaks_exact(dat$RoleAdjustedPost3h, n = 5L)

  ggplot(dat, aes(Matsuda_ISI_Z, RoleAdjustedPost3h)) +
    geom_abline(
      intercept = stat$Intercept, slope = stat$Effect,
      linewidth = 0.55, color = "#3F4A54"
    ) +
    geom_point(
      aes(fill = BodyFat, shape = Role),
      size = 2.25, stroke = 0.45, color = "#30353A", alpha = 0.95
    ) +
    scale_shape_manual(values = role_shapes, drop = FALSE) +
    scale_fill_gradientn(
      colours = body_fat_colors,
      limits = range(clinical$BodyFat, na.rm = TRUE),
      oob = squish,
      name = "Body fat (%)"
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      limits = range(y_breaks),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = stat$DisplayLabel,
      subtitle = paste0(stat$LayerLabel, "\n", stat$Annotation),
      x = NULL,
      y = if (panel_id == "A") "Role-adjusted Post3h molecular value" else NULL
    ) +
    theme_classic(base_family = "Arial", base_size = 7.2) +
    theme(
      plot.title = element_text(face = "bold", size = 7.8, hjust = 0),
      plot.subtitle = element_text(size = 6.2, color = "#4F575D", hjust = 0, lineheight = 1.05),
      axis.title = element_text(face = "bold", size = 6.8, color = "black"),
      axis.text = element_text(size = 6.3, color = "black"),
      axis.line = element_line(linewidth = 0.4, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.4),
      legend.title = element_text(face = "bold", size = 6.5),
      legend.text = element_text(size = 6.2),
      legend.key.height = unit(3.4, "mm"),
      legend.key.width = unit(3.4, "mm"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

panels <- lapply(candidates$Panel, make_panel)
combined_two_by_two <- wrap_plots(panels, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Post-exercise molecular profiles associated with baseline Matsuda ISI",
    subtitle = paste0(
      "Role-adjusted partial scatterplots; formal inference used limma with ",
      "FamilyID duplicateCorrelation"
    ),
    tag_levels = "A",
    theme = theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(face = "bold", size = 10, hjust = 0),
      plot.subtitle = element_text(size = 7, hjust = 0, color = "#4F575D")
    )
  ) &
  theme(legend.position = "right")

combined_one_row <- wrap_plots(panels, nrow = 1, guides = "collect") +
  plot_annotation(
    title = "D  Post3h molecular profiles associated with baseline Matsuda ISI",
    caption = "Baseline Matsuda ISI (rank-normalized SD)",
    theme = theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(face = "bold", size = 9.2, hjust = 0),
      plot.caption = element_text(face = "bold", size = 6.8, hjust = 0.46, margin = margin(t = 2))
    )
  ) &
  theme(legend.position = "right")

ggsave(
  file.path(output_dir, "Fig04D_Post3h_Matsuda_Candidate_Scatterplots_2x2.pdf"),
  combined_two_by_two, width = 160, height = 132, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Fig04D_Post3h_Matsuda_Candidate_Scatterplots_2x2.png"),
  combined_two_by_two, width = 160, height = 132, units = "mm", dpi = 400, bg = "white"
)
ggsave(
  file.path(output_dir, "Fig04D_Post3h_Matsuda_Candidate_Scatterplots_OneRow_190mm.pdf"),
  combined_one_row, width = 190, height = 61, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Fig04D_Post3h_Molecular_Profiles.pdf"),
  combined_one_row, width = 190, height = 61, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Fig04D_Post3h_Molecular_Profiles.png"),
  combined_one_row, width = 190, height = 61, units = "mm", dpi = 500, bg = "white"
)

public_plot_data <- plot_data %>%
  group_by(Panel) %>%
  arrange(Matsuda_ISI_Z, .by_group = TRUE) %>%
  mutate(PlotPointID = sprintf("%s%02d", Panel, row_number())) %>%
  ungroup() %>%
  select(
    Panel, PlotPointID, DisplayLabel, LayerLabel, Role, BodyFat,
    Matsuda_ISI, Matsuda_ISI_Z, Post3hValue, RoleAdjustedPost3h
  )

readme_table <- tribble(
  ~Field, ~Value,
  "Figure purpose", "Visualize Post3h molecular profiles associated with baseline Matsuda ISI.",
  "Population", "Locked primary cohort with observed paired Pre and Post3h platform data.",
  "X axis", "Matsuda ISI rank-normalized across the locked clinical cohort.",
  "Y axis", "Post3h analysis-scale value after removal of the fixed Role effect for visualization.",
  "Line", "Slope equals the formal role-adjusted limma coefficient; FamilyID was handled by duplicateCorrelation in inference.",
  "Color", "Continuous body-fat percentage.",
  "Shape", "Family role.",
  "Multiplicity", "BH correction within dataset x state x predictor.",
  "Interpretation", "Post3h absolute-profile association; not evidence that exercise changed Matsuda ISI or caused molecular change."
)

write_xlsx(
  list(
    Plot_Data = as.data.frame(public_plot_data),
    Model_Statistics = as.data.frame(model_stats %>% select(
      Panel, Dataset, FeatureID, DisplayLabel, LayerLabel, Effect,
      P_Value, BH_FDR, ObservedN, ObservedFamilies, Subjects, Families
    )),
    Figure_Notes = as.data.frame(readme_table)
  ),
  file.path(output_dir, "Table_Fig04D_Post3h_Candidate_Scatterplots.xlsx")
)

saveRDS(
  list(
    AnalysisID = analysis_id,
    PlotDataPrivate = plot_data,
    ModelStatistics = model_stats,
    CandidateRegistry = candidates
  ),
  file.path(output_dir, "Fig04D_Post3h_Candidate_Scatterplots_Source.rds"),
  compress = "gzip"
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_SESSION_INFO.txt"))
writeLines(
  c(
    "# Figure 4D Post3h candidate scatterplots",
    "",
    "- These are candidate visualizations and do not replace the current Figure 4 panel.",
    "- The x axis is rank-normalized baseline Matsuda ISI, matching the formal model predictor.",
    "- The y axis is the observed Post3h analysis-scale value after removing the fixed Role effect for visualization.",
    "- The plotted line uses the formal limma coefficient; inference accounted for FamilyID using duplicateCorrelation.",
    "- Point fill indicates continuous body-fat percentage; point shape indicates family role.",
    "- Post3h association does not establish an exercise-induced molecular change; the corresponding Delta models were not BH-FDR significant.",
    "",
    "## Candidate results",
    paste0(
      "- ", model_stats$DisplayLabel, " (", model_stats$LayerLabel, "): beta = ",
      sprintf("%+.3f", model_stats$Effect), ", P = ",
      vapply(model_stats$P_Value, format_probability, character(1)), ", BH-FDR = ",
      vapply(model_stats$BH_FDR, format_probability, character(1)),
      ", n = ", model_stats$ObservedN, "."
    )
  ),
  file.path(output_dir, "README.md")
)

message("Completed: ", output_dir)
