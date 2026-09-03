# Candidate Figure 6B: matched-cohort clinical associations of pre-exercise
# and Post 3-h osteocalcin measures. Existing Figure 6 artwork is not replaced.

.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.this_script <- if (length(.file_arg)) sub("^--file=", "", .file_arg[[1]]) else NA_character_
.script_dir <- if (!is.na(.this_script)) dirname(normalizePath(.this_script)) else file.path(getwd(), "code")
.config_file <- file.path(.script_dir, "Fig00_Config.R")
if (!file.exists(.config_file)) stop("Cannot locate Fig00_Config.R.")
source(.config_file)
rm(.file_arg, .this_script, .script_dir, .config_file)

p_load(dplyr, ggplot2, lme4, lmerTest, patchwork, readr, scales, tibble, tidyr, writexl)

analysis_id <- "Fig06B_Pre_Post_Clinical_Associations_Candidate_2026-08-30"
output_dir <- path_project("fig06", "analysis", analysis_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

clinical_file <- path_analysis_ready("Metadata", "Clinical_Observed_Primary.csv")
registry_file <- path_analysis_ready("Metadata", "Clinical_Variable_Registry.csv")
longitudinal_file <- file.path(
  path_intermediate("Clinical_Exercise_Source_Audit_2026-08-21"),
  "Restricted", "OC_BCAA_Serotonin_Long_Internal.csv"
)
for (required_file in c(clinical_file, registry_file, longitudinal_file)) {
  if (!file.exists(required_file)) stop("Required input is missing: ", required_file)
}

clinical <- read_csv(clinical_file, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    FamilyID = factor(as.character(FamilyID)),
    Role = factor(Role, levels = c("Daughter", "Mother", "Father"))
  )
if (nrow(clinical) != 81L || n_distinct(clinical$FamilyID) != 27L) {
  stop("Expected the locked 81-participant, 27-family clinical cohort.")
}

registry <- read_csv(registry_file, show_col_types = FALSE, progress = FALSE)

marker_registry <- tribble(
  ~Marker, ~MarkerLabel, ~MarkerShort, ~Transform, ~MarkerOrder, ~Color,
  "cOC", "Carboxylated osteocalcin", "cOC", "log", 1L, "#5AA79D",
  "TotalOC", "Total osteocalcin", "tOC", "log", 2L, "#D47E70",
  "cOC_ratio", "cOC/tOC ratio", "cOC/tOC", "logit", 3L, "#8E88B9"
)

excluded_outcomes <- c(
  "S_TotalOC_84m_G", "S_cOC_84m_G", "cOC_totalOC_ratio_84m_G",
  "LntOC_G", "LncOC_G", "FamilyID", "Membercode", "Role",
  "Clinical_Subject_ID"
)

outcome_registry <- registry %>%
  filter(
    ReleaseRole == "Curated observed or source-derived phenotype",
    SourceDataType == "numeric",
    SourceVariable %in% names(clinical),
    !SourceVariable %in% excluded_outcomes
  ) %>%
  transmute(
    Outcome = SourceVariable,
    Phenotype = if_else(is.na(VariableLabel) | VariableLabel == "", SourceVariable, VariableLabel),
    Domain = if_else(is.na(Domain) | Domain == "", "Other", Domain),
    Unit = if_else(is.na(Unit) | Unit == "", "not documented", Unit),
    DomainOrder = suppressWarnings(as.numeric(DomainOrder)),
    VariableOrder = suppressWarnings(as.numeric(VariableOrder))
  ) %>%
  distinct(Outcome, .keep_all = TRUE) %>%
  arrange(DomainOrder, VariableOrder, Outcome)

longitudinal <- read_csv(longitudinal_file, show_col_types = FALSE, progress = FALSE) %>%
  filter(
    MatchStatus == "Unique locked-cohort match",
    Marker %in% marker_registry$Marker,
    Timepoint %in% c("Pre", "Post1h", "Post3h")
  ) %>%
  transmute(
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Marker = as.character(Marker),
    Timepoint = as.character(Timepoint),
    Value = suppressWarnings(as.numeric(Value))
  ) %>%
  filter(is.finite(Value), Value > 0) %>%
  distinct(Clinical_Subject_ID, Marker, Timepoint, .keep_all = TRUE)

complete_ids <- longitudinal %>%
  count(Clinical_Subject_ID, Marker, name = "ObservedTimepoints") %>%
  filter(ObservedTimepoints == 3L) %>%
  count(Clinical_Subject_ID, name = "CompleteMarkers") %>%
  filter(CompleteMarkers == 3L) %>%
  pull(Clinical_Subject_ID)

state_values <- longitudinal %>%
  filter(Clinical_Subject_ID %in% complete_ids, Timepoint %in% c("Pre", "Post3h")) %>%
  select(Clinical_Subject_ID, Marker, Timepoint, Value) %>%
  pivot_wider(names_from = c(Timepoint, Marker), values_from = Value)

analysis_data <- clinical %>%
  inner_join(state_values, by = "Clinical_Subject_ID") %>%
  droplevels()
if (nrow(analysis_data) != 50L || n_distinct(analysis_data$FamilyID) != 18L) {
  stop("Expected 50 matched participants from 18 families.")
}

rank_normalize <- function(x) {
  result <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  n <- sum(keep)
  if (n < 3L) return(result)
  result[keep] <- qnorm((rank(x[keep], ties.method = "average") - 0.5) / n)
  result
}

z_score <- function(x) {
  result <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  if (sum(keep) < 3L || sd(x[keep]) == 0) return(result)
  result[keep] <- as.numeric(scale(x[keep]))
  result
}

transform_marker <- function(x, method) {
  x <- suppressWarnings(as.numeric(x))
  if (method == "log") return(ifelse(is.finite(x) & x > 0, log(x), NA_real_))
  if (method == "logit") {
    clipped <- pmin(pmax(x, 1e-4), 1 - 1e-4)
    return(ifelse(is.finite(x) & x > 0 & x < 1, qlogis(clipped), NA_real_))
  }
  x
}

fit_one <- function(state, marker, marker_label, marker_short, marker_transform, outcome) {
  predictor_column <- paste0(state, "_", marker)
  dat <- tibble(
    FamilyID = analysis_data$FamilyID,
    Role = analysis_data$Role,
    PredictorRaw = transform_marker(analysis_data[[predictor_column]], marker_transform),
    OutcomeRaw = suppressWarnings(as.numeric(analysis_data[[outcome]]))
  ) %>%
    filter(
      is.finite(PredictorRaw), is.finite(OutcomeRaw),
      !is.na(FamilyID), !is.na(Role)
    ) %>%
    mutate(
      PredictorZ = z_score(PredictorRaw),
      OutcomeZ = rank_normalize(OutcomeRaw)
    ) %>%
    filter(is.finite(PredictorZ), is.finite(OutcomeZ)) %>%
    droplevels()

  empty_result <- function() tibble(
    State = state, Marker = marker, MarkerLabel = marker_label,
    MarkerShort = marker_short, Transform = marker_transform,
    Outcome = outcome, EffectiveN = nrow(dat), Families = n_distinct(dat$FamilyID),
    StandardizedEffect = NA_real_, SE = NA_real_, CI_Low = NA_real_,
    CI_High = NA_real_, P_Value = NA_real_, SingularFit = NA
  )
  if (
    nrow(dat) < 20L || n_distinct(dat$FamilyID) < 10L ||
    n_distinct(dat$Role) < 2L || n_distinct(dat$OutcomeRaw) < 5L
  ) return(empty_result())

  fit <- tryCatch(
    suppressMessages(suppressWarnings(lmerTest::lmer(
      OutcomeZ ~ PredictorZ + Role + (1 | FamilyID),
      data = dat, REML = TRUE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    ))),
    error = function(e) NULL
  )
  if (is.null(fit)) return(empty_result())
  co <- summary(fit)$coefficients["PredictorZ", ]
  tibble(
    State = state, Marker = marker, MarkerLabel = marker_label,
    MarkerShort = marker_short, Transform = marker_transform,
    Outcome = outcome, EffectiveN = nrow(dat), Families = n_distinct(dat$FamilyID),
    StandardizedEffect = unname(co["Estimate"]),
    SE = unname(co["Std. Error"]),
    CI_Low = unname(co["Estimate"] - 1.96 * co["Std. Error"]),
    CI_High = unname(co["Estimate"] + 1.96 * co["Std. Error"]),
    P_Value = unname(co["Pr(>|t|)"]),
    SingularFit = lme4::isSingular(fit, tol = 1e-5)
  )
}

results <- tidyr::crossing(
  State = c("Pre", "Post3h"),
  marker_registry %>% select(Marker, MarkerLabel, MarkerShort, Transform),
  outcome_registry %>% select(Outcome)
) %>%
  purrr::pmap_dfr(function(State, Marker, MarkerLabel, MarkerShort, Transform, Outcome) {
    fit_one(State, Marker, MarkerLabel, MarkerShort, Transform, Outcome)
  }) %>%
  left_join(outcome_registry, by = "Outcome") %>%
  group_by(State, Marker) %>%
  mutate(BH_FDR_WithinStateMarker = p.adjust(P_Value, method = "BH")) %>%
  ungroup() %>%
  group_by(State) %>%
  mutate(BH_FDR_WithinStateJoint = p.adjust(P_Value, method = "BH")) %>%
  ungroup() %>%
  mutate(BH_FDR_AllStatesJoint = p.adjust(P_Value, method = "BH")) %>%
  arrange(State, Marker, BH_FDR_WithinStateMarker, P_Value)

display_registry <- tribble(
  ~Outcome, ~DomainDisplay, ~PhenotypeDisplay, ~DisplayOrder,
  "TRACP5b_84m_G", "Bone resorption", "TRACP5b", 1L,
  "ALPPlus", "Liver enzymes", "Alkaline phosphatase", 2L,
  "Fat_p_G", "Adiposity", "Body-fat percentage", 3L,
  "Total_FM_2010_G", "Adiposity", "Total fat mass", 4L,
  "ArMs_FM_2010_G", "Adiposity", "Arm fat mass", 5L,
  "Trunk_FM_2010_G", "Adiposity", "Trunk fat mass", 6L,
  "Android_FM_2010_G", "Adiposity", "Android fat mass", 7L,
  "SubcutaneousFatKg_G", "Adiposity", "Subcutaneous fat mass", 8L,
  "IntraperitonealFatKg_G", "Adiposity", "Visceral fat mass", 9L,
  "ALT", "Liver enzymes", "ALT", 10L
)

selected_fdr_outcomes <- results %>%
  filter(is.finite(BH_FDR_WithinStateMarker), BH_FDR_WithinStateMarker < 0.05) %>%
  distinct(Outcome) %>%
  pull(Outcome)
if (!setequal(selected_fdr_outcomes, display_registry$Outcome)) {
  stop(
    "The displayed registry must equal the union of FDR-significant clinical outcomes. ",
    "Missing: ", paste(setdiff(selected_fdr_outcomes, display_registry$Outcome), collapse = ", "),
    "; extra: ", paste(setdiff(display_registry$Outcome, selected_fdr_outcomes), collapse = ", ")
  )
}

sig_label <- function(x) case_when(
  is.na(x) ~ "",
  x < 0.001 ~ "***",
  x < 0.01 ~ "**",
  x < 0.05 ~ "*",
  TRUE ~ ""
)

plot_data <- results %>%
  inner_join(display_registry, by = "Outcome") %>%
  left_join(marker_registry %>% select(Marker, MarkerOrder, Color), by = "Marker") %>%
  mutate(
    MarkerShort = factor(MarkerShort, levels = c("cOC", "tOC", "cOC/tOC")),
    DomainDisplay = factor(DomainDisplay, levels = c("Bone resorption", "Adiposity", "Liver enzymes")),
    PhenotypeDisplay = factor(PhenotypeDisplay, levels = rev(display_registry$PhenotypeDisplay)),
    Significance = sig_label(BH_FDR_WithinStateMarker)
  ) %>%
  arrange(State, DisplayOrder, MarkerOrder)
if (nrow(plot_data) != 60L) stop("Expected 60 displayed state-marker-phenotype rows.")

effect_limit <- max(0.8, ceiling(max(abs(c(plot_data$CI_Low, plot_data$CI_High)), na.rm = TRUE) * 10) / 10)
effect_limit <- min(effect_limit, 1.4)
effect_breaks <- seq(-effect_limit, effect_limit, length.out = 5L)

common_theme <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.text = element_text(color = "black", size = 6.2),
    axis.title = element_text(color = "black", size = 6.8, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    strip.background = element_blank(),
    panel.spacing.y = grid::unit(0.7, "mm")
  )

make_state_figure <- function(state_name) {
  dat <- plot_data %>% filter(State == state_name)
  state_title <- if (state_name == "Pre") {
    "Associations of osteocalcin with clinical traits (pre-exercise)"
  } else {
    "Associations of osteocalcin with clinical traits (Post 3 h)"
  }

  heatmap <- ggplot(dat, aes(MarkerShort, PhenotypeDisplay, fill = StandardizedEffect)) +
    geom_tile(color = "black", linewidth = 0.18, width = 0.88, height = 0.88) +
    geom_text(
      aes(label = Significance), color = "black", family = "Arial",
      size = 2.45, fontface = "bold"
    ) +
    facet_grid(rows = vars(DomainDisplay), scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2(
      low = "#3E6E9E", mid = "#FFFFFF", high = "#C8564B",
      midpoint = 0, limits = c(-effect_limit, effect_limit), oob = squish,
      breaks = effect_breaks, name = "Standardized\neffect"
    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    labs(title = "Association matrix", x = NULL, y = NULL) +
    common_theme +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, color = "black", size = 6.0, face = "bold"),
      axis.text.x.top = element_text(color = "black", size = 6.3, face = "bold"),
      axis.text.y = element_text(color = "black", size = 6.0),
      axis.ticks = element_blank(), axis.line = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )

  forest <- ggplot(dat, aes(StandardizedEffect, PhenotypeDisplay, color = MarkerShort)) +
    geom_vline(xintercept = 0, linewidth = 0.22, linetype = "dashed", color = "black") +
    geom_errorbarh(aes(xmin = CI_Low, xmax = CI_High), height = 0, linewidth = 0.34, color = "black") +
    geom_point(size = 1.8) +
    facet_grid(
      rows = vars(DomainDisplay), cols = vars(MarkerShort),
      scales = "free_y", space = "free_y"
    ) +
    scale_color_manual(values = c(tOC = "#D47E70", cOC = "#5AA79D", `cOC/tOC` = "#8E88B9"), guide = "none") +
    scale_x_continuous(
      limits = c(-effect_limit, effect_limit), breaks = effect_breaks,
      labels = label_number(accuracy = 0.1), expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "Standardized effect", y = NULL) +
    common_theme +
    theme(
      strip.text.x = element_text(color = "black", size = 6.3, face = "bold"),
      strip.text.y = element_blank(), axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), axis.line.y = element_blank(),
      axis.text.x = element_text(color = "black", size = 5.5),
      axis.title.x = element_text(color = "black", size = 6.3, face = "bold"),
      plot.margin = margin(2, 1.5, 2, 1.5, unit = "mm")
    )

  legend_gradient <- tibble(x = 0, y = seq(-effect_limit, effect_limit, length.out = 161))
  legend <- ggplot(legend_gradient, aes(x, y, fill = y)) +
    geom_tile(width = 0.28, height = diff(range(effect_breaks)) / 160) +
    scale_fill_gradient2(
      low = "#3E6E9E", mid = "#FFFFFF", high = "#C8564B",
      midpoint = 0, limits = c(-effect_limit, effect_limit), guide = "none"
    ) +
    annotate(
      "text", x = 0.14, y = effect_limit * 1.35,
      label = "Standardized\neffect", family = "Arial",
      size = 2.15, fontface = "bold", hjust = 0.5
    ) +
    annotate(
      "text", x = -0.05, y = -effect_limit * 1.35,
      label = "BH-FDR\n* < 0.05\n** < 0.01\n*** < 0.001",
      family = "Arial", size = 1.95, hjust = 0, vjust = 1, lineheight = 1.15
    ) +
    coord_cartesian(
      xlim = c(-0.12, 0.95),
      ylim = c(-effect_limit * 2.1, effect_limit * 1.55), clip = "off"
    ) +
    theme_void(base_family = "Arial") +
    theme(plot.margin = margin(6, 1, 6, 4, unit = "mm"))

  heatmap + forest + legend +
    plot_layout(widths = c(2.7, 4.45, 0.75)) +
    plot_annotation(
      title = state_title,
      theme = theme(
        plot.title = element_text(
          family = "Arial", color = "black", size = 8.5,
          face = "bold", hjust = 0
        )
      )
    )
}

pre_figure <- make_state_figure("Pre")
post_figure <- make_state_figure("Post3h")
comparison <- wrap_elements(full = pre_figure) / wrap_elements(full = post_figure) +
  plot_layout(heights = c(1, 1))

pre_pdf <- file.path(output_dir, "Fig06B_PreExercise_Clinical_Associations_Candidate.pdf")
post_pdf <- file.path(output_dir, "Fig06B_Post3h_Clinical_Associations_Candidate.pdf")
comparison_pdf <- file.path(output_dir, "Fig06B_Pre_vs_Post3h_Clinical_Associations_Comparison.pdf")
ggsave(pre_pdf, pre_figure, width = 190, height = 94, units = "mm", device = cairo_pdf)
ggsave(post_pdf, post_figure, width = 190, height = 94, units = "mm", device = cairo_pdf)
ggsave(comparison_pdf, comparison, width = 190, height = 178, units = "mm", device = cairo_pdf)

summary_table <- results %>%
  group_by(State, MarkerShort) %>%
  summarise(
    TestedOutcomes = sum(is.finite(P_Value)),
    NominalP05 = sum(P_Value < 0.05, na.rm = TRUE),
    BH_FDR05 = sum(BH_FDR_WithinStateMarker < 0.05, na.rm = TRUE),
    BH_FDR10 = sum(BH_FDR_WithinStateMarker < 0.10, na.rm = TRUE),
    BH_FDR20 = sum(BH_FDR_WithinStateMarker < 0.20, na.rm = TRUE),
    MinimumBH_FDR = min(BH_FDR_WithinStateMarker, na.rm = TRUE),
    .groups = "drop"
  )

methods <- tribble(
  ~Item, ~Description,
  "Analysis cohort", "Matched subset of 50 participants from 18 families with complete Pre, Post 1 h, and Post 3 h tOC, cOC, and cOC/tOC values.",
  "Model", "Rank-normalized clinical phenotype ~ standardized state-specific osteocalcin measure + Role + (1 | FamilyID).",
  "Marker transformations", "Natural-log tOC and cOC; logit cOC/tOC ratio; each transformed marker was standardized within the analyzed state.",
  "Clinical outcomes", paste0(nrow(outcome_registry), " curated numeric clinical phenotypes were screened."),
  "Multiplicity", "BH correction was applied separately for each state and osteocalcin measure across all eligible clinical phenotypes.",
  "Display", "The same ten clinical phenotypes forming the union of BH-FDR-significant outcomes across Pre and Post 3-h states are shown in both figures.",
  "Interpretation", "Post 3-h absolute associations do not test exercise response because the Post value includes stable between-participant differences present before exercise."
)

sheet_data <- list(
  Summary = summary_table,
  Displayed_results = plot_data %>%
    select(
      State, DomainDisplay, PhenotypeDisplay, Outcome, MarkerShort,
      EffectiveN, Families, StandardizedEffect, SE, CI_Low, CI_High,
      P_Value, BH_FDR_WithinStateMarker, BH_FDR_WithinStateJoint,
      BH_FDR_AllStatesJoint, SingularFit
    ),
  Significant_FDR05 = results %>% filter(BH_FDR_WithinStateMarker < 0.05),
  All_results = results,
  Outcome_registry = outcome_registry,
  Methods = methods
)
writexl::write_xlsx(
  sheet_data,
  file.path(output_dir, "Fig06B_Pre_Post_Clinical_Associations_Source_Data.xlsx")
)

writeLines(
  c(
    "# Figure 6B pre/post clinical association candidates",
    "",
    "Pre-exercise and Post 3-h candidates use the same matched participants and the same displayed clinical phenotypes.",
    "",
    "Post 3-h absolute associations describe post-exercise state. They are not change-score or baseline-adjusted response models."
  ),
  file.path(output_dir, "README.md"), useBytes = TRUE
)

writeLines(
  c(
    paste0("Analysis ID: ", analysis_id),
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Participants: ", nrow(analysis_data)),
    paste0("Families: ", n_distinct(analysis_data$FamilyID)),
    paste0("Clinical outcomes screened: ", nrow(outcome_registry)),
    "",
    capture.output(sessionInfo())
  ),
  file.path(output_dir, "RUN_LOG.txt"), useBytes = TRUE
)

message("Created Figure 6B pre/post clinical association candidates in: ", output_dir)
