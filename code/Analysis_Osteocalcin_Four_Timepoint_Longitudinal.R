# Four-time-point longitudinal analysis of circulating osteocalcin.
# Fast is a fasting reference sample; Pre is the formal exercise baseline.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  file.path(getwd(), "Fig00_Config.R"),
  if (!is.na(.local_script)) file.path(dirname(.local_script), "Fig00_Config.R") else NA_character_
))
.config_file <- .config_candidates[file.exists(.config_candidates)][1]
if (is.na(.config_file)) stop("Cannot locate Fig00_Config.R.")
source(.config_file)
.project_root <- PROJECT_ROOT
rm(.local_script, .config_candidates, .config_file)

p_load(digest, dplyr, emmeans, lme4, lmerTest, readr, tibble, tidyr, writexl)

analysis_id <- "Osteocalcin_Four_Timepoint_Longitudinal_2026-08-31"
module_root <- file.path(
  .project_root, "Pre_Post_Response_Analysis", "Results",
  "08_Osteocalcin_Molecular_State_and_Response"
)
output_dir <- file.path(module_root, "05_Four_Timepoint_Longitudinal_Analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source_file <- file.path(
  path_intermediate("Clinical_Exercise_Source_Audit_2026-08-21"),
  "Restricted", "OC_BCAA_Serotonin_Long_Internal.csv"
)
if (!file.exists(source_file)) stop("Longitudinal clinical source table is missing.")

marker_order <- c("cOC", "TotalOC", "cOC_ratio")
marker_labels <- c(
  TotalOC = "Total osteocalcin (tOC)",
  cOC = "Carboxylated osteocalcin (cOC)",
  cOC_ratio = "cOC/tOC ratio"
)
time_order <- c("Fast", "Pre", "Post1h", "Post3h")
time_labels <- c(
  Fast = "Fasting reference",
  Pre = "Pre-exercise",
  Post1h = "Post 1 h",
  Post3h = "Post 3 h"
)
role_order <- c("Daughter", "Mother", "Father")

clinical <- read_csv(source_file, show_col_types = FALSE, progress = FALSE) %>%
  filter(
    MatchStatus == "Unique locked-cohort match",
    Marker %in% marker_order,
    Timepoint %in% time_order
  ) %>%
  transmute(
    FamilyID = factor(as.character(FamilyID)),
    ParticipantID = factor(as.character(Clinical_Subject_ID)),
    Role = factor(Role, levels = role_order),
    Marker,
    Timepoint = factor(Timepoint, levels = time_order),
    Value = suppressWarnings(as.numeric(Value))
  ) %>%
  filter(is.finite(Value), Value > 0, !is.na(Role)) %>%
  distinct(FamilyID, ParticipantID, Marker, Timepoint, .keep_all = TRUE) %>%
  group_by(Marker, ParticipantID) %>%
  filter(n_distinct(Timepoint) == length(time_order)) %>%
  ungroup() %>%
  mutate(Log2Value = log2(Value))

coverage <- clinical %>%
  count(Marker, Timepoint, name = "ObservedN") %>%
  complete(Marker = marker_order, Timepoint = factor(time_order, levels = time_order), fill = list(ObservedN = 0L)) %>%
  mutate(
    MarkerLabel = unname(marker_labels[Marker]),
    TimepointLabel = unname(time_labels[as.character(Timepoint)])
  ) %>%
  arrange(match(Marker, marker_order), Timepoint)

if (any(coverage$ObservedN != 50L)) {
  stop("Expected 50 complete participants per osteocalcin marker and time point.")
}

subject_key <- clinical %>%
  distinct(ParticipantID) %>%
  arrange(ParticipantID) %>%
  mutate(PlotSubject = sprintf("P%03d", row_number()))

plot_data <- clinical %>%
  left_join(subject_key, by = "ParticipantID") %>%
  group_by(Marker, ParticipantID) %>%
  mutate(
    PreValue = Value[Timepoint == "Pre"][1],
    RelativeToPre = Value / PreValue
  ) %>%
  ungroup()

descriptive <- plot_data %>%
  group_by(Marker, Timepoint) %>%
  summarise(
    N = n(),
    ArithmeticMean = mean(Value),
    SD = sd(Value),
    GeometricMean = 2^mean(Log2Value),
    Median = median(Value),
    Q1 = quantile(Value, 0.25),
    Q3 = quantile(Value, 0.75),
    MeanRelativeToPre = mean(RelativeToPre),
    .groups = "drop"
  ) %>%
  mutate(
    MarkerLabel = unname(marker_labels[Marker]),
    TimepointLabel = unname(time_labels[as.character(Timepoint)])
  ) %>%
  arrange(match(Marker, marker_order), Timepoint)

contrast_definitions <- list(
  Pre_minus_Fast = c(-1, 1, 0, 0),
  Post1h_minus_Pre = c(0, -1, 1, 0),
  Post3h_minus_Pre = c(0, -1, 0, 1),
  Post3h_minus_Post1h = c(0, 0, -1, 1)
)

fit_marker <- function(marker) {
  d <- clinical %>% filter(Marker == marker) %>% droplevels()
  fit <- lmerTest::lmer(
    Log2Value ~ Timepoint + Role + (1 | FamilyID) + (1 | ParticipantID),
    data = d, REML = TRUE
  )
  emm <- emmeans::emmeans(fit, ~ Timepoint)
  contrasts <- emmeans::contrast(emm, method = contrast_definitions, adjust = "none") %>%
    summary(infer = c(TRUE, TRUE), adjust = "none") %>%
    as_tibble() %>%
    transmute(
      Marker = marker,
      Contrast = as.character(contrast),
      Log2Estimate = estimate,
      SE,
      Df = df,
      Lower95Log2 = lower.CL,
      Upper95Log2 = upper.CL,
      T = t.ratio,
      P_Value = p.value,
      FoldChange = 2^Log2Estimate,
      PercentChange = 100 * (FoldChange - 1),
      Lower95Percent = 100 * (2^Lower95Log2 - 1),
      Upper95Percent = 100 * (2^Upper95Log2 - 1),
      Participants = n_distinct(d$ParticipantID),
      Families = n_distinct(d$FamilyID)
    )
  global <- as.data.frame(anova(fit, ddf = "Satterthwaite")) %>%
    rownames_to_column("Term") %>%
    filter(Term == "Timepoint") %>%
    transmute(
      Marker = marker,
      NumDf = NumDF,
      DenDf = DenDF,
      F = `F value`,
      P_Value = `Pr(>F)`
    )
  list(fit = fit, contrasts = contrasts, global = global)
}

fits <- lapply(marker_order, fit_marker)
names(fits) <- marker_order

time_effects <- bind_rows(lapply(fits, `[[`, "contrasts")) %>%
  mutate(
    MarkerLabel = unname(marker_labels[Marker]),
    ContrastFamily = case_when(
      Contrast %in% c("Post1h_minus_Pre", "Post3h_minus_Pre") ~ "Primary exercise contrasts",
      Contrast == "Pre_minus_Fast" ~ "Pre-intervention contrast",
      Contrast == "Post3h_minus_Post1h" ~ "Recovery contrast",
      TRUE ~ "Other"
    ),
    BH_FDR_All12 = p.adjust(P_Value, method = "BH")
  ) %>%
  group_by(ContrastFamily) %>%
  mutate(BH_FDR_WithinFamily = p.adjust(P_Value, method = "BH")) %>%
  ungroup() %>%
  arrange(match(Marker, marker_order), match(Contrast, names(contrast_definitions)))

global_tests <- bind_rows(lapply(fits, `[[`, "global")) %>%
  mutate(
    MarkerLabel = unname(marker_labels[Marker]),
    BH_FDR = p.adjust(P_Value, method = "BH")
  ) %>%
  arrange(match(Marker, marker_order))

model_parameters <- tribble(
  ~Item, ~Description,
  "Population", "50 participants with complete Fast, Pre, Post1h, and Post3h tOC, cOC, and cOC/tOC measurements.",
  "Outcome", "Log2-transformed osteocalcin value.",
  "Model", "log2(osteocalcin) ~ Timepoint + Role + (1 | FamilyID) + (1 | ParticipantID).",
  "Time reference", "Fast was retained as a fasting reference sample; Pre was the formal exercise baseline.",
  "Primary exercise contrasts", "Post1h minus Pre and Post3h minus Pre.",
  "Secondary contrasts", "Pre minus Fast and Post3h minus Post1h.",
  "Multiplicity", "BH-FDR is reported across all 12 marker-by-contrast tests and within each prespecified contrast family.",
  "Causal boundary", "Without a time-matched non-exercise control, estimated differences are within-person exercise-associated changes rather than fully identified causal effects.",
  "Ratio caveat", "The cOC/tOC ratio is algebraically derived from cOC and tOC."
)

plot_export <- plot_data %>%
  transmute(
    PlotSubject,
    Marker,
    MarkerLabel = unname(marker_labels[Marker]),
    Timepoint = as.character(Timepoint),
    TimepointLabel = unname(time_labels[as.character(Timepoint)]),
    AbsoluteValue = Value,
    RelativeToPre
  ) %>%
  arrange(match(Marker, marker_order), PlotSubject, match(Timepoint, time_order))

provenance <- tibble(
  Item = c("Source longitudinal table", "Analysis script"),
  ProjectRelativePath = c(
    file.path("Data", "Processed", "Intermediate", "Clinical_Exercise_Source_Audit_2026-08-21", "Restricted", basename(source_file)),
    file.path("code", "Analysis_Osteocalcin_Four_Timepoint_Longitudinal.R")
  ),
  SHA256 = c(
    digest::digest(file = source_file, algo = "sha256"),
    NA_character_
  )
)

write_csv(time_effects, file.path(output_dir, "OC_Four_Timepoint_Time_Effects.csv"))
write_csv(global_tests, file.path(output_dir, "OC_Four_Timepoint_Global_Tests.csv"))
write_csv(descriptive, file.path(output_dir, "OC_Four_Timepoint_Descriptive_Summary.csv"))

write_xlsx(
  list(
    Time_effects = time_effects,
    Global_tests = global_tests,
    Descriptive_summary = descriptive,
    Coverage = coverage,
    Plot_data = plot_export,
    Model_parameters = model_parameters,
    Provenance = provenance
  ),
  file.path(output_dir, "OC_Four_Timepoint_Longitudinal_Analysis.xlsx")
)

saveRDS(
  list(
    time_effects = time_effects,
    global_tests = global_tests,
    descriptive = descriptive,
    coverage = coverage,
    plot_data = plot_export,
    model_parameters = model_parameters,
    fitted_models = lapply(fits, `[[`, "fit")
  ),
  file.path(output_dir, "OC_Four_Timepoint_Longitudinal_Analysis.rds"),
  compress = "xz"
)

writeLines(
  c(
    "# Four-time-point longitudinal osteocalcin analysis",
    "",
    "This module adds the fasting reference sample to the existing Pre/Post1h/Post3h analysis.",
    "Pre remains the formal exercise baseline. Post1h-versus-Pre and Post3h-versus-Pre are the primary exercise contrasts.",
    "Pre-versus-Fast is a pre-intervention contrast and must not be interpreted as an exercise response.",
    "The complete 50-participant cohort is primary; no participant was excluded from the formal model.",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Analysis ID: ", analysis_id)
  ),
  file.path(output_dir, "README.md"), useBytes = TRUE
)

writeLines(
  c(
    paste0("Analysis ID: ", analysis_id),
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Source SHA256: ", digest::digest(file = source_file, algo = "sha256")),
    "",
    capture.output(sessionInfo())
  ),
  file.path(output_dir, "RUN_LOG.txt"), useBytes = TRUE
)

message("Four-time-point osteocalcin analysis written to: ", output_dir)
