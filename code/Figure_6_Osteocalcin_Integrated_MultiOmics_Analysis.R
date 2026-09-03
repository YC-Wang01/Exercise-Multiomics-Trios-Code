# Integrated osteocalcin audit for Figure 6 planning.
# New outputs only: existing Figure 6, Figure S6 and editable artwork are untouched.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  if (!is.na(.local_script)) file.path(dirname(.local_script), "Fig00_Config.R") else NA_character_
))
.config_file <- .config_candidates[file.exists(.config_candidates)][1]
if (is.na(.config_file)) stop("Cannot locate Fig00_Config.R.")
source(.config_file)
rm(.local_script, .config_candidates, .config_file)

p_load(dplyr, emmeans, ggplot2, lme4, lmerTest, patchwork, readr,
       scales, stringr, tibble, tidyr, writexl)

analysis_id <- "Fig06_OC_Integrated_Analysis_2026-08-22"
output_dir <- file.path(PROJECT_ROOT, "fig06", "Temp_OC_Integrated_Analysis_2026-08-22")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sentinel <- file.path(output_dir, "Fig06_OC_Integrated_Results.xlsx")
if (file.exists(sentinel)) {
  message("Replacing validated derived output: ", sentinel)
}

role_levels <- c("Daughter", "Mother", "Father")
marker_order <- c("cOC", "TotalOC", "cOC_ratio")
marker_labels <- c(
  TotalOC = "Total osteocalcin",
  cOC = "Carboxylated osteocalcin",
  cOC_ratio = "cOC/tOC ratio"
)
marker_variables <- c(
  TotalOC = "S_TotalOC_84m_G",
  cOC = "S_cOC_84m_G",
  cOC_ratio = "cOC_totalOC_ratio_84m_G"
)

theme_atm <- function(base_size = 8) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(
      axis.text = element_text(color = "#17242D"),
      axis.title = element_text(color = "#17242D", face = "bold"),
      plot.title = element_text(face = "bold", size = 9, color = "#17242D"),
      strip.text = element_text(face = "bold", size = 7.5, color = "#17242D"),
      strip.background = element_rect(fill = "white", color = "#17242D", linewidth = 0.35),
      panel.border = element_rect(fill = NA, color = "#17242D", linewidth = 0.4),
      axis.line = element_blank(),
      legend.title = element_text(face = "bold", size = 7),
      legend.text = element_text(size = 6.5)
    )
}

rank_normalize <- function(x) {
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  n <- sum(keep)
  if (n >= 3L) out[keep] <- qnorm((rank(x[keep], ties.method = "average") - 0.5) / n)
  out
}

safe_log2 <- function(x) ifelse(is.finite(x) & x > 0, log2(x), NA_real_)

significance <- function(x) {
  case_when(
    is.na(x) ~ "",
    x < 0.001 ~ "***",
    x < 0.01 ~ "**",
    x < 0.05 ~ "*",
    TRUE ~ ""
  )
}

clean_label <- function(x) {
  x %>%
    str_replace_all("_84m_G$|_2010_G$|_G$", "") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\bS TotalOC\\b", "Total osteocalcin") %>%
    str_replace_all("\\bS cOC\\b", "Carboxylated osteocalcin") %>%
    str_replace_all("cOC totalOC ratio", "cOC/tOC ratio") %>%
    str_squish()
}

# -------------------------------------------------------------------------
# 1. Expanded-cohort clinical associations and role interaction tests
# -------------------------------------------------------------------------
phen <- read_csv(
  path_intermediate("Clinical_Observed_SourcePool.csv"),
  show_col_types = FALSE, progress = FALSE
) %>%
  mutate(
    Family_UID = factor(as.character(FamilyID)),
    Role = factor(Role, levels = role_levels)
  )

audit <- read_csv(
  file.path(path_results("Clinical_Association_Atlas_179"),
            "Clinical_Association_Atlas_VariableAudit.csv"),
  show_col_types = FALSE, progress = FALSE
)
registry <- read_csv(
  path_analysis_ready("Metadata", "Clinical_Variable_Registry.csv"),
  show_col_types = FALSE, progress = FALSE
)
label_lookup <- setNames(registry$VariableLabel, registry$SourceVariable)
domain_lookup <- setNames(registry$Domain, registry$SourceVariable)

predictors <- audit %>%
  filter(EligibleContinuous) %>%
  pull(Variable) %>%
  setdiff(unname(marker_variables))

extract_term <- function(fit, term = "x_z") {
  tab <- summary(fit)$coefficients
  if (!term %in% rownames(tab)) return(rep(NA_real_, 6))
  beta <- tab[term, "Estimate"]
  se <- tab[term, "Std. Error"]
  df <- tab[term, "df"]
  p <- tab[term, "Pr(>|t|)"]
  critical <- qt(0.975, df = df)
  c(beta = beta, se = se, p = p,
    ci_low = beta - critical * se, ci_high = beta + critical * se,
    df = df)
}

fit_clinical_pair <- function(marker, predictor) {
  outcome <- unname(marker_variables[marker])
  dat <- phen[, c("Family_UID", "Role", predictor, outcome)]
  names(dat)[3:4] <- c("x", "y")
  dat <- dat[complete.cases(dat), , drop = FALSE]
  dat$Family_UID <- droplevels(dat$Family_UID)
  dat$Role <- droplevels(dat$Role)
  n <- nrow(dat)
  families <- nlevels(dat$Family_UID)
  if (n < 40L || families < 15L || nlevels(dat$Role) < 3L ||
      length(unique(dat$x)) < 5L || length(unique(dat$y)) < 5L) {
    return(tibble(
      Marker = marker, MarkerLabel = unname(marker_labels[marker]),
      Predictor = predictor, PredictorLabel = clean_label(ifelse(
        is.na(label_lookup[predictor]) | !nzchar(label_lookup[predictor]),
        predictor, label_lookup[predictor]
      )), PredictorDomain = unname(domain_lookup[predictor]),
      N = n, Families = families, Effect = NA_real_, SE = NA_real_,
      CI_Low = NA_real_, CI_High = NA_real_, P_Value = NA_real_,
      RoleInteractionP = NA_real_, SingularFit = NA, InteractionSingularFit = NA,
      Status = "Insufficient data"
    ))
  }
  dat$x_z <- rank_normalize(dat$x)
  dat$y_z <- rank_normalize(dat$y)
  base <- tryCatch(
    suppressMessages(suppressWarnings(lmer(
      y_z ~ x_z + Role + (1 | Family_UID), data = dat, REML = FALSE
    ))),
    error = function(e) e
  )
  full <- tryCatch(
    suppressMessages(suppressWarnings(lmer(
      y_z ~ x_z * Role + (1 | Family_UID), data = dat, REML = FALSE
    ))),
    error = function(e) e
  )
  if (inherits(base, "error")) {
    values <- rep(NA_real_, 6)
    base_singular <- NA
    status <- conditionMessage(base)
  } else {
    values <- extract_term(base)
    base_singular <- isSingular(base, tol = 1e-4)
    status <- "OK"
  }
  interaction_p <- NA_real_
  full_singular <- NA
  if (!inherits(base, "error") && !inherits(full, "error")) {
    comparison <- suppressMessages(anova(base, full))
    interaction_p <- comparison$`Pr(>Chisq)`[2]
    full_singular <- isSingular(full, tol = 1e-4)
  }
  tibble(
    Marker = marker,
    MarkerLabel = unname(marker_labels[marker]),
    Predictor = predictor,
    PredictorLabel = clean_label(ifelse(
      is.na(label_lookup[predictor]) | !nzchar(label_lookup[predictor]),
      predictor, label_lookup[predictor]
    )),
    PredictorDomain = unname(domain_lookup[predictor]),
    N = n,
    Families = families,
    Effect = unname(values["beta"]),
    SE = unname(values["se"]),
    CI_Low = unname(values["ci_low"]),
    CI_High = unname(values["ci_high"]),
    P_Value = unname(values["p"]),
    RoleInteractionP = interaction_p,
    SingularFit = base_singular,
    InteractionSingularFit = full_singular,
    Status = status
  )
}

clinical_associations <- bind_rows(lapply(marker_order, function(marker) {
  bind_rows(lapply(predictors, function(predictor) fit_clinical_pair(marker, predictor)))
})) %>%
  mutate(
    BH_FDR_Global = p.adjust(P_Value, method = "BH"),
    RoleInteraction_BH_FDR = p.adjust(RoleInteractionP, method = "BH"),
    Significance = significance(BH_FDR_Global),
    RoleInteractionSignificance = significance(RoleInteraction_BH_FDR)
  ) %>%
  arrange(BH_FDR_Global, P_Value)

# -------------------------------------------------------------------------
# 2. Direct acute changes and time-by-role models
# -------------------------------------------------------------------------
clinical_long <- read_csv(
  file.path(path_intermediate("Clinical_Exercise_Source_Audit_2026-08-21"),
            "Restricted", "OC_BCAA_Serotonin_Long_Internal.csv"),
  show_col_types = FALSE, progress = FALSE
) %>%
  filter(
    MatchStatus == "Unique locked-cohort match",
    Marker %in% marker_order,
    Timepoint %in% c("Pre", "Post1h", "Post3h")
  ) %>%
  mutate(
    FamilyID = factor(as.character(FamilyID)),
    Clinical_Subject_ID = factor(as.character(Clinical_Subject_ID)),
    Role = factor(Role, levels = role_levels),
    Timepoint = factor(Timepoint, levels = c("Pre", "Post1h", "Post3h")),
    Value = suppressWarnings(as.numeric(Value)),
    TransformedValue = safe_log2(Value)
  ) %>%
  filter(is.finite(TransformedValue))

time_model_results <- list()
time_emmeans <- list()
role_contrasts <- list()

for (marker in marker_order) {
  dat <- clinical_long %>% filter(Marker == marker) %>% droplevels()
  complete_subjects <- dat %>%
    count(Clinical_Subject_ID, name = "Timepoints") %>%
    filter(Timepoints == 3L) %>%
    pull(Clinical_Subject_ID)
  dat <- dat %>% filter(Clinical_Subject_ID %in% complete_subjects) %>% droplevels()

  additive <- suppressMessages(suppressWarnings(lmer(
    TransformedValue ~ Timepoint + Role + (1 | FamilyID) + (1 | Clinical_Subject_ID),
    data = dat, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa", check.conv.singular = "ignore")
  )))
  interaction <- suppressMessages(suppressWarnings(lmer(
    TransformedValue ~ Timepoint * Role + (1 | FamilyID) + (1 | Clinical_Subject_ID),
    data = dat, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa", check.conv.singular = "ignore")
  )))
  comparison <- suppressMessages(anova(additive, interaction))
  interaction_p <- comparison$`Pr(>Chisq)`[2]

  emm <- emmeans(interaction, ~ Timepoint * Role)
  emm_table <- as.data.frame(emm) %>%
    as_tibble() %>%
    transmute(
      Marker = marker,
      MarkerLabel = unname(marker_labels[marker]),
      Timepoint,
      Role,
      Estimate_log2 = emmean,
      SE,
      DF = df,
      Lower_log2 = lower.CL,
      Upper_log2 = upper.CL,
      GeometricMean = 2^emmean,
      Lower95 = 2^lower.CL,
      Upper95 = 2^upper.CL,
      Subjects = n_distinct(dat$Clinical_Subject_ID),
      Families = n_distinct(dat$FamilyID)
    )
  time_emmeans[[marker]] <- emm_table

  role_pairs <- as.data.frame(pairs(emmeans(interaction, ~ Role | Timepoint), adjust = "none")) %>%
    as_tibble() %>%
    filter(contrast %in% c("Daughter - Mother", "Mother - Father")) %>%
    transmute(
      Marker = marker,
      MarkerLabel = unname(marker_labels[marker]),
      Timepoint,
      Contrast = contrast,
      Effect_log2 = estimate,
      FoldDifference = 2^estimate,
      SE,
      DF = df,
      Lower95_log2 = estimate - qt(0.975, df) * SE,
      Upper95_log2 = estimate + qt(0.975, df) * SE,
      P_Value = p.value,
      Subjects = n_distinct(dat$Clinical_Subject_ID),
      Families = n_distinct(dat$FamilyID)
    )
  role_contrasts[[marker]] <- role_pairs
  time_model_results[[marker]] <- tibble(
    Marker = marker,
    MarkerLabel = unname(marker_labels[marker]),
    Subjects = n_distinct(dat$Clinical_Subject_ID),
    Families = n_distinct(dat$FamilyID),
    CompleteTrios = dat %>% distinct(FamilyID, Role) %>% count(FamilyID) %>%
      summarise(n = sum(n == 3L)) %>% pull(n),
    TimeByRole_LRT_P = interaction_p,
    AdditiveSingularFit = isSingular(additive, tol = 1e-4),
    InteractionSingularFit = isSingular(interaction, tol = 1e-4)
  )
}

time_model_results <- bind_rows(time_model_results) %>%
  mutate(TimeByRole_BH_FDR = p.adjust(TimeByRole_LRT_P, method = "BH"))
time_emmeans <- bind_rows(time_emmeans)
role_contrasts <- bind_rows(role_contrasts) %>%
  mutate(BH_FDR_Global = p.adjust(P_Value, method = "BH"),
         Significance = significance(BH_FDR_Global))

direct_effects <- read_csv(
  file.path(path_results("Six_Marker_Family_State_Analysis_2026-08-22"),
            "Direct_Time_Effects.csv"),
  show_col_types = FALSE, progress = FALSE
) %>%
  filter(Marker %in% marker_order)

# -------------------------------------------------------------------------
# 3. Family-level state and response similarity
# -------------------------------------------------------------------------
family_dir <- path_results("Six_Marker_Family_State_Analysis_2026-08-22")
family_icc <- read_csv(file.path(family_dir, "Family_ICC_Pre_Post_Delta.csv"),
                       show_col_types = FALSE, progress = FALSE) %>%
  filter(Marker %in% marker_order)
family_pairs <- read_csv(file.path(family_dir, "DM_DF_MF_Pre_Post_Delta.csv"),
                         show_col_types = FALSE, progress = FALSE) %>%
  filter(Marker %in% marker_order)

# -------------------------------------------------------------------------
# 4. Omics models: Pre-pre, Pre-response, Post-post and Delta-delta
# -------------------------------------------------------------------------
three_state <- readRDS(file.path(
  path_results("Six_Markers_Three_State_Omics_2026-08-22"),
  "Six_Markers_Three_State_Omics.rds"
))

three_feature_summary <- three_state$feature_summary %>%
  filter(ClinicalMarker %in% marker_order)
three_feature_hits <- three_state$feature_results %>%
  filter(ClinicalMarker %in% marker_order, BH_FDR < 0.05)
three_pathway_summary <- three_state$pathway_summary %>%
  filter(ClinicalMarker %in% marker_order)
three_pathway_hits <- three_state$pathway_results %>%
  filter(ClinicalMarker %in% marker_order, BH_FDR < 0.05)

tissue_pre_response <- readRDS(file.path(
  PROJECT_ROOT, "fig06", "analysis",
  "Fig06_OC_GO_Hallmark_PreExercise_Exercise_2026-08-18",
  "Fig06_OC_GO_Hallmark_PreExercise_Exercise.rds"
))

tissue_feature_summary <- tissue_pre_response$FeatureSummary %>%
  filter(State == "ExerciseResponse", Marker %in% marker_order) %>%
  transmute(
    Dataset, Layer = "Proteomics", Tissue,
    Model = "Pre_Response", Timepoint = "Post3h",
    ClinicalMarker = Marker, ClinicalMarkerLabel = MarkerLabel,
    Subjects, Families, TestedFeatures, NominalP05, BH_FDR05, BH_FDR10,
    BH_FDR20, MinimumP, MinimumBH_FDR, FamilyCorrelation
  )
tissue_feature_hits <- tissue_pre_response$FeatureResults %>%
  filter(State == "ExerciseResponse", Marker %in% marker_order, BH_FDR < 0.05) %>%
  transmute(
    AnalysisID, Dataset, Layer = "Proteomics", Tissue,
    Model = "Pre_Response", Timepoint = "Post3h",
    ClinicalMarker = Marker, ClinicalMarkerLabel = MarkerLabel,
    ClinicalVariable = "Pre_log2", FeatureID, Effect, ModeratedT,
    P_Value, BH_FDR, ObservedN, Subjects, Families, FamilyCorrelation,
    GeneSymbol = NA_character_
  )
tissue_pathway_summary <- tissue_pre_response$PathwaySummary %>%
  filter(Database == "Hallmark", State == "ExerciseResponse", Marker %in% marker_order) %>%
  transmute(
    Dataset, Layer = "Proteomics", Tissue,
    Model = "Pre_Response", Timepoint = "Post3h",
    ClinicalMarker = Marker, ClinicalMarkerLabel = MarkerLabel,
    Subjects = NA_integer_, Families = NA_integer_, TestedPathways, NominalP05, BH_FDR05, BH_FDR10,
    BH_FDR20, MinimumP, MinimumBH_FDR
  )
tissue_pathway_hits <- tissue_pre_response$PathwayResults %>%
  filter(Database == "Hallmark", State == "ExerciseResponse",
         Marker %in% marker_order, BH_FDR < 0.05) %>%
  transmute(
    AnalysisID, Dataset, Layer = "Proteomics", Tissue,
    Model = "Pre_Response", Timepoint = "Post3h",
    ClinicalMarker = Marker, ClinicalMarkerLabel = MarkerLabel,
    Subjects = NA_integer_, Families = NA_integer_, Database = "MSigDB Hallmark",
    PathwayID, Pathway, Size = NGenes, NES, P_Value, BH_FDR, LeadingEdge
  )

serum_pre_response <- readRDS(file.path(output_dir, "PreResponse_NonMethylation_Results.rds"))
methyl_pre_response <- readRDS(file.path(output_dir, "PreResponse_Methylation_Results.rds"))

omics_feature_summary <- bind_rows(
  three_feature_summary,
  tissue_feature_summary,
  serum_pre_response$feature_summary,
  methyl_pre_response$feature_summary
) %>%
  distinct(Dataset, Model, Timepoint, ClinicalMarker, .keep_all = TRUE) %>%
  mutate(
    ClinicalMarker = factor(ClinicalMarker, levels = marker_order),
    Model = factor(Model, levels = c("Pre_Pre", "Pre_Response", "Post_Post", "Delta_Delta"))
  ) %>%
  arrange(ClinicalMarker, Model, Tissue, Layer, Timepoint)

omics_feature_hits <- bind_rows(
  three_feature_hits,
  tissue_feature_hits,
  serum_pre_response$feature_results %>% filter(BH_FDR < 0.05),
  methyl_pre_response$feature_results %>% filter(BH_FDR < 0.05)
) %>%
  distinct(Dataset, Model, Timepoint, ClinicalMarker, FeatureID, .keep_all = TRUE) %>%
  arrange(Model, ClinicalMarker, BH_FDR, P_Value)

omics_pathway_summary <- bind_rows(
  three_pathway_summary,
  tissue_pathway_summary,
  serum_pre_response$pathway_summary,
  methyl_pre_response$pathway_summary
) %>%
  distinct(Dataset, Model, Timepoint, ClinicalMarker, .keep_all = TRUE) %>%
  mutate(
    ClinicalMarker = factor(ClinicalMarker, levels = marker_order),
    Model = factor(Model, levels = c("Pre_Pre", "Pre_Response", "Post_Post", "Delta_Delta"))
  ) %>%
  arrange(ClinicalMarker, Model, Tissue, Layer, Timepoint)

omics_pathway_hits <- bind_rows(
  three_pathway_hits,
  tissue_pathway_hits,
  serum_pre_response$pathway_results %>% filter(BH_FDR < 0.05),
  methyl_pre_response$pathway_results %>% filter(BH_FDR < 0.05)
) %>%
  distinct(Dataset, Model, Timepoint, ClinicalMarker, PathwayID, .keep_all = TRUE) %>%
  arrange(Model, ClinicalMarker, BH_FDR, P_Value)

# -------------------------------------------------------------------------
# 5. Reproducible working figures
# -------------------------------------------------------------------------
clinical_plot_vars <- clinical_associations %>%
  filter(BH_FDR_Global < 0.05, !is.na(PredictorLabel)) %>%
  group_by(Predictor, PredictorLabel, PredictorDomain) %>%
  summarise(MinFDR = min(BH_FDR_Global), MaxAbsEffect = max(abs(Effect)), .groups = "drop") %>%
  arrange(MinFDR, desc(MaxAbsEffect)) %>%
  slice_head(n = 24)

clinical_plot_data <- clinical_associations %>%
  semi_join(clinical_plot_vars, by = c("Predictor", "PredictorLabel", "PredictorDomain")) %>%
  mutate(
    MarkerLabel = factor(MarkerLabel, levels = unname(marker_labels[marker_order])),
    PredictorLabel = factor(PredictorLabel, levels = rev(clinical_plot_vars$PredictorLabel)),
    Label = significance(BH_FDR_Global)
  )

p_clinical_heatmap <- ggplot(clinical_plot_data,
                             aes(x = MarkerLabel, y = PredictorLabel, fill = Effect)) +
  geom_tile(color = "#17242D", linewidth = 0.25) +
  geom_text(aes(label = Label), family = "Arial", fontface = "bold", size = 2.8) +
  scale_fill_gradient2(low = "#6D95B6", mid = "white", high = "#D7796D",
                       midpoint = 0, limits = c(-1, 1), oob = squish) +
  labs(
    x = NULL, y = NULL, fill = "Standardized\neffect",
    title = "Clinical correlates of osteocalcin measures"
  ) +
  theme_atm() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

ggsave(file.path(output_dir, "Fig06_Temp_Clinical_Associations.pdf"),
       p_clinical_heatmap, width = 150, height = 145, units = "mm", device = cairo_pdf)

p_time <- time_emmeans %>%
  mutate(
    MarkerLabel = factor(MarkerLabel, levels = unname(marker_labels[marker_order])),
    Timepoint = factor(Timepoint, levels = c("Pre", "Post1h", "Post3h"),
                       labels = c("Pre-exercise", "Post 1 h", "Post 3 h"))
  ) %>%
  ggplot(aes(x = Timepoint, y = GeometricMean, color = Role, group = Role)) +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0, linewidth = 0.55) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.0) +
  facet_wrap(~MarkerLabel, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c(Daughter = "#E58B73", Mother = "#4AA4A7", Father = "#527EAD")) +
  labs(x = NULL, y = "Model-estimated geometric mean", color = NULL,
       title = "Osteocalcin trajectories by family role") +
  theme_atm() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")

ggsave(file.path(output_dir, "Fig06_Temp_OC_Time_by_Role.pdf"),
       p_time, width = 190, height = 72, units = "mm", device = cairo_pdf)

p_role <- role_contrasts %>%
  mutate(
    MarkerLabel = factor(MarkerLabel, levels = rev(unname(marker_labels[marker_order]))),
    Timepoint = factor(Timepoint, levels = c("Pre", "Post1h", "Post3h"),
                       labels = c("Pre-exercise", "Post 1 h", "Post 3 h")),
    Contrast = recode(Contrast,
                      "Daughter - Mother" = "Daughter vs mother",
                      "Mother - Father" = "Mother vs father")
  ) %>%
  ggplot(aes(x = Effect_log2, y = MarkerLabel, color = Timepoint)) +
  geom_vline(xintercept = 0, color = "#66757F", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = Lower95_log2, xmax = Upper95_log2), height = 0,
                 position = position_dodge(width = 0.55), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.55), size = 1.8) +
  geom_text(aes(label = Significance), position = position_dodge(width = 0.55),
            hjust = -0.7, color = "#17242D", fontface = "bold", size = 2.5) +
  facet_wrap(~Contrast, nrow = 1) +
  scale_color_manual(values = c("Pre-exercise" = "#7E8D98", "Post 1 h" = "#D78361",
                                "Post 3 h" = "#4B8CA8")) +
  labs(x = "Role contrast in log2 osteocalcin level", y = NULL, color = NULL,
       title = "Family-role contrasts across sampling times") +
  theme_atm() +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "Fig06_Temp_OC_Role_Contrasts.pdf"),
       p_role, width = 180, height = 82, units = "mm", device = cairo_pdf)

state_levels <- c("Pre-exercise", "Post 1 h", "Post 3 h", "Post 1 h response", "Post 3 h response")
icc_plot <- family_icc %>%
  mutate(
    MarkerLabel = factor(MarkerLabel, levels = rev(unname(marker_labels[marker_order]))),
    StateLabel = factor(StateLabel, levels = state_levels),
    Label = ifelse(is.finite(FamilyICC),
                   sprintf("%.2f%s", FamilyICC, significance(PermutationBH_FDR)), "-")
  ) %>%
  ggplot(aes(x = StateLabel, y = MarkerLabel, fill = FamilyICC)) +
  geom_tile(color = "#17242D", linewidth = 0.25) +
  geom_text(aes(label = Label), family = "Arial", size = 2.4) +
  scale_fill_gradient(low = "white", high = "#4B9A8E", limits = c(0, 0.7),
                      na.value = "#EDF0F1", oob = squish) +
  labs(x = NULL, y = NULL, fill = "Family ICC", title = "Role-adjusted family clustering") +
  theme_atm() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

pair_plot <- family_pairs %>%
  mutate(
    MarkerPair = paste(MarkerLabel, Pair, sep = " | "),
    MarkerPair = factor(MarkerPair, levels = rev(unique(MarkerPair))),
    StateLabel = factor(StateLabel, levels = state_levels),
    Label = ifelse(is.finite(SpearmanRho),
                   sprintf("%.2f%s", SpearmanRho, significance(PermutationBH_FDR)), "-")
  ) %>%
  ggplot(aes(x = StateLabel, y = MarkerPair, fill = SpearmanRho)) +
  geom_tile(color = "#17242D", linewidth = 0.22) +
  geom_text(aes(label = Label), family = "Arial", size = 2.1) +
  scale_fill_gradient2(low = "#5E86AD", mid = "white", high = "#D7796D",
                       midpoint = 0, limits = c(-0.9, 0.9), na.value = "#EDF0F1", oob = squish) +
  labs(x = NULL, y = NULL, fill = "Spearman rho",
       title = "Within-family pair similarity") +
  theme_atm() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 6.3),
        legend.position = "right")

family_combined <- icc_plot / pair_plot + plot_layout(heights = c(0.7, 1.3))
ggsave(file.path(output_dir, "Fig06_Temp_OC_Family_State_and_Response.pdf"),
       family_combined, width = 190, height = 175, units = "mm", device = cairo_pdf)

dataset_labels <- c(
  Serum_Proteomics = "Serum proteomics",
  Serum_Metabonomics = "Serum metabolomics",
  Adipose_Microarray = "Adipose transcriptomics",
  Adipose_Proteomics = "Adipose proteomics",
  Adipose_Methylation = "Adipose DNA methylation",
  Muscle_Microarray = "Muscle transcriptomics",
  Muscle_Proteomics = "Muscle proteomics",
  Muscle_Methylation = "Muscle DNA methylation"
)
model_labels <- c(
  Pre_Pre = "Pre-exercise level",
  Pre_Response = "Pre-exercise predictor of response",
  Post_Post = "Post-exercise level",
  Delta_Delta = "Change-to-change association"
)

feature_plot_data <- omics_feature_summary %>%
  mutate(
    DatasetLabel = recode(Dataset, !!!dataset_labels),
    ModelLabel = recode(as.character(Model), !!!model_labels),
    ModelTime = ifelse(Model %in% c("Post_Post", "Delta_Delta") & Timepoint == "Post1h",
                       paste0(ModelLabel, " (1 h)"), ModelLabel),
    ModelTime = ifelse(Model %in% c("Post_Post", "Delta_Delta", "Pre_Response") & Timepoint == "Post3h",
                       paste0(ModelLabel, " (3 h)"), ModelTime),
    MarkerLabel = factor(ClinicalMarkerLabel, levels = unname(marker_labels[marker_order])),
    FillValue = log10(BH_FDR05 + 1),
    Label = as.character(BH_FDR05)
  )

p_feature_counts <- ggplot(feature_plot_data,
                           aes(x = DatasetLabel, y = ModelTime, fill = FillValue)) +
  geom_tile(color = "#17242D", linewidth = 0.2) +
  geom_text(aes(label = Label), size = 2.1, family = "Arial") +
  facet_wrap(~MarkerLabel, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#C65F73", name = "log10(FDR hits + 1)") +
  labs(x = NULL, y = NULL, title = "Feature-level BH-FDR discoveries across osteocalcin models") +
  theme_atm() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6.3), legend.position = "right")

ggsave(file.path(output_dir, "Fig06_Temp_Omics_Feature_FDR_Overview.pdf"),
       p_feature_counts, width = 190, height = 105, units = "mm", device = cairo_pdf)

pathway_plot_data <- omics_pathway_summary %>%
  mutate(
    DatasetLabel = recode(Dataset, !!!dataset_labels),
    ModelLabel = recode(as.character(Model), !!!model_labels),
    ModelTime = ifelse(Model %in% c("Post_Post", "Delta_Delta") & Timepoint == "Post1h",
                       paste0(ModelLabel, " (1 h)"), ModelLabel),
    ModelTime = ifelse(Model %in% c("Post_Post", "Delta_Delta", "Pre_Response") & Timepoint == "Post3h",
                       paste0(ModelLabel, " (3 h)"), ModelTime),
    MarkerLabel = factor(ClinicalMarkerLabel, levels = unname(marker_labels[marker_order])),
    FillValue = log10(BH_FDR05 + 1),
    Label = as.character(BH_FDR05)
  )

p_pathway_counts <- ggplot(pathway_plot_data,
                           aes(x = DatasetLabel, y = ModelTime, fill = FillValue)) +
  geom_tile(color = "#17242D", linewidth = 0.2) +
  geom_text(aes(label = Label), size = 2.1, family = "Arial") +
  facet_wrap(~MarkerLabel, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#4C8FAE", name = "log10(FDR pathways + 1)") +
  labs(x = NULL, y = NULL, title = "Hallmark BH-FDR discoveries across osteocalcin models") +
  theme_atm() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6.3), legend.position = "right")

ggsave(file.path(output_dir, "Fig06_Temp_Hallmark_FDR_Overview.pdf"),
       p_pathway_counts, width = 190, height = 105, units = "mm", device = cairo_pdf)

directional_pathways <- omics_pathway_hits %>%
  filter(is.finite(NES)) %>%
  group_by(Pathway) %>%
  mutate(PathwayMinFDR = min(BH_FDR, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(rank(PathwayMinFDR, ties.method = "min") <= 18) %>%
  mutate(
    Pathway = factor(Pathway, levels = rev(unique(Pathway[order(PathwayMinFDR)]))),
    ModelLabel = recode(as.character(Model), !!!model_labels),
    ModelTime = ifelse(Timepoint == "Post1h", paste0(ModelLabel, " (1 h)"), ModelLabel),
    ModelTime = ifelse(Timepoint == "Post3h" & Model %in% c("Post_Post", "Delta_Delta", "Pre_Response"),
                       paste0(ModelLabel, " (3 h)"), ModelTime),
    TissueLayer = paste(Tissue, Layer),
    MarkerLabel = factor(ClinicalMarkerLabel, levels = unname(marker_labels[marker_order]))
  )

p_pathways <- ggplot(directional_pathways,
                     aes(x = interaction(ModelTime, TissueLayer, sep = " | "), y = Pathway,
                         color = NES, size = -log10(pmax(BH_FDR, 1e-50)))) +
  geom_point(alpha = 0.85) +
  facet_wrap(~MarkerLabel, nrow = 1, scales = "free_x") +
  scale_color_gradient2(low = "#5E86AD", mid = "white", high = "#D7796D", midpoint = 0) +
  scale_size_continuous(range = c(1.2, 4.5), name = expression(-log[10](BH-FDR))) +
  labs(x = NULL, y = NULL, color = "NES", title = "Top directional Hallmark associations") +
  theme_atm() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.8),
        axis.text.y = element_text(size = 6.2), legend.position = "right")

ggsave(file.path(output_dir, "Fig06_Temp_Top_Hallmark_Associations.pdf"),
       p_pathways, width = 190, height = 135, units = "mm", device = cairo_pdf)

# -------------------------------------------------------------------------
# 6. Workbook, provenance and manuscript-ready summary
# -------------------------------------------------------------------------
write_xlsx(
  list(
    Clinical_associations = clinical_associations,
    Clinical_FDR05 = clinical_associations %>% filter(BH_FDR_Global < 0.05),
    Clinical_role_interactions = clinical_associations %>%
      filter(is.finite(RoleInteractionP)) %>%
      arrange(RoleInteraction_BH_FDR, RoleInteractionP),
    Direct_time_effects = direct_effects,
    Time_by_role_tests = time_model_results,
    Time_by_role_estimates = time_emmeans,
    Role_contrasts = role_contrasts,
    Family_ICC = family_icc,
    Family_pair_similarity = family_pairs,
    Omics_feature_summary = omics_feature_summary,
    Omics_feature_FDR05 = omics_feature_hits,
    Hallmark_summary = omics_pathway_summary,
    Hallmark_FDR05 = omics_pathway_hits
  ),
  sentinel
)

write_csv(clinical_associations,
          file.path(output_dir, "Clinical_OC_Associations_and_Role_Interactions.csv"))
write_csv(time_model_results,
          file.path(output_dir, "OC_Time_by_Role_Tests.csv"))
write_csv(role_contrasts,
          file.path(output_dir, "OC_Role_Contrasts.csv"))
write_csv(omics_feature_summary,
          file.path(output_dir, "OC_Omics_Feature_FDR_Summary.csv"))
write_csv(omics_feature_hits,
          file.path(output_dir, "OC_Omics_Feature_FDR05.csv"))
write_csv(omics_pathway_summary,
          file.path(output_dir, "OC_Hallmark_FDR_Summary.csv"))
write_csv(omics_pathway_hits,
          file.path(output_dir, "OC_Hallmark_FDR05.csv"))

clinical_hits <- clinical_associations %>% filter(BH_FDR_Global < 0.05)
role_interaction_hits <- clinical_associations %>% filter(RoleInteraction_BH_FDR < 0.05)
time_interaction_hits <- time_model_results %>% filter(TimeByRole_BH_FDR < 0.05)
family_icc_hits <- family_icc %>% filter(PermutationBH_FDR < 0.05)
family_pair_hits <- family_pairs %>% filter(PermutationBH_FDR < 0.05)

summary_lines <- c(
  "# Integrated osteocalcin analysis for Figure 6", "",
  "## Analysis questions", "",
  "1. Which clinical phenotypes are associated with pre-exercise tOC, cOC and cOC/tOC ratio?",
  "2. Do acute changes differ by family role?",
  "3. Are pre-exercise, post-exercise or acute-response states similar within families?",
  "4. Which features and Hallmark pathways are associated in Pre-pre, Pre-response, Post-post and Delta-delta models?", "",
  "## Models", "",
  "- Clinical: rank-normalized osteocalcin outcome ~ rank-normalized phenotype + Role + (1 | FamilyID).",
  "- Clinical role interaction: the above model versus a phenotype-by-Role interaction model; BH correction across all tested interactions.",
  "- Time by role: log2(osteocalcin) ~ Timepoint * Role + (1 | FamilyID) + (1 | ParticipantID).",
  "- Omics: limma moderated models with Role and FamilyID duplicateCorrelation. Pre-response means omics delta ~ pre-exercise osteocalcin.",
  "- Pathways: full-ranked Hallmark GSEA on moderated t statistics; methylation uses CpG-count-adjusted unsigned methylglm.", "",
  "## Clinical associations", "",
  sprintf("The expanded clinical cohort yielded %d BH-FDR < 0.05 osteocalcin-phenotype associations across %d tested marker-predictor models.",
          nrow(clinical_hits), sum(is.finite(clinical_associations$P_Value))),
  sprintf("There were %d phenotype-by-role interactions at BH-FDR < 0.05.", nrow(role_interaction_hits)), "",
  "## Acute time and family structure", "",
  sprintf("There were %d osteocalcin time-by-role interactions at BH-FDR < 0.05.", nrow(time_interaction_hits)),
  sprintf("Family ICC analysis identified %d state-marker combinations at permutation BH-FDR < 0.05.", nrow(family_icc_hits)),
  sprintf("DM/DF/MF analysis identified %d state-marker-pair combinations at permutation BH-FDR < 0.05.", nrow(family_pair_hits)),
  "Stable similarity of a pre- or post-exercise level must not be described as familial similarity of the acute response.", "",
  "## Omics", "",
  sprintf("Across the four model classes, %d feature rows and %d Hallmark rows passed BH-FDR < 0.05.",
          nrow(omics_feature_hits), nrow(omics_pathway_hits)),
  "Positive GSEA NES indicates enrichment toward positive predictor coefficients; it does not prove pathway activation.",
  "Methylation pathway results are unsigned and provisional because assay provenance and the known design-version block remain unresolved.", "",
  "## Claim boundary", "",
  "These analyses test association and effect modification. They do not establish that exercise changed osteocalcin, that osteocalcin caused an omics response, or that an enriched pathway was activated or inhibited."
)
writeLines(summary_lines, file.path(output_dir, "RESULTS_SUMMARY.md"), useBytes = TRUE)

manifest <- tibble(
  AnalysisID = analysis_id,
  Generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  RVersion = R.version.string,
  OutputDirectory = output_dir,
  RawDataModified = FALSE,
  ExistingFigureFilesOverwritten = FALSE,
  PrimaryThreshold = "BH-FDR < 0.05",
  SuggestiveThreshold = "0.05 <= BH-FDR < 0.10"
)
write_csv(manifest, file.path(output_dir, "RUN_MANIFEST.csv"))
capture.output(sessionInfo(), file = file.path(output_dir, "SESSION_INFO.txt"))

message("Completed integrated osteocalcin audit: ", output_dir)
