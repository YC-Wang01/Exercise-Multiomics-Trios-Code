# Osteocalcin state and response associations with family, role, and sex code.
# Restricted and analysis-ready inputs are read-only. Outputs are aggregate.

.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.this_script <- if (length(.file_arg)) sub("^--file=", "", .file_arg[[1]]) else NA_character_
.script_dir <- if (!is.na(.this_script)) dirname(normalizePath(.this_script)) else file.path(getwd(), "code")
.config_file <- file.path(.script_dir, "Fig00_Config.R")
if (!file.exists(.config_file)) stop("Cannot locate Fig00_Config.R.")
source(.config_file)
rm(.file_arg, .this_script, .config_file)

p_load(dplyr, emmeans, lme4, lmerTest, readr, tibble, tidyr, writexl)

analysis_id <- "Fig06_OC_State_Family_Sex_Role_Audit_2026-08-30"
output_dir <- path_project("fig06", "analysis", analysis_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

longitudinal_file <- file.path(
  path_intermediate("Clinical_Exercise_Source_Audit_2026-08-21"),
  "Restricted", "OC_BCAA_Serotonin_Long_Internal.csv"
)
phenotype_file <- path_analysis_ready("Metadata", "Clinical_Observed_Primary.csv")
family_icc_file <- path_project(
  "Pre_Post_Response_Analysis", "Results", "08_Osteocalcin_Molecular_State_and_Response",
  "03_Family_Structure", "OC_Family_ICC_Pre_Post_Delta.csv"
)
pair_file <- path_project(
  "Pre_Post_Response_Analysis", "Results", "08_Osteocalcin_Molecular_State_and_Response",
  "03_Family_Structure", "OC_DM_DF_MF_Pre_Post_Delta.csv"
)
for (required_file in c(longitudinal_file, phenotype_file, family_icc_file, pair_file)) {
  if (!file.exists(required_file)) stop("Required input is missing: ", required_file)
}

marker_registry <- tribble(
  ~Marker, ~MarkerLabel,
  "TotalOC", "Total osteocalcin (tOC)",
  "cOC", "Carboxylated osteocalcin (cOC)",
  "cOC_ratio", "cOC/tOC ratio"
)
role_levels <- c("Daughter", "Mother", "Father")

phenotypes <- read_csv(phenotype_file, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    GenderCode = as.character(Gender),
    Sex = factor(
      case_when(GenderCode == "1" ~ "Female", GenderCode == "2" ~ "Male", TRUE ~ NA_character_),
      levels = c("Female", "Male")
    )
  )

longitudinal <- read_csv(longitudinal_file, show_col_types = FALSE, progress = FALSE) %>%
  filter(
    MatchStatus == "Unique locked-cohort match",
    Marker %in% marker_registry$Marker,
    Timepoint %in% c("Pre", "Post1h", "Post3h")
  ) %>%
  transmute(
    FamilyID = factor(as.character(FamilyID)),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = factor(Role, levels = role_levels),
    Marker = as.character(Marker),
    Timepoint = as.character(Timepoint),
    Value = suppressWarnings(as.numeric(Value))
  ) %>%
  filter(is.finite(Value), Value > 0) %>%
  distinct(FamilyID, Clinical_Subject_ID, Role, Marker, Timepoint, .keep_all = TRUE)

complete_ids <- longitudinal %>%
  count(Clinical_Subject_ID, Marker, name = "ObservedTimepoints") %>%
  filter(ObservedTimepoints == 3L) %>%
  count(Clinical_Subject_ID, name = "CompleteMarkers") %>%
  filter(CompleteMarkers == 3L) %>%
  pull(Clinical_Subject_ID)

state_data <- longitudinal %>%
  filter(Clinical_Subject_ID %in% complete_ids) %>%
  select(FamilyID, Clinical_Subject_ID, Role, Marker, Timepoint, Value) %>%
  pivot_wider(names_from = Timepoint, values_from = Value) %>%
  filter(if_all(c(Pre, Post1h, Post3h), ~ is.finite(.x) & .x > 0)) %>%
  mutate(
    Pre = log2(Pre),
    Post1h = log2(Post1h),
    Post3h = log2(Post3h),
    Delta1h = Post1h - Pre,
    Delta3h = Post3h - Pre
  ) %>%
  left_join(phenotypes, by = "Clinical_Subject_ID") %>%
  left_join(marker_registry, by = "Marker") %>%
  droplevels()

if (n_distinct(state_data$Clinical_Subject_ID) != 50L) stop("Expected 50 complete participants.")
if (any(is.na(state_data$Sex))) stop("Sex code is missing in the complete cohort.")

state_registry <- tribble(
  ~State, ~StateLabel,
  "Pre", "Pre-exercise",
  "Post1h", "Post 1 h",
  "Post3h", "Post 3 h",
  "Delta1h", "Post 1 h - Pre",
  "Delta3h", "Post 3 h - Pre"
)

standardize <- function(x) as.numeric(scale(x))

fit_fixed_effect <- function(data, predictor, state_name) {
  dat <- data %>%
    transmute(
      FamilyID, Role, Sex,
      OutcomeZ = standardize(.data[[state_name]])
    ) %>%
    filter(is.finite(OutcomeZ), !is.na(FamilyID), !is.na(Role), !is.na(Sex)) %>%
    droplevels()

  formula_full <- if (predictor == "Role") {
    OutcomeZ ~ Role + (1 | FamilyID)
  } else {
    OutcomeZ ~ Sex + (1 | FamilyID)
  }
  formula_null <- OutcomeZ ~ 1 + (1 | FamilyID)

  full_fit <- suppressMessages(suppressWarnings(lmer(
    formula_full, data = dat, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa", check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))
  )))
  null_fit <- suppressMessages(suppressWarnings(lmer(
    formula_null, data = dat, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa", check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))
  )))
  lrt <- anova(null_fit, full_fit)

  coefficient_table <- as.data.frame(summary(full_fit)$coefficients) %>%
    rownames_to_column("Term") %>%
    as_tibble() %>%
    filter(Term != "(Intercept)") %>%
    transmute(
      Term,
      StandardizedEffect = Estimate,
      SE = `Std. Error`,
      CI_Lower = Estimate - 1.96 * `Std. Error`,
      CI_Upper = Estimate + 1.96 * `Std. Error`,
      CoefficientP = `Pr(>|t|)`
    )

  omnibus <- tibble(
    Predictor = predictor,
    State = state_name,
    Subjects = nrow(dat),
    Families = n_distinct(dat$FamilyID),
    OmnibusChiSq = lrt$Chisq[[2]],
    OmnibusDF = lrt$`Chi Df`[[2]],
    OmnibusP = lrt$`Pr(>Chisq)`[[2]],
    SingularFit = isSingular(full_fit, tol = 1e-4)
  )

  pairwise <- if (predictor == "Role") {
    as.data.frame(emmeans::contrast(emmeans::emmeans(full_fit, ~ Role), method = "pairwise", adjust = "none")) %>%
      as_tibble() %>%
      transmute(
        Contrast = as.character(contrast), StandardizedEffect = estimate, SE,
        DF = df, P_Value = p.value
      )
  } else {
    as.data.frame(emmeans::contrast(emmeans::emmeans(full_fit, ~ Sex), method = "pairwise", adjust = "none")) %>%
      as_tibble() %>%
      transmute(
        Contrast = as.character(contrast), StandardizedEffect = estimate, SE,
        DF = df, P_Value = p.value
      )
  }

  list(omnibus = omnibus, coefficients = coefficient_table, pairwise = pairwise)
}

fits <- list()
for (marker_name in marker_registry$Marker) {
  marker_data <- state_data %>% filter(Marker == marker_name)
  for (state_name in state_registry$State) {
    for (predictor in c("Role", "Sex")) {
      result <- fit_fixed_effect(marker_data, predictor, state_name)
      key <- paste(marker_name, state_name, predictor, sep = "::")
      fits[[key]] <- list(
        marker = marker_name,
        state = state_name,
        predictor = predictor,
        result = result
      )
    }
  }
}

omnibus_results <- bind_rows(lapply(fits, function(item) {
  item$result$omnibus %>% mutate(Marker = item$marker, .before = 1)
})) %>%
  left_join(marker_registry, by = "Marker") %>%
  left_join(state_registry, by = "State") %>%
  group_by(Predictor, State) %>%
  mutate(BH_FDR_WithinState = p.adjust(OmnibusP, method = "BH")) %>%
  ungroup() %>%
  group_by(Predictor) %>%
  mutate(BH_FDR_Global = p.adjust(OmnibusP, method = "BH")) %>%
  ungroup() %>%
  arrange(Predictor, match(State, state_registry$State), BH_FDR_WithinState)

coefficient_results <- bind_rows(lapply(fits, function(item) {
  item$result$coefficients %>%
    mutate(Marker = item$marker, State = item$state, Predictor = item$predictor, .before = 1)
})) %>%
  left_join(marker_registry, by = "Marker") %>%
  left_join(state_registry, by = "State") %>%
  group_by(Predictor, State) %>%
  mutate(BH_FDR_WithinState = p.adjust(CoefficientP, method = "BH")) %>%
  ungroup()

pairwise_results <- bind_rows(lapply(fits, function(item) {
  item$result$pairwise %>%
    mutate(Marker = item$marker, State = item$state, Predictor = item$predictor, .before = 1)
})) %>%
  left_join(marker_registry, by = "Marker") %>%
  left_join(state_registry, by = "State") %>%
  group_by(Predictor, State) %>%
  mutate(BH_FDR_WithinState = p.adjust(P_Value, method = "BH")) %>%
  ungroup()

family_icc <- read_csv(family_icc_file, show_col_types = FALSE, progress = FALSE) %>%
  filter(State %in% c("Pre_log2", "Post1h_log2", "Post3h_log2", "Delta1h", "Delta3h"))

family_pairs <- read_csv(pair_file, show_col_types = FALSE, progress = FALSE) %>%
  filter(State %in% c("Pre_log2", "Post1h_log2", "Post3h_log2", "Delta1h", "Delta3h"))

role_sex_crosstab <- state_data %>%
  distinct(Clinical_Subject_ID, Role, Sex) %>%
  count(Role, Sex, name = "Participants")

primary_states <- c("Pre", "Post3h", "Delta3h")
primary_omnibus <- omnibus_results %>% filter(State %in% primary_states)
primary_pairwise <- pairwise_results %>% filter(State %in% primary_states)
primary_family_icc <- family_icc %>% filter(State %in% c("Pre_log2", "Post3h_log2", "Delta3h"))
primary_family_pairs <- family_pairs %>% filter(State %in% c("Pre_log2", "Post3h_log2", "Delta3h"))

methods <- tribble(
  ~Item, ~Description,
  "Population", "Primary complete-case cohort of 50 participants from 18 families; no visually influential trajectory was excluded from inferential models.",
  "States", "Pre and Post values were log2-transformed; change was log2(Post/Pre). Primary reporting uses Pre, Post 3 h, and Post 3 h minus Pre; Post 1 h is retained in the complete workbook.",
  "Role model", "Standardized OC state or change ~ Role + (1|FamilyID). The omnibus P value compares models with and without Role using maximum likelihood.",
  "Sex-coded model", "Standardized OC state or change ~ Sex + (1|FamilyID). Gender code 1 was mapped to female and code 2 to male according to established project coding.",
  "Collinearity", "Sex and Role are perfectly collinear in this cohort: daughters and mothers are female, whereas fathers are male. Role and sex-coded models are therefore separate descriptions and cannot estimate independent effects.",
  "Family ICC", "Role-adjusted family ICC and role-stratified family-label permutation P values were imported from the validated osteocalcin family-state analysis.",
  "Family pairs", "DM, DF, and MF Spearman correlations, permutation P values, and BH-FDR values were imported from the validated family-state analysis.",
  "Multiplicity", "For Role and sex-coded omnibus tests, BH correction was applied across the three OC markers within each state; global BH-FDR across states is retained as sensitivity information."
)

writexl::write_xlsx(
  list(
    Primary_omnibus = primary_omnibus,
    Primary_pairwise = primary_pairwise,
    Primary_family_ICC = primary_family_icc,
    Primary_family_pairs = primary_family_pairs,
    All_omnibus = omnibus_results,
    All_coefficients = coefficient_results,
    All_pairwise = pairwise_results,
    All_family_ICC = family_icc,
    All_family_pairs = family_pairs,
    Role_sex_crosstab = role_sex_crosstab,
    Methods = methods
  ),
  file.path(output_dir, "Fig06_OC_State_Family_Sex_Role_Audit.xlsx")
)

write_csv(primary_omnibus, file.path(output_dir, "Primary_Role_Sex_Omnibus.csv"))
write_csv(primary_family_icc, file.path(output_dir, "Primary_Family_ICC.csv"))
write_csv(primary_family_pairs, file.path(output_dir, "Primary_Family_Pair_Similarity.csv"))

writeLines(
  c(
    "# Osteocalcin state, family, role, and sex-code audit",
    "",
    "The primary comparison covers Pre, Post 3 h, and log2(Post 3 h/Pre) for tOC, cOC, and cOC/tOC.",
    "Role and sex-coded models are reported separately because sex and family role are perfectly collinear in this family design.",
    "The complete 50-participant cohort is primary; the display-only outlier exclusions are not used for inference."
  ),
  file.path(output_dir, "README.md"), useBytes = TRUE
)

writeLines(
  c(
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Participants: ", n_distinct(state_data$Clinical_Subject_ID)),
    paste0("Families: ", n_distinct(state_data$FamilyID)),
    "",
    capture.output(sessionInfo())
  ),
  file.path(output_dir, "RUN_LOG.txt"), useBytes = TRUE
)

message("Completed osteocalcin family, role, and sex-code audit: ", output_dir)
