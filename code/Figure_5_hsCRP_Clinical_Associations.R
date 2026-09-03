# Figure 5A: clinical correlates of hsCRP.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.local_config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  file.path(getwd(), "Fig00_Config.R"),
  if (!is.na(.local_script)) file.path(dirname(.local_script), "Fig00_Config.R") else NA_character_
))
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (is.na(.local_config)) stop("Cannot find Fig00_Config.R.")
source(.local_config)
rm(.local_script, .local_config_candidates, .local_config)

p_load(digest, dplyr, ggplot2, lme4, lmerTest, patchwork, readr, scales,
       tibble, tidyr, writexl)

set.seed(20260816)

out_dir <- path_project("fig05", "analysis", "Fig05AB_HalfWidth_Candidates_v4_2026-08-29")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clinical_file <- path_intermediate("Clinical_Observed_SourcePool.csv")
registry_file <- path_analysis_ready("Metadata", "Clinical_Variable_Registry.csv")

output_files <- c(
  file.path(out_dir, "Fig05A_hsCRP_Clinical_Associations_HalfWidth_Candidate.pdf"),
  file.path(out_dir, "QA_Fig05A_hsCRP_Clinical_Associations_HalfWidth_Candidate.png"),
  file.path(out_dir, "Fig05A_hsCRP_Clinical_Associations_HalfWidth_SourceData.csv"),
  file.path(out_dir, "Fig05A_hsCRP_Clinical_Associations_HalfWidth_SourceData.xlsx"),
  file.path(out_dir, "PROVENANCE_Fig05A_HalfWidth.csv"),
  file.path(out_dir, "README_Fig05A_HalfWidth.md"),
  file.path(out_dir, "SESSION_INFO_Fig05A_HalfWidth.txt")
)
if (any(file.exists(output_files))) {
  message(
    "Overwriting current Figure 5A half-width candidate output(s): ",
    paste(basename(output_files[file.exists(output_files)]), collapse = ", ")
  )
}

rank_normalize <- function(x) {
  output <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  n <- sum(keep)
  if (n < 3L) return(output)
  ranks <- rank(x[keep], ties.method = "average")
  output[keep] <- qnorm((ranks - 0.5) / n)
  output
}

z_score <- function(x) {
  output <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  if (sum(keep) < 3L || stats::sd(x[keep]) == 0) return(output)
  output[keep] <- as.numeric(scale(x[keep]))
  output
}

clinical <- readr::read_csv(clinical_file, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    Family_Model_ID = factor(FamilyID),
    Role = factor(Role, levels = c("Daughter", "Mother", "Father"))
  )
if (nrow(clinical) != 179L || n_distinct(clinical$Family_Model_ID) != 99L) {
  stop("Expected the total 179-participant, 99-family clinical source pool.")
}

registry <- readr::read_csv(registry_file, show_col_types = FALSE) %>%
  filter(ReleaseRole == "Curated observed or source-derived phenotype") %>%
  transmute(
    Outcome = SourceVariable,
    OutcomeLabel = if_else(
      is.na(VariableLabel) | VariableLabel == "", SourceVariable, VariableLabel
    ),
    Unit = if_else(is.na(Unit) | Unit == "", "not documented", Unit),
    RegistryDomain = if_else(is.na(Domain) | Domain == "", "Other", Domain)
  ) %>%
  distinct(Outcome, .keep_all = TRUE)

# hsCRP is the sole predictor. Leptin, adiponectin and WBC remain eligible
# outcomes in the multiplicity family even when they are not all displayed.
excluded_outcomes <- c("Gender", "CRP_84m_G", "LnCRP_G")
candidate_outcomes <- intersect(registry$Outcome, names(clinical)) %>%
  setdiff(excluded_outcomes)

predictor_raw <- suppressWarnings(as.numeric(clinical$CRP_84m_G))
predictor_log <- ifelse(
  is.finite(predictor_raw) & predictor_raw > 0,
  log(predictor_raw), NA_real_
)

fit_one <- function(outcome) {
  outcome_raw <- suppressWarnings(as.numeric(clinical[[outcome]]))
  dat <- tibble(
    Family_Model_ID = clinical$Family_Model_ID,
    Role = clinical$Role,
    Predictor = predictor_log,
    Outcome = outcome_raw
  ) %>%
    filter(
      is.finite(Predictor), is.finite(Outcome),
      !is.na(Role), !is.na(Family_Model_ID)
    ) %>%
    mutate(Predictor_Z = z_score(Predictor), Outcome_Z = rank_normalize(Outcome)) %>%
    droplevels()

  if (nrow(dat) < 30L || n_distinct(dat$Family_Model_ID) < 10L ||
      n_distinct(dat$Outcome) < 5L || n_distinct(dat$Predictor) < 5L ||
      n_distinct(dat$Role) < 2L) return(NULL)

  model <- tryCatch(
    suppressWarnings(lmerTest::lmer(
      Outcome_Z ~ Predictor_Z + Role + (1 | Family_Model_ID),
      data = dat, REML = TRUE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    )),
    error = function(e) NULL
  )
  if (is.null(model)) return(NULL)
  coef_table <- summary(model)$coefficients
  if (!"Predictor_Z" %in% rownames(coef_table)) return(NULL)

  estimate <- coef_table["Predictor_Z", "Estimate"]
  standard_error <- coef_table["Predictor_Z", "Std. Error"]
  tibble(
    Outcome = outcome,
    EffectiveN = nrow(dat),
    Families = n_distinct(dat$Family_Model_ID),
    StandardizedEffect = estimate,
    SE = standard_error,
    CI_Low = estimate - 1.96 * standard_error,
    CI_High = estimate + 1.96 * standard_error,
    P_Value = coef_table["Predictor_Z", "Pr(>|t|)"],
    SingularFit = lme4::isSingular(model, tol = 1e-5)
  )
}

all_results <- bind_rows(lapply(candidate_outcomes, fit_one)) %>%
  left_join(registry, by = "Outcome") %>%
  mutate(BH_FDR = p.adjust(P_Value, method = "BH")) %>%
  arrange(BH_FDR, P_Value)

display_registry <- tribble(
  ~Outcome, ~Domain, ~ClinicalPhenotype, ~DisplayOrder,
  "Fat_p_G", "Adiposity", "Body-fat percentage", 1,
  "IntraperitonealFatKg_G", "Adiposity", "Visceral fat mass", 2,
  "Liverfat_G", "Ectopic fat", "Liver fat", 3,
  "Leptin_84m_G", "Adipokines", "Leptin", 4,
  "TRIGLY_G", "Serum lipids", "Triglycerides", 5,
  "HDLPlus_G", "Serum lipids", "HDL cholesterol", 6,
  "GGT", "Liver", "GGT", 7,
  "ALT", "Liver", "ALT", 8,
  "ICL_noutleir_G", "Muscle lipid", "Intramyocellular lipid", 9,
  "S_cOC_84m_G", "Bone turnover", "Carboxylated osteocalcin", 10,
  "cOC_totalOC_ratio_84m_G", "Bone turnover", "cOC/total OC ratio", 11
)

display_data <- all_results %>%
  inner_join(display_registry, by = "Outcome") %>%
  arrange(DisplayOrder) %>%
  mutate(
    Domain = factor(
      Domain,
      levels = c("Adiposity", "Ectopic fat", "Adipokines", "Serum lipids",
                 "Liver", "Muscle lipid", "Bone turnover", "Fitness")
    ),
    ClinicalPhenotype = factor(
      ClinicalPhenotype,
      levels = rev(display_registry$ClinicalPhenotype)
    ),
    Significance = case_when(
      BH_FDR < 0.001 ~ "***",
      BH_FDR < 0.01 ~ "**",
      BH_FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Significance = if_else(
      SingularFit & Significance != "", paste0(Significance, "\u2020"), Significance
    ),
    StarX = pmin(CI_High + 0.035, 1.03)
  )
if (nrow(display_data) != nrow(display_registry)) {
  missing <- setdiff(display_registry$Outcome, display_data$Outcome)
  stop("Missing displayed hsCRP outcome(s): ", paste(missing, collapse = ", "))
}

common_theme <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.text = element_text(colour = "#000000", size = 6.5),
    axis.title = element_text(colour = "#000000", size = 7),
    axis.line = element_line(colour = "#000000", linewidth = 0.25),
    axis.ticks = element_line(colour = "#000000", linewidth = 0.25),
    plot.title = element_text(
      colour = "#000000", size = 8, face = "bold", hjust = 0.5,
      margin = margin(b = 1.2, unit = "mm")
    ),
    plot.title.position = "plot",
    strip.background = element_blank(),
    panel.spacing.y = grid::unit(0.8, "mm")
  )

heatmap_plot <- ggplot(
  display_data,
  aes("hsCRP", ClinicalPhenotype, fill = StandardizedEffect)
) +
  geom_tile(colour = "#000000", linewidth = 0.18, width = 0.82) +
  geom_text(
    aes(label = Significance), colour = "#000000", family = "Arial",
    size = 2.45, fontface = "bold"
  ) +
  facet_grid(
    rows = vars(Domain), scales = "free_y", space = "free_y", switch = "y"
  ) +
  scale_fill_gradient2(
    low = "#3E6E9E", mid = "#FFFFFF", high = "#C8564B",
    midpoint = 0, limits = c(-0.9, 0.9), oob = scales::squish,
    breaks = c(-0.9, -0.45, 0, 0.45, 0.9), guide = "none"
  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  labs(title = "Association", x = NULL, y = NULL) +
  common_theme +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0, hjust = 1, colour = "#000000", size = 6.0, face = "bold"
    ),
    axis.text.x.top = element_text(colour = "#000000", size = 6.2, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(2, 0.2, 2, 1, unit = "mm")
  )

forest_plot <- ggplot(display_data, aes(StandardizedEffect, ClinicalPhenotype)) +
  geom_vline(xintercept = 0, linewidth = 0.22, linetype = "dashed", colour = "#000000") +
  geom_errorbarh(
    aes(xmin = CI_Low, xmax = CI_High), height = 0,
    linewidth = 0.34, colour = "#000000"
  ) +
  geom_point(size = 1.9, shape = 16, colour = "#C86B5A") +
  geom_text(
    aes(x = StarX, label = Significance), hjust = 0,
    colour = "#000000", family = "Arial", size = 2.45, fontface = "bold"
  ) +
  facet_grid(rows = vars(Domain), scales = "free_y", space = "free_y") +
  scale_y_discrete(position = "right") +
  scale_x_continuous(
    limits = c(-1.1, 1.1), breaks = c(-1, 0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(title = "Effect and 95% CI", x = "Standardized effect", y = NULL) +
  common_theme +
  theme(
    strip.text.y = element_blank(),
    axis.text.y.right = element_text(colour = "#000000", size = 6.0, hjust = 0),
    axis.ticks.y.right = element_blank(),
    axis.text.x = element_text(colour = "#000000", size = 5.8),
    axis.title.x = element_text(colour = "#000000", size = 6.2, face = "bold"),
    plot.title = element_text(
      colour = "#000000", size = 7.5, face = "bold", hjust = 0.5
    ),
    plot.margin = margin(2, 1.0, 2, 0.2, unit = "mm")
  )

legend_gradient <- tibble(x = 0, y = seq(-0.9, 0.9, length.out = 181))
legend_plot <- ggplot(legend_gradient, aes(x, y, fill = y)) +
  geom_tile(width = 0.18, height = 0.011) +
  scale_fill_gradient2(
    low = "#3E6E9E", mid = "#FFFFFF", high = "#C8564B",
    midpoint = 0, limits = c(-0.9, 0.9), guide = "none"
  ) +
  annotate(
    "text", x = 0.20, y = 1.28, label = "Effect\n(beta)",
    family = "Arial", size = 1.30, fontface = "bold", hjust = 0.5
  ) +
  annotate(
    "segment", x = 0.09, xend = 0.15,
    y = c(-0.9, -0.45, 0, 0.45, 0.9),
    yend = c(-0.9, -0.45, 0, 0.45, 0.9),
    linewidth = 0.22
  ) +
  annotate(
    "text", x = 0.18, y = c(-0.9, -0.45, 0, 0.45, 0.9),
    label = c("-0.9", "-0.45", "0", "0.45", "0.9"),
    family = "Arial", size = 1.18, hjust = 0
  ) +
  annotate(
    "text", x = -0.07, y = -1.12,
    label = "BH-FDR\n*<.05\n**<.01\n***<.001",
    family = "Arial", size = 1.22, hjust = 0, vjust = 1
  ) +
  coord_cartesian(xlim = c(-0.10, 0.75), ylim = c(-1.55, 1.48), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(4, 0, 4, 0, unit = "mm"))

combined <- heatmap_plot + forest_plot + legend_plot +
  patchwork::plot_layout(widths = c(1.35, 2.85, 0.72)) +
  patchwork::plot_annotation(
    title = "Clinical markers associated with hsCRP",
    theme = theme(
      plot.title = element_text(
        family = "Arial", colour = "#000000", size = 8.4,
        face = "bold", hjust = 0
      )
    )
  )

ggsave(output_files[1], combined, width = 68, height = 72, units = "mm", device = cairo_pdf)
ggsave(output_files[2], combined, width = 68, height = 72, units = "mm", dpi = 300, bg = "white")

source_export <- display_data %>%
  transmute(
    Domain = as.character(Domain),
    ClinicalPhenotype = as.character(ClinicalPhenotype),
    SourceVariable = Outcome,
    EffectiveN,
    Families,
    StandardizedEffect,
    SE,
    CI_Low,
    CI_High,
    P_Value,
    BH_FDR,
    SingularFit,
    MultiplicityFamilySize = nrow(all_results)
  )
readr::write_csv(source_export, output_files[3])
writexl::write_xlsx(
  list(
    Displayed_associations = source_export,
    Complete_hsCRP_screen = all_results,
    Analysis_parameters = tibble(
      Parameter = c("Predictor", "Predictor transform", "Outcome transform", "Model", "Multiplicity"),
      Value = c(
        "hsCRP (CRP_84m_G)", "Natural log followed by z score",
        "Rank-based inverse normal transform",
        "Outcome_Z ~ hsCRP_Z + Role + (1 | FamilyID)",
        paste0("Benjamini-Hochberg across ", nrow(all_results), " estimable curated outcomes")
      )
    )
  ),
  output_files[4]
)

provenance <- tibble(
  Item = c(
    "Clinical input", "Clinical input SHA256", "Registry input", "Registry input SHA256",
    "Cohort", "Model", "Multiplicity family", "Displayed outcomes", "Figure size"
  ),
  Value = c(
    clinical_file,
    digest::digest(file = clinical_file, algo = "sha256", serialize = FALSE),
    registry_file,
    digest::digest(file = registry_file, algo = "sha256", serialize = FALSE),
    "Expanded observed clinical source pool: 179 participants from 99 families",
    "rank-normalized phenotype ~ log-z hsCRP + Role + (1 | FamilyID)",
    paste0(nrow(all_results), " estimable curated outcomes; BH correction applied once across the complete hsCRP screen"),
    paste(display_registry$ClinicalPhenotype, collapse = "; "),
    "68 x 78 mm"
  )
)
readr::write_csv(provenance, output_files[5])

writeLines(c(
  "# Figure 5A hsCRP candidate",
  "",
  "- hsCRP is the sole predictor; WBC, leptin and adiponectin are not parallel predictors.",
  "- The display is restricted to the 11 clinical outcomes reaching BH-FDR < 0.05 in the current Figure 5B association profile.",
  "- Adiponectin and VO2max remain in the complete 106-outcome multiplicity family but are not displayed because they did not reach BH-FDR < 0.05.",
  "- Platelet count and the Hematology display block were removed.",
  "- Triglycerides and HDL cholesterol are grouped under Serum lipids.",
  "- Model: `rank-normalized phenotype ~ log-z hsCRP + Role + (1 | FamilyID)`.",
  paste0("- BH-FDR was calculated across all ", nrow(all_results), " estimable curated outcomes, not only the displayed subset."),
  "- Cohort: expanded observed clinical source pool (179 participants, 99 families). This is not the locked 81-participant manuscript cohort.",
  "- The original four-marker candidate and protected Adobe Illustrator artwork were not modified."
), output_files[6])
capture.output(sessionInfo(), file = output_files[7])

message("Figure 5A hsCRP half-width candidate updated: ", out_dir)
print(source_export)
