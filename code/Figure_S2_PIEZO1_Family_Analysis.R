# Four-panel PIEZO1-linked adipose DNA-methylation technical-audit summary.
# Existing model outputs are visualized without rerunning statistical models.

library(dplyr)
library(ggplot2)
library(patchwork)
library(readxl)
library(tibble)
library(writexl)

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  file.path(getwd(), "Fig00_Config.R"),
  if (!is.na(.local_script)) {
    file.path(dirname(normalizePath(.local_script, winslash = "/", mustWork = FALSE)), "Fig00_Config.R")
  } else {
    NA_character_
  }
))
.config_file <- .config_candidates[file.exists(.config_candidates)][1]
if (is.na(.config_file)) stop("Cannot locate Fig00_Config.R.")
source(.config_file)
root <- PROJECT_ROOT
rm(.local_script, .config_candidates, .config_file)
family_dir <- file.path(root, "fig02", "Family_Similarity_Figures")
audit_dir <- file.path(root, "fig02", "PIEZO1_Family_HighLow_Omics_Audit")
result_dir <- file.path(audit_dir, "Full_Results")
output_dir <- file.path(audit_dir, "Figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
lean_color <- "#45A7AE"
obese_color <- "#F28E45"
dm_color <- "#C76D7A"
df_color <- "#4A86C5"
positive_color <- "#D48791"
negative_color <- "#79A6C5"

source_file <- file.path(
  family_dir, "Fig02G_Five_Molecular_Features_SourceData.xlsx"
)
feature_summary <- readxl::read_excel(source_file, "FeatureSummary") %>%
  filter(Label == "PIEZO1")
participant_values <- readxl::read_excel(source_file, "ParticipantValues") %>%
  filter(Label == "PIEZO1") %>%
  mutate(
    Adiposity = factor(Adiposity, levels = c("Lean", "Obese")),
    Role = factor(Role, levels = c("Daughter", "Mother", "Father"))
  )
family_pairs <- readxl::read_excel(source_file, "FamilyPairValues") %>%
  filter(Label == "PIEZO1", Pair %in% c("DM", "DF")) %>%
  mutate(Pair = factor(Pair, levels = c("DM", "DF")))

if (nrow(feature_summary) != 1L || !nrow(participant_values) || !nrow(family_pairs)) {
  stop("PIEZO1 source-data extraction failed.")
}

base_theme <- theme_classic(base_family = font_family, base_size = 5.4) +
  theme(
    axis.title = element_text(face = "bold", size = 5.4),
    axis.text = element_text(size = 4.8, color = "#111111"),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, color = "#222222", linewidth = 0.35),
    plot.title = element_text(face = "bold", size = 6.4, hjust = 0),
    plot.subtitle = element_text(size = 4.6, color = "#333333", hjust = 0),
    plot.margin = margin(3, 3, 3, 3),
    legend.title = element_text(face = "bold", size = 5.0),
    legend.text = element_text(size = 4.7),
    legend.key.height = unit(2.6, "mm"),
    legend.key.width = unit(3.2, "mm")
  )

panel_a <- ggplot(
  participant_values,
  aes(x = BodyFatPercent, y = Value, color = Adiposity, shape = Role)
) +
  geom_smooth(
    aes(group = 1), method = "lm", se = TRUE,
    color = "#555555", fill = "#D7D7D7", linewidth = 0.45
  ) +
  geom_point(size = 1.55, stroke = 0.55, fill = "white") +
  scale_color_manual(values = c(Lean = lean_color, Obese = obese_color)) +
  scale_shape_manual(values = c(Daughter = 2, Mother = 1, Father = 0)) +
  labs(
    title = "A  Adiposity association",
    subtitle = sprintf(
      "Effect %.2f; BH-FDR %.3f",
      feature_summary$BodyFatEffect, feature_summary$BodyFatBH_FDR
    ),
    x = "Body fat (%)", y = "PIEZO1-linked DNAm\n(M value)",
    color = "Adiposity", shape = "Role"
  ) +
  coord_fixed(ratio = diff(range(participant_values$BodyFatPercent, na.rm = TRUE)) /
    diff(range(participant_values$Value, na.rm = TRUE)) * 0.95) +
  base_theme +
  theme(legend.position = "none")

panel_b <- ggplot(
  family_pairs,
  aes(x = ParentValue, y = DaughterValue, color = Pair)
) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5) +
  geom_point(size = 1.6, alpha = 0.92) +
  scale_color_manual(values = c(DM = dm_color, DF = df_color)) +
  labs(
    title = "B  Within-family similarity",
    subtitle = sprintf(
      "DM: r %.2f, q %.3f; DF: r %.2f, q %.3f",
      feature_summary$DM_PearsonR, feature_summary$DM_BH_FDR,
      feature_summary$DF_PearsonR, feature_summary$DF_BH_FDR
    ),
    x = "Parental M value", y = "Daughter M value", color = "Family pair"
  ) +
  coord_fixed() +
  base_theme +
  theme(legend.position = "none")

metabolites <- readRDS(file.path(
  result_dir, "Serum_Metabonomics_PIEZO1_Family_HighLow.rds"
)) %>%
  filter(BH_FDR < 0.05) %>%
  arrange(AdjustedEffectHighMinusLow) %>%
  mutate(
    ShortLabel = factor(FeatureID, levels = FeatureID),
    Direction = "Higher"
  )

proteins_all <- readRDS(file.path(
  result_dir, "Adipose_Proteomics_PIEZO1_Family_HighLow.rds"
)) %>%
  filter(BH_FDR < 0.10) %>%
  mutate(Direction = if_else(AdjustedEffectHighMinusLow > 0, "Higher", "Lower"))

functional_proteins <- c(
  "MGST3", "SDHC", "UQCC3", "HADH", "DNAJC11",
  "SULT1A1", "PRDX4", "NQO1", "HMOX1", "CFD"
)
proteins <- proteins_all %>%
  filter(FeatureLabel %in% functional_proteins) %>%
  arrange(AdjustedEffectHighMinusLow) %>%
  mutate(
    ShortLabel = if_else(
      FeatureLabel == "P20231;Q15661", "TPSB2/TPSAB1", FeatureLabel
    ),
    ShortLabel = factor(ShortLabel, levels = ShortLabel)
  )

effect_scale <- scale_fill_manual(
  values = c(Higher = positive_color, Lower = negative_color),
  breaks = c("Higher", "Lower"),
  labels = c("Higher in higher-M-value families", "Lower in higher-M-value families"),
  name = "Effect direction"
)

panel_c <- ggplot(
  metabolites,
  aes(x = AdjustedEffectHighMinusLow, y = ShortLabel, fill = Direction)
) +
  geom_vline(xintercept = 0, color = "#555555", linewidth = 0.28) +
  geom_col(width = 0.68, color = "#333333", linewidth = 0.15) +
  effect_scale +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(
    title = "C  Serum LDL-related features",
    subtitle = "Serum metabolomics; BH-FDR < 0.05",
    x = "Adjusted effect", y = NULL, fill = "Effect direction"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 4.25),
    aspect.ratio = 1
  )

panel_d <- ggplot(
  proteins,
  aes(x = AdjustedEffectHighMinusLow, y = ShortLabel, fill = Direction)
) +
  geom_vline(xintercept = 0, color = "#555555", linewidth = 0.28) +
  geom_col(width = 0.68, color = "#333333", linewidth = 0.15) +
  effect_scale +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.04))) +
  labs(
    title = "D  Adipose proteins",
    subtitle = "Function-prioritized; BH-FDR < 0.10",
    x = "Adjusted effect", y = NULL, fill = "Effect direction"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 4.35),
    aspect.ratio = 1
  )

legend_panel <- ggplot() +
  annotate("text", x = 2, y = 8, label = "Adiposity", hjust = 0,
           family = font_family, fontface = "bold", size = 1.65) +
  annotate("point", x = 19, y = 8, color = lean_color, shape = 1, size = 1.8, stroke = 0.55) +
  annotate("text", x = 23, y = 8, label = "Lean", hjust = 0,
           family = font_family, size = 1.55) +
  annotate("point", x = 35, y = 8, color = obese_color, shape = 1, size = 1.8, stroke = 0.55) +
  annotate("text", x = 39, y = 8, label = "Obese", hjust = 0,
           family = font_family, size = 1.55) +
  annotate("text", x = 55, y = 8, label = "Role", hjust = 0,
           family = font_family, fontface = "bold", size = 1.65) +
  annotate("point", x = 66, y = 8, color = "#222222", shape = 2, size = 1.8, stroke = 0.55) +
  annotate("text", x = 70, y = 8, label = "Daughter", hjust = 0,
           family = font_family, size = 1.55) +
  annotate("point", x = 84, y = 8, color = "#222222", shape = 1, size = 1.8, stroke = 0.55) +
  annotate("text", x = 88, y = 8, label = "Mother", hjust = 0,
           family = font_family, size = 1.55) +
  annotate("point", x = 104, y = 8, color = "#222222", shape = 0, size = 1.8, stroke = 0.55) +
  annotate("text", x = 108, y = 8, label = "Father", hjust = 0,
           family = font_family, size = 1.55) +
  annotate("text", x = 2, y = 3.5, label = "Family pair", hjust = 0,
           family = font_family, fontface = "bold", size = 1.65) +
  annotate("point", x = 24, y = 3.5, color = dm_color, shape = 16, size = 1.8) +
  annotate("text", x = 28, y = 3.5, label = "DM", hjust = 0,
           family = font_family, size = 1.55) +
  annotate("point", x = 39, y = 3.5, color = df_color, shape = 16, size = 1.8) +
  annotate("text", x = 43, y = 3.5, label = "DF", hjust = 0,
           family = font_family, size = 1.55) +
  annotate("text", x = 55, y = 3.5, label = "Effect", hjust = 0,
           family = font_family, fontface = "bold", size = 1.65) +
  annotate("point", x = 69, y = 3.5, color = positive_color, shape = 15, size = 2.2) +
  annotate("text", x = 74, y = 3.5, label = "Higher", hjust = 0,
           family = font_family, size = 1.55) +
  annotate("point", x = 91, y = 3.5, color = negative_color, shape = 15, size = 2.2) +
  annotate("text", x = 96, y = 3.5, label = "Lower", hjust = 0,
           family = font_family, size = 1.55) +
  coord_cartesian(xlim = c(0, 122), ylim = c(0, 10), clip = "off") +
  theme_void()

combined <- (panel_a | panel_b | panel_c | panel_d) / legend_panel +
  plot_layout(heights = c(1, 0.16), widths = c(1, 1, 1, 1)) +
  plot_annotation(
    title = "PIEZO1-linked adipose DNA methylation and cross-platform profiles",
    caption = paste(
      "Technical audit: the higher/lower M-value family grouping coincides",
      "with a documented adipose DNA-methylation design-version block."
    ),
    theme = theme(
      plot.title = element_text(
        family = font_family, face = "bold", size = 7.2, hjust = 0
      ),
      plot.caption = element_text(
        family = font_family, size = 4.5, color = "#555555", hjust = 0
      )
    )
  )

output_pdf <- file.path(
  output_dir, "Fig02S_PIEZO1_FourPanel_Summary_Candidate.pdf"
)
ggsave(
  output_pdf, combined,
  width = 190, height = 62, units = "mm", device = cairo_pdf
)

writexl::write_xlsx(
  list(
    `Panel A participants` = participant_values,
    `Panel B family pairs` = family_pairs,
    `Panel C metabolites` = metabolites,
    `Panel D proteins` = proteins,
    Parameters = tibble::tribble(
      ~Parameter, ~Value,
      "Figure size", "190 mm wide x 62 mm high; four panels in one row",
      "Panel A model", "PIEZO1-linked adipose DNAm ~ standardized body-fat percentage + Role; FamilyID duplicateCorrelation block",
      "Panel B", "Daughter-mother and daughter-father Pearson correlations; BH correction within family-pair context",
      "Panel C", "All serum metabolites with BH-FDR < 0.05",
      "Panel D", "Function-prioritized adipose proteins among BH-FDR < 0.10: mitochondrial bioenergetics, redox/lipid processing, oxidative stress, and complement",
      "Effect direction", "Higher-M-value families minus lower-M-value families",
      "Critical limitation", "Family group coincides with the documented adipose methylation design-version block"
    )
  ),
  file.path(output_dir, "Fig02S_PIEZO1_FourPanel_Summary_SourceData.xlsx")
)

message("Wrote: ", output_pdf)
