# Compact Figure S6B: parental-role contrasts in osteocalcin states.

.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.this_script <- if (length(.file_arg)) sub("^--file=", "", .file_arg[[1]]) else NA_character_
.script_dir <- if (!is.na(.this_script)) dirname(normalizePath(.this_script)) else file.path(getwd(), "code")
.config_file <- file.path(.script_dir, "Fig00_Config.R")
if (!file.exists(.config_file)) stop("Cannot locate Fig00_Config.R.")
source(.config_file)
rm(.file_arg, .this_script, .script_dir, .config_file)

p_load(dplyr, ggplot2, readxl, writexl)

output_dir <- path_project(
  "fig06", "analysis", "Fig06_OC_Role_Forest_Candidate_2026-08-30"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

audit_workbook <- path_project(
  "fig06", "analysis", "Fig06_OC_State_Family_Sex_Role_Audit_2026-08-30",
  "Fig06_OC_State_Family_Sex_Role_Audit.xlsx"
)
if (!file.exists(audit_workbook)) stop("Required audit workbook is missing: ", audit_workbook)

marker_colors <- c(tOC = "#D47E70", cOC = "#5AA79D", `cOC/tOC` = "#8E88B9")
primary_states <- c("Pre", "Post3h", "Delta3h")
state_labels <- c(
  Pre = "Pre-exercise", Post3h = "Post 3 h", Delta3h = "Post 3 h - Pre"
)

fdr_stars <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "",
    x < 0.001 ~ "***",
    x < 0.01 ~ "**",
    x < 0.05 ~ "*",
    TRUE ~ ""
  )
}

pairwise <- read_excel(audit_workbook, sheet = "All_pairwise") %>%
  filter(Predictor == "Role", State %in% primary_states) %>%
  mutate(
    CI_Lower = StandardizedEffect - qt(0.975, DF) * SE,
    CI_Upper = StandardizedEffect + qt(0.975, DF) * SE,
    StateLabel = state_labels[State],
    ContrastCode = recode(
      Contrast,
      "Daughter - Mother" = "DM",
      "Daughter - Father" = "DF",
      "Mother - Father" = "MF"
    ),
    ContrastCode = factor(ContrastCode, levels = c("DM", "DF", "MF")),
    RowLabel = paste(StateLabel, ContrastCode, sep = "  |  "),
    RowLabel = factor(RowLabel, levels = rev(unique(RowLabel))),
    MarkerShort = factor(
      recode(Marker, TotalOC = "tOC", cOC = "cOC", cOC_ratio = "cOC/tOC"),
      levels = c("cOC", "tOC", "cOC/tOC")
    ),
    Stars = fdr_stars(BH_FDR_WithinState),
    StarX = CI_Upper + 0.05
  )

forest_plot <- ggplot(pairwise, aes(x = StandardizedEffect, y = RowLabel, color = MarkerShort)) +
  geom_vline(xintercept = 0, linewidth = 0.32, linetype = "dashed", color = "#444444") +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0, linewidth = 0.50) +
  geom_point(size = 1.95) +
  geom_text(
    aes(x = StarX, label = Stars), color = "black", hjust = 0,
    size = 2.45, fontface = "bold", show.legend = FALSE
  ) +
  facet_grid(. ~ MarkerShort) +
  scale_color_manual(values = marker_colors, guide = "none") +
  scale_x_continuous(
    limits = c(-1.50, 1.90), breaks = c(-1.5, -0.75, 0, 0.75, 1.5),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Parental role-associated variation in osteocalcin profiles",
    x = "Standardized pairwise difference (95% CI)", y = NULL,
    caption = "DM, daughter-mother; DF, daughter-father; MF, mother-father. Stars denote pairwise BH-FDR: * < 0.05, ** < 0.01, *** < 0.001."
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "Arial", base_size = 6.5) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 5.8),
    axis.text.y = element_text(size = 5.45),
    axis.title.x = element_text(color = "black", face = "bold", size = 6.3),
    plot.title = element_text(face = "bold", size = 8.0, hjust = 0),
    plot.caption = element_text(size = 5.3, hjust = 0, margin = margin(t = 2.5, unit = "pt")),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.32),
    strip.text = element_text(face = "bold", size = 6.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.34),
    panel.spacing.x = grid::unit(2.0, "mm"),
    plot.margin = margin(2.5, 4.5, 2.0, 3.0, unit = "mm")
  )

ggsave(
  file.path(output_dir, "FigS06C_Role_Associated_Osteocalcin_Forest.pdf"),
  forest_plot, width = 190, height = 58, units = "mm", device = cairo_pdf
)

supp_b_dir <- path_project("fig06", "figS06B")
dir.create(supp_b_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path(supp_b_dir, "FigS06B_Parental_Role_Associated_Osteocalcin_Profiles.pdf"),
  forest_plot, width = 190, height = 58, units = "mm", device = cairo_pdf
)

writexl::write_xlsx(
  list(
    Role_pairwise_contrasts = pairwise %>%
      select(
        Marker, MarkerShort, State, StateLabel, Contrast, ContrastCode,
        StandardizedEffect, SE, DF, CI_Lower, CI_Upper, P_Value,
        BH_FDR_WithinState
      )
  ),
  file.path(output_dir, "FigS06C_Role_Associated_Osteocalcin_Forest_Source_Data.xlsx")
)

writeLines(
  c(
    "# Figure S6B",
    "",
    "Compact forest plot of family-role pairwise contrasts for osteocalcin states and responses.",
    "DM denotes daughter-mother, DF daughter-father, and MF mother-father.",
    "Effects and 95% confidence intervals come from family-random-intercept mixed models.",
    "Stars denote pairwise BH-FDR thresholds."
  ),
  file.path(output_dir, "README.md"), useBytes = TRUE
)

message("Created compact Figure S6C forest plot: ", output_dir)
