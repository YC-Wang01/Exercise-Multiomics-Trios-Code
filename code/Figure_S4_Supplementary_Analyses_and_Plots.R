# Assemble the current Figure S4 panel series from verified Figure 4 outputs.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  file.path(getwd(), "Fig00_Config.R"),
  if (!is.na(.local_script)) file.path(
    dirname(normalizePath(.local_script, winslash = "/", mustWork = FALSE)),
    "Fig00_Config.R"
  ) else NA_character_
))
.config_file <- .config_candidates[file.exists(.config_candidates)][1]
if (is.na(.config_file)) stop("Cannot locate Fig00_Config.R.")
source(.config_file)
rm(.local_script, .config_candidates, .config_file)

p_load(dplyr, ggplot2, openxlsx, patchwork, readr)

analysis_id <- "FigS04_Current_Series_2026-08-28"
figure_root <- path_project("fig04")
pathway_dir <- file.path(figure_root, "fig04B")
feature_dir <- file.path(figure_root, "fig04D")
output_dir <- file.path(figure_root, "figS4", "Current_Series")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  S4A_pdf = file.path(
    pathway_dir, "FigS04B_All_Hallmark_Pre_Post_Delta_FDR_Tiers.pdf"
  ),
  S4A_png = file.path(
    pathway_dir, "FigS04B_All_Hallmark_Pre_Post_Delta_FDR_Tiers.png"
  ),
  S4A_source = file.path(
    pathway_dir, "Table_Fig04B_Hallmark_Pathway_Associations.xlsx"
  ),
  S4B_pdf = file.path(pathway_dir, "FigS04B_Two_Core_GSEA_Profiles.pdf"),
  S4B_png = file.path(pathway_dir, "FigS04B_Two_Core_GSEA_Profiles.png"),
  S4B_source = file.path(pathway_dir, "Fig04B_Core_GSEA_Profile_SourceData.xlsx"),
  S4C_source = file.path(feature_dir, "Fig04D_Feature_Level_SourceData.csv"),
  S4C_links = file.path(
    feature_dir, "Table_Fig04D_Feature_Level_Insulin_Associations.xlsx"
  )
)

if (!all(file.exists(inputs))) {
  stop(
    "Missing required input: ",
    paste(inputs[!file.exists(inputs)], collapse = "; ")
  )
}

copy_verified <- function(source, destination) {
  copied <- file.copy(source, destination, overwrite = TRUE, copy.date = TRUE)
  if (!isTRUE(copied)) stop("Failed to copy: ", source)
  if (!identical(unname(tools::md5sum(source)), unname(tools::md5sum(destination)))) {
    stop("Checksum mismatch after copy: ", destination)
  }
}

outputs <- c(
  S4A_pdf = file.path(
    output_dir, "FigS04A_Complete_Hallmark_Pathway_Associations.pdf"
  ),
  S4A_png = file.path(
    output_dir, "FigS04A_Complete_Hallmark_Pathway_Associations.png"
  ),
  S4A_source = file.path(output_dir, "FigS04A_SourceData.xlsx"),
  S4B_pdf = file.path(
    output_dir, "FigS04B_Heme_Metabolism_and_Myogenesis_GSEA_Profiles.pdf"
  ),
  S4B_png = file.path(
    output_dir, "FigS04B_Heme_Metabolism_and_Myogenesis_GSEA_Profiles.png"
  ),
  S4B_source = file.path(output_dir, "FigS04B_SourceData.xlsx")
)

copy_verified(inputs[["S4A_pdf"]], outputs[["S4A_pdf"]])
copy_verified(inputs[["S4A_png"]], outputs[["S4A_png"]])
copy_verified(inputs[["S4A_source"]], outputs[["S4A_source"]])
copy_verified(inputs[["S4B_pdf"]], outputs[["S4B_pdf"]])
copy_verified(inputs[["S4B_png"]], outputs[["S4B_png"]])
copy_verified(inputs[["S4B_source"]], outputs[["S4B_source"]])

candidate_spec <- tibble::tribble(
  ~Dataset, ~FeatureLabel, ~DisplayFeature,
  "Adipose_Proteomics", "GLIPR2", "GLIPR2 (adipose)",
  "Adipose_Proteomics", "CD276", "CD276 (adipose)",
  "Adipose_Proteomics", "AUH", "AUH (adipose)",
  "Muscle_Proteomics", "SUCLG2", "SUCLG2 (muscle)"
)

display_results <- read_csv(inputs[["S4C_source"]], show_col_types = FALSE) %>%
  inner_join(candidate_spec, by = c("Dataset", "FeatureLabel"), suffix = c("", ".selected")) %>%
  mutate(
    DisplayFeature = coalesce(DisplayFeature.selected, DisplayFeature),
    PredictorShort = factor(PredictorShort, levels = c("HOMA-IR", "Matsuda ISI")),
    StateShort = factor(
      StateShort, levels = c("Pre-exercise", "Post3h", "Response")
    ),
    DisplayFeature = factor(
      DisplayFeature, levels = rev(candidate_spec$DisplayFeature)
    ),
    FDRTier = case_when(
      BH_FDR < 0.05 ~ "***",
      BH_FDR < 0.10 ~ "**",
      BH_FDR < 0.20 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(-any_of("DisplayFeature.selected"))

expected_rows <- nrow(candidate_spec) * 6L
if (nrow(display_results) != expected_rows) {
  stop(
    "Selected molecular model grid is incomplete: expected ", expected_rows,
    " rows, found ", nrow(display_results), "."
  )
}

pathway_links <- read.xlsx(inputs[["S4C_links"]], sheet = "Pathway_links") %>%
  filter(as.character(DisplayFeature) %in% candidate_spec$DisplayFeature) %>%
  mutate(
    PlotPathwayLink = dplyr::recode(
      PlotPathwayLink,
      "Il2 stat5 signaling" = "IL-2/STAT5 signaling"
    ),
    DisplayFeature = factor(
      as.character(DisplayFeature), levels = rev(candidate_spec$DisplayFeature)
    )
  )

effect_limit <- ceiling(max(abs(display_results$Effect), na.rm = TRUE) * 10) / 10

p_heat <- ggplot(
  display_results,
  aes(x = StateShort, y = DisplayFeature, fill = Effect)
) +
  geom_tile(colour = "black", linewidth = 0.25) +
  geom_text(
    aes(label = FDRTier), family = "Arial", fontface = "bold",
    size = 3.1, colour = "black"
  ) +
  facet_grid(. ~ PredictorShort) +
  scale_fill_gradient2(
    low = "#7799BE", mid = "white", high = "#D88E89", midpoint = 0,
    limits = c(-effect_limit, effect_limit), name = "Coefficient"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_family = "Arial", base_size = 7) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.30),
    panel.spacing.x = grid::unit(0.8, "mm"),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.30),
    strip.text = element_text(size = 7.2, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 6.3, colour = "black"),
    axis.text.y = element_text(size = 6.5, colour = "black"),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.7),
    legend.key.height = grid::unit(12, "mm"),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  )

link_colors <- c(
  "Fatty-acid metabolism" = "#E7B55D",
  "Oxidative phosphorylation" = "#7FA9C4",
  "Other Hallmark pathway" = "#B8A8CF",
  "None" = "#E1E1E1"
)

p_links <- ggplot(
  pathway_links,
  aes(x = 1, y = DisplayFeature, fill = CorePathwayLink)
) +
  geom_tile(width = 0.22, colour = "black", linewidth = 0.25) +
  geom_text(
    aes(x = 1.18, label = PlotPathwayLink),
    hjust = 0, family = "Arial", size = 2.3, colour = "black"
  ) +
  scale_fill_manual(values = link_colors, guide = "none") +
  scale_x_continuous(limits = c(0.86, 2.65), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Hallmark leading-edge relation", x = NULL, y = NULL) +
  theme_void(base_family = "Arial", base_size = 7) +
  theme(
    plot.title = element_text(size = 6.5, face = "bold", hjust = 0),
    plot.margin = margin(1, 1, 1, 0.5, unit = "mm")
  )

combined <- p_heat + p_links +
  plot_layout(widths = c(1.8, 0.95)) +
  plot_annotation(
    title = "Selected molecular associations",
    caption = paste0(
      "Stars denote BH-FDR tiers: *** < 0.05, ** < 0.10, and * < 0.20. ",
      "Response denotes Post3h minus Pre."
    ),
    theme = theme(
      plot.title = element_text(
        family = "Arial", size = 8.5, face = "bold", hjust = 0.5
      ),
      plot.caption = element_text(
        family = "Arial", size = 5.5, hjust = 0, colour = "black"
      ),
      plot.margin = margin(1, 1, 1, 1, unit = "mm")
    )
  )

s4c_pdf <- file.path(
  output_dir, "FigS04C_Representative_Molecular_Associations.pdf"
)
s4c_png <- sub("\\.pdf$", ".png", s4c_pdf)
s4c_source <- file.path(output_dir, "FigS04C_SourceData.xlsx")

ggsave(s4c_pdf, combined, width = 190, height = 58, units = "mm", device = cairo_pdf)
ggsave(s4c_png, combined, width = 190, height = 58, units = "mm", dpi = 300)

write.xlsx(
  list(
    Displayed_models = display_results %>%
      mutate(across(where(is.factor), as.character)),
    Hallmark_links = pathway_links %>%
      mutate(across(where(is.factor), as.character)),
    Selection_notes = data.frame(
      Item = c(
        "GLIPR2", "CD276", "AUH", "SUCLG2", "Model", "Multiplicity"
      ),
      Description = c(
        "Post3h adipose association with Matsuda ISI at BH-FDR < 0.05; significant Hallmark leading-edge member",
        "Post3h adipose association with Matsuda ISI at BH-FDR < 0.05; significant Hallmark leading-edge member",
        "Post3h adipose association with Matsuda ISI at BH-FDR < 0.10; fatty-acid-metabolism leading-edge protein",
        "Post3h muscle association with Matsuda ISI at BH-FDR < 0.20; fatty-acid-metabolism leading-edge protein",
        "Feature ~ rank-normalized insulin metric + Role; FamilyID duplicateCorrelation",
        "BH correction within dataset x state x predictor x model tier"
      )
    )
  ),
  s4c_source,
  overwrite = TRUE
)

readme_file <- file.path(output_dir, "README.md")
writeLines(c(
  "# Current Figure S4 series",
  "",
  paste0("Analysis ID: `", analysis_id, "`"),
  "",
  "## Panels",
  "",
  "- `FigS04A_Complete_Hallmark_Pathway_Associations.pdf`: all Hallmark pathways with BH-FDR < 0.20 in at least one displayed tissue, predictor, or state model.",
  "- `FigS04B_Heme_Metabolism_and_Myogenesis_GSEA_Profiles.pdf`: one-row GSEA profiles for the two supplementary pathways, heme metabolism and myogenesis.",
  "- `FigS04C_Representative_Molecular_Associations.pdf`: role-adjusted scatterplots for TNMD, AUH, LDHD, and SUCLG2 against baseline Matsuda ISI.",
  "",
  "## Interpretation boundary",
  "",
  "TNMD is a Pre-exercise adipose-transcriptomic association. AUH, LDHD, and SUCLG2 are Post3h absolute-profile associations and must not be described as exercise-induced molecular changes because their response models did not reach BH-FDR < 0.20.",
  "",
  "## Provenance",
  "",
  "S4A and S4B are checksum-verified copies of the current vector source panels generated in `fig04/fig04B`. S4C is generated by `code/Figure_S4_Molecular_Association_Scatterplots.R` from the current analysis-ready matrices and formal model table."
), readme_file)

log_file <- file.path(output_dir, "run_log.txt")
writeLines(c(
  paste0("Analysis ID: ", analysis_id),
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("Output directory: ", output_dir),
  "",
  "Inputs:",
  paste0(names(inputs), ": ", unname(inputs)),
  "",
  "Input MD5:",
  paste0(names(inputs), ": ", unname(tools::md5sum(inputs))),
  "",
  "Outputs:",
  paste0(
    c(names(outputs), "S4C_pdf", "S4C_png", "S4C_source", "README"),
    ": ",
    c(unname(outputs), s4c_pdf, s4c_png, s4c_source, readme_file)
  ),
  "",
  capture.output(sessionInfo())
), log_file)

s4c_scatter_script <- file.path(
  dirname(CONFIG_FILE), "Figure_S4_Molecular_Association_Scatterplots.R"
)
if (!file.exists(s4c_scatter_script)) {
  stop("Missing current Figure S4C scatterplot script: ", s4c_scatter_script)
}
sys.source(s4c_scatter_script, envir = new.env(parent = globalenv()))

message("Completed current Figure S4 series in: ", output_dir)
