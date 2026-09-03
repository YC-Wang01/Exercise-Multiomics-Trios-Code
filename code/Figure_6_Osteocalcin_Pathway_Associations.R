# Candidate Figure 6C: longitudinal osteocalcin-associated proteomic pathways.

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

p_load(digest, dplyr, ggplot2, readr, scales, tibble, tidyr, writexl)

analysis_id <- "Fig06C_OC_Pre_Post_Delta_Hallmark_Bubble_Candidate_2026-08-29"
figure_dir <- path_project("fig06", "fig06C")
analysis_dir <- path_project("fig06", "analysis", analysis_id)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

source_dir <- path_results("Six_Markers_Three_State_Omics_2026-08-22")
hallmark_file <- file.path(source_dir, "Hallmark_Selected.csv")
summary_file <- file.path(source_dir, "Combined_Hallmark_FDR_Summary.csv")
if (!file.exists(hallmark_file)) stop("Hallmark result table is missing.")
if (!file.exists(summary_file)) stop("Hallmark summary table is missing.")

marker_order <- c("cOC", "TotalOC")
marker_labels <- c(TotalOC = "Total osteocalcin (tOC)", cOC = "Carboxylated osteocalcin (cOC)")
model_order <- c("Pre_Pre", "Post_Post", "Delta_Delta")
state_labels <- c(
  Pre_Pre = "Pre-exercise",
  Post_Post = "Post 3 h",
  Delta_Delta = "Post3h - Pre"
)
tissue_order <- c("Adipose", "Muscle")
tissue_labels <- c(Adipose = "Adipose tissue", Muscle = "Skeletal muscle")

pathway_registry <- tribble(
  ~FunctionalGroup, ~Pathway, ~PathwayOrder,
  "Energy metabolism", "Oxidative phosphorylation", 1L,
  "Energy metabolism", "Fatty acid metabolism", 2L,
  "Energy metabolism", "Glycolysis", 3L,
  "Energy metabolism", "Cholesterol homeostasis", 4L,
  "Heme and hemostasis", "Heme metabolism", 5L,
  "Heme and hemostasis", "Coagulation", 6L,
  "Heme and hemostasis", "Complement", 7L,
  "Growth and remodeling", "Adipogenesis", 8L,
  "Growth and remodeling", "Epithelial mesenchymal transition", 9L,
  "Growth and remodeling", "G2m checkpoint", 10L
)

hallmark <- read_csv(hallmark_file, show_col_types = FALSE, progress = FALSE) %>%
  filter(
    Dataset %in% c("Adipose_Proteomics", "Muscle_Proteomics"),
    ClinicalMarker %in% marker_order,
    Model %in% model_order,
    Pathway %in% pathway_registry$Pathway
  ) %>%
  transmute(
    Dataset,
    Tissue,
    Model,
    Timepoint,
    ClinicalMarker,
    ClinicalMarkerLabel,
    Subjects,
    Families,
    Pathway,
    Size,
    NES = suppressWarnings(as.numeric(NES)),
    P_Value = suppressWarnings(as.numeric(P_Value)),
    BH_FDR = suppressWarnings(as.numeric(BH_FDR)),
    LeadingEdge
  ) %>%
  distinct(Dataset, Tissue, Model, ClinicalMarker, Pathway, .keep_all = TRUE)

plot_grid <- tidyr::crossing(
  Tissue = tissue_order,
  ClinicalMarker = marker_order,
  Model = model_order,
  Pathway = pathway_registry$Pathway
) %>%
  left_join(
    hallmark,
    by = c("Tissue", "ClinicalMarker", "Model", "Pathway")
  ) %>%
  left_join(pathway_registry, by = "Pathway") %>%
  mutate(
    MarkerLabel = factor(
      unname(marker_labels[ClinicalMarker]),
      levels = unname(marker_labels[marker_order])
    ),
    TissueLabel = factor(
      unname(tissue_labels[Tissue]),
      levels = unname(tissue_labels[tissue_order])
    ),
    StateLabel = factor(
      unname(state_labels[Model]),
      levels = unname(state_labels[model_order])
    ),
    Pathway = factor(Pathway, levels = rev(pathway_registry$Pathway)),
    FDRTier = case_when(
      is.finite(BH_FDR) & BH_FDR < 0.05 ~ "***",
      is.finite(BH_FDR) & BH_FDR < 0.10 ~ "**",
      is.finite(BH_FDR) & BH_FDR < 0.20 ~ "*",
      TRUE ~ ""
    ),
    DisplayNES = if_else(is.finite(BH_FDR) & BH_FDR < 0.20, NES, NA_real_),
    DisplayStrength = if_else(
      is.finite(BH_FDR) & BH_FDR < 0.20,
      pmin(-log10(pmax(BH_FDR, 1e-12)), 10),
      NA_real_
    )
  )

displayed <- plot_grid %>%
  filter(is.finite(DisplayNES)) %>%
  arrange(TissueLabel, MarkerLabel, StateLabel, BH_FDR)

if (nrow(displayed) == 0L) stop("No Hallmark associations met BH-FDR < 0.20.")

nes_limit <- max(3, ceiling(max(abs(displayed$NES), na.rm = TRUE) * 2) / 2)

plot_bubble <- ggplot(plot_grid, aes(x = StateLabel, y = Pathway)) +
  geom_point(
    aes(fill = DisplayNES, size = DisplayStrength),
    shape = 21,
    colour = "#303030",
    stroke = 0.25,
    na.rm = TRUE
  ) +
  geom_text(
    data = displayed,
    aes(label = FDRTier),
    family = "Arial",
    fontface = "bold",
    colour = "black",
    size = 2.15,
    vjust = 0.40,
    na.rm = TRUE
  ) +
  facet_grid(
    rows = vars(TissueLabel),
    cols = vars(MarkerLabel),
    switch = "y"
  ) +
  scale_fill_gradient2(
    low = "#82A9C6",
    mid = "white",
    high = "#D9919E",
    midpoint = 0,
    limits = c(-nes_limit, nes_limit),
    oob = scales::squish,
    name = "GSEA NES"
  ) +
  scale_size_continuous(
    range = c(2.0, 6.0),
    limits = c(0.7, 10),
    breaks = c(1, 2, 5, 10),
    name = expression(-log[10]("BH-FDR"))
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(23, "mm"),
      barwidth = grid::unit(3.8, "mm")
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      keyheight = grid::unit(3.4, "mm"),
      keywidth = grid::unit(5.2, "mm")
    )
  ) +
  labs(
    title = "Osteocalcin-associated proteomic programs across molecular states",
    x = NULL,
    y = NULL,
    caption = paste0(
      "Stars denote BH-FDR tiers: ***BH-FDR < 0.05, ",
      "**BH-FDR < 0.10, and *BH-FDR < 0.20."
    )
  ) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
    panel.spacing.x = grid::unit(0.8, "mm"),
    panel.spacing.y = grid::unit(0.8, "mm"),
    axis.text.x = element_text(size = 6.0, colour = "black", lineheight = 0.9),
    axis.text.y = element_text(size = 6.0, colour = "black"),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.35),
    strip.text.x = element_text(size = 6.5, face = "bold", colour = "black"),
    strip.text.y.left = element_text(size = 6.5, face = "bold", colour = "black", angle = 90),
    strip.placement = "outside",
    plot.title = element_text(size = 8.5, face = "bold", hjust = 0),
    plot.caption = element_text(size = 5.5, hjust = 0, colour = "black"),
    legend.position = "right",
    legend.justification = "center",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.y = grid::unit(1.0, "mm"),
    legend.title = element_text(size = 5.4, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 5.0, colour = "black"),
    plot.margin = margin(2.5, 2.5, 2.5, 2.5, unit = "mm")
  )

figure_pdf <- file.path(figure_dir, "Fig06C_OC_Pre_Post_Delta_Hallmark_Bubble_Candidate.pdf")
preview_png <- file.path(analysis_dir, "Fig06C_OC_Pre_Post_Delta_Hallmark_Bubble_Candidate.png")
ggsave(figure_pdf, plot_bubble, width = 190, height = 96, units = "mm", device = cairo_pdf)
ggsave(preview_png, plot_bubble, width = 190, height = 96, units = "mm", dpi = 350, bg = "white")

summary_export <- read_csv(summary_file, show_col_types = FALSE, progress = FALSE) %>%
  filter(
    Dataset %in% c("Adipose_Proteomics", "Muscle_Proteomics"),
    ClinicalMarker %in% marker_order,
    Model %in% model_order
  ) %>%
  arrange(ClinicalMarker, Model, Dataset)

source_checksums <- tibble(
  Source = c("Selected Hallmark GSEA results", "Hallmark model summary"),
  ProjectRelativePath = c(
    file.path("results", "Six_Markers_Three_State_Omics_2026-08-22", basename(hallmark_file)),
    file.path("results", "Six_Markers_Three_State_Omics_2026-08-22", basename(summary_file))
  ),
  SHA256 = c(
    digest::digest(file = hallmark_file, algo = "sha256"),
    digest::digest(file = summary_file, algo = "sha256")
  )
)

workbook_file <- file.path(analysis_dir, "Fig06C_OC_Pre_Post_Delta_Hallmark_Source_Data.xlsx")
writexl::write_xlsx(
  list(
    Displayed_pathways = displayed %>%
      transmute(
        Tissue,
        Marker = as.character(MarkerLabel),
        MolecularState = as.character(StateLabel),
        Pathway = as.character(Pathway),
        NES,
        P_Value,
        BH_FDR,
        FDRTier,
        Subjects,
        Families,
        LeadingEdge
      ),
    Model_summary = summary_export,
    Pathway_registry = pathway_registry,
    Provenance = source_checksums
  ),
  workbook_file
)

readme <- c(
  "# Candidate Figure 6C: longitudinal osteocalcin-associated proteomic pathways",
  "",
  "## Design",
  "",
  "- Pre-exercise: Pre omics abundance associated with Pre osteocalcin.",
  "- Post 3 h: Post3h omics abundance associated with Post3h osteocalcin.",
  "- Post3h - Pre: molecular change associated with osteocalcin change.",
  "- Each feature model included family role; FamilyID was used as the limma duplicateCorrelation block.",
  "- Hallmark GSEA ranked all detected proteins by the moderated t statistic for the osteocalcin term.",
  "- BH correction was performed within each dataset, marker, molecular state, and time point.",
  "",
  "## Display",
  "",
  "- Fill color: GSEA normalized enrichment score (NES).",
  "- Bubble size: -log10(BH-FDR), capped at 10.",
  "- Stars: *** BH-FDR <0.05; ** <0.10; * <0.20.",
  "- Cells without BH-FDR <0.20 are blank.",
  "",
  "## Scope boundary",
  "",
  "The main candidate uses tissue proteomics because it supports directional GSEA in all three states. Microarray data are Pre-only. Methylation methylglm enrichment is non-directional and should be shown separately if retained. Post-Post is a contemporaneous association; only the change-change model evaluates covariation between osteocalcin and molecular responses.",
  "",
  "## Reproduction",
  "",
  "Run: Rscript code/Fig06C_OC_Pre_Post_Delta_Hallmark_Bubble_Candidate.R",
  paste0("Analysis ID: ", analysis_id)
)
writeLines(readme, file.path(analysis_dir, "README.md"), useBytes = TRUE)

writeLines(
  c(
    paste0("Analysis ID: ", analysis_id),
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Figure: ", figure_pdf),
    paste0("Workbook: ", workbook_file),
    "",
    capture.output(sessionInfo())
  ),
  file.path(analysis_dir, "RUN_LOG.txt"),
  useBytes = TRUE
)

message("Created candidate Figure 6C: ", figure_pdf)
