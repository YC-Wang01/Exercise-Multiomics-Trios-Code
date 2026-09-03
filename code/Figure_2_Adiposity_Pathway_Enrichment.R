# Candidate Figure 2 serum and tissue bubble plots.
# Existing Figure 2 files are not overwritten.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.local_config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  file.path(getwd(), "Fig00_Config.R"),
  if (!is.na(.local_script)) {
    file.path(dirname(normalizePath(.local_script, winslash = "/", mustWork = FALSE)), "Fig00_Config.R")
  } else NA_character_
))
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (is.na(.local_config)) stop("Cannot find CellMetabolism_Transfer code/Fig00_Config.R.")
source(.local_config)
rm(.local_script, .local_config_candidates, .local_config)

p_load(dplyr, ggplot2, patchwork, readxl, readr, stringr, tibble, tidyr)

analysis_id <- "Fig02_Serum_Tissue_Pathway_Bubbles_Candidate_2026-08-24"
output_dir <- path_project("fig02", "Pathway_Bubble_Candidates")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

font_title_pt <- 9
font_section_pt <- 7
font_body_pt <- 7
font_small_pt <- 6
pt_to_mm <- 1 / ggplot2::.pt

direction_colors <- c(
  "Higher in lean" = "#4C90C6",
  "Higher in obese" = "#E6816F"
)

make_bubble_plot <- function(
    display,
    selection,
    column_levels,
    column_labels,
    title,
    right_header,
    width,
    height,
    legend_position = "right",
    column_positions = NULL) {
  selection <- selection %>%
    mutate(
      Row = rev(seq_len(n())),
      ItemLabel = stringr::str_wrap(Item, width = 43)
    )

  if (is.null(column_positions)) {
    column_positions <- seq_along(column_levels)
  }
  if (length(column_positions) != length(column_levels)) {
    stop("column_positions must match column_levels.")
  }
  names(column_positions) <- column_levels

  display <- tidyr::crossing(selection, Column = column_levels) %>%
    left_join(display, by = c("Item", "Column")) %>%
    mutate(
      Column = factor(Column, levels = column_levels),
      ColumnX = unname(column_positions[as.character(Column)]),
      FDRStrength = if_else(
        is.finite(BH_FDR),
        pmin(-log10(pmax(BH_FDR, 1e-12)), 8),
        NA_real_
      )
    )

  block_rectangles <- selection %>%
    group_by(Block) %>%
    summarise(
      YMin = min(Row) - 0.5,
      YMax = max(Row) + 0.5,
      Y = mean(c(YMin, YMax)),
      .groups = "drop"
    )

  column_step <- if (length(column_positions) > 1) {
    min(diff(sort(unique(column_positions))))
  } else {
    1
  }
  panel_xmin <- min(column_positions) - column_step / 2
  panel_xmax <- max(column_positions) + column_step / 2
  right_x <- panel_xmax + 0.10
  x_limit_right <- right_x + 4.65
  block_label_x <- panel_xmin - 0.15
  x_limit_left <- panel_xmin - 1.40

  ggplot(display, aes(x = ColumnX, y = Row)) +
    geom_rect(
      data = block_rectangles,
      aes(xmin = panel_xmin, xmax = panel_xmax, ymin = YMin, ymax = YMax),
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.38
    ) +
    geom_point(
      data = display %>% filter(is.finite(BH_FDR)),
      aes(size = FDRStrength, fill = Direction),
      shape = 21,
      color = "black",
      stroke = 0.32,
      alpha = 0.98
    ) +
    geom_text(
      data = selection,
      aes(x = right_x, y = Row, label = ItemLabel),
      inherit.aes = FALSE,
      hjust = 0,
      family = "Arial",
      size = font_body_pt * pt_to_mm,
      lineheight = 0.90,
      color = "black"
    ) +
    geom_text(
      data = block_rectangles,
      aes(x = block_label_x, y = Y, label = Block),
      inherit.aes = FALSE,
      hjust = 1,
      family = "Arial",
      fontface = "bold",
      size = font_section_pt * pt_to_mm,
      lineheight = 0.90,
      color = "black"
    ) +
    annotate(
      "text",
      x = right_x,
      y = nrow(selection) + 0.72,
      label = right_header,
      hjust = 0,
      family = "Arial",
      fontface = "bold",
      size = font_section_pt * pt_to_mm,
      color = "black"
    ) +
    scale_x_continuous(
      breaks = column_positions,
      labels = column_labels,
      limits = c(x_limit_left, x_limit_right),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = c(0.35, nrow(selection) + 0.88),
      expand = expansion(mult = 0)
    ) +
    scale_fill_manual(values = direction_colors, name = "Association direction") +
    scale_size_continuous(
      range = c(2.0, 5.8),
      limits = c(0, 8),
      breaks = c(2, 4, 6, 8),
      labels = c("2", "4", "6", ">=8"),
      name = expression(-log[10]~"BH-FDR")
    ) +
    coord_cartesian(clip = "off") +
    labs(title = title, x = NULL, y = NULL) +
    guides(
      fill = guide_legend(order = 1, title.position = "top", title.hjust = 0.5),
      size = guide_legend(order = 2, title.position = "top", title.hjust = 0.5)
    ) +
    theme_void(base_family = "Arial", base_size = font_body_pt) +
    theme(
      plot.title = element_text(
        face = "bold", size = font_title_pt, color = "black",
        hjust = 0.5, margin = margin(b = 6)
      ),
      axis.text.x.bottom = element_text(
        face = "bold", size = font_body_pt, color = "black",
        angle = 35, hjust = 1, vjust = 1, margin = margin(t = 4)
      ),
      legend.position = legend_position,
      legend.box = "vertical",
      legend.box.just = "center",
      legend.title = element_text(face = "bold", size = font_section_pt, color = "black"),
      legend.text = element_text(size = font_small_pt, color = "black"),
      plot.margin = margin(8, 18, 8, 100)
    )
}

write_plot <- function(plot, stem, width, height) {
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 450,
    bg = "white"
  )
}

# Serum: protein candidates are feature-level results, whereas the HDL entry is
# a correlation-adjusted biochemical module. They are deliberately separated.
statistics_object <- readRDS(path_project("fig02", "fig02A", "Fig02A_Full_Statistics.rds"))
serum_functional <- readRDS(path_project(
  "fig02", "Serum_Functional_Analysis",
  "Serum_PreExercise_Functional_Analysis.rds"
))

serum_protein <- statistics_object$Results$Serum_Proteomics %>%
  filter(FeatureLabel %in% c("SHBG", "CUL4B")) %>%
  transmute(
    Block = "Protein\ncandidates",
    Item = FeatureLabel,
    Column = "Proteomics",
    Direction = if_else(Effect < 0, "Higher in lean", "Higher in obese"),
    Effect,
    P_Value,
    BH_FDR,
    EvidenceLevel = if_else(BH_FDR < 0.05, "BH-FDR < 0.05", "Exploratory BH-FDR < 0.20"),
    AnalysisUnit = "Protein feature"
  )

serum_module <- serum_functional$MetaboliteModuleCAMERA %>%
  filter(SetID == "HDL-related lipids", BH_FDR < 0.10) %>%
  transmute(
    Block = "Lipoprotein\nmodule",
    Item = Pathway,
    Column = "Metabolomics",
    Direction,
    Effect = Statistic,
    P_Value,
    BH_FDR,
    EvidenceLevel = "Suggestive BH-FDR < 0.10",
    AnalysisUnit = "Correlation-adjusted biochemical module"
  )

# Expand the serum candidate with representative feature-level NMR results.
# These rows are biochemical measures, not canonical pathway-enrichment hits.
serum_metabolite <- statistics_object$Results$Serum_Metabonomics %>%
  filter(FeatureLabel %in% c(
    "HDL2_C", "HDL_C", "L_HDL_P",
    "VLDL_D", "L_VLDL_P", "M_VLDL_P", "VLDL_TG"
  )) %>%
  mutate(
    Block = case_when(
      FeatureLabel %in% c("HDL2_C", "HDL_C", "L_HDL_P") ~ "HDL-associated\nmeasures",
      TRUE ~ "VLDL-associated\nmeasures"
    ),
    Item = recode(
      FeatureLabel,
      HDL2_C = "Total cholesterol in HDL2",
      HDL_C = "Total cholesterol in HDL",
      L_HDL_P = "Large HDL particle concentration",
      VLDL_D = "Mean VLDL particle diameter",
      L_VLDL_P = "Large VLDL particle concentration",
      M_VLDL_P = "Medium VLDL particle concentration",
      VLDL_TG = "Triglycerides in VLDL"
    )
  ) %>%
  transmute(
    Block,
    Item,
    Column = "Metabolomics",
    Direction = if_else(Effect < 0, "Higher in lean", "Higher in obese"),
    Effect,
    P_Value,
    BH_FDR,
    EvidenceLevel = "Feature-level BH-FDR < 0.10",
    AnalysisUnit = "NMR metabolite or lipoprotein feature"
  ) %>%
  filter(BH_FDR < 0.10)

serum_display <- bind_rows(serum_protein, serum_module)
serum_selection <- tibble::tribble(
  ~Block, ~Item,
  "Protein\ncandidates", "SHBG",
  "Protein\ncandidates", "CUL4B",
  "Lipoprotein\nmodule", "HDL-related lipids"
)

serum_plot <- make_bubble_plot(
  serum_display,
  serum_selection,
  column_levels = c("Proteomics", "Metabolomics"),
  column_labels = c("Proteomics", "Metabolomics"),
  title = "Pre-exercise serum associations with adiposity",
  right_header = "Feature or biochemical module",
  width = 7.6,
  height = 3.0
)
write_plot(serum_plot, "Fig02_Serum_Molecular_Module_Bubble_Candidate", 7.6, 3.0)
readr::write_csv(serum_display, file.path(output_dir, "Fig02_Serum_Molecular_Module_Bubble_SourceData.csv"))

serum_expanded_display <- bind_rows(serum_protein, serum_module, serum_metabolite)
serum_expanded_selection <- tibble::tribble(
  ~Block, ~Item,
  "Protein\ncandidates", "SHBG",
  "Protein\ncandidates", "CUL4B",
  "Coordinated\nbiochemical module", "HDL-related lipids",
  "HDL-associated\nmeasures", "Total cholesterol in HDL2",
  "HDL-associated\nmeasures", "Total cholesterol in HDL",
  "HDL-associated\nmeasures", "Large HDL particle concentration",
  "VLDL-associated\nmeasures", "Mean VLDL particle diameter",
  "VLDL-associated\nmeasures", "Large VLDL particle concentration",
  "VLDL-associated\nmeasures", "Medium VLDL particle concentration",
  "VLDL-associated\nmeasures", "Triglycerides in VLDL"
)

serum_expanded_plot <- make_bubble_plot(
  serum_expanded_display,
  serum_expanded_selection,
  column_levels = c("Proteomics", "Metabolomics"),
  column_labels = c("Proteomics", "Metabolomics"),
  title = "Pre-exercise serum molecular and biochemical associations with adiposity",
  right_header = "Feature or biochemical module",
  width = 8.6,
  height = 4.8
)
write_plot(
  serum_expanded_plot,
  "Fig02_Serum_Molecular_Module_Bubble_Expanded_Candidate",
  8.6,
  4.8
)
readr::write_csv(
  serum_expanded_display,
  file.path(output_dir, "Fig02_Serum_Molecular_Module_Bubble_Expanded_SourceData.csv")
)

adipose_selection <- tibble::tribble(
  ~Block, ~Item,
  "Mitochondrial\nbioenergetics", "Aerobic respiration",
  "Mitochondrial\nbioenergetics", "Mitochondrial respiratory chain complex I assembly",
  "Mitochondrial\nbioenergetics", "Proton motive force-driven mitochondrial ATP synthesis",
  "Mitochondrial\nbioenergetics", "Tricarboxylic acid cycle",
  "Protein synthesis and\nRNA processing", "Cytoplasmic translation",
  "Protein synthesis and\nRNA processing", "Mitochondrial translation",
  "Protein synthesis and\nRNA processing", "mRNA splicing, via spliceosome",
  "Lipid\nmetabolism", "Fatty acid beta-oxidation",
  "Lipid\nmetabolism", "Cholesterol efflux",
  "Immune and\nhemostatic responses", "Complement activation",
  "Immune and\nhemostatic responses", "Acute-phase response",
  "Immune and\nhemostatic responses", "Adaptive immune response",
  "Tissue remodeling\nand signaling", "Cell adhesion mediated by integrin",
  "Tissue remodeling\nand signaling", "Positive regulation of angiogenesis",
  "Tissue remodeling\nand signaling", "Steroid hormone receptor signaling pathway"
)

muscle_selection <- tibble::tribble(
  ~Block, ~Item,
  "Mitochondrial\nbioenergetics", "Proton motive force-driven mitochondrial ATP synthesis",
  "Mitochondrial\nbioenergetics", "Aerobic respiration",
  "Mitochondrial\nbioenergetics", "Mitochondrial electron transport, NADH to ubiquinone",
  "Mitochondrial\nbioenergetics", "Mitochondrial respiratory chain complex I assembly",
  "Mitochondrial\nbioenergetics", "Tricarboxylic acid cycle",
  "Mitochondrial\nbioenergetics", "Proton transmembrane transport",
  "Protein synthesis\nand proteostasis", "Cytoplasmic translation",
  "Protein synthesis\nand proteostasis", "Mitochondrial translation",
  "Protein synthesis\nand proteostasis", "Protein K11-linked ubiquitination",
  "Developmental patterning\n(exploratory DNAm)", "Endocardial cushion formation",
  "Developmental patterning\n(exploratory DNAm)", "Secondary heart field specification",
  "Developmental patterning\n(exploratory DNAm)", "Endoderm formation"
)

read_tissue <- function(tissue) {
  source_file <- if (tissue == "Adipose") {
    path_project("fig02", "SourceData", "Fig02D_Adipose_Pathway_All.xlsx")
  } else {
    path_project("fig02", "SourceData", "FigS02D_Muscle_Pathway_All.xlsx")
  }
  readxl::read_xlsx(source_file, sheet = "Pathway_All") %>%
    mutate(
      ItemSource = as.character(Pathway),
      ItemKey = stringr::str_to_lower(ItemSource),
      Column = recode(
        as.character(Layer),
        "Proteomics" = "Proteomics",
        "Transcriptomics" = "Microarray",
        "DNA methylation" = "DNA methylation"
      ),
      Direction = case_when(
        Direction == "Down" ~ "Higher in lean",
        Direction == "Up" ~ "Higher in obese",
        TRUE ~ NA_character_
      )
    ) %>%
    select(ItemSource, ItemKey, Column, GO_ID, Direction, DirectionStatistic, P_Value, BH_FDR, CoreGenes)
}

prepare_tissue_display <- function(tissue, selection) {
  selection_keyed <- selection %>%
    mutate(ItemKey = stringr::str_to_lower(Item))
  read_tissue(tissue) %>%
    filter(ItemKey %in% selection_keyed$ItemKey, BH_FDR < 0.05) %>%
    left_join(selection_keyed, by = "ItemKey") %>%
    mutate(
      Tissue = tissue,
      AnalysisUnit = "GO Biological Process",
      EvidenceLevel = "BH-FDR < 0.05"
    )
}

adipose_display <- prepare_tissue_display("Adipose", adipose_selection)
muscle_display <- prepare_tissue_display("Muscle", muscle_selection)
tissue_columns <- c("Proteomics", "Microarray", "DNA methylation")
tissue_labels <- c("Proteomics", "Microarray", "DNA methylation")

adipose_plot <- make_bubble_plot(
  adipose_display,
  adipose_selection,
  tissue_columns,
  tissue_labels,
  title = "Pre-exercise adipose pathway associations with adiposity",
  right_header = "GO Biological Process",
  width = 9.0,
  height = 5.8
)

muscle_plot <- make_bubble_plot(
  muscle_display,
  muscle_selection,
  tissue_columns,
  tissue_labels,
  title = "Pre-exercise muscle pathway associations with adiposity",
  right_header = "GO Biological Process",
  width = 9.0,
  height = 5.1
)

write_plot(adipose_plot, "Fig02_Adipose_Pathway_Bubble_Candidate", 9.0, 5.8)
write_plot(muscle_plot, "Fig02_Muscle_Pathway_Bubble_Candidate", 9.0, 5.1)

combined_plot <- adipose_plot + muscle_plot +
  patchwork::plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "right")
write_plot(combined_plot, "Fig02_Tissue_Pathway_Bubbles_Combined_Candidate", 15.0, 5.8)

# Six-column tissue matrix: a pathway occupies one row across both tissues.
# This view is intended to expose cross-tissue overlap without duplicating rows.
tissue_six_selection <- tibble::tribble(
  ~Block, ~Item,
  "Mitochondrial\nbioenergetics", "Aerobic respiration",
  "Mitochondrial\nbioenergetics", "Mitochondrial respiratory chain complex I assembly",
  "Mitochondrial\nbioenergetics", "Mitochondrial electron transport, NADH to ubiquinone",
  "Mitochondrial\nbioenergetics", "Proton motive force-driven mitochondrial ATP synthesis",
  "Mitochondrial\nbioenergetics", "Tricarboxylic acid cycle",
  "Mitochondrial\nbioenergetics", "Proton transmembrane transport",
  "Protein and RNA\nhomeostasis", "Cytoplasmic translation",
  "Protein and RNA\nhomeostasis", "Mitochondrial translation",
  "Protein and RNA\nhomeostasis", "mRNA splicing, via spliceosome",
  "Protein and RNA\nhomeostasis", "Protein K11-linked ubiquitination",
  "Lipid\nmetabolism", "Fatty acid beta-oxidation",
  "Lipid\nmetabolism", "Cholesterol efflux",
  "Immune and\nacute-phase responses", "Complement activation",
  "Immune and\nacute-phase responses", "Acute-phase response",
  "Immune and\nacute-phase responses", "Adaptive immune response",
  "Immune and\nacute-phase responses", "Positive regulation of phagocytosis, engulfment",
  "Immune and\nacute-phase responses", "Chemotaxis",
  "Immune and\nacute-phase responses", "Inflammatory response",
  "Tissue remodeling\nand signaling", "Cell adhesion mediated by integrin",
  "Tissue remodeling\nand signaling", "Positive regulation of angiogenesis",
  "Tissue remodeling\nand signaling", "Steroid hormone receptor signaling pathway"
)

tissue_six_keyed <- tissue_six_selection %>%
  mutate(ItemKey = stringr::str_to_lower(Item))

tissue_six_display <- bind_rows(
  read_tissue("Adipose") %>% mutate(Tissue = "Adipose"),
  read_tissue("Muscle") %>% mutate(Tissue = "Muscle")
) %>%
  filter(ItemKey %in% tissue_six_keyed$ItemKey, BH_FDR < 0.05) %>%
  left_join(tissue_six_keyed, by = "ItemKey") %>%
  mutate(
    Column = paste(Tissue, Column, sep = "__"),
    AnalysisUnit = "GO Biological Process",
    EvidenceLevel = "BH-FDR < 0.05"
  )

tissue_six_columns <- c(
  "Adipose__Proteomics", "Adipose__Microarray", "Adipose__DNA methylation",
  "Muscle__Proteomics", "Muscle__Microarray", "Muscle__DNA methylation"
)
tissue_six_labels <- rep(c("Proteomics", "Microarray", "DNA methylation"), 2)

tissue_six_plot <- make_bubble_plot(
  tissue_six_display,
  tissue_six_selection,
  tissue_six_columns,
  tissue_six_labels,
  title = NULL,
  right_header = "GO Biological Process",
  width = 11.5,
  height = 6.2,
  legend_position = "none",
  column_positions = c(1.00, 1.82, 2.64, 3.72, 4.54, 5.36)
) +
  annotate(
    "text", x = 1.82, y = nrow(tissue_six_selection) + 0.72,
    label = "Adipose tissue", family = "Arial", fontface = "bold",
    size = font_section_pt * pt_to_mm
  ) +
  annotate(
    "text", x = 4.54, y = nrow(tissue_six_selection) + 0.72,
    label = "Skeletal muscle", family = "Arial", fontface = "bold",
    size = font_section_pt * pt_to_mm
  )

write_plot(
  tissue_six_plot,
  "Fig02_Tissue_Pathway_Bubble_SixColumn_Compact_Candidate",
  11.5,
  6.2
)

legend_plot <- ggplot(
  tibble(
    x = 1,
    y = 1,
    Direction = factor(
      c("Higher in lean", "Higher in obese", "Higher in lean", "Higher in obese"),
      levels = names(direction_colors)
    ),
    FDRStrength = c(2, 4, 6, 8)
  ),
  aes(x = x, y = y, fill = Direction, size = FDRStrength)
) +
  geom_point(shape = 21, color = "black", stroke = 0.32, alpha = 0) +
  scale_fill_manual(values = direction_colors, name = "Association direction") +
  scale_size_continuous(
    range = c(2.0, 5.8), limits = c(0, 8), breaks = c(2, 4, 6, 8),
    labels = c("2", "4", "6", ">=8"), name = expression(-log[10]~"BH-FDR")
  ) +
  guides(
    fill = guide_legend(
      order = 1, title.position = "top", title.hjust = 0.5,
      override.aes = list(alpha = 1, size = 3)
    ),
    size = guide_legend(
      order = 2, title.position = "top", title.hjust = 0.5,
      override.aes = list(alpha = 1, fill = "white")
    )
  ) +
  theme_void(base_family = "sans") +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = font_section_pt, color = "black"),
    legend.text = element_text(size = font_small_pt, color = "black")
  )

legend_grob <- ggplotGrob(legend_plot)
legend_grob <- gtable::gtable_filter(legend_grob, "guide-box-right", trim = TRUE)
cairo_pdf(
  file.path(output_dir, "Fig02_Tissue_Pathway_Bubble_Legend.pdf"),
  width = 2.25, height = 2.65
)
grid::grid.newpage()
grid::grid.draw(legend_grob)
dev.off()
png(
  file.path(output_dir, "Fig02_Tissue_Pathway_Bubble_Legend.png"),
  width = 2.25, height = 2.65, units = "in", res = 450,
  type = "cairo", bg = "white"
)
grid::grid.newpage()
grid::grid.draw(legend_grob)
dev.off()
readr::write_csv(
  tissue_six_display,
  file.path(output_dir, "Fig02_Tissue_Pathway_Bubble_SixColumn_SourceData.csv")
)

tissue_all_pathways <- bind_rows(
  read_tissue("Adipose") %>% mutate(Tissue = "Adipose"),
  read_tissue("Muscle") %>% mutate(Tissue = "Muscle")
)
pathway_thresholds <- c(0.05, 0.10, 0.20)
threshold_counts <- bind_rows(lapply(pathway_thresholds, function(q_cutoff) {
  tissue_all_pathways %>%
    group_by(Tissue, Column) %>%
    summarise(
      BH_FDR_Threshold = q_cutoff,
      SignificantPathways = sum(BH_FDR < q_cutoff, na.rm = TRUE),
      .groups = "drop"
    )
}))
threshold_overlap <- bind_rows(lapply(pathway_thresholds, function(q_cutoff) {
  adipose_ids <- tissue_all_pathways %>%
    filter(Tissue == "Adipose", BH_FDR < q_cutoff) %>%
    distinct(GO_ID) %>%
    pull(GO_ID)
  muscle_ids <- tissue_all_pathways %>%
    filter(Tissue == "Muscle", BH_FDR < q_cutoff) %>%
    distinct(GO_ID) %>%
    pull(GO_ID)
  tibble(
    BH_FDR_Threshold = q_cutoff,
    ExactCrossTissueGOOverlap = length(intersect(adipose_ids, muscle_ids))
  )
}))
readr::write_csv(
  threshold_counts,
  file.path(output_dir, "Fig02_Tissue_Pathway_Threshold_Sensitivity.csv")
)
readr::write_csv(
  threshold_overlap,
  file.path(output_dir, "Fig02_Tissue_Pathway_Overlap_Sensitivity.csv")
)

readr::write_csv(adipose_display, file.path(output_dir, "Fig02_Adipose_Pathway_Bubble_SourceData.csv"))
readr::write_csv(muscle_display, file.path(output_dir, "Fig02_Muscle_Pathway_Bubble_SourceData.csv"))

readme <- c(
  "# Figure 2 pathway bubble candidates",
  "",
  paste0("Analysis ID: `", analysis_id, "`."),
  "",
  "- The compact six-column candidate is updated in place; formal assembled Figure 2 artwork is not modified.",
  "- Tissue columns are ordered as Proteomics, Microarray, and DNA methylation.",
  "- The six-column tissue candidate uses one shared pathway row set across adipose and skeletal muscle.",
  "- Its six-column matrix width is compressed while preserving the pathway order, statistics, bubble encoding, and output canvas size.",
  "- Tissue bubbles show GO Biological Process results with BH-FDR < 0.05.",
  "- Developmental muscle DNA-methylation terms are excluded from the compact adiposity-focused candidate but remain in the complete source tables.",
  "- Adipose DNA-methylation terms for phagocytosis, chemotaxis, and inflammatory response are included because each meets the formal two-sided methylglm BH-FDR < 0.05 threshold.",
  "- Directional sensitivity-only methylation findings are not included.",
  "- For transcriptomics and proteomics, color denotes the Obese-minus-Lean pathway direction.",
  "- For DNA methylation, color summarizes the median signed CpG statistic and does not imply gene expression or pathway activation.",
  "- Serum protein entries are feature-level candidates: SHBG meets BH-FDR < 0.05 and CUL4B meets exploratory BH-FDR < 0.20.",
  "- The compact serum metabolomics entry is a correlation-adjusted biochemical module meeting suggestive BH-FDR < 0.10.",
  "- The expanded serum candidate additionally shows selected feature-level HDL and VLDL measures meeting BH-FDR < 0.10; these are not canonical pathway-enrichment results.",
  "- Bubble size represents -log10(BH-FDR)."
)
writeLines(readme, file.path(output_dir, "README.md"))

writeLines(
  c(
    paste0("AnalysisID: ", analysis_id),
    paste0("R version: ", R.version.string),
    paste0("ggplot2: ", packageVersion("ggplot2")),
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  ),
  file.path(output_dir, "RUN_LOG.txt")
)
