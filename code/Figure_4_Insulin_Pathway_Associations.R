# Figure 4B candidate panels: pre-exercise, Post3h, and exercise-change
# Hallmark pathway associations with HOMA-IR and Matsuda ISI.

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

p_load(dplyr, ggplot2, ggrepel, openxlsx, patchwork, readr, scales, tidyr)

analysis_id <- "Fig04B_Pre_Post_Delta_HOMA_ISI_Pathway_Comparison_2026-08-28"
input_file <- path_project(
  "fig04", "analysis", "Fig04S04_Complete_Table_Analysis_2026-08-27",
  "molecular_pathway", "Hallmark_GSEA_CAMERA_All.csv"
)
output_dir <- path_project("fig04", "fig04B")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Missing formal pathway table: ", input_file)

state_levels <- c("PreExercise", "Post3hAbsolute", "ExerciseChange")
state_labels <- c(
  PreExercise = "Pre-exercise",
  Post3hAbsolute = "Post3h molecular profile",
  ExerciseChange = "Exercise response (Post3h - Pre)"
)
predictor_levels <- c("HOMA_IR", "Matsuda_ISI")
predictor_labels <- c(HOMA_IR = "HOMA-IR", Matsuda_ISI = "Matsuda ISI")

raw <- read_csv(input_file, show_col_types = FALSE) %>%
  filter(
    State %in% state_levels,
    ModelTier == "Primary",
    Predictor %in% predictor_levels,
    Dataset %in% c("Adipose_Proteomics", "Muscle_Proteomics")
  ) %>%
  mutate(
    State = factor(State, levels = state_levels),
    StateLabel = factor(
      unname(state_labels[as.character(State)]),
      levels = unname(state_labels[state_levels])
    ),
    Predictor = factor(Predictor, levels = predictor_levels),
    PredictorLabel = factor(
      unname(predictor_labels[as.character(Predictor)]),
      levels = unname(predictor_labels[predictor_levels])
    ),
    TissueLabel = recode(
      Dataset,
      Adipose_Proteomics = "Adipose tissue",
      Muscle_Proteomics = "Skeletal muscle"
    )
  )

build_cross_tissue <- function(state, predictor) {
  adipose <- raw %>%
    filter(State == state, Predictor == predictor, Dataset == "Adipose_Proteomics") %>%
    transmute(
      HallmarkID, Pathway,
      NES_Adipose = NES,
      P_Adipose = P_Value,
      BH_FDR_Adipose = BH_FDR,
      CAMERA_BH_FDR_Adipose = CAMERA_BH_FDR,
      Subjects_Adipose = Subjects,
      Families_Adipose = Families
    )

  muscle <- raw %>%
    filter(State == state, Predictor == predictor, Dataset == "Muscle_Proteomics") %>%
    transmute(
      HallmarkID, Pathway,
      NES_Muscle = NES,
      P_Muscle = P_Value,
      BH_FDR_Muscle = BH_FDR,
      CAMERA_BH_FDR_Muscle = CAMERA_BH_FDR,
      Subjects_Muscle = Subjects,
      Families_Muscle = Families
    )

  inner_join(adipose, muscle, by = c("HallmarkID", "Pathway")) %>%
    mutate(
      State = state,
      StateLabel = unname(state_labels[state]),
      Predictor = predictor,
      PredictorLabel = unname(predictor_labels[predictor]),
      Significant_Adipose = BH_FDR_Adipose < 0.05,
      Significant_Muscle = BH_FDR_Muscle < 0.05,
      Significance = case_when(
        Significant_Adipose & Significant_Muscle ~ "Both tissues",
        Significant_Adipose ~ "Adipose only",
        Significant_Muscle ~ "Muscle only",
        TRUE ~ "Neither tissue"
      ),
      Display_BH_FDR = case_when(
        Significant_Adipose & Significant_Muscle ~
          pmax(BH_FDR_Adipose, BH_FDR_Muscle),
        Significant_Adipose ~ BH_FDR_Adipose,
        Significant_Muscle ~ BH_FDR_Muscle,
        TRUE ~ pmin(BH_FDR_Adipose, BH_FDR_Muscle)
      ),
      FDR_Strength = pmin(-log10(pmax(Display_BH_FDR, 1e-10)), 10),
      CAMERA_Concordant =
        CAMERA_BH_FDR_Adipose < 0.05 | CAMERA_BH_FDR_Muscle < 0.05,
      DisplayPathway = recode(
        Pathway,
        `Fatty acid metabolism` = "Fatty-acid metabolism",
        `Epithelial mesenchymal transition` = "Epithelial-mesenchymal transition"
      )
    )
}

cross_tissue_all <- bind_rows(lapply(state_levels, function(state) {
  bind_rows(lapply(predictor_levels, function(predictor) {
    build_cross_tissue(state, predictor)
  }))
})) %>%
  mutate(
    StateLabel = factor(StateLabel, levels = unname(state_labels[state_levels])),
    PredictorLabel = factor(
      PredictorLabel, levels = unname(predictor_labels[predictor_levels])
    ),
    Significance = factor(
      Significance,
      levels = c("Both tissues", "Adipose only", "Muscle only", "Neither tissue")
    )
  )

cross_tissue_sig <- cross_tissue_all %>%
  filter(Significance != "Neither tissue") %>%
  mutate(
    Significance = factor(
      as.character(Significance),
      levels = c("Both tissues", "Adipose only", "Muscle only")
    )
  )

core_pathways <- c(
  "Oxidative phosphorylation",
  "Fatty acid metabolism",
  "Adipogenesis",
  "Heme metabolism",
  "Myogenesis",
  "Coagulation",
  "Complement",
  "Inflammatory response",
  "Glycolysis",
  "Cholesterol homeostasis",
  "Peroxisome",
  "Epithelial mesenchymal transition"
)

label_pathways_2d <- c(
  "Oxidative phosphorylation",
  "Fatty acid metabolism",
  "Adipogenesis",
  "Heme metabolism",
  "Myogenesis",
  "Coagulation",
  "Complement"
)

significance_colors <- c(
  "Both tissues" = "#B7A4D6",
  "Adipose only" = "#E5A1B7",
  "Muscle only" = "#91BBD3"
)

theme_atm <- theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 0.3),
    axis.text = element_text(size = 6.2, colour = "black"),
    axis.title = element_text(size = 7.1, face = "bold", colour = "black"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.margin = margin(2, 2, 2, 2, unit = "mm"),
    legend.position = "right",
    legend.justification = "center",
    legend.box = "vertical",
    legend.title = element_text(size = 6.1, face = "bold"),
    legend.text = element_text(size = 5.7),
    legend.key.height = unit(3.8, "mm")
  )

make_2d_panel <- function(state, predictor, show_y_title = TRUE, show_labels = TRUE) {
  d <- cross_tissue_sig %>%
    filter(State == state, Predictor == predictor)
  label_data <- d %>% filter(Pathway %in% label_pathways_2d)

  p <- ggplot(d, aes(x = NES_Muscle, y = NES_Adipose)) +
    geom_hline(
      yintercept = 0, linetype = "dashed", colour = "#777777", linewidth = 0.25
    ) +
    geom_vline(
      xintercept = 0, linetype = "dashed", colour = "#777777", linewidth = 0.25
    ) +
    geom_point(
      aes(fill = Significance, size = FDR_Strength),
      shape = 21, colour = "black", stroke = 0.25, alpha = 0.95
    ) +
    scale_fill_manual(
      values = significance_colors,
      limits = c("Both tissues", "Adipose only", "Muscle only"),
      breaks = c("Both tissues", "Adipose only", "Muscle only"),
      drop = FALSE,
      name = "GSEA BH-FDR < 0.05"
    ) +
    scale_size_continuous(
      range = c(2.0, 5.2), limits = c(1.3, 10), breaks = c(2, 5, 10),
      name = expression(-log[10]("displayed BH-FDR"))
    ) +
    scale_x_continuous(limits = c(-4, 4), breaks = c(-4, -2, 0, 2, 4)) +
    scale_y_continuous(limits = c(-4, 4), breaks = c(-4, -2, 0, 2, 4)) +
    coord_fixed(clip = "off") +
    labs(
      title = unname(predictor_labels[predictor]),
      x = "Skeletal-muscle GSEA NES",
      y = if (show_y_title) "Adipose-tissue GSEA NES" else NULL
    ) +
    theme_atm

  if (show_labels && nrow(label_data) > 0) {
    p <- p + geom_text_repel(
      data = label_data,
      aes(label = DisplayPathway),
      family = "Arial", size = 2.05, colour = "black",
      min.segment.length = 0, segment.colour = "#777777", segment.size = 0.18,
      box.padding = 0.28, point.padding = 0.28, max.overlaps = Inf,
      seed = 260828, force = 2.0, max.time = 2
    )
  }
  p
}

save_state_2d <- function(state) {
  state_file <- recode(
    state,
    PreExercise = "PreExercise",
    Post3hAbsolute = "Post3hAbsolute",
    ExerciseChange = "ExerciseResponse"
  )
  title_text <- unname(state_labels[state])
  p_homa <- make_2d_panel(state, "HOMA_IR", TRUE, TRUE) +
    theme(legend.position = "none")
  p_isi <- make_2d_panel(state, "Matsuda_ISI", FALSE, TRUE) +
    theme(legend.position = "right", legend.justification = "center")
  combined <- p_homa + p_isi +
    plot_layout(nrow = 1, widths = c(1, 1.28)) +
    plot_annotation(
      title = paste0("Tissue-proteomic pathway associations: ", title_text),
      theme = theme(
        plot.title = element_text(
          family = "Arial", face = "bold", size = 9,
          colour = "black", hjust = 0, margin = margin(0, 0, 2, 0, unit = "mm")
        )
      )
    )

  pdf_file <- file.path(
    output_dir,
    paste0("Fig04B_2D_", state_file, "_HOMA_ISI_Pathway_Landscape.pdf")
  )
  png_file <- sub("\\.pdf$", ".png", pdf_file)
  ggsave(pdf_file, combined, width = 190, height = 91, units = "mm", device = cairo_pdf)
  ggsave(png_file, combined, width = 190, height = 91, units = "mm", dpi = 300)
  c(pdf_file, png_file)
}

state_outputs <- unlist(lapply(state_levels, save_state_2d))

overview_panels <- list()
for (state in state_levels) {
  right_legend <- if (state == "Post3hAbsolute") "right" else "none"
  overview_panels[[length(overview_panels) + 1]] <-
    make_2d_panel(state, "HOMA_IR", TRUE, TRUE) +
    labs(subtitle = unname(state_labels[state])) +
    theme(
      plot.subtitle = element_text(size = 6.6, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  overview_panels[[length(overview_panels) + 1]] <-
    make_2d_panel(state, "Matsuda_ISI", FALSE, TRUE) +
    labs(subtitle = unname(state_labels[state])) +
    theme(
      plot.subtitle = element_text(size = 6.6, face = "bold", hjust = 0.5),
      legend.position = right_legend,
      legend.justification = "center"
    )
}

overview <- wrap_plots(overview_panels, ncol = 2, widths = c(1, 1.28)) +
  plot_annotation(
    title = "Pre-exercise, Post3h, and exercise-response pathway associations",
    theme = theme(
      plot.title = element_text(
        family = "Arial", face = "bold", size = 9,
        colour = "black", hjust = 0, margin = margin(0, 0, 2, 0, unit = "mm")
      )
    )
  )

overview_pdf <- file.path(
  output_dir, "Fig04B_2D_Pre_Post_Delta_HOMA_ISI_Pathway_Overview.pdf"
)
overview_png <- sub("\\.pdf$", ".png", overview_pdf)
ggsave(overview_pdf, overview, width = 190, height = 245, units = "mm", device = cairo_pdf)
ggsave(overview_png, overview, width = 190, height = 245, units = "mm", dpi = 300)

bubble_base <- raw %>%
  mutate(
    StateShort = recode(
      as.character(State),
      PreExercise = "Pre-exercise",
      Post3hAbsolute = "Post3h profile",
      ExerciseChange = "Post3h - Pre"
    ),
    ModelColumn = factor(
      paste(as.character(PredictorLabel), StateShort, sep = "\n"),
      levels = c(
        "HOMA-IR\nPre-exercise", "HOMA-IR\nPost3h profile",
        "HOMA-IR\nPost3h - Pre", "Matsuda ISI\nPre-exercise",
        "Matsuda ISI\nPost3h profile", "Matsuda ISI\nPost3h - Pre"
      )
    ),
    TissueLabel = factor(
      TissueLabel, levels = c("Adipose tissue", "Skeletal muscle")
    ),
    Significant = BH_FDR < 0.05,
    BubbleSize = pmin(-log10(pmax(BH_FDR, 1e-10)), 10)
  )

significant_union <- bubble_base %>%
  filter(Significant) %>%
  distinct(Pathway) %>%
  pull(Pathway)

make_bubble <- function(pathways, title_text, output_stub, height_mm) {
  pathway_order <- pathways
  d <- bubble_base %>%
    filter(Pathway %in% pathway_order) %>%
    mutate(
      Pathway = factor(Pathway, levels = rev(pathway_order)),
      DisplayNES = if_else(Significant, NES, NA_real_),
      DisplaySize = if_else(Significant, BubbleSize, NA_real_)
    )

  p <- ggplot(d, aes(x = ModelColumn, y = Pathway)) +
    geom_point(
      aes(fill = DisplayNES, size = DisplaySize),
      shape = 21, colour = "black", stroke = 0.2, na.rm = TRUE
    ) +
    facet_grid(rows = vars(TissueLabel), scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "#83A9C7", mid = "white", high = "#D9909E",
      midpoint = 0, limits = c(-4, 4), oob = squish,
      name = "GSEA NES"
    ) +
    scale_size_continuous(
      range = c(1.6, 5.3), limits = c(1.3, 10), breaks = c(2, 5, 10),
      name = expression(-log[10]("BH-FDR"))
    ) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
      panel.grid.major.x = element_line(colour = "#E2E2E2", linewidth = 0.22),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
      axis.text.x = element_text(size = 5.8, colour = "black", angle = 35, hjust = 1),
      axis.text.y = element_text(size = 5.9, colour = "black"),
      axis.ticks = element_blank(),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.35),
      strip.text = element_text(size = 6.5, face = "bold", colour = "black"),
      plot.title = element_text(size = 8.5, face = "bold", hjust = 0),
      legend.position = "right",
      legend.justification = "center",
      legend.title = element_text(size = 6.1, face = "bold"),
      legend.text = element_text(size = 5.7),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )

  pdf_file <- file.path(output_dir, paste0(output_stub, ".pdf"))
  png_file <- sub("\\.pdf$", ".png", pdf_file)
  ggsave(pdf_file, p, width = 190, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(png_file, p, width = 190, height = height_mm, units = "mm", dpi = 300)
  c(pdf_file, png_file)
}

core_order <- core_pathways[core_pathways %in% significant_union]
bubble_core_outputs <- make_bubble(
  core_order,
  "Core Hallmark pathway associations across molecular states",
  "Fig04B_Bubble_Core_Hallmark_Pre_Post_Delta_HOMA_ISI",
  112
)
bubble_all_outputs <- make_bubble(
  sort(significant_union),
  "Hallmark pathways with BH-FDR < 0.05 in at least one model",
  "Fig04B_Bubble_AllSignificant_Hallmark_Pre_Post_Delta_HOMA_ISI",
  190
)

# Tiered-FDR candidates requested for the main and supplementary figures.
main_core_pathways <- c(
  "Oxidative phosphorylation",
  "Fatty acid metabolism",
  "Adipogenesis",
  "Myogenesis",
  "Heme metabolism",
  "Coagulation",
  "Complement",
  "Inflammatory response"
)

tiered_base <- bubble_base %>%
  mutate(
    StateAxis = factor(
      StateShort,
      levels = c("Pre-exercise", "Post3h profile", "Post3h - Pre"),
      labels = c("Pre-exercise", "Post3h", "Response")
    ),
    FDRTier = case_when(
      BH_FDR < 0.05 ~ "***",
      BH_FDR < 0.10 ~ "**",
      BH_FDR < 0.20 ~ "*",
      TRUE ~ ""
    ),
    DisplayNES20 = if_else(BH_FDR < 0.20, NES, NA_real_),
    DisplaySize20 = if_else(
      BH_FDR < 0.20,
      pmin(-log10(pmax(BH_FDR, 1e-10)), 10),
      NA_real_
    )
  )

make_tiered_bubble <- function(
    pathways, title_text, output_stub, width_mm, height_mm,
    tissue_layout = c("rows", "columns")) {
  tissue_layout <- match.arg(tissue_layout)
  pathway_order <- pathways
  d <- tiered_base %>%
    filter(Pathway %in% pathway_order) %>%
    mutate(Pathway = factor(Pathway, levels = rev(pathway_order)))

  tissue_facet <- if (tissue_layout == "rows") {
    facet_grid(
      rows = vars(TissueLabel), cols = vars(PredictorLabel),
      scales = "free_y", space = "free_y"
    )
  } else {
    facet_grid(cols = vars(TissueLabel, PredictorLabel))
  }

  p <- ggplot(d, aes(x = StateAxis, y = Pathway)) +
    geom_point(
      aes(fill = DisplayNES20, size = DisplaySize20),
      shape = 21, colour = "#333333", stroke = 0.22, na.rm = TRUE
    ) +
    geom_text(
      data = d %>% filter(FDRTier != ""),
      aes(label = FDRTier),
      family = "Arial", fontface = "bold", colour = "black",
      size = 1.85, vjust = 0.42, na.rm = TRUE
    ) +
    tissue_facet +
    scale_fill_gradient2(
      low = "#83A9C7", mid = "white", high = "#D9909E",
      midpoint = 0, limits = c(-4, 4), oob = squish,
      name = "GSEA NES"
    ) +
    scale_size_continuous(
      range = c(2.3, 6.1), limits = c(0.7, 10), breaks = c(1, 2, 5, 10),
      name = expression(-log[10]("BH-FDR"))
    ) +
    guides(
      fill = guide_colourbar(
        title.position = "top", title.hjust = 0.5,
        barheight = unit(21, "mm"), barwidth = unit(3.5, "mm")
      ),
      size = guide_legend(
        title.position = "top", title.hjust = 0.5,
        keyheight = unit(3.4, "mm"), keywidth = unit(5.5, "mm")
      )
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = NULL,
      caption = paste0(
        "Stars denote BH-FDR tiers: ***BH-FDR < 0.05, ",
        "**BH-FDR < 0.10, and *BH-FDR < 0.20."
      )
    ) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.28),
      panel.spacing.x = unit(0.7, "mm"),
      panel.spacing.y = unit(0.7, "mm"),
      axis.text.x = element_text(size = 5.5, colour = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 6.1, colour = "black"),
      axis.ticks = element_blank(),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.28),
      strip.text = element_text(size = 6.1, face = "bold", colour = "black"),
      plot.title = element_text(size = 8.5, face = "bold", hjust = 0),
      plot.caption = element_text(size = 5.7, hjust = 0, colour = "black"),
      legend.position = "right",
      legend.justification = "center",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.spacing.y = unit(1.2, "mm"),
      legend.title = element_text(size = 5.4, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 5.0),
      legend.key.height = unit(3.4, "mm"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )

  pdf_file <- file.path(output_dir, paste0(output_stub, ".pdf"))
  png_file <- sub("\\.pdf$", ".png", pdf_file)
  ggsave(pdf_file, p, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(png_file, p, width = width_mm, height = height_mm, units = "mm", dpi = 300)
  c(pdf_file, png_file)
}

fdr20_union <- tiered_base %>%
  filter(BH_FDR < 0.20) %>%
  distinct(Pathway) %>%
  arrange(Pathway) %>%
  pull(Pathway)

tiered_main_outputs <- make_tiered_bubble(
  main_core_pathways,
  "Core Hallmark pathway associations across molecular states",
  "Fig04B_Main_Core_Hallmark_Pre_Post_Delta_FDR_Tiers",
  125,
  93,
  "rows"
)

tiered_supp_outputs <- make_tiered_bubble(
  fdr20_union,
  "Hallmark pathways with BH-FDR < 0.20 in at least one model",
  "FigS04B_All_Hallmark_Pre_Post_Delta_FDR_Tiers",
  190,
  150,
  "columns"
)

supplement_file <- file.path(
  output_dir, "Table_Fig04B_Hallmark_Pathway_Associations.xlsx"
)

format_table_data <- function(data) {
  data %>%
    transmute(
      Tissue = as.character(TissueLabel),
      Predictor = as.character(PredictorLabel),
      Molecular_state = as.character(StateLabel),
      Hallmark_pathway = Pathway,
      NES = NES,
      P_value = P_Value,
      BH_FDR = BH_FDR,
      FDR_tier = case_when(
        BH_FDR < 0.05 ~ "BH-FDR < 0.05",
        BH_FDR < 0.10 ~ "0.05 <= BH-FDR < 0.10",
        BH_FDR < 0.20 ~ "0.10 <= BH-FDR < 0.20",
        TRUE ~ "BH-FDR >= 0.20"
      ),
      CAMERA_BH_FDR = CAMERA_BH_FDR,
      Effective_n = Subjects,
      Families = Families
    ) %>%
    arrange(Tissue, Predictor, Molecular_state, BH_FDR, Hallmark_pathway)
}

core_table <- raw %>%
  filter(Pathway %in% main_core_pathways) %>%
  format_table_data()
all_table <- format_table_data(raw)

counts <- raw %>%
  group_by(StateLabel, TissueLabel, PredictorLabel) %>%
  summarise(
    Subjects = first(Subjects),
    Families = first(Families),
    TestedPathways = n(),
    GSEA_FDR05 = sum(BH_FDR < 0.05, na.rm = TRUE),
    GSEA_FDR10 = sum(BH_FDR < 0.10, na.rm = TRUE),
    GSEA_FDR20 = sum(BH_FDR < 0.20, na.rm = TRUE),
    CAMERA_FDR05 = sum(CAMERA_BH_FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

model_summary <- counts %>%
  transmute(
    Molecular_state = as.character(StateLabel),
    Tissue = as.character(TissueLabel),
    Predictor = as.character(PredictorLabel),
    Effective_n = Subjects,
    Families = Families,
    Hallmark_pathways_tested = TestedPathways,
    GSEA_BH_FDR_lt_0.05 = GSEA_FDR05,
    GSEA_BH_FDR_lt_0.10 = GSEA_FDR10,
    GSEA_BH_FDR_lt_0.20 = GSEA_FDR20,
    CAMERA_BH_FDR_lt_0.05 = CAMERA_FDR05
  ) %>%
  arrange(Predictor, Molecular_state, Tissue)

wb <- createWorkbook(creator = "WYC")
header_style <- createStyle(
  fontName = "Arial", fontSize = 9, textDecoration = "bold",
  halign = "center", valign = "center",
  border = "Bottom", borderStyle = "medium"
)
body_style <- createStyle(fontName = "Arial", fontSize = 9, valign = "center")
bottom_style <- createStyle(
  fontName = "Arial", fontSize = 9, valign = "center",
  border = "Bottom", borderStyle = "medium"
)
scientific_style <- createStyle(fontName = "Arial", fontSize = 9, numFmt = "0.00E+00")
decimal_style <- createStyle(fontName = "Arial", fontSize = 9, numFmt = "0.000")

write_three_line_sheet <- function(sheet, data) {
  addWorksheet(wb, sheet, gridLines = FALSE)
  writeData(wb, sheet, data, headerStyle = header_style)
  if (nrow(data) > 0L) {
    addStyle(wb, sheet, body_style, rows = 2:(nrow(data) + 1L),
             cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
    addStyle(wb, sheet, bottom_style, rows = nrow(data) + 1L,
             cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
  }
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = "auto")
  setColWidths(wb, sheet, cols = which(names(data) == "Hallmark_pathway"), widths = 34)
  setRowHeights(wb, sheet, rows = 1, heights = 24)
  if ("NES" %in% names(data)) {
    addStyle(wb, sheet, decimal_style, rows = 2:(nrow(data) + 1L),
             cols = which(names(data) == "NES"), gridExpand = TRUE, stack = TRUE)
  }
  p_cols <- which(names(data) %in% c("P_value", "BH_FDR", "CAMERA_BH_FDR"))
  if (length(p_cols) > 0L) {
    addStyle(wb, sheet, scientific_style, rows = 2:(nrow(data) + 1L),
             cols = p_cols, gridExpand = TRUE, stack = TRUE)
  }
}

write_three_line_sheet("Core pathways", core_table)
write_three_line_sheet("All Hallmark", all_table)
write_three_line_sheet("Model summary", model_summary)

notes_start <- nrow(model_summary) + 4L
notes <- data.frame(
  Item = c(
    "Primary model", "Family correlation", "Pathway method", "Multiple testing",
    "Molecular states", "NES interpretation", "Figure symbols"
  ),
  Description = c(
    "Protein abundance or Post3h-minus-Pre response ~ rank-normalized baseline HOMA-IR or Matsuda ISI + family role.",
    "FamilyID was handled using limma duplicateCorrelation.",
    "Hallmark GSEA used the limma moderated t-statistic; CAMERA was retained as a correlation-aware sensitivity analysis.",
    "Benjamini-Hochberg correction was applied within each tissue, molecular state, predictor, and pathway method.",
    "Pre-exercise and Post3h are absolute molecular profiles; exercise response is Post3h minus Pre.",
    "Positive NES indicates enrichment toward positive phenotype-protein associations; it does not by itself establish pathway activation.",
    "*** BH-FDR < 0.05; ** 0.05 <= BH-FDR < 0.10; * 0.10 <= BH-FDR < 0.20."
  ),
  stringsAsFactors = FALSE
)
writeData(wb, "Model summary", notes, startRow = notes_start, headerStyle = header_style)
addStyle(wb, "Model summary", body_style,
         rows = (notes_start + 1L):(notes_start + nrow(notes)), cols = 1:2,
         gridExpand = TRUE, stack = TRUE)
addStyle(wb, "Model summary", bottom_style,
         rows = notes_start + nrow(notes), cols = 1:2,
         gridExpand = TRUE, stack = TRUE)
setColWidths(wb, "Model summary", cols = 1, widths = 24)
setColWidths(wb, "Model summary", cols = 2, widths = 105)
setRowHeights(wb, "Model summary", rows = (notes_start + 1L):(notes_start + nrow(notes)), heights = 30)
addStyle(
  wb, "Model summary", createStyle(wrapText = TRUE, valign = "top"),
  rows = (notes_start + 1L):(notes_start + nrow(notes)), cols = 1:2,
  gridExpand = TRUE, stack = TRUE
)
saveWorkbook(wb, supplement_file, overwrite = TRUE)

source_file <- file.path(
  output_dir, "Fig04B_Pre_Post_Delta_HOMA_ISI_Pathway_SourceData.csv"
)
counts_file <- file.path(
  output_dir, "Fig04B_Pre_Post_Delta_HOMA_ISI_Pathway_Counts.csv"
)
results_file <- file.path(
  output_dir, "Fig04B_Pre_Post_Delta_HOMA_ISI_Pathway_Results.md"
)
log_file <- file.path(output_dir, "run_log_pre_post_delta.txt")

write_csv(cross_tissue_all, source_file, na = "")

write_csv(counts, counts_file)

results_lines <- c(
  "# Figure 4B pathway comparison",
  "",
  "## Model",
  "",
  "Primary role-adjusted tissue-proteomic models used rank-normalized baseline HOMA-IR or Matsuda ISI, with FamilyID handled by limma duplicateCorrelation. Hallmark GSEA used the moderated t-statistic; CAMERA was retained as a correlation-aware sensitivity analysis.",
  "",
  "## Interpretation",
  "",
  "Pre-exercise oxidative-phosphorylation and fatty-acid-metabolism rankings were negative for HOMA-IR and positive for Matsuda ISI in both adipose tissue and skeletal muscle. Exercise-response models identified coordinated pathway-level associations despite no individual protein reaching BH-FDR < 0.20. Post3h absolute-profile associations were strongest for Matsuda ISI but may retain substantial pre-exercise information and therefore should not be interpreted as exercise-induced changes.",
  "",
  "## Figure use",
  "",
  "The two-dimensional panels emphasize cross-tissue concordance and tissue specificity. The bubble matrices provide the complete state-by-predictor comparison. Blank cells in the tiered-FDR figures indicate BH-FDR >= 0.20, not missing measurements.",
  "",
  "## Candidate pathways",
  "",
  "Oxidative phosphorylation; fatty-acid metabolism; adipogenesis; heme metabolism; myogenesis; coagulation; complement; inflammatory response; glycolysis; cholesterol homeostasis; peroxisome; epithelial-mesenchymal transition.",
  "",
  "## Molecular boundary",
  "",
  "GLIPR2 and CD276 were associated with Matsuda ISI in Post3h adipose abundance (BH-FDR = 0.0437), whereas no exercise-change protein reached BH-FDR < 0.20. GSDME was associated with HOMA-IR in the body-fat-adjusted pre-exercise muscle model (BH-FDR = 0.0309). These features should be shown with their Pre, Post3h, and Post3h-minus-Pre coefficients to distinguish absolute-profile associations from exercise-response associations."
)
writeLines(results_lines, results_file)

log_lines <- c(
  paste0("Analysis ID: ", analysis_id),
  paste0("Input: ", input_file),
  "Models: primary role-adjusted Hallmark GSEA; FamilyID duplicateCorrelation.",
  "Two-dimensional panels include pathways with GSEA BH-FDR < 0.05 in at least one tissue.",
  "Strict bubble plots show points where GSEA BH-FDR < 0.05; tiered-FDR plots show points where BH-FDR < 0.20 and distinguish the <0.05, <0.10, and <0.20 tiers.",
  paste0("Outputs: ", paste(c(
    state_outputs, overview_pdf, overview_png,
    bubble_core_outputs, bubble_all_outputs,
    tiered_main_outputs, tiered_supp_outputs,
    supplement_file, source_file, counts_file, results_file
  ), collapse = "; ")),
  "",
  capture.output(sessionInfo())
)
writeLines(log_lines, log_file)

message("Completed Figure 4B pathway comparison outputs in: ", output_dir)
