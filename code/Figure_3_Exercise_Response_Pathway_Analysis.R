# Final Figure 3D-E Hallmark pathway panels and pathway-member audit tables.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.local_config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  file.path(getwd(), "Fig00_Config.R"),
  if (!is.na(.local_script)) {
    file.path(dirname(normalizePath(.local_script, winslash = "/", mustWork = FALSE)),
              "Fig00_Config.R")
  } else NA_character_
))
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (is.na(.local_config)) stop("Cannot find CellMetabolism_Transfer code/Fig00_Config.R.")
source(.local_config)
rm(.local_script, .local_config_candidates, .local_config)

p_load(dplyr, ggplot2, readxl, stringr, tibble, writexl, msigdbr)

analysis_dir <- path_project("fig03", "Figure3")
input_file <- file.path(analysis_dir, "Fig03_Exercise_Hallmark_Pathway_Results.xlsx")
gene_stats_file <- path_project(
  "fig03", "Figure3",
  "Fig03DE_Proteomics_CAMERA_Source.xlsx"
)
output_dir <- path_project("fig03", "Figure3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Missing Hallmark results: ", input_file)
if (!file.exists(gene_stats_file)) stop("Missing protein statistics: ", gene_stats_file)

all_results <- readxl::read_excel(input_file, sheet = "All_Pathways")
gene_stats <- readxl::read_excel(gene_stats_file, sheet = "Gene_Statistics")

display_map <- tibble::tribble(
  ~Tissue, ~Pathway, ~Block, ~Block_order, ~Pathway_order,
  "Adipose", "Adipogenesis", "Adipocyte and lipid\nmetabolism", 1, 1,
  "Adipose", "Fatty acid metabolism", "Adipocyte and lipid\nmetabolism", 1, 2,
  "Adipose", "Oxidative phosphorylation", "Mitochondrial and\nheme metabolism", 2, 3,
  "Adipose", "Heme metabolism", "Mitochondrial and\nheme metabolism", 2, 4,
  "Muscle", "Oxidative phosphorylation", "Mitochondrial\nenergetics", 1, 1,
  "Muscle", "Coagulation", "Immune and\nhemostatic response", 2, 2,
  "Muscle", "Complement", "Immune and\nhemostatic response", 2, 3,
  "Muscle", "KRAS signaling up", "Growth and tissue\norganization", 3, 4,
  "Muscle", "Apical junction", "Growth and tissue\norganization", 3, 5,
  "Muscle", "Myogenesis", "Myogenic and endocrine\nprograms", 4, 6,
  "Muscle", "Estrogen response early", "Myogenic and endocrine\nprograms", 4, 7
)

size_limits <- c(1.5, 8.0)
size_breaks <- c(2, 4, 6, 8)

plot_data <- all_results %>%
  filter(
    Tissue %in% c("Adipose", "Muscle"),
    Layer %in% c("Proteomics", "DNA methylation"),
    BH_FDR < 0.05
  ) %>%
  inner_join(display_map, by = c("Tissue", "Pathway")) %>%
  mutate(
    Omics = recode(Layer, Proteomics = "Proteomics",
                   `DNA methylation` = "DNA methylation"),
    Evidence = -log10(pmax(BH_FDR, .Machine$double.xmin)),
    Evidence_display = pmin(pmax(Evidence, size_limits[1]), size_limits[2]),
    Evidence_type = case_when(
      Layer == "DNA methylation" ~ "DNA methylation enrichment",
      Direction == "Positive" ~ "Increased after exercise",
      Direction == "Negative" ~ "Decreased after exercise",
      TRUE ~ "Undetermined"
    ),
    Pathway_label = case_when(
      Pathway == "KRAS signaling up" ~ "KRAS signaling (up)",
      Pathway == "Estrogen response early" ~ "Early estrogen response",
      TRUE ~ str_to_sentence(Pathway)
    )
  )

evidence_colors <- c(
  "Decreased after exercise" = "#4F81AD",
  "Increased after exercise" = "#D96B5F",
  "DNA methylation enrichment" = "#8F79B5"
)

make_panel <- function(tissue, title, omics_levels, show_legend, stem,
                       width, height) {
  panel_data <- plot_data %>%
    filter(Tissue == tissue) %>%
    arrange(Pathway_order) %>%
    mutate(
      y = rev(seq_len(n())),
      x = match(Omics, omics_levels)
    )

  n_rows <- nrow(panel_data)
  n_cols <- length(omics_levels)
  blocks <- panel_data %>%
    group_by(Block, Block_order) %>%
    summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5,
              ymid = mean(range(y)), .groups = "drop") %>%
    arrange(Block_order)
  ordered_blocks <- blocks %>% arrange(ymin)
  separators <- if (nrow(ordered_blocks) > 1L) {
    ordered_blocks$ymax[seq_len(nrow(ordered_blocks) - 1L)]
  } else numeric()

  frame_left <- 0.55
  frame_right <- n_cols + 0.45
  label_x <- frame_right + 0.32
  x_max <- if (show_legend) label_x + 3.75 else label_x + 3.30

  p <- ggplot(panel_data, aes(x = x, y = y)) +
    annotate(
      "rect", xmin = frame_left, xmax = frame_right,
      ymin = 0.5, ymax = n_rows + 0.5,
      fill = "white", colour = "black", linewidth = 0.32
    ) +
    geom_segment(
      data = tibble(y = separators),
      aes(x = frame_left, xend = frame_right, y = y, yend = y),
      inherit.aes = FALSE, colour = "black", linewidth = 0.26
    ) +
    geom_segment(
      aes(x = frame_right, xend = frame_right + 0.18, y = y, yend = y),
      colour = "black", linewidth = 0.26
    ) +
    geom_point(
      aes(size = Evidence_display, fill = Evidence_type),
      shape = 21, colour = "black", stroke = 0.30
    ) +
    geom_text(
      aes(x = label_x, label = str_wrap(Pathway_label, 35)),
      hjust = 0, family = "Arial", size = 2.45, colour = "black"
    ) +
    geom_text(
      data = blocks,
      aes(x = frame_left - 0.25, y = ymid, label = Block),
      inherit.aes = FALSE, hjust = 1, family = "Arial", size = 2.45,
      colour = "black", lineheight = 0.92
    ) +
    geom_text(
      data = tibble(x = seq_along(omics_levels), label = omics_levels),
      aes(x = x, y = 0.27, label = label),
      inherit.aes = FALSE, angle = 45, hjust = 1, vjust = 1,
      family = "Arial", fontface = "bold", size = 2.45, colour = "black"
    ) +
    scale_fill_manual(values = evidence_colors, drop = FALSE) +
    scale_size_continuous(
      limits = size_limits, breaks = size_breaks, range = c(2.2, 6.2)
    ) +
    scale_x_continuous(
      breaks = NULL,
      limits = c(frame_left - 1.90, x_max), expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = c(-0.52, n_rows + 0.67), expand = expansion(mult = 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = title, x = NULL, y = NULL,
      fill = "Pathway evidence",
      size = expression(-log[10]("BH-FDR"))
    ) +
    theme_void(base_family = "Arial", base_size = 7.5) +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0,
                                colour = "black", margin = margin(b = 2.2, unit = "mm")),
      axis.text.x = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      legend.justification = "center",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(1.5, "mm"),
      legend.title = element_text(face = "bold", size = 7, colour = "black"),
      legend.text = element_text(size = 6.4, colour = "black"),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(3.5, "mm"),
      plot.margin = margin(3, 4, 10, 4, unit = "mm")
    ) +
    guides(
      fill = guide_legend(order = 1, override.aes = list(size = 3.8)),
      size = guide_legend(order = 2)
    )

  pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
  png_file <- file.path(output_dir, paste0(stem, ".png"))
  ggsave(pdf_file, p, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(png_file, p, width = width, height = height, units = "in", dpi = 500,
         bg = "white")
  invisible(list(pdf = pdf_file, png = png_file))
}

make_panel(
  tissue = "Adipose",
  title = "Adipose Hallmark pathway responses",
  omics_levels = "Proteomics",
  show_legend = FALSE,
  stem = "Fig03D_Adipose_Hallmark_Pathways",
  width = 4.85,
  height = 2.85
)
make_panel(
  tissue = "Muscle",
  title = "Muscle Hallmark pathway responses",
  omics_levels = c("Proteomics", "DNA methylation"),
  show_legend = TRUE,
  stem = "Fig03D_Muscle_Hallmark_Pathways",
  width = 7.10,
  height = 3.80
)

# Trace protein-level drivers of every significant proteomic Hallmark result.
hallmark_sets <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
  transmute(Pathway_key = gs_name, Gene = gene_symbol) %>%
  distinct()

proteomic_pathways <- plot_data %>%
  filter(Layer == "Proteomics") %>%
  mutate(
    Pathway_key = paste0(
      "HALLMARK_", str_replace_all(str_to_upper(Pathway), "[^A-Z0-9]+", "_")
    )
  ) %>%
  select(Tissue, Pathway, Pathway_key, Pathway_BH_FDR = BH_FDR,
         Pathway_direction = Direction)

protein_members <- proteomic_pathways %>%
  inner_join(hallmark_sets, by = "Pathway_key", relationship = "many-to-many") %>%
  inner_join(gene_stats, by = c("Tissue", "Gene"), relationship = "many-to-many") %>%
  group_by(Tissue, Pathway) %>%
  arrange(desc(abs(ModeratedT)), .by_group = TRUE) %>%
  mutate(
    Driver_rank = row_number(),
    Nominal_P_lt_0_05 = FeatureP < 0.05,
    Individual_BH_FDR_lt_0_05 = FeatureBH_FDR < 0.05
  ) %>%
  ungroup()

member_summary <- protein_members %>%
  group_by(Tissue, Pathway, Pathway_direction, Pathway_BH_FDR) %>%
  summarise(
    Measured_proteins = n_distinct(Gene),
    Negative_effect = n_distinct(Gene[Effect < 0]),
    Positive_effect = n_distinct(Gene[Effect > 0]),
    Nominal_P_lt_0_05 = n_distinct(Gene[FeatureP < 0.05]),
    Individual_BH_FDR_lt_0_05 = n_distinct(Gene[FeatureBH_FDR < 0.05]),
    .groups = "drop"
  ) %>%
  arrange(Tissue, Pathway_BH_FDR)

top_drivers <- protein_members %>%
  group_by(Tissue, Pathway) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  arrange(Tissue, Pathway, Driver_rank)

writexl::write_xlsx(
  list(
    Figure_Pathways = plot_data %>%
      arrange(Tissue, Pathway_order) %>%
      select(Tissue, Layer, Block, Pathway, Direction, NGenes,
             P_Value, BH_FDR, Method, Evidence_type),
    Protein_Member_Summary = member_summary,
    Top10_Protein_Drivers = top_drivers,
    All_Protein_Members = protein_members,
    Notes = tibble::tribble(
      ~Field, ~Value,
      "Contrast", "Paired Post3h - Pre exercise response",
      "Model", "Delta(Post3h - Pre) ~ centered adiposity + centered role; FamilyID handled by limma duplicateCorrelation",
      "Collection", "MSigDB Hallmark 2026.1.Hs",
      "Pathway display", "BH-FDR < 0.05 within each tissue-layer analysis",
      "Proteomics effect", "Adjusted log2 relative-abundance change for Post3h - Pre",
      "Protein significance", "Feature-level BH-FDR is reported separately and is not inferred from pathway significance",
      "DNA methylation", "Purple denotes methylglm enrichment without directional inference",
      "Network recommendation", "Use a bipartite pathway-protein membership network; protein fill encodes adjusted log2 change, pathway outline encodes BH-FDR, and edges encode curated Hallmark membership"
    )
  ),
  file.path(output_dir, "Fig03DE_Hallmark_Pathways_and_Protein_Members.xlsx")
)

writeLines(
  c(
    "Figure 3D-E Hallmark BH-FDR final panels",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Hallmark input: ", normalizePath(input_file, winslash = "/")),
    paste0("Protein statistics: ", normalizePath(gene_stats_file, winslash = "/")),
    "Pathway threshold: BH-FDR < 0.05 within each tissue-layer analysis",
    "Protein colours/directions: signed CAMERA results",
    "DNA methylation: unsigned methylglm enrichment shown in neutral purple",
    "Panel D omits the empty adipose DNA-methylation column",
    "Only Panel E carries the shared legend",
    paste0("R: ", R.version.string)
  ),
  file.path(output_dir, "Fig03DE_Hallmark_RunLog.txt")
)

message("Saved final Figure 3D-E files in: ", normalizePath(output_dir, winslash = "/"))
