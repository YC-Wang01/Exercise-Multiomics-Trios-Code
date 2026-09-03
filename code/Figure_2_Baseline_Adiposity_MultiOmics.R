# Figure 2 Ed2: modality-specific BH-FDR thresholds adapted from the acute
# exercise reference figure. This script reuses the locked Figure 2 statistics
# and writes only to fig02.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.local_config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  file.path(getwd(), "Fig00_Config.R"),
  if (!is.na(.local_script)) {
    file.path(dirname(normalizePath(.local_script, winslash = "/", mustWork = FALSE)), "Fig00_Config.R")
  } else {
    NA_character_
  }
))
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (is.na(.local_config)) stop("Cannot find code/Fig00_Config.R.")
source(.local_config)
rm(.local_script, .local_config_candidates, .local_config)

p_load(cowplot, digest, dplyr, ggplot2, ggrepel, patchwork, png, ragg,
       readr, tibble, tidyr, writexl)
select <- dplyr::select
set.seed(20260822)

analysis_id <- "Fig02_Ed2_NominalPVolcano_12FeatureBar_2026-08-22"
nominal_p_threshold <- 0.05
bar_feature_total <- 12L
bar_platform_cap <- 3L
volcano_y_cap <- 8
output_dir <- path_project("fig02")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

statistics_file <- path_project("fig02", "fig02A", "Fig02A_Full_Statistics.rds")
annotation_file <- path_project(
  "fig02", "fig02A", "Fig02A_Methylation_Label_Annotation.csv.gz"
)
statistics <- readRDS(statistics_file)

expected_model <- "Feature ~ Adiposity + Role; FamilyID duplicateCorrelation block"
if (!identical(statistics$Model, expected_model)) {
  stop("The locked Figure 2 model does not match the requested categorical model.")
}

# The reference figure does not include DNA methylation. A conventional
# within-platform BH-FDR <0.05 is therefore prespecified for methylation.
registry <- tibble::tribble(
  ~Dataset, ~Tissue, ~Omics, ~PanelTitle, ~Order, ~FDRThreshold, ~XAxisTitle, ~XLimit, ~XBreaks,
  "Serum_Metabonomics", "Serum", "Metabolomics", "Serum metabolomics", 1L, 0.10,
  "Adjusted difference\n(Obese - Lean)", 6.5, list(c(-6, -3, 0, 3, 6)),
  "Adipose_Proteomics", "Adipose", "Proteomics", "Adipose proteomics", 2L, 0.20,
  "Adjusted log2FC\n(Obese - Lean)", 3.4, list(c(-3, -1.5, 0, 1.5, 3)),
  "Adipose_Microarray", "Adipose", "Transcriptomics", "Adipose transcriptomics", 3L, 0.05,
  "Adjusted log2FC\n(Obese - Lean)", 3.4, list(c(-3, -1.5, 0, 1.5, 3)),
  "Adipose_Methylation", "Adipose", "DNA methylation", "Adipose DNA methylation", 4L, 0.05,
  "Adjusted M-value\n(Obese - Lean)", 3.4, list(c(-3, -1.5, 0, 1.5, 3)),
  "Serum_Proteomics", "Serum", "Proteomics", "Serum proteomics", 5L, 0.20,
  "Adjusted log2FC\n(Obese - Lean)", 1.05, list(c(-0.8, -0.4, 0, 0.4, 0.8)),
  "Muscle_Proteomics", "Muscle", "Proteomics", "Muscle proteomics", 6L, 0.20,
  "Adjusted log2FC\n(Obese - Lean)", 2.4, list(c(-2, -1, 0, 1, 2)),
  "Muscle_Microarray", "Muscle", "Transcriptomics", "Muscle transcriptomics", 7L, 0.05,
  "Adjusted log2FC\n(Obese - Lean)", 2.4, list(c(-2, -1, 0, 1, 2)),
  "Muscle_Methylation", "Muscle", "DNA methylation", "Muscle DNA methylation", 8L, 0.05,
  "Adjusted M-value\n(Obese - Lean)", 2.4, list(c(-2, -1, 0, 1, 2))
)

direction_colors <- c(
  "Not nominally significant" = "#D7DEE2",
  "Nominal P, higher in Lean" = "#3FA0A8",
  "Nominal P, higher in Obese" = "#F59646"
)

platform_colors <- c(
  "Serum metabolomics" = "#C99768",
  "Adipose proteomics" = "#D7A34B",
  "Adipose transcriptomics" = "#78A982",
  "Adipose DNA methylation" = "#68A6A8",
  "Serum proteomics" = "#D9827B",
  "Muscle proteomics" = "#6D9BB2",
  "Muscle transcriptomics" = "#718CB7",
  "Muscle DNA methylation" = "#9A7EAA"
)

methylation_annotation <- readr::read_csv(
  annotation_file, show_col_types = FALSE, progress = FALSE
) %>%
  select(Dataset, FeatureID, GeneLabel, UCSC_RefGene_Group) %>%
  distinct(Dataset, FeatureID, .keep_all = TRUE)

prepare_dataset <- function(dataset) {
  spec <- registry %>% filter(Dataset == dataset)
  data <- tibble::as_tibble(statistics$Results[[dataset]]) %>%
    transmute(
      Dataset = dataset,
      Tissue = spec$Tissue[[1]],
      Omics = spec$Omics[[1]],
      PanelTitle = spec$PanelTitle[[1]],
      PlatformOrder = spec$Order[[1]],
      FDRThreshold = spec$FDRThreshold[[1]],
      FeatureID = as.character(FeatureID),
      OriginalLabel = as.character(FeatureLabel),
      Effect = as.numeric(Effect),
      ModeratedT = as.numeric(ModeratedT),
      P_Value = as.numeric(P_Value),
      BH_FDR = as.numeric(BH_FDR),
      Samples = as.integer(Samples),
      Families = as.integer(Families),
      LeanN = as.integer(LeanN),
      ObeseN = as.integer(ObeseN),
      FamilyCorrelation = as.numeric(FamilyCorrelation)
    ) %>%
    filter(is.finite(Effect), is.finite(P_Value), is.finite(BH_FDR))

  if (grepl("Methylation$", dataset)) {
    data <- data %>%
      left_join(
        methylation_annotation %>% filter(Dataset == dataset),
        by = c("Dataset", "FeatureID")
      ) %>%
      mutate(
        DisplayLabel = if_else(!is.na(GeneLabel) & nzchar(GeneLabel), GeneLabel, FeatureID),
        AnnotatedLabel = !is.na(GeneLabel) & nzchar(GeneLabel)
      )
  } else {
    data <- data %>%
      mutate(
        GeneLabel = NA_character_,
        UCSC_RefGene_Group = NA_character_,
        DisplayLabel = if_else(
          !is.na(OriginalLabel) & nzchar(OriginalLabel), OriginalLabel, FeatureID
        ),
        DisplayLabel = if (dataset == "Serum_Metabonomics") DisplayLabel else toupper(DisplayLabel),
        AnnotatedLabel = TRUE
      )
  }

  data %>%
    mutate(
      PlotLabel = iconv(DisplayLabel, from = "", to = "ASCII//TRANSLIT", sub = ""),
      FDRPass = BH_FDR < FDRThreshold,
      Direction = if_else(Effect < 0, "Higher in Lean", "Higher in Obese"),
      PlotY = pmin(-log10(pmax(P_Value, 1e-300)), volcano_y_cap),
      DisplayClass = case_when(
        P_Value < nominal_p_threshold & Effect < 0 ~ "Nominal P, higher in Lean",
        P_Value < nominal_p_threshold & Effect > 0 ~ "Nominal P, higher in Obese",
        TRUE ~ "Not nominally significant"
      ),
      DisplayClass = factor(DisplayClass, levels = names(direction_colors)),
      SignedEvidence = sign(Effect) * -log10(pmax(BH_FDR, 1e-300))
    )
}

prepared <- setNames(lapply(registry$Dataset, prepare_dataset), registry$Dataset)
all_results <- bind_rows(prepared)

summarise_platforms <- function(metric = c("P", "FDR")) {
  metric <- match.arg(metric)
  bind_rows(lapply(registry$Dataset, function(dataset) {
    data <- prepared[[dataset]]
    pass <- if (metric == "P") data$P_Value < nominal_p_threshold else data$FDRPass
    tibble(
      Dataset = dataset,
      FeaturesTested = nrow(data),
      Threshold = if (metric == "P") nominal_p_threshold else unique(data$FDRThreshold),
      HigherLean = sum(pass & data$Effect < 0, na.rm = TRUE),
      HigherObese = sum(pass & data$Effect > 0, na.rm = TRUE)
    )
  })) %>%
    left_join(registry %>% select(Dataset, Tissue, Omics, PanelTitle, Order), by = "Dataset") %>%
    mutate(
      HigherLeanPct = 100 * HigherLean / FeaturesTested,
      HigherObesePct = 100 * HigherObese / FeaturesTested,
      Total = HigherLean + HigherObese
    )
}

p_summary <- summarise_platforms("P")
fdr_summary <- summarise_platforms("FDR")
fdr_positive <- fdr_summary %>% filter(Total > 0)

make_landscape <- function(summary, x_title) {
  long <- summary %>%
    pivot_longer(
      c(HigherLeanPct, HigherObesePct),
      names_to = "DirectionCode", values_to = "Percent"
    ) %>%
    mutate(
      Direction = if_else(DirectionCode == "HigherLeanPct", "Higher in Lean", "Higher in Obese"),
      PlotValue = if_else(Direction == "Higher in Lean", -Percent, Percent),
      PanelTitle = factor(PanelTitle, levels = rev(summary$PanelTitle[order(summary$Order)]))
    )
  limit <- max(abs(long$PlotValue), na.rm = TRUE) * 1.12
  if (!is.finite(limit) || limit <= 0) limit <- 1
  ggplot(long, aes(PlotValue, PanelTitle, fill = Direction)) +
    geom_col(width = 0.58, color = "black", linewidth = 0.28) +
    geom_vline(xintercept = 0, linewidth = 0.40, color = "black") +
    scale_fill_manual(
      values = c("Higher in Lean" = "#3FA0A8", "Higher in Obese" = "#F59646"),
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(-limit, limit), labels = function(x) abs(x),
      breaks = scales::pretty_breaks(n = 5), expand = expansion(mult = 0)
    ) +
    labs(x = x_title, y = NULL) +
    theme_classic(base_family = "Arial", base_size = 6) +
    theme(
      axis.title.x = element_text(face = "bold", size = 6.2, color = "black"),
      axis.text = element_text(size = 5.8, color = "black"),
      legend.text = element_text(size = 5.8, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.25),
      panel.grid = element_blank(),
      plot.margin = margin(3, 3, 3, 3)
    )
}

p_landscape <- make_landscape(
  p_summary,
  "Proportion of tested features with nominal P < 0.05 (%)"
) +
  labs(title = "Pre-exercise multi-omic associations with adiposity") +
  theme(
    plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
    legend.position = "none"
  )

fdr_inset <- make_landscape(
  fdr_positive,
  "Features meeting modality-specific\nBH-FDR threshold (%)"
) +
  labs(title = "BH-FDR-supported platforms") +
  theme(
    plot.title = element_text(face = "bold", size = 6.3, hjust = 0.5),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA)
  )

shared_legend <- cowplot::get_legend(
  make_landscape(p_summary, "") +
    theme(legend.position = "bottom", legend.direction = "horizontal")
)

right_column <- cowplot::plot_grid(
  NULL, fdr_inset, shared_legend, NULL,
  ncol = 1, rel_heights = c(0.12, 0.70, 0.14, 0.04)
)

panel_a <- cowplot::plot_grid(
  p_landscape, right_column,
  nrow = 1, rel_widths = c(0.70, 0.30), align = "h", axis = "tb"
)

panel_a_file <- file.path(
  output_dir, "Fig02A_EightOmics_NominalP_with_ModalityBH_FDR_Inset.pdf"
)
ggsave(panel_a_file, panel_a, width = 190, height = 122, units = "mm", device = cairo_pdf)

select_volcano_labels <- function(data) {
  label_data <- data %>%
    filter(
      P_Value < nominal_p_threshold,
      AnnotatedLabel, !is.na(PlotLabel), nzchar(PlotLabel)
    )
  if (grepl("Microarray$", unique(data$Dataset))) {
    label_data <- label_data %>%
      arrange(P_Value, BH_FDR, desc(abs(Effect))) %>%
      distinct(PlotLabel, .keep_all = TRUE)
  }
  label_data %>%
    group_by(Direction) %>%
    arrange(P_Value, BH_FDR, desc(abs(Effect)), .by_group = TRUE) %>%
    slice_head(n = 2L) %>%
    ungroup()
}

make_volcano <- function(dataset, show_y_title = FALSE) {
  data <- prepared[[dataset]]
  labels <- select_volcano_labels(data)
  spec <- registry %>% filter(Dataset == dataset)
  x_limit <- spec$XLimit[[1]]
  x_breaks <- unlist(spec$XBreaks[[1]], use.names = FALSE)

  raster_file <- tempfile(paste0(dataset, "_points_"), fileext = ".png")
  ragg::agg_png(
    raster_file, width = 1400, height = 1400, units = "px", res = 600,
    background = "transparent"
  )
  par(mar = rep(0, 4), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(-x_limit, x_limit), ylim = c(0, volcano_y_cap))
  for (display_class in levels(data$DisplayClass)) {
    sub <- data[data$DisplayClass == display_class, , drop = FALSE]
    if (!nrow(sub)) next
    nominal_class <- display_class != "Not nominally significant"
    point_color <- grDevices::adjustcolor(
      direction_colors[[display_class]], alpha.f = if (nominal_class) 0.88 else 0.46
    )
    graphics::points(
      sub$Effect, sub$PlotY, pch = 16,
      cex = if (nominal_class) 0.28 else 0.20, col = point_color
    )
  }
  grDevices::dev.off()
  raster_points <- png::readPNG(raster_file)
  unlink(raster_file)

  legend_data <- tibble(
    Effect = 0, PlotY = 0,
    DisplayClass = factor(names(direction_colors), levels = names(direction_colors))
  )

  ggplot() +
    annotation_raster(
      raster_points, xmin = -x_limit, xmax = x_limit,
      ymin = 0, ymax = volcano_y_cap, interpolate = FALSE
    ) +
    geom_hline(
      yintercept = -log10(nominal_p_threshold), linetype = "22",
      linewidth = 0.28, color = "#5F6A70"
    ) +
    geom_vline(xintercept = 0, linewidth = 0.28, color = "#4C5960") +
    geom_point(
      data = legend_data, aes(Effect, PlotY, color = DisplayClass),
      alpha = 0, size = 0.1, show.legend = TRUE
    ) +
    ggrepel::geom_text_repel(
      data = labels, aes(Effect, PlotY, label = PlotLabel),
      family = "sans", size = 1.7, color = "black",
      segment.color = "#68747C", segment.size = 0.20,
      min.segment.length = 0, box.padding = 0.16, point.padding = 0.08,
      max.overlaps = Inf, seed = 20260822
    ) +
    scale_color_manual(values = direction_colors, name = NULL, drop = FALSE) +
    scale_x_continuous(
      limits = c(-x_limit, x_limit), breaks = x_breaks,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = c(0, volcano_y_cap), breaks = c(0, 2, 4, 6, 8),
      expand = expansion(mult = c(0, 0.01))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = spec$PanelTitle[[1]], x = spec$XAxisTitle[[1]],
      y = if (show_y_title) expression(-log[10]~"nominal P") else NULL
    ) +
    theme_classic(base_family = "sans", base_size = 6) +
    theme(
      plot.title = element_text(face = "bold", size = 6.7, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 5.7, color = "black"),
      axis.text = element_text(size = 5.3, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
      axis.line = element_blank(), panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.23),
      aspect.ratio = 1, legend.position = "none",
      plot.margin = margin(3, 3, 3, 3)
    )
}

volcano_plots <- lapply(seq_len(nrow(registry)), function(index) {
  make_volcano(registry$Dataset[[index]], show_y_title = index %in% c(1L, 5L))
})

volcano_legend <- cowplot::get_legend(
  ggplot(
    tibble(
      x = seq_along(direction_colors), y = 1,
      DisplayClass = factor(names(direction_colors), levels = names(direction_colors))
    ),
    aes(x, y, color = DisplayClass)
  ) +
    geom_point(size = 2.2) +
    scale_color_manual(values = direction_colors, name = NULL) +
    theme_void(base_family = "sans", base_size = 6) +
    theme(
      legend.position = "bottom", legend.direction = "horizontal",
      legend.text = element_text(size = 5.5, color = "black")
    )
)

volcano_grid <- patchwork::wrap_plots(volcano_plots, ncol = 4, nrow = 2)
panel_b <- cowplot::plot_grid(
  volcano_grid,
  volcano_legend,
  ncol = 1,
  rel_heights = c(0.92, 0.08)
)

panel_b_file <- file.path(output_dir, "Fig02B_EightOmics_NominalP_Volcanoes.pdf")
pdf_device <- function(filename, width, height, ...) {
  grDevices::pdf(
    file = filename, width = width, height = height,
    family = "Helvetica", useDingbats = FALSE
  )
}
ggsave(
  panel_b_file, panel_b, width = 190, height = 104, units = "mm",
  device = pdf_device
)

bar_pool <- bind_rows(lapply(registry$Dataset, function(dataset) {
  data <- prepared[[dataset]] %>% filter(FDRPass)
  if (grepl("Microarray$|Methylation$", dataset)) {
    data <- data %>%
      arrange(BH_FDR, P_Value, desc(abs(Effect))) %>%
      distinct(PlotLabel, .keep_all = TRUE)
  }
  if (grepl("Methylation$", dataset)) {
    data <- data %>% filter(AnnotatedLabel)
  }
  data
}))

bar_allocation <- bar_pool %>%
  count(Dataset, PanelTitle, PlatformOrder, name = "CandidateCount") %>%
  arrange(CandidateCount, PlatformOrder) %>%
  mutate(
    Quota = if_else(CandidateCount <= 2L, CandidateCount, 1L),
    MaximumQuota = pmin(CandidateCount, bar_platform_cap)
  )

remaining_slots <- bar_feature_total - sum(bar_allocation$Quota)
while (remaining_slots > 0L) {
  eligible <- which(bar_allocation$Quota < bar_allocation$MaximumQuota)
  if (!length(eligible)) break
  for (index in eligible) {
    if (remaining_slots <= 0L) break
    bar_allocation$Quota[[index]] <- bar_allocation$Quota[[index]] + 1L
    remaining_slots <- remaining_slots - 1L
  }
}
if (sum(bar_allocation$Quota) != bar_feature_total) {
  stop("Unable to allocate exactly ", bar_feature_total, " Figure 2C features.")
}

select_balanced_bar_features <- function(data, quota) {
  ranked <- data %>% arrange(BH_FDR, P_Value, desc(abs(Effect)))
  if (quota <= 1L) return(slice_head(ranked, n = quota))

  direction_seeds <- ranked %>%
    group_by(Direction) %>%
    slice_head(n = 1L) %>%
    ungroup() %>%
    arrange(BH_FDR, P_Value, desc(abs(Effect)))
  if (nrow(direction_seeds) >= quota) return(slice_head(direction_seeds, n = quota))

  seed_ids <- direction_seeds$FeatureID
  bind_rows(
    direction_seeds,
    ranked %>% filter(!FeatureID %in% seed_ids) %>%
      slice_head(n = quota - nrow(direction_seeds))
  )
}

bar_candidates <- bind_rows(lapply(seq_len(nrow(bar_allocation)), function(index) {
  dataset <- bar_allocation$Dataset[[index]]
  quota <- bar_allocation$Quota[[index]]
  select_balanced_bar_features(bar_pool %>% filter(Dataset == dataset), quota)
})) %>%
  mutate(
    Platform = factor(PanelTitle, levels = names(platform_colors)),
    DirectionOrder = if_else(Direction == "Higher in Lean", 1L, 2L),
    RowID = paste(Dataset, Direction, FeatureID, sep = "__")
  ) %>%
  arrange(PlatformOrder, DirectionOrder, BH_FDR) %>%
  mutate(RowID = factor(RowID, levels = rev(unique(RowID))))

bar_labels <- setNames(bar_candidates$PlotLabel, as.character(bar_candidates$RowID))
bar_limit <- ceiling(max(abs(bar_candidates$SignedEvidence), na.rm = TRUE))
bar_limit <- max(bar_limit, 3L)
bar_breaks <- seq(-bar_limit, bar_limit, by = 1)

panel_c <- ggplot(
  bar_candidates,
  aes(x = SignedEvidence, y = RowID, fill = Platform)
) +
  geom_vline(xintercept = 0, linewidth = 0.32, color = "#4E5A60") +
  geom_col(width = 0.64, color = "black", linewidth = 0.27) +
  scale_fill_manual(values = platform_colors, drop = TRUE, name = "Omics platform") +
  scale_x_continuous(
    limits = c(-bar_limit, bar_limit), breaks = bar_breaks,
    expand = expansion(mult = c(0.015, 0.015))
  ) +
  scale_y_discrete(labels = bar_labels) +
  labs(
    title = "Strongest pre-exercise adiposity associations meeting modality-specific BH-FDR thresholds",
    subtitle = "Higher in Lean  <-  12 features; scarce platforms prioritized; maximum three per platform  ->  Higher in Obese",
    x = expression("Signed statistical evidence (" * -log[10] * " BH-FDR)"),
    y = NULL,
    caption = paste0(
      "Thresholds: transcripts and DNA methylation <0.05; metabolites <0.10; proteins <0.20. ",
      "BH-FDR was calculated separately within each omics platform."
    )
  ) +
  theme_classic(base_family = "Arial", base_size = 7) +
  theme(
    plot.title = element_text(face = "bold", size = 8.5, hjust = 0.5),
    plot.subtitle = element_text(size = 6.1, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 6.8, color = "black"),
    axis.text.x = element_text(size = 5.8, color = "black"),
    axis.text.y = element_text(size = 5.8, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.36),
    axis.line = element_blank(), panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.24),
    legend.position = "right", legend.justification = "center",
    legend.title = element_text(face = "bold", size = 6.2),
    legend.text = element_text(size = 5.7),
    legend.key.height = grid::unit(3.7, "mm"),
    legend.key.width = grid::unit(4.0, "mm"),
    plot.caption = element_text(size = 5.2, hjust = 0.5, color = "#343A3D"),
    plot.margin = margin(4, 5, 4, 5)
  )

panel_c_height <- max(92, 7.0 * nrow(bar_candidates))
panel_c_file <- file.path(output_dir, "Fig02C_ModalityThreshold_FDR_Features_Bar.pdf")
ggsave(
  panel_c_file, panel_c, width = 180, height = panel_c_height,
  units = "mm", device = cairo_pdf, limitsize = FALSE
)

# A compact assembly candidate is supplied for placement review. Individual
# panels remain the authoritative editing inputs.
composite <- patchwork::wrap_elements(panel_a) /
  patchwork::wrap_elements(panel_c) /
  patchwork::wrap_elements(panel_b) +
  patchwork::plot_layout(heights = c(1.02, 1.08, 1.00))
composite_file <- file.path(output_dir, "Fig02_Ed2_ABC_Layout_Candidate.pdf")
ggsave(
  composite_file, composite, width = 190, height = 270,
  units = "mm", device = cairo_pdf, limitsize = FALSE
)

label_table <- bind_rows(lapply(registry$Dataset, function(dataset) {
  select_volcano_labels(prepared[[dataset]]) %>%
    transmute(
      Dataset, PanelTitle, FeatureID, DisplayLabel, GeneLabel,
      Effect, Direction, P_Value, BH_FDR, FDRThreshold,
      Samples, Families, LeanN, ObeseN
    )
}))

parameter_table <- tibble(
  Field = c(
    "Analysis ID", "Source statistics", "Source statistics SHA256", "Model",
    "Adiposity definition", "Contrast", "Nominal P threshold",
    "Transcript threshold", "Protein threshold", "Metabolite threshold",
    "Complex-lipid threshold", "DNA methylation threshold", "FDR scope",
    "Volcano y axis", "Volcano horizontal line", "Volcano point color", "Bar selection"
  ),
  Value = c(
    analysis_id, statistics_file,
    digest::digest(file = statistics_file, algo = "sha256"),
    statistics$Model,
    "Obese if body-fat percentage >30%; Lean otherwise",
    "Obese minus Lean",
    "P <0.05",
    "Within-platform BH-FDR <0.05",
    "Within-platform BH-FDR <0.20",
    "Within-platform BH-FDR <0.10",
    "Reference threshold is BH-FDR <0.20; no separate complex-lipid platform is present in the locked eight-omics object",
    "Within-platform BH-FDR <0.05; prespecified because the reference figure did not include DNA methylation",
    "Benjamini-Hochberg correction separately within each omics platform",
    "-log10(nominal P), capped at 8 for display only",
    "Nominal P =0.05",
    "Nominal P <0.05 with direction-specific colors; BH-FDR does not alter volcano color",
    "Exactly 12 threshold-positive features; platforms with only one or two candidates retained first; maximum three per platform; directions balanced where available; transcript and methylation labels deduplicated by gene"
  )
)

source_workbook <- file.path(output_dir, "Fig02_Ed2_SourceData.xlsx")
writexl::write_xlsx(
  list(
    PanelA_NominalP = p_summary,
    PanelA_ModalityFDR = fdr_summary,
    PanelB_Labels = label_table,
    PanelC_Selected = bar_candidates %>%
      mutate(Platform = as.character(Platform), RowID = as.character(RowID)) %>%
      select(
        Dataset, PanelTitle, Tissue, Omics, FeatureID, DisplayLabel,
        GeneLabel, UCSC_RefGene_Group, Effect, Direction, P_Value, BH_FDR,
        FDRThreshold, SignedEvidence, Samples, Families, LeanN, ObeseN
      ),
    PanelC_Allocation = bar_allocation %>%
      arrange(PlatformOrder) %>%
      select(Dataset, PanelTitle, CandidateCount, Quota, MaximumQuota),
    ThresholdPositive_All = all_results %>%
      filter(FDRPass) %>%
      select(
        Dataset, PanelTitle, Tissue, Omics, FeatureID, DisplayLabel,
        GeneLabel, UCSC_RefGene_Group, Effect, ModeratedT, Direction,
        P_Value, BH_FDR, FDRThreshold, Samples, Families, LeanN, ObeseN,
        FamilyCorrelation
      ),
    Parameters = parameter_table
  ),
  source_workbook
)

readme_file <- file.path(output_dir, "README.md")
writeLines(c(
  "# Figure 2 Ed2",
  "",
  "This directory was generated independently from the current fig02 release. No current Figure 2 PDF, PNG, AI, table, or final-submission file was overwritten.",
  "",
  "## Panels",
  "",
  "- `Fig02A_EightOmics_NominalP_with_ModalityBH_FDR_Inset.pdf`: nominal-P landscape with a modality-specific BH-FDR inset.",
  "- `Fig02B_EightOmics_NominalP_Volcanoes.pdf`: eight nominal-P volcano plots in the requested 2 x 4 order.",
  "- `Fig02C_ModalityThreshold_FDR_Features_Bar.pdf`: exactly 12 threshold-positive features, prioritizing platforms with only one or two candidates and limiting each platform to at most three.",
  "- `Fig02_Ed2_ABC_Layout_Candidate.pdf`: vertical layout candidate; use individual PDFs for final Illustrator assembly.",
  "",
  "## Thresholds",
  "",
  "- Transcriptomics: within-platform BH-FDR <0.05.",
  "- Proteomics: within-platform BH-FDR <0.20.",
  "- Serum metabolomics: within-platform BH-FDR <0.10.",
  "- DNA methylation: within-platform BH-FDR <0.05 because the reference figure did not define a methylation threshold.",
  "- The locked analysis contains no separate complex-lipid platform; the reference BH-FDR <0.20 complex-lipid threshold was therefore not applied to selected serum metabolite features.",
  "",
  "## Interpretation",
  "",
  "Volcano height and point color are based on nominal P; teal and orange denote P <0.05 with higher abundance in Lean and Obese participants, respectively. BH-FDR does not alter volcano color. Panel C contains exactly 12 threshold-positive features and no filler features."
), readme_file)

log_file <- file.path(output_dir, "RUN_LOG.txt")
writeLines(c(
  paste0("AnalysisID=", analysis_id),
  paste0("Completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("SourceStatistics=", statistics_file),
  paste0("Model=", statistics$Model),
  "FDRScope=Benjamini-Hochberg separately within each omics platform",
  "Thresholds=transcripts 0.05; proteins 0.20; metabolites 0.10; DNA methylation 0.05",
  "VolcanoColor=Nominal P <0.05 and effect direction only; BH-FDR does not alter volcano color",
  "PanelCSelection=Exactly 12 threshold-positive features; scarce platforms retained first; maximum three per platform; directions balanced where available",
  "ComplexLipids=reference threshold 0.20 not applied because no separate complex-lipid platform is present",
  paste0("PanelA=", panel_a_file),
  paste0("PanelB=", panel_b_file),
  paste0("PanelC=", panel_c_file),
  paste0("Composite=", composite_file),
  paste0("SourceWorkbook=", source_workbook),
  capture.output(sessionInfo())
), log_file)

message("Completed Figure 2 Ed2 outputs in: ", output_dir)
