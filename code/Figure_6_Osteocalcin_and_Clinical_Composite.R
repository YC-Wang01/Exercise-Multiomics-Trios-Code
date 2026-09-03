# Figure 6A-B in the reviewed one-row geometry, updated only by adding the
# fasting osteocalcin measurement to the left trajectory panel.

.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.this_script <- if (length(.file_arg)) sub("^--file=", "", .file_arg[[1]]) else NA_character_
.script_dir <- if (!is.na(.this_script)) dirname(normalizePath(.this_script)) else file.path(getwd(), "code")
.config_file <- file.path(.script_dir, "Fig00_Config.R")
if (!file.exists(.config_file)) stop("Cannot locate Fig00_Config.R.")
source(.config_file)
rm(.file_arg, .this_script, .config_file)

p_load(
  dplyr, emmeans, ggplot2, lme4, lmerTest, patchwork, readr, readxl,
  scales, tibble, tidyr, writexl
)

analysis_id <- "Fig06AB_FourTimepoint_OriginalLayout_2026-08-31"
main_dir <- path_project("fig06", "fig06AB")
analysis_dir <- path_project("fig06", "analysis", analysis_id)
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

four_timepoint_rds <- path_project(
  "Pre_Post_Response_Analysis", "Results",
  "08_Osteocalcin_Molecular_State_and_Response",
  "05_Four_Timepoint_Longitudinal_Analysis",
  "OC_Four_Timepoint_Longitudinal_Analysis.rds"
)
clinical_results_file <- path_project(
  "fig06", "analysis", "Fig06B_Pre_Post_Clinical_Associations_Candidate_2026-08-30",
  "Fig06B_Pre_Post_Clinical_Associations_Source_Data.xlsx"
)
if (!all(file.exists(c(four_timepoint_rds, clinical_results_file)))) stop("Required input is missing.")

markers <- c(
  "cOC", "TotalOC", "cOC_ratio", "BCAA", "Serotonin", "Leptin",
  "Serotonin_Leptin_Ratio"
)
marker_labels <- c(
  cOC = "Carboxylated osteocalcin (cOC)",
  TotalOC = "Total osteocalcin (tOC)",
  cOC_ratio = "cOC/tOC ratio",
  BCAA = "Branched-chain amino acids",
  Serotonin = "Serotonin",
  Leptin = "Leptin",
  Serotonin_Leptin_Ratio = "Serotonin/leptin ratio"
)
marker_colors <- c(cOC = "#5AA79D", TotalOC = "#D47E70", cOC_ratio = "#8E88B9")
marker_short_labels <- c(cOC = "cOC", TotalOC = "tOC", cOC_ratio = "cOC/tOC ratio")
time_levels <- c("Fast", "Pre", "Post1h", "Post3h")
time_labels <- c("Fast", "Pre", "1 h", "3 h")

four_timepoint <- readRDS(four_timepoint_rds)
all_plot_data <- four_timepoint$plot_data %>%
  transmute(
    Participant = as.character(PlotSubject),
    Marker = as.character(Marker),
    Timepoint = factor(as.character(Timepoint), levels = time_levels, labels = time_labels),
    RelativeToPre = as.numeric(RelativeToPre),
    AbsoluteValue = as.numeric(AbsoluteValue)
  ) %>%
  group_by(Participant, Marker) %>%
  mutate(
    FastValue = AbsoluteValue[Timepoint == "Fast"][1],
    AbsoluteChangeFromFast = AbsoluteValue - FastValue,
    RelativeToFast = AbsoluteValue / FastValue
  ) %>%
  ungroup()
if (n_distinct(all_plot_data$Participant) != 50L) stop("Expected n = 50 in the four-timepoint source data.")

# Display-only sensitivity exclusions. The mixed-model estimates and BH-FDR
# annotations below remain the complete-cohort (n = 50) primary results.
absolute_display_outlier <- all_plot_data %>%
  filter(Marker %in% c("TotalOC", "cOC")) %>%
  arrange(desc(AbsoluteValue)) %>%
  slice(1) %>%
  pull(Participant)

relative_display_outlier <- all_plot_data %>%
  filter(Marker == "cOC", Timepoint %in% c("1 h", "3 h")) %>%
  arrange(desc(RelativeToPre)) %>%
  slice(1) %>%
  pull(Participant)

absolute_plot_data <- all_plot_data %>%
  filter(Participant != absolute_display_outlier)
relative_plot_data <- all_plot_data %>%
  filter(Participant != relative_display_outlier)

if (n_distinct(absolute_plot_data$Participant) != 49L) stop("Expected n = 49 in the absolute display.")
if (n_distinct(relative_plot_data$Participant) != 49L) stop("Expected n = 49 in the relative display.")

relative_effects <- four_timepoint$time_effects %>%
  filter(Contrast %in% c("Pre_minus_Fast", "Post1h_minus_Pre", "Post3h_minus_Pre")) %>%
  transmute(
    Marker = as.character(Marker),
    Contrast = as.character(Contrast),
    BH_FDR = as.numeric(BH_FDR_WithinFamily)
  )

format_fdr <- function(x) {
  if (x >= 0.01) return(formatC(x, format = "f", digits = 3))
  formatC(x, format = "e", digits = 1)
}

fixed_nice_breaks <- function(limits, n = 4L) {
  limits <- range(limits, finite = TRUE)
  if (diff(limits) <= 0) return(rep(limits[[1]], n))
  raw_step <- diff(limits) / (n - 1L)
  magnitude <- 10^floor(log10(raw_step))
  candidate_steps <- c(1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 7.5, 10) * magnitude
  for (step in candidate_steps) {
    start <- floor(limits[[1]] / step) * step
    candidate <- start + seq.int(0, n - 1L) * step
    if (candidate[[1]] <= limits[[1]] && candidate[[n]] >= limits[[2]]) {
      return(signif(candidate, 8))
    }
  }
  step <- 10 * magnitude
  start <- floor(limits[[1]] / step) * step
  signif(start + seq.int(0, n - 1L) * step, 8)
}

absolute_shared_breaks <- fixed_nice_breaks(
  c(0, 1.12 * max(
    absolute_plot_data$AbsoluteValue[absolute_plot_data$Marker %in% c("TotalOC", "cOC")],
    na.rm = TRUE
  )),
  n = 4L
)
absolute_ratio_breaks <- fixed_nice_breaks(
  c(0, 1.04 * max(
    absolute_plot_data$AbsoluteValue[absolute_plot_data$Marker == "cOC_ratio"],
    na.rm = TRUE
  )),
  n = 4L
)
relative_shared_breaks <- c(0.5, 1.0, 1.5, 2.0)
relative_ratio_breaks <- c(0.5, 1.0, 1.5, 2.0)

trajectory_theme <- theme_classic(base_family = "Arial", base_size = 6.5) +
  theme(
    plot.title = element_text(size = 7.0, face = "bold", hjust = 0.5, margin = margin(b = 1.5)),
    axis.title.y = element_text(size = 6.2, face = "bold", margin = margin(r = 2)),
    axis.text.x = element_text(size = 5.2, color = "black"),
    axis.text.y = element_text(size = 5.2, color = "black"),
    axis.ticks = element_line(linewidth = 0.28, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.38),
    plot.margin = margin(2.2, 2.0, 1.7, 2.0)
  )

make_trajectory_panel <- function(
  data, effects, marker_name, relative = TRUE, show_y = FALSE,
  axis_breaks
) {
  value_column <- if (relative) "RelativeToFast" else "AbsoluteValue"
  dat <- data %>% filter(Marker == marker_name)
  axis_limits <- range(axis_breaks)
  axis_accuracy <- if (max(abs(axis_breaks)) <= 2) 0.1 else 1
  stat <- effects %>% filter(Marker == marker_name) %>%
    mutate(
      DisplayTime = recode(
        Contrast,
        Pre_minus_Fast = 2,
        Post1h_minus_Pre = 3,
        Post3h_minus_Pre = 4
      )
    )
  ann <- tibble(
    x = stat$DisplayTime,
    y = axis_limits[[2]] - 0.03 * diff(axis_limits),
    Label = paste0("FDR\n", vapply(stat$BH_FDR, format_fdr, character(1)))
  )
  y_title <- if (relative) {
    if (show_y) "Relative level (Fast = 1)" else NULL
  } else if (marker_name == "cOC_ratio") {
    "cOC/tOC ratio"
  } else if (show_y) {
    "Concentration (ng/mL)"
  } else NULL
  ggplot(dat, aes(x = Timepoint, y = .data[[value_column]], group = Participant)) +
    {if (relative) geom_hline(yintercept = 1, linewidth = 0.24, color = "#6E6E6E", linetype = "dashed") else NULL} +
    geom_vline(xintercept = 2.5, linewidth = 0.22, color = "#777777", linetype = "dashed") +
    geom_line(linewidth = 0.18, alpha = 0.25, color = "#909090") +
    geom_point(size = 0.48, alpha = 0.34, color = "#858585") +
    stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 0.72, color = marker_colors[[marker_name]]) +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.55, stroke = 0.42, fill = "white", color = marker_colors[[marker_name]]) +
    geom_text(
      data = ann, aes(x = x, y = y, label = Label), inherit.aes = FALSE,
      family = "Arial", size = 1.45, lineheight = 0.86, color = "black",
      vjust = 1
    ) +
    scale_y_continuous(
      limits = axis_limits, breaks = axis_breaks,
      labels = label_number(accuracy = axis_accuracy),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(title = unname(marker_short_labels[[marker_name]]), x = NULL, y = y_title) +
    trajectory_theme
}

assemble_trajectory <- function(data, effects, relative, title) {
  shared_breaks <- if (relative) relative_shared_breaks else absolute_shared_breaks
  ratio_breaks <- if (relative) relative_ratio_breaks else absolute_ratio_breaks
  make_trajectory_panel(data, effects, "cOC", relative, TRUE, shared_breaks) +
    make_trajectory_panel(data, effects, "TotalOC", relative, FALSE, shared_breaks) +
    make_trajectory_panel(data, effects, "cOC_ratio", relative, FALSE, ratio_breaks) +
    plot_layout(widths = c(1.05, 1, 1)) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(family = "Arial", size = 8.0, face = "bold", hjust = 0.5, margin = margin(b = 0.8)),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm")
      )
    )
}

relative_trajectory <- assemble_trajectory(
  relative_plot_data, relative_effects, TRUE,
  "Relative osteocalcin profiles following acute exercise"
)
absolute_trajectory <- assemble_trajectory(
  absolute_plot_data, relative_effects, FALSE,
  "Effect of acute exercise on circulating osteocalcin profiles"
)

clinical_display <- readxl::read_xlsx(clinical_results_file, sheet = "Displayed_results") %>%
  mutate(
    MarkerShort = factor(MarkerShort, levels = c("cOC", "tOC", "cOC/tOC")),
    DomainDisplay = factor(DomainDisplay, levels = c("Bone resorption", "Adiposity", "Liver enzymes"))
  )

phenotype_order <- c(
  "TRACP5b", "Alkaline phosphatase", "Body-fat percentage", "Total fat mass",
  "Arm fat mass", "Trunk fat mass", "Android fat mass", "Subcutaneous fat mass",
  "Visceral fat mass", "ALT"
)
clinical_display <- clinical_display %>%
  mutate(
    PhenotypeDisplay = factor(PhenotypeDisplay, levels = rev(phenotype_order)),
    Significance = case_when(
      BH_FDR_WithinStateMarker < 0.001 ~ "***",
      BH_FDR_WithinStateMarker < 0.01 ~ "**",
      BH_FDR_WithinStateMarker < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

effect_limit <- max(0.8, ceiling(max(abs(c(clinical_display$CI_Low, clinical_display$CI_High)), na.rm = TRUE) * 10) / 10)
effect_limit <- min(effect_limit, 1.4)
effect_breaks <- seq(-effect_limit, effect_limit, length.out = 3L)

clinical_theme <- theme_classic(base_size = 6.1, base_family = "Arial") +
  theme(
    axis.text = element_text(color = "black", size = 5.3),
    axis.title = element_text(color = "black", size = 5.8, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.23),
    axis.ticks = element_line(color = "black", linewidth = 0.23),
    strip.background = element_blank(),
    panel.spacing.y = grid::unit(0.35, "mm")
  )

make_clinical_figure <- function(state_name) {
  dat <- clinical_display %>% filter(State == state_name)
  state_title <- if (state_name == "Pre") {
    "Associations of osteocalcin with clinical traits (pre-exercise)"
  } else {
    "Associations of osteocalcin with clinical traits (Post 3 h)"
  }

  heatmap <- ggplot(dat, aes(MarkerShort, PhenotypeDisplay, fill = StandardizedEffect)) +
    geom_tile(color = "black", linewidth = 0.16, width = 0.86, height = 0.86) +
    geom_text(aes(label = Significance), color = "black", family = "Arial", size = 2.05, fontface = "bold") +
    facet_grid(rows = vars(DomainDisplay), scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2(
      low = "#3E6E9E", mid = "#FFFFFF", high = "#C8564B", midpoint = 0,
      limits = c(-effect_limit, effect_limit), oob = squish, guide = "none"
    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    labs(title = "Association matrix", x = NULL, y = NULL) +
    clinical_theme +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, color = "black", size = 5.0, face = "bold"),
      axis.text.x.top = element_text(color = "black", size = 5.4, face = "bold"),
      axis.text.y = element_text(color = "black", size = 4.9),
      axis.ticks = element_blank(), axis.line = element_blank(),
      plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5, margin = margin(b = 1)),
      plot.margin = margin(1.2, 1.0, 1.2, 1.0, unit = "mm")
    )

  forest <- ggplot(dat, aes(StandardizedEffect, PhenotypeDisplay, color = MarkerShort)) +
    geom_vline(xintercept = 0, linewidth = 0.20, linetype = "dashed", color = "black") +
    geom_errorbarh(aes(xmin = CI_Low, xmax = CI_High), height = 0, linewidth = 0.28, color = "black") +
    geom_point(size = 1.35) +
    facet_grid(rows = vars(DomainDisplay), cols = vars(MarkerShort), scales = "free_y", space = "free_y") +
    scale_color_manual(values = c(tOC = "#D47E70", cOC = "#5AA79D", `cOC/tOC` = "#8E88B9"), guide = "none") +
    scale_x_continuous(
      limits = c(-effect_limit, effect_limit), breaks = effect_breaks,
      labels = label_number(accuracy = 0.1), expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "Standardized effect", y = NULL) +
    clinical_theme +
    theme(
      strip.text.x = element_text(color = "black", size = 5.2, face = "bold"),
      strip.text.y = element_blank(), axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), axis.line.y = element_blank(),
      axis.text.x = element_text(color = "black", size = 4.8),
      axis.title.x = element_text(color = "black", size = 5.4, face = "bold"),
      plot.margin = margin(1.2, 0.8, 1.2, 0.8, unit = "mm")
    )

  # A high-resolution interpolated raster avoids the hairline seams that PDF
  # viewers can show between adjacent vector tiles. Width and height use fixed
  # physical units so 6B and S6C retain the same 2 x 15 mm colourbar.
  gradient_colours <- grDevices::colorRampPalette(
    c("#C8564B", "#FFFFFF", "#3E6E9E")
  )(512)
  gradient_raster <- as.raster(matrix(gradient_colours, ncol = 1))
  bar_x <- grid::unit(0.30, "npc")
  bar_y <- grid::unit(0.62, "npc")

  legend_grob <- grid::grobTree(
    grid::textGrob(
      "Effect", x = bar_x, y = bar_y + grid::unit(10.5, "mm"),
      gp = grid::gpar(fontfamily = "Arial", fontsize = 5.3, fontface = "bold")
    ),
    grid::rasterGrob(
      gradient_raster, x = bar_x, y = bar_y,
      width = grid::unit(2, "mm"), height = grid::unit(15, "mm"),
      interpolate = TRUE
    ),
    grid::textGrob(
      "+1.0", x = bar_x + grid::unit(2.4, "mm"),
      y = bar_y + grid::unit(7.5, "mm"), just = c("left", "center"),
      gp = grid::gpar(fontfamily = "Arial", fontsize = 4.7)
    ),
    grid::textGrob(
      "0", x = bar_x + grid::unit(2.4, "mm"), y = bar_y,
      just = c("left", "center"),
      gp = grid::gpar(fontfamily = "Arial", fontsize = 4.7)
    ),
    grid::textGrob(
      "-1.0", x = bar_x + grid::unit(2.4, "mm"),
      y = bar_y - grid::unit(7.5, "mm"), just = c("left", "center"),
      gp = grid::gpar(fontfamily = "Arial", fontsize = 4.7)
    ),
    grid::textGrob(
      "BH-FDR\n* < 0.05\n** < 0.01\n*** < 0.001",
      x = grid::unit(0.12, "npc"), y = bar_y - grid::unit(11.0, "mm"),
      just = c("left", "top"),
      gp = grid::gpar(fontfamily = "Arial", fontsize = 4.45, lineheight = 1.08)
    )
  )
  legend <- wrap_elements(full = legend_grob)

  heatmap + forest + legend +
    plot_layout(widths = c(2.35, 4.3, 0.7)) +
    plot_annotation(
      title = state_title,
      theme = theme(
        plot.title = element_text(family = "Arial", color = "black", size = 8.0, face = "bold", hjust = 0.02),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm")
      )
    )
}

pre_clinical <- make_clinical_figure("Pre")

compact_left_panel <- function(plot_object) {
  plot_spacer() /
    wrap_elements(full = plot_object) /
    plot_spacer() +
    plot_layout(heights = c(0.07, 0.86, 0.07))
}

relative_trajectory_compact <- compact_left_panel(relative_trajectory)
absolute_trajectory_compact <- compact_left_panel(absolute_trajectory)

main_one_row <- wrap_elements(full = absolute_trajectory_compact) +
  wrap_elements(full = pre_clinical) +
  plot_layout(widths = c(1.00, 1.00)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 9, color = "black"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm")
    )
  )

main_pdf <- file.path(
  main_dir,
  "Fig06AB_Absolute_OC_Trajectories_and_PreExercise_Clinical_Associations_190mm.pdf"
)

ggsave(main_pdf, main_one_row, width = 190, height = 69, units = "mm", device = cairo_pdf)

main_a_pdf <- file.path(main_dir, "Fig06A_Absolute_OC_Trajectories_190mm.pdf")
ggsave(
  main_a_pdf, absolute_trajectory_compact,
  width = 95, height = 69, units = "mm", device = cairo_pdf
)

supp_a_dir <- path_project("fig06", "figS06A")
dir.create(supp_a_dir, recursive = TRUE, showWarnings = FALSE)
supp_a_pdf <- file.path(supp_a_dir, "FigS06A_Relative_OC_Trajectories_190mm.pdf")
ggsave(
  supp_a_pdf, relative_trajectory_compact,
  width = 95, height = 69, units = "mm", device = cairo_pdf
)

export_plot_data <- bind_rows(
  absolute_plot_data %>% mutate(Display = "Figure 6A absolute concentration"),
  relative_plot_data %>% mutate(Display = "Figure S6A Fast-normalized relative value")
) %>%
  transmute(
    Display, Participant, Marker, Timepoint, AbsoluteValue, FastValue,
    AbsoluteChangeFromFast, RelativeToFast, RelativeToPre
  ) %>%
  arrange(Marker, Participant, Timepoint)

display_exclusions <- tibble(
  Display = c(
    "Figure 6A absolute concentration",
    "Figure S6A Fast-normalized relative value"
  ),
  Participant = c(absolute_display_outlier, relative_display_outlier),
  Reason = c(
    "Largest observed tOC/cOC absolute concentration",
    "Largest observed post-exercise cOC relative-to-Pre value"
  ),
  Inference = "Display-only exclusion; model estimates and BH-FDR use the complete n = 50 cohort"
)

methods <- tribble(
  ~Item, ~Description,
  "Main one-row figure", "Paired absolute concentrations at Fast, Pre, Post 1 h, and Post 3 h after one display-only high-value exclusion (n = 49), shown beside Pre-exercise clinical associations (n = 50).",
  "Supplementary Figure S6A", "Paired concentrations normalized to each participant's Fast value (Fast = 1) after one display-only influential cOC exclusion (n = 49).",
  "Time-effect model", "log2(marker) ~ Timepoint + Role + (1|FamilyID) + (1|ParticipantID), fitted with REML.",
  "Multiplicity", "BH correction was applied within each prespecified osteocalcin contrast family across cOC, tOC, and cOC/tOC ratio.",
  "Clinical model", "Rank-normalized clinical phenotype ~ standardized state-specific osteocalcin measure + Role + (1|FamilyID).",
  "Interpretation", "Fast is a fasting reference and Pre is the formal exercise baseline. Display exclusions do not alter the complete-cohort mixed-model estimates or BH-FDR annotations."
)

writexl::write_xlsx(
  list(
    Display_exclusions = display_exclusions,
    Time_effects = relative_effects,
    Plot_data = export_plot_data,
    Methods = methods
  ),
  file.path(analysis_dir, "Fig06AB_FourTimepoint_OriginalLayout_Source_Data.xlsx")
)

writeLines(
  c(
    "# Figure 6A-B and Figure S6A four-timepoint display update",
    "",
    "The 190 x 69 mm geometry, equal panel widths, clinical panel, borders, typography, and legend are inherited from the reviewed one-row layout.",
    "Figure 6A shows absolute concentrations after one display-only high-value exclusion; tOC and cOC use one shared y-axis, and the ratio uses an independent scale with the same number of tick labels.",
    "Figure S6A shows concentrations normalized to Fast (Fast = 1).",
    "FDR labels correspond to Pre minus Fast, Post 1 h minus Pre, and Post 3 h minus Pre, respectively; Pre remains the formal exercise baseline for inferential exercise contrasts.",
    "Each trajectory panel uses one prespecified display-only exclusion recorded in the source workbook; inference remains based on all 50 participants.",
    "No clinical association result in Figure 6B was changed."
  ),
  file.path(analysis_dir, "README.md"),
  useBytes = TRUE
)

writeLines(
  c(
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Main PDF: ", main_pdf),
    paste0("Figure 6A PDF: ", main_a_pdf),
    paste0("Figure S6A PDF: ", supp_a_pdf),
    paste0("Relative display participants: ", n_distinct(relative_plot_data$Participant)),
    paste0("Absolute display participants: ", n_distinct(absolute_plot_data$Participant)),
    "",
    capture.output(sessionInfo())
  ),
  file.path(analysis_dir, "RUN_LOG.txt"),
  useBytes = TRUE
)

message("Updated Figure 6A-B in the reviewed one-row layout.")
