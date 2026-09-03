# Figure 4A: clinical-insulin heatmap and four aligned association forests.

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
if (is.na(.local_config)) stop("Cannot find CellMetabolism_Transfer code/Fig00_Config.R.")
source(.local_config)
rm(.local_script, .local_config_candidates, .local_config)

p_load(readr, dplyr, tidyr, tibble, stringr, ggplot2, patchwork, writexl)

input_file <- path_project(
  "fig04", "fig04A", "Fig04A_Clinical_Insulin_SourceData.csv"
)
if (!file.exists(input_file)) stop("Missing Figure 4A source data: ", input_file)

output_dir <- path_project("fig04", "analysis", "Fig04A_YAxis12_Heatmap_Four_Forests")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metric_order <- c("Fasting insulin", "Fasting glucose", "HOMA-IR", "Matsuda ISI")
block_order <- c(
  "General and central adiposity",
  "Regional fat distribution",
  "Abdominal and ectopic fat",
  "Adipokines and circulating lipids",
  "Fitness and energy metabolism"
)
block_labels <- c(
  "General and central adiposity" = "Adiposity",
  "Regional fat distribution" = "Fat distribution",
  "Abdominal and ectopic fat" = "Ectopic fat",
  "Adipokines and circulating lipids" = "Lipid metabolism",
  "Fitness and energy metabolism" = "Fitness"
)

figure_data <- readr::read_csv(input_file, show_col_types = FALSE)
required_columns <- c(
  "FunctionalBlock", "Outcome", "InsulinMetric", "N", "Families",
  "Effect", "CI_Low", "CI_High", "P_Value", "BH_FDR", "SingularFit",
  "SelectionReason"
)
missing_columns <- setdiff(required_columns, names(figure_data))
if (length(missing_columns) > 0L) {
  stop("Figure 4A source data lacks: ", paste(missing_columns, collapse = ", "))
}
outcome_order <- unique(figure_data$Outcome)
if (length(outcome_order) != 12L || nrow(figure_data) != 48L ||
    anyDuplicated(figure_data[c("Outcome", "InsulinMetric")])) {
  stop("Expected 12 outcomes crossed with four insulin-related metrics.")
}
if (!setequal(unique(figure_data$InsulinMetric), metric_order)) {
  stop("Unexpected insulin-related metric labels in Figure 4A source data.")
}

selected <- figure_data %>%
  group_by(FunctionalBlock, Outcome, SelectionReason) %>%
  summarise(
    Insulin_Significant_Count = sum(BH_FDR < 0.05, na.rm = TRUE),
    Significant_Insulin_Metrics = paste(
      InsulinMetric[is.finite(BH_FDR) & BH_FDR < 0.05], collapse = "; "
    ),
    Minimum_Insulin_BH_FDR = min(BH_FDR, na.rm = TRUE),
    Any_Singular_Insulin_Model = any(SingularFit %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(DisplayOrder = match(Outcome, outcome_order)) %>%
  arrange(DisplayOrder)

figure_data <- figure_data %>%
  mutate(DisplayOrder = match(Outcome, outcome_order)) %>%
  mutate(
    FunctionalBlock = factor(FunctionalBlock, levels = block_order),
    Outcome = factor(Outcome, levels = rev(selected$Outcome)),
    InsulinMetric = factor(InsulinMetric, levels = metric_order),
    Significance = case_when(
      BH_FDR < 0.001 ~ "***",
      BH_FDR < 0.01 ~ "**",
      BH_FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    DisplayOutcome = if_else(as.character(Outcome) == "VO2max", "VO2max", as.character(Outcome)),
    StarX = if_else(CI_High > 0.72, 0.84, CI_High + 0.035),
    StarHjust = if_else(CI_High > 0.72, 1, 0)
  )

metric_colors <- c(
  "Fasting insulin" = "#C86B5A",
  "Fasting glucose" = "#4E9F8E",
  "HOMA-IR" = "#7A5AA6",
  "Matsuda ISI" = "#2B78A6"
)

common_theme <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.text = element_text(color = "black", size = 6.5),
    axis.title = element_text(color = "black", size = 7),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    plot.title = element_text(color = "black", size = 8, face = "bold", hjust = 0.5),
    strip.background = element_blank(),
    panel.spacing.y = unit(0.8, "mm")
  )

heatmap_plot <- ggplot(figure_data, aes(InsulinMetric, Outcome, fill = Effect)) +
  geom_tile(color = "black", linewidth = 0.18) +
  geom_text(aes(label = Significance), color = "black", size = 2.45, fontface = "bold") +
  facet_grid(
    rows = vars(FunctionalBlock), scales = "free_y", space = "free_y", switch = "y",
    labeller = as_labeller(block_labels)
  ) +
  scale_fill_gradient2(
    low = "#3E6E9E", mid = "white", high = "#C8564B",
    midpoint = 0, limits = c(-0.6, 0.6), oob = scales::squish,
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
    guide = "none"
  ) +
  scale_x_discrete(labels = c(
    "Fasting insulin" = "Insulin",
    "Fasting glucose" = "Glucose",
    "HOMA-IR" = "HOMA",
    "Matsuda ISI" = "Matsuda"
  ), position = "top") +
  scale_y_discrete(labels = c("VO2max" = "VO2max"), position = "right") +
  labs(title = "Association matrix", x = NULL, y = NULL) +
  common_theme +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, color = "black", size = 6.0, face = "bold"),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, vjust = 0, color = "black", size = 5.8),
    axis.text.y = element_text(color = "black", size = 6.2),
    axis.ticks = element_blank(),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )

make_forest <- function(metric) {
  dat <- figure_data %>% filter(InsulinMetric == metric)
  metric_title <- dplyr::recode(
    metric,
    "Fasting insulin" = "Fasting\ninsulin",
    "Fasting glucose" = "Fasting\nglucose",
    "HOMA-IR" = "HOMA-IR",
    "Matsuda ISI" = "Matsuda ISI"
  )
  ggplot(dat, aes(Effect, Outcome)) +
    geom_vline(xintercept = 0, linewidth = 0.22, linetype = "dashed", color = "black") +
    geom_errorbarh(aes(xmin = CI_Low, xmax = CI_High), height = 0, linewidth = 0.34, color = "black") +
    geom_point(size = 1.8, shape = 16, color = unname(metric_colors[metric])) +
    geom_text(
      aes(x = StarX, label = Significance, hjust = StarHjust),
      color = "black", size = 2.45, fontface = "bold"
    ) +
    facet_grid(rows = vars(FunctionalBlock), scales = "free_y", space = "free_y") +
    scale_x_continuous(
      limits = c(-0.9, 0.9), breaks = c(-0.8, 0, 0.8),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(title = metric_title, x = NULL, y = NULL) +
    common_theme +
    theme(
      strip.text.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(color = "black", size = 5.8),
      axis.title.x = element_text(color = "black", size = 6.2),
      plot.title = element_text(color = "black", size = 7.5, face = "bold", hjust = 0.5, lineheight = 0.9),
      plot.margin = margin(2, 1.5, 2, 1.5, unit = "mm")
    )
}

forest_plots <- lapply(metric_order, make_forest)

legend_gradient <- tibble(
  x = 0,
  y = seq(-0.6, 0.6, length.out = 121)
)
legend_plot <- ggplot(legend_gradient, aes(x, y, fill = y)) +
  geom_tile(width = 0.28, height = 0.011) +
  scale_fill_gradient2(
    low = "#3E6E9E", mid = "white", high = "#C8564B",
    midpoint = 0, limits = c(-0.6, 0.6), guide = "none"
  ) +
  annotate("text", x = 0.14, y = 1.00, label = "Standardized\neffect",
           family = "Arial", size = 2.35, fontface = "bold", hjust = 0.5, color = "black") +
  annotate("segment", x = 0.15, xend = 0.30, y = c(-0.6, -0.3, 0, 0.3, 0.6),
           yend = c(-0.6, -0.3, 0, 0.3, 0.6), linewidth = 0.22, color = "black") +
  annotate("text", x = 0.40, y = c(-0.6, -0.3, 0, 0.3, 0.6),
           label = c("-0.6", "-0.3", "0", "0.3", "0.6"),
           family = "Arial", size = 2.1, hjust = 0, color = "black") +
  annotate("text", x = -0.05, y = -0.84,
           label = "BH-FDR\n* < 0.05\n** < 0.01\n*** < 0.001",
           family = "Arial", size = 2.1, hjust = 0, vjust = 1, lineheight = 1.15, color = "black") +
  coord_cartesian(xlim = c(-0.12, 1.05), ylim = c(-1.35, 1.15), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(8, 1, 8, 1, unit = "mm"))

combined <- heatmap_plot + forest_plots[[1]] + forest_plots[[2]] + forest_plots[[3]] + forest_plots[[4]] + legend_plot +
  plot_layout(widths = c(2.25, 1.20, 1.20, 1.20, 1.20, 0.90)) +
  plot_annotation(
    caption = "Points and horizontal lines denote standardized effects and 95% confidence intervals.",
    theme = theme(
      plot.caption = element_text(family = "Arial", color = "black", size = 5.8, hjust = 0)
    )
  )

ggsave(
  file.path(output_dir, "Fig04A_Clinical_Insulin_Heatmap_Four_Forests.pdf"),
  combined, width = 174, height = 88, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Fig04A_Clinical_Insulin_Heatmap_Four_Forests.png"),
  combined, width = 174, height = 88, units = "mm", dpi = 500, bg = "white"
)

ggsave(file.path(output_dir, "Fig04A_Heatmap.pdf"), heatmap_plot,
       width = 78, height = 92, units = "mm", device = cairo_pdf)
for (i in seq_along(metric_order)) {
  safe_name <- str_replace_all(metric_order[i], "[^A-Za-z0-9]+", "_")
  ggsave(
    file.path(output_dir, paste0("Fig04A_Forest_", safe_name, ".pdf")),
    forest_plots[[i]], width = 42, height = 92, units = "mm", device = cairo_pdf
  )
}

source_table <- figure_data %>%
  mutate(
    FunctionalBlock = as.character(FunctionalBlock),
    Outcome = as.character(Outcome),
    InsulinMetric = as.character(InsulinMetric)
  ) %>%
  arrange(DisplayOrder, InsulinMetric) %>%
  select(
    FunctionalBlock, Outcome, InsulinMetric, N, Families,
    Effect, CI_Low, CI_High, P_Value, BH_FDR, Significance, SingularFit,
    SelectionReason
  )

selection_table <- selected %>%
  select(
    DisplayOrder, FunctionalBlock, Outcome, SelectionReason,
    Insulin_Significant_Count, Significant_Insulin_Metrics,
    Minimum_Insulin_BH_FDR, Any_Singular_Insulin_Model
  )

writexl::write_xlsx(
  list(
    `Figure source data` = source_table,
    `Y-axis selection` = selection_table
  ),
  file.path(output_dir, "Table_Fig04A_Clinical_Insulin_Associations.xlsx")
)

readr::write_csv(source_table, file.path(output_dir, "Fig04A_Clinical_Insulin_SourceData.csv"))

capture.output(
  {
    cat("Figure 4A clinical-insulin heatmap and four forests\n")
    cat("Locked figure source:", input_file, "\n")
    cat("Rows: 12 evidence-ranked clinical phenotypes grouped into five functional blocks.\n")
    cat("Columns: fasting insulin, fasting glucose, HOMA-IR, and Matsuda ISI.\n")
    cat("Model effects and BH-FDR are inherited from the 23-outcome x 4-metric analysis (92 tests).\n")
    cat("Figure size: 174 x 92 mm; PDF vector and 500-dpi PNG.\n\n")
    print(sessionInfo())
  },
  file = file.path(output_dir, "run_log.txt")
)

message("Completed Figure 4A heatmap, four forests, composite, and source workbook.")
