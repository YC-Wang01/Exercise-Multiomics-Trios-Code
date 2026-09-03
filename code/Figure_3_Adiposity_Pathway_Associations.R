# Re-render Figure 3E as a portrait heatmap from the locked source-data CSV.
# Statistical values are not recalculated by this script.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

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
project_dir <- PROJECT_ROOT
rm(.local_script, .config_candidates, .config_file)

figure_dir <- file.path(project_dir, "fig03", "Figure3")
working_dir <- file.path(project_dir, "fig03", "Temp_Pathway_Response_Figures_2026-08-26")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(working_dir, recursive = TRUE, showWarnings = FALSE)

source_file <- file.path(
  figure_dir,
  "Fig03E_Adiposity_Pathway_Response_Associations_SourceData.csv"
)
data <- read_csv(source_file, show_col_types = FALSE)

required_columns <- c(
  "Tissue", "Pathway", "Phenotype", "ClusterOrder",
  "StandardizedEffect", "BH_FDR_WithinPhenotype"
)
missing_columns <- setdiff(required_columns, names(data))
if (length(missing_columns)) {
  stop("Missing source columns: ", paste(missing_columns, collapse = ", "))
}

pathway_order <- tibble::tribble(
  ~Tissue, ~Pathway, ~DisplayOrder,
  "Adipose", "Adipogenesis", 1,
  "Adipose", "Fatty acid metabolism", 2,
  "Adipose", "Oxidative phosphorylation", 3,
  "Adipose", "Heme metabolism", 4,
  "Muscle", "Oxidative phosphorylation", 5,
  "Muscle", "Coagulation", 6,
  "Muscle", "Complement", 7,
  "Muscle", "KRAS signaling up", 8,
  "Muscle", "Apical junction", 9
)

phenotype_order <- data %>%
  distinct(Phenotype, ClusterOrder) %>%
  arrange(ClusterOrder) %>%
  pull(Phenotype)

plot_data <- data %>%
  inner_join(pathway_order, by = c("Tissue", "Pathway")) %>%
  mutate(
    XPos = DisplayOrder,
    YPos = length(phenotype_order) + 1L - ClusterOrder,
    FDRMark = case_when(
      BH_FDR_WithinPhenotype < 0.05 ~ "***",
      BH_FDR_WithinPhenotype < 0.10 ~ "**",
      BH_FDR_WithinPhenotype < 0.20 ~ "*",
      TRUE ~ ""
    )
  )

if (nrow(plot_data) != 9L * length(phenotype_order)) {
  stop("Unexpected Figure 3E matrix dimensions: ", nrow(plot_data), " rows.")
}

pt_to_mm <- 1 / ggplot2::.pt
pathway_labels <- pathway_order$Pathway[order(pathway_order$DisplayOrder)]
phenotype_labels <- phenotype_order

figure_3e <- ggplot(plot_data, aes(XPos, YPos)) +
  geom_tile(
    aes(fill = StandardizedEffect),
    width = 1, height = 1, color = "white", linewidth = 0.42
  ) +
  geom_text(
    aes(label = FDRMark),
    family = "Arial", fontface = "bold", size = 6.5 * pt_to_mm,
    color = "black", vjust = 0.62
  ) +
  annotate(
    "rect", xmin = 0.5, xmax = 4.5, ymin = 17.55, ymax = 18.35,
    fill = "white", color = "black", linewidth = 0.45
  ) +
  annotate(
    "rect", xmin = 4.5, xmax = 9.5, ymin = 17.55, ymax = 18.35,
    fill = "white", color = "black", linewidth = 0.45
  ) +
  annotate(
    "text", x = 2.5, y = 17.95, label = "Adipose",
    family = "Arial", fontface = "bold", size = 7 * pt_to_mm
  ) +
  annotate(
    "text", x = 7, y = 17.95, label = "Muscle",
    family = "Arial", fontface = "bold", size = 7 * pt_to_mm
  ) +
  scale_x_continuous(
    breaks = seq_along(pathway_labels), labels = pathway_labels,
    position = "top", expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = seq_along(phenotype_labels), labels = rev(phenotype_labels),
    expand = expansion(mult = 0), position = "left"
  ) +
  scale_fill_gradient2(
    low = "#7FA9C7", mid = "white", high = "#D9897F",
    midpoint = 0, limits = c(-1.1, 1.1), oob = scales::squish,
    name = "Standardized effect"
  ) +
  coord_fixed(
    xlim = c(0.5, 9.5), ylim = c(0.5, 18.35),
    ratio = 1, clip = "off", expand = FALSE
  ) +
  labs(
    title = "Adiposity associations with pathway responses",
    x = NULL, y = NULL,
    caption = "Within-phenotype BH-FDR: * <0.20, ** <0.10, *** <0.05"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = unit(23, "mm"), barheight = unit(2.8, "mm")
    )
  ) +
  theme_void(base_family = "Arial", base_size = 6.5) +
  theme(
    plot.title = element_text(
      family = "Arial", face = "bold", size = 9,
      hjust = 0.5, margin = margin(b = 2)
    ),
    axis.text.x.top = element_text(
      family = "Arial", face = "bold", size = 5.5,
      angle = 45, hjust = 0, vjust = 0, color = "black",
      margin = margin(b = 1.5)
    ),
    axis.text.y.left = element_text(
      family = "Arial", size = 5.5, color = "black",
      hjust = 1, margin = margin(r = 2)
    ),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.45),
    legend.position = "bottom",
    legend.title = element_text(family = "Arial", face = "bold", size = 7),
    legend.text = element_text(family = "Arial", size = 5.5),
    plot.caption = element_text(
      family = "Arial", size = 5.5, hjust = 0.5,
      color = "black", margin = margin(t = 1.5)
    ),
    plot.margin = margin(2, 4, 2, 4)
  )

output_pdf <- file.path(figure_dir, "Fig03E_Adiposity_Associated_Pathway_Responses.pdf")
output_png <- file.path(figure_dir, "Fig03E_Adiposity_Associated_Pathway_Responses.png")
working_pdf <- file.path(
  working_dir,
  "Fig03_Adiposity_Pathway_Response_Associations_Clustered_FDR20_Candidate.pdf"
)

ggsave(
  output_pdf, figure_3e,
  width = 102 / 25.4, height = 108 / 25.4,
  device = cairo_pdf, bg = "white"
)
ggsave(
  output_png, figure_3e,
  width = 102 / 25.4, height = 108 / 25.4,
  dpi = 450, bg = "white"
)
ggsave(
  working_pdf, figure_3e,
  width = 102 / 25.4, height = 108 / 25.4,
  device = cairo_pdf, bg = "white"
)

writeLines(
  c(
    "Figure 3E portrait re-render",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Source: ", source_file),
    "Statistics were not recalculated.",
    "Rows: adiposity phenotypes in the locked clustered order.",
    "Columns: pathways grouped by tissue.",
    "Color: standardized association effect.",
    "Stars: within-phenotype BH-FDR (* <0.20, ** <0.10, *** <0.05)."
  ),
  file.path(figure_dir, "Fig03E_Vertical_Render_RunLog.txt"),
  useBytes = TRUE
)

message("Wrote: ", output_pdf)
