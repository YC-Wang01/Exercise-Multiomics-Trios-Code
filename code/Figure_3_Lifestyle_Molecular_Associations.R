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
if (is.na(.local_config)) stop("Cannot find code/Fig00_Config.R.")
source(.local_config)
rm(.local_script, .local_config_candidates, .local_config)

p_load(readr, dplyr, pheatmap, ComplexHeatmap, circlize, grid, ggplot2, patchwork)

transfer_dir <- PROJECT_ROOT
input_file <- file.path(
  transfer_dir,
  "fig03/Figure3",
  "Fig03A_Clinical_Trait_Contributions_SourceData.csv"
)
output_dir <- file.path(
  transfer_dir,
  "fig03/Temp_Fig03F_FDR_2026-08-25",
  "Cell_ExactStyle_DisplayedBH_FDR"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cell_size <- 12
metadata_order <- c(
  "PA(Duration)", "PA(Frequence)", "IPA",
  "BMI", "Fat%", "WC", "HOMA-IR", "HFC",
  "IEE", "Protein", "Fats", "Carbohydrate", "Sucrose", "Alcohol"
)

panel_order <- c(
  "Serum Metabonomics",
  "Serum Proteomics",
  "Adipose",
  "Muscle"
)

file_stubs <- c(
  "Serum Metabonomics" = "Serum_Metabonomics_Heatmap_BH_FDR",
  "Serum Proteomics" = "Serum_Proteomics_Heatmap_BH_FDR",
  "Adipose" = "Adipose_Heatmap_BH_FDR",
  "Muscle" = "Muscle_Heatmap_BH_FDR"
)

panel_titles <- c(
  "Serum Metabonomics" = "Serum_Metabonomics",
  "Serum Proteomics" = "Serum_Proteomics",
  "Adipose" = "Adipose",
  "Muscle" = "Muscle"
)

plot_data <- read_csv(input_file, show_col_types = FALSE)

make_matrix <- function(panel_data, value_column) {
  feature_order <- unique(panel_data$Feature)
  result <- matrix(
    NA_real_,
    nrow = length(feature_order),
    ncol = length(metadata_order),
    dimnames = list(feature_order, metadata_order)
  )
  for (i in seq_len(nrow(panel_data))) {
    result[panel_data$Feature[i], panel_data$Metadata[i]] <- panel_data[[value_column]][i]
  }
  result
}

summary_rows <- list()

for (panel_name in panel_order) {
  panel_data <- plot_data %>% filter(Panel == panel_name)
  cor_mat <- make_matrix(panel_data, "Spearman_R")
  fdr_mat <- make_matrix(panel_data, "BH_FDR_Displayed_Audited")

  star_mat <- matrix("", nrow = nrow(fdr_mat), ncol = ncol(fdr_mat), dimnames = dimnames(fdr_mat))
  star_mat[is.finite(fdr_mat) & fdr_mat < 0.05] <- "*"
  star_mat[is.finite(fdr_mat) & fdr_mat < 0.01] <- "**"

  col_anno <- data.frame(
    Category = c(rep("PA", 3), rep("Adiposity", 5), rep("Dietary intake", 6)),
    row.names = metadata_order,
    stringsAsFactors = FALSE
  )
  anno_colors <- list(
    Category = c(
      "PA" = "#FFD700",
      "Adiposity" = "#4DBBD5",
      "Dietary intake" = "#00A087"
    )
  )

  pdf_w <- (ncol(cor_mat) * cell_size / 72) + 3.5
  pdf_h <- (nrow(cor_mat) * cell_size / 72) + 2.5

  pdf(
    file.path(output_dir, paste0(unname(file_stubs[panel_name]), ".pdf")),
    width = pdf_w,
    height = pdf_h
  )
  pheatmap::pheatmap(
    cor_mat,
    main = unname(panel_titles[panel_name]),
    annotation_col = col_anno,
    annotation_colors = anno_colors,
    annotation_row = NULL,
    color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    treeheight_row = 0,
    display_numbers = star_mat,
    number_color = "black",
    fontsize_number = 10,
    border_color = "white",
    cellwidth = cell_size,
    cellheight = cell_size,
    gaps_col = c(3, 8, 14),
    legend = FALSE,
    annotation_legend = FALSE,
    fontsize = 8,
    angle_col = 45
  )
  dev.off()

  summary_rows[[panel_name]] <- tibble(
    Panel = panel_name,
    DisplayedCells = length(fdr_mat),
    ValidCells = sum(is.finite(fdr_mat)),
    BH_FDR_LT_005 = sum(fdr_mat < 0.05, na.rm = TRUE),
    BH_FDR_LT_001 = sum(fdr_mat < 0.01, na.rm = TRUE),
    MinimumBH_FDR = min(fdr_mat, na.rm = TRUE)
  )
}

pdf(file.path(output_dir, "Fig3F_Standalone_Legend_BH_FDR.pdf"), width = 3, height = 4)
col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("#3C5488", "white", "#E64B35"))
heatmap_legend <- Legend(
  col_fun = col_fun,
  title = "Correlation",
  direction = "horizontal"
)
anno_legend <- Legend(
  labels = c("PA", "Adiposity", "Dietary intake"),
  title = "Category",
  legend_gp = gpar(fill = c("#FFD700", "#4DBBD5", "#00A087"))
)
fdr_legend <- Legend(
  labels = c("BH-FDR < 0.05", "BH-FDR < 0.01"),
  title = "Statistical significance",
  type = "points",
  pch = c("*", "**"),
  size = unit(c(4, 4), "mm"),
  legend_gp = gpar(col = "black")
)
draw(packLegend(heatmap_legend, anno_legend, fdr_legend, direction = "vertical"))
dev.off()

make_compact_panel <- function(panel_name) {
  panel_data <- plot_data %>% filter(Panel == panel_name)
  feature_order <- unique(panel_data$Feature)
  n_features <- length(feature_order)
  panel_data <- panel_data %>%
    mutate(
      x = match(Metadata, metadata_order),
      y = n_features - match(Feature, feature_order) + 1,
      Symbol = case_when(
        BH_FDR_Displayed_Audited < 0.01 ~ "**",
        BH_FDR_Displayed_Audited < 0.05 ~ "*",
        TRUE ~ ""
      )
    )

  category_segments <- tibble(
    x = c(1, 4, 9),
    xend = c(3, 8, 14),
    y = n_features + 0.72,
    Category = factor(
      c("PA", "Adiposity", "Dietary intake"),
      levels = c("PA", "Adiposity", "Dietary intake")
    )
  )

  ggplot(panel_data, aes(x = x, y = y)) +
    geom_tile(aes(fill = Spearman_R), color = "white", linewidth = 0.16) +
    geom_text(aes(label = Symbol), size = 1.65, fontface = "bold", family = "Arial") +
    geom_segment(
      data = category_segments,
      aes(x = x, xend = xend, y = y, yend = y, color = Category),
      inherit.aes = FALSE,
      linewidth = 3.6,
      lineend = "butt"
    ) +
    scale_fill_gradient2(
      low = "#3C5488",
      mid = "white",
      high = "#E64B35",
      midpoint = 0,
      limits = c(-0.8, 0.8),
      oob = scales::squish
    ) +
    scale_color_manual(
      values = c("PA" = "#FFD700", "Adiposity" = "#4DBBD5", "Dietary intake" = "#00A087")
    ) +
    scale_x_continuous(
      breaks = seq_along(metadata_order),
      labels = metadata_order,
      limits = c(0.5, 14.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = n_features:1,
      labels = feature_order,
      limits = c(0.5, n_features + 1.05),
      expand = c(0, 0),
      position = "right"
    ) +
    coord_fixed(clip = "off") +
    labs(title = unname(panel_titles[panel_name]), x = NULL, y = NULL) +
    theme_void(base_family = "Arial") +
    theme(
      axis.text.x = element_text(
        size = 4.0,
        angle = 45,
        hjust = 1,
        vjust = 1,
        color = "black",
        margin = margin(t = 0.6, unit = "mm")
      ),
      axis.text.y = element_text(
        size = 4.2,
        color = "black",
        margin = margin(l = 0.7, unit = "mm")
      ),
      plot.title = element_text(face = "bold", size = 6.2, hjust = 0.5, margin = margin(b = 1.4, unit = "mm")),
      plot.margin = margin(t = 2, r = 11, b = 5, l = 7, unit = "mm"),
      legend.position = "none"
    )
}

legend_gradient <- tibble(
  x = seq(-0.5, 0.5, length.out = 101),
  y = 8.9,
  rho = seq(-0.5, 0.5, length.out = 101)
)

compact_legend <- ggplot() +
  geom_tile(
    data = legend_gradient,
    aes(x = x, y = y, fill = rho),
    width = 0.011,
    height = 0.48
  ) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0) +
  annotate("text", x = -0.5, y = 9.65, label = "Correlation", hjust = 0, size = 2.0, fontface = "bold", family = "Arial") +
  annotate("text", x = c(-0.5, 0, 0.5), y = 8.25, label = c("-0.5", "0", "0.5"), size = 1.55, family = "Arial") +
  annotate("text", x = -0.5, y = 7.3, label = "Category", hjust = 0, size = 2.0, fontface = "bold", family = "Arial") +
  annotate("rect", xmin = -0.5, xmax = -0.38, ymin = c(6.4, 5.6, 4.8), ymax = c(6.75, 5.95, 5.15), fill = c("#FFD700", "#4DBBD5", "#00A087")) +
  annotate("text", x = -0.32, y = c(6.58, 5.78, 4.98), label = c("PA", "Adiposity", "Dietary intake"), hjust = 0, size = 1.55, family = "Arial") +
  annotate("text", x = -0.5, y = 3.75, label = "BH-FDR", hjust = 0, size = 2.0, fontface = "bold", family = "Arial") +
  annotate("text", x = -0.46, y = c(2.9, 2.1), label = c("*", "**"), hjust = 0.5, size = 2.2, fontface = "bold", family = "Arial") +
  annotate("text", x = -0.32, y = c(2.9, 2.1), label = c("< 0.05", "< 0.01"), hjust = 0, size = 1.55, family = "Arial") +
  coord_cartesian(xlim = c(-0.55, 0.65), ylim = c(1.4, 10), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(t = 2, r = 2, b = 4, l = 1, unit = "mm"), legend.position = "none")

compact_panels <- lapply(panel_order, make_compact_panel)
combined_figure <- wrap_plots(
  c(compact_panels, list(compact_legend)),
  nrow = 1,
  widths = c(1, 1, 1, 1, 0.68)
)

ggsave(
  file.path(output_dir, "Fig03F_Cell_ExactStyle_FourPanel_BH_FDR_190mm.pdf"),
  combined_figure,
  width = 190,
  height = 44,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)

write_csv(
  bind_rows(summary_rows),
  file.path(output_dir, "Fig03F_Cell_ExactStyle_BH_FDR_Summary.csv")
)

message("Exact-style BH-FDR heatmaps written to: ", output_dir)
