suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
  library(stringr)
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
figure_dir <- file.path(project_dir, "fig03", "FigureS3")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

go_display_file <- file.path(
  figure_dir,
  "FigS03B_GO_Pathway_Response_DisplayData.csv"
)
family_file <- file.path(
  figure_dir,
  "FigS03F_Family_Module_Stability_SourceData.xlsx"
)
family_all_file <- file.path(
  figure_dir,
  "FigS03G_Family_Metabolic_Response_Heatmap_SourceData.xlsx"
)

stopifnot(
  file.exists(go_display_file),
  file.exists(family_file),
  file.exists(family_all_file)
)

font_family <- "Arial"
positive_color <- "#D16A5D"
negative_color <- "#4E82A8"
pair_colors <- c(
  "Daughter-mother" = "#C94F45",
  "Daughter-father" = "#2F6FA3",
  "Mother-father" = "#5B8E7D"
)
pair_shapes <- c(
  "Daughter-mother" = 16,
  "Daughter-father" = 17,
  "Mother-father" = 15
)
cluster_colors <- c(
  "Cluster 1" = "#D65A31",
  "Cluster 2" = "#3B82A0",
  "Cluster 3" = "#5B8E7D",
  "Cluster 4" = "#C45B7A"
)
cluster_labels <- c(
  "Cluster 1" = "C1: VLDL/triglycerides",
  "Cluster 2" = "C2: Glycolysis/amino acids",
  "Cluster 3" = "C3: Ketone/BCAA/fatty acids",
  "Cluster 4" = "C4: HDL remodeling"
)

panel_theme <- theme_classic(base_family = font_family, base_size = 6.4) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.32),
    axis.text = element_text(color = "black", size = 5.8),
    axis.title = element_text(color = "black", size = 6.4, face = "bold"),
    axis.ticks = element_line(linewidth = 0.28, color = "black"),
    plot.title = element_text(color = "black", size = 8.2, face = "bold"),
    plot.subtitle = element_text(color = "black", size = 5.5, margin = margin(b = 2.2)),
    legend.title = element_text(color = "black", size = 5.4),
    legend.text = element_text(color = "black", size = 5.1),
    plot.margin = margin(2.2, 2.2, 2.2, 2.2)
  )

# -----------------------------------------------------------------------------
# S3B: GO biological-process responses, FDR-only display.
# -----------------------------------------------------------------------------
go_display <- read.csv(go_display_file, check.names = FALSE) %>%
  arrange(DisplayOrder) %>%
  mutate(
    ColumnText = paste(Tissue, Layer),
    Column = factor(
      ColumnText,
      levels = c(
        "Adipose Proteomics",
        "Adipose DNA methylation",
        "Muscle Proteomics",
        "Muscle DNA methylation"
      ),
      labels = c(
        "Adipose\nproteomics",
        "Adipose\nDNA methylation",
        "Muscle\nproteomics",
        "Muscle\nDNA methylation"
      )
    ),
    DirectionLabel = if_else(
      Direction == "Increased",
      "Positive-ranked response",
      "Negative-ranked response"
    ),
    SignificanceStrength = pmin(-log10(pmax(BH_FDR, 1e-12)), 10),
    PathwayLabelText = case_when(
      GO_ID == "GO:0050911" ~ "sensory chemical-stimulus detection",
      GO_ID == "GO:0042776" ~ "mitochondrial ATP synthesis",
      GO_ID == "GO:0032981" ~ "respiratory-chain complex I assembly",
      GO_ID == "GO:0006958" ~ "classical complement activation",
      TRUE ~ as.character(Pathway)
    ),
    PathwayLabelText = str_wrap(PathwayLabelText, width = 34),
    PathwayLabel = factor(
      PathwayLabelText,
      levels = rev(unique(PathwayLabelText))
    )
  )

if (any(go_display$BH_FDR >= 0.20, na.rm = TRUE)) {
  stop("S3B display contains a pathway with BH-FDR >= 0.20.")
}

plot_b <- ggplot(go_display, aes(Column, PathwayLabel)) +
  geom_point(
    aes(size = SignificanceStrength, fill = DirectionLabel),
    shape = 21, color = "black", stroke = 0.32
  ) +
  scale_fill_manual(
    values = c(
      "Negative-ranked response" = negative_color,
      "Positive-ranked response" = positive_color
    ),
    name = "Ranked response"
  ) +
  scale_size_continuous(
    range = c(1.7, 5.0), limits = c(0, 10),
    breaks = c(1, 4, 8, 10), labels = c("1", "4", "8", ">=10"),
    name = expression(-log[10]~"(BH-FDR)")
  ) +
  scale_y_discrete(expand = expansion(add = c(0.45, 0.45))) +
  labs(
    title = "B  GO biological-process responses",
    subtitle = "Post3h - Pre; displayed pathways meet BH-FDR < 0.20",
    x = NULL, y = NULL
  ) +
  panel_theme +
  theme(
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 5.6, face = "bold"),
    axis.text.y = element_text(size = 5.35, lineheight = 0.88),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = grid::unit(0.7, "mm"),
    legend.spacing.y = grid::unit(0.2, "mm")
  ) +
  guides(
    size = guide_legend(order = 1, title.position = "top", nrow = 1),
    fill = guide_legend(order = 2, title.position = "top", nrow = 1,
                        override.aes = list(size = 2.5))
  )

# -----------------------------------------------------------------------------
# S3C: pair-specific within-family similarity. No ICC is repeated here.
# -----------------------------------------------------------------------------
pair_summary <- read.xlsx(
  family_all_file,
  sheet = "Similarity Summary",
  check.names = FALSE
) %>%
  mutate(
    Scope = factor(
      Scope,
      levels = rev(c("Overall response", paste0("Cluster ", 1:4)))
    ),
    PairType = factor(
      PairType,
      levels = c("Daughter-mother", "Daughter-father", "Mother-father")
    )
  )

scope_labels <- c(
  "Overall response" = "Overall response",
  cluster_labels
)

plot_c <- ggplot(
  pair_summary,
  aes(DeltaRho, Scope, color = PairType, shape = PairType)
) +
  geom_vline(xintercept = 0, linewidth = 0.30, color = "black", linetype = "22") +
  geom_errorbarh(
    aes(xmin = Lower95, xmax = Upper95),
    height = 0, linewidth = 0.44,
    position = position_dodge(width = 0.52)
  ) +
  geom_point(size = 1.75, stroke = 0.35, position = position_dodge(width = 0.52)) +
  scale_color_manual(
    values = pair_colors,
    labels = c("D-M", "D-F", "M-F"),
    name = "Family pair"
  ) +
  scale_shape_manual(
    values = pair_shapes,
    labels = c("D-M", "D-F", "M-F"),
    name = "Family pair"
  ) +
  scale_x_continuous(
    limits = c(-0.42, 0.42),
    breaks = c(-0.4, -0.2, 0, 0.2, 0.4)
  ) +
  scale_y_discrete(labels = scope_labels) +
  labs(
    title = "C  Within-family response similarity",
    subtitle = "Family pairs versus role-matched unrelated pairs; all BH-FDR >= 0.40",
    x = expression(Delta*" Spearman "*rho), y = NULL,
    color = "Family pair", shape = "Family pair"
  ) +
  panel_theme +
  theme(
    axis.text.y = element_text(size = 5.0),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.key.width = grid::unit(2.4, "mm"),
    legend.spacing.x = grid::unit(0.5, "mm")
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      override.aes = list(linewidth = 0.5)
    ),
    shape = guide_legend(nrow = 1)
  )

# -----------------------------------------------------------------------------
# S3D: family-level ICC with bootstrap confidence and LOFO stability ranges.
# -----------------------------------------------------------------------------
icc_primary <- read.xlsx(
  family_file,
  sheet = "Primary Family ICC",
  check.names = FALSE
) %>%
  select(
    Cluster, DominantTime, Subjects, Families, ICC,
    Lower95, Upper95, PermutationP, BH_FDR
  )
icc_lofo <- read.xlsx(
  family_file,
  sheet = "LOFO ICC Summary",
  check.names = FALSE
) %>%
  select(Cluster, LOFOMin, LOFOMax)

icc_display <- icc_primary %>%
  left_join(icc_lofo, by = "Cluster") %>%
  mutate(
    Cluster = factor(Cluster, levels = rev(paste0("Cluster ", 1:4))),
    FDRLabel = if_else(BH_FDR < 0.05, "*", "")
  )

plot_d <- ggplot(icc_display, aes(ICC, Cluster, color = Cluster)) +
  geom_vline(xintercept = 0, linewidth = 0.30, color = "black", linetype = "22") +
  geom_errorbarh(
    aes(xmin = pmax(0, Lower95), xmax = pmin(0.70, Upper95)),
    height = 0, linewidth = 0.34, color = "#9A9A9A"
  ) +
  geom_errorbarh(
    aes(xmin = LOFOMin, xmax = LOFOMax),
    height = 0, linewidth = 1.05
  ) +
  geom_point(size = 2.15, stroke = 0.35) +
  geom_text(
    aes(x = pmin(0.68, pmax(Upper95, LOFOMax) + 0.035), label = FDRLabel),
    color = "black", family = font_family, fontface = "bold", size = 2.8
  ) +
  scale_color_manual(values = cluster_colors, guide = "none") +
  scale_y_discrete(labels = cluster_labels) +
  scale_x_continuous(
    limits = c(0, 0.70), breaks = c(0, 0.2, 0.4, 0.6),
    expand = expansion(mult = c(0, 0.015))
  ) +
  labs(
    title = "D  Family-level response clustering",
    subtitle = "ICC (point), 95% CI (gray), LOFO range (color); * FDR < 0.05",
    x = "Family-level ICC", y = NULL
  ) +
  panel_theme +
  theme(axis.text.y = element_text(size = 5.25))

# -----------------------------------------------------------------------------
# S3E: anonymous family-by-role module-score heatmap.
# -----------------------------------------------------------------------------
module_scores <- read.xlsx(
  family_file,
  sheet = "Anonymous Module Scores",
  check.names = FALSE
) %>%
  mutate(
    FamilyNumber = as.integer(str_extract(DisplayFamily, "[0-9]+")),
    DisplayFamily = factor(
      DisplayFamily,
      levels = unique(DisplayFamily[order(FamilyNumber)])
    ),
    Role = factor(
      Role,
      levels = c("Daughter", "Mother", "Father"),
      labels = c("D", "M", "F")
    ),
    Cluster = factor(Cluster, levels = paste0("Cluster ", 1:4))
  )

heatmap_cluster_labels <- c(
  "Cluster 1" = "C1: VLDL/\ntriglycerides",
  "Cluster 2" = "C2: Glycolysis/\namino acids",
  "Cluster 3" = "C3: Ketone/BCAA/\nfatty acids",
  "Cluster 4" = "C4: HDL\nremodeling"
)

plot_e <- ggplot(module_scores, aes(Role, DisplayFamily, fill = ModuleScore)) +
  geom_tile(color = "black", linewidth = 0.22) +
  facet_grid(. ~ Cluster, scales = "free_x", space = "free_x",
             labeller = as_labeller(heatmap_cluster_labels)) +
  scale_fill_gradient2(
    low = "#4E82A8", mid = "white", high = "#D16A5D",
    midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish,
    na.value = "#8A8A8A", name = "Standardized\nmodule score",
    guide = guide_colorbar(
      barheight = grid::unit(34, "mm"),
      barwidth = grid::unit(3.2, "mm"),
      ticks = TRUE, frame.colour = "black"
    )
  ) +
  labs(
    title = "E  Family metabolic-response modules",
    subtitle = "D, daughter; M, mother; F, father; gray, unavailable",
    x = NULL, y = "Family"
  ) +
  theme_minimal(base_family = font_family, base_size = 6.0) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.32),
    axis.text.x = element_text(color = "black", size = 5.0),
    axis.text.y = element_text(color = "black", size = 4.8),
    axis.title.y = element_text(color = "black", size = 6.2, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 5.2, face = "bold"),
    strip.text.x = element_text(margin = margin(b = 1.0)),
    plot.title = element_text(color = "black", size = 8.2, face = "bold"),
    plot.subtitle = element_text(color = "black", size = 5.5, margin = margin(b = 2.2)),
    legend.position = "right",
    legend.title = element_text(color = "black", size = 5.2),
    legend.text = element_text(color = "black", size = 4.9),
    plot.margin = margin(2.2, 2.2, 2.2, 2.2)
  )

# Natural panel sizes correspond to their final shares in the 190-mm assembly.
ggsave(
  file.path(figure_dir, "FigS03B_GO_Pathway_Response_BH_FDR20.pdf"),
  plot_b, width = 118 / 25.4, height = 76 / 25.4,
  device = grDevices::cairo_pdf
)
ggsave(
  file.path(figure_dir, "FigS03C_Within_Family_Response_Similarity.pdf"),
  plot_c, width = 72 / 25.4, height = 76 / 25.4,
  device = grDevices::cairo_pdf
)
ggsave(
  file.path(figure_dir, "FigS03D_Family_Module_Stability.pdf"),
  plot_d, width = 72 / 25.4, height = 82 / 25.4,
  device = grDevices::cairo_pdf
)
ggsave(
  file.path(figure_dir, "FigS03E_Family_Metabolic_Response_Heatmap.pdf"),
  plot_e, width = 118 / 25.4, height = 82 / 25.4,
  device = grDevices::cairo_pdf
)

top_row <- plot_b + plot_c + plot_layout(widths = c(1.64, 1.00))
bottom_row <- plot_d + plot_e + plot_layout(widths = c(1.00, 1.64))
combined <- top_row / bottom_row +
  plot_layout(heights = c(1.00, 1.08)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

combined_pdf <- file.path(
  figure_dir,
  "FigS03BCDE_TwoRow_190mm.pdf"
)
ggsave(
  combined_pdf, combined,
  width = 190 / 25.4, height = 160 / 25.4,
  device = grDevices::cairo_pdf
)
ggsave(
  file.path(figure_dir, "FigS03BCDE_TwoRow_190mm.png"),
  combined,
  width = 190 / 25.4, height = 160 / 25.4,
  dpi = 600, bg = "white"
)

write.xlsx(
  list(
    `S3B GO pathways` = as.data.frame(go_display),
    `S3C similarity` = as.data.frame(pair_summary),
    `S3D family ICC` = as.data.frame(icc_display),
    `S3E module scores` = as.data.frame(module_scores)
  ),
  file.path(figure_dir, "FigS03BCDE_Display_SourceData.xlsx"),
  overwrite = TRUE
)

writeLines(
  c(
    "AnalysisID: FigS03BCDE_Final_Panels_2026-08-27",
    "S3B displays GO biological-process pathways meeting BH-FDR < 0.20.",
    "S3C compares within-family response-profile similarity with role-matched unrelated pairs using family-label permutation; BH-FDR is across 15 scope-by-pair tests.",
    "S3D shows family-level ICC estimates, bootstrap 95% confidence intervals, and leave-one-family-out ranges.",
    "S3E displays standardized serum-metabolic response-module scores using anonymous family identifiers.",
    "The family analyses use 99 serum metabolites with acute-response BH-FDR < 0.05 and adjust response profiles for body-fat percentage and family role.",
    "Family-associated clustering does not estimate heritability and does not distinguish genetic from shared-environment effects.",
    capture.output(sessionInfo())
  ),
  file.path(figure_dir, "FigS03BCDE_RunLog.txt")
)

message("Created: ", combined_pdf)
