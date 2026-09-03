# Figure S6G-H: selected Hallmark running-enrichment profiles for osteocalcin.
#
# This script visualizes two prespecified results from the frozen Figure 6
# body-fat-adjusted proteomic analysis. It does not refit statistical models.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  file.path(getwd(), "Fig00_Config.R"),
  if (!is.na(.local_script)) file.path(dirname(normalizePath(
    .local_script, winslash = "/", mustWork = FALSE
  )), "Fig00_Config.R") else NA_character_
))
.config_file <- .config_candidates[file.exists(.config_candidates)][1]
if (is.na(.config_file)) stop("Cannot locate Fig00_Config.R.")
source(.config_file)
rm(.local_script, .config_candidates, .config_file)

p_load(dplyr, ggplot2, patchwork, stringr, tibble, tidyr, writexl)

analysis_id <- "Fig06_OC_BodyFatAdjusted_Hallmark_2026-08-18"
input_file <- path_project(
  "fig06", "analysis", analysis_id,
  "Fig06_OC_BodyFatAdjusted_Hallmark.rds"
)
output_dir <- path_project("fig06", "figS06GH")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Missing frozen Figure 6 analysis: ", input_file)

analysis <- readRDS(input_file)
if (is.null(analysis$FeatureResults) || is.null(analysis$PathwayResults)) {
  stop("Figure 6 analysis object lacks feature or pathway results.")
}

panel_spec <- tribble(
  ~Panel, ~PathwayID, ~Pathway, ~Title, ~FileStem,
  "G",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "Oxidative phosphorylation",
  "Oxidative phosphorylation across osteocalcin models",
  "FigS06G_Oxidative_Phosphorylation_GSEA_Profiles",
  "H",
  "HALLMARK_HEME_METABOLISM", "Heme metabolism",
  "Heme metabolism across osteocalcin models",
  "FigS06H_Heme_Metabolism_GSEA_Profiles"
)

model_spec <- tribble(
  ~Marker, ~MarkerLabel, ~State, ~StateLabel,
  "TotalOC", "tOC", "PreExercise", "Pre-exercise",
  "TotalOC", "tOC", "ExerciseResponse", "Exercise response",
  "cOC", "cOC", "PreExercise", "Pre-exercise",
  "cOC", "cOC", "ExerciseResponse", "Exercise response"
)

hallmark_sets <- load_hallmark_gene_sets()

missing_sets <- setdiff(panel_spec$PathwayID, names(hallmark_sets))
if (length(missing_sets) > 0L) {
  stop("Missing Hallmark set(s): ", paste(missing_sets, collapse = ", "))
}

running_score <- function(stats, pathway_genes) {
  stats <- sort(stats[is.finite(stats)], decreasing = TRUE)
  hit <- names(stats) %in% pathway_genes
  nh <- sum(hit)
  if (nh == 0L || nh == length(stats)) {
    stop("Invalid pathway membership for running enrichment score.")
  }
  hit_weights <- abs(stats) * hit
  hit_increment <- hit_weights / sum(hit_weights)
  miss_increment <- (!hit) / (length(stats) - nh)
  tibble(
    Rank = seq_along(stats),
    RankPercent = 100 * Rank / length(stats),
    RunningES = cumsum(hit_increment - miss_increment),
    IsHit = hit,
    GeneSymbol = names(stats),
    ModeratedT = unname(stats)
  )
}

format_fdr <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < 1e-4) return(formatC(x, format = "e", digits = 2))
  formatC(x, format = "f", digits = 4)
}

significance_symbol <- function(x) {
  case_when(
    !is.finite(x) ~ "",
    x < 0.001 ~ "***",
    x < 0.01 ~ "**",
    x < 0.05 ~ "*",
    TRUE ~ ""
  )
}

prepare_ranks <- function(marker, state, tissue) {
  selected <- analysis$FeatureResults %>%
    filter(Marker == marker, State == state, Tissue == tissue) %>%
    filter(!is.na(GeneSymbol), nzchar(GeneSymbol), is.finite(ModeratedT)) %>%
    mutate(
      GeneSymbol = str_trim(str_split_i(GeneSymbol, ";", 1L))
    ) %>%
    filter(nzchar(GeneSymbol)) %>%
    group_by(GeneSymbol) %>%
    arrange(desc(ObservedN), P_Value, FeatureID, .by_group = TRUE) %>%
    slice_head(n = 1L) %>%
    ungroup()
  ranks <- selected$ModeratedT
  names(ranks) <- selected$GeneSymbol
  sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
}

profile_list <- list()
summary_list <- list()
index <- 1L
for (i in seq_len(nrow(panel_spec))) {
  spec <- panel_spec[i, ]
  for (j in seq_len(nrow(model_spec))) {
    model <- model_spec[j, ]
    for (tissue in c("Adipose", "Muscle")) {
      ranks <- prepare_ranks(model$Marker, model$State, tissue)
      profile <- running_score(ranks, hallmark_sets[[spec$PathwayID]])
      pathway_row <- analysis$PathwayResults %>%
        filter(
          Database == "Hallmark", Marker == model$Marker,
          State == model$State, Tissue == tissue,
          PathwayID == spec$PathwayID
        ) %>%
        slice_head(n = 1L)
      if (nrow(pathway_row) != 1L) {
        stop(
          "Missing pathway result for panel ", spec$Panel, " / ",
          model$Marker, " / ", model$State, " / ", tissue
        )
      }
      profile_list[[index]] <- profile %>%
        mutate(
          Panel = spec$Panel,
          Marker = model$Marker,
          MarkerLabel = model$MarkerLabel,
          State = model$State,
          StateLabel = model$StateLabel,
          PathwayID = spec$PathwayID,
          Pathway = spec$Pathway,
          Tissue = tissue,
          NES = pathway_row$NES,
          BH_FDR = pathway_row$BH_FDR
        )
      summary_list[[index]] <- pathway_row %>%
        transmute(
          Panel = spec$Panel,
          Marker,
          MarkerLabel = model$MarkerLabel,
          State,
          StateLabel = model$StateLabel,
          Tissue, PathwayID, Pathway,
          NES, BH_FDR, NGenes,
          Significant = BH_FDR < 0.05
        )
      index <- index + 1L
    }
  }
}

profile_data <- bind_rows(profile_list) %>%
  mutate(
    MarkerLabel = factor(MarkerLabel, levels = c("cOC", "tOC")),
    StateLabel = factor(StateLabel, levels = c("Pre-exercise", "Exercise response")),
    Tissue = factor(Tissue, levels = c("Adipose", "Muscle"))
  )
summary_data <- bind_rows(summary_list) %>%
  mutate(
    MarkerLabel = factor(MarkerLabel, levels = c("cOC", "tOC")),
    StateLabel = factor(StateLabel, levels = c("Pre-exercise", "Exercise response")),
    Tissue = factor(Tissue, levels = c("Adipose", "Muscle")),
    Significance = vapply(BH_FDR, significance_symbol, character(1)),
    Annotation = paste0(
      Tissue, ": NES ", sprintf("%.2f", NES),
      "; BH-FDR ", vapply(BH_FDR, format_fdr, character(1)),
      if_else(nzchar(Significance), paste0(" ", Significance), "")
    )
  )

tissue_colors <- c("Adipose" = "#D98C9B", "Muscle" = "#6FA7C5")

theme_atm <- theme_classic(base_family = "Arial", base_size = 7) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = 6),
    axis.title = element_text(colour = "black", size = 7, face = "bold"),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.ticks = element_line(colour = "black", linewidth = 0.25),
    plot.title = element_text(size = 8.2, face = "bold", hjust = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 5.8),
    plot.margin = margin(3, 3, 3, 3, unit = "mm")
  )

make_profile_plot <- function(panel_code) {
  spec <- panel_spec %>% filter(Panel == panel_code)
  data <- profile_data %>% filter(Panel == panel_code)
  annotations <- summary_data %>%
    filter(Panel == panel_code) %>%
    arrange(MarkerLabel, StateLabel, Tissue) %>%
    mutate(VJust = if_else(Tissue == "Adipose", 1.25, 2.75))

  curve_min <- min(data$RunningES, na.rm = TRUE)
  curve_max <- max(data$RunningES, na.rm = TRUE)
  curve_span <- max(curve_max - curve_min, 0.35)
  track_height <- 0.052 * curve_span
  track_base <- curve_min - 2.25 * track_height

  hit_tracks <- data %>%
    filter(IsHit) %>%
    mutate(
      TrackY = if_else(
        Tissue == "Adipose",
        track_base + track_height,
        track_base
      ),
      TrackYEnd = TrackY + track_height
    )

  ggplot(data, aes(RankPercent, RunningES, colour = Tissue)) +
    geom_hline(yintercept = 0, colour = "#666666", linewidth = 0.22) +
    geom_line(linewidth = 0.55) +
    geom_segment(
      data = hit_tracks,
      aes(
        x = RankPercent, xend = RankPercent,
        y = TrackY, yend = TrackYEnd,
        colour = Tissue
      ),
      inherit.aes = FALSE,
      linewidth = 0.28, alpha = 0.95, lineend = "butt"
    ) +
    geom_text(
      data = annotations,
      aes(
        x = 98, y = Inf, label = Annotation,
        colour = Tissue, vjust = VJust
      ),
      inherit.aes = FALSE,
      hjust = 1, size = 1.65, show.legend = FALSE
    ) +
    facet_grid(MarkerLabel ~ StateLabel) +
    scale_colour_manual(values = tissue_colors, breaks = c("Adipose", "Muscle")) +
    scale_x_continuous(
      breaks = c(0, 50, 100),
      limits = c(0, 100),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      limits = c(track_base, curve_max + 0.10 * curve_span),
      expand = expansion(mult = c(0.005, 0.01))
    ) +
    labs(
      title = spec$Title,
      x = "Ranked protein list (%)",
      y = "Running enrichment score"
    ) +
    theme_atm +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(
        fill = "white", colour = "black", linewidth = 0.25
      ),
      strip.text = element_text(size = 6.3, face = "bold"),
      panel.spacing.x = grid::unit(2.2, "mm"),
      panel.spacing.y = grid::unit(0.8, "mm"),
      legend.position = "right",
      legend.justification = "center"
    )
}

p_g <- make_profile_plot("G")
p_h <- make_profile_plot("H")

ggsave(
  file.path(output_dir, paste0(panel_spec$FileStem[panel_spec$Panel == "G"], ".pdf")),
  p_g, width = 190, height = 90, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, paste0(panel_spec$FileStem[panel_spec$Panel == "H"], ".pdf")),
  p_h, width = 190, height = 90, units = "mm", device = cairo_pdf
)

combined <- p_g / p_h +
  plot_layout(ncol = 1, guides = "collect") +
  plot_annotation(tag_levels = list(c("G", "H"))) &
  theme(legend.position = "right")

ggsave(
  file.path(output_dir, "FigS06GH_Selected_Hallmark_GSEA_Profiles.pdf"),
  combined, width = 190, height = 180, units = "mm", device = cairo_pdf
)

write_xlsx(
  list(
    GSEA_Summary = summary_data %>%
      mutate(Tissue = as.character(Tissue)) %>%
      select(-Annotation, -Significance),
    Running_Profile_Data = profile_data %>%
      mutate(Tissue = as.character(Tissue))
  ),
  file.path(output_dir, "FigS06GH_Hallmark_GSEA_Source_Data.xlsx")
)

writeLines(c(
  "# Figure S6G-H selected Hallmark GSEA profiles",
  "",
  paste0("Input: fig06/analysis/", analysis_id, "/Fig06_OC_BodyFatAdjusted_Hallmark.rds"),
  "S6G: oxidative phosphorylation across tOC and cOC pre-exercise and exercise-response models.",
  "S6H: heme metabolism across tOC and cOC pre-exercise and exercise-response models.",
  "Protein lists are ranked by the limma moderated t statistic from the body-fat-adjusted, role-adjusted, family-blocked model.",
  "Curves display the running enrichment score; the two compact tracks show Hallmark-member positions in the adipose and muscle ranked lists.",
  "Reported pathway statistics are the frozen fgseaMultilevel NES and BH-FDR values.",
  "These profiles visualize pathway-level rankings and do not imply feature-level BH-FDR significance."
), file.path(output_dir, "README.md"))

writeLines(capture.output(sessionInfo()), file.path(output_dir, "SESSION_INFO.txt"))

message("Figure S6G-H outputs written to: ", output_dir)
