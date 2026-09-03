# Compact Figure 4B Hallmark GSEA running-enrichment profiles.

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

p_load(dplyr, ggplot2, openxlsx, patchwork, scales, tidyr)

analysis_id <- "Fig04B_Core_Hallmark_GSEA_Profiles_2026-08-28"
analysis_rds <- path_project(
  "fig04", "analysis", "Fig04S04_Complete_Table_Analysis_2026-08-27",
  "molecular_pathway", "Fig04S04_Molecular_Pathway_Tables.rds"
)
hallmark_file <- path_project(
  "fig04", "figS4", "Analysis_Inputs", "FigS04_Hallmark_GSEA_Analysis.rds"
)
output_dir <- path_project("fig04", "fig04B")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(analysis_rds, hallmark_file)
if (!all(file.exists(required_files))) {
  stop("Missing required input: ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}

analysis_object <- readRDS(analysis_rds)
feature_results <- analysis_object$FeatureAssociations %>%
  filter(
    ModelTier == "Primary",
    Dataset %in% c("Adipose_Proteomics", "Muscle_Proteomics"),
    State %in% c("PreExercise", "Post3hAbsolute", "ExerciseChange"),
    Predictor %in% c("HOMA_IR", "Matsuda_ISI")
  )
pathway_results <- analysis_object$PathwayResults %>%
  filter(
    ModelTier == "Primary",
    Dataset %in% c("Adipose_Proteomics", "Muscle_Proteomics"),
    State %in% c("PreExercise", "Post3hAbsolute", "ExerciseChange"),
    Predictor %in% c("HOMA_IR", "Matsuda_ISI")
  )
hallmark_sets <- readRDS(hallmark_file)$HallmarkSets

selected_pathways <- c(
  "Oxidative phosphorylation" = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "Fatty acid metabolism" = "HALLMARK_FATTY_ACID_METABOLISM",
  "Heme metabolism" = "HALLMARK_HEME_METABOLISM",
  "Myogenesis" = "HALLMARK_MYOGENESIS"
)

running_score <- function(stats, pathway_genes) {
  stats <- sort(stats[is.finite(stats)], decreasing = TRUE)
  hit <- names(stats) %in% pathway_genes
  nh <- sum(hit)
  if (nh == 0L || nh == length(stats)) stop("Invalid pathway membership.")
  hit_weights <- abs(stats) * hit
  tibble(
    Rank = seq_along(stats),
    RankPercent = 100 * Rank / length(stats),
    RunningES = cumsum(
      hit_weights / sum(hit_weights) - (!hit) / (length(stats) - nh)
    ),
    IsHit = hit,
    Gene = names(stats)
  )
}

model_keys <- feature_results %>%
  distinct(Dataset, Tissue, State, StateLabel, Predictor, PredictorLabel) %>%
  arrange(Dataset, Predictor, State)

profile_list <- vector("list", nrow(model_keys) * length(selected_pathways))
k <- 1L
for (i in seq_len(nrow(model_keys))) {
  spec <- model_keys[i, ]
  stats_table <- feature_results %>%
    semi_join(spec, by = c(
      "Dataset", "Tissue", "State", "StateLabel", "Predictor", "PredictorLabel"
    )) %>%
    filter(!is.na(FeatureLabel), nzchar(FeatureLabel), is.finite(ModeratedT)) %>%
    arrange(desc(ObservedN), P_Value, FeatureID) %>%
    distinct(FeatureLabel, .keep_all = TRUE)
  ranks <- setNames(stats_table$ModeratedT, stats_table$FeatureLabel)

  for (pathway_name in names(selected_pathways)) {
    hallmark_id <- selected_pathways[[pathway_name]]
    pathway_row <- pathway_results %>%
      filter(
        Dataset == spec$Dataset,
        State == spec$State,
        Predictor == spec$Predictor,
        HallmarkID == hallmark_id
      ) %>%
      slice_head(n = 1)
    if (nrow(pathway_row) != 1L) stop("Missing pathway statistic: ", pathway_name)

    profile_list[[k]] <- running_score(ranks, hallmark_sets[[hallmark_id]]) %>%
      mutate(
        Pathway = pathway_name,
        HallmarkID = hallmark_id,
        Dataset = spec$Dataset,
        Tissue = spec$Tissue,
        State = spec$State,
        Predictor = spec$Predictor,
        PredictorLabel = spec$PredictorLabel,
        NES = pathway_row$NES,
        P_Value = pathway_row$P_Value,
        BH_FDR = pathway_row$BH_FDR,
        Subjects = pathway_row$Subjects,
        Families = pathway_row$Families
      )
    k <- k + 1L
  }
}

profiles <- bind_rows(profile_list) %>%
  mutate(
    Tissue = factor(Tissue, levels = c("Adipose tissue", "Skeletal muscle")),
    PredictorLabel = factor(PredictorLabel, levels = c("HOMA-IR", "Matsuda ISI")),
    StateShort = factor(
      State,
      levels = c("PreExercise", "Post3hAbsolute", "ExerciseChange"),
      labels = c("Pre-exercise", "Post3h", "Response")
    ),
    FDRTier = case_when(
      BH_FDR < 0.05 ~ "***",
      BH_FDR < 0.10 ~ "**",
      BH_FDR < 0.20 ~ "*",
      TRUE ~ ""
    )
  )

summary_table <- profiles %>%
  group_by(
    Pathway, HallmarkID, Dataset, Tissue, State, StateShort,
    Predictor, PredictorLabel
  ) %>%
  summarise(
    NES = first(NES), P_Value = first(P_Value), BH_FDR = first(BH_FDR),
    FDRTier = first(FDRTier), Subjects = first(Subjects), Families = first(Families),
    .groups = "drop"
  )

tissue_colors <- c("Adipose tissue" = "#D491A0", "Skeletal muscle" = "#74A2BC")

format_fdr <- function(x) {
  sub("e([+-])0+", "e\\1", formatC(x, format = "e", digits = 1))
}

make_profile_plot <- function(pathway_name) {
  d <- profiles %>% filter(Pathway == pathway_name)
  s <- summary_table %>% filter(Pathway == pathway_name)
  curve_min <- min(d$RunningES, na.rm = TRUE)
  curve_max <- max(d$RunningES, na.rm = TRUE)
  curve_span <- max(curve_max - curve_min, 0.4)
  track_height <- 0.038 * curve_span
  track_base <- curve_min - 2.4 * track_height

  hits <- d %>%
    filter(IsHit) %>%
    mutate(
      TrackY = if_else(Tissue == "Adipose tissue", track_base + track_height, track_base),
      TrackYEnd = TrackY + 0.82 * track_height
    )

  annotation_data <- s %>%
    mutate(
      Label = paste0(
        if_else(Tissue == "Adipose tissue", "A: ", "M: "),
        "NES ", sprintf("%.2f", NES), "; FDR ", format_fdr(BH_FDR)
      ),
      LabelY = if_else(
        Tissue == "Adipose tissue",
        curve_max + 0.055 * curve_span,
        curve_max - 0.045 * curve_span
      )
    )

  p <- ggplot(d, aes(RankPercent, RunningES, colour = Tissue)) +
    geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.20) +
    geom_line(linewidth = 0.42) +
    geom_segment(
      data = hits,
      aes(
        x = RankPercent, xend = RankPercent,
        y = TrackY, yend = TrackYEnd, colour = Tissue
      ),
      inherit.aes = FALSE, linewidth = 0.22, alpha = 0.9, lineend = "butt"
    ) +
    geom_text(
      data = annotation_data,
      aes(x = 50, y = LabelY, label = Label, colour = Tissue),
      inherit.aes = FALSE, hjust = 0.5, size = 1.12, show.legend = FALSE
    ) +
    facet_grid(PredictorLabel ~ StateShort) +
    scale_colour_manual(values = tissue_colors, name = NULL) +
    scale_x_continuous(
      breaks = c(0, 100), labels = c("0", "100"),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      limits = c(track_base, curve_max + 0.16 * curve_span),
      breaks = pretty(c(curve_min, curve_max), n = 3),
      expand = expansion(mult = c(0.002, 0.002))
    ) +
    labs(
      title = pathway_name,
      x = "Ranked protein list (%)",
      y = "Running enrichment score",
      caption = NULL
    ) +
    theme_bw(base_family = "Arial", base_size = 6) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
      panel.spacing = unit(0.6, "mm"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),
      strip.text = element_text(size = 4.8, face = "bold", colour = "black"),
      axis.text = element_text(size = 4.3, colour = "black"),
      axis.title = element_text(size = 4.9, face = "bold", colour = "black"),
      axis.ticks = element_line(linewidth = 0.22, colour = "black"),
      plot.title = element_text(size = 6.6, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text = element_text(size = 4.4),
      legend.key.width = unit(5, "mm"),
      legend.key.height = unit(2.5, "mm"),
      plot.margin = margin(1.0, 1.0, 1.0, 1.0, unit = "mm")
    )

  p
}

save_profile <- function(plot, output_stub, width_mm = 65, height_mm = 48) {

  pdf_file <- file.path(output_dir, paste0(output_stub, ".pdf"))
  png_file <- sub("\\.pdf$", ".png", pdf_file)
  ggsave(pdf_file, plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(png_file, plot, width = width_mm, height = height_mm, units = "mm", dpi = 300)
  c(pdf_file, png_file)
}

p_oxphos <- make_profile_plot("Oxidative phosphorylation")
p_fatty_acid <- make_profile_plot("Fatty acid metabolism")
p_heme <- make_profile_plot("Heme metabolism")
p_myogenesis <- make_profile_plot("Myogenesis")

outputs <- c(
  save_profile(p_oxphos, "Fig04B_CoreProfile_Oxidative_Phosphorylation"),
  save_profile(p_fatty_acid, "Fig04B_CoreProfile_Fatty_Acid_Metabolism"),
  save_profile(p_heme, "Fig04B_CoreProfile_Heme_Metabolism"),
  save_profile(p_myogenesis, "Fig04B_CoreProfile_Myogenesis")
)

main_profiles <- (
  (p_oxphos + theme(axis.title.x = element_blank(), legend.position = "none")) /
    p_fatty_acid
) +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")

main_pdf <- file.path(output_dir, "Fig04B_Main_Two_Core_GSEA_Profiles.pdf")
main_png <- sub("\\.pdf$", ".png", main_pdf)
ggsave(main_pdf, main_profiles, width = 65, height = 93, units = "mm", device = cairo_pdf)
ggsave(main_png, main_profiles, width = 65, height = 93, units = "mm", dpi = 300)

supp_profiles <- (
  p_heme +
    (p_myogenesis + theme(axis.title.y = element_blank()))
) +
  plot_layout(nrow = 1, widths = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")

supp_pdf <- file.path(output_dir, "FigS04B_Two_Core_GSEA_Profiles.pdf")
supp_png <- sub("\\.pdf$", ".png", supp_pdf)
ggsave(supp_pdf, supp_profiles, width = 190, height = 70, units = "mm", device = cairo_pdf)
ggsave(supp_png, supp_profiles, width = 190, height = 70, units = "mm", dpi = 300)

outputs <- c(outputs, main_pdf, main_png, supp_pdf, supp_png)

source_file <- file.path(output_dir, "Fig04B_Core_GSEA_Profile_SourceData.xlsx")
write.xlsx(
  list(
    Profile_summary = summary_table %>% mutate(across(where(is.factor), as.character)),
    Running_scores = profiles %>% mutate(across(where(is.factor), as.character))
  ),
  source_file,
  overwrite = TRUE
)

log_file <- file.path(output_dir, "run_log_core_gsea_profiles.txt")
writeLines(c(
  paste0("Analysis ID: ", analysis_id),
  paste0("Analysis input: ", analysis_rds),
  paste0("Hallmark input: ", hallmark_file),
  "Ranking statistic: limma moderated t statistic from the formal primary role-adjusted, family-blocked models.",
  "States: Pre-exercise abundance, Post3h absolute abundance, and Post3h-minus-Pre response.",
  "Profiles: weighted running enrichment scores; curve colors denote tissue.",
  paste0("Outputs: ", paste(c(outputs, source_file), collapse = "; ")),
  "",
  capture.output(sessionInfo())
), log_file)

message("Completed compact Figure 4B GSEA profiles in: ", output_dir)
