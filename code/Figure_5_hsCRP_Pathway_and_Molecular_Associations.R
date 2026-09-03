# Build the current Figure 5C-E and Figure S5 panels from the unified
# pre-exercise, Post3h, and Post3h-minus-Pre hsCRP association audit.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  if (!is.na(.local_script)) file.path(dirname(.local_script), "Fig00_Config.R") else NA_character_
))
.config_file <- .config_candidates[file.exists(.config_candidates)][1]
if (is.na(.config_file)) stop("Cannot locate Fig00_Config.R.")
source(.config_file)
rm(.local_script, .config_candidates, .config_file)

p_load(
  dplyr, ggplot2, openxlsx, patchwork, readr, scales, stringr,
  tibble, tidyr
)

analysis_id <- "Fig05_Current_hsCRP_Hallmark_and_Features_2026-08-29"
formal_rds <- path_project(
  "results", "Unified_Pre_Post_Response_Association_Audit_2026-08-28",
  "Unified_NonMethylation_Pre_Post_Response_Audit.rds"
)
hallmark_rds <- path_project(
  "fig04", "figS4", "Analysis_Inputs", "FigS04_Hallmark_GSEA_Analysis.rds"
)
clinical_file <- path_analysis_ready("Metadata", "Clinical_Observed_Primary.csv")
sample_sheet_file <- path_analysis_ready("Metadata", "Project_Sample_Sheet.csv")

main_dir <- path_project("fig05", "fig5")
supp_dir <- path_project("fig05", "figS5")
table_dir <- path_project("fig05", "Suptable")
analysis_dir <- path_project(
  "fig05", "analysis", "Fig05_Current_hsCRP_Hallmark_2026-08-29"
)
for (directory in c(main_dir, supp_dir, table_dir, analysis_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

required_files <- c(formal_rds, hallmark_rds, clinical_file, sample_sheet_file)
if (!all(file.exists(required_files))) {
  stop("Missing required input: ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}

formal <- readRDS(formal_rds)
hallmark_object <- readRDS(hallmark_rds)
feature_results <- formal$feature_results
pathway_results <- formal$pathway_results
hallmark_sets <- hallmark_object$HallmarkSets

role_levels <- c("Daughter", "Mother", "Father")
role_shapes <- c("Daughter" = 24, "Mother" = 21, "Father" = 22)
body_fat_colors <- c("#6D9FC4", "#F4F1EC", "#D77A72")
tissue_colors <- c("Adipose tissue" = "#D491A0", "Skeletal muscle" = "#74A2BC")
state_levels <- c("PreExercise", "Post3hAbsolute", "ExerciseChange")
state_labels <- c(
  "PreExercise" = "Pre-exercise",
  "Post3hAbsolute" = "Post3h",
  "ExerciseChange" = "Response"
)

rank_normalize <- function(x) {
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  n <- sum(keep)
  if (n >= 3L) out[keep] <- qnorm((rank(x[keep], ties.method = "average") - 0.5) / n)
  out
}

nice_breaks_exact <- function(values, n = 5L) {
  values <- values[is.finite(values)]
  if (length(values) < 2L || diff(range(values)) == 0) return(pretty(values, n = n))
  value_range <- range(values)
  exponent <- floor(log10(diff(value_range) / (n - 1L)))
  candidate_steps <- sort(unique(as.vector(outer(
    c(1, 1.5, 2, 2.5, 5),
    10 ^ seq(exponent - 1L, exponent + 2L),
    `*`
  ))))
  for (step in candidate_steps) {
    start <- floor(value_range[1] / step) * step
    breaks <- start + seq.int(0L, n - 1L) * step
    if (breaks[n] >= value_range[2]) return(breaks)
  }
  pretty(values, n = n)
}

format_probability <- function(x, digits = 2L) {
  if (!is.finite(x)) return("NA")
  if (x < 0.001) {
    return(sub("e([+-])0+", "e\\1", formatC(x, format = "e", digits = digits)))
  }
  formatC(x, format = "f", digits = 3)
}

display_pathway <- function(x) {
  recode(
    x,
    "Mtorc1 signaling" = "mTORC1 signaling",
    "Myc targets v1" = "MYC targets v1",
    "Myc targets v2" = "MYC targets v2",
    "Dna repair" = "DNA repair",
    "Kras signaling dn" = "KRAS signaling down",
    "Kras signaling up" = "KRAS signaling up",
    "Pi3k akt mtor signaling" = "PI3K-AKT-mTOR signaling",
    "Tnfa signaling via nfkb" = "TNFA signaling via NFKB",
    "Uv response dn" = "UV response down",
    "Epithelial mesenchymal transition" = "Epithelial-mesenchymal transition",
    .default = x
  )
}

hscrp_pathways <- pathway_results %>%
  filter(
    Predictor == "hsCRP",
    Dataset %in% c("Adipose_Proteomics", "Muscle_Proteomics"),
    State %in% state_levels
  ) %>%
  mutate(
    Tissue = recode(
      Tissue,
      "Adipose" = "Adipose tissue",
      "Skeletal muscle" = "Skeletal muscle"
    ),
    StateShort = factor(state_labels[State], levels = unname(state_labels)),
    DisplayPathway = display_pathway(Pathway),
    FDRTier = case_when(
      BH_FDR < 0.05 ~ "***",
      BH_FDR < 0.10 ~ "**",
      BH_FDR < 0.20 ~ "*",
      TRUE ~ ""
    ),
    DisplayNES = if_else(BH_FDR < 0.20, NES, NA_real_),
    FDRStrength = if_else(
      BH_FDR < 0.20,
      pmin(-log10(pmax(BH_FDR, 1e-12)), 10),
      NA_real_
    )
  )

core_pathways <- c(
  "Oxidative phosphorylation",
  "Fatty acid metabolism",
  "Adipogenesis",
  "mTORC1 signaling",
  "Epithelial-mesenchymal transition",
  "Coagulation",
  "Complement",
  "MYC targets v1"
)

expanded_pathways <- hscrp_pathways %>%
  group_by(DisplayPathway) %>%
  summarise(MinFDR = min(BH_FDR, na.rm = TRUE), .groups = "drop") %>%
  filter(MinFDR < 0.05) %>%
  arrange(MinFDR) %>%
  pull(DisplayPathway)

make_bubble_plot <- function(
    pathways, title, tissue_layout = c("rows", "columns"), show_caption = TRUE,
    combined_header = FALSE) {
  tissue_layout <- match.arg(tissue_layout)
  complete_grid <- tidyr::crossing(
    Tissue = factor(c("Adipose tissue", "Skeletal muscle"),
                    levels = c("Adipose tissue", "Skeletal muscle")),
    DisplayPathway = factor(pathways, levels = rev(pathways)),
    StateShort = factor(unname(state_labels), levels = unname(state_labels))
  )
  d <- hscrp_pathways %>%
    filter(DisplayPathway %in% pathways) %>%
    mutate(
      Tissue = factor(Tissue, levels = c("Adipose tissue", "Skeletal muscle")),
      DisplayPathway = factor(DisplayPathway, levels = rev(pathways)),
      PredictorLabel = factor("hsCRP", levels = "hsCRP")
    ) %>%
    right_join(complete_grid, by = c("Tissue", "DisplayPathway", "StateShort"))

  d <- d %>% mutate(
    PredictorLabel = factor("hsCRP", levels = "hsCRP"),
    FacetLabel = factor(
      paste0(as.character(Tissue), " (hsCRP)"),
      levels = c("Adipose tissue (hsCRP)", "Skeletal muscle (hsCRP)")
    )
  )
  tissue_facet <- if (tissue_layout == "rows") {
    facet_grid(
      rows = vars(Tissue), cols = vars(PredictorLabel),
      scales = "free_y", space = "free_y"
    )
  } else if (combined_header) {
    facet_grid(cols = vars(FacetLabel))
  } else {
    facet_grid(cols = vars(Tissue, PredictorLabel))
  }

  ggplot(d, aes(StateShort, DisplayPathway)) +
    geom_point(
      aes(fill = DisplayNES, size = FDRStrength),
      shape = 21, colour = "#333333", stroke = 0.22, na.rm = TRUE
    ) +
    geom_text(
      data = d %>% filter(FDRTier != ""),
      aes(label = FDRTier), family = "Arial", fontface = "bold",
      colour = "black", size = 1.85, vjust = 0.42, na.rm = TRUE
    ) +
    tissue_facet +
    scale_fill_gradient2(
      low = "#83A9C7", mid = "white", high = "#D9909E",
      midpoint = 0, limits = c(-4, 4), oob = squish,
      name = "GSEA NES"
    ) +
    scale_size_continuous(
      range = c(2.3, 6.1), limits = c(0.7, 10),
      breaks = c(1, 2, 5, 10), name = expression(-log[10]("BH-FDR"))
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
      title = str_wrap(title, width = 40),
      x = NULL,
      y = NULL,
      caption = if (show_caption) paste0(
        "Stars denote BH-FDR tiers: ***BH-FDR < 0.05, ",
        "**BH-FDR < 0.10, and *BH-FDR < 0.20."
      ) else NULL
    ) +
    theme_bw(base_family = "Arial", base_size = 7) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.28),
      panel.spacing.x = unit(0.7, "mm"),
      panel.spacing.y = unit(0.7, "mm"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.28),
      strip.text = element_text(size = 6.1, face = "bold", colour = "black"),
      axis.text.x = element_text(size = 5.5, colour = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 6.1, colour = "black"),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 8.0, hjust = 0, lineheight = 0.95),
      plot.caption = element_text(size = 5.7, hjust = 0, colour = "black"),
      legend.position = "right",
      legend.justification = "center",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.spacing.y = unit(1.2, "mm"),
      legend.title = element_text(face = "bold", size = 5.4, hjust = 0.5),
      legend.text = element_text(size = 5.0),
      legend.key.height = unit(3.4, "mm"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )
}

main_bubble <- make_bubble_plot(
  core_pathways,
  "Core Hallmark pathway associations across molecular states",
  tissue_layout = "rows"
) +
  labs(caption = "BH-FDR: *** < 0.05; ** < 0.10;\n* < 0.20.") +
  guides(
    fill = guide_colourbar(
      title.position = "top", title.hjust = 0.5, order = 1,
      direction = "vertical", barheight = unit(12, "mm"),
      barwidth = unit(2.6, "mm")
    ),
    size = guide_legend(
      title.position = "top", title.hjust = 0.5, order = 2,
      direction = "vertical", ncol = 1,
      keyheight = unit(2.5, "mm"), keywidth = unit(4.0, "mm")
    )
  ) +
  theme(
    legend.position = "right", legend.justification = "center",
    legend.box = "vertical", legend.box.just = "center",
    legend.title = element_text(face = "bold", size = 5.0, hjust = 0.5),
    legend.text = element_text(size = 4.7)
  )
supp_bubble <- make_bubble_plot(
  expanded_pathways,
  "Hallmark pathways with BH-FDR < 0.20 in at least one hsCRP model",
  tissue_layout = "columns",
  combined_header = TRUE
)

ggsave(
  file.path(main_dir, "Fig05C_Core_Hallmark_Pathway_Associations.pdf"),
  main_bubble, width = 82, height = 78, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(analysis_dir, "QA_Fig05C_Core_Hallmark_Pathway_Associations.png"),
  main_bubble, width = 82, height = 78, units = "mm", dpi = 300, bg = "white"
)
supp_bubble_centered <- plot_spacer() + supp_bubble + plot_spacer() +
  plot_layout(widths = c(14, 162, 14))
ggsave(
  file.path(supp_dir, "FigS05A_Expanded_Hallmark_Pathway_Associations.pdf"),
  supp_bubble_centered, width = 190, height = 88, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(analysis_dir, "QA_FigS05A_Expanded_Hallmark_Pathway_Associations.png"),
  supp_bubble_centered, width = 190, height = 88, units = "mm", dpi = 300, bg = "white"
)

# Running-enrichment profiles -------------------------------------------------
hscrp_features <- feature_results %>%
  filter(
    Predictor == "hsCRP",
    Dataset %in% c("Adipose_Proteomics", "Muscle_Proteomics"),
    State %in% state_levels,
    !is.na(GeneSymbol), nzchar(GeneSymbol), is.finite(ModeratedT)
  )

model_keys <- hscrp_features %>%
  distinct(Dataset, Tissue, State)

running_score <- function(stats, pathway_genes) {
  stats <- sort(stats[is.finite(stats)], decreasing = TRUE)
  hit <- names(stats) %in% pathway_genes
  n_hit <- sum(hit)
  if (n_hit == 0L || n_hit == length(stats)) return(tibble())
  weights <- abs(stats) * hit
  tibble(
    Rank = seq_along(stats),
    RankPercent = 100 * Rank / length(stats),
    RunningES = cumsum(weights / sum(weights) - (!hit) / (length(stats) - n_hit)),
    IsHit = hit,
    Feature = names(stats)
  )
}

profile_specs <- tribble(
  ~DisplayPathway, ~HallmarkID, ~FigureSet,
  "Oxidative phosphorylation", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "Main",
  "Epithelial-mesenchymal transition", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "Main",
  "Fatty acid metabolism", "HALLMARK_FATTY_ACID_METABOLISM", "Supplementary",
  "Adipogenesis", "HALLMARK_ADIPOGENESIS", "Supplementary",
  "Coagulation", "HALLMARK_COAGULATION", "Supplementary",
  "MYC targets v1", "HALLMARK_MYC_TARGETS_V1", "Supplementary"
)

profile_list <- list()
k <- 1L
for (i in seq_len(nrow(model_keys))) {
  key <- model_keys[i, ]
  ranks_table <- hscrp_features %>%
    filter(Dataset == key$Dataset, State == key$State) %>%
    mutate(GeneSymbol = str_split(GeneSymbol, ";", simplify = TRUE)[, 1]) %>%
    filter(!is.na(GeneSymbol), nzchar(GeneSymbol)) %>%
    arrange(P_Value, desc(ObservedN), FeatureID) %>%
    distinct(GeneSymbol, .keep_all = TRUE)
  ranks <- setNames(ranks_table$ModeratedT, ranks_table$GeneSymbol)
  for (j in seq_len(nrow(profile_specs))) {
    spec <- profile_specs[j, ]
    path_row <- hscrp_pathways %>%
      filter(
        Dataset == key$Dataset,
        State == key$State,
        DisplayPathway == spec$DisplayPathway
      ) %>%
      slice_head(n = 1L)
    if (nrow(path_row) != 1L || !spec$HallmarkID %in% names(hallmark_sets)) next
    curve <- running_score(ranks, hallmark_sets[[spec$HallmarkID]])
    if (!nrow(curve)) next
    profile_list[[k]] <- curve %>%
      mutate(
        DisplayPathway = spec$DisplayPathway,
        HallmarkID = spec$HallmarkID,
        FigureSet = spec$FigureSet,
        Dataset = key$Dataset,
        Tissue = if_else(key$Dataset == "Adipose_Proteomics", "Adipose tissue", "Skeletal muscle"),
        State = key$State,
        StateShort = factor(state_labels[key$State], levels = unname(state_labels)),
        PredictorLabel = factor("hsCRP", levels = "hsCRP"),
        NES = path_row$NES,
        P_Value = path_row$P_Value,
        BH_FDR = path_row$BH_FDR,
        Subjects = path_row$Subjects,
        Families = path_row$Families
      )
    k <- k + 1L
  }
}
profile_data <- bind_rows(profile_list) %>%
  mutate(Tissue = factor(Tissue, levels = c("Adipose tissue", "Skeletal muscle")))

profile_summary <- profile_data %>%
  distinct(
    DisplayPathway, HallmarkID, FigureSet, Dataset, Tissue, State, StateShort,
    PredictorLabel,
    NES, P_Value, BH_FDR, Subjects, Families
  )

make_profile_plot <- function(pathway_name) {
  d <- profile_data %>% filter(DisplayPathway == pathway_name)
  s <- profile_summary %>% filter(DisplayPathway == pathway_name)
  curve_min <- min(d$RunningES, na.rm = TRUE)
  curve_max <- max(d$RunningES, na.rm = TRUE)
  curve_span <- max(curve_max - curve_min, 0.4)
  track_height <- 0.045 * curve_span
  track_base <- curve_min - 2.3 * track_height
  y_breaks <- nice_breaks_exact(d$RunningES, n = 3L)
  y_upper <- max(curve_max + 0.16 * curve_span, max(y_breaks))

  hits <- d %>%
    filter(IsHit) %>%
    mutate(
      TrackY = if_else(Tissue == "Adipose tissue", track_base + track_height, track_base),
      TrackYEnd = TrackY + 0.82 * track_height
    )
  annotations <- s %>%
    mutate(
      Label = paste0(
        if_else(Tissue == "Adipose tissue", "A: ", "M: "),
        "NES ", sprintf("%.2f", NES), "; FDR ",
        vapply(BH_FDR, format_probability, character(1))
      ),
      LabelY = if_else(
        Tissue == "Adipose tissue",
        curve_max + 0.055 * curve_span,
        curve_max - 0.045 * curve_span
      )
    )

  ggplot(d, aes(RankPercent, RunningES, colour = Tissue)) +
    geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.20) +
    geom_line(linewidth = 0.42) +
    geom_segment(
      data = hits,
      aes(
        x = RankPercent, xend = RankPercent,
        y = TrackY, yend = TrackYEnd, colour = Tissue
      ),
      inherit.aes = FALSE, linewidth = 0.22, alpha = 0.90, lineend = "butt"
    ) +
    geom_text(
      data = annotations,
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
      limits = c(track_base, y_upper),
      breaks = y_breaks,
      expand = expansion(mult = c(0.002, 0.002))
    ) +
    labs(
      title = pathway_name,
      x = "Ranked protein list (%)",
      y = "Running enrichment score"
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
      plot.margin = margin(1, 1, 1, 1, unit = "mm")
    )
}

p_oxphos <- make_profile_plot("Oxidative phosphorylation")
p_emt <- make_profile_plot("Epithelial-mesenchymal transition")
p_fatty <- make_profile_plot("Fatty acid metabolism")
p_adipogenesis <- make_profile_plot("Adipogenesis")
p_coag <- make_profile_plot("Coagulation")
p_myc <- make_profile_plot("MYC targets v1")

main_profiles <- (p_oxphos + theme(axis.title.x = element_blank(), legend.position = "none")) /
  p_emt + plot_layout(heights = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")
supp_profiles <- wrap_plots(
  p_fatty + theme(axis.title.x = element_blank(), legend.position = "none"),
  p_adipogenesis + theme(
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    legend.position = "none"
  ),
  p_coag,
  p_myc + theme(axis.title.y = element_blank()),
  ncol = 2, guides = "collect"
) & theme(legend.position = "bottom")

ggsave(
  file.path(main_dir, "Fig05D_Core_GSEA_Profiles.pdf"),
  main_profiles, width = 118, height = 78, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(analysis_dir, "QA_Fig05D_Core_GSEA_Profiles.png"),
  main_profiles, width = 118, height = 78, units = "mm", dpi = 300, bg = "white"
)
ggsave(
  file.path(supp_dir, "FigS05C_Additional_GSEA_Profiles.pdf"),
  supp_profiles, width = 190, height = 93, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(analysis_dir, "QA_FigS05C_Additional_GSEA_Profiles.png"),
  supp_profiles, width = 190, height = 93, units = "mm", dpi = 300, bg = "white"
)

# Exact 190-mm one-row composition for direct Figure 5C-D placement.
main_cd <- main_bubble + main_profiles +
  plot_layout(widths = c(82, 108))
ggsave(
  file.path(main_dir, "Fig05CD_Core_Hallmark_and_GSEA_Profiles.pdf"),
  main_cd, width = 190, height = 78, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(analysis_dir, "QA_Fig05CD_Core_Hallmark_and_GSEA_Profiles.png"),
  main_cd, width = 190, height = 78, units = "mm", dpi = 300, bg = "white"
)

# Feature-level panels --------------------------------------------------------
clinical <- read_csv(clinical_file, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = factor(Role, levels = role_levels),
    BodyFat = suppressWarnings(as.numeric(Fat_p_G)),
    hsCRP = suppressWarnings(as.numeric(CRP_84m_G))
  ) %>%
  distinct(Clinical_Subject_ID, .keep_all = TRUE) %>%
  mutate(hsCRP_Z = rank_normalize(if_else(hsCRP > 0, hsCRP, NA_real_)))

sample_sheet <- read_csv(sample_sheet_file, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    SampleID = tolower(as.character(SampleID)),
    Role = factor(Role, levels = role_levels),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort)
  ) %>%
  filter(AnalysisSet == "Primary", Is_Primary_Cohort)

candidates <- tribble(
  ~Panel, ~Dataset, ~FeatureID, ~FeatureLabel, ~State, ~StateLabel, ~YAxisLabel,
  "A", "Adipose_Proteomics", "Q93062", "RBPMS", "ExerciseChange", "Response",
  "Role-adjusted molecular value",
  "B", "Adipose_Proteomics", "Q9Y6D6", "ARFGEF1", "ExerciseChange", "Response",
  "Role-adjusted molecular value",
  "C", "Adipose_Proteomics", "C4AMC7;Q6VEQ5", "WASH3P/WASH2P", "ExerciseChange", "Response",
  "Role-adjusted molecular value"
)

adipose_matrix <- read_analysis_matrix(path_analysis_ready("Analysis_Adipose_Proteomics.csv"))

extract_candidate <- function(spec) {
  formal_row <- feature_results %>%
    filter(
      Predictor == "hsCRP", Dataset == spec$Dataset,
      State == spec$State, FeatureID == spec$FeatureID
    ) %>%
    slice_head(n = 1L)
  if (nrow(formal_row) != 1L) stop("Missing formal feature result: ", spec$FeatureLabel)
  if (!spec$FeatureID %in% rownames(adipose_matrix)) stop("Missing matrix feature: ", spec$FeatureID)

  available <- sample_sheet %>%
    filter(Dataset == spec$Dataset, Timepoint %in% c("Pre", "Post3h"),
           SampleID %in% colnames(adipose_matrix)) %>%
    inner_join(clinical, by = c("FamilyID", "Clinical_Subject_ID", "Role")) %>%
    distinct(Clinical_Subject_ID, Timepoint, .keep_all = TRUE)

  if (spec$State %in% c("PreExercise", "Post3hAbsolute")) {
    target_timepoint <- if (spec$State == "PreExercise") "Pre" else "Post3h"
    metadata <- available %>% filter(Timepoint == target_timepoint) %>%
      arrange(FamilyID, Role, Clinical_Subject_ID)
    outcome <- as.numeric(adipose_matrix[spec$FeatureID, metadata$SampleID])
  } else {
    metadata <- available %>%
      select(FamilyID, Clinical_Subject_ID, Role, BodyFat, hsCRP, hsCRP_Z, Timepoint, SampleID) %>%
      pivot_wider(names_from = Timepoint, values_from = SampleID) %>%
      filter(!is.na(Pre), !is.na(Post3h)) %>%
      arrange(FamilyID, Role, Clinical_Subject_ID)
    pre <- as.numeric(adipose_matrix[spec$FeatureID, metadata$Pre])
    post <- as.numeric(adipose_matrix[spec$FeatureID, metadata$Post3h])
    outcome <- post - pre
  }

  keep <- is.finite(outcome) & is.finite(metadata$hsCRP_Z) &
    is.finite(metadata$BodyFat) & !is.na(metadata$Role)
  metadata <- droplevels(metadata[keep, , drop = FALSE])
  outcome <- outcome[keep]
  role_fit <- lm(outcome ~ Role, data = metadata)
  role_adjusted <- residuals(role_fit) + mean(outcome)

  stats <- formal_row %>%
    transmute(
      Panel = spec$Panel,
      FeatureLabel = spec$FeatureLabel,
      State = spec$State,
      StateLabel = spec$StateLabel,
      Effect = Effect,
      ModeratedT = ModeratedT,
      P_Value = P_Value,
      BH_FDR = BH_FDR,
      Subjects = Subjects,
      Families = Families,
      Intercept = mean(role_adjusted) - Effect * mean(metadata$hsCRP_Z),
      Annotation = paste0(
        intToUtf8(946), " = ", sprintf("%+.3f", Effect),
        "\nP = ", format_probability(P_Value),
        "\nBH-FDR = ", format_probability(BH_FDR)
      )
    )

  data <- metadata %>%
    mutate(
      Panel = spec$Panel,
      FeatureLabel = spec$FeatureLabel,
      State = spec$State,
      StateLabel = spec$StateLabel,
      YAxisLabel = spec$YAxisLabel,
      MolecularValue = outcome,
      RoleAdjustedValue = role_adjusted
    )
  list(data = data, stats = stats)
}

extracted <- lapply(seq_len(nrow(candidates)), function(i) extract_candidate(candidates[i, ]))
plot_data <- bind_rows(lapply(extracted, `[[`, "data"))
model_stats <- bind_rows(lapply(extracted, `[[`, "stats"))

make_feature_panel <- function(panel_id) {
  d <- plot_data %>% filter(Panel == panel_id)
  s <- model_stats %>% filter(Panel == panel_id)
  y_breaks <- nice_breaks_exact(d$RoleAdjustedValue, n = 5L)
  ggplot(d, aes(hsCRP_Z, RoleAdjustedValue)) +
    geom_abline(
      intercept = s$Intercept, slope = s$Effect,
      linewidth = 0.55, colour = "#3F4A54"
    ) +
    geom_point(
      aes(fill = BodyFat, shape = Role),
      size = 2.25, stroke = 0.45, colour = "#30353A", alpha = 0.95
    ) +
    scale_shape_manual(values = role_shapes, drop = FALSE) +
    scale_fill_gradientn(
      colours = body_fat_colors,
      limits = range(clinical$BodyFat, na.rm = TRUE),
      oob = squish,
      name = "Body fat (%)"
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      limits = range(y_breaks),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = s$FeatureLabel,
      subtitle = paste0("Adipose proteomics, ", s$StateLabel, "\n", s$Annotation),
      x = NULL,
      y = if (panel_id == "A") d$YAxisLabel[1] else NULL
    ) +
    theme_classic(base_family = "Arial", base_size = 7.2) +
    theme(
      plot.title = element_text(face = "bold", size = 7.8, hjust = 0, colour = "black"),
      plot.subtitle = element_text(size = 6.2, colour = "black", hjust = 0, lineheight = 1.05),
      axis.title = element_text(face = "bold", size = 6.8, colour = "black"),
      axis.text = element_text(size = 6.3, colour = "black"),
      axis.line = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.4),
      legend.title = element_text(face = "bold", size = 6.5),
      legend.text = element_text(size = 6.2),
      legend.key.height = unit(3.4, "mm"),
      legend.key.width = unit(3.4, "mm"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

feature_panels <- lapply(candidates$Panel, make_feature_panel)
feature_figure <- wrap_plots(feature_panels, nrow = 1, guides = "collect") +
  plot_annotation(
    title = "E  Selected 3-h post-exercise adipose protein responses associated with baseline hsCRP",
    caption = "Baseline hsCRP (rank-normalized SD)",
    theme = theme(
      text = element_text(family = "Arial", colour = "black"),
      plot.title = element_text(face = "bold", size = 9.2, hjust = 0),
      plot.caption = element_text(face = "bold", size = 6.8, hjust = 0.46, margin = margin(t = 2))
    )
  ) & theme(legend.position = "right")

feature_figure_centered <- plot_spacer() + feature_figure + plot_spacer() +
  plot_layout(widths = c(32, 126, 32))

ggsave(
  file.path(main_dir, "Fig05E_hsCRP_Molecular_Associations.pdf"),
  feature_figure_centered, width = 190, height = 54, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(analysis_dir, "QA_Fig05E_hsCRP_Molecular_Associations.png"),
  feature_figure_centered, width = 190, height = 54, units = "mm", dpi = 500, bg = "white"
)

# Source data and supplementary table ----------------------------------------
feature_candidate_ids <- c("Q9H2H8", "Q93062", "Q9Y6D6", "C4AMC7;Q6VEQ5")
feature_candidate_results <- feature_results %>%
  filter(
    Predictor == "hsCRP", Dataset == "Adipose_Proteomics",
    FeatureID %in% feature_candidate_ids,
    State %in% state_levels
  ) %>%
  mutate(
    FeatureLabel = case_when(
      FeatureID == "Q9H2H8" ~ "PPIL3",
      FeatureID == "Q93062" ~ "RBPMS",
      FeatureID == "Q9Y6D6" ~ "ARFGEF1",
      FeatureID == "C4AMC7;Q6VEQ5" ~ "WASH3P/WASH2P",
      TRUE ~ coalesce(GeneSymbol, FeatureID)
    ),
    DisplayedInFigure =
      FeatureID %in% c("Q93062", "Q9Y6D6", "C4AMC7;Q6VEQ5") &
      State == "ExerciseChange",
    Interpretation = case_when(
      FeatureID == "Q9H2H8" & State == "PreExercise" ~ "Strict pre-exercise feature-level association",
      FeatureID %in% c("Q93062", "Q9Y6D6") & State == "ExerciseChange" ~
        "Exploratory response association (BH-FDR < 0.20)",
      FeatureID == "C4AMC7;Q6VEQ5" & State == "ExerciseChange" ~
        "Exploratory shared-peptide/pseudogene protein group (BH-FDR < 0.20)",
      TRUE ~ "Context model for the same candidate"
    ),
    FigureLocation = case_when(
      FeatureID %in% c("Q9H2H8", "Q93062") ~ "Figure S5B",
      FeatureID != "Q9H2H8" & State == "ExerciseChange" ~ "Figure 5E",
      TRUE ~ "Table S10"
    )
  ) %>%
  arrange(match(FeatureID, feature_candidate_ids), match(State, state_levels)) %>%
  select(
    FeatureLabel, FeatureID, Tissue, State, StateLabel, Contrast, Effect, ModeratedT,
    P_Value, BH_FDR, Subjects, Families, FamilyCorrelation, DisplayedInFigure,
    Interpretation, FigureLocation
  )

# Figure S5B: two proteins across pre-exercise, Post3h, and response states --
s5b_candidates <- tribble(
  ~Panel, ~Dataset, ~FeatureID, ~FeatureLabel, ~State, ~StateLabel, ~YAxisLabel,
  "A", "Adipose_Proteomics", "Q9H2H8", "PPIL3", "PreExercise", "Pre-exercise", "Role-adjusted protein value",
  "B", "Adipose_Proteomics", "Q9H2H8", "PPIL3", "Post3hAbsolute", "Post3h", "Role-adjusted protein value",
  "C", "Adipose_Proteomics", "Q9H2H8", "PPIL3", "ExerciseChange", "Response", "Role-adjusted protein change",
  "D", "Adipose_Proteomics", "Q93062", "RBPMS", "PreExercise", "Pre-exercise", "Role-adjusted protein value",
  "E", "Adipose_Proteomics", "Q93062", "RBPMS", "Post3hAbsolute", "Post3h", "Role-adjusted protein value",
  "F", "Adipose_Proteomics", "Q93062", "RBPMS", "ExerciseChange", "Response", "Role-adjusted protein change"
)

s5b_extracted <- lapply(
  seq_len(nrow(s5b_candidates)),
  function(i) extract_candidate(s5b_candidates[i, ])
)
s5b_plot_data <- bind_rows(lapply(s5b_extracted, `[[`, "data"))
s5b_model_stats <- bind_rows(lapply(s5b_extracted, `[[`, "stats"))
s5b_x_limit <- ceiling(max(abs(s5b_plot_data$hsCRP_Z), na.rm = TRUE) * 2) / 2

make_s5b_panel <- function(panel_id) {
  d <- s5b_plot_data %>% filter(Panel == panel_id)
  s <- s5b_model_stats %>% filter(Panel == panel_id)
  y_breaks <- nice_breaks_exact(d$RoleAdjustedValue, n = 4L)
  ggplot(d, aes(hsCRP_Z, RoleAdjustedValue)) +
    geom_abline(
      intercept = s$Intercept, slope = s$Effect,
      linewidth = 0.48, colour = "#3F4A54"
    ) +
    geom_point(
      aes(fill = BodyFat, shape = Role),
      size = 1.75, stroke = 0.35, colour = "#30353A", alpha = 0.93
    ) +
    scale_shape_manual(values = role_shapes, drop = FALSE) +
    scale_fill_gradientn(
      colours = body_fat_colors,
      limits = range(clinical$BodyFat, na.rm = TRUE),
      oob = squish,
      name = "Body fat (%)"
    ) +
    scale_x_continuous(
      breaks = c(-2, 0, 2), limits = c(-s5b_x_limit, s5b_x_limit),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = y_breaks, limits = range(y_breaks),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = paste0(s$FeatureLabel, " (", s$StateLabel, ")"),
      subtitle = paste0(
        intToUtf8(946), " = ", sprintf("%+.2f", s$Effect),
        "\nP = ", format_probability(s$P_Value),
        "; FDR = ", format_probability(s$BH_FDR)
      ),
      x = NULL,
      y = if (panel_id %in% c("A", "D")) d$YAxisLabel[1] else NULL
    ) +
    theme_classic(base_family = "Arial", base_size = 6.5) +
    theme(
      plot.title = element_text(face = "bold", size = 6.9, hjust = 0, colour = "black"),
      plot.subtitle = element_text(size = 5.1, colour = "black", hjust = 0, lineheight = 1.02),
      axis.title = element_text(face = "bold", size = 5.8, colour = "black"),
      axis.text = element_text(size = 5.5, colour = "black"),
      axis.line = element_line(linewidth = 0.34, colour = "black"),
      axis.ticks = element_line(linewidth = 0.34, colour = "black"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.34),
      legend.title = element_text(face = "bold", size = 5.7),
      legend.text = element_text(size = 5.3),
      legend.key.height = unit(3.1, "mm"),
      legend.key.width = unit(3.1, "mm"),
      plot.margin = margin(2, 1.5, 2, 1.5, unit = "mm")
    )
}

s5b_panels <- lapply(s5b_candidates$Panel, make_s5b_panel)
s5b_plot <- wrap_plots(s5b_panels, ncol = 3, nrow = 2, guides = "collect") +
  plot_annotation(
    title = "Selected adipose-protein associations with baseline hsCRP across molecular states",
    caption = "Baseline hsCRP (rank-normalized SD)",
    theme = theme(
      text = element_text(family = "Arial", colour = "black"),
      plot.title = element_text(face = "bold", size = 8.5, hjust = 0),
      plot.caption = element_text(face = "bold", size = 6.0, hjust = 0.45, margin = margin(t = 1.5))
    )
  ) & theme(legend.position = "right", legend.justification = "center")

s5b_plot_centered <- plot_spacer() + s5b_plot + plot_spacer() +
  plot_layout(widths = c(20, 150, 20))

ggsave(
  file.path(supp_dir, "FigS05B_Selected_hsCRP_Protein_Associations.pdf"),
  s5b_plot_centered, width = 190, height = 102, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(analysis_dir, "QA_FigS05B_Selected_hsCRP_Protein_Associations.png"),
  s5b_plot_centered, width = 190, height = 102, units = "mm", dpi = 450, bg = "white"
)
write_csv(
  s5b_plot_data %>%
    select(
      Panel, FeatureLabel, State, StateLabel, FamilyID, Clinical_Subject_ID,
      Role, BodyFat, hsCRP, hsCRP_Z, MolecularValue, RoleAdjustedValue
    ),
  file.path(analysis_dir, "FigS05B_Selected_hsCRP_Protein_Associations_SourceData.csv")
)
write_csv(
  s5b_model_stats,
  file.path(analysis_dir, "FigS05B_Selected_hsCRP_Protein_Associations_Statistics.csv")
)

pathway_table <- hscrp_pathways %>%
  mutate(
    FigureLocation = case_when(
      DisplayPathway %in% c("Oxidative phosphorylation", "Epithelial-mesenchymal transition") ~
        "Figures 5C and 5D; Figure S5A",
      DisplayPathway %in% c(
        "Fatty acid metabolism", "Adipogenesis", "Coagulation", "MYC targets v1"
      ) ~ "Figure 5C; Figures S5A and S5B",
      DisplayPathway %in% core_pathways ~ "Figure 5C; Figure S5A",
      DisplayPathway %in% expanded_pathways ~ "Figure S5A",
      TRUE ~ "Table only"
    ),
    FDRTier = case_when(
      BH_FDR < 0.05 ~ "BH-FDR < 0.05",
      BH_FDR < 0.10 ~ "0.05 <= BH-FDR < 0.10",
      BH_FDR < 0.20 ~ "0.10 <= BH-FDR < 0.20",
      TRUE ~ "BH-FDR >= 0.20"
    )
  ) %>%
  arrange(Tissue, State, BH_FDR, DisplayPathway) %>%
  select(
    DisplayPathway, PathwayID, Tissue, StateLabel, Contrast, Size, NES, P_Value,
    BH_FDR, FDRTier, Subjects, Families, LeadingEdge, FigureLocation
  )

sample_method_table <- bind_rows(
  hscrp_pathways %>%
    distinct(Dataset, Tissue, State, StateLabel, Contrast, Subjects, Families) %>%
    mutate(
      RecordType = "Effective sample size",
      Item = paste(Dataset, State, sep = " / "),
      Value = paste0(Subjects, " participants; ", Families, " families")
    ) %>%
    select(RecordType, Item, Value),
  tribble(
    ~RecordType, ~Item, ~Value,
    "Method", "Predictor", "Baseline hsCRP; positive values rank-normalized",
    "Method", "Feature model", "Feature ~ hsCRP + Role; FamilyID duplicateCorrelation",
    "Method", "Response model", "Post3h-minus-Pre feature ~ hsCRP + Role; FamilyID duplicateCorrelation",
    "Method", "Pathway analysis", "MSigDB Hallmark GSEA ranked by limma moderated t statistic",
    "Method", "Multiplicity", "Benjamini-Hochberg correction within dataset, state, and predictor",
    "Method", "Primary threshold", "BH-FDR < 0.05",
    "Method", "Exploratory tiers", "BH-FDR < 0.10 and BH-FDR < 0.20",
    "Interpretation", "NES direction", "Direction within the hsCRP-associated ranked feature list; not pathway activation",
    "Interpretation", "Post3h state", "Absolute molecular profile; not an exercise-response contrast"
  )
)

write_csv(
  hscrp_pathways,
  file.path(analysis_dir, "Fig05C_Pathway_SourceData.csv")
)
write_csv(
  plot_data %>% mutate(Role = as.character(Role)),
  file.path(analysis_dir, "Fig05E_Molecular_SourceData.csv")
)
write.xlsx(
  list(
    Pathway_summary = profile_summary %>% mutate(across(where(is.factor), as.character)),
    Running_scores = profile_data %>% mutate(across(where(is.factor), as.character))
  ),
  file.path(analysis_dir, "Fig05D_GSEA_Profile_SourceData.xlsx"),
  overwrite = TRUE
)

apply_three_line_style <- function(workbook, sheet, data, title, note) {
  addWorksheet(workbook, sheet, gridLines = FALSE)
  writeData(workbook, sheet, title, startRow = 1, startCol = 1)
  mergeCells(workbook, sheet, cols = 1:max(1, ncol(data)), rows = 1)
  addStyle(
    workbook, sheet,
    createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold"),
    rows = 1, cols = 1, gridExpand = TRUE
  )
  writeData(workbook, sheet, data, startRow = 3, startCol = 1, withFilter = FALSE)
  header_style <- createStyle(
    fontName = "Arial", fontSize = 9, textDecoration = "bold",
    border = c("top", "bottom"), borderStyle = "thin",
    halign = "center", valign = "center", wrapText = TRUE
  )
  body_style <- createStyle(fontName = "Arial", fontSize = 8.5, valign = "center")
  bottom_style <- createStyle(border = "bottom", borderStyle = "thin")
  addStyle(workbook, sheet, header_style, rows = 3, cols = seq_len(ncol(data)), gridExpand = TRUE)
  if (nrow(data) > 0L) {
    addStyle(workbook, sheet, body_style, rows = 4:(3 + nrow(data)),
             cols = seq_len(ncol(data)), gridExpand = TRUE)
    addStyle(workbook, sheet, bottom_style, rows = 3 + nrow(data),
             cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
  }
  note_row <- 5 + nrow(data)
  writeData(workbook, sheet, note, startRow = note_row, startCol = 1)
  mergeCells(workbook, sheet, cols = 1:max(1, ncol(data)), rows = note_row)
  addStyle(
    workbook, sheet,
    createStyle(
      fontName = "Arial", fontSize = 8,
      textDecoration = "italic", wrapText = TRUE
    ),
    rows = note_row, cols = 1, gridExpand = TRUE
  )
  freezePane(workbook, sheet, firstActiveRow = 4)
  pageSetup(
    workbook, sheet, orientation = "landscape", paperSize = 9,
    fitToWidth = 1, fitToHeight = 0
  )
  setColWidths(workbook, sheet, cols = seq_len(ncol(data)), widths = "auto")
  setRowHeights(workbook, sheet, rows = note_row, heights = 34)
}

table_s10 <- file.path(
  table_dir,
  "Table S10 hsCRP molecular and Hallmark associations (Figure 5C-E and Figure S5A-C).xlsx"
)
wb <- createWorkbook(creator = "WYC")
apply_three_line_style(
  wb, "Hallmark pathways", pathway_table,
  "Table S10. Molecular and Hallmark pathway associations with baseline hsCRP.",
  paste0(
    "GSEA used all mapped proteins ranked by the moderated t statistic from role-adjusted, ",
    "family-blocked limma models. BH-FDR was calculated within each dataset-state analysis. ",
    "Positive and negative NES values indicate opposite ends of the hsCRP-associated ranking ",
    "and do not establish pathway activation or inhibition."
  )
)
apply_three_line_style(
  wb, "Feature candidates", feature_candidate_results,
  "Feature-level hsCRP association candidates across molecular states.",
  paste0(
    "PPIL3 met BH-FDR < 0.05 in the pre-exercise adipose-proteomic model. ",
    "RBPMS, ARFGEF1, and WASH3P/WASH2P met BH-FDR < 0.20 only in the response model. ",
    "WASH3P/WASH2P is displayed and reported explicitly as a non-unique measured protein group."
  )
)
apply_three_line_style(
  wb, "Methods and sample sizes", sample_method_table,
  "Analysis specifications and effective sample sizes.",
  paste0(
    "Pre-exercise hsCRP was the predictor in all molecular models. Post3h denotes the absolute ",
    "post-exercise molecular profile; response denotes Post3h minus Pre."
  )
)
saveWorkbook(wb, table_s10, overwrite = TRUE)

# Preserve the current Figure 5A table as Table S9 without altering its content.
table_s9_source <- path_project("fig05", "Suptable", "Figure5A_hsCRP_Clinical_Associations.xlsx")
table_s9_target <- file.path(
  table_dir,
  "Table S9 hsCRP clinical associations (Figure 5A).xlsx"
)
if (file.exists(table_s9_source) && !file.exists(table_s9_target)) {
  file.copy(table_s9_source, table_s9_target, overwrite = FALSE)
}

readme_lines <- c(
  "# Figure 5 current molecular update",
  "",
  paste0("Analysis ID: `", analysis_id, "`"),
  "",
  "## Protected panels",
  "",
  "- Figure 5A and Figure 5B were not modified.",
  "- Files under `fig05/editable/` and their same-basename exports were not modified.",
  "",
  "## Current panel plan",
  "",
  "- Figure 5C: core Hallmark pathway associations with baseline hsCRP across pre-exercise, Post3h, and response states.",
  "- Figure 5D: oxidative-phosphorylation and epithelial-mesenchymal-transition GSEA profiles.",
  "- Figure 5C-D: an exact 190-mm, 125:65 one-row composition is provided for direct assembly.",
  "- Figure 5E: PPIL3, RBPMS, ARFGEF1, and WASH3P/WASH2P feature-level associations.",
  "- Figure S5A: expanded Hallmark pathway matrix.",
  "- Figure S5B: PPIL3 and RBPMS associations with baseline hsCRP at pre-exercise, Post3h, and Post3h-minus-Pre states.",
  "- Figure S5C: fatty-acid-metabolism, adipogenesis, coagulation, and MYC-targets-v1 GSEA profiles.",
  "",
  "## Statistical interpretation",
  "",
  "- Primary feature threshold: BH-FDR < 0.05.",
  "- BH-FDR < 0.10 and BH-FDR < 0.20 are displayed as progressively weaker exploratory tiers.",
  "- GSEA NES describes enrichment within the hsCRP-associated molecular ranking; it is not direct evidence of pathway activation or inhibition.",
  "- Post3h is an absolute molecular profile. Response is Post3h minus Pre.",
  "",
  "## Main result",
  "",
  paste0(
    "Higher baseline hsCRP aligned with positively ranked adipose oxidative and lipid-metabolic ",
    "programs before exercise. In skeletal muscle, higher baseline hsCRP aligned with a negative ",
    "oxidative-phosphorylation response and positive remodeling- and coagulation-related responses. ",
    "Feature-level evidence was sparse: PPIL3 met BH-FDR < 0.05 before exercise, whereas RBPMS ",
    "and ARFGEF1 were exploratory response candidates at BH-FDR < 0.20."
  ),
  "",
  "## Tables",
  "",
  "- Table S9: hsCRP clinical associations supporting Figure 5A.",
  "- Table S10: all hsCRP feature and Hallmark pathway results supporting Figure 5C-E and Figure S5A-B."
)
writeLines(readme_lines, file.path(path_project("fig05"), "FIGURE05_RELEASE_NOTES.md"))

writeLines(c(
  paste0("Analysis ID: ", analysis_id),
  paste0("Formal input: ", formal_rds),
  paste0("Hallmark input: ", hallmark_rds),
  "Figure 5A and Figure 5B were not modified.",
  "No editable artwork or same-basename manual export was modified.",
  "",
  capture.output(sessionInfo())
), file.path(analysis_dir, "RUN_LOG.txt"))

message("Completed current Figure 5 molecular update in: ", path_project("fig05"))
