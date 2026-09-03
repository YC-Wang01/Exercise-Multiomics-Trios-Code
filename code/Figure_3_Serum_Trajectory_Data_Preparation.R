# Four-timepoint data preparation for Figure 3B and Figure S3A.
#
# Cluster assignments remain defined by the verified Pre/Post1h/Post3h
# exercise-response analyses. Four-timepoint trajectories are referenced to Fast.

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0L) {
  normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), "code", "Figure_3_Serum_Trajectory_Data_Preparation.R"),
                winslash = "/", mustWork = TRUE)
}
.local_config <- file.path(dirname(script_file), "Fig00_Config.R")
source(.local_config)
project <- PROJECT_ROOT

time_order <- c("Fast", "Pre", "Post1h", "Post3h")
time_labels <- c(Fast = "Fast", Pre = "Pre", Post1h = "1 h", Post3h = "3 h")
cluster_colors <- c("1" = "#0B8994", "2" = "#D86542", "3" = "#6C5A98", "4" = "#4F8A5B")
change_low <- "#4C78A8"
change_mid <- "#FFFFFF"
change_high <- "#D66A6A"

metadata_dir <- path_analysis_ready("Metadata")
sample_sheet <- read_csv(
  file.path(metadata_dir, "Project_Sample_Sheet.csv"), show_col_types = FALSE
) %>%
  mutate(
    SampleID = tolower(as.character(SampleID)),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    FamilyID = as.character(FamilyID),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort),
    RetainInPrimaryInput = as.logical(RetainInPrimaryInput),
    Timepoint = as.character(Timepoint)
  ) %>%
  filter(
    AnalysisSet == "Primary", Is_Primary_Cohort,
    is.na(RetainInPrimaryInput) | RetainInPrimaryInput,
    Timepoint %in% time_order
  )

met_source <- file.path(
  project, "fig03", "Figure3", "Fig03B_ThreePoint_Cluster_Definitions_Source.xlsx"
)
protein_source <- file.path(
  project, "fig03", "FigureS3", "FigS03ABC_Serum_Proteomics_SourceData.xlsx"
)
stopifnot(file.exists(met_source), file.exists(protein_source))

met_clusters <- read.xlsx(met_source, sheet = "Feature_Clusters") %>%
  filter(Dataset == "Serum_Metabonomics") %>%
  transmute(
    FeatureID = as.character(FeatureID),
    FeatureLabel = as.character(FeatureLabel),
    Cluster = as.integer(Cluster),
    SelectionMetric = as.numeric(MinimumFDR)
  )
met_modules <- read.xlsx(met_source, sheet = "Module_Definitions") %>%
  transmute(
    Cluster = as.integer(Cluster), FeatureID = as.character(FeatureID),
    Module = as.character(Module), DisplayName = as.character(DisplayName)
  )
met_display <- met_modules %>%
  left_join(select(met_clusters, FeatureID, FeatureLabel, SelectionMetric), by = "FeatureID") %>%
  mutate(FeatureLabel = coalesce(DisplayName, FeatureLabel, FeatureID)) %>%
  group_by(Cluster) %>%
  arrange(SelectionMetric, FeatureID, .by_group = TRUE) %>%
  slice_head(n = 4L) %>%
  ungroup()

protein_clusters <- read.xlsx(protein_source, sheet = "Selected_Proteins") %>%
  transmute(
    FeatureID = as.character(FeatureID),
    FeatureLabel = as.character(FeatureLabel),
    Cluster = as.integer(Cluster),
    Module = as.character(Module),
    SelectionMetric = as.numeric(MinimumP)
  )
protein_display <- protein_clusters %>%
  group_by(Cluster) %>%
  arrange(SelectionMetric, FeatureLabel, .by_group = TRUE) %>%
  slice_head(n = 4L) %>%
  ungroup()

dataset_specs <- list(
  Serum_Metabonomics = list(
    input = path_analysis_ready("Analysis_Serum_Metabonomics.csv"),
    clusters = met_clusters,
    display = met_display,
    panel = "Serum metabolite response clusters",
    prefix = "Fig03B_Serum_Metabolite_Clusters"
  ),
  Serum_Proteomics = list(
    input = path_analysis_ready("Analysis_Serum_Proteomics.csv"),
    clusters = protein_clusters,
    display = protein_display,
    panel = "Serum protein response clusters",
    prefix = "FigS03A_Serum_Protein_Clusters"
  )
)

mean_ci <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(tibble(Mean = NA_real_, Lower95 = NA_real_, Upper95 = NA_real_, N = 0L))
  se <- if (length(x) > 1L) sd(x) / sqrt(length(x)) else NA_real_
  tibble(
    Mean = mean(x), Lower95 = mean(x) - 1.96 * se,
    Upper95 = mean(x) + 1.96 * se, N = length(x)
  )
}

prepare_four_timepoint <- function(dataset_name, spec) {
  matrix <- read_analysis_matrix(spec$input)
  meta <- sample_sheet %>%
    filter(Dataset == dataset_name, SampleID %in% colnames(matrix)) %>%
    distinct(Clinical_Subject_ID, Timepoint, .keep_all = TRUE)

  complete_subjects <- meta %>%
    distinct(Clinical_Subject_ID, Timepoint) %>%
    count(Clinical_Subject_ID, name = "Timepoints") %>%
    filter(Timepoints == length(time_order)) %>%
    pull(Clinical_Subject_ID)
  meta <- meta %>%
    filter(Clinical_Subject_ID %in% complete_subjects) %>%
    mutate(Timepoint = factor(Timepoint, levels = time_order)) %>%
    arrange(Clinical_Subject_ID, Timepoint)
  if (length(complete_subjects) < 10L) stop("Too few complete four-timepoint subjects for ", dataset_name)

  cluster_key <- spec$clusters %>%
    filter(FeatureID %in% rownames(matrix), Cluster %in% 1:4) %>%
    distinct(FeatureID, .keep_all = TRUE)
  if (nrow(cluster_key) == 0L) stop("No cluster features matched ", dataset_name)

  mat <- matrix[cluster_key$FeatureID, meta$SampleID, drop = FALSE]
  row_mean <- rowMeans(mat, na.rm = TRUE)
  row_sd <- apply(mat, 1, sd, na.rm = TRUE)
  row_sd[!is.finite(row_sd) | row_sd == 0] <- 1
  zmat <- sweep(sweep(mat, 1, row_mean, "-"), 1, row_sd, "/")

  absolute <- as.data.frame(zmat) %>%
    rownames_to_column("FeatureID") %>%
    pivot_longer(-FeatureID, names_to = "SampleID", values_to = "AbsoluteZ") %>%
    left_join(select(meta, SampleID, Clinical_Subject_ID, Timepoint), by = "SampleID") %>%
    left_join(select(cluster_key, FeatureID, FeatureLabel, Cluster), by = "FeatureID")

  fast_values <- absolute %>%
    filter(Timepoint == "Fast") %>%
    select(Clinical_Subject_ID, FeatureID, FastZ = AbsoluteZ)
  values <- absolute %>%
    left_join(fast_values, by = c("Clinical_Subject_ID", "FeatureID")) %>%
    mutate(RelativeToFast = AbsoluteZ - FastZ)

  summaries <- bind_rows(
    values %>%
      group_by(FeatureID, FeatureLabel, Cluster, Timepoint) %>%
      group_modify(~mean_ci(.x$AbsoluteZ)) %>%
      ungroup() %>% mutate(DisplayType = "Absolute", Value = Mean),
    values %>%
      group_by(FeatureID, FeatureLabel, Cluster, Timepoint) %>%
      group_modify(~mean_ci(.x$RelativeToFast)) %>%
      ungroup() %>% mutate(DisplayType = "Relative", Value = Mean)
  )

  list(
    values = values,
    summaries = summaries,
    clusters = cluster_key,
    display = spec$display %>% filter(FeatureID %in% cluster_key$FeatureID),
    complete_subjects = complete_subjects
  )
}

make_segments <- function(data, value_col = "Value") {
  data %>%
    mutate(TimeIndex = as.integer(Timepoint)) %>%
    arrange(FeatureID, TimeIndex) %>%
    group_by(FeatureID) %>%
    mutate(
      Xend = lead(TimeIndex), Yend = lead(.data[[value_col]]),
      Segment = case_when(
        TimeIndex == 1L ~ "Pre-intervention",
        TimeIndex %in% c(2L, 3L) ~ "Exercise and recovery",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Xend)) %>%
    ungroup()
}

make_dataset_plot <- function(dataset_name, spec, prepared, display_type) {
  summary_data <- prepared$summaries %>%
    filter(DisplayType == display_type) %>%
    mutate(Timepoint = factor(Timepoint, levels = time_order))
  display_key <- prepared$display %>%
    distinct(Cluster, FeatureID, FeatureLabel, .keep_all = TRUE)
  display_summary <- summary_data %>%
    semi_join(display_key, by = c("Cluster", "FeatureID")) %>%
    select(-FeatureLabel) %>%
    left_join(
      display_key %>% select(Cluster, FeatureID, FeatureLabel),
      by = c("Cluster", "FeatureID")
    )

  heat_limit <- quantile(abs(display_summary$Value), 0.95, na.rm = TRUE)
  heat_limit <- ceiling(heat_limit * 10) / 10
  if (!is.finite(heat_limit) || heat_limit <= 0) heat_limit <- 0.5

  make_cluster <- function(cluster_id) {
    cluster_all <- summary_data %>% filter(Cluster == cluster_id)
    cluster_display <- display_summary %>% filter(Cluster == cluster_id)
    key <- display_key %>%
      filter(Cluster == cluster_id) %>%
      arrange(SelectionMetric, FeatureLabel)
    feature_levels <- key$FeatureLabel
    feature_palette <- setNames(
      grDevices::hcl.colors(max(1L, length(feature_levels)), palette = "Dark 3"),
      feature_levels
    )

    centroid <- cluster_all %>%
      group_by(Timepoint) %>%
      summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
      mutate(FeatureID = "Cluster centroid", TimeIndex = as.integer(Timepoint))
    all_segments <- make_segments(cluster_all)
    centroid_segments <- make_segments(centroid)

    profile <- ggplot() +
      geom_hline(yintercept = 0, linewidth = 0.20, colour = "black") +
      geom_segment(
        data = all_segments,
        aes(x = TimeIndex, y = Value, xend = Xend, yend = Yend, linetype = Segment),
        colour = cluster_colors[as.character(cluster_id)], linewidth = 0.20, alpha = 0.18
      ) +
      geom_segment(
        data = centroid_segments,
        aes(x = TimeIndex, y = Value, xend = Xend, yend = Yend, linetype = Segment),
        colour = cluster_colors[as.character(cluster_id)], linewidth = 0.80
      ) +
      geom_point(
        data = centroid, aes(x = TimeIndex, y = Value),
        colour = cluster_colors[as.character(cluster_id)], size = 0.95
      ) +
      scale_linetype_manual(values = c("Pre-intervention" = "22", "Exercise and recovery" = "solid"), guide = "none") +
      scale_x_continuous(breaks = 1:4, labels = unname(time_labels), limits = c(0.85, 4.15)) +
      labs(
        x = NULL,
        y = if (cluster_id %in% c(1L, 3L)) {
          if (display_type == "Absolute") "Standardized\nabundance" else "Difference\nfrom Fast (SD)"
        } else NULL
      ) +
      theme_classic(base_family = "Arial", base_size = 5.0) +
      theme(
        axis.title = element_text(size = 5.0, colour = "black"),
        axis.text = element_text(size = 4.5, colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.30),
        axis.line = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(0.7, 0.4, 0.7, 0.7)
      )

    line_data <- cluster_display %>%
      mutate(FeatureLabel = factor(FeatureLabel, levels = feature_levels))
    line_segments <- make_segments(line_data) %>%
      mutate(FeatureLabel = factor(FeatureLabel, levels = feature_levels))
    feature_lines <- ggplot() +
      geom_hline(yintercept = 0, linewidth = 0.20, colour = "black") +
      geom_segment(
        data = line_segments,
        aes(x = TimeIndex, y = Value, xend = Xend, yend = Yend,
            colour = FeatureLabel, linetype = Segment),
        linewidth = 0.48
      ) +
      geom_point(
        data = line_data,
        aes(x = as.integer(Timepoint), y = Value, colour = FeatureLabel), size = 0.82
      ) +
      scale_colour_manual(values = feature_palette, guide = "none") +
      scale_linetype_manual(values = c("Pre-intervention" = "22", "Exercise and recovery" = "solid"), guide = "none") +
      scale_x_continuous(breaks = 1:4, labels = unname(time_labels), limits = c(0.85, 4.15)) +
      labs(x = NULL, y = NULL) +
      theme_classic(base_family = "Arial", base_size = 5.0) +
      theme(
        axis.text = element_text(size = 4.5, colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.30),
        axis.line = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(0.7, -1.0, 0.7, 0.3)
      )

    tile_data <- cluster_display %>%
      mutate(
        FeatureLabel = factor(FeatureLabel, levels = rev(feature_levels)),
        TimeIndex = as.integer(Timepoint)
      )
    label_data <- tile_data %>%
      distinct(FeatureLabel) %>%
      mutate(
        Label = str_wrap(as.character(FeatureLabel), width = 18),
        FeatureColour = feature_palette[as.character(FeatureLabel)]
      )
    heatmap <- ggplot(tile_data, aes(x = TimeIndex, y = FeatureLabel, fill = Value)) +
      geom_tile(colour = "black", linewidth = 0.12, width = 0.98, height = 0.98) +
      geom_point(
        data = label_data,
        aes(x = -0.20, y = FeatureLabel, colour = FeatureColour),
        inherit.aes = FALSE, size = 2.10
      ) +
      geom_text(
        data = label_data,
        aes(x = 4.55, y = FeatureLabel, label = Label),
        inherit.aes = FALSE, hjust = 0, size = 1.65,
        family = "Arial", colour = "black", lineheight = 0.84
      ) +
      scale_colour_identity() +
      scale_fill_gradient2(
        low = change_low, mid = change_mid, high = change_high,
        midpoint = 0, limits = c(-heat_limit, heat_limit),
        oob = scales::squish, guide = "none"
      ) +
      scale_x_continuous(breaks = 1:4, labels = unname(time_labels), limits = c(-0.7, 7.20)) +
      scale_y_discrete(drop = FALSE) +
      coord_cartesian(clip = "off") +
      labs(x = NULL, y = NULL) +
      theme_classic(base_family = "Arial", base_size = 4.8) +
      theme(
        axis.text.x = element_text(size = 4.0, angle = 0, hjust = 0.5, colour = "black"),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        panel.border = element_blank(), axis.line = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(0.7, 0.4, 0.7, -1.4)
      )

    aligned <- cowplot::align_plots(profile, feature_lines, heatmap, align = "v", axis = "tb")
    body <- cowplot::plot_grid(
      plotlist = aligned, nrow = 1,
      rel_widths = c(0.24, 0.28, 0.48), align = "h", axis = "tb"
    )
    module_label <- key$Module[which(!is.na(key$Module) & nzchar(key$Module))[1]]
    if (length(module_label) == 0L || is.na(module_label)) module_label <- "Temporal response module"
    title <- ggdraw() + draw_label(
      paste0("Cluster ", cluster_id, "\n", module_label),
      x = 0.5, y = 0.48, fontfamily = "Arial", fontface = "bold",
      size = 5.4, colour = "black", lineheight = 0.90
    )
    wrap_elements(full = title) / wrap_elements(full = body) +
      plot_layout(heights = c(0.17, 0.83))
  }

  panels <- lapply(1:4, make_cluster)
  main <- wrap_plots(panels, ncol = 2)

  legend_title <- if (display_type == "Absolute") "Standardized abundance" else "Difference from Fast (SD)"
  legend_source <- ggplot(
    tibble(x = 1, y = 1, z = 0), aes(x = x, y = y, fill = z)
  ) +
    geom_tile() +
    scale_fill_gradient2(
      low = change_low, mid = change_mid, high = change_high,
      midpoint = 0, limits = c(-heat_limit, heat_limit),
      name = legend_title,
      guide = guide_colourbar(
        direction = "vertical", title.position = "top", title.hjust = 0.5,
        barwidth = grid::unit(2, "mm"), barheight = grid::unit(15, "mm"),
        ticks = FALSE, frame.colour = NA
      )
    ) +
    theme_void(base_family = "Arial") +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 4.7, face = "bold", colour = "black", hjust = 0.5),
      legend.text = element_text(size = 4.3, colour = "black"),
      legend.title.position = "top",
      legend.key.height = grid::unit(15, "mm"),
      legend.key.width = grid::unit(2, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0, 0, 0, 0)
    )
  colour_legend <- cowplot::get_legend(legend_source)
  line_legend <- ggdraw() +
    draw_line(x = c(0.05, 0.28), y = c(0.68, 0.68), colour = "black", linewidth = 0.55, linetype = "22") +
    draw_label("Fast-Pre", x = 0.34, y = 0.68, hjust = 0, fontfamily = "Arial", size = 4.3) +
    draw_line(x = c(0.05, 0.28), y = c(0.28, 0.28), colour = "black", linewidth = 0.55) +
    draw_label("Pre-3 h", x = 0.34, y = 0.28, hjust = 0, fontfamily = "Arial", size = 4.3)
  right_legend <- cowplot::plot_grid(
    NULL, colour_legend, line_legend, NULL,
    ncol = 1, rel_heights = c(0.20, 0.34, 0.16, 0.30)
  )

  combined <- wrap_elements(full = main) | wrap_elements(full = right_legend)
  combined <- combined +
    plot_layout(widths = c(0.91, 0.09)) +
    plot_annotation(
      title = paste0(spec$panel, " (four-timepoint display)"),
      subtitle = if (display_type == "Absolute") {
        "Original three-point cluster assignments; standardized abundance shown at fasting, Pre, 1 h, and 3 h"
      } else {
        "Original three-point cluster assignments; within-participant differences are referenced to Fast"
      },
      theme = theme(
        text = element_text(family = "Arial", colour = "black"),
        plot.title = element_text(size = 8.2, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 5.4, hjust = 0, colour = "black")
      )
    )
  combined
}

if (identical(Sys.getenv("FIG03_FOURPOINT_LEGACY_DISPLAY"), "1")) {
out_dir <- file.path(project, "fig03", "analysis", "Legacy_FourTimepoint_Display")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prepared_all <- list()
source_sheets <- list()
summary_rows <- list()

for (dataset_name in names(dataset_specs)) {
  spec <- dataset_specs[[dataset_name]]
  prepared <- prepare_four_timepoint(dataset_name, spec)
  prepared_all[[dataset_name]] <- prepared
  summary_rows[[dataset_name]] <- tibble(
    Dataset = dataset_name,
    CompleteFourTimepointSubjects = length(prepared$complete_subjects),
    ClusterFeatures = nrow(prepared$clusters),
    DisplayFeatures = nrow(prepared$display)
  )

  for (display_type in c("Absolute", "Relative")) {
    plot <- make_dataset_plot(dataset_name, spec, prepared, display_type)
    filename <- paste0(spec$prefix, "_FourTimepoint_", display_type, "_Legacy.pdf")
    ggsave(
      file.path(out_dir, filename), plot,
      width = 190, height = 80, units = "mm", device = cairo_pdf
    )
  }

  source_sheets[[paste0(dataset_name, "_Summary")]] <- prepared$summaries %>%
    mutate(Timepoint = as.character(Timepoint))
  source_sheets[[paste0(dataset_name, "_Values")]] <- prepared$values %>%
    mutate(Timepoint = as.character(Timepoint))
}

analysis_summary <- bind_rows(summary_rows)
source_sheets <- c(list(Analysis_Summary = analysis_summary), source_sheets)
write.xlsx(
  source_sheets,
  file.path(out_dir, "Fig03B_S03A_FourTimepoint_Display_SourceData.xlsx"),
  overwrite = TRUE
)

writeLines(
  c(
    "# Legacy Figure 3B and Figure S3A four-timepoint displays",
    "",
    "- Cluster assignments were not recalculated.",
    "- Fast is the reference state for displayed four-timepoint differences.",
    "- Pre remains the formal exercise baseline for the original three-point cluster assignments.",
    "- The dashed Fast-to-Pre segment is a pre-intervention interval.",
    "- Solid Pre-to-Post segments describe exercise and recovery.",
    "- Absolute panels show feature-wise standardized abundance.",
    "- Relative panels show within-participant differences from Fast in feature SD units.",
    "- Only participants with all four serum timepoints were retained for these displays.",
    "",
    paste(capture.output(print(analysis_summary)), collapse = "\n")
  ),
  file.path(out_dir, "README.md"), useBytes = TRUE
)
capture.output(sessionInfo(), file = file.path(out_dir, "R_SESSION_INFO.txt"))
warning_log <- warnings()
if (is.null(warning_log)) {
  writeLines("No warnings.", file.path(out_dir, "WARNINGS.txt"))
} else {
  capture.output(warning_log, file = file.path(out_dir, "WARNINGS.txt"))
}
message("Legacy four-timepoint displays written to: ", out_dir)
}
