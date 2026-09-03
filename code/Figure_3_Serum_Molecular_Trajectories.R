# Overall and representative-feature four-timepoint displays for Figure 3B and Figure S3A.
# Cluster assignments remain defined by the verified Pre/Post1h/Post3h analysis;
# displayed four-timepoint differences are referenced to Fast.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0L) {
  normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), "code", "Figure_3_Serum_Molecular_Trajectories.R"),
                winslash = "/", mustWork = TRUE)
}
Sys.setenv(FIG03_FOURPOINT_SOURCE_ONLY = "1")
source(file.path(dirname(script_file), "Figure_3_Serum_Trajectory_Data_Preparation.R"))
Sys.unsetenv("FIG03_FOURPOINT_SOURCE_ONLY")

output_map <- list(
  Serum_Metabonomics = list(
    pdf = file.path(project, "fig03", "Figure3", "Fig03B_Serum_Metabolite_Response_Clusters.pdf"),
    source = file.path(project, "fig03", "Figure3", "Fig03B_Serum_Metabolite_Response_Clusters_SourceData.xlsx"),
    runlog = file.path(project, "fig03", "Figure3", "Fig03B_Serum_Metabolite_Response_Clusters_RunLog.txt")
  ),
  Serum_Proteomics = list(
    pdf = file.path(project, "fig03", "FigureS3", "FigS03A_Serum_Protein_Response_Clusters_Integrated.pdf"),
    source = file.path(project, "fig03", "FigureS3", "FigS03A_Serum_Protein_Response_Clusters_SourceData.xlsx"),
    runlog = file.path(project, "fig03", "FigureS3", "FigS03A_Serum_Protein_Response_Clusters_RunLog.txt")
  )
)

make_overall_representative_plot <- function(dataset_name, spec, prepared) {
  summary_data <- prepared$summaries %>%
    filter(DisplayType == "Relative") %>%
    mutate(Timepoint = factor(Timepoint, levels = time_order))

  representative_key <- prepared$display %>%
    filter(FeatureID %in% prepared$clusters$FeatureID) %>%
    group_by(Cluster) %>%
    arrange(SelectionMetric, FeatureLabel, .by_group = TRUE) %>%
    slice_head(n = 3L) %>%
    ungroup() %>%
    distinct(Cluster, FeatureID, .keep_all = TRUE)

  representative_summary <- summary_data %>%
    semi_join(representative_key, by = c("Cluster", "FeatureID")) %>%
    select(-FeatureLabel) %>%
    left_join(
      representative_key %>% select(Cluster, FeatureID, FeatureLabel),
      by = c("Cluster", "FeatureID")
    )

  make_cluster_panel <- function(cluster_id) {
    cluster_all <- summary_data %>% filter(Cluster == cluster_id)
    cluster_representatives <- representative_summary %>% filter(Cluster == cluster_id)
    key <- representative_key %>%
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

    overall <- ggplot() +
      geom_hline(yintercept = 0, linewidth = 0.20, colour = "black") +
      geom_segment(
        data = all_segments,
        aes(x = TimeIndex, y = Value, xend = Xend, yend = Yend, linetype = Segment),
        colour = cluster_colors[as.character(cluster_id)], linewidth = 0.22, alpha = 0.18
      ) +
      geom_segment(
        data = centroid_segments,
        aes(x = TimeIndex, y = Value, xend = Xend, yend = Yend, linetype = Segment),
        colour = cluster_colors[as.character(cluster_id)], linewidth = 0.90
      ) +
      geom_point(
        data = centroid, aes(x = TimeIndex, y = Value),
        colour = cluster_colors[as.character(cluster_id)], size = 1.05
      ) +
      scale_linetype_manual(
        values = c("Pre-intervention" = "22", "Exercise and recovery" = "solid"),
        guide = "none"
      ) +
      scale_x_continuous(
        breaks = 1:4, labels = unname(time_labels), limits = c(0.85, 4.15)
      ) +
      labs(
        x = NULL,
        y = if (cluster_id == 1L) "Difference from Fast (SD)" else NULL,
        title = "Cluster-level response"
      ) +
      theme_classic(base_family = "Arial", base_size = 5.0) +
      theme(
        plot.title = element_text(size = 5.0, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 5.0, colour = "black"),
        axis.text = element_text(size = 4.5, colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.30),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.line = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(0.8, 0.6, 0.8, 0.8)
      )

    representative_data <- cluster_representatives %>%
      mutate(
        FeatureLabel = factor(FeatureLabel, levels = feature_levels),
        TimeIndex = as.integer(Timepoint)
      )
    representative_segments <- make_segments(representative_data) %>%
      mutate(FeatureLabel = factor(FeatureLabel, levels = feature_levels))
    endpoint_labels <- representative_data %>%
      filter(Timepoint == "Post3h") %>%
      mutate(Label = str_wrap(as.character(FeatureLabel), width = 18))

    representatives <- ggplot() +
      geom_hline(yintercept = 0, linewidth = 0.20, colour = "black") +
      geom_segment(
        data = representative_segments,
        aes(
          x = TimeIndex, y = Value, xend = Xend, yend = Yend,
          colour = FeatureLabel, linetype = Segment
        ),
        linewidth = 0.62
      ) +
      geom_point(
        data = representative_data,
        aes(x = TimeIndex, y = Value, colour = FeatureLabel), size = 1.00
      ) +
      ggrepel::geom_text_repel(
        data = endpoint_labels,
        aes(x = TimeIndex, y = Value, label = Label, colour = FeatureLabel),
        direction = "y", hjust = 0, nudge_x = 0.16,
        min.segment.length = 0, segment.colour = "#777777", segment.size = 0.18,
        box.padding = 0.10, point.padding = 0.05, max.overlaps = Inf,
        size = 1.65, family = "Arial", lineheight = 0.85, show.legend = FALSE
      ) +
      scale_colour_manual(values = feature_palette, guide = "none") +
      scale_linetype_manual(
        values = c("Pre-intervention" = "22", "Exercise and recovery" = "solid"),
        guide = "none"
      ) +
      scale_x_continuous(
        breaks = 1:4, labels = unname(time_labels), limits = c(0.85, 4.95)
      ) +
      coord_cartesian(clip = "off") +
      labs(x = NULL, y = NULL, title = "Representative molecular trajectories") +
      theme_classic(base_family = "Arial", base_size = 5.0) +
      theme(
        plot.title = element_text(size = 4.5, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 4.0, colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.border = element_rect(
          colour = "black", fill = NA, linewidth = 0.30, linetype = "22"
        ),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.line = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(0.8, 4.5, 0.8, 0.5)
      )

    module_label <- key$Module[which(!is.na(key$Module) & nzchar(key$Module))[1]]
    if (length(module_label) == 0L || is.na(module_label)) {
      module_label <- "Temporal response module"
    }
    cluster_title <- ggdraw() + draw_label(
      paste0("Cluster ", cluster_id, " (", module_label, ")"),
      x = 0.5, y = 0.48, fontfamily = "Arial", fontface = "bold",
      size = 5.2, colour = "black"
    )
    body_base <- cowplot::plot_grid(
      overall, representatives, nrow = 1,
      rel_widths = c(0.36, 0.64), align = "h", axis = "tb"
    )
    body <- ggdraw() +
      draw_line(
        x = c(0.345, 0.382), y = c(0.70, 0.80),
        colour = "black", linewidth = 0.32, linetype = "22"
      ) +
      draw_line(
        x = c(0.345, 0.382), y = c(0.30, 0.20),
        colour = "black", linewidth = 0.32, linetype = "22"
      ) +
      draw_plot(body_base, x = 0, y = 0, width = 1, height = 1)
    wrap_elements(full = cluster_title) / wrap_elements(full = body) +
      plot_layout(heights = c(0.13, 0.87))
  }

  panels <- lapply(1:4, make_cluster_panel)
  main <- wrap_plots(panels, ncol = 2)
  main +
    plot_annotation(
      title = paste0(spec$panel, " and representative molecular trajectories"),
      theme = theme(
        text = element_text(family = "Arial", colour = "black"),
        plot.title = element_text(size = 8.2, face = "bold", hjust = 0)
      )
    )
}

summary_rows <- list()

for (dataset_name in names(dataset_specs)) {
  spec <- dataset_specs[[dataset_name]]
  output <- output_map[[dataset_name]]
  prepared <- prepare_four_timepoint(dataset_name, spec)
  representative_key <- prepared$display %>%
    group_by(Cluster) %>%
    arrange(SelectionMetric, FeatureLabel, .by_group = TRUE) %>%
    slice_head(n = 3L) %>%
    ungroup()
  plot <- make_overall_representative_plot(dataset_name, spec, prepared)
  ggsave(
    output$pdf, plot,
    width = 190, height = 78, units = "mm", device = cairo_pdf
  )

  dataset_summary <- tibble(
    Dataset = dataset_name,
    CompleteFourTimepointSubjects = length(prepared$complete_subjects),
    ClusterFeatures = nrow(prepared$clusters),
    RepresentativeFeatures = nrow(representative_key)
  )
  summary_rows[[dataset_name]] <- dataset_summary
  relative_summary <- prepared$summaries %>%
    filter(DisplayType == "Relative") %>%
    mutate(Timepoint = as.character(Timepoint))
  writexl::write_xlsx(
    x = list(
      Analysis_Summary = dataset_summary,
      Representative_Features = representative_key,
      Relative_Trajectories = relative_summary
    ),
    path = output$source
  )
  writeLines(
    c(
      "Figure 3 four-timepoint response-cluster analysis",
      "",
      paste0("Dataset: ", dataset_name),
      paste0("Complete four-timepoint participants: ", length(prepared$complete_subjects)),
      paste0("Cluster features: ", nrow(prepared$clusters)),
      paste0("Representative features: ", nrow(representative_key)),
      "Clusters were defined from Pre, Post1h, and Post3h.",
      "Fast is the display reference and did not contribute to clustering.",
      "Displayed values are within-participant differences from Fast in feature SD units.",
      "No inferential significance is encoded in these descriptive panels."
    ),
    output$runlog, useBytes = TRUE
  )
  capture.output(
    sessionInfo(),
    file = sub("_RunLog\\.txt$", "_SessionInfo.txt", output$runlog)
  )
}

analysis_summary <- bind_rows(summary_rows)
message("Final Figure 3B and Figure S3A written to their formal locations.")
