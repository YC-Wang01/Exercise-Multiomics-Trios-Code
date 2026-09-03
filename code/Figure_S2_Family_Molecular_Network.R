# Candidate pair-specific familial molecular similarity networks for Figure S2.
# The primary panel excludes DNA methylation because baseline between-family
# methylation inference remains confounded by the unresolved design-version block.

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
if (is.na(.local_config)) stop("Cannot find CellMetabolism_Transfer code/Fig00_Config.R.")
source(.local_config)
rm(.local_script, .local_config_candidates, .local_config)

p_load(
  digest, dplyr, ggforce, ggplot2, ggraph, ggrepel, graphlayouts, igraph,
  patchwork, readr, scales, stringr, tibble, tidyr, writexl
)

analysis_id <- "Fig02S_Family_Molecular_Network_Stress_2026-08-25"
family_fdr <- 0.05
minimum_fraction_of_platform_max_pairs <- 0.80
layout_seed <- 121L
primary_label_cap <- 12L
hybrid_label_cap <- 16L
exploratory_label_cap <- 14L

input_file <- path_project(
  "fig02", "Family_Similarity_Audit",
  "All_FDR_Significant_Family_Features.csv"
)
output_dir <- path_project(
  "fig02", "Family_Molecular_Network_2026-08-24"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Missing family-similarity audit input: ", input_file)

pair_levels <- c("DM", "DF", "MF")
pair_labels <- c(
  DM = "Daughter-mother",
  DF = "Daughter-father",
  MF = "Mother-father"
)
pair_colors <- c(DM = "#C45A78", DF = "#3E7CB1", MF = "#5B8E6F")
tissue_colors <- c(Serum = "#D56A6A", Adipose = "#D9A441", Muscle = "#527EAF")
tissue_dark <- c(Serum = "#A94747", Adipose = "#A67821", Muscle = "#315F92")
omics_shapes <- c(
  Metabolomics = 21,
  Proteomics = 22,
  Transcriptomics = 24,
  `DNA methylation` = 23
)

audit <- readr::read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    Gene = str_split_i(coalesce(GeneSymbol, ""), fixed(";"), 1),
    Gene = na_if(Gene, ""),
    DisplayLabel = case_when(
      !is.na(Gene) ~ Gene,
      TRUE ~ as.character(FeatureID)
    )
  )

make_pair_edges <- function(data, pair) {
  q_col <- paste0(pair, "_SimilarityBH_FDR")
  p_col <- paste0(pair, "_SimilarityP")
  r_col <- paste0(pair, "_Correlation")
  n_col <- paste0(pair, "_NPairs")
  stable_col <- paste0(pair, "_StabilityFlag")

  data %>%
    transmute(
      Dataset, Tissue, Omics, FeatureID, Gene, DisplayLabel, BroadFunction,
      Pair = pair,
      NPairs = suppressWarnings(as.numeric(.data[[n_col]])),
      Correlation = suppressWarnings(as.numeric(.data[[r_col]])),
      P_Value = suppressWarnings(as.numeric(.data[[p_col]])),
      BH_FDR = suppressWarnings(as.numeric(.data[[q_col]])),
      StabilityFlag = as.character(.data[[stable_col]])
    )
}

all_edges <- bind_rows(lapply(pair_levels, function(pair) make_pair_edges(audit, pair))) %>%
  mutate(
    Pair = factor(Pair, levels = pair_levels),
    Tissue = factor(Tissue, levels = c("Serum", "Adipose", "Muscle")),
    Omics = factor(
      Omics,
      levels = c("Metabolomics", "Proteomics", "Transcriptomics", "DNA methylation")
    ),
    NodeID = paste(Dataset, FeatureID, sep = "__")
  )

platform_pair_max_n <- all_edges %>%
  filter(is.finite(NPairs)) %>%
  group_by(Dataset, Pair) %>%
  summarise(MaxPairs = max(NPairs, na.rm = TRUE), .groups = "drop")

eligible_edges <- all_edges %>%
  left_join(platform_pair_max_n, by = c("Dataset", "Pair")) %>%
  filter(
    is.finite(BH_FDR), BH_FDR < family_fdr,
    is.finite(NPairs)
  ) %>%
  mutate(
    Pair = factor(Pair, levels = pair_levels),
    CoverageOK = NPairs >= ceiling(minimum_fraction_of_platform_max_pairs * MaxPairs),
    LOFOStable = StabilityFlag == "Stable by prespecified LOFO diagnostics" & CoverageOK,
    StabilityClass = if_else(LOFOStable, "LOFO-stable", "FDR-only"),
    EvidenceClass = if_else(
      Omics == "DNA methylation",
      "Exploratory DNAm",
      "Primary non-DNAm"
    )
  )

# The primary graph includes every non-methylation feature with BH-FDR < 0.05
# in at least one family pair. The exploratory graph includes all molecular
# layers at the same threshold. Multi-pair nodes are visually prioritized.
# PIEZO1-linked DNAm is intentionally retained in the exploratory graph but
# remains subject to the documented design-version-block limitation.
primary_edges <- eligible_edges %>%
  filter(Omics != "DNA methylation")

hybrid_edges <- eligible_edges %>%
  group_by(NodeID) %>%
  filter(
    as.character(first(Omics)) != "DNA methylation" |
      n_distinct(Pair) >= 2
  ) %>%
  ungroup() %>%
  distinct(Pair, NodeID, .keep_all = TRUE)

exploratory_edges <- eligible_edges %>%
  distinct(Pair, NodeID, .keep_all = TRUE)

build_plot_data <- function(edges) {
  nodes <- edges %>%
    group_by(NodeID, Dataset, Tissue, Omics, FeatureID, Gene, DisplayLabel, EvidenceClass) %>%
    summarise(
      MinBH_FDR = min(BH_FDR, na.rm = TRUE),
      MaxAbsCorrelation = max(abs(Correlation), na.rm = TRUE),
      PairCount = n_distinct(Pair),
      .groups = "drop"
    ) %>%
    mutate(
      Tissue = factor(Tissue, levels = c("Serum", "Adipose", "Muscle")),
      Omics = factor(
        Omics,
        levels = c("Metabolomics", "Proteomics", "Transcriptomics", "DNA methylation")
      ),
      Label = if_else(
        Omics == "DNA methylation",
        paste0(DisplayLabel, "-linked DNAm"),
        DisplayLabel
      ),
      NodeSignificance = pmin(pmax(-log10(MinBH_FDR), 1.3), 4.0)
    ) %>%
    arrange(Tissue, Omics, Label)

  y_cursor <- 0
  node_positions <- list()
  tissue_bands <- list()
  for (tissue_name in c("Serum", "Adipose", "Muscle")) {
    tissue_nodes <- nodes %>% filter(as.character(Tissue) == tissue_name)
    if (nrow(tissue_nodes) == 0) next
    y_values <- y_cursor - seq(0, nrow(tissue_nodes) - 1)
    tissue_nodes$y <- y_values
    node_positions[[tissue_name]] <- tissue_nodes
    tissue_bands[[tissue_name]] <- tibble(
      Tissue = tissue_name,
      ymin = min(y_values) - 0.45,
      ymax = max(y_values) + 0.45,
      ymid = mean(range(y_values))
    )
    y_cursor <- min(y_values) - 1.8
  }
  nodes <- bind_rows(node_positions)
  bands <- bind_rows(tissue_bands)

  y_shift <- mean(range(nodes$y))
  nodes <- nodes %>% mutate(y = y - y_shift)
  bands <- bands %>% mutate(
    ymin = ymin - y_shift,
    ymax = ymax - y_shift,
    ymid = ymid - y_shift
  )

  all_y <- range(nodes$y)
  pair_y <- c(
    DM = all_y[2] - 0.20 * diff(all_y),
    DF = mean(all_y),
    MF = all_y[1] + 0.20 * diff(all_y)
  )
  hubs <- tibble(
    Pair = factor(pair_levels, levels = pair_levels),
    PairLabel = unname(pair_labels[pair_levels]),
    x = 0.24,
    y = unname(pair_y[pair_levels])
  )

  plot_edges <- edges %>%
    select(Pair, NodeID, Correlation, P_Value, BH_FDR, NPairs, StabilityFlag, EvidenceClass) %>%
    left_join(nodes %>% select(NodeID, node_y = y), by = "NodeID") %>%
    left_join(hubs %>% select(Pair, hub_y = y), by = "Pair") %>%
    mutate(
      x = 0.35,
      xend = 1.08,
      EdgeAlpha = if_else(EvidenceClass == "Exploratory DNAm", 0.62, 0.78)
    )

  list(nodes = nodes, bands = bands, hubs = hubs, edges = plot_edges)
}

plot_family_network <- function(edges, title, subtitle) {
  pd <- build_plot_data(edges)
  label_x <- 1.17

  p <- ggplot() +
    geom_rect(
      data = pd$bands,
      aes(xmin = 1.00, xmax = 1.85, ymin = ymin, ymax = ymax, fill = Tissue),
      alpha = 0.055, colour = NA
    ) +
    geom_curve(
      data = pd$edges,
      aes(
        x = x, y = hub_y, xend = xend, yend = node_y,
        colour = Pair, alpha = EdgeAlpha,
        linetype = EvidenceClass
      ),
      curvature = 0.12,
      linewidth = 0.58,
      lineend = "round"
    ) +
    geom_label(
      data = pd$hubs,
      aes(x = x, y = y, label = PairLabel, fill = Pair),
      colour = "white", size = 2.55, fontface = "bold",
      label.size = 0, label.padding = unit(0.18, "lines"),
      show.legend = FALSE
    ) +
    geom_point(
      data = pd$nodes,
      aes(
        x = 1.08, y = y, shape = Omics, fill = Tissue,
        colour = Tissue, size = NodeSignificance
      ),
      stroke = 0.65
    ) +
    geom_text(
      data = pd$nodes,
      aes(x = label_x, y = y, label = Label),
      hjust = 0, size = 2.25, colour = "black"
    ) +
    geom_text(
      data = pd$bands,
      aes(x = 1.015, y = ymax + 0.12, label = Tissue, colour = Tissue),
      hjust = 0, vjust = 0, size = 2.8, fontface = "bold",
      show.legend = FALSE
    ) +
    scale_colour_manual(values = c(pair_colors, tissue_dark), breaks = pair_levels) +
    scale_fill_manual(values = c(pair_colors, tissue_colors)) +
    scale_shape_manual(values = omics_shapes) +
    scale_size_continuous(
      name = expression(-log[10]("BH-FDR")),
      range = c(2.0, 4.4), limits = c(1.3, 4.0),
      breaks = c(2, 3, 4)
    ) +
    scale_alpha_identity() +
    scale_linetype_manual(
      name = "Evidence",
      values = c(`Primary non-DNAm` = "solid", `Exploratory DNAm` = "22"),
      drop = FALSE
    ) +
    guides(
      colour = guide_legend(
        title = "Family pair", override.aes = list(linewidth = 1.1, alpha = 1)
      ),
      fill = "none",
      shape = guide_legend(title = "Molecular layer", override.aes = list(size = 3)),
      size = guide_legend(order = 3),
      linetype = guide_legend(order = 4)
    ) +
    coord_cartesian(xlim = c(-0.04, 1.86), clip = "off") +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6.5, colour = "black", margin = margin(b = 5)),
      legend.position = "right",
      legend.title = element_text(size = 6.2, face = "bold"),
      legend.text = element_text(size = 5.8),
      legend.key.height = unit(3.4, "mm"),
      legend.key.width = unit(5.0, "mm"),
      plot.margin = margin(5, 18, 5, 5)
    )

  list(plot = p, data = pd)
}

plot_family_network_stress <- function(edges) {
  molecule_nodes <- edges %>%
    group_by(NodeID, Dataset, Tissue, Omics, FeatureID, Gene, DisplayLabel, EvidenceClass) %>%
    summarise(
      MinBH_FDR = min(BH_FDR, na.rm = TRUE),
      MaxAbsCorrelation = max(abs(Correlation), na.rm = TRUE),
      PairCount = n_distinct(Pair),
      .groups = "drop"
    ) %>%
    mutate(
      NodeType = "Molecule",
      Label = if_else(
        Omics == "DNA methylation",
        paste0(DisplayLabel, "-linked DNAm"),
        DisplayLabel
      ),
      NodeSignificance = pmin(pmax(-log10(MinBH_FDR), 1.3), 4.0)
    )

  hub_nodes <- tibble(
    NodeID = paste0("PairHub_", pair_levels),
    Dataset = NA_character_,
    Tissue = NA_character_,
    Omics = NA_character_,
    FeatureID = NA_character_,
    Gene = NA_character_,
    DisplayLabel = pair_levels,
    EvidenceClass = "Family-pair hub",
    MinBH_FDR = NA_real_,
    MaxAbsCorrelation = NA_real_,
    PairCount = NA_integer_,
    NodeType = "Family pair",
    Label = pair_levels,
    NodeSignificance = NA_real_,
    Pair = factor(pair_levels, levels = pair_levels)
  )

  graph_edges <- edges %>%
    transmute(
      from = paste0("PairHub_", as.character(Pair)),
      to = NodeID,
      NodeID = to,
      Pair = factor(Pair, levels = pair_levels),
      Correlation, P_Value, BH_FDR, NPairs, StabilityFlag, EvidenceClass,
      EdgeStrength = abs(Correlation)
    )

  graph_vertices <- bind_rows(
    molecule_nodes %>% mutate(Pair = factor(NA_character_, levels = pair_levels)),
    hub_nodes
  )

  graph <- igraph::graph_from_data_frame(
    graph_edges,
    directed = FALSE,
    vertices = graph_vertices %>% rename(name = NodeID)
  )
  set.seed(layout_seed)
  layout <- ggraph::create_layout(graph, layout = "stress")
  layout_table <- as_tibble(layout) %>% mutate(NodeID = as.character(name))

  molecule_layout <- layout_table %>%
    filter(NodeType == "Molecule") %>%
    mutate(
      Tissue = factor(Tissue, levels = c("Serum", "Adipose", "Muscle")),
      Omics = factor(
        Omics,
        levels = c("Metabolomics", "Proteomics", "Transcriptomics", "DNA methylation")
      )
    )
  hub_layout <- layout_table %>%
    filter(NodeType == "Family pair") %>%
    mutate(Pair = factor(Pair, levels = pair_levels))

  hull_members <- bind_rows(
    graph_edges %>% transmute(Pair, NodeID = to),
    graph_edges %>% distinct(Pair, from) %>% transmute(Pair, NodeID = from)
  ) %>%
    distinct(Pair, NodeID) %>%
    left_join(layout_table %>% select(NodeID, x, y), by = "NodeID") %>%
    mutate(Pair = factor(Pair, levels = pair_levels))

  current_label_cap <- if (any(molecule_layout$Omics == "DNA methylation")) {
    exploratory_label_cap
  } else {
    primary_label_cap
  }

  label_table <- molecule_layout %>%
    arrange(MinBH_FDR, desc(MaxAbsCorrelation), Label) %>%
    slice_head(n = current_label_cap)

  group_colors <- c(
    Serum = "#E9B85F",
    Adipose = "#E98773",
    Muscle = "#72A8CC"
  )
  modality_shapes <- c(
    Metabolomics = 23,
    Proteomics = 21,
    Transcriptomics = 24,
    `DNA methylation` = 25
  )

  p <- ggraph::ggraph(layout) +
    ggraph::geom_edge_link(
      aes(
        edge_colour = Pair,
        edge_alpha = EdgeStrength,
        edge_width = EdgeStrength,
        edge_linetype = EvidenceClass
      ),
      lineend = "round"
    ) +
    ggforce::geom_mark_hull(
      data = hull_members,
      aes(
        x = x, y = y, group = Pair,
        fill = Pair
      ),
      concavity = 1.6,
      expand = grid::unit(1.8, "mm"),
      radius = grid::unit(1.8, "mm"),
      alpha = 0.065,
      colour = "#777777",
      linetype = "dashed",
      linewidth = 0.25,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = molecule_layout,
      aes(
        x = x, y = y, size = NodeSignificance,
        fill = Tissue, shape = Omics
      ),
      colour = "white", stroke = 0.50,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = hub_layout,
      aes(x = x, y = y, fill = Pair),
      shape = 21, size = 5.6, colour = "white", stroke = 0.70,
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    geom_text(
      data = hub_layout,
      aes(x = x, y = y, label = DisplayLabel),
      family = "Arial", fontface = "bold", size = 5.8 / ggplot2::.pt,
      colour = "white", inherit.aes = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = label_table,
      aes(x = x, y = y, label = Label),
      family = "Arial", fontface = "bold", size = 5.6 / ggplot2::.pt,
      colour = "black", box.padding = 0.20, point.padding = 0.10,
      min.segment.length = Inf, segment.colour = NA,
      max.overlaps = Inf, seed = layout_seed,
      force = 0.40, force_pull = 0.60,
      inherit.aes = FALSE
    ) +
    ggraph::scale_edge_colour_manual(
      values = pair_colors,
      labels = pair_labels,
      name = "Family pair"
    ) +
    ggraph::scale_edge_alpha_continuous(range = c(0.18, 0.50), guide = "none") +
    ggraph::scale_edge_width_continuous(range = c(0.22, 0.72), guide = "none") +
    ggraph::scale_edge_linetype_manual(
      values = c(`Primary non-DNAm` = "solid", `Exploratory DNAm` = "22"),
      name = "Evidence",
      drop = TRUE
    ) +
    scale_fill_manual(
      values = c(group_colors, pair_colors),
      breaks = names(group_colors),
      name = "Tissue"
    ) +
    scale_shape_manual(
      values = modality_shapes,
      name = "Molecular layer",
      drop = TRUE
    ) +
    scale_size_continuous(
      name = expression(-log[10]("BH-FDR")),
      range = c(2.0, 5.1), limits = c(1.3, 4.0), breaks = c(2, 3, 4)
    ) +
    coord_equal(clip = "off") +
    theme_void(base_family = "Arial", base_size = 7) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.title = element_text(face = "bold", size = 6.5, colour = "black"),
      legend.text = element_text(size = 5.8, colour = "black"),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.spacing.y = grid::unit(1.1, "mm"),
      legend.box.margin = margin(0, 0, 1.5, 0, unit = "mm"),
      plot.margin = margin(3, 3, 4, 3, unit = "mm")
    ) +
    guides(
      fill = guide_legend(
        order = 1, nrow = 1, byrow = TRUE,
        override.aes = list(shape = 21, size = 2.8, alpha = 1)
      ),
      shape = guide_legend(
        order = 2, nrow = 1, byrow = TRUE,
        override.aes = list(size = 2.8, fill = "#D9D9D9")
      ),
      edge_colour = guide_legend(
        order = 3, nrow = 1,
        override.aes = list(edge_width = 0.70, edge_alpha = 0.85)
      ),
      edge_linetype = guide_legend(
        order = 4, nrow = 1,
        override.aes = list(edge_width = 0.65, edge_alpha = 0.85)
      ),
      size = guide_legend(order = 5, nrow = 1)
    )

  list(
    plot = p,
    data = list(
      nodes = molecule_layout,
      hubs = hub_layout,
      edges = graph_edges,
      hulls = hull_members
    )
  )
}

plot_family_network_triangle <- function(edges, label_cap) {
  triangle_pair_levels <- c("DF", "DM", "MF")
  signature_layout <- tibble::tribble(
    ~PairSignature, ~CenterX, ~CenterY, ~Radius, ~AngleOffset,
    "DF",          -1.18, -0.68, 0.55, 2.70,
    "DM",           1.18, -0.68, 0.55, 0.45,
    "MF",           0.00,  1.18, 0.55, 1.55,
    "DF+DM",        0.00, -0.63, 0.26, 0.20,
    "DM+MF",        0.56,  0.25, 0.26, 1.10,
    "DF+MF",       -0.56,  0.25, 0.26, 2.10,
    "DF+DM+MF",     0.00,  0.12, 0.34, 0.70
  )
  hub_positions <- tibble::tribble(
    ~Pair, ~HubX, ~HubY,
    "DF", -1.38, -0.88,
    "DM",  1.38, -0.88,
    "MF",  0.00,  1.48
  ) %>%
    mutate(Pair = factor(Pair, levels = pair_levels))

  molecule_nodes <- edges %>%
    group_by(NodeID, Dataset, Tissue, Omics, FeatureID, Gene, DisplayLabel, EvidenceClass) %>%
    summarise(
      MinBH_FDR = min(BH_FDR, na.rm = TRUE),
      MaxAbsCorrelation = max(abs(Correlation), na.rm = TRUE),
      PairCount = n_distinct(Pair),
      PairSignature = paste(
        triangle_pair_levels[triangle_pair_levels %in% as.character(Pair)],
        collapse = "+"
      ),
      AnyLOFOStable = any(LOFOStable),
      .groups = "drop"
    ) %>%
    mutate(
      Tissue = factor(Tissue, levels = c("Serum", "Adipose", "Muscle")),
      Omics = factor(
        Omics,
        levels = c("Metabolomics", "Proteomics", "Transcriptomics", "DNA methylation")
      ),
      Label = case_when(
        Omics == "DNA methylation" & !is.na(Gene) ~ Gene,
        Omics == "DNA methylation" ~ NA_character_,
        TRUE ~ DisplayLabel
      ),
      PairSupport = factor(
        PairCount, levels = c(1, 2, 3), labels = c("1 pair", "2 pairs", "3 pairs")
      ),
      SupportClass = if_else(PairCount >= 2, "2-3 pairs", "1 pair")
    ) %>%
    left_join(signature_layout, by = "PairSignature") %>%
    group_by(PairSignature) %>%
    arrange(desc(PairCount), MinBH_FDR, desc(MaxAbsCorrelation), NodeID, .by_group = TRUE) %>%
    mutate(
      PositionIndex = row_number(),
      SignatureN = n(),
      Theta = AngleOffset + (PositionIndex - 1) * pi * (3 - sqrt(5)),
      RadialDistance = Radius * sqrt((PositionIndex - 0.35) / SignatureN),
      x = CenterX + RadialDistance * cos(Theta),
      y = CenterY + RadialDistance * sin(Theta)
    ) %>%
    ungroup()

  edge_plot <- edges %>%
    left_join(
      molecule_nodes %>% select(NodeID, NodeX = x, NodeY = y, PairCount),
      by = "NodeID"
    ) %>%
    left_join(hub_positions, by = "Pair") %>%
    mutate(
      Pair = factor(Pair, levels = pair_levels),
      StabilityClass = factor(StabilityClass, levels = c("LOFO-stable", "FDR-only")),
      EdgeSign = factor(if_else(Correlation >= 0, "Positive", "Negative"),
                        levels = c("Negative", "Positive")),
      EdgeWidth = scales::rescale(abs(Correlation), to = c(0.08, 0.34), from = c(0.45, 1.00)),
      EdgeAlpha = if_else(LOFOStable, 0.14, 0.035)
    ) %>%
    arrange(LOFOStable, PairCount, BH_FDR)

  unique_label_nodes <- molecule_nodes %>%
    filter(!is.na(Label), Label != "") %>%
    arrange(desc(PairCount), MinBH_FDR, desc(MaxAbsCorrelation), Label) %>%
    distinct(Label, .keep_all = TRUE)
  mandatory_labels <- unique_label_nodes %>%
    filter(Label %in% c("PIEZO1", "PQLC3")) %>%
    mutate(LabelPriority = 0L)
  multipair_by_layer <- unique_label_nodes %>%
    filter(PairCount >= 2) %>%
    group_by(Tissue, Omics) %>%
    slice_min(MinBH_FDR, n = 4, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(LabelPriority = 1L)
  multipair_global <- unique_label_nodes %>%
    filter(PairCount >= 2) %>%
    slice_min(MinBH_FDR, n = label_cap, with_ties = FALSE) %>%
    mutate(LabelPriority = 2L)
  single_by_layer <- unique_label_nodes %>%
    filter(PairCount == 1) %>%
    group_by(Tissue, Omics) %>%
    slice_min(MinBH_FDR, n = 2, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(LabelPriority = 3L)
  label_table <- bind_rows(
    mandatory_labels, multipair_by_layer, multipair_global, single_by_layer
  ) %>%
    arrange(LabelPriority, desc(PairCount), MinBH_FDR, desc(MaxAbsCorrelation), Label) %>%
    distinct(NodeID, .keep_all = TRUE) %>%
    slice_head(n = label_cap)

  molecule_nodes <- molecule_nodes %>%
    mutate(IsLabeled = NodeID %in% label_table$NodeID)

  pair_hulls <- bind_rows(
    molecule_nodes %>%
      filter(PairSignature %in% triangle_pair_levels) %>%
      transmute(Pair = factor(PairSignature, levels = pair_levels), x, y),
    hub_positions %>% transmute(Pair, x = HubX, y = HubY)
  )

  triangle_outline <- bind_rows(
    hub_positions %>% filter(as.character(Pair) == "DF"),
    hub_positions %>% filter(as.character(Pair) == "MF"),
    hub_positions %>% filter(as.character(Pair) == "DM"),
    hub_positions %>% filter(as.character(Pair) == "DF")
  )
  group_colors <- c(
    Serum = "#E9B85F",
    Adipose = "#E98773",
    Muscle = "#72A8CC"
  )
  modality_shapes <- c(
    Metabolomics = 23,
    Proteomics = 21,
    Transcriptomics = 24,
    `DNA methylation` = 25
  )
  edge_colors <- c(Negative = "#5D91B8", Positive = "#D96B63")
  ordinary_node_size <- if (nrow(molecule_nodes) > 500) {
    0.62
  } else if (nrow(molecule_nodes) > 100) {
    0.90
  } else {
    1.55
  }
  multipair_node_size <- if (nrow(molecule_nodes) > 500) {
    1.45
  } else if (nrow(molecule_nodes) > 100) {
    1.75
  } else {
    2.15
  }

  p <- ggplot() +
    geom_path(
      data = triangle_outline,
      aes(x = HubX, y = HubY),
      colour = "#8A8A8A", linewidth = 0.28, linetype = "22",
      inherit.aes = FALSE
    ) +
    ggforce::geom_mark_hull(
      data = pair_hulls,
      aes(x = x, y = y, group = Pair),
      concavity = 1.5,
      expand = grid::unit(1.2, "mm"),
      radius = grid::unit(1.2, "mm"),
      fill = NA, colour = "#777777", linetype = "22", linewidth = 0.23,
      show.legend = FALSE, inherit.aes = FALSE
    ) +
    geom_segment(
      data = edge_plot,
      aes(
        x = HubX, y = HubY, xend = NodeX, yend = NodeY,
        colour = EdgeSign,
        linewidth = EdgeWidth
      ),
      alpha = if (nrow(molecule_nodes) > 500) 0.055 else 0.20,
      lineend = "round", inherit.aes = FALSE
    ) +
    geom_point(
      data = molecule_nodes %>% filter(!IsLabeled, PairCount == 1),
      aes(
        x = x, y = y, fill = Tissue, shape = Omics
      ),
      size = ordinary_node_size, alpha = 0.48,
      colour = "white", stroke = 0.18, inherit.aes = FALSE
    ) +
    geom_point(
      data = molecule_nodes %>% filter(!IsLabeled, PairCount >= 2),
      aes(x = x, y = y, fill = Tissue, shape = Omics),
      size = multipair_node_size, alpha = 0.82,
      colour = "white", stroke = 0.28, inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_point(
      data = molecule_nodes %>% filter(IsLabeled),
      aes(x = x, y = y, fill = Tissue, shape = Omics),
      size = 3.35, alpha = 1,
      colour = "white", stroke = 0.48, inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_point(
      data = hub_positions,
      aes(x = HubX, y = HubY),
      shape = 21, size = 5.8, fill = pair_colors[as.character(hub_positions$Pair)],
      colour = "white", stroke = 0.70, inherit.aes = FALSE
    ) +
    geom_text(
      data = hub_positions,
      aes(x = HubX, y = HubY, label = Pair),
      family = "Arial", fontface = "bold", size = 5.8 / ggplot2::.pt,
      colour = "white", inherit.aes = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = label_table,
      aes(x = x, y = y, label = Label),
      family = "Arial", fontface = "bold", size = 5.6 / ggplot2::.pt,
      colour = "black", box.padding = 0.25, point.padding = 0.10,
      min.segment.length = Inf, segment.colour = NA,
      max.overlaps = Inf, seed = layout_seed,
      force = 1.10, force_pull = 0.20,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(values = edge_colors, name = "Association direction") +
    scale_fill_manual(values = group_colors, name = "Tissue", drop = TRUE) +
    scale_shape_manual(values = modality_shapes, name = "Molecular layer", drop = TRUE) +
    scale_linewidth_identity() +
    guides(
      fill = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 21, size = 2.8, alpha = 1)),
      shape = guide_legend(order = 2, nrow = 1, override.aes = list(size = 2.8, fill = "#D9D9D9", alpha = 1)),
      colour = guide_legend(order = 3, nrow = 1, override.aes = list(linewidth = 0.70, alpha = 0.85))
    ) +
    coord_equal(xlim = c(-1.95, 1.95), ylim = c(-1.48, 2.00), clip = "off") +
    theme_void(base_family = "Arial", base_size = 7) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.title = element_text(face = "bold", size = 6.4, colour = "black"),
      legend.text = element_text(size = 5.7, colour = "black"),
      legend.key.height = grid::unit(3.3, "mm"),
      legend.spacing.y = grid::unit(0.9, "mm"),
      legend.box.margin = margin(0, 0, 1.5, 0, unit = "mm"),
      plot.margin = margin(3, 3, 4, 3, unit = "mm")
    )

  list(
    plot = p,
    data = list(
      nodes = molecule_nodes,
      hubs = hub_positions,
      edges = edge_plot,
      hulls = triangle_outline,
      labels = label_table
    )
  )
}

plot_family_network_stress_allfdr <- function(edges, label_cap = 40L) {
  molecule_nodes <- edges %>%
    group_by(NodeID, Dataset, Tissue, Omics, FeatureID, Gene, DisplayLabel, EvidenceClass) %>%
    summarise(
      MinBH_FDR = min(BH_FDR, na.rm = TRUE),
      MaxAbsCorrelation = max(abs(Correlation), na.rm = TRUE),
      PairCount = n_distinct(Pair),
      AnyLOFOStable = any(LOFOStable),
      .groups = "drop"
    ) %>%
    mutate(
      NodeType = "Molecule",
      Tissue = factor(Tissue, levels = c("Serum", "Adipose", "Muscle")),
      Omics = factor(
        Omics,
        levels = c("Metabolomics", "Proteomics", "Transcriptomics", "DNA methylation")
      ),
      Label = case_when(
        Omics == "DNA methylation" & !is.na(Gene) ~ Gene,
        Omics == "DNA methylation" ~ NA_character_,
        TRUE ~ DisplayLabel
      )
    )

  hub_nodes <- tibble(
    NodeID = paste0("PairHub_", pair_levels),
    Dataset = NA_character_, Tissue = NA_character_, Omics = NA_character_,
    FeatureID = NA_character_, Gene = NA_character_, DisplayLabel = pair_levels,
    EvidenceClass = "Family-pair hub", MinBH_FDR = NA_real_,
    MaxAbsCorrelation = NA_real_, PairCount = NA_integer_,
    AnyLOFOStable = NA, NodeType = "Family pair", Label = pair_levels,
    Pair = factor(pair_levels, levels = pair_levels)
  )

  graph_edges <- edges %>%
    transmute(
      from = paste0("PairHub_", as.character(Pair)),
      to = NodeID,
      NodeID = to,
      Pair = factor(Pair, levels = pair_levels),
      Correlation, P_Value, BH_FDR, NPairs, StabilityFlag,
      CoverageOK, LOFOStable, StabilityClass, EvidenceClass,
      EdgeSign = factor(if_else(Correlation >= 0, "Positive", "Negative"),
                        levels = c("Negative", "Positive")),
      EdgeWeight = abs(Correlation)
    )

  graph_vertices <- bind_rows(
    molecule_nodes %>% mutate(Pair = factor(NA_character_, levels = pair_levels)),
    hub_nodes
  )
  graph <- igraph::graph_from_data_frame(
    graph_edges, directed = FALSE,
    vertices = graph_vertices %>% rename(name = NodeID)
  )
  set.seed(layout_seed)
  layout <- ggraph::create_layout(graph, layout = "stress")
  layout_table <- as_tibble(layout) %>% mutate(NodeID = as.character(name))

  molecule_layout <- layout_table %>%
    filter(NodeType == "Molecule") %>%
    mutate(
      Tissue = factor(Tissue, levels = c("Serum", "Adipose", "Muscle")),
      Omics = factor(
        Omics,
        levels = c("Metabolomics", "Proteomics", "Transcriptomics", "DNA methylation")
      )
    )
  hub_layout <- layout_table %>%
    filter(NodeType == "Family pair") %>%
    mutate(Pair = factor(Pair, levels = pair_levels))

  unique_labels <- molecule_layout %>%
    filter(!is.na(Label), Label != "") %>%
    arrange(desc(PairCount), MinBH_FDR, desc(MaxAbsCorrelation), Label) %>%
    distinct(Label, .keep_all = TRUE)
  mandatory_labels <- unique_labels %>%
    filter(Label %in% c("PIEZO1", "PQLC3")) %>%
    mutate(LabelPriority = 0L)
  protein_labels <- molecule_layout %>%
    filter(Omics == "Proteomics", !is.na(Label), Label != "") %>%
    arrange(desc(PairCount), MinBH_FDR, desc(MaxAbsCorrelation), Label) %>%
    distinct(Label, .keep_all = TRUE) %>%
    mutate(LabelPriority = 1L)
  transcript_labels <- unique_labels %>%
    filter(Omics == "Transcriptomics") %>%
    mutate(LabelPriority = 2L)
  multipair_dnam_labels <- unique_labels %>%
    filter(Omics == "DNA methylation", PairCount >= 2) %>%
    slice_min(MinBH_FDR, n = 8, with_ties = FALSE) %>%
    mutate(LabelPriority = 3L)
  metabolite_labels <- unique_labels %>%
    filter(Omics == "Metabolomics") %>%
    slice_min(MinBH_FDR, n = 4, with_ties = FALSE) %>%
    mutate(LabelPriority = 4L)
  required_labels <- bind_rows(mandatory_labels, protein_labels) %>%
    arrange(LabelPriority, desc(PairCount), MinBH_FDR, desc(MaxAbsCorrelation), Label) %>%
    distinct(NodeID, .keep_all = TRUE)
  optional_labels <- bind_rows(
    transcript_labels, multipair_dnam_labels, metabolite_labels
  ) %>%
    filter(
      !NodeID %in% required_labels$NodeID,
      !Label %in% required_labels$Label
    ) %>%
    arrange(LabelPriority, desc(PairCount), MinBH_FDR, desc(MaxAbsCorrelation), Label) %>%
    distinct(NodeID, .keep_all = TRUE) %>%
    slice_head(n = max(0L, label_cap - nrow(required_labels)))
  label_table <- bind_rows(required_labels, optional_labels) %>%
    arrange(LabelPriority, desc(PairCount), MinBH_FDR, desc(MaxAbsCorrelation), Label)

  molecule_layout <- molecule_layout %>%
    mutate(IsLabeled = NodeID %in% label_table$NodeID)
  label_table <- molecule_layout %>% filter(IsLabeled)

  group_colors <- c(Serum = "#E9B85F", Adipose = "#E98773", Muscle = "#72A8CC")
  modality_shapes <- c(
    Metabolomics = 23, Proteomics = 21,
    Transcriptomics = 24, `DNA methylation` = 25
  )

  p <- ggraph::ggraph(layout) +
    ggraph::geom_edge_link(
      aes(edge_width = EdgeWeight, edge_alpha = EdgeWeight, edge_colour = Pair),
      lineend = "round"
    ) +
    geom_point(
      data = molecule_layout %>%
        filter(!IsLabeled, Omics == "DNA methylation", PairCount == 1),
      aes(x = x, y = y, fill = Tissue, shape = Omics),
      size = 0.46, alpha = 0.34, colour = "white", stroke = 0.10,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = molecule_layout %>%
        filter(!IsLabeled, !(Omics == "DNA methylation" & PairCount == 1), PairCount == 1),
      aes(x = x, y = y, fill = Tissue, shape = Omics),
      size = 0.82, alpha = 0.55, colour = "white", stroke = 0.16,
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    geom_point(
      data = molecule_layout %>% filter(!IsLabeled, PairCount >= 2),
      aes(x = x, y = y, fill = Tissue, shape = Omics),
      size = 1.35, alpha = 0.82, colour = "white", stroke = 0.25,
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    geom_point(
      data = molecule_layout %>% filter(IsLabeled),
      aes(x = x, y = y, fill = Tissue, shape = Omics),
      size = 3.10, alpha = 1, colour = "white", stroke = 0.48,
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    geom_point(
      data = hub_layout,
      aes(x = x, y = y),
      shape = 21, size = 5.6,
      fill = pair_colors[as.character(hub_layout$Pair)],
      colour = "white", stroke = 0.70,
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    geom_text(
      data = hub_layout,
      aes(x = x, y = y, label = DisplayLabel),
      family = "Arial", fontface = "bold", size = 5.8 / ggplot2::.pt,
      colour = "white", inherit.aes = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = label_table,
      aes(x = x, y = y, label = Label),
      family = "Arial", fontface = "bold", size = 5.5 / ggplot2::.pt,
      colour = "black", box.padding = 0.20, point.padding = 0.08,
      min.segment.length = Inf, segment.colour = NA,
      max.overlaps = Inf, seed = layout_seed,
      force = 0.75, force_pull = 0.45,
      inherit.aes = FALSE
    ) +
    ggraph::scale_edge_colour_manual(
      values = pair_colors,
      breaks = pair_levels,
      labels = unname(pair_labels[pair_levels]),
      name = "Family pair"
    ) +
    ggraph::scale_edge_alpha_continuous(range = c(0.025, 0.12), guide = "none") +
    ggraph::scale_edge_width_continuous(range = c(0.07, 0.34), guide = "none") +
    scale_fill_manual(values = group_colors, name = "Tissue", drop = TRUE) +
    scale_shape_manual(values = modality_shapes, name = "Molecular layer", drop = TRUE) +
    coord_equal(clip = "off") +
    theme_void(base_family = "Arial", base_size = 7) +
    theme(
      legend.position = "bottom", legend.box = "vertical", legend.box.just = "center",
      legend.title = element_text(face = "bold", size = 6.4, colour = "black"),
      legend.text = element_text(size = 5.7, colour = "black"),
      legend.key.height = grid::unit(3.3, "mm"),
      legend.spacing.y = grid::unit(0.9, "mm"),
      legend.box.margin = margin(0, 0, 1.5, 0, unit = "mm"),
      plot.margin = margin(3, 3, 4, 3, unit = "mm")
    ) +
    guides(
      fill = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 21, size = 2.8, alpha = 1)),
      shape = guide_legend(order = 2, nrow = 1, override.aes = list(size = 2.8, fill = "#D9D9D9", alpha = 1)),
      edge_colour = guide_legend(order = 3, nrow = 1, override.aes = list(edge_width = 0.70, edge_alpha = 0.85))
    )

  list(
    plot = p,
    data = list(
      nodes = molecule_layout,
      hubs = hub_layout,
      edges = graph_edges,
      hulls = tibble(),
      labels = label_table
    )
  )
}

primary_result <- plot_family_network_triangle(primary_edges, primary_label_cap)

hybrid_result <- plot_family_network_triangle(hybrid_edges, hybrid_label_cap)

exploratory_result <- plot_family_network_stress_allfdr(exploratory_edges, label_cap = 40L)

dataset_key <- tibble::tribble(
  ~Dataset, ~Tissue, ~Omics, ~DatasetLabel,
  "Serum_Metabonomics", "Serum", "Metabolomics", "Serum metabolomics",
  "Serum_Proteomics", "Serum", "Proteomics", "Serum proteomics",
  "Adipose_Proteomics", "Adipose", "Proteomics", "Adipose proteomics",
  "Adipose_Microarray", "Adipose", "Transcriptomics", "Adipose transcriptomics",
  "Adipose_Methylation", "Adipose", "DNA methylation", "Adipose DNAm",
  "Muscle_Proteomics", "Muscle", "Proteomics", "Muscle proteomics",
  "Muscle_Microarray", "Muscle", "Transcriptomics", "Muscle transcriptomics",
  "Muscle_Methylation", "Muscle", "DNA methylation", "Muscle DNAm"
)

count_observed <- all_edges %>%
  mutate(
    Significant = is.finite(BH_FDR) & BH_FDR < family_fdr,
    Stable = Significant & StabilityFlag == "Stable by prespecified LOFO diagnostics"
  ) %>%
  group_by(Dataset, Pair) %>%
  summarise(
    FDRSignificant = sum(Significant, na.rm = TRUE),
    LOFOStable = sum(Stable, na.rm = TRUE),
    .groups = "drop"
  )

count_table <- tidyr::crossing(
  dataset_key,
  Pair = factor(pair_levels, levels = pair_levels)
) %>%
  left_join(count_observed, by = c("Dataset", "Pair")) %>%
  mutate(
    FDRSignificant = coalesce(FDRSignificant, 0L),
    LOFOStable = coalesce(LOFOStable, 0L),
    DatasetLabel = factor(
      DatasetLabel,
      levels = rev(c(
        "Serum metabolomics", "Serum proteomics",
        "Adipose proteomics", "Adipose transcriptomics", "Adipose DNAm",
        "Muscle proteomics", "Muscle transcriptomics", "Muscle DNAm"
      ))
    )
  )

count_plot <- ggplot(count_table, aes(x = Pair, y = DatasetLabel, fill = log10(LOFOStable + 1))) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = LOFOStable), size = 2.55, colour = "black") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#DDEBE7", "#95CEBE", "#527EAF"),
    name = expression(log[10]("stable count" + 1))
  ) +
  scale_x_discrete(labels = pair_labels) +
  labs(
    title = "Stable pair-specific molecular similarity signals",
    subtitle = "Counts satisfy BH-FDR < 0.05 and prespecified leave-one-family-out diagnostics",
    x = NULL, y = NULL
  ) +
  theme_classic(base_family = "Arial", base_size = 7) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 6.5, face = "bold"),
    axis.text.y = element_text(size = 6.2, colour = "black"),
    plot.title = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 6.5),
    legend.title = element_text(size = 6.2),
    legend.text = element_text(size = 5.8),
    legend.position = "right"
  )

ggsave(
  file.path(output_dir, "FigS02_Family_Molecular_Network_Primary_NonDNAm_Candidate.pdf"),
  primary_result$plot, width = 155, height = 140, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "FigS02_Family_Molecular_Network_Hybrid_MultiPairDNAm_Candidate.pdf"),
  hybrid_result$plot, width = 155, height = 140, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "FigS02_Family_Molecular_Network_Exploratory_DNAm_Candidate.pdf"),
  exploratory_result$plot, width = 155, height = 140, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "FigS02_Family_Molecular_Similarity_Evidence_Counts.pdf"),
  count_plot, width = 125, height = 92, units = "mm", device = cairo_pdf
)

selection_parameters <- tibble(
  Parameter = c(
    "Analysis ID", "Time point", "Family-pair BH-FDR threshold",
    "Display inclusion", "LOFO reporting", "Minimum pair coverage for LOFO-stable edge",
    "Primary network", "Hybrid network", "Exploratory network", "Prespecified DNAm candidate",
    "Interpretation boundary"
  ),
  Value = c(
    analysis_id, "Pre-exercise", as.character(family_fdr),
    "At least one family pair with BH-FDR < 0.05",
    "Recorded in the source workbook; not used as a visual edge encoding",
    paste0(minimum_fraction_of_platform_max_pairs * 100, "% of platform-pair maximum"),
    "All non-methylation features meeting display inclusion",
    "All non-methylation features plus DNAm features significant in at least two family pairs",
    "All molecular layers meeting display inclusion; multi-pair features prioritized",
    "PIEZO1-linked DNAm included and labeled; interpretation remains exploratory because of documented design-version block confounding",
    "Pair-specific significance is not a direct difference between family-pair correlations"
  )
)

input_manifest <- tibble(
  InputRole = "Family similarity audit",
  InputFile = file.path(
    "fig02", "Family_Similarity_Audit",
    "All_FDR_Significant_Family_Features.csv"
  ),
  SHA256 = digest::digest(file = input_file, algo = "sha256", serialize = FALSE)
)

source_workbook <- list(
  PrimaryNodes = primary_result$data$nodes,
  PrimaryPairEdges = primary_result$data$edges %>%
    select(NodeID, Pair, Correlation, P_Value, BH_FDR, NPairs, StabilityFlag,
           CoverageOK, LOFOStable, StabilityClass, EvidenceClass),
  HybridNodes = hybrid_result$data$nodes,
  HybridPairEdges = hybrid_result$data$edges %>%
    select(NodeID, Pair, Correlation, P_Value, BH_FDR, NPairs, StabilityFlag,
           CoverageOK, LOFOStable, StabilityClass, EvidenceClass),
  HybridLabels = hybrid_result$data$labels,
  ExploratoryNodes = exploratory_result$data$nodes,
  ExploratoryPairEdges = exploratory_result$data$edges %>%
    select(NodeID, Pair, Correlation, P_Value, BH_FDR, NPairs, StabilityFlag,
           CoverageOK, LOFOStable, StabilityClass, EvidenceClass),
  PrimaryLabels = primary_result$data$labels,
  ExploratoryLabels = exploratory_result$data$labels,
  CompleteEvidenceCounts = count_table,
  SelectionParameters = selection_parameters,
  InputManifest = input_manifest
)
writexl::write_xlsx(
  source_workbook,
  file.path(output_dir, "FigS02_Family_Molecular_Network_Source_Data.xlsx")
)

run_summary <- c(
  paste0("Analysis ID: ", analysis_id),
  paste0("Input SHA256: ", input_manifest$SHA256),
  paste0("Primary selected feature-pair edges: ", nrow(primary_edges)),
  paste0("Primary unique molecular nodes: ", n_distinct(primary_edges$NodeID)),
  paste0("Hybrid selected feature-pair edges: ", nrow(hybrid_edges)),
  paste0("Hybrid unique molecular nodes: ", n_distinct(hybrid_edges$NodeID)),
  paste0("Hybrid DNA-methylation nodes: ", n_distinct(hybrid_edges$NodeID[as.character(hybrid_edges$Omics) == "DNA methylation"])),
  paste0("Exploratory selected feature-pair edges: ", nrow(exploratory_edges)),
  paste0("Exploratory unique molecular nodes: ", n_distinct(exploratory_edges$NodeID)),
  paste0("Primary nodes significant in >=2 family pairs: ", sum(primary_result$data$nodes$PairCount >= 2)),
  paste0("Exploratory nodes significant in >=2 family pairs: ", sum(exploratory_result$data$nodes$PairCount >= 2)),
  paste0("Exploratory nodes significant in all 3 family pairs: ", sum(exploratory_result$data$nodes$PairCount == 3)),
  "Primary graph excludes DNA methylation.",
  "Exploratory graph includes every feature with BH-FDR < 0.05 in at least one family pair.",
  "Exploratory graph edge colour denotes DM, DF or MF family-pair similarity; correlation sign remains available in the source workbook.",
  "PIEZO1-linked DNAm is included and labeled but remains exploratory because of documented design-version block confounding.",
  "Exploratory stress-layout position is for visualization and is not a biological distance metric.",
  "Pair-specific significance must not be interpreted as a direct DM-versus-DF-versus-MF difference."
)
writeLines(run_summary, file.path(output_dir, "RunLog.txt"), useBytes = TRUE)

message(paste(run_summary, collapse = "\n"))
