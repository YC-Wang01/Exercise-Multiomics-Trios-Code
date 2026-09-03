# Candidate Figure 5B: hsCRP-centred adipose-proteomic association networks.
#
# Direct hsCRP-protein spokes use the current family-aware limma results.
# Protein-protein edges are exploratory role-adjusted Spearman correlations
# with BH correction and leave-one-family-out sign-stability assessment.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script) || is.na(.local_script)) {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  .local_script <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else NA_character_
}
.script_dir <- if (!is.na(.local_script)) {
  dirname(normalizePath(.local_script, winslash = "/", mustWork = TRUE))
} else {
  file.path(getwd(), "code")
}
source(file.path(.script_dir, "Fig00_Config.R"))

p_load(
  dplyr, tidyr, tibble, stringr, readr, purrr, ggplot2, igraph,
  ggrepel, patchwork, scales, openxlsx
)

set.seed(20260829L)

analysis_id <- "Fig05B_hsCRP_Adipose_Association_Network_HalfWidth_2026-08-29"
output_dir <- path_project("fig05", "analysis", analysis_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

formal_file <- path_project(
  "results", "Unified_Pre_Post_Response_Association_Audit_2026-08-28",
  "Unified_NonMethylation_Pre_Post_Response_Audit.rds"
)
clinical_file <- path_analysis_ready("Metadata", "Clinical_Observed_Primary.csv")
sample_sheet_file <- path_analysis_ready("Metadata", "Project_Sample_Sheet.csv")
required_files <- c(formal_file, clinical_file, sample_sheet_file)
if (!all(file.exists(required_files))) {
  stop("Missing required input: ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}

formal <- readRDS(formal_file)
feature_results <- formal$feature_results

clinical <- read_csv(clinical_file, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = factor(Role, levels = c("Daughter", "Mother", "Father")),
    hsCRP = suppressWarnings(as.numeric(CRP_84m_G))
  ) %>%
  distinct(Clinical_Subject_ID, .keep_all = TRUE) %>%
  mutate(
    hsCRP_Z = if_else(
      is.finite(hsCRP) & hsCRP > 0,
      qnorm((rank(hsCRP, na.last = "keep", ties.method = "average") - 0.5) /
              sum(is.finite(hsCRP) & hsCRP > 0)),
      NA_real_
    )
  )

sample_sheet <- read_csv(sample_sheet_file, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = factor(Role, levels = c("Daughter", "Mother", "Father")),
    SampleID = tolower(as.character(SampleID)),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort)
  ) %>%
  filter(AnalysisSet == "Primary", Is_Primary_Cohort)

state_levels <- c("PreExercise", "ExerciseChange")
state_labels <- c(
  "PreExercise" = "Pre-exercise",
  "ExerciseChange" = "Exercise response"
)
dataset_labels <- c(
  "Adipose_Proteomics" = "Adipose"
)
tissue_colors <- c(
  "Adipose" = "#E98773",
  "Clinical" = "#79B8A9"
)
edge_colors <- c("Positive" = "#D96B63", "Negative" = "#5D91B8")
state_colors <- c("Pre-exercise" = "#6F9FBF", "Exercise response" = "#D77A72")

clean_label <- function(gene, feature_id) {
  gene <- as.character(gene)
  invalid <- is.na(gene) | gene == "" | gene == "NA" |
    str_detect(gene, "^[0-9.+-]+$")
  out <- if_else(invalid, as.character(feature_id), gene)
  str_replace_all(out, ";", "/")
}

# Select a bounded nominal context set while retaining every BH-FDR < 0.20 node.
candidate_results <- feature_results %>%
  filter(
    Predictor == "hsCRP",
    Dataset %in% names(dataset_labels),
    State %in% state_levels
  ) %>%
  mutate(
    TissueDisplay = unname(dataset_labels[Dataset]),
    DisplayLabel = clean_label(GeneSymbol, FeatureID),
    FDRNode = BH_FDR < 0.20
  ) %>%
  group_by(State, Dataset) %>%
  arrange(P_Value, .by_group = TRUE) %>%
  filter(P_Value < 0.01 | FDRNode) %>%
  slice_head(n = 16L) %>%
  ungroup()

fdr_nodes_all <- feature_results %>%
  filter(
    Predictor == "hsCRP", Dataset %in% names(dataset_labels),
    State %in% state_levels, BH_FDR < 0.20
  ) %>%
  mutate(
    TissueDisplay = unname(dataset_labels[Dataset]),
    DisplayLabel = clean_label(GeneSymbol, FeatureID),
    FDRNode = TRUE
  )
candidate_results <- bind_rows(candidate_results, fdr_nodes_all) %>%
  distinct(State, Dataset, FeatureID, .keep_all = TRUE)

read_state_values <- function(dataset, state, features) {
  matrix_file <- switch(
    dataset,
    "Adipose_Proteomics" = path_analysis_ready("Analysis_Adipose_Proteomics.csv"),
    "Muscle_Proteomics" = path_analysis_ready("Analysis_Muscle_Proteomics.csv"),
    stop("Unsupported dataset: ", dataset)
  )
  matrix_data <- read_analysis_matrix(matrix_file)
  features <- intersect(features, rownames(matrix_data))

  available <- sample_sheet %>%
    filter(
      Dataset == dataset, Timepoint %in% c("Pre", "Post3h"),
      SampleID %in% colnames(matrix_data)
    ) %>%
    inner_join(clinical, by = c("FamilyID", "Clinical_Subject_ID", "Role")) %>%
    distinct(Clinical_Subject_ID, Timepoint, .keep_all = TRUE)

  if (state == "PreExercise") {
    metadata <- available %>%
      filter(Timepoint == "Pre") %>%
      arrange(FamilyID, Role, Clinical_Subject_ID)
    values <- matrix_data[features, metadata$SampleID, drop = FALSE]
  } else {
    metadata <- available %>%
      select(FamilyID, Clinical_Subject_ID, Role, hsCRP_Z, Timepoint, SampleID) %>%
      pivot_wider(names_from = Timepoint, values_from = SampleID) %>%
      filter(!is.na(Pre), !is.na(Post3h)) %>%
      arrange(FamilyID, Role, Clinical_Subject_ID)
    values <- matrix_data[features, metadata$Post3h, drop = FALSE] -
      matrix_data[features, metadata$Pre, drop = FALSE]
  }

  role_adjusted <- t(vapply(seq_len(nrow(values)), function(i) {
    value <- as.numeric(values[i, ])
    keep <- is.finite(value) & !is.na(metadata$Role)
    output <- rep(NA_real_, length(value))
    if (sum(keep) >= 8L && nlevels(droplevels(metadata$Role[keep])) >= 2L) {
      output[keep] <- residuals(lm(value[keep] ~ droplevels(metadata$Role[keep])))
    }
    output
  }, numeric(ncol(values))))
  rownames(role_adjusted) <- rownames(values)
  colnames(role_adjusted) <- metadata$Clinical_Subject_ID
  list(values = role_adjusted, metadata = metadata)
}

correlation_tests <- list()
value_objects <- list()
for (state in state_levels) {
  for (dataset in names(dataset_labels)) {
    specs <- candidate_results %>% filter(State == state, Dataset == dataset)
    object <- read_state_values(dataset, state, specs$FeatureID)
    value_objects[[paste(state, dataset, sep = "__")]] <- object
    features <- rownames(object$values)
    if (length(features) < 2L) next
    pairs <- combn(features, 2, simplify = FALSE)
    correlation_tests[[paste(state, dataset, sep = "__")]] <- map_dfr(pairs, function(pair) {
      x <- as.numeric(object$values[pair[1], ])
      y <- as.numeric(object$values[pair[2], ])
      keep <- is.finite(x) & is.finite(y)
      n <- sum(keep)
      if (n < 10L) return(tibble())
      test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
      rho <- unname(test$estimate)
      p_value <- test$p.value
      families <- as.character(object$metadata$FamilyID[keep])
      lofo_rho <- map_dbl(unique(families), function(family) {
        subset <- families != family
        if (sum(subset) < 8L) return(NA_real_)
        suppressWarnings(cor(x[keep][subset], y[keep][subset], method = "spearman"))
      })
      stability <- mean(sign(lofo_rho[is.finite(lofo_rho)]) == sign(rho), na.rm = TRUE)
      tibble(
        State = state, Dataset = dataset, NodeA = pair[1], NodeB = pair[2],
        N = n, Rho = rho, P_Value = p_value, LOFO_Sign_Stability = stability
      )
    })
  }
}

all_correlations <- bind_rows(correlation_tests) %>%
  group_by(State) %>%
  mutate(BH_FDR = p.adjust(P_Value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    Direction = if_else(Rho >= 0, "Positive", "Negative"),
    Evidence = case_when(
      BH_FDR < 0.05 & abs(Rho) >= 0.55 & LOFO_Sign_Stability >= 0.80 ~ "BH-FDR < 0.05",
      BH_FDR < 0.20 & abs(Rho) >= 0.55 & LOFO_Sign_Stability >= 0.80 ~ "BH-FDR 0.05-0.20",
      P_Value < 0.01 & abs(Rho) >= 0.65 & LOFO_Sign_Stability >= 0.80 ~ "Nominal P < 0.01",
      TRUE ~ "Not displayed"
    )
  )

fdr_correlations <- all_correlations %>%
  filter(Evidence %in% c("BH-FDR < 0.05", "BH-FDR 0.05-0.20"))

nominal_correlations <- all_correlations %>%
  filter(Evidence == "Nominal P < 0.01") %>%
  group_by(State, Dataset) %>%
  arrange(P_Value, desc(abs(Rho)), .by_group = TRUE) %>%
  slice_head(n = 32L) %>%
  ungroup()

display_correlations <- bind_rows(fdr_correlations, nominal_correlations)

node_table <- candidate_results %>%
  transmute(
    State, Dataset, Tissue = TissueDisplay, FeatureID,
    NodeID = paste(Dataset, FeatureID, sep = "::"),
    DisplayLabel, Effect, P_Value, BH_FDR, Subjects, Families,
    FDRNode,
    EvidenceTier = case_when(
      BH_FDR < 0.05 ~ "BH-FDR < 0.05",
      BH_FDR < 0.10 ~ "BH-FDR 0.05-0.10",
      BH_FDR < 0.20 ~ "BH-FDR 0.10-0.20",
      TRUE ~ "Nominal P < 0.01"
    )
  )

display_correlations <- display_correlations %>%
  mutate(
    NodeA_ID = paste(Dataset, NodeA, sep = "::"),
    NodeB_ID = paste(Dataset, NodeB, sep = "::")
  )

make_network_panel <- function(state) {
  nodes <- node_table %>% filter(State == state)
  internal <- display_correlations %>% filter(State == state)

  coordinates <- map_dfr(names(dataset_labels), function(dataset) {
    tissue_nodes <- nodes %>% filter(Dataset == dataset)
    tissue_edges <- internal %>% filter(Dataset == dataset)
    graph <- graph_from_data_frame(
      tissue_edges %>% select(from = NodeA_ID, to = NodeB_ID),
      directed = FALSE,
      vertices = tissue_nodes %>% transmute(name = NodeID)
    )
    if (vcount(graph) == 1L) {
      layout <- matrix(c(0, 0), ncol = 2)
    } else if (ecount(graph) == 0L) {
      layout <- layout_in_circle(graph)
    } else {
      set.seed(20260829L + match(state, state_levels) * 100L + match(dataset, names(dataset_labels)))
      layout <- layout_with_fr(graph, niter = 2500, grid = "nogrid")
    }
    layout <- norm_coords(layout, xmin = -0.72, xmax = 0.72, ymin = -0.72, ymax = 0.72)
    tibble(NodeID = V(graph)$name, x = layout[, 1], y = layout[, 2] - 0.10)
  })

  hs_node <- tibble(
    State = state, Dataset = "Clinical", Tissue = "Clinical", FeatureID = "hsCRP",
    NodeID = "Clinical::hsCRP", DisplayLabel = "hsCRP", Effect = NA_real_,
    P_Value = NA_real_, BH_FDR = NA_real_, Subjects = NA_integer_, Families = NA_integer_,
    FDRNode = TRUE, EvidenceTier = "Anchor", x = 0, y = 1.10
  )
  nodes <- nodes %>%
    left_join(coordinates, by = "NodeID") %>%
    bind_rows(hs_node)

  degree_source <- internal %>%
    select(NodeA_ID, NodeB_ID) %>%
    pivot_longer(everything(), values_to = "NodeID") %>%
    count(NodeID, name = "Degree")
  nodes <- nodes %>%
    left_join(degree_source, by = "NodeID") %>%
    mutate(
      Degree = coalesce(Degree, 0L),
      LabelNode = FDRNode | (Tissue != "Clinical" & min_rank(desc(Degree)) <= 3L),
      NodeLabel = if_else(LabelNode, str_wrap(DisplayLabel, width = 12), ""),
      NodeSize = case_when(
        Tissue == "Clinical" ~ 5.3,
        FDRNode ~ 4.6,
        TRUE ~ rescale(Degree, to = c(1.85, 3.30))
      ),
      NodeAlpha = if_else(FDRNode | Tissue == "Clinical", 1, 0.62)
    )

  internal_plot <- internal %>%
    left_join(nodes %>% select(NodeA_ID = NodeID, x, y), by = "NodeA_ID") %>%
    left_join(nodes %>% select(NodeB_ID = NodeID, xend = x, yend = y), by = "NodeB_ID") %>%
    mutate(
      EdgeAlpha = case_when(
        Evidence == "BH-FDR < 0.05" ~ 0.76,
        Evidence == "BH-FDR 0.05-0.20" ~ 0.48,
        TRUE ~ 0.16
      ),
      EdgeWidth = rescale(abs(Rho), to = c(0.25, 0.85))
    )

  spokes <- nodes %>%
    filter(Tissue != "Clinical") %>%
    mutate(
      xend = 0, yend = 1.10,
      Direction = if_else(Effect >= 0, "Positive", "Negative"),
      EdgeAlpha = if_else(BH_FDR < 0.20, 0.78, 0.10),
      EdgeWidth = if_else(BH_FDR < 0.20, 0.88, 0.24)
    )

  ggplot() +
    geom_segment(
      data = spokes,
      aes(x, y, xend = xend, yend = yend, colour = Direction),
      linewidth = spokes$EdgeWidth, alpha = spokes$EdgeAlpha, lineend = "round"
    ) +
    geom_segment(
      data = internal_plot,
      aes(x, y, xend = xend, yend = yend, colour = Direction),
      linewidth = internal_plot$EdgeWidth, alpha = internal_plot$EdgeAlpha,
      lineend = "round"
    ) +
    geom_point(
      data = nodes,
      aes(x, y, fill = Tissue, size = NodeSize, alpha = NodeAlpha),
      shape = 21, colour = "white", stroke = 0.45
    ) +
    geom_point(
      data = nodes %>% filter(FDRNode, Tissue != "Clinical"),
      aes(x, y, size = NodeSize), shape = 21, fill = NA,
      colour = "#30353A", stroke = 0.60
    ) +
    geom_text_repel(
      data = nodes %>% filter(NodeLabel != "", Tissue != "Clinical"),
      aes(x, y, label = NodeLabel), family = "Arial", fontface = "bold",
      size = 1.78, min.segment.length = 0, segment.size = 0.18,
      segment.color = "#777777", box.padding = 0.15, point.padding = 0.08,
      max.overlaps = Inf, seed = 20260829L
    ) +
    geom_text(
      data = nodes %>% filter(Tissue == "Clinical"),
      aes(x, y, label = NodeLabel), family = "Arial", fontface = "bold",
      size = 1.95, nudge_y = 0.13
    ) +
    scale_colour_manual(values = edge_colors) +
    scale_fill_manual(values = tissue_colors) +
    scale_size_identity() +
    scale_alpha_identity() +
    coord_equal(xlim = c(-0.94, 0.94), ylim = c(-1.02, 1.31), clip = "off") +
    labs(title = unname(state_labels[state])) +
    theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = 6.6, hjust = 0.5),
      plot.margin = margin(1, 1, 1, 1), legend.position = "none"
    )
}

bar_data <- node_table %>%
  filter(FDRNode) %>%
  mutate(
    StateLabel = unname(state_labels[State]),
    StateLabel = factor(StateLabel, levels = unname(state_labels[state_levels])),
    StateAxis = factor(if_else(State == "PreExercise", "Pre", "Response"), levels = c("Pre", "Response")),
    MoleculeLabel = if_else(DisplayLabel == "WASH3P/WASH2P", "WASH3P/\nWASH2P", DisplayLabel),
    Stars = case_when(
      BH_FDR < 0.05 ~ "***",
      BH_FDR < 0.10 ~ "**",
      BH_FDR < 0.20 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(StateAxis, desc(Effect)) %>%
  group_by(StateAxis) %>%
  mutate(MoleculeLabel = factor(MoleculeLabel, levels = unique(MoleculeLabel))) %>%
  ungroup()

bar_limit <- 2
bar_plot <- ggplot(bar_data, aes(MoleculeLabel, Effect, fill = StateLabel)) +
  geom_hline(yintercept = 0, linewidth = 0.30, colour = "#444444") +
  geom_col(width = 0.60) +
  geom_text(
    aes(
      y = if_else(Effect >= 0, Effect + 0.08, Effect - 0.08),
      label = Stars, vjust = if_else(Effect >= 0, 0, 1)
    ),
    family = "Arial", fontface = "bold", size = 2.20
  ) +
  scale_fill_manual(values = state_colors) +
  scale_y_continuous(
    limits = c(-2, 2), breaks = c(-2, 0, 2),
    expand = expansion(mult = 0)
  ) +
  facet_grid(
    cols = vars(StateAxis), scales = "free_x", space = "free_x"
  ) +
  labs(
    title = "hsCRP-associated proteins",
    subtitle = NULL,
    x = NULL, y = expression("Effect (" * beta * ")"), fill = NULL
  ) +
  theme_classic(base_family = "Arial", base_size = 6.3) +
  theme(
    plot.title = element_text(face = "bold", size = 5.8, hjust = 0.5),
    axis.title.y = element_text(face = "bold", size = 5.6),
    axis.text.x = element_text(size = 4.2, colour = "black", angle = 60, hjust = 1),
    axis.text.y = element_text(size = 5.0, colour = "black"),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.35),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 4.9, face = "bold", colour = "black"),
    legend.position = "none",
    legend.text = element_text(size = 5.0),
    legend.key.width = unit(4.5, "mm"),
    plot.margin = margin(7, 1, 4, 1)
  )

bar_block <- wrap_plots(
  plot_spacer(), bar_plot, plot_spacer(),
  ncol = 1, heights = c(0.14, 0.72, 0.14)
)

legend_plot <- ggplot() +
  annotate("point", x = 0.30, y = c(4.0, 3.2),
           shape = 21, size = 2.25, stroke = 0.38,
           fill = unname(tissue_colors[c("Adipose", "Clinical")]), colour = "white") +
  annotate("text", x = 0.60, y = c(4.0, 3.2),
           label = c("Adipose protein", "hsCRP"), hjust = 0,
           family = "Arial", size = 1.55) +
  annotate("segment", x = 0.10, xend = 0.50,
           y = c(2.4, 1.6), yend = c(2.4, 1.6),
           colour = unname(edge_colors[c("Positive", "Negative")]),
           linewidth = 0.62, lineend = "round") +
  annotate("text", x = 0.60, y = c(2.4, 1.6),
           label = c("Positive", "Negative"), hjust = 0,
           family = "Arial", size = 1.55) +
  coord_cartesian(xlim = c(0, 2.9), ylim = c(1.2, 4.4), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(2, 0, 2, 1, unit = "mm"))

main_block <- wrap_plots(
  make_network_panel("PreExercise"),
  make_network_panel("ExerciseChange"),
  bar_block,
  nrow = 1, widths = c(1.22, 1.22, 0.60)
)

final_plot <- wrap_plots(main_block, legend_plot, nrow = 1, widths = c(8.45, 1.55)) +
  plot_annotation(
    title = "hsCRP-associated adipose-proteomic networks",
    theme = theme(
      text = element_text(family = "Arial", colour = "black"),
      plot.title = element_text(face = "bold", size = 8.4, hjust = 0, margin = margin(b = 0)),
      plot.margin = margin(1, 1, 0, 1)
    )
  )

pdf_file <- file.path(output_dir, "Fig05B_hsCRP_Adipose_Association_Network_HalfWidth_Candidate.pdf")
png_file <- file.path(output_dir, "Fig05B_hsCRP_Adipose_Association_Network_HalfWidth_Candidate.png")
ggsave(pdf_file, final_plot, width = 122, height = 72, units = "mm", device = cairo_pdf)
ggsave(png_file, final_plot, width = 122, height = 72, units = "mm", dpi = 450, bg = "white")

summary_table <- candidate_results %>%
  group_by(State, Dataset) %>%
  summarise(
    Displayed_candidates = n(), NominalP01 = sum(P_Value < 0.01),
    BH_FDR05 = sum(BH_FDR < 0.05), BH_FDR10 = sum(BH_FDR < 0.10),
    BH_FDR20 = sum(BH_FDR < 0.20), .groups = "drop"
  )

workbook_file <- file.path(output_dir, "Fig05B_hsCRP_Adipose_Association_Network_HalfWidth_SourceData.xlsx")
write.xlsx(
  list(
    Displayed_nodes = node_table %>% arrange(State, Dataset, BH_FDR, P_Value),
    Displayed_internal_edges = display_correlations %>% arrange(State, Evidence, BH_FDR, P_Value),
    All_tested_internal_edges = all_correlations %>% arrange(State, BH_FDR, P_Value),
    FDR_bar_data = bar_data %>% mutate(MoleculeLabel = as.character(MoleculeLabel)),
    Display_summary = summary_table,
    Methods = tibble(
      Item = c(
        "Direct association model", "Context-node threshold", "FDR-node threshold",
        "Internal-edge model", "Internal-edge FDR family", "Stability",
        "Interpretation"
      ),
      Value = c(
        "Protein state ~ baseline hsCRP + Role; FamilyID duplicateCorrelation",
        "Nominal P < 0.01; bounded to the top 16 proteins per tissue and state",
        "BH-FDR < 0.20 within dataset, state and predictor",
        "Spearman correlation of Role-adjusted protein values within tissue and state",
        "BH correction across all selected protein pairs within each state",
        "Leave-one-family-out sign stability >= 0.80 for displayed internal edges",
        "Exploratory association network; layout distance and edges do not imply causality"
      )
    )
  ),
  workbook_file, overwrite = TRUE
)

run_log <- c(
  paste0("Analysis ID: ", analysis_id),
  paste0("Formal input: ", normalizePath(formal_file, winslash = "/")),
  paste0("Output PDF: ", normalizePath(pdf_file, winslash = "/")),
  paste0("Output workbook: ", normalizePath(workbook_file, winslash = "/")),
  "Figure 5B candidate only; no editable artwork or assembled figure was modified.",
  capture.output(sessionInfo())
)
writeLines(run_log, file.path(output_dir, "RUN_LOG.txt"))

message("Figure 5B candidate written to: ", output_dir)
