# Revised Figure 3C: adiposity-stratified multi-omic exercise-response networks.
#
# Analysis-ready inputs are read only. The script writes a new candidate folder
# and does not replace assembled Figure 3 artwork or Adobe Illustrator files.

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
  dplyr, tidyr, tibble, stringr, readr, ggplot2, ggraph, igraph,
  ggrepel, patchwork, scales, openxlsx, graphlayouts
)

set.seed(20260826L)

analysis_id <- "Fig03C_Adiposity_Response_Network_Revised_2026-08-26"
output_dir <- path_project("fig03", "Figure3C_Adiposity_Response_Network_Revised_2026-08-26")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

minimum_pair_n <- 6L
nominal_p_threshold <- 0.05
absolute_rho_threshold <- 0.50
display_edge_cap <- 180L

clinical_spec <- tibble::tribble(
  ~Variable, ~DisplayLabel, ~Order,
  "Fat_p_G", "Body-fat percentage", 1L,
  "Leptin_84m_G", "Leptin", 2L,
  "TRIGLY_G", "Triglycerides", 3L,
  "RetroperitonealFatKg_G", "Retroperitoneal fat", 4L,
  "IntraperitonealFatKg_G", "Visceral fat", 5L,
  "Liverfat_G", "Liver fat", 6L
)

dataset_spec <- tibble::tribble(
  ~Dataset, ~MatrixFile, ~Tissue, ~Modality, ~FeatureCap,
  "Serum_Metabonomics", "Analysis_Serum_Metabonomics.csv", "Serum", "Metabolite", 20L,
  "Serum_Proteomics", "Analysis_Serum_Proteomics.csv", "Serum", "Protein", 20L,
  "Adipose_Proteomics", "Analysis_Adipose_Proteomics.csv", "Adipose", "Protein", 15L,
  "Muscle_Proteomics", "Analysis_Muscle_Proteomics.csv", "Muscle", "Protein", 15L
)

clinical_file <- path_analysis_ready("Metadata", "Clinical_Observed_Primary.csv")
sample_sheet_file <- path_analysis_ready("Metadata", "Project_Sample_Sheet.csv")
metabolite_dictionary_file <- path_analysis_ready(
  "Metadata", "Feature_Dictionary_Serum_Metabonomics.csv"
)
serum_protein_dictionary_file <- path_analysis_ready(
  "Metadata", "Feature_Dictionary_Serum_Proteomics.csv"
)
adipose_protein_annotation_file <- path_analysis_ready(
  "Metadata", "Feature_Annotation_Adipose_Proteomics.csv"
)
muscle_protein_annotation_file <- path_analysis_ready(
  "Metadata", "Feature_Annotation_Muscle_Proteomics.csv"
)

required_files <- c(
  clinical_file,
  sample_sheet_file,
  metabolite_dictionary_file,
  serum_protein_dictionary_file,
  adipose_protein_annotation_file,
  muscle_protein_annotation_file,
  file.path(DIR_ANALYSIS_READY, dataset_spec$MatrixFile)
)
if (any(!file.exists(required_files))) {
  stop("Missing input(s): ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}

clinical_raw <- read_csv(clinical_file, show_col_types = FALSE)
missing_clinical <- setdiff(c("FamilyID", "Clinical_Subject_ID", "Role", clinical_spec$Variable), names(clinical_raw))
if (length(missing_clinical)) {
  stop("Clinical input is missing: ", paste(missing_clinical, collapse = ", "))
}

clinical <- clinical_raw %>%
  transmute(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = as.character(Role),
    across(all_of(clinical_spec$Variable), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  mutate(
    AdiposityGroup = classify_body_fat(Fat_p_G),
    AdiposityGroup = factor(AdiposityGroup, levels = c("Lean", "Obese"))
  ) %>%
  distinct(Clinical_Subject_ID, .keep_all = TRUE)

sample_sheet <- read_csv(sample_sheet_file, show_col_types = FALSE) %>%
  mutate(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    SampleID = tolower(as.character(SampleID)),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort),
    RetainInPrimaryInput = as.logical(RetainInPrimaryInput)
  ) %>%
  filter(
    AnalysisSet == "Primary", Is_Primary_Cohort,
    is.na(RetainInPrimaryInput) | RetainInPrimaryInput
  ) %>%
  inner_join(
    clinical %>% select(FamilyID, Clinical_Subject_ID, Role, AdiposityGroup),
    by = c("FamilyID", "Clinical_Subject_ID", "Role")
  )

metabolite_dictionary <- read_csv(metabolite_dictionary_file, show_col_types = FALSE) %>%
  transmute(
    FeatureID = as.character(FeatureID),
    LongName = na_if(str_squish(as.character(LongName)), "")
  ) %>%
  distinct(FeatureID, .keep_all = TRUE)

protein_dictionary <- bind_rows(
  read_csv(serum_protein_dictionary_file, show_col_types = FALSE) %>%
    transmute(
      Dataset = "Serum_Proteomics",
      FeatureID = as.character(FeatureID),
      GeneSymbol = na_if(str_squish(as.character(GeneSymbol)), "")
    ),
  read_csv(adipose_protein_annotation_file, show_col_types = FALSE) %>%
    transmute(
      Dataset = "Adipose_Proteomics",
      FeatureID = as.character(FeatureID),
      GeneSymbol = coalesce(
        na_if(str_squish(as.character(IntegrationGene)), ""),
        na_if(str_squish(as.character(PG.Genes)), "")
      )
    ),
  read_csv(muscle_protein_annotation_file, show_col_types = FALSE) %>%
    transmute(
      Dataset = "Muscle_Proteomics",
      FeatureID = as.character(FeatureID),
      GeneSymbol = coalesce(
        na_if(str_squish(as.character(IntegrationGene)), ""),
        na_if(str_squish(as.character(PG.Genes)), "")
      )
    )
) %>%
  mutate(GeneSymbol = str_replace(GeneSymbol, ";.*$", "")) %>%
  distinct(Dataset, FeatureID, .keep_all = TRUE)

make_display_label <- function(feature_id, dataset) {
  if (dataset == "Serum_Metabonomics") {
    long_name <- metabolite_dictionary$LongName[match(feature_id, metabolite_dictionary$FeatureID)]
    label <- ifelse(!is.na(long_name) & nchar(long_name) <= 28L, long_name, feature_id)
  } else {
    gene_symbol <- protein_dictionary$GeneSymbol[
      match(paste(dataset, feature_id), paste(protein_dictionary$Dataset, protein_dictionary$FeatureID))
    ]
    label <- ifelse(!is.na(gene_symbol) & nzchar(gene_symbol), gene_symbol, feature_id)
  }
  label <- str_replace(label, "\\.\\.\\..*$", "")
  label <- str_replace_all(label, "_", " ")
  str_squish(label)
}

prepare_response <- function(dataset, matrix_file) {
  matrix <- read_analysis_matrix(path_analysis_ready(matrix_file))
  metadata <- sample_sheet %>%
    filter(
      Dataset == dataset,
      Timepoint %in% c("Pre", "Post3h"),
      SampleID %in% colnames(matrix)
    ) %>%
    distinct(Clinical_Subject_ID, Timepoint, .keep_all = TRUE) %>%
    select(FamilyID, Clinical_Subject_ID, Role, AdiposityGroup, Timepoint, SampleID) %>%
    pivot_wider(names_from = Timepoint, values_from = SampleID) %>%
    filter(!is.na(Pre), !is.na(Post3h)) %>%
    arrange(FamilyID, Role, Clinical_Subject_ID)

  response <- matrix[, metadata$Post3h, drop = FALSE] - matrix[, metadata$Pre, drop = FALSE]
  colnames(response) <- metadata$Clinical_Subject_ID
  storage.mode(response) <- "double"
  list(response = response, metadata = metadata)
}

response_objects <- lapply(seq_len(nrow(dataset_spec)), function(i) {
  prepare_response(dataset_spec$Dataset[i], dataset_spec$MatrixFile[i])
})
names(response_objects) <- dataset_spec$Dataset

# Shared feature selection prevents group-specific node selection from being
# mistaken for a biological network difference.
selected_features <- lapply(seq_len(nrow(dataset_spec)), function(i) {
  spec <- dataset_spec[i, ]
  response <- response_objects[[spec$Dataset]]$response
  feature_sd <- apply(response, 1L, stats::sd, na.rm = TRUE)
  feature_sd[!is.finite(feature_sd)] <- -Inf
  tibble(
    Dataset = spec$Dataset,
    FeatureID = names(sort(feature_sd, decreasing = TRUE))[seq_len(min(spec$FeatureCap, nrow(response)))],
    FeatureSD = sort(feature_sd, decreasing = TRUE)[seq_len(min(spec$FeatureCap, nrow(response)))],
    Tissue = spec$Tissue,
    Modality = spec$Modality
  )
}) %>% bind_rows() %>%
  mutate(
    NodeID = paste(Dataset, FeatureID, sep = "::"),
    DisplayLabel = mapply(make_display_label, FeatureID, Dataset, USE.NAMES = FALSE)
  )

molecular_values <- lapply(seq_len(nrow(selected_features)), function(i) {
  row <- selected_features[i, ]
  response <- response_objects[[row$Dataset]]$response
  tibble(
    NodeID = row$NodeID,
    Clinical_Subject_ID = colnames(response),
    Value = as.numeric(response[row$FeatureID, ])
  )
}) %>% bind_rows()

clinical_nodes <- clinical_spec %>%
  transmute(
    Dataset = "Clinical",
    FeatureID = Variable,
    FeatureSD = NA_real_,
    Tissue = "Clinical",
    Modality = "Clinical",
    NodeID = paste("Clinical", Variable, sep = "::"),
    DisplayLabel,
    CoreNode = TRUE,
    ClinicalOrder = Order
  )

clinical_values <- lapply(clinical_spec$Variable, function(variable) {
  tibble(
    NodeID = paste("Clinical", variable, sep = "::"),
    Clinical_Subject_ID = clinical$Clinical_Subject_ID,
    Value = clinical[[variable]]
  )
}) %>% bind_rows()

node_table <- bind_rows(
  selected_features %>%
    mutate(CoreNode = FALSE, ClinicalOrder = NA_integer_) %>%
    select(Dataset, FeatureID, FeatureSD, Tissue, Modality, NodeID, DisplayLabel, CoreNode, ClinicalOrder),
  clinical_nodes
) %>%
  mutate(
    Tissue = factor(Tissue, levels = c("Clinical", "Serum", "Adipose", "Muscle")),
    Modality = factor(Modality, levels = c("Clinical", "Metabolite", "Protein"))
  )

value_table <- bind_rows(molecular_values, clinical_values) %>%
  left_join(
    clinical %>% select(Clinical_Subject_ID, FamilyID, Role, AdiposityGroup),
    by = "Clinical_Subject_ID"
  ) %>%
  filter(!is.na(AdiposityGroup))

node_vectors <- split(value_table, value_table$NodeID)
node_ids <- node_table$NodeID
node_pairs <- t(combn(node_ids, 2L)) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  setNames(c("NodeA", "NodeB"))

spearman_test <- function(x, y) {
  observed <- is.finite(x) & is.finite(y)
  x <- x[observed]
  y <- y[observed]
  n <- length(x)
  if (n < minimum_pair_n || sd(x) == 0 || sd(y) == 0) {
    return(c(N = n, Rho = NA_real_, PValue = NA_real_))
  }
  has_ties <- anyDuplicated(x) > 0L || anyDuplicated(y) > 0L
  test <- suppressWarnings(stats::cor.test(
    x, y, method = "spearman", exact = n <= 9L && !has_ties
  ))
  p_value <- unname(test$p.value)
  if (!is.finite(p_value) || p_value <= 0) {
    # Conservative finite lower bound for small no-tie samples.
    p_value <- if (n <= 10L) min(1, 2 / factorial(n)) else .Machine$double.xmin
  }
  c(N = n, Rho = unname(test$estimate), PValue = p_value)
}

test_group_network <- function(group_name) {
  group_subjects <- clinical %>%
    filter(AdiposityGroup == group_name) %>%
    pull(Clinical_Subject_ID)

  tests <- lapply(seq_len(nrow(node_pairs)), function(i) {
    node_a <- node_pairs$NodeA[i]
    node_b <- node_pairs$NodeB[i]
    values_a <- node_vectors[[node_a]] %>%
      filter(Clinical_Subject_ID %in% group_subjects) %>%
      select(Clinical_Subject_ID, ValueA = Value)
    values_b <- node_vectors[[node_b]] %>%
      filter(Clinical_Subject_ID %in% group_subjects) %>%
      select(Clinical_Subject_ID, ValueB = Value)
    merged <- inner_join(values_a, values_b, by = "Clinical_Subject_ID")
    result <- spearman_test(merged$ValueA, merged$ValueB)
    tibble(
      Group = group_name,
      NodeA = node_a,
      NodeB = node_b,
      N = as.integer(result[["N"]]),
      Rho = result[["Rho"]],
      PValue = result[["PValue"]]
    )
  }) %>% bind_rows() %>%
    mutate(
      BH_FDR = p.adjust(PValue, method = "BH"),
      Direction = if_else(Rho >= 0, "Positive", "Negative"),
      NominalEdge = is.finite(PValue) & PValue < nominal_p_threshold &
        is.finite(Rho) & abs(Rho) >= absolute_rho_threshold,
      FDR05 = is.finite(BH_FDR) & BH_FDR < 0.05,
      FDR10 = is.finite(BH_FDR) & BH_FDR < 0.10,
      FDR20 = is.finite(BH_FDR) & BH_FDR < 0.20
    ) %>%
    left_join(node_table %>% select(NodeA = NodeID, NodeA_Tissue = Tissue, NodeA_Label = DisplayLabel), by = "NodeA") %>%
    left_join(node_table %>% select(NodeB = NodeID, NodeB_Tissue = Tissue, NodeB_Label = DisplayLabel), by = "NodeB") %>%
    mutate(
      NodeA_Tissue = as.character(NodeA_Tissue),
      NodeB_Tissue = as.character(NodeB_Tissue),
      TissuePair = if_else(
        NodeA_Tissue <= NodeB_Tissue,
        paste(NodeA_Tissue, NodeB_Tissue, sep = " - "),
        paste(NodeB_Tissue, NodeA_Tissue, sep = " - ")
      ),
      ClinicalClinical = NodeA_Tissue == "Clinical" & NodeB_Tissue == "Clinical",
      CoreRelated = xor(NodeA_Tissue == "Clinical", NodeB_Tissue == "Clinical"),
      DisplayEligible = NominalEdge & !ClinicalClinical
    )
  tests
}

all_tests <- bind_rows(test_group_network("Lean"), test_group_network("Obese"))

select_display_edges <- function(tests) {
  nominal <- tests %>% filter(DisplayEligible)
  if (!nrow(nominal)) return(nominal)

  # Retain all FDR-supported and clinical-to-molecular core edges, then add strongest
  # nominal edges until the visual cap is reached.
  required <- nominal %>%
    filter(FDR20 | CoreRelated) %>%
    arrange(BH_FDR, PValue, desc(abs(Rho)))
  remaining <- nominal %>%
    anti_join(required %>% select(NodeA, NodeB), by = c("NodeA", "NodeB")) %>%
    arrange(PValue, desc(abs(Rho)), desc(N))
  bind_rows(required, remaining) %>%
    distinct(NodeA, NodeB, .keep_all = TRUE) %>%
    slice_head(n = display_edge_cap)
}

display_edges <- all_tests %>%
  group_split(Group) %>%
  lapply(select_display_edges) %>%
  bind_rows()

select_fdr05_display_edges <- function(tests, within_source_global_cap = 12L) {
  strict_edges <- tests %>%
    filter(
      !ClinicalClinical,
      FDR05,
      is.finite(Rho),
      abs(Rho) >= absolute_rho_threshold
    ) %>%
    arrange(BH_FDR, PValue, desc(abs(Rho)))
  if (!nrow(strict_edges)) return(strict_edges)

  cross_source <- strict_edges %>%
    filter(NodeA_Tissue != NodeB_Tissue)
  within_source <- strict_edges %>%
    filter(NodeA_Tissue == NodeB_Tissue)

  # Preserve each connected node's strongest within-source relationship, then
  # add a small global set of the strongest within-source edges. This affects
  # visualization only; every tested and FDR-supported edge remains in the
  # source-data workbook.
  strongest_by_node_a <- within_source %>%
    group_by(NodeA) %>%
    slice_min(order_by = BH_FDR, n = 1, with_ties = FALSE) %>%
    ungroup()
  strongest_by_node_b <- within_source %>%
    group_by(NodeB) %>%
    slice_min(order_by = BH_FDR, n = 1, with_ties = FALSE) %>%
    ungroup()
  strongest_global <- within_source %>%
    slice_head(n = within_source_global_cap)

  bind_rows(
    cross_source,
    strongest_by_node_a,
    strongest_by_node_b,
    strongest_global
  ) %>%
    distinct(NodeA, NodeB, .keep_all = TRUE) %>%
    mutate(
      DisplayBasis = case_when(
        CoreRelated ~ "Clinical-to-molecular",
        NodeA_Tissue != NodeB_Tissue ~ "Cross-source",
        TRUE ~ "Within-source anchor"
      )
    ) %>%
    arrange(BH_FDR, PValue, desc(abs(Rho)))
}

fdr05_display_edges <- all_tests %>%
  group_split(Group) %>%
  lapply(select_fdr05_display_edges) %>%
  bind_rows()

retain_explanatory_components <- function(edges) {
  if (!nrow(edges)) return(edges)
  component_graph <- graph_from_data_frame(
    edges %>% select(from = NodeA, to = NodeB),
    directed = FALSE
  )
  component_info <- components(component_graph)
  explanatory_nodes <- edges %>%
    filter(CoreRelated | NodeA_Tissue != NodeB_Tissue) %>%
    select(NodeA, NodeB) %>%
    unlist(use.names = FALSE) %>%
    unique()
  retained_components <- union(
    which(component_info$csize >= 3),
    unname(component_info$membership[intersect(explanatory_nodes, names(component_info$membership))])
  )
  retained_nodes <- names(component_info$membership)[
    component_info$membership %in% retained_components
  ]
  edges %>%
    filter(NodeA %in% retained_nodes, NodeB %in% retained_nodes)
}

fdr05_plot_edges <- fdr05_display_edges %>%
  group_split(Group) %>%
  lapply(retain_explanatory_components) %>%
  bind_rows()

tissue_colors <- c(
  "Adipose" = "#E98773",
  "Clinical" = "#79B8A9",
  "Muscle" = "#72A8CC",
  "Serum" = "#E9B85F"
)
edge_colors <- c("Positive" = "#D96B63", "Negative" = "#5D91B8")
modality_shapes <- c("Clinical" = 22, "Metabolite" = 23, "Protein" = 21)

make_panel <- function(group_name, evidence = c("Nominal", "FDR20", "FDR05Optimized")) {
  evidence <- match.arg(evidence)
  group_edges <- if (evidence == "Nominal") {
    display_edges %>% filter(Group == group_name)
  } else if (evidence == "FDR05Optimized") {
    fdr05_plot_edges %>% filter(Group == group_name)
  } else {
    all_tests %>%
      filter(Group == group_name, !ClinicalClinical, FDR20) %>%
      arrange(BH_FDR, PValue, desc(abs(Rho)))
  }
  if (!nrow(group_edges)) stop("No edges available for ", group_name, " / ", evidence)

  edge_frame <- group_edges %>%
    transmute(
      from = NodeA, to = NodeB, N, Rho, PValue, BH_FDR, Direction,
      FDRClass = case_when(
        FDR05 ~ "BH-FDR < 0.05",
        FDR10 ~ "BH-FDR < 0.10",
        FDR20 ~ "BH-FDR < 0.20",
        TRUE ~ "Nominal P < 0.05"
      ),
      EdgeWidth = rescale(abs(Rho), to = c(0.22, 0.95))
    )
  connected_nodes <- union(edge_frame$from, edge_frame$to)
  graph <- graph_from_data_frame(
    edge_frame,
    directed = FALSE,
    vertices = node_table %>%
      filter(NodeID %in% connected_nodes) %>%
      transmute(
        name = NodeID, Dataset, FeatureID, Tissue, Modality,
        DisplayLabel, CoreNode, ClinicalOrder
      )
  )
  message(
    "Panel ", group_name, " / ", evidence,
    ": ", gsize(graph), " edges; edge attributes = ",
    paste(edge_attr_names(graph), collapse = ", ")
  )
  V(graph)$Degree <- degree(graph)
  evidence_seed <- case_when(
    evidence == "FDR20" ~ 100L,
    evidence == "FDR05Optimized" ~ 200L,
    TRUE ~ 0L
  )
  key_edge_nodes <- group_edges %>%
    filter(CoreRelated) %>%
    transmute(NodeA, NodeB) %>%
    unlist(use.names = FALSE) %>%
    unique()
  set.seed(if_else(group_name == "Lean", 202608261L, 202608262L) + evidence_seed)
  layout <- if (evidence == "FDR05Optimized") {
    create_layout(graph, layout = "fr", niter = 2500, grid = "nogrid")
  } else {
    create_layout(graph, layout = "stress")
  }
  layout <- layout %>%
    group_by(Tissue) %>%
    mutate(TissueDegreeRank = min_rank(desc(Degree))) %>%
    ungroup() %>%
    mutate(
      LabelNode = case_when(
        CoreNode ~ Degree >= 2,
        name %in% key_edge_nodes ~ TRUE,
        Tissue == "Serum" ~ TissueDegreeRank <= 4,
        Tissue %in% c("Adipose", "Muscle") ~ TissueDegreeRank <= 3,
        TRUE ~ FALSE
      ),
      NodeSize = if_else(CoreNode, 4.8, rescale(Degree, to = c(1.7, 4.0))),
      NodeLabel = if_else(LabelNode, str_wrap(DisplayLabel, 15), "")
    )
  node_plot_data <- as_tibble(layout)
  node_coordinates <- node_plot_data %>%
    select(name, x, y)
  edge_plot_data <- edge_frame %>%
    left_join(
      node_coordinates %>% rename(from = name, x = x, y = y),
      by = "from"
    ) %>%
    left_join(
      node_coordinates %>% rename(to = name, xend = x, yend = y),
      by = "to"
    )

  ggplot(node_plot_data, aes(x = x, y = y)) +
    geom_segment(
      data = edge_plot_data,
      aes(x = x, y = y, xend = xend, yend = yend,
          color = Direction, linewidth = EdgeWidth),
      alpha = if_else(evidence == "FDR05Optimized", 0.52, 0.68),
      lineend = "round", show.legend = TRUE,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = edge_colors, name = "Association") +
    scale_linewidth_identity() +
    geom_point(
      aes(fill = Tissue, shape = Modality, size = NodeSize),
      color = "white", stroke = 0.48
    ) +
    scale_fill_manual(values = tissue_colors, name = "Tissue / source") +
    scale_shape_manual(values = modality_shapes, name = "Molecular layer") +
    scale_size_identity() +
    geom_text_repel(
      aes(label = NodeLabel),
      family = "Arial", size = 1.95, fontface = "bold",
      force = 4, max.time = 5, box.padding = 0.34, point.padding = 0.15,
      max.overlaps = Inf, min.segment.length = 0,
      segment.size = 0.16, segment.color = "#888888"
    ) +
    labs(
      title = group_name,
      subtitle = if (evidence == "Nominal") {
        paste0(nrow(group_edges), " displayed nominal edges; ",
               sum(group_edges$FDR05), " at BH-FDR < 0.05")
      } else if (evidence == "FDR05Optimized") {
        NULL
      } else {
        paste0(nrow(group_edges), " edges at BH-FDR < 0.20")
      }
    ) +
    theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
      plot.subtitle = element_text(size = 5.7, color = "#555555", hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2)
    )
}

lean_plot <- make_panel("Lean", "Nominal")
obese_plot <- make_panel("Obese", "Nominal")
lean_fdr20_plot <- make_panel("Lean", "FDR20")
obese_fdr20_plot <- make_panel("Obese", "FDR20")
lean_fdr05_optimized_plot <- make_panel("Lean", "FDR05Optimized")
obese_fdr05_optimized_plot <- make_panel("Obese", "FDR05Optimized")

combined_plot <- lean_plot + obese_plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Adiposity-stratified multi-omic exercise-response networks",
    subtitle = paste0(
      "Post3h - Pre molecular responses; ",
      "displayed edges require abs(Spearman rho) >= ", absolute_rho_threshold,
      " and nominal P < ", nominal_p_threshold
    )
  ) &
  theme(
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(face = "bold", size = 9.2, hjust = 0),
    plot.subtitle = element_text(size = 5.8, color = "#555555", hjust = 0)
  )

fdr20_combined_plot <- lean_fdr20_plot + obese_fdr20_plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "BH-FDR-supported multi-omic exercise-response networks",
    subtitle = "Post3h - Pre molecular responses; edges require BH-FDR < 0.20 within each adiposity group"
  ) &
  theme(
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(face = "bold", size = 9.2, hjust = 0),
    plot.subtitle = element_text(size = 5.8, color = "#555555", hjust = 0)
  )

fdr05_optimized_plot <- lean_fdr05_optimized_plot + obese_fdr05_optimized_plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Adiposity-stratified molecular response networks",
    subtitle = paste0(
      "Post3h - Pre responses; displayed edges require BH-FDR < 0.05 and absolute Spearman rho >= ",
      absolute_rho_threshold
    )
  ) &
  theme(
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(face = "bold", size = 9.2, hjust = 0),
    plot.subtitle = element_text(size = 5.8, color = "#555555", hjust = 0)
  )

legend_plot <- ggplot() +
  annotate("point", x = c(0.08, 0.28, 0.48, 0.68), y = 0.73,
           shape = 21, size = 4.0, stroke = 0.45, color = "white",
           fill = unname(tissue_colors[c("Clinical", "Serum", "Adipose", "Muscle")])) +
  annotate("text", x = c(0.12, 0.32, 0.52, 0.72), y = 0.73,
           label = c("Clinical", "Serum", "Adipose", "Muscle"),
           hjust = 0, family = "Arial", size = 2.7) +
  annotate("point", x = c(0.08, 0.28, 0.48), y = 0.43,
           shape = unname(modality_shapes[c("Clinical", "Metabolite", "Protein")]),
           size = 3.7, stroke = 0.45, color = "#555555", fill = "white") +
  annotate("text", x = c(0.12, 0.32, 0.52), y = 0.43,
           label = c("Clinical", "Metabolite", "Protein"),
           hjust = 0, family = "Arial", size = 2.7) +
  annotate("segment", x = c(0.08, 0.28), xend = c(0.15, 0.35),
           y = 0.15, yend = 0.15,
           color = unname(edge_colors[c("Positive", "Negative")]),
           linewidth = 0.9) +
  annotate("text", x = c(0.17, 0.37), y = 0.15,
           label = c("Positive", "Negative"), hjust = 0,
           family = "Arial", size = 2.7) +
  annotate("text", x = 0.57, y = 0.15,
           label = "Edge width: |rho|; all edges BH-FDR < 0.05", hjust = 0,
           family = "Arial", size = 2.7) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(1, 1, 1, 1))

ggsave(
  file.path(output_dir, "Fig03C_Adiposity_Response_Networks_Revised.pdf"),
  combined_plot, width = 190, height = 91, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Fig03C_Adiposity_Response_Networks_BH_FDR20.pdf"),
  fdr20_combined_plot, width = 190, height = 91, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Fig03C_Adiposity_Response_Networks_BH_FDR05_Optimized.pdf"),
  fdr05_optimized_plot, width = 190, height = 91, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Fig03C_Adiposity_Response_Networks_BH_FDR05_Optimized.png"),
  fdr05_optimized_plot, width = 190, height = 91, units = "mm", dpi = 400, bg = "white"
)
ggsave(
  file.path(output_dir, "Fig03C_Adiposity_Response_Networks_BH_FDR20.png"),
  fdr20_combined_plot, width = 190, height = 91, units = "mm", dpi = 400, bg = "white"
)
ggsave(
  file.path(output_dir, "Fig03C_Adiposity_Response_Networks_Revised.png"),
  combined_plot, width = 190, height = 91, units = "mm", dpi = 400, bg = "white"
)
ggsave(
  file.path(output_dir, "Fig03C_Lean_Response_Network_BH_FDR20.pdf"),
  lean_fdr20_plot, width = 92, height = 86, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Fig03C_Obese_Response_Network_BH_FDR20.pdf"),
  obese_fdr20_plot, width = 92, height = 86, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Fig03C_Lean_Response_Network_BH_FDR05_Optimized.pdf"),
  lean_fdr05_optimized_plot, width = 92, height = 86, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Fig03C_Obese_Response_Network_BH_FDR05_Optimized.pdf"),
  obese_fdr05_optimized_plot, width = 92, height = 86, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Fig03C_Lean_Response_Network_Revised.pdf"),
  lean_plot, width = 92, height = 86, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Fig03C_Obese_Response_Network_Revised.pdf"),
  obese_plot, width = 92, height = 86, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Fig03C_Response_Network_Legend.pdf"),
  legend_plot, width = 118, height = 27, units = "mm", device = cairo_pdf
)

summary_table <- all_tests %>%
  group_by(Group) %>%
  summarise(
    Participants = n_distinct(clinical$Clinical_Subject_ID[clinical$AdiposityGroup == first(Group)]),
    CandidateNodes = nrow(node_table),
    TestedPairs = sum(is.finite(PValue)),
    NominalP05 = sum(is.finite(PValue) & PValue < 0.05),
    NominalP05_Rho05 = sum(NominalEdge),
    BH_FDR05 = sum(FDR05),
    BH_FDR10 = sum(FDR10),
    BH_FDR20 = sum(FDR20),
    Molecular_BH_FDR05 = sum(FDR05 & !ClinicalClinical),
    Molecular_BH_FDR10 = sum(FDR10 & !ClinicalClinical),
    Molecular_BH_FDR20 = sum(FDR20 & !ClinicalClinical),
    CoreRelated_Nominal = sum(NominalEdge & CoreRelated),
    CoreRelated_FDR05 = sum(FDR05 & CoreRelated),
    DisplayedEdges = sum(display_edges$Group == first(Group)),
    DisplayedFDR05OptimizedEdges = sum(fdr05_plot_edges$Group == first(Group)),
    .groups = "drop"
  )

sample_availability <- bind_rows(lapply(seq_len(nrow(dataset_spec)), function(i) {
  spec <- dataset_spec[i, ]
  metadata <- response_objects[[spec$Dataset]]$metadata
  metadata %>%
    count(AdiposityGroup, name = "PairedN") %>%
    mutate(Dataset = spec$Dataset, Tissue = spec$Tissue, Modality = spec$Modality)
})) %>%
  bind_rows(
    lapply(clinical_spec$Variable, function(variable) {
      clinical %>%
        group_by(AdiposityGroup) %>%
        summarise(PairedN = sum(is.finite(.data[[variable]])), .groups = "drop") %>%
        mutate(Dataset = "Clinical", Tissue = "Clinical", Modality = variable)
    }) %>% bind_rows()
  ) %>%
  select(Dataset, Tissue, Modality, AdiposityGroup, PairedN)

edge_type_summary <- all_tests %>%
  group_by(Group, TissuePair) %>%
  summarise(
    TestedPairs = sum(is.finite(PValue)),
    NominalP05 = sum(is.finite(PValue) & PValue < 0.05),
    NominalP05_Rho05 = sum(NominalEdge),
    BH_FDR05 = sum(FDR05),
    BH_FDR10 = sum(FDR10),
    BH_FDR20 = sum(FDR20),
    MedianN = median(N[is.finite(PValue)]),
    MinimumN = min(N[is.finite(PValue)]),
    MaximumN = max(N[is.finite(PValue)]),
    .groups = "drop"
  )

core_edge_table <- all_tests %>%
  filter(CoreRelated, is.finite(PValue)) %>%
  arrange(Group, BH_FDR, PValue, desc(abs(Rho)))

write_csv(summary_table, file.path(output_dir, "Fig03C_Response_Network_Statistics_Summary.csv"))

workbook <- createWorkbook()
sheet_data <- list(
  "Summary" = summary_table,
  "Edge_Type_Summary" = edge_type_summary,
  "Sample_Availability" = sample_availability,
  "Nodes" = node_table %>% mutate(across(where(is.factor), as.character)),
  "Displayed_Edges" = display_edges %>% arrange(Group, BH_FDR, PValue),
  "FDR05_Display_Edges" = fdr05_plot_edges %>%
    arrange(Group, BH_FDR, PValue),
  "FDR05_All_Edges" = all_tests %>%
    filter(!ClinicalClinical, FDR05, abs(Rho) >= absolute_rho_threshold) %>%
    arrange(Group, BH_FDR, PValue),
  "FDR20_Edges" = all_tests %>%
    filter(!ClinicalClinical, FDR20) %>%
    arrange(Group, BH_FDR, PValue),
  "Core_Node_Edges" = core_edge_table,
  "All_Tested_Pairs" = all_tests %>% arrange(Group, BH_FDR, PValue)
)
for (sheet_name in names(sheet_data)) {
  addWorksheet(workbook, sheet_name)
  writeData(workbook, sheet_name, sheet_data[[sheet_name]], withFilter = TRUE)
  freezePane(workbook, sheet_name, firstRow = TRUE)
  setColWidths(workbook, sheet_name, cols = seq_len(ncol(sheet_data[[sheet_name]])), widths = "auto")
}
saveWorkbook(
  workbook,
  file.path(output_dir, "Fig03C_Adiposity_Response_Networks_SourceData.xlsx"),
  overwrite = TRUE
)

readme_lines <- c(
  "# Revised Figure 3C adiposity-stratified response networks",
  "",
  "## Scope",
  "",
  "This is a candidate revision and does not replace the assembled Figure 3 artwork.",
  "",
  "## Analysis",
  "",
  "- Molecular response: Post3h minus Pre for serum metabolomics, serum proteomics, adipose proteomics, and skeletal-muscle proteomics.",
  "- Adiposity groups: Obese if pre-exercise body-fat percentage >30%; otherwise Lean.",
  "- Shared molecular node set: top response-SD features across the full paired cohort (20 serum metabolites, 20 serum proteins, 15 adipose proteins, and 15 muscle proteins).",
  "- Fixed clinical core nodes: body-fat percentage, leptin, triglycerides, retroperitoneal fat, visceral fat, and liver fat.",
  "- Association: pairwise Spearman correlation within each adiposity group, using available paired observations for each edge.",
  paste0("- Display rule: nominal P < ", nominal_p_threshold, " and absolute rho >= ", absolute_rho_threshold, "."),
  "- Multiplicity: Benjamini-Hochberg correction across all finite candidate-node pairs separately within Lean and Obese groups.",
  "- Each panel retains only nodes with at least one displayed edge in that group and uses a Figure 2-style stress layout.",
  "- The optimized main candidate requires BH-FDR <0.05 and absolute rho >=0.50. All cross-source edges are retained; within-source edges are thinned for display by retaining each connected node's strongest relationship plus the 12 strongest within-source edges per group.",
  "- Clinical labels are shown only for nodes with at least two displayed relationships. Molecular labels prioritize clinical-linked features and high-degree nodes within each source.",
  "- The additional FDR panels retain all molecular or clinical-to-molecular edges with BH-FDR <0.20, without a separate rho cutoff.",
  "",
  "## Visual changes",
  "",
  "- Figure 2 tissue/source colors are used exactly: serum #E9B85F, adipose #E98773, muscle #72A8CC, and clinical #79B8A9.",
  "- Background tissue hulls and boxes were removed.",
  "- Unconnected nodes are not drawn.",
  "- Node shapes encode molecular layer; edge color encodes correlation direction; edge width encodes absolute Spearman rho.",
  "",
  "## Interpretation boundary",
  "",
  "This network is exploratory. An edge indicates within-group correlation, not causation or direct molecular interaction. An edge present in only one panel is not evidence of rewiring unless a formal between-group correlation-difference test is significant. Group-specific sample sizes and pairwise missingness are retained in the source-data workbook."
)
writeLines(readme_lines, file.path(output_dir, "README.md"), useBytes = TRUE)

session_lines <- c(
  paste0("Analysis ID: ", analysis_id),
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R version: ", R.version.string),
  "",
  capture.output(sessionInfo())
)
writeLines(session_lines, file.path(output_dir, "SESSION_INFO.txt"), useBytes = TRUE)

print(summary_table)
message("Outputs written to: ", output_dir)
