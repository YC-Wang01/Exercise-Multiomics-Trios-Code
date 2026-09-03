# Candidate Figure 2E integrating a pre-exercise multi-omic network with
# adiposity-associated proteins and a PQLC3-linked within-family inset.

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
  patchwork, readr, readxl, scales, stringr, tibble, tidyr, writexl
)

analysis_id <- "Fig02E_MultiOmics_PQLC3_Integrated_2026-08-24"
layout_seed <- 121L
minimum_pair_n <- 10L
strict_edge_fdr <- 0.05
label_cap <- 28L

output_dir <- path_project(
  "fig02", "MultiOmics_PQLC3_Integrated_2026-08-24"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

node_source <- path_project(
  "fig02", "Clinical_Centered_MultiOmics_CellStyle_Candidate_2026-08-24",
  "Fig02E_ClinicalCentered_MultiOmics_CellStyle_Source_Data.xlsx"
)
family_source <- path_project(
  "fig02", "Family_Similarity_Audit",
  "All_FDR_Significant_Family_Features.csv"
)
pqlc3_panel_source <- path_project(
  "fig02", "Family_Similarity_Figures",
  "Representative_Family_Adiposity_Features_SourceData.xlsx"
)
pqlc3_contrast_source <- path_project(
  "fig02", "Family_Similarity_Audit", "PQLC3_DM_DF_Direct_Contrast.csv"
)
reference_network_source <- path_project(
  "fig02", "Cell_Stress_Layout_Comparison_2026-08-24",
  "Fig02E_Stress_Layout_Comparison_Source_Data.xlsx"
)
sample_sheet_file <- path_analysis_ready("Metadata", "Project_Sample_Sheet.csv")
clinical_file <- path_analysis_ready("Metadata", "Clinical_Observed_Primary.csv")

dataset_files <- c(
  Serum_Metabonomics = path_analysis_ready("Analysis_Serum_Metabonomics.csv"),
  Serum_Proteomics = path_analysis_ready("Analysis_Serum_Proteomics.csv"),
  Adipose_Proteomics = path_analysis_ready("Analysis_Adipose_Proteomics.csv"),
  Adipose_Microarray = path_analysis_ready("Analysis_Adipose_Microarray.csv"),
  Adipose_Methylation = path_analysis_ready("Analysis_Adipose_Methylation_MValue.rds"),
  Muscle_Proteomics = path_analysis_ready("Analysis_Muscle_Proteomics.csv"),
  Muscle_Microarray = path_analysis_ready("Analysis_Muscle_Microarray.csv"),
  Muscle_Methylation = path_analysis_ready("Analysis_Muscle_Methylation_MValue.rds")
)

input_files <- c(
  CurrentNodeSource = node_source,
  FamilyFeatureAudit = family_source,
  PQLC3PanelSource = pqlc3_panel_source,
  PQLC3DirectContrast = pqlc3_contrast_source,
  ReferenceNetworkSource = reference_network_source,
  ProjectSampleSheet = sample_sheet_file,
  ClinicalPhenotypes = clinical_file,
  dataset_files
)
if (any(!file.exists(input_files))) {
  stop(
    "Missing integration input(s): ",
    paste(names(input_files)[!file.exists(input_files)], collapse = ", ")
  )
}

relative_path <- function(path) {
  root <- paste0(normalizePath(WORKSPACE_ROOT, winslash = "/"), "/")
  sub(paste0("^", stringr::fixed(root)), "", normalizePath(path, winslash = "/"))
}

input_manifest <- tibble(
  InputRole = names(input_files),
  InputFile = vapply(input_files, relative_path, character(1)),
  SHA256 = vapply(
    input_files,
    function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE),
    character(1)
  )
)

current_nodes <- readxl::read_excel(node_source, sheet = "DisplayedNodes") %>%
  as_tibble() %>%
  transmute(
    NodeID = as.character(NodeID),
    Dataset = as.character(Dataset),
    Tissue = as.character(Tissue),
    Omics = as.character(Omics),
    Modality = as.character(Modality),
    FeatureID = as.character(FeatureID),
    DisplayLabel = as.character(DisplayLabel),
    Effect = suppressWarnings(as.numeric(Effect)),
    P_Value = suppressWarnings(as.numeric(P_Value)),
    BH_FDR = suppressWarnings(as.numeric(BH_FDR)),
    SelectionThreshold = suppressWarnings(as.numeric(SelectionThreshold)),
    Direction = as.character(Direction),
    Group = as.character(Group),
    StrictFDR = as.logical(StrictFDR),
    EvidenceClass = as.character(EvidenceClass),
    SerumModule = as.character(SerumModule)
  ) %>%
  mutate(
    FamilyStrict = FALSE,
    FamilyPair = NA_character_,
    FamilyR = NA_real_,
    FamilyP = NA_real_,
    FamilyBH_FDR = NA_real_,
    FamilyNPairs = NA_integer_,
    FamilyLOFO = NA_character_
  )

family_audit <- readr::read_csv(family_source, show_col_types = FALSE) %>%
  mutate(across(everything(), ~ .x))

family_spec <- tibble::tribble(
  ~Dataset, ~FeatureID, ~DisplayLabel, ~FamilyPair, ~Modality, ~Group,
  "Serum_Proteomics", "Q96JM4;Q9HAE3", "LRRIQ1", "DF", "Protein", "Serum factors",
  "Serum_Proteomics", "P22352", "GPX3", "DF", "Protein", "Serum factors",
  "Muscle_Proteomics", "P24310", "COX7A1", "DM", "Protein", "Muscle",
  "Adipose_Methylation", "cg06510534", "PQLC3", "DM", "DNA methylation", "Adipose"
)

family_rows <- family_spec %>%
  left_join(
    family_audit,
    by = c("Dataset", "FeatureID")
  ) %>%
  rowwise() %>%
  mutate(
    NodeID = if_else(
      Dataset == "Serum_Proteomics",
      paste0("Prot_", str_replace_all(FeatureID, ";", "_")),
      paste0(Dataset, "__", FeatureID)
    ),
    Effect = as.numeric(AdiposityGroupEffect),
    P_Value = as.numeric(AdiposityGroupP),
    BH_FDR = as.numeric(AdiposityGroupBH_FDR),
    SelectionThreshold = if_else(Modality == "Protein", 0.20, 0.05),
    Direction = if_else(Effect < 0, "Higher in lean", "Higher in obese"),
    StrictFDR = is.finite(BH_FDR) & BH_FDR < 0.05,
    EvidenceClass = "Parent-daughter correlation BH-FDR < 0.05",
    SerumModule = NA_character_,
    FamilyStrict = TRUE,
    FamilyR = if_else(FamilyPair == "DM", as.numeric(DM_Correlation), as.numeric(DF_Correlation)),
    FamilyP = if_else(FamilyPair == "DM", as.numeric(DM_SimilarityP), as.numeric(DF_SimilarityP)),
    FamilyBH_FDR = if_else(
      FamilyPair == "DM",
      as.numeric(DM_SimilarityBH_FDR),
      as.numeric(DF_SimilarityBH_FDR)
    ),
    FamilyNPairs = if_else(FamilyPair == "DM", as.integer(DM_NPairs), as.integer(DF_NPairs)),
    FamilyLOFO = if_else(FamilyPair == "DM", as.character(DM_StabilityFlag), as.character(DF_StabilityFlag))
  ) %>%
  ungroup() %>%
  select(
    NodeID, Dataset, Tissue, Omics, Modality, FeatureID, DisplayLabel,
    Effect, P_Value, BH_FDR, SelectionThreshold, Direction, Group,
    StrictFDR, EvidenceClass, SerumModule, FamilyStrict, FamilyPair,
    FamilyR, FamilyP, FamilyBH_FDR, FamilyNPairs, FamilyLOFO
  )

if (any(!is.finite(family_rows$FamilyBH_FDR)) || any(family_rows$FamilyBH_FDR >= 0.05)) {
  stop("The selected within-family rows do not all satisfy BH-FDR < 0.05.")
}

pqlc3_dm <- readxl::read_excel(pqlc3_panel_source, sheet = "Daughter-mother pairs") %>%
  filter(FeatureID == "cg06510534") %>%
  transmute(
    FamilyID = as.character(FamilyID), Pair = "Daughter-mother",
    ParentValue = as.numeric(ParentValue), DaughterValue = as.numeric(DaughterValue)
  )
pqlc3_df <- readxl::read_excel(pqlc3_panel_source, sheet = "Daughter-father pairs") %>%
  filter(FeatureID == "cg06510534") %>%
  transmute(
    FamilyID = as.character(FamilyID), Pair = "Daughter-father",
    ParentValue = as.numeric(ParentValue), DaughterValue = as.numeric(DaughterValue)
  )
pqlc3_summary <- readxl::read_excel(pqlc3_panel_source, sheet = "Candidate summary") %>%
  filter(FeatureID == "cg06510534")
pqlc3_contrast <- readr::read_csv(pqlc3_contrast_source, show_col_types = FALSE)
if (nrow(pqlc3_dm) != 10L || nrow(pqlc3_df) != 10L ||
    nrow(pqlc3_summary) != 1L || nrow(pqlc3_contrast) != 1L) {
  stop("Unexpected PQLC3 panel source dimensions.")
}
pqlc3_pairs <- bind_rows(pqlc3_dm, pqlc3_df) %>%
  mutate(Pair = factor(Pair, levels = c("Daughter-mother", "Daughter-father")))

node_table <- bind_rows(current_nodes, family_rows) %>%
  arrange(desc(FamilyStrict), Dataset, BH_FDR, DisplayLabel) %>%
  distinct(NodeID, .keep_all = TRUE)

sample_sheet <- readr::read_csv(sample_sheet_file, show_col_types = FALSE) %>%
  mutate(
    Dataset = as.character(Dataset),
    SampleKey = tolower(as.character(SampleID)),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Timepoint = as.character(Timepoint),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort)
  ) %>%
  filter(Timepoint == "Pre", Is_Primary_Cohort)

clinical_variables <- node_table %>%
  filter(Dataset == "Clinical") %>%
  pull(FeatureID)

clinical <- readr::read_csv(clinical_file, show_col_types = FALSE) %>%
  transmute(
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    FamilyID = as.character(FamilyID),
    Role = as.character(Role),
    across(all_of(clinical_variables), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  distinct(Clinical_Subject_ID, .keep_all = TRUE)

read_node_values <- function(dataset, nodes) {
  mat <- read_analysis_matrix(dataset_files[[dataset]])
  row_index <- match(nodes$FeatureID, rownames(mat))
  missing_features <- nodes$FeatureID[is.na(row_index)]
  if (length(missing_features)) {
    warning(dataset, " missing feature(s): ", paste(missing_features, collapse = ", "))
  }
  nodes <- nodes[!is.na(row_index), , drop = FALSE]
  row_index <- row_index[!is.na(row_index)]
  mat <- mat[row_index, , drop = FALSE]
  rownames(mat) <- nodes$NodeID

  mapping <- sample_sheet %>%
    filter(Dataset == dataset) %>%
    distinct(SampleKey, Clinical_Subject_ID)
  subject_ids <- mapping$Clinical_Subject_ID[match(tolower(colnames(mat)), mapping$SampleKey)]
  keep <- !is.na(subject_ids) & nzchar(subject_ids)
  mat <- mat[, keep, drop = FALSE]
  subject_ids <- subject_ids[keep]
  if (anyDuplicated(subject_ids)) {
    unique_subjects <- unique(subject_ids)
    mat <- vapply(unique_subjects, function(id) {
      rowMeans(mat[, subject_ids == id, drop = FALSE], na.rm = TRUE)
    }, numeric(nrow(mat)))
    rownames(mat) <- nodes$NodeID
    colnames(mat) <- unique_subjects
  } else {
    colnames(mat) <- subject_ids
  }
  list(nodes = nodes, matrix = mat)
}

omics_datasets <- setdiff(unique(node_table$Dataset), "Clinical")
matrix_pieces <- lapply(omics_datasets, function(dataset) {
  read_node_values(dataset, node_table %>% filter(Dataset == dataset))
})
names(matrix_pieces) <- omics_datasets

subjects <- sort(unique(c(
  clinical$Clinical_Subject_ID,
  unlist(lapply(matrix_pieces, function(x) colnames(x$matrix)))
)))

network_matrix <- matrix(
  NA_real_, nrow = nrow(node_table), ncol = length(subjects),
  dimnames = list(node_table$NodeID, subjects)
)

for (variable in clinical_variables) {
  node_id <- node_table$NodeID[node_table$Dataset == "Clinical" & node_table$FeatureID == variable][1]
  network_matrix[node_id, match(clinical$Clinical_Subject_ID, subjects)] <- clinical[[variable]]
}
for (piece in matrix_pieces) {
  network_matrix[rownames(piece$matrix), match(colnames(piece$matrix), subjects)] <- piece$matrix
}

metadata <- clinical %>%
  select(Clinical_Subject_ID, FamilyID, Role) %>%
  right_join(tibble(Clinical_Subject_ID = subjects), by = "Clinical_Subject_ID")

pair_list <- utils::combn(node_table$NodeID, 2L, simplify = FALSE)
edge_tests <- bind_rows(lapply(pair_list, function(pair) {
  from_meta <- node_table %>% filter(NodeID == pair[[1]])
  to_meta <- node_table %>% filter(NodeID == pair[[2]])
  if (from_meta$Dataset[[1]] == to_meta$Dataset[[1]]) return(NULL)
  if (from_meta$Dataset[[1]] == "Clinical" && to_meta$Dataset[[1]] == "Clinical") return(NULL)

  x <- as.numeric(network_matrix[pair[[1]], ])
  y <- as.numeric(network_matrix[pair[[2]], ])
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  family_n <- n_distinct(metadata$FamilyID[ok & !is.na(metadata$FamilyID)])
  if (n < minimum_pair_n || length(unique(x[ok])) < 3L || length(unique(y[ok])) < 3L) {
    return(tibble(
      From = pair[[1]], To = pair[[2]], N = n, Families = family_n,
      Rho = NA_real_, P_Value = NA_real_
    ))
  }
  test <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(
    From = pair[[1]], To = pair[[2]], N = n, Families = family_n,
    Rho = unname(test$estimate), P_Value = test$p.value
  )
})) %>%
  mutate(
    BH_FDR = p.adjust(P_Value, method = "BH"),
    EdgeSign = if_else(Rho >= 0, "Positive", "Negative"),
    StrictPass = is.finite(BH_FDR) & BH_FDR < strict_edge_fdr
  )

strict_edges <- edge_tests %>% filter(StrictPass)

node_attributes <- node_table %>%
  select(NodeID, Group, Modality, FamilyStrict)

edge_endpoints <- strict_edges %>%
  mutate(EdgeID = row_number()) %>%
  pivot_longer(c(From, To), names_to = "EndpointSide", values_to = "Endpoint") %>%
  left_join(node_attributes, by = c("Endpoint" = "NodeID")) %>%
  mutate(
    EdgeCap = case_when(
      FamilyStrict ~ 6L,
      Group == "Clinical phenotype" ~ 7L,
      Modality == "Protein" ~ 4L,
      TRUE ~ 3L
    )
  ) %>%
  group_by(Endpoint) %>%
  arrange(BH_FDR, desc(abs(Rho)), .by_group = TRUE) %>%
  mutate(EndpointRank = row_number()) %>%
  filter(EndpointRank <= first(EdgeCap)) %>%
  ungroup()

display_edge_ids <- sort(unique(edge_endpoints$EdgeID))
display_edges <- strict_edges %>%
  mutate(EdgeID = row_number()) %>%
  filter(EdgeID %in% display_edge_ids) %>%
  select(-EdgeID)

edge_node_ids <- unique(c(display_edges$From, display_edges$To))
required_node_ids <- node_table %>%
  filter(
    StrictFDR | FamilyStrict |
      (Modality == "Protein" & is.finite(BH_FDR) & BH_FDR <= SelectionThreshold)
  ) %>%
  pull(NodeID) %>%
  unique()
display_nodes <- node_table %>%
  filter(NodeID %in% union(edge_node_ids, required_node_ids))

node_inclusion_audit <- node_table %>%
  mutate(
    RequiredForDisplay = NodeID %in% required_node_ids,
    HasDisplayedEdge = NodeID %in% edge_node_ids,
    IncludedInNetwork = NodeID %in% display_nodes$NodeID
  ) %>%
  filter(RequiredForDisplay | DisplayLabel == "PQLC3") %>%
  arrange(desc(DisplayLabel == "PQLC3"), Dataset, BH_FDR, DisplayLabel)

if (!all(node_inclusion_audit$IncludedInNetwork) ||
    !any(node_inclusion_audit$DisplayLabel == "PQLC3")) {
  stop("Required strict-FDR/PQLC3 nodes were not retained in the network.")
}

graph <- igraph::graph_from_data_frame(
  display_edges %>%
    transmute(from = From, to = To, Rho, BH_FDR, EdgeSign, EdgeWeight = abs(Rho)),
  directed = FALSE,
  vertices = display_nodes %>%
    transmute(
      name = NodeID, Group, Modality, DisplayLabel, SerumModule,
      StrictFDR, NodeBH_FDR = BH_FDR, FamilyStrict, FamilyPair,
      FamilyR, FamilyBH_FDR
    )
)

igraph::V(graph)$Degree <- igraph::degree(graph)
set.seed(layout_seed)
layout <- ggraph::create_layout(graph, layout = "stress")
layout_table <- as_tibble(layout) %>% mutate(NodeID = as.character(name))

label_table <- layout_table %>%
  filter(is.na(SerumModule) | !nzchar(SerumModule)) %>%
  mutate(
    LabelPriority = case_when(
      FamilyStrict ~ 1L,
      Group == "Clinical phenotype" ~ 2L,
      Modality == "Protein" ~ 3L,
      StrictFDR ~ 4L,
      TRUE ~ 5L
    ),
    StableFDR = if_else(is.finite(NodeBH_FDR), NodeBH_FDR, Inf)
  ) %>%
  arrange(LabelPriority, desc(Degree), StableFDR, DisplayLabel) %>%
  slice_head(n = label_cap)

label_sides <- tibble::tribble(
  ~DisplayLabel, ~LabelSide,
  "PQLC3", "above",
  "GPX3", "right",
  "LRRIQ1", "above",
  "COX7A1", "left",
  "Visceral fat", "below",
  "Triglycerides", "left",
  "Retroperitoneal fat", "left",
  "Body-fat percentage", "above",
  "Leptin", "left",
  "Liver fat", "below",
  "TLCD4", "left",
  "CUL4B", "below",
  "C9", "below",
  "CLMP", "left",
  "SAA4", "above",
  "SERPINA3", "above",
  "PCMTD1", "left",
  "GSDMB", "above",
  "CENPV", "below",
  "IGDCC4", "above",
  "SENP2", "left",
  "CSNK2A1", "below",
  "SVIP", "above",
  "OSBPL3", "below",
  "TBC1D20", "below",
  "OTOG", "left",
  "SCUBE1", "above",
  "IRAK3", "right"
)

label_table <- label_table %>%
  left_join(label_sides, by = "DisplayLabel") %>%
  mutate(
    LabelOffset = 0.11 + pmin(Degree, 20) * 0.002,
    LabelX = case_when(
      LabelSide == "left" ~ x - LabelOffset,
      LabelSide == "right" ~ x + LabelOffset,
      TRUE ~ x
    ),
    LabelY = case_when(
      LabelSide == "above" ~ y + LabelOffset,
      LabelSide == "below" ~ y - LabelOffset,
      TRUE ~ y
    ),
    LabelHJust = case_when(
      LabelSide == "left" ~ 1,
      LabelSide == "right" ~ 0,
      TRUE ~ 0.5
    ),
    LabelVJust = case_when(
      LabelSide == "above" ~ 0,
      LabelSide == "below" ~ 1,
      TRUE ~ 0.5
    )
  )

if (any(is.na(label_table$LabelSide))) {
  stop(
    "Missing fixed label position for: ",
    paste(label_table$DisplayLabel[is.na(label_table$LabelSide)], collapse = ", ")
  )
}

group_colors <- c(
  "Adipose" = "#E98773",
  "Clinical phenotype" = "#79B8A9",
  "Muscle" = "#72A8CC",
  "Serum factors" = "#E9B85F"
)
module_colors <- c(
  "HDL remodeling" = "#C6E4E5",
  "VLDL/TG composition" = "#F7D6B6"
)
edge_colors <- c("Positive" = "#D96B63", "Negative" = "#5D91B8")
modality_shapes <- c(
  "Clinical phenotype" = 22,
  "Metabolite" = 23,
  "Protein" = 21,
  "Transcript" = 24,
  "DNA methylation" = 25
)

network_plot <- ggraph::ggraph(layout) +
  ggraph::geom_edge_link(
    aes(edge_width = EdgeWeight, edge_alpha = EdgeWeight, edge_color = EdgeSign),
    lineend = "round"
  ) +
  ggforce::geom_mark_hull(
    data = layout_table,
    aes(
      x = x, y = y, group = SerumModule,
      filter = !is.na(SerumModule) & nzchar(SerumModule),
      fill = SerumModule, label = SerumModule
    ),
    concavity = 1.5,
    expand = grid::unit(2.0, "mm"),
    radius = grid::unit(2.0, "mm"),
    alpha = 0.14,
    color = "#777777",
    linetype = "dashed",
    linewidth = 0.25,
    label.fontsize = 6.2,
    label.fontface = "bold",
    label.fill = "white",
    label.buffer = grid::unit(2.0, "mm"),
    con.colour = "#777777",
    con.cap = grid::unit(0, "mm"),
    show.legend = FALSE
  ) +
  ggraph::geom_node_point(
    aes(size = Degree, fill = Group, shape = Modality),
    color = "white", stroke = 0.50
  ) +
  geom_text(
    data = label_table,
    aes(
      x = LabelX, y = LabelY, label = DisplayLabel,
      hjust = LabelHJust, vjust = LabelVJust
    ),
    inherit.aes = FALSE,
    family = "Arial", fontface = "bold", size = 6.0 / ggplot2::.pt,
    color = "black"
  ) +
  ggraph::scale_edge_color_manual(values = edge_colors, name = "Association direction") +
  ggraph::scale_edge_alpha_continuous(range = c(0.10, 0.44), guide = "none") +
  ggraph::scale_edge_width_continuous(range = c(0.18, 0.68), guide = "none") +
  scale_fill_manual(
    values = c(group_colors, module_colors),
    breaks = names(group_colors),
    name = "Tissue / clinical"
  ) +
  scale_shape_manual(
    values = modality_shapes,
    breaks = names(modality_shapes),
    labels = c("Clinical", "Metabolite", "Protein", "Transcript", "DNAm"),
    name = "Omics layer"
  ) +
  scale_size_continuous(range = c(2.1, 6.0), guide = "none") +
  coord_equal(clip = "off") +
  theme_void(base_family = "Arial", base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.title = element_text(face = "bold", size = 6.5, color = "black"),
    legend.text = element_text(size = 5.8, color = "black"),
    legend.key.height = grid::unit(3.5, "mm"),
    legend.spacing.y = grid::unit(1.2, "mm"),
    legend.box.margin = margin(0, 0, 1.5, 0, unit = "mm"),
    plot.margin = margin(2, 2, 4, 2, unit = "mm")
  ) +
  guides(
    fill = guide_legend(order = 1, nrow = 1, byrow = TRUE,
                        override.aes = list(shape = 21, size = 2.8, alpha = 1)),
    shape = guide_legend(order = 2, nrow = 1, byrow = TRUE,
                         override.aes = list(size = 2.8, fill = "#D9D9D9")),
    edge_color = guide_legend(order = 3, nrow = 1,
                              override.aes = list(edge_width = 0.65, edge_alpha = 0.8))
  )

shared_family_protein_labels <- c("LRRIQ1", "GPX3", "COX7A1")
adiposity_priority_protein_labels <- c("PCMTD1", "SERPINA3", "SHBG", "TLCD4", "CLMP")
protein_bar_labels <- c(shared_family_protein_labels, adiposity_priority_protein_labels)

adiposity_proteins <- node_table %>%
  filter(Modality == "Protein", DisplayLabel %in% protein_bar_labels) %>%
  mutate(
    SharedFamilyNetwork = DisplayLabel %in% shared_family_protein_labels,
    EvidenceTier = case_when(
      BH_FDR < 0.05 ~ "Adiposity BH-FDR < 0.05",
      BH_FDR < 0.20 ~ "Adiposity BH-FDR 0.05-0.20",
      TRUE ~ "Family-network core; adiposity BH-FDR >= 0.20"
    ),
    PointClass = if_else(SharedFamilyNetwork, "Also family-associated", "Adiposity-selected"),
    TissueLabel = paste0(DisplayLabel, " (", tolower(Tissue), ")"),
    PlotLabel = factor(TissueLabel, levels = TissueLabel[order(Effect)]),
    QLabel = paste0("q=", formatC(BH_FDR, digits = 2, format = "g"))
  ) %>%
  arrange(Effect)

if (nrow(adiposity_proteins) != 8L ||
    !all(shared_family_protein_labels %in% adiposity_proteins$DisplayLabel)) {
  stop("The updated Figure 2E protein bar must contain eight proteins, including all shared family-network proteins.")
}

adiposity_plot <- ggplot(adiposity_proteins, aes(y = PlotLabel, x = Effect)) +
  geom_vline(xintercept = 0, linewidth = 0.28, color = "black") +
  geom_segment(aes(x = 0, xend = Effect, yend = PlotLabel, color = Group,
                   alpha = EvidenceTier), linewidth = 1.8, lineend = "round") +
  geom_point(aes(fill = Group, alpha = EvidenceTier, shape = PointClass),
             color = "black", stroke = 0.30, size = 2.6) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_alpha_manual(
    values = c(
      "Adiposity BH-FDR < 0.05" = 1,
      "Adiposity BH-FDR 0.05-0.20" = 0.46,
      "Family-network core; adiposity BH-FDR >= 0.20" = 0.28
    ),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c("Adiposity-selected" = 21, "Also family-associated" = 23),
    guide = "none"
  ) +
  scale_x_continuous(limits = c(-1.70, 1.10), breaks = c(-1.5, -0.5, 0.5, 1.0)) +
  labs(
    title = "Adiposity effects",
    subtitle = "Diamonds: family-associated\nShade: adiposity BH-FDR",
    x = "Obese - lean difference\n(log2 relative abundance)", y = NULL
  ) +
  theme_classic(base_family = "Arial", base_size = 6.5) +
  theme(
    plot.title = element_text(face = "bold", size = 8, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 5.2, hjust = 0.5, color = "black"),
    axis.title.x = element_text(face = "bold", size = 6.5, color = "black"),
    axis.text = element_text(size = 6.2, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.30),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.30, color = "black"),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )

pqlc3_pair_colors <- c(
  "Daughter-mother" = "#B65C6B",
  "Daughter-father" = "#3A78B8"
)
pqlc3_pair_shapes <- c("Daughter-mother" = 16, "Daughter-father" = 17)
format_q <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < 0.001) return(format(x, scientific = TRUE, digits = 2))
  formatC(x, format = "f", digits = 3)
}
pqlc3_annotation <- paste0(
  "DM: r = ", sprintf("%.2f", pqlc3_summary$DM_PearsonR[[1]]),
  ", FDR = ", format_q(pqlc3_summary$DM_BH_FDR[[1]]), "\n",
  "DF: r = ", sprintf("%.2f", pqlc3_summary$DF_PearsonR[[1]]),
  ", FDR = ", format_q(pqlc3_summary$DF_BH_FDR[[1]]), "\n",
  "DM-DF: Delta r = ", sprintf("%.2f", pqlc3_contrast$DeltaR_DM_minus_DF[[1]]),
  ", exact P = ", format_q(pqlc3_contrast$ExactSwapPermutationP[[1]])
)

pqlc3_plot <- ggplot(
  pqlc3_pairs,
  aes(x = ParentValue, y = DaughterValue, colour = Pair, shape = Pair)
) +
  geom_smooth(
    data = pqlc3_pairs,
    aes(
      x = ParentValue, y = DaughterValue, colour = Pair,
      fill = Pair, group = Pair
    ),
    inherit.aes = FALSE, method = "lm", formula = y ~ x,
    se = TRUE, alpha = 0.10, linewidth = 0.60, show.legend = FALSE
  ) +
  geom_point(size = 2.0, stroke = 0.30) +
  scale_colour_manual(
    values = pqlc3_pair_colors,
    labels = c("Daughter-mother" = "DM", "Daughter-father" = "DF"),
    name = "Family pair"
  ) +
  scale_fill_manual(values = pqlc3_pair_colors) +
  scale_shape_manual(
    values = pqlc3_pair_shapes,
    labels = c("Daughter-mother" = "DM", "Daughter-father" = "DF"),
    name = "Family pair"
  ) +
  labs(
    title = "PQLC3 DNAm similarity",
    subtitle = pqlc3_annotation,
    x = "Parental M-value", y = "Daughter M-value"
  ) +
  coord_equal() +
  theme_classic(base_family = "Arial", base_size = 6.5) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.30),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold", size = 8, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 5.2, hjust = 0, color = "black", lineheight = 1.05),
    axis.title = element_text(face = "bold", size = 6.2, color = "black"),
    axis.text = element_text(size = 5.8, color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 5.8),
    legend.text = element_text(size = 5.4),
    legend.key.height = grid::unit(3.0, "mm"),
    legend.box.margin = margin(0, 0, 1.5, 0, unit = "mm"),
    plot.margin = margin(2, 2, 4, 2, unit = "mm")
  ) +
  guides(
    colour = guide_legend(order = 1, nrow = 1),
    shape = guide_legend(order = 1, nrow = 1)
  )

side_plots <- adiposity_plot / pqlc3_plot +
  patchwork::plot_layout(heights = c(1.25, 1))

integrated_plot <- network_plot + side_plots +
  patchwork::plot_layout(widths = c(2.25, 1))

file_integrated <- file.path(
  output_dir, "Fig02E_MultiOmics_PQLC3_Integrated_Candidate.pdf"
)
file_network <- file.path(
  output_dir, "Fig02E_MultiOmics_PQLC3_Network_Candidate.pdf"
)
file_side <- file.path(
  output_dir, "Fig02E_ProteinBar_and_PQLC3_SidePlots_Candidate.pdf"
)

ggsave(file_integrated, integrated_plot, width = 190, height = 128, units = "mm",
       device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".png", file_integrated), integrated_plot,
       width = 190, height = 128, units = "mm", dpi = 300, bg = "white")
ggsave(file_network, network_plot, width = 142, height = 128, units = "mm",
       device = grDevices::cairo_pdf)
ggsave(file_side, side_plots, width = 62, height = 128, units = "mm",
       device = grDevices::cairo_pdf)

family_connectivity <- family_rows %>%
  transmute(
    NodeID, DisplayLabel, Dataset, FamilyPair, FamilyNPairs, FamilyR,
    FamilyP, FamilyBH_FDR, FamilyLOFO,
    IncludedInNetwork = NodeID %in% layout_table$NodeID,
    NetworkDegree = if_else(
      IncludedInNetwork,
      as.integer(layout_table$Degree[match(NodeID, layout_table$NodeID)]),
      0L
    )
  )

edge_export <- igraph::as_data_frame(graph, what = "edges") %>%
  as_tibble() %>%
  rename(From = from, To = to) %>%
  arrange(BH_FDR, desc(abs(Rho)))

reference_nodes <- readxl::read_excel(reference_network_source, sheet = "A_Nodes")
reference_edges <- readxl::read_excel(reference_network_source, sheet = "A_Edges") %>%
  mutate(EdgeKey = if_else(From < To, paste(From, To, sep = "||"), paste(To, From, sep = "||")))
expanded_edge_keys <- edge_export %>%
  mutate(EdgeKey = if_else(From < To, paste(From, To, sep = "||"), paste(To, From, sep = "||")))
family_node_ids <- family_rows$NodeID
comparison_summary <- tibble(
  Metric = c(
    "Displayed nodes", "Displayed edges", "Shared displayed edges",
    "Reference-only displayed edges", "Expanded-only displayed edges",
    "Expanded-only edges incident to family proteins"
  ),
  Reference = c(nrow(reference_nodes), nrow(reference_edges), NA, NA, NA, NA),
  Expanded = c(
    nrow(layout_table), nrow(edge_export),
    sum(reference_edges$EdgeKey %in% expanded_edge_keys$EdgeKey),
    sum(!reference_edges$EdgeKey %in% expanded_edge_keys$EdgeKey),
    sum(!expanded_edge_keys$EdgeKey %in% reference_edges$EdgeKey),
    sum(
      !expanded_edge_keys$EdgeKey %in% reference_edges$EdgeKey &
        (expanded_edge_keys$From %in% family_node_ids | expanded_edge_keys$To %in% family_node_ids)
    )
  )
)

writexl::write_xlsx(
  list(
    AdiposityProteinBars = adiposity_proteins %>%
      select(Dataset, Tissue, FeatureID, DisplayLabel, Effect, P_Value, BH_FDR,
             SelectionThreshold, Direction, EvidenceTier, SharedFamilyNetwork),
    FamilyProteinCandidates = family_rows %>%
      filter(Modality == "Protein") %>%
      select(Dataset, Tissue, FeatureID, DisplayLabel, FamilyPair, FamilyNPairs,
             FamilyR, FamilyP, FamilyBH_FDR, FamilyLOFO),
    PQLC3FamilyPairs = pqlc3_pairs,
    PQLC3DirectContrast = pqlc3_contrast,
    RequiredNodeAudit = node_inclusion_audit,
    FamilyNetworkConnectivity = family_connectivity,
    ComparisonSummary = comparison_summary,
    DisplayedNodes = layout_table %>% arrange(Group, desc(Degree), DisplayLabel),
    DisplayedEdges = edge_export,
    AllTestedEdges = edge_tests %>% arrange(BH_FDR, P_Value),
    AllCandidateNodes = node_table %>% arrange(Group, Modality, BH_FDR, DisplayLabel),
    AnalysisParameters = tibble(
      Parameter = c(
        "Analysis ID", "Time point", "Adiposity model",
        "Family feature evidence", "Network association",
        "Network edge multiplicity", "Network edge threshold",
        "Minimum complete pairs", "Layout", "Random seed",
        "Required node rule", "Node shape", "Interpretation boundary"
      ),
      Value = c(
        analysis_id, "Pre-exercise", "Feature ~ adiposity group + Role; FamilyID block",
        "Pearson correlation within DM or DF; BH within dataset and pair type",
        "Spearman correlation across measured participants",
        "BH across all tested cross-dataset candidate pairs",
        "BH-FDR < 0.05", as.character(minimum_pair_n),
        "Stress layout", as.character(layout_seed),
        "All strict-FDR selected nodes, adiposity-associated proteins, within-family candidates, and PQLC3",
        "Measurement layer",
        "Network edges and family correlations are distinct statistical quantities"
      )
    ),
    InputManifest = input_manifest
  ),
  file.path(output_dir, "Fig02E_MultiOmics_PQLC3_Integration_Source_Data.xlsx")
)

readme <- c(
  "# Figure 2E pre-exercise multi-omic and PQLC3 integration candidate",
  "",
  paste0("Analysis ID: `", analysis_id, "`."),
  "",
  "## Scope",
  "",
  "- The candidate outputs in this directory were regenerated from the recorded source data and code.",
  "- The network uses pre-exercise observations and recomputed cross-dataset Spearman correlations.",
  "- Every displayed network edge satisfies global BH-FDR < 0.05 in the expanded candidate set.",
  "- The upper side panel shows eight proteins, prioritizing COX7A1, LRRIQ1 and GPX3 because the same protein measurements also reach family-pair BH-FDR < 0.05; PQLC3 is retained in the lower panel.",
  "- PCMTD1, SERPINA3, SHBG, TLCD4 and CLMP complete the protein bar using the strongest and directionally balanced adiposity-associated protein candidates.",
  "- The lower side panel combines the audited PQLC3 daughter-mother and daughter-father correlations.",
  "- Node colour denotes tissue/module and node shape denotes measurement layer.",
  "- Within-family correlation is not encoded as a network edge.",
  "- Network positions and distances are graphical artifacts.",
  "",
  "## Comparison with the reference stress-layout candidate",
  "",
  paste0("- Displayed nodes: ", nrow(reference_nodes), " -> ", nrow(layout_table), "."),
  paste0("- Displayed edges: ", nrow(reference_edges), " -> ", nrow(edge_export), "."),
  paste0(
    "- Shared displayed edges: ",
    sum(reference_edges$EdgeKey %in% expanded_edge_keys$EdgeKey), "."
  ),
  paste0(
    "- Expanded-only edges incident to family proteins: ",
    sum(
      !expanded_edge_keys$EdgeKey %in% reference_edges$EdgeKey &
        (expanded_edge_keys$From %in% family_node_ids | expanded_edge_keys$To %in% family_node_ids)
    ), "."
  ),
  "",
  "## Outputs",
  "",
  "- `Fig02E_MultiOmics_PQLC3_Integrated_Candidate.pdf`",
  "- `Fig02E_MultiOmics_PQLC3_Network_Candidate.pdf`",
  "- `Fig02E_ProteinBar_and_PQLC3_SidePlots_Candidate.pdf`",
  "- `Fig02E_MultiOmics_PQLC3_Integration_Source_Data.xlsx`"
)
writeLines(readme, file.path(output_dir, "README.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "SESSION_INFO.txt"), useBytes = TRUE)

message("Candidate nodes: ", nrow(node_table))
message("Strict expanded edges: ", nrow(strict_edges))
message("Displayed nodes: ", nrow(layout_table))
message("Displayed edges: ", nrow(edge_export))
message(
  "Family proteins connected: ",
  paste(family_connectivity$DisplayLabel[family_connectivity$IncludedInNetwork], collapse = ", ")
)
