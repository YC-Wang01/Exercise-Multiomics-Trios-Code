.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(.file_arg)) stop("Cannot resolve the current script path.")
.script_file <- normalizePath(sub("^--file=", "", .file_arg[[1]]), winslash = "/")
.local_config <- file.path(dirname(.script_file), "Fig00_Config.R")
source(.local_config)

p_load(
  dplyr, tidyr, tibble, readr, stringr, purrr,
  ggplot2, igraph, ggraph, ggforce, scales, patchwork
)

set.seed(20260831)

output_dir <- path_results("Figure_6_Osteocalcin_Modular_Network")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

feature_file <- path_project(
  "Pre_Post_Response_Analysis", "Results",
  "08_Osteocalcin_Molecular_State_and_Response", "01_Feature_Level",
  "OC_Matched_State_Feature_Candidates.csv.gz"
)
sample_file <- path_analysis_ready("Metadata", "Project_Sample_Sheet.csv")

fdr_node_threshold <- 0.10
fdr_edge_threshold <- 0.10
absolute_rho_threshold <- 0.45
minimum_pair_n <- 12L
minimum_lofo_stability <- 0.75

state_order <- c("Pre", "Post", "Delta")
state_labels <- c(
  Pre = "Pre-exercise",
  Post = "Post-exercise",
  Delta = "Response"
)
marker_order <- c("cOC", "tOC", "cOC/tOC ratio")
marker_labels <- c("cOC", "tOC", "cOC/tOC ratio")
tissue_colors <- c(Serum = "#4F9A8D", Adipose = "#DDAA3C", Muscle = "#806FAF")
association_colors <- c(Positive = "#C85A52", Negative = "#3E6E9E")
module_audit <- list()

matrix_registry <- tribble(
  ~Dataset, ~InputFile,
  "Serum_Proteomics", "Analysis_Serum_Proteomics.csv",
  "Serum_Metabonomics", "Analysis_Serum_Metabonomics.csv",
  "Adipose_Microarray", "Analysis_Adipose_Microarray.csv",
  "Muscle_Microarray", "Analysis_Muscle_Microarray.csv",
  "Adipose_Proteomics", "Analysis_Adipose_Proteomics.csv",
  "Muscle_Proteomics", "Analysis_Muscle_Proteomics.csv",
  "Adipose_Methylation", "Analysis_Adipose_Methylation_MValue.rds",
  "Muscle_Methylation", "Analysis_Muscle_Methylation_MValue.rds"
)

feature_display_labels <- c(
  AcAce = "Acetoacetate", Ace = "Acetate", Ala = "Alanine",
  ApoBtoApoA1 = "ApoB/ApoA1 ratio", Cit = "Citrate",
  Gln = "Glutamine", Gly = "Glycine", His = "Histidine",
  Ile = "Isoleucine", Lac = "Lactate", Leu = "Leucine",
  Pyr = "Pyruvate", Serum_C = "Serum cholesterol",
  Serum_TG = "Serum triglycerides", TGtoPG = "TG/PG ratio",
  TotFA = "Total fatty acids", TotPG = "Total phosphoglycerides",
  Tyr = "Tyrosine", bOHBut = "3-hydroxybutyrate",
  IDL_PL = "IDL phospholipids", Glc = "Glucose",
  VLDL_TG = "VLDL triglycerides", FAw6toFA = "Omega-6/total FA",
  FAw79StoFA = "PUFA/total FA", FAw79S = "PUFA"
)

classify_feature <- function(feature_id, display, dataset, tissue) {
  case_when(
    str_detect(feature_id, "(^|_)(VLDL|LDL|HDL|IDL)(_|$)|ApoB|Serum_C|Serum_TG|TGtoPG|TotPG|EstC|^PC$") ~
      "Lipoprotein remodeling",
    feature_id %in% c("Ala", "Gln", "Gly", "His", "Ile", "Leu", "Tyr") ~
      "Amino-acid metabolism",
    feature_id %in% c("Cit", "Lac", "Pyr", "Glc") ~ "Central carbon metabolism",
    feature_id %in% c("AcAce", "Ace", "bOHBut") ~ "Ketone metabolism",
    str_detect(feature_id, "FA|DHA|^LA$|MUFA|PUFA|CH2|DBin|BISto|otPUFA") ~
      "Fatty-acid composition",
    str_detect(display, regex("C4B|CFH|SERPINA|COMPLEMENT", ignore_case = TRUE)) ~
      "Complement and acute-phase proteins",
    str_detect(dataset, "Methylation") ~ paste(tissue, "DNA methylation"),
    str_detect(dataset, "Microarray") ~ paste(tissue, "transcripts"),
    str_detect(dataset, "Proteomics") ~ paste(tissue, "proteins"),
    TRUE ~ "Other circulating features"
  )
}

sample_sheet <- read_csv(sample_file, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    SampleID = tolower(as.character(SampleID)),
    Role = as.character(Role),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort)
  ) %>%
  filter(AnalysisSet == "Primary", Is_Primary_Cohort) %>%
  distinct(Dataset, Timepoint, Clinical_Subject_ID, .keep_all = TRUE)

associations <- read_csv(feature_file, show_col_types = FALSE, progress = FALSE) %>%
  filter(
    Model %in% c("Pre_Pre", "Post_Post", "Delta_Delta"),
    ClinicalMarker %in% c("TotalOC", "cOC", "cOC_ratio"),
    BH_FDR < fdr_node_threshold,
    FeatureID != "cg05101930"
  ) %>%
  mutate(
    State = recode(Model, Pre_Pre = "Pre", Post_Post = "Post", Delta_Delta = "Delta"),
    ClinicalMarker = recode(
      ClinicalMarker,
      TotalOC = "tOC", cOC_ratio = "cOC/tOC ratio", .default = ClinicalMarker
    ),
    TissueClass = case_when(
      str_starts(Dataset, "Serum") ~ "Serum",
      str_starts(Dataset, "Adipose") ~ "Adipose",
      TRUE ~ "Muscle"
    ),
    Direction = case_when(
      is.na(Effect) & is.na(Effect_MValue) ~ "Unavailable",
      coalesce(Effect, Effect_MValue) >= 0 ~ "Positive",
      TRUE ~ "Negative"
    ),
    AssociationEffect = coalesce(Effect, Effect_MValue),
    Display = case_when(
      !is.na(GeneSymbol) & nzchar(GeneSymbol) ~ str_to_upper(GeneSymbol),
      FeatureID %in% names(feature_display_labels) ~ unname(feature_display_labels[FeatureID]),
      TRUE ~ FeatureID
    ),
    NodeID = paste(State, Timepoint, Dataset, FeatureID, sep = "|"),
    FeatureClass = classify_feature(FeatureID, Display, Dataset, TissueClass)
  )

node_table <- associations %>%
  group_by(NodeID, State, Timepoint, Dataset, TissueClass, FeatureID, Display, FeatureClass) %>%
  summarise(
    MinimumOCFDR = min(BH_FDR, na.rm = TRUE),
    MarkerDegree = n_distinct(ClinicalMarker),
    Markers = paste(marker_order[marker_order %in% unique(ClinicalMarker)], collapse = " + "),
    .groups = "drop"
  )

extract_node_profiles <- function(dataset, nodes) {
  input_file <- matrix_registry$InputFile[match(dataset, matrix_registry$Dataset)]
  if (is.na(input_file)) stop("No analysis matrix registered for ", dataset)
  matrix <- read_analysis_matrix(path_analysis_ready(input_file))
  profiles <- vector("list", nrow(nodes))

  for (i in seq_len(nrow(nodes))) {
    node <- nodes[i, ]
    if (!node$FeatureID %in% rownames(matrix)) next
    target_time <- as.character(node$Timepoint)

    if (node$State == "Delta") {
      metadata <- sample_sheet %>%
        filter(
          Dataset == dataset, Timepoint %in% c("Pre", target_time),
          SampleID %in% colnames(matrix)
        ) %>%
        select(FamilyID, Clinical_Subject_ID, Role, Timepoint, SampleID) %>%
        pivot_wider(names_from = Timepoint, values_from = SampleID) %>%
        filter(!is.na(Pre), !is.na(.data[[target_time]]))
      values <- matrix[node$FeatureID, metadata[[target_time]]] - matrix[node$FeatureID, metadata$Pre]
    } else {
      metadata <- sample_sheet %>%
        filter(
          Dataset == dataset, Timepoint == target_time,
          SampleID %in% colnames(matrix)
        ) %>%
        select(FamilyID, Clinical_Subject_ID, Role, SampleID)
      values <- matrix[node$FeatureID, metadata$SampleID]
    }

    frame <- metadata %>%
      transmute(
        NodeID = node$NodeID,
        Clinical_Subject_ID, FamilyID, Role,
        Value = as.numeric(values)
      ) %>%
      filter(is.finite(Value), !is.na(Role), !is.na(FamilyID))

    if (nrow(frame) >= 8L && sd(frame$Value) > 0) {
      ranked <- rank(frame$Value, ties.method = "average")
      design <- model.matrix(~ factor(Role), data = frame)
      frame$Residual <- residuals(lm.fit(design, ranked))
      profiles[[i]] <- frame
    }
  }
  rm(matrix); gc(verbose = FALSE)
  bind_rows(profiles)
}

profile_data <- map_dfr(
  unique(node_table$Dataset),
  ~ extract_node_profiles(.x, filter(node_table, Dataset == .x))
)

pairwise_state_correlations <- function(state) {
  state_nodes <- node_table %>% filter(State == state, NodeID %in% unique(profile_data$NodeID))
  if (nrow(state_nodes) < 2L) return(tibble())
  pairs <- t(combn(state_nodes$NodeID, 2L))

  results <- lapply(seq_len(nrow(pairs)), function(i) {
    a <- profile_data %>%
      filter(NodeID == pairs[i, 1]) %>%
      select(Clinical_Subject_ID, FamilyID, ResidualA = Residual)
    b <- profile_data %>%
      filter(NodeID == pairs[i, 2]) %>%
      select(Clinical_Subject_ID, ResidualB = Residual)
    common <- inner_join(a, b, by = "Clinical_Subject_ID") %>%
      filter(is.finite(ResidualA), is.finite(ResidualB))
    n <- nrow(common)
    if (n < minimum_pair_n || sd(common$ResidualA) == 0 || sd(common$ResidualB) == 0) {
      return(tibble(
        State = state, From = pairs[i, 1], To = pairs[i, 2], N = n,
        Rho = NA_real_, PValue = NA_real_, LOFOSignStability = NA_real_
      ))
    }
    test <- suppressWarnings(cor.test(
      common$ResidualA, common$ResidualB,
      method = "spearman", exact = FALSE
    ))
    rho <- unname(test$estimate)
    families <- unique(common$FamilyID)
    lofo <- vapply(families, function(family_id) {
      keep <- common$FamilyID != family_id
      if (sum(keep) < minimum_pair_n - 2L) return(NA_real_)
      suppressWarnings(cor(
        common$ResidualA[keep], common$ResidualB[keep],
        method = "spearman"
      ))
    }, numeric(1))
    lofo <- lofo[is.finite(lofo)]
    tibble(
      State = state, From = pairs[i, 1], To = pairs[i, 2], N = n,
      Rho = rho, PValue = test$p.value,
      LOFOSignStability = if (length(lofo)) mean(sign(lofo) == sign(rho)) else NA_real_
    )
  })

  bind_rows(results) %>%
    mutate(
      BH_FDR = p.adjust(PValue, method = "BH"),
      Direction = if_else(Rho >= 0, "Positive", "Negative"),
      Qualifying = is.finite(BH_FDR) & BH_FDR < fdr_edge_threshold &
        abs(Rho) >= absolute_rho_threshold &
        LOFOSignStability >= minimum_lofo_stability
    )
}

correlation_results <- map_dfr(state_order, pairwise_state_correlations)
qualifying_edges <- correlation_results %>% filter(Qualifying)

select_display_edges <- function(edges, max_edges = 75L) {
  if (!nrow(edges)) return(edges)
  ranked <- edges %>% mutate(EdgeID = row_number(), Strength = abs(Rho))
  by_from <- ranked %>% group_by(From) %>% slice_max(Strength, n = 2, with_ties = FALSE) %>% ungroup()
  by_to <- ranked %>% group_by(To) %>% slice_max(Strength, n = 2, with_ties = FALSE) %>% ungroup()
  bind_rows(by_from, by_to, ranked %>% slice_max(Strength, n = max_edges, with_ties = FALSE)) %>%
    distinct(EdgeID, .keep_all = TRUE) %>%
    arrange(BH_FDR, desc(Strength)) %>%
    slice_head(n = max_edges) %>%
    select(-EdgeID, -Strength)
}

build_state_plot <- function(state) {
  edge_all <- qualifying_edges %>% filter(State == state)
  if (!nrow(edge_all)) {
    return(ggplot() +
      annotate(
        "text", x = 0.5, y = 0.5,
        label = "No stable molecular correlation edge\npassed BH-FDR < 0.10",
        family = "Arial", size = 2.4
      ) +
      labs(title = state_labels[state]) +
      theme_void(base_family = "Arial") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 7.2)))
  }

  connected_ids <- union(edge_all$From, edge_all$To)
  shared_core_ids <- node_table %>%
    filter(State == state, State == "Post", Display == "C4B") %>%
    pull(NodeID)
  nodes <- node_table %>%
    filter(State == state, NodeID %in% union(connected_ids, shared_core_ids))
  edge_display <- select_display_edges(edge_all, max_edges = if_else(state == "Post", 95L, 70L))
  graph <- graph_from_data_frame(
    edge_all %>% transmute(from = From, to = To, Rho, BH_FDR, weight = abs(Rho)),
    directed = FALSE,
    vertices = nodes %>% rename(name = NodeID)
  )
  community <- cluster_louvain(graph, weights = E(graph)$weight)
  community_map <- tibble(
    NodeID = names(membership(community)),
    Module = as.integer(membership(community))
  )

  graph_display <- graph_from_data_frame(
    edge_display %>% transmute(from = From, to = To, Rho, Direction, BH_FDR, weight = abs(Rho)),
    directed = FALSE,
    vertices = nodes %>%
      left_join(community_map, by = "NodeID") %>%
      rename(name = NodeID)
  )
  set.seed(20260831 + match(state, state_order))
  layout <- create_layout(graph_display, layout = "fr", weights = E(graph_display)$weight, niter = 4000)
  layout <- layout %>%
    mutate(
      x = rescale(x, to = c(-1.52, 1.52)),
      y = rescale(y, to = c(-0.95, 0.95)),
      Degree = degree(graph_display)[name],
      IsUnclusteredCore = Degree == 0,
      x = if_else(IsUnclusteredCore, -1.28, x),
      y = if_else(IsUnclusteredCore, 0.88, y),
      LabelPriority = 100 * as.integer(MinimumOCFDR < 0.05) +
        12 * MarkerDegree + Degree,
      LabelNode = case_when(
        State == "Pre" ~ Display %in% c(
          "C4B", "AK9", "APOF", "DEFA10P"
        ),
        State == "Post" ~ Display %in% c(
          "C4B", "Omega-6/total FA", "PUFA/total FA", "M_VLDL_P", "PALLD"
        ),
        State == "Delta" ~ Display %in% c(
          "IDL phospholipids", "ADAM23", "LRRIQ1", "BCHE", "CMIP"
        ),
        TRUE ~ FALSE
      ),
      Label = case_when(
        Display == "C4B" & State == "Post" ~ "C4B (3 h)",
        LabelNode ~ Display,
        TRUE ~ ""
      )
    )

  module_composition <- as_tibble(layout) %>%
    count(Module, FeatureClass, name = "ClassN") %>%
    group_by(Module) %>%
    mutate(ModuleN = sum(ClassN), ClassFraction = ClassN / ModuleN) %>%
    arrange(desc(ClassN), FeatureClass, .by_group = TRUE) %>%
    mutate(ClassRank = row_number()) %>%
    ungroup()

  module_labels <- module_composition %>%
    filter(ClassRank <= 2L) %>%
    group_by(Module) %>%
    summarise(
      PrimaryClass = first(FeatureClass),
      SecondaryClass = if_else(
        n() >= 2L && nth(ClassN, 2L) >= 2L && nth(ClassFraction, 2L) >= 0.10,
        nth(FeatureClass, 2L), NA_character_
      ),
      ModuleN = first(ModuleN),
      .groups = "drop"
    ) %>%
    mutate(
      ModuleLabel = case_when(
        PrimaryClass == "Lipoprotein remodeling" &
          SecondaryClass == "Fatty-acid composition" ~
          "Lipid-remodeling module",
        PrimaryClass == "Lipoprotein remodeling" &
          SecondaryClass == "Serum proteins" ~
          "Lipoprotein and serum-protein module",
        PrimaryClass == "Lipoprotein remodeling" &
          SecondaryClass == "Complement and acute-phase proteins" ~
          "Lipoprotein and complement module",
        PrimaryClass == "Lipoprotein remodeling" &
          SecondaryClass == "Central carbon metabolism" ~
          "Lipoprotein and central-carbon module",
        PrimaryClass == "Lipoprotein remodeling" &
          is.na(SecondaryClass) ~ "Lipoprotein module",
        PrimaryClass == "Fatty-acid composition" &
          SecondaryClass == "Lipoprotein remodeling" ~
          "Fatty-acid module",
        PrimaryClass == "Serum proteins" &
          SecondaryClass == "Adipose transcripts" ~
          "Cross-tissue molecular module",
        PrimaryClass == "Serum proteins" & is.na(SecondaryClass) ~
          "Serum-protein module",
        !is.na(SecondaryClass) ~ paste0(PrimaryClass, " and ", SecondaryClass),
        TRUE ~ PrimaryClass
      )
    ) %>%
    group_by(ModuleLabel) %>%
    mutate(
      DuplicateIndex = row_number(),
      DuplicateN = n(),
      ModuleLabel = if_else(
        DuplicateN > 1L,
        paste0(ModuleLabel, " ", as.roman(DuplicateIndex)),
        ModuleLabel
      )
    ) %>%
    ungroup()

  module_audit[[state]] <<- module_composition %>%
    left_join(
      module_labels %>% select(Module, ModuleLabel, PrimaryClass, SecondaryClass),
      by = "Module"
    ) %>%
    mutate(State = state_labels[state], .before = 1)
  layout <- layout %>% left_join(module_labels, by = "Module")
  hull_data <- layout %>%
    group_by(Module) %>% mutate(ModuleN = n()) %>% ungroup() %>%
    filter(ModuleN >= 3L)
  module_label_positions <- hull_data %>%
    group_by(Module, ModuleLabel) %>%
    summarise(
      LabelX = pmin(pmax(median(x), -1.28), 1.28),
      LabelY = as.numeric(quantile(y, probs = 0.18, names = FALSE)) + 0.04,
      .groups = "drop"
    ) %>%
    mutate(
      LabelX = case_when(
        state == "Pre" & ModuleLabel == "Cross-tissue molecular module" ~ -1.12,
        state == "Pre" & str_detect(ModuleLabel, "Fatty-acid") ~ -0.20,
        state == "Pre" & ModuleLabel == "Lipoprotein module" ~ 0.72,
        state == "Post" & ModuleLabel == "Lipid-remodeling module" ~ -0.82,
        state == "Post" & ModuleLabel == "Lipoprotein module" ~ 0.62,
        state == "Delta" & ModuleLabel == "Serum-protein module" ~ -0.72,
        state == "Delta" & ModuleLabel == "Lipoprotein module" ~ 0.68,
        TRUE ~ LabelX
      ),
      LabelY = case_when(
        state == "Pre" & ModuleLabel == "Cross-tissue molecular module" ~ -0.08,
        state == "Pre" & str_detect(ModuleLabel, "Fatty-acid") ~ -0.92,
        state == "Pre" & ModuleLabel == "Lipoprotein module" ~ -0.68,
        state == "Post" & ModuleLabel == "Lipid-remodeling module" ~ -0.62,
        state == "Post" & ModuleLabel == "Lipoprotein module" ~ -0.58,
        state == "Delta" & ModuleLabel == "Serum-protein module" ~ -0.58,
        state == "Delta" & ModuleLabel == "Lipoprotein module" ~ -0.58,
        TRUE ~ LabelY
      )
    )

  marker_data <- associations %>%
    filter(State == state, NodeID %in% layout$name) %>%
    distinct(NodeID, ClinicalMarker, Direction, BH_FDR) %>%
    mutate(
      MarkerX = case_when(
        ClinicalMarker == "cOC" ~ 0,
        ClinicalMarker == "tOC" ~ -1.02,
        ClinicalMarker == "cOC/tOC ratio" ~ 1.02,
        TRUE ~ 0
      ),
      MarkerY = case_when(
        ClinicalMarker == "cOC" ~ 1.48,
        ClinicalMarker %in% c("tOC", "cOC/tOC ratio") ~ -1.42,
        TRUE ~ -1.42
      ),
      State = state_labels[state]
    ) %>%
    left_join(layout %>% as_tibble() %>% select(NodeID = name, NodeX = x, NodeY = y), by = "NodeID")
  active_markers <- marker_data %>% distinct(ClinicalMarker, MarkerX, MarkerY)

  ggraph(layout) +
    ggforce::geom_mark_hull(
      data = hull_data,
      aes(x = x, y = y, group = Module),
      inherit.aes = FALSE,
      fill = "#DCE5E8", color = "#707070", alpha = 0.16,
      linewidth = 0.30, linetype = "22",
      expand = grid::unit(1.2, "mm"),
      show.legend = FALSE
    ) +
    geom_edge_link(
      aes(edge_width = abs(Rho)),
      colour = "#858585", alpha = 0.32,
      lineend = "round", show.legend = FALSE
    ) +
    scale_edge_width(range = c(0.18, 0.65), guide = "none") +
    geom_segment(
      data = marker_data,
      aes(
        x = MarkerX, y = MarkerY, xend = NodeX, yend = NodeY,
        colour = Direction,
        linewidth = if_else(BH_FDR < 0.05, 0.28, 0.12),
        alpha = if_else(BH_FDR < 0.05, 0.46, 0.10)
      ),
      inherit.aes = FALSE, lineend = "round"
    ) +
    scale_colour_manual(values = association_colors) +
    scale_linewidth_identity() + scale_alpha_identity() +
    geom_node_point(
      aes(fill = TissueClass), shape = 21,
      size = ifelse(layout$MinimumOCFDR < 0.05, 1.55, 1.05),
      stroke = ifelse(layout$MinimumOCFDR < 0.05, 0.42, 0.16),
      colour = "black"
    ) +
    scale_fill_manual(values = tissue_colors) +
    geom_node_text(
      aes(label = str_wrap(Label, 15)), repel = TRUE,
      size = 1.62, family = "Arial", colour = "black",
      box.padding = 0.14, point.padding = 0.05,
      segment.size = 0.14, max.overlaps = Inf
    ) +
    geom_label(
      data = module_label_positions,
      aes(x = LabelX, y = LabelY, label = str_wrap(ModuleLabel, width = 15)),
      inherit.aes = FALSE, family = "Arial", fontface = "bold",
      size = 1.58, lineheight = 0.92, colour = "black",
      fill = alpha("white", 0.82), label.size = 0,
      label.padding = grid::unit(0.35, "mm")
    ) +
    geom_point(
      data = active_markers,
      aes(x = MarkerX, y = MarkerY), inherit.aes = FALSE,
      shape = 21, size = 5.0, stroke = 0.32,
      fill = "#F4F4F4", colour = "black"
    ) +
    geom_text(
      data = active_markers,
      aes(x = MarkerX, y = MarkerY, label = recode(ClinicalMarker, `cOC/tOC ratio` = "cOC/tOC")),
      inherit.aes = FALSE, family = "Arial", fontface = "bold",
      size = 1.55, colour = "black"
    ) +
    coord_cartesian(xlim = c(-1.82, 1.82), ylim = c(-1.68, 1.70), clip = "off") +
    labs(title = state_labels[state]) +
    theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 7.2, colour = "black"),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2)
    )
}

state_plots <- map(state_order, build_state_plot)

bar_exact <- associations %>%
  filter(
    (State == "Pre" & Display %in% c("C4B", "AK9", "APOF", "DEFA10P")) |
      (State == "Post" & Display %in% c("C4B", "PALLD")) |
      (State == "Delta" & Display %in% c("IDL phospholipids", "BCHE"))
  ) %>%
  group_by(State, Display, ClinicalMarker) %>%
  slice_min(BH_FDR, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(
    !(State == "Pre" & Display == "C4B" & ClinicalMarker != "tOC"),
    !(State == "Pre" & Display %in% c("AK9", "APOF") & ClinicalMarker != "cOC/tOC ratio"),
    !(State == "Pre" & Display == "DEFA10P" & ClinicalMarker != "cOC/tOC ratio"),
    !(State == "Post" & Display == "PALLD" & ClinicalMarker != "cOC/tOC ratio"),
    !(State == "Delta" & ClinicalMarker != "cOC")
  )

bar_post_lipids <- associations %>%
  filter(
    State == "Post",
    Display %in% c("M_VLDL_P", "Omega-6/total FA", "PUFA/total FA")
  ) %>%
  group_by(State, Display) %>%
  slice_min(BH_FDR, n = 1, with_ties = FALSE) %>%
  ungroup()

bar_data <- bind_rows(bar_exact, bar_post_lipids) %>%
  mutate(
    Display = recode(Display, M_VLDL_P = "M-VLDL-P"),
    BarDisplay = case_when(
      Display == "DEFA10P" & str_detect(Dataset, "Methylation") ~ "DEFA10P-DNAm",
      Display == "PALLD" & str_detect(Dataset, "Methylation") ~ "PALLD-DNAm",
      str_detect(Dataset, "Methylation") ~ paste0(Display, "-DNAm"),
      TRUE ~ Display
    ),
    FeatureLabel = paste0(
      BarDisplay, " (", recode(ClinicalMarker, `cOC/tOC ratio` = "ratio"), ")"
    ),
    StateLabel = factor(
      State, levels = state_order,
      labels = unname(state_labels[state_order])
    ),
    EvidenceStrength = -log10(BH_FDR),
    FDRLabel = sprintf("%.3f", BH_FDR),
    FeatureLabel = factor(FeatureLabel, levels = rev(unique(FeatureLabel)))
  )

bar_plot <- ggplot(bar_data, aes(x = EvidenceStrength, y = FeatureLabel, fill = Direction)) +
  geom_col(width = 0.58, colour = "black", linewidth = 0.18) +
  geom_text(
    aes(label = FDRLabel), hjust = -0.08,
    family = "Arial", size = 1.72, colour = "black"
  ) +
  facet_grid(StateLabel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = association_colors, guide = "none") +
  scale_x_continuous(
    limits = c(0, max(bar_data$EvidenceStrength) * 1.34),
    breaks = c(0, 1, 2), expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(position = "right") +
  labs(
    title = "Core molecular evidence",
    x = expression(-log[10]("BH-FDR")), y = NULL
  ) +
  theme_classic(base_family = "Arial", base_size = 6.2) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 7.2, colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.30),
    axis.line = element_blank(),
    axis.title.x = element_text(face = "bold", size = 5.8, colour = "black"),
    axis.text = element_text(size = 4.8, colour = "black"),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0, face = "bold", size = 4.8,
      colour = "black", hjust = 1
    ),
    panel.spacing.y = grid::unit(0.45, "mm"),
    plot.margin = margin(12.2, 1, 12.2, 1, unit = "mm")
  )

legend_plot <- ggplot() +
  annotate("text", x = 0.08, y = 0.90, label = "Tissue", hjust = 0,
           family = "Arial", fontface = "bold", size = 2.15) +
  annotate("point", x = 0.13, y = c(0.80, 0.70, 0.60), shape = 21, size = 2.8,
           fill = unname(tissue_colors), colour = "black", stroke = 0.25) +
  annotate("text", x = 0.28, y = c(0.80, 0.70, 0.60),
           label = names(tissue_colors), hjust = 0, family = "Arial", size = 2.0) +
  annotate("text", x = 0.08, y = 0.48, label = "Edges", hjust = 0,
           family = "Arial", fontface = "bold", size = 2.15) +
  annotate("segment", x = 0.08, xend = 0.27, y = 0.39, yend = 0.39,
           colour = "#858585", linewidth = 0.8) +
  annotate("text", x = 0.30, y = 0.39, label = "Stable molecular\ncorrelation",
           hjust = 0, vjust = 0.5, family = "Arial", size = 1.9, lineheight = 0.92) +
  annotate("segment", x = 0.08, xend = 0.27, y = c(0.27, 0.17), yend = c(0.27, 0.17),
           colour = unname(association_colors), linewidth = 0.8) +
  annotate("text", x = 0.30, y = c(0.27, 0.17),
           label = c("Positive OC\nassociation", "Negative OC\nassociation"),
           hjust = 0, vjust = 0.5, family = "Arial", size = 1.9, lineheight = 0.92) +
  annotate("point", x = 0.13, y = 0.06, shape = 21, size = 3.1,
           fill = "white", colour = "black", stroke = 0.7) +
  annotate("text", x = 0.28, y = 0.06, label = "OC association\nBH-FDR < 0.05",
           hjust = 0, vjust = 0.5, family = "Arial", size = 1.9, lineheight = 0.92) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(5, 1, 5, 1, unit = "mm"))

combined <- wrap_plots(state_plots, nrow = 1) / legend_plot +
  plot_layout(heights = c(1, 0.12))

combined_four_panel <-
  state_plots[[1]] + state_plots[[2]] + state_plots[[3]] + bar_plot + legend_plot +
  plot_layout(widths = c(1.03, 1.03, 1.03, 0.90, 0.72))

ggsave(
  file.path(output_dir, "Figure_6D_Osteocalcin_Modular_Molecular_Network.pdf"),
  combined, width = 190, height = 78, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Figure_6D_Osteocalcin_Modular_Network_with_FDR_Bars_190mm.pdf"),
  combined_four_panel, width = 190, height = 72, units = "mm", device = cairo_pdf
)

write_csv(node_table, file.path(output_dir, "Fig06D_OC_Modular_Network_Node_Audit.csv"))
write_csv(
  bind_rows(module_audit),
  file.path(output_dir, "Fig06D_OC_Modular_Network_Module_Composition.csv")
)
write_csv(correlation_results, file.path(output_dir, "Fig06D_OC_Molecular_Correlation_All_Tests.csv"))
write_csv(qualifying_edges, file.path(output_dir, "Fig06D_OC_Molecular_Correlation_Qualifying_Edges.csv"))
write_csv(
  bar_data %>%
    transmute(
      State, Timepoint, Feature = Display, ClinicalMarker,
      Direction, Effect = AssociationEffect, P_Value, BH_FDR
    ),
  file.path(output_dir, "Fig06D_OC_Core_Molecular_FDR_Bar_Data.csv")
)
write_csv(
  correlation_results %>%
    group_by(State) %>%
    summarise(
      TestedPairs = sum(is.finite(PValue)),
      NominalP05 = sum(PValue < 0.05, na.rm = TRUE),
      BH_FDR05 = sum(BH_FDR < 0.05, na.rm = TRUE),
      BH_FDR10 = sum(BH_FDR < 0.10, na.rm = TRUE),
      QualifyingStableEdges = sum(Qualifying, na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(output_dir, "Fig06D_OC_Modular_Network_Edge_Counts.csv")
)
writeLines(
  c(
    "# Figure 6D modular molecular network",
    "",
    "Nodes were preselected by an osteocalcin-feature association at BH-FDR < 0.10.",
    "Molecular edges are Role-residualized Spearman correlations tested within each state.",
    paste0("Displayed molecular edges require BH-FDR < ", fdr_edge_threshold,
           ", absolute rho >= ", absolute_rho_threshold,
           ", pairwise N >= ", minimum_pair_n,
           ", and leave-one-family-out sign stability >= ", minimum_lofo_stability, "."),
    "Communities were detected with the Louvain algorithm using absolute rho as the edge weight.",
    "Module labels report the two most abundant annotation classes among connected nodes and are descriptive.",
    "OC-feature edges retain the formal Role-adjusted, FamilyID-blocked association models.",
    "Network positions are layout artifacts and do not encode biological distance or causality."
  ),
  file.path(output_dir, "README.md"), useBytes = TRUE
)

capture.output(sessionInfo(), file = file.path(output_dir, "R_SESSION_INFO.txt"))
message("Modular network written to: ", output_dir)
