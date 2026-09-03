# Unified pre-exercise, Post3h absolute, and Post3h-minus-Pre association audit.
# Raw/source data are read only. Outputs are written to a new dated result folder.

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

p_load(digest, dplyr, fgsea, limma, matrixStats, msigdbr, readr, stringr,
       tibble, tidyr, writexl)

analysis_id <- "Unified_Pre_Post_Response_Association_Audit_2026-08-28"
output_dir <- path_results(analysis_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
minimum_subjects <- 10L
maximum_correlation_features <- 20000L
minimum_gene_set_size <- 10L
maximum_gene_set_size <- 500L
role_levels <- c("Daughter", "Mother", "Father")
set.seed(20260828)

rank_normalize <- function(x) {
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  n <- sum(keep)
  if (n >= 3L) out[keep] <- qnorm((rank(x[keep], ties.method = "average") - 0.5) / n)
  out
}

clinical_file <- path_analysis_ready("Metadata", "Clinical_Observed_Primary.csv")
sample_sheet_file <- path_analysis_ready("Metadata", "Project_Sample_Sheet.csv")

clinical <- read_csv(clinical_file, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = factor(Role, levels = role_levels),
    BodyFat = suppressWarnings(as.numeric(Fat_p_G)),
    HOMA_IR = suppressWarnings(as.numeric(HOMA_IR_84m_G)),
    Matsuda_ISI = suppressWarnings(as.numeric(ISI_G)),
    hsCRP = suppressWarnings(as.numeric(CRP_84m_G)),
    TotalOC = suppressWarnings(as.numeric(S_TotalOC_84m_G)),
    cOC = suppressWarnings(as.numeric(S_cOC_84m_G)),
    cOC_ratio = suppressWarnings(as.numeric(cOC_totalOC_ratio_84m_G))
  ) %>%
  distinct(Clinical_Subject_ID, .keep_all = TRUE) %>%
  mutate(
    AdiposityGroup = factor(if_else(BodyFat > 30, "Obese", "Lean"), levels = c("Lean", "Obese")),
    BodyFat_Z = rank_normalize(BodyFat),
    HOMA_IR_Z = rank_normalize(if_else(HOMA_IR > 0, HOMA_IR, NA_real_)),
    Matsuda_ISI_Z = rank_normalize(if_else(Matsuda_ISI > 0, Matsuda_ISI, NA_real_)),
    hsCRP_Z = rank_normalize(if_else(hsCRP > 0, hsCRP, NA_real_)),
    TotalOC_Z = rank_normalize(if_else(TotalOC > 0, TotalOC, NA_real_)),
    cOC_Z = rank_normalize(if_else(cOC > 0, cOC, NA_real_)),
    cOC_ratio_Z = rank_normalize(if_else(cOC_ratio > 0 & cOC_ratio < 1, cOC_ratio, NA_real_))
  )

if (nrow(clinical) != 81L || n_distinct(clinical$FamilyID) != 27L) {
  stop("Clinical input does not match the locked 81-participant, 27-family cohort.")
}

sample_sheet <- read_csv(sample_sheet_file, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    SampleID = tolower(as.character(SampleID)),
    Role = factor(Role, levels = role_levels),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort)
  ) %>%
  filter(AnalysisSet == "Primary", Is_Primary_Cohort)

predictors <- tribble(
  ~PredictorFamily, ~Predictor, ~PredictorLabel, ~PredictorColumn, ~PredictorType, ~EffectDefinition,
  "Adiposity", "AdiposityGroup", "Adiposity group", "AdiposityGroup", "Categorical", "Obese minus Lean",
  "Adiposity", "BodyFat", "Body-fat percentage", "BodyFat_Z", "Continuous", "Per rank-normalized SD",
  "Insulin", "HOMA_IR", "HOMA-IR", "HOMA_IR_Z", "Continuous", "Per rank-normalized SD",
  "Insulin", "Matsuda_ISI", "Matsuda ISI", "Matsuda_ISI_Z", "Continuous", "Per rank-normalized SD",
  "Inflammation", "hsCRP", "hsCRP", "hsCRP_Z", "Continuous", "Per rank-normalized SD",
  "Osteocalcin", "TotalOC", "Total osteocalcin", "TotalOC_Z", "Continuous", "Per rank-normalized SD",
  "Osteocalcin", "cOC", "Carboxylated osteocalcin", "cOC_Z", "Continuous", "Per rank-normalized SD",
  "Osteocalcin", "cOC_ratio", "cOC-to-total osteocalcin ratio", "cOC_ratio_Z", "Continuous", "Per rank-normalized SD"
)

datasets <- tribble(
  ~Dataset, ~Layer, ~Tissue, ~InputFile, ~PathwayEligible, ~AvailableStates,
  "Serum_Metabonomics", "Metabolomics", "Serum", "Analysis_Serum_Metabonomics.csv", FALSE, "PreExercise;Post3hAbsolute;ExerciseChange",
  "Serum_Proteomics", "Proteomics", "Serum", "Analysis_Serum_Proteomics.csv", TRUE, "PreExercise;Post3hAbsolute;ExerciseChange",
  "Adipose_Microarray", "Transcriptomics", "Adipose", "Analysis_Adipose_Microarray.csv", TRUE, "PreExercise",
  "Adipose_Proteomics", "Proteomics", "Adipose", "Analysis_Adipose_Proteomics.csv", TRUE, "PreExercise;Post3hAbsolute;ExerciseChange",
  "Muscle_Microarray", "Transcriptomics", "Skeletal muscle", "Analysis_Muscle_Microarray.csv", TRUE, "PreExercise",
  "Muscle_Proteomics", "Proteomics", "Skeletal muscle", "Analysis_Muscle_Proteomics.csv", TRUE, "PreExercise;Post3hAbsolute;ExerciseChange"
)

states <- tribble(
  ~State, ~StateLabel, ~Contrast,
  "PreExercise", "Pre-exercise molecular profile", "Pre",
  "Post3hAbsolute", "Post-exercise molecular profile", "Post3h",
  "ExerciseChange", "Exercise-induced molecular change", "Post3h minus Pre"
)

read_annotation <- function(dataset) {
  metadata_dir <- path_analysis_ready("Metadata")
  if (dataset == "Serum_Proteomics") {
    return(read_csv(file.path(metadata_dir, "Feature_Dictionary_Serum_Proteomics.csv"),
                    show_col_types = FALSE, progress = FALSE) %>%
             transmute(FeatureID = as.character(FeatureID), GeneSymbol = as.character(GeneSymbol)) %>%
             distinct(FeatureID, .keep_all = TRUE))
  }
  if (grepl("Microarray$", dataset)) {
    return(read_csv(file.path(metadata_dir, paste0("Feature_Annotation_", dataset, ".csv")),
                    show_col_types = FALSE, progress = FALSE) %>%
             transmute(FeatureID = as.character(ProbeSetID), GeneSymbol = as.character(SYMBOL)) %>%
             distinct(FeatureID, .keep_all = TRUE))
  }
  read_csv(file.path(metadata_dir, paste0("Feature_Annotation_", dataset, ".csv")),
           show_col_types = FALSE, progress = FALSE) %>%
    transmute(
      FeatureID = as.character(FeatureID),
      GeneSymbol = case_when(
        "IntegrationGene" %in% names(.) & !is.na(IntegrationGene) & nzchar(IntegrationGene) ~ IntegrationGene,
        "PG.Genes" %in% names(.) & !is.na(PG.Genes) & nzchar(PG.Genes) ~ PG.Genes,
        TRUE ~ NA_character_
      )
    ) %>% distinct(FeatureID, .keep_all = TRUE)
}

estimate_family_correlation <- function(values, design, family_id) {
  use <- values
  if (nrow(use) > maximum_correlation_features) {
    variances <- matrixStats::rowVars(use, na.rm = TRUE)
    use <- use[head(order(variances, decreasing = TRUE, na.last = NA), maximum_correlation_features), , drop = FALSE]
  }
  value <- suppressWarnings(duplicateCorrelation(use, design, block = family_id)$consensus.correlation)
  if (!is.finite(value)) value <- 0
  value
}

prepare_state <- function(dataset, state, matrix) {
  available <- sample_sheet %>%
    filter(Dataset == dataset, Timepoint %in% c("Pre", "Post3h"), SampleID %in% colnames(matrix)) %>%
    inner_join(clinical, by = c("FamilyID", "Clinical_Subject_ID", "Role")) %>%
    distinct(Clinical_Subject_ID, Timepoint, .keep_all = TRUE)
  if (state == "PreExercise") {
    metadata <- available %>% filter(Timepoint == "Pre") %>% arrange(FamilyID, Role, Clinical_Subject_ID)
    values <- matrix[, metadata$SampleID, drop = FALSE]
  } else {
    metadata <- available %>%
      select(FamilyID, Clinical_Subject_ID, Role, all_of(predictors$PredictorColumn), Timepoint, SampleID) %>%
      pivot_wider(names_from = Timepoint, values_from = SampleID) %>%
      filter(!is.na(Pre), !is.na(Post3h)) %>% arrange(FamilyID, Role, Clinical_Subject_ID)
    pre <- matrix[, metadata$Pre, drop = FALSE]
    post <- matrix[, metadata$Post3h, drop = FALSE]
    values <- if (state == "Post3hAbsolute") post else post - pre
  }
  colnames(values) <- metadata$Clinical_Subject_ID
  list(values = values, metadata = metadata)
}

fit_feature_model <- function(prepared, dataset_spec, state_spec, predictor_spec) {
  metadata <- prepared$metadata
  predictor_column <- predictor_spec$PredictorColumn
  if (predictor_spec$PredictorType == "Categorical") {
    complete <- !is.na(metadata[[predictor_column]])
  } else {
    complete <- is.finite(metadata[[predictor_column]])
  }
  complete <- complete & !is.na(metadata$Role) & !is.na(metadata$FamilyID)
  metadata <- metadata[complete, , drop = FALSE] %>% droplevels()
  values <- prepared$values[, complete, drop = FALSE]
  if (predictor_spec$PredictorType == "Categorical" && nlevels(metadata[[predictor_column]]) < 2L) return(NULL)
  if (nrow(metadata) < minimum_subjects || n_distinct(metadata$FamilyID) < 5L) return(NULL)

  observed <- rowSums(is.finite(values))
  variable <- matrixStats::rowSds(values, na.rm = TRUE) > 0
  keep <- observed >= minimum_subjects & is.finite(variable) & variable
  values <- values[keep, , drop = FALSE]
  observed <- observed[keep]
  design <- model.matrix(reformulate(c(predictor_column, "Role")), data = metadata)
  if (qr(design)$rank != ncol(design)) return(NULL)
  coefficient <- if (predictor_spec$PredictorType == "Categorical") {
    paste0(predictor_column, "Obese")
  } else predictor_column
  if (!coefficient %in% colnames(design)) return(NULL)
  family_correlation <- estimate_family_correlation(values, design, metadata$FamilyID)
  fit <- lmFit(values, design, block = metadata$FamilyID, correlation = family_correlation) %>%
    eBayes(trend = TRUE, robust = TRUE)
  stats <- topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
  tibble(
    AnalysisFramework = "Baseline phenotype association",
    Dataset = dataset_spec$Dataset, Layer = dataset_spec$Layer, Tissue = dataset_spec$Tissue,
    State = state_spec$State, StateLabel = state_spec$StateLabel, Contrast = state_spec$Contrast,
    PredictorFamily = predictor_spec$PredictorFamily, Predictor = predictor_spec$Predictor,
    PredictorLabel = predictor_spec$PredictorLabel, EffectDefinition = predictor_spec$EffectDefinition,
    FeatureID = rownames(stats), Effect = stats$logFC, ModeratedT = stats$t,
    P_Value = stats$P.Value, BH_FDR = p.adjust(stats$P.Value, method = "BH"),
    ObservedN = unname(observed[rownames(stats)]), Subjects = nrow(metadata),
    Families = n_distinct(metadata$FamilyID), DaughterN = sum(metadata$Role == "Daughter"),
    MotherN = sum(metadata$Role == "Mother"), FatherN = sum(metadata$Role == "Father"),
    FamilyCorrelation = family_correlation
  )
}

hallmark_raw <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  filter(!is.na(gs_name), !is.na(gene_symbol), nzchar(gene_symbol)) %>% distinct(gs_name, gene_symbol)
hallmark_sets <- lapply(split(hallmark_raw$gene_symbol, hallmark_raw$gs_name), unique)

run_gsea <- function(feature_results) {
  ranked <- feature_results %>%
    filter(!is.na(GeneSymbol), nzchar(GeneSymbol), is.finite(ModeratedT), is.finite(P_Value)) %>%
    mutate(GeneSymbol = str_split(GeneSymbol, ";", simplify = TRUE)[, 1]) %>%
    group_by(GeneSymbol) %>% arrange(P_Value, desc(ObservedN), FeatureID, .by_group = TRUE) %>%
    slice_head(n = 1L) %>% ungroup()
  stats <- ranked$ModeratedT
  names(stats) <- ranked$GeneSymbol
  stats <- sort(stats[is.finite(stats) & !duplicated(names(stats))], decreasing = TRUE)
  if (length(stats) < 50L) return(NULL)
  mapped_sets <- lapply(hallmark_sets, function(x) intersect(unique(x), names(stats)))
  eligible <- lengths(mapped_sets) >= minimum_gene_set_size & lengths(mapped_sets) <= maximum_gene_set_size
  mapped_sets <- mapped_sets[eligible]
  if (!length(mapped_sets)) return(NULL)
  suppressWarnings(fgseaMultilevel(pathways = mapped_sets, stats = stats,
                                   minSize = minimum_gene_set_size,
                                   maxSize = maximum_gene_set_size,
                                   eps = 0, nproc = 1)) %>%
    as_tibble() %>% transmute(
      Database = "MSigDB Hallmark", PathwayID = pathway,
      Pathway = pathway %>% str_remove("^HALLMARK_") %>% str_replace_all("_", " ") %>%
        str_to_lower() %>% str_to_sentence(),
      Size = size, NES = NES, P_Value = pval, BH_FDR = p.adjust(pval, method = "BH"),
      LeadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)),
      Method = "Ranked Hallmark GSEA"
    )
}

feature_list <- list()
pathway_list <- list()
availability <- list()
counter <- 0L
for (i in seq_len(nrow(datasets))) {
  dataset_spec <- datasets[i, ]
  input_path <- path_analysis_ready(dataset_spec$InputFile)
  message("Loading ", dataset_spec$Dataset)
  matrix <- read_analysis_matrix(input_path)
  annotation <- if (dataset_spec$PathwayEligible) read_annotation(dataset_spec$Dataset) else NULL
  state_names <- str_split(dataset_spec$AvailableStates, ";", simplify = TRUE)
  state_names <- state_names[state_names != ""]
  for (state_name in state_names) {
    state_spec <- states %>% filter(State == state_name)
    prepared <- prepare_state(dataset_spec$Dataset, state_name, matrix)
    availability[[paste(dataset_spec$Dataset, state_name, sep = "__")]] <- tibble(
      Dataset = dataset_spec$Dataset, Layer = dataset_spec$Layer, Tissue = dataset_spec$Tissue,
      State = state_name, MatrixSubjects = nrow(prepared$metadata),
      MatrixFamilies = n_distinct(prepared$metadata$FamilyID)
    )
    for (j in seq_len(nrow(predictors))) {
      predictor_spec <- predictors[j, ]
      message("Fitting ", dataset_spec$Dataset, " / ", state_name, " / ", predictor_spec$Predictor)
      result <- fit_feature_model(prepared, dataset_spec, state_spec, predictor_spec)
      if (is.null(result)) next
      if (dataset_spec$PathwayEligible) result <- result %>% left_join(annotation, by = "FeatureID")
      else result <- result %>% mutate(GeneSymbol = NA_character_)
      counter <- counter + 1L
      feature_list[[counter]] <- result
      if (dataset_spec$PathwayEligible) {
        pathway <- run_gsea(result)
        if (!is.null(pathway)) {
          pathway_list[[counter]] <- pathway %>% mutate(
            AnalysisFramework = "Baseline phenotype association",
            Dataset = dataset_spec$Dataset, Layer = dataset_spec$Layer, Tissue = dataset_spec$Tissue,
            State = state_spec$State, StateLabel = state_spec$StateLabel, Contrast = state_spec$Contrast,
            PredictorFamily = predictor_spec$PredictorFamily, Predictor = predictor_spec$Predictor,
            PredictorLabel = predictor_spec$PredictorLabel, Subjects = first(result$Subjects),
            Families = first(result$Families), .before = 1
          )
        }
      }
    }
    rm(prepared); invisible(gc())
  }
  rm(matrix, annotation); invisible(gc())
}

feature_all <- bind_rows(feature_list)
pathway_all <- bind_rows(pathway_list)

feature_summary <- feature_all %>%
  group_by(AnalysisFramework, PredictorFamily, Predictor, PredictorLabel, Dataset, Layer, Tissue,
           State, StateLabel, Contrast, EffectDefinition) %>%
  summarise(
    Subjects = first(Subjects), Families = first(Families), TestedFeatures = n(),
    NominalP05 = sum(P_Value < 0.05, na.rm = TRUE),
    BH_FDR05 = sum(BH_FDR < 0.05, na.rm = TRUE),
    BH_FDR10 = sum(BH_FDR < 0.10, na.rm = TRUE),
    BH_FDR20 = sum(BH_FDR < 0.20, na.rm = TRUE),
    MinimumP = min(P_Value, na.rm = TRUE), MinimumBH_FDR = min(BH_FDR, na.rm = TRUE),
    FamilyCorrelation = first(FamilyCorrelation), .groups = "drop"
  )

pathway_summary <- pathway_all %>%
  group_by(AnalysisFramework, PredictorFamily, Predictor, PredictorLabel, Dataset, Layer, Tissue,
           State, StateLabel, Contrast, Database, Method) %>%
  summarise(
    Subjects = first(Subjects), Families = first(Families), TestedPathways = n(),
    NominalP05 = sum(P_Value < 0.05, na.rm = TRUE),
    BH_FDR05 = sum(BH_FDR < 0.05, na.rm = TRUE),
    BH_FDR10 = sum(BH_FDR < 0.10, na.rm = TRUE),
    BH_FDR20 = sum(BH_FDR < 0.20, na.rm = TRUE),
    MinimumP = min(P_Value, na.rm = TRUE), MinimumBH_FDR = min(BH_FDR, na.rm = TRUE),
    .groups = "drop"
  )

feature_selected <- feature_all %>% filter(P_Value < 0.05 | BH_FDR < 0.20) %>%
  arrange(PredictorFamily, Predictor, Dataset, State, BH_FDR, P_Value)
pathway_selected <- pathway_all %>% filter(P_Value < 0.05 | BH_FDR < 0.20) %>%
  arrange(PredictorFamily, Predictor, Dataset, State, BH_FDR, P_Value)

write_csv(feature_summary, file.path(output_dir, "NonMethylation_Feature_Significance_Summary.csv"))
write_csv(pathway_summary, file.path(output_dir, "NonMethylation_Pathway_Significance_Summary.csv"))
write_csv(feature_selected, file.path(output_dir, "NonMethylation_Feature_Selected.csv"))
write_csv(pathway_selected, file.path(output_dir, "NonMethylation_Pathway_Selected.csv"))
write_csv(bind_rows(availability), file.path(output_dir, "NonMethylation_Sample_Availability.csv"))
write_csv(predictors, file.path(output_dir, "Predictor_Definitions.csv"))
write_xlsx(list(
  Feature_summary = feature_summary,
  Pathway_summary = pathway_summary,
  Feature_selected = feature_selected,
  Pathway_selected = pathway_selected,
  Sample_availability = bind_rows(availability),
  Predictor_definitions = predictors
), file.path(output_dir, "Unified_NonMethylation_Pre_Post_Response_Audit.xlsx"))
saveRDS(list(feature_summary = feature_summary, pathway_summary = pathway_summary,
             feature_results = feature_all, pathway_results = pathway_all,
             predictors = predictors, availability = bind_rows(availability)),
        file.path(output_dir, "Unified_NonMethylation_Pre_Post_Response_Audit.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_SESSION_INFO_NonMethylation.txt"))
message("Completed: ", output_dir)
