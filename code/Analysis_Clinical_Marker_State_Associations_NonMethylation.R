# Six clinical exercise markers across pre-pre, post-post and delta-delta
# omics models. This is a new analysis release and does not overwrite prior runs.

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

p_load(dplyr, fgsea, limma, matrixStats, msigdbr, readr, stringr, tibble, tidyr, writexl)

analysis_id <- "Six_Markers_Three_State_Omics_2026-08-22"
output_dir <- path_results(analysis_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
minimum_subjects <- 10L
minimum_gene_set_size <- 10L
maximum_gene_set_size <- 500L
role_levels <- c("Daughter", "Mother", "Father")
requested_markers <- c("cOC", "TotalOC", "cOC_ratio", "BCAA", "Serotonin", "Leptin")

rank_normalize <- function(x) {
  out <- rep(NA_real_, length(x)); keep <- is.finite(x); n <- sum(keep)
  if (n >= 3L) out[keep] <- qnorm((rank(x[keep], ties.method = "average") - 0.5) / n)
  out
}
safe_log2 <- function(x) ifelse(is.finite(x) & x > 0, log2(x), NA_real_)

clinical <- read_csv(
  file.path(path_intermediate("Clinical_Exercise_Source_Audit_2026-08-21"),
            "Restricted", "OC_BCAA_Serotonin_Long_Internal.csv"),
  show_col_types = FALSE, progress = FALSE
) %>%
  filter(MatchStatus == "Unique locked-cohort match", Marker %in% requested_markers) %>%
  mutate(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = factor(Role, levels = role_levels),
    Value = suppressWarnings(as.numeric(Value))
  ) %>%
  select(FamilyID, Clinical_Subject_ID, Role, Marker, MarkerLabel, Timepoint, Value) %>%
  distinct(FamilyID, Clinical_Subject_ID, Marker, Timepoint, .keep_all = TRUE) %>%
  pivot_wider(names_from = Timepoint, values_from = Value) %>%
  mutate(
    Pre_log2 = safe_log2(Pre), Post1h_log2 = safe_log2(Post1h), Post3h_log2 = safe_log2(Post3h),
    Delta1h = Post1h_log2 - Pre_log2, Delta3h = Post3h_log2 - Pre_log2
  )
marker_registry <- clinical %>% distinct(Marker, MarkerLabel) %>% arrange(match(Marker, requested_markers))

sample_sheet <- read_csv(
  path_analysis_ready("Metadata", "Project_Sample_Sheet.csv"),
  show_col_types = FALSE, progress = FALSE
) %>%
  mutate(
    FamilyID = as.character(FamilyID), Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    SampleID = tolower(as.character(SampleID)), Role = factor(Role, levels = role_levels),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort)
  ) %>% filter(AnalysisSet == "Primary", Is_Primary_Cohort)

registry <- tribble(
  ~Dataset, ~Layer, ~Tissue, ~InputFile, ~Timepoint, ~PathwayEligible,
  "Serum_Proteomics", "Proteomics", "Serum", "Analysis_Serum_Proteomics.csv", "Pre", TRUE,
  "Serum_Proteomics", "Proteomics", "Serum", "Analysis_Serum_Proteomics.csv", "Post1h", TRUE,
  "Serum_Proteomics", "Proteomics", "Serum", "Analysis_Serum_Proteomics.csv", "Post3h", TRUE,
  "Serum_Metabonomics", "Metabolomics", "Serum", "Analysis_Serum_Metabonomics.csv", "Pre", FALSE,
  "Serum_Metabonomics", "Metabolomics", "Serum", "Analysis_Serum_Metabonomics.csv", "Post1h", FALSE,
  "Serum_Metabonomics", "Metabolomics", "Serum", "Analysis_Serum_Metabonomics.csv", "Post3h", FALSE,
  "Adipose_Microarray", "Transcriptomics", "Adipose", "Analysis_Adipose_Microarray.csv", "Pre", TRUE,
  "Muscle_Microarray", "Transcriptomics", "Muscle", "Analysis_Muscle_Microarray.csv", "Pre", TRUE,
  "Adipose_Proteomics", "Proteomics", "Adipose", "Analysis_Adipose_Proteomics.csv", "Pre", TRUE,
  "Adipose_Proteomics", "Proteomics", "Adipose", "Analysis_Adipose_Proteomics.csv", "Post3h", TRUE,
  "Muscle_Proteomics", "Proteomics", "Muscle", "Analysis_Muscle_Proteomics.csv", "Pre", TRUE,
  "Muscle_Proteomics", "Proteomics", "Muscle", "Analysis_Muscle_Proteomics.csv", "Post3h", TRUE
)

read_annotation <- function(dataset) {
  if (dataset == "Serum_Proteomics") {
    return(read_csv(path_analysis_ready("Metadata", "Feature_Dictionary_Serum_Proteomics.csv"),
                    show_col_types = FALSE, progress = FALSE) %>%
             transmute(FeatureID = as.character(FeatureID), GeneSymbol = as.character(GeneSymbol)) %>%
             distinct(FeatureID, .keep_all = TRUE))
  }
  if (grepl("Microarray$", dataset)) {
    return(read_csv(path_analysis_ready("Metadata", paste0("Feature_Annotation_", dataset, ".csv")),
                    show_col_types = FALSE, progress = FALSE) %>%
             transmute(FeatureID = as.character(ProbeSetID), GeneSymbol = as.character(SYMBOL)) %>%
             distinct(FeatureID, .keep_all = TRUE))
  }
  read_csv(path_analysis_ready("Metadata", paste0("Feature_Annotation_", dataset, ".csv")),
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
  if (nrow(use) > 20000L) {
    variances <- matrixStats::rowVars(use, na.rm = TRUE)
    use <- use[head(order(variances, decreasing = TRUE, na.last = NA), 20000L), , drop = FALSE]
  }
  value <- suppressWarnings(duplicateCorrelation(use, design, block = family_id)$consensus.correlation)
  if (!is.finite(value)) value <- 0
  value
}

hallmark_raw <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  filter(!is.na(gs_name), !is.na(gene_symbol), nzchar(gene_symbol)) %>% distinct(gs_name, gene_symbol)
hallmark_sets <- lapply(split(hallmark_raw$gene_symbol, hallmark_raw$gs_name), unique)

run_gsea <- function(feature_results, annotation) {
  ranked <- feature_results
  if (!"GeneSymbol" %in% names(ranked)) ranked <- ranked %>% left_join(annotation, by = "FeatureID")
  ranked <- ranked %>%
    filter(!is.na(GeneSymbol), nzchar(GeneSymbol), is.finite(ModeratedT)) %>%
    mutate(GeneSymbol = str_split(GeneSymbol, ";", simplify = TRUE)[, 1]) %>%
    group_by(GeneSymbol) %>% slice_max(abs(ModeratedT), n = 1, with_ties = FALSE) %>% ungroup()
  stats <- ranked$ModeratedT; names(stats) <- ranked$GeneSymbol; stats <- sort(stats, decreasing = TRUE)
  if (length(stats) < 50L) return(NULL)
  suppressWarnings(fgseaMultilevel(
    pathways = hallmark_sets, stats = stats, minSize = minimum_gene_set_size,
    maxSize = maximum_gene_set_size, eps = 1e-20, nproc = 1
  )) %>% as_tibble() %>% transmute(
    PathwayID = pathway,
    Pathway = pathway %>% str_remove("^HALLMARK_") %>% str_replace_all("_", " ") %>%
      str_to_lower() %>% str_to_sentence(),
    Size = size, NES = NES, P_Value = pval, BH_FDR = padj,
    LeadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  )
}

prepare_model <- function(dataset, timepoint, matrix, model, marker) {
  post_time <- if (timepoint == "Post1h") "Post1h" else "Post3h"
  clinical_column <- switch(model,
    Pre_Pre = "Pre_log2",
    Post_Post = paste0(post_time, "_log2"),
    Delta_Delta = if (post_time == "Post1h") "Delta1h" else "Delta3h"
  )
  if (model == "Delta_Delta") {
    metadata <- sample_sheet %>%
      filter(Dataset == dataset, Timepoint %in% c("Pre", post_time), SampleID %in% colnames(matrix)) %>%
      distinct(Clinical_Subject_ID, Timepoint, .keep_all = TRUE) %>%
      select(FamilyID, Clinical_Subject_ID, Role, Timepoint, SampleID) %>%
      pivot_wider(names_from = Timepoint, values_from = SampleID) %>%
      filter(!is.na(Pre), !is.na(.data[[post_time]])) %>% arrange(FamilyID, Role, Clinical_Subject_ID)
    values <- matrix[, metadata[[post_time]], drop = FALSE] - matrix[, metadata$Pre, drop = FALSE]
  } else {
    target_time <- if (model == "Pre_Pre") "Pre" else post_time
    metadata <- sample_sheet %>%
      filter(Dataset == dataset, Timepoint == target_time, SampleID %in% colnames(matrix)) %>%
      distinct(Clinical_Subject_ID, .keep_all = TRUE) %>%
      select(FamilyID, Clinical_Subject_ID, Role, SampleID) %>% arrange(FamilyID, Role, Clinical_Subject_ID)
    values <- matrix[, metadata$SampleID, drop = FALSE]
  }
  colnames(values) <- metadata$Clinical_Subject_ID
  predictor <- clinical %>% filter(Marker == marker) %>%
    transmute(Clinical_Subject_ID, Predictor = .data[[clinical_column]])
  metadata <- metadata %>% inner_join(predictor, by = "Clinical_Subject_ID") %>%
    filter(is.finite(Predictor), !is.na(Role), !is.na(FamilyID)) %>%
    mutate(Predictor_Z = rank_normalize(Predictor)) %>%
    arrange(match(Clinical_Subject_ID, colnames(values))) %>% droplevels()
  values <- values[, metadata$Clinical_Subject_ID, drop = FALSE]
  list(values = values, metadata = metadata, clinical_column = clinical_column)
}

feature_all <- list(); feature_summary <- list(); pathway_all <- list()
loaded_matrices <- list()
for (i in seq_len(nrow(registry))) {
  spec <- registry[i, ]
  if (is.null(loaded_matrices[[spec$Dataset]])) {
    loaded_matrices[[spec$Dataset]] <- read_analysis_matrix(path_analysis_ready(spec$InputFile))
  }
  matrix <- loaded_matrices[[spec$Dataset]]
  annotation <- if (spec$PathwayEligible) read_annotation(spec$Dataset) else NULL
  models <- if (spec$Timepoint == "Pre") "Pre_Pre" else c("Post_Post", "Delta_Delta")
  for (model in models) for (j in seq_len(nrow(marker_registry))) {
    marker <- marker_registry$Marker[j]
    prepared <- prepare_model(spec$Dataset, spec$Timepoint, matrix, model, marker)
    metadata <- prepared$metadata; values <- prepared$values
    if (nrow(metadata) < minimum_subjects || n_distinct(metadata$FamilyID) < 5L) next
    observed <- rowSums(is.finite(values))
    variable <- matrixStats::rowSds(values, na.rm = TRUE) > 0
    keep <- observed >= minimum_subjects & is.finite(variable) & variable
    values <- values[keep, , drop = FALSE]; observed <- observed[keep]
    design <- model.matrix(~ Predictor_Z + Role, data = metadata)
    if (qr(design)$rank != ncol(design)) next
    family_correlation <- estimate_family_correlation(values, design, metadata$FamilyID)
    fit <- lmFit(values, design, block = metadata$FamilyID, correlation = family_correlation) %>%
      eBayes(trend = TRUE, robust = TRUE)
    stats <- topTable(fit, coef = "Predictor_Z", number = Inf, sort.by = "none")
    result <- tibble(
      AnalysisID = analysis_id, Dataset = spec$Dataset, Layer = spec$Layer, Tissue = spec$Tissue,
      Model = model, Timepoint = spec$Timepoint, ClinicalMarker = marker,
      ClinicalMarkerLabel = marker_registry$MarkerLabel[j], ClinicalVariable = prepared$clinical_column,
      FeatureID = rownames(stats), Effect = stats$logFC, ModeratedT = stats$t,
      P_Value = stats$P.Value, BH_FDR = stats$adj.P.Val,
      ObservedN = unname(observed[rownames(stats)]), Subjects = nrow(metadata),
      Families = n_distinct(metadata$FamilyID), FamilyCorrelation = family_correlation
    )
    if (spec$PathwayEligible) result <- result %>% left_join(annotation, by = "FeatureID")
    else result <- result %>% mutate(GeneSymbol = NA_character_)
    key <- paste(spec$Dataset, spec$Timepoint, model, marker, sep = "__")
    feature_all[[key]] <- result
    feature_summary[[key]] <- result %>% summarise(
      Dataset = first(Dataset), Layer = first(Layer), Tissue = first(Tissue), Model = first(Model),
      Timepoint = first(Timepoint), ClinicalMarker = first(ClinicalMarker),
      ClinicalMarkerLabel = first(ClinicalMarkerLabel), Subjects = first(Subjects), Families = first(Families),
      TestedFeatures = n(), NominalP05 = sum(P_Value < 0.05, na.rm = TRUE),
      BH_FDR05 = sum(BH_FDR < 0.05, na.rm = TRUE), BH_FDR10 = sum(BH_FDR < 0.10, na.rm = TRUE),
      BH_FDR20 = sum(BH_FDR < 0.20, na.rm = TRUE), MinimumP = min(P_Value, na.rm = TRUE),
      MinimumBH_FDR = min(BH_FDR, na.rm = TRUE), FamilyCorrelation = first(FamilyCorrelation)
    )
    if (spec$PathwayEligible) {
      pathway <- run_gsea(result, annotation)
      if (!is.null(pathway)) pathway_all[[key]] <- pathway %>% mutate(
        AnalysisID = analysis_id, Dataset = spec$Dataset, Layer = spec$Layer, Tissue = spec$Tissue,
        Model = model, Timepoint = spec$Timepoint, ClinicalMarker = marker,
        ClinicalMarkerLabel = marker_registry$MarkerLabel[j], Subjects = nrow(metadata),
        Families = n_distinct(metadata$FamilyID), Database = "MSigDB Hallmark", .before = 1
      )
    }
    rm(fit, stats, result, values, prepared); gc(verbose = FALSE)
  }
}

feature_results <- bind_rows(feature_all) %>% arrange(Dataset, Model, Timepoint, ClinicalMarker, BH_FDR)
feature_summary_table <- bind_rows(feature_summary) %>% arrange(Dataset, Model, Timepoint, ClinicalMarker)
pathway_results <- bind_rows(pathway_all) %>% arrange(Dataset, Model, Timepoint, ClinicalMarker, BH_FDR)
pathway_summary_table <- pathway_results %>%
  group_by(Dataset, Layer, Tissue, Model, Timepoint, ClinicalMarker, ClinicalMarkerLabel, Subjects, Families) %>%
  summarise(TestedPathways = n(), NominalP05 = sum(P_Value < 0.05, na.rm = TRUE),
            BH_FDR05 = sum(BH_FDR < 0.05, na.rm = TRUE), BH_FDR10 = sum(BH_FDR < 0.10, na.rm = TRUE),
            BH_FDR20 = sum(BH_FDR < 0.20, na.rm = TRUE), MinimumP = min(P_Value, na.rm = TRUE),
            MinimumBH_FDR = min(BH_FDR, na.rm = TRUE), .groups = "drop")

feature_selected <- feature_results %>% filter(P_Value < 0.05 | BH_FDR < 0.20)
pathway_selected <- pathway_results %>% filter(P_Value < 0.05 | BH_FDR < 0.20)
write_csv(feature_summary_table, file.path(output_dir, "Feature_FDR_Summary.csv"))
write_csv(pathway_summary_table, file.path(output_dir, "Hallmark_FDR_Summary.csv"))
write_csv(feature_selected, file.path(output_dir, "Feature_Selected.csv"))
write_csv(pathway_selected, file.path(output_dir, "Hallmark_Selected.csv"))
write_xlsx(list(Feature_summary = feature_summary_table, Hallmark_summary = pathway_summary_table,
                Feature_selected = feature_selected, Hallmark_selected = pathway_selected),
           file.path(output_dir, "Six_Markers_Three_State_Omics.xlsx"))
saveRDS(list(feature_summary = feature_summary_table, pathway_summary = pathway_summary_table,
             feature_results = feature_results, pathway_results = pathway_results),
        file.path(output_dir, "Six_Markers_Three_State_Omics.rds"), compress = "xz")
message("Completed: ", output_dir)
