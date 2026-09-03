# Methylation component of the unified pre/Post3h/response association audit.
# Uses formal EPIC M-value matrices and CpG-count-adjusted Hallmark enrichment.

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

methylation_library <- Sys.getenv(
  "ATM_METHYLATION_R_LIBRARY",
  unset = path_project("R_Library_Methylation")
)
if (dir.exists(methylation_library)) {
  .libPaths(unique(c(methylation_library, DIR_R_LIBRARY, .libPaths())))
}

p_load(dplyr, IlluminaHumanMethylationEPICanno.ilm10b4.hg19, limma,
       matrixStats, methylGSA, msigdbr, readr, stringr, tibble, tidyr, writexl)

analysis_id <- "Unified_Pre_Post_Response_Association_Audit_2026-08-28"
output_dir <- path_results(analysis_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
minimum_subjects <- 10L
maximum_correlation_features <- 20000L
minimum_set_size <- 10L
maximum_set_size <- 500L
role_levels <- c("Daughter", "Mother", "Father")

rank_normalize <- function(x) {
  out <- rep(NA_real_, length(x)); keep <- is.finite(x); n <- sum(keep)
  if (n >= 3L) out[keep] <- qnorm((rank(x[keep], ties.method = "average") - 0.5) / n)
  out
}

clinical <- read_csv(path_analysis_ready("Metadata", "Clinical_Observed_Primary.csv"),
                     show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    FamilyID = as.character(FamilyID), Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = factor(Role, levels = role_levels),
    BodyFat = suppressWarnings(as.numeric(Fat_p_G)),
    HOMA_IR = suppressWarnings(as.numeric(HOMA_IR_84m_G)),
    Matsuda_ISI = suppressWarnings(as.numeric(ISI_G)),
    hsCRP = suppressWarnings(as.numeric(CRP_84m_G)),
    TotalOC = suppressWarnings(as.numeric(S_TotalOC_84m_G)),
    cOC = suppressWarnings(as.numeric(S_cOC_84m_G)),
    cOC_ratio = suppressWarnings(as.numeric(cOC_totalOC_ratio_84m_G))
  ) %>% distinct(Clinical_Subject_ID, .keep_all = TRUE) %>%
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

sample_sheet <- read_csv(path_analysis_ready("Metadata", "Project_Sample_Sheet.csv"),
                         show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    FamilyID = as.character(FamilyID), Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    SampleID = tolower(as.character(SampleID)), Role = factor(Role, levels = role_levels),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort),
    RetainInPrimaryInput = as.logical(RetainInPrimaryInput)
  ) %>%
  filter(AnalysisSet == "Primary", Is_Primary_Cohort,
         is.na(RetainInPrimaryInput) | RetainInPrimaryInput)

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
  ~Dataset, ~Layer, ~Tissue, ~InputFile,
  "Adipose_Methylation", "DNA methylation", "Adipose", "Analysis_Adipose_Methylation_MValue.rds",
  "Muscle_Methylation", "DNA methylation", "Skeletal muscle", "Analysis_Muscle_Methylation_MValue.rds"
)
dataset_filter <- Sys.getenv("ATM_METHYLATION_DATASET", unset = "")
if (nzchar(dataset_filter)) {
  datasets <- datasets %>% filter(Dataset == dataset_filter)
  if (!nrow(datasets)) stop("Unknown ATM_METHYLATION_DATASET: ", dataset_filter)
}
output_suffix <- if (nzchar(dataset_filter)) paste0("_", dataset_filter) else ""
states <- tribble(
  ~State, ~StateLabel, ~Contrast,
  "PreExercise", "Pre-exercise molecular profile", "Pre",
  "Post3hAbsolute", "Post-exercise molecular profile", "Post3h",
  "ExerciseChange", "Exercise-induced molecular change", "Post3h minus Pre"
)

hallmark_raw <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  filter(!is.na(gs_name), !is.na(gene_symbol), nzchar(gene_symbol)) %>% distinct(gs_name, gene_symbol)
hallmark_sets <- lapply(split(hallmark_raw$gene_symbol, hallmark_raw$gs_name), unique)
cpg_annotation <- methylGSA:::getAnnot("EPIC", "all") %>% as.data.frame() %>% as_tibble() %>%
  transmute(FeatureID = as.character(Name), GeneSymbol = as.character(UCSC_RefGene_Name),
            GeneRegion = as.character(UCSC_RefGene_Group)) %>% distinct(FeatureID, .keep_all = TRUE)

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

summary_list <- list(); selected_list <- list(); pathway_list <- list(); availability <- list()
counter <- 0L
for (i in seq_len(nrow(datasets))) {
  spec <- datasets[i, ]
  message("Loading ", spec$Dataset)
  matrix <- read_analysis_matrix(path_analysis_ready(spec$InputFile))
  for (s in seq_len(nrow(states))) {
    state_spec <- states[s, ]
    prepared <- prepare_state(spec$Dataset, state_spec$State, matrix)
    availability[[paste(spec$Dataset, state_spec$State, sep = "__")]] <- tibble(
      Dataset = spec$Dataset, Layer = spec$Layer, Tissue = spec$Tissue,
      State = state_spec$State, MatrixSubjects = nrow(prepared$metadata),
      MatrixFamilies = n_distinct(prepared$metadata$FamilyID)
    )
    for (j in seq_len(nrow(predictors))) {
      pred <- predictors[j, ]
      metadata <- prepared$metadata
      pred_col <- pred$PredictorColumn
      complete <- if (pred$PredictorType == "Categorical") !is.na(metadata[[pred_col]]) else is.finite(metadata[[pred_col]])
      complete <- complete & !is.na(metadata$Role) & !is.na(metadata$FamilyID)
      metadata <- metadata[complete, , drop = FALSE] %>% droplevels()
      values <- prepared$values[, complete, drop = FALSE]
      if (nrow(metadata) < minimum_subjects || n_distinct(metadata$FamilyID) < 5L) next
      if (pred$PredictorType == "Categorical" && nlevels(metadata[[pred_col]]) < 2L) next
      observed <- rowSums(is.finite(values))
      variable <- matrixStats::rowSds(values, na.rm = TRUE) > 0
      keep <- observed >= minimum_subjects & is.finite(variable) & variable
      values <- values[keep, , drop = FALSE]; observed <- observed[keep]
      design <- model.matrix(reformulate(c(pred_col, "Role")), data = metadata)
      if (qr(design)$rank != ncol(design)) next
      coefficient <- if (pred$PredictorType == "Categorical") paste0(pred_col, "Obese") else pred_col
      if (!coefficient %in% colnames(design)) next
      message("Fitting ", spec$Dataset, " / ", state_spec$State, " / ", pred$Predictor)
      family_correlation <- estimate_family_correlation(values, design, metadata$FamilyID)
      fit <- lmFit(values, design, block = metadata$FamilyID, correlation = family_correlation) %>%
        eBayes(trend = TRUE, robust = TRUE)
      p_values <- fit$p.value[, coefficient]
      bh <- p.adjust(p_values, method = "BH")
      result <- tibble(
        AnalysisFramework = "Baseline phenotype association", Dataset = spec$Dataset,
        Layer = spec$Layer, Tissue = spec$Tissue, State = state_spec$State,
        StateLabel = state_spec$StateLabel, Contrast = state_spec$Contrast,
        PredictorFamily = pred$PredictorFamily, Predictor = pred$Predictor,
        PredictorLabel = pred$PredictorLabel, EffectDefinition = pred$EffectDefinition,
        FeatureID = rownames(values), Effect = unname(fit$coefficients[, coefficient]),
        ModeratedT = unname(fit$t[, coefficient]), P_Value = unname(p_values), BH_FDR = unname(bh),
        ObservedN = unname(observed[rownames(values)]), Subjects = nrow(metadata),
        Families = n_distinct(metadata$FamilyID), DaughterN = sum(metadata$Role == "Daughter"),
        MotherN = sum(metadata$Role == "Mother"), FatherN = sum(metadata$Role == "Father"),
        FamilyCorrelation = family_correlation
      )
      key <- paste(spec$Dataset, state_spec$State, pred$Predictor, sep = "__")
      summary_list[[key]] <- result %>% summarise(
        AnalysisFramework = dplyr::first(AnalysisFramework), PredictorFamily = dplyr::first(PredictorFamily),
        Predictor = dplyr::first(Predictor), PredictorLabel = dplyr::first(PredictorLabel), Dataset = dplyr::first(Dataset),
        Layer = dplyr::first(Layer), Tissue = dplyr::first(Tissue), State = dplyr::first(State), StateLabel = dplyr::first(StateLabel),
        Contrast = dplyr::first(Contrast), EffectDefinition = dplyr::first(EffectDefinition), Subjects = dplyr::first(Subjects),
        Families = dplyr::first(Families), TestedFeatures = n(), NominalP05 = sum(P_Value < 0.05, na.rm = TRUE),
        BH_FDR05 = sum(BH_FDR < 0.05, na.rm = TRUE), BH_FDR10 = sum(BH_FDR < 0.10, na.rm = TRUE),
        BH_FDR20 = sum(BH_FDR < 0.20, na.rm = TRUE), MinimumP = min(P_Value, na.rm = TRUE),
        MinimumBH_FDR = min(BH_FDR, na.rm = TRUE), FamilyCorrelation = dplyr::first(FamilyCorrelation)
      )
      selected_list[[key]] <- result %>% filter(P_Value < 1e-4 | BH_FDR < 0.20) %>%
        left_join(cpg_annotation, by = "FeatureID")
      named_p <- pmax(result$P_Value, .Machine$double.xmin); names(named_p) <- result$FeatureID
      methyl_fit <- suppressWarnings(methylglm(
        cpg.pval = named_p, array.type = "EPIC", group = "all", GS.list = hallmark_sets,
        GS.idtype = "SYMBOL", GS.type = "GO", minsize = minimum_set_size,
        maxsize = maximum_set_size, parallel = FALSE
      )) %>% as.data.frame()
      pathway_list[[key]] <- tibble(
        AnalysisFramework = "Baseline phenotype association", Dataset = spec$Dataset,
        Layer = spec$Layer, Tissue = spec$Tissue, State = state_spec$State,
        StateLabel = state_spec$StateLabel, Contrast = state_spec$Contrast,
        PredictorFamily = pred$PredictorFamily, Predictor = pred$Predictor,
        PredictorLabel = pred$PredictorLabel, Database = "MSigDB Hallmark",
        Method = "CpG-count-adjusted methylglm", PathwayID = rownames(methyl_fit),
        Pathway = rownames(methyl_fit) %>% str_remove("^HALLMARK_") %>%
          str_replace_all("_", " ") %>% str_to_lower() %>% str_to_sentence(),
        Size = methyl_fit$Size, P_Value = methyl_fit$pvalue,
        BH_FDR = p.adjust(methyl_fit$pvalue, method = "BH"),
        Subjects = nrow(metadata), Families = n_distinct(metadata$FamilyID)
      )
      rm(fit, result, methyl_fit, values); invisible(gc())
      counter <- counter + 1L
    }
    rm(prepared); invisible(gc())
  }
  rm(matrix); invisible(gc())
}

feature_summary <- bind_rows(summary_list) %>%
  arrange(PredictorFamily, Predictor, Dataset, State)
feature_selected <- bind_rows(selected_list) %>%
  arrange(PredictorFamily, Predictor, Dataset, State, BH_FDR, P_Value)
pathway_all <- bind_rows(pathway_list)
pathway_summary <- pathway_all %>%
  group_by(AnalysisFramework, PredictorFamily, Predictor, PredictorLabel, Dataset, Layer, Tissue,
           State, StateLabel, Contrast, Database, Method) %>%
  summarise(
    Subjects = dplyr::first(Subjects), Families = dplyr::first(Families), TestedPathways = n(),
    NominalP05 = sum(P_Value < 0.05, na.rm = TRUE),
    BH_FDR05 = sum(BH_FDR < 0.05, na.rm = TRUE),
    BH_FDR10 = sum(BH_FDR < 0.10, na.rm = TRUE),
    BH_FDR20 = sum(BH_FDR < 0.20, na.rm = TRUE),
    MinimumP = min(P_Value, na.rm = TRUE), MinimumBH_FDR = min(BH_FDR, na.rm = TRUE),
    .groups = "drop"
  )
pathway_selected <- pathway_all %>% filter(P_Value < 0.05 | BH_FDR < 0.20) %>%
  arrange(PredictorFamily, Predictor, Dataset, State, BH_FDR, P_Value)

write_csv(feature_summary, file.path(output_dir, paste0("Methylation_Feature_Significance_Summary", output_suffix, ".csv")))
write_csv(pathway_summary, file.path(output_dir, paste0("Methylation_Pathway_Significance_Summary", output_suffix, ".csv")))
write_csv(feature_selected, file.path(output_dir, paste0("Methylation_Feature_Selected", output_suffix, ".csv")))
write_csv(pathway_selected, file.path(output_dir, paste0("Methylation_Pathway_Selected", output_suffix, ".csv")))
write_csv(bind_rows(availability), file.path(output_dir, paste0("Methylation_Sample_Availability", output_suffix, ".csv")))
write_xlsx(list(
  Feature_summary = feature_summary, Pathway_summary = pathway_summary,
  Feature_selected = feature_selected, Pathway_selected = pathway_selected,
  Sample_availability = bind_rows(availability)
), file.path(output_dir, paste0("Unified_Methylation_Pre_Post_Response_Audit", output_suffix, ".xlsx")))
saveRDS(list(feature_summary = feature_summary, pathway_summary = pathway_summary,
             feature_selected = feature_selected, pathway_results = pathway_all,
             availability = bind_rows(availability)),
        file.path(output_dir, paste0("Unified_Methylation_Pre_Post_Response_Audit", output_suffix, ".rds")), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(output_dir, paste0("R_SESSION_INFO_Methylation", output_suffix, ".txt")))
message("Completed methylation audit: ", output_dir)
