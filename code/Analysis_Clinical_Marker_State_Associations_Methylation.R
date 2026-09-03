# Methylation pre-pre and post-post models for six clinical exercise markers.
# Delta-delta results are reused from the validated 2026-08-21 release.

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

p_load(dplyr, IlluminaHumanMethylationEPICanno.ilm10b4.hg19, limma, matrixStats,
       methylGSA, msigdbr, readr, stringr, tibble, tidyr, writexl)

analysis_id <- "Six_Markers_Three_State_Methylation_2026-08-22"
output_dir <- path_results("Six_Markers_Three_State_Omics_2026-08-22")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
minimum_subjects <- 10L
maximum_correlation_features <- 20000L
minimum_set_size <- 10L
maximum_set_size <- 500L
role_levels <- c("Daughter", "Mother", "Father")
requested_markers <- c("cOC", "TotalOC", "cOC_ratio", "BCAA", "Serotonin", "Leptin")

rank_normalize <- function(x) {
  out <- rep(NA_real_, length(x)); keep <- is.finite(x); n <- sum(keep)
  if (n >= 3L) out[keep] <- qnorm((rank(x[keep], ties.method = "average") - 0.5) / n)
  out
}
safe_log2 <- function(x) ifelse(is.finite(x) & x > 0, log2(x), NA_real_)
format_pathway <- function(x) x %>% str_remove("^HALLMARK_") %>% str_replace_all("_", " ") %>%
  str_to_lower() %>% str_to_sentence()

clinical <- read_csv(
  file.path(path_intermediate("Clinical_Exercise_Source_Audit_2026-08-21"),
            "Restricted", "OC_BCAA_Serotonin_Long_Internal.csv"),
  show_col_types = FALSE, progress = FALSE
) %>%
  filter(MatchStatus == "Unique locked-cohort match", Marker %in% requested_markers) %>%
  mutate(FamilyID = as.character(FamilyID), Clinical_Subject_ID = as.character(Clinical_Subject_ID),
         Role = factor(Role, levels = role_levels), Value = suppressWarnings(as.numeric(Value))) %>%
  select(FamilyID, Clinical_Subject_ID, Role, Marker, MarkerLabel, Timepoint, Value) %>%
  distinct(FamilyID, Clinical_Subject_ID, Marker, Timepoint, .keep_all = TRUE) %>%
  pivot_wider(names_from = Timepoint, values_from = Value) %>%
  mutate(Pre_log2 = safe_log2(Pre), Post3h_log2 = safe_log2(Post3h))
marker_registry <- clinical %>% distinct(Marker, MarkerLabel) %>% arrange(match(Marker, requested_markers))

sample_sheet <- read_csv(path_analysis_ready("Metadata", "Project_Sample_Sheet.csv"),
                         show_col_types = FALSE, progress = FALSE) %>%
  mutate(FamilyID = as.character(FamilyID), Clinical_Subject_ID = as.character(Clinical_Subject_ID),
         SampleID = tolower(as.character(SampleID)), Role = factor(Role, levels = role_levels),
         Is_Primary_Cohort = as.logical(Is_Primary_Cohort),
         RetainInPrimaryInput = as.logical(RetainInPrimaryInput)) %>%
  filter(AnalysisSet == "Primary", Is_Primary_Cohort,
         is.na(RetainInPrimaryInput) | RetainInPrimaryInput)

registry <- tribble(
  ~Tissue, ~Dataset, ~InputFile,
  "Adipose", "Adipose_Methylation", "Analysis_Adipose_Methylation_MValue.rds",
  "Muscle", "Muscle_Methylation", "Analysis_Muscle_Methylation_MValue.rds"
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

summary_out <- list(); selected_out <- list(); pathway_out <- list()
for (i in seq_len(nrow(registry))) {
  spec <- registry[i, ]
  matrix <- read_analysis_matrix(path_analysis_ready(spec$InputFile))
  for (model in c("Pre_Pre", "Post_Post")) {
    target_time <- if (model == "Pre_Pre") "Pre" else "Post3h"
    clinical_column <- if (model == "Pre_Pre") "Pre_log2" else "Post3h_log2"
    state_metadata <- sample_sheet %>%
      filter(Dataset == spec$Dataset, Timepoint == target_time, SampleID %in% colnames(matrix)) %>%
      distinct(Clinical_Subject_ID, .keep_all = TRUE) %>%
      select(FamilyID, Clinical_Subject_ID, Role, SampleID) %>% arrange(FamilyID, Role, Clinical_Subject_ID)
    for (j in seq_len(nrow(marker_registry))) {
      marker <- marker_registry$Marker[j]
      predictor <- clinical %>% filter(Marker == marker) %>%
        transmute(Clinical_Subject_ID, Predictor = .data[[clinical_column]])
      metadata <- state_metadata %>% inner_join(predictor, by = "Clinical_Subject_ID") %>%
        filter(is.finite(Predictor), !is.na(Role), !is.na(FamilyID)) %>%
        mutate(Predictor_Z = rank_normalize(Predictor)) %>% droplevels()
      values <- matrix[, metadata$SampleID, drop = FALSE]
      colnames(values) <- metadata$Clinical_Subject_ID
      if (nrow(metadata) < minimum_subjects || n_distinct(metadata$FamilyID) < 5L) next
      observed <- rowSums(is.finite(values)); variable <- matrixStats::rowSds(values, na.rm = TRUE) > 0
      keep <- observed >= minimum_subjects & is.finite(variable) & variable
      values <- values[keep, , drop = FALSE]; observed <- observed[keep]
      design <- model.matrix(~ Predictor_Z + Role, data = metadata)
      if (qr(design)$rank != ncol(design)) next
      family_correlation <- estimate_family_correlation(values, design, metadata$FamilyID)
      fit <- lmFit(values, design, block = metadata$FamilyID, correlation = family_correlation) %>%
        eBayes(trend = TRUE, robust = TRUE)
      p_values <- fit$p.value[, "Predictor_Z"]; bh <- p.adjust(p_values, method = "BH")
      result <- tibble(
        AnalysisID = analysis_id, Dataset = spec$Dataset, Tissue = spec$Tissue,
        Model = model, Timepoint = target_time, ClinicalMarker = marker,
        ClinicalMarkerLabel = marker_registry$MarkerLabel[j], FeatureID = rownames(values),
        Effect_MValue = unname(fit$coefficients[, "Predictor_Z"]),
        ModeratedT = unname(fit$t[, "Predictor_Z"]), P_Value = unname(p_values), BH_FDR = unname(bh),
        ObservedN = unname(observed[rownames(values)]), Subjects = nrow(metadata),
        Families = n_distinct(metadata$FamilyID), FamilyCorrelation = family_correlation
      )
      key <- paste(spec$Dataset, model, marker, sep = "__")
      summary_out[[key]] <- result %>% summarise(
        Dataset = dplyr::first(Dataset), Tissue = dplyr::first(Tissue),
        Model = dplyr::first(Model), Timepoint = dplyr::first(Timepoint),
        ClinicalMarker = dplyr::first(ClinicalMarker),
        ClinicalMarkerLabel = dplyr::first(ClinicalMarkerLabel),
        Subjects = dplyr::first(Subjects), Families = dplyr::first(Families), TestedCpGs = n(),
        NominalP05 = sum(P_Value < 0.05, na.rm = TRUE), BH_FDR05 = sum(BH_FDR < 0.05, na.rm = TRUE),
        BH_FDR10 = sum(BH_FDR < 0.10, na.rm = TRUE), BH_FDR20 = sum(BH_FDR < 0.20, na.rm = TRUE),
        MinimumP = min(P_Value, na.rm = TRUE), MinimumBH_FDR = min(BH_FDR, na.rm = TRUE),
        FamilyCorrelation = dplyr::first(FamilyCorrelation)
      )
      selected_out[[key]] <- result %>% filter(P_Value < 1e-4 | BH_FDR < 0.20) %>%
        left_join(cpg_annotation, by = "FeatureID")
      named_p <- pmax(result$P_Value, .Machine$double.xmin); names(named_p) <- result$FeatureID
      methyl_fit <- suppressWarnings(methylglm(
        cpg.pval = named_p, array.type = "EPIC", group = "all", GS.list = hallmark_sets,
        GS.idtype = "SYMBOL", GS.type = "GO", minsize = minimum_set_size,
        maxsize = maximum_set_size, parallel = FALSE
      )) %>% as.data.frame()
      pathway_out[[key]] <- tibble(
        AnalysisID = analysis_id, Dataset = spec$Dataset, Tissue = spec$Tissue,
        Model = model, Timepoint = target_time, ClinicalMarker = marker,
        ClinicalMarkerLabel = marker_registry$MarkerLabel[j], PathwayID = rownames(methyl_fit),
        Pathway = format_pathway(rownames(methyl_fit)), Size = methyl_fit$Size,
        P_Value = methyl_fit$pvalue, BH_FDR = p.adjust(methyl_fit$pvalue, "BH"),
        Subjects = nrow(metadata), Families = n_distinct(metadata$FamilyID),
        Method = "methylglm CpG-count-adjusted Hallmark enrichment"
      )
      rm(fit, result, methyl_fit, values); gc(verbose = FALSE)
    }
  }
  rm(matrix); gc(verbose = FALSE)
}

state_summary <- bind_rows(summary_out) %>% arrange(Dataset, Model, ClinicalMarker)
state_selected <- bind_rows(selected_out) %>% arrange(BH_FDR, P_Value)
state_pathways <- bind_rows(pathway_out) %>% arrange(Dataset, Model, ClinicalMarker, BH_FDR)
state_pathway_summary <- state_pathways %>%
  group_by(Dataset, Tissue, Model, Timepoint, ClinicalMarker, ClinicalMarkerLabel, Subjects, Families) %>%
  summarise(TestedPathways = n(), NominalP05 = sum(P_Value < 0.05), BH_FDR05 = sum(BH_FDR < 0.05),
            BH_FDR10 = sum(BH_FDR < 0.10), BH_FDR20 = sum(BH_FDR < 0.20),
            MinimumP = min(P_Value), MinimumBH_FDR = min(BH_FDR), .groups = "drop")

previous_dir <- path_results("Clinical_Exercise_Methylation_Change_2026-08-21")
delta_summary <- read_csv(file.path(previous_dir, "Methylation_Change_Feature_FDR_Summary.csv"), show_col_types = FALSE) %>%
  filter(ClinicalMarker %in% requested_markers) %>% mutate(Model = "Delta_Delta", Timepoint = "Post3h", BH_FDR20 = NA_integer_)
delta_selected <- read_csv(file.path(previous_dir, "Methylation_Change_Feature_Selected.csv"), show_col_types = FALSE) %>%
  filter(ClinicalMarker %in% requested_markers) %>% mutate(Model = "Delta_Delta", Timepoint = "Post3h")
delta_pathway_summary <- read_csv(file.path(previous_dir, "Methylation_Change_Hallmark_FDR_Summary.csv"), show_col_types = FALSE) %>%
  filter(ClinicalMarker %in% requested_markers) %>% mutate(Model = "Delta_Delta", Timepoint = "Post3h", NominalP05 = NA_integer_, BH_FDR20 = NA_integer_)
delta_pathways <- read_csv(file.path(previous_dir, "Methylation_Change_Hallmark_Selected.csv"), show_col_types = FALSE) %>%
  filter(ClinicalMarker %in% requested_markers) %>% mutate(Model = "Delta_Delta", Timepoint = "Post3h")

write_csv(state_summary, file.path(output_dir, "Methylation_State_Feature_FDR_Summary.csv"))
write_csv(state_selected, file.path(output_dir, "Methylation_State_Feature_Selected.csv"))
write_csv(state_pathway_summary, file.path(output_dir, "Methylation_State_Hallmark_FDR_Summary.csv"))
write_csv(state_pathways %>% filter(P_Value < 0.05 | BH_FDR < 0.20),
          file.path(output_dir, "Methylation_State_Hallmark_Selected.csv"))
write_xlsx(list(State_feature_summary = state_summary, Delta_feature_summary = delta_summary,
                State_feature_selected = state_selected, Delta_feature_selected = delta_selected,
                State_pathway_summary = state_pathway_summary, Delta_pathway_summary = delta_pathway_summary,
                State_pathway_selected = state_pathways %>% filter(P_Value < 0.05 | BH_FDR < 0.20),
                Delta_pathway_selected = delta_pathways),
           file.path(output_dir, "Six_Markers_Three_State_Methylation.xlsx"))
message("Completed methylation state models: ", output_dir)
