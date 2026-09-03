# GO biological-process GSEA for matched osteocalcin molecular states and changes.
#
# This script reads the validated full non-methylation statistics from the
# six-marker audit and writes an osteocalcin-only release. Source data and
# upstream analysis results are never modified.

.local_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.local_script)) .local_script <- NA_character_
.config_candidates <- unique(c(
  file.path(getwd(), "code", "Fig00_Config.R"),
  file.path(getwd(), "Fig00_Config.R"),
  if (!is.na(.local_script)) file.path(dirname(.local_script), "Fig00_Config.R") else NA_character_
))
.config_file <- .config_candidates[file.exists(.config_candidates)][1]
if (is.na(.config_file)) stop("Cannot locate Fig00_Config.R.")
source(.config_file)
.project_root <- PROJECT_ROOT
rm(.local_script, .config_candidates, .config_file)

p_load(dplyr, fgsea, msigdbr, readr, stringr, tibble, writexl)

set.seed(20260830)

source_dir <- path_results("Six_Markers_Three_State_Omics_2026-08-22")
source_rds <- file.path(source_dir, "Six_Markers_Three_State_Omics.rds")
module_root <- file.path(
  .project_root, "Pre_Post_Response_Analysis", "Results",
  "08_Osteocalcin_Molecular_State_and_Response"
)
feature_dir <- file.path(module_root, "01_Feature_Level")
pathway_dir <- file.path(module_root, "02_Pathway_Enrichment")
dir.create(feature_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)

markers <- c("cOC", "TotalOC", "cOC_ratio")
minimum_gene_set_size <- 10L
maximum_gene_set_size <- 500L

format_state <- function(model, timepoint) {
  dplyr::case_when(
    model == "Pre_Pre" ~ "Pre OC vs Pre omics",
    model == "Post_Post" & timepoint == "Post1h" ~ "Post1h OC vs Post1h omics",
    model == "Post_Post" & timepoint == "Post3h" ~ "Post3h OC vs Post3h omics",
    model == "Delta_Delta" & timepoint == "Post1h" ~ "Delta(Post1h-Pre) OC vs Delta(Post1h-Pre) omics",
    model == "Delta_Delta" & timepoint == "Post3h" ~ "Delta(Post3h-Pre) OC vs Delta(Post3h-Pre) omics",
    TRUE ~ paste(model, timepoint)
  )
}

source_object <- readRDS(source_rds)
feature_results <- source_object$feature_results %>%
  filter(ClinicalMarker %in% markers) %>%
  mutate(
    Comparison = format_state(Model, Timepoint),
    FDR_Tier = case_when(
      BH_FDR < 0.05 ~ "BH-FDR < 0.05",
      BH_FDR < 0.10 ~ "0.05 <= BH-FDR < 0.10",
      BH_FDR < 0.20 ~ "0.10 <= BH-FDR < 0.20",
      TRUE ~ "BH-FDR >= 0.20"
    )
  )

feature_summary <- feature_results %>%
  group_by(Dataset, Layer, Tissue, Model, Timepoint, Comparison,
           ClinicalMarker, ClinicalMarkerLabel, Subjects, Families,
           FamilyCorrelation) %>%
  summarise(
    TestedFeatures = n(),
    NominalP05 = sum(P_Value < 0.05, na.rm = TRUE),
    BH_FDR05 = sum(BH_FDR < 0.05, na.rm = TRUE),
    BH_FDR10 = sum(BH_FDR < 0.10, na.rm = TRUE),
    BH_FDR20 = sum(BH_FDR < 0.20, na.rm = TRUE),
    MinimumP = min(P_Value, na.rm = TRUE),
    MinimumBH_FDR = min(BH_FDR, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(ClinicalMarker, Model, Timepoint, Dataset)

feature_candidates <- feature_results %>%
  filter(P_Value < 0.05 | BH_FDR < 0.20) %>%
  arrange(ClinicalMarker, Model, Timepoint, Dataset, BH_FDR, P_Value)

saveRDS(
  list(
    feature_results = feature_results,
    feature_summary = feature_summary,
    feature_candidates = feature_candidates,
    source_rds = source_rds
  ),
  file.path(feature_dir, "OC_NonMethylation_Full_Statistics.rds"),
  compress = "xz"
)
write_csv(feature_summary, file.path(feature_dir, "OC_NonMethylation_Feature_Summary.csv"))
write_csv(feature_candidates, file.path(feature_dir, "OC_NonMethylation_Feature_Candidates.csv.gz"))

go_raw <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
  filter(!is.na(gs_name), !is.na(gene_symbol), nzchar(gene_symbol)) %>%
  distinct(gs_name, gene_symbol)
go_sets <- lapply(split(go_raw$gene_symbol, go_raw$gs_name), unique)

run_go_gsea <- function(data) {
  ranked <- data %>%
    filter(!is.na(GeneSymbol), nzchar(GeneSymbol), is.finite(ModeratedT)) %>%
    mutate(GeneSymbol = str_split(GeneSymbol, ";", simplify = TRUE)[, 1]) %>%
    group_by(GeneSymbol) %>%
    slice_max(abs(ModeratedT), n = 1, with_ties = FALSE) %>%
    ungroup()
  stats <- ranked$ModeratedT
  names(stats) <- ranked$GeneSymbol
  stats <- sort(stats, decreasing = TRUE)
  if (length(stats) < 50L) return(tibble())
  suppressWarnings(
    fgseaMultilevel(
      pathways = go_sets,
      stats = stats,
      minSize = minimum_gene_set_size,
      maxSize = maximum_gene_set_size,
      eps = 1e-20,
      nproc = 1
    )
  ) %>%
    as_tibble() %>%
    transmute(
      PathwayID = pathway,
      Pathway = pathway %>%
        str_remove("^GOBP_") %>%
        str_replace_all("_", " ") %>%
        str_to_lower() %>%
        str_to_sentence(),
      Size = size,
      NES = NES,
      P_Value = pval,
      BH_FDR = padj,
      LeadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
    )
}

group_keys <- feature_results %>%
  filter(!is.na(GeneSymbol), nzchar(GeneSymbol)) %>%
  distinct(Dataset, Layer, Tissue, Model, Timepoint, Comparison,
           ClinicalMarker, ClinicalMarkerLabel, Subjects, Families)

go_results <- vector("list", nrow(group_keys))
for (i in seq_len(nrow(group_keys))) {
  key <- group_keys[i, ]
  subset <- feature_results %>%
    filter(
      Dataset == key$Dataset,
      Model == key$Model,
      Timepoint == key$Timepoint,
      ClinicalMarker == key$ClinicalMarker
    )
  result <- run_go_gsea(subset)
  if (nrow(result)) {
    go_results[[i]] <- result %>%
      mutate(
        Dataset = key$Dataset,
        Layer = key$Layer,
        Tissue = key$Tissue,
        Model = key$Model,
        Timepoint = key$Timepoint,
        Comparison = key$Comparison,
        ClinicalMarker = key$ClinicalMarker,
        ClinicalMarkerLabel = key$ClinicalMarkerLabel,
        Subjects = key$Subjects,
        Families = key$Families,
        Database = "MSigDB C5 GO biological process",
        Method = "Full-ranked GSEA using limma moderated t statistic",
        .before = 1
      )
  }
  message(sprintf("GO-BP GSEA %d/%d: %s | %s | %s | %s",
                  i, nrow(group_keys), key$ClinicalMarker, key$Dataset,
                  key$Model, key$Timepoint))
}

go_results <- bind_rows(go_results) %>%
  mutate(
    FDR_Tier = case_when(
      BH_FDR < 0.05 ~ "BH-FDR < 0.05",
      BH_FDR < 0.10 ~ "0.05 <= BH-FDR < 0.10",
      BH_FDR < 0.20 ~ "0.10 <= BH-FDR < 0.20",
      TRUE ~ "BH-FDR >= 0.20"
    )
  ) %>%
  arrange(ClinicalMarker, Model, Timepoint, Dataset, BH_FDR, P_Value)

go_summary <- go_results %>%
  group_by(Dataset, Layer, Tissue, Model, Timepoint, Comparison,
           ClinicalMarker, ClinicalMarkerLabel, Subjects, Families,
           Database, Method) %>%
  summarise(
    TestedPathways = n(),
    NominalP05 = sum(P_Value < 0.05, na.rm = TRUE),
    BH_FDR05 = sum(BH_FDR < 0.05, na.rm = TRUE),
    BH_FDR10 = sum(BH_FDR < 0.10, na.rm = TRUE),
    BH_FDR20 = sum(BH_FDR < 0.20, na.rm = TRUE),
    MinimumP = min(P_Value, na.rm = TRUE),
    MinimumBH_FDR = min(BH_FDR, na.rm = TRUE),
    .groups = "drop"
  )

go_selected <- go_results %>% filter(P_Value < 0.05 | BH_FDR < 0.20)
write_csv(go_results, file.path(pathway_dir, "OC_NonMethylation_GO_BP_GSEA_All.csv.gz"))
write_csv(go_selected, file.path(pathway_dir, "OC_NonMethylation_GO_BP_GSEA_Selected.csv"))
write_csv(go_summary, file.path(pathway_dir, "OC_NonMethylation_GO_BP_GSEA_Summary.csv"))
write_xlsx(
  list(
    Feature_summary = feature_summary,
    Feature_FDR05 = feature_results %>% filter(BH_FDR < 0.05),
    GO_BP_summary = go_summary,
    GO_BP_FDR05 = go_results %>% filter(BH_FDR < 0.05),
    GO_BP_FDR10 = go_results %>% filter(BH_FDR < 0.10),
    GO_BP_FDR20 = go_results %>% filter(BH_FDR < 0.20)
  ),
  file.path(pathway_dir, "OC_NonMethylation_GO_BP_GSEA.xlsx")
)

writeLines(c(
  paste0("Run time: ", format(Sys.time(), tz = "Asia/Shanghai")),
  paste0("Source RDS: ", source_rds),
  "Markers: TotalOC, cOC, cOC_ratio",
  "Feature model: omic state/change ~ rank-normalized matched-state OC + Role; FamilyID duplicateCorrelation block",
  "GO-BP method: full-ranked fgseaMultilevel using limma moderated t statistics",
  "BH correction: within dataset x model x timepoint x marker",
  "BH-FDR <0.05 is significant; <0.10 suggestive; <0.20 exploratory."
), file.path(pathway_dir, "RUN_LOG_GO_BP.txt"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(pathway_dir, "R_SESSION_INFO_GO_BP.txt"))
message("Osteocalcin GO-BP GSEA release written to: ", module_root)
