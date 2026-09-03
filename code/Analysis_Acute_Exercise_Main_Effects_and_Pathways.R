# Consolidated acute-exercise main-effect and pathway-enrichment analysis.
# Raw and analysis-ready inputs are read only. All outputs are written to the
# Pre_Post_Response_Analysis archive.

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

p_load(digest, dplyr, fgsea, limma, matrixStats, msigdbr, readr, readxl,
       stringr, tibble, tidyr, writexl)

set.seed(20260828)
archive_root <- file.path(.project_root, "Pre_Post_Response_Analysis")
feature_dir <- file.path(archive_root, "Results", "05_Exercise_Main_Effect")
pathway_dir <- file.path(archive_root, "Results", "06_Pathway_Enrichment")
dir.create(feature_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)

minimum_pairs <- 10L
maximum_correlation_features <- 20000L
minimum_gene_set_size <- 10L
maximum_gene_set_size <- 500L
role_levels <- c("Daughter", "Mother", "Father")

clinical <- read_csv(path_analysis_ready("Metadata", "Clinical_Observed_Primary.csv"),
                     show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = factor(Role, levels = role_levels),
    BodyFat = suppressWarnings(as.numeric(Fat_p_G)),
    Adiposity = factor(classify_body_fat(BodyFat), levels = c("Lean", "Obese"))
  ) %>%
  distinct(Clinical_Subject_ID, .keep_all = TRUE)

sample_sheet <- read_csv(path_analysis_ready("Metadata", "Project_Sample_Sheet.csv"),
                         show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    FamilyID = as.character(FamilyID),
    Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    SampleID = tolower(as.character(SampleID)),
    Role = factor(Role, levels = role_levels),
    Is_Primary_Cohort = as.logical(Is_Primary_Cohort),
    RetainInPrimaryInput = as.logical(RetainInPrimaryInput)
  ) %>%
  filter(
    AnalysisSet == "Primary", Is_Primary_Cohort,
    is.na(RetainInPrimaryInput) | RetainInPrimaryInput
  )

datasets <- tribble(
  ~Dataset, ~Layer, ~Tissue, ~InputFile, ~AnnotationFile, ~LabelColumn, ~Contrasts,
  "Serum_Metabonomics", "Metabolomics", "Serum", "Analysis_Serum_Metabonomics.csv",
  "Feature_Dictionary_Serum_Metabonomics.csv", "FeatureID", "Post1h;Post3h",
  "Serum_Proteomics", "Proteomics", "Serum", "Analysis_Serum_Proteomics.csv",
  "Feature_Dictionary_Serum_Proteomics.csv", "GeneSymbol", "Post1h;Post3h",
  "Adipose_Proteomics", "Proteomics", "Adipose", "Analysis_Adipose_Proteomics.csv",
  "Feature_Annotation_Adipose_Proteomics.csv", "IntegrationGene", "Post3h",
  "Muscle_Proteomics", "Proteomics", "Skeletal muscle", "Analysis_Muscle_Proteomics.csv",
  "Feature_Annotation_Muscle_Proteomics.csv", "IntegrationGene", "Post3h"
)

read_feature_annotation <- function(spec) {
  path <- path_analysis_ready("Metadata", spec$AnnotationFile)
  tab <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  label <- as.character(tab[[spec$LabelColumn]])
  if (spec$Layer == "Proteomics" && "PG.Genes" %in% names(tab)) {
    fallback <- as.character(tab$PG.Genes)
    replace <- is.na(label) | !nzchar(trimws(label))
    label[replace] <- fallback[replace]
  }
  label[is.na(label) | !nzchar(trimws(label))] <- as.character(tab$FeatureID)[
    is.na(label) | !nzchar(trimws(label))
  ]
  tibble(FeatureID = as.character(tab$FeatureID), FeatureLabel = label,
         GeneSymbol = if (spec$Layer == "Proteomics") label else NA_character_) %>%
    distinct(FeatureID, .keep_all = TRUE)
}

prepare_delta <- function(spec, post_timepoint) {
  matrix <- read_analysis_matrix(path_analysis_ready(spec$InputFile))
  available <- sample_sheet %>%
    filter(
      Dataset == spec$Dataset,
      Timepoint %in% c("Pre", post_timepoint),
      SampleID %in% colnames(matrix)
    ) %>%
    inner_join(clinical, by = c("FamilyID", "Clinical_Subject_ID", "Role")) %>%
    distinct(Clinical_Subject_ID, Timepoint, .keep_all = TRUE) %>%
    select(FamilyID, Clinical_Subject_ID, Role, BodyFat, Adiposity, Timepoint, SampleID) %>%
    pivot_wider(names_from = Timepoint, values_from = SampleID) %>%
    filter(!is.na(Pre), !is.na(.data[[post_timepoint]])) %>%
    arrange(FamilyID, Role, Clinical_Subject_ID)
  pre <- matrix[, available$Pre, drop = FALSE]
  post <- matrix[, available[[post_timepoint]], drop = FALSE]
  delta <- post - pre
  colnames(delta) <- available$Clinical_Subject_ID
  list(delta = delta, metadata = available,
       average = rowMeans((pre + post) / 2, na.rm = TRUE), total = nrow(matrix))
}

fit_main_effect <- function(prepared, spec, post_timepoint) {
  metadata <- prepared$metadata
  adiposity <- as.numeric(metadata$Adiposity == "Obese")
  role_mother <- as.numeric(metadata$Role == "Mother")
  role_father <- as.numeric(metadata$Role == "Father")
  design <- cbind(
    Overall = 1,
    AdiposityObese = adiposity - mean(adiposity),
    RoleMother = role_mother - mean(role_mother),
    RoleFather = role_father - mean(role_father)
  )
  observed <- rowSums(is.finite(prepared$delta))
  variable <- matrixStats::rowSds(prepared$delta, na.rm = TRUE) > 0
  keep <- observed >= minimum_pairs & is.finite(variable) & variable
  values <- prepared$delta[keep, , drop = FALSE]
  observed <- observed[keep]
  average <- prepared$average[keep]
  correlation_input <- values
  if (nrow(correlation_input) > maximum_correlation_features) {
    variance <- matrixStats::rowVars(correlation_input, na.rm = TRUE)
    index <- head(order(variance, decreasing = TRUE, na.last = NA), maximum_correlation_features)
    correlation_input <- correlation_input[index, , drop = FALSE]
  }
  family_correlation <- suppressWarnings(
    duplicateCorrelation(correlation_input, design, block = metadata$FamilyID)$consensus.correlation
  )
  if (!is.finite(family_correlation)) family_correlation <- 0
  fit <- lmFit(values, design, block = metadata$FamilyID,
               correlation = family_correlation) %>%
    eBayes(trend = TRUE, robust = TRUE)
  stats <- topTable(fit, coef = "Overall", number = Inf, sort.by = "none")
  tibble(
    Analysis = "Acute exercise main effect",
    Dataset = spec$Dataset, Layer = spec$Layer, Tissue = spec$Tissue,
    Contrast = paste0(post_timepoint, " minus Pre"),
    FeatureID = rownames(stats), Effect = stats$logFC,
    ModeratedT = stats$t, P_Value = stats$P.Value,
    BH_FDR = p.adjust(stats$P.Value, method = "BH"),
    AverageAbundance = unname(average[rownames(stats)]),
    ObservedN = unname(observed[rownames(stats)]),
    PairedSubjects = nrow(metadata), Families = n_distinct(metadata$FamilyID),
    DaughterN = sum(metadata$Role == "Daughter"),
    MotherN = sum(metadata$Role == "Mother"), FatherN = sum(metadata$Role == "Father"),
    FamilyCorrelation = family_correlation,
    EffectDefinition = if_else(
      spec$Layer == "Proteomics",
      "Adjusted mean paired log2-relative-abundance difference",
      "Adjusted mean paired difference on the supplied analysis scale"
    )
  )
}

feature_results <- list()
availability <- list()
counter <- 0L
for (i in seq_len(nrow(datasets))) {
  spec <- datasets[i, ]
  annotation <- read_feature_annotation(spec)
  for (post_timepoint in str_split(spec$Contrasts, ";", simplify = TRUE)) {
    message("Fitting ", spec$Dataset, " / ", post_timepoint, " minus Pre")
    prepared <- prepare_delta(spec, post_timepoint)
    result <- fit_main_effect(prepared, spec, post_timepoint) %>%
      left_join(annotation, by = "FeatureID")
    counter <- counter + 1L
    feature_results[[counter]] <- result
    availability[[counter]] <- tibble(
      Dataset = spec$Dataset, Layer = spec$Layer, Tissue = spec$Tissue,
      Contrast = paste0(post_timepoint, " minus Pre"),
      PairedSubjects = nrow(prepared$metadata),
      Families = n_distinct(prepared$metadata$FamilyID),
      DaughterN = sum(prepared$metadata$Role == "Daughter"),
      MotherN = sum(prepared$metadata$Role == "Mother"),
      FatherN = sum(prepared$metadata$Role == "Father"),
      TestedFeatures = nrow(result)
    )
    rm(prepared, result); invisible(gc())
  }
}

# The formal tissue methylation response fits are reused from the locked Figure 3 cache.
cache_dir <- file.path(.project_root, "R_Cache", "Fig03DE")
for (dataset in c("Adipose_Methylation", "Muscle_Methylation")) {
  cache_file <- file.path(cache_dir, paste0("FeatureResults_", dataset, ".rds"))
  if (!file.exists(cache_file)) stop("Missing frozen methylation response cache: ", cache_file)
  cache <- readRDS(cache_file)
  result <- as_tibble(cache$Results) %>%
    transmute(
      Analysis = "Acute exercise main effect", Dataset,
      Layer = "DNA methylation",
      Tissue = if_else(Dataset == "Adipose_Methylation", "Adipose", "Skeletal muscle"),
      Contrast = "Post3h minus Pre", FeatureID, Effect, ModeratedT,
      P_Value, BH_FDR, AverageAbundance = AverageExpression,
      ObservedN, PairedSubjects, Families,
      DaughterN = NA_integer_, MotherN = NA_integer_, FatherN = NA_integer_,
      FamilyCorrelation,
      EffectDefinition = "Adjusted mean paired DNA-methylation M-value difference",
      FeatureLabel = as.character(FeatureLabel), GeneSymbol = NA_character_
    )
  counter <- counter + 1L
  feature_results[[counter]] <- result
  availability[[counter]] <- result %>% summarise(
    Dataset = first(Dataset), Layer = first(Layer), Tissue = first(Tissue),
    Contrast = first(Contrast), PairedSubjects = first(PairedSubjects),
    Families = first(Families), DaughterN = NA_integer_, MotherN = NA_integer_,
    FatherN = NA_integer_, TestedFeatures = n()
  )
}

feature_all <- bind_rows(feature_results)
availability_all <- bind_rows(availability) %>%
  bind_rows(tribble(
    ~Dataset, ~Layer, ~Tissue, ~Contrast, ~PairedSubjects, ~Families,
    ~DaughterN, ~MotherN, ~FatherN, ~TestedFeatures,
    "Adipose_Microarray", "Transcriptomics", "Adipose", "Post3h minus Pre unavailable",
    NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, 0L,
    "Muscle_Microarray", "Transcriptomics", "Skeletal muscle", "Post3h minus Pre unavailable",
    NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, 0L
  ))

feature_summary <- feature_all %>%
  group_by(Dataset, Layer, Tissue, Contrast) %>%
  summarise(
    PairedSubjects = first(PairedSubjects), Families = first(Families),
    TestedFeatures = n(), NominalP05 = sum(P_Value < 0.05, na.rm = TRUE),
    BH_FDR05 = sum(BH_FDR < 0.05, na.rm = TRUE),
    BH_FDR10 = sum(BH_FDR < 0.10, na.rm = TRUE),
    BH_FDR20 = sum(BH_FDR < 0.20, na.rm = TRUE),
    MinimumP = min(P_Value, na.rm = TRUE),
    MinimumBH_FDR = min(BH_FDR, na.rm = TRUE),
    MedianAbsoluteEffect = median(abs(Effect), na.rm = TRUE),
    MaximumAbsoluteEffect = max(abs(Effect), na.rm = TRUE),
    FamilyCorrelation = first(FamilyCorrelation), .groups = "drop"
  )

feature_selected <- feature_all %>%
  filter(P_Value < 0.05 | BH_FDR < 0.20) %>%
  arrange(Dataset, Contrast, BH_FDR, P_Value)

saveRDS(feature_all, file.path(feature_dir, "Exercise_Main_Effect_All_Features.rds"), compress = "xz")
write_csv(feature_summary, file.path(feature_dir, "Exercise_Main_Effect_Significance_Summary.csv"))
write_csv(feature_selected, file.path(feature_dir, "Exercise_Main_Effect_Selected_Features.csv.gz"))
write_csv(availability_all, file.path(feature_dir, "Exercise_Main_Effect_Sample_Availability.csv"))
write_xlsx(list(
  Significance_summary = feature_summary,
  Selected_features = feature_selected %>% filter(BH_FDR < 0.20 | P_Value < 1e-4),
  Sample_availability = availability_all,
  Methods = tribble(
    ~Field, ~Value,
    "Primary contrast", "Within-participant Post minus Pre molecular response",
    "Model", "Delta ~ centered adiposity group + centered Role; FamilyID duplicateCorrelation block",
    "Overall coefficient", "Role- and adiposity-adjusted grand mean response",
    "Multiplicity", "Benjamini-Hochberg within dataset and time contrast",
    "Primary threshold", "BH-FDR < 0.05",
    "Suggestive threshold", "0.05 <= BH-FDR < 0.10",
    "Exploratory threshold", "0.10 <= BH-FDR < 0.20",
    "Microarrays", "Post-exercise arrays unavailable; no longitudinal test"
  )
), file.path(feature_dir, "Exercise_Main_Effect_Audit.xlsx"))

# Full-ranked GSEA for protein-mapped platforms. Hallmark and GO:BP are gene-set
# collections; GSEA is the enrichment algorithm.
hallmark_raw <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  filter(!is.na(gs_name), !is.na(gene_symbol), nzchar(gene_symbol)) %>%
  distinct(gs_name, gene_symbol, db_version)
go_raw <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
  filter(!is.na(gs_name), !is.na(gene_symbol), nzchar(gene_symbol)) %>%
  distinct(gs_name, gene_symbol, db_version)
gene_sets <- list(
  Hallmark = split(hallmark_raw$gene_symbol, hallmark_raw$gs_name),
  GO_BP = split(go_raw$gene_symbol, go_raw$gs_name)
)
database_versions <- tibble(
  Database = c("MSigDB Hallmark", "MSigDB GO biological process"),
  Version = c(first(hallmark_raw$db_version), first(go_raw$db_version))
)

run_gsea <- function(results, collection_name, sets) {
  ranked <- results %>%
    filter(!is.na(GeneSymbol), nzchar(GeneSymbol), is.finite(ModeratedT)) %>%
    mutate(GeneSymbol = str_split(GeneSymbol, ";", simplify = TRUE)[, 1]) %>%
    group_by(GeneSymbol) %>%
    arrange(P_Value, desc(abs(ModeratedT)), FeatureID, .by_group = TRUE) %>%
    slice_head(n = 1L) %>% ungroup()
  stats <- ranked$ModeratedT
  names(stats) <- ranked$GeneSymbol
  stats <- sort(stats[is.finite(stats) & !duplicated(names(stats))], decreasing = TRUE)
  mapped <- lapply(sets, function(x) intersect(unique(x), names(stats)))
  mapped <- mapped[lengths(mapped) >= minimum_gene_set_size & lengths(mapped) <= maximum_gene_set_size]
  if (!length(mapped)) return(NULL)
  suppressWarnings(fgseaMultilevel(
    pathways = mapped, stats = stats, minSize = minimum_gene_set_size,
    maxSize = maximum_gene_set_size, eps = 0, nproc = 1
  )) %>%
    as_tibble() %>%
    transmute(
      Database = if_else(collection_name == "Hallmark", "MSigDB Hallmark",
                         "MSigDB GO biological process"),
      Method = "Full-ranked GSEA using limma moderated t statistics",
      PathwayID = pathway,
      Pathway = pathway %>% str_remove("^HALLMARK_") %>%
        str_remove("^GOBP_") %>% str_replace_all("_", " ") %>%
        str_to_lower() %>% str_to_sentence(),
      Size = size, NES = NES, P_Value = pval,
      BH_FDR = p.adjust(pval, method = "BH"),
      LeadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
    )
}

gsea_results <- list(); gsea_counter <- 0L
gsea_keys <- feature_all %>% filter(Layer == "Proteomics") %>% distinct(Dataset, Contrast)
for (key_index in seq_len(nrow(gsea_keys))) {
  dataset <- gsea_keys$Dataset[key_index]
  contrast <- gsea_keys$Contrast[key_index]
  subset <- feature_all %>% filter(Dataset == dataset, Contrast == contrast)
  for (collection in names(gene_sets)) {
    message("Running ", collection, " GSEA: ", dataset, " / ", contrast)
    enriched <- run_gsea(subset, collection, gene_sets[[collection]])
    if (!is.null(enriched)) {
      gsea_counter <- gsea_counter + 1L
      gsea_results[[gsea_counter]] <- enriched %>%
        mutate(
          Dataset = dataset, Layer = first(subset$Layer), Tissue = first(subset$Tissue),
          Contrast = contrast, PairedSubjects = first(subset$PairedSubjects),
          Families = first(subset$Families), .before = 1
        )
    }
  }
}
gsea_all <- bind_rows(gsea_results)
gsea_summary <- gsea_all %>%
  group_by(Dataset, Layer, Tissue, Contrast, Database, Method) %>%
  summarise(
    PairedSubjects = first(PairedSubjects), Families = first(Families),
    TestedPathways = n(), NominalP05 = sum(P_Value < 0.05, na.rm = TRUE),
    BH_FDR05 = sum(BH_FDR < 0.05, na.rm = TRUE),
    BH_FDR10 = sum(BH_FDR < 0.10, na.rm = TRUE),
    BH_FDR20 = sum(BH_FDR < 0.20, na.rm = TRUE),
    MinimumP = min(P_Value, na.rm = TRUE),
    MinimumBH_FDR = min(BH_FDR, na.rm = TRUE), .groups = "drop"
  )
gsea_selected <- gsea_all %>%
  filter(P_Value < 0.05 | BH_FDR < 0.20) %>%
  arrange(Database, Dataset, Contrast, BH_FDR, P_Value)

# Preserve the frozen CpG-count-adjusted methylation pathway analyses used for
# Figure 3. They are not relabelled as GSEA.
go_source <- file.path(.project_root, "fig03", "FigureS3", "FigS03D_GO_Pathway_Response_SourceData.xlsx")
hallmark_candidates <- list.files(
  file.path(.project_root, "0_RecycleBin"),
  pattern = "^Fig03_Exercise_Hallmark_Pathway_Results\\.xlsx$",
  recursive = TRUE, full.names = TRUE
)
if (!file.exists(go_source)) stop("Missing frozen GO-BP pathway source: ", go_source)
if (!length(hallmark_candidates)) stop("Missing frozen Hallmark pathway source.")
hallmark_source <- hallmark_candidates[1]
methyl_go <- read_excel(go_source, sheet = "All_GO_Pathways") %>%
  filter(Layer == "DNA methylation") %>%
  mutate(Database = "GO biological process", SourceStatus = "Frozen Figure 3 formal output")
methyl_hallmark <- read_excel(hallmark_source, sheet = "All_Pathways") %>%
  filter(Layer == "DNA methylation") %>%
  mutate(Database = "MSigDB Hallmark", SourceStatus = "Frozen Figure 3 formal output")

write_csv(gsea_summary, file.path(pathway_dir, "Proteomics_GSEA_Significance_Summary.csv"))
write_csv(gsea_selected, file.path(pathway_dir, "Proteomics_GSEA_Selected_Pathways.csv"))
write_csv(gsea_all, file.path(pathway_dir, "Proteomics_GSEA_All_Pathways.csv.gz"))
write_csv(methyl_go, file.path(pathway_dir, "Methylation_GO_BP_methylglm_All_Pathways.csv.gz"))
write_csv(methyl_hallmark, file.path(pathway_dir, "Methylation_Hallmark_methylglm_All_Pathways.csv"))
write_xlsx(list(
  GSEA_summary = gsea_summary,
  GSEA_FDR20 = gsea_selected %>% filter(BH_FDR < 0.20),
  GSEA_nominal_P05 = gsea_selected %>% filter(P_Value < 0.05),
  Methylation_Hallmark = methyl_hallmark,
  Methylation_GO_BP_FDR20 = methyl_go %>% filter(BH_FDR < 0.20),
  Database_versions = database_versions,
  Methods = tribble(
    ~Field, ~Value,
    "GSEA input", "All measured gene-mapped proteins ranked by the limma moderated t statistic",
    "GSEA collections", "MSigDB Hallmark and C5 GO biological process",
    "GSEA filtering", "No feature-level P-value prefiltering",
    "GSEA multiplicity", "BH within dataset, time contrast, and gene-set collection",
    "DNA methylation", "CpG-count-adjusted methylglm; not ordinary GSEA",
    "GO-BP role", "Secondary detailed pathway annotation because GO-BP terms are numerous and redundant",
    "Hallmark role", "Primary concise pathway collection",
    "Serum metabolomics", "Gene-based GSEA not performed because validated gene mappings are unavailable"
  )
), file.path(pathway_dir, "Pathway_Enrichment_GO_BP_Hallmark_GSEA.xlsx"))

writeLines(capture.output(sessionInfo()), file.path(feature_dir, "R_SESSION_INFO.txt"))
writeLines(capture.output(sessionInfo()), file.path(pathway_dir, "R_SESSION_INFO.txt"))
message("Completed exercise main-effect and pathway analyses in: ", archive_root)
