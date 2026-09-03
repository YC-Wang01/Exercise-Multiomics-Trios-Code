# Family aggregation of individual Hallmark proteomic response scores.
# GSVA provides subject-level pathway scores; ICC and within-family profile
# similarity test family-associated response structure. GO-BP is not used here
# because its high redundancy would create an unstable family-level test space.

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

p_load(dplyr, GSVA, matrixStats, msigdbr, readr, readxl, stringr, tibble, tidyr, writexl)
select <- dplyr::select
filter <- dplyr::filter
count <- dplyr::count

set.seed(20260828)
archive_root <- file.path(.project_root, "Pre_Post_Response_Analysis")
output_dir <- file.path(archive_root, "Results", "07_Family_Pathway_Response")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

role_levels <- c("Daughter", "Mother", "Father")
n_permutations <- 1000L
n_bootstrap <- 1000L
minimum_gene_set_size <- 10L
maximum_gene_set_size <- 500L

rank_normalize <- function(x) {
  out <- rep(NA_real_, length(x)); keep <- is.finite(x); n <- sum(keep)
  if (n >= 3L) out[keep] <- qnorm((rank(x[keep], ties.method = "average") - 0.5) / n)
  out
}

clinical <- read_csv(path_analysis_ready("Metadata", "Clinical_Observed_Primary.csv"),
                     show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    FamilyID = as.character(FamilyID), Clinical_Subject_ID = as.character(Clinical_Subject_ID),
    Role = factor(Role, levels = role_levels), BodyFat = suppressWarnings(as.numeric(Fat_p_G))
  ) %>% distinct(Clinical_Subject_ID, .keep_all = TRUE) %>%
  mutate(BodyFat_Z = rank_normalize(BodyFat))

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

datasets <- tribble(
  ~Dataset, ~Tissue, ~InputFile, ~AnnotationFile, ~GeneColumn,
  "Serum_Proteomics", "Serum", "Analysis_Serum_Proteomics.csv",
  "Feature_Dictionary_Serum_Proteomics.csv", "GeneSymbol",
  "Adipose_Proteomics", "Adipose", "Analysis_Adipose_Proteomics.csv",
  "Feature_Annotation_Adipose_Proteomics.csv", "IntegrationGene",
  "Muscle_Proteomics", "Skeletal muscle", "Analysis_Muscle_Proteomics.csv",
  "Feature_Annotation_Muscle_Proteomics.csv", "IntegrationGene"
)

hallmark_raw <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  filter(!is.na(gs_name), !is.na(gene_symbol), nzchar(gene_symbol)) %>%
  distinct(gs_name, gene_symbol, db_version)
hallmark_sets <- lapply(split(hallmark_raw$gene_symbol, hallmark_raw$gs_name), unique)

aggregate_to_gene <- function(matrix, spec) {
  annotation <- read_csv(path_analysis_ready("Metadata", spec$AnnotationFile),
                         show_col_types = FALSE, progress = FALSE)
  genes <- as.character(annotation[[spec$GeneColumn]])
  if ("PG.Genes" %in% names(annotation)) {
    fallback <- as.character(annotation$PG.Genes)
    replace <- is.na(genes) | !nzchar(trimws(genes))
    genes[replace] <- fallback[replace]
  }
  key <- tibble(FeatureID = as.character(annotation$FeatureID), GeneSymbol = genes) %>%
    mutate(GeneSymbol = str_split(GeneSymbol, ";", simplify = TRUE)[, 1]) %>%
    filter(!is.na(GeneSymbol), nzchar(GeneSymbol)) %>% distinct(FeatureID, .keep_all = TRUE)
  index <- match(rownames(matrix), key$FeatureID)
  keep <- !is.na(index)
  matrix <- matrix[keep, , drop = FALSE]
  genes <- key$GeneSymbol[index[keep]]
  summed <- rowsum(matrix, group = genes, reorder = FALSE, na.rm = TRUE)
  counts <- rowsum(is.finite(matrix) * 1, group = genes, reorder = FALSE, na.rm = TRUE)
  summed / pmax(counts, 1)
}

prepare_scores <- function(spec) {
  matrix <- read_analysis_matrix(path_analysis_ready(spec$InputFile))
  metadata <- sample_sheet %>%
    filter(Dataset == spec$Dataset, Timepoint %in% c("Pre", "Post3h"),
           SampleID %in% colnames(matrix)) %>%
    inner_join(clinical, by = c("FamilyID", "Clinical_Subject_ID", "Role")) %>%
    distinct(Clinical_Subject_ID, Timepoint, .keep_all = TRUE) %>%
    select(FamilyID, Clinical_Subject_ID, Role, BodyFat, BodyFat_Z, Timepoint, SampleID) %>%
    pivot_wider(names_from = Timepoint, values_from = SampleID) %>%
    filter(!is.na(Pre), !is.na(Post3h)) %>% arrange(FamilyID, Role, Clinical_Subject_ID)
  sample_ids <- c(metadata$Pre, metadata$Post3h)
  gene_matrix <- aggregate_to_gene(matrix[, sample_ids, drop = FALSE], spec)
  gene_matrix <- gene_matrix[matrixStats::rowSds(gene_matrix, na.rm = TRUE) > 0, , drop = FALSE]
  mapped <- lapply(hallmark_sets, function(x) intersect(x, rownames(gene_matrix)))
  mapped <- mapped[lengths(mapped) >= minimum_gene_set_size & lengths(mapped) <= maximum_gene_set_size]
  parameter <- GSVA::gsvaParam(gene_matrix, mapped, kcdf = "Gaussian",
                               minSize = minimum_gene_set_size,
                               maxSize = maximum_gene_set_size)
  scores <- GSVA::gsva(parameter, verbose = FALSE)
  pre_scores <- scores[, match(metadata$Pre, colnames(scores)), drop = FALSE]
  post_scores <- scores[, match(metadata$Post3h, colnames(scores)), drop = FALSE]
  delta <- t(post_scores - pre_scores)
  rownames(delta) <- metadata$Clinical_Subject_ID
  colnames(delta) <- str_remove(rownames(scores), "^HALLMARK_") %>%
    str_replace_all("_", " ") %>% str_to_lower() %>% str_to_sentence()
  list(delta = delta, metadata = metadata, measured_genes = nrow(gene_matrix),
       scored_pathways = ncol(delta))
}

residualize_pathway <- function(score, metadata) {
  data <- data.frame(Score = score, Role = metadata$Role, BodyFat_Z = metadata$BodyFat_Z)
  keep <- complete.cases(data)
  residual <- rep(NA_real_, nrow(data))
  if (sum(keep) >= 8L) residual[keep] <- residuals(lm(Score ~ Role + BodyFat_Z, data = data[keep, ]))
  residual
}

balanced_icc <- function(values, family, role) {
  data <- tibble(Value = values, FamilyID = family, Role = role) %>% filter(is.finite(Value))
  counts <- data %>% distinct(FamilyID, Role) %>% count(FamilyID, name = "Roles")
  complete <- counts %>% filter(Roles == 3L) %>% pull(FamilyID)
  data <- data %>% filter(FamilyID %in% complete)
  if (n_distinct(data$FamilyID) < 5L) return(NA_real_)
  grand <- mean(data$Value)
  family_means <- data %>% group_by(FamilyID) %>% summarise(Mean = mean(Value), .groups = "drop")
  k <- 3; n_family <- nrow(family_means)
  ms_between <- k * sum((family_means$Mean - grand)^2) / (n_family - 1)
  joined <- data %>% left_join(family_means, by = "FamilyID")
  ms_within <- sum((joined$Value - joined$Mean)^2) / (n_family * (k - 1))
  (ms_between - ms_within) / (ms_between + (k - 1) * ms_within)
}

permuted_family <- function(metadata) {
  output <- as.character(metadata$FamilyID)
  for (role in role_levels) {
    index <- which(metadata$Role == role)
    output[index] <- sample(output[index])
  }
  output
}

icc_one <- function(pathway, scores, metadata) {
  residual <- residualize_pathway(scores, metadata)
  complete_families <- metadata %>% distinct(FamilyID, Role) %>% count(FamilyID, name = "Roles") %>%
    filter(Roles == 3L) %>% pull(FamilyID)
  if (length(complete_families) < 5L) return(NULL)
  trio_matrix <- vapply(role_levels, function(role) {
    role_meta <- metadata %>% filter(Role == role, FamilyID %in% complete_families) %>% arrange(FamilyID)
    residual[match(role_meta$Clinical_Subject_ID, metadata$Clinical_Subject_ID)]
  }, numeric(length(complete_families)))
  rownames(trio_matrix) <- sort(complete_families)
  colnames(trio_matrix) <- role_levels
  if (any(!is.finite(trio_matrix))) return(NULL)
  matrix_icc <- function(value_matrix) {
    k <- ncol(value_matrix); n_family <- nrow(value_matrix)
    grand <- mean(value_matrix)
    family_means <- rowMeans(value_matrix)
    ms_between <- k * sum((family_means - grand)^2) / (n_family - 1)
    ms_within <- sum((value_matrix - family_means)^2) / (n_family * (k - 1))
    (ms_between - ms_within) / (ms_between + (k - 1) * ms_within)
  }
  raw_icc <- matrix_icc(trio_matrix)
  permutation <- replicate(n_permutations, {
    permuted <- trio_matrix
    for (column in seq_len(ncol(permuted))) permuted[, column] <- sample(permuted[, column])
    matrix_icc(permuted)
  })
  p_value <- (1 + sum(permutation >= raw_icc)) / (n_permutations + 1)
  bootstrap <- replicate(n_bootstrap, {
    sampled <- sample.int(nrow(trio_matrix), nrow(trio_matrix), replace = TRUE)
    matrix_icc(trio_matrix[sampled, , drop = FALSE])
  })
  tibble(
    Pathway = pathway, Subjects = sum(is.finite(residual)),
    CompleteTrioFamilies = length(complete_families),
    ICC_Raw = raw_icc, FamilyICC = pmax(raw_icc, 0),
    Lower95 = pmax(quantile(bootstrap, 0.025, na.rm = TRUE), 0),
    Upper95 = pmax(quantile(bootstrap, 0.975, na.rm = TRUE), 0),
    PermutationP = p_value
  )
}

pair_specs <- tribble(
  ~PairType, ~RoleA, ~RoleB,
  "Daughter-mother", "Daughter", "Mother",
  "Daughter-father", "Daughter", "Father",
  "Mother-father", "Mother", "Father"
)

pairwise_pathway_correlations <- function(delta, metadata, spec) {
  rows <- list(); counter <- 0L
  for (i in seq_len(nrow(pair_specs))) {
    pair <- pair_specs[i, ]
    a <- metadata %>% filter(Role == pair$RoleA)
    b <- metadata %>% filter(Role == pair$RoleB)
    common <- intersect(a$FamilyID, b$FamilyID)
    a <- a %>% filter(FamilyID %in% common) %>% arrange(FamilyID)
    b <- b %>% filter(FamilyID %in% common) %>% arrange(FamilyID)
    for (pathway in colnames(delta)) {
      x <- delta[match(a$Clinical_Subject_ID, rownames(delta)), pathway]
      y <- delta[match(b$Clinical_Subject_ID, rownames(delta)), pathway]
      keep <- is.finite(x) & is.finite(y)
      test <- if (sum(keep) >= 5L && sd(x[keep]) > 0 && sd(y[keep]) > 0) {
        suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
      } else NULL
      counter <- counter + 1L
      rows[[counter]] <- tibble(
        Dataset = spec$Dataset, Tissue = spec$Tissue, PairType = pair$PairType,
        Pathway = pathway, FamilyPairs = sum(keep),
        SpearmanRho = if (is.null(test)) NA_real_ else unname(test$estimate),
        P_Value = if (is.null(test)) NA_real_ else test$p.value
      )
    }
  }
  bind_rows(rows) %>% group_by(Dataset, PairType) %>%
    mutate(BH_FDR = p.adjust(P_Value, method = "BH")) %>% ungroup()
}

profile_similarity_test <- function(delta, metadata, spec, pair) {
  a <- metadata %>% filter(Role == pair$RoleA)
  b <- metadata %>% filter(Role == pair$RoleB)
  common <- intersect(a$FamilyID, b$FamilyID)
  a <- a %>% filter(FamilyID %in% common) %>% arrange(FamilyID)
  b <- b %>% filter(FamilyID %in% common) %>% arrange(FamilyID)
  matrix_a <- delta[match(a$Clinical_Subject_ID, rownames(delta)), , drop = FALSE]
  matrix_b <- delta[match(b$Clinical_Subject_ID, rownames(delta)), , drop = FALSE]
  similarity <- outer(seq_len(nrow(matrix_a)), seq_len(nrow(matrix_b)), Vectorize(function(i, j) {
    keep <- is.finite(matrix_a[i, ]) & is.finite(matrix_b[j, ])
    if (sum(keep) < 10L) return(NA_real_)
    suppressWarnings(cor(matrix_a[i, keep], matrix_b[j, keep], method = "spearman"))
  }))
  familial <- diag(similarity)
  unrelated_by_row <- vapply(seq_len(nrow(similarity)), function(i) {
    mean(similarity[i, -i], na.rm = TRUE)
  }, numeric(1))
  contribution <- familial - unrelated_by_row
  estimate <- mean(contribution, na.rm = TRUE)
  permutation <- replicate(n_permutations, {
    perm <- sample.int(nrow(similarity))
    mean(similarity[cbind(seq_len(nrow(similarity)), perm)] - unrelated_by_row, na.rm = TRUE)
  })
  bootstrap <- replicate(n_bootstrap, mean(sample(contribution, replace = TRUE), na.rm = TRUE))
  tibble(
    Dataset = spec$Dataset, Tissue = spec$Tissue, PairType = pair$PairType,
    FamilyPairs = length(familial), PathwayDimensions = ncol(delta),
    MeanFamilialRho = mean(familial, na.rm = TRUE),
    MeanRoleMatchedUnrelatedRho = mean(similarity[row(similarity) != col(similarity)], na.rm = TRUE),
    FamilialMinusUnrelatedRho = estimate,
    Lower95 = quantile(bootstrap, 0.025, na.rm = TRUE),
    Upper95 = quantile(bootstrap, 0.975, na.rm = TRUE),
    PermutationP = (1 + sum(abs(permutation) >= abs(estimate))) / (n_permutations + 1)
  )
}

icc_results <- list(); pathway_pair_results <- list(); profile_results <- list()
score_tables <- list(); availability <- list()
for (i in seq_len(nrow(datasets))) {
  spec <- datasets[i, ]
  message("Calculating Hallmark GSVA response scores: ", spec$Dataset)
  prepared <- prepare_scores(spec)
  score_tables[[spec$Dataset]] <- as_tibble(prepared$delta, rownames = "Clinical_Subject_ID") %>%
    left_join(prepared$metadata %>% select(Clinical_Subject_ID, FamilyID, Role, BodyFat, BodyFat_Z),
              by = "Clinical_Subject_ID")
  availability[[spec$Dataset]] <- tibble(
    Dataset = spec$Dataset, Tissue = spec$Tissue,
    PairedSubjects = nrow(prepared$metadata), Families = n_distinct(prepared$metadata$FamilyID),
    CompleteTrioFamilies = prepared$metadata %>% distinct(FamilyID, Role) %>%
      count(FamilyID, name = "Roles") %>% summarise(N = sum(Roles == 3L)) %>% pull(N),
    MeasuredGenes = prepared$measured_genes, ScoredHallmarkPathways = prepared$scored_pathways
  )
  icc <- bind_rows(lapply(colnames(prepared$delta), function(pathway) {
    icc_one(pathway, prepared$delta[, pathway], prepared$metadata)
  }))
  if (nrow(icc)) {
    icc <- icc %>%
      mutate(Dataset = spec$Dataset, Tissue = spec$Tissue, .before = 1) %>%
      mutate(BH_FDR = p.adjust(PermutationP, method = "BH"))
  } else {
    icc <- tibble(
      Dataset = character(), Tissue = character(), Pathway = character(),
      Subjects = integer(), CompleteTrioFamilies = integer(), ICC_Raw = double(),
      FamilyICC = double(), Lower95 = double(), Upper95 = double(),
      PermutationP = double(), BH_FDR = double()
    )
  }
  icc_results[[spec$Dataset]] <- icc
  pathway_pair_results[[spec$Dataset]] <- pairwise_pathway_correlations(
    prepared$delta, prepared$metadata, spec
  )
  profile_results[[spec$Dataset]] <- bind_rows(lapply(seq_len(nrow(pair_specs)), function(j) {
    profile_similarity_test(prepared$delta, prepared$metadata, spec, pair_specs[j, ])
  })) %>% mutate(BH_FDR = p.adjust(PermutationP, method = "BH"))
  rm(prepared); invisible(gc())
}

icc_all <- bind_rows(icc_results) %>% arrange(Dataset, BH_FDR, PermutationP)
pathway_pairs_all <- bind_rows(pathway_pair_results) %>% arrange(Dataset, PairType, BH_FDR, P_Value)
profiles_all <- bind_rows(profile_results) %>% arrange(Dataset, BH_FDR, PermutationP)
availability_all <- bind_rows(availability)

# Reuse the frozen serum-metabolite module family analysis, including LOFO checks.
serum_family_source <- file.path(.project_root, "fig03", "FigureS3",
                                 "FigS03F_Family_Module_Stability_SourceData.xlsx")
if (!file.exists(serum_family_source)) stop("Missing frozen serum family source: ", serum_family_source)
serum_icc <- read_excel(serum_family_source, sheet = "Primary Family ICC")
serum_pairs <- read_excel(serum_family_source, sheet = "Primary Pair Tests")
serum_lofo_icc <- read_excel(serum_family_source, sheet = "LOFO ICC Summary")
serum_lofo_pairs <- read_excel(serum_family_source, sheet = "LOFO Pair Summary")

write_csv(icc_all, file.path(output_dir, "Proteomic_Hallmark_Response_Family_ICC.csv"))
write_csv(pathway_pairs_all, file.path(output_dir, "Proteomic_Hallmark_Response_DM_DF_MF_Pathway_Correlations.csv"))
write_csv(profiles_all, file.path(output_dir, "Proteomic_Hallmark_Response_Whole_Profile_Similarity.csv"))
write_csv(availability_all, file.path(output_dir, "Family_Pathway_Response_Sample_Availability.csv"))
saveRDS(score_tables, file.path(output_dir, "Proteomic_Hallmark_Response_GSVA_Scores.rds"), compress = "xz")
write_xlsx(list(
  Proteomic_Hallmark_ICC = icc_all,
  Proteomic_pathway_pairs = pathway_pairs_all,
  Proteomic_profile_similarity = profiles_all,
  Serum_metabolite_ICC = serum_icc,
  Serum_metabolite_pairs = serum_pairs,
  Serum_metabolite_LOFO_ICC = serum_lofo_icc,
  Serum_metabolite_LOFO_pairs = serum_lofo_pairs,
  Sample_availability = availability_all,
  Methods = tribble(
    ~Field, ~Value,
    "Proteomic pathway scores", "GSVA on gene-aggregated Pre and Post3h proteomic matrices; response score = Post3h minus Pre",
    "Gene aggregation", "Mean abundance across protein features mapped to the same primary gene symbol",
    "ICC covariates", "Pathway response residualized for Role and rank-normalized body-fat percentage",
    "Family ICC", "Balanced one-way ICC estimated in complete daughter-mother-father trios",
    "ICC inference", "One-sided within-role family-label permutation; BH within dataset across Hallmark pathways",
    "Pairwise pathway tests", "Spearman correlation across families for daughter-mother, daughter-father, and mother-father pairs",
    "Whole-profile similarity", "Matched-family minus role-matched unrelated Spearman similarity across Hallmark response profiles",
    "Serum metabolomics", "Frozen four-module Figure S3 family analysis with permutation, bootstrap, and leave-one-family-out checks",
    "GO-BP", "Not used for family ICC because thousands of redundant GO-BP terms create an unstable low-power test space",
    "DNA methylation", "Not scored by GSVA; CpG-count-adjusted methylglm is retained in the pathway-enrichment module"
  )
), file.path(output_dir, "Family_Pathway_Response_GSVA_ICC.xlsx"))

writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_SESSION_INFO.txt"))
message("Completed family pathway-response analysis: ", output_dir)
