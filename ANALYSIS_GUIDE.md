# Analysis guide

## Purpose

This document maps the public scripts to the analyses supporting Figures 1-6
and Figures S1-S6. The repository documents the statistical models, pathway
analyses, network construction, and scientific plotting workflow. Exact
numerical reproduction requires authorized access to the controlled study
inputs.

## Configuration

Run scripts from the repository root. Set `ATM_PROJECT_ROOT` to an authorized
analysis workspace before sourcing the shared configuration:

```r
Sys.setenv(ATM_PROJECT_ROOT = "/path/to/CellMetabolism_Transfer")
source("code/Fig00_Config.R")
```

The configuration resolves controlled source and derived-data locations
outside the Git repository. Raw and source-of-truth files must remain
immutable and outside version control.

## Computational environment

The release was prepared and statically checked with R 4.5.1 on Windows. The
resolved R package set is recorded in `renv.lock`. Principal packages include
`limma`, `lme4`, `lmerTest`, `emmeans`, `fgsea`, `GSVA`, `msigdbr`,
`methylGSA`, `ComplexHeatmap`, `igraph`, `ggraph`, `ggplot2`, `dplyr`,
`tidyr`, `readr`, and `readxl`.

All 42 public R scripts parse successfully with R 4.5.1. A complete numerical
rerun is not included because participant-level inputs remain controlled.

## Shared analyses

Run only the shared modules required for the target figure, in the following
order:

1. `code/Fig00_Config.R`
2. `code/Analysis_Acute_Exercise_Main_Effects_and_Pathways.R`
3. `code/Analysis_Clinical_Phenotype_Omics_Associations_NonMethylation.R`
4. `code/Analysis_Clinical_Phenotype_Omics_Associations_Methylation.R`
5. `code/Analysis_Family_Pathway_Response_GSVA_ICC.R`
6. `code/Analysis_Clinical_Marker_State_Associations_NonMethylation.R`
7. `code/Analysis_Clinical_Marker_State_Associations_Methylation.R`
8. `code/Analysis_Osteocalcin_Four_Timepoint_Longitudinal.R`
9. `code/Analysis_Osteocalcin_GO_BP_GSEA.R`

These modules fit the shared acute-response, clinical-association,
methylation, pathway, family-similarity, and osteocalcin models reused by
multiple figure workflows.

## Figure 1

1. `code/Figure_1_Baseline_Molecular_Variability.R`
2. `code/Figure_S1_MultiOmics_PCA_and_Clinical_Associations.R`

These scripts summarize baseline molecular variability, principal-component
structure, and clinical cohort comparisons.

## Figure 2

1. `code/Figure_2_Baseline_Adiposity_Analysis.R`
2. `code/Figure_2_Baseline_Adiposity_MultiOmics.R`
3. `code/Figure_2_Adiposity_Pathway_Enrichment.R`
4. `code/Figure_2_Adiposity_MultiOmics_Network.R`
5. `code/Figure_S2_Family_Molecular_Network.R`
6. `code/Figure_S2_PIEZO1_Family_Analysis.R`

These scripts fit obese-versus-lean models across the eight omics platforms,
perform adiposity-associated pathway analyses, construct the adiposity and
family-similarity networks, and evaluate PIEZO1-linked methylation results.

## Figure 3

1. `code/Figure_3_Lifestyle_Molecular_Associations.R`
2. `code/Figure_3_Serum_Response_and_Cluster_Analysis.R`
3. `code/Figure_3_Serum_Trajectory_Data_Preparation.R`
4. `code/Figure_3_Serum_Molecular_Trajectories.R`
5. `code/Figure_3_Adiposity_Response_Network.R`
6. `code/Figure_3_Exercise_Response_Pathway_Analysis.R`
7. `code/Figure_3_Adiposity_Pathway_Associations.R`
8. `code/Figure_S3_Pathway_and_Family_Response_Analyses.R`

These scripts estimate acute serum responses, derive and display trajectory
clusters, relate molecular responses to lifestyle and adiposity, and generate
the response-network and pathway analyses.

## Figure 4

1. `code/Figure_4_Clinical_Insulin_Associations.R`
2. `code/Figure_4_Insulin_Pathway_Associations.R`
3. `code/Figure_4_Core_Hallmark_GSEA_Profiles.R`
4. `code/Figure_4_Post3h_Molecular_Associations.R`
5. `code/Figure_S4_Molecular_Association_Scatterplots.R`
6. `code/Figure_S4_Supplementary_Analyses_and_Plots.R`

These scripts test clinical associations with insulin-related phenotypes,
perform pre-exercise, Post3h, and response pathway analyses for HOMA-IR and
Matsuda ISI, and display selected molecular associations.

## Figure 5

1. `code/Figure_5_hsCRP_Clinical_Associations.R`
2. `code/Figure_5_hsCRP_Adipose_Network.R`
3. `code/Figure_5_hsCRP_Pathway_and_Molecular_Associations.R`

These scripts test clinical, adipose-proteomic, pathway, and molecular
associations with pre-exercise hsCRP.

## Figure 6

1. `code/Figure_6_Osteocalcin_Integrated_MultiOmics_Analysis.R`
2. `code/Figure_6_Osteocalcin_Clinical_Associations.R`
3. `code/Figure_6_Osteocalcin_and_Clinical_Composite.R`
4. `code/Figure_6_Osteocalcin_Pathway_Associations.R`
5. `code/Figure_6_Osteocalcin_Modular_Network.R`
6. `code/Figure_S6_Osteocalcin_Family_Role_Sex_Analysis.R`
7. `code/Figure_S6_Parental_Role_Associations.R`
8. `code/Figure_S6_Core_Hallmark_GSEA_Profiles.R`

These scripts analyze four-timepoint osteocalcin trajectories, clinical and
multi-omic associations, Hallmark pathway enrichment, modular molecular
networks, and family-, role-, and sex-associated patterns.

## Statistical interpretation

Benjamini-Hochberg correction is applied within the prespecified analysis
families defined by the scripts and manuscript. Exploratory thresholds remain
explicitly labeled. Network geometry is descriptive and does not establish a
physical interaction or causal relationship. Positive or negative GSEA
normalized enrichment scores describe enrichment within a ranked association
statistic and do not, by themselves, demonstrate pathway activation or
inhibition.
