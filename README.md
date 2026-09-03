# Exercise Multi-omics in Family Trios

This code-only repository contains the curated statistical analysis and
scientific plotting scripts supporting Figures 1-6 and Figures S1-S6 of the
associated Cell Metabolism manuscript.

## Release status

This repository contains no raw data, analysis-ready participant-level
matrices, direct identifiers, identifier crosswalks, or private
correspondence. Full numerical reproduction requires authorized access to the
controlled clinical and omics inputs. See `CODE_AND_DATA_AVAILABILITY.md` for
the public-code scope and controlled-data access conditions.

## Repository structure

```text
code/                         Scientific analysis and plotting scripts
ANALYSIS_GUIDE.md             Script map, execution order, and environment
CODE_AND_DATA_AVAILABILITY.md Public-code scope and controlled-data statement
renv.lock                     Resolved R dependency record
```

`ANALYSIS_GUIDE.md` provides the script-to-analysis mapping and figure-specific
execution order. Scripts used only for spreadsheet styling, submission-table
formatting, delivery-file assembly, document composition, or internal release
auditing are intentionally excluded.

## Configuration

Run scripts from the repository root. Set `ATM_PROJECT_ROOT` to an authorized
local analysis workspace before execution:

```r
Sys.setenv(ATM_PROJECT_ROOT = "/path/to/CellMetabolism_Transfer")
source("code/Fig00_Config.R")
```

The configuration resolves controlled inputs outside this repository and
writes derived results to the authorized project workspace. Raw inputs must
remain outside version control and must not be modified.

## Execution

1. Restore the packages recorded in `renv.lock`, or install compatible package
   versions listed in `ANALYSIS_GUIDE.md`.
2. Set `ATM_PROJECT_ROOT` in a clean R session.
3. Read the shared-analysis order in `ANALYSIS_GUIDE.md`.
4. Follow the relevant figure-specific section in `ANALYSIS_GUIDE.md`.
5. Compare generated statistics and plots with the manuscript release.

Scripts remain modular because analysis, network construction, pathway
enrichment and plotting have distinct inputs and validation requirements. All
figure scripts are stored directly in `code/`; public filenames identify the
manuscript figure and scientific task. Historical output names inside scripts
are retained where needed for provenance and comparison with submitted
artifacts.

## Statistical scope

The repository retains the models and display thresholds used for the reported
analyses. Exploratory thresholds and network analyses must remain identified as
exploratory. Network geometry is descriptive and does not establish physical
interaction or causality. A positive or negative GSEA normalized enrichment
score denotes enrichment within the ranked association statistic and does not,
by itself, establish pathway activation or inhibition.

## Data availability

Raw and individual-level source data are not distributed because they contain
potentially sensitive participant information and are governed by consent,
ethics, institutional, legal and data-use restrictions. Qualified researchers
may request access from the corresponding author, subject to the applicable
approvals and agreements.

## Release validation

All 42 public R scripts were parsed successfully with R 4.5.1 and were checked
for machine-specific paths and restricted-data references before release. A
full numerical rerun requires authorized controlled inputs and is not possible
from the public repository alone.

## Citation and license

Citation metadata are provided in `CITATION.cff`. The code is distributed under
the MIT License; see `LICENSE`. Repository metadata and the public-release
decision should be confirmed by the corresponding authors before the first
remote push.
