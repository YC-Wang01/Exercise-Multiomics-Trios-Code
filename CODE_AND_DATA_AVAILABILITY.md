# Code and data availability

## Public code

This repository contains the curated R analysis and scientific plotting code
supporting Figures 1-6 and Figures S1-S6 of the associated manuscript. The
planned permanent GitHub location is:

`https://github.com/YC-Wang01/Exercise-Multiomics-Trios-Code`

The repository is organized around the scientific analyses reported in the
manuscript. `ANALYSIS_GUIDE.md` provides the figure-specific execution order
and links each script to its analytical purpose.

## Data access

Individual-level clinical and omics data are not distributed in this public
repository because they contain potentially sensitive human-participant
information and are governed by consent, ethics, institutional, legal, and
data-use restrictions. Qualified researchers may request controlled access
from the corresponding author. Requests are subject to the applicable
institutional approvals and data-use agreements.

No raw data, participant-level analysis matrices, participant identifiers,
family crosswalks, or private correspondence are included in this repository.
The expected controlled-data boundary and configuration are described in
`ANALYSIS_GUIDE.md` and `code/Fig00_Config.R`.

## Computational environment

The code was curated and statically checked with R 4.5.1. `renv.lock` records
the resolved R dependency set available during release preparation. All public
R scripts parse successfully, and the documented execution order is provided
in `ANALYSIS_GUIDE.md`.

A full numerical rerun is not possible from the public repository alone
because the participant-level inputs remain controlled. The released code is
therefore intended to document the analysis logic, statistical models,
multiple-testing procedures, pathway analyses, network construction, and
figure-generation workflow. Exact numerical reproduction requires authorized
access to the controlled inputs and a compatible software environment.

## Scope exclusions

Scripts used solely for spreadsheet styling, Word or PDF assembly,
supplementary-table delivery formatting, and internal quality-control
packaging are not included. These utilities do not define the reported
statistical analyses. Figure 7 is not included because the final manuscript
uses a graphical abstract rather than a quantitative Figure 7 analysis.

## Release integrity

Before release, the repository was screened for raw or participant-level data,
direct identifiers, identifier crosswalks, machine-specific absolute paths,
credentials, private keys, and inappropriate large files. No findings were
detected. Internal release-audit utilities and reports are retained locally
and are not part of the public scientific-code package.
