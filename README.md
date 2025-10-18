# JPP-RPMGRS: Reproducible Analysis Pipeline

This repository contains a reproducible R pipeline for comparing two rating methods across conditions (RPM vs. GRS). It covers reliability (ICC), construct validity (CFA), leniency (elevation and discrepancy), and convergent validity (correlations), and produces cleaned data and publication-ready tables.

## Contents

- Project overview and requirements
- Data expectations
- Quick start and how to run
- Repository structure
- Outputs produced
- Methods per step
- Troubleshooting and reproducibility notes

---

## Overview

- Entry point: `main.R` orchestrates the pipeline (Steps 02–06) and logs runtime.
- Steps:
  1) Load and clean data; derive composites and discrepancies
  2) ICC(1) by variable/cluster/condition
  3) Confirmatory Factor Analyses (by rater set and condition)
  4) Leniency analyses (elevation and discrepancy) with achieved power
  5) Correlations by condition, Steiger Z-tests, and power

## Requirements

- R 4.2+ (4.3+ recommended)
- Toolchain for compiled packages (needed by `lme4`, `nlme`, `glmmTMB`, `lavaan`):
  - macOS: Xcode Command Line Tools
  - Windows: Rtools
  - Linux: build-essential/compilers (e.g., `g++`, `make`)
- Internet access to install CRAN packages on first run

R packages (installed automatically when running scripts):

```
here, haven, dplyr, tidyr, lme4, lmerTest, glmmTMB, nlme, lavaan,
correlation, pwr, knitr, tibble, rlang
```

> Note: Packages are installed on demand via helper utilities; see `R/01_helpers.R`.

## Data Expectations

Place the raw SPSS dataset at `data/data.sav`.

Required columns:

- IDs/clusters:
  - `GroupID` (renamed to `Group` in cleaning)
  - `Experimenter`
  - `Condition` coded as 1/2 (recoded to `RPM`/`GRS`, reference = `RPM`)
- Rating dimensions with prefixes (used for composites and CFA):
  - `Self_*`, `Peer_*`, `Super_*`
  - CFA indicators used explicitly:
    - `Super_Organization`, `Super_Physical`, `Super_Visual`, `Super_Vocal`
    - `Peer_Organization`, `Peer_Physical`, `Peer_Visual`, `Peer_Vocal`

Derived during cleaning (`R/02_load_clean.R`):

- Composite means: `Self_Mean`, `Peer_Mean`, `Super_Mean`
- Discrepancies: `Peer_Discrep = Peer_Mean - Super_Mean`, `Self_Discrep = Self_Mean - Super_Mean`

## Quick Start

1) Open the R project `JPP Submission.Rproj` (recommended) or set working directory to the repo root.
2) Ensure the input file exists: `data/data.sav`.
3) Create output folders (first run only):

```r
dir.create("output", showWarnings = FALSE)
dir.create(file.path("output", "tables"), showWarnings = FALSE, recursive = TRUE)
```

4) Run the full pipeline:

```r
source("main.R")
```

Or from a shell:

```sh
Rscript main.R
```

## Repository Structure

```
.
├── main.R                  # Master script to run all analyses
├── R/
│   ├── 01_helpers.R        # Helpers: package loading, ICC, model & power utils, display
│   ├── 02_load_clean.R     # Load SPSS, clean, derive composites/discrepancies
│   ├── 03_icc_analysis.R   # ICC(1) by variable/cluster/condition
│   ├── 04_cfa_analysis.R   # Single-factor CFA per rater set & condition
│   ├── 05_leniency.R       # Elevation/discrepancy analyses and achieved power
│   └── 06_correlations.R   # Correlations by condition; Steiger Z comparisons and power
├── data/
│   └── data.sav            # Raw SPSS data (input)
├── output/                 # Generated outputs (created at runtime)
│   └── tables/             # CSV tables saved by steps 03–06
└── JPP Submission.Rproj    # RStudio project file
```

## Outputs

- Clean data: `output/cleaned_data.rds`
- Tables (CSV) under `output/tables/`:
  - `03_icc_summary.csv`
  - `04_cfa_fit_summary.csv`
  - `04_cfa_loadings_summary.csv`
  - `05_elevation_leniency_summary.csv`
  - `05_discrepancy_leniency_summary.csv`
  - `05_elevation_power_summary.csv`
  - `05_discrepancy_power_summary.csv`
  - `06_correlation_summary.csv`

## Methods Summary

- Step 02 (Load/Clean): Reads SPSS, recodes `Condition` (1/2 → RPM/GRS), converts clusters to factors, derives composites and discrepancies, saves `cleaned_data.rds`.
- Step 03 (ICC): Computes ICC(1) for `Self_Mean`, `Peer_Mean`, `Super_Mean` across clusters (`Group`, `Experimenter`) separately by `Condition`.
- Step 04 (CFA): Single-factor CFA for Supervisor and Peer indicators under `RPM` and `GRS`. Chooses a cluster for multilevel CFA when max ICC > 0.05.
- Step 05 (Leniency):
  - Elevation: tests effect of `Condition` on composite means, selecting random effects via ICC-guided candidates (engines: `glmmTMB`, `nlme`, fallback `lm`).
  - Discrepancy: intercept-only tests of `Self_Discrep` and `Peer_Discrep` within each condition.
  - Achieved power: computed for elevation and discrepancy outcomes.
- Step 06 (Correlations): Multilevel correlations for Self–Peer, Self–Super, Peer–Super by condition; Steiger Z test for differences between conditions; achieved power per correlation.

## Reproducibility Notes

- Paths use `here::here()` for portability; no `setwd()` needed.
- Scripts attempt package installation automatically on first run.
- Some models use iterative optimizers; minor numerical variation is expected across runs/platforms.
- Missing data are handled via `na.omit`/`drop_na` or `lavaan`'s `missing = "pairwise"`.

## Troubleshooting

- Missing output folders: Create `output/` and `output/tables/` before running (see Quick Start).
- Package compilation errors (`glmmTMB`, `nlme`, `lavaan`): ensure compilers/toolchains are installed for your OS.
- Data column mismatches: Verify variable names and `Condition` coding match “Data Expectations”.
- Convergence/singularity: The helpers reduce random-effects complexity automatically; check console messages to see which model/engine was used.
- Not found: `%||%` operator: if you encounter an error about `%||%`, ensure `rlang` is installed and attached, or replace with a standard null-coalescing check in code.

## Citation and License

TBD

## Maintainers

Justin R Feeney, PhD, SHRM-SCP
Associate Professor of Management
Department of Business Administration and Economics
jfeeney@saintmarys.edu | 574-284-4488
