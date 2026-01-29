# IPCW_RWD

Inverse probability of censoring weighting (IPCW) and inverse probability of treatment weighting (IPTW) workflow for OHDSI studies using real-world data. This repository implements a full pipeline from cohort generation through weighted outcome models and bootstrap inference, illustrated on the **LEGEND-HTN** study (thiazide diuretics vs. ACE inhibitors, outcome: acute myocardial infarction).

## Overview

The workflow combines:

- **IPTW** — propensity score–based weights to balance confounding.
- **IPCW** — weights from a Cox model for censoring to address informative censoring.
- **Combined weights** — product of (truncated) IPTW and stabilized IPCW for confounding and censoring adjustment.
- **m-of-n bootstrap** — for hazard ratio (HR), restricted mean survival time (RMST), and survival curve confidence intervals.

## Study setup (LEGEND-HTN)

- **Target exposure:** Thiazide diuretics (Atlas ID 1788868)
- **Comparator:** ACE inhibitors (Atlas ID 1788867)
- **Outcome:** Acute myocardial infarction (AMI, Atlas ID 1788866)

Cohort definitions live in the `cohorts/` folder. You can skip cohort generation if you already have LEGEND-HTN CohortMethod data.

## Workflow

Scripts are intended to be run in order. A single entry point is `workflow_for_all_codes.R` (see **How to run** below).

### 1. Cohort and data

| Script | Purpose |
|--------|--------|
| `generate_cohort_script.R` | Builds LEGEND-HTN cohorts and exports `CohortMethodData` (target, comparator, outcome). Uses covariate exclusion for thiazides and ACEi. Produces `results/cohortMethodData_t1788868_c1788867_o1788866.zip` (and an “all covar” version for IPCW). |

### 2. Propensity scores and IPTW

| Script | Purpose |
|--------|--------|
| `find_ps_legend_script.R` | Fits propensity score model (LSPS), creates IPTW, and fits unadjusted and IPTW-weighted Cox outcome models. Saves `results/ps_study.rds`. |

### 3. Censoring model (IPCW)

| Script | Purpose |
|--------|--------|
| `censoring_vars_legend_script.R` | Builds censoring outcome from CohortMethod data and fits a penalized Cox model for censoring risk (Cyclops). Requires the “all covar” CohortMethodData. Saves `results/Cox_censoring.rds`. |
| `find_kZ_script.R` | Computes IPCW weights (stabilized and unstabilized) using the censoring model, aligns with PS cohort, and outputs long-format survival data with weights (e.g. `survival_weights_endObsDate.csv`). |

### 4. Combined weights and final models

| Script | Purpose |
|--------|--------|
| `weights_script.R` | Truncates IPCW and IPTW (e.g. at 1st/99th percentiles), forms combined weight, and fits Cox models: unadjusted, IPCW-only, IPTW-only, and combined (IPTW + IPCW). |

### 5. Bootstrap and inference

| Script | Purpose |
|--------|--------|
| `bootstrap_m_of_n.R` | m-of-n bootstrap: resamples, re-runs PS → censoring model → IPCW → weights → Cox models. Produces bootstrap HR, RMST, and survival curve summaries in `bootstrap_results/`. |
| `hr_point_est_with_bootstrap_ci.R` | Combines point-estimate HR (from `weights_script.R`) with bootstrap HR CIs. |
| `rmst_point_est_with_bootstrap_ci.R` | Point-estimate RMST with bootstrap CIs. |
| `plot_hr_bootstrap_results.R` | Plots HR point estimates and bootstrap CIs. |
| `plot_rmst_bootstrap_results.R` | Plots RMST and bootstrap CIs. |
| `plot_bootstrap_survival_curves.R` | Plots bootstrap-based survival curves. |

### Other

| Script | Purpose |
|--------|--------|
| `study_info.R` | Study population, follow-up distribution, attrition, balance (Table 1), and covariate summaries. |
| `workflow_for_all_codes.R` | Sources the main analysis scripts in sequence (cohort → PS → censoring → IPCW → weights). Bootstrap and plotting are run separately. |

## Prerequisites

- **R** (e.g. 4.4.x; see `renv.lock`).
- **Database:** Cohort generation expects a configured database connection (e.g. SQL Server) and schema names set in `generate_cohort_script.R`. If you already have CohortMethod data, you can start from `find_ps_legend_script.R` and use the “all covar” zip for IPCW scripts.
- **OHDSI packages:** e.g. `DatabaseConnector`, `CohortMethod`, `CohortGenerator`, `ROhdsiWebApi`, `Cyclops`, plus `dplyr`, `survival`, `data.table`, etc. The project uses **renv** for reproducible dependencies.

## Setup

```r
# Restore packages from lockfile
renv::restore()
```

Configure database and paths in `generate_cohort_script.R` if you are generating cohorts from scratch.

## How to run

**Full analysis (cohort → PS → censoring → IPCW → weights):**

```r
source("workflow_for_all_codes.R")
```

Set `MAX_THREADS` in that file as needed. Ensure `results/` exists and that you have either:

- Generated data via `generate_cohort_script.R`, or  
- Placed the required `cohortMethodData_*.zip` files in `results/`.

**Bootstrap (after the main workflow):**

```r
source("bootstrap_m_of_n.R")
```

Then run the point-estimate + CI and plotting scripts as needed:

```r
source("hr_point_est_with_bootstrap_ci.R")
source("rmst_point_est_with_bootstrap_ci.R")
source("plot_hr_bootstrap_results.R")
source("plot_rmst_bootstrap_results.R")
source("plot_bootstrap_survival_curves.R")
```

## Outputs

- **results/**  
  CohortMethod data zips, `ps_study.rds`, `Cox_censoring.rds`, survival/weights CSV, Table 1 and covariate CSVs.

- **bootstrap_results/**  
  Bootstrap HR and RMST summaries, survival curve outputs, and any other files written by `bootstrap_m_of_n.R` and the plotting scripts.

## License

See repository and OHDSI package licenses for terms of use.
