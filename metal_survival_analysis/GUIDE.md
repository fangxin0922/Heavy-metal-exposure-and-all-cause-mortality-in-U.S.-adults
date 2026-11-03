# Metal survival analysis

This is the November 3, 2026 code edition of the heavy-metal censored quantile regression (CQR) pipeline. It includes input reconstruction, multiple imputation, primary and sensitivity models, comparisons with archived fits, supplemental Cox models, and R figures.

The scripts require the existing study archives and NHANES files listed in [INPUTS.md](INPUTS.md). This repository is a code release, not a standalone raw-data download or a complete data archive. No participant-level data or saved model objects are included. The publication copy has been checked for syntax and preserved statistical expressions; it has not been refitted after relocation.

## Scripts

| File | Purpose |
|---|---|
| [run_analysis.py](run_analysis.py) | Schedule R jobs, record status, and validate resumable outputs. |
| [settings.R](R/settings.R) | Set paths, packages, seeds, and analysis constants. |
| [data_tools.R](R/data_tools.R) | Read source files and validate participant keys and values. |
| [exposure_records.R](R/exposure_records.R) | Rebuild metal exposures, weights, and urinary standardization. |
| [survival_records.R](R/survival_records.R) | Read public mortality linkage and follow-up time. |
| [adjustment_data.R](R/adjustment_data.R) | Assemble covariates and audited baseline corrections. |
| [cohort_builder.R](R/cohort_builder.R) | Join records, define cohorts, and compare inputs with the archive. |
| [model_catalog.R](R/model_catalog.R) | Define the 23 scenarios and adjustment formulas. |
| [impute_covariates.R](R/impute_covariates.R) | Create five MICE datasets for one scenario. |
| [fit_cqr.R](R/fit_cqr.R) | Fit one CQR model and obtain package confidence intervals. |
| [compare_runs.R](R/compare_runs.R) | Compare current estimates with saved historical results. |
| [survey_cox.R](R/survey_cox.R) | Produce weighted descriptions and supplemental Cox models. |
| [draw_figures.R](R/draw_figures.R) | Plot saved model results in R. |
| [replay_baseline.R](R/replay_baseline.R) | Refit 20 archived primary inputs as point-estimate controls. |
| [check_saved_models.R](review/check_saved_models.R) | Validate hashes and keys in completed output objects. |
| [plot_run_comparison.R](review/plot_run_comparison.R) | Draw the primary comparison with labeled availability counts. |

## Environment and execution

Use a Unix-like environment with Python 3 and `Rscript` on `PATH`; the runner uses `fcntl` and the control replay uses forked R workers. R requires `dplyr`, `tidyr`, `haven`, `survival`, `survey`, `mice`, and `quantreg`, plus their dependencies. The completed source run used R 4.6.1 on Linux, with `quantreg` 6.1, `mice` 3.19.0, `survey` 4.5, `survival` 3.8-6, `haven` 2.5.5, `tidyr` 1.3.2, and `dplyr` 1.2.1. These describe that run; they are not an installation lockfile.

All commands below run from the repository root. This command displays options without reading data or fitting models:

```sh
python3 metal_survival_analysis/run_analysis.py --help
```

R searches its standard library paths, `metal_survival_analysis/library/`, and the legacy `replication_10_24/library/`. Set `HM_R_LIBRARY` to one additional package-library directory if needed. After supplying every required input in [INPUTS.md](INPUTS.md), run:

```sh
python3 metal_survival_analysis/run_analysis.py --workers 16 --prepare-workers 4
Rscript --vanilla metal_survival_analysis/review/check_saved_models.R
Rscript --vanilla metal_survival_analysis/review/plot_run_comparison.R
```

Choose worker counts suitable for available memory and CPUs. `--scope primary` restricts fresh fitting to `primary_Pb`, `primary_Hg`, `egfr60_Cd`, and `egfr60_As`; the driver still runs the historical control and comparison stages. `--points-only` is a diagnostic option and does not calculate intervals. The full artifact-check script expects the complete 23-scenario run.

After the four input stages have completed, one scenario can be prepared and fitted directly:

```sh
Rscript --vanilla metal_survival_analysis/R/impute_covariates.R egfr60_Cd
Rscript --vanilla metal_survival_analysis/R/fit_cqr.R egfr60_Cd 1 400 fresh
```

All calculations and figures use R; Python schedules jobs. Source and input hashes guard reuse of saved outputs. Start with an empty output state for a new run rather than copying earlier fits into this folder.

## Models and interpretation

The primary analysis calls `quantreg::crq(method = "PengHuang")`; Portnoy models are a sensitivity analysis. It retains the package fitting grid, quantiles 0.10/0.25/0.50, original covariate units, and weights divided by their cohort mean. MICE uses five imputations, 20 iterations, and seed 4242. Estimable models request 400 package jackknife resamples through `summary(..., se = "boot", bmethod = "jack")`. Unavailable quantiles remain unavailable.

The four primary cohorts are blood Pb and Hg, plus urinary Cd and As restricted to observed baseline eGFR ≥60 before imputation. Other scenarios preserve the historical no-BMI, age ≥55, Portnoy, joint-metal, EPA+DHA, and spline comparisons. The complete specification is in [model_catalog.R](R/model_catalog.R).

CQR coefficients describe differences in conditional survival-time quantiles, in years, per exposure doubling. They are not mean years of life lost or causal effects. The CQR confidence intervals are separate for each imputation, without Rubin pooling or NHANES stratified-PSU variance estimation. In the source run's `quantreg` implementation, the Peng–Huang summary bootstrap omits original sampling weights from resampled fits. Survey-weighted Cox models and unweighted Cox spline checks are separate supplemental analyses.

The original population-screening builder and some questionnaire/dietary builders are unavailable. Those fields retain the archived values and missingness described in [input_lineage.md](review/input_lineage.md). Blood logarithms are independently checked before their archived numeric representation is retained; see [numeric_precision.md](review/numeric_precision.md). Raw NHANES downloads alone cannot reproduce these exact archived inputs.

## Outputs and provenance

Execution creates `RUN_MANIFEST.json`, `RUN_STATUS.json`, `RUN_REPORT.md`, and the `data/`, `models/`, `fits/`, `logs/`, `audit/`, `results/`, and `figures/` directories beneath `metal_survival_analysis/`. These are generated outputs, not bundled results.

This package derives from the local `metal_survival_analysis/` presentation edition, itself derived from the completed `fresh_reference_pipeline_20261103/` run. The publication adaptation only replaces one personal package-library path with `HM_R_LIBRARY` and updates documentation. See [REVISION_NOTES.md](REVISION_NOTES.md) and [SOURCES.md](SOURCES.md) for the exact provenance and methodological attribution.
