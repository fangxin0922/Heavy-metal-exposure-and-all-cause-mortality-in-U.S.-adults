# Script presentation revision

The November 3 completed run is preserved. This separate code edition changes script names, internal path references, comments and report text. It does not change analysis formulas, exposure definitions, cohort selection, imputation settings, model calls, seeds or interval calculations.

Comments now state each operation's purpose or a necessary numerical/data caveat. Reference attribution is collected in [SOURCES.md](SOURCES.md), with its license notice retained.

## Filename map

| Previous file | Revised file |
|---|---|
| `run_pipeline.py` | `run_analysis.py` |
| `00_setup.R` | `settings.R` |
| `01_04_helpers.R` | `data_tools.R` |
| `01_dataset_creation_chemicals.R` | `exposure_records.R` |
| `02_dataset_creation_mortality.R` | `survival_records.R` |
| `03_dataset_creation_covariates.R` | `adjustment_data.R` |
| `04_data_cleaning.R` | `cohort_builder.R` |
| `05_model_specifications.R` | `model_catalog.R` |
| `06_covariate_imputation.R` | `impute_covariates.R` |
| `07_censored_quantile_analysis.R` | `fit_cqr.R` |
| `08_compare_results.R` | `compare_runs.R` |
| `09_descriptive_and_cox.R` | `survey_cox.R` |
| `10_figures.R` | `draw_figures.R` |
| `11_replay_control.R` | `replay_baseline.R` |
| `audit/validate_completed_artifacts.R` | `review/check_saved_models.R` |
| `audit/plot_primary_comparison.R` | `review/plot_run_comparison.R` |
| Generated `RESULTS.md` | Generated `RUN_REPORT.md` |

The report wording now explicitly includes jointly unavailable cells in agreement counts and distinguishes resolving a summary stage from executing resamples. These are presentation changes; the underlying counts are unchanged.

## Publication adaptation and checks

The publication copy derives from the November 3, 2026 `metal_survival_analysis/` edition. One executable expression changes: the personal package-library fallback in `R/settings.R` becomes `Sys.getenv("HM_R_LIBRARY", unset = "")`. The standard and project-relative library search paths remain. Documentation now describes the external input requirements and uses links within this publication package.

All 15 R scripts were parsed without evaluating the analysis. Their parsed expressions match the source edition after normalizing that one library-path expression. The Python runner is byte-identical to the source edition, passes syntax compilation, and displays `--help` without starting a job. All local R helper references and runner script names resolve. The metadata CSV is copied exactly from its source.

No imputation or model fitting was performed for publication. No source-edition scripts were overwritten. Code-equivalence checks establish the scope of these edits; they do not guarantee identical numerical results under a different R, package, BLAS, or operating-system environment.
