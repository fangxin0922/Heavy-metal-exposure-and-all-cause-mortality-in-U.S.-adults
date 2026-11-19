# NHANES III replication

Reproduces the 1988–1994 replication with mortality follow-up through 2019: blood Pb (9,388 participants; 6,146 deaths) and urinary Cd (7,924; 4,823). Only Cd requires observed eGFR ≥60. Alcohol categories use the final harmonised thresholds.

From the repository root, with the R environment described in [the main analysis guide](../metal_survival_analysis/README.md):

```sh
export HM_R_LIBRARY=/path/to/R/library
python3 nhanes3_analysis/run_analysis.py --workers 8
```

Use Python 3.9 or later. `--source-dir /path/to/original/raw` reuses the original CDC files. Public files are verified against `sources.json`; selected fields are read from the official SAS column definitions. No historical study RDS or previous imputations are needed for this pipeline.

Preparation fits the urinary-creatinine standardisation model, applies cohort eligibility, and runs five MICE imputations with 50 iterations. Imputation seeds are 2026102904 for Pb and 2026102905 for Cd. Sampling weights are extracted from the full MEC survey design before cohort selection, retaining the original package representation; direct raw weights can differ at machine precision and alter a sensitive CQR path. Public `quantreg::crq(method="PengHuang")` uses these weights after mean normalisation; its `summary()` requests 400 delete-d jackknife resamples. Results remain separate by imputation, with unavailable quantiles retained. These intervals are neither MI-pooled nor NHANES-design-adjusted.

Outputs: `results/results_by_imputation.csv`, R-generated forest plots under `figures/`, and preparation/fit diagnostics. `--points-only` omits intervals. `--fit-only` uses already completed model preparations; fit output directories must still be empty. Set `NH3_OUTPUT_DIR` for a separate output location.
