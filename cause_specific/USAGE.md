# Cause-specific mortality

Cancer is underlying-cause code 2; cardiovascular mortality combines heart disease (1) and stroke (5). The full primary metal cohorts are retained. Other-cause deaths are censored at their observed death times. The pipeline reuses the primary five completed datasets, formulas and sampling weights, changing only the event indicator.

After preparing the four primary models in `metal_survival_analysis/models/`:

```sh
export HM_R_LIBRARY=/path/to/R/library
python3 cause_specific/build_quantreg.py
python3 cause_specific/run_analysis.py --workers 8
```

Use Python 3.12+, the main analysis's R packages, a Fortran compiler, and BLAS/LAPACK development libraries. The builder downloads the pinned [quantreg source from CRAN](https://cran.r-project.org/package=quantreg), checking both its current and archived locations. `build_quantreg.py --archive /path/to/quantreg_6.1.tar.gz` can use an existing source archive. The archive and native source are checksum-verified.

The isolated quantreg 6.1 build contains the manuscript's memory repair: initialise the `wn` work array over `n1` separately from the `u` loop over `n`. This avoids an out-of-bounds write for sparse outcomes. It does not change the estimator, covariates, weights or quantile grid. The patch is in `native/memory_safety.patch`; the installed package is confined to `cause_specific/library/` and verified before each fit. The main and NHANES III pipelines continue using the unmodified package. The quantreg-derived patch is covered by [GPL-2 or later](../licenses/GPL-2.txt).

Each job requests 400 public-package jackknife resamples, with seed `20261102 + endpoint_index*100000 + scenario_index*1000 + imputation`. Endpoints are ordered cancer, cardiovascular; scenarios are Pb, Hg, Cd, As. Output is imputation-specific, with no MI pooling or design adjustment and no substitution for unavailable quantiles. Interpretation depends on conditional independent censoring of competing events.

Results and cohort counts go to `cause_specific/results/`; R forest plots go to `figures/`. `--points-only` omits intervals; `--fit-only` uses an existing endpoint map. Set `HM_ANALYSIS_DIR` to a different prepared main-analysis directory and `CAUSE_OUTPUT_DIR` for separate outputs. Use empty output directories for a fresh run.
