# Biomarker effect modification with CQR

This code accompanies the completed November 3, 2026 analysis of blood lead, blood mercury, urinary cadmium and urinary arsenic with serum vitamin C, serum alpha-carotene and RBC folate. It fits censored quantile regression (CQR), with conditional survival time in years.

The completed run included 12 original fits and 4,800 fixed survey replicate fits. Its independent audit passed 54,439 checks and reproduced the previous fitted coefficient matrices exactly. **Seven of 12 primary interaction intervals remain provisional, five are unavailable, and no available primary test meets BH q < 0.05.** Exact rerun agreement does not resolve the documented sensitivity to the numerical solver.

## Scripts

| Script | Purpose |
| --- | --- |
| [01_prepare_biomarkers.R](R/01_prepare_biomarkers.R) | Read laboratory records, calibrate assays, select complete cases and native survey weights |
| [02_prepare_cqr_models.R](R/02_prepare_cqr_models.R) | Prepare continuous interactions, fixed quartiles, model matrices and survey replicate weights |
| [03_fit_cqr.R](R/03_fit_cqr.R) | Fit original models and fixed replicate blocks using the checked Peng–Huang implementation |
| [04_summarize_emm.R](full_run_2026-11-03/R/04_summarize_emm.R) | Summarize the fresh fits with survey covariance, contrasts, intervals and multiplicity adjustment |
| [05_plot_emm.R](full_run_2026-11-03/R/05_plot_emm.R) | Make the R figure, retaining provisional and unavailable interval labels |
| [run_models.py](full_run_2026-11-03/run_models.py) | Schedule bounded, restartable R jobs |

Stages 01–03 come from the revised preparation/fitting directory. Stages 04–05 and the runner come from its completed `full_run_2026-11-03` directory. Earlier reporting copies read an older fit directory and are deliberately omitted.

The supporting R files and supplied native source are included under `replication_through_2014/`. The required aggregate numerical comparison is included at `nutrient_biomarker_emm_run_2026-11-02/solver_diagnostic_v2/aggregate/path_equivalence.csv`.

## Input contract

This is a code publication. Participant data, prepared matrices, replicate-weight banks, fitted objects and raw NHANES files are not included. Run the commands from the repository root after supplying these upstream study artifacts in their original layout:

- `replication_through_2014/manuscript_completion/data/analytic_corrected.rds`: the observed corrected analysis frame, including the existing exposure transforms, outcomes and covariates.
- `replication_through_2014/remote_execution/models/{primary_Pb,primary_Hg,egfr60_Cd,egfr60_As}/{specification.rds,replicate_design.rds}`: parent eligibility/covariates, row identifiers and the fixed full-component survey draws.
- `replication_10_24/revised/raw/`: native NHANES laboratory and demographic files used for period-specific weights.
- `nutrient_biomarker_feasibility_2026-11-02/raw/`: native biomarker files listed in stage 01.
- `nutrient_biomarker_emm_run_2026-11-02/raw/`: CBC files `LAB25.xpt`, `L25_B.xpt`, `L25_C.xpt` and `CBC_D.xpt` for the early folate calibration.

Stage 01 defines all required laboratory filenames and fields. Stages 01–02 generate `data/*.rds`, `models/*/prepared.rds` and `aggregate/model_preparation.csv`; they stop if prepared outputs already exist. The public main-analysis scripts do not automatically generate all of these historical parent artifacts.

Checkpoint guards validate paths, source hashes, runtime and original-fit links. The publication copy retains those guards. Personal library paths have been replaced with an optional environment variable, so this copy must not be presented as a byte-identical historical checkpoint environment. No models were refitted for publication.

## Run

The audited numerical environment was Linux with single-threaded OpenBLAS, R 4.6.1, quantreg 6.1, survival 3.8-6, survey 4.5 and tidyverse 2.0.0. Preparation also uses haven. The plotting script requires ggplot2 3.5.2, the version used successfully for the audited figure. R packages are not installed by these scripts.

`HM_REMOTE_R_LIBRARY` optionally selects the project R package library. `HM_R_LIBRARY` optionally supplies an additional library and selects the plotting library; omit it when the required packages are already on the default search path.

Compile the supplied audited Fortran source on the execution machine. A Fortran compiler and BLAS/LAPACK development libraries are required. The native code and compilation notice are described [here](../../replication_through_2014/remote_execution/validation/instrumented_ph/README.md).

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

# Run this command inside the native source directory.
cd replication_through_2014/remote_execution/validation/instrumented_ph/src
R CMD SHLIB ph_diagnostic.f -o ph_diagnostic.so -llapack -lblas
cd ../../../../..

# Run the remaining commands from the repository root.
Rscript nutrient_biomarker_emm_run_2026-11-02/project_style_v3/R/01_prepare_biomarkers.R
Rscript nutrient_biomarker_emm_run_2026-11-02/project_style_v3/R/02_prepare_cqr_models.R
python3 nutrient_biomarker_emm_run_2026-11-02/project_style_v3/full_run_2026-11-03/run_models.py original 12
python3 nutrient_biomarker_emm_run_2026-11-02/project_style_v3/full_run_2026-11-03/run_models.py replicates 24
Rscript nutrient_biomarker_emm_run_2026-11-02/project_style_v3/full_run_2026-11-03/R/04_summarize_emm.R
Rscript nutrient_biomarker_emm_run_2026-11-02/project_style_v3/full_run_2026-11-03/R/05_plot_emm.R
```

All original jobs must finish before the replicate stage starts. The runner permits at most 24 single-threaded R workers. A direct worker call takes a pair and inclusive replicate limits, for example `Rscript nutrient_biomarker_emm_run_2026-11-02/project_style_v3/R/03_fit_cqr.R Pb_vitamin_c 0 0`. Original replicate 0 must be requested separately from replicates 1–400.

## Model and interpretation

The model is `Surv(followup_yrs, event) ~ metal_log2 * biomarker_iqr + parent covariates`, with the existing age-squared term. The continuous biomarker is scaled as `(value − weighted P25)/(weighted P75 − weighted P25)`. The interaction therefore estimates the difference in the metal-doubling time association at biomarker P75 versus P25. These are evaluation levels of a continuous variable, not separately fitted participant subgroups.

The primary survival quantile is 0.10. Quantiles 0.25 and 0.50 are descriptive secondary analyses. Each of the 12 pairs uses its own observed complete-case sample; there is no imputation. Vitamin C uses 2003–2006 measurements. Alpha-carotene uses 2001–2006 for Pb/Cd and 2003–2006 for Hg/As. RBC folate uses 1999–2014 for Pb/Cd and 2003–2014 for Hg/As. Urinary metal models retain observed baseline eGFR ≥60.

The EMM implementation uses a checked native Peng–Huang solver, its existing explicit grid and fixed numerical retry rules. It is distinct from the public `quantreg::crq()` package-summary calculation used by the main analysis. It does not use Cox regression or the reference project's Cox mixture interaction estimator.

Inference uses 400 inherited survey draws, positive affine rescaling of Rao–Wu multipliers with kappa=2, and a fourfold variance-scale correction. `survey::svrVar()` uses replicate-mean centering and stored domain degrees of freedom. Each contrast requires a finite original estimate and all 400 finite replicate contrasts. The rescaling check preserves linear-total variance; it does not establish exact variance preservation for nonlinear CQR. BH adjustment retains the full planned 12-test primary family, including unavailable tests.

The prior alternative-solver comparison passed only 16 of 24 control comparisons. Solver sensitivity also affected controls for As–vitamin C and As–RBC folate, despite their complete intervals. All available inference remains provisional and exploratory. Missing intervals are not null findings, and a positive interaction alone does not establish nutrient or supplement protection. Complete-case selection, assay conversions, residual confounding and censoring assumptions remain additional limitations.

## Sources

Workflow organization and compatible R package use were informed by [Plasticizer-mixture-and-mortality](https://github.com/YZhang9719/Plasticizer-mixture-and-mortality/tree/da34a8bdeff9d665da38f7f07c853aaaf71bc5c3), pinned to commit `da34a8bdeff9d665da38f7f07c853aaaf71bc5c3` (MIT). Its Cox mixture estimand and imputation strategy were not adopted.

The native implementation derives from `quantreg` 6.1 and retains its GPL (version 2 or later) provenance. See the [native source notice](../../replication_through_2014/remote_execution/validation/instrumented_ph/README.md) and repository license notices. This is separate from the reference workflow's MIT license.
