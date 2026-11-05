# Heavy-metal exposure and all-cause mortality in U.S. adults

Analysis code for NHANES heavy-metal exposure and conditional survival time, with a separate nutrient-biomarker effect-measure modification (EMM) analysis.

| Folder | Analysis |
|---|---|
| [metal_survival_analysis](metal_survival_analysis) | Input preparation, MICE imputation, CQR models, supplemental survey Cox models, result tables and figures. |
| [biomarker_emm](biomarker_emm) | Biomarker preparation, continuous interactions, checked Peng–Huang CQR, survey replicate inference and figures. |

All analyses and figures use R. Python schedules the model fits. Participant data and fitted objects are not included. Reproduction requires the study input frames and fixed survey-draw banks described below; this repository alone cannot reconstruct all required inputs from raw NHANES downloads.

## Main analysis

Run from the repository root with Python 3 and R. The analysis environment uses R 4.6.1, `dplyr` 1.2.1, `tidyr` 1.3.2, `haven` 2.5.5, `survival` 3.8-6, `survey` 4.5, `mice` 3.19.0 and `quantreg` 6.1. Set `HM_R_LIBRARY` if using a separate R package library. Use a Unix-like environment and empty output folders.

```sh
python3 metal_survival_analysis/run_analysis.py --workers 16 --prepare-workers 4
```

The runner executes the input scripts, `impute_covariates.R`, `fit_cqr.R`, `summarize_results.R`, `survey_cox.R` and `draw_figures.R`. `--scope primary` selects the four primary cohorts; the default runs all 23 scenarios.

### Urinary sensitivity cohort correction

The Cd/As Portnoy, age ≥55 and BMI-omission sensitivities, and both joint Pb/Cd specifications, require **finite observed baseline eGFR ≥60 mL/min/1.73 m² before imputation**. Missing observed eGFR and values below 60 are excluded. These eight corrected scenarios use distinct names (`portnoy_egfr60_Cd/As`, `age55plus_egfr60_Cd/As`, `no_bmi_egfr60_Cd/As` and `joint_egfr60_Pb/Cd`) so earlier unrestricted outputs are not reused or relabelled. All runs require empty output folders.

The primary Cd/As cohorts are `egfr60_Cd` and `egfr60_As`. Their preparations finish before dependent sensitivities begin. The restricted Portnoy models reuse the exact five primary completed datasets, sampling weights and formulas, changing only the estimator. Age ≥55 and BMI-omission sensitivities receive five fresh imputations after eligibility selection. Continuous eGFR adjustment is retained in Cd/As formulas.

Both joint specifications use the same Cd-eligible cohort with finite Pb and Cd log exposures and Cd component sampling weights. They are additive mutual-adjustment models, not Pb×Cd product interactions. The Pb-target formula retains its original covariates; the Cd-target formula additionally adjusts for continuous eGFR. Both apply the renal eligibility restriction.

The explicit unrestricted urinary comparisons remain `primary_Cd` and `primary_As`, despite their historical names; they are not the designated primary analyses. The separate `rcs_Cd` branch in this repository also retains its historical unrestricted cohort and does not reproduce the separately corrected renal-restricted cadmium spline. Blood Pb/Hg sensitivities and biomarker EMM retain their existing definitions.

This update corrects analysis code. Earlier unrestricted sensitivity estimates are not corrected results; new estimates and intervals require a fresh run. Participant data, fitted objects and newly computed estimates are not included in this update.

Supply these study inputs under `data/metal_survival/`, or set `HM_DATA_DIR` to their parent directory:

- `analytic.rds`: the complete study-supplied reference frame covering NHANES cycles 1–10, before the baseline corrections, with the specified column types and factor levels. The outcome models select cycles 1–8; the full reference frame is needed for creatinine standardization.
- `raw/`: demographic, metal-laboratory and public linked-mortality files.
- `creatinine/`: `L16_B.xpt`, `ALB_CR_I.xpt` and `ALB_CR_J.xpt`.
- `baseline/`: `L40_B.xpt`, `BMX_B.xpt` and `L13_B.xpt`.
- `diet/`: `DRXTOT_B.xpt` and `DR1TOT_C.xpt` through `DR1TOT_H.xpt`.

Exact filenames and fields are specified in the input scripts. Questionnaire and dietary fields, participant selection, categorical levels and some numeric representations must be supplied in the study frame; their complete raw-data derivation is not implemented here. The included weight-specification CSV is metadata, not participant data.

Main CQR uses five imputations and 400 package jackknife resamples where estimable. Intervals are reported separately by imputation; they are not Rubin-pooled or NHANES-design-adjusted. The Peng–Huang package resampling implementation used here omits the original sampling weights from resampled fits.

## Biomarker EMM

Use Linux with single-threaded OpenBLAS 0.3.26, LAPACK 3.12.0, and the packages above plus `tidyverse` 2.0.0. The separate figure script requires ggplot2 3.5.2 and uses `readr` 2.2.0. A Fortran compiler and BLAS/LAPACK development libraries are needed for the supplied native solver. Package, compiler and numerical-library differences can affect CQR estimates; matching package names alone does not guarantee identical results.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
(cd biomarker_emm/native && R CMD SHLIB ph_diagnostic.f -o ph_diagnostic.so -llapack -lblas)
Rscript biomarker_emm/R/prepare_biomarkers.R
Rscript biomarker_emm/R/prepare_models.R
python3 biomarker_emm/run_models.py original 12
python3 biomarker_emm/run_models.py replicates 24
Rscript biomarker_emm/R/summarize_results.R
Rscript biomarker_emm/R/plot_interactions.R
```

Supply the following under `data/biomarker_emm/`, or set `EMM_DATA_DIR` to their parent directory:

- `analytic.rds`: the observed corrected study frame.
- `parent_models/{primary_Pb,primary_Hg,egfr60_Cd,egfr60_As}/{specification.rds,replicate_design.rds}`: parent model specifications and the fixed 400 survey draws.
- `raw/labs/`, `raw/biomarkers/` and `raw/cbc/`: the XPT files listed in `prepare_biomarkers.R`.

Outputs go to `biomarker_emm/output/`, configurable through `EMM_OUTPUT_DIR`. The `original` stage fits the full sample; complete it before the replicate stage. The main pipeline does not create the required parent survey-draw banks. These fixed banks and both study input frames are prerequisites for reproducing the reported analysis.

EMM uses observed complete cases, continuous biomarker interactions and survey replicate variance. Available intervals remain provisional because numerical solver checks were not fully satisfied. An unavailable interval is not a null finding. CQR effects are conditional survival-time differences in years, not mean years of life lost or causal effects.

## Attribution

Data-preparation and supplemental-model conventions draw on [Plasticizer-mixture-and-mortality](https://github.com/YZhang9719/Plasticizer-mixture-and-mortality); its [MIT notice](licenses/reference-MIT.txt) is retained. The EMM native solver and R interface derive from [quantreg](https://github.com/cran/quantreg), by Roger Koenker and contributors, including Stephen Portnoy for censored QR. Study modifications add 33 native diagnostic fields and path-validity checks. The quantreg-derived code retains [GPL version 2 or later](licenses/GPL-2.txt).
