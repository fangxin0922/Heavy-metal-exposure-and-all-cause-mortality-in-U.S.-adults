# Heavy-metal exposure and all-cause mortality in U.S. adults

Analysis code for NHANES heavy-metal exposure and conditional survival time, with a separate nutrient-biomarker effect-measure modification (EMM) analysis.

| Folder | Analysis |
|---|---|
| [metal_survival_analysis](metal_survival_analysis) | Input preparation, MICE imputation, CQR models, supplemental survey Cox models, result tables and figures. |
| [biomarker_emm](biomarker_emm) | Biomarker preparation, continuous interactions, checked Peng–Huang CQR, survey replicate inference and figures. |

All analyses and figures use R. Python schedules the model fits. Participant data and fitted objects are not included.

## Main analysis

Run from the repository root with Python 3 and R packages `dplyr`, `tidyr`, `haven`, `survival`, `survey`, `mice` and `quantreg`. Set `HM_R_LIBRARY` if using a separate R package library. Use a Unix-like environment and empty output folders.

```sh
python3 metal_survival_analysis/run_analysis.py --workers 16 --prepare-workers 4
```

The runner executes the input scripts, `impute_covariates.R`, `fit_cqr.R`, `summarize_results.R`, `survey_cox.R` and `draw_figures.R`. `--scope primary` selects the four primary cohorts; the default runs all 23 scenarios.

Supply these study inputs under `data/metal_survival/`, or set `HM_DATA_DIR` to their parent directory:

- `analytic.rds`: the complete inherited 1999–2018 study frame, before the baseline corrections, with its original column types and factor levels.
- `raw/`: demographic, metal-laboratory and public linked-mortality files.
- `creatinine/`: `L16_B.xpt`, `ALB_CR_I.xpt` and `ALB_CR_J.xpt`.
- `baseline/`: `L40_B.xpt`, `BMX_B.xpt` and `L13_B.xpt`.
- `diet/`: `DRXTOT_B.xpt` and `DR1TOT_C.xpt` through `DR1TOT_H.xpt`.

Exact filenames and fields are specified in the input scripts. Some questionnaire and dietary fields are inherited from the study frame, so raw NHANES downloads alone are insufficient. The included weight-specification CSV is metadata, not participant data.

Main CQR uses five imputations and 400 package jackknife resamples where estimable. Intervals are reported separately by imputation; they are not Rubin-pooled or NHANES-design-adjusted. The Peng–Huang package resampling implementation used here omits the original sampling weights from resampled fits.

## Biomarker EMM

Use Linux with single-threaded OpenBLAS, quantreg 6.1, and the packages above plus `tidyverse`. The separate figure script requires ggplot2 3.5.2. A Fortran compiler and BLAS/LAPACK development libraries are needed for the supplied native solver.

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

Outputs go to `biomarker_emm/output/`, configurable through `EMM_OUTPUT_DIR`. Complete the original fits before the replicate stage. The main pipeline does not create the required parent survey-draw banks.

EMM uses observed complete cases, continuous biomarker interactions and survey replicate variance. Available intervals remain provisional because numerical solver checks were not fully satisfied. An unavailable interval is not a null finding. CQR effects are conditional survival-time differences in years, not mean years of life lost or causal effects.

## Attribution

Data-preparation and supplemental-model conventions draw on [Plasticizer-mixture-and-mortality](https://github.com/YZhang9719/Plasticizer-mixture-and-mortality); its [MIT notice](licenses/reference-MIT.txt) is retained. The EMM native solver and R interface derive from [quantreg](https://github.com/cran/quantreg), by Roger Koenker and contributors, including Stephen Portnoy for censored QR. Study modifications add 33 native diagnostic fields and path-validity checks. The quantreg-derived code retains [GPL version 2 or later](licenses/GPL-2.txt).
