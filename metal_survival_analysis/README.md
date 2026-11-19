# Main analysis

Input preparation, MICE imputation, CQR models, supplemental survey Cox models, result tables, and figures. All analyses and figures use R; Python schedules model fits.

Additional manuscript analyses: [NHANES III](../nhanes3_analysis/USAGE.md), [cause-specific mortality](../cause_specific/USAGE.md), and [numerical comparison](../validation/PROTOCOL.md).

## Run

Run from the repository root with Python 3 and R. The analysis environment uses R 4.6.1, `dplyr` 1.2.1, `tidyr` 1.3.2, `haven` 2.5.5, `survival` 3.8-6, `survey` 4.5, `mice` 3.19.0 and `quantreg` 6.1. Set `HM_R_LIBRARY` if using a separate R package library. Use a Unix-like environment and empty output folders.

```sh
python3 metal_survival_analysis/run_analysis.py --workers 16 --prepare-workers 4
```

The runner prepares inputs, imputes covariates, fits the models, aggregates results and creates R figures. `--scope primary` selects the four primary cohorts; the default runs all 23 scenarios.

Figure 1 is rebuilt from native demographic, mortality and metal files, with final cohort membership checked against the prepared analysis. Counts are saved under `results/participant_flow/` and R-generated PDF/PNG diagrams under `figures/`. Regenerate them with `Rscript metal_survival_analysis/R/participant_flow.R` after primary preparation.

Observed baseline summaries for Table 1 are saved under `results/baseline/`: cohort counts, weighted quartiles, unweighted category counts, weighted percentages and missingness. These use the primary preimputation frames and sampling weights. Run `Rscript metal_survival_analysis/R/baseline_characteristics.R` to regenerate them after primary preparation.

`results/supplementary_descriptives/` contains the inherited activity-quartile ranges, observed Hg fish-source counts, and renal-exclusion counts for restricted sensitivities. Regenerate them with `Rscript metal_survival_analysis/R/supplementary_descriptives.R` after primary preparation. Activity ranks are checked on their original reference population before model eligibility restrictions.

The full run also reconstructs Figure 2 from the fitted Pb and corrected Cd splines. `results/spline_predictions/` contains predictions for all five imputations at the observed median/mode covariate profile, over the cohort's P5–P95 exposure range. The two spline PDFs under `figures/` show the curves and imputation-specific pointwise intervals from the saved package covariance. Unavailable curves remain unavailable. To regenerate these outputs from completed fits, run `Rscript metal_survival_analysis/R/spline_predictions.R`.

## Inputs

Use the [input preparation scripts](../scripts/README.md) to stage the reference archive and retrieve the native files. Supply these study inputs under `data/metal_survival/`, or set `HM_DATA_DIR` to their parent directory:

- `analytic.rds`: the complete study-supplied reference frame covering NHANES cycles 1–10, before the baseline corrections, with the specified column types and factor levels. The outcome models select cycles 1–8; the full reference frame is needed for creatinine standardization.
- `raw/`: demographic, metal-laboratory and public linked-mortality files.
- `creatinine/`: `L16_B.xpt`, `ALB_CR_I.xpt` and `ALB_CR_J.xpt`.
- `baseline/`: `L40_B.xpt`, `BMX_B.xpt` and `L13_B.xpt`.
- `diet/`: `DRXTOT_B.xpt` and `DR1TOT_C.xpt` through `DR1TOT_H.xpt`.

Exact filenames and fields are specified in the input scripts. Questionnaire and dietary fields, participant selection, categorical levels and some numeric representations must be supplied in the study frame; their complete raw-data derivation is not implemented here. The included weight-specification CSV is metadata, not participant data.

## Inference

Main CQR uses five imputations and 400 package jackknife resamples where estimable. Intervals are reported separately by imputation; they are not Rubin-pooled or NHANES-design-adjusted. The Peng–Huang package resampling implementation used here omits the original sampling weights from resampled fits.

## Urinary sensitivity cohort correction

The Cd/As Portnoy, age ≥55 and BMI-omission sensitivities, both joint Pb/Cd specifications, and the cadmium spline require **finite observed baseline eGFR ≥60 mL/min/1.73 m² before imputation**. Missing observed eGFR and values below 60 are excluded. These nine corrected scenarios use distinct names (`portnoy_egfr60_Cd/As`, `age55plus_egfr60_Cd/As`, `no_bmi_egfr60_Cd/As`, `joint_egfr60_Pb/Cd`, and `rcs_egfr60_Cd`) so earlier unrestricted outputs are not reused or relabelled. All runs require empty output folders.

The primary Cd/As cohorts are `egfr60_Cd` and `egfr60_As`. Their preparations finish before dependent sensitivities begin. The restricted Portnoy models reuse the exact five primary completed datasets, sampling weights and formulas, changing only the estimator. Age ≥55 and BMI-omission sensitivities receive five fresh imputations after eligibility selection. Continuous eGFR adjustment is retained in Cd/As formulas.

Both joint specifications use the same Cd-eligible cohort with finite Pb and Cd log exposures and Cd component sampling weights. They are additive mutual-adjustment models, not Pb×Cd product interactions. The Pb-target formula retains its original covariates; the Cd-target formula additionally adjusts for continuous eGFR. Both apply the renal eligibility restriction.

The explicit unrestricted urinary comparisons remain `primary_Cd` and `primary_As`, despite their historical names; they are not the designated primary analyses. The cadmium spline is now `rcs_egfr60_Cd`; the old `rcs_Cd` name is rejected. Its ordered participants and sampling weights must match `egfr60_Cd`. Knots are the unweighted 5th, 35th, 65th and 95th exposure percentiles in that restricted cohort. Five fresh imputations retain observed eGFR and spline terms; the corrected amendment summary seeds are 20262105–20262109. Blood Pb/Hg sensitivities and biomarker EMM retain their existing definitions.

This update corrects analysis code. Earlier unrestricted sensitivity estimates are not corrected results; new estimates and intervals require a fresh run. Participant data, fitted objects and newly computed estimates are not included in this update.

## Attribution

Data-preparation and supplemental-model conventions draw on [Plasticizer-mixture-and-mortality](https://github.com/YZhang9719/Plasticizer-mixture-and-mortality); its [MIT notice](../licenses/reference-MIT.txt) is retained.
