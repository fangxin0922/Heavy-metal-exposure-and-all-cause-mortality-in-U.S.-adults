# Heavy-metal exposure and all-cause mortality in U.S. adults

R code for NHANES analyses of heavy-metal exposure and conditional survival time.

| Analysis | Contents |
|---|---|
| [Main analysis](metal_survival_analysis/README.md) | Cohort preparation, imputation, censored quantile regression (CQR), survey Cox models, and figures. |
| [Biomarker EMM](biomarker_emm/README.md) | Nutrient-biomarker effect-measure modification and survey replicate inference. |

## Getting started

Requires R, Python 3, and the study input files. Participant data and fitted objects are not included; raw NHANES downloads alone are insufficient. Biomarker EMM also requires fixed parent survey-draw banks.

Follow the linked guides for dependencies, inputs, and commands. For the main analysis, run from the repository root with empty output folders:

```sh
python3 metal_survival_analysis/run_analysis.py --workers 16 --prepare-workers 4
```

The default runs all 23 scenarios; add `--scope primary` for the four primary cohorts.

## Interpretation

- CQR estimates conditional survival-time differences in years, not mean years of life lost or causal effects.
- Main CQR intervals are imputation-specific, without Rubin pooling or NHANES design adjustment.
- Biomarker EMM intervals remain provisional because numerical solver checks were not fully satisfied.

See the [main analysis guide](metal_survival_analysis/README.md#urinary-sensitivity-cohort-correction) for renal eligibility restrictions and historical unrestricted comparisons.

## Attribution

Adapted components draw on [Plasticizer-mixture-and-mortality](https://github.com/YZhang9719/Plasticizer-mixture-and-mortality) and [quantreg](https://github.com/cran/quantreg). Their [MIT](licenses/reference-MIT.txt) and [GPL-2-or-later](licenses/GPL-2.txt) notices are retained.
