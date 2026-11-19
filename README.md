# Heavy-metal exposure and all-cause mortality in U.S. adults

R code for NHANES analyses of heavy-metal exposure and conditional survival time.

| Analysis | Contents |
|---|---|
| [Main analysis](metal_survival_analysis/README.md) | Cohort preparation, imputation, censored quantile regression (CQR), survey Cox models, and figures. |
| [Biomarker EMM](biomarker_emm/README.md) | Nutrient-biomarker effect-measure modification and survey replicate inference. |

## Getting started

Requires R, Python 3, and the retained study reference archive. [Input preparation](scripts/README.md) stages that archive, retrieves public files, and generates biomarker parent inputs. Participant data and fitted objects are not included. The original reference archive still cannot be reconstructed from raw NHANES alone.

Follow the linked guides for dependencies, inputs, and commands. For the main analysis, run from the repository root with empty output folders:

```sh
python3 metal_survival_analysis/run_analysis.py --workers 16 --prepare-workers 4
```

The default runs all 23 scenarios; add `--scope primary` for the four primary cohorts.

## Attribution

Adapted components draw on [Plasticizer-mixture-and-mortality](https://github.com/YZhang9719/Plasticizer-mixture-and-mortality) and [quantreg](https://github.com/cran/quantreg). Their [MIT](licenses/reference-MIT.txt) and [GPL-2-or-later](licenses/GPL-2.txt) notices are retained.
