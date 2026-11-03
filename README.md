# Heavy-metal exposure and all-cause mortality in U.S. adults

Analysis code for NHANES heavy-metal exposure and conditional survival time, with a separate nutrient-biomarker effect-measure modification (EMM) analysis. 

## Analysis scripts

| Analysis | Start here | Main method |
|---|---|---|
| Metal exposure and survival | [Pipeline guide and script index](metal_survival_analysis/GUIDE.md) | `quantreg::crq()` censored quantile regression (CQR), five MICE imputations, and package resampling intervals. |
| Nutrient-biomarker EMM | [EMM guide and script index](nutrient_biomarker_emm_run_2026-11-02/project_style_v3/README.md) | Observed complete cases, continuous metal–biomarker interactions, checked native Peng–Huang CQR, and 400 survey replicate fits per model. |

The main pipeline also includes survey-weighted descriptive statistics and supplemental Cox models. All figures are generated in R. 


