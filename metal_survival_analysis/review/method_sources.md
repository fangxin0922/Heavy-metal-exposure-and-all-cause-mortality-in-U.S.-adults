# Function and method map

The workflow reference is [YZhang9719/Plasticizer-mixture-and-mortality at da34a8b](https://github.com/YZhang9719/Plasticizer-mixture-and-mortality/tree/da34a8bdeff9d665da38f7f07c853aaaf71bc5c3). Attribution concerns the workflow and the functions actually reused; the two studies do not use the same primary estimator.

| Operation | This pipeline | Relationship to the reference |
|---|---|---|
| Native laboratory and demographic input | `haven::read_xpt`, keyed `dplyr` joins, explicit variable mappings | Uses the same package approach, with metal-specific files and weights. |
| Mortality input | Fixed-width parsing of public 2019 linkage | Same source format; this pipeline parses the published positions with base R. |
| Urinary standardization | `stats::lm` and `stats::predict` | Shared covariate-adjusted standardization functions; the study's original covariates and reference population are preserved. |
| Missing covariates | `mice::mice` and `mice::complete` | Study-specific five-imputation method. The reference instead uses `missForest`. |
| Primary survival model | `quantreg::crq(method = "PengHuang")` | Study-specific censored quantile regression, absent from the reference. |
| Alternative CQR fit | `quantreg::crq(method = "Portnoy")` | Study-specific sensitivity model. |
| CQR uncertainty | `summary` with `se = "boot"`, `bmethod = "jack"`, `R = 400` | Package intervals per imputation; not reference Cox inference and not MI-pooled survey intervals. |
| Survey descriptions and Cox model | `survey::svydesign`, `survey::svycoxph` | Shared direct package calls with the study's strata, PSUs, weights, and covariates. |
| Cox shape checks | `survival::coxph` with `pspline` | Shared penalized-spline approach; supplemental unweighted checks. |
| Plotting | R graphics | All figures are created in R from saved results. |

Exposure definitions, cohort restrictions, imputation strategy, and estimands remain study-specific. The reference's phthalate/BPA mixture, `qgcomp` models, attributable fractions, and nutrient interactions are not relabeled as heavy-metal CQR results. Historical panel intervals produced by a different solver/bootstrap pipeline are compared only at the point-estimate level.
