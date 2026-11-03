# Heavy-metal exposure and all-cause mortality in U.S. adults

Analysis code for NHANES heavy-metal exposure and conditional survival time, with a separate nutrient-biomarker effect-measure modification (EMM) analysis. This source release contains the revised main pipeline and the scripts used for the latest completed EMM run on November 3, 2026.

## Analysis scripts

| Analysis | Start here | Main method |
|---|---|---|
| Metal exposure and survival | [Pipeline guide and script index](metal_survival_analysis/GUIDE.md) | `quantreg::crq()` censored quantile regression (CQR), five MICE imputations, and package resampling intervals. |
| Nutrient-biomarker EMM | [EMM guide and script index](nutrient_biomarker_emm_run_2026-11-02/project_style_v3/README.md) | Observed complete cases, continuous metal–biomarker interactions, checked native Peng–Huang CQR, and 400 survey replicate fits per model. |

The main pipeline also includes survey-weighted descriptive statistics and supplemental Cox models. All figures are generated in R. Python only schedules the R jobs.

## Running the code

Clone the repository and run the commands in each analysis guide from the repository root. First review the [software environment](ENVIRONMENT.md) and the input requirements in those guides.

**Archived study inputs are required.** This is a code release, not a standalone NHANES downloader. It does not include participant-level records, imputed datasets, fitted model objects, survey replicate banks, compiled libraries, or historical comparator datasets. The two small CSV inputs included in the repository contain survey-weight specifications and aggregate solver diagnostics.

The publication copy makes library locations configurable while preserving the analysis calculations and numerical checks. It has been checked for syntax, dependency paths, and agreement with the archived source. The publication copy has not been used for another numerical run. Historical checkpoints cannot be treated as new completed fits after changing source paths or runtime libraries.

## Interpretation

CQR coefficients describe differences in conditional survival-time quantiles, in years. They are not estimates of mean years of life lost or causal effects.

The main CQR intervals are specific to each imputation; they are neither Rubin-pooled nor adjusted for the NHANES survey design. The Peng–Huang package-summary resampling implementation used by this run omits the original sampling weights from resampled fits. The separate EMM workflow uses survey replication and its documented numerical acceptance rules. Its available intervals remain provisional because solver-sensitivity checks were not fully satisfied; an unavailable interval is not a null finding.

## Sources

The main workflow draws on [YZhang9719/Plasticizer-mixture-and-mortality](https://github.com/YZhang9719/Plasticizer-mixture-and-mortality/tree/da34a8bdeff9d665da38f7f07c853aaaf71bc5c3). See [method attribution](metal_survival_analysis/SOURCES.md) and [third-party notices](THIRD_PARTY_NOTICES.md). The EMM native solver derives from quantreg 6.1 and carries separate GPL attribution.

The [source manifest](SOURCE_MANIFEST.csv) records the SHA-256 checksum of each published file other than the manifest itself. [Publication provenance](provenance/README.md) identifies the source editions and the small library-path adaptations.
