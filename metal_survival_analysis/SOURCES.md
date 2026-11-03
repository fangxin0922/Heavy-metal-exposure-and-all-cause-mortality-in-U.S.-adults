# Sources and attribution

This publication copy derives from the project's `metal_survival_analysis/` presentation edition dated November 3, 2026. That edition reorganized the completed `fresh_reference_pipeline_20261103/` heavy-metal CQR pipeline. These names identify the local source editions; their data and results are not bundled in this repository.

The external workflow reference is YZhang9719's [Plasticizer-mixture-and-mortality](https://github.com/YZhang9719/Plasticizer-mixture-and-mortality/tree/da34a8bdeff9d665da38f7f07c853aaaf71bc5c3), pinned to commit `da34a8bdeff9d665da38f7f07c853aaaf71bc5c3`. Its MIT copyright and permission notice is retained in [licenses/reference-MIT.txt](licenses/reference-MIT.txt).

The project adopts the reference's explicit data-building stages, `haven` reads and keyed joins, `lm`/`predict` urinary standardization, and supplemental `survey::svydesign`, `survey::svycoxph`, and `survival::coxph(pspline())` models. The main analysis retains the project's `quantreg::crq` CQR models and MICE imputation. The reference's `missForest` imputation and Cox mixture models are not substituted for those study methods.

The [function map](review/method_sources.md) distinguishes shared package functions from the study-specific model choices. The [input audit](review/input_lineage.md) and [blood-log audit](review/numeric_precision.md) record limitations in the inherited source data. These limitations also remain visible in the executable input stages.

The publication adaptation replaces the personal `/home/<user>/R/library` fallback with the optional `HM_R_LIBRARY` environment variable. It does not change formulas, cohorts, seeds, model calls, interval calculations, or output paths. No numerical analysis was repeated for this code publication.
