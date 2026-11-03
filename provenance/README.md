# Publication provenance

This November 3, 2026 source release combines the revised `metal_survival_analysis/` code edition with the latest completed biomarker EMM execution. It does not overwrite the archived study scripts or results.

The main code edition reorganized the completed `fresh_reference_pipeline_20261103/` run with descriptive filenames and concise comments. Its [revision notes](../metal_survival_analysis/REVISION_NOTES.md) map the former names. All 15 published main-analysis R scripts preserve the source edition's parsed expressions after normalization of the package-library path. Its Python runner is byte-identical and its `--help` command was checked without running an analysis.

The EMM release combines preparation/fitting stages 01–03 from `project_style_v3/` with reporting stages 04–05 and the runner from `project_style_v3/full_run_2026-11-03/`. All nine EMM R scripts parse, its Python runner passes syntax checks, and the native Fortran source is byte-identical to the audited 33-field source. Required R helper paths resolve within this release. External data and compiled-library requirements remain documented in the [EMM guide](../nutrient_biomarker_emm_run_2026-11-02/project_style_v3/README.md).

Publication changes to executable files are limited to:

- `metal_survival_analysis/R/settings.R`: replace a personal library directory with optional `HM_R_LIBRARY`.
- `replication_through_2014/remote_execution/R/runtime.R`: replace a personal fallback library directory with optional `HM_R_LIBRARY`, retaining the project-library setting and runtime checks.
- EMM `05_plot_emm.R`: select the plotting library through optional `HM_R_LIBRARY`, retaining the ggplot2 3.5.2 version check.

No model formulas, cohorts, seeds, fit calls, interval calculations, acceptance criteria or hash guards were changed for publication. No imputation or model fitting was performed for this release. Documentation was adapted for the repository, and the EMM runner's output-directory placeholder was retained.

[source-map.csv](source-map.csv) records the relative archived source path and original/published SHA-256 hashes for the 29 executable and required-metadata source files. These archived paths describe the study workspace; they are not links to an additional public dataset. The repository-wide [manifest](../SOURCE_MANIFEST.csv) additionally covers documentation, licenses and recorded sessions.
