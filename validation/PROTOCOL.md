# Manuscript comparison

Supply local reference CSVs exported from the manuscript's archived analysis outputs. They are comparison targets and are never substituted for newly fitted estimates. Manuscript reference files and participant data are not bundled with this code.

Run after all four pipelines finish:

```sh
python3 validation/compare_results.py --reference-dir /path/to/manuscript/reference
Rscript validation/plot_comparison.R validation/output/cell_comparison.csv
```

Use `--run-root /path/to/isolated/run` to compare an external run, or `--families main nhanes3 cause biomarker` to select completed families. Cohort and event counts, analysis keys, and uniqueness must match exactly.

The reference filenames are `main_cqr.csv`, `nhanes3_cqr.csv`, `cause_cqr.csv`, and `biomarker.csv`. Use the same column names as each pipeline's aggregate output. CQR comparisons require scenario, imputation, term, tau, n, deaths, estimate, lower and upper, plus endpoint for cause-specific results. Biomarker comparisons require pair, tau, n, deaths, biomarker_low, biomarker_high, valid_replicates, the P25/P75 effect columns and their intervals, and the interaction estimate and interval. Optional p_value and q_primary columns support significance checks. The comparison reader does not use participant-level files.

The comparison allows numerical variation. It records absolute differences, signs, interval changes, and quantile availability. A difference greater than the larger of 0.1 years and 20% of the reference magnitude is flagged for review. This is a screening rule, not a prespecified statistical equivalence margin or proof of scientific agreement. Near-zero sign changes, confidence intervals crossing zero, and changed availability are reported separately and must be reviewed rather than silently counted as agreement. No estimates are rounded or changed to force a match.

Biomarker comparisons include the P25 effect, P75 effect and their interaction contrast, plus changes in the number of finite interaction replicates. Interval, P-value and BH-value availability changes are separate from significance comparisons among values available in both runs. The report distinguishes missing intervals from null effects and records the largest difference between comparable interval bounds. Manuscript CQR intervals use public-package resampling; biomarker intervals use the separate survey-replicate procedure.

For Figure 2, compare the prediction CSVs separately; spline-basis coefficients alone do not describe the exposure contrast:

```sh
Rscript validation/compare_splines.R /path/to/reference_predictions.csv metal_survival_analysis/results/spline_predictions/figure2_predictions_by_imputation.csv validation/output
```

This checks the complete exposure grid and cohort counts, records prediction and interval differences, counts missing interval bounds separately, and compares P95-minus-P5 exposure contrasts using the same 0.1-year/20% review rule. The R overlay preserves missing curves and shows each imputation separately.
