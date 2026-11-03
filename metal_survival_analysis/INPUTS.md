# Required study inputs

Paths below are relative to the repository root. Except for the weight-specification CSV, these files are not distributed in this code release. They must be the exact archived objects and native source files expected by the validation assertions; generic replacement datasets are insufficient. No script downloads or synthesizes the missing archives.

## Included metadata

[`replication_10_24/weight_specification.csv`](../replication_10_24/weight_specification.csv) is a cycle-by-metal metadata table containing laboratory filenames, weight variables, and public source URLs. It has no participant records. The table retains its original 1999–2018 metadata; `exposure_records.R` selects cycles 1–8 and recalculates pooling multipliers for the through-2014 analysis. The table's original `K`, `pooling_multiplier`, and status labels must not be interpreted as the final analysis weights.

## Input construction

| Required path | Use |
|---|---|
| `replication_10_24/revised/data/analytic_revised.rds` | Pre-correction 1999–2018 analytic archive. Supplies the retained population, inherited covariates, factor levels, urinary-thallium field, and archived blood-log representations. The full 1999–2018 archive is needed for creatinine standardization, even though outcome models end in 2014. |
| `replication_10_24/revised/raw/DEMO.xpt`, `DEMO_B.xpt` through `DEMO_H.xpt` | Eight demographic files: participant keys, demographics, survey design IDs, and applicable examination weights. |
| `replication_10_24/revised/raw/LAB06.xpt`, `L06_B.xpt`, `L06BMT_C.xpt`, `PBCD_D.xpt` through `PBCD_H.xpt` | Blood lead and mercury laboratory files; early mercury is retained for descriptions only. |
| `replication_10_24/revised/raw/LAB06HM.xpt`, `L06HM_B.xpt`, `L06HM_C.xpt`, `UHM_D.xpt` through `UHM_G.xpt`, `UM_H.xpt` | Urinary cadmium and applicable subsample weights. |
| `replication_10_24/revised/raw/L06UAS_C.xpt`, `UAS_D.xpt` through `UAS_H.xpt` | Arsenic species (`URXUAS3`, `URXUAS5`, `URXUMMA`, `URXUDMA`) and subsample weights. |
| `replication_10_24/revised/raw/NHANES_1999_2000_MORT_2019_PUBLIC.dat` through `NHANES_2013_2014_MORT_2019_PUBLIC.dat` | Eight biennial public linked-mortality files with 2019 follow-up. |
| `replication_through_2014/manuscript_audit/raw_creatinine_check/L16_B.xpt`, `ALB_CR_I.xpt`, `ALB_CR_J.xpt` | Native urinary creatinine used to restore the 2001–2002 baseline value and check the 2015–2018 reference cycles. |
| `replication_through_2014/manuscript_audit/raw_baseline_covariate_check/L40_B.xpt`, `BMX_B.xpt`, `L13_B.xpt` | Native 2001–2002 serum creatinine, BMI, and total cholesterol. |
| `replication_through_2014/manuscript_audit/raw_fish_zero_check/DRXTOT_B.xpt`, `DR1TOT_C.xpt` through `DR1TOT_H.xpt` | Native fish-detail refusal coding used in the audited correction. |
| `replication_through_2014/manuscript_completion/data/analytic_corrected.rds` | Required validation comparator after reconstructed inputs have been saved; does not supply rebuilt exposure or covariate values. |
| `replication_through_2014/data/full_designs.rds` | Required validation comparator for survey-design IDs and weight vectors. |

A comma-separated filename list in a row uses the same directory for every listed file. All `.xpt` filenames use the lowercase extension shown because the scripts read that exact spelling.

The inherited RDS must retain its original column types and factor levels. Among its required fields are `SEQN`, `SDDSRVYR`, demographic factors, BMI, calibrated serum creatinine, cholesterol, smoking pack-years, alcohol and physical-activity categories, fish variables, EPA/DHA fields, urinary creatinine/thallium, and the archived blood-log columns. The scripts also preserve additional inherited columns. The complete object is required, not only the named examples; original builders for several of these fields are unavailable. See [input_lineage.md](review/input_lineage.md) and [numeric_precision.md](review/numeric_precision.md).

## Control and comparison stages

The driver always attempts the historical control replay and comparison, including with `--scope primary` or `--points-only`.

| Required path | Use |
|---|---|
| `replication_through_2014/remote_execution/models/<scenario>/specification.rds` | Frozen primary control formulas and participant order. |
| `replication_through_2014/remote_execution/models/<scenario>/imputations.rds` | Five completed primary control datasets. |
| `replication_through_2014/remote_execution/models/<scenario>/replicate_design.rds` | Control participant order and sampling weights; control fitting is point-only. |
| `replication_through_2014/official_crq_builtin_20261029/results_by_imputation.csv` | Direct-package primary comparator used by reporting and comparison plots. |
| `replication_through_2014/remote_execution/results/panel_results.csv` | Historical panel comparator used once fresh fits exist; only point comparisons are made because inference differs. |

For the three control RDS paths, `<scenario>` is each of `primary_Pb`, `primary_Hg`, `egfr60_Cd`, and `egfr60_As` (12 files total). An optional `replication_through_2014/remote_execution/models/<scenario>/preimputation_frame.rds` is checked when present for any of the 23 scenarios. It is not used to fill the new imputed data.

## Internal outputs

`chemicals.rds`, `mortality.rds`, `covariates.rds`, `analytic_fresh.rds`, `full_designs.rds`, the scenario specifications, and the new MICE datasets are produced inside `metal_survival_analysis/` by the documented stage order. They are not additional downloaded inputs. All R helper scripts needed by this package are included under `R/`; no external R source file is required.
