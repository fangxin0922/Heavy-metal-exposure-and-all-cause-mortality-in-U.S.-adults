# Prepare study inputs

Run these commands from the repository root. Set `HM_R_LIBRARY` for a separate R library, `HM_DATA_DIR` for main inputs, and `EMM_DATA_DIR` for biomarker inputs. Defaults are `data/metal_survival/` and `data/biomarker_emm/`.

## 1. Supply the original reference archive

```sh
Rscript scripts/stage_reference_input.R /path/to/reference_analytic.rds
```

This validates and copies the retained **pre-correction, 1999–2018 study archive** without changing values, row order, factor levels or numeric representations. It must contain all ten cycles; the urinary-creatinine reference model needs cycles 9–10 even though outcome models use cycles 1–8. A corrected through-2014 dataset cannot replace this input.

**Remaining prerequisite:** the complete original archive-construction, questionnaire harmonization, physical-activity/MET, fish food-code and dietary recall-day selection builders have not been recovered. These scripts do not reconstruct that source archive from raw NHANES, infer missing rules, or substitute imputed values. The archive's original values and missingness are required. `input_manifest.R` lists the required columns and validates their structure. The separate `metal_survival_analysis/R/participant_flow.R` script reconstructs the current Figure 1 selection counts from native releases.

## 2. Retrieve public native files

```sh
Rscript scripts/fetch_public_inputs.R --biomarkers
```

The manifest covers demographic and laboratory XPT files, baseline corrections, dietary files, the **2019 public linked-mortality release**, and biomarker/CBC files. Omit `--biomarkers` for the main analysis only. Files are downloaded from CDC, parsed before use, and recorded with checksums. Existing files are validated and retained.

To use files already downloaded instead, supply their parent directory:

```sh
Rscript scripts/fetch_public_inputs.R --biomarkers /path/to/native_files
```

The search is recursive and case-insensitive. Conflicting copies of the same filename cause an error; use a curated directory in that case. Input data and manifests under `data/` are ignored by Git.

## 3. Build the corrected data

The normal main-analysis runner executes these stages automatically. To prepare biomarker inputs without running the main outcome models, run only:

```sh
Rscript metal_survival_analysis/R/exposure_records.R
Rscript metal_survival_analysis/R/survival_records.R
Rscript metal_survival_analysis/R/adjustment_data.R
Rscript metal_survival_analysis/R/cohort_builder.R
```

These stages rebuild the available native measurements, apply documented corrections, and write `metal_survival_analysis/data/analytic.rds` and `full_designs.rds`. Existing cohort-count checks remain active. Use a clean checkout/output directory for a new full main-analysis run; the main runner rejects populated output folders.

## 4. Generate biomarker parent inputs

```sh
Rscript biomarker_emm/R/prepare_parent_inputs.R
```

This writes the observed corrected `analytic.rds` and each primary scenario's `specification.rds` and `replicate_design.rds` under `EMM_DATA_DIR`. It draws 400 Rao–Wu subbootstrap replicates, with seed 20260101 and an explicit RNG, from each full component **before** selecting the analysis domain. It checks the weight/variance calculation against `survey::withReplicates()` and records source hashes and runtime versions. Existing parent inputs are never overwritten.

The biomarker analysis retains complete-case selection and applies its positive weight rescaling later in `prepare_models.R`. The [biomarker guide](../biomarker_emm/README.md) gives the remaining commands.

## Validation

```sh
Rscript --vanilla tests/test_preparation.R
python3 tests/test_runner.py
```

These synthetic checks cover the renal boundary, restricted-cohort knots, fresh spline imputation, summary seeds, full-component survey draws, generated parent inputs, dependency scheduling, and overwrite guards. They do not refit the study's outcome models.
