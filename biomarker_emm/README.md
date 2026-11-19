# Biomarker effect-measure modification

[Repository overview](../README.md)

Nutrient-biomarker preparation, continuous interactions, checked Peng–Huang CQR, survey replicate inference, and figures.

## Requirements

Use Linux with single-threaded OpenBLAS 0.3.26, LAPACK 3.12.0, and the R environment and packages in the [main analysis guide](../metal_survival_analysis/README.md#run) plus `tidyverse` 2.0.0. The separate figure script requires ggplot2 3.5.2 and uses `readr` 2.2.0. A Fortran compiler and BLAS/LAPACK development libraries are needed for the supplied native solver. Package, compiler and numerical-library differences can affect CQR estimates; matching package names alone does not guarantee identical results.

## Inputs

Supply the following under `data/biomarker_emm/`, or set `EMM_DATA_DIR` to their parent directory:

- `analytic.rds`: the observed corrected study frame.
- `parent_models/{primary_Pb,primary_Hg,egfr60_Cd,egfr60_As}/{specification.rds,replicate_design.rds}`: parent model specifications and the fixed 400 survey draws.
- `raw/labs/`, `raw/biomarkers/` and `raw/cbc/`: the XPT files listed in `prepare_biomarkers.R`.

Outputs go to `biomarker_emm/output/`, configurable through `EMM_OUTPUT_DIR`. The `original` stage fits the full sample; complete it before the replicate stage. The main pipeline does not create the required parent survey-draw banks. These fixed banks and both study input frames are prerequisites for reproducing the reported analysis.

## Run

Run from the repository root, completing the `original` stage before `replicates`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
(cd biomarker_emm/native && R CMD SHLIB ph_diagnostic.f -o ph_diagnostic.so -llapack -lblas)
Rscript biomarker_emm/R/prepare_biomarkers.R
Rscript biomarker_emm/R/prepare_models.R
python3 biomarker_emm/run_models.py original 12
python3 biomarker_emm/run_models.py replicates 24
Rscript biomarker_emm/R/summarize_results.R
Rscript biomarker_emm/R/plot_interactions.R
```

## Interpretation

EMM uses observed complete cases, continuous biomarker interactions and survey replicate variance. Available intervals remain provisional because numerical solver checks were not fully satisfied. An unavailable interval is not a null finding. CQR effects are conditional survival-time differences in years, not mean years of life lost or causal effects.

## Attribution

The native solver and R interface derive from [quantreg](https://github.com/cran/quantreg), by Roger Koenker and contributors, including Stephen Portnoy for censored QR. Study modifications add 33 native diagnostic fields and path-validity checks. The quantreg-derived code retains [GPL version 2 or later](../licenses/GPL-2.txt).
