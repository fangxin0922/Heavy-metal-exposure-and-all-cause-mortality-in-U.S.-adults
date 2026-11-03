# Software environment

Run all commands from the repository root. The completed analyses used Linux with R 4.6.1, Ubuntu 24.04.5, OpenBLAS 0.3.26 and LAPACK 3.12.0. EMM fitting explicitly checks Linux, OpenBLAS and single-threaded execution. The Python runners use the standard library; the main runner uses POSIX `fcntl`.

| Package | Recorded version | Use |
|---|---|---|
| quantreg | 6.1 | CQR fitting and coefficient extraction |
| survival | 3.8-6 | Survival responses and supplemental Cox models |
| survey | 4.5 | Survey designs, summaries and replicate variance |
| mice | 3.19.0 | Main-analysis covariate imputation |
| haven | 2.5.5 | NHANES XPT input |
| dplyr | 1.2.1 | Data preparation |
| tidyr | 1.3.2 | Data preparation |
| tidyverse | 2.0.0 | EMM preparation and summary |
| readr | 2.2.0 | EMM input/output |
| Matrix | 1.7-5 | Matrix dependency |
| SparseM | 1.84-2 | CQR dependency |
| MatrixModels | 0.5-4 | CQR dependency |
| ggplot2 | 3.5.2 | EMM figure, in a separate R process |

The EMM summary session loaded ggplot2 4.0.3 through tidyverse. Its separate figure script requires 3.5.2 because the other environment had a recorded S3/S7 dispatch conflict. Follow the EMM guide's plotting-library setting. The main pipeline draws its figures with base R graphics.

Install these packages and their dependencies in a compatible R environment. This version table records the completed run; it is not an `renv` lockfile or a guarantee that a different operating system, BLAS, package build or compiler will produce bitwise-identical estimates. Full recorded sessions are available for the [main control](provenance/main-session.txt), [EMM summary](provenance/emm-summary-session.txt) and [EMM figure](provenance/emm-plot-session.txt).

The publication scripts accept an optional `HM_R_LIBRARY` directory in place of a personal library path. Set it before launching R or Python if your packages are outside the standard R library search path. Existing project-library settings remain documented in the analysis guides.

Before fitting, set:

```sh
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
```

For EMM, compile the supplied audited Fortran source with the R-compatible compiler and BLAS/LAPACK libraries, following its guide. The repository does not distribute a compiled `.so`. Compile the supplied `ph_diagnostic.f` directly; the earlier instrumentation generator is not part of this source release.

Both analyses record source and runtime hashes and reject incompatible saved checkpoints. Preserve those checks when reproducing an analysis.
