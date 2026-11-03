# Third-party sources and notices

## Reference workflow

The main analysis draws on YZhang9719's `Plasticizer-mixture-and-mortality`, pinned to commit [`da34a8bdeff9d665da38f7f07c853aaaf71bc5c3`](https://github.com/YZhang9719/Plasticizer-mixture-and-mortality/tree/da34a8bdeff9d665da38f7f07c853aaaf71bc5c3). Its copyright and MIT permission notice are retained in [reference-MIT.txt](metal_survival_analysis/licenses/reference-MIT.txt). [SOURCES.md](metal_survival_analysis/SOURCES.md) explains which preparation and supplemental-model functions informed this study.

## quantreg native solver and R interface

The supplied `replication_through_2014/remote_execution/validation/instrumented_ph/src/ph_diagnostic.f` derives from quantreg 6.1's [`crqfnb.f`](https://github.com/cran/quantreg/blob/6.1/src/crqfnb.f) and [`rqfnb.f`](https://github.com/cran/quantreg/blob/6.1/src/rqfnb.f). The R native-call interface in `replication_through_2014/R/revised_cqr_solver.R` follows quantreg's `crq.fit.pen` conventions; `checked_ph_solver.R` adds the study's validity checks.

quantreg 6.1 is maintained and authored by Roger Koenker, with Stephen Portnoy credited for contributions to censored QR code and other contributors listed in the upstream [DESCRIPTION](https://github.com/cran/quantreg/blob/6.1/DESCRIPTION). Its declared license is **GPL (>= 2)**. The quantreg-derived source retains that license; it is not covered by the reference workflow's MIT notice. Copies of [GPL version 2](licenses/GPL-2.txt) and [GPL version 3](licenses/GPL-3.txt) are included.

The study modifications rename native entry points to `crqdiag`, `rqdiag`, `lpdiag` and `stepdiag`, expose 33 diagnostic fields, and retain native status information for the R acceptance checks. The audited source used for the November 3, 2026 EMM run is supplied directly. It retains the native tolerance and iteration limit. The accepted-path checks, fixed retry policy and missing-result handling are described in the EMM guide and executable helpers.

The upstream Ratfor-generated source comments are retained. The dated study modifications are distributed with their source; no compiled native library is included.
