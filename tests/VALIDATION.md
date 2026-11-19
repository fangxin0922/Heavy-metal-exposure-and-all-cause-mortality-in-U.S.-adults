# Validation — 19 November 2026

- `tests/test_preparation.R` passed on synthetic data: eGFR boundary exclusions, fresh spline imputations, restricted-cohort knots, corrected summary seeds, survey draws generated before domain selection, generated parent inputs, and overwrite guards.
- All four `tests/test_runner.py` tests passed, including blocking the corrected spline when primary Cd preparation fails.
- The reference staging script and all four main input stages completed using the retained local archive and native files. All 24 checked outcome-model input columns agreed with the corrected reference within `1e-12`.
- The corrected cadmium spline selected exactly the primary Cd participants: **7,061 participants, 1,510 deaths**.
- All four generated parent banks reproduced the saved historical full-component sampling weights and the 400 replicate multipliers for the selected rows with **zero difference**. The variance scales also matched. Cohort sizes were Pb 24,021; Hg 18,258; Cd 7,061; As 5,552.
- All 76 native input targets were parsed and validated; the biomarker downloader was exercised against CDC.

Local validation used R 4.6.1, survey 4.5, mice 3.19.0 and quantreg 6.1. Participant data, fitted objects and local validation artifacts are excluded from the repository.

These checks validate input preparation and cohort selection. They do not constitute a new full outcome-model run or recover the missing original archive-construction rules. The retained pre-correction reference archive remains required.
