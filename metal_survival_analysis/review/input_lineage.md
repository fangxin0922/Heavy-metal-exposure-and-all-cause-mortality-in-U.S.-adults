Reconstruction starts with the retained, pre-correction 1999–2018 analytic archive. Its population-wide screening builder is unavailable; no new claim of a full raw-NHANES cohort is made.
Native DEMO restores age, sex, race, education, marital status, PIR and design identifiers; factors are rebuilt with the stored reference levels.
Native 2001–2002 baseline serum creatinine, BMI and cholesterol are read directly. Other cycles retain supplied baseline BMI/cholesterol/calibrated serum creatinine. Race-free CKD-EPI 2021 eGFR is recalculated.
Pack-years, alcohol, physical-activity categories and fish recall gram totals are inherited with their original missingness. The original questionnaire harmonization, MET builder and food-code lookup are unavailable.
Fish keeps the stored screener-first/recall-second hierarchy and the exact /4.345 and /85*7 conversions. Native DRD370V invalidates the documented species-detail refusal. Observed zeros are not automatically reclassified as missing.
EPA+DHA sensitivity fields DRXTP205 and DRXTP226 retain inherited recall-day selection; their original day-selection builder is unavailable.
The old corrected analytic dataset is read only by stage 04 for validation. It supplies no rebuilt analysis values.
