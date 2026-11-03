# Rebuild available covariates from native files. Retain inherited questionnaire
# and dietary fields, including missingness, where original builders are unavailable.
source("metal_survival_analysis/R/settings.R")
source(file.path(RUN_DIR, "R/data_tools.R"))

inherited <- read_inherited() %>% dplyr::filter(SDDSRVYR <= 8)
demo_fields <- c("RIDAGEYR", "RIAGENDR", "RIDRETH1", "DMDEDUC2", "DMDMARTL",
                 "INDFMPIR", "SDMVSTRA", "SDMVPSU")
demo_files <- file.path(raw_directory, paste0(c("DEMO", paste0("DEMO_", LETTERS[2:8])), ".xpt"))
demographics <- lapply(seq_along(demo_files), function(cycle) {
  read_native(demo_files[cycle]) %>%
    dplyr::select(SEQN, dplyr::all_of(demo_fields)) %>%
    dplyr::mutate(SDDSRVYR=cycle)
}) %>% dplyr::bind_rows()
covariates <- inherited %>%
  dplyr::select(-dplyr::all_of(demo_fields)) %>%
  dplyr::left_join(demographics, by=c("SEQN", "SDDSRVYR")) %>%
  dplyr::mutate(
    RIAGENDR_f=factor(RIAGENDR, levels=1:2, labels=levels(inherited$RIAGENDR_f)),
    RIDRETH1_f=factor(RIDRETH1, levels=1:5, labels=levels(inherited$RIDRETH1_f)),
    DMDEDUC2_f=factor(DMDEDUC2, levels=1:5, labels=levels(inherited$DMDEDUC2_f)),
    marital_3cat=factor(dplyr::case_when(
      DMDMARTL %in% c(1, 6) ~ "married_or_cohabiting",
      DMDMARTL %in% 2:4 ~ "widowed_divorced_separated",
      DMDMARTL == 5 ~ "never_married", TRUE ~ NA_character_),
      levels=levels(inherited$marital_3cat)))

baseline_directory <- file.path(PROJECT_ROOT, "replication_through_2014/manuscript_audit/raw_baseline_covariate_check")
baseline_files <- file.path(baseline_directory, c("L40_B.xpt", "BMX_B.xpt", "L13_B.xpt"))
serum <- read_native(baseline_files[1])
bmi <- read_native(baseline_files[2])
cholesterol <- read_native(baseline_files[3])
i <- which(covariates$SDDSRVYR == 2)
covariates$LBXSCR[i] <- serum$LBDSCR[match(covariates$SEQN[i], serum$SEQN)]
covariates$LBXSCR_cal[i] <- covariates$LBXSCR[i]
covariates$BMXBMI[i] <- bmi$BMXBMI[match(covariates$SEQN[i], bmi$SEQN)]
covariates$LBXTC[i] <- cholesterol$LBXTC[match(covariates$SEQN[i], cholesterol$SEQN)]
female <- covariates$RIAGENDR == 2
kappa <- ifelse(female, .7, .9)
alpha <- ifelse(female, -.241, -.302)
covariates$egfr_ckdepi_2021 <- 142 * pmin(covariates$LBXSCR_cal/kappa, 1)^alpha *
  pmax(covariates$LBXSCR_cal/kappa, 1)^(-1.2) * .9938^covariates$RIDAGEYR *
  ifelse(female, 1.012, 1)

# Correct one false fish-screener zero using the native refusal code.
# Retain audited fish totals and missingness; the original food-code map is unavailable.
fish_files <- file.path(PROJECT_ROOT, "replication_through_2014/manuscript_audit/raw_fish_zero_check",
                        c("DRXTOT_B.xpt", paste0("DR1TOT_", LETTERS[3:8], ".xpt")))
covariates$fish_detail_refusal_code_native <- NA_real_
for (k in seq_along(fish_files)) {
  native <- read_native(fish_files[k])
  i <- which(covariates$SDDSRVYR == k + 1)
  covariates$fish_detail_refusal_code_native[i] <- native$DRD370V[match(covariates$SEQN[i], native$SEQN)]
}
covariates$fish_detail_refused_native <- ifelse(is.na(covariates$fish_detail_refusal_code_native),
                                              NA, covariates$fish_detail_refusal_code_native == 1)
refused <- !is.na(covariates$fish_detail_refused_native) & covariates$fish_detail_refused_native
stopifnot(identical(as.numeric(covariates$SEQN[refused]), 65758),
          all(is.na(covariates$fish_servings_wk_t2[refused])))
for (field in c("dark_fish_30d_anchored", "fish_servings_wk_t1", "fish_servings_wk"))
  covariates[[field]][refused] <- NA_real_
covariates$fish_source[refused] <- "missing"
covariates <- droplevels(covariates)
saveRDS(covariates, data_file("covariates.rds"))
saveRDS(demographics, data_file("demographics.rds"))
record_input_sources(c(inherited_file, demo_files, baseline_files, fish_files), "03_covariates")
writeLines(c(
  "Reconstruction starts with the retained, pre-correction 1999–2018 analytic archive. Its population-wide screening builder is unavailable; no new claim of a full raw-NHANES cohort is made.",
  "Native DEMO restores age, sex, race, education, marital status, PIR and design identifiers; factors are rebuilt with the stored reference levels.",
  "Native 2001–2002 baseline serum creatinine, BMI and cholesterol are read directly. Other cycles retain supplied baseline BMI/cholesterol/calibrated serum creatinine. Race-free CKD-EPI 2021 eGFR is recalculated.",
  "Pack-years, alcohol, physical-activity categories and fish recall gram totals are inherited with their original missingness. The original questionnaire harmonization, MET builder and food-code lookup are unavailable.",
  "Fish keeps the stored screener-first/recall-second hierarchy and the exact /4.345 and /85*7 conversions. Native DRD370V invalidates the documented species-detail refusal. Observed zeros are not automatically reclassified as missing.",
  "EPA+DHA sensitivity fields DRXTP205 and DRXTP226 retain inherited recall-day selection; their original day-selection builder is unavailable.",
  "The old corrected analytic dataset is read only by stage 04 for validation. It supplies no rebuilt analysis values."
), audit_file("INHERITED_VARIABLES.md"))
message("Covariate reconstruction complete: ", nrow(covariates), " retained archive rows.")
