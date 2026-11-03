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

baseline_directory <- file.path(INPUT_DIR, "baseline")
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
# Retain supplied fish totals and missingness; the original food-code map is unavailable.
fish_files <- file.path(INPUT_DIR, "diet",
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
message("Covariate reconstruction complete: ", nrow(covariates), " retained archive rows.")
