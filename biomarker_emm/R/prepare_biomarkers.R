# Prepare serum and RBC biomarker cohorts for the CQR interaction analysis.

source("biomarker_emm/R/settings.R")
suppressPackageStartupMessages({
  library(haven)
  library(tidyverse)
})

stopifnot(Sys.info()[["sysname"]] == "Linux")
dir.create(file.path(OUTPUT_DIR, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "validation"), showWarnings = FALSE)

SCENARIOS <- c(Pb = "primary_Pb", Hg = "primary_Hg", Cd = "egfr60_Cd", As = "egfr60_As")
METALS <- c(Pb = "LBXBPB", Hg = "LBXTHG", Cd = "URXUCD", As = "VNSUMARS")
UNITS <- c(vitamin_c = "mg/dL", alpha_carotene = "ug/dL", rbc_folate = "nmol/L")
source_labels <- list()

# Biomarker files and published assay bridges ----

registry <- data.frame(
  biomarker = c(rep("vitamin_c", 2), rep("alpha_carotene", 3), rep("rbc_folate", 8)),
  cycle = c(3:4, 2:4, 1:8),
  file = c("L06VIT_C", "VIC_D", "L06VIT_B", "L45VIT_C", "VITAEC_D",
           "LAB06", "L06_B", "L06NB_C", "FOLATE_D", "FOLATE_E", "FOLATE_F", "FOLATE_G", "FOLATE_H"),
  variable = c(rep("LBXVIC", 2), rep("LBXALC", 3), rep("LBDRBFSI", 4),
               rep("LBXRBFSI", 2), rep("LBDRFOSI", 2)))
registry$year <- 1997 + 2 * registry$cycle
registry$units <- UNITS[registry$biomarker]
registry$url <- sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%s/DataFiles/%s.htm",
                        registry$year, registry$file)
CBC_FILES <- c("LAB25", "L25_B", "L25_C", "CBC_D")

read_native <- function(path, fields, cycle) {
  x <- read_xpt(path)
  stopifnot(!anyDuplicated(x$SEQN), all(fields %in% names(x)))
  for (field in fields) {
    source_labels[[length(source_labels) + 1L]] <<- data.frame(
      file = basename(path), cycle, variable = field,
      label = as.character(attr(x[[field]], "label")))
  }
  return(x)
}

adjust_folate <- function(rbc, serum, hct) {
  # CDC FOLATE_F: Bio-Rad 1999-2006 -> microbiologic-assay equivalent.
  # Serum must exceed 1 nmol/L so log10(serum)^(-1/2) is defined.
  valid <- is.finite(rbc) & rbc > 0 & is.finite(serum) & serum > 1 &
    is.finite(hct) & hct > 0 & hct < 100
  whole_blood <- whole_adjusted <- serum_adjusted <- rbc_adjusted <- rep(NA_real_, length(rbc))
  fraction <- hct[valid] / 100
  whole_blood[valid] <- rbc[valid] * fraction + serum[valid] * (1 - fraction)
  whole_adjusted[valid] <- 10^(0.2204 + 1.017 * log10(whole_blood[valid]))
  x <- log10(serum[valid])
  serum_adjusted[valid] <- 10^(0.0188 * x^3 - 2.7109 * x^(-1/2) + 3.8276)
  rbc_adjusted[valid] <- (whole_adjusted[valid] - serum_adjusted[valid] * (1 - fraction)) / fraction
  valid <- valid & is.finite(rbc_adjusted) & rbc_adjusted > 0
  rbc_adjusted[!valid] <- NA_real_
  return(data.frame(biomarker_value = rbc_adjusted, calibration_valid = valid,
    whole_blood_folate_raw_nmol_l = whole_blood,
    whole_blood_folate_adjusted_nmol_l = whole_adjusted,
    serum_folate_adjusted_nmol_l = serum_adjusted))
}

# Check against the CDC worked example before processing participant records.
example <- adjust_folate(378.3, 30.4, 41.4)
expected_rbc <- (316.323 - 45.9955 * (1 - 0.414)) / 0.414
stopifnot(abs(example$whole_blood_folate_raw_nmol_l - 174.4306) < 1e-7,
          abs(example$whole_blood_folate_adjusted_nmol_l - 316.323) < 0.002,
          abs(example$serum_folate_adjusted_nmol_l - 45.9955) < 0.002,
          abs(example$biomarker_value - expected_rbc) < 0.005)
calibration_checks <- data.frame(
  check = c("CDC whole-blood example", "CDC adjusted whole-blood example", "CDC adjusted serum example", "CDC derived RBC example"),
  calculated = c(example$whole_blood_folate_raw_nmol_l, example$whole_blood_folate_adjusted_nmol_l,
                 example$serum_folate_adjusted_nmol_l, example$biomarker_value),
  expected = c(174.4306, 316.323, 45.9955, expected_rbc))

marker_parts <- list()
calibration_summary <- list()
for (i in seq_len(nrow(registry))) {
  item <- registry[i, ]
  fields <- item$variable
  if (item$biomarker == "rbc_folate" && item$cycle <= 4) fields <- c(fields, "LBDFOLSI")
  x <- read_native(file.path(MARKER_DIR, paste0(item$file, ".xpt")), fields, item$cycle)
  marker <- data.frame(SEQN = as.numeric(x$SEQN), biomarker_cycle = item$cycle,
    biomarker_raw = as.numeric(x[[item$variable]]), biomarker_value = as.numeric(x[[item$variable]]),
    biomarker_source_file = paste0(item$file, ".xpt"), biomarker_source_variable = item$variable,
    calibration_method = "released value", calibration_valid = TRUE,
    serum_folate_raw_nmol_l = NA_real_, hematocrit_pct = NA_real_,
    whole_blood_folate_raw_nmol_l = NA_real_, whole_blood_folate_adjusted_nmol_l = NA_real_,
    serum_folate_adjusted_nmol_l = NA_real_, hematocrit_source_file = NA_character_)
  if (item$biomarker == "alpha_carotene" && item$cycle == 3) {
    marker$biomarker_value <- 0.9837 * marker$biomarker_raw + 0.2137
    marker$calibration_method <- "Craft to CDC: combined Deming 0.9837*x+0.2137 ug/dL"
  }
  if (item$biomarker == "rbc_folate" && item$cycle <= 4) {
    hct_file <- paste0(CBC_FILES[item$cycle], ".xpt")
    cbc <- read_native(file.path(SUPPORT_DIR, hct_file), "LBXHCT", item$cycle)
    marker$hematocrit_pct <- as.numeric(cbc$LBXHCT[match(marker$SEQN, cbc$SEQN)])
    marker$serum_folate_raw_nmol_l <- as.numeric(x$LBDFOLSI)
    marker$hematocrit_source_file <- hct_file
    adjusted <- adjust_folate(marker$biomarker_raw, marker$serum_folate_raw_nmol_l, marker$hematocrit_pct)
    for (field in names(adjusted)) marker[[field]] <- adjusted[[field]]
    marker$calibration_method <- "CDC four-step Bio-Rad to microbiologic-assay equivalent"
  }
  if (item$biomarker == "rbc_folate" && item$cycle == 7) {
    marker$calibration_method <- "released LBDRFOSI; February 2019 CDC revision already applied"
  }
  marker$calibration_valid <- marker$calibration_valid & is.finite(marker$biomarker_value) &
    marker$biomarker_value > 0 & is.finite(marker$biomarker_raw) & marker$biomarker_raw > 0
  marker$biomarker_value[!marker$calibration_valid] <- NA_real_
  observed <- is.finite(marker$biomarker_raw) & marker$biomarker_raw > 0
  calibrated <- marker$calibration_valid
  calibration_summary[[length(calibration_summary) + 1L]] <- data.frame(
    biomarker = item$biomarker, cycle = item$cycle, file = item$file,
    native_n = nrow(marker), native_observed_n = sum(observed), calibrated_n = sum(calibrated),
    observed_but_not_calibrated_n = sum(observed & !calibrated),
    raw_median = median(marker$biomarker_raw[observed]),
    calibrated_median = median(marker$biomarker_value[calibrated]),
    calibrated_min = min(marker$biomarker_value[calibrated]),
    calibrated_max = max(marker$biomarker_value[calibrated]),
    calibration = unique(marker$calibration_method))
  marker_parts[[paste(item$biomarker, item$cycle, sep = "_")]] <- marker
}

# Native weights for the shortened vitamin C and alpha-carotene windows ----

native_weights <- list()
for (cycle in 2:4) {
  file <- paste0("DEMO_", LETTERS[cycle], ".xpt")
  x <- read_native(file.path(LAB_DIR, file), "WTMEC2YR", cycle)
  native_weights[[paste0("blood_", cycle)]] <- data.frame(SEQN = as.numeric(x$SEQN), weight = as.numeric(x$WTMEC2YR))
  file <- c("L06HM_B", "L06HM_C", "UHM_D")[cycle - 1]
  variable <- if (cycle == 2) "WTSHM2YR" else "WTSA2YR"
  x <- read_native(file.path(LAB_DIR, paste0(file, ".xpt")), variable, cycle)
  native_weights[[paste0("Cd_", cycle)]] <- data.frame(SEQN = as.numeric(x$SEQN), weight = as.numeric(x[[variable]]))
  if (cycle >= 3) {
    file <- if (cycle == 3) "L06UAS_C" else "UAS_D"
    x <- read_native(file.path(LAB_DIR, paste0(file, ".xpt")), "WTSA2YR", cycle)
    native_weights[[paste0("As_", cycle)]] <- data.frame(SEQN = as.numeric(x$SEQN), weight = as.numeric(x$WTSA2YR))
  }
}

# Restore observed parent records; require only analysis-variable completeness ----

archive <- readRDS(INPUT_FILE)
stopifnot(!anyDuplicated(archive$SEQN))
archive$cycle_f <- factor(archive$SDDSRVYR)
flow <- by_cycle <- missingness <- weight_checks <- pair_registry <- list()
for (metal in names(SCENARIOS)) {
  spec_file <- file.path(PARENT_DIR, SCENARIOS[[metal]], "specification.rds")
  bank_file <- file.path(PARENT_DIR, SCENARIOS[[metal]], "replicate_design.rds")
  spec <- readRDS(spec_file)
  bank <- readRDS(bank_file)
  j <- match(bank$row_ids, archive$SEQN)
  stopifnot(!anyNA(j), identical(as.numeric(bank$row_ids), as.numeric(spec$row_ids)))
  parent <- archive[j, ]
  stopifnot(nrow(parent) == spec$n, sum(parent$event) == spec$events,
            all(parent$RIDAGEYR >= 40), all(parent$SDDSRVYR <= 8), all(parent$ELIGSTAT == 1))
  if (metal %in% c("Cd", "As")) stopifnot(all(parent$egfr_ckdepi_2021 >= 60))
  parent$parent_sampling_weight <- as.numeric(bank$sampling_weights)
  parent$metal_log2 <- parent[[paste0("log_primary_", METALS[[metal]])]]

  for (biomarker in names(UNITS)) {
    marker_cycles <- registry$cycle[registry$biomarker == biomarker]
    cycles <- intersect(marker_cycles, sort(unique(parent$SDDSRVYR)))
    d <- parent %>% filter(SDDSRVYR %in% cycles)
    marker <- do.call(rbind, marker_parts[paste(biomarker, cycles, sep = "_")])
    stopifnot(!anyDuplicated(marker$SEQN))
    j <- match(d$SEQN, marker$SEQN)
    d$has_native_record <- !is.na(j)
    for (field in setdiff(names(marker), "SEQN")) d[[field]] <- marker[[field]][j]
    stopifnot(all(d$SDDSRVYR[d$has_native_record] == d$biomarker_cycle[d$has_native_record]))
    d$analysis_weight <- d$parent_sampling_weight
    d$weight_source <- "parent pooled component weight"
    d$weight_divisor <- if (metal %in% c("Pb", "Cd")) 8 else 6
    if (biomarker != "rbc_folate") {
      for (cycle in cycles) {
        selected <- d$SDDSRVYR == cycle
        w <- native_weights[[paste0(if (metal %in% c("Pb", "Hg")) "blood" else metal, "_", cycle)]]
        d$analysis_weight[selected] <- w$weight[match(d$SEQN[selected], w$SEQN)] / length(cycles)
        d$weight_source[selected] <- if (metal %in% c("Pb", "Hg")) "WTMEC2YR" else if (cycle == 2) "WTSHM2YR" else "WTSA2YR"
        d$weight_divisor[selected] <- length(cycles)
      }
    } else {
      d$weight_source <- if (metal %in% c("Pb", "Hg")) {
        ifelse(d$SDDSRVYR <= 2, "WTMEC4YR*2/8", ifelse(d$SDDSRVYR == 8, "WTSH2YR/K", "WTMEC2YR/K"))
      } else if (metal == "Cd") {
        ifelse(d$SDDSRVYR <= 2, "WTSHM4YR*2/8", "WTSA2YR/8")
      } else "WTSA2YR/6"
    }
    covariates <- spec$covariates
    complete_case_variables <- unique(c("followup_yrs", "event", "metal_log2", "biomarker_value",
      "analysis_weight", "SDMVSTRA", "SDMVPSU", covariates))
    stopifnot(all(complete_case_variables %in% names(d)))
    # Infinite or invalid analysis values are missing; ancillary fields never drive drop_na.
    for (field in complete_case_variables) if (is.numeric(d[[field]])) d[[field]][!is.finite(d[[field]])] <- NA_real_
    d$analysis_weight[d$analysis_weight <= 0 & !is.na(d$analysis_weight)] <- NA_real_
    d$followup_yrs[d$followup_yrs <= 0 & !is.na(d$followup_yrs)] <- NA_real_
    for (field in complete_case_variables) missingness[[length(missingness) + 1L]] <- data.frame(
      metal, biomarker, variable = field, n = nrow(d), missing_n = sum(is.na(d[[field]])))
    complete <- d %>% drop_na(all_of(complete_case_variables))
    for (field in covariates) if (is.factor(complete[[field]])) complete[[field]] <- droplevels(complete[[field]])
    stopifnot(!anyDuplicated(complete$SEQN), all(complete$event %in% 0:1),
              all(complete$biomarker_value > 0), all(complete$analysis_weight > 0))
    flow[[length(flow) + 1L]] <- data.frame(metal, biomarker,
      cycles = paste(cycles, collapse = ";"), parent_n = nrow(parent), eligible_cycles_n = nrow(d),
      native_record_n = sum(d$has_native_record), native_observed_n = sum(is.finite(d$biomarker_raw) & d$biomarker_raw > 0),
      calibrated_n = sum(is.finite(d$biomarker_value)), positive_weight_n = sum(is.finite(d$analysis_weight)),
      complete_case_n = nrow(complete), deaths = sum(complete$event),
      excluded_after_calibration_n = sum(is.finite(d$biomarker_value)) - nrow(complete))
    for (cycle in cycles) {
      z <- complete[complete$SDDSRVYR == cycle, ]
      b <- d[d$SDDSRVYR == cycle, ]
      by_cycle[[length(by_cycle) + 1L]] <- data.frame(metal, biomarker, cycle,
        eligible_n = nrow(b), native_observed_n = sum(is.finite(b$biomarker_raw) & b$biomarker_raw > 0),
        calibrated_n = sum(is.finite(b$biomarker_value)), complete_case_n = nrow(z), deaths = sum(z$event))
      ratio <- z$analysis_weight / z$parent_sampling_weight
      weight_checks[[length(weight_checks) + 1L]] <- data.frame(metal, biomarker, cycle,
        n = nrow(z), source = paste(unique(z$weight_source), collapse = ";"),
        divisor = unique(z$weight_divisor), ratio_to_parent_min = min(ratio), ratio_to_parent_max = max(ratio),
        weight_min = min(z$analysis_weight), weight_max = max(z$analysis_weight), weight_sum = sum(z$analysis_weight))
    }
    keep <- unique(c("SEQN", "SDDSRVYR", complete_case_variables, setdiff(names(marker), "SEQN"),
      "has_native_record", "parent_sampling_weight", "weight_source", "weight_divisor"))
    prepared <- list(data = complete[keep], covariates = covariates,
      complete_case_variables = complete_case_variables, metal = metal, biomarker = biomarker,
      units = UNITS[[biomarker]], parent_scenario = SCENARIOS[[metal]], cycles = cycles,
      calibration = registry[registry$biomarker == biomarker & registry$cycle %in% cycles, ],
      calibration_note = if (biomarker == "rbc_folate") "1999-2006 calibrated with observed serum folate and hematocrit; 2011-2012 revised release used once" else if (biomarker == "alpha_carotene") "2003-2004 Craft values mapped to CDC scale with combined Deming equation" else "Released serum vitamin C mg/dL",
      missing_data_method = "tidyr::drop_na on explicit analysis variables; no imputation",
      parent_variance_metadata = bank[c("scale", "rscales", "mse", "R", "seed", "complete_data_df", "df_method")])
    filename <- paste0(metal, "_", biomarker, ".rds")
    stopifnot(!file.exists(file.path(OUTPUT_DIR, "data", filename)))
    saveRDS(prepared, file.path(OUTPUT_DIR, "data", filename))
    pair_registry[[length(pair_registry) + 1L]] <- data.frame(metal, biomarker, file = filename,
      n = nrow(complete), deaths = sum(complete$event), units = UNITS[[biomarker]],
      first_cycle = min(cycles), last_cycle = max(cycles), parent_scenario = SCENARIOS[[metal]])
  }
}

# Cohort and calibration checks ----

write.csv(registry, file.path(OUTPUT_DIR, "biomarker_registry.csv"), row.names = FALSE)
write.csv(calibration_checks, file.path(OUTPUT_DIR, "validation", "calibration_example_checks.csv"), row.names = FALSE)
for (item in c("flow", "by_cycle", "missingness", "weight_checks", "pair_registry", "calibration_summary", "source_labels")) {
  write.csv(do.call(rbind, get(item)), file.path(OUTPUT_DIR, "validation", paste0(item, ".csv")), row.names = FALSE)
}
