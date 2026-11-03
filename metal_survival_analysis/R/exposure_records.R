# Rebuild metal concentrations, survey weights and urinary standardization.
# Run from the repository root.
source("metal_survival_analysis/R/settings.R")
source(file.path(RUN_DIR, "R/data_tools.R"))

inherited <- read_inherited()
chemicals <- inherited %>%
  dplyr::filter(SDDSRVYR <= 8) %>%
  dplyr::select(SEQN, SDDSRVYR, URXUCR, URXUTL)
specification_file <- file.path(RUN_DIR, "inputs/weight_specification.csv")
specification <- read.csv(specification_file) %>%
  dplyr::filter(SDDSRVYR <= 8) %>%
  dplyr::mutate(metal_code=unname(c(blood_lead="LBXBPB", blood_mercury="LBXTHG",
                                   urinary_cadmium="URXUCD", urinary_arsenic_sum="VNSUMARS")[metal]),
                K=ifelse(metal_code %in% c("LBXBPB", "URXUCD"), 8L, 6L),
                pooling_multiplier=ifelse(SDDSRVYR <= 2, 2/K, 1/K))
demographic_files <- file.path(raw_directory, paste0(c("DEMO", paste0("DEMO_", LETTERS[2:8])), ".xpt"))
demographics <- lapply(demographic_files, read_native)
full_designs <- setNames(vector("list", 4), metal_names)

for (metal in metal_names) {
  rows <- specification %>% dplyr::filter(metal_code == .env$metal)
  stopifnot(nrow(rows) > 0)
  chemicals[[metal]] <- NA_real_
  chemicals[[paste0("weight_", metal)]] <- NA_real_
  chemicals[[paste0("weight_raw_", metal)]] <- NA_real_
  chemicals[[paste0("weight_source_", metal)]] <- NA_character_
  components <- vector("list", nrow(rows))
  for (k in seq_len(nrow(rows))) {
    row <- rows[k, ]
    path <- file.path(raw_directory, paste0(row$lab_file, ".xpt"))
    lab <- read_native(path)
    demo <- demographics[[row$SDDSRVYR]]
    stopifnot(all(lab$SEQN %in% demo$SEQN))
    if (row$raw_weight %in% names(demo)) {
      component <- demo %>%
        dplyr::transmute(SEQN, SDMVSTRA, SDMVPSU, RIDAGEYR,
                         w_native=.data[[row$raw_weight]]) %>%
        dplyr::filter(is.finite(w_native), w_native >= 0) %>%
        dplyr::mutate(lab_present=SEQN %in% lab$SEQN)
      weight_source <- basename(demographic_files[row$SDDSRVYR])
    } else {
      component <- lab %>%
        dplyr::transmute(SEQN, w_native=.data[[row$raw_weight]]) %>%
        dplyr::filter(is.finite(w_native), w_native >= 0) %>%
        dplyr::left_join(demo %>% dplyr::select(SEQN, SDMVSTRA, SDMVPSU, RIDAGEYR), by="SEQN") %>%
        dplyr::mutate(lab_present=TRUE)
      weight_source <- paste0(row$lab_file, ".xpt")
    }
    component <- component %>% dplyr::mutate(
      SDDSRVYR=row$SDDSRVYR, w_metal=w_native * row$pooling_multiplier,
      weight_variable=row$raw_weight, weight_source_file=sub("[.]xpt$", "", weight_source),
      lab_source_file=row$lab_file, pooling_multiplier=row$pooling_multiplier,
      K=row$K, in_analytic=SEQN %in% chemicals$SEQN)
    stopifnot(!anyNA(component[c("SEQN", "SDMVSTRA", "SDMVPSU")]))
    components[[k]] <- component
    i <- which(chemicals$SDDSRVYR == row$SDDSRVYR)
    j <- match(chemicals$SEQN[i], component$SEQN)
    chemicals[[paste0("weight_", metal)]][i] <- component$w_metal[j]
    chemicals[[paste0("weight_raw_", metal)]][i] <- component$w_native[j]
    chemicals[[paste0("weight_source_", metal)]][i] <- row$raw_weight
    j <- match(chemicals$SEQN[i], lab$SEQN)
    if (metal == "VNSUMARS") {
      species <- c("URXUAS3", "URXUAS5", "URXUMMA", "URXUDMA")
      values <- as.matrix(lab[j, species])
      concentration <- rowSums(values)
      concentration[rowSums(is.finite(values)) != 4] <- NA_real_
      for (v in species) {
        name <- paste0("arsenic_native_", v)
        if (!name %in% names(chemicals)) chemicals[[name]] <- NA_real_
        chemicals[[name]][i] <- lab[[v]][j]
      }
    } else {
      field <- if (metal == "URXUCD" && row$SDDSRVYR <= 2) "URDUCD" else metal
      concentration <- lab[[field]][j]
    }
    chemicals[[metal]][i] <- concentration
  }
  full_designs[[metal]] <- dplyr::bind_rows(components)
  stopifnot(!anyDuplicated(full_designs[[metal]]$SEQN))
}

# Retain early mercury for descriptions only. Models start in 2003
# because the earlier laboratory sample had different coverage.
for (cycle in 1:2) {
  path <- file.path(raw_directory, paste0(c("LAB06", "L06_B")[cycle], ".xpt"))
  lab <- read_native(path)
  i <- which(chemicals$SDDSRVYR == cycle)
  chemicals$LBXTHG[i] <- lab$LBXTHG[match(chemicals$SEQN[i], lab$SEQN)]
}

# Fit the original creatinine model in its 1999–2018 cohort.
# Restore the 2001–2002 baseline urinary measurement before fitting.
reference <- inherited
reference$cycle_f <- factor(reference$SDDSRVYR)
baseline_file <- file.path(INPUT_DIR, "creatinine/L16_B.xpt")
baseline <- read_native(baseline_file)
i <- which(reference$SDDSRVYR == 2)
j <- match(reference$SEQN[i], baseline$SEQN)
stopifnot(!anyNA(j))
reference$URXUCR[i] <- baseline$URXUCR[j]
later_files <- file.path(INPUT_DIR, "creatinine",
                        c("ALB_CR_I.xpt", "ALB_CR_J.xpt"))
for (k in 1:2) {
  lab <- read_native(later_files[k])
  i <- which(reference$SDDSRVYR == k + 8)
  j <- match(reference$SEQN[i], lab$SEQN)
  stopifnot(!anyNA(j), isTRUE(all.equal(reference$URXUCR[i], lab$URXUCR[j], check.attributes=FALSE)))
  reference$URXUCR[i] <- lab$URXUCR[j]
}
creatinine_model <- stats::lm(log(URXUCR) ~ RIDAGEYR + I(RIDAGEYR^2) +
                              RIAGENDR_f + RIDRETH1_f + BMXBMI + cycle_f,
                            data=reference, na.action=stats::na.exclude)
stopifnot(stats::nobs(creatinine_model) == 28555L)
reference$URXUCR_pred <- as.numeric(exp(stats::predict(creatinine_model, newdata=reference)))
j <- match(chemicals$SEQN, reference$SEQN)
stopifnot(!anyNA(j))
chemicals <- chemicals %>% dplyr::mutate(
  URXUCR=reference$URXUCR[j], URXUCR_pred=reference$URXUCR_pred[j],
  log_URXUCR=log(URXUCR), URXUCR_cs_ratio=URXUCR/URXUCR_pred,
  log_primary_LBXBPB=safe_log2(LBXBPB), log_primary_LBXTHG=safe_log2(LBXTHG))
# Validate blood logs, then retain the supplied numerical representation.
# Machine-precision differences can affect CQR estimates.
blood_log_checks <- list()
for (metal in c("LBXBPB", "LBXTHG")) {
  field <- paste0("log_primary_", metal)
  calculated <- chemicals[[field]]
  archived <- inherited[[field]][j]
  check <- compare_column(calculated, archived, field, tolerance=1e-12)
  check$values_differing_at_machine_precision <- sum(calculated != archived, na.rm=TRUE)
  stopifnot(check$pass)
  chemicals[[paste0(field, "_recalculated")]] <- calculated
  chemicals[[field]] <- archived
  blood_log_checks[[metal]] <- check
}
write.csv(dplyr::bind_rows(blood_log_checks), diagnostic_file("blood_log_numeric_representation.csv"), row.names=FALSE)
chemicals$log_VNSUMARS <- safe_log2(chemicals$VNSUMARS)
chemicals$z_log_VNSUMARS_barr <- (chemicals$log_VNSUMARS - mean(chemicals$log_VNSUMARS, na.rm=TRUE)) /
  stats::sd(chemicals$log_VNSUMARS, na.rm=TRUE)
for (metal in c("URXUCD", "VNSUMARS", "URXUTL")) {
  chemicals[[paste0(metal, "_cs")]] <- chemicals[[metal]] / chemicals$URXUCR_cs_ratio
  chemicals[[paste0("log_", metal, "_cs")]] <- safe_log2(chemicals[[paste0(metal, "_cs")]])
  chemicals[[paste0("log_primary_", metal)]] <- chemicals[[paste0("log_", metal, "_cs")]]
  x <- chemicals[[paste0("log_", metal, "_cs")]]
  chemicals[[paste0("z_log_", metal)]] <- (x - mean(x, na.rm=TRUE)) / stats::sd(x, na.rm=TRUE)
}
saveRDS(chemicals, data_file("chemicals.rds"))
saveRDS(full_designs, data_file("full_designs.rds"))
saveRDS(creatinine_model, data_file("creatinine_model.rds"))
saveRDS(reference %>% dplyr::select(SEQN, SDDSRVYR, URXUCR, URXUCR_pred),
        data_file("creatinine_reference.rds"))
write.csv(specification, diagnostic_file("native_component_weight_specification.csv"), row.names=FALSE)
message("Chemical reconstruction complete: ", nrow(chemicals), " retained archive rows.")
