# Summarize the saved CQR fits using survey replicate covariance.
# This script reads completed checkpoints; it does not fit or replace models.

source("replication_through_2014/remote_execution/R/runtime.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(survey)
})

BASE_DIR <- "nutrient_biomarker_emm_run_2026-11-02"
FIT_DIR <- file.path(BASE_DIR, "project_style_v3")
OUTPUT_DIR <- file.path(FIT_DIR, "full_run_2026-11-03")
N_REPLICATES <- 400L
PRIMARY_TAU <- .10
N_PRIMARY_TESTS <- 12L
TERMS <- c("metal_log2", "biomarker_iqr", "metal_log2:biomarker_iqr")
dir.create(file.path(OUTPUT_DIR, "aggregate"), showWarnings = FALSE)

# Contrasts on the existing survival-time scale ----

biomarker_contrasts <- function(coefficients) {
  c(metal_low = unname(coefficients["metal_log2"]),
    metal_high = unname(coefficients["metal_log2"] + coefficients["metal_log2:biomarker_iqr"]),
    ic = unname(coefficients["metal_log2:biomarker_iqr"]))
}

preparation <- read.csv(file.path(FIT_DIR, "aggregate", "model_preparation.csv"))
results <- diagnostics <- precision <- provenance <- list()

for (pair in preparation$pair) {
  model_dir <- file.path(FIT_DIR, "models", pair)
  model_file <- file.path(model_dir, "prepared.rds")
  spec <- readRDS(model_file)
  design <- spec$replicate_design
  stopifnot(design$R == N_REPLICATES, identical(design$mse, FALSE),
    length(design$rscales) == N_REPLICATES, all(is.finite(design$rscales)),
    all(design$rscales >= 0), is.finite(design$scale), design$scale > 0,
    design$rescale_divisor == 2, design$complete_data_df > 0,
    identical(spec$terms, TERMS), identical(spec$taus, c(.10, .25, .50)))
  stopifnot(identical(unname(tools::md5sum(names(spec$source_hashes))),
    unname(spec$source_hashes)))

  fit_files <- file.path(model_dir, "fits", sprintf("rep%03d.rds", 0:N_REPLICATES))
  stopifnot(all(file.exists(fit_files)))
  fits <- lapply(fit_files, readRDS)
  point <- fits[[1]]
  point_hash <- tools::md5sum(fit_files[1])
  stopifnot(point$replicate == 0L, point$pair == pair,
    identical(unname(tools::md5sum(names(point$hashes))), unname(point$hashes)))
  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    stopifnot(fit$replicate == i - 1L, fit$pair == pair,
      identical(rownames(fit$coefficients), TERMS),
      identical(colnames(fit$coefficients), sprintf("tau=%.2f", spec$taus)),
      identical(fit$hashes, point$hashes), identical(fit$runtime, point$runtime))
    if (i > 1L) stopifnot(identical(fit$original_hash, point_hash))
    dg <- fit$base$diagnostics
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      pair, replicate = fit$replicate, base_failure = dg$failure_class,
      base_finite_targets = sum(is.finite(fit$base$coefficients)),
      effective_finite_targets = sum(is.finite(fit$coefficients)),
      retries = length(fit$attempts))
  }
  fit_hashes <- tools::md5sum(c(model_file, fit_files))
  provenance[[pair]] <- data.frame(file = names(fit_hashes), md5 = unname(fit_hashes))

  for (j in seq_along(spec$taus)) {
    estimate <- biomarker_contrasts(point$coefficients[, j])
    draws <- t(vapply(fits[-1], function(fit) {
      biomarker_contrasts(fit$coefficients[, j])
    }, numeric(3)))
    colnames(draws) <- names(estimate)

    # Each contrast requires its own finite point and all 400 finite draws.
    # An available interaction does not imply available simple-slope intervals.
    finite_n <- colSums(is.finite(draws))
    available <- is.finite(estimate) & finite_n == N_REPLICATES
    covariance <- matrix(NA_real_, 3, 3,
      dimnames = list(names(estimate), names(estimate)))
    if (any(available)) {
      covariance[available, available] <- survey::svrVar(
        draws[, available, drop = FALSE], scale = design$scale,
        rscales = design$rscales, mse = FALSE, coef = estimate[available],
        na.action = "na.fail")
    }
    se <- sqrt(diag(covariance))
    df <- design$complete_data_df
    lower <- estimate - qt(.975, df) * se
    upper <- estimate + qt(.975, df) * se
    p_value <- if (is.finite(se["ic"]) && se["ic"] > 0) {
      2 * pt(-abs(estimate["ic"] / se["ic"]), df)
    } else NA_real_
    status <- case_when(
      !is.finite(estimate["ic"]) ~ "original_quantile_unavailable",
      finite_n["ic"] < N_REPLICATES ~ "incomplete_valid_replicates",
      !is.finite(p_value) ~ "variance_unavailable",
      TRUE ~ "complete")

    results[[length(results) + 1L]] <- data.frame(
      pair, metal = spec$metal, biomarker = spec$biomarker, tau = spec$taus[j],
      n = spec$n, deaths = spec$deaths, units = spec$units,
      biomarker_low = spec$biomarker_low, biomarker_high = spec$biomarker_high,
      metal_low_years = estimate["metal_low"],
      metal_low_lower = lower["metal_low"], metal_low_upper = upper["metal_low"],
      metal_high_years = estimate["metal_high"],
      metal_high_lower = lower["metal_high"], metal_high_upper = upper["metal_high"],
      ic_years = estimate["ic"], ic_se = se["ic"],
      ic_lower_years = lower["ic"], ic_upper_years = upper["ic"],
      p_value, df, replicates_planned = N_REPLICATES,
      valid_replicates = finite_n["ic"],
      low_valid_replicates = finite_n["metal_low"],
      high_valid_replicates = finite_n["metal_high"], status,
      row.names = NULL)
    for (b in c(100L, 200L, 400L)) {
      values <- draws[seq_len(b), "ic"]
      precision[[length(precision) + 1L]] <- data.frame(
        pair, tau = spec$taus[j], replicates = b, valid = sum(is.finite(values)),
        ic_se = if (all(is.finite(values))) design$rescale_divisor * sd(values) else NA_real_)
    }
  }
}

# Planned multiplicity family and numerical qualifications ----

results <- bind_rows(results) %>%
  mutate(ic_days = ic_years * 365.25,
    ic_lower_days = ic_lower_years * 365.25,
    ic_upper_days = ic_upper_years * 365.25,
    q_primary = NA_real_)
primary_rows <- which(results$tau == PRIMARY_TAU)
stopifnot(length(primary_rows) == N_PRIMARY_TESTS,
  !anyDuplicated(results$pair[primary_rows]))
tested <- primary_rows[is.finite(results$p_value[primary_rows])]
results$q_primary[tested] <- p.adjust(results$p_value[tested], "BH", n = N_PRIMARY_TESTS)

path_file <- file.path(BASE_DIR, "solver_diagnostic_v2", "aggregate", "path_equivalence.csv")
path_checks <- read.csv(path_file)
# The probe examined one successful replicate per pair, not every saved fit.
control_checks <- path_checks %>% filter(kind == "successful_control") %>%
  group_by(pair) %>%
  summarise(successful_control_matches = sum(control_equivalence %in% c("True", "TRUE")),
    successful_control_comparisons = n(), .groups = "drop")
stopifnot(nrow(control_checks) == N_PRIMARY_TESTS,
  all(control_checks$successful_control_comparisons == 2L))
results <- results %>% left_join(control_checks, by = "pair") %>%
  mutate(full_solver_equivalence_validated = FALSE,
    reporting_status = if_else(status == "complete", "provisional_exploratory", "unavailable"),
    successful_control_sensitivity = successful_control_matches < successful_control_comparisons)

write.csv(results, file.path(OUTPUT_DIR, "aggregate", "results.csv"), row.names = FALSE)
write.csv(bind_rows(diagnostics), file.path(OUTPUT_DIR, "aggregate", "fit_diagnostics.csv"), row.names = FALSE)
write.csv(bind_rows(precision), file.path(OUTPUT_DIR, "aggregate", "bootstrap_precision.csv"), row.names = FALSE)
write.csv(bind_rows(provenance), file.path(OUTPUT_DIR, "aggregate", "checkpoint_provenance.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "aggregate", "summary_session.txt"))

# Readable primary table ----

labels <- c(vitamin_c = "Serum vitamin C", alpha_carotene = "Serum alpha-carotene",
  rbc_folate = "RBC folate")
primary <- results %>% filter(tau == PRIMARY_TAU)
number <- function(x, digits = 3) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}
lines <- c("# Biomarker effect modification: revised scripts and audited reporting", "",
  "The numbered R scripts follow the reference project's practical workflow and use compatible packages. These results were calculated from the newly executed revised CQR fitting script, including all 12 original fits and 4,800 fixed survey replicate fits. The preceding numerical limitations remain relevant.", "",
  sprintf("%d/12 primary interactions have provisional intervals; %d meet BH q<0.05. %d tests remain unavailable. No alternative-solver coefficient was substituted.",
    sum(primary$status == "complete"), sum(primary$q_primary < .05, na.rm = TRUE),
    sum(primary$status != "complete")), "",
  sprintf("All available intervals remain provisional and conditional on the native implementation. The earlier alternative-solver gate passed only %d/%d successful-control comparisons. Controls for arsenic–vitamin C and arsenic–RBC folate also disagreed, so the concern extends to complete analyses. These results are not ready for definitive manuscript inference.",
    sum(control_checks$successful_control_matches),
    sum(control_checks$successful_control_comparisons)), "",
  "| Metal | Biomarker | N / deaths | Interaction years (95% CI) | P | BH q | Valid draws | Reporting status |",
  "| --- | --- | ---: | --- | ---: | ---: | ---: | --- |")
for (i in seq_len(nrow(primary))) {
  z <- primary[i, ]
  lines <- c(lines, sprintf("| %s | %s | %s / %s | %s (%s, %s) | %s | %s | %s/400 | %s |",
    z$metal, labels[z$biomarker], z$n, z$deaths, number(z$ic_years),
    number(z$ic_lower_years), number(z$ic_upper_years), number(z$p_value, 4),
    number(z$q_primary, 4), z$valid_replicates, z$reporting_status))
}
lines <- c(lines, "",
  "The contrast is the per-doubling metal association at biomarker P75 minus that at P25 at the 10th conditional survival-time percentile. A positive contrast suggests attenuation only if the low-status metal association is adverse. This is not RERI, an individual gain in lifespan, or evidence of supplement protection.", "",
  "Each metal–biomarker pair uses its own observed complete-case sample and measured cycles. Serum vitamin C is restricted to 2003–2006. The calibrated continuous biomarkers, model covariates, native weights, fixed quartiles, grid, checked native solver and saved 400 survey draws are retained. Variance uses survey::svrVar with the documented kappa=2 positive-rescaling correction, replicate-mean centering and stored survey degrees of freedom. Each contrast independently requires a finite point and 400/400 finite draws. BH uses the planned 12 primary tests; secondary quantiles remain descriptive.", "",
  "Complete-case selection, residual confounding, simultaneous exposure/biomarker measurement, assay conversion uncertainty, a linear interaction specification, censoring assumptions and numerical path sensitivity remain limitations. Missing intervals are not null findings. Following dietary analyses makes this an exploratory extension.", "",
  "Outputs: [all results](aggregate/results.csv), [R figure](figures/primary_interaction.pdf), [workflow and package choices](README.md), [numerical diagnosis](../../solver_diagnostic_v2/RESULTS.md), and [independent numerical review](../../solver_diagnostic_v2/audit/INDEPENDENT_REVIEW.md).")
writeLines(lines, file.path(OUTPUT_DIR, "RESULTS.md"))
