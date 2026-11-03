# Prepare complete-case biomarker CQR models and survey replicate weights.

source("replication_through_2014/remote_execution/R/runtime.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survey)
  library(quantreg)
})
source("replication_through_2014/R/revised_cqr_solver.R")

BASE_DIR <- "nutrient_biomarker_emm_run_2026-11-02"
OUTPUT_DIR <- file.path(BASE_DIR, "project_style_v3")
KAPPA <- 2
PARENT_DIR <- "replication_through_2014/remote_execution/models"
SCENARIOS <- c(Pb = "primary_Pb", Hg = "primary_Hg", Cd = "egfr60_Cd", As = "egfr60_As")
BIOMARKERS <- c("vitamin_c", "alpha_carotene", "rbc_folate")
TAUS <- c(.10, .25, .50)
TERMS <- c("metal_log2", "biomarker_iqr", "metal_log2:biomarker_iqr")
dir.create(file.path(OUTPUT_DIR, "models"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "aggregate"), showWarnings = FALSE)

# Fixed reporting levels ----

weighted_quantile <- function(x, w, probs) {
  i <- order(x)
  cumulative <- cumsum(w[i]) / sum(w)
  vapply(probs, function(p) x[i[which(cumulative >= p)[1]]], numeric(1))
}

# Model, sample and survey weights ----

rows <- list()
for (metal in names(SCENARIOS)) {
  parent_path <- file.path(PARENT_DIR, SCENARIOS[[metal]], "replicate_design.rds")
  bank <- readRDS(parent_path)
  stopifnot(bank$R == 400L, identical(bank$mse, FALSE), all(bank$sampling_weights > 0))
  for (biomarker in BIOMARKERS) {
    pair <- paste(metal, biomarker, sep = "_")
    input_path <- file.path(OUTPUT_DIR, "data", paste0(pair, ".rds"))
    input <- readRDS(input_path)
    d <- droplevels(as.data.frame(input$data))
    covariates <- input$covariates
    required <- unique(c("SEQN", "SDDSRVYR", "SDMVSTRA", "SDMVPSU", "followup_yrs",
                         "event", "metal_log2", "biomarker_value", "analysis_weight", covariates))
    stopifnot(all(required %in% names(d)), !anyNA(d[required]), !anyDuplicated(d$SEQN),
              all(is.finite(d$analysis_weight) & d$analysis_weight > 0),
              all(is.finite(d$biomarker_value) & d$biomarker_value >= 0),
              all(is.finite(d$metal_log2)), all(d$event %in% 0:1), all(d$followup_yrs > 0))
    if (biomarker == "vitamin_c") stopifnot(all(d$SDDSRVYR %in% 3:4))
    if (metal %in% c("Cd", "As")) stopifnot(all(d$egfr_ckdepi_2021 >= 60))
    levels <- weighted_quantile(d$biomarker_value, d$analysis_weight, c(.25, .75))
    stopifnot(all(is.finite(levels)), diff(levels) > 0)
    # The product coefficient is the P75-minus-P25 interaction in years.
    d <- d %>% mutate(biomarker_iqr = (biomarker_value - levels[1]) / diff(levels))
    formula <- as.formula(paste("survival::Surv(followup_yrs, event) ~",
      "metal_log2 * biomarker_iqr + RIDAGEYR + I(RIDAGEYR^2) +",
      paste(setdiff(covariates, "RIDAGEYR"), collapse = " + ")))
    prepared <- revised_cqr_prepare(formula, d, taus = TAUS, standardize = TRUE)
    stopifnot(prepared$matrix_rank == prepared$p, all(TERMS %in% prepared$coefficient_names))

    # Preserve full-component PSU draws, then select the biomarker/CC domain.
    # New native weights are needed for alpha-carotene's 2001-start window.
    j <- match(d$SEQN, bank$row_ids)
    stopifnot(!anyNA(j))
    multipliers <- bank$analysis_weights[j, , drop = FALSE] / as.numeric(bank$sampling_weights[j])
    replicate_weights <- multipliers * d$analysis_weight
    stopifnot(all(is.finite(replicate_weights)), all(replicate_weights >= 0),
              all(colSums(replicate_weights) > 0))
    psu_key <- interaction(d$SDMVSTRA, d$SDMVPSU, drop = TRUE)
    first <- !duplicated(psu_key)
    within_psu_error <- max(abs(multipliers - multipliers[match(psu_key, psu_key), , drop = FALSE]))
    stopifnot(within_psu_error < 1e-10)
    rank_df <- qr(t(multipliers[first, , drop = FALSE]), tol = 1e-5)$rank - 1L
    psu_strata_df <- sum(first) - length(unique(d$SDMVSTRA))
    inference_df <- min(rank_df, psu_strata_df)
    stopifnot(inference_df > 0)
    rp <- list(row_ids = d$SEQN, sampling_weights = d$analysis_weight,
      analysis_weights = replicate_weights, scale = bank$scale, rscales = bank$rscales,
      mse = bank$mse, R = bank$R, seed = bank$seed, complete_data_df = inference_df,
      rank_df = rank_df, psu_strata_df = psu_strata_df,
      parent_bank = parent_path, full_component_n = length(bank$full_row_ids))

    # Preserve the documented kappa=2 bootstrap amendment and its arithmetic order.
    # This preserves linear-total variance, not exact nonlinear CQR variance.
    original_weights <- rp$analysis_weights
    original_scale <- rp$scale
    multipliers <- original_weights / as.numeric(rp$sampling_weights)
    positive_multipliers <- 1 + (multipliers - 1) / KAPPA
    rp$analysis_weights <- positive_multipliers * as.numeric(rp$sampling_weights)
    rp$scale <- original_scale * KAPPA^2
    rp$rescale_divisor <- KAPPA
    rp$rescale_method <- "Positive affine rescaling of Rao-Wu bootstrap multipliers; global variance scale multiplied by kappa squared"
    stopifnot(all(rp$analysis_weights > 0),
      qr(t(positive_multipliers[first, , drop = FALSE]), tol = 1e-5)$rank - 1L == rank_df)
    audit_values <- as.matrix(d[c("biomarker_value", "metal_log2", "followup_yrs")])
    old_totals <- t(crossprod(audit_values, original_weights))
    new_totals <- t(crossprod(audit_values, rp$analysis_weights))
    old_variance <- svrVar(old_totals, original_scale, rp$rscales, mse = FALSE, na.action = "na.fail")
    new_variance <- svrVar(new_totals, rp$scale, rp$rscales, mse = FALSE, na.action = "na.fail")
    stopifnot(isTRUE(all.equal(old_variance, new_variance, tolerance = 1e-10,
      check.attributes = FALSE)))
    replicate_weights <- rp$analysis_weights

    # A simple statistic verifies the weight/variance chain independently.
    rep_design <- svrepdesign(data = d, weights = ~analysis_weight,
      repweights = replicate_weights, type = "other", combined.weights = TRUE,
      scale = rp$scale, rscales = rp$rscales, mse = FALSE, degf = inference_df)
    means <- colSums(replicate_weights * d$biomarker_value) / colSums(replicate_weights)
    direct_variance <- as.numeric(svrVar(means, rp$scale, rp$rscales, mse = FALSE, na.action = "na.fail"))
    check <- withReplicates(rep_design, function(w, data) weighted.mean(data$biomarker_value, w))
    stopifnot(isTRUE(all.equal(direct_variance, as.numeric(vcov(check)), tolerance = 1e-10)))
    dest <- file.path(OUTPUT_DIR, "models", pair)
    dir.create(dest, showWarnings = FALSE)
    target <- file.path(dest, "prepared.rds")
    stopifnot(!file.exists(target))
    spec <- list(pair = pair, metal = metal, biomarker = biomarker,
      units = input$units, n = nrow(d), deaths = sum(d$event), covariates = covariates,
      formula = formula, prepared = prepared, replicate_design = rp, terms = TERMS,
      biomarker_low = levels[1], biomarker_high = levels[2], taus = TAUS,
      primary_tau = .10, weight_constant = mean(d$analysis_weight),
      source_hashes = tools::md5sum(c(input_path, parent_path,
        file.path(OUTPUT_DIR, "R", "02_prepare_cqr_models.R"),
        "replication_through_2014/R/revised_cqr_solver.R")),
      inference_amendment = "Fixed kappa=2 positive rescaled Rao-Wu bootstrap, adopted after numerical failure diagnostics; original model unchanged",
      missing_data = "Complete cases; no imputation; separate sample for each biomarker and metal",
      contrast = "Fixed weighted P25 to P75 from the complete-case sample",
      df_rule = "Minimum of domain replicate-weight rank minus one and represented PSUs minus strata")
    saveRDS(spec, target)
    rows[[length(rows) + 1L]] <- data.frame(pair, metal, biomarker, n = spec$n, deaths = spec$deaths,
      units = input$units, biomarker_low = levels[1], biomarker_high = levels[2],
      cycles = paste(sort(unique(d$SDDSRVYR)), collapse = ";"),
      p = prepared$p, matrix_rank = prepared$matrix_rank, condition_number = prepared$condition_number,
      strata = length(unique(d$SDMVSTRA)), psus = sum(first), rank_df, psu_strata_df, inference_df,
      rescale_divisor = KAPPA, R = rp$R, multiplier_within_psu_error = within_psu_error,
      replicate_variance_check = direct_variance, survey_variance_check = as.numeric(vcov(check)))
  }
}
write.csv(do.call(rbind, rows), file.path(OUTPUT_DIR, "aggregate", "model_preparation.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "aggregate", "preparation_session.txt"))
