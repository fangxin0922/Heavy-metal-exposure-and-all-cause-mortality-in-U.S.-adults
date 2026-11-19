# Reconstruct manuscript spline predictions from newly saved fits; never refit.
source("metal_survival_analysis/R/model_catalog.R")
OUT <- file.path(RUN_DIR, "results", "spline_predictions")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
predictions <- profiles <- knots_out <- rugs <- checks <- coefficients_out <- list()
source_paths <- character()
for (scenario in c("rcs_Pb", "rcs_egfr60_Cd")) {
  model_dir <- file.path(RUN_DIR, "models", scenario)
  spec_path <- file.path(model_dir, "specification.rds")
  frame_path <- file.path(model_dir, "preimputation_frame.rds")
  spec <- readRDS(spec_path)
  frame <- readRDS(frame_path)
  stopifnot(identical(as.numeric(frame$SEQN), as.numeric(spec$row_ids)),
            nrow(frame) == spec$n, sum(frame$event) == spec$events,
            identical(spec$method, "PengHuang"), identical(spec$terms, c("rcs1", "rcs2", "rcs3")))
  if (scenario == "rcs_egfr60_Cd") {
    primary <- readRDS(file.path(RUN_DIR, "models/egfr60_Cd/specification.rds"))
    stopifnot(spec$n == 7061L, spec$events == 1510L,
      all(is.finite(frame$egfr_ckdepi_2021)), all(frame$egfr_ckdepi_2021 >= 60),
      identical(as.numeric(frame$SEQN), as.numeric(primary$row_ids)),
      identical(as.numeric(spec$sampling_weights), as.numeric(primary$sampling_weights)),
      "egfr_ckdepi_2021" %in% spec$covariates,
      identical(spec$knots, unname(quantile(frame$log_primary_URXUCD, c(.05,.35,.65,.95), type=7))))
  }
  exposure_term <- paste0("log_primary_", spec$metal)
  profile <- frame[1L, spec$covariates, drop = FALSE]
  for (variable in spec$covariates) {
    values <- frame[[variable]]
    if (is.factor(values)) {
      # Ties follow the saved factor-level order; no integer median of factors.
      counts <- table(values)
      selected <- names(counts)[which.max(counts)]
      profile[[variable]] <- factor(selected, levels = levels(values), ordered = is.ordered(values))
      rule <- "Unweighted observed modal category; ties use saved factor-level order"
    } else {
      selected <- median(values, na.rm = TRUE)
      stopifnot(is.finite(selected))
      profile[[variable]] <- selected
      rule <- "Unweighted median of observed preimputation covariate values"
    }
    profiles[[length(profiles) + 1L]] <- data.frame(scenario, variable,
      value = as.character(selected), rule, observed_n = sum(!is.na(values)))
  }
  concentration_knots <- 2^spec$knots
  concentration <- sort(unique(c(seq(concentration_knots[1L], concentration_knots[4L],
                                     length.out = 301L), concentration_knots)))
  prediction_data <- profile[rep(1L, length(concentration)), , drop = FALSE]
  prediction_data[[exposure_term]] <- log2(concentration)
  prediction_data[, spec$terms] <- restricted_cubic_basis(log2(concentration), spec$knots)
  X <- model.matrix(delete.response(terms(spec$formula)), prediction_data)
  observed <- 2^frame[[exposure_term]]
  stopifnot(all(is.finite(observed)), all(observed > 0))
  h <- hist(observed[observed >= min(concentration) & observed <= max(concentration)],
    breaks = seq(min(concentration), max(concentration), length.out = 201L), plot = FALSE)
  rugs[[scenario]] <- data.frame(scenario, concentration = h$mids, count = h$counts)
  knots_out[[scenario]] <- data.frame(scenario, percentile = c(.05, .35, .65, .95),
    log2_exposure = spec$knots, concentration = concentration_knots, n = spec$n, deaths = spec$events)
  source_paths <- c(source_paths, spec_path, frame_path, file.path(model_dir, "imputations.rds"))
  for (imputation in seq_len(M)) {
    fit_dir <- file.path(RUN_DIR, "fits", scenario, sprintf("imp%02d", imputation))
    point_path <- file.path(fit_dir, "point.rds")
    result_path <- file.path(fit_dir, "result.rds")
    point <- readRDS(point_path)
    result <- if (file.exists(result_path)) readRDS(result_path) else NULL
    stopifnot(!nzchar(point$fit_error),
      identical(point$coefficients, as.matrix(coef(point$fit, taus = TAUS))),
      setequal(colnames(X), rownames(point$coefficients)))
    if (!is.null(result)) stopifnot(!nzchar(result$summary_error), result$requested_R == 400L)
    Z <- X[, rownames(point$coefficients), drop = FALSE]
    for (j in seq_along(TAUS)) {
      tau <- TAUS[j]
      beta <- point$coefficients[, j]
      valid <- all(is.finite(beta))
      stopifnot(valid || all(!is.finite(beta)))
      estimate <- se <- lower <- upper <- rep(NA_real_, nrow(Z))
      missing_draws <- NA_integer_
      covariance_symmetric <- covariance_positive <- NA
      prediction_arithmetic_error <- NA_real_
      if (valid) estimate <- drop(Z %*% beta)
      if (valid && !is.null(result)) {
        matched <- Filter(function(x) abs(x$tau - tau) < 1e-10, result$package_summary)
        stopifnot(length(matched) >= 1L)
        summary <- matched[[1L]]
        stopifnot(identical(rownames(summary$coefficients), names(beta)),
                  isTRUE(all.equal(unname(beta), unname(summary$coefficients[, "Value"]), tolerance = 0)))
        V <- summary$cov
        stopifnot(identical(dim(V), c(length(beta), length(beta))), all(is.finite(V)))
        covariance_symmetric <- max(abs(V - t(V))) < 1e-10
        covariance_positive <- min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) >= -1e-10
        stopifnot(covariance_symmetric, covariance_positive)
        estimate <- drop(Z %*% beta)
        variance <- rowSums((Z %*% V) * Z)
        stopifnot(all(variance >= -1e-10))
        se <- sqrt(pmax(variance, 0))
        lower <- estimate - qnorm(.975) * se
        upper <- estimate + qnorm(.975) * se
        # Independent scalar check at both endpoints and a midpoint.
        sample_rows <- c(1L, ceiling(nrow(Z) / 2L), nrow(Z))
        scalar <- vapply(sample_rows, function(k) sum(Z[k, ] * beta), numeric(1))
        prediction_arithmetic_error <- max(abs(scalar - estimate[sample_rows]))
        stopifnot(prediction_arithmetic_error < 1e-10,
          all(is.finite(c(estimate, se, lower, upper))),
          all(lower <= estimate & estimate <= upper))
        missing_draws <- summary$NAs
      }
      predictions[[length(predictions) + 1L]] <- data.frame(scenario, imputation, tau,
        concentration, log2_exposure = log2(concentration), estimate, std_error = se,
        lower, upper, n = spec$n, deaths = spec$events,
        package_missing_draws = missing_draws, package_attempts = if (valid && !is.null(result)) 400L else NA_integer_,
        status = if (valid) "available" else "percentile_unavailable",
        interval_method = "Approximate pointwise normal interval: X beta +/- 1.96 sqrt(X V X'); saved package bootstrap covariance; no MI pooling")
      checks[[length(checks) + 1L]] <- data.frame(scenario, imputation, tau,
        available = valid, grid_rows = nrow(Z), covariance_symmetric, covariance_positive,
        prediction_arithmetic_error, package_missing_draws = missing_draws)
      coefficients_out[[length(coefficients_out) + 1L]] <- data.frame(scenario, imputation, tau,
        term = names(beta), estimate = as.numeric(beta))
    }
    source_paths <- c(source_paths, point_path, if (file.exists(result_path)) result_path)
    rm(point, result)
  }
}
predictions <- bind_rows(predictions)
checks <- bind_rows(checks)
stopifnot(nrow(checks) == 30L,
  !anyDuplicated(predictions[, c("scenario", "imputation", "tau", "concentration")]))
write_table(predictions, file.path(OUT, "figure2_predictions_by_imputation.csv"))
write_table(bind_rows(profiles), file.path(OUT, "figure2_reference_profiles.csv"))
write_table(bind_rows(knots_out), file.path(OUT, "figure2_knots_and_cohorts.csv"))
write_table(bind_rows(rugs), file.path(OUT, "figure2_binned_rug.csv"))
write_table(bind_rows(coefficients_out), file.path(OUT, "figure2_full_coefficients.csv"))
write_table(checks, file.path(OUT, "figure2_validation.csv"))
source_paths <- unique(c(source_paths, "metal_survival_analysis/R/model_catalog.R",
  "metal_survival_analysis/R/spline_predictions.R"))
write_table(data.frame(file = source_paths, md5 = unname(tools::md5sum(source_paths))),
  file.path(OUT, "figure2_source_hashes.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "extraction_session_info.txt"))
cat("Validated", sum(checks$available), "available and", sum(!checks$available),
    "unavailable imputation/quantile curves.\n")
source("metal_survival_analysis/R/plot_cqr_splines.R")
