# Supplementary survey Cox models estimate hazard ratios.
# Run after impute_covariates.R has prepared all primary models.
source("metal_survival_analysis/R/model_catalog.R")

full_designs <- readRDS(file.path(RUN_DIR, "data", "full_designs.rds"))
descriptive_rows <- list()
cox_rows <- list()
spline_rows <- list()
spline_tests <- list()
diagnostic_rows <- list()
result_index <- 0L

capture_analysis <- function(expression) {
  recorded_warnings <- character()
  answer <- tryCatch(withCallingHandlers(
    expression,
    warning = function(w) {
      recorded_warnings <<- c(recorded_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  list(value = answer, warnings = paste(unique(recorded_warnings), collapse = " | "),
       error = if (inherits(answer, "error")) conditionMessage(answer) else "")
}

make_completed_design <- function(full, completed, specification) {
  stopifnot(!anyDuplicated(full$SEQN), !anyDuplicated(completed$SEQN),
            all(completed$SEQN %in% full$SEQN),
            setequal(completed$SEQN, specification$row_ids),
            all(is.finite(full$w_metal)), all(full$w_metal >= 0))
  # Build the full PSU/stratum design before selecting the analysis cohort.
  design_fields <- c("SEQN", "SDDSRVYR", "SDMVSTRA", "SDMVPSU", "w_metal")
  full <- full[, design_fields]
  overlap <- intersect(setdiff(design_fields, "SEQN"), names(completed))
  matched <- match(completed$SEQN, full$SEQN)
  stopifnot(all(full$w_metal[matched] > 0))
  for (variable in overlap) {
    # Use original design weights; fitting weights may have been rescaled.
    if (variable != "w_metal") {
      stopifnot(isTRUE(all.equal(as.vector(completed[[variable]]),
                                as.vector(full[[variable]][matched]),
                                check.attributes = FALSE)))
    }
  }
  completed$.cohort <- TRUE
  joined <- dplyr::left_join(full,
                            completed[, setdiff(names(completed), overlap)],
                            by = "SEQN")
  joined$.cohort[is.na(joined$.cohort)] <- FALSE
  design <- survey::svydesign(data = joined, ids = ~SDMVPSU,
                             strata = ~SDMVSTRA, weights = ~w_metal,
                             nest = TRUE)
  design <- subset(design, .cohort)
  stopifnot(nrow(design$variables) == nrow(completed))
  design
}

summarise_variable <- function(design, variable, scenario, imputation) {
  values <- design$variables[[variable]]
  categories <- if (is.factor(values)) levels(values) else if (is.character(values)) {
    sort(unique(values[!is.na(values)]))
  } else NULL
  n_missing <- sum(is.na(values))
  if (is.null(categories)) {
    fit <- survey::svymean(reformulate(variable), design, na.rm = TRUE)
    interval <- confint(fit)
    return(data.frame(scenario, imputation, variable, level = "",
                      statistic = "survey_mean", unweighted_n = sum(!is.na(values)),
                      missing_n = n_missing, estimate = unname(coef(fit)[1]),
                      std_error = unname(survey::SE(fit)[1]),
                      lower = interval[1, 1], upper = interval[1, 2]))
  }
  bind_rows(lapply(categories, function(level) {
    design$variables$.summary_indicator <- as.numeric(values == level)
    fit <- survey::svymean(~.summary_indicator, design, na.rm = TRUE)
    interval <- confint(fit)
    data.frame(scenario, imputation, variable, level, statistic = "survey_proportion",
               unweighted_n = sum(values == level, na.rm = TRUE),
               missing_n = n_missing, estimate = unname(coef(fit)[1]),
               std_error = unname(survey::SE(fit)[1]),
               lower = interval[1, 1], upper = interval[1, 2])
  }))
}

for (scenario in PRIMARY) {
  model_dir <- file.path(RUN_DIR, "models", scenario)
  specification <- readRDS(file.path(model_dir, "specification.rds"))
  imputations <- readRDS(file.path(model_dir, "imputations.rds"))$data
  stopifnot(length(imputations) == M, length(specification$term) == 1L)
  full <- full_designs[[specification$design_metal]]
  stopifnot(!is.null(full))
  for (imputation in seq_along(imputations)) {
    completed <- imputations[[imputation]]
    term <- specification$term
    stopifnot(all(c("SEQN", term, "followup_yrs", "event") %in% names(completed)),
              all(all.vars(specification$formula) %in% names(completed)),
              all(complete.cases(completed[, all.vars(specification$formula)])))
    design <- make_completed_design(full, completed, specification)
    result_index <- result_index + 1L

    variables <- unique(c(term, specification$covariates, "followup_yrs", "event"))
    descriptive_rows[[result_index]] <- bind_rows(lapply(variables, function(variable) {
      summarise_variable(design, variable, scenario, imputation)
    }))

    # Fit survey Cox models with the same covariates and exposure scale.
    result <- capture_analysis(survey::svycoxph(specification$formula, design = design))
    output <- data.frame(scenario, imputation, term,
                         n = nrow(completed), deaths = sum(completed$event),
                         log_hr = NA_real_, std_error = NA_real_, HR = NA_real_,
                         lower = NA_real_, upper = NA_real_, p_value = NA_real_,
                         status = "error", warnings = result$warnings, error = result$error,
                         contrast = paste("One-unit increase in", term),
                         inference = "NHANES PSU/stratum survey Cox; imputation-specific")
    if (!nzchar(result$error)) {
      fit <- result$value
      estimate <- unname(coef(fit)[term])
      standard_error <- sqrt(unname(vcov(fit)[term, term]))
      output$log_hr <- estimate
      output$std_error <- standard_error
      output$HR <- exp(estimate)
      output$lower <- exp(estimate - 1.96 * standard_error)
      output$upper <- exp(estimate + 1.96 * standard_error)
      output$p_value <- 2 * pnorm(-abs(estimate / standard_error))
      output$status <- if (all(is.finite(c(estimate, standard_error)))) "ok" else "nonfinite"
      save_result(list(formula = specification$formula, coefficients = coef(fit),
                       covariance = vcov(fit), summary = summary(fit),
                       warnings = result$warnings),
                  file.path(model_dir, paste0("cox_imputation_", imputation, ".rds")))
    }
    cox_rows[[result_index]] <- output

    # Explore exposure shape with unweighted Cox splines.
    # These curves do not provide survey-adjusted or CQR inference.
    labels <- attr(terms(specification$formula), "term.labels")
    spline_term <- paste0("pspline(", term, ")")
    spline_formula <- as.formula(paste("Surv(followup_yrs, event) ~",
                                       paste(c(spline_term, setdiff(labels, term)), collapse = " + ")))
    spline_result <- capture_analysis(survival::coxph(spline_formula, data = completed))
    diagnostic_rows[[result_index]] <- data.frame(
      scenario, imputation, method = "unweighted_cox_pspline",
      status = if (nzchar(spline_result$error)) "error" else "ok",
      warnings = spline_result$warnings, error = spline_result$error)
    if (!nzchar(spline_result$error)) {
      spline_fit <- spline_result$value
      curve_result <- capture_analysis(termplot(spline_fit, se = TRUE, plot = FALSE,
                                                terms = spline_term))
      if (nzchar(curve_result$error)) {
        diagnostic_rows[[result_index]]$status <- "curve_error"
        diagnostic_rows[[result_index]]$error <- curve_result$error
      } else {
        curve <- curve_result$value[[1]]
        spline_rows[[result_index]] <- data.frame(
          scenario, imputation, term, exposure = curve$x,
          log_hazard_contribution = curve$y, std_error = curve$se,
          lower = curve$y - 1.96 * curve$se, upper = curve$y + 1.96 * curve$se,
          method = "Unweighted Cox pspline; centered partial log-hazard contribution")
      }
      coefficients <- as.data.frame(summary(spline_fit)$coefficients)
      coefficients$coefficient <- rownames(coefficients)
      rownames(coefficients) <- NULL
      spline_tests[[result_index]] <- cbind(data.frame(scenario, imputation), coefficients)
    }
    message("Supplementary analyses finished: ", scenario, ", imputation ", imputation)
  }
}

write_table(bind_rows(descriptive_rows), file.path(RUN_DIR, "results", "descriptive.csv"))
write_table(bind_rows(cox_rows), file.path(RUN_DIR, "results", "cox_by_imputation.csv"))
write_table(bind_rows(spline_rows), file.path(RUN_DIR, "results", "cox_pspline_curves.csv"))
write_table(bind_rows(spline_tests), file.path(RUN_DIR, "results", "cox_pspline_tests.csv"))
write_table(bind_rows(diagnostic_rows), file.path(RUN_DIR, "results", "cox_pspline_status.csv"))
writeLines(c(
  "Descriptives use survey::svymean for means/proportions and design-based standard errors.",
  "These are explicit survey equivalents of the reference's TableOne summaries; tableone is not required.",
  "unweighted_n is a raw observed count; estimates are weighted means or proportions, not weighted counts.",
  "All five imputations are reported separately. No MI pooling is performed.",
  "The full design is constructed before selecting each model cohort.",
  "Cox HRs are supplemental estimates per one unit of the saved exposure term; CQR estimates remain in years.",
  "Cox HR confidence intervals follow the reference's coefficient +/- 1.96*SE calculation.",
  "Spline diagnostics follow the reference's unweighted coxph(pspline()) method and are exploratory.",
  "Spline curve bands describe centered term contributions; they are not CQR, survey, or MI-pooled intervals.",
  "Quantile g-computation is not run: it changes the estimand and needs a joint-cohort weighting specification."
), file.path(RUN_DIR, "results", "DESCRIPTIVE_COX_METHODS.txt"))
