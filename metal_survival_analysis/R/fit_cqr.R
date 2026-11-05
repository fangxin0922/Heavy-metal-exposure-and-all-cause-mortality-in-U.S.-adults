# Fit one scenario and imputation with quantreg::crq().
# Usage: Rscript metal_survival_analysis/R/fit_cqr.R scenario imputation [400|0]
source("metal_survival_analysis/R/model_catalog.R")
args <- commandArgs(trailingOnly = TRUE)
scenario <- args[1]
imputation <- as.integer(args[2])
resamples <- if (length(args) >= 3L) as.integer(args[3]) else BOOTSTRAPS
stopifnot(scenario %in% SCENARIOS$scenario, imputation %in% seq_len(M),
          resamples %in% c(0L, BOOTSTRAPS))
destination <- file.path(RUN_DIR, "fits", scenario, sprintf("imp%02d", imputation))
if (dir.exists(destination) && length(list.files(destination, all.files = TRUE, no.. = TRUE))) {
  stop("Fit output already exists; use an empty output directory for a new run.")
}
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
status <- function(message) {
  cat(format(Sys.time()), scenario, imputation, message, "\n")
  flush.console()
  writeLines(c(message, as.character(Sys.time())), file.path(destination, "STATUS.txt"))
}
input_dir <- file.path(RUN_DIR, "models", scenario)
specification <- readRDS(file.path(input_dir, "specification.rds"))
imputations <- readRDS(file.path(input_dir, "imputations.rds"))
analysis <- imputations$data[[imputation]]
rm(imputations)
catalog <- SCENARIOS[match(scenario, SCENARIOS$scenario), ]
stopifnot(identical(specification$scenario, scenario))
if (catalog$require_egfr60) {
  stopifnot(isTRUE(specification$require_egfr60),
            length(specification$observed_egfr) == nrow(analysis),
            all(is.finite(specification$observed_egfr)),
            all(specification$observed_egfr >= 60))
  if ("egfr_ckdepi_2021" %in% names(analysis)) {
    stopifnot(identical(analysis$egfr_ckdepi_2021, specification$observed_egfr))
  }
}
weights <- specification$sampling_weights
analysis$.analysis_weight <- as.numeric(weights) / mean(weights)
stopifnot(length(weights) == nrow(analysis),
          identical(as.numeric(specification$row_ids), as.numeric(analysis$SEQN)),
          nrow(analysis) == specification$n, sum(analysis$event) == specification$events,
          all(is.finite(analysis$.analysis_weight)), all(analysis$.analysis_weight > 0))
formula <- specification$formula
environment(formula) <- .GlobalEnv
method <- specification$method
if (is.null(method)) method <- SCENARIOS$method[match(scenario, SCENARIOS$scenario)]

# Guard against native Peng-Huang workspace overflow when n > 9 * events.
# This protects memory; it is not a cohort exclusion rule.
if (method == "PengHuang" && nrow(analysis) > 9 * sum(analysis$event)) {
  stop("Native Peng-Huang workspace risk for this outcome; model not launched")
}
stopifnot(identical(environment(quantreg::crq), asNamespace("quantreg")))
warnings_seen <- character()
capture_warning <- function(w) {
  warnings_seen <<- c(warnings_seen, conditionMessage(w))
  invokeRestart("muffleWarning")
}
point_path <- file.path(destination, "point.rds")
status("fitting point estimates")
started <- proc.time()["elapsed"]
fit_error <- ""
fit <- tryCatch(withCallingHandlers(
  quantreg::crq(formula, data = analysis, weights = .analysis_weight,
                method = method, na.action = na.fail),
  warning = capture_warning), error = function(e) {fit_error <<- conditionMessage(e); NULL})
coefficients <- if (is.null(fit)) NULL else as.matrix(coef(fit, taus = TAUS))
point <- list(fit = fit, coefficients = coefficients, fit_error = fit_error,
              fit_seconds = unname(proc.time()["elapsed"] - started),
              warnings = unique(warnings_seen))
save_result(point, point_path)
targets <- specification$terms
if (is.null(targets)) targets <- specification$term
result_table <- expand.grid(term = targets, tau = TAUS, stringsAsFactors = FALSE) %>%
  mutate(scenario = scenario, imputation = imputation,
         method = method, n = specification$n, deaths = specification$events,
         period = specification$period, estimate = NA_real_, lower = NA_real_,
         upper = NA_real_, std_error = NA_real_, p_value = NA_real_,
         package_missing_draws = NA_integer_, package_Brep = NA_integer_,
         status = "percentile_unavailable")
if (!is.null(coefficients)) {
  for (i in seq_len(nrow(result_table))) {
    result_table$estimate[i] <- coefficients[result_table$term[i], match(result_table$tau[i], TAUS)]
  }
}
result_table$status[is.finite(result_table$estimate)] <- "point_only"
if (nzchar(point$fit_error)) result_table$status <- "package_fit_error"
write_table(result_table, file.path(destination, "point_estimates.csv"))
if (resamples == 0L) {
  status("point estimates complete")
  quit(save = "no", status = if (nzchar(point$fit_error)) 1L else 0L)
}

available <- if (is.null(coefficients)) numeric() else TAUS[apply(coefficients, 2, function(x) any(!is.na(x)))]
# Use fixed seeds for imputation-specific summaries.
seed <- if (scenario %in% PRIMARY) 20261029L + match(scenario, PRIMARY) * 1000L + imputation else
  20261103L + match(scenario, SCENARIOS$scenario) * 1000L + imputation
set.seed(seed)
summary_error <- ""
package_summary <- NULL
started <- proc.time()["elapsed"]
if (length(available)) {
  status("running public summary: 400 jackknife resamples")
  package_summary <- tryCatch(withCallingHandlers(
    summary(fit, taus = available, se = "boot", bmethod = "jack", R = resamples, alpha = .05),
    warning = capture_warning), error = function(e) {summary_error <<- conditionMessage(e); NULL})
}
if (!is.null(package_summary)) {
  # A single-quantile summary duplicates entries; both update the same row.
  for (item in package_summary) for (term in targets) {
    index <- which(abs(result_table$tau - item$tau) < 1e-10 & result_table$term == term)
    row <- item$coefficients[term, ]
    stopifnot(length(index) == 1L,
              isTRUE(all.equal(unname(row["Value"]), result_table$estimate[index], tolerance = 1e-12)))
    result_table[index, c("lower", "upper", "std_error", "p_value")] <-
      as.list(unname(row[c("Lower Bd", "Upper Bd", "Std Error", "Pr(>|t|)")]))
    result_table$package_missing_draws[index] <- item$NAs
    result_table$package_Brep[index] <- item$Brep
    result_table$status[index] <- if (all(is.finite(row[c("Lower Bd", "Upper Bd")])))
      "package_CI_available" else "package_CI_unavailable"
  }
}
if (nzchar(summary_error)) result_table$status[is.finite(result_table$estimate)] <- "package_summary_error"
result <- list(table = result_table, package_summary = package_summary,
  seed = seed, requested_R = resamples, fit_error = point$fit_error,
  summary_error = summary_error, warnings = unique(c(point$warnings, warnings_seen)),
  fit_seconds = point$fit_seconds, summary_seconds = unname(proc.time()["elapsed"] - started))
save_result(result, file.path(destination, "result.rds"))
write_table(result_table, file.path(destination, "results.csv"))
status(if (nzchar(summary_error) || nzchar(point$fit_error)) "finished with package error" else "complete")
if (nzchar(summary_error) || nzchar(point$fit_error)) quit(save = "no", status = 1L)
