# Compare current fits with frozen historical results.
source("metal_survival_analysis/R/model_catalog.R")
read_fit_tables <- function(branch) {
  paths <- list.files(file.path(RUN_DIR, "fits", branch), pattern = "(point_estimates[.]csv|result[.]rds)$",
                      recursive = TRUE, full.names = TRUE)
  bind_rows(lapply(unique(dirname(paths)), function(directory) {
    completed <- file.path(directory, "result.rds")
    table <- if (file.exists(completed)) readRDS(completed)$table else read.csv(file.path(directory, "point_estimates.csv"))
    table$summary_attempted <- file.exists(completed)
    table
  }))
}
previous_path <- file.path(PROJECT_ROOT,
  "replication_through_2014/official_crq_builtin_20261029/results_by_imputation.csv")
previous <- read.csv(previous_path)
previous <- previous %>% select(scenario, imputation, term, tau, n, deaths,
  estimate, lower, upper, std_error, p_value, status) %>%
  rename_with(~paste0("previous_", .x), -all_of(c("scenario", "imputation", "term", "tau")))
same_number <- function(x, y, tolerance = 1e-8) {
  (is.na(x) & is.na(y)) | (is.finite(x) & is.finite(y) & abs(x - y) <= tolerance)
}
fresh <- read_fit_tables("fresh")
control <- read_fit_tables("control")
compare_primary <- function(current, name) {
  if (!nrow(current)) return(data.frame())
  comparison <- previous %>% left_join(current, by = c("scenario", "imputation", "term", "tau")) %>%
    mutate(executed = !is.na(branch),
           estimate_difference = estimate - previous_estimate,
           fit_success = executed & !status %in% c("package_fit_error"),
           same_availability = fit_success & (is.finite(estimate) == is.finite(previous_estimate)),
           same_estimate = fit_success & same_number(estimate, previous_estimate),
           same_lower = executed & !is.na(package_Brep) & same_number(lower, previous_lower),
           same_upper = executed & !is.na(package_Brep) & same_number(upper, previous_upper))
  write_table(comparison, file.path(RUN_DIR, "results", paste0(name, ".csv")))
  comparison
}
control_comparison <- compare_primary(control, "control_vs_previous")
fresh_comparison <- compare_primary(fresh, "fresh_vs_previous")
if (nrow(fresh)) write_table(fresh, file.path(RUN_DIR, "results/results_by_imputation.csv"))
if (nrow(control)) write_table(control, file.path(RUN_DIR, "results/control_by_imputation.csv"))

# Historical panel intervals used survey bootstrap and MI pooling.
# Compare points only: their intervals use different methods.
if (nrow(fresh)) {
  historical <- read.csv(file.path(PROJECT_ROOT, "replication_through_2014/remote_execution/results/panel_results.csv"))
  means <- fresh %>% group_by(scenario, term, tau) %>%
    summarise(imputations = n_distinct(imputation), finite = sum(is.finite(estimate)),
      fresh_mean = if (n_distinct(imputation) == M && all(is.finite(estimate))) mean(estimate) else NA_real_,
      .groups = "drop")
  manuscript <- historical %>% select(scenario, term, tau, previous_pooled_estimate = estimate,
    previous_status = status, previous_ci_lower = ci_lower, previous_ci_upper = ci_upper) %>%
    left_join(means, by = c("scenario", "term", "tau")) %>%
    mutate(difference = fresh_mean - previous_pooled_estimate,
      comparison = "Point comparison only: historical solver/grid and inference differ")
  write_table(manuscript, file.path(RUN_DIR, "results/manuscript_point_comparison.csv"))
}
expected <- expand.grid(scenario = SCENARIOS$scenario, imputation = seq_len(M), stringsAsFactors = FALSE)
requested_path <- file.path(RUN_DIR, "audit/requested_jobs.csv")
if (file.exists(requested_path)) expected <- read.csv(requested_path)
if (nrow(fresh)) {
  jobs <- fresh %>% group_by(scenario, imputation) %>%
    summarise(point_complete = TRUE,
      summary_finished = all(summary_attempted) & !any(status %in% c("package_fit_error", "package_summary_error")),
      package_error = any(status %in% c("package_fit_error", "package_summary_error")),
      finite_points = sum(is.finite(estimate)), package_CIs = sum(status == "package_CI_available"),
      .groups = "drop")
  expected <- expected %>% left_join(jobs, by = c("scenario", "imputation"))
} else {
  expected$point_complete <- FALSE
  expected$summary_finished <- FALSE
}
write_table(expected, file.path(RUN_DIR, "results/job_status.csv"))
format_range <- function(x) {
  if (!any(is.finite(x))) return("unavailable")
  paste(sprintf("%.3f", range(x, na.rm = TRUE)), collapse = " to ")
}
lines <- c("# Metal survival analysis", "",
  paste("Updated:", Sys.time()), "",
  "The seven-file reference is pinned to commit `da34a8bdeff9d665da38f7f07c853aaaf71bc5c3`.",
  "Inputs are rebuilt from available NHANES files. CQR models retain the study specifications, with separate survey Cox analyses.", "",
  sprintf("Fresh point jobs: %d / %d. Jobs with the summary stage resolved, including unavailable quantiles: %d / %d.",
    sum(expected$point_complete %in% TRUE), nrow(expected),
    sum(expected$summary_finished %in% TRUE), nrow(expected)),
  if ("resamples" %in% names(expected) && all(expected$resamples == 0L)) "This is an explicitly requested points-only diagnostic run; bootstrap summaries were not requested." else
    "Estimable models request 400 package jackknife resamples; models with no available requested quantile skip this step.",
  sprintf("Jobs with recorded package errors: %d.", if ("package_error" %in% names(expected)) sum(expected$package_error %in% TRUE) else 0L), "",
  "## Does the rewrite reproduce the previous result?", "")
if (nrow(control_comparison)) {
  lines <- c(lines, sprintf("Control refits on the previous frozen imputations: %d/60 cells executed; %d/60 match previous estimates, including matching unavailable cells (absolute tolerance 1e-8 years).",
      sum(control_comparison$executed), sum(control_comparison$same_estimate)),
    "This control isolates the modeling-code rewrite. It does not substitute old imputations for the fresh run.", "")
} else lines <- c(lines, "Control refits pending.", "")
if (nrow(fresh_comparison)) {
  lines <- c(lines, sprintf("Fresh analysis: %d/60 primary cells executed; %d/60 agree with the previous entries, including jointly unavailable cells, and %d/60 match their availability.",
    sum(fresh_comparison$executed), sum(fresh_comparison$same_estimate), sum(fresh_comparison$same_availability)),
    sprintf("Both interval bounds match in %d/60 primary cells with completed summaries.",
      sum(fresh_comparison$same_lower & fresh_comparison$same_upper)), "",
    "| Primary scenario | Quantile | Previous point range across MI | Fresh point range across MI | Fresh finite / executed |",
    "|---|---:|---|---|---:|")
  for (scenario in PRIMARY) for (tau in TAUS) {
    part <- fresh_comparison %>% filter(.data$scenario == .env$scenario, abs(.data$tau - .env$tau) < 1e-10)
    lines <- c(lines, sprintf("| %s | %.2f | %s | %s | %d / %d |", scenario, tau,
      format_range(part$previous_estimate), format_range(part$estimate),
      sum(is.finite(part$estimate)), sum(part$executed)))
  }
} else lines <- c(lines, "Fresh model results pending.")
lines <- c(lines, "", "## Interpretation", "",
  "CQR coefficients are conditional survival-time differences in years per exposure doubling. Negative values indicate shorter conditional survival time. RCS basis coefficients are not per-doubling effects.",
  "Package CQR intervals are imputation-specific, not Rubin-pooled or NHANES-design-adjusted. In the installed quantreg implementation, summary resampling does not carry original sampling weights into Peng-Huang bootstraps. The report preserves those package outputs and documents this limitation.",
  "Cox hazard ratios are a different estimand. They are supplemental comparisons; their numerical equality with CQR coefficients is neither expected nor tested.",
  "Missing quantiles remain unavailable. No failed result is filled from another imputation or from historical fits.",
  "Input reconstruction preserves documented inherited behavioral and dietary fields where original raw builders are unavailable; see audit/INHERITED_VARIABLES.md.", "",
  "Full rows: results_by_imputation.csv, control_vs_previous.csv, fresh_vs_previous.csv and manuscript_point_comparison.csv.")
writeLines(lines, file.path(RUN_DIR, "RUN_REPORT.md"))
