# Combine current CQR results without pooling imputations.
source("metal_survival_analysis/R/model_catalog.R")
args <- commandArgs(trailingOnly = TRUE)
scenarios <- if (length(args)) args else SCENARIOS$scenario
stopifnot(all(scenarios %in% SCENARIOS$scenario), !anyDuplicated(scenarios))

results <- bind_rows(lapply(scenarios, function(scenario) {
  bind_rows(lapply(seq_len(M), function(imputation) {
    directory <- file.path(RUN_DIR, "fits", scenario, sprintf("imp%02d", imputation))
    completed <- file.path(directory, "result.rds")
    point_file <- file.path(directory, "point_estimates.csv")
    stopifnot(file.exists(completed) || file.exists(point_file))
    table <- if (file.exists(completed)) readRDS(completed)$table else read.csv(point_file)
    stopifnot(all(table$scenario == scenario), all(table$imputation == imputation),
              !anyDuplicated(table[c("term", "tau")]))
    table$summary_attempted <- file.exists(completed)
    table
  }))
}))
write_table(results, file.path(RUN_DIR, "results", "results_by_imputation.csv"))

jobs <- results %>% group_by(scenario, imputation) %>%
  summarise(summary_attempted = all(summary_attempted),
    package_error = any(status %in% c("package_fit_error", "package_summary_error")),
    finite_points = sum(is.finite(estimate)),
    available_intervals = sum(status == "package_CI_available"), .groups = "drop")
write_table(jobs, file.path(RUN_DIR, "results", "cqr_status.csv"))

ranges <- results %>% group_by(scenario, term, tau) %>%
  summarise(imputations = n_distinct(imputation), finite_points = sum(is.finite(estimate)),
    minimum = if (any(is.finite(estimate))) min(estimate[is.finite(estimate)]) else NA_real_,
    maximum = if (any(is.finite(estimate))) max(estimate[is.finite(estimate)]) else NA_real_,
    available_intervals = sum(status == "package_CI_available"), .groups = "drop")
write_table(ranges, file.path(RUN_DIR, "results", "cqr_imputation_ranges.csv"))
