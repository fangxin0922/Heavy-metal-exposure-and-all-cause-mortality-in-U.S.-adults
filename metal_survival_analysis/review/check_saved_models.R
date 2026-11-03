# Check saved model hashes and result keys.
source("metal_survival_analysis/R/model_catalog.R")
hashes_match <- function(recorded) {
  length(recorded) > 0L && all(file.exists(names(recorded))) &&
    identical(unname(tools::md5sum(names(recorded))), unname(recorded))
}
model_checks <- lapply(SCENARIOS$scenario, function(scenario) {
  directory <- file.path(RUN_DIR, "models", scenario)
  manifest <- readRDS(file.path(directory, "fresh_artifact_manifest.rds"))
  specification <- readRDS(file.path(directory, "specification.rds"))
  data.frame(scenario = scenario,
    artifact_hashes_match = hashes_match(manifest$hashes),
    input_hashes_match = hashes_match(specification$input_hashes),
    R_matches = identical(manifest$R, R.version.string),
    mice_matches = identical(manifest$mice, as.character(packageVersion("mice"))))
})
model_checks <- bind_rows(model_checks)
write_table(model_checks, file.path(RUN_DIR, "audit/final_imputation_integrity.csv"))
stopifnot(nrow(model_checks) == 23L,
          all(as.matrix(model_checks[, -1, drop = FALSE])))

fit_checks <- list()
for (branch in c("fresh", "control")) {
  scenarios <- if (branch == "fresh") SCENARIOS$scenario else PRIMARY
  for (scenario in scenarios) for (imputation in seq_len(M)) {
    directory <- file.path(RUN_DIR, "fits", branch, scenario, sprintf("imp%02d", imputation))
    point <- readRDS(file.path(directory, "point.rds"))
    checked <- data.frame(branch = branch, scenario = scenario, imputation = imputation,
      point_hashes_match = hashes_match(point$input_hashes),
      point_fit_ok = !nzchar(point$fit_error),
      summary_requested = branch == "fresh", summary_hashes_match = NA,
      summary_fit_ok = NA, requested_400 = NA, table_keys_unique = NA,
      expected_rows = NA_integer_, observed_rows = NA_integer_,
      point_result_input_identity = NA, package_summary_present = NA,
      finite_points = NA_integer_)
    if (branch == "fresh") {
      result <- readRDS(file.path(directory, "result.rds"))
      specification <- readRDS(file.path(RUN_DIR, "models", scenario, "specification.rds"))
      targets <- specification$terms
      if (is.null(targets)) targets <- specification$term
      checked$summary_hashes_match <- hashes_match(result$input_hashes)
      checked$summary_fit_ok <- !nzchar(result$fit_error) && !nzchar(result$summary_error)
      checked$requested_400 <- identical(result$requested_R, BOOTSTRAPS)
      checked$table_keys_unique <- !anyDuplicated(result$table[, c("scenario", "imputation", "term", "tau")])
      checked$expected_rows <- length(targets) * length(TAUS)
      checked$observed_rows <- nrow(result$table)
      checked$point_result_input_identity <- identical(point$input_hashes, result$input_hashes)
      checked$package_summary_present <- !is.null(result$package_summary)
      checked$finite_points <- sum(is.finite(result$table$estimate))
      rm(result)
    }
    fit_checks[[length(fit_checks) + 1L]] <- checked
    rm(point)
    invisible(gc(FALSE))
  }
}
fit_checks <- bind_rows(fit_checks)
write_table(fit_checks, file.path(RUN_DIR, "audit/final_fit_integrity.csv"))
fresh_checks <- fit_checks[fit_checks$branch == "fresh", ]
stopifnot(nrow(fit_checks) == 135L, nrow(fresh_checks) == 115L,
  all(fit_checks$point_hashes_match), all(fit_checks$point_fit_ok),
  all(fresh_checks$summary_hashes_match), all(fresh_checks$summary_fit_ok),
  all(fresh_checks$requested_400), all(fresh_checks$table_keys_unique),
  all(fresh_checks$expected_rows == fresh_checks$observed_rows),
  all(fresh_checks$point_result_input_identity))
writeLines(c("PASS: 23 fresh imputation manifests, 115 fresh result objects and 20 control point fits.",
  sprintf("Actual package summaries: %d; result objects without a summary: %d.",
    sum(fresh_checks$package_summary_present), sum(!fresh_checks$package_summary_present)),
  "Saved input/artifact hashes match; no recorded fit/summary errors; expected table rows and keys verified.",
  "This verifies execution/provenance, not statistical equivalence or stability.", as.character(Sys.time())),
  file.path(RUN_DIR, "audit/FINAL_ARTIFACT_VALIDATION.txt"))
