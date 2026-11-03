# Fit full-sample or survey-replicate CQR models without package summary().

source("biomarker_emm/R/settings.R")
suppressPackageStartupMessages({
  library(quantreg)
  library(survival)
})
source("biomarker_emm/R/quantile_solver.R")
source("biomarker_emm/R/numerical_retries.R")
source("biomarker_emm/R/checked_solver.R")

# Select one model and a fixed range of replicate IDs ----
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
pair <- args[1]
replicates <- seq.int(as.integer(args[2]), as.integer(args[3]))
stopifnot(length(args) == 3L, all(replicates %in% 0:400),
  Sys.getenv("OPENBLAS_NUM_THREADS") == "1", Sys.getenv("OMP_NUM_THREADS") == "1")
if (0L %in% replicates && length(replicates) > 1L) {
  stop("Run the original fit as 0 0 first, then request bootstrap IDs 1 through 400.")
}
dest <- file.path(OUTPUT_DIR, "models", pair)
prepared_path <- file.path(dest, "prepared.rds")
spec <- readRDS(prepared_path)
stopifnot(identical(unname(tools::md5sum(names(spec$source_hashes))), unname(spec$source_hashes)))
dir.create(file.path(dest, "fits"), showWarnings = FALSE)
code_files <- c("biomarker_emm/R/fit_models.R",
  "biomarker_emm/R/quantile_solver.R",
  "biomarker_emm/R/numerical_retries.R",
  "biomarker_emm/R/checked_solver.R",
  "biomarker_emm/R/settings.R", ph_native_files())
stopifnot(all(file.exists(code_files)))
hashes <- tools::md5sum(c(prepared_path, code_files))
runtime <- runtime_identity(require_linux = TRUE)
eligible <- spec$taus
original_hash <- character()
if (any(replicates > 0)) {
  original_path <- file.path(dest, "fits", "rep000.rds")
  original <- readRDS(original_path)
  stopifnot(identical(original$hashes, hashes), identical(original$runtime, runtime))
  eligible <- spec$taus[colSums(!is.finite(original$coefficients)) == 0]
  original_hash <- tools::md5sum(original_path)
}

# Replicate 0 fits the full sample; each other ID uses its fixed PSU draw.
# Validate each native solution against the numerical acceptance criteria.

for (b in replicates) {
  path <- file.path(dest, "fits", sprintf("rep%03d.rds", b))
  if (file.exists(path)) {
    old <- readRDS(path)
    stopifnot(identical(old$hashes, hashes), old$replicate == b,
              identical(old$original_hash, original_hash), identical(old$runtime, runtime))
    next
  }
  w <- if (b == 0L) spec$replicate_design$sampling_weights else spec$replicate_design$analysis_weights[, b]
  base <- checked_ph_fit(spec$prepared, w, coefficient_names = spec$terms,
    weight_constant = spec$weight_constant, keep_path = b == 0L)
  attempts <- list()
  effective <- base$coefficients
  decisions <- NULL
  if (any(!is.finite(base$coefficients[, spec$taus %in% eligible, drop = FALSE]))) {
    for (mult in retry_policy()$multipliers) {
      attempts[[paste0("mult", mult)]] <- checked_ph_fit(spec$prepared, w,
        coefficient_names = spec$terms, weight_constant = spec$weight_constant / mult,
        keep_path = b == 0L)
    }
    retry <- retry_apply(base$coefficients, attempts, eligible)
    effective <- retry$effective_coefficients
    decisions <- retry$policy_decisions
  }
  result <- list(pair = pair, replicate = b, coefficients = effective, base = base,
    attempts = attempts, retry_decisions = decisions, eligible_taus = eligible,
    hashes = hashes, original_hash = original_hash, runtime = runtime,
    retry_policy = retry_policy())
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(result, temporary)
  stopifnot(file.rename(temporary, path))
}
