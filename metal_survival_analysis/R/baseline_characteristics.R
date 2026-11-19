# Run from the repository root.
# Aggregate observed baseline values in the four current primary cohorts.
source("metal_survival_analysis/R/model_catalog.R")
OUT <- file.path(RUN_DIR, "results", "baseline")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
analytic_path <- file.path(RUN_DIR, "data/analytic.rds")
analytic <- readRDS(analytic_path)
stopifnot(!anyDuplicated(analytic$SEQN))
rows <- coverage <- cohorts <- list()
sources <- analytic_path
continuous <- c("RIDAGEYR", "age_squared", "INDFMPIR", "BMXBMI", "pack_years",
                "LBXTC", "egfr_ckdepi_2021", "fish_servings_wk", "exposure", "followup_yrs")
categorical <- c("RIAGENDR_f", "RIDRETH1_f", "DMDEDUC2_f", "marital_3cat", "alq_5cat",
                 "pa_q4_cyc", "cycle_f")
expected_n <- c(24021L, 18258L, 7061L, 5552L)
expected_events <- c(6885L, 4358L, 1510L, 960L)
for (i in seq_along(PRIMARY)) {
  scenario <- PRIMARY[i]
  folder <- file.path(RUN_DIR, "models", scenario)
  paths <- file.path(folder, c("specification.rds", "preimputation_frame.rds"))
  spec <- readRDS(paths[1])
  frame <- readRDS(paths[2])
  stopifnot(identical(as.numeric(frame$SEQN), as.numeric(spec$row_ids)),
            nrow(frame) == expected_n[i], sum(frame$event) == expected_events[i])
  index <- match(frame$SEQN, analytic$SEQN)
  stopifnot(!anyNA(index))
  # Describe the union of model covariates in each cohort, including covariates
  # used only in other metals' models. Preserve saved preimputation values.
  for (v in c("egfr_ckdepi_2021", "fish_servings_wk")) {
    if (!v %in% names(frame)) frame[[v]] <- analytic[[v]][index]
    else stopifnot(isTRUE(all.equal(frame[[v]], analytic[[v]][index])))
  }
  frame$age_squared <- frame$RIDAGEYR^2
  frame$exposure <- 2^frame[[spec$term]]
  # Match the exact sampling-weight representation consumed by the fresh fit.
  stopifnot(length(spec$sampling_weights) == nrow(frame),
            isTRUE(all.equal(frame$w_metal, spec$sampling_weights, tolerance = 1e-12)))
  frame$w_metal <- spec$sampling_weights
  frame$descriptive_stratum <- interaction(frame$SDDSRVYR, frame$SDMVSTRA, drop = TRUE)
  stopifnot(all(is.finite(frame$w_metal)), all(frame$w_metal > 0),
            all(spec$covariates %in% c(continuous, categorical)))
  if (scenario %in% c("egfr60_Cd", "egfr60_As")) {
    stopifnot(all(is.finite(frame$egfr_ckdepi_2021)), all(frame$egfr_ckdepi_2021 >= 60))
  }
  design <- survey::svydesign(ids = ~SDMVPSU, strata = ~descriptive_stratum,
                             weights = ~w_metal, data = frame, nest = TRUE)
  append_row <- function(variable, level, type, n, percent = NA_real_,
                         q25 = NA_real_, median = NA_real_, q75 = NA_real_) {
    rows[[length(rows) + 1L]] <<- data.frame(scenario, metal = names(METALS)[i], variable, level,
      type, n, weighted_percent = percent, q25, median, q75,
      modeled_covariate = variable %in% c(spec$covariates, "age_squared"))
  }
  for (v in continuous) {
    x <- frame[[v]]
    observed <- is.finite(x)
    stopifnot(all(observed | is.na(x)), any(observed))
    qs <- as.numeric(coef(survey::svyquantile(reformulate(v), design,
      quantiles = c(.25, .5, .75), ci = FALSE, na.rm = TRUE, qrule = "math")))
    stopifnot(length(qs) == 3L, all(is.finite(qs)), all(diff(qs) >= 0))
    # Independent inverse weighted-ECDF check of the survey point quantiles.
    xx <- x[observed]; ww <- frame$w_metal[observed]; ord <- order(xx)
    cumulative <- cumsum(ww[ord]) / sum(ww)
    check <- vapply(c(.25, .5, .75), function(p) xx[ord][which(cumulative >= p)[1L]], numeric(1))
    stopifnot(max(abs(check - qs)) < 1e-8)
    append_row(v, "Observed", "continuous", sum(observed), q25 = qs[1], median = qs[2], q75 = qs[3])
    append_row(v, "Missing", "missing", sum(!observed),
               100 * sum(frame$w_metal[!observed]) / sum(frame$w_metal))
  }
  for (v in categorical) {
    x <- as.character(frame[[v]])
    levels_v <- levels(frame[[v]])
    if (v == "cycle_f") levels_v <- as.character(1:8)
    total_n <- 0L; total_pct <- 0
    # Category percentages include missing observations in the denominator.
    for (lev in c(levels_v, "Missing")) {
      selected <- if (lev == "Missing") is.na(x) else !is.na(x) & x == lev
      count <- sum(selected)
      pct <- 100 * sum(frame$w_metal[selected]) / sum(frame$w_metal)
      append_row(v, lev, if (lev == "Missing") "missing" else "categorical", count, pct)
      total_n <- total_n + count; total_pct <- total_pct + pct
    }
    stopifnot(total_n == nrow(frame), abs(total_pct - 100) < 1e-8)
  }
  coverage[[i]] <- data.frame(scenario, variable = c(spec$covariates, "age_squared"),
                             included = TRUE)
  cohorts[[i]] <- data.frame(scenario, metal = names(METALS)[i], n = spec$n, deaths = spec$events,
    period = spec$period, weight_sum = sum(frame$w_metal),
    observed_egfr_restriction = if (i >= 3) ">=60" else "None")
  sources <- c(sources, paths)
}
write_table(bind_rows(rows), file.path(OUT, "table1_baseline_long.csv"))
write_table(bind_rows(coverage), file.path(OUT, "table1_covariate_coverage.csv"))
write_table(bind_rows(cohorts), file.path(OUT, "table1_cohorts.csv"))
sources <- unique(c(sources, "metal_survival_analysis/R/baseline_characteristics.R"))
write_table(data.frame(file=sources,md5=unname(tools::md5sum(sources))),
  file.path(OUT, "table1_source_hashes.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "baseline_session_info.txt"))
writeLines(c("PASS: all four primary cohort counts and ordered membership verified; input hashes recorded.",
  "PASS: all model covariates, including the derived age-squared term, are represented.",
  "PASS: all categorical counts sum to cohort N and weighted percentages sum to 100 before rounding.",
  "PASS: all weighted quartiles agree with an independent inverse weighted-ECDF calculation.",
  "Observed preimputation values only; no model refits or new imputations."),
  file.path(OUT, "baseline_validation.txt"))
cat("Baseline summaries complete for four primary cohorts.\n")
