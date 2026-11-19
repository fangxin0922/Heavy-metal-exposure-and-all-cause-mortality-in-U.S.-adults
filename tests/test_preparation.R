# Synthetic integration checks: no participant data or outcome fits are needed.
library_path <- Sys.getenv("HM_R_LIBRARY", "")
if (nzchar(library_path)) .libPaths(c(library_path, .libPaths()))
root <- normalizePath(".")
source("scripts/input_manifest.R")
source("biomarker_emm/R/replicate_weights.R")
stopifnot(nrow(public_input_manifest(TRUE)) == 76L)

# Compare domain draws with an independently generated full-design draw bank.
d <- expand.grid(person = 1:10, psu = 1:2, stratum = 1:6)
d$SEQN <- seq_len(nrow(d)); d$w_metal <- 1 + d$person / 10
d$SDMVSTRA <- d$stratum; d$SDMVPSU <- d$psu
d$analysis_domain <- d$person <= 5 & !(d$stratum == 1 & d$psu == 2)
set.seed(987); old_seed <- .Random.seed; old_kind <- RNGkind()
bank <- build_parent_replicates(d)
stopifnot(identical(old_seed, .Random.seed), identical(old_kind, RNGkind()),
          identical(bank$analysis_weights, build_parent_replicates(d)$analysis_weights))
design <- survey::svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
  weights = ~w_metal, data = d, nest = TRUE)
set.seed(20260101)
full <- survey::as.svrepdesign(design, type = "subbootstrap", replicates = 400, mse = FALSE)
expected <- as.matrix(weights(full, type = "analysis"))[d$analysis_domain, , drop = FALSE]
stopifnot(isTRUE(all.equal(unname(bank$analysis_weights), unname(expected), tolerance = 0)),
          bank$design_checks$n_psus_without_domain == 1,
          isTRUE(bank$design_checks$withReplicates_agreement))
bad <- d[!(d$stratum == 1 & d$psu == 2), ]
stopifnot(inherits(try(build_parent_replicates(bad), silent = TRUE), "try-error"))

# Run the real preparation entry points in an isolated temporary checkout.
sandbox <- tempfile("hm-preparation-")
dir.create(sandbox)
for (folder in c("metal_survival_analysis", "biomarker_emm")) {
  dir.create(file.path(sandbox, folder, "R"), recursive = TRUE)
  stopifnot(all(file.copy(list.files(file.path(folder, "R"), full.names = TRUE),
    file.path(sandbox, folder, "R"))))
}
Sys.setenv(EMM_DATA_DIR = file.path(sandbox, "data/biomarker_emm"))
setwd(sandbox)
source("metal_survival_analysis/R/model_catalog.R")
pick <- function(name) SCENARIOS[SCENARIOS$scenario == name, ]
stopifnot(nrow(SCENARIOS) == 23, !"rcs_Cd" %in% SCENARIOS$scenario,
  identical(vapply(1:5, function(i) cqr_summary_seed("rcs_egfr60_Cd", i), integer(1)),
            20262105:20262109))
boundary <- data.frame(SEQN = 1:8, SDDSRVYR = 1,
  domain_URXUCD = c(rep(TRUE, 7), FALSE),
  egfr_ckdepi_2021 = c(60, 59.99, NA, NaN, Inf, -Inf, 90, 90))
stopifnot(identical(select_analysis_cohort(boundary, pick("rcs_egfr60_Cd"))$SEQN, c(1L, 7L)))

set.seed(7521)
n <- 640L
fixture <- data.frame(SEQN = seq_len(n), SDDSRVYR = rep(1:8, each = 80),
  SDMVSTRA = rep(1:16, each = 40), SDMVPSU = rep(rep(1:2, each = 20), 16),
  RIDAGEYR = sample(40:85, n, TRUE), followup_yrs = runif(n, 1, 20), event = rbinom(n, 1, .4),
  ELIGSTAT = 1, log_primary_LBXBPB = rnorm(n), log_primary_LBXTHG = rnorm(n),
  log_primary_URXUCD = rnorm(n), log_primary_VNSUMARS = rnorm(n),
  domain_LBXBPB = TRUE, domain_LBXTHG = rep(1:8, each = 80) >= 3,
  domain_URXUCD = TRUE, domain_VNSUMARS = rep(1:8, each = 80) >= 3,
  egfr_ckdepi_2021 = sample(c(45, 59.9, 60, 65, 90, NA, Inf), n, TRUE),
  INDFMPIR = runif(n, .2, 5), BMXBMI = runif(n, 18, 40), pack_years = runif(n, 0, 40),
  LBXTC = rnorm(n, 200, 25), fish_servings_wk = runif(n), weight_LBXBPB = runif(n, 10, 100))
for (v in c("RIAGENDR_f", "RIDRETH1_f", "DMDEDUC2_f", "marital_3cat", "alq_5cat", "pa_q4_cyc")) {
  fixture[[v]] <- factor(sample(c("A", "B"), n, TRUE))
}
fixture$INDFMPIR[sample(n, 30)] <- NA_real_
fixture$BMXBMI[sample(n, 30)] <- NA_real_
designs <- lapply(names(METALS), function(metal) {
  fd <- fixture[c("SEQN", "SDDSRVYR", "SDMVSTRA", "SDMVPSU")]
  fd$w_metal <- fixture$weight_LBXBPB
  if (metal %in% c("Hg", "As")) fd <- fd[fd$SDDSRVYR >= 3, ]
  fd
})
names(designs) <- unname(METALS)
saveRDS(fixture, file.path(RUN_DIR, "data/analytic.rds"))
saveRDS(designs, file.path(RUN_DIR, "data/full_designs.rds"))
run_r <- function(script, args = character(), success = TRUE) {
  log <- tempfile(tmpdir = sandbox, fileext = ".log")
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script, args), stdout = log, stderr = log)
  if ((status == 0L) != success) stop(paste(readLines(log), collapse = "\n"))
}
run_r("metal_survival_analysis/R/impute_covariates.R", "egfr60_Cd")
run_r("metal_survival_analysis/R/impute_covariates.R", "rcs_egfr60_Cd")
model <- function(name, file) readRDS(file.path(RUN_DIR, "models", name, paste0(file, ".rds")))
primary <- model("egfr60_Cd", "specification")
spline <- model("rcs_egfr60_Cd", "specification")
frame <- model("rcs_egfr60_Cd", "preimputation_frame")
imps <- model("rcs_egfr60_Cd", "imputations")
expected <- fixture[is.finite(fixture$egfr_ckdepi_2021) & fixture$egfr_ckdepi_2021 >= 60, ]
stopifnot(identical(spline$row_ids, primary$row_ids),
  identical(spline$sampling_weights, primary$sampling_weights),
  identical(spline$observed_egfr, primary$observed_egfr),
  identical(spline$knots, quantile(expected$log_primary_URXUCD, c(.05, .35, .65, .95), names = FALSE)),
  identical(imps$provenance, "Fresh MICE after observed cohort eligibility"), length(imps$data) == 5)
for (d in imps$data) stopifnot(!anyNA(d), identical(d$egfr_ckdepi_2021, frame$egfr_ckdepi_2021),
  identical(d[c("rcs1", "rcs2", "rcs3")], frame[c("rcs1", "rcs2", "rcs3")]))
run_r("metal_survival_analysis/R/impute_covariates.R", "rcs_Cd", success = FALSE)
run_r("metal_survival_analysis/R/impute_covariates.R", "rcs_egfr60_Cd", success = FALSE)

run_r("biomarker_emm/R/prepare_parent_inputs.R")
emm <- Sys.getenv("EMM_DATA_DIR")
stopifnot(identical(readRDS(file.path(emm, "analytic.rds")), fixture))
for (name in PRIMARY) {
  spec <- readRDS(file.path(emm, "parent_models", name, "specification.rds"))
  rp <- readRDS(file.path(emm, "parent_models", name, "replicate_design.rds"))
  expected <- select_analysis_cohort(fixture, pick(name))
  stopifnot(setequal(rp$row_ids, expected$SEQN), identical(spec$row_ids, rp$row_ids),
    spec$n == nrow(expected), spec$events == sum(expected$event), rp$R == 400,
    identical(rp$mse, FALSE), all(rp$sampling_weights > 0), all(colSums(rp$analysis_weights) > 0))
}
run_r("biomarker_emm/R/prepare_parent_inputs.R", success = FALSE)
setwd(root)
unlink(sandbox, recursive = TRUE)
cat("PASS: renal eligibility, fresh spline preparation, preserved seeds, full-design survey draws, parent inputs and overwrite guards.\n")
