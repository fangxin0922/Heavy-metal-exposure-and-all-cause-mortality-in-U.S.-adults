# Build observed EMM inputs and fixed survey draws from the main input stages.
# Run from the repository root after cohort_builder.R; model fitting is not needed.
source("metal_survival_analysis/R/model_catalog.R")
source("biomarker_emm/R/settings.R")
source("biomarker_emm/R/replicate_weights.R")

prepare_parent_inputs <- function() {
  sources <- file.path(RUN_DIR, "data", c("analytic.rds", "full_designs.rds"))
  if (!all(file.exists(sources))) {
    stop("Run the four main input stages (through cohort_builder.R) first.")
  }
  if (file.exists(INPUT_FILE) || file.exists(PARENT_DIR) ||
      file.exists(file.path(DATA_DIR, "parent_input_manifest.rds"))) {
    stop("EMM parent inputs already exist; use a new EMM_DATA_DIR.")
  }
  analytic <- readRDS(sources[1])
  designs <- readRDS(sources[2])
  stopifnot(!anyNA(analytic$SEQN), !anyDuplicated(analytic$SEQN))
  dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
  staging <- tempfile(".parent_inputs_", tmpdir = DATA_DIR)
  dir.create(file.path(staging, "parent_models"), recursive = TRUE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  checks <- list()
  for (name in PRIMARY) {
    scenario <- SCENARIOS[match(name, SCENARIOS$scenario), ]
    cohort <- select_analysis_cohort(analytic, scenario)
    component <- as.data.frame(designs[[scenario$design_metal]])
    stopifnot(nrow(cohort) > 0L, nrow(component) > 0L,
              all(cohort$SEQN %in% component$SEQN),
              all(component$SDDSRVYR >= scenario$first_cycle & component$SDDSRVYR <= 8))
    # Keep every component PSU, including those without eligible participants.
    component$analysis_domain <- component$SEQN %in% cohort$SEQN
    bank <- build_parent_replicates(component, R = 400L, seed = 20260101L)
    cohort <- cohort[match(bank$row_ids, cohort$SEQN), , drop = FALSE]
    j <- match(cohort$SEQN, component$SEQN)
    stopifnot(!anyNA(j), all(bank$sampling_weights > 0),
              identical(as.numeric(cohort$SEQN), as.numeric(bank$row_ids)),
              identical(cohort$SDMVSTRA, component$SDMVSTRA[j]),
              identical(cohort$SDMVPSU, component$SDMVPSU[j]))
    covariates <- model_covariates(scenario$metal, scenario$panel)
    spec <- list(scenario = name, metal = scenario$metal,
      design_metal = scenario$design_metal, panel = scenario$panel,
      n = nrow(cohort), events = sum(cohort$event), row_ids = bank$row_ids,
      covariates = covariates, require_egfr60 = scenario$require_egfr60,
      observed_egfr = if (scenario$require_egfr60) cohort$egfr_ckdepi_2021 else NULL,
      sampling_weights = bank$sampling_weights,
      input_hashes = tools::md5sum(sources),
      provenance = "Observed corrected data; no imputation or fitted results")
    destination <- file.path(staging, "parent_models", name)
    dir.create(destination)
    saveRDS(spec, file.path(destination, "specification.rds"))
    saveRDS(bank, file.path(destination, "replicate_design.rds"))
    checks[[name]] <- data.frame(scenario = name, n = spec$n, events = spec$events,
      full_component_n = length(bank$full_row_ids), replicates = bank$R,
      seed = bank$seed, minimum_observed_egfr = if (scenario$require_egfr60)
        min(cohort$egfr_ckdepi_2021) else NA_real_)
    message("Prepared ", name, ": ", spec$n, " participants, ", bank$R, " draws")
  }
  # Retain the exact corrected observed frame, including missing covariates.
  stopifnot(file.copy(sources[1], file.path(staging, "analytic.rds")))
  files <- list.files(staging, recursive = TRUE, full.names = TRUE)
  hashes <- tools::md5sum(files)
  names(hashes) <- substring(files, nchar(staging) + 2L)
  manifest <- list(schema = "generated_emm_parent_inputs_v1", hashes = hashes,
    source_hashes = tools::md5sum(c(sources,
      "biomarker_emm/R/prepare_parent_inputs.R", "biomarker_emm/R/replicate_weights.R")),
    cohorts = do.call(rbind, checks), R = R.version.string,
    survey = as.character(packageVersion("survey")),
    note = "Historical bank identity requires the same full-component rows, weights and runtime.")
  stopifnot(file.rename(file.path(staging, "parent_models"), PARENT_DIR),
            file.rename(file.path(staging, "analytic.rds"), INPUT_FILE))
  saveRDS(manifest, file.path(DATA_DIR, "parent_input_manifest.rds"))
  message("Prepared EMM inputs in ", DATA_DIR)
}
prepare_parent_inputs()
