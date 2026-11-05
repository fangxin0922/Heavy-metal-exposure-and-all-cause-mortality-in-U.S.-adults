# Impute covariates for one scenario using the study's MICE specification.
source("metal_survival_analysis/R/model_catalog.R")
args <- commandArgs(trailingOnly = TRUE)
requested <- if (length(args)) args[1] else PRIMARY[1]
stopifnot(requested %in% SCENARIOS$scenario)
scenario <- SCENARIOS %>% filter(.data$scenario == requested)
destination <- file.path(RUN_DIR, "models", requested)
if (dir.exists(destination) && length(list.files(destination, all.files = TRUE, no.. = TRUE))) {
  stop("Imputation output already exists; use an empty output directory for a new run.")
}
dir.create(destination, recursive = TRUE, showWarnings = FALSE)

input <- readRDS(file.path(RUN_DIR, "data/analytic.rds"))
full_designs <- readRDS(file.path(RUN_DIR, "data/full_designs.rds"))
metal <- scenario$metal
exposure <- paste0("log_primary_", metal)
terms <- exposure
covariates <- model_covariates(metal, scenario$panel)
analysis <- select_analysis_cohort(input, scenario)
if (scenario$panel == "joint") {
  terms <- c(exposure, setdiff(c("log_primary_LBXBPB", "log_primary_URXUCD"), exposure))
}
if (scenario$panel == "epa_dha") analysis <- analysis %>% mutate(epa_dha_mg = (DRXTP205 + DRXTP226) * 1000)
analysis <- droplevels(analysis)

# Match sampling weights from the full component design.
component <- full_designs[[scenario$design_metal]]
index <- match(analysis$SEQN, component$SEQN)
stopifnot(!anyNA(index), !anyDuplicated(component$SEQN))
analysis$w_metal <- component$w_metal[index]
# Preserve the numerical precision of the supplied Pb weights.
if (metal == "LBXBPB" && scenario$panel %in% c("primary", "portnoy")) {
  analysis$w_metal <- analysis$weight_LBXBPB
}
knots <- NULL
if (scenario$panel == "rcs") {
  knots <- quantile(analysis[[exposure]], c(.05, .35, .65, .95), names = FALSE, type = 7)
  basis <- restricted_cubic_basis(analysis[[exposure]], knots)
  analysis <- bind_cols(analysis, as.data.frame(basis))
  terms <- colnames(basis)
}
columns <- unique(c("SEQN", "SDDSRVYR", "SDMVSTRA", "SDMVPSU", "w_metal",
                    "followup_yrs", "event", exposure, terms, covariates))
frame <- as.data.frame(analysis[, columns])
rownames(frame) <- NULL
stopifnot(!anyDuplicated(frame$SEQN), all(frame$w_metal > 0),
          all(is.finite(frame[[exposure]])), all(frame$event %in% 0:1),
          all(frame$followup_yrs > 0))
formula <- analysis_formula(terms, covariates)
primary_scenario <- NULL
if (scenario$require_egfr60 && !scenario$primary) {
  primary_label <- if (scenario$panel == "joint") "Cd" else scenario$label
  primary_scenario <- paste0("egfr60_", primary_label)
  primary_path <- file.path(RUN_DIR, "models", primary_scenario)
  primary_spec <- readRDS(file.path(primary_path, "specification.rds"))
  primary_index <- match(frame$SEQN, primary_spec$row_ids)
  stopifnot(isTRUE(primary_spec$require_egfr60), !anyNA(primary_index),
            identical(as.numeric(frame$w_metal), as.numeric(primary_spec$sampling_weights[primary_index])),
            identical(as.numeric(analysis$egfr_ckdepi_2021), as.numeric(primary_spec$observed_egfr[primary_index])))
  if (scenario$panel %in% c("portnoy", "no_bmi")) {
    stopifnot(identical(as.numeric(frame$SEQN), as.numeric(primary_spec$row_ids)))
  }
}
specification <- list(scenario = requested, panel = scenario$panel, metal = metal,
  design_metal = scenario$design_metal, method = scenario$method,
  period = scenario$period, primary = scenario$primary, n = nrow(frame),
  events = sum(frame$event), row_ids = frame$SEQN, formula = formula,
  term = exposure, terms = terms, covariates = covariates, knots = knots,
  M = M, iterations = MAXIT, seed = MI_SEED,
  previous_scenario = scenario$previous_scenario,
  require_egfr60 = scenario$require_egfr60,
  observed_egfr = if (scenario$require_egfr60) analysis$egfr_ckdepi_2021 else NULL,
  primary_scenario = primary_scenario,
  sampling_weights = component$w_metal[index])
save_result(frame, file.path(destination, "preimputation_frame.rds"))
save_result(specification, file.path(destination, "specification.rds"))

fixed <- unique(c("SEQN", "SDDSRVYR", "SDMVSTRA", "SDMVPSU", "w_metal",
                  "followup_yrs", "event", exposure, terms, "cycle_f"))
if (scenario$require_egfr60) fixed <- union(fixed, intersect("egfr_ckdepi_2021", names(frame)))
if (scenario$panel == "portnoy" && scenario$require_egfr60) {
  # Change the estimator while retaining the exact primary completed datasets.
  primary_frame <- readRDS(file.path(primary_path, "preimputation_frame.rds"))
  stopifnot(identical(frame, primary_frame),
            identical(deparse(formula), deparse(primary_spec$formula)),
            identical(primary_spec$M, M), identical(primary_spec$iterations, MAXIT),
            identical(primary_spec$seed, MI_SEED))
  primary_imputations <- readRDS(file.path(primary_path, "imputations.rds"))
  stopifnot(identical(primary_imputations$row_ids, frame$SEQN))
  imputed <- readRDS(file.path(primary_path, "mice_object.rds"))
  completed <- primary_imputations$data
  method <- primary_imputations$method
  predictors <- primary_imputations$predictor_matrix
  provenance <- paste("Exact completed datasets reused from", primary_scenario)
} else {
  initial <- mice(frame, m = 1, maxit = 0, printFlag = FALSE, remove.collinear = FALSE)
  method <- initial$method
  predictors <- initial$predictorMatrix
  method[intersect(c("alq_5cat", "marital_3cat", "pa_q4_cyc"), names(method))] <- "polyreg"
  method[fixed] <- ""
  predictors[fixed, ] <- 0
  predictors[, c("SEQN", "SDDSRVYR", "SDMVSTRA", "SDMVPSU")] <- 0
  if (scenario$panel == "rcs") predictors[, exposure] <- 0
  cat(format(Sys.time()), "Starting MICE:", requested, "N=", nrow(frame), "\n")
  imputed <- mice(frame, m = M, maxit = MAXIT, method = method,
                  predictorMatrix = predictors, seed = MI_SEED,
                  remove.collinear = FALSE, printFlag = TRUE)
  completed <- lapply(seq_len(M), function(i) as.data.frame(complete(imputed, i)))
  provenance <- "Fresh MICE after observed cohort eligibility"
}
stopifnot(length(completed) == M)
for (dataset in completed) {
  stopifnot(!anyNA(dataset), all(vapply(fixed, function(variable) {
    identical(dataset[[variable]], frame[[variable]])
  }, logical(1))))
  if (scenario$require_egfr60 && !scenario$primary) {
    stopifnot(all(vapply(names(frame), function(variable) {
      observed <- !is.na(frame[[variable]])
      identical(dataset[[variable]][observed], frame[[variable]][observed])
    }, logical(1))))
    design_matrix <- model.matrix(delete.response(stats::terms(formula)), dataset)
    stopifnot(qr(design_matrix)$rank == ncol(design_matrix))
  }
}
save_result(imputed, file.path(destination, "mice_object.rds"))
save_result(list(data = completed, row_ids = frame$SEQN, method = method,
                 predictor_matrix = predictors, seed = MI_SEED, iterations = MAXIT,
                 provenance = provenance),
            file.path(destination, "imputations.rds"))
events <- imputed$loggedEvents
if (is.null(events)) events <- data.frame(it = integer(), im = integer(), dep = character(), meth = character(), out = character())
write_table(events, file.path(destination, "mice_logged_events.csv"))
pdf(file.path(destination, "mice_trace.pdf"), width = 11, height = 8.5)
print(plot(imputed))
dev.off()
