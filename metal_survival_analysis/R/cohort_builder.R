# Join inputs by participant and cycle, then define analysis cohorts.
# Use the previous corrected dataset only to validate the completed build.
source("metal_survival_analysis/R/settings.R")
source(file.path(RUN_DIR, "R/data_tools.R"))

chemicals <- readRDS(data_file("chemicals.rds"))
mortality <- readRDS(data_file("mortality.rds"))
covariates <- readRDS(data_file("covariates.rds"))
keys <- c("SEQN", "SDDSRVYR")
replacement_fields <- union(setdiff(names(chemicals), keys), setdiff(names(mortality), keys))
analytic <- covariates %>%
  dplyr::select(-dplyr::any_of(replacement_fields)) %>%
  dplyr::left_join(chemicals, by=keys) %>%
  dplyr::left_join(mortality, by=keys)
# Keep pre-correction values from the source archive as audit columns.
historical <- read_inherited()
j <- match(analytic$SEQN, historical$SEQN)
for (field in c("URXUCR", "LBXSCR"))
  analytic[[paste0(field, "_original_averaged")]] <- historical[[field]][j]
for (field in c("URXUCR_pred", "URXUCR_cs_ratio", "URXUCD_cs", "VNSUMARS_cs", "URXUTL_cs",
                "log_primary_URXUCD", "log_primary_VNSUMARS", "log_primary_URXUTL",
                "LBXSCR_cal", "egfr_ckdepi_2021", "fish_servings_wk", "fish_source"))
  analytic[[paste0(field, "_original")]] <- historical[[field]][j]
analytic$VNSUMARS_original_two_species <- historical$VNSUMARS[j]
stopifnot(nrow(analytic) == nrow(covariates), !anyDuplicated(analytic$SEQN),
          all(analytic$RIDAGEYR >= 40), all(analytic$ELIGSTAT == 1),
          all(analytic$event %in% 0:1), all(analytic$followup_yrs > 0),
          !anyNA(analytic[c("SDMVSTRA", "SDMVPSU")]))
full_designs <- readRDS(data_file("full_designs.rds"))
for (metal in metal_names) {
  first_cycle <- if (metal %in% c("LBXTHG", "VNSUMARS")) 3L else 1L
  exposure <- paste0("log_primary_", metal)
  weight <- paste0("weight_", metal)
  domain <- paste0("domain_", metal)
  analytic[[domain]] <- analytic$SDDSRVYR >= first_cycle & analytic$SDDSRVYR <= 8 &
    is.finite(analytic[[exposure]]) & is.finite(analytic[[weight]]) &
    analytic[[weight]] > 0 & analytic$RIDAGEYR >= 40 &
    stats::complete.cases(analytic[c("followup_yrs", "event", "SDMVSTRA", "SDMVPSU")])
  full_designs[[metal]]$domain_primary <- full_designs[[metal]]$SEQN %in% analytic$SEQN[analytic[[domain]]]
}
analytic$primary_Pb <- analytic$domain_LBXBPB
analytic$primary_Hg <- analytic$domain_LBXTHG
analytic$primary_Cd <- analytic$domain_URXUCD & is.finite(analytic$egfr_ckdepi_2021) & analytic$egfr_ckdepi_2021 >= 60
analytic$primary_As <- analytic$domain_VNSUMARS & is.finite(analytic$egfr_ckdepi_2021) & analytic$egfr_ckdepi_2021 >= 60
union <- rowSums(is.finite(as.matrix(analytic[metal_names]))) > 0
analytic <- droplevels(analytic)
saveRDS(analytic, data_file("analytic_corrected.rds"))
saveRDS(analytic, data_file("analytic_fresh.rds"))
saveRDS(analytic[union, ], data_file("analytic_four_metals.rds"))
saveRDS(full_designs, data_file("full_designs.rds"))

cohorts <- dplyr::bind_rows(lapply(c("primary_Pb", "primary_Hg", "primary_Cd", "primary_As"), function(scenario) {
  selected <- analytic[[scenario]]
  data.frame(scenario=scenario, n=sum(selected), deaths=sum(analytic$event[selected]))
}))
write.csv(cohorts, audit_file("fresh_primary_cohorts.csv"), row.names=FALSE)
stopifnot(identical(cohorts$n, c(24021L, 18258L, 7061L, 5552L)),
          all(cohorts$deaths == c(6885, 4358, 1510, 960)), sum(union) == 24338L)

# Read previous results only for validation after saving the new inputs.
comparator_file <- file.path(PROJECT_ROOT, "replication_through_2014/manuscript_completion/data/analytic_corrected.rds")
previous <- readRDS(comparator_file)
stopifnot(identical(as.numeric(analytic$SEQN), as.numeric(previous$SEQN)),
          identical(as.numeric(analytic$SDDSRVYR), as.numeric(previous$SDDSRVYR)))
model_fields <- unique(c(keys, "SDMVSTRA", "SDMVPSU", "ELIGSTAT", "MORTSTAT", "PERMTH_EXM",
  "event", "followup_yrs", metal_names, paste0("log_primary_", metal_names),
  paste0("weight_", metal_names), paste0("domain_", metal_names),
  "URXUCR", "URXUCR_pred", "URXUCR_cs_ratio", "URXUCD_cs", "VNSUMARS_cs",
  "log_VNSUMARS", "z_log_VNSUMARS_barr", "z_log_URXUCD", "z_log_VNSUMARS", "z_log_URXUTL",
  "RIDAGEYR", "RIAGENDR_f", "RIDRETH1_f", "DMDEDUC2_f", "marital_3cat", "INDFMPIR",
  "BMXBMI", "pack_years", "alq_5cat", "pa_q4_cyc", "LBXTC", "LBXSCR_cal",
  "egfr_ckdepi_2021", "fish_servings_wk", "fish_source", "DRXTP205", "DRXTP226"))
comparison <- dplyr::bind_rows(lapply(model_fields, function(field)
  compare_column(analytic[[field]], previous[[field]], field)))
write.csv(comparison, audit_file("fresh_vs_previous_inputs.csv"), row.names=FALSE)
stopifnot(all(comparison$pass), all(comparison$factor_levels_identity))
design_comparator <- file.path(PROJECT_ROOT, "replication_through_2014/data/full_designs.rds")
old_designs <- readRDS(design_comparator)
design_checks <- dplyr::bind_rows(lapply(metal_names, function(metal) {
  dplyr::bind_rows(lapply(c("SEQN", "SDDSRVYR", "w_native", "w_metal", "SDMVSTRA", "SDMVPSU"), function(field) {
    check <- compare_column(full_designs[[metal]][[field]], old_designs[[metal]][[field]], field)
    check$metal <- metal
    check
  }))
}))
write.csv(design_checks, audit_file("fresh_vs_previous_full_designs.csv"), row.names=FALSE)
stopifnot(all(design_checks$exact_value_identity))
record_input_sources(c(data_file(c("chemicals.rds", "mortality.rds", "covariates.rds")),
                       comparator_file, design_comparator), "04_assembly_validation")
writeLines(c(paste("Validation completed:", Sys.time()),
             paste(nrow(comparison), "model-relevant fields matched in values/missingness at tolerance 1e-10."),
             paste(sum(comparison$exact_value_identity), "fields matched exactly in values."),
             "All categorical level sets and the four primary cohort counts matched.",
             "All four full-design ID, cycle, native weight, pooled weight, stratum and PSU vectors matched exactly."),
           audit_file("INPUT_RECONSTRUCTION_COMPLETE.txt"))
message("Fresh input reconstruction validated against ", nrow(comparison), " model-relevant fields.")
