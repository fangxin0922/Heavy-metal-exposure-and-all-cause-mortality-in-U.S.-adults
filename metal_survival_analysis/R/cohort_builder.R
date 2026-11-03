# Join inputs by participant and cycle, then define analysis cohorts.
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
saveRDS(analytic, data_file("analytic.rds"))
saveRDS(analytic[union, ], data_file("analytic_four_metals.rds"))
saveRDS(full_designs, data_file("full_designs.rds"))

cohorts <- dplyr::bind_rows(lapply(c("primary_Pb", "primary_Hg", "primary_Cd", "primary_As"), function(scenario) {
  selected <- analytic[[scenario]]
  data.frame(scenario=scenario, n=sum(selected), deaths=sum(analytic$event[selected]))
}))
write.csv(cohorts, diagnostic_file("primary_cohorts.csv"), row.names=FALSE)
stopifnot(identical(cohorts$n, c(24021L, 18258L, 7061L, 5552L)),
          all(cohorts$deaths == c(6885, 4358, 1510, 960)), sum(union) == 24338L)
