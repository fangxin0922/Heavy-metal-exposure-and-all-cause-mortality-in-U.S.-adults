# Descriptive appendix checks; no outcome fitting or changes to stored covariates.
source('metal_survival_analysis/R/model_catalog.R')
out <- file.path(RUN_DIR,'results/supplementary_descriptives')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
input <- readRDS(file.path(RUN_DIR,'data/analytic.rds'))
archive <- readRDS(file.path(INPUT_DIR,'analytic.rds'))
a <- archive[archive$SDDSRVYR<=8,,drop=FALSE]
stopifnot(nrow(a)==24339L,!anyDuplicated(a$SEQN))
# Reconstruct the inherited, within-cycle rank groups, breaking ties by row order.
rows <- list();quartile <- rep(NA_integer_,nrow(a))
for(cycle in sort(unique(a$SDDSRVYR))) {
  ix <- which(a$SDDSRVYR==cycle & is.finite(a$pa_mets_wkly_all))
  ix <- ix[order(a$pa_mets_wkly_all[ix],ix)]
  n <- length(ix);sizes <- rep(n%/%4,4)
  if(n%%4>0)sizes[seq_len(n%%4)] <- sizes[seq_len(n%%4)]+1L
  quartile[ix] <- rep(1:4,sizes)
  for(q in 1:4) {
    values <- a$pa_mets_wkly_all[ix[quartile[ix]==q]]
    rows[[length(rows)+1L]] <- data.frame(cycle,quartile=q,n=length(values),
      minimum=min(values),maximum=max(values))
  }
}
stored <- as.integer(gsub('Q','',as.character(a$pa_q4_cyc)))
stopifnot(identical(is.na(quartile),is.na(stored)),all(quartile==stored,na.rm=TRUE),
  sum(!is.na(quartile))==23845L)
write_table(bind_rows(rows),file.path(out,'activity_ranges.csv'))

spec <- readRDS(file.path(RUN_DIR,'models/primary_Hg/specification.rds'))
ix <- match(spec$row_ids,input$SEQN);stopifnot(!anyNA(ix))
h <- input[ix,,drop=FALSE]
sources <- c('screener_species','screener_zero','recall_2day','recall_1day','missing')
tiers <- c('Tier 1','Tier 1','Tier 2','Tier 2','Missing')
stopifnot(nrow(h)==18258L,sum(h$event)==4358L,!anyNA(h$fish_source),
  all(as.character(h$fish_source)%in%sources))
t1 <- !is.na(h$fish_servings_wk_t1)
t2 <- !t1 & !is.na(h$fish_servings_wk_t2)
missing <- !t1 & is.na(h$fish_servings_wk_t2)
stopifnot(identical(t1,as.character(h$fish_source)%in%sources[1:2]),
  identical(t2,as.character(h$fish_source)%in%sources[3:4]),
  identical(missing,as.character(h$fish_source)=='missing'),
  identical(is.na(h$fish_servings_wk),missing),
  all(h$fish_servings_wk[t1]==h$fish_servings_wk_t1[t1]),
  all(h$fish_servings_wk[t2]==h$fish_servings_wk_t2[t2]))
counts <- as.integer(table(factor(h$fish_source,levels=sources)))
write_table(data.frame(source=sources,tier=tiers,n=counts,denominator=nrow(h),
  percent_unweighted=100*counts/nrow(h)),file.path(out,'fish_source_counts.csv'))
tier_names <- c('Tier 1','Tier 2','Missing')
counts <- vapply(tier_names,function(z)sum(counts[tiers==z]),integer(1))
stopifnot(identical(unname(counts),c(14757L,2505L,996L)))
write_table(data.frame(tier=tier_names,n=unname(counts),denominator=nrow(h),
  percent_unweighted=100*unname(counts)/nrow(h)),file.path(out,'fish_tier_counts.csv'))

rows <- list()
for(sc in SCENARIOS$scenario[SCENARIOS$require_egfr60 & !SCENARIOS$primary]) {
  scenario <- SCENARIOS[SCENARIOS$scenario==sc,]
  unrestricted <- scenario;unrestricted$require_egfr60 <- FALSE
  before <- select_analysis_cohort(input,unrestricted)
  after <- select_analysis_cohort(input,scenario)
  rows[[sc]] <- data.frame(scenario=sc,n=nrow(after),deaths=sum(after$event),
    excluded_missing_egfr=sum(!is.finite(before$egfr_ckdepi_2021)),
    excluded_egfr_below60=sum(is.finite(before$egfr_ckdepi_2021)&before$egfr_ckdepi_2021<60))
}
write_table(bind_rows(rows),file.path(out,'restricted_sensitivity_cohorts.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
