source('cause_specific/R/settings.R')
stopifnot(!file.exists(file.path(base,'data/manifest.rds')))
a <- readRDS(file.path(parent,'data/analytic.rds'))
expected <- data.frame(scenario=scenarios,n=c(24021,18258,7061,5552),
  all_cause=c(6885,4358,1510,960),cancer=c(1503,987,391,264),cardiovascular=c(2233,1410,448,283))
maps <- counts <- list();inputs <- file.path(parent,'data/analytic.rds')
for(sc in scenarios) {
  files <- file.path(parent,'models',sc,c('specification.rds','imputations.rds','preimputation_frame.rds'))
  spec <- readRDS(files[1]);imp <- readRDS(files[2]);frame <- readRDS(files[3])
  index <- match(spec$row_ids,a$SEQN);stopifnot(!anyNA(index),length(imp$data)==5L)
  observed <- a[index,,drop=FALSE];mp <- endpoint_map(observed)
  stopifnot(spec$scenario==sc,spec$M==5L,spec$iterations==20L,spec$seed==4242L,
    identical(as.numeric(frame$SEQN),mp$SEQN),identical(as.numeric(frame$event),as.numeric(mp$event_allcause)),
    identical(frame$followup_yrs,mp$followup_yrs),length(spec$sampling_weights)==nrow(mp),
    all(is.finite(spec$sampling_weights)),all(spec$sampling_weights>0))
  for(d in imp$data) {
    stopifnot(!anyNA(d),identical(names(d),names(frame)),
      identical(as.numeric(d$SEQN),mp$SEQN),identical(as.numeric(d$event),as.numeric(mp$event_allcause)),
      identical(d$followup_yrs,mp$followup_yrs))
    for(v in names(frame)) {ok <- !is.na(frame[[v]]);stopifnot(identical(d[[v]][ok],frame[[v]][ok]))}
    if(grepl('egfr60',sc))stopifnot(all(is.finite(d$egfr_ckdepi_2021)),all(d$egfr_ckdepi_2021>=60))
  }
  want <- expected[expected$scenario==sc,]
  stopifnot(nrow(mp)==want$n,sum(mp$event_allcause)==want$all_cause,
    sum(mp$cancer)==want$cancer,sum(mp$cardiovascular)==want$cardiovascular)
  maps[[sc]] <- mp
  counts[[sc]] <- data.frame(scenario=sc,n=nrow(mp),all_cause_deaths=sum(mp$event_allcause),
    cancer_deaths=sum(mp$cancer),
    heart_disease_deaths=sum(mp$event_allcause==1L & mp$UCOD_LEADING %in% 1L),
    stroke_deaths=sum(mp$event_allcause==1L & mp$UCOD_LEADING %in% 5L),
    cardiovascular_deaths=sum(mp$cardiovascular),
    other_deaths=sum(mp$event_allcause-mp$cancer-mp$cardiovascular),alive=sum(mp$event_allcause==0L))
  inputs <- c(inputs,files)
}
saveRDS(maps,file.path(base,'data/endpoint_maps.rds'))
saveRDS(list(input_hashes=tools::md5sum(inputs),new_imputation=FALSE),file.path(base,'data/manifest.rds'))
write.csv(do.call(rbind,counts),file.path(base,'results/cohort_counts.csv'),row.names=FALSE)
