source('cause_specific/R/settings.R')
source('cause_specific/R/runtime.R')
source('scripts/cqr_fit.R')
args <- commandArgs(TRUE)
ep <- args[1];sc <- args[2];i <- as.integer(args[3]);R <- if(length(args)>=4)as.integer(args[4])else 400L
stopifnot(ep %in% endpoints,sc %in% scenarios,i %in% 1:5)
manifest <- readRDS(file.path(base,'data/manifest.rds'))
stopifnot(identical(unname(tools::md5sum(names(manifest$input_hashes))),unname(manifest$input_hashes)))
spec <- readRDS(file.path(parent,'models',sc,'specification.rds'))
d <- readRDS(file.path(parent,'models',sc,'imputations.rds'))$data[[i]]
mp <- readRDS(file.path(base,'data/endpoint_maps.rds'))[[sc]]
stopifnot(identical(as.numeric(d$SEQN),mp$SEQN),identical(d$followup_yrs,mp$followup_yrs),
  identical(as.numeric(d$event),as.numeric(mp$event_allcause)))
all_deaths <- sum(d$event);d$event <- mp[[ep]]
fit_cqr_job(d,spec$formula,spec$sampling_weights,spec$term,
  20261102L+match(ep,endpoints)*100000L+match(sc,scenarios)*1000L+i,
  list(endpoint=ep,scenario=sc,imputation=i,period=spec$period,n=nrow(d),deaths=sum(d$event),
       all_cause_deaths=all_deaths,competing_deaths=all_deaths-sum(d$event)),
  file.path(base,'fits',ep,sc,sprintf('imp%02d',i)),R,allow_sparse=TRUE)
