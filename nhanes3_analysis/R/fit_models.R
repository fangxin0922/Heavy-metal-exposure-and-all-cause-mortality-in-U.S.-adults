source('nhanes3_analysis/R/settings.R')
source('scripts/cqr_fit.R')
args <- commandArgs(TRUE)
metal <- args[1];i <- as.integer(args[2]);R <- if(length(args)>=3)as.integer(args[3])else 400L
stopifnot(metal %in% scenarios,i %in% 1:5)
model <- file.path(base,'models',metal)
spec <- readRDS(file.path(model,'specification.rds'))
d <- readRDS(file.path(model,'imputations.rds'))$data[[i]]
stopifnot(identical(as.numeric(d$SEQN),as.numeric(spec$row_ids)),nrow(d)==spec$n,sum(d$event)==spec$events)
if(metal=='cadmium')stopifnot(all(is.finite(d$egfr)),all(d$egfr>=60))
if(metal=='lead')stopifnot(!'egfr' %in% all.vars(spec$formula))
fit_cqr_job(d,spec$formula,spec$sampling_weights,spec$term,
  20261030L+match(metal,scenarios)*1000L+i,
  list(scenario=metal,imputation=i,period=spec$period,n=spec$n,deaths=spec$events),
  file.path(base,'fits',metal,sprintf('imp%02d',i)),R)
