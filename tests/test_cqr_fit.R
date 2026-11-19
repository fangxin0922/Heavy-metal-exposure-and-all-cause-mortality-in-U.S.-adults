library_path <- Sys.getenv('HM_R_LIBRARY','')
if(nzchar(library_path)) .libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages({library(quantreg);library(survival)})
source('scripts/cqr_fit.R')
set.seed(91)
n <- 160L
x <- rnorm(n)
failure <- rexp(n,exp(.2*x));censor <- rexp(n,.3)
d <- data.frame(time=pmin(failure,censor),event=as.integer(failure<=censor),x=x)
destination <- tempfile('cqr_interval_test_')
fit_cqr_job(d,Surv(time,event)~x,rep(1,n),'x',123L,
            list(scenario='synthetic',imputation=1L,n=n,deaths=sum(d$event)),destination,400L)
r <- readRDS(file.path(destination,'result.rds'))
stopifnot(nrow(r$table)==3L,r$requested_R==400L,
          all(r$table$status=='package_CI_available'),all(is.finite(r$table$p_value)))
# The saved fit also retains enough context for public model.frame reconstruction.
point <- readRDS(file.path(destination,'point.rds'))
stopifnot(nrow(model.frame(point$fit))==n)
cat('Shared fitting helper completed all 400 resamples with reconstructible model context.\n')
