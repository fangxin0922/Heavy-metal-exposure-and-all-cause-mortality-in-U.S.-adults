# A missing interval must be reported even when all point predictions still agree.
tmp <- tempfile('spline_comparison_')
dir.create(tmp)
x <- expand.grid(scenario='test',imputation=1:5,tau=c(.1,.25,.5),concentration=c(1,2,3))
x$n <- 100L;x$deaths <- 40L
x$estimate <- 5+x$concentration;x$lower <- x$estimate-1;x$upper <- x$estimate+1
y <- x
y$lower[1:2] <- NA_real_;y$upper[2] <- NA_real_
write.csv(x,file.path(tmp,'reference.csv'),row.names=FALSE)
write.csv(y,file.path(tmp,'fresh.csv'),row.names=FALSE)
status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',
  'validation/compare_splines.R',shQuote(file.path(tmp,'reference.csv')),
  shQuote(file.path(tmp,'fresh.csv')),shQuote(tmp)))
stopifnot(status==0L)
z <- read.csv(file.path(tmp,'spline_comparison.csv'))
stopifnot(nrow(z)==15L,all(z$maximum_prediction_difference_years==0),
  sum(z$interval_availability_changes)==2L,
  sum(z$lower_availability_changes)==2L,sum(z$upper_availability_changes)==1L,
  sum(z$reference_interval_rows)==45L,sum(z$fresh_interval_rows)==43L)
unlink(tmp,recursive=TRUE)
cat('Spline comparison availability checks passed.\n')
