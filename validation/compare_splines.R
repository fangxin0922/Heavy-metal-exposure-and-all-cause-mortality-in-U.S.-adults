# Usage: Rscript validation/compare_splines.R reference.csv fresh.csv output_directory
# The reference CSV is supplied locally and is never used to calculate fitted results.
args <- commandArgs(TRUE);stopifnot(length(args)==3L)
a <- read.csv(args[1]);b <- read.csv(args[2]);out <- args[3]
keys <- c('scenario','imputation','tau','concentration')
stopifnot(!anyDuplicated(a[keys]),!anyDuplicated(b[keys]))
z <- merge(a,b,by=keys,suffixes=c('_reference','_fresh'))
stopifnot(nrow(z)==nrow(a),nrow(z)==nrow(b),
  all(z$n_reference==z$n_fresh),all(z$deaths_reference==z$deaths_fresh))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
groups <- split(z,interaction(z$scenario,z$imputation,z$tau,drop=TRUE))
maximum_difference <- function(a,b) {
  ok <- is.finite(a)&is.finite(b)
  if(any(ok))max(abs(a[ok]-b[ok]))else NA_real_
}
report <- do.call(rbind,lapply(groups,function(x) {
  x <- x[order(x$concentration),]
  old <- x$estimate_reference;new <- x$estimate_fresh
  old_interval <- is.finite(x$lower_reference)&is.finite(x$upper_reference)
  new_interval <- is.finite(x$lower_fresh)&is.finite(x$upper_fresh)
  stopifnot(all(is.finite(old))||all(!is.finite(old)),
            all(is.finite(new))||all(!is.finite(new)))
  d0 <- tail(old,1)-old[1];d1 <- tail(new,1)-new[1]
  threshold <- if(is.finite(d0))max(.1,.2*abs(d0))else NA_real_
  data.frame(scenario=x$scenario[1],imputation=x$imputation[1],tau=x$tau[1],
    grid_rows=nrow(x),reference_available=all(is.finite(old)),
    fresh_available=all(is.finite(new)),
    reference_interval_rows=sum(old_interval),fresh_interval_rows=sum(new_interval),
    interval_availability_changes=sum(old_interval!=new_interval),
    lower_availability_changes=sum(is.finite(x$lower_reference)!=is.finite(x$lower_fresh)),
    upper_availability_changes=sum(is.finite(x$upper_reference)!=is.finite(x$upper_fresh)),
    maximum_prediction_difference_years=maximum_difference(old,new),
    reference_P95_minus_P5=d0,fresh_P95_minus_P5=d1,
    contrast_difference_years=d1-d0,contrast_review_threshold_years=threshold,
    contrast_magnitude_flag=abs(d1-d0)>threshold,contrast_sign_changed=d0*d1<0,
    maximum_lower_difference_years=maximum_difference(x$lower_reference,x$lower_fresh),
    maximum_upper_difference_years=maximum_difference(x$upper_reference,x$upper_fresh))
}))
write.csv(report,file.path(out,'spline_comparison.csv'),row.names=FALSE,na='')
pdf(file.path(out,'spline_comparison.pdf'),width=12,height=13)
for(sc in unique(z$scenario)) {
  part <- z[z$scenario==sc,]
  par(mfrow=c(5,3),mar=c(3,3,2,1),oma=c(3,2,3,0))
  for(mi in 1:5)for(tau in c(.1,.25,.5)) {
    x <- part[part$imputation==mi & abs(part$tau-tau)<1e-10,]
    x <- x[order(x$concentration),]
    yy <- c(x$estimate_reference,x$estimate_fresh);yy <- yy[is.finite(yy)]
    limits <- if(length(yy))range(yy)else c(0,1)
    if(diff(limits)==0)limits <- limits+c(-.1,.1)
    plot(NA,xlim=range(x$concentration),ylim=limits,xlab='',ylab='',
      main=paste('Imputation',mi,'; tau =',tau))
    if(any(is.finite(x$estimate_reference)))lines(x$concentration,x$estimate_reference,col='grey50',lwd=2)
    if(any(is.finite(x$estimate_fresh)))lines(x$concentration,x$estimate_fresh,col='#2166AC',lwd=1.5,lty=2)
    if(!length(yy))text(mean(range(x$concentration)),.5,'Unavailable in both',cex=.8)
  }
  mtext(paste(sc,'- manuscript reference (grey) and fresh run (blue dashed)'),side=3,outer=TRUE,line=1)
  mtext('Concentration',side=1,outer=TRUE,line=1)
  mtext('Predicted conditional survival time (years)',side=2,outer=TRUE,line=.5)
}
dev.off()
cat('Compared',nrow(report),'imputation/quantile curves;',
  sum(report$reference_available!=report$fresh_available),'availability changes;',
  sum(report$interval_availability_changes),'grid-row interval availability changes;',
  sum(report$contrast_magnitude_flag,na.rm=TRUE),'contrast magnitude flags.\n')
