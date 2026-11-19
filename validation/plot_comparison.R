# Usage: Rscript validation/plot_comparison.R /path/to/cell_comparison.csv
args <- commandArgs(TRUE)
stopifnot(length(args)==1L,file.exists(args[1]))
d <- read.csv(args[1],stringsAsFactors=FALSE)
path <- file.path(dirname(args[1]),'point_estimate_comparison.pdf')
pdf(path,width=10,height=8)
par(mfrow=c(2,2),mar=c(4,4,3,1),oma=c(2,0,1,0))
labels <- c(main='Main and sensitivity CQR',nhanes3='NHANES III',cause='Cause-specific CQR',biomarker='Biomarker interaction contrasts')
for(family in names(labels)) {
  z <- d[d$family==family & is.finite(d$reference) & is.finite(d$rerun),]
  # Spline basis coefficients are not exposure contrasts; compare their curves separately.
  if(family=='main')z <- z[!grepl('^rcs_',z$cell),]
  if(family=='biomarker')z <- z[z$measure=='ic_years',]
  if(!nrow(z)) {plot.new();title(labels[[family]]);text(.5,.5,'No paired finite estimates');next}
  limits <- range(c(z$reference,z$rerun));padding <- max(diff(limits)*.05,.02)
  limits <- limits+c(-padding,padding)
  flag <- !is.na(z$magnitude_flag)&tolower(as.character(z$magnitude_flag))=='true'
  plot(z$reference,z$rerun,xlim=limits,ylim=limits,xlab='Manuscript reference (years)',
       ylab='Fresh estimate (years)',main=labels[[family]],pch=19,cex=.7,
       col=ifelse(flag,'#D55E00','#0072B2'))
  abline(a=0,b=1,lty=2,col='grey45')
  legend('topleft',legend=sprintf('%d paired estimates; %d magnitude flags',nrow(z),sum(flag)),bty='n',cex=.8)
}
mtext('Each point is a separate imputation or interaction contrast. Unavailable estimates and interval changes are reported in the CSV.',
      side=1,outer=TRUE,line=.5,cex=.75)
dev.off()
cat('Saved',path,'\n')
