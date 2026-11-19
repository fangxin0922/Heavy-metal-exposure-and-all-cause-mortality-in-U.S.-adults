# Aggregate imputation-specific output and draw R forest plots.
kind <- commandArgs(TRUE)[1]
stopifnot(kind %in% c('nhanes3','cause'))
base <- if(kind=='nhanes3')Sys.getenv('NH3_OUTPUT_DIR','nhanes3_analysis')else Sys.getenv('CAUSE_OUTPUT_DIR','cause_specific')
scenarios <- if(kind=='nhanes3')c('lead','cadmium')else c('primary_Pb','primary_Hg','egfr60_Cd','egfr60_As')
endpoints <- if(kind=='nhanes3')'allcause'else c('cancer','cardiovascular')
paths <- unlist(lapply(endpoints,function(ep)unlist(lapply(scenarios,function(sc)
  file.path(base,'fits',if(kind=='cause')ep else '',sc,sprintf('imp%02d',1:5),'result.rds')))))
stopifnot(all(file.exists(paths)))
jobs <- lapply(paths,readRDS)
stopifnot(length(unique(vapply(jobs,`[[`,integer(1),'requested_R')))==1L)
tab <- do.call(rbind,lapply(jobs,`[[`,'table'))
stopifnot(nrow(tab)==if(kind=='nhanes3')30L else 120L)
dir.create(file.path(base,'results'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(base,'figures'),recursive=TRUE,showWarnings=FALSE)
write.csv(tab,file.path(base,'results/results_by_imputation.csv'),row.names=FALSE,na='')
if(kind=='nhanes3')tab$endpoint <- 'allcause'
for(ep in endpoints) {
  pdf(file.path(base,'figures',paste0(ep,'_cqr_by_imputation.pdf')),width=10,height=3.8*length(scenarios))
  par(mfrow=c(length(scenarios),3),mar=c(4,4,2,1),oma=c(0,0,2,0))
  for(sc in scenarios)for(tau_value in c(.1,.25,.5)) {
    z <- tab[tab$endpoint==ep & tab$scenario==sc & abs(tab$tau-tau_value)<1e-8,]
    finite <- c(0,z$estimate,z$lower,z$upper);finite <- finite[is.finite(finite)]
    limits <- range(finite);if(diff(limits)==0)limits <- limits+c(-1,1)
    plot(NA,xlim=limits,ylim=c(.5,5.5),xlab='Survival-time difference (years)',ylab='Imputation',
         main=paste(sc,'; tau =',tau_value),yaxt='n');axis(2,at=1:5);abline(v=0,lty=3,col='grey60')
    ok <- is.finite(z$lower)&is.finite(z$upper)
    segments(z$lower[ok],z$imputation[ok],z$upper[ok],z$imputation[ok],col='#22577A')
    points(z$estimate,z$imputation,pch=19,col='#22577A')
    if(any(!is.finite(z$estimate)))text(mean(limits),z$imputation[!is.finite(z$estimate)],'Unavailable',cex=.7)
  }
  heading <- c(allcause='NHANES III: all-cause mortality',cancer='Cancer mortality',
               cardiovascular='Cardiovascular mortality')[[ep]]
  mtext(paste(heading,'- CQR estimates by imputation'),side=3,outer=TRUE,line=.4,font=2)
  dev.off()
}
