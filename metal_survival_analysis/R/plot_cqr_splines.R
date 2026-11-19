# Plot aggregate spline predictions in R, separately for each imputation.
root <- file.path('metal_survival_analysis','results','spline_predictions')
p <- read.csv(file.path(root,'figure2_predictions_by_imputation.csv'))
k <- read.csv(file.path(root,'figure2_knots_and_cohorts.csv'))
r <- read.csv(file.path(root,'figure2_binned_rug.csv'))
scenarios <- c('rcs_Pb','rcs_egfr60_Cd')
titles <- c('Blood lead (ug/dL)','Urinary cadmium (ug/L; eGFR >=60)')
colors <- c('#8EBAD9','#3182BD','#08306B')
finite <- c(p$estimate,p$lower,p$upper)
limits <- range(finite[is.finite(finite)])
stopifnot(all(is.finite(limits)))
draw <- function(imputations, bands=FALSE) {
  par(mfrow=c(1,2),mar=c(4.5,4.5,4,1),oma=c(4,0,2,0))
  for(j in seq_along(scenarios)) {
    sc <- scenarios[j];z <- p[p$scenario==sc & p$imputation %in% imputations,]
    knots <- k[k$scenario==sc,]
    plot(NA,xlim=range(z$concentration),ylim=limits,xlab=titles[j],
      ylab='Predicted conditional survival time (years)',
      main=paste0('N = ',unique(z$n),'; deaths = ',unique(z$deaths)))
    abline(v=knots$concentration,lty=3,col='grey75')
    rugdata <- r[r$scenario==sc & r$count>0,]
    rug(rugdata$concentration,col=adjustcolor('grey30',alpha.f=.25),ticksize=.025)
    unavailable <- character()
    for(mi in imputations)for(t in seq_along(c(.1,.25,.5))) {
      q <- z[z$imputation==mi & abs(z$tau-c(.1,.25,.5)[t])<1e-10,]
      q <- q[order(q$concentration),]
      if(!any(is.finite(q$estimate))) {
        unavailable <- c(unavailable,paste0('I',mi,': ',c(.1,.25,.5)[t]));next
      }
      if(bands && all(is.finite(c(q$lower,q$upper))))
        polygon(c(q$concentration,rev(q$concentration)),c(q$lower,rev(q$upper)),
          col=adjustcolor(colors[t],alpha.f=.15),border=NA)
      lines(q$concentration,q$estimate,col=colors[t],lty=if(bands)1 else mi,lwd=1.5)
    }
    legend('topright',legend=c('10th percentile','25th percentile','50th percentile'),
      col=colors,lty=1,lwd=2,bty='n',cex=.75)
    if(length(unavailable))mtext(paste('Unavailable:',paste(unavailable,collapse='; ')),
      side=3,line=.3,cex=.6)
  }
  mtext(if(bands)paste('Imputation',imputations)else'Five imputations (line types 1-5)',
    outer=TRUE,side=3,cex=1.1)
  mtext('Profile: observed numeric medians and categorical modes. Dotted lines: P5/P35/P65/P95 knots.',
    outer=TRUE,side=1,line=1,cex=.7)
  if(bands)mtext('Bands: pointwise 95% normal intervals from package resampling covariance; no MI pooling or survey adjustment.',
    outer=TRUE,side=1,line=2.3,cex=.65)
}
pdf(file.path('metal_survival_analysis','figures','cqr_spline_overview.pdf'),width=12,height=7)
draw(1:5);dev.off()
png(file.path('metal_survival_analysis','figures','cqr_spline_overview.png'),
  width=12,height=7,units='in',res=180)
draw(1:5);dev.off()
pdf(file.path('metal_survival_analysis','figures','cqr_splines_by_imputation.pdf'),width=12,height=7)
for(mi in 1:5)draw(mi,TRUE)
dev.off()
