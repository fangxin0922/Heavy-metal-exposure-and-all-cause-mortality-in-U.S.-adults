source('nhanes3_analysis/R/settings.R')
source('nhanes3_analysis/R/sampling_weights.R')
metal <- commandArgs(TRUE)[1]
stopifnot(metal %in% scenarios)
dest <- file.path(base,'models',metal)
stopifnot(!file.exists(file.path(dest,'imputations.rds')))
d <- readRDS(file.path(dest,'cohort.rds'))
covariates <- c('age','age2','sex','race','education','marital','pir','bmi',
  'packyears','alcohol','activity','phase','cholesterol',if(metal=='cadmium')'egfr')
formula <- as.formula(paste('Surv(time,event) ~ exposure +',paste(covariates,collapse=' + ')))
z <- d[c('time','event','exposure',covariates)]
z$log_weight <- log(d$WTPFEX6)
z$stratum_aux <- factor(d$SDPSTRA6)
z$psu_aux <- factor(d$SDPPSU6)
method <- make.method(z);method[] <- ''
for(v in names(z)[colSums(is.na(z))>0])
  method[v] <- if(is.factor(z[[v]])) if(nlevels(z[[v]])==2)'logreg'else'polyreg' else 'pmm'
predictors <- make.predictorMatrix(z);predictors[method=='',] <- 0
seed <- if(metal=='lead')2026102904L else 2026102905L
if(metal=='cadmium')stopifnot(all(is.finite(z$egfr)),all(z$egfr>=60),method['egfr']=='')
mi <- mice(z,m=5,maxit=50,method=method,predictorMatrix=predictors,seed=seed,printFlag=TRUE)
completed <- lapply(1:5,function(i) {
  x <- complete(mi,i);x$SEQN <- d$SEQN
  stopifnot(!anyNA(x))
  for(v in names(z)) {observed <- !is.na(z[[v]]);stopifnot(identical(x[[v]][observed],z[[v]][observed]))}
  x
})
saveRDS(z,file.path(dest,'preimputation_frame.rds'))
saveRDS(mi,file.path(dest,'mids.rds'))
saveRDS(list(data=completed,seed=seed,iterations=50L,method=method,predictor_matrix=predictors),
        file.path(dest,'imputations.rds'))
sampling_weights <- iii_sampling_weights(readRDS(file.path(base,'data/full_mec.rds')),d$SEQN)
saveRDS(list(scenario=metal,formula=formula,row_ids=d$SEQN,term='exposure',
  n=nrow(d),events=sum(d$event),sampling_weights=sampling_weights,
  weight_source='Full MEC survey-design sampling weights',period='1988-1994',M=5L,seed=seed),
  file.path(dest,'specification.rds'))
write.csv(data.frame(variable=names(z),n_missing=colSums(is.na(z))),
          file.path(dest,'missingness_before_mi.csv'),row.names=FALSE)
pdf(file.path(dest,'mice_trace.pdf'),width=11,height=8.5);print(plot(mi));dev.off()
writeLines(capture.output(sessionInfo()),file.path(dest,'session_info.txt'))
