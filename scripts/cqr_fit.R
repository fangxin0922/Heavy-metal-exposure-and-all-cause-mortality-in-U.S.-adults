# Shared extraction of the public Peng-Huang fit and imputation-specific intervals.
# The caller selects the correct quantreg library before sourcing this file.
fit_cqr_job <- function(analysis, formula, weights, term, seed, metadata, destination,
                        resamples=400L, allow_sparse=FALSE) {
  stopifnot(resamples %in% c(0L,400L),length(weights)==nrow(analysis),
    all(is.finite(weights)),all(weights>0),!anyNA(analysis),
    identical(environment(quantreg::crq),asNamespace('quantreg')))
  if(dir.exists(destination) && length(list.files(destination,all.files=TRUE,no..=TRUE)))
    stop('Use an empty fit directory: ',destination)
  dir.create(destination,recursive=TRUE,showWarnings=FALSE)
  response <- model.response(model.frame(formula,data=analysis,na.action=na.fail))
  if(!allow_sparse && nrow(analysis)>9*sum(response[,2]))stop('Sparse outcome requires the isolated memory repair.')
  analysis$.analysis_weight <- as.numeric(weights)/mean(weights)
  # summary.crq reconstructs the original model frame in the formula environment.
  # Retain only the job's data/formula there, not a transient caller's frame.
  evaluation_environment <- list2env(list(analysis=analysis),parent=.GlobalEnv)
  environment(formula) <- evaluation_environment
  evaluation_environment$formula <- formula
  warnings <- character()
  record <- function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}
  stage <- function(x) {writeLines(c(x,as.character(Sys.time())),file.path(destination,'STATUS.txt'));cat(x,'\n');flush.console()}
  taus <- c(.1,.25,.5)
  stage('Fitting point estimates')
  fit <- withCallingHandlers(quantreg::crq(formula,data=analysis,weights=.analysis_weight,
    method='PengHuang',na.action=na.fail),warning=record)
  coefficients <- as.matrix(coef(fit,taus=taus))
  stopifnot(term %in% rownames(coefficients),ncol(coefficients)==3L)
  saveRDS(list(fit=fit,coefficients=coefficients),file.path(destination,'point.rds'))
  tab <- cbind(as.data.frame(metadata,stringsAsFactors=FALSE)[rep(1,3),,drop=FALSE],
    data.frame(term=term,tau=taus,estimate=as.numeric(coefficients[term,]),lower=NA_real_,
      upper=NA_real_,std_error=NA_real_,p_value=NA_real_,package_missing_draws=NA_integer_,
      package_Brep=NA_integer_,status=ifelse(is.finite(coefficients[term,]),'point_only','percentile_unavailable')))
  available <- taus[apply(coefficients,2,function(z)any(is.finite(z)))]
  package_summary <- NULL
  set.seed(seed)
  if(resamples && length(available)) {
    stage('Running 400 package jackknife resamples')
    package_summary <- withCallingHandlers(summary(fit,taus=available,se='boot',
      bmethod='jack',R=resamples,alpha=.05),warning=record)
    for(s in package_summary) {
      k <- which(abs(tab$tau-s$tau)<1e-10)
      z <- s$coefficients[term,]
      stopifnot(length(k)==1L,isTRUE(all.equal(unname(z['Value']),tab$estimate[k],tolerance=1e-12)))
      tab[k,c('lower','upper','std_error','p_value')] <- as.list(unname(z[c('Lower Bd','Upper Bd','Std Error','Pr(>|t|)')]))
      tab$package_missing_draws[k] <- s$NAs;tab$package_Brep[k] <- s$Brep
      tab$status[k] <- if(all(is.finite(z[c('Lower Bd','Upper Bd')])))'package_CI_available'else'package_CI_unavailable'
    }
  }
  saveRDS(list(table=tab,package_summary=package_summary,seed=seed,requested_R=resamples,
    warnings=unique(warnings),package_path=find.package('quantreg')),file.path(destination,'result.rds'))
  write.csv(tab,file.path(destination,'results.csv'),row.names=FALSE,na='')
  writeLines(unique(warnings),file.path(destination,'warnings.txt'))
  writeLines(capture.output(sessionInfo()),file.path(destination,'session_info.txt'))
  stage('complete')
}
