# Prepare and fit Peng-Huang censored quantile models.
# Native-call conventions follow quantreg (GPL >= 2); see README.md.

cqr_grid <- function(n, taus=c(.10, .25, .50), upper=NULL) {
  stopifnot(length(n)==1L, is.finite(n), n>2,
            all(is.finite(taus)), all(taus>0 & taus<1))
  step <- min(.01, 1/(2*n^.7))
  if (is.null(upper)) upper <- min(1-1/n, max(taus)+max(.01, 4*step))
  g <- seq(1/n, 1-1/n, by=step)
  g <- g[g<=upper]
  if (length(g)<4L || max(g)<=max(taus)+2*step)
    stop("Grid must extend several steps beyond the largest requested tau")
  g
}

cqr_prepare <- function(formula, data, grid=NULL,
                                taus=c(.10, .25, .50), standardize=TRUE) {
  if (!requireNamespace("quantreg", quietly=TRUE)) stop("quantreg is required")
  if (!requireNamespace("survival", quietly=TRUE)) stop("survival is required")
  taus <- sort(unique(as.numeric(taus)))
  if (!length(taus) || any(!is.finite(taus)) || any(taus<=0 | taus>=1))
    stop("taus must lie strictly between zero and one")
  mf <- stats::model.frame(formula, data=as.data.frame(data), na.action=stats::na.fail,
                           drop.unused.levels=TRUE)
  X <- stats::model.matrix(attr(mf,"terms"), mf)
  Y <- stats::model.response(mf)
  if (!inherits(Y,"Surv") || attr(Y,"type")!="right")
    stop("A right-censored survival::Surv response is required")
  if (any(!is.finite(X)) || any(!is.finite(Y)) || any(!Y[,2] %in% c(0,1)))
    stop("Nonfinite model inputs or an invalid event indicator")
  if (is.null(grid)) grid <- cqr_grid(nrow(X), taus)
  if (!is.numeric(grid) || length(grid)<4 || any(!is.finite(grid)) ||
      any(diff(grid)<=0) || min(grid)<=0 || max(grid)>=1 ||
      min(grid)>=min(taus) || max(grid)<=max(taus))
    stop("An increasing explicit grid strictly inside (0,1) must bracket all taus")
  center <- setNames(rep(0,ncol(X)),colnames(X))
  scale <- setNames(rep(1,ncol(X)),colnames(X))
  int <- match("(Intercept)", colnames(X))
  if (standardize) {
    # Center only when the formula includes an intercept, preserving its span.
    nonconstant <- vapply(seq_len(ncol(X)), function(j) diff(range(X[,j]))>0, logical(1))
    if (!is.na(int)) center[nonconstant] <- colMeans(X[,nonconstant,drop=FALSE])
    scale[nonconstant] <- apply(X[,nonconstant,drop=FALSE],2,stats::sd)
    if (any(!is.finite(scale)) || any(scale<=0)) stop("Invalid matrix scales")
  }
  Xs <- sweep(sweep(X,2,center,"-"),2,scale,"/")
  if (!is.na(int)) Xs[,int] <- X[,int]
  structure(list(X=Xs, y=as.numeric(Y[,1]), event=as.numeric(Y[,2]),
    coefficient_names=colnames(X), center=center, scale=scale, intercept=int,
    taus=taus, grid=grid, n=nrow(X), p=ncol(X), standardize=standardize,
    formula=formula, contrasts=attr(X,"contrasts"),
    matrix_rank=qr(Xs,tol=1e-10)$rank,
    condition_number=tryCatch(kappa(Xs,exact=FALSE),error=function(e)NA_real_),
    solver="Peng-Huang"), class="cqr_prepared")
}

ph_native_unchecked <- function(x,y,cen,weights,grid) {
  # Source-equivalent call, with information retained rather than discarded.
  loadNamespace("quantreg")
  p <- ncol(x); n <- length(y); m <- length(grid)
  x <- x*weights; y <- y*weights
  cc <- as.logical(cen); y1 <- y[cc]; x1 <- x[cc,,drop=FALSE]
  n1 <- length(y1)
  z <- .Fortran("crqfnb", as.integer(n), as.integer(p),
    a1=as.double(t(x1)), c1=as.double(-y1), n1=as.integer(n1),
    as.double(x), as.double(y), as.double(cen), B=double(p*m),
    g=as.double(grid), m=as.integer(m), rep(1,p), double(n),
    rep(1,n), rep(1,n), wn=double(n1*9), wp=double((p+3)*p),
    info=integer(1), PACKAGE="quantreg")
  J <- max(0L,min(m, z$m-1L))
  B <- -matrix(z$B,p,m)[,seq_len(J),drop=FALSE]
  rownames(B) <- colnames(x)
  list(beta=B, tau=grid[seq_len(J)], native_info=z$info, native_m=z$m)
}

cqr_backtransform <- function(prepared, beta) {
  beta <- as.matrix(beta)
  ans <- sweep(beta,1,prepared$scale,"/")
  if (!is.na(prepared$intercept)) {
    j <- prepared$intercept
    ans[j,] <- beta[j,] - colSums(sweep(ans,1,prepared$center,"*"))
  }
  rownames(ans) <- prepared$coefficient_names
  ans
}

cqr_fit <- function(prepared, weights, coefficient_names=NULL,
                            weight_constant=NULL, keep_path=FALSE) {
  stopifnot(inherits(prepared,"cqr_prepared"))
  started <- proc.time()[["elapsed"]]
  weights <- as.numeric(weights)
  if (is.null(coefficient_names)) coefficient_names <- prepared$coefficient_names
  if (any(!coefficient_names %in% prepared$coefficient_names))
    stop("Unknown requested coefficient name")
  cf <- matrix(NA_real_,length(coefficient_names),length(prepared$taus),
    dimnames=list(coefficient_names,sprintf("tau=%.2f",prepared$taus)))
  dg <- list(
    n_original=prepared$n,n_positive=0L,n_zero=0L,p=prepared$p,rank=NA_integer_,
    event_rows=NA_integer_,weight_constant=NA_real_,condition_number=NA_real_,
    standardize=prepared$standardize,grid_n=length(prepared$grid),
    grid_min=min(prepared$grid),grid_max=max(prepared$grid),
    native_info=NA_integer_,native_m=NA_integer_,tau_max=NA_real_,
    last_finite_tau=NA_real_,finite_prefix_tau=NA_real_,first_nonfinite_tau=NA_real_,
    failure_class="not_run",all_targets_finite=FALSE,elapsed=NA_real_,error="",warnings="")
  finish <- function(path=NULL) {
    dg$elapsed <- proc.time()[["elapsed"]]-started
    dg$all_targets_finite <- all(is.finite(cf))
    ans <- list(coefficients=cf,diagnostics=dg)
    if (keep_path) ans$path <- path
    ans
  }
  if (length(weights)!=prepared$n || any(!is.finite(weights)) || any(weights<0)) {
    dg$failure_class <- "invalid_weights"; return(finish())
  }
  take <- weights>0
  dg$n_positive <- sum(take); dg$n_zero <- sum(!take)
  if (is.null(weight_constant)) weight_constant <- mean(weights[take])
  if (length(weight_constant)!=1L || !is.finite(weight_constant) || weight_constant<=0) {
    dg$failure_class <- "invalid_weight_constant"; return(finish())
  }
  dg$weight_constant <- weight_constant
  X <- prepared$X[take,,drop=FALSE]
  y <- prepared$y[take]; event <- prepared$event[take]
  if (nrow(X)<prepared$p) {dg$failure_class <- "insufficient_positive_rows";return(finish())}
  dg$rank <- qr(X,tol=1e-10)$rank
  dg$event_rows <- sum(event==1)
  dg$condition_number <- tryCatch(kappa(X,exact=FALSE),error=function(e)NA_real_)
  if (dg$rank<prepared$p) {dg$failure_class <- "rank_deficient";return(finish())}
  if (dg$event_rows<prepared$p) {dg$failure_class <- "insufficient_event_rows";return(finish())}
  # The native LP uses uncensored-row X; check its rank before calling Fortran.
  if (qr(X[event==1,,drop=FALSE],tol=1e-10)$rank<prepared$p) {
    dg$failure_class <- "event_matrix_rank_deficient";return(finish())
  }
  warn <- character()
  z <- tryCatch(withCallingHandlers(
    ph_native_unchecked(X,y,event,weights[take]/weight_constant,prepared$grid),
    warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart("muffleWarning")}),
    error=function(e)e)
  dg$warnings <- paste(unique(warn),collapse=" | ")
  if (inherits(z,"error")) {dg$failure_class<-"native_call_error";dg$error<-conditionMessage(z);return(finish())}
  dg$native_info <- z$native_info; dg$native_m <- z$native_m
  if (!length(z$tau)) {
    dg$failure_class <- if(z$native_info>0)"native_cholesky_failure_no_solution" else "no_solution"
    return(finish(z))
  }
  beta <- cqr_backtransform(prepared,z$beta)
  finite <- colSums(!is.finite(beta))==0L
  dg$tau_max <- max(z$tau)
  if (any(finite)) dg$last_finite_tau <- max(z$tau[finite])
  first_bad <- which(!finite)[1]
  prefix_n <- if(is.na(first_bad))length(finite)else first_bad-1L
  if (!is.na(first_bad)) dg$first_nonfinite_tau <- z$tau[first_bad]
  if (prefix_n>0L) dg$finite_prefix_tau <- z$tau[prefix_n]
  if (prefix_n>=2L) {
    # Use quantreg's coefficient interpolation exactly, after excluding a
    # nonfinite suffix. Its last row is Qhat, omitted by coef.crq; dummy zeros
    # suffice because only regression coefficients are requested here.
    sol <- rbind(tau=z$tau[seq_len(prefix_n)],beta[,seq_len(prefix_n),drop=FALSE],Qhat=0)
    ff <- structure(list(sol=sol,ctype="right"),class="crq")
    extracted <- getFromNamespace("coef.crq","quantreg")(ff,taus=prepared$taus)
    cf[,] <- extracted[coefficient_names,,drop=FALSE]
  }
  dg$failure_class <- if(any(!finite)) "nonfinite_path" else if(z$native_info>0)
    "native_cholesky_failure" else if(z$native_info<0) "native_argument_error" else
    if(all(is.finite(cf))) "none" else "requested_tau_not_reached"
  finish(list(tau=z$tau,coefficients=beta,native_info=z$native_info))
}
