# Separately versioned, mathematically equivalent numerical retries.
# The frozen revised_cqr_solver.R and base checkpoints remain untouched.
revised_retry_policy <- function() list(
  multipliers=c(.1,10,100), tolerance=.01, min_finite=2L,
  selection="first_finite_in_fixed_order",
  anchor_gate="all_attempts_all_base_finite_targets")

revised_retry_fingerprint <- function(object) {
  tf <- tempfile(); on.exit(unlink(tf))
  saveRDS(object,tf,version=2)
  unname(tools::md5sum(tf))
}

revised_retry_atomic_save <- function(object,path) {
  temporary <- paste0(path,".tmp-",Sys.getpid())
  saveRDS(object,temporary)
  stopifnot(file.rename(temporary,path))
}

revised_retry_apply <- function(base, attempts, eligible_taus,
                                policy=revised_retry_policy()) {
  stopifnot(identical(policy,revised_retry_policy()),is.matrix(base),
    !is.null(rownames(base)),!is.null(colnames(base)),
    identical(names(attempts),paste0("mult",policy$multipliers)))
  taus <- as.numeric(sub("^tau=","",colnames(base)))
  stopifnot(all(is.finite(taus)),all(eligible_taus %in% taus))
  alternate <- lapply(attempts,function(a) {
    stopifnot(is.matrix(a$coefficients),identical(dim(a$coefficients),dim(base)),
      identical(dimnames(a$coefficients),dimnames(base)))
    a$coefficients
  })
  # An anchor conflict anywhere blocks every fill in this case. An attempt
  # cannot be removed selectively to make the remaining candidates agree.
  anchor_differences <- unlist(lapply(alternate,function(a) {
    overlap <- is.finite(base)&is.finite(a)
    abs(a[overlap]-base[overlap])
  }),use.names=FALSE)
  anchor_max <- if(length(anchor_differences)) max(anchor_differences) else NA_real_
  anchors_pass <- !length(anchor_differences) || anchor_max<=policy$tolerance
  effective <- base; decisions <- list()
  for(j in seq_len(ncol(base))) for(i in seq_len(nrow(base))) {
    values <- vapply(alternate,function(a)a[i,j],numeric(1))
    finite <- which(is.finite(values)); n_finite <- length(finite)
    lo <- if(n_finite)min(values[finite])else NA_real_
    hi <- if(n_finite)max(values[finite])else NA_real_
    span <- hi-lo; eligible <- taus[j] %in% eligible_taus
    accepted <- FALSE; chosen <- NA_real_
    reason <- if(is.finite(base[i,j])) "base_finite_preserved" else if(!eligible)
      "target_ineligible" else if(!anchors_pass) "base_anchor_conflict" else if(
        n_finite<policy$min_finite) "fewer_than_two_finite_alternatives" else if(
          span>policy$tolerance) "alternative_range_exceeds_tolerance" else "recovered"
    if(identical(reason,"recovered")) {
      chosen <- policy$multipliers[finite[1L]]
      effective[i,j] <- values[finite[1L]]; accepted <- TRUE
    }
    decisions[[length(decisions)+1L]] <- data.frame(
      term=rownames(base)[i],tau=taus[j],base_value=base[i,j],eligible=eligible,
      n_finite_alternatives=n_finite,alternative_min=lo,alternative_max=hi,
      alternative_range=span,anchor_max_abs_difference=anchor_max,
      anchors_pass=anchors_pass,accepted=accepted,selected_multiplier=chosen,
      effective_value=effective[i,j],reason=reason,stringsAsFactors=FALSE)
  }
  list(effective_coefficients=effective,policy_decisions=do.call(rbind,decisions))
}

revised_retry_validate_manifest <- function(dest) {
  manifest <- readRDS(file.path(dest,"retry_manifest.rds"))
  core <- manifest; core$policy_fingerprint <- NULL
  base_path <- file.path(dest,"runner_manifest.rds")
  base <- readRDS(base_path)
  stopifnot(identical(manifest$schema_version,"revised-retry-1.0"),
    identical(manifest$policy,revised_retry_policy()),
    identical(manifest$policy_fingerprint,revised_retry_fingerprint(core)),
    identical(unname(tools::md5sum(base_path)),manifest$base_manifest_md5),
    identical(base$input_fingerprint,manifest$input_fingerprint),
    identical(tools::md5sum(names(manifest$source_hashes)),manifest$source_hashes),
    identical(tools::md5sum(names(base$input_hashes)),base$input_hashes))
  manifest
}

revised_retry_validate_checkpoint <- function(path,base_path,manifest) {
  cp <- readRDS(path); base <- readRDS(base_path)
  stopifnot(identical(cp$schema_version,"revised-retry-1.0"),
    identical(cp$scenario,manifest$scenario),
    identical(cp$input_fingerprint,manifest$input_fingerprint),
    identical(cp$policy_fingerprint,manifest$policy_fingerprint),
    identical(cp$policy,manifest$policy),
    identical(cp$base_checkpoint_md5,unname(tools::md5sum(base_path))),
    identical(cp$original_coefficients,base$result$coefficients),
    identical(cp$imputation,base$imputation),identical(cp$replicate,base$replicate))
  recomputed <- revised_retry_apply(cp$original_coefficients,cp$attempts,
                                   cp$eligible_taus,cp$policy)
  stopifnot(identical(cp$effective_coefficients,recomputed$effective_coefficients),
            identical(cp$policy_decisions,recomputed$policy_decisions))
  cp
}

revised_retry_original_evidence <- function(dest,manifest,require_retries=TRUE) {
  values <- list(); evidence <- character()
  for(m in seq_len(manifest$M)) {
    stem <- sprintf("imp%02d_rep000.rds",m)
    bp <- file.path(dest,"fits",stem)
    stopifnot(file.exists(bp))
    base <- readRDS(bp)
    stopifnot(identical(base$input_fingerprint,manifest$input_fingerprint),
      base$imputation==m,base$replicate==0L,base$scenario==manifest$scenario)
    rp <- file.path(dest,"retries",stem)
    values[[m]] <- base$result$coefficients
    evidence <- c(evidence,tools::md5sum(bp))
    if(file.exists(rp)) {
      cp <- revised_retry_validate_checkpoint(rp,bp,manifest)
      stopifnot(identical(cp$eligible_taus,manifest$taus))
      values[[m]] <- cp$effective_coefficients
      evidence <- c(evidence,tools::md5sum(rp))
    } else if(require_retries && any(!is.finite(values[[m]]))) {
      stop("Original numerical retries have not finished: ",manifest$scenario," MI ",m)
    }
  }
  okay <- Reduce(`&`,lapply(values,function(a)apply(is.finite(a),2,all)))
  list(eligible_taus=manifest$taus[okay],original_evidence_hashes=evidence,
       effective_originals=values)
}

revised_retry_job <- function(project_root,scenario,imputation,replicates) {
  setwd(project_root)
  .libPaths(c(normalizePath("replication_10_24/library"),.libPaths()))
  suppressPackageStartupMessages({library(quantreg);library(survival)})
  source("replication_through_2014/R/revised_cqr_solver.R")
  dest <- file.path("replication_through_2014/models",scenario)
  manifest <- revised_retry_validate_manifest(dest)
  spec <- readRDS(file.path(dest,"specification.rds"))
  imp <- readRDS(file.path(dest,"imputations.rds"))
  rp <- readRDS(file.path(dest,"replicate_design.rds"))
  info <- sessionInfo()
  stopifnot(grepl("vecLib|Accelerate",info$BLAS),
    identical(Sys.getenv("VECLIB_MAXIMUM_THREADS"),"1"),
    identical(as.numeric(imp$row_ids),as.numeric(rp$row_ids)),
    identical(as.numeric(imp$data[[imputation]]$SEQN),as.numeric(rp$row_ids)))
  prepared <- revised_cqr_prepare(spec$formula,imp$data[[imputation]],
    taus=manifest$taus,standardize=TRUE)
  stopifnot(prepared$n==length(rp$row_ids))
  constant <- mean(rp$sampling_weights)
  original <- if(any(replicates>0L))revised_retry_original_evidence(dest,manifest)else NULL
  for(b in replicates) {
    stem <- sprintf("imp%02d_rep%03d.rds",imputation,b)
    base_path <- file.path(dest,"fits",stem)
    target <- file.path(dest,"retries",stem)
    base <- readRDS(base_path)
    stopifnot(base$scenario==scenario,base$imputation==imputation,base$replicate==b,
      identical(base$input_fingerprint,manifest$input_fingerprint))
    base_hash <- unname(tools::md5sum(base_path))
    eligible <- if(b==0L)manifest$taus else original$eligible_taus
    evidence <- if(b==0L)character() else original$original_evidence_hashes
    if(file.exists(target)) {
      cp <- revised_retry_validate_checkpoint(target,base_path,manifest)
      stopifnot(identical(cp$eligible_taus,eligible),
        identical(cp$original_evidence_hashes,evidence))
      next
    }
    needed <- colnames(base$result$coefficients) %in% sprintf("tau=%.2f",eligible)
    if(!any(!is.finite(base$result$coefficients[,needed,drop=FALSE])))next
    w <- if(b==0L)rp$sampling_weights else rp$analysis_weights[,b]
    attempts <- list()
    for(mult in manifest$policy$multipliers) {
      label <- paste0("mult",mult)
      attempt_path <- file.path(dest,"retries","attempts",
        sprintf("imp%02d_rep%03d_%s.rds",imputation,b,label))
      if(file.exists(attempt_path)) {
        at <- readRDS(attempt_path)
        stopifnot(identical(at$base_checkpoint_md5,base_hash),
          identical(at$policy_fingerprint,manifest$policy_fingerprint),
          identical(at$multiplier,mult),at$scenario==scenario,
          at$imputation==imputation,at$replicate==b)
      } else {
        message(Sys.time()," retry ",scenario," MI ",imputation," replicate ",b,
                " weight multiplier ",mult)
        started <- Sys.time()
        result <- tryCatch(revised_cqr_fit(prepared,w,coefficient_names=spec$term,
          weight_constant=constant/mult,keep_path=(b==0L)),error=function(e) {
            empty <- base$result$coefficients;empty[] <- NA_real_
            list(coefficients=empty,diagnostics=list(failure_class="retry_R_error",
              error=conditionMessage(e),elapsed=as.numeric(difftime(Sys.time(),started,units="secs"))))
          })
        at <- list(scenario=scenario,imputation=imputation,replicate=b,
          input_fingerprint=manifest$input_fingerprint,base_checkpoint_md5=base_hash,
          policy_fingerprint=manifest$policy_fingerprint,multiplier=mult,
          result=result,completed_at=as.character(Sys.time()),backend=info$BLAS)
        revised_retry_atomic_save(at,attempt_path)
      }
      attempts[[label]] <- at$result
    }
    evaluated <- revised_retry_apply(base$result$coefficients,attempts,eligible,manifest$policy)
    stopifnot(identical(unname(tools::md5sum(base_path)),base_hash))
    cp <- list(schema_version="revised-retry-1.0",scenario=scenario,
      imputation=imputation,replicate=b,input_fingerprint=manifest$input_fingerprint,
      base_checkpoint_md5=base_hash,policy_fingerprint=manifest$policy_fingerprint,
      policy=manifest$policy,eligible_taus=eligible,original_evidence_hashes=evidence,
      original_coefficients=base$result$coefficients,
      effective_coefficients=evaluated$effective_coefficients,
      attempts=attempts,policy_decisions=evaluated$policy_decisions,
      completed_at=as.character(Sys.time()),backend=info$BLAS)
    revised_retry_atomic_save(cp,target)
    message("Saved ",scenario,"/",stem," recovered ",sum(cp$policy_decisions$accepted))
  }
  list(scenario=scenario,imputation=imputation,replicates=replicates,
    session=info,finished=as.character(Sys.time()))
}
