# Recover missing targets only when fixed numerical retries agree.

retry_policy <- function() list(
  multipliers=c(.1,10,100), tolerance=.01, min_finite=2L,
  selection="first_finite_in_fixed_order",
  anchor_gate="all_attempts_all_base_finite_targets")

retry_apply <- function(base, attempts, eligible_taus,
                                policy=retry_policy()) {
  stopifnot(identical(policy,retry_policy()),is.matrix(base),
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

