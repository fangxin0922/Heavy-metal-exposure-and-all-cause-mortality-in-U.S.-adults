# Retain the sampling-weight representation produced by the original survey design.
# survey stores inverse inclusion probabilities; direct raw weights can differ at
# machine precision, which matters for the numerically sensitive CQR path.
iii_sampling_weights <- function(full, row_ids) {
  stopifnot(!anyDuplicated(full$SEQN),!anyDuplicated(row_ids),
    all(is.finite(full$WTPFEX6)),all(full$WTPFEX6>0))
  design <- survey::svydesign(ids=~SDPPSU6,strata=~SDPSTRA6,weights=~WTPFEX6,
    data=full,nest=TRUE)
  index <- match(row_ids,full$SEQN)
  stopifnot(!anyNA(index))
  weights <- as.numeric(stats::weights(design))[index]
  stopifnot(all(is.finite(weights)),all(weights>0),
    isTRUE(all.equal(weights,as.numeric(full$WTPFEX6[index]),tolerance=1e-12)))
  weights
}
