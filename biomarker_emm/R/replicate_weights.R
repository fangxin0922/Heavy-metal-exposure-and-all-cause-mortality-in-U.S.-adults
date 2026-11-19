# Generate the original Rao-Wu full-component draws before selecting a domain.
# Adapted from the study's recorded revised_inference.R; no saved banks required.
.replicate_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

build_parent_replicates <- function(full_data, id_col = "SEQN",
    weight_col = "w_metal", domain_col = "analysis_domain",
    strata_col = "SDMVSTRA", psu_col = "SDMVPSU", R = 400L,
    seed = 20260101L, retain_design = FALSE) {
  .replicate_assert(requireNamespace("survey", quietly = TRUE), "survey is required")
  needed <- c(id_col, weight_col, domain_col, strata_col, psu_col)
  .replicate_assert(all(needed %in% names(full_data)), "Missing full-component design columns")
  .replicate_assert(length(R) == 1L && is.finite(R) && R >= 2 && R == as.integer(R), "Invalid replicate count")
  .replicate_assert(length(seed) == 1L && is.finite(seed), "A finite replicate seed is required")
  ids <- full_data[[id_col]]
  w <- full_data[[weight_col]]
  domain <- full_data[[domain_col]]
  .replicate_assert(!anyNA(ids) && !anyDuplicated(ids), "Participant IDs must be unique and nonmissing")
  .replicate_assert(is.numeric(w) && all(is.finite(w)) && all(w >= 0) && any(w > 0), "Component weights must be finite, nonnegative, and not all zero")
  .replicate_assert(is.logical(domain) && !anyNA(domain) && any(domain & w > 0), "Domain must be logical, nonmissing, and contain positive weight")
  .replicate_assert(!anyNA(full_data[c(strata_col, psu_col)]), "Stratum/PSU IDs must be nonmissing")
  dd <- data.frame(.id = ids, .weight = w, .stratum = full_data[[strata_col]],
                   .psu = full_data[[psu_col]], .domain = domain,
                   .audit_value = seq_along(ids) / length(ids))
  psu_key <- interaction(dd$.stratum, dd$.psu, drop = TRUE)
  npsu <- tapply(as.character(psu_key), dd$.stratum, function(x) length(unique(x)))
  .replicate_assert(all(npsu >= 2L), "Full component has a stratum with fewer than two PSUs; investigate before resampling")
  des <- survey::svydesign(ids = ~.psu, strata = ~.stratum, weights = ~.weight,
                           data = dd, nest = TRUE)
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed))
  rep_full <- survey::as.svrepdesign(des, type = "subbootstrap", replicates = R, mse = FALSE)
  # Replicates are generated BEFORE dropping any off-domain participant or PSU.
  rep_domain <- rep_full[domain, ]
  aw <- as.matrix(stats::weights(rep_domain, type = "analysis"))
  rownames(aw) <- as.character(ids[domain])
  colnames(aw) <- sprintf("rep%03d", seq_len(R))
  sw <- as.numeric(stats::weights(rep_domain, type = "sampling"))
  names(sw) <- as.character(ids[domain])
  .replicate_assert(identical(dim(aw), c(sum(domain), as.integer(R))), "Unexpected replicate-weight dimensions")
  .replicate_assert(all(is.finite(aw)) && all(aw >= 0), "Invalid analysis replicate weights")
  .replicate_assert(all(colSums(aw) > 0), "An entire domain has zero weight in a replicate")
  .replicate_assert(isTRUE(all.equal(sw, stats::setNames(w[domain], as.character(ids[domain])), check.attributes = TRUE)), "Sampling weights changed during design construction")
  # Independent package-level check of the complete design-to-variance chain.
  y <- dd$.audit_value[domain]
  q <- sum(sw * y) / sum(sw)
  qr <- colSums(aw * y) / colSums(aw)
  v <- survey::svrVar(qr, rep_domain$scale, rep_domain$rscales,
                      na.action = "na.fail", mse = FALSE, coef = q)
  wr <- survey::withReplicates(rep_domain,
          function(weights, data) stats::weighted.mean(data$.audit_value, weights),
          return.replicates = TRUE)
  .replicate_assert(isTRUE(all.equal(as.numeric(q), as.numeric(stats::coef(wr)), tolerance = 1e-11)) &&
                 isTRUE(all.equal(as.numeric(v), as.numeric(stats::vcov(wr)), tolerance = 1e-11)),
                 "Replicate inference does not agree with survey::withReplicates")
  traditional_df <- function(rows) length(unique(psu_key[rows])) - length(unique(dd$.stratum[rows]))
  df_domain <- as.numeric(survey::degf(rep_domain))
  .replicate_assert(is.finite(df_domain) && df_domain > 0, "Nonpositive domain survey degrees of freedom")
  out <- list(row_ids = ids[domain], full_row_ids = ids, domain_index = which(domain),
    sampling_weights = sw, analysis_weights = aw, scale = rep_domain$scale,
    rscales = rep_domain$rscales, mse = FALSE, R = as.integer(R), seed = as.integer(seed),
    rng_kind = RNGkind(), R_version = R.version.string,
    survey_version = as.character(utils::packageVersion("survey")),
    complete_data_df = df_domain, df_method = "rank of domain analysis replicate weights minus one",
    design_checks = data.frame(n_full = nrow(dd), n_domain = sum(domain),
      n_positive_full = sum(w > 0), n_positive_domain = sum(domain & w > 0),
      n_strata_full = length(unique(dd$.stratum)), n_psus_full = length(unique(psu_key)),
      n_psus_without_domain = sum(!levels(psu_key) %in% unique(as.character(psu_key[domain]))),
      full_replicate_df = as.numeric(survey::degf(rep_full)), domain_replicate_df = df_domain,
      full_psu_minus_strata_df = traditional_df(w > 0),
      domain_psu_minus_strata_df = traditional_df(domain & w > 0),
      weighted_mean_check = q, weighted_mean_variance_check = as.numeric(v),
      withReplicates_agreement = TRUE))
  if (retain_design) {
    out$full_design <- des
    out$full_rep_design <- rep_full
    out$domain_rep_design <- rep_domain
  }
  out
}
