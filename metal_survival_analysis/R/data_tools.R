# Input readers and numerical checks.
raw_directory <- file.path(INPUT_DIR, "raw")
inherited_file <- file.path(INPUT_DIR, "analytic.rds")
metal_names <- c("LBXBPB", "LBXTHG", "URXUCD", "VNSUMARS")
data_file <- function(name) file.path(RUN_DIR, "data", name)
diagnostic_file <- function(name) file.path(RUN_DIR, "diagnostics", name)

read_native <- function(path) {
  stopifnot(file.exists(path))
  x <- as.data.frame(haven::read_xpt(path))
  x[] <- lapply(x, function(z) if (is.numeric(z)) as.numeric(z) else z)
  stopifnot("SEQN" %in% names(x), !anyDuplicated(x$SEQN), !anyNA(x$SEQN))
  x
}

read_inherited <- function() {
  x <- readRDS(inherited_file)
  stopifnot(!anyDuplicated(x$SEQN), !anyNA(x$SEQN))
  x
}

safe_log2 <- function(x) {
  answer <- rep(NA_real_, length(x))
  observed <- is.finite(x) & x > 0
  answer[observed] <- log2(x[observed])
  answer
}

compare_column <- function(a, b, variable, tolerance=1e-10) {
  stopifnot(length(a) == length(b))
  missing_mismatch <- sum(is.na(a) != is.na(b))
  observed <- !is.na(a) & !is.na(b)
  if (is.numeric(a) && is.numeric(b)) {
    difference <- abs(as.numeric(a[observed]) - as.numeric(b[observed]))
    largest <- if (length(difference)) max(difference) else 0
    n_different <- sum(difference > tolerance)
  } else {
    largest <- NA_real_
    n_different <- sum(as.character(a[observed]) != as.character(b[observed]))
  }
  data.frame(variable=variable, n=length(a), missing_mismatch=missing_mismatch,
             unequal_values=n_different, max_absolute_difference=largest,
             exact_value_identity=identical(as.vector(a), as.vector(b)),
             factor_levels_identity=if (is.factor(a) || is.factor(b))
               identical(levels(a), levels(b)) else TRUE,
             pass=missing_mismatch == 0 && n_different == 0)
}
