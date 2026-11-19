# Stage the retained reference archive without changing values or factor levels.
# This does not reconstruct the unavailable original screening/harmonization code.
source("scripts/input_manifest.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript scripts/stage_reference_input.R /path/to/reference_analytic.rds")
source_file <- normalizePath(args[1], mustWork = TRUE)
target_dir <- Sys.getenv("HM_DATA_DIR", "data/metal_survival")
target <- file.path(target_dir, "analytic.rds")
validate_reference_input(readRDS(source_file))
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(target)) {
  if (!identical(unname(tools::md5sum(source_file)), unname(tools::md5sum(target)))) {
    stop("A different reference input already exists: ", target)
  }
} else {
  stopifnot(file.copy(source_file, target, overwrite = FALSE))
}
stopifnot(identical(unname(tools::md5sum(source_file)), unname(tools::md5sum(target))))
message("Validated and staged the supplied reference archive: ", target)
