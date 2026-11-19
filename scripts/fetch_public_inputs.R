# Download native files, or copy byte-identical files from a local raw-data tree.
# Usage: Rscript scripts/fetch_public_inputs.R [--biomarkers] [local_raw_directory]
source("scripts/input_manifest.R")
library_path <- Sys.getenv("HM_R_LIBRARY", "")
if (nzchar(library_path)) .libPaths(c(library_path, .libPaths()))
stopifnot(requireNamespace("haven", quietly = TRUE))
args <- commandArgs(trailingOnly = TRUE)
include_biomarkers <- "--biomarkers" %in% args
args <- setdiff(args, "--biomarkers")
if (length(args) > 1L) stop("Supply at most one local raw-data directory.")
manifest <- public_input_manifest(include_biomarkers)
native <- if (length(args)) list.files(normalizePath(args[1], mustWork = TRUE),
  recursive = TRUE, full.names = TRUE, pattern = "[.](xpt|dat)$", ignore.case = TRUE) else character()
options(timeout = max(300, getOption("timeout")))

validate_native <- function(path, target) {
  if (grepl("[.]xpt$", target, ignore.case = TRUE)) {
    d <- haven::read_xpt(path)
    stopifnot(nrow(d) > 0, "SEQN" %in% names(d), !anyNA(d$SEQN), !anyDuplicated(d$SEQN))
  } else {
    lines <- readLines(path, warn = FALSE)
    ids <- suppressWarnings(as.numeric(trimws(substr(lines, 1, 6))))
    # Ineligible records can omit trailing blank follow-up columns.
    stopifnot(length(lines) > 1000, all(nchar(lines) >= 46), !anyNA(ids), !anyDuplicated(ids))
  }
}

fetch_one <- function(target, url) {
  if (file.exists(target)) {
    validate_native(target, target)
    return("existing validated file")
  }
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(".download_", tmpdir = dirname(target))
  on.exit(unlink(temporary), add = TRUE)
  if (length(args)) {
    candidates <- native[tolower(basename(native)) == tolower(basename(target))]
    if (!length(candidates)) stop("Native file not found: ", basename(target))
    if (length(unique(unname(tools::md5sum(candidates)))) != 1L) {
      stop("Conflicting local copies of ", basename(target), "; supply a curated raw-data directory.")
    }
    stopifnot(file.copy(candidates[1], temporary))
  } else {
    status <- utils::download.file(url, temporary, mode = "wb", method = "libcurl", quiet = TRUE)
    if (status != 0L) stop("Download failed: ", url)
  }
  validate_native(temporary, target)
  stopifnot(file.rename(temporary, target))
  if (length(args)) "local native file" else "CDC download"
}
manifest$source <- NA_character_
for (i in seq_len(nrow(manifest))) {
  manifest$source[i] <- fetch_one(manifest$target[i], manifest$url[i])
  message("Ready: ", manifest$target[i])
}
manifest$md5 <- unname(tools::md5sum(manifest$target))
destination <- file.path(Sys.getenv("HM_DATA_DIR", "data/metal_survival"), "public_input_manifest.csv")
write.csv(manifest, destination, row.names = FALSE)
message("Validated ", nrow(manifest), " native input files.")
