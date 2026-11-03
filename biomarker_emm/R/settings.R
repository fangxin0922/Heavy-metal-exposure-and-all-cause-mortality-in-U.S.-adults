# Paths and runtime checks for the biomarker analysis.
DATA_DIR <- Sys.getenv("EMM_DATA_DIR", "data/biomarker_emm")
OUTPUT_DIR <- Sys.getenv("EMM_OUTPUT_DIR", "biomarker_emm/output")
INPUT_FILE <- file.path(DATA_DIR, "analytic.rds")
PARENT_DIR <- file.path(DATA_DIR, "parent_models")
LAB_DIR <- file.path(DATA_DIR, "raw", "labs")
MARKER_DIR <- file.path(DATA_DIR, "raw", "biomarkers")
SUPPORT_DIR <- file.path(DATA_DIR, "raw", "cbc")
LIBRARY_DIR <- Sys.getenv("HM_R_LIBRARY", "")
if (nzchar(LIBRARY_DIR)) .libPaths(c(LIBRARY_DIR, .libPaths()))

# Prevent saved fits from being combined across incompatible runtimes.
runtime_identity <- function(require_linux = FALSE) {
  info <- sessionInfo()
  blas <- normalizePath(info$BLAS, mustWork = TRUE)
  lapack <- normalizePath(La_library(), mustWork = TRUE)
  if (require_linux) stopifnot(Sys.info()[["sysname"]] == "Linux",
    grepl("openblas", blas, ignore.case = TRUE), !nzchar(Sys.getenv("LD_PRELOAD")),
    Sys.getenv("OPENBLAS_NUM_THREADS") == "1", Sys.getenv("OMP_NUM_THREADS") == "1")
  packages <- c("quantreg", "survival", "SparseM", "MatrixModels", "Matrix", "stats", "methods")
  package_files <- unlist(lapply(packages, function(package) {
    path <- find.package(package)
    c(file.path(path, "DESCRIPTION"),
      list.files(file.path(path, "R"), pattern = "[.](rdb|rdx)$", full.names = TRUE),
      list.files(file.path(path, "libs"), recursive = TRUE, full.names = TRUE))
  }))
  files <- c(blas, lapack, file.path(R.home("bin"), "exec", "R"),
    list.files(R.home("lib"), pattern = "^libR([.]|blas|lapack)", full.names = TRUE),
    package_files)
  files <- unique(normalizePath(files[!dir.exists(files)], mustWork = TRUE))
  list(platform = R.version$platform, os = unname(Sys.info()[c("sysname", "release", "machine")]),
    threads = Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_CORETYPE", "LD_PRELOAD")),
    files = tools::md5sum(files))
}
