# Shared paths, packages and analysis settings.
PROJECT_ROOT <- normalizePath(".", mustWork = TRUE)
RUN_DIR <- file.path(PROJECT_ROOT, "metal_survival_analysis")
.libPaths(c(file.path(RUN_DIR, "library"),
            file.path(PROJECT_ROOT, "replication_10_24/library"),
            Sys.getenv("HM_R_LIBRARY", unset = ""), .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(haven)
  library(survival)
  library(survey)
  library(mice)
  library(quantreg)
})
for (folder in c("data", "models", "results", "logs", "audit", "figures")) {
  dir.create(file.path(RUN_DIR, folder), recursive = TRUE, showWarnings = FALSE)
}
options(stringsAsFactors = FALSE, survey.lonely.psu = "fail")
M <- 5L
MAXIT <- 20L
BOOTSTRAPS <- 400L
TAUS <- c(0.10, 0.25, 0.50)
MI_SEED <- 4242L
REFERENCE_COMMIT <- "da34a8bdeff9d665da38f7f07c853aaaf71bc5c3"

save_result <- function(object, path) {
  temporary <- paste0(path, ".partial_", Sys.getpid())
  saveRDS(object, temporary)
  stopifnot(file.rename(temporary, path))
}

write_table <- function(object, path) {
  write.csv(object, path, row.names = FALSE, na = "")
}

record_sources <- function(paths, destination) {
  stopifnot(all(file.exists(paths)))
  hashes <- tools::md5sum(paths)
  write_table(data.frame(path = names(hashes), md5 = unname(hashes)), destination)
  invisible(hashes)
}
