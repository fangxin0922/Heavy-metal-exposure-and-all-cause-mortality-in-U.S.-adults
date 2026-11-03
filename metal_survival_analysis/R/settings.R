# Shared paths, packages and analysis settings.
PROJECT_ROOT <- normalizePath(".", mustWork = TRUE)
RUN_DIR <- file.path(PROJECT_ROOT, "metal_survival_analysis")
INPUT_DIR <- Sys.getenv("HM_DATA_DIR", unset = file.path(PROJECT_ROOT, "data/metal_survival"))
.libPaths(c(file.path(RUN_DIR, "library"),
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
for (folder in c("data", "models", "fits", "results", "logs", "diagnostics", "figures")) {
  dir.create(file.path(RUN_DIR, folder), recursive = TRUE, showWarnings = FALSE)
}
options(stringsAsFactors = FALSE, survey.lonely.psu = "fail")
M <- 5L
MAXIT <- 20L
BOOTSTRAPS <- 400L
TAUS <- c(0.10, 0.25, 0.50)
MI_SEED <- 4242L

save_result <- function(object, path) {
  temporary <- paste0(path, ".partial_", Sys.getpid())
  saveRDS(object, temporary)
  stopifnot(file.rename(temporary, path))
}

write_table <- function(object, path) {
  write.csv(object, path, row.names = FALSE, na = "")
}
