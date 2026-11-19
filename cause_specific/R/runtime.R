# The sparse-event repair is isolated from the main and NHANES III fits.
repair <- readRDS('cause_specific/library/repair_manifest.rds')
stopifnot(repair$repair=='quantreg-6.1-crqfnb-n1-initialisation',
  identical(unname(tools::md5sum(names(repair$hashes))),unname(repair$hashes)))
.libPaths(c(normalizePath('cause_specific/library'),.libPaths()))
if('quantreg' %in% loadedNamespaces())stop('Start a fresh R process for the repaired package.')
suppressPackageStartupMessages({library(quantreg);library(survival)})
stopifnot(identical(normalizePath(find.package('quantreg')),repair$package_path),
          as.character(packageVersion('quantreg'))=='6.1')
