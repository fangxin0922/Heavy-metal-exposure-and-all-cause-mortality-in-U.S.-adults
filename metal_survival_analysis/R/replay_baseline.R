# Check the fitting code with 20 point-only refits of frozen historical data.
# Keep these control estimates separate from the fresh analysis.
source("metal_survival_analysis/R/model_catalog.R")
jobs <- expand.grid(scenario = PRIMARY, imputation = seq_len(M), stringsAsFactors = FALSE)
writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "audit/control_session_info.txt"))
codes <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, ]
  logfile <- file.path(RUN_DIR, "logs", paste0("control_", job$scenario, "_", job$imputation, ".log"))
  system2(file.path(R.home("bin"), "Rscript"),
          c("--vanilla", "metal_survival_analysis/R/fit_cqr.R",
            job$scenario, job$imputation, 0, "control"), stdout = logfile, stderr = logfile)
}, mc.cores = 8L)
stopifnot(all(unlist(codes) == 0L))
writeLines(c("20 control point fits complete", as.character(Sys.time())),
           file.path(RUN_DIR, "audit/CONTROL_COMPLETE.txt"))
