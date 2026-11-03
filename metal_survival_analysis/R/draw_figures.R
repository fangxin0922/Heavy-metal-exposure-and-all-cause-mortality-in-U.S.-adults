# Plot saved CQR and Cox estimates.
source("metal_survival_analysis/R/model_catalog.R")

colors <- c(primary_Pb = "#255F85", primary_Hg = "#3E8C70",
            egfr60_Cd = "#B66137", egfr60_As = "#80569A")
labels <- c(primary_Pb = "Blood lead", primary_Hg = "Blood mercury",
            egfr60_Cd = "Urinary cadmium (eGFR >=60)",
            egfr60_As = "Urinary arsenic (eGFR >=60)")

save_figure <- function(name, draw, width = 10, height = 7) {
  for (extension in c("pdf", "png")) {
    path <- file.path(RUN_DIR, "figures", paste0(name, ".", extension))
    if (extension == "pdf") grDevices::pdf(path, width = width, height = height,
                                           useDingbats = FALSE)
    else grDevices::png(path, width = width, height = height, units = "in", res = 180)
    tryCatch(draw(), finally = grDevices::dev.off())
  }
}

cox_path <- file.path(RUN_DIR, "results", "cox_by_imputation.csv")
if (file.exists(cox_path)) {
  cox <- read.csv(cox_path)
  usable <- with(cox, status == "ok" & is.finite(HR) & is.finite(lower) &
                   is.finite(upper) & lower > 0 & upper > 0)
  cox <- cox[usable & cox$scenario %in% PRIMARY, ]
  if (nrow(cox)) save_figure("cox_forest_by_imputation", function() {
    par(mfrow = c(2, 2), mar = c(4.2, 5, 3, 1), oma = c(2, 0, 2, 0))
    limits <- range(c(1, cox$lower, cox$upper))
    limits <- exp(log(limits) + c(-1, 1) * max(diff(log(limits)) * 0.08, 0.02))
    for (scenario in PRIMARY) {
      result <- cox[cox$scenario == scenario, ]
      plot(NA, NA, xlim = limits, ylim = c(0.5, 5.5), log = "x", yaxt = "n",
           xlab = "Hazard ratio (95% CI)", ylab = "", main = labels[[scenario]])
      axis(2, at = 1:5, labels = paste("MI", 5:1), las = 1)
      abline(v = 1, col = "grey65", lty = 2)
      if (nrow(result)) {
        y <- 6 - result$imputation
        segments(result$lower, y, result$upper, y, col = colors[[scenario]], lwd = 2)
        points(result$HR, y, col = colors[[scenario]], pch = 19)
      }
    }
    mtext("Supplementary survey Cox estimates, reported separately for five imputations", outer = TRUE, side = 3)
    mtext("One-unit increase in the saved exposure term; these effects are not years lost", outer = TRUE, side = 1, cex = 0.85)
  })
}

spline_path <- file.path(RUN_DIR, "results", "cox_pspline_curves.csv")
if (file.exists(spline_path)) {
  curves <- read.csv(spline_path)
  if (nrow(curves)) {
    for (mi in sort(unique(curves$imputation))) {
      save_figure(paste0("cox_pspline_imputation_", mi), function() {
        par(mfrow = c(2, 2), mar = c(4.5, 5, 3, 1), oma = c(2, 0, 2, 0))
        for (scenario in PRIMARY) {
          curve <- curves[curves$scenario == scenario & curves$imputation == mi, ]
          curve <- curve[order(curve$exposure), ]
          finite <- with(curve, is.finite(exposure) & is.finite(lower) & is.finite(upper))
          curve <- curve[finite, ]
          if (!nrow(curve)) {
            plot.new(); title(labels[[scenario]]); text(0.5, 0.5, "Curve unavailable")
            next
          }
          plot(curve$exposure, curve$log_hazard_contribution, type = "n",
               ylim = range(c(curve$lower, curve$upper)),
               xlab = unique(curve$term), ylab = "Centered partial log hazard",
               main = labels[[scenario]])
          polygon(c(curve$exposure, rev(curve$exposure)),
                  c(curve$lower, rev(curve$upper)), border = NA,
                  col = adjustcolor(colors[[scenario]], alpha.f = 0.18))
          lines(curve$exposure, curve$log_hazard_contribution, col = colors[[scenario]], lwd = 2)
          abline(h = 0, lty = 2, col = "grey65")
        }
        mtext(paste("Exploratory unweighted Cox penalized splines - imputation", mi), outer = TRUE, side = 3)
        mtext("Shading: centered term contribution +/- 1.96 SE; no survey adjustment or MI pooling", outer = TRUE, side = 1, cex = 0.8)
      })
    }
  }
}

# Compare available primary estimates with the baseline run.
fresh_path <- file.path(RUN_DIR, "results", "results_by_imputation.csv")
previous_path <- file.path(PROJECT_ROOT, "replication_through_2014",
                           "official_crq_builtin_20261029", "results_by_imputation.csv")
if (file.exists(fresh_path) && file.exists(previous_path)) {
  fresh <- read.csv(fresh_path)
  previous <- read.csv(previous_path)
  if ("branch" %in% names(fresh)) fresh <- fresh[fresh$branch == "fresh", ]
  required <- c("scenario", "imputation", "term", "tau", "estimate")
  stopifnot(all(required %in% names(fresh)), all(required %in% names(previous)))
  comparison <- merge(previous[, required], fresh[, required],
                      by = c("scenario", "imputation", "term", "tau"),
                      suffixes = c("_previous", "_fresh"), all = TRUE)
  comparison <- comparison[comparison$scenario %in% PRIMARY, ]
  good <- with(comparison, is.finite(estimate_previous) & is.finite(estimate_fresh))
  finite_comparison <- comparison[good, ]
  if (nrow(finite_comparison)) save_figure("cqr_previous_vs_fresh", function() {
    par(mfrow = c(1, length(TAUS)), mar = c(5, 5, 3, 1), oma = c(4, 0, 2, 0))
    for (tau in TAUS) {
      part <- finite_comparison[abs(finite_comparison$tau - tau) < 1e-10, ]
      if (!nrow(part)) {
        plot.new(); title(paste("Quantile", tau)); text(0.5, 0.5, "No comparable estimates")
        next
      }
      extent <- range(c(part$estimate_previous, part$estimate_fresh), finite = TRUE)
      span <- max(diff(extent), 0.05)
      extent <- extent + c(-0.1, 0.1) * span
      plot(part$estimate_previous, part$estimate_fresh,
           xlim = extent, ylim = extent, asp = 1,
           pch = part$imputation + 14, col = colors[part$scenario],
           xlab = "Previous CQR estimate (years)", ylab = "Fresh CQR estimate (years)",
           main = paste("Quantile", tau))
      abline(a = 0, b = 1, col = "grey55", lty = 2)
    }
    mtext("CQR refactor agreement: each symbol is one metal and imputation", side = 3, outer = TRUE)
    par(xpd = NA)
    legend("bottom", inset = c(0, -0.65), legend = unname(labels),
           col = unname(colors), pch = 19, bty = "n", ncol = 2, cex = 0.8)
    mtext(paste(sum(!good), "matched or unmatched cells with unavailable estimates omitted"),
          side = 1, outer = TRUE, line = 3, cex = 0.8)
  }, width = 13, height = 5)
}

message("Available result figures were saved in ", file.path(RUN_DIR, "figures"))
