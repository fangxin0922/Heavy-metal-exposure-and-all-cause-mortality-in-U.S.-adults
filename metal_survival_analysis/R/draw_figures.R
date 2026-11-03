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

# Show each CQR estimate and its imputation-specific interval.
cqr_path <- file.path(RUN_DIR, "results", "results_by_imputation.csv")
if (file.exists(cqr_path)) {
  cqr <- read.csv(cqr_path)
  cqr <- cqr[cqr$scenario %in% PRIMARY, ]
  for (tau in TAUS) {
    save_figure(paste0("cqr_quantile_", format(tau, trim = TRUE)), function() {
      par(mfrow = c(2, 2), mar = c(4.2, 5, 3, 1), oma = c(2, 0, 2, 0))
      for (scenario in PRIMARY) {
        part <- cqr[cqr$scenario == scenario & abs(cqr$tau - tau) < 1e-10, ]
        limits <- c(0, part$estimate, part$lower, part$upper)
        limits <- range(limits[is.finite(limits)])
        limits <- limits + c(-1, 1) * max(diff(limits) * 0.08, 0.05)
        plot(NA, NA, xlim = limits, ylim = c(0.5, 5.5), yaxt = "n",
             xlab = "Conditional survival-time difference (years)", ylab = "",
             main = labels[[scenario]])
        axis(2, at = 1:5, labels = paste("MI", 5:1), las = 1)
        abline(v = 0, col = "grey65", lty = 2)
        y <- 6 - part$imputation
        interval <- is.finite(part$lower) & is.finite(part$upper)
        segments(part$lower[interval], y[interval], part$upper[interval], y[interval],
                 col = colors[[scenario]], lwd = 2)
        observed <- is.finite(part$estimate)
        points(part$estimate[observed], y[observed], col = colors[[scenario]], pch = 19)
        if (any(!observed)) text(mean(limits), y[!observed], "Unavailable", cex = 0.75, col = "grey45")
      }
      mtext(paste("Censored quantile regression: quantile", tau), outer = TRUE, side = 3)
      mtext("Per exposure doubling; package intervals are not MI-pooled or survey-design-adjusted",
            outer = TRUE, side = 1, cex = 0.8)
    })
  }
}
message("Figures saved in ", file.path(RUN_DIR, "figures"))
