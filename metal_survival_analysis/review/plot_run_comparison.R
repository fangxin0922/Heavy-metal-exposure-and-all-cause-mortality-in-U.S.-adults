# Plot primary estimates with legends and availability counts.
project_root <- normalizePath(".", mustWork = TRUE)
run_dir <- file.path(project_root, "metal_survival_analysis")
fresh_path <- file.path(run_dir, "results", "results_by_imputation.csv")
previous_path <- file.path(project_root, "replication_through_2014",
                           "official_crq_builtin_20261029", "results_by_imputation.csv")
primary <- c("primary_Pb", "primary_Hg", "egfr60_Cd", "egfr60_As")
taus <- c(0.10, 0.25, 0.50)
colors <- c(primary_Pb = "#255F85", primary_Hg = "#3E8C70",
            egfr60_Cd = "#B66137", egfr60_As = "#80569A")
metal_labels <- c("Lead (blood)", "Mercury (blood)",
                  "Cadmium (urine; eGFR >=60)", "Arsenic (urine; eGFR >=60)")
fresh <- read.csv(fresh_path)
previous <- read.csv(previous_path)
if ("branch" %in% names(fresh)) fresh <- fresh[fresh$branch == "fresh", ]
columns <- c("scenario", "imputation", "term", "tau", "estimate")
keys <- setdiff(columns, "estimate")
stopifnot(all(columns %in% names(fresh)), all(columns %in% names(previous)))
fresh <- fresh[fresh$scenario %in% primary, columns]
previous <- previous[previous$scenario %in% primary, columns]
stopifnot(!anyDuplicated(fresh[keys]), !anyDuplicated(previous[keys]))
comparison <- merge(previous, fresh, by = keys, all = TRUE,
                     suffixes = c("_previous", "_fresh"))
stopifnot(nrow(comparison) == 60L,
          all(comparison$imputation %in% 1:5), all(comparison$tau %in% taus))
old_available <- is.finite(comparison$estimate_previous)
new_available <- is.finite(comparison$estimate_fresh)
paired <- old_available & new_available
plotted <- comparison[paired, ]
stopifnot(nrow(plotted) > 0L)
extent <- range(c(plotted$estimate_previous, plotted$estimate_fresh))
extent <- extent + c(-0.08, 0.08) * max(diff(extent), 0.05)
missing_note <- sprintf(
  "%d/%d cells plotted. Omitted: %d previous-only, %d fresh-only, %d unavailable in both runs.",
  sum(paired), nrow(comparison), sum(old_available & !new_available),
  sum(!old_available & new_available), sum(!old_available & !new_available))

draw_comparison <- function() {
  layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE),
         heights = c(4, 1.25))
  par(oma = c(0.3, 0, 2.3, 0), mar = c(4.3, 4.4, 2.4, 1))
  for (tau in taus) {
    part <- plotted[abs(plotted$tau - tau) < 1e-10, ]
    plot(part$estimate_previous, part$estimate_fresh,
         xlim = extent, ylim = extent, asp = 1,
         pch = part$imputation + 14L, col = colors[part$scenario], cex = 1.05,
         xlab = "Previous CQR estimate (years)", ylab = "Fresh CQR estimate (years)",
         main = sprintf("Quantile %.2f", tau))
    abline(a = 0, b = 1, col = "grey55", lty = 2)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  legend(x = 0.5, y = 0.98, xjust = 0.5, yjust = 1,
         legend = metal_labels, col = unname(colors), pch = 19,
         ncol = 4, bty = "n", cex = 0.9)
  legend(x = 0.5, y = 0.66, xjust = 0.5, yjust = 1,
         legend = paste("Imputation", 1:5), col = "grey30", pch = 15:19,
         ncol = 5, bty = "n", cex = 0.85)
  text(0.5, 0.24, missing_note, cex = 0.9)
  text(0.5, 0.05,
       "Dashed line: exact agreement. Estimates are conditional survival-time differences per exposure doubling.",
       cex = 0.85)
  mtext("Primary CQR estimates: previous versus fresh run", side = 3,
        outer = TRUE, line = 0.5, font = 2, cex = 1.1)
}

for (extension in c("pdf", "png")) {
  destination <- file.path(run_dir, "figures", paste0("cqr_previous_vs_fresh_labeled.", extension))
  if (extension == "pdf") {
    grDevices::pdf(destination, width = 13, height = 6.2, useDingbats = FALSE)
  } else {
    grDevices::png(destination, width = 13, height = 6.2, units = "in", res = 180)
  }
  tryCatch(draw_comparison(), finally = grDevices::dev.off())
  message("Saved presentation supplement: ", destination)
}
message(missing_note)
