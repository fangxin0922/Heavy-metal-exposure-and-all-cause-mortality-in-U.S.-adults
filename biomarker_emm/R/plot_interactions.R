# Plot the primary CQR interaction estimates in years.

source("biomarker_emm/R/settings.R")
# Set HM_R_LIBRARY to the ggplot2 3.5.2 library if it is outside the default search path.
# This script runs in a separate process and does not load or fit CQR models.
PLOT_LIBRARY <- Sys.getenv("HM_R_LIBRARY", unset = "")
suppressPackageStartupMessages({
  library(ggplot2, lib.loc = if (nzchar(PLOT_LIBRARY)) PLOT_LIBRARY else NULL)
  library(dplyr)
  library(readr)
})
stopifnot(as.character(packageVersion("ggplot2")) == "3.5.2")

dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
results <- read_csv(file.path(OUTPUT_DIR, "aggregate", "results.csv"),
  show_col_types = FALSE)
labels <- c(vitamin_c = "Serum vitamin C", alpha_carotene = "Serum alpha-carotene",
  rbc_folate = "RBC folate")

# Open points retain estimates whose intervals are unavailable ----

primary <- results %>% filter(tau == .10) %>%
  mutate(marker = factor(labels[biomarker], levels = unname(labels)),
    metal = factor(metal, levels = c("As", "Cd", "Hg", "Pb")),
    interval = if_else(status == "complete", "Provisional 95% CI", "CI unavailable"))
stopifnot(nrow(primary) == 12L, !anyDuplicated(primary$pair))

figure <- ggplot(primary, aes(x = ic_years, y = metal)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_segment(data = filter(primary, status == "complete"),
    aes(x = ic_lower_years, xend = ic_upper_years, yend = metal),
    linewidth = .7, colour = "#22577A") +
  geom_point(aes(shape = interval), size = 2.8, stroke = .8, colour = "#22577A") +
  scale_shape_manual(values = c("Provisional 95% CI" = 16, "CI unavailable" = 1),
    breaks = c("Provisional 95% CI", "CI unavailable")) +
  facet_wrap(~marker, nrow = 1, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = .15)) +
  labs(x = "P75 minus P25 interaction (years)", y = NULL, shape = NULL,
    title = "Biomarker effect modification: provisional CQR results",
    subtitle = "10th conditional survival-time percentile; a doubling of metal concentration",
    caption = paste(
      "Complete cases; fixed biomarker quartiles. All intervals remain provisional because solver equivalence failed.",
      "Sensitivity also affects complete analyses of arsenic with vitamin C and RBC folate.",
      "An unavailable interval is not a null finding. A positive interaction alone does not establish protection.",
      sep = "\n")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 9),
    plot.title = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "figures", "primary_interaction.pdf"),
  figure, width = 12, height = 5.3, units = "in")
ggsave(file.path(OUTPUT_DIR, "figures", "primary_interaction.png"),
  figure, width = 12, height = 5.3, units = "in", dpi = 180)
