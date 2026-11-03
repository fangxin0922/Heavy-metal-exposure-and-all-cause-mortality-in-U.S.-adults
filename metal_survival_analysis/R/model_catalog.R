# Model formulas and analysis scenarios.
source("metal_survival_analysis/R/settings.R")

METALS <- c(Pb = "LBXBPB", Hg = "LBXTHG", Cd = "URXUCD", As = "VNSUMARS")
PRIMARY <- c("primary_Pb", "primary_Hg", "egfr60_Cd", "egfr60_As")
PANELS <- list(primary = names(METALS), no_bmi = names(METALS),
               age55plus = names(METALS), egfr60 = c("Cd", "As"),
               portnoy = names(METALS), joint = c("Pb", "Cd"),
               epa_dha = "Hg", rcs = c("Pb", "Cd"))
SCENARIOS <- bind_rows(lapply(names(PANELS), function(panel) {
  data.frame(scenario = paste(panel, PANELS[[panel]], sep = "_"),
             panel = panel, label = PANELS[[panel]])
})) %>%
  mutate(metal = unname(METALS[label]),
         design_metal = if_else(panel == "joint", "URXUCD", metal),
         first_cycle = if_else(label %in% c("Hg", "As"), 3L, 1L),
         period = if_else(first_cycle == 3L, "2003-2014", "1999-2014"),
         method = if_else(panel == "portnoy", "Portnoy", "PengHuang"),
         primary = scenario %in% PRIMARY)

model_covariates <- function(metal, panel) {
  variables <- c("RIDAGEYR", "RIAGENDR_f", "RIDRETH1_f", "DMDEDUC2_f",
                 "marital_3cat", "INDFMPIR", "BMXBMI", "pack_years",
                 "alq_5cat", "pa_q4_cyc", "cycle_f", "LBXTC")
  if (metal %in% c("URXUCD", "VNSUMARS")) variables <- c(variables, "egfr_ckdepi_2021")
  if (metal == "LBXTHG") {
    variables <- c(variables, if (panel == "epa_dha") "epa_dha_mg" else "fish_servings_wk")
  }
  if (panel == "no_bmi") variables <- setdiff(variables, "BMXBMI")
  variables
}

# Keep the historical CQR spline basis for comparison.
# Cox shape checks use a separate penalized spline.
restricted_cubic_basis <- function(x, knots) {
  stopifnot(length(knots) == 4L, all(diff(knots) > 0))
  positive_cube <- function(z) pmax(z, 0)^3
  nonlinear <- sapply(1:2, function(j) {
    (positive_cube(x - knots[j]) -
       positive_cube(x - knots[3]) * (knots[4] - knots[j]) / (knots[4] - knots[3]) +
       positive_cube(x - knots[4]) * (knots[3] - knots[j]) / (knots[4] - knots[3])) /
      (knots[4] - knots[1])^2
  })
  basis <- cbind(x, matrix(nonlinear, nrow = length(x), ncol = 2))
  colnames(basis) <- c("rcs1", "rcs2", "rcs3")
  basis
}

analysis_formula <- function(exposures, covariates) {
  as.formula(paste("Surv(followup_yrs, event) ~",
                  paste(c(exposures, "RIDAGEYR", "I(RIDAGEYR^2)",
                          setdiff(covariates, "RIDAGEYR")), collapse = " + ")))
}
