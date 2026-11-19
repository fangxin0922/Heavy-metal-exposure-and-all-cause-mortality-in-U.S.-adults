# Native public files consumed by the existing analysis scripts.
public_input_manifest <- function(include_biomarkers = FALSE) {
  hm <- Sys.getenv("HM_DATA_DIR", "data/metal_survival")
  emm <- Sys.getenv("EMM_DATA_DIR", "data/biomarker_emm")
  xpt <- function(root, folder, files, cycles) {
    data.frame(target = file.path(root, folder, paste0(files, ".xpt")),
      url = sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/%s.xpt",
                    1997L + 2L * cycles, files), stringsAsFactors = FALSE)
  }
  weights <- read.csv("metal_survival_analysis/inputs/weight_specification.csv")
  weights <- unique(weights[weights$SDDSRVYR <= 8, c("lab_file", "SDDSRVYR")])
  result <- rbind(
    xpt(hm, "raw", c("DEMO", paste0("DEMO_", LETTERS[2:8])), 1:8),
    xpt(hm, "raw", weights$lab_file, weights$SDDSRVYR),
    xpt(hm, "creatinine", c("L16_B", "ALB_CR_I", "ALB_CR_J"), c(2, 9, 10)),
    xpt(hm, "baseline", c("L40_B", "BMX_B", "L13_B"), rep(2, 3)),
    xpt(hm, "diet", c("DRXTOT_B", paste0("DR1TOT_", LETTERS[3:8])), 2:8))
  mortality <- sprintf("NHANES_%d_%d_MORT_2019_PUBLIC.dat", seq(1999, 2013, 2), seq(2000, 2014, 2))
  result <- rbind(result, data.frame(target = file.path(hm, "raw", mortality),
    url = paste0("https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/", mortality)))
  if (include_biomarkers) {
    result <- rbind(result,
      xpt(emm, "raw/labs", paste0("DEMO_", LETTERS[2:4]), 2:4),
      xpt(emm, "raw/labs", c("L06HM_B", "L06HM_C", "UHM_D", "L06UAS_C", "UAS_D"), c(2:4, 3:4)),
      xpt(emm, "raw/biomarkers", c("L06VIT_C", "VIC_D", "L06VIT_B", "L45VIT_C", "VITAEC_D",
        "LAB06", "L06_B", "L06NB_C", paste0("FOLATE_", LETTERS[4:8])), c(3:4, 2:4, 1:8)),
      xpt(emm, "raw/cbc", c("LAB25", "L25_B", "L25_C", "CBC_D"), 1:4))
  }
  result <- unique(result)
  stopifnot(!anyDuplicated(result$target))
  result
}

validate_reference_input <- function(x) {
  required <- c("SEQN", "SDDSRVYR", "RIDAGEYR", "RIAGENDR", "RIDRETH1", "DMDEDUC2",
    "DMDMARTL", "INDFMPIR", "SDMVSTRA", "SDMVPSU", "URXUCR", "URXUTL", "BMXBMI",
    "LBXSCR", "LBXSCR_cal", "LBXTC", "DRXTP205", "DRXTP226", "pack_years",
    "alq_5cat", "pa_q4_cyc", "pa_mets_wkly_all", "fish_servings_wk", "fish_servings_wk_t1",
    "fish_servings_wk_t2", "dark_fish_30d_anchored", "fish_source",
    "log_primary_LBXBPB", "log_primary_LBXTHG", "RIAGENDR_f", "RIDRETH1_f",
    "DMDEDUC2_f", "marital_3cat")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Reference archive is missing: ", paste(missing, collapse = ", "))
  stopifnot(is.data.frame(x), nrow(x) > 0,
    is.numeric(x$SEQN), is.numeric(x$SDDSRVYR), is.numeric(x$RIDAGEYR),
    is.numeric(x$pa_mets_wkly_all),
    !anyNA(x$SEQN), !anyDuplicated(x$SEQN),
    !anyNA(x$SDDSRVYR), identical(sort(unique(as.integer(x$SDDSRVYR))), 1:10),
    all(is.finite(x$RIDAGEYR) & x$RIDAGEYR >= 40),
    all(vapply(x[c("RIAGENDR_f", "RIDRETH1_f", "DMDEDUC2_f", "marital_3cat",
                   "alq_5cat", "pa_q4_cyc")], is.factor, logical(1))))
  invisible(TRUE)
}
