# Read public 2019 mortality records using the CDC field positions.
source("metal_survival_analysis/R/settings.R")
source(file.path(RUN_DIR, "R/data_tools.R"))

files <- file.path(raw_directory, paste0("NHANES_", seq(1999, 2013, 2), "_",
                                       seq(2000, 2014, 2), "_MORT_2019_PUBLIC.dat"))
mortality <- lapply(seq_along(files), function(cycle) {
  lines <- readLines(files[cycle], warn=FALSE)
  stopifnot(length(lines) > 1000, all(nchar(lines) >= 46))
  field <- function(first, last) {
    x <- trimws(substr(lines, first, last))
    x[x %in% c("", ".")] <- NA_character_
    as.integer(x)
  }
  data.frame(SEQN=as.numeric(field(1, 6)), SDDSRVYR=cycle, ELIGSTAT=field(15, 15),
             MORTSTAT=field(16, 16), UCOD_LEADING=field(17, 19),
             DIABETES=field(20, 20), HYPERTEN=field(21, 21),
             PERMTH_INT=field(43, 45), PERMTH_EXM=field(46, 48))
}) %>% dplyr::bind_rows() %>%
  dplyr::mutate(event=MORTSTAT, followup_yrs=PERMTH_EXM/12)
stopifnot(!anyDuplicated(mortality$SEQN), all(mortality$ELIGSTAT %in% 1:3))
saveRDS(mortality, data_file("mortality.rds"))
message("Mortality reconstruction complete: ", nrow(mortality), " public records.")
