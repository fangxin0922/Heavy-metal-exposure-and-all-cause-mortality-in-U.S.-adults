# Keep competing deaths in the cohort, censoring them at their observed death time.
endpoint_map <- function(data) {
  cause <- suppressWarnings(as.integer(as.character(data$UCOD_LEADING)))
  event <- as.integer(data$event)
  stopifnot(!anyDuplicated(data$SEQN),all(event %in% 0:1),
    all(cause[event==1L] %in% 1:10),all(is.finite(data$followup_yrs)),all(data$followup_yrs>0))
  data.frame(SEQN=as.numeric(data$SEQN),followup_yrs=data$followup_yrs,event_allcause=event,
    UCOD_LEADING=cause,cancer=as.integer(event==1L & cause %in% 2L),
    cardiovascular=as.integer(event==1L & cause %in% c(1L,5L)))
}
