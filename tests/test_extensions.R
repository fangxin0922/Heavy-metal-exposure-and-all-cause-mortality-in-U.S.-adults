# Boundary checks for the scientific distinctions in the added pipelines.
source('cause_specific/R/endpoint_data.R')
x <- data.frame(SEQN=1:7,followup_yrs=c(1,2,3,4,5,6,7),
  event=c(1,1,1,1,0,0,1),UCOD_LEADING=c(1,2,5,10,NA,2,3))
m <- endpoint_map(x)
stopifnot(identical(m$cancer,c(0L,1L,0L,0L,0L,0L,0L)),
  identical(m$cardiovascular,c(1L,0L,1L,0L,0L,0L,0L)),
  identical(m$followup_yrs,x$followup_yrs),nrow(m)==nrow(x),
  sum(m$event_allcause-m$cancer)==4L)
x$UCOD_LEADING[2] <- NA
stopifnot(inherits(try(endpoint_map(x),silent=TRUE),'try-error'))
source('nhanes3_analysis/R/alcohol.R')
z <- data.frame(MAPE1=rep(1,8),MAPE2=rep(1,8),MAPE3S=rep(c(364,365,364,365),2),
 MAPE4=c(2,2,3,3,1,1,2,2),sex=rep(c('Male','Female'),each=4))
stopifnot(identical(as.character(iii_alcohol(z)),rep(c('Light','Moderate','Moderate','Heavy'),2)))
z <- data.frame(MAPE1=c(2,1,NA,1),MAPE2=c(NA,2,1,1),MAPE3S=c(NA,NA,365,888),
 MAPE4=c(NA,NA,3,2),sex='Male')
stopifnot(identical(as.character(iii_alcohol(z)),c('Under12Lifetime','Under12PastYear',NA_character_,NA_character_)))
cat('Endpoint coding, competing-death retention, invalid-cause rejection, and alcohol boundaries passed.\n')
