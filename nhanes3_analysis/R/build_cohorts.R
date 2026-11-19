source('nhanes3_analysis/R/settings.R')
source('nhanes3_analysis/R/alcohol.R')
stopifnot(!file.exists(file.path(base,'data/full_mec.rds')))
read <- function(n) read.csv(file.path(base,'data',paste0(n,'_selected.csv')),na.strings=c('','.','NA'))
x <- read('exam')
for(n in c('lab','adult','mortality')) {
 y<-read(n);stopifnot(!anyDuplicated(y$SEQN),!anyDuplicated(x$SEQN))
 z<-match(x$SEQN,y$SEQN);x<-cbind(x,y[z,setdiff(names(y),'SEQN'),drop=FALSE])
}
clean <- function(z,lo,hi) ifelse(is.finite(z)&z>=lo&z<=hi,z,NA_real_)
# Codebook ranges remove special codes without inventing broad sentinel rules.
x$age<-clean(x$HSAGEIR,0,90);x$age2<-x$age^2
x$sex<-factor(x$HSSEX,levels=c(1,2),labels=c('Male','Female'))
r<-ifelse(x$DMAETHNR==1,'MexicanAmerican',ifelse(x$DMAETHNR==2,'OtherHispanic',
 ifelse(x$DMAETHNR==3 & x$DMARACER==1,'NHWhite',ifelse(x$DMAETHNR==3 & x$DMARACER==2,'NHBlack',
 ifelse(x$DMAETHNR==3 & x$DMARACER==3,'Other',NA)))))
x$race<-factor(r,levels=c('NHWhite','NHBlack','MexicanAmerican','OtherHispanic','Other'))
e<-clean(x$HFA8R,0,17)
x$education<-factor(ifelse(e<12,'LessThan12',ifelse(e==12,'12Years','MoreThan12')),levels=c('MoreThan12','12Years','LessThan12'))
x$marital<-factor(ifelse(x$HFA12%in%1:3,'MarriedCohabiting',ifelse(x$HFA12%in%4:6,'PreviouslyMarried',ifelse(x$HFA12==7,'NeverMarried',NA))))
x$pir<-pmin(clean(x$DMPPIR,0,11.889),5)
x$bmi<-clean(x$BMPBMI,7.3,79.6)
x$cholesterol<-clean(x$TCP,59,702)
x$phase<-factor(x$SDPPHASE,levels=1:2,labels=c('1988_1991','1991_1994'))
start<-clean(x$HAR2,5,78);quit<-clean(x$HAR11R,8,90)
intensity<-ifelse(x$HAR3==1,x$HAR4S,ifelse(x$HAR3==2,x$HAR12S,NA_real_))
intensity<-ifelse(intensity==777,.5,clean(intensity,1,240))
duration<-ifelse(x$HAR3==1,x$age-start,ifelse(x$HAR3==2,quit-start,NA_real_))
duration[!is.na(quit)&x$HAR3%in%2&quit>x$age]<-NA_real_
x$packyears<-ifelse(x$HAR1==2,0,ifelse(x$HAR1==1,intensity/20*clean(duration,0,90),NA_real_))
weekly<-clean(x$MAPE3S,1,365)*clean(x$MAPE4,1,90)/365.25*7
limit<-ifelse(x$sex=='Female',7,14)
a<-ifelse(x$MAPE1==2,'Under12Lifetime',ifelse(x$MAPE1==1 & x$MAPE2==2,'Under12PastYear',
 ifelse(x$MAPE1==1 & x$MAPE2==1,ifelse(weekly<=3,'Light',ifelse(weekly<=limit,'Moderate','Heavy')),NA)))
x$alcohol<-factor(a,levels=c('Under12Lifetime','Under12PastYear','Light','Moderate','Heavy'))
walk<-clean(x$HAT1S,0,100);acts<-as.matrix(x[paste0('HAT',seq(2,18,2))]);acts[!acts%in%c(1,2)]<-NA_real_
active<-(is.finite(walk)&walk>0)|rowSums(acts==1,na.rm=TRUE)>0
inactive<-is.finite(walk)&walk==0&rowSums(acts==2,na.rm=TRUE)==ncol(acts)
x$activity<-factor(ifelse(active,'Active',ifelse(inactive,'Inactive',NA)),levels=c('Inactive','Active'))
x$pb<-clean(x$PBP,.7,71.8);x$cd<-clean(x$UDP,.01,16.65);x$ucr<-clean(x$URP,7.9,682.1)
scr<-.960*clean(x$CEP,.3,13.9)-.184
k<-ifelse(x$sex=='Female',.7,.9);alpha<-ifelse(x$sex=='Female',-.241,-.302)
x$egfr<-142*pmin(scr/k,1)^alpha*pmax(scr/k,1)^(-1.2)*.9938^x$age*ifelse(x$sex=='Female',1.012,1)
x$event<-x$MORTSTAT;x$time<-x$PERMTH_EXM/12
stopifnot(nrow(x)==31311L,!anyDuplicated(x$SEQN))
# Full MEC design, including all ages and those without a laboratory measurement.
full<-x[x$DMPSTAT%in%2 & is.finite(x$WTPFEX6)&x$WTPFEX6>0,,drop=FALSE]
stopifnot(!anyNA(full[c('SEQN','WTPFEX6','SDPSTRA6','SDPPSU6')]),
 length(unique(full$SDPSTRA6))==49L,length(unique(paste(full$SDPSTRA6,full$SDPPSU6)))==98L)
eligible<-is.finite(full$age)&full$age>=40&full$ELIGSTAT%in%1&full$event%in%c(0,1)&
 is.finite(full$time)&full$time>0&!(full$MAPF12R%in%1)
crvars<-c('age','age2','sex','race','bmi','phase','ucr')
reference<-eligible&complete.cases(full[crvars])&full$ucr>0
stopifnot(!anyNA(reference),sum(reference)>1000)
crfit<-lm(log(ucr)~age+age2+sex+race+bmi+phase,data=full[reference,])
full$predicted_ucr<-exp(predict(crfit,newdata=full,na.action=na.pass))
full$cd_adjusted<-full$cd*full$predicted_ucr/full$ucr
saveRDS(crfit,file.path(base,'audit','urine_creatinine_reference.rds'))
write.csv(data.frame(term=names(coef(crfit)),coefficient=coef(crfit)),file.path(base,'audit','urine_creatinine_reference.csv'),row.names=FALSE)
write.csv(data.frame(reference_n=sum(reference),R_squared=summary(crfit)$r.squared),file.path(base,'audit','urine_creatinine_reference_summary.csv'),row.names=FALSE)
flow<-data.frame(stage=c('Examination release','Full MEC positive weight','Age40plus','Mortality eligible','Valid positive exam followup','Not known pregnant'),
 n=c(nrow(x),nrow(full),sum(full$age>=40,na.rm=TRUE),sum(full$age>=40&full$ELIGSTAT%in%1,na.rm=TRUE),
 sum(full$age>=40&full$ELIGSTAT%in%1&full$event%in%c(0,1)&is.finite(full$time)&full$time>0,na.rm=TRUE),sum(eligible)))
write.csv(flow,file.path(base,'audit','cohort_flow.csv'),row.names=FALSE)
# Full raw/range audit to expose import/scaling/sentinel mistakes.
rawvars<-c('HSAGEIR','DMPPIR','WTPFEX6','BMPBMI','PBP','UDP','URP','CEP','TCP','MAPE3S','MAPE4')
write.csv(do.call(rbind,lapply(rawvars,function(v){z<-x[[v]];data.frame(variable=v,observed=sum(!is.na(z)),minimum=min(z,na.rm=TRUE),median=median(z,na.rm=TRUE),maximum=max(z,na.rm=TRUE))})),file.path(base,'audit','raw_range_audit.csv'),row.names=FALSE)
saveRDS(full,file.path(base,'data','full_mec.rds'))
covars<-c('age','age2','sex','race','education','marital','pir','bmi','packyears','alcohol','activity','phase','cholesterol')

# Use the final harmonised alcohol categories before model preparation.
full$alcohol <- iii_alcohol(full)
saveRDS(full,file.path(base,'data/full_mec.rds'))
for (metal in c('lead','cadmium')) {
  exposure <- if(metal=='lead')full$pb else full$cd_adjusted
  keep <- eligible & is.finite(exposure) & exposure>0
  if(metal=='cadmium')keep <- keep & is.finite(full$egfr) & full$egfr>=60
  keep[is.na(keep)] <- FALSE
  d <- full[keep,,drop=FALSE]
  d$exposure <- log2(if(metal=='lead')d$pb else d$cd_adjusted)
  stopifnot(nrow(d)==if(metal=='lead')9388L else 7924L,
            sum(d$event)==if(metal=='lead')6146L else 4823L)
  dest <- file.path(base,'models',metal)
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  saveRDS(d,file.path(dest,'cohort.rds'))
}
writeLines(capture.output(sessionInfo()),file.path(base,'audit/build_session.txt'))
