# Match primary intensity thresholds; retain truthful survey-specific low-intake labels.
iii_alcohol<-function(d) {
 clean<-function(x,lo,hi)ifelse(is.finite(x)&x>=lo&x<=hi,x,NA_real_)
 days<-clean(d$MAPE3S,1,365);quantity<-clean(d$MAPE4,1,90)
 daily<-days*quantity/365
 sex<-as.character(d$sex);female<-sex=='Female';valid_sex<-sex%in%c('Male','Female')
 category<-rep(NA_character_,nrow(d));level<-c('Under12Lifetime','Under12PastYear','Light','Moderate','Heavy')
 category[d$MAPE1%in%2]<-level[1]
 category[d$MAPE1%in%1 & d$MAPE2%in%2]<-level[2]
 current<-d$MAPE1%in%1 & d$MAPE2%in%1 & is.finite(daily)&valid_sex
 low<-ifelse(female,1,2);high<-ifelse(female,2,3)
 category[current]<-ifelse(daily[current]<low[current],'Light',ifelse(daily[current]<high[current],'Moderate','Heavy'))
 factor(category,levels=level)
}
