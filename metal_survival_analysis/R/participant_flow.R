# Reconstruct Figure 1 counts from native releases and the prepared model cohorts.
source('metal_survival_analysis/R/model_catalog.R')
source('metal_survival_analysis/R/data_tools.R')
out <- file.path(RUN_DIR,'results/participant_flow')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
demo_paths <- file.path(raw_directory,paste0(c('DEMO',paste0('DEMO_',LETTERS[2:8])),'.xpt'))
demo <- bind_rows(lapply(seq_along(demo_paths),function(cycle) {
  d <- read_native(demo_paths[cycle])
  data.frame(SEQN=d$SEQN,SDDSRVYR=cycle,RIDAGEYR=d$RIDAGEYR,RIDEXPRG=d$RIDEXPRG)
}))
mortality <- readRDS(data_file('mortality.rds'))
stopifnot(!anyDuplicated(demo$SEQN),!anyDuplicated(mortality$SEQN))
d <- left_join(demo,mortality,by=c('SEQN','SDDSRVYR'))
stopifnot(!anyNA(d$ELIGSTAT))
specification <- read.csv(file.path(RUN_DIR,'inputs/weight_specification.csv'))
specification <- specification[specification$SDDSRVYR<=8,]
codes <- c(blood_lead='LBXBPB',blood_mercury='LBXTHG',urinary_cadmium='URXUCD',urinary_arsenic_sum='VNSUMARS')
for(v in unname(codes))d[[v]] <- NA_real_
lab_paths <- character()
for(i in seq_len(nrow(specification))) {
  s <- specification[i,];metal <- unname(codes[s$metal])
  path <- file.path(raw_directory,paste0(s$lab_file,'.xpt'))
  lab <- read_native(path);lab_paths <- c(lab_paths,path)
  ix <- which(d$SDDSRVYR==s$SDDSRVYR);j <- match(d$SEQN[ix],lab$SEQN)
  if(metal=='VNSUMARS') {
    species <- as.matrix(lab[j,c('URXUAS3','URXUAS5','URXUMMA','URXUDMA')])
    values <- rowSums(species);values[rowSums(is.finite(species))!=4L] <- NA_real_
  } else {
    field <- if(metal=='URXUCD' && s$SDDSRVYR<=2)'URDUCD'else metal
    values <- lab[[field]][j]
  }
  d[[metal]][ix] <- values
}
# Early mercury measurements count toward having any measured study metal,
# although the comparable primary mercury analysis starts in 2003.
for(cycle in 1:2) {
  path <- file.path(raw_directory,paste0(c('LAB06','L06_B')[cycle],'.xpt'))
  lab <- read_native(path);lab_paths <- c(lab_paths,path)
  ix <- which(d$SDDSRVYR==cycle)
  d$LBXTHG[ix] <- lab$LBXTHG[match(d$SEQN[ix],lab$SEQN)]
}
adult <- is.finite(d$RIDAGEYR)&d$RIDAGEYR>=40
eligible <- adult & d$ELIGSTAT%in%1
followed <- eligible & is.finite(d$followup_yrs)&d$followup_yrs>0
nonpregnant <- followed & !(d$RIDEXPRG%in%1)
measured <- rowSums(is.finite(as.matrix(d[unname(codes)])))>0
included <- nonpregnant & measured
stopifnot(all(d$event[included]%in%0:1))
counts <- c(all_participants=nrow(d),age_below40=sum(!adult),age40plus=sum(adult),
  linkage_ineligible=sum(adult&!eligible),followup_missing_or_nonpositive=sum(eligible&!followed),
  linked_positive_followup=sum(followed),known_pregnancy=sum(followed&!nonpregnant),
  after_pregnancy_exclusion=sum(nonpregnant),no_measured_study_metal=sum(nonpregnant&!measured),
  overall_n=sum(included),overall_deaths=sum(d$event[included]),overall_censored=sum(d$event[included]==0))
expected <- c(82091,53487,28604,43,1489,27072,23,27049,2711,24338,6992,17346)
stopifnot(identical(as.numeric(counts),expected))
analytic <- readRDS(data_file('analytic_four_metals.rds'))
stopifnot(setequal(d$SEQN[included],analytic$SEQN),sum(analytic$event)==counts['overall_deaths'])
write_table(data.frame(stage=names(counts),n=unname(counts)),file.path(out,'selection_counts.csv'))
model_counts <- bind_rows(lapply(PRIMARY,function(sc) {
  spec <- readRDS(file.path(RUN_DIR,'models',sc,'specification.rds'))
  scenario <- SCENARIOS[SCENARIOS$scenario==sc,];scenario$require_egfr60 <- FALSE
  before <- select_analysis_cohort(analytic,scenario)
  renal <- sc%in%c('egfr60_Cd','egfr60_As')
  data.frame(scenario=sc,before_renal_n=nrow(before),not_in_before_renal_cohort=nrow(analytic)-nrow(before),
    excluded_low_egfr=if(renal)sum(is.finite(before$egfr_ckdepi_2021)&before$egfr_ckdepi_2021<60)else 0L,
    excluded_missing_egfr=if(renal)sum(!is.finite(before$egfr_ckdepi_2021))else 0L,
    n=spec$n,deaths=spec$events,censored=spec$n-spec$events)
}))
stopifnot(identical(as.numeric(model_counts$before_renal_n),c(24021,18258,8466,6669)),
  identical(as.numeric(model_counts$excluded_low_egfr),c(0,0,973,779)),
  identical(as.numeric(model_counts$excluded_missing_egfr),c(0,0,432,338)),
  all(model_counts$before_renal_n-model_counts$excluded_low_egfr-model_counts$excluded_missing_egfr==model_counts$n))
write_table(model_counts,file.path(out,'model_counts.csv'))
paths <- unique(c(demo_paths,lab_paths,data_file('mortality.rds'),data_file('analytic_four_metals.rds'),
  file.path(RUN_DIR,'models',PRIMARY,'specification.rds')))
write_table(data.frame(file=paths,md5=unname(tools::md5sum(paths))),file.path(out,'source_hashes.csv'))
cat('All native participant-selection counts and final cohort membership verified.\n')

# All graphics are generated in R; exclusions are unweighted sequential counts.
draw_flow <- function() {
  grid::grid.newpage()
  box <- function(x,y,w,h,label,fill='#F5F8FA',size=11) {
    grid::grid.roundrect(x,y,w,h,r=grid::unit(.015,'npc'),gp=grid::gpar(fill=fill,col='#8296A5'))
    grid::grid.text(label,x,y,gp=grid::gpar(fontsize=size,col='#243746'))
  }
  arrow <- function(x,y,xend,yend)grid::grid.lines(c(x,xend),c(y,yend),
    arrow=grid::arrow(length=grid::unit(.1,'inches'),type='closed'),gp=grid::gpar(col='#738793'))
  fmt <- function(x)format(x,big.mark=',',scientific=FALSE,trim=TRUE)
  ys <- c(.94,.79,.64,.49,.34)
  labels <- c(sprintf('NHANES 1999–2014 participants\nn = %s',fmt(counts['all_participants'])),
    sprintf('Adults aged ≥40 years\nn = %s',fmt(counts['age40plus'])),
    sprintf('Eligible linkage and positive examination follow-up\nn = %s',fmt(counts['linked_positive_followup'])),
    sprintf('After excluding known pregnancies\nn = %s',fmt(counts['after_pregnancy_exclusion'])),
    sprintf('At least one study metal measured\nn = %s\n%s deaths; %s censored',fmt(counts['overall_n']),fmt(counts['overall_deaths']),fmt(counts['overall_censored'])))
  for(i in seq_along(ys)) {
    box(.29,ys[i],.54,.105,labels[i],if(i==5)'#E5F0F4'else'#F5F8FA')
    if(i<5)arrow(.29,ys[i]-.0525,.29,ys[i+1]+.0525)
  }
  excluded <- c(sprintf('Excluded: %s\nAge <40 years',fmt(counts['age_below40'])),
    sprintf('Excluded: %s\nLinkage ineligible: %s\nMissing/nonpositive follow-up: %s',fmt(counts['linkage_ineligible']+counts['followup_missing_or_nonpositive']),fmt(counts['linkage_ineligible']),fmt(counts['followup_missing_or_nonpositive'])),
    sprintf('Excluded: %s\nKnown pregnancy',fmt(counts['known_pregnancy'])),
    sprintf('Excluded: %s\nNo measured study metal',fmt(counts['no_measured_study_metal'])))
  for(i in 1:4) {
    y <- mean(ys[i:(i+1)])
    arrow(.29,y,.62,y);box(.8,y,.36,.095,excluded[i],size=10)
  }
  xs <- c(.13,.38,.63,.88)
  titles <- c('Blood lead','Blood mercury','Urinary cadmium','Urinary arsenic sum')
  arrow(.29,.2875,.29,.255)
  grid::grid.lines(c(xs[1],xs[4]),c(.255,.255),gp=grid::gpar(col='#738793'))
  for(i in 1:4) {
    z <- model_counts[i,]
    arrow(xs[i],.255,xs[i],.225)
    extra <- if(i>=3)sprintf('Before eGFR screen: %s\nNot included: %s\neGFR <60: %s; missing: %s\nPrimary analysis: eGFR ≥60',
      fmt(z$before_renal_n),fmt(z$not_in_before_renal_cohort),fmt(z$excluded_low_egfr),fmt(z$excluded_missing_egfr))else
      sprintf('Not included: %s\nPrimary analysis',fmt(z$not_in_before_renal_cohort))
    label <- sprintf('%s\n%s\n%s\nn = %s\nDeaths: %s\nCensored: %s',titles[i],
      if(i%in%c(1,3))'1999–2014'else'2003–2014',extra,fmt(z$n),fmt(z$deaths),fmt(z$censored))
    box(xs[i],.135,.23,.18,label,'#E5F0F4',size=10)
  }
  grid::grid.text('Cohorts overlap. Counts are unweighted. Missing covariates are imputed after cohort selection.',
    .5,.025,gp=grid::gpar(fontsize=9,col='#526878'))
}
grDevices::cairo_pdf(file.path(RUN_DIR,'figures/participant_flow.pdf'),width=11,height=14)
draw_flow();dev.off()
png(file.path(RUN_DIR,'figures/participant_flow.png'),width=1650,height=2100,res=150,type='cairo')
draw_flow();dev.off()
