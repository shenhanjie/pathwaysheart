# Create IPW Datasets
# load R packages
require(LtAtStructuR)
require(data.table)
require(lubridate)
require(future) # optional (for parallel processing)
plan(multiprocess) # optional (for parallel processing)


create_dat_tab = function(d, x, eoftype, covar, trt){
  dat1 = d[,covar]
  
  for (i in 1:dim(dat1)[2]){
    dat1[,i] = as.character(dat1[,i])
  }
  
  dat1$ID = d$cvd_studyid
  dat1$IndexDate = d$index_date
  
  dat1$EOFDate = d[,"eofdate"]
  dat1$EOFtype = d[,eoftype]
  dat1$EOFtype = as.character(dat1$EOFtype)
  dat1$fakeDx = 0
  dat1 = as.data.table(dat1)
  
  
  dat2 = data.frame(ID = d[,"cvd_studyid"])
  trt_strat = paste0(trt, "_start_dt")
  trt_end = paste0(trt, "_end_dt")
  dat2$startA = d[,trt_strat]
  dat2$endA = d[,trt_end]
  dat2 = dat2[!is.na(dat2$startA),]
  dat2 = as.data.table(dat2)
  
  id = d[,"cvd_studyid"]
  id_long = rep(id,3)
  fake_dt = c(d[,"index_date"] + 10, d[,"index_date"] + 20, d[,"index_date"] + 30)
  fakeDx = c(rep(1,length(id)), rep(2,length(id)), rep(3,length(id)))
  dat3 = data.frame(ID = id_long, fakeDate = fake_dt, fakeDx = fakeDx)
  dat3 = dat3[order(dat3$ID),]
  dat3 = as.data.table(dat3)
  
  covar2 = c(covar, "fakeDx")
  list_cohort0 = list("categorical"=TRUE,
                      "impute"=NA,
                      "impute_default_level"=NA)
  list_cohort = rep(list(list_cohort0), length(covar))
  names(list_cohort) = covar
  
  cohort = setCohort(dat1, "ID", "IndexDate", "EOFDate", "EOFtype",
                     x, covar2, list_cohort) 
  
  exposure = setExposure(dat2, "ID", "startA", "endA")
  covariate = setCovariate(dat3, "sporadic", "ID", "fakeDate", "fakeDx",
                           categorical = FALSE)
  
  LtAt.specification = cohort + exposure + covariate 
  LtAt.data = construct(LtAt.specification, time_unit = 365, first_exp_rule = 1,
                        exp_threshold = 0.75)
  
  
  return(list(dat1 = dat1, dat2 = dat2, dat3 = dat3, LtAt.data = LtAt.data))
}


# Individual Chemotherapy
covar_chemo = c("ajcc_stage", "diab_bl", "htn_bl", "dyslipid_bl", "bmicat1", "smok", "charlson", "edu_cat", "income_cat",
                "agegrp", "race", "horm_yn","rad_tx_yn", "menop")

# Taxane
phs.aim2.data$eofdate = phs.aim2.data$death_date
phs.aim2.data$eofdate[is.na(phs.aim2.data$eofdate)] = phs.aim2.data$end_of_study

phs.aim2.data$eoftype_ihd = "None"
phs.aim2.data$eoftype_ihd[!is.na(phs.aim2.data$death_date)] = "Death"
phs.aim2.data$eoftype_ihd[phs.aim2.data$ischemic_heart_disease_grp_inc == 1] = "Ischemic Heart Disease"

# Ischemic Heart Disease
dat_ihd_taxane = create_dat_tab(phs.aim2.data, "Ischemic Heart Disease", "eoftype_ihd", covar_chemo, "taxane")




