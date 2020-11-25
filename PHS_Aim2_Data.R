##### Pathways Heart Study Aim 2 Data Generation ######
### Hanjie Shen, updated on 11/25/2020

require(tidyverse)
require(haven)
require(ggplot2)
require(dplyr)
require(plyr)
require(naniar)
require(reshape2)
require(doParallel)
require(lsr)
require(LtAtStructuR)
require(data.table)
require(lubridate)
require(future) # optional (for parallel processing)
plan(multiprocess) # optional (for parallel processing)


# setup parallel computing
nodes = detectCores()
cl = makeCluster(nodes)
registerDoParallel(cl)

################################################################################
# Load data
setwd("~/Desktop/Pathways Heart Study/Data/full data")
# setwd("C:/Users/shenh/Desktop/Pathways Heart Study/Data/full data")

# cases data
cases = as.data.frame(read_sas('cases.sas7bdat'))
cases_rad = as.data.frame(read_sas('aim1_cases_rad_tx_02mar2020.sas7bdat'))
cases$rad_tx_yn = cases_rad$rad_tx_yn[match(cases_rad$CVD_studyid,cases$CVD_studyid)]
cases$laterality_tx = cases_rad$LATERALITY[match(cases_rad$CVD_studyid,cases$CVD_studyid)]
cases$rad_tx_start_dt = cases_rad$rad_tx_start_dt[match(cases_rad$CVD_studyid,cases$CVD_studyid)]
cases = cases[,!names(cases) %in% c("rad_yn","rad_date","laterality")]

# controls data 
controls = as.data.frame(read_sas('controls.sas7bdat'))
controls$cvd_caseid = controls$case_id

# risk factor data 
bmi = as.data.frame(read_sas('baseline_bmi.sas7bdat'))
bp = as.data.frame(read_sas('baseline_bp.sas7bdat'))
lipid = as.data.frame(read_sas('dyslipidemia.sas7bdat'))
diab = as.data.frame(read_sas('diabetes.sas7bdat'))
smok = as.data.frame(read_sas('smoking.sas7bdat'))
hyper = as.data.frame(read_sas('phase_hypertension.sas7bdat'))

# income and education data
educ_income = as.data.frame(read_sas('pw_educ_income.sas7bdat'))

# menopause, parity, census data 
menop = as.data.frame(read_sas('menopause.sas7bdat'))
census = as.data.frame(read_sas('census.sas7bdat'))


# censoring events 
censor = as.data.frame(read_sas('censoring.sas7bdat'))
censor_phs = as.data.frame(read_sas('censoring_pathways.sas7bdat'))
censor_phs_new = as.data.frame(read_sas('cvd_censoring_19nov2020.sas7bdat'))


# lab data  
labs = as.data.frame(read_sas('baseline_labs.sas7bdat'))

# chemotherapy data
chemo = as.data.frame(read_sas('aim2_chemo_data_08apr2020.sas7bdat'))

# hormonal data
horm = as.data.frame(read_sas('aim2_horm_rx_08may2020.sas7bdat'))
names(horm) = tolower(names(horm))

# EF dataset
ef = as.data.frame(read_sas('aim2_ef.sas7bdat'))

# Charlson dataset
charlson = as.data.frame(read_sas('kwan_num0090_01jun20.sas7bdat'))

# Labs
labs$result = as.numeric(labs$RESULT_C)
names(labs) = tolower(names(labs))
labs_w = dcast(data=labs,cvd_studyid~test_type,
                value.var = 'result', mean, na.rm=T)


# CVD outcome data 
cvd = as.data.frame(read_sas('cvd_events.sas7bdat'))

# clean CVD data
names(cvd) = tolower(names(cvd))
# recode cvd group
cvd$cvd_group[cvd$cvd_condition=='Heart Failure'] = 'Heart Failure'
cvd$cvd_group[cvd$cvd_condition=='Cardiomyopathy'] = 'Cardiomyopathy'
cvd$cvd_group[cvd$cvd_condition=='Transient Ischemic Attack (TIA)'] = 'TIA'
cvd$cvd_group[cvd$cvd_condition %in% c('Acute Ischemic stroke',
                                       'Acute myocardial infarction',
                                       'Intracerebral hemorrhage',
                                       'Other cerebrovascular disease',
                                       'Retinal vascular occlusion',
                                       'Subarachnoid hemorrhage')] = 'Stroke'

#### define true incidence and prevalence
cvd2 = dcast(data=cvd, cvd_studyid ~ cvd_group,
             value.var = 'daysto_cvd_event', min,na.rm=T)

cvd2[,-1] = lapply(cvd2[,-1], function(x) {
  x[x==Inf] = NA
  x
})

names(cvd2) = gsub('/| ','_',names(cvd2))
cvd_grp = paste0(names(cvd2)[-1],'_grp')
names(cvd2)[-1] = paste0(cvd_grp,'_dt')

cvd2[,cvd_grp] = lapply(cvd2[,-1], function(x) {
  y = ifelse(!is.na(x),1,0)
  y[y==1 & x>0] = 2
  y = factor(y, levels=0:2, labels = c('No','Prevalent','Incident'))
  y
})

cvd2list = lapply(cvd_grp, function(x){
  dates = eval(parse(text=paste0("dcast(data=cvd2, cvd_studyid~",x,", value.var = '",x,"_dt')")))
  prev = ifelse(cvd2[,x]=='Prevalent',1,0)
  inc = ifelse(cvd2[,x]=='Incident',1,0)
  sum = data.frame(cbind(prev,inc,dates))
  sum = sum[,c('cvd_studyid','prev','inc','Incident')]
  names(sum) = c('cvd_studyid',paste0(x,c('_prev','_inc','_incdt')))
  sum
})

#### define recurrence
cvd2r = dcast(data=cvd[which(cvd$daysto_cvd_event>=0),], 
              cvd_studyid~cvd_group, value.var = 'daysto_cvd_event', min,na.rm=T)

cvd2r[,-1] = lapply(cvd2r[,-1], function(x) {
  x[x==Inf] = NA
  x
})

names(cvd2r) = gsub('/| ','_',names(cvd2r))
cvd_grpr = paste0(names(cvd2r)[-1],'_grp_rec')
names(cvd2r)[-1] = paste0(cvd_grpr,'dt')

cvd2r[,cvd_grpr] = lapply(cvd2r[,-1], function(x) {
  y = ifelse(!is.na(x),1,0)
  y
})

#### merge all together
cvd2list[[12]] = cvd2
cvd2list[[13]] = cvd2r
cvd2_final = Reduce(function(...) merge(...,by='cvd_studyid', all.x=T), cvd2list)


#### define true incidence and prevalence
cvd3 = dcast(data=cvd, cvd_studyid~cvd_condition, 
             value.var = 'daysto_cvd_event', min,na.rm=T)


names(cvd3) = gsub('/| ','_',names(cvd3))
names(cvd3) = gsub('&_|\\(|\\)','',names(cvd3))
names(cvd3)[17] = "Percutaneous_coronary_intervention_PCI"

cvd3[,-1] = lapply(cvd3[,-1], function(x) {
  x[x==Inf] = NA
  x
})

cvd_cond = names(cvd3)[-1]
names(cvd3)[-1] = paste0(names(cvd3)[-1],'_dt')
cvd3[,cvd_cond] = lapply(cvd3[,-1],  function(x) {
  y = ifelse(!is.na(x),1,0)
  y[y==1 & x>0] = 2
  y = factor(y, levels=0:2, labels = c('No','Prevalent','Incident'))
  y
})

cvd3list = lapply(cvd_cond, function(x){
  dates = eval(parse(text=paste0("dcast(data=cvd3, cvd_studyid~",x,", value.var = '",x,"_dt')")))
  prev = ifelse(cvd3[,x]=='Prevalent',1,0)
  inc = ifelse(cvd3[,x]=='Incident',1,0)
  sum = data.frame(cbind(prev,inc,dates))
  sum = sum[,c('cvd_studyid','prev','inc','Incident')]
  names(sum) = c('cvd_studyid',paste0(x,c('_prev','_inc','_incdt')))
  sum
})

#### define recurrence
cvd3r = dcast(data=cvd[which(cvd$daysto_cvd_event>=0),], 
              cvd_studyid~cvd_condition, value.var = 'daysto_cvd_event', min,na.rm=T)

cvd3r[,-1] = lapply(cvd3r[,-1], function(x) {
  x[x==Inf] = NA
  x
})

names(cvd3r) = gsub('/| ','_',names(cvd3r))
names(cvd3r) = gsub('&_|\\(|\\)','',names(cvd3r))
names(cvd3r)[17] = "Percutaneous_coronary_intervention_PCI"
cvd_condr = paste0(names(cvd3r)[-1],'_rec')
names(cvd3r)[-1] = paste0(cvd_condr,'dt')

cvd3r[,cvd_condr] = lapply(cvd3r[,-1], function(x) {
  y = ifelse(!is.na(x),1,0)
  y
})


cvd3list[[23]] = cvd3
cvd3list[[24]] = cvd3r
cvd3_final = Reduce(function(...) merge(...,by='cvd_studyid', all.x=T), cvd3list)


################################################################################################################################

####################################################################
################## Create Aim 2 Analytic Data ######################
####################################################################

datalist = list(chemo,
                 cases[,-1],
                 bmi[,c("CVD_studyid","BMI","BMI_Category_II","daysto_bmi")],
                 bp[,c("CVD_studyid","systolic","diastolic","daysto_bp","hypertension")],
                 labs_w[,c("cvd_studyid","GLU_F","GTT75_PRE","HDL","HGBA1C","LDL_CLC_NS","TOT_CHOLES","TRIGL_NS")],
                 diab[,c(6,7,1)],lipid[,c(4,5,2)], hyper[,c(6,2,7)],
                 smok[,c(5,2)],menop[,c(1,3)],
                 census[,c(4,5,6,7,8,9,10,11,12,13,14,15,16,2)],
                 cvd2_final,
                 cvd3_final,
                 censor[,c(1,9,10,11,12,13)],
                 censor_phs)

## make all var names lower case
datalist = lapply(datalist,function(x){
  names(x) = tolower(names(x))
  x
})

## merge all data by cvd_studyid
phs.aim2.data = Reduce(function(...) merge(...,by='cvd_studyid', all.x=T), datalist)

phs.aim2.data$ischemic_heart_disease_event = cvd$ischemic_heart_disease[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]
phs.aim2.data$heart_failure_event = cvd$heart_failure[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]
phs.aim2.data$cardiomyopathy_event = cvd$cardiomyopathy[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]
phs.aim2.data$stroke_event = cvd$stroke[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]
phs.aim2.data$tia_event = cvd$tia[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]
phs.aim2.data$valvular_disease_event = cvd$valvular_disease[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]
phs.aim2.data$myocarditis_pericarditis_event = cvd$myocarditis_pericarditis[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]
phs.aim2.data$carotid_disease_event = cvd$carotid_disease[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]
phs.aim2.data$cardiac_arrest_event = cvd$cardiac_arrest[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]
phs.aim2.data$arrhythmia_event = cvd$arrhythmia[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]
phs.aim2.data$vte_event = cvd$vte[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]

phs.aim2.data$cvd_condition = cvd$cvd_condition[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]
phs.aim2.data$cvd_group = cvd$cvd_group[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]
phs.aim2.data$daysto_cvd_event = cvd$daysto_cvd_event[match(phs.aim2.data$cvd_studyid, cvd$cvd_studyid)]

phs.aim2.data$daysto_cvd_dt = phs.aim2.data$index_date + phs.aim2.data$daysto_cvd_event


# format and create new variables
# age groups
phs.aim2.data$agegrp = cut(phs.aim2.data$dxage, c(0, 40, 50, 60, 70, 101), right=F,
                           labels = c('<40 yo','40-49','50-59','60-69','70+ yo'))

# race
phs.aim2.data$race = factor(phs.aim2.data$raceethn1, 
                            labels=c("WHITE","BLACK","ASIAN","HISPANIC","PI","AI-AN"))

# postive nodes
phs.aim2.data$nodal = factor(phs.aim2.data$nodal, labels = c('Positive','Negative'))

# tumor markers
phs.aim2.data[,c('er','pr','her2')] = lapply(phs.aim2.data[,c('er','pr','her2')], function(x)
  factor(x, levels = c(0,1,2,3,8,9),
         labels = c('Not done','Positive','Negative',
                    'Borderline','Ordered, N/A',
                    'Unknown')))

# erpr status
phs.aim2.data$erpr = factor(phs.aim2.data$erpr, labels = c('ER+/PR+','ER+/PR-','ER-/PR+','ER-/PR-',
                                                           'UNKNOWN'))

# more tumor markers and treatment received
phs.aim2.data[,c('tri_neg','chemo_yn','rad_tx_yn','horm_yn')] = 
  lapply(phs.aim2.data[,c('tri_neg','chemo_yn','rad_tx_yn','horm_yn')], function(x)
    factor(x, levels = c(0,1,9),
           labels = c('No','Yes','Other')))

# AJCC stage
phs.aim2.data$ajcc_stage = factor(phs.aim2.data$ajcc_stage,
                                  labels = c('Stage I','Stage II', 'Stage III','Stage IV'))


# BMI category
phs.aim2.data$bmi1 = as.numeric(gsub('[<|>|=]','',phs.aim2.data$bmi))
phs.aim2.data$bmicat = cut(phs.aim2.data$bmi1, c(0,18.5,25,30,35,Inf),right=F,
                           labels = c('Underweight','Normal','Overweight','Obese I','Obese II+'))
phs.aim2.data$bmicat = as.character(phs.aim2.data$bmicat)
phs.aim2.data$bmicat[phs.aim2.data$raceethn1=='ASIAN' & (phs.aim2.data$bmi1 >= 23 & phs.aim2.data$bmi1 < 27.5)] = 'Overweight'
phs.aim2.data$bmicat[phs.aim2.data$raceethn1=='ASIAN' & (phs.aim2.data$bmi1 >= 27.5 & phs.aim2.data$bmi1 < 35)] = 'Obese I'
phs.aim2.data$bmicat[is.na(phs.aim2.data$bmicat)] = 'Unknown'
phs.aim2.data$bmicat = factor(phs.aim2.data$bmicat,
                              levels = c('Underweight','Normal','Overweight','Obese I','Obese II+','Unknown'))

# smoking status
phs.aim2.data$smok = phs.aim2.data$smoke_status_6m
phs.aim2.data$smok[phs.aim2.data$smok %in% c('Q')] = 'Quited'
phs.aim2.data$smok[phs.aim2.data$smok %in% c('Y')] = 'Current smoker'
phs.aim2.data$smok[phs.aim2.data$smok %in% c('N')] = 'Never smoker'
phs.aim2.data$smok[(phs.aim2.data$smok %in% c(''))|is.na(phs.aim2.data$smok)] = 'Unknown'

phs.aim2.data$smok = factor(phs.aim2.data$smok, levels=c('Never smoker','Current smoker','Quited','Unknown'))

# menopausal status
phs.aim2.data$menop = factor(phs.aim2.data$bl_meno_status)
phs.aim2.data$menop[is.na(phs.aim2.data$menop) & phs.aim2.data$age > 51] = '1'
phs.aim2.data$menop[is.na(phs.aim2.data$menop) & phs.aim2.data$age <= 51] = '0'

# baseline diabetes                                                 
phs.aim2.data$diab_bl = ifelse(phs.aim2.data$daysto_diabetes<=0, 1, 0)
phs.aim2.data$diab_bl[is.na(phs.aim2.data$diab_bl)] = 0
phs.aim2.data$diab_bl = factor(phs.aim2.data$diab_bl,levels=c(1,0),labels = c('Yes','No'))

# baseline dyslipidemia
phs.aim2.data$dyslipid_bl = ifelse(phs.aim2.data$daysto_dyslipidemia<=0, 1, 0)
phs.aim2.data$dyslipid_bl[is.na(phs.aim2.data$dyslipid_bl)] = 0
phs.aim2.data$dyslipid_bl = factor(phs.aim2.data$dyslipid_bl,levels=c(1,0),labels = c('Yes','No'))

# baseline hypertension
phs.aim2.data$htn_bl = ifelse(phs.aim2.data$daysto_phase_htn<=0, 1, 0)
phs.aim2.data$htn_bl[is.na(phs.aim2.data$htn_bl)] = 0
phs.aim2.data$htn_bl = factor(phs.aim2.data$htn_bl,levels=c(1,0),labels = c('Yes','No'))

# education
phs.aim2.data$edu_cat = educ_income$BL_EDUCLVL_5CAT[match(phs.aim2.data$cvd_studyid,educ_income$CVD_studyid)]
phs.aim2.data$edu_cat = as.factor(phs.aim2.data$edu_cat)

# income
phs.aim2.data$income_cat = educ_income$BL_INC_CAT2[match(phs.aim2.data$cvd_studyid,educ_income$CVD_studyid)]
phs.aim2.data$income_cat = as.factor(phs.aim2.data$income_cat)



### Define CVD outcomes
## set all NA for cvd outcomes as 0
phs.aim2.data[,tolower(c(cvd_grp,cvd_cond))] = lapply(phs.aim2.data[,tolower(c(cvd_grp,cvd_cond))], function(x) {
  x[is.na(x)] = 'No'
  x
})

phs.aim2.data[,grep('inc$|prev$|rec$',names(phs.aim2.data))] = 
  lapply(phs.aim2.data[,grep('inc$|prev$|rec$',names(phs.aim2.data))], function(x) {
    x[is.na(x)] = 0
    x
  })


# define true incindence
# combined event of ischemic heart disease, stroke/TIA, cardiomyopathy/heart failure
phs.aim2.data$cvdcombo_grp_prev = ifelse(phs.aim2.data$ischemic_heart_disease_grp_prev==1 | 
                                           phs.aim2.data$stroke_grp_prev==1 | 
                                           phs.aim2.data$cardiomyopathy_grp_prev==1|
                                           phs.aim2.data$heart_failure_grp_prev==1, 1, 0)


phs.aim2.data$cvdcombo_grp_inc = ifelse(phs.aim2.data$ischemic_heart_disease_grp_inc==1 | 
                                          phs.aim2.data$stroke_grp_inc==1 |  
                                          phs.aim2.data$cardiomyopathy_grp_inc==1|
                                          phs.aim2.data$heart_failure_grp_inc==1, 1, 0)


phs.aim2.data$cvdcombo_grp_incdt = apply(phs.aim2.data[,c("ischemic_heart_disease_grp_incdt",
                                                          "stroke_grp_incdt",
                                                          "cardiomyopathy_grp_incdt",
                                                          "heart_failure_grp_incdt")],1,min,na.rm=T)

phs.aim2.data$heart_failure_cardiomyopathy_grp_prev = ifelse(phs.aim2.data$cardiomyopathy_grp_prev==1|
                                           phs.aim2.data$heart_failure_grp_prev==1, 1, 0)

phs.aim2.data$heart_failure_cardiomyopathy_grp_inc =  ifelse(phs.aim2.data$cardiomyopathy_grp_inc==1|
                                                               phs.aim2.data$heart_failure_grp_inc==1, 1, 0)

phs.aim2.data$heart_failure_cardiomyopathy_grp_incdt = apply(phs.aim2.data[,c("cardiomyopathy_grp_incdt",
                                                          "heart_failure_grp_incdt")],1,min,na.rm=T)

# define censoring time
phs.aim2.data$censor_dt = apply(phs.aim2.data[,c("daysto_death","daysto_enr_end","daysto_end_study")],1,min,na.rm=T)


# calculate time of follow up
phs.aim2.data$ischemic_heart_disease_grp_inc_fu = phs.aim2.data$ischemic_heart_disease_grp_incdt
phs.aim2.data$ischemic_heart_disease_grp_inc_fu[which(phs.aim2.data$ischemic_heart_disease_grp_inc==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$ischemic_heart_disease_grp_inc==0)] 

phs.aim2.data$stroke_grp_inc_fu = phs.aim2.data$stroke_grp_incdt
phs.aim2.data$stroke_grp_inc_fu[which(phs.aim2.data$stroke_grp_inc==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$stroke_grp_inc==0)]

phs.aim2.data$cardiomyopathy_grp_inc_fu = phs.aim2.data$cardiomyopathy_grp_incdt
phs.aim2.data$cardiomyopathy_grp_inc_fu[which(phs.aim2.data$cardiomyopathy_grp_inc==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$cardiomyopathy_grp_inc==0)] 

phs.aim2.data$heart_failure_grp_inc_fu = phs.aim2.data$heart_failure_grp_incdt
phs.aim2.data$heart_failure_grp_inc_fu[which(phs.aim2.data$heart_failure_grp_inc==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$heart_failure_grp_inc==0)] 

phs.aim2.data$heart_failure_cardiomyopathy_grp_inc_fu = phs.aim2.data$heart_failure_cardiomyopathy_grp_incdt
phs.aim2.data$heart_failure_cardiomyopathy_grp_inc_fu[which(phs.aim2.data$heart_failure_cardiomyopathy_grp_inc==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$heart_failure_cardiomyopathy_grp_inc==0)] 

phs.aim2.data$cvdcombo_grp_inc_fu = phs.aim2.data$cvdcombo_grp_incdt
phs.aim2.data$cvdcombo_grp_inc_fu[which(phs.aim2.data$cvdcombo_grp_inc==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$cvdcombo_grp_inc==0)] 

phs.aim2.data$arrhythmia_grp_inc_fu = phs.aim2.data$arrhythmia_grp_incdt
phs.aim2.data$arrhythmia_grp_inc_fu[which(phs.aim2.data$arrhythmia_grp_inc==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$arrhythmia_grp_inc==0)] 

phs.aim2.data$cardiac_arrest_grp_inc_fu = phs.aim2.data$cardiac_arrest_grp_incdt
phs.aim2.data$cardiac_arrest_grp_inc_fu[which(phs.aim2.data$cardiac_arrest_grp_inc==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$cardiac_arrest_grp_inc==0)] 

phs.aim2.data$carotid_disease_grp_inc_fu = phs.aim2.data$carotid_disease_grp_incdt
phs.aim2.data$carotid_disease_grp_inc_fu[which(phs.aim2.data$carotid_disease_grp_inc==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$carotid_disease_grp_inc==0)] 

phs.aim2.data$myocarditis_pericarditis_grp_inc_fu = phs.aim2.data$myocarditis_pericarditis_grp_incdt
phs.aim2.data$myocarditis_pericarditis_grp_inc_fu[which(phs.aim2.data$myocarditis_pericarditis_grp_inc==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$myocarditis_pericarditis_grp_inc==0)] 

phs.aim2.data$tia_grp_inc_fu = phs.aim2.data$tia_grp_incdt
phs.aim2.data$tia_grp_inc_fu[which(phs.aim2.data$tia_grp_inc==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$tia_grp_inc==0)] 

phs.aim2.data$valvular_disease_grp_inc_fu = phs.aim2.data$valvular_disease_grp_incdt
phs.aim2.data$valvular_disease_grp_inc_fu[which(phs.aim2.data$valvular_disease_grp_inc==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$valvular_disease_grp_inc==0)] 

phs.aim2.data$venous_thromboembolic_disease_grp_inc_fu = phs.aim2.data$venous_thromboembolic_disease_grp_incdt
phs.aim2.data$venous_thromboembolic_disease_grp_inc_fu[which(phs.aim2.data$venous_thromboembolic_disease_grp_inc==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$venous_thromboembolic_disease_grp_inc==0)] 

# define any new onset = true incindence + recurrence
# combined event of ischemic heart disease, stroke/TIA, cardiomyopathy/heart failure
phs.aim2.data$cvdcombo_grp_rec = ifelse(phs.aim2.data$ischemic_heart_disease_grp_rec==1 | 
                                          phs.aim2.data$stroke_grp_rec==1 |  
                                          phs.aim2.data$cardiomyopathy_grp_rec==1|
                                          phs.aim2.data$heart_failure_grp_rec==1, 1, 0)


phs.aim2.data$cvdcombo_grp_recdt = apply(phs.aim2.data[,c("ischemic_heart_disease_grp_recdt",
                                                          "stroke_grp_recdt",
                                                          "cardiomyopathy_grp_recdt",
                                                          "heart_failure_grp_recdt")],1,min,na.rm=T)

# calculate time of follow up
phs.aim2.data$ischemic_heart_disease_grp_rec_fu = phs.aim2.data$ischemic_heart_disease_grp_recdt
phs.aim2.data$ischemic_heart_disease_grp_rec_fu[which(phs.aim2.data$ischemic_heart_disease_grp_rec==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$ischemic_heart_disease_grp_rec==0)] 

phs.aim2.data$stroke_grp_rec_fu = phs.aim2.data$stroke_grp_recdt
phs.aim2.data$stroke_grp_rec_fu[which(phs.aim2.data$stroke_grp_rec==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$stroke_grp_rec==0)]

phs.aim2.data$cardiomyopathy_grp_rec_fu = phs.aim2.data$cardiomyopathy_grp_recdt
phs.aim2.data$cardiomyopathy_grp_rec_fu[which(phs.aim2.data$cardiomyopathy_grp_rec==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$cardiomyopathy_grp_rec==0)] 

phs.aim2.data$heart_failure_grp_rec_fu = phs.aim2.data$heart_failure_grp_recdt
phs.aim2.data$heart_failure_grp_rec_fu[which(phs.aim2.data$heart_failure_grp_rec==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$heart_failure_grp_rec==0)] 

phs.aim2.data$cvdcombo_grp_rec_fu = phs.aim2.data$cvdcombo_grp_recdt
phs.aim2.data$cvdcombo_grp_rec_fu[which(phs.aim2.data$cvdcombo_grp_rec==0)] = 
  phs.aim2.data$censor_dt[which(phs.aim2.data$cvdcombo_grp_rec==0)] 


### Creat the combination of prevalent CVDs
phs.aim2.data$prevcvd = ifelse(phs.aim2.data$arrhythmia_grp_prev == 1 | phs.aim2.data$cardiomyopathy_grp_prev == 1 | 
                                 phs.aim2.data$heart_failure_grp_prev == 1 | phs.aim2.data$ischemic_heart_disease_grp_prev == 1 | 
                                 phs.aim2.data$myocarditis_pericarditis_grp_prev == 1 | phs.aim2.data$stroke_grp_prev == 1 | 
                                 phs.aim2.data$tia_grp_prev == 1 | phs.aim2.data$valvular_disease_grp_prev == 1 | 
                                 phs.aim2.data$venous_thromboembolic_disease_grp_prev == 1,1,0)



### Define CVD risk factor outcomes
# define incident diabetes, hypertension and dyslipidemia
phs.aim2.data$cvdrf_diab = ifelse(phs.aim2.data$daysto_diabetes>0,1,0)
phs.aim2.data$cvdrf_diab[is.na(phs.aim2.data$cvdrf_diab)] = 0
phs.aim2.data$cvdrf_htn = ifelse(phs.aim2.data$daysto_phase_htn>0,1,0)
phs.aim2.data$cvdrf_htn[is.na(phs.aim2.data$cvdrf_htn)] = 0
phs.aim2.data$cvdrf_dyslipid = ifelse(phs.aim2.data$daysto_dyslipidemia>0,1,0)
phs.aim2.data$cvdrf_dyslipid[is.na(phs.aim2.data$cvdrf_dyslipid)] = 0

# combined event of diabetes, hypertension or dyslipidemia
phs.aim2.data$cvdrfcombo = ifelse(phs.aim2.data$cvdrf_diab==1 | 
                                    phs.aim2.data$cvdrf_htn==1 | 
                                    phs.aim2.data$cvdrf_dyslipid==1, 1, 0)

phs.aim2.data$cvdrfcombo_dt = apply(phs.aim2.data[,c("daysto_diabetes","daysto_phase_htn",
                                                     "daysto_dyslipidemia")],1,function(x){
                                                       min(x[x>0], na.rm=T)
                                                     })
phs.aim2.data$cvdrfcombo_dt[phs.aim2.data$cvdrfcombo_dt==Inf] = NA

# define censoring time
phs.aim2.data$censor_dt = apply(phs.aim2.data[,c("daysto_death","daysto_enr_end","daysto_end_study")],1,min,na.rm=T)

# calculate time of follow up
phs.aim2.data$cvdrf_diab_fu = phs.aim2.data$daysto_diabetes
phs.aim2.data$cvdrf_diab_fu[which(phs.aim2.data$cvdrf_diab==0)] = phs.aim2.data$censor_dt[which(phs.aim2.data$cvdrf_diab==0)] 

phs.aim2.data$cvdrf_htn_fu = phs.aim2.data$daysto_phase_htn
phs.aim2.data$cvdrf_htn_fu[which(phs.aim2.data$cvdrf_htn==0)] = phs.aim2.data$censor_dt[which(phs.aim2.data$cvdrf_htn==0)] 

phs.aim2.data$cvdrf_dyslipid_fu = phs.aim2.data$daysto_dyslipidemia
phs.aim2.data$cvdrf_dyslipid_fu[which(phs.aim2.data$cvdrf_dyslipid==0)] = phs.aim2.data$censor_dt[which(phs.aim2.data$cvdrf_dyslipid==0)] 

phs.aim2.data$cvdrfcombo_fu = phs.aim2.data$cvdrfcombo_dt
phs.aim2.data$cvdrfcombo_fu[which(phs.aim2.data$cvdrfcombo==0)] = phs.aim2.data$censor_dt[which(phs.aim2.data$cvdrfcombo==0)] 

# relevel bmi var
phs.aim2.data$bmicat1 = factor(phs.aim2.data$bmicat, levels = c('Normal','Underweight','Overweight','Obese I','Obese II+','Unknown'))

# relevel baseline diabetes 
phs.aim2.data$diab_bl1 = factor(phs.aim2.data$diab_bl, levels = c('No','Yes'))

# relevel dyslipidemia
phs.aim2.data$dyslipid_bl1 = factor(phs.aim2.data$dyslipid_bl, levels=c('No/Unknown','Yes'))


# Create Chemotherapy Variables
phs.aim2.data$anthra = phs.aim2.data$anthra_yn
phs.aim2.data$anthra[phs.aim2.data$drug_comb == 9] = 0
phs.aim2.data$anthra_full = phs.aim2.data$anthra_full_yn
phs.aim2.data$anthra_full[phs.aim2.data$drug_comb == 9] = 0
phs.aim2.data$anthra_part = phs.aim2.data$anthra_part_yn
phs.aim2.data$anthra_part[phs.aim2.data$drug_comb == 9] = 0
phs.aim2.data$anthra = as.factor(phs.aim2.data$anthra)

phs.aim2.data$tras = phs.aim2.data$tras_yn
phs.aim2.data$tras[phs.aim2.data$drug_comb == 9] = 0
phs.aim2.data$tras_full = phs.aim2.data$tras_full_yn
phs.aim2.data$tras_full[phs.aim2.data$drug_comb == 9] = 0
phs.aim2.data$tras_part = phs.aim2.data$tras_part_yn
phs.aim2.data$tras_part[phs.aim2.data$drug_comb == 9] = 0
phs.aim2.data$tras = as.factor(phs.aim2.data$tras)

phs.aim2.data$taxane = phs.aim2.data$taxane_yn
phs.aim2.data$taxane[phs.aim2.data$drug_comb == 9] = 0
phs.aim2.data$taxane_full = phs.aim2.data$taxane_full_yn
phs.aim2.data$taxane_full[phs.aim2.data$drug_comb == 9] = 0
phs.aim2.data$taxane_part = phs.aim2.data$taxane_part_yn
phs.aim2.data$taxane_part[phs.aim2.data$drug_comb == 9] = 0
phs.aim2.data$taxane = as.factor(phs.aim2.data$taxane)

#	Anthracycline without trastuzumab, without taxane 
phs.aim2.data$comb1 = ifelse(phs.aim2.data$drug_comb == 1,1,0)

# Anthracycline without trastuzumab, with taxane 
phs.aim2.data$comb2 = ifelse(phs.aim2.data$drug_comb == 2,1,0)

# Anthracycline with trastuzumab, without taxane 
phs.aim2.data$comb3 = ifelse(phs.aim2.data$drug_comb == 3,1,0)

# Anthracycline with trastuzumab, with taxane 
phs.aim2.data$comb4 = ifelse(phs.aim2.data$drug_comb == 4,1,0)

# Trastuzumab without anthracycline, without taxane (aka Trastuzumab alone)
phs.aim2.data$comb5 = ifelse(phs.aim2.data$drug_comb == 5,1,0)

# Trastuzumab without anthracycline, with taxane 
phs.aim2.data$comb6 = ifelse(phs.aim2.data$drug_comb == 6,1,0)

# Other chemotherapy, with taxane 
phs.aim2.data$comb7 = ifelse(phs.aim2.data$drug_comb == 7,1,0)

# Other chemotherapy, without taxane 
phs.aim2.data$comb8 = ifelse(phs.aim2.data$drug_comb == 8,1,0)

# No chemotherapy
phs.aim2.data$comb9 = ifelse(phs.aim2.data$drug_comb == 9,1,0)

#	Anthracycline with trastuzumab
phs.aim2.data$com10 = ifelse(phs.aim2.data$anthra == 1 & phs.aim2.data$tras == 1,1,0)

# Hormonal Therapy variables
phs.aim2.data$tmt_type = horm$tmt_type[match(phs.aim2.data$cvd_studyid, horm$cvd_studyid)]
phs.aim2.data$ai = NA
phs.aim2.data$ai[phs.aim2.data$tmt_type == "AI"] = 1
phs.aim2.data$ai[phs.aim2.data$tmt_type == "TAM"] = 0
phs.aim2.data$ai[is.na(phs.aim2.data$ai)] = 0
phs.aim2.data$ai = as.factor(phs.aim2.data$ai)
phs.aim2.data$tam = NA
phs.aim2.data$tam[phs.aim2.data$tmt_type == "TAM"] = 1
phs.aim2.data$tam[phs.aim2.data$tmt_type == "AI"] = 0
phs.aim2.data$tam[is.na(phs.aim2.data$tam)] = 0
phs.aim2.data$tam = as.factor(phs.aim2.data$tam)

phs.aim2.data$horm_any = ifelse(phs.aim2.data$ai == 1 | phs.aim2.data$tam == 1,1,0)
phs.aim2.data$horm_any = as.factor(phs.aim2.data$horm_any)

# Radiation Therapy variable
phs.aim2.data$rad = as.character(phs.aim2.data$rad_tx_yn)
phs.aim2.data$rad[phs.aim2.data$rad == "Yes"] = 1
phs.aim2.data$rad[phs.aim2.data$rad == "No"] = 0
phs.aim2.data$rad = as.factor(phs.aim2.data$rad)

# Any anthracycline with trastuzumab, with endocrine therapy, with radiation therapy
phs.aim2.data$comb11 = ifelse(phs.aim2.data$anthra == 1 & phs.aim2.data$tras == 1 & phs.aim2.data$horm_any == 1 & phs.aim2.data$rad == 1,1,0)

# LVEF variables
phs.aim2.data$final_EF = ef$final_EF[match(phs.aim2.data$cvd_studyid,ef$CVD_studyid)]
phs.aim2.data$lvef_ind = ifelse(!is.na(phs.aim2.data$final_EF), 1, 0)
phs.aim2.data$lvef_ind = as.factor(phs.aim2.data$lvef_ind)
phs.aim2.data$lvef = phs.aim2.data$final_EF
phs.aim2.data$lvef[is.na(phs.aim2.data$lvef)] = mean(phs.aim2.data$final_EF, na.rm = T)


# Charlson Score
phs.aim2.data$charlson = charlson$CHRLSON[match(phs.aim2.data$cvd_studyid, charlson$CVD_studyid)]
phs.aim2.data$charlson = as.factor(phs.aim2.data$charlson)



# all
a2_ihd = which(phs.aim2.data$ischemic_heart_disease_grp_prev==0)
a2_hf_cm = which(phs.aim2.data$heart_failure_cardiomyopathy_grp_prev==0)
a2_stroke = which(phs.aim2.data$stroke_grp_prev==0)
a2_combo = which(phs.aim2.data$cvdcombo_grp_prev==0)
a2_arrhythmia = which(phs.aim2.data$arrhythmia_grp_prev==0)
a2_cardiac = which(phs.aim2.data$cardiac_arrest_grp_prev==0)
a2_carotid = which(phs.aim2.data$carotid_disease_grp_prev==0)
a2_myocarditis = which(phs.aim2.data$myocarditis_pericarditis_grp_prev==0)
a2_tia = which(phs.aim2.data$tia_grp_prev==0)
a2_valvular = which(phs.aim2.data$valvular_disease_grp_prev==0)
a2_dvt = which(phs.aim2.data$venous_thromboembolic_disease_grp_prev==0)
a2_diab = which(phs.aim2.data$diab_bl=='No')
a2_htn = which(phs.aim2.data$htn_bl=='No')
a2_dyslipid = which(phs.aim2.data$dyslipid_bl=='No')
a2_rfcombo = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$htn_bl=='No' & phs.aim2.data$dyslipid_bl=='No')

# excluding 110 prevalent cvd at baseline
a2_ihd_cvd = which(phs.aim2.data$ischemic_heart_disease_grp_prev==0 & phs.aim2.data$prevcvd != 1)
a2_hf_cm_cvd = which(phs.aim2.data$heart_failure_cardiomyopathy_grp_prev==0 & phs.aim2.data$prevcvd != 1)
a2_stroke_cvd = which(phs.aim2.data$stroke_grp_prev==0 & phs.aim2.data$prevcvd != 1)
a2_combo_cvd = which(phs.aim2.data$cvdcombo_grp_prev==0 & phs.aim2.data$prevcvd != 1)
a2_arrhythmia_cvd = which(phs.aim2.data$arrhythmia_grp_prev==0 & phs.aim2.data$prevcvd != 1)
a2_cardiac_cvd = which(phs.aim2.data$cardiac_arrest_grp_prev==0 & phs.aim2.data$prevcvd != 1)
a2_carotid_cvd = which(phs.aim2.data$carotid_disease_grp_prev==0 & phs.aim2.data$prevcvd != 1)
a2_myocarditis_cvd = which(phs.aim2.data$myocarditis_pericarditis_grp_prev==0 & phs.aim2.data$prevcvd != 1)
a2_tia_cvd = which(phs.aim2.data$tia_grp_prev==0 & phs.aim2.data$prevcvd != 1)
a2_valvular_cvd = which(phs.aim2.data$valvular_disease_grp_prev==0 & phs.aim2.data$prevcvd != 1)
a2_dvt_cvd = which(phs.aim2.data$venous_thromboembolic_disease_grp_prev==0 & phs.aim2.data$prevcvd != 1)
a2_diab_cvd = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$prevcvd != 1)
a2_htn_cvd = which(phs.aim2.data$htn_bl=='No' & phs.aim2.data$prevcvd != 1)
a2_dyslipid_cvd= which(phs.aim2.data$dyslipid_bl=='No' & phs.aim2.data$prevcvd != 1)
a2_rfcombo_cvd = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$htn_bl=='No' & phs.aim2.data$dyslipid_bl=='No' & phs.aim2.data$prevcvd != 1)



# left sided radiation therapy
id_l = phs.aim2.data$cvd_studyid[which(phs.aim2.data$laterality==2)]
a2_l = which(phs.aim2.data$cvd_studyid %in% id_l)
a2_ihd_l = which(phs.aim2.data$ischemic_heart_disease_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l)
a2_hf_cm_l = which(phs.aim2.data$heart_failure_cardiomyopathy_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l)
a2_stroke_l = which(phs.aim2.data$stroke_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l)
a2_combo_l = which(phs.aim2.data$cvdcombo_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l)
a2_arrhythmia_l = which(phs.aim2.data$arrhythmia_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l)
a2_cardiac_l = which(phs.aim2.data$cardiac_arrest_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l)
a2_carotid_l = which(phs.aim2.data$carotid_disease_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l)
a2_myocarditis_l = which(phs.aim2.data$myocarditis_pericarditis_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l)
a2_tia_l = which(phs.aim2.data$tia_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l)
a2_valvular_l = which(phs.aim2.data$valvular_disease_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l)
a2_dvt_l = which(phs.aim2.data$venous_thromboembolic_disease_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l)
a2_diab_l = which(phs.aim2.data$diab_bl=='No'& phs.aim2.data$cvd_studyid %in% id_l)
a2_htn_l = which(phs.aim2.data$htn_bl=='No'& phs.aim2.data$cvd_studyid %in% id_l)
a2_dyslipid_l = which(phs.aim2.data$dyslipid_bl=='No'& phs.aim2.data$cvd_studyid %in% id_l)
a2_rfcombo_l = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$htn_bl=='No' & phs.aim2.data$dyslipid_bl=='No' & phs.aim2.data$cvd_studyid %in% id_l)


# left sided radiation therapy excluding 110 prevalent cvd at baseline
id_l_cvd = phs.aim2.data$cvd_studyid[which(phs.aim2.data$laterality==2)]
a2_l_cvd = which(phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_ihd_l_cvd = which(phs.aim2.data$ischemic_heart_disease_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_hf_cm_l_cvd = which(phs.aim2.data$heart_failure_cardiomyopathy_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_stroke_l_cvd = which(phs.aim2.data$stroke_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_combo_l_cvd = which(phs.aim2.data$cvdcombo_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_arrhythmia_l_cvd = which(phs.aim2.data$arrhythmia_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_cardiac_l_cvd = which(phs.aim2.data$cardiac_arrest_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_carotid_l_cvd = which(phs.aim2.data$carotid_disease_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_myocarditis_l_cvd = which(phs.aim2.data$myocarditis_pericarditis_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_tia_l_cvd = which(phs.aim2.data$tia_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_valvular_l_cvd = which(phs.aim2.data$valvular_disease_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_dvt_l_cvd = which(phs.aim2.data$venous_thromboembolic_disease_grp_prev==0 & phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1) 
a2_diab_l_cvd = which(phs.aim2.data$diab_bl=='No'& phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_htn_l_cvd = which(phs.aim2.data$htn_bl=='No'& phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_dyslipid_l_cvd = which(phs.aim2.data$dyslipid_bl=='No'& phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)
a2_rfcombo_l_cvd = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$htn_bl=='No' & phs.aim2.data$dyslipid_bl=='No' & phs.aim2.data$cvd_studyid %in% id_l & phs.aim2.data$prevcvd != 1)



# LVEF > 50 
a2_ihd_ef = which(phs.aim2.data$ischemic_heart_disease_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_hf_cm_ef = which(phs.aim2.data$heart_failure_cardiomyopathy_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_stroke_ef = which(phs.aim2.data$stroke_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_combo_ef = which(phs.aim2.data$cvdcombo_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_arrhythmia_ef = which(phs.aim2.data$arrhythmia_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_cardiac_ef = which(phs.aim2.data$cardiac_arrest_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_carotid_ef = which(phs.aim2.data$carotid_disease_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_myocarditis_ef = which(phs.aim2.data$myocarditis_pericarditis_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_tia_ef = which(phs.aim2.data$tia_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_valvular_ef = which(phs.aim2.data$valvular_disease_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_dvt_ef = which(phs.aim2.data$venous_thromboembolic_disease_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_diab_ef = which(phs.aim2.data$diab_bl=='No' & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_htn_ef = which(phs.aim2.data$htn_bl=='No' & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_dyslipid_ef = which(phs.aim2.data$dyslipid_bl=='No' & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))
a2_rfcombo_ef = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$htn_bl=='No' & phs.aim2.data$dyslipid_bl=='No' & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)))


# LVEF > 50 excluding 110 prevalent cvd at baseline
a2_ihd_ef_cvd = which(phs.aim2.data$ischemic_heart_disease_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_hf_cm_ef_cvd = which(phs.aim2.data$heart_failure_cardiomyopathy_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_stroke_ef_cvd = which(phs.aim2.data$stroke_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_combo_ef_cvd = which(phs.aim2.data$cvdcombo_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_arrhythmia_ef_cvd = which(phs.aim2.data$arrhythmia_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_cardiac_ef_cvd = which(phs.aim2.data$cardiac_arrest_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_carotid_ef_cvd = which(phs.aim2.data$carotid_disease_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_myocarditis_ef_cvd = which(phs.aim2.data$myocarditis_pericarditis_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_tia_ef_cvd = which(phs.aim2.data$tia_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_valvular_ef_cvd = which(phs.aim2.data$valvular_disease_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_dvt_ef_cvd = which(phs.aim2.data$venous_thromboembolic_disease_grp_prev==0 & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_diab_ef_cvd = which(phs.aim2.data$diab_bl=='No' & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_htn_ef_cvd = which(phs.aim2.data$htn_bl=='No' & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_dyslipid_ef_cvd = which(phs.aim2.data$dyslipid_bl=='No' & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)
a2_rfcombo_ef_cvd = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$htn_bl=='No' & phs.aim2.data$dyslipid_bl=='No' & (phs.aim2.data$final_EF > 50 | is.na(phs.aim2.data$final_EF)) & phs.aim2.data$prevcvd != 1)


# post menopausal
a2_ihd_post_menop = which(phs.aim2.data$ischemic_heart_disease_grp_prev==0 & phs.aim2.data$menop == 1)
a2_hf_cm_post_menop = which(phs.aim2.data$heart_failure_cardiomyopathy_grp_prev==0 & phs.aim2.data$menop == 1)
a2_stroke_post_menop = which(phs.aim2.data$stroke_grp_prev==0 & phs.aim2.data$menop == 1)
a2_combo_post_menop = which(phs.aim2.data$cvdcombo_grp_prev==0 & phs.aim2.data$menop == 1)
a2_arrhythmia_post_menop = which(phs.aim2.data$arrhythmia_grp_prev==0 & phs.aim2.data$menop == 1)
a2_cardiac_post_menop = which(phs.aim2.data$cardiac_arrest_grp_prev==0 & phs.aim2.data$menop == 1)
a2_carotid_post_menop = which(phs.aim2.data$carotid_disease_grp_prev==0 & phs.aim2.data$menop == 1)
a2_myocarditis_post_menop = which(phs.aim2.data$myocarditis_pericarditis_grp_prev==0 & phs.aim2.data$menop == 1)
a2_tia_post_menop = which(phs.aim2.data$tia_grp_prev==0 & phs.aim2.data$menop == 1)
a2_valvular_post_menop = which(phs.aim2.data$valvular_disease_grp_prev==0 & phs.aim2.data$menop == 1)
a2_dvt_post_menop = which(phs.aim2.data$venous_thromboembolic_disease_grp_prev==0 & phs.aim2.data$menop == 1)
a2_diab_post_menop = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$menop == 1)
a2_htn_post_menop = which(phs.aim2.data$htn_bl=='No' & phs.aim2.data$menop == 1)
a2_dyslipid_post_menop = which(phs.aim2.data$dyslipid_bl=='No' & phs.aim2.data$menop == 1)
a2_rfcombo_post_menop = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$htn_bl=='No' & phs.aim2.data$dyslipid_bl=='No' & phs.aim2.data$menop == 1)

# post menopausal excluding 110 prevalent cvd at baseline
a2_ihd_post_menop_cvd = which(phs.aim2.data$ischemic_heart_disease_grp_prev==0 & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_hf_cm_post_menop_cvd = which(phs.aim2.data$heart_failure_cardiomyopathy_grp_prev==0 & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_stroke_post_menop_cvd = which(phs.aim2.data$stroke_grp_prev==0 & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_combo_post_menop_cvd = which(phs.aim2.data$cvdcombo_grp_prev==0 & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_arrhythmia_post_menop_cvd = which(phs.aim2.data$arrhythmia_grp_prev==0 & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_cardiac_post_menop_cvd = which(phs.aim2.data$cardiac_arrest_grp_prev==0 & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_carotid_post_menop_cvd = which(phs.aim2.data$carotid_disease_grp_prev==0 & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_myocarditis_post_menop_cvd = which(phs.aim2.data$myocarditis_pericarditis_grp_prev==0 & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_tia_post_menop_cvd = which(phs.aim2.data$tia_grp_prev==0 & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_valvular_post_menop_cvd = which(phs.aim2.data$valvular_disease_grp_prev==0 & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_dvt_post_menop_cvd = which(phs.aim2.data$venous_thromboembolic_disease_grp_prev==0 & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_diab_post_menop_cvd = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_htn_post_menop_cvd = which(phs.aim2.data$htn_bl=='No' & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_dyslipid_post_menop_cvd = which(phs.aim2.data$dyslipid_bl=='No' & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)
a2_rfcombo_post_menop_cvd = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$htn_bl=='No' & phs.aim2.data$dyslipid_bl=='No' & phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1)



# pre menopausal
a2_ihd_pre_menop = which(phs.aim2.data$ischemic_heart_disease_grp_prev==0 & phs.aim2.data$menop == 0)
a2_hf_cm_pre_menop = which(phs.aim2.data$heart_failure_cardiomyopathy_grp_prev==0 & phs.aim2.data$menop == 0)
a2_stroke_pre_menop = which(phs.aim2.data$stroke_grp_prev==0 & phs.aim2.data$menop == 0)
a2_combo_pre_menop = which(phs.aim2.data$cvdcombo_grp_prev==0 & phs.aim2.data$menop == 0)
a2_arrhythmia_pre_menop = which(phs.aim2.data$arrhythmia_grp_prev==0 & phs.aim2.data$menop == 0)
a2_cardiac_pre_menop = which(phs.aim2.data$cardiac_arrest_grp_prev==0 & phs.aim2.data$menop == 0)
a2_carotid_pre_menop = which(phs.aim2.data$carotid_disease_grp_prev==0 & phs.aim2.data$menop == 0)
a2_myocarditis_pre_menop = which(phs.aim2.data$myocarditis_pericarditis_grp_prev==0 & phs.aim2.data$menop == 0)
a2_tia_pre_menop = which(phs.aim2.data$tia_grp_prev==0 & phs.aim2.data$menop == 0)
a2_valvular_pre_menop = which(phs.aim2.data$valvular_disease_grp_prev==0 & phs.aim2.data$menop == 0)
a2_dvt_pre_menop = which(phs.aim2.data$venous_thromboembolic_disease_grp_prev==0 & phs.aim2.data$menop == 0)
a2_diab_pre_menop = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$menop == 0)
a2_htn_pre_menop = which(phs.aim2.data$htn_bl=='No' & phs.aim2.data$menop == 0)
a2_dyslipid_pre_menop = which(phs.aim2.data$dyslipid_bl=='No' & phs.aim2.data$menop == 0)
a2_rfcombo_pre_menop = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$htn_bl=='No' & phs.aim2.data$dyslipid_bl=='No' & phs.aim2.data$menop == 0)

# pre menopausal excluding 110 prevalent cvd at baseline
a2_ihd_pre_menop_cvd = which(phs.aim2.data$ischemic_heart_disease_grp_prev==0 & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_hf_cm_pre_menop_cvd = which(phs.aim2.data$heart_failure_cardiomyopathy_grp_prev==0 & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_stroke_pre_menop_cvd = which(phs.aim2.data$stroke_grp_prev==0 & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_combo_pre_menop_cvd = which(phs.aim2.data$cvdcombo_grp_prev==0 & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_arrhythmia_pre_menop_cvd = which(phs.aim2.data$arrhythmia_grp_prev==0 & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_cardiac_pre_menop_cvd = which(phs.aim2.data$cardiac_arrest_grp_prev==0 & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_carotid_pre_menop_cvd = which(phs.aim2.data$carotid_disease_grp_prev==0 & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_myocarditis_pre_menop_cvd = which(phs.aim2.data$myocarditis_pericarditis_grp_prev==0 & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_tia_pre_menop_cvd = which(phs.aim2.data$tia_grp_prev==0 & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_valvular_pre_menop_cvd = which(phs.aim2.data$valvular_disease_grp_prev==0 & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_dvt_pre_menop_cvd = which(phs.aim2.data$venous_thromboembolic_disease_grp_prev==0 & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_diab_pre_menop_cvd = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_htn_pre_menop_cvd = which(phs.aim2.data$htn_bl=='No' & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_dyslipid_pre_menop_cvd = which(phs.aim2.data$dyslipid_bl=='No' & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)
a2_rfcombo_pre_menop_cvd = which(phs.aim2.data$diab_bl=='No' & phs.aim2.data$htn_bl=='No' & phs.aim2.data$dyslipid_bl=='No' & phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1)


# Create some subsets of Aim 2 dataset
a2_ef = phs.aim2.data[phs.aim2.data$final_EF > 50 & !is.na(phs.aim2.data$final_EF),]
a2_post_menop = phs.aim2.data[phs.aim2.data$menop == 1,]
a2_pre_menop = phs.aim2.data[phs.aim2.data$menop == 0,]
a2_rad_l = phs.aim2.data[phs.aim2.data$laterality==2,]

# Create some subsets of Aim 2 dataset excluding lvef <= 50
phs.aim2.data.lvefgr50 = phs.aim2.data[phs.aim2.data$final_EF > 50,]
a2_ef_lvefgr50 = phs.aim2.data[phs.aim2.data$final_EF > 50,]

# Create some subsets of Aim 2 dataset excluding 110 prevalent CVD at baseline
phs.aim2.data.cvd = phs.aim2.data[phs.aim2.data$prevcvd != 1,]
a2_ef_cvd = phs.aim2.data[phs.aim2.data$final_EF > 50 & phs.aim2.data$prevcvd != 1,]
a2_post_menop_cvd = phs.aim2.data[phs.aim2.data$menop == 1 & phs.aim2.data$prevcvd != 1,]
a2_pre_menop_cvd = phs.aim2.data[phs.aim2.data$menop == 0 & phs.aim2.data$prevcvd != 1,]
a2_rad_l_cvd = phs.aim2.data[phs.aim2.data$laterality==2 & phs.aim2.data$prevcvd != 1,]


##########################################################################################################################
#####################################################################
################# Generate IPW Method Datasets  #####################
#####################################################################
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
  dat1 = as.data.table(dat1)
  
  dat2 = data.frame(ID = d[,"cvd_studyid"])
  trt_strat = paste0(trt, "_start_dt")
  trt_end = paste0(trt, "_end_dt")
  dat2$startA = d[,trt_strat]
  dat2$endA = d[,trt_end]
  dat2 = dat2[!is.na(dat2$startA),]
  dat2 = as.data.table(dat2)
  

  list_cohort0 = list("categorical"=TRUE,
                     "impute"=NA,
                     "impute_default_level"=NA)
  list_cohort = rep(list(list_cohort0), length(covar))
  names(list_cohort) = covar
  
  cohort = setCohort(dat1, "ID", "IndexDate", "EOFDate", "EOFtype",
                     x, covar, list_cohort) 
  
  exposure = setExposure(dat2, "ID", "startA", "endA")
  
  LtAt.specification = cohort + exposure
  LtAt.data = construct(LtAt.specification, time_unit = 365, first_exp_rule = 1,
                         exp_threshold = 0.75)
  
  
  return(list(dat1 = dat1, dat2 = dat2, LtAt.data = LtAt.data))
}


# Individual Chemotherapy
covar_chemo = c("ajcc_stage", "diab_bl", "htn_bl", "dyslipid_bl", "bmicat1", "smok", "charlson", "edu_cat", "income_cat",
                "agegrp", "race", "horm_yn","rad_tx_yn", "menop")

phs.aim2.data$eofdate = pmin(phs.aim2.data$death_date, phs.aim2.data$end_of_study, phs.aim2.data$daysto_cvd_dt, phs.aim2.data$enr_end, na.rm = T)
phs.aim2.data$eofdate[phs.aim2.data$eofdate < phs.aim2.data$index_date] = phs.aim2.data$index_date[phs.aim2.data$eofdate < phs.aim2.data$index_date]

phs.aim2.data$eoftype = NA
phs.aim2.data$eoftype[!is.na(phs.aim2.data$death_date)] = "Death"
phs.aim2.data$eoftype[phs.aim2.data$eofdate == "2018-12-31"] = "End of Study"
phs.aim2.data$eoftype[phs.aim2.data$ischemic_heart_disease_event == 1] = "Ischemic Heart Disease"
phs.aim2.data$eoftype[phs.aim2.data$heart_failure_event == 1] = "Heart Failure"
phs.aim2.data$eoftype[phs.aim2.data$cardiomyopathy_event == 1] = "Cardiomyopathy"
phs.aim2.data$eoftype[phs.aim2.data$stroke_event == 1] = "Stroke"
phs.aim2.data$eoftype[phs.aim2.data$arrhythmia_event == 1] = "Arrhythmia"
phs.aim2.data$eoftype[phs.aim2.data$cardiac_arrest_event == 1] = "Cardiac Arrest"
phs.aim2.data$eoftype[phs.aim2.data$carotid_disease_event == 1] = "Carotid Disease"
phs.aim2.data$eoftype[phs.aim2.data$myocarditis_pericarditis_event == 1] = "Myocarditis Pericarditis"
phs.aim2.data$eoftype[phs.aim2.data$tia_event == 1] = "TIA"
phs.aim2.data$eoftype[phs.aim2.data$valvular_disease_event == 1] = "Valvular Disease"
phs.aim2.data$eoftype[phs.aim2.data$vte_event == 1] = "VTE"
phs.aim2.data$eoftype[is.na(phs.aim2.data$eoftype)] = "Disenrollment"

# Taxane
# Primary Outcomes
dat_ihd_taxane = create_dat_tab(phs.aim2.data, "Ischemic Heart Disease", "eoftype", covar_chemo, "taxane")
dat_hf_taxane = create_dat_tab(phs.aim2.data, "Heart Failure", "eoftype", covar_chemo, "taxane")
dat_cardio_taxane = create_dat_tab(phs.aim2.data, "Cardiomyopathy", "eoftype", covar_chemo, "taxane")
dat_stroke_taxane = create_dat_tab(phs.aim2.data, "Stroke", "eoftype", covar_chemo, "taxane")













