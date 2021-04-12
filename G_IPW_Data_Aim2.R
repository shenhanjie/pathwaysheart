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




phs.aim2.data$eofdate = pmin(phs.aim2.data$death_date, phs.aim2.data$end_of_study, phs.aim2.data$daysto_cvd_dt, phs.aim2.data$enr_end, na.rm = T)
phs.aim2.data$eofdate[phs.aim2.data$eofdate < phs.aim2.data$index_date] = phs.aim2.data$index_date[phs.aim2.data$eofdate < phs.aim2.data$index_date]

phs.aim2.data$eoftype = NA
phs.aim2.data$eoftype[!is.na(phs.aim2.data$death_date)] = "Death"
phs.aim2.data$eoftype[phs.aim2.data$eofdate == "2018-12-31"] = "End of Study"
phs.aim2.data$eoftype[phs.aim2.data$ischemic_heart_disease_event == 1] = "Ischemic Heart Disease"
phs.aim2.data$eoftype[phs.aim2.data$heart_failure_event == 1 | phs.aim2.data$cardiomyopathy_event == 1] = "Heart Failure/Cardiomyopathy"
#phs.aim2.data$eoftype[phs.aim2.data$cardiomyopathy_event == 1] = "Cardiomyopathy"
phs.aim2.data$eoftype[phs.aim2.data$stroke_event == 1] = "Stroke"
phs.aim2.data$eoftype[phs.aim2.data$arrhythmia_event == 1] = "Arrhythmia"
phs.aim2.data$eoftype[phs.aim2.data$cardiac_arrest_event == 1] = "Cardiac Arrest"
phs.aim2.data$eoftype[phs.aim2.data$carotid_disease_event == 1] = "Carotid Disease"
phs.aim2.data$eoftype[phs.aim2.data$myocarditis_pericarditis_event == 1] = "Myocarditis Pericarditis"
phs.aim2.data$eoftype[phs.aim2.data$tia_event == 1] = "TIA"
phs.aim2.data$eoftype[phs.aim2.data$valvular_disease_event == 1] = "Valvular Disease"
phs.aim2.data$eoftype[phs.aim2.data$vte_event == 1] = "VTE"
phs.aim2.data$eoftype[is.na(phs.aim2.data$eoftype)] = "Disenrollment"



# Chemotherapy
covar_chemo = c("ajcc_stage", "diab_bl", "htn_bl", "dyslipid_bl", "bmicat1", "smok", "charlson", "edu_cat", "income_cat",
                "agegrp", "race", "horm_yn","rad_tx_yn", "menop")

# Anthracycline
# Primary Outcomes
dat_ihd_anthra = create_dat_tab(phs.aim2.data[a2_ihd,], "Ischemic Heart Disease", "eoftype", covar_chemo, "anthra")
dat_hf_cm_anthra = create_dat_tab(phs.aim2.data[a2_hf_cm,], "Heart Failure/Cardiomyopathy", "eoftype", covar_chemo, "anthra")
dat_stroke_anthra = create_dat_tab(phs.aim2.data[a2_stroke,], "Stroke", "eoftype", covar_chemo, "anthra")
#write.csv(dat_ihd_anthra$LtAt.data,"dat_ihd_anthra.csv")
#write.csv(dat_hf_cm_anthra$LtAt.data,"dat_hf_cm_anthra.csv")
#write.csv(dat_stroke_anthra$LtAt.data,"dat_stroke_anthra.csv")


# Trastuzumab
dat_ihd_tras = create_dat_tab(phs.aim2.data[a2_ihd,], "Ischemic Heart Disease", "eoftype", covar_chemo, "tras")
dat_hf_cm_tras = create_dat_tab(phs.aim2.data[a2_hf_cm,], "Heart Failure/Cardiomyopathy", "eoftype", covar_chemo, "tras")
dat_stroke_tras = create_dat_tab(phs.aim2.data[a2_stroke,], "Stroke", "eoftype", covar_chemo, "tras")
#write.csv(dat_ihd_tras$LtAt.data,"dat_ihd_tras.csv")
#write.csv(dat_hf_cm_tras$LtAt.data,"dat_hf_cm_tras.csv")
#write.csv(dat_stroke_tras$LtAt.data,"dat_stroke_tras.csv")

# Taxane
# Primary Outcomes
dat_ihd_taxane = create_dat_tab(phs.aim2.data[a2_ihd,], "Ischemic Heart Disease", "eoftype", covar_chemo, "taxane")
dat_hf_cm_taxane = create_dat_tab(phs.aim2.data[a2_hf_cm,], "Heart Failure/Cardiomyopathy", "eoftype", covar_chemo, "taxane")
dat_stroke_taxane = create_dat_tab(phs.aim2.data[a2_stroke,], "Stroke", "eoftype", covar_chemo, "taxane")
#write.csv(dat_ihd_taxane$LtAt.data,"dat_ihd_taxane.csv")
#write.csv(dat_hf_cm_taxane$LtAt.data,"dat_hf_cm_taxane.csv")
#write.csv(dat_stroke_taxane$LtAt.data,"dat_stroke_taxane.csv")


# Hormonal Therapy
covar_horm = c("ajcc_stage", "diab_bl", "htn_bl", "dyslipid_bl", "bmicat1", "smok", "charlson", "edu_cat", "income_cat",
               "agegrp", "race", "chemo_yn","rad_tx_yn")

# ### Aromatase inhibitor (AI) (post menopausal women)
# Primary Outcomes
dat_ihd_ai = create_dat_tab(phs.aim2.data[a2_ihd_post_menop,], "Ischemic Heart Disease", "eoftype", covar_horm, "ai")
dat_hf_cm_ai = create_dat_tab(phs.aim2.data[a2_hf_cm_post_menop,], "Heart Failure/Cardiomyopathy", "eoftype", covar_horm, "ai")
dat_stroke_ai = create_dat_tab(phs.aim2.data[a2_stroke_post_menop,], "Stroke", "eoftype", covar_horm, "ai")

# Tamoxifen (TAM) (pre menopausal women)
dat_ihd_tam = create_dat_tab(phs.aim2.data[a2_ihd_pre_menop,], "Ischemic Heart Disease", "eoftype", covar_horm, "tam")
dat_hf_cm_tam = create_dat_tab(phs.aim2.data[a2_hf_cm_pre_menop,], "Heart Failure/Cardiomyopathy", "eoftype", covar_horm, "tam")
dat_stroke_tam = create_dat_tab(phs.aim2.data[a2_stroke_pre_menop,], "Stroke", "eoftype", covar_horm, "tam")

### Any hormonal therapy
dat_ihd_horm_any = create_dat_tab(phs.aim2.data[a2_ihd,], "Ischemic Heart Disease", "eoftype", covar_horm, "horm_any")
dat_hf_cm_horm_any = create_dat_tab(phs.aim2.data[a2_hf_cm,], "Heart Failure/Cardiomyopathy", "eoftype", covar_horm, "horm_any")
dat_stroke_horm_any = create_dat_tab(phs.aim2.data[a2_stroke,], "Stroke", "eoftype", covar_horm, "horm_any")


# Radiation Therapy
# Primary Outcomes
dat_ihd_rad = create_dat_tab(phs.aim2.data[a2_ihd,], "Ischemic Heart Disease", "eoftype", covar_chemo, "rad")
dat_hf_cm_rad = create_dat_tab(phs.aim2.data[a2_hf_cm,], "Heart Failure/Cardiomyopathy", "eoftype", covar_chemo, "rad")
dat_stroke_rad = create_dat_tab(phs.aim2.data[a2_stroke,], "Stroke", "eoftype", covar_chemo, "rad")


# Left-sided Radiation Therapy
# Primary Outcomes
dat_ihd_rad_l= create_dat_tab(phs.aim2.data[a2_ihd_l,], "Ischemic Heart Disease", "eoftype", covar_chemo, "rad")
dat_hf_cm_rad_l = create_dat_tab(phs.aim2.data[a2_hf_cm_l,], "Heart Failure/Cardiomyopathy", "eoftype", covar_chemo, "rad")
dat_stroke_rad_l = create_dat_tab(phs.aim2.data[a2_stroke_l,], "Stroke", "eoftype", covar_chemo, "rad")




