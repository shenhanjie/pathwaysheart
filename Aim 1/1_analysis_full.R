# Pathways Heart Study Aim 1 data analysis
# Hanjie Shen, 3/10/2020


library(tidyverse)
library(haven)
library(ggplot2)
library(dplyr)
library(plyr)
library(naniar)
library(reshape2)
library(survival)
library(survminer)
library(cmprsk)

# define function round2() to round number to exact decimal place
round2 <- function(x, k) trimws(format(round(x, k), nsmall=k))


setwd("~/Desktop/Pathways Heart Study/Aim 1 Results")


################################################################################
# compare characteristics between cases and controls

tab1 <- function(d,x){
  if(class(d[,x]) %in% c('numeric')){
    n <- c(paste0(length(which(!is.na(d[,x]))),' (',
                  round(100*length(which(is.na(d[,x])))/nrow(d)),'%)'),
           tapply(d[,x], d$group, function(y) 
             paste0(length(which(!is.na(y))),' (',
                    round(100*length(which(is.na(y)))/length(y)),'%)')))
    m <- c(mean(d[,x], na.rm=T),
           tapply(d[,x], d$group,mean, na.rm=T))
    sd <- c(sd(d[,x], na.rm=T),
            tapply(d[,x], d$group,sd, na.rm=T))
    p <- t.test(d[,x]~d$group)$p.value
    sum <- c(x,n[1],m[1],sd[1],'',n[2],m[2],sd[2],'',n[3],m[3],sd[3],'',p)
  } 
  if(class(d[,x]) %in% c('factor')){
    n <- cbind(table(d[,x]),table(d[,x], d$group))
    pct <- cbind(prop.table(table(d[,x])),
                 prop.table(table(d[,x], d$group),2))
    p <- chisq.test(table(d[,x], d$group))$p.value
    sum <- cbind(x,'',n[,1],pct[,1],'','',n[,2],pct[,2],'','',n[,3],pct[,3],'',p)
  } 
  sum
}

# variables to analyze
varlist <- c('dxage', 'enr_len','cops2','bmi1','agegrp','raceethn1',
             'bmicat','smok','menop','diab_bl','htn_bl','dyslipid_bl',
             "systolic" ,"diastolic","glu_f","hdl","hgba1c",
             "ldl_clc_ns","tot_choles","trigl_ns",
             'medhousincome','houspoverty', "edu1","edu2","edu3","edu4")

# all sample
table1_all <- do.call(rbind,lapply(varlist,function(x) tab1(a1,x)))
table1_all <- data.frame(cbind(row.names(table1_all), table1_all))

# cases receiving chemo
table1_chemo <- do.call(rbind,lapply(varlist,function(x) tab1(a1[a1_chemo,],x)))
table1_chemo <- data.frame(cbind(row.names(table1_chemo), table1_chemo))

# cases receiving hormonal therapy
table1_horm <- do.call(rbind,lapply(varlist,function(x) tab1(a1[a1_horm,],x)))
table1_horm <- data.frame(cbind(row.names(table1_horm), table1_horm))

# cases receiving radiation
table1_rad <- do.call(rbind,lapply(varlist,function(x) tab1(a1[a1_rad,],x)))
table1_rad <- data.frame(cbind(row.names(table1_rad), table1_rad))

# export tables
write_csv(table1_all,'table1_all.csv')
write_csv(table1_chemo,'table1_chemo.csv')
write_csv(table1_horm,'table1_horm.csv')
write_csv(table1_rad,'table1_rad.csv')







################################################################################
# Compare cvd prevalence between cases and controls

# get cvd prevalence variables
cvdvars <- grep('_prev$',names(a1),value = T)
cvdvars <- cvdvars[c(1:3,5,4,6:11,c(16,17,18,21,19,13:15,20,24,26,27,28,23,12,22,25,29:33))]
# function to create
tab2 <- function(d,x){
  d[,x] <- factor(d[,x],levels=c(0,1),labels = c('No','Yes'))
  n <- cbind(table(d[,x]),table(d[,x], d$group))
  pct <- cbind(prop.table(table(d[,x])),
               prop.table(table(d[,x], d$group),2))
  p <- chisq.test(table(d[,x], d$group))$p.value
  sum <- cbind(x,'',n[,1],pct[,1],'',n[,2],pct[,2],'',n[,3],pct[,3],'',p)
  sum
}

## compare prevalent CVD by age group
cvdp_age <- data.frame(t(sapply(cvdvars, function(x) 
  prop.table(table(factor(a1[,x],levels=0:1), a1$agegrp),2)[2,])))
write_csv(cvdp_age, 'cvd_prev_agegrp.csv')

## compare incidence between all cases and controls
table2_prev_cvd <- do.call(rbind,lapply(cvdvars ,function(x) tab2(a1,x)))
table2_prev_cvd <- data.frame(cbind(row.names(table2_prev_cvd), table2_prev_cvd))

## compare incidence between cases receiving chemo and controls
table2_prev_cvd_chemo <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_chemo,],x)))
table2_prev_cvd_chemo <- data.frame(cbind(row.names(table2_prev_cvd_chemo), table2_prev_cvd_chemo))

## compare incidence between cases receiving horm and controls
table2_prev_cvd_horm <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_horm,],x)))
table2_prev_cvd_horm <- data.frame(cbind(row.names(table2_prev_cvd_horm), table2_prev_cvd_horm))

## compare incidence between cases receiving rad and controls
table2_prev_cvd_rad <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_rad,],x)))
table2_prev_cvd_rad <- data.frame(cbind(row.names(table2_prev_cvd_rad), table2_prev_cvd_rad))

# export table
write_csv(table2_prev_cvd[table2_prev_cvd$V1=='Yes',], 'table2_prev_cvd.csv')
write_csv(table2_prev_cvd_chemo[table2_prev_cvd_chemo$V1=='Yes',], 'table2_prev_cvd_chemo.csv')
write_csv(table2_prev_cvd_horm[table2_prev_cvd_horm$V1=='Yes',], 'table2_prev_cvd_horm.csv')
write_csv(table2_prev_cvd_rad[table2_prev_cvd_rad$V1=='Yes',], 'table2_prev_cvd_rad.csv')









################################################################################
# Compare cvd frequencies between cases and controls

## True incidence
cvdvars <- grep('_inc$',names(a1),value = T)
cvdvars <- cvdvars[c(1:3,5,4,6:11,c(16,17,18,21,19,13:15,20,24,26,27,28,23,12,22,25,29:33))]


## compare incidence between all cases and controls
table2_cvd <- do.call(rbind,lapply(cvdvars ,function(x) tab2(a1,x)))
table2_cvd <- data.frame(cbind(row.names(table2_cvd), table2_cvd))

## compare incidence between cases receiving chemo and controls
table2_cvd_chemo <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_chemo,],x)))
table2_cvd_chemo <- data.frame(cbind(row.names(table2_cvd_chemo), table2_cvd_chemo))

## compare incidence between cases receiving horm and controls
table2_cvd_horm <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_horm,],x)))
table2_cvd_horm <- data.frame(cbind(row.names(table2_cvd_horm), table2_cvd_horm))

## compare incidence between cases receiving rad and controls
table2_cvd_rad <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_rad,],x)))
table2_cvd_rad <- data.frame(cbind(row.names(table2_cvd_rad), table2_cvd_rad))


write_csv(table2_cvd[table2_cvd$V1=='Yes',], 'table2_cvd.csv')
write_csv(table2_cvd_chemo[table2_cvd_chemo$V1=='Yes',], 'table2_cvd_chemo.csv')
write_csv(table2_cvd_horm[table2_cvd_horm$V1=='Yes',], 'table2_cvd_horm.csv')
write_csv(table2_cvd_rad[table2_cvd_rad$V1=='Yes',], 'table2_cvd_rad.csv')



## Any new onset = true incidence + recurrence
cvdvars <- grep('_rec$',names(a1),value = T)
cvdvars <- cvdvars[c(1:3,5,4,6:11,c(16,17,18,21,19,13:15,20,24,26,27,28,23,12,22,25,29:33))]

## compare incidence between all cases and controls
table2_rec_cvd <- do.call(rbind,lapply(cvdvars ,function(x) tab2(a1,x)))
table2_rec_cvd <- data.frame(cbind(row.names(table2_rec_cvd), table2_rec_cvd))

## compare incidence between cases receiving chemo and controls
table2_rec_cvd_chemo <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_chemo,],x)))
table2_rec_cvd_chemo <- data.frame(cbind(row.names(table2_rec_cvd_chemo), table2_rec_cvd_chemo))

## compare incidence between cases receiving horm and controls
table2_rec_cvd_horm <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_horm,],x)))
table2_rec_cvd_horm <- data.frame(cbind(row.names(table2_rec_cvd_horm), table2_rec_cvd_horm))

## compare incidence between cases receiving rad and controls
table2_rec_cvd_rad <- do.call(rbind,lapply(cvdvars,function(x) tab2(a1[a1_rad,],x)))
table2_rec_cvd_rad <- data.frame(cbind(row.names(table2_rec_cvd_rad), table2_rec_cvd_rad))


write_csv(table2_rec_cvd[table2_rec_cvd$V1=='Yes',], 'table2_rec_cvd.csv')
write_csv(table2_rec_cvd_chemo[table2_rec_cvd_chemo$V1=='Yes',], 'table2_rec_cvd_chemo.csv')
write_csv(table2_rec_cvd_horm[table2_rec_cvd_horm$V1=='Yes',], 'table2_rec_cvd_horm.csv')
write_csv(table2_rec_cvd_rad[table2_rec_cvd_rad$V1=='Yes',], 'table2_rec_cvd_rad.csv')








################################################################################
# Compare cvd risk factor frequencies between cases and controls

## cvd rf variables
cvdrfvars <- c("cvdrf_diab","cvdrf_htn","cvdrf_dyslipid","cvdrfcombo")

## compare incidence between all cases and controls
table2_cvdrf <- do.call(rbind,lapply(cvdrfvars ,function(x) tab2(a1,x)))
table2_cvdrf <- data.frame(cbind(row.names(table2_cvdrf), table2_cvdrf))

## compare incidence between cases receiving chemo and controls
table2_cvdrf_chemo <- do.call(rbind,lapply(cvdrfvars,function(x) tab2(a1[a1_chemo,],x)))
table2_cvdrf_chemo <- data.frame(cbind(row.names(table2_cvdrf_chemo), table2_cvdrf_chemo))

## compare incidence between cases receiving horm and controls
table2_cvdrf_horm <- do.call(rbind,lapply(cvdrfvars,function(x) tab2(a1[a1_horm,],x)))
table2_cvdrf_horm <- data.frame(cbind(row.names(table2_cvdrf_horm), table2_cvdrf_horm))

## compare incidence between cases receiving rad and controls
table2_cvdrf_rad <- do.call(rbind,lapply(cvdrfvars,function(x) tab2(a1[a1_rad,],x)))
table2_cvdrf_rad <- data.frame(cbind(row.names(table2_cvdrf_rad), table2_cvdrf_rad))


write_csv(table2_cvdrf[table2_cvdrf$V1=='Yes',], 'table2_cvdrf.csv')
write_csv(table2_cvdrf_chemo[table2_cvdrf_chemo$V1=='Yes',], 'table2_cvdrf_chemo.csv')
write_csv(table2_cvdrf_horm[table2_cvdrf_horm$V1=='Yes',], 'table2_cvdrf_horm.csv')
write_csv(table2_cvdrf_rad[table2_cvdrf_rad$V1=='Yes',], 'table2_cvdrf_rad.csv')








################################################################################
# Aim 1a - log rank test at 5 and 10 year follow up

# for all cases v controls

# true incidence

# function to calculate time specific incidence rate
inc_rate <- function(index, outcomevar, timevar){
  surv_dat <- a1[index,c(outcomevar,timevar,'group')]
  surv_dat$time <- round(surv_dat[,2]/30)
  surv_object <- Surv(time = surv_dat$time, 
                      event = surv_dat[,1])
  fit1 <- survfit(surv_object ~ group, data = surv_dat)
  logrank.24 <- survdiff(Surv(time = surv_dat$time[surv_dat$time <= 24], 
                              event = surv_dat[surv_dat$time <= 24,1]) ~ group, 
                         data = surv_dat[surv_dat$time <= 24,])
  logrank.60 <- survdiff(Surv(time = surv_dat$time[surv_dat$time <= 60], 
                              event = surv_dat[surv_dat$time <= 60,1]) ~ group, 
                         data = surv_dat[surv_dat$time <= 60,])
  logrank.120 <- survdiff(Surv(time = surv_dat$time[surv_dat$time <= 120], 
                              event = surv_dat[surv_dat$time <= 120,1]) ~ group, 
                         data = surv_dat[surv_dat$time <= 120,])
  pvalue.24 <- round(1 - pchisq(logrank.24$chisq, length(logrank.24$n) - 1),2)
  pvalue.60 <- round(1 - pchisq(logrank.60$chisq, length(logrank.60$n) - 1),2)
  pvalue.120 <- round(1 - pchisq(logrank.120$chisq, length(logrank.120$n) - 1),2)
  pvalue = c(pvalue.24, pvalue.60, pvalue.120)
  sum1 <- summary(fit1, times=c(24,60,120))
  sum2 <- 100*(1-cbind(sum1$surv,sum1$upper,sum1$lower))
  sum3 <- cbind(round2(sum2[,1],1),
                paste0("(",round2(sum2[,2],1),", ",round2(sum2[,3],1),")"))
  sum4 <- cbind(outcomevar,sum3[1:3,],pvalue,sum3[4:6,],round2(sum2[1:3,2],1),round2(sum2[1:3,3],1),round2(sum2[4:6,2],1),round2(sum2[4:6,3],1))
  sum4
}


inc_all <- data.frame(rbind(inc_rate(a1_ihd,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu'),
                            '',inc_rate(a1_stroke,'stroke_grp_inc','stroke_grp_inc_fu'),
                            '',inc_rate(a1_hf,'heart_failure_grp_inc','heart_failure_grp_inc_fu'),
                            '',inc_rate(a1_cm,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu'),
                            '',inc_rate(a1_combo,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu'),
                            '',inc_rate(a1_diab,'cvdrf_diab','cvdrf_diab_fu'),
                            '',inc_rate(a1_htn,'cvdrf_htn','cvdrf_htn_fu'),
                            '',inc_rate(a1_dyslipid,'cvdrf_dyslipid','cvdrf_dyslipid_fu'),
                            '',inc_rate(a1_rfcombo,'cvdrfcombo','cvdrfcombo_fu')))


inc_chemo <- data.frame(rbind(inc_rate(a1_ihd_chemo,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu'),
                              '',inc_rate(a1_stroke_chemo,'stroke_grp_inc','stroke_grp_inc_fu'),
                              '',inc_rate(a1_hf_chemo,'heart_failure_grp_inc','heart_failure_grp_inc_fu'),
                              '',inc_rate(a1_cm_chemo,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu'),
                              '',inc_rate(a1_combo_chemo,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu'),
                              '',inc_rate(a1_diab_chemo,'cvdrf_diab','cvdrf_diab_fu'),
                              '',inc_rate(a1_htn_chemo,'cvdrf_htn','cvdrf_htn_fu'),
                              '',inc_rate(a1_dyslipid_chemo,'cvdrf_dyslipid','cvdrf_dyslipid_fu'),
                              '',inc_rate(a1_rfcombo_chemo,'cvdrfcombo','cvdrfcombo_fu')))


inc_horm <- data.frame(rbind(inc_rate(a1_ihd_horm,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu'),
                              '',inc_rate(a1_stroke_horm,'stroke_grp_inc','stroke_grp_inc_fu'),
                              '',inc_rate(a1_hf_horm,'heart_failure_grp_inc','heart_failure_grp_inc_fu'),
                              '',inc_rate(a1_cm_horm,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu'),
                              '',inc_rate(a1_combo_horm,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu'),
                              '',inc_rate(a1_diab_horm,'cvdrf_diab','cvdrf_diab_fu'),
                              '',inc_rate(a1_htn_horm,'cvdrf_htn','cvdrf_htn_fu'),
                              '',inc_rate(a1_dyslipid_horm,'cvdrf_dyslipid','cvdrf_dyslipid_fu'),
                              '',inc_rate(a1_rfcombo_horm,'cvdrfcombo','cvdrfcombo_fu')))


inc_rad <- data.frame(rbind(inc_rate(a1_ihd_rad,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu'),
                              '',inc_rate(a1_stroke_rad,'stroke_grp_inc','stroke_grp_inc_fu'),
                              '',inc_rate(a1_hf_rad,'heart_failure_grp_inc','heart_failure_grp_inc_fu'),
                              '',inc_rate(a1_cm_rad,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu'),
                              '',inc_rate(a1_combo_rad,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu'),
                              '',inc_rate(a1_diab_rad,'cvdrf_diab','cvdrf_diab_fu'),
                              '',inc_rate(a1_htn_rad,'cvdrf_htn','cvdrf_htn_fu'),
                              '',inc_rate(a1_dyslipid_rad,'cvdrf_dyslipid','cvdrf_dyslipid_fu'),
                              '',inc_rate(a1_rfcombo_rad,'cvdrfcombo','cvdrfcombo_fu')))

inc_rad_l <- data.frame(rbind(inc_rate(a1_ihd_rad_l,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu'),
                            '',inc_rate(a1_stroke_rad_l,'stroke_grp_inc','stroke_grp_inc_fu'),
                            '',inc_rate(a1_hf_rad_l,'heart_failure_grp_inc','heart_failure_grp_inc_fu'),
                            '',inc_rate(a1_cm_rad_l,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu'),
                            '',inc_rate(a1_combo_rad_l,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu'),
                            '',inc_rate(a1_diab_rad_l,'cvdrf_diab','cvdrf_diab_fu'),
                            '',inc_rate(a1_htn_rad_l,'cvdrf_htn','cvdrf_htn_fu'),
                            '',inc_rate(a1_dyslipid_rad_l,'cvdrf_dyslipid','cvdrf_dyslipid_fu'),
                            '',inc_rate(a1_rfcombo_rad_l,'cvdrfcombo','cvdrfcombo_fu')))

inc_combined <- cbind(inc_all, inc_chemo, inc_horm, inc_rad, inc_rad_l)

write_csv(inc_combined, 'incidence_all.csv')


# Generate plots
# all
cvd = c(rep("Ischemic Heart Disease",6),rep("Stroke",6),rep("Heart Failure",6),rep("Cardiomyopathy",6),rep("Any",6))
cvd = factor(cvd, levels = c("Ischemic Heart Disease","Stroke","Heart Failure","Cardiomyopathy","Any"))
cvd.dat.all = data.frame(cvd = cvd)
cvd.dat.all$group = rep(c(rep("Case",3),rep("Control",3)),5)
inc_all2 = inc_all[!apply(inc_all == "", 1, all),]

cvd.dat.all$inc = NA
cvd.dat.all$inc_l = NA
cvd.dat.all$inc_u = NA
cvd.dat.all$p = NA

cvd.dat.all$inc[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_all2$V2[1:15]))
cvd.dat.all$inc_l[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_all2$V7[1:15]))
cvd.dat.all$inc_u[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_all2$V8[1:15]))
cvd.dat.all$p[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_all2$pvalue[1:15]))

cvd.dat.all$inc[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_all2$V5[1:15]))
cvd.dat.all$inc_l[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_all2$V9[1:15]))
cvd.dat.all$inc_u[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_all2$V10[1:15]))
cvd.dat.all$p[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_all2$pvalue[1:15]))

cvd.dat.all$year = rep(c("2 Years","5 Years","10 Years"),10)
cvd.dat.all$year = factor(cvd.dat.all$year, levels = c("2 Years","5 Years","10 Years"))


cvd.dat.all$inc2 = NA
for (i in 1:dim(cvd.dat.all)[1]){
  if (cvd.dat.all$p[i] < 0.05){
    cvd.dat.all$inc2[i] = paste0(format(round(cvd.dat.all$inc[i],1),nsmall=1), "*")
  }
  else{
    cvd.dat.all$inc2[i] = format(round(cvd.dat.all$inc[i],1),nsmall=1)
  }
  
}


png("Figure_inc_all.png", width = 14, height=8, res=300, units = 'in')
p1 <- ggplot(cvd.dat.all, aes(x = group,
                         y = inc, 
                         fill=year, 
                         color=year)) +
  geom_bar(stat="identity", color="white", 
           position=position_dodge()) + facet_wrap(~cvd, nrow = 1) + 
  geom_errorbar(aes(ymin  =inc_l, ymax = inc_u), width = .2, color = "gold",position=position_dodge(.9)) +
  geom_text(aes(group, inc+0.2, label = inc2, colour = NULL), position=position_dodge(width=0.9), data = cvd.dat.all, size = 4) +
  labs(y = "Cumulative Incidence Rate (%)", 
       x = "") + theme(legend.position = "bottom",legend.title=element_blank(),
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            axis.line.x = element_line(colour = "black"),
                            axis.line.y = element_line(colour = "black"),
                            panel.background = element_blank(),
                       axis.text.x = element_text(size = 12, face='bold'),
                       axis.text.y = element_text(size = 12, face='bold'),
                       axis.title.y = element_text(size = 12, face='bold'),
                       title = element_text(size=16, face='bold'),
                       legend.text = element_text(size = 12, face='bold'),
                       legend.margin=margin(-20, 0, 0, 0)) + ggtitle("A")
p1
dev.off()

# chemo
cvd = c(rep("Ischemic Heart Disease",6),rep("Stroke",6),rep("Heart Failure",6),rep("Cardiomyopathy",6),rep("Any",6))
cvd = factor(cvd, levels = c("Ischemic Heart Disease","Stroke","Heart Failure","Cardiomyopathy","Any"))
cvd.dat.chemo = data.frame(cvd = cvd)
cvd.dat.chemo$group = rep(c(rep("Case",3),rep("Control",3)),5)
inc_chemo2 = inc_chemo[!apply(inc_chemo == "", 1, all),]

cvd.dat.chemo$inc = NA
cvd.dat.chemo$inc_l = NA
cvd.dat.chemo$inc_u = NA
cvd.dat.chemo$p = NA

cvd.dat.chemo$inc[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_chemo2$V2[1:15]))
cvd.dat.chemo$inc_l[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_chemo2$V7[1:15]))
cvd.dat.chemo$inc_u[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_chemo2$V8[1:15]))
cvd.dat.chemo$p[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_chemo2$pvalue[1:15]))

cvd.dat.chemo$inc[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_chemo2$V5[1:15]))
cvd.dat.chemo$inc_l[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_chemo2$V9[1:15]))
cvd.dat.chemo$inc_u[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_chemo2$V10[1:15]))
cvd.dat.chemo$p[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_chemo2$pvalue[1:15]))

cvd.dat.chemo$year = rep(c("2 Years","5 Years","10 Years"),10)
cvd.dat.chemo$year = factor(cvd.dat.chemo$year, levels = c("2 Years","5 Years","10 Years"))

cvd.dat.chemo$inc2 = NA
for (i in 1:dim(cvd.dat.chemo)[1]){
  if (cvd.dat.chemo$p[i] < 0.05){
    cvd.dat.chemo$inc2[i] = paste0(format(round(cvd.dat.chemo$inc[i],1),nsmall=1), "*")
  }
  else{
    cvd.dat.chemo$inc2[i] = format(round(cvd.dat.chemo$inc[i],1),nsmall=1)
  }
  
}


png("Figure_inc_chemo.png", width = 14, height=8, res=300, units = 'in')
p2 <- ggplot(cvd.dat.chemo, aes(x = group,
                        y = inc, 
                        fill=year, 
                        color=year)) +
  geom_bar(stat="identity", color="white", 
           position=position_dodge()) + facet_wrap(~cvd, nrow = 1) + 
  geom_errorbar(aes(ymin  =inc_l, ymax = inc_u), width = .2, color = "gold",position=position_dodge(.9)) +
  geom_text(aes(group, inc+0.2, label = inc2, colour = NULL), position=position_dodge(width=0.9), data = cvd.dat.chemo, size = 4) +
  labs(y = "Cumulative Incidence Rate (%)", 
       x = "") + theme(legend.position = "bottom",legend.title=element_blank(),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.line.x = element_line(colour = "black"),
                       axis.line.y = element_line(colour = "black"),
                       panel.background = element_blank(),
                       axis.text.x = element_text(size = 12, face='bold'),
                       axis.text.y = element_text(size = 12, face='bold'),
                       axis.title.y = element_text(size = 12, face='bold'),
                       title = element_text(size=16, face='bold'),
                       legend.text = element_text(size = 12, face='bold'),
                       legend.margin=margin(-20, 0, 0, 0)) + ggtitle("B")
p2
dev.off()

# left-sided rad
cvd = c(rep("Ischemic Heart Disease",6),rep("Stroke",6),rep("Heart Failure",6),rep("Cardiomyopathy",6),rep("Any",6))
cvd = factor(cvd, levels = c("Ischemic Heart Disease","Stroke","Heart Failure","Cardiomyopathy","Any"))
cvd.dat.rad_l = data.frame(cvd = cvd)
cvd.dat.rad_l$group = rep(c(rep("Case",3),rep("Control",3)),5)
inc_rad_l2 = inc_rad_l[!apply(inc_rad_l == "", 1, all),]

cvd.dat.rad_l$inc = NA
cvd.dat.rad_l$inc_l = NA
cvd.dat.rad_l$inc_u = NA
cvd.dat.rad_l$p = NA

cvd.dat.rad_l$inc[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_rad_l2$V2[1:15]))
cvd.dat.rad_l$inc_l[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_rad_l2$V7[1:15]))
cvd.dat.rad_l$inc_u[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_rad_l2$V8[1:15]))
cvd.dat.rad_l$p[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_rad_l2$pvalue[1:15]))

cvd.dat.rad_l$inc[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_rad_l2$V5[1:15]))
cvd.dat.rad_l$inc_l[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_rad_l2$V9[1:15]))
cvd.dat.rad_l$inc_u[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_rad_l2$V10[1:15]))
cvd.dat.rad_l$p[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_rad_l2$pvalue[1:15]))


cvd.dat.rad_l$year = rep(c("2 Years","5 Years","10 Years"),10)
cvd.dat.rad_l$year = factor(cvd.dat.rad_l$year, levels = c("2 Years","5 Years","10 Years"))

cvd.dat.rad_l$inc2 = NA
for (i in 1:dim(cvd.dat.rad_l)[1]){
  if (cvd.dat.rad_l$p[i] < 0.05){
    cvd.dat.rad_l$inc2[i] = paste0(format(round(cvd.dat.rad_l$inc[i],1),nsmall=1), "*")
  }
  else{
    cvd.dat.rad_l$inc2[i] = format(round(cvd.dat.rad_l$inc[i],1),nsmall=1)
  }
  
}


png("Figure_inc_rad_l.png", width = 14, height=8, res=300, units = 'in')
p3 <- ggplot(cvd.dat.rad_l, aes(x = group,
                          y = inc, 
                          fill=year, 
                          color=year)) +
  geom_bar(stat="identity", color="white", 
           position=position_dodge()) + facet_wrap(~cvd, nrow = 1) + 
  geom_errorbar(aes(ymin  =inc_l, ymax = inc_u), width = .2, color = "gold",position=position_dodge(.9)) +
  geom_text(aes(group, inc+0.15, label = inc2, colour = NULL), position=position_dodge(width=0.9), data = cvd.dat.rad_l, size = 4) +
  labs(y = "Cumulative Incidence Rate (%)", 
       x = "") + theme(legend.position = "bottom",legend.title=element_blank(),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.line.x = element_line(colour = "black"),
                       axis.line.y = element_line(colour = "black"),
                       panel.background = element_blank(),
                       axis.text.x = element_text(size = 12, face='bold'),
                       axis.text.y = element_text(size = 12, face='bold'),
                       axis.title.y = element_text(size = 12, face='bold'),
                       title = element_text(size=16, face='bold'),
                       legend.text = element_text(size = 12, face='bold'),
                       legend.margin=margin(-20, 0, 0, 0)) + ggtitle("C")
p3
dev.off()

# horm
cvd = c(rep("Ischemic Heart Disease",6),rep("Stroke",6),rep("Heart Failure",6),rep("Cardiomyopathy",6),rep("Any",6))
cvd = factor(cvd, levels = c("Ischemic Heart Disease","Stroke","Heart Failure","Cardiomyopathy","Any"))
cvd.dat.horm = data.frame(cvd = cvd)
cvd.dat.horm$group = rep(c(rep("Case",3),rep("Control",3)),5)
inc_horm2 = inc_horm[!apply(inc_horm == "", 1, all),]

cvd.dat.horm$inc = NA
cvd.dat.horm$inc_l = NA
cvd.dat.horm$inc_u = NA
cvd.dat.horm$p = NA

cvd.dat.horm$inc[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_horm2$V2[1:15]))
cvd.dat.horm$inc_l[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_horm2$V7[1:15]))
cvd.dat.horm$inc_u[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_horm2$V8[1:15]))
cvd.dat.horm$p[c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_horm2$pvalue[1:15]))

cvd.dat.horm$inc[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_horm2$V5[1:15]))
cvd.dat.horm$inc_l[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_horm2$V9[1:15]))
cvd.dat.horm$inc_u[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_horm2$V10[1:15]))
cvd.dat.horm$p[-c(1:3,7:9,13:15,19:21,25:27)] = as.numeric(as.character(inc_horm2$pvalue[1:15]))


cvd.dat.horm$year = rep(c("2 Years","5 Years","10 Years"),10)
cvd.dat.horm$year = factor(cvd.dat.horm$year, levels = c("2 Years","5 Years","10 Years"))

cvd.dat.horm$inc2 = NA
for (i in 1:dim(cvd.dat.horm)[1]){
  if (cvd.dat.horm$p[i] < 0.05){
    cvd.dat.horm$inc2[i] = paste0(format(round(cvd.dat.horm$inc[i],1),nsmall=1), "*")
  }
  else{
    cvd.dat.horm$inc2[i] = format(round(cvd.dat.horm$inc[i],1),nsmall=1)
  }
  
}

png("Figure_inc_horm.png", width = 14, height=8, res=300, units = 'in')
p4 <- ggplot(cvd.dat.horm, aes(x = group,
                        y = inc, 
                        fill=year, 
                        color=year)) +
  geom_bar(stat="identity", color="white", 
           position=position_dodge()) + facet_wrap(~cvd, nrow = 1) + 
  geom_errorbar(aes(ymin  =inc_l, ymax = inc_u), width = .2, color = "gold",position=position_dodge(.9)) +
  geom_text(aes(group, inc+0.15, label = inc2, colour = NULL), position=position_dodge(width=0.9), data = cvd.dat.horm, size = 4) +
  labs(y = "Cumulative Incidence Rate (%)", 
       x = "") + theme(legend.position = "bottom",legend.title=element_blank(),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.line.x = element_line(colour = "black"),
                       axis.line.y = element_line(colour = "black"),
                       panel.background = element_blank(),
                       axis.text.x = element_text(size = 12, face='bold'),
                       axis.text.y = element_text(size = 12, face='bold'),
                       axis.title.y = element_text(size = 12, face='bold'),
                       title = element_text(size=16, face='bold'),
                       legend.text = element_text(size = 12, face='bold'),
                       legend.margin=margin(-20, 0, 0, 0)) + ggtitle("D")
p4
dev.off()


# Figure Panel 
png("Figure_1.png", width = 20, height=10, res=300, units = 'in')
ggarrange(p1, p2, p3, p4,
           ncol = 2, nrow = 2)
dev.off()

################################################################################
# Aim 1a - K-M analysis figures


# for all cases v controls

# true incidence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_ihd,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1[a1_ihd,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_ihd,],fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(a1[a1_stroke,]$stroke_grp_inc_fu/30), 
                    event = a1[a1_stroke,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_stroke,],fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(a1[a1_hf,]$heart_failure_grp_inc_fu/30), 
                    event = a1[a1_hf,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_hf,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_hf,],fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(a1[a1_cm,]$cardiomyopathy_grp_inc_fu/30), 
                    event = a1[a1_cm,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_cm,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_cm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(a1[a1_combo,]$cvdcombo_grp_inc_fu/30), 
                    event = a1[a1_combo,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_combo,])
splots[[5]] <- ggsurvplot(fit1, data = a1[a1_combo,],fun = "event", size=1, risk.table = T,pval = T,
           censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig1 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig1.png", fig1,width = 13, height = 18)


# any new onset=true incidence+recurrence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_ihd,]$ischemic_heart_disease_grp_rec_fu/30), 
                    event = a1[a1_ihd,]$ischemic_heart_disease_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_ihd,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(a1[a1_stroke,]$stroke_grp_rec_fu/30), 
                    event = a1[a1_stroke,]$stroke_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_stroke,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(a1[a1_hf,]$heart_failure_grp_rec_fu/30), 
                    event = a1[a1_hf,]$heart_failure_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_hf,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_hf,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(a1[a1_cm,]$cardiomyopathy_grp_rec_fu/30), 
                    event = a1[a1_cm,]$cardiomyopathy_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_cm,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_cm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(a1[a1_combo,]$cvdcombo_grp_rec_fu/30), 
                    event = a1[a1_combo,]$cvdcombo_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_combo,])
splots[[5]] <- ggsurvplot(fit1, data = a1[a1_combo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig2 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig2.png", fig2,width = 13, height = 18)








################################################################################
# Aim 1b - K-M analysis figures

###
# for chemo cases v controls

# true incidence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_ihd_chemo,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1[a1_ihd_chemo,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd_chemo,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_ihd_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(a1[a1_stroke_chemo,]$stroke_grp_inc_fu/30), 
                    event = a1[a1_stroke_chemo,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke_chemo,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_stroke_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(a1[a1_hf_chemo,]$heart_failure_grp_inc_fu/30), 
                    event = a1[a1_hf_chemo,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_hf_chemo,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_hf_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(a1[a1_cm_chemo,]$cardiomyopathy_grp_inc_fu/30), 
                    event = a1[a1_cm_chemo,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_cm_chemo,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_cm_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(a1[a1_combo_chemo,]$cvdcombo_grp_inc_fu/30), 
                    event = a1[a1_combo_chemo,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_combo_chemo,])
splots[[5]] <- ggsurvplot(fit1, data = a1[a1_combo_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig3 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig3.png", fig3,width = 13, height = 18)


# any new onset=true incidence+recurrence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_ihd_chemo,]$ischemic_heart_disease_grp_rec_fu/30), 
                    event = a1[a1_ihd_chemo,]$ischemic_heart_disease_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd_chemo,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_ihd_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(a1[a1_stroke_chemo,]$stroke_grp_rec_fu/30), 
                    event = a1[a1_stroke_chemo,]$stroke_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke_chemo,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_stroke_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(a1[a1_hf_chemo,]$heart_failure_grp_rec_fu/30), 
                    event = a1[a1_hf_chemo,]$heart_failure_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_hf_chemo,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_hf_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(a1[a1_cm_chemo,]$cardiomyopathy_grp_rec_fu/30), 
                    event = a1[a1_cm_chemo,]$cardiomyopathy_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_cm_chemo,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_cm_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(a1[a1_combo_chemo,]$cvdcombo_grp_rec_fu/30), 
                    event = a1[a1_combo_chemo,]$cvdcombo_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_combo_chemo,])
splots[[5]] <- ggsurvplot(fit1, data = a1[a1_combo_chemo,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig4 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig4.png", fig4,width = 13, height = 18)





###
# for hormonal therapy cases v controls

# true incidence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_ihd_horm,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1[a1_ihd_horm,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd_horm,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_ihd_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(a1[a1_stroke_horm,]$stroke_grp_inc_fu/30), 
                    event = a1[a1_stroke_horm,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke_horm,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_stroke_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(a1[a1_hf_horm,]$heart_failure_grp_inc_fu/30), 
                    event = a1[a1_hf_horm,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_hf_horm,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_hf_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(a1[a1_cm_horm,]$cardiomyopathy_grp_inc_fu/30), 
                    event = a1[a1_cm_horm,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_cm_horm,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_cm_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(a1[a1_combo_horm,]$cvdcombo_grp_inc_fu/30), 
                    event = a1[a1_combo_horm,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_combo_horm,])
splots[[5]] <- ggsurvplot(fit1, data = a1[a1_combo_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig5 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig5.png", fig5, width = 13, height = 18)


# any new onset=true incidence+recurrence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_ihd_horm,]$ischemic_heart_disease_grp_rec_fu/30), 
                    event = a1[a1_ihd_horm,]$ischemic_heart_disease_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd_horm,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_ihd_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(a1[a1_stroke_horm,]$stroke_grp_rec_fu/30), 
                    event = a1[a1_stroke_horm,]$stroke_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke_horm,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_stroke_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(a1[a1_hf_horm,]$heart_failure_grp_rec_fu/30), 
                    event = a1[a1_hf_horm,]$heart_failure_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_hf_horm,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_hf_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(a1[a1_cm_horm,]$cardiomyopathy_grp_rec_fu/30), 
                    event = a1[a1_cm_horm,]$cardiomyopathy_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_cm_horm,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_cm_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(a1[a1_combo_horm,]$cvdcombo_grp_rec_fu/30), 
                    event = a1[a1_combo_horm,]$cvdcombo_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_combo_horm,])
splots[[5]] <- ggsurvplot(fit1, data = a1[a1_combo_horm,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig6 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig6.png", fig6, width = 13, height = 18)





###
# for radiation cases v controls

# true incidence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_ihd_rad,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1[a1_ihd_rad,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd_rad,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_ihd_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(a1[a1_stroke_rad,]$stroke_grp_inc_fu/30), 
                    event = a1[a1_stroke_rad,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke_rad,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_stroke_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(a1[a1_hf_rad,]$heart_failure_grp_inc_fu/30), 
                    event = a1[a1_hf_rad,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_hf_rad,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_hf_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(a1[a1_cm_rad,]$cardiomyopathy_grp_inc_fu/30), 
                    event = a1[a1_cm_rad,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_cm_rad,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_cm_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(a1[a1_combo_rad,]$cvdcombo_grp_inc_fu/30), 
                    event = a1[a1_combo_rad,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_combo_rad,])
splots[[5]] <- ggsurvplot(fit1, data = a1[a1_combo_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig7 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig7.png", fig7, width = 13, height = 18)


# any new onset=true incidence+recurrence
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_ihd_rad,]$ischemic_heart_disease_grp_rec_fu/30), 
                    event = a1[a1_ihd_rad,]$ischemic_heart_disease_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd_rad,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_ihd_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(a1[a1_stroke_rad,]$stroke_grp_rec_fu/30), 
                    event = a1[a1_stroke_rad,]$stroke_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke_rad,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_stroke_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(a1[a1_hf_rad,]$heart_failure_grp_rec_fu/30), 
                    event = a1[a1_hf_rad,]$heart_failure_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_hf_rad,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_hf_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(a1[a1_cm_rad,]$cardiomyopathy_grp_rec_fu/30), 
                    event = a1[a1_cm_rad,]$cardiomyopathy_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_cm_rad,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_cm_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(a1[a1_combo_rad,]$cvdcombo_grp_rec_fu/30), 
                    event = a1[a1_combo_rad,]$cvdcombo_grp_rec)
fit1 <- survfit(surv_object ~ group, data = a1[a1_combo_rad,])
splots[[5]] <- ggsurvplot(fit1, data = a1[a1_combo_rad,],fun = "event", size=1, risk.table = T,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig8 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig8.png", fig8, width = 13, height = 18)




################################################################################
# Aim 1a - Cox model

# function to summary cox model results
coxtab <- function(d,x,covar){
  covars <- paste(covar, collapse = '+')
  surv_object <- Surv(time = floor(d[,paste0(x,'_fu')]/30), event = d[,x])
  fit.coxph1 <- coxph(surv_object ~ group1, data = d)
  fit.coxph2 <- eval(parse(text=(paste0("coxph(surv_object ~ group1+bmicat1+menop+smok+
                        diab_bl+htn_bl+dyslipid_bl+medhousincome_q+edu1_q+",covars,", data = d)"))))
  sum1 <- c(summary(fit.coxph1)$coef[c(2,5)],summary(fit.coxph1)$conf.int[3:4])
  sum2 <- cbind(summary(fit.coxph2)$coef[,c(2,5)],summary(fit.coxph2)$conf.int[,3:4])
  sum3 <- rbind(sum1, sum2)[,c(1,3,4,2)]
  sum4 <- data.frame(cbind(var=row.names(sum3), sum3))
  sum4$var <- factor(sum4$var, levels=row.names(sum3))
  sum4
}



coxtab_htn <- function(d,x,covar){
  covars <- paste(covar, collapse = '+')
  surv_object <- Surv(time = floor(d[,paste0(x,'_fu')]/30), event = d[,x])
  fit.coxph1 <- coxph(surv_object ~ group1, data = d)
  fit.coxph2 <- eval(parse(text=(paste0("coxph(surv_object ~ group1+bmicat1+menop+smok+
                        diab_bl+dyslipid_bl+medhousincome_q+edu1_q+",covars,", data = d)"))))
  sum1 <- c(summary(fit.coxph1)$coef[c(2,5)],summary(fit.coxph1)$conf.int[3:4])
  sum2 <- cbind(summary(fit.coxph2)$coef[,c(2,5)],summary(fit.coxph2)$conf.int[,3:4])
  sum3 <- rbind(sum1, sum2)[,c(1,3,4,2)]
  sum4 <- data.frame(cbind(var=row.names(sum3), sum3))
  sum4$var <- factor(sum4$var, levels=row.names(sum3))
  sum4
}

coxtab_diab <- function(d,x,covar){
  covars <- paste(covar, collapse = '+')
  surv_object <- Surv(time = floor(d[,paste0(x,'_fu')]/30), event = d[,x])
  fit.coxph1 <- coxph(surv_object ~ group1, data = d)
  fit.coxph2 <- eval(parse(text=(paste0("coxph(surv_object ~ group1+bmicat1+menop+smok+
                        htn_bl+dyslipid_bl+medhousincome_q+edu1_q+",covars,", data = d)"))))
  sum1 <- c(summary(fit.coxph1)$coef[c(2,5)],summary(fit.coxph1)$conf.int[3:4])
  sum2 <- cbind(summary(fit.coxph2)$coef[,c(2,5)],summary(fit.coxph2)$conf.int[,3:4])
  sum3 <- rbind(sum1, sum2)[,c(1,3,4,2)]
  sum4 <- data.frame(cbind(var=row.names(sum3), sum3))
  sum4$var <- factor(sum4$var, levels=row.names(sum3))
  sum4
}

coxtab_dys <- function(d,x,covar){
  covars <- paste(covar, collapse = '+')
  surv_object <- Surv(time = floor(d[,paste0(x,'_fu')]/30), event = d[,x])
  fit.coxph1 <- coxph(surv_object ~ group1, data = d)
  fit.coxph2 <- eval(parse(text=(paste0("coxph(surv_object ~ group1+bmicat1+menop+smok+
                        diab_bl+htn_bl+medhousincome_q+edu1_q+",covars,", data = d)"))))
  sum1 <- c(summary(fit.coxph1)$coef[c(2,5)],summary(fit.coxph1)$conf.int[3:4])
  sum2 <- cbind(summary(fit.coxph2)$coef[,c(2,5)],summary(fit.coxph2)$conf.int[,3:4])
  sum3 <- rbind(sum1, sum2)[,c(1,3,4,2)]
  sum4 <- data.frame(cbind(var=row.names(sum3), sum3))
  sum4$var <- factor(sum4$var, levels=row.names(sum3))
  sum4
}

# list of prevalent cvd
prevcvd <- c("arrhythmia_grp_prev","cardiomyopathy_grp_prev","heart_failure_grp_prev",
            "ischemic_heart_disease_grp_prev","myocarditis_pericarditis_grp_prev",    
            "stroke_grp_prev","tia_grp_prev","valvular_disease_grp_prev","venous_thromboembolic_disease_grp_prev")

###
# Pure incidence
cox_isch <- coxtab(a1[a1_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[-4])
cox_stroke <- coxtab(a1[a1_stroke,],'stroke_grp_inc',prevcvd[-6])
cox_hf <- coxtab(a1[a1_hf,],'heart_failure_grp_inc',prevcvd[-3])
cox_cm <- coxtab(a1[a1_cm,],'cardiomyopathy_grp_inc',prevcvd[-2])
cox_hf_cm <- coxtab(a1[a1_hf_cm,],'hfcm_grp_inc',prevcvd[c(-2,-3)])
cox_cvdcombo <- coxtab(a1[a1_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)])
cox_arrhythmia <- coxtab(a1[a1_arrhythmia,],'arrhythmia_grp_inc',prevcvd[c(-1)])
cox_cardiac <- coxtab(a1[a1_cardiac,],'cardiac_arrest_grp_inc',prevcvd[c(-1)])
cox_carotid <- coxtab(a1[a1_carotid,],'carotid_disease_grp_inc',prevcvd[c(-1)])
cox_myocarditis <- coxtab(a1[a1_myocarditis,],'myocarditis_pericarditis_grp_inc',prevcvd[c(-1)])
cox_tia <- coxtab(a1[a1_tia,],'tia_grp_inc',prevcvd[c(-1)])
cox_valvular <- coxtab(a1[a1_valvular,],'valvular_disease_grp_inc',prevcvd[c(-1)])
cox_dvt <- coxtab(a1[a1_dvt,],'venous_thromboembolic_disease_grp_inc',prevcvd[c(-1)])
cox_diab <- coxtab_diab(a1[a1_diab,],'cvdrf_diab',prevcvd)
cox_htn <- coxtab_htn(a1[a1_htn,],'cvdrf_htn',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo,],'cvdrfcombo',prevcvd)
cox_allcvd <- coxtab(a1[a1_allcvd,],'cvdrfcombo',prevcvd[c(-1,-2,-3,-4,-6)])
cox1 <- rbind(cox_isch[1:2,],cox_hf[1:2,], cox_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,],"",
              cox_arrhythmia[1:2,], cox_cardiac[1:2,], cox_carotid[1:2,], 
              cox_myocarditis[1:2,], cox_tia[1:2,], cox_valvular[1:2,], cox_dvt[1:2,], cox_allcvd[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round2(cox1[,2],2),paste0("(",round2(cox1[,3],2),", ", round2(cox1[,4],2),")"))

# # Any new onset
# cox_isch <- coxtab(a1,'ischemic_heart_disease_grp_rec',prevcvd)
# cox_stroke <- coxtab(a1,'stroke_tia_grp_rec',prevcvd)
# cox_chf <- coxtab(a1,'cardiomyopathy_heart_failure_grp_rec',prevcvd)
# cox_cvdcombo <- coxtab(a1,'cvdcombo_grp_rec',prevcvd)
# cox1rec <- Reduce(function(...) merge(...,by='var', all=T),
#                list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
# cox1rec$var <- factor(cox1rec$var, levels=levels(cox_isch$var))
# cox1rec[,-1] <- lapply(cox1rec[,-1], as.character)
# cox1rec <- cox1rec[order(cox1rec$var),]
# cox1rec <- rbind(cox1rec[1,],'',cox1rec[2,],'',cox1rec[3:7,],cox1rec[8,],'',cox1rec[9:11,],cox1rec[12,],'',
#               cox1rec[13:14,],'',cox1rec[15:16,],'',cox1rec[17:20,],'',cox1rec[21:24,],'',
#               cox1rec[25:27,],'',cox1rec[28:34,])
# cox1rec <- cbind(cox1rec[,1:5],'',cox1rec[,6:9],'',cox1rec[,10:13],'',cox1rec[,14:17])
# 







################################################################################
# Aim 1b - Cox model

###
# chemo cases
# Pure incidence
cox_isch <- coxtab(a1[a1_ihd_chemo,],'ischemic_heart_disease_grp_inc',prevcvd[-4])
cox_stroke <- coxtab(a1[a1_stroke_chemo,],'stroke_grp_inc',prevcvd[-6])
cox_hf <- coxtab(a1[a1_hf_chemo,],'heart_failure_grp_inc',prevcvd[-3])
cox_cm <- coxtab(a1[a1_cm_chemo,],'cardiomyopathy_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_chemo,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)])
cox_arrhythmia <- coxtab(a1[a1_arrhythmia_chemo,],'arrhythmia_grp_inc',prevcvd[c(-1)])
cox_cardiac <- coxtab(a1[a1_cardiac_chemo,],'cardiac_arrest_grp_inc',prevcvd[c(-1)])
cox_carotid <- coxtab(a1[a1_carotid_chemo,],'carotid_disease_grp_inc',prevcvd[c(-1)])
cox_myocarditis <- coxtab(a1[a1_myocarditis_chemo,],'myocarditis_pericarditis_grp_inc',prevcvd[c(-1)])
cox_tia <- coxtab(a1[a1_tia_chemo,],'tia_grp_inc',prevcvd[c(-1)])
cox_valvular <- coxtab(a1[a1_valvular_chemo,],'valvular_disease_grp_inc',prevcvd[c(-1)])
cox_dvt <- coxtab(a1[a1_dvt_chemo,],'venous_thromboembolic_disease_grp_inc',prevcvd[c(-1)])
cox_diab <- coxtab_diab(a1[a1_diab_chemo,],'cvdrf_diab',prevcvd)
cox_htn <- coxtab_htn(a1[a1_htn_chemo,],'cvdrf_htn',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_chemo,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_chemo,],'cvdrfcombo',prevcvd)
cox2 <- rbind(cox_isch[1:2,],cox_hf[1:2,], cox_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,],"",
              cox_arrhythmia[1:2,], cox_cardiac[1:2,], cox_carotid[1:2,], 
              cox_myocarditis[1:2,], cox_tia[1:2,], cox_valvular[1:2,], cox_dvt[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round2(cox2[,2],2),paste0("(",round2(cox2[,3],2),", ", round2(cox2[,4],2),")"))

# # Any new onset
# cox_isch <- coxtab(a1[a1_chemo,],'ischemic_heart_disease_grp_rec',prevcvd)
# cox_stroke <- coxtab(a1[a1_chemo,],'stroke_tia_grp_rec',prevcvd)
# cox_chf <- coxtab(a1[a1_chemo,],'cardiomyopathy_heart_failure_grp_rec',prevcvd)
# cox_cvdcombo <- coxtab(a1[a1_chemo,],'cvdcombo_grp_rec',prevcvd)
# cox2rec <- Reduce(function(...) merge(...,by='var', all=T),
#                   list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
# cox2rec$var <- factor(cox2rec$var, levels=levels(cox_isch$var))
# cox2rec[,-1] <- lapply(cox2rec[,-1], as.character)
# cox2rec <- cox2rec[order(cox2rec$var),]
# cox2rec <- rbind(cox2rec[1,],'',cox2rec[2,],'',cox2rec[3:7,],cox2rec[8,],'',cox2rec[9:11,],cox2rec[12,],'',
#                  cox2rec[13:14,],'',cox2rec[15:16,],'',cox2rec[17:20,],'',cox2rec[21:24,],'',
#                  cox2rec[25:27,],'',cox2rec[28:34,])
# cox2rec <- cbind(cox2rec[,1:5],'',cox2rec[,6:9],'',cox2rec[,10:13],'',cox2rec[,14:17])



# horm
# Pure incidence
cox_isch <- coxtab(a1[a1_ihd_horm,],'ischemic_heart_disease_grp_inc',prevcvd[-4])
cox_stroke <- coxtab(a1[a1_stroke_horm,],'stroke_grp_inc',prevcvd[-6])
cox_hf <- coxtab(a1[a1_hf_horm,],'heart_failure_grp_inc',prevcvd[-3])
cox_cm <- coxtab(a1[a1_cm_horm,],'cardiomyopathy_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_horm,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)])
cox_arrhythmia <- coxtab(a1[a1_arrhythmia_horm,],'arrhythmia_grp_inc',prevcvd[c(-1)])
cox_cardiac <- coxtab(a1[a1_cardiac_horm,],'cardiac_arrest_grp_inc',prevcvd[c(-1)])
cox_carotid <- coxtab(a1[a1_carotid_horm,],'carotid_disease_grp_inc',prevcvd[c(-1)])
cox_myocarditis <- coxtab(a1[a1_myocarditis_horm,],'myocarditis_pericarditis_grp_inc',prevcvd[c(-1)])
cox_tia <- coxtab(a1[a1_tia_horm,],'tia_grp_inc',prevcvd[c(-1)])
cox_valvular <- coxtab(a1[a1_valvular_horm,],'valvular_disease_grp_inc',prevcvd[c(-1)])
cox_dvt <- coxtab(a1[a1_dvt_horm,],'venous_thromboembolic_disease_grp_inc',prevcvd[c(-1)])
cox_diab <- coxtab_diab(a1[a1_diab_horm,],'cvdrf_diab',prevcvd)
cox_htn <- coxtab_htn(a1[a1_htn_horm,],'cvdrf_htn',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_horm,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_horm,],'cvdrfcombo',prevcvd)
cox3 <- rbind(cox_isch[1:2,],cox_hf[1:2,], cox_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,],"",
              cox_arrhythmia[1:2,], cox_cardiac[1:2,], cox_carotid[1:2,], 
              cox_myocarditis[1:2,], cox_tia[1:2,], cox_valvular[1:2,], cox_dvt[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
cox3 <- cbind(round2(cox3[,2],2),paste0("(",round2(cox3[,3],2),", ", round2(cox3[,4],2),")"))

# # Any new onset
# cox_isch <- coxtab(a1[a1_horm,],'ischemic_heart_disease_grp_rec',prevcvd)
# cox_stroke <- coxtab(a1[a1_horm,],'stroke_tia_grp_rec',prevcvd)
# cox_chf <- coxtab(a1[a1_horm,],'cardiomyopathy_heart_failure_grp_rec',prevcvd)
# cox_cvdcombo <- coxtab(a1[a1_horm,],'cvdcombo_grp_rec',prevcvd)
# cox3rec <- Reduce(function(...) merge(...,by='var', all=T),
#                   list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
# cox3rec$var <- factor(cox3rec$var, levels=levels(cox_isch$var))
# cox3rec[,-1] <- lapply(cox3rec[,-1], as.character)
# cox3rec <- cox3rec[order(cox3rec$var),]
# cox3rec <- rbind(cox3rec[1,],'',cox3rec[2,],'',cox3rec[3:7,],cox3rec[8,],'',cox3rec[9:11,],cox3rec[12,],'',
#                  cox3rec[13:14,],'',cox3rec[15:16,],'',cox3rec[17:20,],'',cox3rec[21:24,],'',
#                  cox3rec[25:27,],'',cox3rec[28:34,])
# cox3rec <- cbind(cox3rec[,1:5],'',cox3rec[,6:9],'',cox3rec[,10:13],'',cox3rec[,14:17])


# rad
# Pure incidence
cox_isch <- coxtab(a1[a1_ihd_rad,],'ischemic_heart_disease_grp_inc',prevcvd[-4])
cox_stroke <- coxtab(a1[a1_stroke_rad,],'stroke_grp_inc',prevcvd[-6])
cox_hf <- coxtab(a1[a1_hf_rad,],'heart_failure_grp_inc',prevcvd[-3])
cox_cm <- coxtab(a1[a1_cm_rad,],'cardiomyopathy_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_rad,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)])
cox_arrhythmia <- coxtab(a1[a1_arrhythmia_rad,],'arrhythmia_grp_inc',prevcvd[c(-1)])
cox_cardiac <- coxtab(a1[a1_cardiac_rad,],'cardiac_arrest_grp_inc',prevcvd[c(-1)])
cox_carotid <- coxtab(a1[a1_carotid_rad,],'carotid_disease_grp_inc',prevcvd[c(-1)])
cox_myocarditis <- coxtab(a1[a1_myocarditis_rad,],'myocarditis_pericarditis_grp_inc',prevcvd[c(-1)])
cox_tia <- coxtab(a1[a1_tia_rad,],'tia_grp_inc',prevcvd[c(-1)])
cox_valvular <- coxtab(a1[a1_valvular_rad,],'valvular_disease_grp_inc',prevcvd[c(-1)])
cox_dvt <- coxtab(a1[a1_dvt_rad,],'venous_thromboembolic_disease_grp_inc',prevcvd[c(-1)])
cox_diab <- coxtab_diab(a1[a1_diab_rad,],'cvdrf_diab',prevcvd)
cox_htn <- coxtab_htn(a1[a1_htn_rad,],'cvdrf_htn',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_rad,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_rad,],'cvdrfcombo',prevcvd)
cox4 <- rbind(cox_isch[1:2,],cox_hf[1:2,], cox_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,],"",
              cox_arrhythmia[1:2,], cox_cardiac[1:2,], cox_carotid[1:2,], 
              cox_myocarditis[1:2,], cox_tia[1:2,], cox_valvular[1:2,], cox_dvt[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
cox4 <- cbind(round2(cox4[,2],2),paste0("(",round2(cox4[,3],2),", ", round2(cox4[,4],2),")"))



# left side BC, receive RT
# Pure incidence
cox_isch <- coxtab(a1[a1_ihd_rad_l,],'ischemic_heart_disease_grp_inc',prevcvd[-4])
cox_stroke <- coxtab(a1[a1_stroke_rad_l,],'stroke_grp_inc',prevcvd[-6])
cox_hf <- coxtab(a1[a1_hf_rad_l,],'heart_failure_grp_inc',prevcvd[-3])
cox_cm <- coxtab(a1[a1_cm_rad_l,],'cardiomyopathy_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_rad_l,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)])
cox_arrhythmia <- coxtab(a1[a1_arrhythmia_rad_l,],'arrhythmia_grp_inc',prevcvd[c(-1)])
cox_cardiac <- coxtab(a1[a1_cardiac_rad_l,],'cardiac_arrest_grp_inc',prevcvd[c(-1)])
cox_carotid <- coxtab(a1[a1_carotid_rad_l,],'carotid_disease_grp_inc',prevcvd[c(-1)])
cox_myocarditis <- coxtab(a1[a1_myocarditis_rad_l,],'myocarditis_pericarditis_grp_inc',prevcvd[c(-1)])
cox_tia <- coxtab(a1[a1_tia_rad_l,],'tia_grp_inc',prevcvd[c(-1)])
cox_valvular <- coxtab(a1[a1_valvular_rad_l,],'valvular_disease_grp_inc',prevcvd[c(-1)])
cox_dvt <- coxtab(a1[a1_dvt_rad_l,],'venous_thromboembolic_disease_grp_inc',prevcvd[c(-1)])
cox_diab <- coxtab_diab(a1[a1_diab_rad_l,],'cvdrf_diab',prevcvd)
cox_htn <- coxtab_htn(a1[a1_htn_rad_l,],'cvdrf_htn',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_rad_l,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_rad_l,],'cvdrfcombo',prevcvd)
cox5 <- rbind(cox_isch[1:2,],cox_hf[1:2,], cox_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,],"",
              cox_arrhythmia[1:2,], cox_cardiac[1:2,], cox_carotid[1:2,], 
              cox_myocarditis[1:2,], cox_tia[1:2,], cox_valvular[1:2,], cox_dvt[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
cox5 <- cbind(round2(cox5[,2],2),paste0("(",round2(cox5[,3],2),", ", round2(cox5[,4],2),")"))

# # Any new onset
# cox_isch <- coxtab(a1[a1_rad,],'ischemic_heart_disease_grp_rec',prevcvd)
# cox_stroke <- coxtab(a1[a1_rad,],'stroke_tia_grp_rec',prevcvd)
# cox_chf <- coxtab(a1[a1_rad,],'cardiomyopathy_heart_failure_grp_rec',prevcvd)
# cox_cvdcombo <- coxtab(a1[a1_rad,],'cvdcombo_grp_rec',prevcvd)
# cox4rec <- Reduce(function(...) merge(...,by='var', all=T),
#                   list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
# cox4rec$var <- factor(cox4rec$var, levels=levels(cox_isch$var))
# cox4rec[,-1] <- lapply(cox4rec[,-1], as.character)
# cox4rec <- cox4rec[order(cox4rec$var),]
# cox4rec <- rbind(cox4rec[1,],'',cox4rec[2,],'',cox4rec[3:7,],cox4rec[8,],'',cox4rec[9:11,],cox4rec[12,],'',
#                  cox4rec[13:14,],'',cox4rec[15:16,],'',cox4rec[17:20,],'',cox4rec[21:24,],'',
#                  cox4rec[25:27,],'',cox4rec[28:34,])
# cox4rec <- cbind(cox4rec[,1:5],'',cox4rec[,6:9],'',cox4rec[,10:13],'',cox4rec[,14:17])


# # left side tumor, rad
# # Pure incidence
# cox_isch <- coxtab(a1[a1_ihd_rad_l,],'ischemic_heart_disease_grp_inc',prevcvd[-3])
# cox_stroke <- coxtab(a1[a1_stroke_rad_l,],'stroke_tia_grp_inc',prevcvd[-5])
# cox_chf <- coxtab(a1[a1_chf_rad_l,],'cardiomyopathy_heart_failure_grp_inc',prevcvd[-2])
# cox_cvdcombo <- coxtab(a1[a1_combo_rad_l,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-5)])
# cox5 <- Reduce(function(...) merge(...,by='var', all=T),
#                list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
# cox5[,-1] <- lapply(cox5[,-1], as.character)
# cox5 <- rbind(cox5[1,],'',cox5[2,],'',cox5[3:7,],cox5[8,],'',cox5[9:11,],cox5[12,],'',
#               cox5[13:14,],'',cox5[15:16,],'',cox5[17:20,],'',cox5[21:24,],'',
#               cox5[25:27,],'',cox5[28:34,])
# cox5 <- cbind(cox5[,1:5],'',cox5[,6:9],'',cox5[,10:13],'',cox5[,14:17])

# # Any new onset
# cox_isch <- coxtab(a1[a1_rad_l,],'ischemic_heart_disease_grp_rec',prevcvd)
# cox_stroke <- coxtab(a1[a1_rad_l,],'stroke_tia_grp_rec',prevcvd)
# cox_chf <- coxtab(a1[a1_rad_l,],'cardiomyopathy_heart_failure_grp_rec',prevcvd)
# cox_cvdcombo <- coxtab(a1[a1_rad_l,],'cvdcombo_grp_rec',prevcvd)
# cox5rec <- Reduce(function(...) merge(...,by='var', all=T),
#                   list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
# cox5rec$var <- factor(cox5rec$var, levels=levels(cox_isch$var))
# cox5rec[,-1] <- lapply(cox5rec[,-1], as.character)
# cox5rec <- cox5rec[order(cox5rec$var),]
# cox5rec <- rbind(cox5rec[1,],'',cox5rec[2,],'',cox5rec[3:7,],cox5rec[8,],'',cox5rec[9:11,],cox5rec[12,],'',
#                  cox5rec[13:14,],'',cox5rec[15:16,],'',cox5rec[17:20,],'',cox5rec[21:24,],'',
#                  cox5rec[25:27,],'',cox5rec[28:34,])
# cox5rec <- cbind(cox5rec[,1:5],'',cox5rec[,6:9],'',cox5rec[,10:13],'',cox5rec[,14:17])
# 

# # right side tumor, rad
# # Pure incidence
# cox_isch <- coxtab(a1[a1_ihd_rad_r,],'ischemic_heart_disease_grp_inc',prevcvd[-3])
# cox_stroke <- coxtab(a1[a1_stroke_rad_r,],'stroke_tia_grp_inc',prevcvd[-5])
# cox_chf <- coxtab(a1[a1_chf_rad_r,],'cardiomyopathy_heart_failure_grp_inc',prevcvd[-2])
# cox_cvdcombo <- coxtab(a1[a1_combo_rad_r,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-5)])
# cox6 <- Reduce(function(...) merge(...,by='var', all=T),
#                list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
# cox6[,-1] <- lapply(cox6[,-1], as.character)
# cox6 <- rbind(cox6[1,],'',cox6[2,],'',cox6[3:7,],cox6[8,],'',cox6[9:11,],cox6[12,],'',
#               cox6[13:14,],'',cox6[15:16,],'',cox6[17:20,],'',cox6[21:24,],'',
#               cox6[25:27,],'',cox6[28:34,])
# cox6 <- cbind(cox6[,1:5],'',cox6[,6:9],'',cox6[,10:13],'',cox6[,14:17])

# # Any new onset
# cox_isch <- coxtab(a1[a1_rad_r,],'ischemic_heart_disease_grp_rec',prevcvd)
# cox_stroke <- coxtab(a1[a1_rad_r,],'stroke_tia_grp_rec',prevcvd)
# cox_chf <- coxtab(a1[a1_rad_r,],'cardiomyopathy_heart_failure_grp_rec',prevcvd)
# cox_cvdcombo <- coxtab(a1[a1_rad_r,],'cvdcombo_grp_rec',prevcvd)
# cox6rec <- Reduce(function(...) merge(...,by='var', all=T),
#                   list(cox_isch, cox_stroke, cox_chf, cox_cvdcombo))
# cox6rec$var <- factor(cox6rec$var, levels=levels(cox_isch$var))
# cox6rec[,-1] <- lapply(cox6rec[,-1], as.character)
# cox6rec <- cox6rec[order(cox6rec$var),]
# cox6rec <- rbind(cox6rec[1,],'',cox6rec[2,],'',cox6rec[3:7,],cox6rec[8,],'',cox6rec[9:11,],cox6rec[12,],'',
#                  cox6rec[13:14,],'',cox6rec[15:16,],'',cox6rec[17:20,],'',cox6rec[21:24,],'',
#                  cox6rec[25:27,],'',cox6rec[28:34,])
# cox6rec <- cbind(cox6rec[,1:5],'',cox6rec[,6:9],'',cox6rec[,10:13],'',cox6rec[,14:17])


# export tables
write_csv(data.frame(cox1),'cox_all_incid.csv')
write_csv(data.frame(cox2),'cox_chemo_incid.csv')
#write_csv(cox2rec,'cox_chemo_newonset.csv')
write_csv(data.frame(cox3),'cox_horm_incid.csv')
#write_csv(cox3rec,'cox_horm_newonset.csv')
write_csv(data.frame(cox4),'cox_rad_incid.csv')
#write_csv(cox4rec,'cox_rad_newonset.csv')
#write_csv(cox5,'cox_rad_left_incid.csv')
#write_csv(cox5rec,'cox_rad_left_newonset.csv')
#write_csv(cox6,'cox_rad_right_incid.csv')
#write_csv(cox6rec,'cox_rad_right_newonset.csv')

cox_all <- data.frame(cbind(cox1, '', cox2,'', cox3,'', cox4,'',cox5))
write_csv(cox_all,'cox_all.csv')


# New Table 5
# chemo
cox_isch <- coxtab(a1[a1_ihd_chemo_only,],'ischemic_heart_disease_grp_inc',prevcvd[-4])
cox_stroke <- coxtab(a1[a1_stroke_chemo_only,],'stroke_grp_inc',prevcvd[-6])
cox_hf <- coxtab(a1[a1_hf_chemo_only,],'heart_failure_grp_inc',prevcvd[-3])
cox_cm <- coxtab(a1[a1_cm_chemo_only,],'cardiomyopathy_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_chemo_only,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)])
cox_arrhythmia <- coxtab(a1[a1_arrhythmia_chemo_only,],'arrhythmia_grp_inc',prevcvd[c(-1)])
cox_cardiac <- coxtab(a1[a1_cardiac_chemo_only,],'cardiac_arrest_grp_inc',prevcvd[c(-1)])
cox_carotid <- coxtab(a1[a1_carotid_chemo_only,],'carotid_disease_grp_inc',prevcvd[c(-1)])
cox_myocarditis <- coxtab(a1[a1_myocarditis_chemo_only,],'myocarditis_pericarditis_grp_inc',prevcvd[c(-1)])
cox_tia <- coxtab(a1[a1_tia_chemo_only,],'tia_grp_inc',prevcvd[c(-1)])
cox_valvular <- coxtab(a1[a1_valvular_chemo_only,],'valvular_disease_grp_inc',prevcvd[c(-1)])
cox_dvt <- coxtab(a1[a1_dvt_chemo_only,],'venous_thromboembolic_disease_grp_inc',prevcvd[c(-1)])
cox_diab <- coxtab_diab(a1[a1_diab_chemo_only,],'cvdrf_diab',prevcvd)
cox_htn <- coxtab_htn(a1[a1_htn_chemo_only,],'cvdrf_htn',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_chemo_only,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_chemo_only,],'cvdrfcombo',prevcvd)
cox_allcvd <- coxtab(a1[a1_allcvd_chemo_only,],'allcvd_inc',prevcvd[c(-1,-2,-3,-4,-6)])
cox_hf_cm <- coxtab(a1[a1_hf_cm_chemo_only,],'hfcm_grp_inc',prevcvd[c(-2,-3)])


cox1 <- rbind(cox_isch[1:2,],cox_hf[1:2,], cox_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],"",
              cox_arrhythmia[1:2,], cox_cardiac[1:2,], cox_carotid[1:2,], 
              cox_myocarditis[1:2,], cox_tia[1:2,], cox_valvular[1:2,], cox_dvt[1:2,],"",
              cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,], cox_allcvd[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round2(cox1[,2],2),paste0("(",round2(cox1[,3],2),", ", round2(cox1[,4],2),")"))



# lest-sided radiation 
cox_isch <- coxtab(a1[a1_ihd_rad_l_only,],'ischemic_heart_disease_grp_inc',prevcvd[-4])
cox_stroke <- coxtab(a1[a1_stroke_rad_l_only,],'stroke_grp_inc',prevcvd[-6])
cox_hf <- coxtab(a1[a1_hf_rad_l_only,],'heart_failure_grp_inc',prevcvd[-3])
cox_cm <- coxtab(a1[a1_cm_rad_l_only,],'cardiomyopathy_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_rad_l_only,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)])
cox_arrhythmia <- coxtab(a1[a1_arrhythmia_rad_l_only,],'arrhythmia_grp_inc',prevcvd[c(-1)])
cox_cardiac <- coxtab(a1[a1_cardiac_rad_l_only,],'cardiac_arrest_grp_inc',prevcvd[c(-1)])
cox_carotid <- coxtab(a1[a1_carotid_rad_l_only,],'carotid_disease_grp_inc',prevcvd[c(-1)])
cox_myocarditis <- coxtab(a1[a1_myocarditis_rad_l_only,],'myocarditis_pericarditis_grp_inc',prevcvd[c(-1)])
cox_tia <- coxtab(a1[a1_tia_rad_l_only,],'tia_grp_inc',prevcvd[c(-1)])
cox_valvular <- coxtab(a1[a1_valvular_rad_l_only,],'valvular_disease_grp_inc',prevcvd[c(-1)])
cox_dvt <- coxtab(a1[a1_dvt_rad_l_only,],'venous_thromboembolic_disease_grp_inc',prevcvd[c(-1)])
cox_diab <- coxtab_diab(a1[a1_diab_rad_l_only,],'cvdrf_diab',prevcvd)
cox_htn <- coxtab_htn(a1[a1_htn_rad_l_only,],'cvdrf_htn',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_rad_l_only,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_rad_l_only,],'cvdrfcombo',prevcvd)
cox_allcvd <- coxtab(a1[a1_allcvd_rad_l_only,],'allcvd_inc',prevcvd[c(-1,-2,-3,-4,-6)])
cox_hf_cm <- coxtab(a1[a1_hf_cm_rad_l_only,],'hfcm_grp_inc',prevcvd[c(-2,-3)])


cox2 <- rbind(cox_isch[1:2,],cox_hf[1:2,], cox_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],"",
              cox_arrhythmia[1:2,], cox_cardiac[1:2,], cox_carotid[1:2,], 
              cox_myocarditis[1:2,], cox_tia[1:2,], cox_valvular[1:2,], cox_dvt[1:2,], "",
              cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,], cox_allcvd[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round2(cox2[,2],2),paste0("(",round2(cox2[,3],2),", ", round2(cox2[,4],2),")"))


# hormonal 
cox_isch <- coxtab(a1[a1_ihd_horm_only,],'ischemic_heart_disease_grp_inc',prevcvd[-4])
cox_stroke <- coxtab(a1[a1_stroke_horm_only,],'stroke_grp_inc',prevcvd[-6])
cox_hf <- coxtab(a1[a1_hf_horm_only,],'heart_failure_grp_inc',prevcvd[-3])
cox_cm <- coxtab(a1[a1_cm_horm_only,],'cardiomyopathy_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_horm_only,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)])
cox_arrhythmia <- coxtab(a1[a1_arrhythmia_horm_only,],'arrhythmia_grp_inc',prevcvd[c(-1)])
cox_cardiac <- coxtab(a1[a1_cardiac_horm_only,],'cardiac_arrest_grp_inc',prevcvd[c(-1)])
cox_carotid <- coxtab(a1[a1_carotid_horm_only,],'carotid_disease_grp_inc',prevcvd[c(-1)])
cox_myocarditis <- coxtab(a1[a1_myocarditis_horm_only,],'myocarditis_pericarditis_grp_inc',prevcvd[c(-1)])
cox_tia <- coxtab(a1[a1_tia_horm_only,],'tia_grp_inc',prevcvd[c(-1)])
cox_valvular <- coxtab(a1[a1_valvular_horm_only,],'valvular_disease_grp_inc',prevcvd[c(-1)])
cox_dvt <- coxtab(a1[a1_dvt_horm_only,],'venous_thromboembolic_disease_grp_inc',prevcvd[c(-1)])
cox_diab <- coxtab_diab(a1[a1_diab_horm_only,],'cvdrf_diab',prevcvd)
cox_htn <- coxtab_htn(a1[a1_htn_horm_only,],'cvdrf_htn',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_horm_only,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_horm_only,],'cvdrfcombo',prevcvd)
cox_allcvd <- coxtab(a1[a1_allcvd_horm_only,],'allcvd_inc',prevcvd[c(-1,-2,-3,-4,-6)])
cox_hf_cm <- coxtab(a1[a1_hf_cm_horm_only,],'hfcm_grp_inc',prevcvd[c(-2,-3)])


cox3 <- rbind(cox_isch[1:2,],cox_hf[1:2,], cox_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],"",
              cox_arrhythmia[1:2,], cox_cardiac[1:2,], cox_carotid[1:2,], 
              cox_myocarditis[1:2,], cox_tia[1:2,], cox_valvular[1:2,], cox_dvt[1:2,], "",
              cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,],cox_allcvd[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
cox3 <- cbind(round2(cox3[,2],2),paste0("(",round2(cox3[,3],2),", ", round2(cox3[,4],2),")"))




# chemo and lest-sided radiation
cox_isch <- coxtab(a1[a1_ihd_chemo_rad_l,],'ischemic_heart_disease_grp_inc',prevcvd[-4])
cox_stroke <- coxtab(a1[a1_stroke_chemo_rad_l,],'stroke_grp_inc',prevcvd[-6])
cox_hf <- coxtab(a1[a1_hf_chemo_rad_l,],'heart_failure_grp_inc',prevcvd[-3])
cox_cm <- coxtab(a1[a1_cm_chemo_rad_l,],'cardiomyopathy_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_chemo_rad_l,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)])
cox_arrhythmia <- coxtab(a1[a1_arrhythmia_chemo_rad_l,],'arrhythmia_grp_inc',prevcvd[c(-1)])
cox_cardiac <- coxtab(a1[a1_cardiac_chemo_rad_l,],'cardiac_arrest_grp_inc',prevcvd[c(-1)])
cox_carotid <- coxtab(a1[a1_carotid_chemo_rad_l,],'carotid_disease_grp_inc',prevcvd[c(-1)])
cox_myocarditis <- coxtab(a1[a1_myocarditis_chemo_rad_l,],'myocarditis_pericarditis_grp_inc',prevcvd[c(-1)])
cox_tia <- coxtab(a1[a1_tia_chemo_rad_l,],'tia_grp_inc',prevcvd[c(-1)])
cox_valvular <- coxtab(a1[a1_valvular_chemo_rad_l,],'valvular_disease_grp_inc',prevcvd[c(-1)])
cox_dvt <- coxtab(a1[a1_dvt_chemo_rad_l,],'venous_thromboembolic_disease_grp_inc',prevcvd[c(-1)])
cox_diab <- coxtab_diab(a1[a1_diab_chemo_rad_l,],'cvdrf_diab',prevcvd)
cox_htn <- coxtab_htn(a1[a1_htn_chemo_rad_l,],'cvdrf_htn',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_chemo_rad_l,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_chemo_rad_l,],'cvdrfcombo',prevcvd)
cox_allcvd <- coxtab(a1[a1_allcvd_chemo_rad_l,],'allcvd_inc',prevcvd[c(-1,-2,-3,-4,-6)])
cox_hf_cm <- coxtab(a1[a1_hf_cm_chemo_rad_l,],'hfcm_grp_inc',prevcvd[c(-2,-3)])

cox4 <- rbind(cox_isch[1:2,],cox_hf[1:2,], cox_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],"",
              cox_arrhythmia[1:2,], cox_cardiac[1:2,], cox_carotid[1:2,], 
              cox_myocarditis[1:2,], cox_tia[1:2,], cox_valvular[1:2,], cox_dvt[1:2,], "",
              cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,], cox_allcvd[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
cox4 <- cbind(round2(cox4[,2],2),paste0("(",round2(cox4[,3],2),", ", round2(cox4[,4],2),")"))



# chemo and horm
cox_isch <- coxtab(a1[a1_ihd_chemo_horm,],'ischemic_heart_disease_grp_inc',prevcvd[-4])
cox_stroke <- coxtab(a1[a1_stroke_chemo_horm,],'stroke_grp_inc',prevcvd[-6])
cox_hf <- coxtab(a1[a1_hf_chemo_horm,],'heart_failure_grp_inc',prevcvd[-3])
cox_cm <- coxtab(a1[a1_cm_chemo_horm,],'cardiomyopathy_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_chemo_horm,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)])
cox_arrhythmia <- coxtab(a1[a1_arrhythmia_chemo_horm,],'arrhythmia_grp_inc',prevcvd[c(-1)])
cox_cardiac <- coxtab(a1[a1_cardiac_chemo_horm,],'cardiac_arrest_grp_inc',prevcvd[c(-1)])
cox_carotid <- coxtab(a1[a1_carotid_chemo_horm,],'carotid_disease_grp_inc',prevcvd[c(-1)])
cox_myocarditis <- coxtab(a1[a1_myocarditis_chemo_horm,],'myocarditis_pericarditis_grp_inc',prevcvd[c(-1)])
cox_tia <- coxtab(a1[a1_tia_chemo_horm,],'tia_grp_inc',prevcvd[c(-1)])
cox_valvular <- coxtab(a1[a1_valvular_chemo_horm,],'valvular_disease_grp_inc',prevcvd[c(-1)])
cox_dvt <- coxtab(a1[a1_dvt_chemo_horm,],'venous_thromboembolic_disease_grp_inc',prevcvd[c(-1)])
cox_diab <- coxtab_diab(a1[a1_diab_chemo_horm,],'cvdrf_diab',prevcvd)
cox_htn <- coxtab_htn(a1[a1_htn_chemo_horm,],'cvdrf_htn',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_chemo_horm,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_chemo_horm,],'cvdrfcombo',prevcvd)
cox_allcvd <- coxtab(a1[a1_allcvd_chemo_horm,],'allcvd_inc',prevcvd[c(-1,-2,-3,-4,-6)])
cox_hf_cm <- coxtab(a1[a1_hf_cm_chemo_horm,],'hfcm_grp_inc',prevcvd[c(-2,-3)])


cox5 <- rbind(cox_isch[1:2,],cox_hf[1:2,], cox_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],"",
              cox_arrhythmia[1:2,], cox_cardiac[1:2,], cox_carotid[1:2,], 
              cox_myocarditis[1:2,], cox_tia[1:2,], cox_valvular[1:2,], cox_dvt[1:2,], "",
              cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,],cox_allcvd[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
cox5 <- cbind(round2(cox5[,2],2),paste0("(",round2(cox5[,3],2),", ", round2(cox5[,4],2),")"))



# left-sided radiation and horm
cox_isch <- coxtab(a1[a1_ihd_rad_l_horm,],'ischemic_heart_disease_grp_inc',prevcvd[-4])
cox_stroke <- coxtab(a1[a1_stroke_rad_l_horm,],'stroke_grp_inc',prevcvd[-6])
cox_hf <- coxtab(a1[a1_hf_rad_l_horm,],'heart_failure_grp_inc',prevcvd[-3])
cox_cm <- coxtab(a1[a1_cm_rad_l_horm,],'cardiomyopathy_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_rad_l_horm,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)])
cox_arrhythmia <- coxtab(a1[a1_arrhythmia_rad_l_horm,],'arrhythmia_grp_inc',prevcvd[c(-1)])
cox_cardiac <- coxtab(a1[a1_cardiac_rad_l_horm,],'cardiac_arrest_grp_inc',prevcvd[c(-1)])
cox_carotid <- coxtab(a1[a1_carotid_rad_l_horm,],'carotid_disease_grp_inc',prevcvd[c(-1)])
cox_myocarditis <- coxtab(a1[a1_myocarditis_rad_l_horm,],'myocarditis_pericarditis_grp_inc',prevcvd[c(-1)])
cox_tia <- coxtab(a1[a1_tia_rad_l_horm,],'tia_grp_inc',prevcvd[c(-1)])
cox_valvular <- coxtab(a1[a1_valvular_rad_l_horm,],'valvular_disease_grp_inc',prevcvd[c(-1)])
cox_dvt <- coxtab(a1[a1_dvt_rad_l_horm,],'venous_thromboembolic_disease_grp_inc',prevcvd[c(-1)])
cox_diab <- coxtab_diab(a1[a1_diab_rad_l_horm,],'cvdrf_diab',prevcvd)
cox_htn <- coxtab_htn(a1[a1_htn_rad_l_horm,],'cvdrf_htn',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_rad_l_horm,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_rad_l_horm,],'cvdrfcombo',prevcvd)
cox_allcvd <- coxtab(a1[a1_allcvd_rad_l_horm,],'allcvd_inc',prevcvd[c(-1,-2,-3,-4,-6)])
cox_hf_cm <- coxtab(a1[a1_hf_cm_rad_l_horm,],'hfcm_grp_inc',prevcvd[c(-2,-3)])


cox6 <- rbind(cox_isch[1:2,],cox_hf[1:2,], cox_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],"",
              cox_arrhythmia[1:2,], cox_cardiac[1:2,], cox_carotid[1:2,], 
              cox_myocarditis[1:2,], cox_tia[1:2,], cox_valvular[1:2,], cox_dvt[1:2,], "",
              cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,], cox_allcvd[1:2,])
cox6[,-1] <- lapply(cox6[,-1], function(x) as.numeric(as.character(x)))
cox6 <- cbind(round2(cox6[,2],2),paste0("(",round2(cox6[,3],2),", ", round2(cox6[,4],2),")"))

# chemo, left-sided radiation and horm
cox_isch <- coxtab(a1[a1_ihd_chemo_horm_rad_l,],'ischemic_heart_disease_grp_inc',prevcvd[-4])
cox_stroke <- coxtab(a1[a1_stroke_chemo_horm_rad_l,],'stroke_grp_inc',prevcvd[-6])
cox_hf <- coxtab(a1[a1_hf_chemo_horm_rad_l,],'heart_failure_grp_inc',prevcvd[-3])
cox_cm <- coxtab(a1[a1_cm_chemo_horm_rad_l,],'cardiomyopathy_grp_inc',prevcvd[-2])
cox_cvdcombo <- coxtab(a1[a1_combo_chemo_horm_rad_l,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)])
cox_arrhythmia <- coxtab(a1[a1_arrhythmia_chemo_horm_rad_l,],'arrhythmia_grp_inc',prevcvd[c(-1)])
cox_cardiac <- coxtab(a1[a1_cardiac_chemo_horm_rad_l,],'cardiac_arrest_grp_inc',prevcvd[c(-1)])
cox_carotid <- coxtab(a1[a1_carotid_chemo_horm_rad_l,],'carotid_disease_grp_inc',prevcvd[c(-1)])
cox_myocarditis <- coxtab(a1[a1_myocarditis_chemo_horm_rad_l,],'myocarditis_pericarditis_grp_inc',prevcvd[c(-1)])
cox_tia <- coxtab(a1[a1_tia_chemo_horm_rad_l,],'tia_grp_inc',prevcvd[c(-1)])
cox_valvular <- coxtab(a1[a1_valvular_chemo_horm_rad_l,],'valvular_disease_grp_inc',prevcvd[c(-1)])
cox_dvt <- coxtab(a1[a1_dvt_chemo_horm_rad_l,],'venous_thromboembolic_disease_grp_inc',prevcvd[c(-1)])
cox_diab <- coxtab_diab(a1[a1_diab_chemo_horm_rad_l,],'cvdrf_diab',prevcvd)
cox_htn <- coxtab_htn(a1[a1_htn_chemo_horm_rad_l,],'cvdrf_htn',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_chemo_horm_rad_l,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_chemo_horm_rad_l,],'cvdrfcombo',prevcvd)
cox_allcvd <- coxtab(a1[a1_allcvd_chemo_horm_rad_l,],'allcvd_inc',prevcvd[c(-1,-2,-3,-4,-6)])
cox_hf_cm <- coxtab(a1[a1_hf_cm_chemo_horm_rad_l,],'hfcm_grp_inc',prevcvd[c(-2,-3)])


cox7 <- rbind(cox_isch[1:2,],cox_hf[1:2,], cox_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],"",
              cox_arrhythmia[1:2,], cox_cardiac[1:2,], cox_carotid[1:2,], 
              cox_myocarditis[1:2,], cox_tia[1:2,], cox_valvular[1:2,], cox_dvt[1:2,], "",
              cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,], cox_allcvd[1:2,])
cox7[,-1] <- lapply(cox7[,-1], function(x) as.numeric(as.character(x)))
cox7 <- cbind(round2(cox7[,2],2),paste0("(",round2(cox7[,3],2),", ", round2(cox7[,4],2),")"))


cox_tab5 <- data.frame(cbind(cox1, '', cox2,'', cox3,'', cox4,'',cox5, '', cox6, '', cox7))
write_csv(cox_tab5,'cox_tab5.csv')

coxns =  data.frame(t(rbind(c(length(a1_chemo_only),sapply(list(a1_ihd_chemo_only, 
                                                                a1_hf_chemo_only,
                                                                a1_cm_chemo_only,
                                                                a1_stroke_chemo_only,
                                                                a1_combo2_chemo_only,
                                                                a1_arrhythmia_chemo_only,
                                                                a1_cardiac_chemo_only,
                                                                a1_carotid_chemo_only,
                                                                a1_myocarditis_chemo_only,
                                                                a1_tia_chemo_only,
                                                                a1_valvular_chemo_only,
                                                                a1_dvt_chemo_only,
                                                                a1_htn_chemo_only,
                                                                a1_diab_chemo_only,
                                                                a1_dyslipid_chemo_only,
                                                                a1_rfcombo_chemo_only,
                                                                a1_allcvd_chemo_only),length)),
                            c(length(a1_chemo_only_cases),sapply(list(a1_ihd_chemo_only_cases, 
                                                                a1_hf_chemo_only_cases,
                                                                a1_cm_chemo_only_cases,
                                                                a1_stroke_chemo_only_cases,
                                                                a1_combo2_chemo_only_cases,
                                                                a1_arrhythmia_chemo_only_cases,
                                                                a1_cardiac_chemo_only_cases,
                                                                a1_carotid_chemo_only_cases,
                                                                a1_myocarditis_chemo_only_cases,
                                                                a1_tia_chemo_only_cases,
                                                                a1_valvular_chemo_only_cases,
                                                                a1_dvt_chemo_only_cases,
                                                                a1_htn_chemo_only_cases,
                                                                a1_diab_chemo_only_cases,
                                                                a1_dyslipid_chemo_only_cases,
                                                                a1_rfcombo_chemo_only_cases,
                                                                a1_allcvd_chemo_only_cases),length)),
                            c(length(a1_chemo_only_controls),sapply(list(a1_ihd_chemo_only_controls, 
                                                                a1_hf_chemo_only_controls,
                                                                a1_cm_chemo_only_controls,
                                                                a1_stroke_chemo_only_controls,
                                                                a1_combo2_chemo_only_controls,
                                                                a1_arrhythmia_chemo_only_controls,
                                                                a1_cardiac_chemo_only_controls,
                                                                a1_carotid_chemo_only_controls,
                                                                a1_myocarditis_chemo_only_controls,
                                                                a1_tia_chemo_only_controls,
                                                                a1_valvular_chemo_only_controls,
                                                                a1_dvt_chemo_only_controls,
                                                                a1_htn_chemo_only_controls,
                                                                a1_diab_chemo_only_controls,
                                                                a1_dyslipid_chemo_only_controls,
                                                                a1_rfcombo_chemo_only_controls,
                                                                a1_allcvd_chemo_only_controls),length)),
                            c(length(a1_rad_l_only),sapply(list(a1_ihd_rad_l_only, 
                                                                a1_hf_rad_l_only,
                                                                a1_cm_rad_l_only,
                                                                a1_stroke_rad_l_only,
                                                                a1_combo2_rad_l_only,
                                                                a1_arrhythmia_rad_l_only,
                                                                a1_cardiac_rad_l_only,
                                                                a1_carotid_rad_l_only,
                                                                a1_myocarditis_rad_l_only,
                                                                a1_tia_rad_l_only,
                                                                a1_valvular_rad_l_only,
                                                                a1_dvt_rad_l_only,
                                                                a1_htn_rad_l_only,
                                                                a1_diab_rad_l_only,
                                                                a1_dyslipid_rad_l_only,
                                                                a1_rfcombo_rad_l_only,
                                                                a1_allcvd_rad_l_only),length)),
                            c(length(a1_rad_l_only_cases),sapply(list(a1_ihd_rad_l_only_cases, 
                                                                      a1_hf_rad_l_only_cases,
                                                                      a1_cm_rad_l_only_cases,
                                                                      a1_stroke_rad_l_only_cases,
                                                                      a1_combo2_rad_l_only_cases,
                                                                      a1_arrhythmia_rad_l_only_cases,
                                                                      a1_cardiac_rad_l_only_cases,
                                                                      a1_carotid_rad_l_only_cases,
                                                                      a1_myocarditis_rad_l_only_cases,
                                                                      a1_tia_rad_l_only_cases,
                                                                      a1_valvular_rad_l_only_cases,
                                                                      a1_dvt_rad_l_only_cases,
                                                                      a1_htn_rad_l_only_cases,
                                                                      a1_diab_rad_l_only_cases,
                                                                      a1_dyslipid_rad_l_only_cases,
                                                                      a1_rfcombo_rad_l_only_cases,
                                                                      a1_allcvd_rad_l_only_cases),length)),
                            c(length(a1_rad_l_only_controls),sapply(list(a1_ihd_rad_l_only_controls, 
                                                                         a1_hf_rad_l_only_controls,
                                                                         a1_cm_rad_l_only_controls,
                                                                         a1_stroke_rad_l_only_controls,
                                                                         a1_combo2_rad_l_only_controls,
                                                                         a1_arrhythmia_rad_l_only_controls,
                                                                         a1_cardiac_rad_l_only_controls,
                                                                         a1_carotid_rad_l_only_controls,
                                                                         a1_myocarditis_rad_l_only_controls,
                                                                         a1_tia_rad_l_only_controls,
                                                                         a1_valvular_rad_l_only_controls,
                                                                         a1_dvt_rad_l_only_controls,
                                                                         a1_htn_rad_l_only_controls,
                                                                         a1_diab_rad_l_only_controls,
                                                                         a1_dyslipid_rad_l_only_controls,
                                                                         a1_rfcombo_rad_l_only_controls,
                                                                         a1_allcvd_rad_l_only_controls),length)),
                            c(length(a1_horm_only),sapply(list(a1_ihd_horm_only, 
                                                               a1_hf_horm_only,
                                                               a1_cm_horm_only,
                                                               a1_stroke_horm_only,
                                                               a1_combo2_horm_only,
                                                               a1_arrhythmia_horm_only,
                                                               a1_cardiac_horm_only,
                                                               a1_carotid_horm_only,
                                                               a1_myocarditis_horm_only,
                                                               a1_tia_horm_only,
                                                               a1_valvular_horm_only,
                                                               a1_dvt_horm_only,
                                                               a1_htn_horm_only,
                                                               a1_diab_horm_only,
                                                               a1_dyslipid_horm_only,
                                                               a1_rfcombo_horm_only,
                                                               a1_allcvd_horm_only),length)),
                            c(length(a1_horm_only_cases),sapply(list(a1_ihd_horm_only_cases, 
                                                                     a1_hf_horm_only_cases,
                                                                     a1_cm_horm_only_cases,
                                                                     a1_stroke_horm_only_cases,
                                                                     a1_combo2_horm_only_cases,
                                                                     a1_arrhythmia_horm_only_cases,
                                                                     a1_cardiac_horm_only_cases,
                                                                     a1_carotid_horm_only_cases,
                                                                     a1_myocarditis_horm_only_cases,
                                                                     a1_tia_horm_only_cases,
                                                                     a1_valvular_horm_only_cases,
                                                                     a1_dvt_horm_only_cases,
                                                                     a1_htn_horm_only_cases,
                                                                     a1_diab_horm_only_cases,
                                                                     a1_dyslipid_horm_only_cases,
                                                                     a1_rfcombo_horm_only_cases,
                                                                     a1_allcvd_horm_only_cases),length)),
                            c(length(a1_horm_only_controls),sapply(list(a1_ihd_horm_only_controls, 
                                                                        a1_hf_horm_only_controls,
                                                                        a1_cm_horm_only_controls,
                                                                        a1_stroke_horm_only_controls,
                                                                        a1_combo2_horm_only_controls,
                                                                        a1_arrhythmia_horm_only_controls,
                                                                        a1_cardiac_horm_only_controls,
                                                                        a1_carotid_horm_only_controls,
                                                                        a1_myocarditis_horm_only_controls,
                                                                        a1_tia_horm_only_controls,
                                                                        a1_valvular_horm_only_controls,
                                                                        a1_dvt_horm_only_controls,
                                                                        a1_htn_horm_only_controls,
                                                                        a1_diab_horm_only_controls,
                                                                        a1_dyslipid_horm_only_controls,
                                                                        a1_rfcombo_horm_only_controls,
                                                                        a1_allcvd_horm_only_controls),length)),
                            c(length(a1_chemo_rad_l),sapply(list(a1_ihd_chemo_rad_l, 
                                                                 a1_hf_chemo_rad_l,
                                                                 a1_cm_chemo_rad_l,
                                                                 a1_stroke_chemo_rad_l,
                                                                 a1_combo2_chemo_rad_l,
                                                                 a1_arrhythmia_chemo_rad_l,
                                                                 a1_cardiac_chemo_rad_l,
                                                                 a1_carotid_chemo_rad_l,
                                                                 a1_myocarditis_chemo_rad_l,
                                                                 a1_tia_chemo_rad_l,
                                                                 a1_valvular_chemo_rad_l,
                                                                 a1_dvt_chemo_rad_l,
                                                                 a1_htn_chemo_rad_l,
                                                                 a1_diab_chemo_rad_l,
                                                                 a1_dyslipid_chemo_rad_l,
                                                                 a1_rfcombo_chemo_rad_l,
                                                                 a1_allcvd_chemo_rad_l),length)),
                            c(length(a1_chemo_rad_l_cases),sapply(list(a1_ihd_chemo_rad_l_cases, 
                                                                       a1_hf_chemo_rad_l_cases,
                                                                       a1_cm_chemo_rad_l_cases,
                                                                       a1_stroke_chemo_rad_l_cases,
                                                                       a1_combo2_chemo_rad_l_cases,
                                                                       a1_arrhythmia_chemo_rad_l_cases,
                                                                       a1_cardiac_chemo_rad_l_cases,
                                                                       a1_carotid_chemo_rad_l_cases,
                                                                       a1_myocarditis_chemo_rad_l_cases,
                                                                       a1_tia_chemo_rad_l_cases,
                                                                       a1_valvular_chemo_rad_l_cases,
                                                                       a1_dvt_chemo_rad_l_cases,
                                                                       a1_htn_chemo_rad_l_cases,
                                                                       a1_diab_chemo_rad_l_cases,
                                                                       a1_dyslipid_chemo_rad_l_cases,
                                                                       a1_rfcombo_chemo_rad_l_cases,
                                                                       a1_allcvd_chemo_rad_l_cases),length)),
                            c(length(a1_chemo_rad_l_controls),sapply(list(a1_ihd_chemo_rad_l_controls, 
                                                                          a1_hf_chemo_rad_l_controls,
                                                                          a1_cm_chemo_rad_l_controls,
                                                                          a1_stroke_chemo_rad_l_controls,
                                                                          a1_combo2_chemo_rad_l_controls,
                                                                          a1_arrhythmia_chemo_rad_l_controls,
                                                                          a1_cardiac_chemo_rad_l_controls,
                                                                          a1_carotid_chemo_rad_l_controls,
                                                                          a1_myocarditis_chemo_rad_l_controls,
                                                                          a1_tia_chemo_rad_l_controls,
                                                                          a1_valvular_chemo_rad_l_controls,
                                                                          a1_dvt_chemo_rad_l_controls,
                                                                          a1_htn_chemo_rad_l_controls,
                                                                          a1_diab_chemo_rad_l_controls,
                                                                          a1_dyslipid_chemo_rad_l_controls,
                                                                          a1_rfcombo_chemo_rad_l_controls,
                                                                          a1_allcvd_chemo_rad_l_controls),length)),
                            c(length(a1_chemo_horm),sapply(list(a1_ihd_chemo_horm, 
                                                                a1_hf_chemo_horm,
                                                                a1_cm_chemo_horm,
                                                                a1_stroke_chemo_horm,
                                                                a1_combo2_chemo_horm,
                                                                a1_arrhythmia_chemo_horm,
                                                                a1_cardiac_chemo_horm,
                                                                a1_carotid_chemo_horm,
                                                                a1_myocarditis_chemo_horm,
                                                                a1_tia_chemo_horm,
                                                                a1_valvular_chemo_horm,
                                                                a1_dvt_chemo_horm,
                                                                a1_htn_chemo_horm,
                                                                a1_diab_chemo_horm,
                                                                a1_dyslipid_chemo_horm,
                                                                a1_rfcombo_chemo_horm,
                                                                a1_allcvd_chemo_horm),length)),
                            c(length(a1_chemo_horm_cases),sapply(list(a1_ihd_chemo_horm_cases, 
                                                                      a1_hf_chemo_horm_cases,
                                                                      a1_cm_chemo_horm_cases,
                                                                      a1_stroke_chemo_horm_cases,
                                                                      a1_combo2_chemo_horm_cases,
                                                                      a1_arrhythmia_chemo_horm_cases,
                                                                      a1_cardiac_chemo_horm_cases,
                                                                      a1_carotid_chemo_horm_cases,
                                                                      a1_myocarditis_chemo_horm_cases,
                                                                      a1_tia_chemo_horm_cases,
                                                                      a1_valvular_chemo_horm_cases,
                                                                      a1_dvt_chemo_horm_cases,
                                                                      a1_htn_chemo_horm_cases,
                                                                      a1_diab_chemo_horm_cases,
                                                                      a1_dyslipid_chemo_horm_cases,
                                                                      a1_rfcombo_chemo_horm_cases,
                                                                      a1_allcvd_chemo_horm_cases),length)),
                            c(length(a1_chemo_horm_controls),sapply(list(a1_ihd_chemo_horm_controls, 
                                                                         a1_hf_chemo_horm_controls,
                                                                         a1_cm_chemo_horm_controls,
                                                                         a1_stroke_chemo_horm_controls,
                                                                         a1_combo2_chemo_horm_controls,
                                                                         a1_arrhythmia_chemo_horm_controls,
                                                                         a1_cardiac_chemo_horm_controls,
                                                                         a1_carotid_chemo_horm_controls,
                                                                         a1_myocarditis_chemo_horm_controls,
                                                                         a1_tia_chemo_horm_controls,
                                                                         a1_valvular_chemo_horm_controls,
                                                                         a1_dvt_chemo_horm_controls,
                                                                         a1_htn_chemo_horm_controls,
                                                                         a1_diab_chemo_horm_controls,
                                                                         a1_dyslipid_chemo_horm_controls,
                                                                         a1_rfcombo_chemo_horm_controls,
                                                                         a1_allcvd_chemo_horm_controls),length)),
                            c(length(a1_rad_l_horm),sapply(list(a1_ihd_rad_l_horm, 
                                                                a1_hf_rad_l_horm,
                                                                a1_cm_rad_l_horm,
                                                                a1_stroke_rad_l_horm,
                                                                a1_combo2_rad_l_horm,
                                                                a1_arrhythmia_rad_l_horm,
                                                                a1_cardiac_rad_l_horm,
                                                                a1_carotid_rad_l_horm,
                                                                a1_myocarditis_rad_l_horm,
                                                                a1_tia_rad_l_horm,
                                                                a1_valvular_rad_l_horm,
                                                                a1_dvt_rad_l_horm,
                                                                a1_htn_rad_l_horm,
                                                                a1_diab_rad_l_horm,
                                                                a1_dyslipid_rad_l_horm,
                                                                a1_rfcombo_rad_l_horm,
                                                                a1_allcvd_rad_l_horm),length)),
                            c(length(a1_rad_l_horm_cases),sapply(list(a1_ihd_rad_l_horm_cases, 
                                                                      a1_hf_rad_l_horm_cases,
                                                                      a1_cm_rad_l_horm_cases,
                                                                      a1_stroke_rad_l_horm_cases,
                                                                      a1_combo2_rad_l_horm_cases,
                                                                      a1_arrhythmia_rad_l_horm_cases,
                                                                      a1_cardiac_rad_l_horm_cases,
                                                                      a1_carotid_rad_l_horm_cases,
                                                                      a1_myocarditis_rad_l_horm_cases,
                                                                      a1_tia_rad_l_horm_cases,
                                                                      a1_valvular_rad_l_horm_cases,
                                                                      a1_dvt_rad_l_horm_cases,
                                                                      a1_htn_rad_l_horm_cases,
                                                                      a1_diab_rad_l_horm_cases,
                                                                      a1_dyslipid_rad_l_horm_cases,
                                                                      a1_rfcombo_rad_l_horm_cases,
                                                                      a1_allcvd_rad_l_horm_cases),length)),
                            c(length(a1_rad_l_horm_controls),sapply(list(a1_ihd_rad_l_horm_controls, 
                                                                         a1_hf_rad_l_horm_controls,
                                                                         a1_cm_rad_l_horm_controls,
                                                                         a1_stroke_rad_l_horm_controls,
                                                                         a1_combo2_rad_l_horm_controls,
                                                                         a1_arrhythmia_rad_l_horm_controls,
                                                                         a1_cardiac_rad_l_horm_controls,
                                                                         a1_carotid_rad_l_horm_controls,
                                                                         a1_myocarditis_rad_l_horm_controls,
                                                                         a1_tia_rad_l_horm_controls,
                                                                         a1_valvular_rad_l_horm_controls,
                                                                         a1_dvt_rad_l_horm_controls,
                                                                         a1_htn_rad_l_horm_controls,
                                                                         a1_diab_rad_l_horm_controls,
                                                                         a1_dyslipid_rad_l_horm_controls,
                                                                         a1_rfcombo_rad_l_horm_controls,
                                                                         a1_allcvd_rad_l_horm_controls),length)),
                            c(length(a1_chemo_horm_rad_l),sapply(list(a1_ihd_chemo_horm_rad_l, 
                                                                      a1_hf_chemo_horm_rad_l,
                                                                      a1_cm_chemo_horm_rad_l,
                                                                      a1_stroke_chemo_horm_rad_l,
                                                                      a1_combo2_chemo_horm_rad_l,
                                                                      a1_arrhythmia_chemo_horm_rad_l,
                                                                      a1_cardiac_chemo_horm_rad_l,
                                                                      a1_carotid_chemo_horm_rad_l,
                                                                      a1_myocarditis_chemo_horm_rad_l,
                                                                      a1_tia_chemo_horm_rad_l,
                                                                      a1_valvular_chemo_horm_rad_l,
                                                                      a1_dvt_chemo_horm_rad_l,
                                                                      a1_htn_chemo_horm_rad_l,
                                                                      a1_diab_chemo_horm_rad_l,
                                                                      a1_dyslipid_chemo_horm_rad_l,
                                                                      a1_rfcombo_chemo_horm_rad_l,
                                                                      a1_allcvd_chemo_horm_rad_l),length)),
                            c(length(a1_chemo_horm_rad_l_cases),sapply(list(a1_ihd_chemo_horm_rad_l_cases, 
                                                                            a1_hf_chemo_horm_rad_l_cases,
                                                                            a1_cm_chemo_horm_rad_l_cases,
                                                                            a1_stroke_chemo_horm_rad_l_cases,
                                                                            a1_combo2_chemo_horm_rad_l_cases,
                                                                            a1_arrhythmia_chemo_horm_rad_l_cases,
                                                                            a1_cardiac_chemo_horm_rad_l_cases,
                                                                            a1_carotid_chemo_horm_rad_l_cases,
                                                                            a1_myocarditis_chemo_horm_rad_l_cases,
                                                                            a1_tia_chemo_horm_rad_l_cases,
                                                                            a1_valvular_chemo_horm_rad_l_cases,
                                                                            a1_dvt_chemo_horm_rad_l_cases,
                                                                            a1_htn_chemo_horm_rad_l_cases,
                                                                            a1_diab_chemo_horm_rad_l_cases,
                                                                            a1_dyslipid_chemo_horm_rad_l_cases,
                                                                            a1_rfcombo_chemo_horm_rad_l_cases,
                                                                            a1_allcvd_chemo_horm_rad_l_cases),length)),
                            c(length(a1_chemo_horm_rad_l_controls),sapply(list(a1_ihd_chemo_horm_rad_l_controls, 
                                                                               a1_hf_chemo_horm_rad_l_controls,
                                                                               a1_cm_chemo_horm_rad_l_controls,
                                                                               a1_stroke_chemo_horm_rad_l_controls,
                                                                               a1_combo2_chemo_horm_rad_l_controls,
                                                                               a1_arrhythmia_chemo_horm_rad_l_controls,
                                                                               a1_cardiac_chemo_horm_rad_l_controls,
                                                                               a1_carotid_chemo_horm_rad_l_controls,
                                                                               a1_myocarditis_chemo_horm_rad_l_controls,
                                                                               a1_tia_chemo_horm_rad_l_controls,
                                                                               a1_valvular_chemo_horm_rad_l_controls,
                                                                               a1_dvt_chemo_horm_rad_l_controls,
                                                                               a1_htn_chemo_horm_rad_l_controls,
                                                                               a1_diab_chemo_horm_rad_l_controls,
                                                                               a1_dyslipid_chemo_horm_rad_l_controls,
                                                                               a1_rfcombo_chemo_horm_rad_l_controls,
                                                                               a1_allcvd_chemo_horm_rad_l_controls),length))
                            
                            )))

coxns[,c(1,4,7,10,13,16,19)] <- lapply(coxns[,c(1,4,7,10,13,16,19)],function(x) paste0('Total: n = ',x))
coxns[,c(2,5,8,11,14,17,20)] <- lapply(coxns[,c(2,5,8,11,14,17,20)],function(x) paste0('Cases: n = ',x))
coxns[,c(3,6,9,12,15,18,21)] <- lapply(coxns[,c(3,6,9,12,15,18,21)],function(x) paste0('Controls: n = ',x))

names(coxns) <- c('chemo','cases','controls','left rad','cases','controls','horm','cases','controls','chemo and left rad','cases','controls',
                  'chemo and horm','cases','controls', 'left rad and horm','cases','controls', 'all therapy','cases','controls')
coxns$sample <- c('all','ihd','hf','cm','stroke','all parimary','arrhythmia','cardiac','carotid','myocarditis','tia','valvular','dvt','htn','diabetes', 'dyslipid','cvdrf combo','all cvd')

coxns_all = coxns[,c(1,4,7,10,13,16,19)]
coxns_cases_controls = coxns[,c(2,3,5,6,8,9,11,12,14,15,17,18,20,21)]


write_csv(coxns_all, 'sample_size_cox_tab5_1.csv')
write_csv(coxns_cases_controls, 'sample_size_cox_tab5_2.csv')




# CVD Death Analysis (Table 6)

# function to summary cox model results
coxtab <- function(d,x,covar = NULL){
  
  if (is.null(covar)){
    surv_object <- Surv(time = floor(d[,paste0(x,'_fu')]/30), event = d[,x])
    fit.coxph1 <- coxph(surv_object ~ group1, data = d)
    fit.coxph2 <- eval(parse(text=(paste0("coxph(surv_object ~ group1+bmicat1+menop+smok+
                        diab_bl+htn_bl+dyslipid_bl+medhousincome_q+edu1_q",", data = d)"))))
    sum1 <- c(summary(fit.coxph1)$coef[c(2,5)],summary(fit.coxph1)$conf.int[3:4])
    sum2 <- cbind(summary(fit.coxph2)$coef[,c(2,5)],summary(fit.coxph2)$conf.int[,3:4])
    sum3 <- rbind(sum1, sum2)[,c(1,3,4,2)]
    sum4 <- data.frame(cbind(var=row.names(sum3), sum3))
    sum4$var <- factor(sum4$var, levels=row.names(sum3))
   
  }
  if (!is.null(covar)){
    covars <- paste(covar, collapse = '+')
    surv_object <- Surv(time = floor(d[,paste0(x,'_fu')]/30), event = d[,x])
    fit.coxph1 <- coxph(surv_object ~ group1, data = d)
    fit.coxph2 <- eval(parse(text=(paste0("coxph(surv_object ~ group1+bmicat1+menop+smok+
                        diab_bl+htn_bl+dyslipid_bl+medhousincome_q+edu1_q+",covars,", data = d)"))))
    sum1 <- c(summary(fit.coxph1)$coef[c(2,5)],summary(fit.coxph1)$conf.int[3:4])
    sum2 <- cbind(summary(fit.coxph2)$coef[,c(2,5)],summary(fit.coxph2)$conf.int[,3:4])
    sum3 <- rbind(sum1, sum2)[,c(1,3,4,2)]
    sum4 <- data.frame(cbind(var=row.names(sum3), sum3))
    sum4$var <- factor(sum4$var, levels=row.names(sum3))
    
  }
  sum4
}


prevcvd <- c("arrhythmia_grp_prev","cardiomyopathy_grp_prev","heart_failure_grp_prev",
             "ischemic_heart_disease_grp_prev","myocarditis_pericarditis_grp_prev",    
             "stroke_grp_prev","tia_grp_prev","valvular_disease_grp_prev","venous_thromboembolic_disease_grp_prev")


# All
cox_primary_cvd <- coxtab(a1[a1_combo,],'primary_cvd_new', prevcvd[-c(2,3,4,6)])
cox_primary_general_cvd <- coxtab(a1[a1_combo,],'primary_general_cvd', prevcvd[-c(2,3,4,6)])
cox_primary_cvd_death <- coxtab(a1[a1_combo,],'primary_cvd_death', prevcvd[-c(2,3,4,6)])
cox_secondary_cvd <- coxtab(a1[a1_combo2,],'secondary_cvd_new', prevcvd[-c(1,5,7,8,9)])
cox_any_cvd <- coxtab(a1[a1_allcombo,],'any_cvd')
cox_any_death_cvd <- coxtab(a1[a1_allcombo,],'any_death_cvd')
cox_primary_secondary_cvd_death = coxtab(a1[a1_allcombo,],'primary_secondary_cvd_death')
cox1 <- rbind(cox_primary_cvd[1:2,],cox_primary_general_cvd[1:2,], cox_primary_cvd_death[1:2,], cox_secondary_cvd[1:2,], cox_any_cvd[1:2,],
              cox_any_death_cvd[1:2,], cox_primary_secondary_cvd_death[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round2(cox1[,2],2),paste0("(",round2(cox1[,3],2),", ", round2(cox1[,4],2),")"))

# chemo cases
cox_primary_cvd <- coxtab(a1[a1_combo_chemo,],'primary_cvd_new', prevcvd[-c(2,3,4,6)])
cox_primary_general_cvd <- coxtab(a1[a1_combo_chemo,],'primary_general_cvd', prevcvd[-c(2,3,4,6)])
cox_primary_cvd_death <- coxtab(a1[a1_combo_chemo,],'primary_cvd_death', prevcvd[-c(2,3,4,6)])
cox_secondary_cvd <- coxtab(a1[a1_combo2_chemo,],'secondary_cvd_new', prevcvd[-c(1,5,7,8,9)])
cox_any_cvd <- coxtab(a1[a1_allcombo_chemo,],'any_cvd')
cox_any_death_cvd <- coxtab(a1[a1_allcombo_chemo,],'any_death_cvd')
cox_primary_secondary_cvd_death = coxtab(a1[a1_allcombo_chemo,],'primary_secondary_cvd_death')
cox2 <- rbind(cox_primary_cvd[1:2,],cox_primary_general_cvd[1:2,], cox_primary_cvd_death[1:2,], cox_secondary_cvd[1:2,], cox_any_cvd[1:2,],
              cox_any_death_cvd[1:2,], cox_primary_secondary_cvd_death[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round2(cox2[,2],2),paste0("(",round2(cox2[,3],2),", ", round2(cox2[,4],2),")"))

# horm
cox_primary_cvd <- coxtab(a1[a1_combo_horm,],'primary_cvd_new', prevcvd[-c(2,3,4,6)])
cox_primary_general_cvd <- coxtab(a1[a1_combo_horm,],'primary_general_cvd', prevcvd[-c(2,3,4,6)])
cox_primary_cvd_death <- coxtab(a1[a1_combo_horm,],'primary_cvd_death', prevcvd[-c(2,3,4,6)])
cox_secondary_cvd <- coxtab(a1[a1_combo2_horm,],'secondary_cvd_new', prevcvd[-c(1,5,7,8,9)])
cox_any_cvd <- coxtab(a1[a1_allcombo_horm,],'any_cvd')
cox_any_death_cvd <- coxtab(a1[a1_allcombo_horm,],'any_death_cvd')
cox_primary_secondary_cvd_death = coxtab(a1[a1_allcombo_horm,],'primary_secondary_cvd_death')
cox3 <- rbind(cox_primary_cvd[1:2,],cox_primary_general_cvd[1:2,], cox_primary_cvd_death[1:2,], cox_secondary_cvd[1:2,], cox_any_cvd[1:2,],
              cox_any_death_cvd[1:2,], cox_primary_secondary_cvd_death[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
cox3 <- cbind(round2(cox3[,2],2),paste0("(",round2(cox3[,3],2),", ", round2(cox3[,4],2),")"))


# rad
cox_primary_cvd <- coxtab(a1[a1_combo_rad,],'primary_cvd_new', prevcvd[-c(2,3,4,6)])
cox_primary_general_cvd <- coxtab(a1[a1_combo_rad,],'primary_general_cvd', prevcvd[-c(2,3,4,6)])
cox_primary_cvd_death <- coxtab(a1[a1_combo_rad,],'primary_cvd_death', prevcvd[-c(2,3,4,6)])
cox_secondary_cvd <- coxtab(a1[a1_combo2_rad,],'secondary_cvd_new', prevcvd[-c(1,5,7,8,9)])
cox_any_cvd <- coxtab(a1[a1_allcombo_rad,],'any_cvd')
cox_any_death_cvd <- coxtab(a1[a1_allcombo_rad,],'any_death_cvd')
cox_primary_secondary_cvd_death = coxtab(a1[a1_allcombo_rad,],'primary_secondary_cvd_death')
cox4 <- rbind(cox_primary_cvd[1:2,],cox_primary_general_cvd[1:2,], cox_primary_cvd_death[1:2,], cox_secondary_cvd[1:2,], cox_any_cvd[1:2,],
              cox_any_death_cvd[1:2,], cox_primary_secondary_cvd_death[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
cox4 <- cbind(round2(cox4[,2],2),paste0("(",round2(cox4[,3],2),", ", round2(cox4[,4],2),")"))


# left-sided rad
cox_primary_cvd <- coxtab(a1[a1_combo_rad_l,],'primary_cvd_new', prevcvd[-c(2,3,4,6)])
cox_primary_general_cvd <- coxtab(a1[a1_combo_rad_l,],'primary_general_cvd', prevcvd[-c(2,3,4,6)])
cox_primary_cvd_death <- coxtab(a1[a1_combo_rad_l,],'primary_cvd_death', prevcvd[-c(2,3,4,6)])
cox_secondary_cvd <- coxtab(a1[a1_combo2_rad_l,],'secondary_cvd_new', prevcvd[-c(1,5,7,8,9)])
cox_any_cvd <- coxtab(a1[a1_allcombo_rad_l,],'any_cvd')
cox_any_death_cvd <- coxtab(a1[a1_allcombo_rad_l,],'any_death_cvd')
cox_primary_secondary_cvd_death = coxtab(a1[a1_allcombo_rad_l,],'primary_secondary_cvd_death')
cox5 <- rbind(cox_primary_cvd[1:2,],cox_primary_general_cvd[1:2,], cox_primary_cvd_death[1:2,], cox_secondary_cvd[1:2,], cox_any_cvd[1:2,],
              cox_any_death_cvd[1:2,], cox_primary_secondary_cvd_death[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
cox5 <- cbind(round2(cox5[,2],2),paste0("(",round2(cox5[,3],2),", ", round2(cox5[,4],2),")"))



cox_all <- data.frame(cbind(cox1, '', cox2,'', cox3,'', cox4,'',cox5))
write_csv(cox_all,'cox_all_death.csv')




### Table 6 Sup

# All
cox_primary_cvd <- coxtab(a1[a1_combo,],'primary_cvd', prevcvd[-c(2,3,4,6)])
cox_secondary_cvd <- coxtab(a1[a1_combo2,],'secondary_cvd', prevcvd[-c(1,5,7,8,9)])
cox_general <- coxtab(a1[a1_allcombo,],'general_cvd')
cox_death_cvd <- coxtab(a1[a1_allcombo,],'death_all')
cox1 <- rbind(cox_primary_cvd[1:2,],cox_secondary_cvd[1:2,], cox_general[1:2,],
              cox_death_cvd[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round2(cox1[,2],2),paste0("(",round2(cox1[,3],2),", ", round2(cox1[,4],2),")"))

# chemo only
cox_primary_cvd <- coxtab(a1[a1_combo_chemo_only,],'primary_cvd', prevcvd[-c(2,3,4,6)])
cox_secondary_cvd <- coxtab(a1[a1_combo2_chemo_only,],'secondary_cvd', prevcvd[-c(1,5,7,8,9)])
cox_general <- coxtab(a1[a1_allcombo_chemo_only,],'general_cvd')
cox_death_cvd <- coxtab(a1[a1_allcombo_chemo_only,],'death_all')
cox2 <- rbind(cox_primary_cvd[1:2,],cox_secondary_cvd[1:2,], cox_general[1:2,],
              cox_death_cvd[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round2(cox2[,2],2),paste0("(",round2(cox2[,3],2),", ", round2(cox2[,4],2),")"))

# rad_l only
cox_primary_cvd <- coxtab(a1[a1_combo_rad_l_only,],'primary_cvd', prevcvd[-c(2,3,4,6)])
cox_secondary_cvd <- coxtab(a1[a1_combo2_rad_l_only,],'secondary_cvd', prevcvd[-c(1,5,7,8,9)])
cox_general <- coxtab(a1[a1_allcombo_rad_l_only,],'general_cvd')
cox_death_cvd <- coxtab(a1[a1_allcombo_rad_l_only,],'death_all')
cox3 <- rbind(cox_primary_cvd[1:2,],cox_secondary_cvd[1:2,], cox_general[1:2,],
              cox_death_cvd[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
cox3 <- cbind(round2(cox3[,2],2),paste0("(",round2(cox3[,3],2),", ", round2(cox3[,4],2),")"))


# horm only
cox_primary_cvd <- coxtab(a1[a1_combo_horm_only,],'primary_cvd', prevcvd[-c(2,3,4,6)])
cox_secondary_cvd <- coxtab(a1[a1_combo2_horm_only,],'secondary_cvd', prevcvd[-c(1,5,7,8,9)])
cox_general <- coxtab(a1[a1_allcombo_horm_only,],'general_cvd')
cox_death_cvd <- coxtab(a1[a1_allcombo_horm_only,],'death_all')
cox4 <- rbind(cox_primary_cvd[1:2,],cox_secondary_cvd[1:2,], cox_general[1:2,],
              cox_death_cvd[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
cox4 <- cbind(round2(cox4[,2],2),paste0("(",round2(cox4[,3],2),", ", round2(cox4[,4],2),")"))


# chemo and rad_l
cox_primary_cvd <- coxtab(a1[a1_combo_chemo_rad_l,],'primary_cvd', prevcvd[-c(2,3,4,6)])
cox_secondary_cvd <- coxtab(a1[a1_combo2_chemo_rad_l,],'secondary_cvd', prevcvd[-c(1,5,7,8,9)])
cox_general <- coxtab(a1[a1_allcombo_chemo_rad_l,],'general_cvd')
cox_death_cvd <- coxtab(a1[a1_allcombo_chemo_rad_l,],'death_all')
cox5 <- rbind(cox_primary_cvd[1:2,],cox_secondary_cvd[1:2,], cox_general[1:2,],
              cox_death_cvd[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
cox5 <- cbind(round2(cox5[,2],2),paste0("(",round2(cox5[,3],2),", ", round2(cox5[,4],2),")"))


# chemo and horm
cox_primary_cvd <- coxtab(a1[a1_combo_chemo_horm,],'primary_cvd', prevcvd[-c(2,3,4,6)])
cox_secondary_cvd <- coxtab(a1[a1_combo2_chemo_horm,],'secondary_cvd', prevcvd[-c(1,5,7,8,9)])
cox_general <- coxtab(a1[a1_allcombo_chemo_horm,],'general_cvd')
cox_death_cvd <- coxtab(a1[a1_allcombo_chemo_horm,],'death_all')
cox6 <- rbind(cox_primary_cvd[1:2,],cox_secondary_cvd[1:2,], cox_general[1:2,],
              cox_death_cvd[1:2,])
cox6[,-1] <- lapply(cox6[,-1], function(x) as.numeric(as.character(x)))
cox6 <- cbind(round2(cox6[,2],2),paste0("(",round2(cox6[,3],2),", ", round2(cox6[,4],2),")"))

# rad_l and horm
cox_primary_cvd <- coxtab(a1[a1_combo_rad_l_horm,],'primary_cvd', prevcvd[-c(2,3,4,6)])
cox_secondary_cvd <- coxtab(a1[a1_combo2_rad_l_horm,],'secondary_cvd', prevcvd[-c(1,5,7,8,9)])
cox_general <- coxtab(a1[a1_allcombo_rad_l_horm,],'general_cvd')
cox_death_cvd <- coxtab(a1[a1_allcombo_rad_l_horm,],'death_all')
cox7 <- rbind(cox_primary_cvd[1:2,],cox_secondary_cvd[1:2,], cox_general[1:2,],
              cox_death_cvd[1:2,])
cox7[,-1] <- lapply(cox7[,-1], function(x) as.numeric(as.character(x)))
cox7 <- cbind(round2(cox7[,2],2),paste0("(",round2(cox7[,3],2),", ", round2(cox7[,4],2),")"))

# chemo, rad_l, and horm
cox_primary_cvd <- coxtab(a1[a1_combo_chemo_horm_rad_l,],'primary_cvd', prevcvd[-c(2,3,4,6)])
cox_secondary_cvd <- coxtab(a1[a1_combo2_chemo_horm_rad_l,],'secondary_cvd', prevcvd[-c(1,5,7,8,9)])
cox_general <- coxtab(a1[a1_allcombo_chemo_horm_rad_l,],'general_cvd')
cox_death_cvd <- coxtab(a1[a1_allcombo_chemo_horm_rad_l,],'death_all')
cox8 <- rbind(cox_primary_cvd[1:2,],cox_secondary_cvd[1:2,], cox_general[1:2,],
              cox_death_cvd[1:2,])
cox8[,-1] <- lapply(cox8[,-1], function(x) as.numeric(as.character(x)))
cox8 <- cbind(round2(cox8[,2],2),paste0("(",round2(cox8[,3],2),", ", round2(cox8[,4],2),")"))



cox_all <- data.frame(cbind(cox1, '', cox2,'', cox3,'', cox4,'',cox5,'',cox6,'',cox7,'',cox8))
write_csv(cox_all,'cox_all_death_sup.csv')


table(a1[a1_combo,]$primary_cvd, a1[a1_combo,]$group1)
table(a1[a1_combo_chemo_only,]$primary_cvd, a1[a1_combo_chemo_only,]$group1)
table(a1[a1_combo_rad_l_only,]$primary_cvd, a1[a1_combo_rad_l_only,]$group1)
table(a1[a1_combo_horm_only,]$primary_cvd, a1[a1_combo_horm_only,]$group1)
table(a1[a1_combo_chemo_rad_l,]$primary_cvd, a1[a1_combo_chemo_rad_l,]$group1)
table(a1[a1_combo_chemo_horm,]$primary_cvd, a1[a1_combo_chemo_horm,]$group1)
table(a1[a1_combo_rad_l_horm,]$primary_cvd, a1[a1_combo_rad_l_horm,]$group1)
table(a1[a1_combo_chemo_horm_rad_l,]$primary_cvd, a1[a1_combo_chemo_horm_rad_l,]$group1)



table(a1[a1_combo2,]$group1)
table(a1[a1_combo2_chemo,]$group1)
table(a1[a1_combo2_horm,]$group1)
table(a1[a1_combo2_rad,]$group1)
table(a1[a1_combo2_rad_l,]$group1)

table(a1[a1_allcombo,]$group1)
table(a1[a1_allcombo_chemo,]$group1)
table(a1[a1_allcombo_horm,]$group1)
table(a1[a1_allcombo_rad,]$group1)
table(a1[a1_allcombo_rad_l,]$group1)


table(a1[a1_combo2,]$secondary_cvd, a1[a1_combo2,]$group1)
table(a1[a1_combo2_chemo,]$secondary_cvd, a1[a1_combo2_chemo,]$group1)
table(a1[a1_combo2_horm,]$secondary_cvd, a1[a1_combo2_horm,]$group1)
table(a1[a1_combo2_rad,]$secondary_cvd, a1[a1_combo2_rad,]$group1)
table(a1[a1_combo2_rad_l,]$secondary_cvd, a1[a1_combo2_rad_l,]$group1)

table(a1[a1_allcombo,]$general_cvd, a1[a1_allcombo,]$group1)
table(a1[a1_allcombo_chemo,]$general_cvd, a1[a1_allcombo_chemo,]$group1)
table(a1[a1_allcombo_horm,]$general_cvd, a1[a1_allcombo_horm,]$group1)
table(a1[a1_allcombo_rad,]$general_cvd, a1[a1_allcombo_rad,]$group1)
table(a1[a1_allcombo_rad_l,]$general_cvd, a1[a1_allcombo_rad_l,]$group1)

table(a1[a1_allcombo,]$death_all, a1[a1_allcombo,]$group1)
table(a1[a1_allcombo_chemo,]$death_all, a1[a1_allcombo_chemo,]$group1)
table(a1[a1_allcombo_horm,]$death_all, a1[a1_allcombo_horm,]$group1)
table(a1[a1_allcombo_rad,]$death_all, a1[a1_allcombo_rad,]$group1)
table(a1[a1_allcombo_rad_l,]$death_all, a1[a1_allcombo_rad_l,]$group1)



################################################################################
# Risk factor analysis
# chemo
cox_htn <- coxtab_htn(a1[a1_htn_chemo_only,],'cvdrf_htn',prevcvd)
cox_diab <- coxtab_diab(a1[a1_diab_chemo_only,],'cvdrf_diab',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_chemo_only,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_chemo_only,],'cvdrfcombo',prevcvd)
cox1 <- rbind(cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round2(cox1[,2],2),paste0("(",round2(cox1[,3],2),", ", round2(cox1[,4],2),")"))

# lest-sided radiation 
cox_htn <- coxtab_htn(a1[a1_htn_rad_l_only,],'cvdrf_htn',prevcvd)
cox_diab <- coxtab_diab(a1[a1_diab_rad_l_only,],'cvdrf_diab',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_rad_l_only,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_rad_l_only,],'cvdrfcombo',prevcvd)
cox2 <- rbind(cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round2(cox2[,2],2),paste0("(",round2(cox2[,3],2),", ", round2(cox2[,4],2),")"))

# hormonal 
cox_htn <- coxtab_htn(a1[a1_htn_horm_only,],'cvdrf_htn',prevcvd)
cox_diab <- coxtab_diab(a1[a1_diab_horm_only,],'cvdrf_diab',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_horm_only,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_horm_only,],'cvdrfcombo',prevcvd)
cox3 <- rbind(cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
cox3 <- cbind(round2(cox3[,2],2),paste0("(",round2(cox3[,3],2),", ", round2(cox3[,4],2),")"))


# chemo and lest-sided radiation
cox_htn <- coxtab_htn(a1[a1_htn_chemo_rad_l,],'cvdrf_htn',prevcvd)
cox_diab <- coxtab_diab(a1[a1_diab_chemo_rad_l,],'cvdrf_diab',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_chemo_rad_l,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_chemo_rad_l,],'cvdrfcombo',prevcvd)
cox4 <- rbind(cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
cox4 <- cbind(round2(cox4[,2],2),paste0("(",round2(cox4[,3],2),", ", round2(cox4[,4],2),")"))


# chemo and horm
cox_htn <- coxtab_htn(a1[a1_htn_chemo_horm,],'cvdrf_htn',prevcvd)
cox_diab <- coxtab_diab(a1[a1_diab_chemo_horm,],'cvdrf_diab',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_chemo_horm,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_chemo_horm,],'cvdrfcombo',prevcvd)
cox5 <- rbind(cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
cox5 <- cbind(round2(cox5[,2],2),paste0("(",round2(cox5[,3],2),", ", round2(cox5[,4],2),")"))

# left-sided radiation and horm
cox_htn <- coxtab_htn(a1[a1_htn_rad_l_horm,],'cvdrf_htn',prevcvd)
cox_diab <- coxtab_diab(a1[a1_diab_rad_l_horm,],'cvdrf_diab',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_rad_l_horm,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_rad_l_horm,],'cvdrfcombo',prevcvd)
cox6 <- rbind(cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,])
cox6[,-1] <- lapply(cox6[,-1], function(x) as.numeric(as.character(x)))
cox6 <- cbind(round2(cox6[,2],2),paste0("(",round2(cox6[,3],2),", ", round2(cox6[,4],2),")"))

# chemo, left-sided radiation and horm
cox_htn <- coxtab_htn(a1[a1_htn_chemo_horm_rad_l,],'cvdrf_htn',prevcvd)
cox_diab <- coxtab_diab(a1[a1_diab_chemo_horm_rad_l,],'cvdrf_diab',prevcvd)
cox_dyslipid <- coxtab_dys(a1[a1_dyslipid_chemo_horm_rad_l,],'cvdrf_dyslipid',prevcvd)
cox_rfcombo <- coxtab(a1[a1_rfcombo_chemo_horm_rad_l,],'cvdrfcombo',prevcvd)
cox7 <- rbind(cox_htn[1:2,],cox_diab[1:2,],cox_dyslipid[1:2,],
              cox_rfcombo[1:2,])
cox7[,-1] <- lapply(cox7[,-1], function(x) as.numeric(as.character(x)))
cox7 <- cbind(round2(cox7[,2],2),paste0("(",round2(cox7[,3],2),", ", round2(cox7[,4],2),")"))


cox_risk <- data.frame(cbind(cox1, '', cox2,'', cox3,'', cox4,'',cox5, '', cox6, '', cox7))
write_csv(cox_risk,'cox_risk.csv')

coxns =  data.frame(t(rbind(c(length(a1_chemo_only),sapply(list(a1_htn_chemo_only,a1_diab_chemo_only,a1_dyslipid_chemo_only,a1_rfcombo_chemo_only),length)),
                            c(length(a1_chemo_only_cases),sapply(list(a1_htn_chemo_only_cases,a1_diab_chemo_only_cases,a1_dyslipid_chemo_only_cases,a1_rfcombo_chemo_only_cases),length)),
                            c(length(a1_chemo_only_controls),sapply(list(a1_htn_chemo_only_controls,a1_diab_chemo_only_controls,a1_dyslipid_chemo_only_controls,a1_rfcombo_chemo_only_controls),length)),
                            c(length(a1_rad_l_only),sapply(list(a1_htn_rad_l_only,a1_diab_rad_l_only,a1_dyslipid_rad_l_only,a1_rfcombo_rad_l_only),length)),
                            c(length(a1_rad_l_only_cases),sapply(list(a1_htn_rad_l_only_cases,a1_diab_rad_l_only_cases,a1_dyslipid_rad_l_only_cases,a1_rfcombo_rad_l_only_cases),length)),
                            c(length(a1_rad_l_only_controls),sapply(list(a1_htn_rad_l_only_controls,a1_diab_rad_l_only_controls,a1_dyslipid_rad_l_only_controls,a1_rfcombo_rad_l_only_controls),length)),
                            c(length(a1_horm_only),sapply(list(a1_htn_horm_only,a1_diab_horm_only,a1_dyslipid_horm_only,a1_rfcombo_horm_only),length)),
                            c(length(a1_horm_only_cases),sapply(list(a1_htn_horm_only_cases,a1_diab_horm_only_cases,a1_dyslipid_horm_only_cases,a1_rfcombo_horm_only_cases),length)),
                            c(length(a1_horm_only_controls),sapply(list(a1_htn_horm_only_controls,a1_diab_horm_only_controls,a1_dyslipid_horm_only_controls,a1_rfcombo_horm_only_controls),length)),
                            c(length(a1_chemo_rad_l),sapply(list(a1_htn_chemo_rad_l,a1_diab_chemo_rad_l,a1_dyslipid_chemo_rad_l,a1_rfcombo_chemo_rad_l),length)),
                            c(length(a1_chemo_rad_l_cases),sapply(list(a1_htn_chemo_rad_l_cases,a1_diab_chemo_rad_l_cases,a1_dyslipid_chemo_rad_l_cases,a1_rfcombo_chemo_rad_l_cases),length)),
                            c(length(a1_chemo_rad_l_controls),sapply(list(a1_htn_chemo_rad_l_controls,a1_diab_chemo_rad_l_controls,a1_dyslipid_chemo_rad_l_controls,a1_rfcombo_chemo_rad_l_controls),length)),
                            c(length(a1_chemo_horm),sapply(list(a1_htn_chemo_horm,a1_diab_chemo_horm,a1_dyslipid_chemo_horm,a1_rfcombo_chemo_horm),length)),
                            c(length(a1_chemo_horm_cases),sapply(list(a1_htn_chemo_horm_cases,a1_diab_chemo_horm_cases,a1_dyslipid_chemo_horm_cases,a1_rfcombo_chemo_horm_cases),length)),
                            c(length(a1_chemo_horm_controls),sapply(list(a1_htn_chemo_horm_controls,a1_diab_chemo_horm_controls,a1_dyslipid_chemo_horm_controls,a1_rfcombo_chemo_horm_controls),length)),
                            c(length(a1_rad_l_horm),sapply(list(a1_htn_rad_l_horm,a1_diab_rad_l_horm,a1_dyslipid_rad_l_horm,a1_rfcombo_rad_l_horm),length)),
                            c(length(a1_rad_l_horm_cases),sapply(list(a1_htn_rad_l_horm_cases,a1_diab_rad_l_horm_cases,a1_dyslipid_rad_l_horm_cases,a1_rfcombo_rad_l_horm_cases),length)),
                            c(length(a1_rad_l_horm_controls),sapply(list(a1_htn_rad_l_horm_controls,a1_diab_rad_l_horm_controls,a1_dyslipid_rad_l_horm_controls,a1_rfcombo_rad_l_horm_controls),length)),
                            c(length(a1_chemo_horm_rad_l),sapply(list(a1_htn_chemo_horm_rad_l,a1_diab_chemo_horm_rad_l,a1_dyslipid_chemo_horm_rad_l,a1_rfcombo_chemo_horm_rad_l),length)),
                            c(length(a1_chemo_horm_rad_l_cases),sapply(list(a1_htn_chemo_horm_rad_l_cases,a1_diab_chemo_horm_rad_l_cases,a1_dyslipid_chemo_horm_rad_l_cases,a1_rfcombo_chemo_horm_rad_l_cases),length)),
                            c(length(a1_chemo_horm_rad_l_controls),sapply(list(a1_htn_chemo_horm_rad_l_controls,a1_diab_chemo_horm_rad_l_controls,a1_dyslipid_chemo_horm_rad_l_controls,a1_rfcombo_chemo_horm_rad_l_controls),length)))))

coxns[,1:21] <- lapply(coxns[,1:21],function(x) paste0('n = ',x))
names(coxns) <- c('chemo','cases','controls','left rad','cases','controls','horm','cases','controls','chemo and left rad','cases','controls',
                  'chemo and horm','cases','controls', 'left rad and horm','cases','controls', 'all therapy','cases','controls')
coxns$sample <- c('all','htn','diabetes', 'dyslipid','cvdrf combo')

write_csv(coxns, 'sample_size_cox_risk.csv')

################################################################################
# Sample size for each analysis

tab1ns <- data.frame(rbind(c(nrow(a1),unlist(table(a1$group))),
                c(nrow(a1[a1_chemo,]),unlist(table(a1[a1_chemo,]$group))),
                c(nrow(a1[a1_horm,]),unlist(table(a1[a1_horm,]$group))),
                c(nrow(a1[a1_rad,]),unlist(table(a1[a1_rad,]$group))),
                c(nrow(a1[a1_rad_l,]),unlist(table(a1[a1_rad_l,]$group)))))
names(tab1ns) <- c('all','cases','controls')
tab1ns$sample <- c('all','chemo','hormonal','radiation','radiation left side')


coxns <- data.frame(t(rbind(c(nrow(a1),
                            sapply(list(a1_ihd,a1_hf,a1_cm,a1_stroke,a1_combo,
                                        a1_htn,a1_diab,a1_dyslipid,a1_rfcombo,
                                        a1_arrhythmia,a1_cardiac,a1_carotid,
                                        a1_myocarditis,a1_tia,a1_valvular,a1_dvt),length)),
                            c(nrow(a1_cases),
                              sapply(list(a1_ihd_cases,a1_hf_cases,a1_cm_cases,a1_stroke_cases,a1_combo_cases,
                                          a1_htn_cases,a1_diab_cases,a1_dyslipid_cases,a1_rfcombo_cases,
                                          a1_arrhythmia_cases,a1_cardiac_cases,a1_carotid_cases,
                                          a1_myocarditis_cases,a1_tia_cases,a1_valvular_cases,a1_dvt_cases),length)),
                            c(nrow(a1_controls),
                              sapply(list(a1_ihd_controls,a1_hf_controls,a1_cm_controls,a1_stroke_controls,a1_combo_controls,
                                          a1_htn_controls,a1_diab_controls,a1_dyslipid_controls,a1_rfcombo_controls,
                                          a1_arrhythmia_controls,a1_cardiac_controls,a1_carotid_controls,
                                          a1_myocarditis_controls,a1_tia_controls,a1_valvular_controls,a1_dvt_controls),length)),
               c(length(a1_chemo),sapply(list(a1_ihd_chemo,a1_hf_chemo,a1_cm_chemo,a1_stroke_chemo,a1_combo_chemo,
                                      a1_htn_chemo,a1_diab_chemo,a1_dyslipid_chemo,a1_rfcombo_chemo,
                                      a1_arrhythmia_chemo,a1_cardiac_chemo,a1_carotid_chemo,
                                      a1_myocarditis_chemo,a1_tia_chemo,a1_valvular_chemo,a1_dvt_chemo),length)),
               c(length(a1_chemo_cases),sapply(list(a1_ihd_chemo_cases,a1_hf_chemo_cases,a1_cm_chemo_cases,a1_stroke_chemo_cases,a1_combo_chemo_cases,
                                              a1_htn_chemo_cases,a1_diab_chemo_cases,a1_dyslipid_chemo_cases,a1_rfcombo_chemo_cases,
                                              a1_arrhythmia_chemo_cases,a1_cardiac_chemo_cases,a1_carotid_chemo_cases,
                                              a1_myocarditis_chemo_cases,a1_tia_chemo_cases,a1_valvular_chemo_cases,a1_dvt_chemo_cases),length)),
               c(length(a1_chemo_controls),sapply(list(a1_ihd_chemo_controls,a1_hf_chemo_controls,a1_cm_chemo_controls,a1_stroke_chemo_controls,a1_combo_chemo_controls,
                                              a1_htn_chemo_controls,a1_diab_chemo_controls,a1_dyslipid_chemo_controls,a1_rfcombo_chemo_controls,
                                              a1_arrhythmia_chemo_controls,a1_cardiac_chemo_controls,a1_carotid_chemo_controls,
                                              a1_myocarditis_chemo_controls,a1_tia_chemo_controls,a1_valvular_chemo_controls,a1_dvt_chemo_controls),length)),
               c(length(a1_horm),sapply(list(a1_ihd_horm,a1_hf_horm,a1_cm_horm,a1_stroke_horm,a1_combo_horm,
                                              a1_htn_horm,a1_diab_horm,a1_dyslipid_horm,a1_rfcombo_horm,
                                              a1_arrhythmia_horm,a1_cardiac_horm,a1_carotid_horm,
                                              a1_myocarditis_horm,a1_tia_horm,a1_valvular_horm,a1_dvt_horm),length)),
               c(length(a1_horm_cases),sapply(list(a1_ihd_horm_cases,a1_hf_horm_cases,a1_cm_horm_cases,a1_stroke_horm_cases,a1_combo_horm_cases,
                                             a1_htn_horm_cases,a1_diab_horm_cases,a1_dyslipid_horm_cases,a1_rfcombo_horm_cases,
                                             a1_arrhythmia_horm_cases,a1_cardiac_horm_cases,a1_carotid_horm_cases,
                                             a1_myocarditis_horm_cases,a1_tia_horm_cases,a1_valvular_horm_cases,a1_dvt_horm_cases),length)),
               c(length(a1_horm_controls),sapply(list(a1_ihd_horm_controls,a1_hf_horm_controls,a1_cm_horm_controls,a1_stroke_horm_controls,a1_combo_horm_controls,
                                             a1_htn_horm_controls,a1_diab_horm_controls,a1_dyslipid_horm_controls,a1_rfcombo_horm_controls,
                                             a1_arrhythmia_horm_controls,a1_cardiac_horm_controls,a1_carotid_horm_controls,
                                             a1_myocarditis_horm_controls,a1_tia_horm_controls,a1_valvular_horm_controls,a1_dvt_horm_controls),length)),
               c(length(a1_rad),sapply(list(a1_ihd_rad,a1_hf_rad,a1_cm_rad,a1_stroke_rad,a1_combo_rad,
                                            a1_htn_rad,a1_diab_rad,a1_dyslipid_rad,a1_rfcombo_rad,
                                            a1_arrhythmia_rad,a1_cardiac_rad,a1_carotid_rad,
                                            a1_myocarditis_rad,a1_tia_rad,a1_valvular_rad,a1_dvt_rad),length)),
               c(length(a1_rad_cases),sapply(list(a1_ihd_rad_cases,a1_hf_rad_cases,a1_cm_rad_cases,a1_stroke_rad_cases,a1_combo_rad_cases,
                                            a1_htn_rad_cases,a1_diab_rad_cases,a1_dyslipid_rad_cases,a1_rfcombo_rad_cases,
                                            a1_arrhythmia_rad_cases,a1_cardiac_rad_cases,a1_carotid_rad_cases,
                                            a1_myocarditis_rad_cases,a1_tia_rad_cases,a1_valvular_rad_cases,a1_dvt_rad_cases),length)),
               c(length(a1_rad_controls),sapply(list(a1_ihd_rad_controls,a1_hf_rad_controls,a1_cm_rad_controls,a1_stroke_rad_controls,a1_combo_rad_controls,
                                            a1_htn_rad_controls,a1_diab_rad_controls,a1_dyslipid_rad_controls,a1_rfcombo_rad_controls,
                                            a1_arrhythmia_rad_controls,a1_cardiac_rad_controls,a1_carotid_rad_controls,
                                            a1_myocarditis_rad_controls,a1_tia_rad_controls,a1_valvular_rad_controls,a1_dvt_rad_controls),length)),
               c(length(a1_rad_l),sapply(list(a1_ihd_rad_l,a1_hf_rad_l,a1_cm_rad_l,a1_stroke_rad_l,a1_combo_rad_l,
                                            a1_htn_rad_l,a1_diab_rad_l,a1_dyslipid_rad_l,a1_rfcombo_rad_l,
                                            a1_arrhythmia_rad_l,a1_cardiac_rad_l,a1_carotid_rad_l,
                                            a1_myocarditis_rad_l,a1_tia_rad_l,a1_valvular_rad_l,a1_dvt_rad_l),length)),
               c(length(a1_rad_l_cases),sapply(list(a1_ihd_rad_l_cases,a1_hf_rad_l_cases,a1_cm_rad_l_cases,a1_stroke_rad_l_cases,a1_combo_rad_l_cases,
                                              a1_htn_rad_l_cases,a1_diab_rad_l_cases,a1_dyslipid_rad_l_cases,a1_rfcombo_rad_l_cases,
                                              a1_arrhythmia_rad_l_cases,a1_cardiac_rad_l_cases,a1_carotid_rad_l_cases,
                                              a1_myocarditis_rad_l_cases,a1_tia_rad_l_cases,a1_valvular_rad_l_cases,a1_dvt_rad_l_cases),length)),
               c(length(a1_rad_l_controls),sapply(list(a1_ihd_rad_l_controls,a1_hf_rad_l_controls,a1_cm_rad_l_controls,a1_stroke_rad_l_controls,a1_combo_rad_l_controls,
                                              a1_htn_rad_l_controls,a1_diab_rad_l_controls,a1_dyslipid_rad_l_controls,a1_rfcombo_rad_l_controls,
                                              a1_arrhythmia_rad_l_controls,a1_cardiac_rad_l_controls,a1_carotid_rad_l_controls,
                                              a1_myocarditis_rad_l_controls,a1_tia_rad_l_controls,a1_valvular_rad_l_controls,a1_dvt_rad_l_controls),length)))))
coxns[,1:15] <- lapply(coxns[,1:15],function(x) paste0('n=',x))
names(coxns) <- c('all','cases','controls','chemo','cases','controls','horm','cases','controls','rad','cases','controls','rad left side','cases','controls')
coxns$sample <- c('all','IHD','HF','CM','stroke','cvd combo','htn','diabetes', 'dyslipid',
                  'cvdrf combo','arrhyth','cardiac','carotid','myocard','tia','valvular','dvt')

# number of events
events <- function(index, outcomevar, timevar){
  surv_dat <- a1[index,c(outcomevar,timevar,'group')]
  surv_dat$time <- round(surv_dat[,2]/30)
  ncases = sum(surv_dat[(surv_dat$time <= 120 & surv_dat$group == "Case"),1])
  ncontrols =  sum(surv_dat[(surv_dat$time <= 120 & surv_dat$group == "Control"),1])
  return(list(n1 = ncases, n2 = ncontrols))
}

riskfactor = c('IHD','HF','CM','stroke','cvd combo','htn','diabetes', 'dyslipid','cvdrf combo')




# all
events_isch <- events(a1_ihd,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
events_stroke <- events(a1_stroke,'stroke_grp_inc','stroke_grp_inc_fu')
events_hf <- events(a1_hf,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
events_cm <- events(a1_cm,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu')
events_cvdcombo <- events(a1_combo,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu')
events_diab <- events(a1_diab,'cvdrf_diab','cvdrf_diab_fu')
events_htn <- events(a1_htn,'cvdrf_htn','cvdrf_htn_fu')
events_dyslipid <- events(a1_dyslipid,'cvdrf_dyslipid','cvdrf_dyslipid_fu')
events_rfcombo <- events(a1_rfcombo,'cvdrfcombo','cvdrfcombo_fu')
case_events = c(events_isch$n1,events_hf$n1,events_cm$n1,events_stroke$n1,events_cvdcombo$n1,
                events_htn$n1,events_diab$n1,events_dyslipid$n1,events_rfcombo$n1)
control_events = c(events_isch$n2,events_hf$n2,events_cm$n2,events_stroke$n2,events_cvdcombo$n2,
                events_htn$n2,events_diab$n2,events_dyslipid$n2,events_rfcombo$n2)

# chemo
events_isch <- events(a1_ihd_chemo,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
events_stroke <- events(a1_stroke_chemo,'stroke_grp_inc','stroke_grp_inc_fu')
events_hf <- events(a1_hf_chemo,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
events_cm <- events(a1_cm_chemo,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu')
events_cvdcombo <- events(a1_combo_chemo,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu')
events_diab <- events(a1_diab_chemo,'cvdrf_diab','cvdrf_diab_fu')
events_htn <- events(a1_htn_chemo,'cvdrf_htn','cvdrf_htn_fu')
events_dyslipid <- events(a1_dyslipid_chemo,'cvdrf_dyslipid','cvdrf_dyslipid_fu')
events_rfcombo <- events(a1_rfcombo_chemo,'cvdrfcombo','cvdrfcombo_fu')
case_events_chemo = c(events_isch$n1,events_hf$n1,events_cm$n1,events_stroke$n1,events_cvdcombo$n1,
                events_htn$n1,events_diab$n1,events_dyslipid$n1,events_rfcombo$n1)
control_events_chemo = c(events_isch$n2,events_hf$n2,events_cm$n2,events_stroke$n2,events_cvdcombo$n2,
                   events_htn$n2,events_diab$n2,events_dyslipid$n2,events_rfcombo$n2)

# horm
events_isch <- events(a1_ihd_horm,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
events_stroke <- events(a1_stroke_horm,'stroke_grp_inc','stroke_grp_inc_fu')
events_hf <- events(a1_hf_horm,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
events_cm <- events(a1_cm_horm,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu')
events_cvdcombo <- events(a1_combo_horm,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu')
events_diab <- events(a1_diab_horm,'cvdrf_diab','cvdrf_diab_fu')
events_htn <- events(a1_htn_horm,'cvdrf_htn','cvdrf_htn_fu')
events_dyslipid <- events(a1_dyslipid_horm,'cvdrf_dyslipid','cvdrf_dyslipid_fu')
events_rfcombo <- events(a1_rfcombo_horm,'cvdrfcombo','cvdrfcombo_fu')
case_events_horm = c(events_isch$n1,events_hf$n1,events_cm$n1,events_stroke$n1,events_cvdcombo$n1,
                      events_htn$n1,events_diab$n1,events_dyslipid$n1,events_rfcombo$n1)
control_events_horm = c(events_isch$n2,events_hf$n2,events_cm$n2,events_stroke$n2,events_cvdcombo$n2,
                         events_htn$n2,events_diab$n2,events_dyslipid$n2,events_rfcombo$n2)

# RT
events_isch <- events(a1_ihd_rad,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
events_stroke <- events(a1_stroke_rad,'stroke_grp_inc','stroke_grp_inc_fu')
events_hf <- events(a1_hf_rad,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
events_cm <- events(a1_cm_rad,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu')
events_cvdcombo <- events(a1_combo_rad,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu')
events_diab <- events(a1_diab_rad,'cvdrf_diab','cvdrf_diab_fu')
events_htn <- events(a1_htn_rad,'cvdrf_htn','cvdrf_htn_fu')
events_dyslipid <- events(a1_dyslipid_rad,'cvdrf_dyslipid','cvdrf_dyslipid_fu')
events_rfcombo <- events(a1_rfcombo_rad,'cvdrfcombo','cvdrfcombo_fu')
case_events_rad = c(events_isch$n1,events_hf$n1,events_cm$n1,events_stroke$n1,events_cvdcombo$n1,
                     events_htn$n1,events_diab$n1,events_dyslipid$n1,events_rfcombo$n1)
control_events_rad = c(events_isch$n2,events_hf$n2,events_cm$n2,events_stroke$n2,events_cvdcombo$n2,
                        events_htn$n2,events_diab$n2,events_dyslipid$n2,events_rfcombo$n2)


# left side RT
events_isch <- events(a1_ihd_rad_l,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
events_stroke <- events(a1_stroke_rad_l,'stroke_grp_inc','stroke_grp_inc_fu')
events_hf <- events(a1_hf_rad_l,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
events_cm <- events(a1_cm_rad_l,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu')
events_cvdcombo <- events(a1_combo_rad_l,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu')
events_diab <- events(a1_diab_rad_l,'cvdrf_diab','cvdrf_diab_fu')
events_htn <- events(a1_htn_rad_l,'cvdrf_htn','cvdrf_htn_fu')
events_dyslipid <- events(a1_dyslipid_rad_l,'cvdrf_dyslipid','cvdrf_dyslipid_fu')
events_rfcombo <- events(a1_rfcombo_rad_l,'cvdrfcombo','cvdrfcombo_fu')
case_events_rad_l = c(events_isch$n1,events_hf$n1,events_cm$n1,events_stroke$n1,events_cvdcombo$n1,
                    events_htn$n1,events_diab$n1,events_dyslipid$n1,events_rfcombo$n1)
control_events_rad_l = c(events_isch$n2,events_hf$n2,events_cm$n2,events_stroke$n2,events_cvdcombo$n2,
                       events_htn$n2,events_diab$n2,events_dyslipid$n2,events_rfcombo$n2)

eventstab = data.frame(riskfactor,case_events,control_events,case_events_chemo,control_events_chemo,
           case_events_horm,control_events_horm,case_events_rad,control_events_rad,
           case_events_rad_l,control_events_rad_l)


# export sample sizes
write_csv(tab1ns, 'sample_size_tab1.csv')
write_csv(coxns, 'sample_size_cox_models.csv')
write_csv(eventstab, 'sample_size_events.csv')


## Cases events for secondary outcomes
# all
events_arrhythmia = events(a1_arrhythmia,'arrhythmia_grp_inc','arrhythmia_grp_inc_fu')
events_cardiac = events(a1_cardiac,'cardiac_arrest_grp_inc','cardiac_arrest_grp_inc_fu')
events_carotid = events(a1_carotid,'carotid_disease_grp_inc','carotid_disease_grp_inc_fu')
events_myocarditis = events(a1_myocarditis,'myocarditis_pericarditis_grp_inc','myocarditis_pericarditis_grp_inc_fu')
events_tia = events(a1_tia,'tia_grp_inc','tia_grp_inc_fu')
events_valvular = events(a1_valvular,'valvular_disease_grp_inc','valvular_disease_grp_inc_fu')
events_venous = events(a1_dvt,'venous_thromboembolic_disease_grp_inc','venous_thromboembolic_disease_grp_inc_fu')

case_events = c(events_arrhythmia$n1,events_cardiac$n1,events_carotid$n1,events_myocarditis$n1,events_tia$n1,
                events_valvular$n1,events_venous$n1)

# chemo
events_arrhythmia = events(a1_arrhythmia_chemo,'arrhythmia_grp_inc','arrhythmia_grp_inc_fu')
events_cardiac = events(a1_cardiac_chemo,'cardiac_arrest_grp_inc','cardiac_arrest_grp_inc_fu')
events_carotid = events(a1_carotid_chemo,'carotid_disease_grp_inc','carotid_disease_grp_inc_fu')
events_myocarditis = events(a1_myocarditis_chemo,'myocarditis_pericarditis_grp_inc','myocarditis_pericarditis_grp_inc_fu')
events_tia = events(a1_tia_chemo,'tia_grp_inc','tia_grp_inc_fu')
events_valvular = events(a1_valvular_chemo,'valvular_disease_grp_inc','valvular_disease_grp_inc_fu')
events_venous = events(a1_dvt_chemo,'venous_thromboembolic_disease_grp_inc','venous_thromboembolic_disease_grp_inc_fu')
case_events_chemo = c(events_arrhythmia$n1,events_cardiac$n1,events_carotid$n1,events_myocarditis$n1,events_tia$n1,
                events_valvular$n1,events_venous$n1)

# left side RT
events_arrhythmia = events(a1_arrhythmia_rad_l,'arrhythmia_grp_inc','arrhythmia_grp_inc_fu')
events_cardiac = events(a1_cardiac_rad_l,'cardiac_arrest_grp_inc','cardiac_arrest_grp_inc_fu')
events_carotid = events(a1_carotid_rad_l,'carotid_disease_grp_inc','carotid_disease_grp_inc_fu')
events_myocarditis = events(a1_myocarditis_rad_l,'myocarditis_pericarditis_grp_inc','myocarditis_pericarditis_grp_inc_fu')
events_tia = events(a1_tia_rad_l,'tia_grp_inc','tia_grp_inc_fu')
events_valvular = events(a1_valvular_rad_l,'valvular_disease_grp_inc','valvular_disease_grp_inc_fu')
events_venous = events(a1_dvt_rad_l,'venous_thromboembolic_disease_grp_inc','venous_thromboembolic_disease_grp_inc_fu')
case_events_rad_l = c(events_arrhythmia$n1,events_cardiac$n1,events_carotid$n1,events_myocarditis$n1,events_tia$n1,
                      events_valvular$n1,events_venous$n1)
# horm
events_arrhythmia = events(a1_arrhythmia_horm,'arrhythmia_grp_inc','arrhythmia_grp_inc_fu')
events_cardiac = events(a1_cardiac_horm,'cardiac_arrest_grp_inc','cardiac_arrest_grp_inc_fu')
events_carotid = events(a1_carotid_horm,'carotid_disease_grp_inc','carotid_disease_grp_inc_fu')
events_myocarditis = events(a1_myocarditis_horm,'myocarditis_pericarditis_grp_inc','myocarditis_pericarditis_grp_inc_fu')
events_tia = events(a1_tia_horm,'tia_grp_inc','tia_grp_inc_fu')
events_valvular = events(a1_valvular_horm,'valvular_disease_grp_inc','valvular_disease_grp_inc_fu')
events_venous = events(a1_dvt_horm,'venous_thromboembolic_disease_grp_inc','venous_thromboembolic_disease_grp_inc_fu')
case_events_horm = c(events_arrhythmia$n1,events_cardiac$n1,events_carotid$n1,events_myocarditis$n1,events_tia$n1,
                      events_valvular$n1,events_venous$n1)


cases_2nd = data.frame(all = case_events, chemo = case_events_chemo, radl = case_events_rad_l, horm = case_events_horm)



####### Events
events <- function(index, outcomevar, timevar){
  surv_dat <- a1[index,c(outcomevar,timevar,'group')]
  ncases = sum(surv_dat[surv_dat$group == "Case",1])
  ncontrols =  sum(surv_dat[surv_dat$group == "Control",1])
  return(list(n1 = ncases, n2 = ncontrols))
}


events_isch <- events(a1_ihd,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
events_stroke <- events(a1_stroke,'stroke_grp_inc','stroke_grp_inc_fu')
events_hf <- events(a1_hf,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
events_cm <- events(a1_cm,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu')
events_cvdcombo <- events(a1_combo,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu')
events_arrhythmia <- events(a1_arrhythmia,'arrhythmia_grp_inc','arrhythmia_grp_inc_fu')
events_cardiac <- events(a1_cardiac,'cardiac_arrest_grp_inc','cardiac_arrest_grp_inc_fu')
events_carotid <- events(a1_carotid,'carotid_disease_grp_inc','carotid_disease_grp_inc_fu')
events_myocarditis <- events(a1_myocarditis,'myocarditis_pericarditis_grp_inc','myocarditis_pericarditis_grp_inc_fu')
events_tia <- events(a1_tia,'tia_grp_inc','tia_grp_inc_fu')
events_valvular <- events(a1_valvular,'valvular_disease_grp_inc','valvular_disease_grp_inc_fu')
events_dvt <- events(a1_dvt,'venous_thromboembolic_disease_grp_inc','venous_thromboembolic_disease_grp_inc_fu')
events_allcvd <- events(a1_allcvd,'allcvd_inc','allcvd_inc_fu')
events_primary_cvd <- events(a1_combo,'primary_cvd', 'primary_cvd_fu')
events_secondary_cvd <- events(a1_combo2,'secondary_cvd', 'secondary_cvd_fu')
events_general <- events(a1_allcombo,'general_cvd','general_cvd_fu')
events_death_cvd <- events(a1_allcombo,'death_all','death_all_fu')


case_events = c(events_isch$n1,events_hf$n1,events_cm$n1,events_stroke$n1,events_cvdcombo$n1,
                events_arrhythmia$n1,events_cardiac$n1, events_carotid$n1, events_myocarditis$n1,
                events_tia$n1, events_valvular$n1, events_dvt$n1, events_allcvd$n1,
                events_primary_cvd$n1, events_secondary_cvd$n1, events_general$n1, events_death_cvd$n1)
control_events = c(events_isch$n2,events_hf$n2,events_cm$n2,events_stroke$n2,events_cvdcombo$n2,
                events_arrhythmia$n2,events_cardiac$n2, events_carotid$n2, events_myocarditis$n2,
                events_tia$n2, events_valvular$n2, events_dvt$n2, events_allcvd$n2,
                events_primary_cvd$n2, events_secondary_cvd$n2, events_general$n2, events_death_cvd$n2)


events_isch <- events(a1_ihd_chemo_only,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
events_stroke <- events(a1_stroke_chemo_only,'stroke_grp_inc','stroke_grp_inc_fu')
events_hf <- events(a1_hf_chemo_only,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
events_cm <- events(a1_cm_chemo_only,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu')
events_cvdcombo <- events(a1_combo_chemo_only,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu')
events_arrhythmia <- events(a1_arrhythmia_chemo_only,'arrhythmia_grp_inc','arrhythmia_grp_inc_fu')
events_cardiac <- events(a1_cardiac_chemo_only,'cardiac_arrest_grp_inc','cardiac_arrest_grp_inc_fu')
events_carotid <- events(a1_carotid_chemo_only,'carotid_disease_grp_inc','carotid_disease_grp_inc_fu')
events_myocarditis <- events(a1_myocarditis_chemo_only,'myocarditis_pericarditis_grp_inc','myocarditis_pericarditis_grp_inc_fu')
events_tia <- events(a1_tia_chemo_only,'tia_grp_inc','tia_grp_inc_fu')
events_valvular <- events(a1_valvular_chemo_only,'valvular_disease_grp_inc','valvular_disease_grp_inc_fu')
events_dvt <- events(a1_dvt_chemo_only,'venous_thromboembolic_disease_grp_inc','venous_thromboembolic_disease_grp_inc_fu')
events_allcvd <- events(a1_allcvd_chemo_only,'allcvd_inc','allcvd_inc_fu')
events_primary_cvd <- events(a1_combo_chemo_only,'primary_cvd', 'primary_cvd_fu')
events_secondary_cvd <- events(a1_combo2_chemo_only,'secondary_cvd', 'secondary_cvd_fu')
events_general <- events(a1_allcombo_chemo_only,'general_cvd','general_cvd_fu')
events_death_cvd <- events(a1_allcombo_chemo_only,'death_all','death_all_fu')


case_events_chemo_only = c(events_isch$n1,events_hf$n1,events_cm$n1,events_stroke$n1,events_cvdcombo$n1,
                events_arrhythmia$n1,events_cardiac$n1, events_carotid$n1, events_myocarditis$n1,
                events_tia$n1, events_valvular$n1, events_dvt$n1, events_allcvd$n1,
                events_primary_cvd$n1, events_secondary_cvd$n1, events_general$n1, events_death_cvd$n1)
control_events_chemo_only = c(events_isch$n2,events_hf$n2,events_cm$n2,events_stroke$n2,events_cvdcombo$n2,
                   events_arrhythmia$n2,events_cardiac$n2, events_carotid$n2, events_myocarditis$n2,
                   events_tia$n2, events_valvular$n2, events_dvt$n2, events_allcvd$n2,
                   events_primary_cvd$n2, events_secondary_cvd$n2, events_general$n2, events_death_cvd$n2)

events_isch <- events(a1_ihd_rad_l_only,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
events_stroke <- events(a1_stroke_rad_l_only,'stroke_grp_inc','stroke_grp_inc_fu')
events_hf <- events(a1_hf_rad_l_only,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
events_cm <- events(a1_cm_rad_l_only,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu')
events_cvdcombo <- events(a1_combo_rad_l_only,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu')
events_arrhythmia <- events(a1_arrhythmia_rad_l_only,'arrhythmia_grp_inc','arrhythmia_grp_inc_fu')
events_cardiac <- events(a1_cardiac_rad_l_only,'cardiac_arrest_grp_inc','cardiac_arrest_grp_inc_fu')
events_carotid <- events(a1_carotid_rad_l_only,'carotid_disease_grp_inc','carotid_disease_grp_inc_fu')
events_myocarditis <- events(a1_myocarditis_rad_l_only,'myocarditis_pericarditis_grp_inc','myocarditis_pericarditis_grp_inc_fu')
events_tia <- events(a1_tia_rad_l_only,'tia_grp_inc','tia_grp_inc_fu')
events_valvular <- events(a1_valvular_rad_l_only,'valvular_disease_grp_inc','valvular_disease_grp_inc_fu')
events_dvt <- events(a1_dvt_rad_l_only,'venous_thromboembolic_disease_grp_inc','venous_thromboembolic_disease_grp_inc_fu')
events_allcvd <- events(a1_allcvd_rad_l_only,'allcvd_inc','allcvd_inc_fu')
events_primary_cvd <- events(a1_combo_rad_l_only,'primary_cvd', 'primary_cvd_fu')
events_secondary_cvd <- events(a1_combo2_rad_l_only,'secondary_cvd', 'secondary_cvd_fu')
events_general <- events(a1_allcombo_rad_l_only,'general_cvd','general_cvd_fu')
events_death_cvd <- events(a1_allcombo_rad_l_only,'death_all','death_all_fu')


case_events_rad_l_only = c(events_isch$n1,events_hf$n1,events_cm$n1,events_stroke$n1,events_cvdcombo$n1,
                events_arrhythmia$n1,events_cardiac$n1, events_carotid$n1, events_myocarditis$n1,
                events_tia$n1, events_valvular$n1, events_dvt$n1, events_allcvd$n1,
                events_primary_cvd$n1, events_secondary_cvd$n1, events_general$n1, events_death_cvd$n1)
control_events_rad_l_only = c(events_isch$n2,events_hf$n2,events_cm$n2,events_stroke$n2,events_cvdcombo$n2,
                   events_arrhythmia$n2,events_cardiac$n2, events_carotid$n2, events_myocarditis$n2,
                   events_tia$n2, events_valvular$n2, events_dvt$n2, events_allcvd$n2,
                   events_primary_cvd$n2, events_secondary_cvd$n2, events_general$n2, events_death_cvd$n2)


events_isch <- events(a1_ihd_horm_only,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
events_stroke <- events(a1_stroke_horm_only,'stroke_grp_inc','stroke_grp_inc_fu')
events_hf <- events(a1_hf_horm_only,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
events_cm <- events(a1_cm_horm_only,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu')
events_cvdcombo <- events(a1_combo_horm_only,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu')
events_arrhythmia <- events(a1_arrhythmia_horm_only,'arrhythmia_grp_inc','arrhythmia_grp_inc_fu')
events_cardiac <- events(a1_cardiac_horm_only,'cardiac_arrest_grp_inc','cardiac_arrest_grp_inc_fu')
events_carotid <- events(a1_carotid_horm_only,'carotid_disease_grp_inc','carotid_disease_grp_inc_fu')
events_myocarditis <- events(a1_myocarditis_horm_only,'myocarditis_pericarditis_grp_inc','myocarditis_pericarditis_grp_inc_fu')
events_tia <- events(a1_tia_horm_only,'tia_grp_inc','tia_grp_inc_fu')
events_valvular <- events(a1_valvular_horm_only,'valvular_disease_grp_inc','valvular_disease_grp_inc_fu')
events_dvt <- events(a1_dvt_horm_only,'venous_thromboembolic_disease_grp_inc','venous_thromboembolic_disease_grp_inc_fu')
events_allcvd <- events(a1_allcvd_horm_only,'allcvd_inc','allcvd_inc_fu')
events_primary_cvd <- events(a1_combo_horm_only,'primary_cvd', 'primary_cvd_fu')
events_secondary_cvd <- events(a1_combo2_horm_only,'secondary_cvd', 'secondary_cvd_fu')
events_general <- events(a1_allcombo_horm_only,'general_cvd','general_cvd_fu')
events_death_cvd <- events(a1_allcombo_horm_only,'death_all','death_all_fu')


case_events_horm_only = c(events_isch$n1,events_hf$n1,events_cm$n1,events_stroke$n1,events_cvdcombo$n1,
                events_arrhythmia$n1,events_cardiac$n1, events_carotid$n1, events_myocarditis$n1,
                events_tia$n1, events_valvular$n1, events_dvt$n1, events_allcvd$n1,
                events_primary_cvd$n1, events_secondary_cvd$n1, events_general$n1, events_death_cvd$n1)
control_events_horm_only = c(events_isch$n2,events_hf$n2,events_cm$n2,events_stroke$n2,events_cvdcombo$n2,
                   events_arrhythmia$n2,events_cardiac$n2, events_carotid$n2, events_myocarditis$n2,
                   events_tia$n2, events_valvular$n2, events_dvt$n2, events_allcvd$n2,
                   events_primary_cvd$n2, events_secondary_cvd$n2, events_general$n2, events_death_cvd$n2)


events_isch <- events(a1_ihd_chemo_rad_l,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
events_stroke <- events(a1_stroke_chemo_rad_l,'stroke_grp_inc','stroke_grp_inc_fu')
events_hf <- events(a1_hf_chemo_rad_l,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
events_cm <- events(a1_cm_chemo_rad_l,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu')
events_cvdcombo <- events(a1_combo_chemo_rad_l,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu')
events_arrhythmia <- events(a1_arrhythmia_chemo_rad_l,'arrhythmia_grp_inc','arrhythmia_grp_inc_fu')
events_cardiac <- events(a1_cardiac_chemo_rad_l,'cardiac_arrest_grp_inc','cardiac_arrest_grp_inc_fu')
events_carotid <- events(a1_carotid_chemo_rad_l,'carotid_disease_grp_inc','carotid_disease_grp_inc_fu')
events_myocarditis <- events(a1_myocarditis_chemo_rad_l,'myocarditis_pericarditis_grp_inc','myocarditis_pericarditis_grp_inc_fu')
events_tia <- events(a1_tia_chemo_rad_l,'tia_grp_inc','tia_grp_inc_fu')
events_valvular <- events(a1_valvular_chemo_rad_l,'valvular_disease_grp_inc','valvular_disease_grp_inc_fu')
events_dvt <- events(a1_dvt_chemo_rad_l,'venous_thromboembolic_disease_grp_inc','venous_thromboembolic_disease_grp_inc_fu')
events_allcvd <- events(a1_allcvd_chemo_rad_l,'allcvd_inc','allcvd_inc_fu')
events_primary_cvd <- events(a1_combo_chemo_rad_l,'primary_cvd', 'primary_cvd_fu')
events_secondary_cvd <- events(a1_combo2_chemo_rad_l,'secondary_cvd', 'secondary_cvd_fu')
events_general <- events(a1_allcombo_chemo_rad_l,'general_cvd','general_cvd_fu')
events_death_cvd <- events(a1_allcombo_chemo_rad_l,'death_all','death_all_fu')


case_events_chemo_rad_l = c(events_isch$n1,events_hf$n1,events_cm$n1,events_stroke$n1,events_cvdcombo$n1,
                events_arrhythmia$n1,events_cardiac$n1, events_carotid$n1, events_myocarditis$n1,
                events_tia$n1, events_valvular$n1, events_dvt$n1, events_allcvd$n1,
                events_primary_cvd$n1, events_secondary_cvd$n1, events_general$n1, events_death_cvd$n1)
control_events_chemo_rad_l = c(events_isch$n2,events_hf$n2,events_cm$n2,events_stroke$n2,events_cvdcombo$n2,
                   events_arrhythmia$n2,events_cardiac$n2, events_carotid$n2, events_myocarditis$n2,
                   events_tia$n2, events_valvular$n2, events_dvt$n2, events_allcvd$n2,
                   events_primary_cvd$n2, events_secondary_cvd$n2, events_general$n2, events_death_cvd$n2)


events_isch <- events(a1_ihd_chemo_horm,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
events_stroke <- events(a1_stroke_chemo_horm,'stroke_grp_inc','stroke_grp_inc_fu')
events_hf <- events(a1_hf_chemo_horm,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
events_cm <- events(a1_cm_chemo_horm,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu')
events_cvdcombo <- events(a1_combo_chemo_horm,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu')
events_arrhythmia <- events(a1_arrhythmia_chemo_horm,'arrhythmia_grp_inc','arrhythmia_grp_inc_fu')
events_cardiac <- events(a1_cardiac_chemo_horm,'cardiac_arrest_grp_inc','cardiac_arrest_grp_inc_fu')
events_carotid <- events(a1_carotid_chemo_horm,'carotid_disease_grp_inc','carotid_disease_grp_inc_fu')
events_myocarditis <- events(a1_myocarditis_chemo_horm,'myocarditis_pericarditis_grp_inc','myocarditis_pericarditis_grp_inc_fu')
events_tia <- events(a1_tia_chemo_horm,'tia_grp_inc','tia_grp_inc_fu')
events_valvular <- events(a1_valvular_chemo_horm,'valvular_disease_grp_inc','valvular_disease_grp_inc_fu')
events_dvt <- events(a1_dvt_chemo_horm,'venous_thromboembolic_disease_grp_inc','venous_thromboembolic_disease_grp_inc_fu')
events_allcvd <- events(a1_allcvd_chemo_horm,'allcvd_inc','allcvd_inc_fu')
events_primary_cvd <- events(a1_combo_chemo_horm,'primary_cvd', 'primary_cvd_fu')
events_secondary_cvd <- events(a1_combo2_chemo_horm,'secondary_cvd', 'secondary_cvd_fu')
events_general <- events(a1_allcombo_chemo_horm,'general_cvd','general_cvd_fu')
events_death_cvd <- events(a1_allcombo_chemo_horm,'death_all','death_all_fu')


case_events_chemo_horm = c(events_isch$n1,events_hf$n1,events_cm$n1,events_stroke$n1,events_cvdcombo$n1,
                events_arrhythmia$n1,events_cardiac$n1, events_carotid$n1, events_myocarditis$n1,
                events_tia$n1, events_valvular$n1, events_dvt$n1, events_allcvd$n1,
                events_primary_cvd$n1, events_secondary_cvd$n1, events_general$n1, events_death_cvd$n1)
control_events_chemo_horm = c(events_isch$n2,events_hf$n2,events_cm$n2,events_stroke$n2,events_cvdcombo$n2,
                   events_arrhythmia$n2,events_cardiac$n2, events_carotid$n2, events_myocarditis$n2,
                   events_tia$n2, events_valvular$n2, events_dvt$n2, events_allcvd$n2,
                   events_primary_cvd$n2, events_secondary_cvd$n2, events_general$n2, events_death_cvd$n2)



events_isch <- events(a1_ihd_rad_l_horm,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
events_stroke <- events(a1_stroke_rad_l_horm,'stroke_grp_inc','stroke_grp_inc_fu')
events_hf <- events(a1_hf_rad_l_horm,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
events_cm <- events(a1_cm_rad_l_horm,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu')
events_cvdcombo <- events(a1_combo_rad_l_horm,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu')
events_arrhythmia <- events(a1_arrhythmia_rad_l_horm,'arrhythmia_grp_inc','arrhythmia_grp_inc_fu')
events_cardiac <- events(a1_cardiac_rad_l_horm,'cardiac_arrest_grp_inc','cardiac_arrest_grp_inc_fu')
events_carotid <- events(a1_carotid_rad_l_horm,'carotid_disease_grp_inc','carotid_disease_grp_inc_fu')
events_myocarditis <- events(a1_myocarditis_rad_l_horm,'myocarditis_pericarditis_grp_inc','myocarditis_pericarditis_grp_inc_fu')
events_tia <- events(a1_tia_rad_l_horm,'tia_grp_inc','tia_grp_inc_fu')
events_valvular <- events(a1_valvular_rad_l_horm,'valvular_disease_grp_inc','valvular_disease_grp_inc_fu')
events_dvt <- events(a1_dvt_rad_l_horm,'venous_thromboembolic_disease_grp_inc','venous_thromboembolic_disease_grp_inc_fu')
events_allcvd <- events(a1_allcvd_rad_l_horm,'allcvd_inc','allcvd_inc_fu')
events_primary_cvd <- events(a1_combo_rad_l_horm,'primary_cvd', 'primary_cvd_fu')
events_secondary_cvd <- events(a1_combo2_rad_l_horm,'secondary_cvd', 'secondary_cvd_fu')
events_general <- events(a1_allcombo_rad_l_horm,'general_cvd','general_cvd_fu')
events_death_cvd <- events(a1_allcombo_rad_l_horm,'death_all','death_all_fu')


case_events_rad_l_horm = c(events_isch$n1,events_hf$n1,events_cm$n1,events_stroke$n1,events_cvdcombo$n1,
                events_arrhythmia$n1,events_cardiac$n1, events_carotid$n1, events_myocarditis$n1,
                events_tia$n1, events_valvular$n1, events_dvt$n1, events_allcvd$n1,
                events_primary_cvd$n1, events_secondary_cvd$n1, events_general$n1, events_death_cvd$n1)
control_events_rad_l_horm = c(events_isch$n2,events_hf$n2,events_cm$n2,events_stroke$n2,events_cvdcombo$n2,
                   events_arrhythmia$n2,events_cardiac$n2, events_carotid$n2, events_myocarditis$n2,
                   events_tia$n2, events_valvular$n2, events_dvt$n2, events_allcvd$n2,
                   events_primary_cvd$n2, events_secondary_cvd$n2, events_general$n2, events_death_cvd$n2)


events_isch <- events(a1_ihd_chemo_horm_rad_l,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
events_stroke <- events(a1_stroke_chemo_horm_rad_l,'stroke_grp_inc','stroke_grp_inc_fu')
events_hf <- events(a1_hf_chemo_horm_rad_l,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
events_cm <- events(a1_cm_chemo_horm_rad_l,'cardiomyopathy_grp_inc','cardiomyopathy_grp_inc_fu')
events_cvdcombo <- events(a1_combo_chemo_horm_rad_l,'cvdcombo_grp_inc','cvdcombo_grp_inc_fu')
events_arrhythmia <- events(a1_arrhythmia_chemo_horm_rad_l,'arrhythmia_grp_inc','arrhythmia_grp_inc_fu')
events_cardiac <- events(a1_cardiac_chemo_horm_rad_l,'cardiac_arrest_grp_inc','cardiac_arrest_grp_inc_fu')
events_carotid <- events(a1_carotid_chemo_horm_rad_l,'carotid_disease_grp_inc','carotid_disease_grp_inc_fu')
events_myocarditis <- events(a1_myocarditis_chemo_horm_rad_l,'myocarditis_pericarditis_grp_inc','myocarditis_pericarditis_grp_inc_fu')
events_tia <- events(a1_tia_chemo_horm_rad_l,'tia_grp_inc','tia_grp_inc_fu')
events_valvular <- events(a1_valvular_chemo_horm_rad_l,'valvular_disease_grp_inc','valvular_disease_grp_inc_fu')
events_dvt <- events(a1_dvt_chemo_horm_rad_l,'venous_thromboembolic_disease_grp_inc','venous_thromboembolic_disease_grp_inc_fu')
events_allcvd <- events(a1_allcvd_chemo_horm_rad_l,'allcvd_inc','allcvd_inc_fu')
events_primary_cvd <- events(a1_combo_chemo_horm_rad_l,'primary_cvd', 'primary_cvd_fu')
events_secondary_cvd <- events(a1_combo2_chemo_horm_rad_l,'secondary_cvd', 'secondary_cvd_fu')
events_general <- events(a1_allcombo_chemo_horm_rad_l,'general_cvd','general_cvd_fu')
events_death_cvd <- events(a1_allcombo_chemo_horm_rad_l,'death_all','death_all_fu')


case_events_chemo_horm_rad_l = c(events_isch$n1,events_hf$n1,events_cm$n1,events_stroke$n1,events_cvdcombo$n1,
                events_arrhythmia$n1,events_cardiac$n1, events_carotid$n1, events_myocarditis$n1,
                events_tia$n1, events_valvular$n1, events_dvt$n1, events_allcvd$n1,
                events_primary_cvd$n1, events_secondary_cvd$n1, events_general$n1, events_death_cvd$n1)
control_events_chemo_horm_rad_l = c(events_isch$n2,events_hf$n2,events_cm$n2,events_stroke$n2,events_cvdcombo$n2,
                   events_arrhythmia$n2,events_cardiac$n2, events_carotid$n2, events_myocarditis$n2,
                   events_tia$n2, events_valvular$n2, events_dvt$n2, events_allcvd$n2,
                   events_primary_cvd$n2, events_secondary_cvd$n2, events_general$n2, events_death_cvd$n2)

events_cases = cbind(case_events, case_events_chemo_only, case_events_rad_l_only, case_events_horm_only, case_events_chemo_rad_l, case_events_chemo_horm, case_events_rad_l_horm, 
                     case_events_chemo_horm_rad_l)

events_controls = cbind(control_events, control_events_chemo_only, control_events_rad_l_only, control_events_horm_only, control_events_chemo_rad_l, control_events_chemo_horm, control_events_rad_l_horm, 
                     control_events_chemo_horm_rad_l)

write.csv(as.data.frame(events_cases),"events_cases.csv")
write.csv(as.data.frame(events_controls),"events_controls.csv")

########################### CIF and Forest Plots ###############################

### Forest Plots
require(forestplot)
# Primary outcomes
tabletext <- cbind(
  c("", "Major CVD events", "Ischemic heart disease", "Entire cohort", 
    "Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy",  "", "Heart failure", "Entire cohort", 
    "Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy", "",
    "Cardiomyopathy", "Entire cohort", 
    "Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy", "",
    "Stroke","Entire cohort", 
    "Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy"),
  c("Events", "(Cases)", "", "308","80","98","174", "", "",
    "543","159","167","276","", "",
    "16","8","6","6","","","623","154","181","358"),
  c("", "HR (95% CI)","","0.91 (0.81, 1.02)","0.84 (0.67, 1.06)","0.97 (0.78, 1.19)","0.93 (0.79, 1.09)","","",
    "1.10 (1.00, 1.20)","1.54 (1.30, 1.83)","1.25 (1.06, 1.47)","1.01 (0.89, 1.15)","","",
    "1.46 (0.89, 2.41)","1.97 (1.00, 3.91)","1.62 (0.72, 3.63)","1.11 (0.49, 2.54)","","",
    "1.00 (0.92, 1.09)","0.97 (0.82, 1.15)","0.96 (0.83, 1.13)","1.02 (0.91, 1.14)"))

temp.data = data.frame(ihd = c(0.91,0.84,0.97,0.93,0.81,0.67,0.78,0.79,1.02,1.06,1.19,1.09),
                       hf = c(1.10,1.54,1.25,1.01,1.00,1.30,1.06,0.89,1.20,1.83,1.47,1.15),
                       card = c(1.46,1.97,1.62,1.11,0.89,1.00,0.72,0.49,2.41,3.91,3.63,2.54),
                       stroke = c(1.00,0.97,0.96,1.02,0.92,0.82,0.83,0.91,1.09,1.15,1.13,1.14))

forestplot(tabletext, 
           txt_gp = fpTxtGp(ticks=gpar(cex=1), xlab=gpar(cex=1)),
           boxsize = 0.5,
           title="Relative risk of developing major CVD events",
           is.summary=c(TRUE, TRUE, TRUE, rep(FALSE,5), TRUE, rep(FALSE,5), TRUE, rep(FALSE,5), TRUE, rep(FALSE,5)),
           mean = c(NA, NA, NA, temp.data$ihd[1:4], NA, NA, temp.data$hf[1:4],NA, NA, temp.data$card[1:4], NA, NA, temp.data$stroke[1:4]),
           lower = c(NA, NA, NA, temp.data$ihd[5:8], NA, NA, temp.data$hf[5:8],NA, NA, temp.data$card[5:8], NA, NA, temp.data$stroke[5:8]),
           upper = c(NA, NA, NA, temp.data$ihd[9:12], NA, NA, temp.data$hf[9:12],NA, NA, temp.data$card[9:12], NA, NA, temp.data$stroke[9:12]),
           clip =c(0, 5),
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           zero = 1,
           xlab="Hazard Ratio (95% CI)")





# Secondary outcomes
tabletext <- cbind(
  c("", "Other CVD events", "Arrhythmia", "Entire cohort", 
    "Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy","", "Cardiac arrest", "Entire cohort", 
    "Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy","",
    "Carotid disease", "Entire cohort", 
    "Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy","",
    "Myocarditis pericarditis","Entire cohort", 
    "Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy","",
    "Transient ischemic attack","Entire cohort","Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy","",
    "Valvular disease","Entire cohort","Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy","",
    "Venous thromboembolism","Entire cohort","Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy"),
  c("Events", "(Cases)", "", 703,195,201,386, "", "",
    85,25,25,47,"", "",3,0,1,2,"","",7,4,3,3,"","",
    35,15,12,14,"","",103,24,31,54,"","",370,202,121,176),
  c("", "HR (95% CI)","","1.04 (0.97, 1.13)","1.15 (0.99, 1.34)","0.98 (0.84, 1.13)","0.99 (0.89, 1.10)","","",
    "1.31 (1.03, 1.66)","1.33 (0.85, 2.07)","1.29 (0.83, 2.01)","1.40 (1.01, 1.95)","","",
    "0.97 (0.22, 2.85)","0.00 (0.00, 0.00)","1.68 (0.34, 8.34)","1.79 (0.46, 6.94)","","",
    "0.81 (0.38, 1.73)","1.36 (0.50, 3.69)","1.05 (0.35, 3.14)","0.74 (0.22, 2,51)","","",
    "0.76 (0.54, 1.09)","1.25 (0.71, 2.21)","0.94 (0.52, 1.70)","0.53 (0.30, 0.92)","","",
    "1.00 (0.82, 1.23)","0.99 (0.65, 1.49)","1.02 (0.70, 1.48)","0.88 (0.66, 1.17)","","",
    "1.99 (1.77, 2.24)","3.34 (2.80, 3.98)","1.96 (1.59, 2.41)","1.71 (1.43, 2.03)"))

temp.data = data.frame(arrh = c(1.04,1.15,0.98,0.99,0.97,0.99,0.84,0.89,1.13,1.34,1.13,1.10),
                       card = c(1.31,1.33,1.29,1.40,1.03,0.85,0.83,1.01,1.66,2.07,2.01,1.95),
                       caro = c(0.97,0.00,1.68,1.79,0.22,0.00,0.34,0.46,2.85,0.00,8.34,6.94),
                       myo = c(0.81,1.36,1.05,0.74,0.38,0.50,0.35,0.22,1.73,3.69,3.14,2.51),
                       tia = c(0.76,1.25,0.94,0.53,0.54,0.71,0.52,0.30,1.09,2.21,1.70,0.92),
                       val = c(1.00,0.99,1.02,0.88,0.82,0.65,0.70,0.66,1.23,1.49,1.48,1.17),
                       ven = c(1.99,3.34,1.96,1.71,1.77,2.80,1.59,1.43,2.24,3.98,2.41,2.03))

forestplot(tabletext, 
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=0.85),ticks=gpar(cex=1), xlab=gpar(cex=1)),
           boxsize = 0.5,
           title="Relative risk of developing other CVD events",
           is.summary=c(TRUE, TRUE, TRUE, rep(FALSE,5), TRUE, rep(FALSE,5), TRUE, rep(FALSE,5), TRUE, rep(FALSE,5),TRUE, rep(FALSE,5),TRUE, rep(FALSE,5),TRUE, rep(FALSE,5)),
           mean = c(NA, NA, NA, temp.data$arrh[1:4], NA, NA, temp.data$card[1:4],NA, NA, temp.data$caro[1:4], NA, NA, temp.data$myo[1:4],NA, NA, temp.data$tia[1:4],NA, NA, temp.data$val[1:4],NA, NA, temp.data$ven[1:4]),
           lower = c(NA, NA, NA, temp.data$arrh[5:8], NA, NA, temp.data$card[5:8],NA, NA, temp.data$caro[5:8], NA, NA, temp.data$myo[5:8],NA, NA, temp.data$tia[5:8],NA, NA, temp.data$val[5:8],NA, NA, temp.data$ven[5:8]),
           upper = c(NA, NA, NA, temp.data$arrh[9:12], NA, NA, temp.data$card[9:12],NA, NA, temp.data$caro[9:12], NA, NA, temp.data$myo[9:12], NA, NA, temp.data$tia[9:12], NA, NA, temp.data$val[9:12], NA, NA, temp.data$ven[9:12]),
           clip =c(0, 5),
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           zero = 1,
           xlab="Hazard Ratio (95% CI)")


# Combined primary and secondary forest plot
require(forestplot)
# Primary outcomes
tabletext <- cbind(
  c("", "Primary CVD events", "Ischemic heart disease", "Entire cohort", 
    "Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy",  "", "Heart failure", "Entire cohort", 
    "Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy", "",
    "Cardiomyopathy", "Entire cohort", 
    "Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy", "",
    "Stroke","Entire cohort", 
    "Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy", "",
    "Secondary CVD events", "Cardiac arrest", "Entire cohort", 
    "Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy", "",
    "Venous thromboembolism","Entire cohort","Received chemotherapy", "Received left-sided radiation therapy", "Received hormonal therapy"),
  c("Events", "(Cases)", "", "308","80","98","174", "", "",
    "543","159","167","276","", "",
    "16","8","6","6","","","623","154","181","358","","","","85","25","25","47","","","370","202","121","176"),
  c("", "HR (95% CI)","","0.91 (0.81, 1.02)","0.84 (0.67, 1.06)","0.97 (0.78, 1.19)","0.93 (0.79, 1.09)","","",
    "1.10 (1.00, 1.20)","1.54 (1.30, 1.83)","1.25 (1.06, 1.47)","1.01 (0.89, 1.15)","","",
    "1.46 (0.89, 2.41)","1.97 (1.00, 3.91)","1.62 (0.72, 3.63)","1.11 (0.49, 2.54)","","",
    "1.00 (0.92, 1.09)","0.97 (0.82, 1.15)","0.96 (0.83, 1.13)","1.02 (0.91, 1.14)","","","",
    "1.31 (1.03, 1.66)","1.33 (0.85, 2.07)","1.29 (0.83, 2.01)","1.40 (1.01, 1.95)","","",
    "1.99 (1.77, 2.24)","3.34 (2.80, 3.98)","1.96 (1.59, 2.41)","1.71 (1.43, 2.03)"))

temp.data = data.frame(ihd = c(0.91,0.84,0.97,0.93,0.81,0.67,0.78,0.79,1.02,1.06,1.19,1.09),
                       hf = c(1.10,1.54,1.25,1.01,1.00,1.30,1.06,0.89,1.20,1.83,1.47,1.15),
                       card = c(1.46,1.97,1.62,1.11,0.89,1.00,0.72,0.49,2.41,3.91,3.63,2.54),
                       stroke = c(1.00,0.97,0.96,1.02,0.92,0.82,0.83,0.91,1.09,1.15,1.13,1.14),
                       cardiac = c(1.31,1.33,1.29,1.40,1.03,0.85,0.83,1.01,1.66,2.07,2.01,1.95),
                       ven = c(1.99,3.34,1.96,1.71,1.77,2.80,1.59,1.43,2.24,3.98,2.41,2.03))

forestplot(tabletext, 
           txt_gp = fpTxtGp(ticks=gpar(cex=1), xlab=gpar(cex=1)),
           boxsize = 0.5,
           title="Relative risk of developing primary and secondary CVD events",
           is.summary=c(TRUE, TRUE, TRUE, rep(FALSE,5), TRUE, rep(FALSE,5), TRUE, rep(FALSE,5), TRUE, rep(FALSE,5),TRUE, TRUE,rep(FALSE,5), TRUE, rep(FALSE,5)),
           mean = c(NA, NA, NA, temp.data$ihd[1:4], NA, NA, temp.data$hf[1:4],NA, NA, temp.data$card[1:4], NA, NA, temp.data$stroke[1:4], NA, NA, NA, temp.data$cardiac[1:4], NA, NA, temp.data$ven[1:4]),
           lower = c(NA, NA, NA, temp.data$ihd[5:8], NA, NA, temp.data$hf[5:8],NA, NA, temp.data$card[5:8], NA, NA, temp.data$stroke[5:8], NA, NA, NA, temp.data$cardiac[5:8], NA, NA, temp.data$ven[5:8]),
           upper = c(NA, NA, NA, temp.data$ihd[9:12], NA, NA, temp.data$hf[9:12],NA, NA, temp.data$card[9:12], NA, NA, temp.data$stroke[9:12], NA, NA, NA, temp.data$cardiac[9:12], NA, NA, temp.data$ven[9:12]),
           clip =c(0, 5),
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           zero = 1,
           xlab="Hazard Ratio (95% CI)")



### CIF curves
### CIFs
splots <- list()
### ischemic heart disease
surv_object <- Surv(time = floor(a1[a1_ihd,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1[a1_ihd,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_ihd,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = F,palette = 'jco',title='Ischemic heart disease (Entire cohort)',xlab='Time/months',
                          ylim = c(0,0.05), xlim = c(0,120),
                          ggtheme = theme_classic2(base_size=8, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_ihd_chemo,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1[a1_ihd_chemo,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd_chemo,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_ihd_chemo,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = F,palette = 'jco',title='Ischemic heart disease (Received chemotherapy)',xlab='Time/months',
                          ylim = c(0,0.04),xlim = c(0,120),
                          ggtheme = theme_classic2(base_size=8, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_ihd_rad_l,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1[a1_ihd_rad_l,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd_rad_l,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_ihd_rad_l,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = F,palette = 'jco',title='Ischemic heart disease (Received left-sided radiation therapy)',xlab='Time/months',
                          ylim = c(0,0.05),xlim = c(0,120),
                          ggtheme = theme_classic2(base_size=8, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_ihd_horm,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = a1[a1_ihd_horm,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_ihd_horm,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_ihd_horm,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = F,palette = 'jco',title='Ischemic heart disease (Received hormonal therapy)',xlab='Time/months',
                          ylim = c(0,0.05),xlim = c(0,120),
                          ggtheme = theme_classic2(base_size=8, base_family = "Arial"))


# heart failure
surv_object <- Surv(time = floor(a1[a1_hf,]$heart_failure_grp_inc_fu/30), 
                    event = a1[a1_hf,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_hf,])
splots[[5]] <- ggsurvplot(fit1, data = a1[a1_hf,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = F,palette = 'jco',title='Heart failure (Entire cohort)',xlab='Time/months',
                          ylim = c(0,0.1),xlim = c(0,120),
                          ggtheme = theme_classic2(base_size=8, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_hf_chemo,]$heart_failure_grp_inc_fu/30), 
                    event = a1[a1_hf_chemo,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_hf_chemo,])
splots[[6]] <- ggsurvplot(fit1, data = a1[a1_hf_chemo,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = F,palette = 'jco',title='Heart failure (Received chemotherapy)',xlab='Time/months',
                          ylim = c(0,0.1),xlim = c(0,120),
                          ggtheme = theme_classic2(base_size=8, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_hf_rad_l,]$heart_failure_grp_inc_fu/30), 
                    event = a1[a1_hf_rad_l,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_hf_rad_l,])
splots[[7]] <- ggsurvplot(fit1, data = a1[a1_hf_rad_l,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = F,palette = 'jco',title='Heart failure (Received left-sided radiation therapy)',xlab='Time/months',
                          ylim = c(0,0.1),xlim = c(0,120),
                          ggtheme = theme_classic2(base_size=8, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_hf_horm,]$heart_failure_grp_inc_fu/30), 
                    event = a1[a1_hf_horm,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_hf_horm,])
splots[[8]] <- ggsurvplot(fit1, data = a1[a1_hf_horm,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = F,palette = 'jco',title='Heart failure (Received hormonal therapy)',xlab='Time/months',
                          ylim = c(0,0.1),xlim = c(0,120),
                          ggtheme = theme_classic2(base_size=8, base_family = "Arial"))



### stroke
surv_object <- Surv(time = floor(a1[a1_stroke,]$stroke_grp_inc_fu/30), 
                    event = a1[a1_stroke,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke,])
splots[[9]] <- ggsurvplot(fit1, data = a1[a1_stroke,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = F,palette = 'jco',title='Stroke (Entire cohort)',xlab='Time/months',
                          ylim = c(0,0.1),xlim = c(0,120),
                          ggtheme = theme_classic2(base_size=8, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_stroke_chemo,]$stroke_grp_inc_fu/30), 
                    event = a1[a1_stroke_chemo,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke_chemo,])
splots[[10]] <- ggsurvplot(fit1, data = a1[a1_stroke_chemo,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = F,palette = 'jco',title='Stroke (Received chemotherapy)',xlab='Time/months',
                          ylim = c(0,0.1),xlim = c(0,120),
                          ggtheme = theme_classic2(base_size=8, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_stroke_rad_l,]$stroke_grp_inc_fu/30), 
                    event = a1[a1_stroke_rad_l,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke_rad_l,])
splots[[11]] <- ggsurvplot(fit1, data = a1[a1_stroke_rad_l,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = F,palette = 'jco',title='Stroke (Received left-sided radiation therapy)',xlab='Time/months',
                          ylim = c(0,0.1),xlim = c(0,120),
                          ggtheme = theme_classic2(base_size=8, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_stroke_horm,]$stroke_grp_inc_fu/30), 
                    event = a1[a1_stroke_horm,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_stroke_horm,])
splots[[12]] <- ggsurvplot(fit1, data = a1[a1_stroke_horm,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = F,palette = 'jco',title='Stroke (Received hormonal therapy)',xlab='Time/months',
                          ylim = c(0,0.1),xlim = c(0,120),
                          ggtheme = theme_classic2(base_size=8, base_family = "Arial"))


# cardiomyopathy
surv_object <- Surv(time = floor(a1[a1_cm,]$cardiomyopathy_grp_inc_fu/30), 
                    event = a1[a1_cm,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_cm,])
splots[[13]] <- ggsurvplot(fit1, data = a1[a1_cm,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = F,palette = 'jco',title='Cardiomyopathy (Entire cohort)',xlab='Time/months',
                          ylim = c(0,0.005),xlim = c(0,120),
                          ggtheme = theme_classic2(base_size=8, base_family = "Arial"))


surv_object <- Surv(time = floor(a1[a1_cm_chemo,]$cardiomyopathy_grp_inc_fu/30), 
                    event = a1[a1_cm_chemo,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_cm_chemo,])
splots[[14]] <- ggsurvplot(fit1, data = a1[a1_cm_chemo,],fun = "event", size=1, risk.table = F,pval = T,
                           censor=F,conf.int = F,palette = 'jco',title='Cardiomyopathy (Received chemotherapy)',xlab='Time/months',
                           ylim = c(0,0.005),xlim = c(0,120),
                           ggtheme = theme_classic2(base_size=8, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_cm_rad_l,]$cardiomyopathy_grp_inc_fu/30), 
                    event = a1[a1_cm_rad_l,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_cm_rad_l,])
splots[[15]] <- ggsurvplot(fit1, data = a1[a1_cm_rad_l,],fun = "event", size=1, risk.table = F,pval = T,
                           censor=F,conf.int = F,palette = 'jco',title='Cardiomyopathy (Received left-sided radiation therapy)',xlab='Time/months',
                           ylim = c(0,0.005),xlim = c(0,120),
                           ggtheme = theme_classic2(base_size=8, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_cm_horm,]$cardiomyopathy_grp_inc_fu/30), 
                    event = a1[a1_cm_horm,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_cm_horm,])
splots[[16]] <- ggsurvplot(fit1, data = a1[a1_cm_horm,],fun = "event", size=1, risk.table = F,pval = T,
                           censor=F,conf.int = F,palette = 'jco',title='Cardiomyopathy (Received hormonal therapy)',xlab='Time/months',
                           ylim = c(0,0.005),xlim = c(0,120),
                           ggtheme = theme_classic2(base_size=8, base_family = "Arial"))



splots <- list()
# thromboembolic disease 
surv_object <- Surv(time = floor(a1[a1_dvt,]$venous_thromboembolic_disease_grp_inc_fu/30), 
                    event = a1[a1_dvt,]$venous_thromboembolic_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_dvt,])
splots[[1]] <- ggsurvplot(fit1, data = a1[a1_dvt,],fun = "event", size=1, risk.table = F,pval = T,
                           censor=F,conf.int = F,palette = 'jco',title='A. Entire cohort',xlab='Time/months',
                           ylim = c(0,0.04),xlim = c(0,120),
                          legend.title = "",
                          legend.labs = c("Case", "Control"),
                           ggtheme = theme_classic2(base_size=10, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_dvt_chemo,]$venous_thromboembolic_disease_grp_inc_fu/30), 
                    event = a1[a1_dvt_chemo,]$venous_thromboembolic_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_dvt_chemo,])
splots[[3]] <- ggsurvplot(fit1, data = a1[a1_dvt_chemo,],fun = "event", size=1, risk.table = F,pval = T,
                           censor=F,conf.int = F,palette = 'jco',title='B. Received chemotherapy',xlab='Time/months',
                           ylim = c(0,0.05),xlim = c(0,120),
                          legend.title = "",
                          legend.labs = c("Case", "Control"),
                           ggtheme = theme_classic2(base_size=10, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_dvt_rad_l,]$venous_thromboembolic_disease_grp_inc_fu/30), 
                    event = a1[a1_dvt_rad_l,]$venous_thromboembolic_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_dvt_rad_l,])
splots[[2]] <- ggsurvplot(fit1, data = a1[a1_dvt_rad_l,],fun = "event", size=1, risk.table = F,pval = T,
                           censor=F,conf.int = F,palette = 'jco',title='C. Received left-sided radiation therapy',xlab='Time/months',
                           ylim = c(0,0.04),xlim = c(0,120),
                          legend.title = "",
                          legend.labs = c("Case", "Control"),
                           ggtheme = theme_classic2(base_size=10, base_family = "Arial"))

surv_object <- Surv(time = floor(a1[a1_dvt_horm,]$venous_thromboembolic_disease_grp_inc_fu/30), 
                    event = a1[a1_dvt_horm,]$venous_thromboembolic_disease_grp_inc)
fit1 <- survfit(surv_object ~ group, data = a1[a1_dvt_horm,])
splots[[4]] <- ggsurvplot(fit1, data = a1[a1_dvt_horm,],fun = "event", size=1, risk.table = F,pval = T,
                           censor=F,conf.int = F,palette = 'jco',title='D. Received hormonal therapy',xlab='Time/months',
                           ylim = c(0,0.04),xlim = c(0,120),
                          legend.title = "",
                          legend.labs = c("Case", "Control"),
                           ggtheme = theme_classic2(base_size=10, base_family = "Arial"))





fig <- arrange_ggsurvplots(splots, ncol = 2, nrow = 2, title = "Cumulative incidence of developing other CVD events (secondary outcomes)")




a1$cvdcombo_grp_prev2 <- ifelse(a1$ischemic_heart_disease_grp_prev==1 | 
                                  a1$stroke_grp_prev==1 |  
                                  a1$arrhythmia_grp_prev==1|
                                  a1$heart_failure_grp_prev==1|
                                  a1$venous_thromboembolic_disease_grp_prev == 1, 1, 0)


a1$cvdcombo_grp_inc2 = ifelse(a1$ischemic_heart_disease_grp_inc==1 | 
                                    a1$stroke_grp_inc==1 |  
                                    a1$arrhythmia_grp_inc==1|
                                    a1$heart_failure_grp_inc==1|
                                    a1$venous_thromboembolic_disease_grp_inc == 1, 1, 0)

a1$cvdcombo_grp_incdt2 <- apply(a1[,c("ischemic_heart_disease_grp_incdt",
                                     "stroke_grp_incdt",
                                     "arrhythmia_grp_incdt",
                                     "heart_failure_grp_incdt",
                                     "venous_thromboembolic_disease_grp_incdt")],1,min,na.rm=T)

a1$cvdcombo_grp_inc_fu2 <- a1$cvdcombo_grp_incdt2
a1$cvdcombo_grp_inc_fu2[which(a1$cvdcombo_grp_inc2==0)] <- 
  a1$censor_dt[which(a1$cvdcombo_grp_inc2==0)] 

a1_combo2 = which(a1$cvdcombo_grp_prev2==0)

pwcard = as.data.frame(read_sas('pw_cardiotox.sas7bdat'))

pwcard$anthra_yn = chemo$anthra_yn[match(pwcard$CVD_studyid, chemo$CVD_studyid)]
pwcard$tras_yn = chemo$tras_yn[match(pwcard$CVD_studyid, chemo$CVD_studyid)]
pwcard$taxane_yn = chemo$taxane_yn[match(pwcard$CVD_studyid, chemo$CVD_studyid)]

pwcard$ischemic_heart_disease_grp_inc = a1$ischemic_heart_disease_grp_inc[match(pwcard$CVD_studyid, a1$cvd_studyid)]
pwcard$heart_failure_grp_inc = a1$heart_failure_grp_inc[match(pwcard$CVD_studyid, a1$cvd_studyid)]
pwcard$stroke_grp_inc = a1$stroke_grp_inc[match(pwcard$CVD_studyid, a1$cvd_studyid)]
pwcard$cardiomyopathy_grp_inc = a1$cardiomyopathy_grp_inc[match(pwcard$CVD_studyid, a1$cvd_studyid)]


pwcard$cvdcombo_grp_inc2 = ifelse(pwcard$ischemic_heart_disease_grp_inc==1 | 
                                   pwcard$stroke_grp_inc==1 |  
                                    pwcard$arrhythmia_grp_inc==1|
                                    pwcard$heart_failure_grp_inc==1|
                                    pwcard$venous_thromboembolic_disease_grp_inc == 1, 1, 0)

pwcard$cvdcombo_grp_inc = a1$cvdcombo_grp_inc[match(pwcard$CVD_studyid, a1$cvd_studyid)]




pwcard$arrhythmia_grp_inc = a1$arrhythmia_grp_inc[match(pwcard$CVD_studyid, a1$cvd_studyid)]
pwcard$cardiac_arrest_grp_inc = a1$cardiac_arrest_grp_inc[match(pwcard$CVD_studyid, a1$cvd_studyid)]
pwcard$carotid_disease_grp_inc = a1$carotid_disease_grp_inc[match(pwcard$CVD_studyid, a1$cvd_studyid)]
pwcard$myocarditis_pericarditis_grp_inc = a1$myocarditis_pericarditis_grp_inc[match(pwcard$CVD_studyid, a1$cvd_studyid)]
pwcard$tia_grp_inc = a1$tia_grp_inc[match(pwcard$CVD_studyid, a1$cvd_studyid)]
pwcard$venous_thromboembolic_disease_grp_inc = a1$venous_thromboembolic_disease_grp_inc[match(pwcard$CVD_studyid, a1$cvd_studyid)]
pwcard$valvular_disease_grp_inc = a1$valvular_disease_grp_inc[match(pwcard$CVD_studyid, a1$cvd_studyid)]
  
  

pwcard$cvdrf_htn = a1$cvdrf_htn[match(pwcard$CVD_studyid, a1$cvd_studyid)]
pwcard$cvdrf_diab = a1$cvdrf_diab[match(pwcard$CVD_studyid, a1$cvd_studyid)]
pwcard$cvdrf_dyslipid = a1$cvdrf_dyslipid[match(pwcard$CVD_studyid, a1$cvd_studyid)]
pwcard$cvdrfcombo = a1$cvdrfcombo[match(pwcard$CVD_studyid, a1$cvd_studyid)]

pwcard$index_date = a1$index_date[match(pwcard$CVD_studyid, a1$cvd_studyid)]
pwcard$end_date = "2020-12-31"
pwcard$num_days = as.Date(pwcard$end_date)-as.Date(pwcard$index_date)
pwcard$num_years = pwcard$num_days/365
pwcard$ind_years = NA
pwcard$ind_years[pwcard$num_years >= 10] = 1
pwcard$ind_years[pwcard$num_years < 10 & pwcard$num_years >= 9] = 2
pwcard$ind_years[pwcard$num_years < 9 & pwcard$num_years >= 8] = 3
pwcard$ind_years[pwcard$num_years < 8 & pwcard$num_years >= 7] = 4

table(pwcard$ind_years)

# function to calculate time specific incidence rate
inc_rate <- function(index, outcomevar, timevar){
  surv_dat <- a1[index,c(outcomevar,timevar,'group')]
  surv_dat$time <- round(surv_dat[,2]/30)
  surv_object <- Surv(time = surv_dat$time, 
                      event = surv_dat[,1])
  fit1 <- survfit(surv_object ~ group, data = surv_dat)
  logrank.24 <- survdiff(Surv(time = surv_dat$time[surv_dat$time <= 24], 
                              event = surv_dat[surv_dat$time <= 24,1]) ~ group, 
                         data = surv_dat[surv_dat$time <= 24,])
  logrank.60 <- survdiff(Surv(time = surv_dat$time[surv_dat$time <= 60], 
                              event = surv_dat[surv_dat$time <= 60,1]) ~ group, 
                         data = surv_dat[surv_dat$time <= 60,])
  logrank.84 <- survdiff(Surv(time = surv_dat$time[surv_dat$time <= 84], 
                              event = surv_dat[surv_dat$time <= 84,1]) ~ group, 
                         data = surv_dat[surv_dat$time <= 84,])
  logrank.96 <- survdiff(Surv(time = surv_dat$time[surv_dat$time <= 96], 
                              event = surv_dat[surv_dat$time <= 96,1]) ~ group, 
                         data = surv_dat[surv_dat$time <= 96,])
  logrank.108 <- survdiff(Surv(time = surv_dat$time[surv_dat$time <= 108], 
                              event = surv_dat[surv_dat$time <= 108,1]) ~ group, 
                         data = surv_dat[surv_dat$time <= 108,])
  logrank.120 <- survdiff(Surv(time = surv_dat$time[surv_dat$time <= 120], 
                               event = surv_dat[surv_dat$time <= 120,1]) ~ group, 
                          data = surv_dat[surv_dat$time <= 120,])
  pvalue.24 <- round(1 - pchisq(logrank.24$chisq, length(logrank.24$n) - 1),2)
  pvalue.60 <- round(1 - pchisq(logrank.60$chisq, length(logrank.60$n) - 1),2)
  pvalue.84 <- round(1 - pchisq(logrank.84$chisq, length(logrank.84$n) - 1),2)
  pvalue.96 <- round(1 - pchisq(logrank.96$chisq, length(logrank.96$n) - 1),2)
  pvalue.108 <- round(1 - pchisq(logrank.108$chisq, length(logrank.108$n) - 1),2)
  pvalue.120 <- round(1 - pchisq(logrank.120$chisq, length(logrank.120$n) - 1),2)
  pvalue = c(pvalue.24, pvalue.60, pvalue.84, pvalue.96, pvalue.108, pvalue.120)
  sum1 <- summary(fit1, times=c(24,60,84,96,108,120))
  sum2 <- 1-cbind(sum1$surv,sum1$upper,sum1$lower)
  sum2 <- cbind(round2(sum2[,1],5),
                paste0("(",round2(sum2[,2],1),", ",round2(sum2[,3],1),")"))
  sum3 <- cbind(outcomevar,sum2[1:6,],pvalue,sum2[7:12,])
  sum3
}




inc = inc_rate(a1_ihd,'ischemic_heart_disease_grp_inc','ischemic_heart_disease_grp_inc_fu')
table(pwcard$ind_years)[1]*as.numeric(inc[6,2]) + table(pwcard$ind_years)[2]*as.numeric(inc[5,2]) + 
  table(pwcard$ind_years)[3]*as.numeric(inc[4,2]) + table(pwcard$ind_years)[4]*as.numeric(inc[3,2])

inc = inc_rate(a1_stroke,'stroke_grp_inc','stroke_grp_inc_fu')
table(pwcard$ind_years)[1]*as.numeric(inc[6,2]) + table(pwcard$ind_years)[2]*as.numeric(inc[5,2]) + 
  table(pwcard$ind_years)[3]*as.numeric(inc[4,2]) + table(pwcard$ind_years)[4]*as.numeric(inc[3,2])

inc = inc_rate(a1_hf,'heart_failure_grp_inc','heart_failure_grp_inc_fu')
table(pwcard$ind_years)[1]*as.numeric(inc[6,2]) + table(pwcard$ind_years)[2]*as.numeric(inc[5,2]) + 
  table(pwcard$ind_years)[3]*as.numeric(inc[4,2]) + table(pwcard$ind_years)[4]*as.numeric(inc[3,2])

inc = inc_rate(a1_arrhythmia,'arrhythmia_grp_inc','arrhythmia_grp_inc_fu')
table(pwcard$ind_years)[1]*as.numeric(inc[6,2]) + table(pwcard$ind_years)[2]*as.numeric(inc[5,2]) + 
  table(pwcard$ind_years)[3]*as.numeric(inc[4,2]) + table(pwcard$ind_years)[4]*as.numeric(inc[3,2])


inc = inc_rate(a1_combo2,'cvdcombo_grp_inc2','cvdcombo_grp_inc_fu2')
table(pwcard$ind_years)[1]*as.numeric(inc[6,2]) + table(pwcard$ind_years)[2]*as.numeric(inc[5,2]) + 
  table(pwcard$ind_years)[3]*as.numeric(inc[4,2]) + table(pwcard$ind_years)[4]*as.numeric(inc[3,2])


inc = inc_rate(a1_htn,'cvdrf_htn','cvdrf_htn_fu')
table(pwcard$ind_years)[1]*as.numeric(inc[6,2]) + table(pwcard$ind_years)[2]*as.numeric(inc[5,2]) + 
  table(pwcard$ind_years)[3]*as.numeric(inc[4,2]) + table(pwcard$ind_years)[4]*as.numeric(inc[3,2])

inc = inc_rate(a1_diab,'cvdrf_diab','cvdrf_diab_fu')
table(pwcard$ind_years)[1]*as.numeric(inc[6,2]) + table(pwcard$ind_years)[2]*as.numeric(inc[5,2]) + 
  table(pwcard$ind_years)[3]*as.numeric(inc[4,2]) + table(pwcard$ind_years)[4]*as.numeric(inc[3,2])

inc = inc_rate(a1_dyslipid,'cvdrf_dyslipid','cvdrf_dyslipid_fu')
table(pwcard$ind_years)[1]*as.numeric(inc[6,2]) + table(pwcard$ind_years)[2]*as.numeric(inc[5,2]) + 
  table(pwcard$ind_years)[3]*as.numeric(inc[4,2]) + table(pwcard$ind_years)[4]*as.numeric(inc[3,2])

inc = inc_rate(a1_rfcombo,'cvdrfcombo','cvdrfcombo_fu')
table(pwcard$ind_years)[1]*as.numeric(inc[6,2]) + table(pwcard$ind_years)[2]*as.numeric(inc[5,2]) + 
  table(pwcard$ind_years)[3]*as.numeric(inc[4,2]) + table(pwcard$ind_years)[4]*as.numeric(inc[3,2])



inc_rate2 <- function(outcomevar, timevar){
  surv_dat <- pwcard[,c(outcomevar,timevar)]
  surv_dat$time <- round(surv_dat[,2]/30)
  surv_object <- Surv(time = surv_dat$time, 
                      event = surv_dat[,1])
  fit1 <- survfit(surv_object ~ 1, data = surv_dat)
  sum1 <- summary(fit1, times=c(24,60,84,96,108,120))
  sum2 <- 1-cbind(sum1$surv,sum1$upper,sum1$lower)
  sum2 <- round2(sum2[,1],5)
  sum2
}
pwcard$venous_thromboembolic_disease_grp_inc_fu = a1$venous_thromboembolic_disease_grp_inc_fu[match(pwcard$CVD_studyid, a1$cvd_studyid)]

inc = inc_rate2('venous_thromboembolic_disease_grp_inc','venous_thromboembolic_disease_grp_inc_fu')
table(pwcard$ind_years)[1]*as.numeric(inc[6]) + table(pwcard$ind_years)[2]*as.numeric(inc[5]) + 
  table(pwcard$ind_years)[3]*as.numeric(inc[4]) + table(pwcard$ind_years)[4]*as.numeric(inc[3])





