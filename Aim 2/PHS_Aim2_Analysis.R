##### Pathways Heart Study Aim 2 Analysis ######
### Hanjie Shen, 8/6/2020

require(tidyverse)
require(haven)
require(ggplot2)
require(dplyr)
require(plyr)
require(naniar)
require(reshape2)
require(survival)
require(survminer)
require(cmprsk)
require(pec)
require(LtAtStructuR)
require(stremr)
require(magrittr)
require(grDevices)



###################################################################################
######################### Part 2: Statistical Analysis ############################
###################################################################################
setwd("/Users/hshen2/Desktop/Pathways Heart Study/Aim 2 Results")
#setwd("C:/Users/shenh/Desktop/Pathways Heart Study")
round2 = function(x, k) trimws(format(round(x, k), nsmall=k))



############ Model 1 & 2: Unadjusted KM Curves and Cox Regression #####################
# Summary of Treatment for Primary CVD Outcomes
sum.trt.tab = matrix(NA, 5,14)
colnames(sum.trt.tab) = c("anathra_y","anathra_n","tras_y","tras_n","tax_y","tax_n","ai_y","ai_n","tam_y","tam_n","rad_y","rad_n","rad_l_y","rad_l_n")
rownames(sum.trt.tab) = c("ihd","hf","card","stroke","any")

sum.trt.tab[1,1:2] = table(phs.aim2.data[a2_ihd,"anthra"])[2:1]
sum.trt.tab[1,3:4] = table(phs.aim2.data[a2_ihd,"tras"])[2:1]
sum.trt.tab[1,5:6] = table(phs.aim2.data[a2_ihd,"taxane"])[2:1]
sum.trt.tab[1,7:8] = table(phs.aim2.data[a2_ihd,"ai"])[2:1]
sum.trt.tab[1,9:10] = table(phs.aim2.data[a2_ihd,"tam"])[2:1]
sum.trt.tab[1,11:12] = table(phs.aim2.data[a2_ihd,"rad"])[2:1]
sum.trt.tab[1,13:14] = table(phs.aim2.data[a2_ihd_l,"rad"])[2:1]

sum.trt.tab[2,1:2] = table(phs.aim2.data[a2_hf,"anthra"])[2:1]
sum.trt.tab[2,3:4] = table(phs.aim2.data[a2_hf,"tras"])[2:1]
sum.trt.tab[2,5:6] = table(phs.aim2.data[a2_hf,"taxane"])[2:1]
sum.trt.tab[2,7:8] = table(phs.aim2.data[a2_hf,"ai"])[2:1]
sum.trt.tab[2,9:10] = table(phs.aim2.data[a2_hf,"tam"])[2:1]
sum.trt.tab[2,11:12] = table(phs.aim2.data[a2_hf,"rad"])[2:1]
sum.trt.tab[2,13:14] = table(phs.aim2.data[a2_hf_l,"rad"])[2:1]

sum.trt.tab[3,1:2] = table(phs.aim2.data[a2_cm,"anthra"])[2:1]
sum.trt.tab[3,3:4] = table(phs.aim2.data[a2_cm,"tras"])[2:1]
sum.trt.tab[3,5:6] = table(phs.aim2.data[a2_cm,"taxane"])[2:1]
sum.trt.tab[3,7:8] = table(phs.aim2.data[a2_cm,"ai"])[2:1]
sum.trt.tab[3,9:10] = table(phs.aim2.data[a2_cm,"tam"])[2:1]
sum.trt.tab[3,11:12] = table(phs.aim2.data[a2_cm,"rad"])[2:1]
sum.trt.tab[3,13:14] = table(phs.aim2.data[a2_cm_l,"rad"])[2:1]

sum.trt.tab[4,1:2] = table(phs.aim2.data[a2_stroke,"anthra"])[2:1]
sum.trt.tab[4,3:4] = table(phs.aim2.data[a2_stroke,"tras"])[2:1]
sum.trt.tab[4,5:6] = table(phs.aim2.data[a2_stroke,"taxane"])[2:1]
sum.trt.tab[4,7:8] = table(phs.aim2.data[a2_stroke,"ai"])[2:1]
sum.trt.tab[4,9:10] = table(phs.aim2.data[a2_stroke,"tam"])[2:1]
sum.trt.tab[4,11:12] = table(phs.aim2.data[a2_stroke,"rad"])[2:1]
sum.trt.tab[4,13:14] = table(phs.aim2.data[a2_stroke_l,"rad"])[2:1]

sum.trt.tab[5,1:2] = table(phs.aim2.data[a2_combo,"anthra"])[2:1]
sum.trt.tab[5,3:4] = table(phs.aim2.data[a2_combo,"tras"])[2:1]
sum.trt.tab[5,5:6] = table(phs.aim2.data[a2_combo,"taxane"])[2:1]
sum.trt.tab[5,7:8] = table(phs.aim2.data[a2_combo,"ai"])[2:1]
sum.trt.tab[5,9:10] = table(phs.aim2.data[a2_combo,"tam"])[2:1]
sum.trt.tab[5,11:12] = table(phs.aim2.data[a2_combo,"rad"])[2:1]
sum.trt.tab[5,13:14] = table(phs.aim2.data[a2_combo_l,"rad"])[2:1]

write.csv(sum.trt.tab,"Table2a.csv")

# LVEF > 50
sum.trt.tab = matrix(NA, 5,4)
colnames(sum.trt.tab) = c("anathra_y","anathra_n","tras_y","tras_n")
rownames(sum.trt.tab) = c("ihd","hf","card","stroke","any")

sum.trt.tab[1,1:2] = table(phs.aim2.data[a2_ihd_ef,"anthra"])[2:1]
sum.trt.tab[1,3:4] = table(phs.aim2.data[a2_ihd_ef,"tras"])[2:1]


sum.trt.tab[2,1:2] = table(phs.aim2.data[a2_hf_ef,"anthra"])[2:1]
sum.trt.tab[2,3:4] = table(phs.aim2.data[a2_hf_ef,"tras"])[2:1]


sum.trt.tab[3,1:2] = table(phs.aim2.data[a2_cm_ef,"anthra"])[2:1]
sum.trt.tab[3,3:4] = table(phs.aim2.data[a2_cm_ef,"tras"])[2:1]


sum.trt.tab[4,1:2] = table(phs.aim2.data[a2_stroke_ef,"anthra"])[2:1]
sum.trt.tab[4,3:4] = table(phs.aim2.data[a2_stroke_ef,"tras"])[2:1]

sum.trt.tab[5,1:2] = table(phs.aim2.data[a2_combo_ef,"anthra"])[2:1]
sum.trt.tab[5,3:4] = table(phs.aim2.data[a2_combo_ef,"tras"])[2:1]


write.csv(sum.trt.tab,"Table2a_LVEFgreater50.csv")


# post menopausal
sum.trt.tab = matrix(NA, 5,2)
colnames(sum.trt.tab) = c("ai_y","ai_n")
rownames(sum.trt.tab) = c("ihd","hf","card","stroke","any")

sum.trt.tab[1,1:2] = table(phs.aim2.data[a2_ihd_post_menop,"ai"])[2:1]
sum.trt.tab[2,1:2] = table(phs.aim2.data[a2_hf_post_menop,"ai"])[2:1]
sum.trt.tab[3,1:2] = table(phs.aim2.data[a2_cm_post_menop,"ai"])[2:1]
sum.trt.tab[4,1:2] = table(phs.aim2.data[a2_stroke_post_menop,"ai"])[2:1]
sum.trt.tab[5,1:2] = table(phs.aim2.data[a2_combo_post_menop,"ai"])[2:1]

write.csv(sum.trt.tab,"Table2a_PostMenop.csv")

# pre menopausal
sum.trt.tab = matrix(NA, 5,2)
colnames(sum.trt.tab) = c("ai_y","ai_n")
rownames(sum.trt.tab) = c("ihd","hf","card","stroke","any")

sum.trt.tab[1,1:2] = table(phs.aim2.data[a2_ihd_pre_menop,"ai"])[2:1]
sum.trt.tab[2,1:2] = table(phs.aim2.data[a2_hf_pre_menop,"ai"])[2:1]
sum.trt.tab[3,1:2] = table(phs.aim2.data[a2_cm_pre_menop,"ai"])[2:1]
sum.trt.tab[4,1:2] = table(phs.aim2.data[a2_stroke_pre_menop,"ai"])[2:1]
sum.trt.tab[5,1:2] = table(phs.aim2.data[a2_combo_pre_menop,"ai"])[2:1]

write.csv(sum.trt.tab,"Table2a_PreMenop.csv")



# Summary of Events of Primary CVD Outcomes by Treatment Group
sum.trt.event.tab = matrix(NA, 5,42)
colnames(sum.trt.event.tab) = c(rep("anathra_y",3),rep("anathra_n",3),rep("tras_y",3),rep("tras_n",3),
                                rep("tax_y",3),rep("tax_n",3),rep("ai_y",3),rep("ai_n",3),rep("tam_y",3),
                                rep("tam_n",3),rep("rad_y",3),rep("rad_n",3),rep("rad_l_y",3),rep("rad_l_n",3))
rownames(sum.trt.event.tab) = c("ihd","hf","card","stroke","any")

sum.trt.event.tab[1,1:2] = table(phs.aim2.data[a2_ihd,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd,"anthra"])[2:1,2]
sum.trt.event.tab[1,4:5] = table(phs.aim2.data[a2_ihd,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd,"anthra"])[2:1,1]
sum.trt.event.tab[1,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_ihd),"anthra"])[2:1]
sum.trt.event.tab[1,7:8] = table(phs.aim2.data[a2_ihd,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd,"tras"])[2:1,2]
sum.trt.event.tab[1,10:11] = table(phs.aim2.data[a2_ihd,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd,"tras"])[2:1,1]
sum.trt.event.tab[1,c(9,12)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_ihd),"tras"])[2:1]
sum.trt.event.tab[1,13:14] = table(phs.aim2.data[a2_ihd,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd,"taxane"])[2:1,2]
sum.trt.event.tab[1,16:17] = table(phs.aim2.data[a2_ihd,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd,"taxane"])[2:1,1]
sum.trt.event.tab[1,c(15,18)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_ihd),"taxane"])[2:1]
sum.trt.event.tab[1,19:20] = table(phs.aim2.data[a2_ihd,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd,"ai"])[2:1,2]
sum.trt.event.tab[1,22:23] = table(phs.aim2.data[a2_ihd,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd,"ai"])[2:1,1]
sum.trt.event.tab[1,c(21,24)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_ihd),"ai"])[2:1]
sum.trt.event.tab[1,25:26] = table(phs.aim2.data[a2_ihd,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd,"tam"])[2:1,2]
sum.trt.event.tab[1,28:29] = table(phs.aim2.data[a2_ihd,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd,"tam"])[2:1,1]
sum.trt.event.tab[1,c(27,30)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_ihd),"tam"])[2:1]
sum.trt.event.tab[1,31:32] = table(phs.aim2.data[a2_ihd,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd,"rad"])[2:1,2]
sum.trt.event.tab[1,34:35] = table(phs.aim2.data[a2_ihd,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd,"rad"])[2:1,1]
sum.trt.event.tab[1,c(33,36)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_ihd),"rad"])[2:1]
sum.trt.event.tab[1,37:38] = table(phs.aim2.data[a2_ihd_l,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd_l,"rad"])[2:1,2]
sum.trt.event.tab[1,40:41] = table(phs.aim2.data[a2_ihd_l,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd_l,"rad"])[2:1,1]
sum.trt.event.tab[1,c(39,42)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_ihd_l),"rad"])[2:1]

sum.trt.event.tab[2,1:2] = table(phs.aim2.data[a2_hf,"heart_failure_grp_inc"], phs.aim2.data[a2_hf,"anthra"])[2:1,2]
sum.trt.event.tab[2,4:5] = table(phs.aim2.data[a2_hf,"heart_failure_grp_inc"], phs.aim2.data[a2_hf,"anthra"])[2:1,1]
sum.trt.event.tab[2,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_hf),"anthra"])[2:1]
sum.trt.event.tab[2,7:8] = table(phs.aim2.data[a2_hf,"heart_failure_grp_inc"], phs.aim2.data[a2_hf,"tras"])[2:1,2]
sum.trt.event.tab[2,10:11] = table(phs.aim2.data[a2_hf,"heart_failure_grp_inc"], phs.aim2.data[a2_hf,"tras"])[2:1,1]
sum.trt.event.tab[2,c(9,12)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_hf),"tras"])[2:1]
sum.trt.event.tab[2,13:14] = table(phs.aim2.data[a2_hf,"heart_failure_grp_inc"], phs.aim2.data[a2_hf,"taxane"])[2:1,2]
sum.trt.event.tab[2,16:17] = table(phs.aim2.data[a2_hf,"heart_failure_grp_inc"], phs.aim2.data[a2_hf,"taxane"])[2:1,1]
sum.trt.event.tab[2,c(15,18)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_hf),"taxane"])[2:1]
sum.trt.event.tab[2,19:20] = table(phs.aim2.data[a2_hf,"heart_failure_grp_inc"], phs.aim2.data[a2_hf,"ai"])[2:1,2]
sum.trt.event.tab[2,22:23] = table(phs.aim2.data[a2_hf,"heart_failure_grp_inc"], phs.aim2.data[a2_hf,"ai"])[2:1,1]
sum.trt.event.tab[2,c(21,24)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_hf),"ai"])[2:1]
sum.trt.event.tab[2,25:26] = table(phs.aim2.data[a2_hf,"heart_failure_grp_inc"], phs.aim2.data[a2_hf,"tam"])[2:1,2]
sum.trt.event.tab[2,28:29] = table(phs.aim2.data[a2_hf,"heart_failure_grp_inc"], phs.aim2.data[a2_hf,"tam"])[2:1,1]
sum.trt.event.tab[2,c(27,30)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_hf),"tam"])[2:1]
sum.trt.event.tab[2,31:32] = table(phs.aim2.data[a2_hf,"heart_failure_grp_inc"], phs.aim2.data[a2_hf,"rad"])[2:1,2]
sum.trt.event.tab[2,34:35] = table(phs.aim2.data[a2_hf,"heart_failure_grp_inc"], phs.aim2.data[a2_hf,"rad"])[2:1,1]
sum.trt.event.tab[2,c(33,36)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_hf),"rad"])[2:1]
sum.trt.event.tab[2,37:38] = table(phs.aim2.data[a2_hf_l,"heart_failure_grp_inc"], phs.aim2.data[a2_hf_l,"rad"])[2:1,2]
sum.trt.event.tab[2,40:41] = table(phs.aim2.data[a2_hf_l,"heart_failure_grp_inc"], phs.aim2.data[a2_hf_l,"rad"])[2:1,1]
sum.trt.event.tab[2,c(39,42)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_hf_l),"rad"])[2:1]

sum.trt.event.tab[3,1:2] = table(phs.aim2.data[a2_cm,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm,"anthra"])[2:1,2]
sum.trt.event.tab[3,4:5] = table(phs.aim2.data[a2_cm,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm,"anthra"])[2:1,1]
sum.trt.event.tab[3,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_cm),"anthra"])[2:1]
sum.trt.event.tab[3,7:8] = table(phs.aim2.data[a2_cm,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm,"tras"])[2:1,2]
sum.trt.event.tab[3,10:11] = table(phs.aim2.data[a2_cm,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm,"tras"])[2:1,1]
sum.trt.event.tab[3,c(9,12)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_cm),"tras"])[2:1]
sum.trt.event.tab[3,13:14] = table(phs.aim2.data[a2_cm,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm,"taxane"])[2:1,2]
sum.trt.event.tab[3,16:17] = table(phs.aim2.data[a2_cm,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm,"taxane"])[2:1,1]
sum.trt.event.tab[3,c(15,18)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_cm),"taxane"])[2:1]
sum.trt.event.tab[3,19:20] = table(phs.aim2.data[a2_cm,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm,"ai"])[2:1,2]
sum.trt.event.tab[3,22:23] = table(phs.aim2.data[a2_cm,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm,"ai"])[2:1,1]
sum.trt.event.tab[3,c(21,24)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_cm),"ai"])[2:1]
sum.trt.event.tab[3,25:26] = table(phs.aim2.data[a2_cm,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm,"tam"])[2:1,2]
sum.trt.event.tab[3,28:29] = table(phs.aim2.data[a2_cm,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm,"tam"])[2:1,1]
sum.trt.event.tab[3,c(27,30)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_cm),"tam"])[2:1]
sum.trt.event.tab[3,31:32] = table(phs.aim2.data[a2_cm,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm,"rad"])[2:1,2]
sum.trt.event.tab[3,34:35] = table(phs.aim2.data[a2_cm,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm,"rad"])[2:1,1]
sum.trt.event.tab[3,c(33,36)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_cm),"rad"])[2:1]
sum.trt.event.tab[3,37:38] = table(phs.aim2.data[a2_cm_l,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm_l,"rad"])[2:1,2]
sum.trt.event.tab[3,40:41] = table(phs.aim2.data[a2_cm_l,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm_l,"rad"])[2:1,1]
sum.trt.event.tab[3,c(39,42)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_cm_l),"rad"])[2:1]

sum.trt.event.tab[4,1:2] = table(phs.aim2.data[a2_stroke,"stroke_grp_inc"], phs.aim2.data[a2_stroke,"anthra"])[2:1,2]
sum.trt.event.tab[4,4:5] = table(phs.aim2.data[a2_stroke,"stroke_grp_inc"], phs.aim2.data[a2_stroke,"anthra"])[2:1,1]
sum.trt.event.tab[4,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_stroke),"anthra"])[2:1]
sum.trt.event.tab[4,7:8] = table(phs.aim2.data[a2_stroke,"stroke_grp_inc"], phs.aim2.data[a2_stroke,"tras"])[2:1,2]
sum.trt.event.tab[4,10:11] = table(phs.aim2.data[a2_stroke,"stroke_grp_inc"], phs.aim2.data[a2_stroke,"tras"])[2:1,1]
sum.trt.event.tab[4,c(9,12)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_stroke),"tras"])[2:1]
sum.trt.event.tab[4,13:14] = table(phs.aim2.data[a2_stroke,"stroke_grp_inc"], phs.aim2.data[a2_stroke,"taxane"])[2:1,2]
sum.trt.event.tab[4,16:17] = table(phs.aim2.data[a2_stroke,"stroke_grp_inc"], phs.aim2.data[a2_stroke,"taxane"])[2:1,1]
sum.trt.event.tab[4,c(15,18)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_stroke),"taxane"])[2:1]
sum.trt.event.tab[4,19:20] = table(phs.aim2.data[a2_stroke,"stroke_grp_inc"], phs.aim2.data[a2_stroke,"ai"])[2:1,2]
sum.trt.event.tab[4,22:23] = table(phs.aim2.data[a2_stroke,"stroke_grp_inc"], phs.aim2.data[a2_stroke,"ai"])[2:1,1]
sum.trt.event.tab[4,c(21,24)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_stroke),"ai"])[2:1]
sum.trt.event.tab[4,25:26] = table(phs.aim2.data[a2_stroke,"stroke_grp_inc"], phs.aim2.data[a2_stroke,"tam"])[2:1,2]
sum.trt.event.tab[4,28:29] = table(phs.aim2.data[a2_stroke,"stroke_grp_inc"], phs.aim2.data[a2_stroke,"tam"])[2:1,1]
sum.trt.event.tab[4,c(27,30)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_stroke),"tam"])[2:1]
sum.trt.event.tab[4,31:32] = table(phs.aim2.data[a2_stroke,"stroke_grp_inc"], phs.aim2.data[a2_stroke,"rad"])[2:1,2]
sum.trt.event.tab[4,34:35] = table(phs.aim2.data[a2_stroke,"stroke_grp_inc"], phs.aim2.data[a2_stroke,"rad"])[2:1,1]
sum.trt.event.tab[4,c(33,36)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_stroke),"rad"])[2:1]
sum.trt.event.tab[4,37:38] = table(phs.aim2.data[a2_stroke_l,"stroke_grp_inc"], phs.aim2.data[a2_stroke_l,"rad"])[2:1,2]
sum.trt.event.tab[4,40:41] = table(phs.aim2.data[a2_stroke_l,"stroke_grp_inc"], phs.aim2.data[a2_stroke_l,"rad"])[2:1,1]
sum.trt.event.tab[4,c(39,42)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_stroke_l),"rad"])[2:1]

sum.trt.event.tab[5,1:2] = table(phs.aim2.data[a2_combo,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo,"anthra"])[2:1,2]
sum.trt.event.tab[5,4:5] = table(phs.aim2.data[a2_combo,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo,"anthra"])[2:1,1]
sum.trt.event.tab[5,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_combo),"anthra"])[2:1]
sum.trt.event.tab[5,7:8] = table(phs.aim2.data[a2_combo,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo,"tras"])[2:1,2]
sum.trt.event.tab[5,10:11] = table(phs.aim2.data[a2_combo,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo,"tras"])[2:1,1]
sum.trt.event.tab[5,c(9,12)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_combo),"tras"])[2:1]
sum.trt.event.tab[5,13:14] = table(phs.aim2.data[a2_combo,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo,"taxane"])[2:1,2]
sum.trt.event.tab[5,16:17] = table(phs.aim2.data[a2_combo,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo,"taxane"])[2:1,1]
sum.trt.event.tab[5,c(15,18)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_combo),"taxane"])[2:1]
sum.trt.event.tab[5,19:20] = table(phs.aim2.data[a2_combo,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo,"ai"])[2:1,2]
sum.trt.event.tab[5,22:23] = table(phs.aim2.data[a2_combo,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo,"ai"])[2:1,1]
sum.trt.event.tab[5,c(21,24)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_combo),"ai"])[2:1]
sum.trt.event.tab[5,25:26] = table(phs.aim2.data[a2_combo,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo,"tam"])[2:1,2]
sum.trt.event.tab[5,28:29] = table(phs.aim2.data[a2_combo,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo,"tam"])[2:1,1]
sum.trt.event.tab[5,c(27,30)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_combo),"tam"])[2:1]
sum.trt.event.tab[5,31:32] = table(phs.aim2.data[a2_combo,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo,"rad"])[2:1,2]
sum.trt.event.tab[5,34:35] = table(phs.aim2.data[a2_combo,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo,"rad"])[2:1,1]
sum.trt.event.tab[5,c(33,36)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_combo),"rad"])[2:1]
sum.trt.event.tab[5,37:38] = table(phs.aim2.data[a2_combo_l,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo_l,"rad"])[2:1,2]
sum.trt.event.tab[5,40:41] = table(phs.aim2.data[a2_combo_l,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo_l,"rad"])[2:1,1]
sum.trt.event.tab[5,c(39,42)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_combo_l),"rad"])[2:1]


write.csv(sum.trt.event.tab,"Table2b_all.csv")


# LVEF > 50
sum.trt.event.tab = matrix(NA, 5,12)
colnames(sum.trt.event.tab) = c(rep("anathra_y",3),rep("anathra_n",3),rep("tras_y",3),rep("tras_n",3))
rownames(sum.trt.event.tab) = c("ihd","hf","card","stroke","any")

sum.trt.event.tab[1,1:2] = table(phs.aim2.data[a2_ihd_ef,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd_ef,"anthra"])[2:1,2]
sum.trt.event.tab[1,4:5] = table(phs.aim2.data[a2_ihd_ef,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd_ef,"anthra"])[2:1,1]
sum.trt.event.tab[1,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_ihd_ef) & !is.na(phs.aim2.data$final_EF) & phs.aim2.data$final_EF > 50,"anthra"])[2:1]
sum.trt.event.tab[1,7:8] = table(phs.aim2.data[a2_ihd_ef,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd_ef,"tras"])[2:1,2]
sum.trt.event.tab[1,10:11] = table(phs.aim2.data[a2_ihd_ef,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd_ef,"tras"])[2:1,1]
sum.trt.event.tab[1,c(9,12)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_ihd) & !is.na(phs.aim2.data$final_EF) & phs.aim2.data$final_EF > 50,"tras"])[2:1]

sum.trt.event.tab[2,1:2] = table(phs.aim2.data[a2_hf_ef,"heart_failure_grp_inc"], phs.aim2.data[a2_hf_ef,"anthra"])[2:1,2]
sum.trt.event.tab[2,4:5] = table(phs.aim2.data[a2_hf_ef,"heart_failure_grp_inc"], phs.aim2.data[a2_hf_ef,"anthra"])[2:1,1]
sum.trt.event.tab[2,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_hf_ef) & !is.na(phs.aim2.data$final_EF) & phs.aim2.data$final_EF > 50,"anthra"])[2:1]
sum.trt.event.tab[2,7:8] = table(phs.aim2.data[a2_hf_ef,"heart_failure_grp_inc"], phs.aim2.data[a2_hf_ef,"tras"])[2:1,2]
sum.trt.event.tab[2,10:11] = table(phs.aim2.data[a2_hf_ef,"heart_failure_grp_inc"], phs.aim2.data[a2_hf_ef,"tras"])[2:1,1]
sum.trt.event.tab[2,c(9,12)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_hf_ef) & !is.na(phs.aim2.data$final_EF) & phs.aim2.data$final_EF > 50,"tras"])[2:1]

sum.trt.event.tab[3,1:2] = table(phs.aim2.data[a2_cm_ef,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm_ef,"anthra"])[2:1,2]
sum.trt.event.tab[3,4:5] = table(phs.aim2.data[a2_cm_ef,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm_ef,"anthra"])[2:1,1]
sum.trt.event.tab[3,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_cm_ef) & !is.na(phs.aim2.data$final_EF) & phs.aim2.data$final_EF > 50,"anthra"])[2:1]
sum.trt.event.tab[3,7:8] = table(phs.aim2.data[a2_cm_ef,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm_ef,"tras"])[2:1,2]
sum.trt.event.tab[3,10:11] = table(phs.aim2.data[a2_cm_ef,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm_ef,"tras"])[2:1,1]
sum.trt.event.tab[3,c(9,12)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_cm_ef) & !is.na(phs.aim2.data$final_EF) & phs.aim2.data$final_EF > 50,"tras"])[2:1]

sum.trt.event.tab[4,1:2] = table(phs.aim2.data[a2_stroke_ef,"stroke_grp_inc"], phs.aim2.data[a2_stroke_ef,"anthra"])[2:1,2]
sum.trt.event.tab[4,4:5] = table(phs.aim2.data[a2_stroke_ef,"stroke_grp_inc"], phs.aim2.data[a2_stroke_ef,"anthra"])[2:1,1]
sum.trt.event.tab[4,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_stroke_ef) & !is.na(phs.aim2.data$final_EF) & phs.aim2.data$final_EF > 50,"anthra"])[2:1]
sum.trt.event.tab[4,7:8] = table(phs.aim2.data[a2_stroke_ef,"stroke_grp_inc"], phs.aim2.data[a2_stroke_ef,"tras"])[2:1,2]
sum.trt.event.tab[4,10:11] = table(phs.aim2.data[a2_stroke_ef,"stroke_grp_inc"], phs.aim2.data[a2_stroke_ef,"tras"])[2:1,1]
sum.trt.event.tab[4,c(9,12)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_stroke_ef) & !is.na(phs.aim2.data$final_EF) & phs.aim2.data$final_EF > 50,"tras"])[2:1]

sum.trt.event.tab[5,1:2] = table(phs.aim2.data[a2_combo_ef,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo_ef,"anthra"])[2:1,2]
sum.trt.event.tab[5,4:5] = table(phs.aim2.data[a2_combo_ef,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo_ef,"anthra"])[2:1,1]
sum.trt.event.tab[5,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_combo_ef) & !is.na(phs.aim2.data$final_EF) & phs.aim2.data$final_EF > 50,"anthra"])[2:1]
sum.trt.event.tab[5,7:8] = table(phs.aim2.data[a2_combo_ef,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo_ef,"tras"])[2:1,2]
sum.trt.event.tab[5,10:11] = table(phs.aim2.data[a2_combo_ef,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo_ef,"tras"])[2:1,1]
sum.trt.event.tab[5,c(9,12)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_combo_ef) & !is.na(phs.aim2.data$final_EF) & phs.aim2.data$final_EF > 50,"tras"])[2:1]

write.csv(sum.trt.event.tab,"Table2b_LVEFgreater50.csv")


# post menopausal
sum.trt.event.tab = matrix(NA, 5,6)
colnames(sum.trt.event.tab) = c(rep("ai_y",3),rep("ai_n",3))
rownames(sum.trt.event.tab) = c("ihd","hf","card","stroke","any")

sum.trt.event.tab[1,1:2] = table(phs.aim2.data[a2_ihd_post_menop,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd_post_menop,"ai"])[2:1,2]
sum.trt.event.tab[1,4:5] = table(phs.aim2.data[a2_ihd_post_menop,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd_post_menop,"ai"])[2:1,1]
sum.trt.event.tab[1,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_ihd_post_menop) & phs.aim2.data$menop == 1,"ai"])[2:1]

sum.trt.event.tab[2,1:2] = table(phs.aim2.data[a2_hf_post_menop,"heart_failure_grp_inc"], phs.aim2.data[a2_hf_post_menop,"ai"])[2:1,2]
sum.trt.event.tab[2,4:5] = table(phs.aim2.data[a2_hf_post_menop,"heart_failure_grp_inc"], phs.aim2.data[a2_hf_post_menop,"ai"])[2:1,1]
sum.trt.event.tab[2,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_hf_post_menop) & phs.aim2.data$menop == 1,"ai"])[2:1]

sum.trt.event.tab[3,1:2] = table(phs.aim2.data[a2_cm_post_menop,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm_post_menop,"ai"])[2:1,2]
sum.trt.event.tab[3,4:5] = table(phs.aim2.data[a2_cm_post_menop,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm_post_menop,"ai"])[2:1,1]
sum.trt.event.tab[3,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_cm_post_menop) & phs.aim2.data$menop == 1,"ai"])[2:1]

sum.trt.event.tab[4,1:2] = table(phs.aim2.data[a2_stroke_post_menop,"stroke_grp_inc"], phs.aim2.data[a2_stroke_post_menop,"ai"])[2:1,2]
sum.trt.event.tab[4,4:5] = table(phs.aim2.data[a2_stroke_post_menop,"stroke_grp_inc"], phs.aim2.data[a2_stroke_post_menop,"ai"])[2:1,1]
sum.trt.event.tab[4,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_stroke_post_menop) & phs.aim2.data$menop == 1,"ai"])[2:1]

sum.trt.event.tab[5,1:2] = table(phs.aim2.data[a2_combo_post_menop,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo_post_menop,"ai"])[2:1,2]
sum.trt.event.tab[5,4:5] = table(phs.aim2.data[a2_combo_post_menop,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo_post_menop,"ai"])[2:1,1]
sum.trt.event.tab[5,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_combo_post_menop) & phs.aim2.data$menop == 1,"ai"])[2:1]

write.csv(sum.trt.event.tab,"Table2b_PostMenop.csv")

# pre menopausal
sum.trt.event.tab = matrix(NA, 5,6)
colnames(sum.trt.event.tab) = c(rep("ai_y",3),rep("ai_n",3))
rownames(sum.trt.event.tab) = c("ihd","hf","card","stroke","any")

sum.trt.event.tab[1,1:2] = table(phs.aim2.data[a2_ihd_pre_menop,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd_pre_menop,"ai"])[2:1,2]
sum.trt.event.tab[1,4:5] = table(phs.aim2.data[a2_ihd_pre_menop,"ischemic_heart_disease_grp_inc"], phs.aim2.data[a2_ihd_pre_menop,"ai"])[2:1,1]
sum.trt.event.tab[1,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_ihd_pre_menop) & phs.aim2.data$menop == 0,"ai"])[2:1]

sum.trt.event.tab[2,1:2] = table(phs.aim2.data[a2_hf_pre_menop,"heart_failure_grp_inc"], phs.aim2.data[a2_hf_pre_menop,"ai"])[2:1,2]
sum.trt.event.tab[2,4:5] = table(phs.aim2.data[a2_hf_pre_menop,"heart_failure_grp_inc"], phs.aim2.data[a2_hf_pre_menop,"ai"])[2:1,1]
sum.trt.event.tab[2,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_hf_pre_menop) & phs.aim2.data$menop == 0,"ai"])[2:1]

sum.trt.event.tab[3,1:2] = table(phs.aim2.data[a2_cm_pre_menop,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm_pre_menop,"ai"])[2:1,2]
sum.trt.event.tab[3,4:5] = table(phs.aim2.data[a2_cm_pre_menop,"cardiomyopathy_grp_inc"], phs.aim2.data[a2_cm_pre_menop,"ai"])[2:1,1]
sum.trt.event.tab[3,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_cm_pre_menop) & phs.aim2.data$menop == 0,"ai"])[2:1]

sum.trt.event.tab[4,1:2] = table(phs.aim2.data[a2_stroke_pre_menop,"stroke_grp_inc"], phs.aim2.data[a2_stroke_pre_menop,"ai"])[2:1,2]
sum.trt.event.tab[4,4:5] = table(phs.aim2.data[a2_stroke_pre_menop,"stroke_grp_inc"], phs.aim2.data[a2_stroke_pre_menop,"ai"])[2:1,1]
sum.trt.event.tab[4,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_stroke_pre_menop) & phs.aim2.data$menop == 0,"ai"])[2:1]

sum.trt.event.tab[5,1:2] = table(phs.aim2.data[a2_combo_pre_menop,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo_pre_menop,"ai"])[2:1,2]
sum.trt.event.tab[5,4:5] = table(phs.aim2.data[a2_combo_pre_menop,"cvdcombo_grp_inc"], phs.aim2.data[a2_combo_pre_menop,"ai"])[2:1,1]
sum.trt.event.tab[5,c(3,6)] =  table(phs.aim2.data[!(1:nrow(phs.aim2.data) %in% a2_combo_pre_menop) & phs.aim2.data$menop == 0,"ai"])[2:1]

write.csv(sum.trt.event.tab,"Table2b_PreMenop.csv")

#### Chemotherapy
### Anthracycline: (yes/no) 
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(phs.aim2.data[a2_ihd_ef,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_ihd_ef,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ anthra, data = phs.aim2.data[a2_ihd_ef,])
splots[[1]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_ihd_ef,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(phs.aim2.data[a2_stroke_ef,]$stroke_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_stroke_ef,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ anthra, data = phs.aim2.data[a2_stroke_ef,])
splots[[2]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_stroke_ef,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(phs.aim2.data[a2_hf_ef,]$heart_failure_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_hf_ef,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ anthra, data = phs.aim2.data[a2_hf_ef,])
splots[[3]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_hf_ef,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(phs.aim2.data[a2_cm_ef,]$cardiomyopathy_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_cm_ef,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ anthra, data = phs.aim2.data[a2_cm_ef,])
splots[[4]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_cm_ef,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_ef,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_ef,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ anthra, data = phs.aim2.data[a2_combo_ef,])
splots[[5]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_ef,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig1 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig1_Anthracycline.png", fig1,width = 13, height = 18)


### Trastuzumab: (yes/no) (LVEF > 50)
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(phs.aim2.data[a2_ihd_ef,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_ihd_ef,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ tras, data = phs.aim2.data[a2_ihd_ef,])
splots[[1]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_ihd_ef,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(phs.aim2.data[a2_stroke_ef,]$stroke_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_stroke_ef,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ tras, data = phs.aim2.data[a2_stroke_ef,])
splots[[2]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_stroke_ef,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(phs.aim2.data[a2_hf_ef,]$heart_failure_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_hf_ef,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ tras, data = phs.aim2.data[a2_hf_ef,])
splots[[3]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_hf_ef,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(phs.aim2.data[a2_cm_ef,]$cardiomyopathy_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_cm_ef,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ tras, data = phs.aim2.data[a2_cm_ef,])
splots[[4]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_cm_ef,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_ef,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_ef,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ tras, data = phs.aim2.data[a2_combo_ef,])
splots[[5]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_ef,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig2 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig2_Trastuzumab.png", fig2,width = 13, height = 18)

### Taxane: (yes/no)
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(phs.aim2.data[a2_ihd,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_ihd,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ taxane, data = phs.aim2.data[a2_ihd,])
splots[[1]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_ihd,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(phs.aim2.data[a2_stroke,]$stroke_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_stroke,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ taxane, data = phs.aim2.data[a2_stroke,])
splots[[2]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_stroke,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(phs.aim2.data[a2_hf,]$heart_failure_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_hf,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ taxane, data = phs.aim2.data[a2_hf,])
splots[[3]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_hf,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(phs.aim2.data[a2_cm,]$cardiomyopathy_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_cm,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ taxane, data = phs.aim2.data[a2_cm,])
splots[[4]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_cm,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ taxane, data = phs.aim2.data[a2_combo,])
splots[[5]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig3 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig3_Taxane.png", fig3,width = 13, height = 18)


#### Hormonal Therapy
### Aromatase inhibitor (AI) (post menopausal women)
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(phs.aim2.data[a2_ihd_post_menop,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_ihd_post_menop,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ ai, data = phs.aim2.data[a2_ihd_post_menop,])
splots[[1]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_ihd_post_menop,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(phs.aim2.data[a2_stroke_post_menop,]$stroke_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_stroke_post_menop,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ ai, data = phs.aim2.data[a2_stroke_post_menop,])
splots[[2]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_stroke_post_menop,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(phs.aim2.data[a2_hf_post_menop,]$heart_failure_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_hf_post_menop,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ ai, data = phs.aim2.data[a2_hf_post_menop,])
splots[[3]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_hf_post_menop,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(phs.aim2.data[a2_cm_post_menop,]$cardiomyopathy_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_cm_post_menop,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ ai, data = phs.aim2.data[a2_cm_post_menop,])
splots[[4]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_cm_post_menop,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_post_menop,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_post_menop,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ ai, data = phs.aim2.data[a2_combo_post_menop,])
splots[[5]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_post_menop,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig4 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig4_Hormonal_Therapy_AI.png", fig4,width = 13, height = 18)


### Tamoxifen (TAM)
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(phs.aim2.data[a2_ihd_pre_menop,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_ihd_pre_menop,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ tam, data = phs.aim2.data[a2_ihd_pre_menop,])
splots[[1]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_ihd_pre_menop,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(phs.aim2.data[a2_stroke_pre_menop,]$stroke_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_stroke_pre_menop,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ tam, data = phs.aim2.data[a2_stroke_pre_menop,])
splots[[2]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_stroke_pre_menop,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(phs.aim2.data[a2_hf_pre_menop,]$heart_failure_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_hf_pre_menop,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ tam, data = phs.aim2.data[a2_hf_pre_menop,])
splots[[3]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_hf_pre_menop,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(phs.aim2.data[a2_cm_pre_menop,]$cardiomyopathy_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_cm_pre_menop,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ tam, data = phs.aim2.data[a2_cm_pre_menop,])
splots[[4]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_cm_pre_menop,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_pre_menop,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_pre_menop,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ tam, data = phs.aim2.data[a2_combo_pre_menop,])
splots[[5]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_pre_menop,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig5 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig5_Hormonal_Therapy_TAM.png", fig5,width = 13, height = 18)


#### Radiation Therapy
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(phs.aim2.data[a2_ihd,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_ihd,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ rad, data = phs.aim2.data[a2_ihd,])
splots[[1]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_ihd,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(phs.aim2.data[a2_stroke,]$stroke_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_stroke,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ rad, data = phs.aim2.data[a2_stroke,])
splots[[2]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_stroke,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(phs.aim2.data[a2_hf,]$heart_failure_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_hf,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ rad, data = phs.aim2.data[a2_hf,])
splots[[3]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_hf,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(phs.aim2.data[a2_cm,]$cardiomyopathy_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_cm,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ rad, data = phs.aim2.data[a2_cm,])
splots[[4]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_cm,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ rad, data = phs.aim2.data[a2_combo,])
splots[[5]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig6 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig6_Radiation_Therapy.png", fig6,width = 13, height = 18)

#### Left Sided Radiation Therapy
splots <- list()
# ischemic heart disease
surv_object <- Surv(time = floor(phs.aim2.data[a2_ihd_l,]$ischemic_heart_disease_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_ihd_l,]$ischemic_heart_disease_grp_inc)
fit1 <- survfit(surv_object ~ rad, data = phs.aim2.data[a2_ihd_l,])
splots[[1]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_ihd_l,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Ischemic heart disease',xlab='Time/months')

# stroke
surv_object <- Surv(time = floor(phs.aim2.data[a2_stroke_l,]$stroke_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_stroke_l,]$stroke_grp_inc)
fit1 <- survfit(surv_object ~ rad, data = phs.aim2.data[a2_stroke_l,])
splots[[2]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_stroke_l,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Stroke',xlab='Time/months')

# heart failure
surv_object <- Surv(time = floor(phs.aim2.data[a2_hf_l,]$heart_failure_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_hf_l,]$heart_failure_grp_inc)
fit1 <- survfit(surv_object ~ rad, data = phs.aim2.data[a2_hf_l,])
splots[[3]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_hf_l,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Heart failure',xlab='Time/months')

# cardiomyopathy
surv_object <- Surv(time = floor(phs.aim2.data[a2_cm_l,]$cardiomyopathy_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_cm_l,]$cardiomyopathy_grp_inc)
fit1 <- survfit(surv_object ~ rad, data = phs.aim2.data[a2_cm_l,])
splots[[4]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_cm_l,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Cardiomyopathy',xlab='Time/months')

# 4 cvd combined
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_l,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_l,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ rad, data = phs.aim2.data[a2_combo_l,])
splots[[5]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_l,],fun = "event", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes',xlab='Time/months')

fig7 <- arrange_ggsurvplots(splots,ncol = 2, nrow = 3, risk.table.height = 0.25)
ggsave("fig7_Radiation_Therapy_Left_Sided.png", fig7,width = 13, height = 18)




## Plots for CVD combined
splots <- list()

### Anthracycline: (yes/no) LVEF > 50
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_ef,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_ef,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ anthra, data = phs.aim2.data[a2_combo_ef,])
splots[[1]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_ef,],fun = "pct", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes (Anthracycline: LVEF > 50)',xlab='Time/months',
                          legend.labs=c("Not Received Treatment ","Received Treatment"),
                          ylim = c(75,100), pval.coord = c(0, 80))


### Trastuzumab: (yes/no) LVEF > 50
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_ef,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_ef,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ tras, data = phs.aim2.data[a2_combo_ef,])
splots[[2]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_ef,],fun = "pct", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes (Trastuzumab: LVEF > 50)',xlab='Time/months',
                          legend.labs=c("Not Received Treatment ","Received Treatment"),
                          ylim = c(75,100), pval.coord = c(0, 80))


### Taxane: (yes/no)
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ taxane, data = phs.aim2.data[a2_combo,])
splots[[3]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo,],fun = "pct", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes (Taxane)',xlab='Time/months',
                          legend.labs=c("Note Received Treatment ","Received Treatment"),
                          ylim = c(75,100), pval.coord = c(0, 80))


### Aromatase inhibitor (AI) (post menopausal women)
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_post_menop,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_post_menop,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ ai, data = phs.aim2.data[a2_combo_post_menop,])
splots[[4]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_post_menop,],fun = "pct", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes (Aromatase inhibitor: post menopausal)',
                          xlab='Time/months',legend.labs=c("Not Received Treatment ","Received Treatment"),
                          ylim = c(75,100), pval.coord = c(0, 80))



### Tamoxifen (TAM)
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_pre_menop,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_pre_menop,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ tam, data = phs.aim2.data[a2_combo_pre_menop,])
splots[[5]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_pre_menop,],fun = "pct", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes (Tamoxifen: pre menopausal)',xlab='Time/months',
                          legend.labs=c("Not Received Treatment ","Received Treatment"),
                          ylim = c(75,100), pval.coord = c(0, 80))

#### Radiation Therapy
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ rad, data = phs.aim2.data[a2_combo,])
splots[[6]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo,],fun = "pct", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes (Radiation Therapy)',xlab='Time/months',
                          legend.labs=c("Not Received Treatment ","Received Treatment"),
                          ylim = c(75,100), pval.coord = c(0, 80))


#### Left Sided Radiation Therapy
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_l,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_l,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ rad, data = phs.aim2.data[a2_combo_l,])
splots[[7]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_l,],fun = "pct", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Any of the 4 outcomes (Left Sided Radiation Therapy)',xlab='Time/months',
                          legend.labs=c("Not Received Treatment ","Received Treatment"),
                          ylim = c(75,100), pval.coord = c(0, 80))


fig_cvdcombo <- arrange_ggsurvplots(splots,ncol = 2, nrow = 4, risk.table.height = 0.25)
ggsave("fig_cvdcombo_surv.png", fig_cvdcombo,width = 13, height = 18)


# slide 1: for retreat
## Plots for CVD combined
splots <- list()

### Anthracycline: (yes/no) LVEF > 50
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_ef,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_ef,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ anthra, data = phs.aim2.data[a2_combo_ef,])
splots[[1]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_ef,],fun = "pct", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Anthracycline: LVEF > 50',xlab='Time/months',
                          legend.labs=c("Not Received Treatment  ","Received Treatment"),
                          ylim = c(75,100), pval.coord = c(0, 80),legend.title="")


### Taxane: (yes/no)
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ taxane, data = phs.aim2.data[a2_combo,])
splots[[2]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo,],fun = "pct", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Taxane',xlab='Time/months',
                          legend.labs=c("Not Received Treatment  ","Received Treatment"),
                          ylim = c(75,100), pval.coord = c(0, 80),legend.title="")


### Trastuzumab: (yes/no) LVEF > 50
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_ef,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_ef,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ tras, data = phs.aim2.data[a2_combo_ef,])
splots[[3]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_ef,],fun = "pct", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Trastuzumab: LVEF > 50',xlab='Time/months',
                          legend.labs=c("Not Received Treatment  ","Received Treatment"),
                          ylim = c(75,100), pval.coord = c(0, 80),legend.title="")


#### Left Sided Radiation Therapy
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_l,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_l,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ rad, data = phs.aim2.data[a2_combo_l,])
splots[[4]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_l,],fun = "pct", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Left-sided Radiation Therapy',xlab='Time/months',
                          legend.labs=c("Not Received Treatment  ","Received Treatment"),
                          ylim = c(75,100), pval.coord = c(0, 80),legend.title="")


fig_cvdcombo <- arrange_ggsurvplots(splots,ncol = 2, nrow = 2, risk.table.height = 0.25)
ggsave("fig_slide1.png", fig_cvdcombo,width = 16, height = 14)

# slide 2: for retreat
## Plots for CVD combined
splots <- list()

### Aromatase inhibitor (AI) (post menopausal women)
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_post_menop,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_post_menop,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ ai, data = phs.aim2.data[a2_combo_post_menop,])
splots[[2]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_post_menop,],fun = "pct", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Aromatase inhibitor: postmenopausal women',
                          xlab='Time/months',legend.labs=c("Not Received Treatment  ","Received Treatment"),
                          ylim = c(75,100), pval.coord = c(0, 80), legend.title="")



### Tamoxifen (TAM)
surv_object <- Surv(time = floor(phs.aim2.data[a2_combo_pre_menop,]$cvdcombo_grp_inc_fu/30), 
                    event = phs.aim2.data[a2_combo_pre_menop,]$cvdcombo_grp_inc)
fit1 <- survfit(surv_object ~ tam, data = phs.aim2.data[a2_combo_pre_menop,])
splots[[1]] <- ggsurvplot(fit1, data = phs.aim2.data[a2_combo_pre_menop,],fun = "pct", size=1, risk.table = F,pval = T,
                          censor=F,conf.int = T,palette = 'jco',title='Tamoxifen: premenopausal women',xlab='Time/months',
                          legend.labs=c("Not Received Treatment  ","Received Treatment"),
                          ylim = c(75,100), pval.coord = c(0, 80),legend.title="")

fig_cvdcombo <- arrange_ggsurvplots(splots,ncol = 2, nrow = 1, risk.table.height = 0.25)
ggsave("fig_slide2.png", fig_cvdcombo,width = 18, height = 8)


############################################ Cox Regression ############################3
# function to summary cox model results
coxtab <- function(d,x,covar,trt){
  covars <- paste(covar, collapse = '+')
  surv_object <- Surv(time = floor(d[,paste0(x,'_fu')]/30), event = d[,x])
  fit.coxph1 <- eval(parse(text=(paste0("coxph(surv_object ~ ", trt, ", data = d)"))))
  fit.coxph2 <- eval(parse(text=(paste0("coxph(surv_object ~ ", trt, " + ajcc_stage + diab_bl +  
                                        htn_bl + dyslipid_bl + bmicat1 + smok + charlson + 
                                       agegrp + race + edu_cat + income_cat + ",covars,", data = d)"))))
  sum1 <- c(summary(fit.coxph1)$coef[c(2,5)],summary(fit.coxph1)$conf.int[3:4])
  sum2 <- cbind(summary(fit.coxph2)$coef[,c(2,5)],summary(fit.coxph2)$conf.int[,3:4])
  sum3 <- rbind(sum1, sum2)[,c(1,3,4,2)]
  sum4 <- data.frame(cbind(var=row.names(sum3), sum3))
  sum4$var <- factor(sum4$var, levels=row.names(sum3))
  sum4
}

# list of prevalent cvd and treatment
prevcvd <- c("prevcvd","chemo_yn","horm_yn","rad_tx_yn", "menop", "lvef_ind", "I(lvef*lvef_ind)")

phs.aim2.data$lvef_ind = as.numeric(as.character(phs.aim2.data$lvef_ind))

### Table 3a1
# Chemotherapy:  Anthracycline: (yes/no) 
cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6)],"anthra")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[c(-2,-6)],"anthra")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6)],"anthra")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6)],"anthra")
cox_htn <- coxtab(phs.aim2.data[a2_htn,],'cvdrf_htn',prevcvd[c(-2,-6)],"anthra")
cox_diab <- coxtab(phs.aim2.data[a2_diab,],'cvdrf_diab',prevcvd[c(-2,-6)],"anthra")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid,],'cvdrf_dyslipid',prevcvd[c(-2,-6)],"anthra")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo,],'cvdrfcombo',prevcvd[c(-2,-6)],"anthra")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,], cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))


# Chemotherapy:  Trastuzumab: (yes/no)
cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6)],"tras")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[c(-2,-6)],"tras")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6)],"tras")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6)],"tras")
cox_htn <- coxtab(phs.aim2.data[a2_htn,],'cvdrf_htn',prevcvd[c(-2,-6)],"tras")
cox_diab <- coxtab(phs.aim2.data[a2_diab,],'cvdrf_diab',prevcvd[c(-2,-6)],"tras")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid,],'cvdrf_dyslipid',prevcvd[c(-2,-6)],"tras")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo,],'cvdrfcombo',prevcvd[c(-2,-6)],"tras")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,], cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))


# Chemotherapy:  Taxane: (yes/no) 	
cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_htn <- coxtab(phs.aim2.data[a2_htn,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"taxane")
cox_diab <- coxtab(phs.aim2.data[a2_diab,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"taxane")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"taxane")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"taxane")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,], cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))


# Chemotherapy:  Cyclophasphamide: (yes/no) 
phs.aim2.data$cyclo = as.factor(phs.aim2.data$cyclo_yn)

cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_htn <- coxtab(phs.aim2.data[a2_htn,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"cyclo")
cox_diab <- coxtab(phs.aim2.data[a2_diab,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"cyclo")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"cyclo")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"cyclo")
cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,], cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))

# Chemotherapy:  Doxorubicin: (yes/no) 
phs.aim2.data$doxo = as.factor(phs.aim2.data$doxorub_yn)

cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6)],"doxo")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[c(-2,-6)],"doxo")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6)],"doxo")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6)],"doxo")
cox_htn <- coxtab(phs.aim2.data[a2_htn,],'cvdrf_htn',prevcvd[c(-2,-6)],"doxo")
cox_diab <- coxtab(phs.aim2.data[a2_diab,],'cvdrf_diab',prevcvd[c(-2,-6)],"doxo")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid,],'cvdrf_dyslipid',prevcvd[c(-2,-6)],"doxo")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo,],'cvdrfcombo',prevcvd[c(-2,-6)],"doxo")
cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))

# Chemotherapy:  Epirubicin: (yes/no) 
phs.aim2.data$epi = as.factor(phs.aim2.data$epirub_yn)

cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6)],"epi")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[c(-2,-6)],"epi")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6)],"epi")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6)],"epi")
cox_htn <- coxtab(phs.aim2.data[a2_htn,],'cvdrf_htn',prevcvd[c(-2,-6)],"epi")
cox_diab <- coxtab(phs.aim2.data[a2_diab,],'cvdrf_diab',prevcvd[c(-2,-6)],"epi")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid,],'cvdrf_dyslipid',prevcvd[c(-2,-6)],"epi")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo,],'cvdrfcombo',prevcvd[c(-2,-6)],"epi")
cox6 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox6[,-1] <- lapply(cox6[,-1], function(x) as.numeric(as.character(x)))
cox6 <- cbind(round(cox6[,2],2),paste0("(",round(cox6[,3],2),", ", round(cox6[,4],2),")"))

# Chemotherapy:  Fluoropyrimidine: (yes/no) 
phs.aim2.data$flu = as.factor(phs.aim2.data$fluoro_yn)

cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"flu")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"flu")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"flu")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"flu")
cox_htn <- coxtab(phs.aim2.data[a2_htn,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"flu")
cox_diab <- coxtab(phs.aim2.data[a2_diab,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"flu")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"flu")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"flu")
cox7 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox7[,-1] <- lapply(cox7[,-1], function(x) as.numeric(as.character(x)))
cox7 <- cbind(round(cox7[,2],2),paste0("(",round(cox7[,3],2),", ", round(cox7[,4],2),")"))



cox_all <- data.frame(cbind(cox1, '', cox2,'', cox3, '', cox4, '', cox5, '', cox6, '', cox7))
write.csv(cox_all,'T3a1_cox_all_chemo.csv')

### Table 3a2, excluding 110 participants with prevalent CVD at baseline
# Chemotherapy:  Anthracycline: (yes/no) 
cox_isch <- coxtab(phs.aim2.data[a2_ihd_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-2,-6)],"anthra")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_cvd,],'stroke_grp_inc',prevcvd[c(-1,-2,-6)],"anthra")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-2,-6)],"anthra")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1,-2,-6)],"anthra")
cox_htn <- coxtab(phs.aim2.data[a2_htn_cvd,],'cvdrf_htn',prevcvd[c(-1,-2,-6)],"anthra")
cox_diab <- coxtab(phs.aim2.data[a2_diab_cvd,],'cvdrf_diab',prevcvd[c(-1,-2,-6)],"anthra")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-2,-6)],"anthra")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_cvd,],'cvdrfcombo',prevcvd[c(-1,-2,-6)],"anthra")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,], cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))


# Chemotherapy:  Trastuzumab: (yes/no)
cox_isch <- coxtab(phs.aim2.data[a2_ihd_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-2,-6)],"tras")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_cvd,],'stroke_grp_inc',prevcvd[c(-1,-2,-6)],"tras")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-2,-6)],"tras")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1,-2,-6)],"tras")
cox_htn <- coxtab(phs.aim2.data[a2_htn_cvd,],'cvdrf_htn',prevcvd[c(-1,-2,-6)],"tras")
cox_diab <- coxtab(phs.aim2.data[a2_diab_cvd,],'cvdrf_diab',prevcvd[c(-1,-2,-6)],"tras")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-2,-6)],"tras")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_cvd,],'cvdrfcombo',prevcvd[c(-1,-2,-6)],"tras")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))


# Chemotherapy:  Taxane: (yes/no) 	
cox_isch <- coxtab(phs.aim2.data[a2_ihd_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_cvd,],'stroke_grp_inc',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_htn <- coxtab(phs.aim2.data[a2_htn_cvd,],'cvdrf_htn',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_diab <- coxtab(phs.aim2.data[a2_diab_cvd,],'cvdrf_diab',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_cvd,],'cvdrfcombo',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))


# Chemotherapy:  Cyclophasphamide: (yes/no) 
phs.aim2.data$cyclo = as.factor(phs.aim2.data$cyclo_yn)

cox_isch <- coxtab(phs.aim2.data[a2_ihd_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1, -2,-6,-7)],"cyclo")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_cvd,],'stroke_grp_inc',prevcvd[c(-1, -2,-6,-7)],"cyclo")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1, -2,-6,-7)],"cyclo")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1, -2,-6,-7)],"cyclo")
cox_htn <- coxtab(phs.aim2.data[a2_htn_cvd,],'cvdrf_htn',prevcvd[c(-1,-2,-6,-7)],"cyclo")
cox_diab <- coxtab(phs.aim2.data[a2_diab_cvd,],'cvdrf_diab',prevcvd[c(--1,-2,-6,-7)],"cyclo")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-2,-6,-7)],"cyclo")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_cvd,],'cvdrfcombo',prevcvd[c(-1,-2,-6,-7)],"cyclo")
cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))

# Chemotherapy:  Doxorubicin: (yes/no) 
phs.aim2.data$doxo = as.factor(phs.aim2.data$doxorub_yn)

cox_isch <- coxtab(phs.aim2.data[a2_ihd_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1, -2,-6)],"doxo")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_cvd,],'stroke_grp_inc',prevcvd[c(-1, -2,-6)],"doxo")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1, -2,-6)],"doxo")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1, -2,-6)],"doxo")
cox_htn <- coxtab(phs.aim2.data[a2_htn_cvd,],'cvdrf_htn',prevcvd[c(-1,-2,-6)],"doxo")
cox_diab <- coxtab(phs.aim2.data[a2_diab_cvd,],'cvdrf_diab',prevcvd[c(-1,-2,-6)],"doxo")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-2,-6)],"doxo")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_cvd,],'cvdrfcombo',prevcvd[c(-1,-2,-6)],"doxo")
cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))

# Chemotherapy:  Epirubicin: (yes/no) 
phs.aim2.data$epi = as.factor(phs.aim2.data$epirub_yn)

cox_isch <- coxtab(phs.aim2.data[a2_ihd_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1, -2,-6)],"epi")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_cvd,],'stroke_grp_inc',prevcvd[c(-1, -2,-6)],"epi")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1, -2,-6)],"epi")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1, -2,-6)],"epi")
cox_htn <- coxtab(phs.aim2.data[a2_htn_cvd,],'cvdrf_htn',prevcvd[c(-1,-2,-6)],"epi")
cox_diab <- coxtab(phs.aim2.data[a2_diab_cvd,],'cvdrf_diab',prevcvd[c(-1,-2,-6)],"epi")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-2,-6)],"epi")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_cvd,],'cvdrfcombo',prevcvd[c(-1,-2,-6)],"epi")
cox6 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox6[,-1] <- lapply(cox6[,-1], function(x) as.numeric(as.character(x)))
cox6 <- cbind(round(cox6[,2],2),paste0("(",round(cox6[,3],2),", ", round(cox6[,4],2),")"))

# Chemotherapy:  Fluoropyrimidine: (yes/no) 
phs.aim2.data$flu = as.factor(phs.aim2.data$fluoro_yn)

cox_isch <- coxtab(phs.aim2.data[a2_ihd_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1, -2,-6,-7)],"flu")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_cvd,],'stroke_grp_inc',prevcvd[c(-1, -2,-6,-7)],"flu")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1, -2,-6,-7)],"flu")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1, -2,-6,-7)],"flu")
cox_htn <- coxtab(phs.aim2.data[a2_htn_cvd,],'cvdrf_htn',prevcvd[c(-1, -2,-6,-7)],"flu")
cox_diab <- coxtab(phs.aim2.data[a2_diab_cvd,],'cvdrf_diab',prevcvd[c(-1, -2,-6,-7)],"flu")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_cvd,],'cvdrf_dyslipid',prevcvd[c(-1, -2,-6,-7)],"flu")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_cvd,],'cvdrfcombo',prevcvd[c(-1, -2,-6,-7)],"flu")
cox7 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox7[,-1] <- lapply(cox7[,-1], function(x) as.numeric(as.character(x)))
cox7 <- cbind(round(cox7[,2],2),paste0("(",round(cox7[,3],2),", ", round(cox7[,4],2),")"))



cox_all <- data.frame(cbind(cox1, '', cox2,'', cox3, '', cox4, '', cox5, '', cox6, '', cox7))
write.csv(cox_all,'T3a2_cox_all_chemo.csv')



t3a = c(table(phs.aim2.data$anthra)[2:1],table(phs.aim2.data$tras)[2:1],table(phs.aim2.data$taxane)[2:1],
        table(phs.aim2.data$cyclo)[2:1],table(phs.aim2.data$doxorub_yn)[2:1],table(phs.aim2.data$epirub_yn)[2:1],table(phs.aim2.data$fluoro_yn)[2:1])
t3a = data.frame(rbind(t3a))
for (i in 1:dim(t3a)[2]){
  t3a[1,i] = paste0("n = ", t3a[1,i])
}


t3b = c(table(phs.aim2.data.cvd$anthra)[2:1],table(phs.aim2.data.cvd$tras)[2:1],table(phs.aim2.data.cvd$taxane)[2:1],
        table(phs.aim2.data.cvd$cyclo_yn)[2:1],table(phs.aim2.data.cvd$doxorub_yn)[2:1],table(phs.aim2.data.cvd$epirub_yn)[2:1],table(phs.aim2.data.cvd$fluoro_yn)[2:1])
t3b = data.frame(rbind(t3b))
for (i in 1:dim(t3b)[2]){
  t3b[1,i] = paste0("n = ", t3b[1,i])

}

t3 = rbind(t3a, t3b)
write.csv(t3,"t3.csv")


# ### Table 3a3
# # Chemotherapy:  Anthracycline: (yes/no) 
# cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-3,-4,-6)],"anthra")
# cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[c(-2,-3,-4,-6)],"anthra")
# cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-3,-4,-6)],"anthra")
# cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)],"anthra")
# cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
# cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))
# 
# 
# # Chemotherapy:  Trastuzumab: (yes/no)
# cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-3,-4,-6)],"tras")
# cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[c(-2,-3,-4,-6)],"tras")
# cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-3,-4,-6)],"tras")
# cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6)],"tras")
# cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
# cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))
# 
# 
# # Chemotherapy:  Taxane: (yes/no) 	
# cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-3,-4,-6,-7)],"taxane")
# cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[c(-2,-3,-4,-6,-7)],"taxane")
# cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-3,-4,-6,-7)],"taxane")
# cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-3,-4,-6,-7)],"taxane")
# cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
# cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))
# 
# cox_all <- data.frame(cbind(cox1, '', cox2,'', cox3))
# write.csv(cox_all,'T3a3_cox_all_chemo.csv')
# 
# 
# ### Table 3a4
# # Chemotherapy:  Anthracycline: (yes/no) 
# cox_isch <- coxtab(phs.aim2.data[a2_ihd_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-2,-3,-4,-6)],"anthra")
# cox_stroke <- coxtab(phs.aim2.data[a2_stroke_cvd,],'stroke_grp_inc',prevcvd[c(-1,-2,-3,-4,-6)],"anthra")
# cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-2,-3,-4,-6)],"anthra")
# cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1,-2,-3,-4,-6)],"anthra")
# cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
# cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))
# 
# 
# # Chemotherapy:  Trastuzumab: (yes/no)
# cox_isch <- coxtab(phs.aim2.data[a2_ihd_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-2,-3,-4,-6)],"tras")
# cox_stroke <- coxtab(phs.aim2.data[a2_stroke_cvd,],'stroke_grp_inc',prevcvd[c(-1,-2,-3,-4,-6)],"tras")
# cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-2,-3,-4,-6)],"tras")
# cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1,-2,-3,-4,-6)],"tras")
# cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
# cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))
# 
# 
# # Chemotherapy:  Taxane: (yes/no) 	
# cox_isch <- coxtab(phs.aim2.data[a2_ihd_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-2,-3,-4,-6,-7)],"taxane")
# cox_stroke <- coxtab(phs.aim2.data[a2_stroke_cvd,],'stroke_grp_inc',prevcvd[c(-1,-2,-3,-4,-6,-7)],"taxane")
# cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-2,-3,-4,-6,-7)],"taxane")
# cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1,-2,-3,-4,-6,-7)],"taxane")
# cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
# cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))
# 
# cox_all <- data.frame(cbind(cox1, '', cox2,'', cox3))
# write.csv(cox_all,'T3a4_cox_all_chemo.csv')



### Table 3b1, excluding 20 participants with LVEF <= 50 for Anthracycline and Trastuzumab
# Chemotherapy:  Anthracycline: (yes/no) 
cox_isch <- coxtab(phs.aim2.data[a2_ihd_ef,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_ef,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_ef,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[-c(-2,-6,-7)],"anthra")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_ef,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_htn <- coxtab(phs.aim2.data[a2_htn_ef,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"anthra")
cox_diab <- coxtab(phs.aim2.data[a2_diab_ef,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"anthra")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_ef,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"anthra")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_ef,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"anthra")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))


# Chemotherapy:  Trastuzumab: (yes/no)
cox_isch <- coxtab(phs.aim2.data[a2_ihd_ef,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_ef,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_ef,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_ef,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_htn <- coxtab(phs.aim2.data[a2_htn_ef,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"tras")
cox_diab <- coxtab(phs.aim2.data[a2_diab_ef,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"tras")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_ef,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"tras")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_ef,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"tras")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))


# Chemotherapy:  Taxane: (yes/no) 	
cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_htn <- coxtab(phs.aim2.data[a2_htn_ef,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"taxane")
cox_diab <- coxtab(phs.aim2.data[a2_diab_ef,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"taxane")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_ef,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"taxane")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_ef,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"taxane")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))

cox_all <- data.frame(cbind(cox1, '', cox2,'', cox3))
write.csv(cox_all,'T3b1_cox_all_chemo.csv')

### Table 3b2, excluding 20 participants with LVEF <= 50 for Anthracycline and Trastuzumab and 110 participants with prevalent CVD at baseline
# Chemotherapy:  Anthracycline: (yes/no) 
cox_isch <- coxtab(phs.aim2.data[a2_ihd_ef_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-2,-6,-7)],"anthra")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_ef_cvd,],'stroke_grp_inc',prevcvd[c(-1,-2,-6,-7)],"anthra")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_ef_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-2,-6,-7)],"anthra")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_ef_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1,-2,-6,-7)],"anthra")
cox_htn <- coxtab(phs.aim2.data[a2_htn_ef_cvd,],'cvdrf_htn',prevcvd[c(-1,-2,-6,-7)],"anthra")
cox_diab <- coxtab(phs.aim2.data[a2_diab_ef_cvd,],'cvdrf_diab',prevcvd[c(--1,-2,-6,-7)],"anthra")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_ef_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-2,-6,-7)],"anthra")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_ef_cvd,],'cvdrfcombo',prevcvd[c(-1,-2,-6,-7)],"anthra")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))


# Chemotherapy:  Trastuzumab: (yes/no)
cox_isch <- coxtab(phs.aim2.data[a2_ihd_ef_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-2,-6,-7)],"tras")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_ef_cvd,],'stroke_grp_inc',prevcvd[c(-1,-2,-6,-7)],"tras")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_ef_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-2,-6,-7)],"tras")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_ef_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1,-2,-6,-7)],"tras")
cox_htn <- coxtab(phs.aim2.data[a2_htn_ef_cvd,],'cvdrf_htn',prevcvd[c(-1,-2,-6,-7)],"tras")
cox_diab <- coxtab(phs.aim2.data[a2_diab_ef_cvd,],'cvdrf_diab',prevcvd[c(-1,-2,-6,-7)],"tras")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_ef_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-2,-6,-7)],"tras")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_ef_cvd,],'cvdrfcombo',prevcvd[c(-1,-2,-6,-7)],"tras")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))


# Chemotherapy:  Taxane: (yes/no) 	
cox_isch <- coxtab(phs.aim2.data.cvd[a2_ihd_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_stroke <- coxtab(phs.aim2.data.cvd[a2_stroke_cvd,],'stroke_grp_inc',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_hf_cm <- coxtab(phs.aim2.data.cvd[a2_hf_cm_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_cvdcombo <- coxtab(phs.aim2.data.cvd[a2_combo_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_htn <- coxtab(phs.aim2.data[a2_htn_ef_cvd,],'cvdrf_htn',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_diab <- coxtab(phs.aim2.data[a2_diab_ef_cvd,],'cvdrf_diab',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_ef_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_ef_cvd,],'cvdrfcombo',prevcvd[c(-1,-2,-6,-7)],"taxane")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))

cox_all <- data.frame(cbind(cox1, '', cox2,'', cox3))
write.csv(cox_all,'T3b2_cox_all_chemo.csv')



### Table 4a
### Radiation Therapy and Left Sided Radiation Therapy
cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-4,-6,-7)],"rad")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[c(-4,-6,-7)],"rad")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-4,-6,-7)],"rad")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-4,-6,-7)],"rad")
cox_htn <- coxtab(phs.aim2.data[a2_htn,],'cvdrf_htn',prevcvd[c(-4,-6,-7)],"rad")
cox_diab <- coxtab(phs.aim2.data[a2_diab,],'cvdrf_diab',prevcvd[c(-4,-6,-7)],"rad")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid,],'cvdrf_dyslipid',prevcvd[c(-4,-6,-7)],"rad")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo,],'cvdrfcombo',prevcvd[c(-4,-6,-7)],"rad")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))

cox_isch <- coxtab(phs.aim2.data[a2_ihd_l,],'ischemic_heart_disease_grp_inc',prevcvd[c(-4,-6,-7)],"rad")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_l,],'stroke_grp_inc',prevcvd[c(-4,-6,-7)],"rad")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_l,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-4,-6,-7)],"rad")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_l,],'cvdcombo_grp_inc',prevcvd[c(-4,-6,-7)],"rad")
cox_htn <- coxtab(phs.aim2.data[a2_htn_l,],'cvdrf_htn',prevcvd[c(-4,-6,-7)],"rad")
cox_diab <- coxtab(phs.aim2.data[a2_diab_l,],'cvdrf_diab',prevcvd[c(-4,-6,-7)],"rad")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_l,],'cvdrf_dyslipid',prevcvd[c(-4,-6,-7)],"rad")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_l,],'cvdrfcombo',prevcvd[c(-4,-6,-7)],"rad")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))

cox_all <- data.frame(cbind(cox1, '', cox2))
write.csv(cox_all,'T4a_cox_all_rad.csv')

### Table 4b, excluding 110 prevalent cvd at baseline
### Radiation Therapy and Left Sided Radiation Therapy
cox_isch <- coxtab(phs.aim2.data[a2_ihd_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_cvd,],'stroke_grp_inc',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_htn <- coxtab(phs.aim2.data[a2_htn_cvd,],'cvdrf_htn',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_diab <- coxtab(phs.aim2.data[a2_diab_cvd,],'cvdrf_diab',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_cvd,],'cvdrfcombo',prevcvd[c(-1,-4,-6,-7)],"rad")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))

cox_isch <- coxtab(phs.aim2.data[a2_ihd_l_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_l_cvd,],'stroke_grp_inc',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_l_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_l_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_htn <- coxtab(phs.aim2.data[a2_htn_l_cvd,],'cvdrf_htn',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_diab <- coxtab(phs.aim2.data[a2_diab_l_cvd,],'cvdrf_diab',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_l_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-4,-6,-7)],"rad")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_l_cvd,],'cvdrfcombo',prevcvd[c(-1,-4,-6,-7)],"rad")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))

cox_all <- data.frame(cbind(cox1, '', cox2))
write.csv(cox_all,'T4b_cox_all_rad.csv')



### Table 5a
#### Hormonal Therapy
### Aromatase inhibitor (AI) (post menopausal women)
cox_isch <- coxtab(phs.aim2.data[a2_ihd_post_menop,],'ischemic_heart_disease_grp_inc',prevcvd[c(-3,-5,-6,-7)],"ai")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_post_menop,],'stroke_grp_inc',prevcvd[c(-3,-5,-6,-7)],"ai")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_post_menop,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-3,-5,-6,-7)],"ai")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_post_menop,],'cvdcombo_grp_inc',prevcvd[c(-3,-5,-6,-7)],"ai")
cox_htn <- coxtab(phs.aim2.data[a2_htn_post_menop,],'cvdrf_htn',prevcvd[c(-3,-5,-6,-7)],"ai")
cox_diab <- coxtab(phs.aim2.data[a2_diab_post_menop,],'cvdrf_diab',prevcvd[c(-3,-5,-6,-7)],"ai")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_post_menop,],'cvdrf_dyslipid',prevcvd[c(-3,-5,-6,-7)],"ai")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_post_menop,],'cvdrfcombo',prevcvd[c(-3,-5,-6,-7)],"ai")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round2(cox1[,2],2),paste0("(",round2(cox1[,3],2),", ", round2(cox1[,4],2),")"))

### Tamoxifen (TAM) (pre menopausal women)
cox_isch <- coxtab(phs.aim2.data[a2_ihd_pre_menop,],'ischemic_heart_disease_grp_inc',prevcvd[c(-3,-5,-6,-7)],"tam")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_pre_menop,],'stroke_grp_inc',prevcvd[c(-3,-5,-6,-7)],"tam")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_pre_menop,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-3,-5,-6,-7)],"tam")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_pre_menop,],'cvdcombo_grp_inc',prevcvd[c(-3,-5,-6,-7)],"tam")
cox_htn <- coxtab(phs.aim2.data[a2_htn_pre_menop,],'cvdrf_htn',prevcvd[c(-3,-5,-6,-7)],"tam")
cox_diab <- coxtab(phs.aim2.data[a2_diab_pre_menop,],'cvdrf_diab',prevcvd[c(-3,-5,-6,-7)],"tam")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_pre_menop,],'cvdrf_dyslipid',prevcvd[c(-3,-5,-6,-7)],"tam")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_pre_menop,],'cvdrfcombo',prevcvd[c(-3,-5,-6,-7)],"tam")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round2(cox2[,2],2),paste0("(",round2(cox2[,3],2),", ", round2(cox2[,4],2),")"))

### Any hormonal therapy
cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-3,-6,-7)],"horm_any")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[c(-3,-6,-7)],"horm_any")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-3,-6,-7)],"horm_any")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[c(-3,-6,-7)],"horm_any")
cox_htn <- coxtab(phs.aim2.data[a2_htn,],'cvdrf_htn',prevcvd[c(-3,-6,-7)],"horm_any")
cox_diab <- coxtab(phs.aim2.data[a2_diab,],'cvdrf_diab',prevcvd[c(-3,-6,-7)],"horm_any")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid,],'cvdrf_dyslipid',prevcvd[c(-3,-6,-7)],"horm_any")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo,],'cvdrfcombo',prevcvd[c(-3,-6,-7)],"horm_any")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
cox3 <- cbind(round2(cox3[,2],2),paste0("(",round2(cox3[,3],2),", ", round2(cox3[,4],2),")"))


cox_all <- data.frame(cbind(cox1, '', cox2, '', cox3))
write.csv(cox_all,'T5a_cox_all_horm.csv')

### Table 5b, excluding 110 prevalent cvd at baseline
#### Hormonal Therapy
### Aromatase inhibitor (AI) (post menopausal women)
cox_isch <- coxtab(phs.aim2.data[a2_ihd_post_menop_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-3,-5,-6,-7)],"ai")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_post_menop_cvd,],'stroke_grp_inc',prevcvd[c(-1,-3,-5,-6,-7)],"ai")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_post_menop_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-3,-5,-6,-7)],"ai")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_post_menop_cvd,],'cvdcombo_grp_inc',prevcvd[c(-3,-5,-6,-7)],"ai")
cox_htn <- coxtab(phs.aim2.data[a2_htn_post_menop_cvd,],'cvdrf_htn',prevcvd[c(-1,-3,-5,-6,-7)],"ai")
cox_diab <- coxtab(phs.aim2.data[a2_diab_post_menop_cvd,],'cvdrf_diab',prevcvd[c(-1,-3,-5,-6,-7)],"ai")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_post_menop_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-3,-5,-6,-7)],"ai")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_post_menop_cvd,],'cvdrfcombo',prevcvd[c(-1,-3,-5,-6,-7)],"ai")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round2(cox1[,2],2),paste0("(",round2(cox1[,3],2),", ", round2(cox1[,4],2),")"))

### Tamoxifen (TAM) (pre menopausal women)
cox_isch <- coxtab(phs.aim2.data[a2_ihd_pre_menop_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-3,-5,-6,-7)],"tam")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_pre_menop_cvd,],'stroke_grp_inc',prevcvd[c(-1,-3,-5,-6,-7)],"tam")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_pre_menop_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-3,-5,-6,-7)],"tam")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_pre_menop_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1,-3,-5,-6,-7)],"tam")
cox_htn <- coxtab(phs.aim2.data[a2_htn_pre_menop_cvd,],'cvdrf_htn',prevcvd[c(-1,-3,-5,-6,-7)],"tam")
cox_diab <- coxtab(phs.aim2.data[a2_diab_pre_menop_cvd,],'cvdrf_diab',prevcvd[c(-1,-3,-5,-6,-7)],"tam")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_pre_menop_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-3,-5,-6,-7)],"tam")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_pre_menop_cvd,],'cvdrfcombo',prevcvd[c(-1,-3,-5,-6,-7)],"tam")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round2(cox2[,2],2),paste0("(",round2(cox2[,3],2),", ", round2(cox2[,4],2),")"))

### Any hormonal therapy
cox_isch <- coxtab(phs.aim2.data[a2_ihd_cvd,],'ischemic_heart_disease_grp_inc',prevcvd[c(-1,-3,-6,-7)],"horm_any")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_cvd,],'stroke_grp_inc',prevcvd[c(-1,-3,-6,-7)],"horm_any")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_cvd,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-1,-3,-6,-7)],"horm_any")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_cvd,],'cvdcombo_grp_inc',prevcvd[c(-1,-3,-6,-7)],"horm_any")
cox_htn <- coxtab(phs.aim2.data[a2_htn_cvd,],'cvdrf_htn',prevcvd[c(-1,-3,-6,-7)],"horm_any")
cox_diab <- coxtab(phs.aim2.data[a2_diab_cvd,],'cvdrf_diab',prevcvd[c(-1,-3,-6,-7)],"horm_any")
cox_dyslipid <- coxtab(phs.aim2.data[a2_dyslipid_cvd,],'cvdrf_dyslipid',prevcvd[c(-1,-3,-6,-7)],"horm_any")
cox_rfcomb <- coxtab(phs.aim2.data[a2_rfcombo_cvd,],'cvdrfcombo',prevcvd[c(-1,-3,-6,-7)],"horm_any")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
cox3 <- cbind(round2(cox3[,2],2),paste0("(",round2(cox3[,3],2),", ", round2(cox3[,4],2),")"))

cox_all <- data.frame(cbind(cox1, '', cox2, '', cox3))
write.csv(cox_all,'T5b_cox_all_horm.csv')



### Table 6a, excluding 20 participants with LVEF <= 50 for Anthracycline (BMI Category)
coxtab2 <- function(d,x,covar,trt){
  covars <- paste(covar, collapse = '+')
  surv_object <- Surv(time = floor(d[,paste0(x,'_fu')]/30), event = d[,x])
  fit.coxph1 <- eval(parse(text=(paste0("coxph(surv_object ~ ", trt, ", data = d)"))))
  fit.coxph2 <- eval(parse(text=(paste0("coxph(surv_object ~ ", trt, " + ajcc_stage + diab_bl +  
                                        htn_bl + dyslipid_bl + smok + charlson + 
                                       agegrp + race + edu_cat + income_cat + ",covars,", data = d)"))))
  sum1 <- c(summary(fit.coxph1)$coef[c(2,5)],summary(fit.coxph1)$conf.int[3:4])
  sum2 <- cbind(summary(fit.coxph2)$coef[,c(2,5)],summary(fit.coxph2)$conf.int[,3:4])
  sum3 <- rbind(sum1, sum2)[,c(1,3,4,2)]
  sum4 <- data.frame(cbind(var=row.names(sum3), sum3))
  sum4$var <- factor(sum4$var, levels=row.names(sum3))
  sum4
}

# Chemotherapy:  Anthracycline: (yes/no) 
cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_underweight,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_underweight,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_underweight,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_underweight,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_underweight,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"anthra")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_underweight,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"anthra")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_underweight,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"anthra")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_underweight,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"anthra")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_normal,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_normal,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_normal,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_normal,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_normal,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"anthra")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_normal,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"anthra")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_normal,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"anthra")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_normal,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"anthra")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
hr_adj_normal_anthra = cox2[c(2,4,6,8),2]
lhr_adj_normal_anthra = cox2[c(2,4,6,8),3]
uhr_adj_normal_anthra = cox2[c(2,4,6,8),4]
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_overweight,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_overweight,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_overweight,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_overweight,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_overweight,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"anthra")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_overweight,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"anthra")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_overweight,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"anthra")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_overweight,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"anthra")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
hr_adj_over_anthra = cox3[c(2,4,6,8),2]
lhr_adj_over_anthra = cox3[c(2,4,6,8),3]
uhr_adj_over_anthra = cox3[c(2,4,6,8),4]
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_obese1,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_obese1,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_obese1,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_obese1,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_obese1,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"anthra")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_obese1,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"anthra")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_obese1,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"anthra")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_obese1,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"anthra")
cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
hr_adj_ob1_anthra = cox4[c(2,4,6,8),2]
lhr_adj_ob1_anthra = cox4[c(2,4,6,8),3]
uhr_adj_ob1_anthra = cox4[c(2,4,6,8),4]
cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_obese2,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_obese2,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_obese2,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_obese2,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"anthra")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_obese2,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"anthra")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_obese2,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"anthra")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_obese2,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"anthra")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_obese2,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"anthra")
cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
hr_adj_ob2_anthra = cox5[c(2,4,6,8),2]
lhr_adj_ob2_anthra = cox5[c(2,4,6,8),3]
uhr_adj_ob2_anthra = cox5[c(2,4,6,8),4]
cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))


cox_all <- data.frame(cbind(cox1, '', cox2, '', cox3, '', cox4, '', cox5))
write.csv(cox_all,'T6a_cox_all.csv')

### Table 6b, excluding 20 participants with LVEF <= 50 for Trastuzumab (BMI Category)
# Chemotherapy:  Trastuzumab: (yes/no) 
cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_underweight,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_underweight,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_underweight,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_underweight,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_underweight,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"tras")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_underweight,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"tras")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_underweight,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"tras")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_underweight,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"tras")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_normal,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_normal,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_normal,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_normal,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_normal,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"tras")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_normal,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"tras")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_normal,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"tras")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_normal,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"tras")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
hr_adj_normal_tras = cox2[c(2,4,6,8),2]
lhr_adj_normal_tras = cox2[c(2,4,6,8),3]
uhr_adj_normal_tras = cox2[c(2,4,6,8),4]
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_overweight,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_overweight,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_overweight,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_overweight,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_overweight,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"tras")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_overweight,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"tras")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_overweight,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"tras")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_overweight,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"tras")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
hr_adj_over_tras = cox3[c(2,4,6,8),2]
lhr_adj_over_tras = cox3[c(2,4,6,8),3]
uhr_adj_over_tras = cox3[c(2,4,6,8),4]
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_obese1,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_obese1,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_obese1,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_obese1,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_obese1,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"tras")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_obese1,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"tras")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_obese1,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"tras")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_obese1,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"tras")
cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
hr_adj_ob1_tras = cox4[c(2,4,6,8),2]
lhr_adj_ob1_tras = cox4[c(2,4,6,8),3]
uhr_adj_ob1_tras = cox4[c(2,4,6,8),4]
cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_obese2,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_obese2,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_obese2,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_obese2,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"tras")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_obese2,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"tras")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_obese2,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"tras")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_obese2,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"tras")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_obese2,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"tras")
cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
hr_adj_ob2_tras = cox5[c(2,4,6,8),2]
lhr_adj_ob2_tras = cox5[c(2,4,6,8),3]
uhr_adj_ob2_tras = cox5[c(2,4,6,8),4]
cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))


cox_all <- data.frame(cbind(cox1, '', cox2, '', cox3, '', cox4, '', cox5))
write.csv(cox_all,'T6b_cox_all.csv')

### Table 6c, excluding 20 participants with LVEF <= 50 for Taxane (BMI Category)
# Chemotherapy:  Taxane: (yes/no) 
cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_underweight,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_underweight,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_underweight,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_underweight,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_underweight,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"taxane")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_underweight,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"taxane")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_underweight,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"taxane")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_underweight,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"taxane")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_normal,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_normal,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_normal,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_normal,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_normal,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"taxane")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_normal,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"taxane")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_normal,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"taxane")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_normal,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"taxane")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
hr_adj_normal_taxane = cox2[c(2,4,6,8),2]
lhr_adj_normal_taxane = cox2[c(2,4,6,8),3]
uhr_adj_normal_taxane = cox2[c(2,4,6,8),4]
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_overweight,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_overweight,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_overweight,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_overweight,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_overweight,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"taxane")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_overweight,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"taxane")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_overweight,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"taxane")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_overweight,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"taxane")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
hr_adj_over_taxane = cox3[c(2,4,6,8),2]
lhr_adj_over_taxane = cox3[c(2,4,6,8),3]
uhr_adj_over_taxane = cox3[c(2,4,6,8),4]
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_obese1,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_obese1,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_obese1,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_obese1,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_obese1,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"taxane")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_obese1,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"taxane")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_obese1,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"taxane")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_obese1,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"taxane")
cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
hr_adj_ob1_taxane = cox4[c(2,4,6,8),2]
lhr_adj_ob1_taxane = cox4[c(2,4,6,8),3]
uhr_adj_ob1_taxane = cox4[c(2,4,6,8),4]
cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_obese2,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_obese2,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_obese2,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_obese2,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"taxane")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_obese2,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"taxane")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_obese2,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"taxane")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_obese2,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"taxane")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_obese2,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"taxane")
cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
hr_adj_ob2_taxane = cox5[c(2,4,6,8),2]
lhr_adj_ob2_taxane = cox5[c(2,4,6,8),3]
uhr_adj_ob2_taxane = cox5[c(2,4,6,8),4]
cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))


cox_all <- data.frame(cbind(cox1, '', cox2, '', cox3, '', cox4, '', cox5))
write.csv(cox_all,'T6c_cox_all.csv')


### Table 6d, excluding 20 participants with LVEF <= 50 for Cyclophasphamide (BMI Category)
# Chemotherapy:  Cyclophasphamide: (yes/no) 
phs.aim2.data$cyclo = as.factor(phs.aim2.data$cyclo_yn)

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_underweight,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_underweight,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_underweight,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_underweight,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_underweight,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"cyclo")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_underweight,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"cyclo")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_underweight,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"cyclo")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_underweight,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"cyclo")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_normal,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_normal,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_normal,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_normal,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_normal,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"cyclo")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_normal,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"cyclo")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_normal,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"cyclo")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_normal,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"cyclo")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
hr_adj_normal_cyclo = cox2[c(2,4,6,8),2]
lhr_adj_normal_cyclo = cox2[c(2,4,6,8),3]
uhr_adj_normal_cyclo = cox2[c(2,4,6,8),4]
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_overweight,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_overweight,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_overweight,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_overweight,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_overweight,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"cyclo")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_overweight,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"cyclo")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_overweight,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"cyclo")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_overweight,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"cyclo")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
hr_adj_over_cyclo = cox3[c(2,4,6,8),2]
lhr_adj_over_cyclo = cox3[c(2,4,6,8),3]
uhr_adj_over_cyclo = cox3[c(2,4,6,8),4]
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_obese1,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_obese1,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_obese1,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_obese1,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_obese1,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"cyclo")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_obese1,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"cyclo")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_obese1,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"cyclo")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_obese1,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"cyclo")
cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
hr_adj_ob1_cyclo = cox4[c(2,4,6,8),2]
lhr_adj_ob1_cyclo = cox4[c(2,4,6,8),3]
uhr_adj_ob1_cyclo = cox4[c(2,4,6,8),4]
cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data[a2_ihd_ef_obese2,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_stroke <- coxtab2(phs.aim2.data[a2_stroke_ef_obese2,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_hf_cm <- coxtab2(phs.aim2.data[a2_hf_cm_ef_obese2,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_cvdcombo <- coxtab2(phs.aim2.data[a2_combo_ef_obese2,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"cyclo")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_obese2,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"cyclo")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_obese2,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"cyclo")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_obese2,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"cyclo")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_obese2,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"cyclo")
cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
hr_adj_ob2_cyclo = cox5[c(2,4,6,8),2]
lhr_adj_ob2_cyclo = cox5[c(2,4,6,8),3]
uhr_adj_ob2_cyclo = cox5[c(2,4,6,8),4]
cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))


cox_all <- data.frame(cbind(cox1, '', cox2, '', cox3, '', cox4, '', cox5))
write.csv(cox_all,'T6d_cox_all.csv')


### Table 6e, excluding 20 participants with LVEF <= 50 and 19 participants receiving Epirubicin for Doxorubicin (BMI Category)
# Chemotherapy:  Doxorubicin: (yes/no) 
phs.aim2.data.doxo$doxo = as.factor(phs.aim2.data.doxo$doxorub_yn)

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_underweight_doxo,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_underweight_doxo,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_underweight_doxo,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_underweight_doxo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_underweight_doxo,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"doxo")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_underweight_doxo,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"doxo")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_underweight_doxo,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"doxo")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_underweight_doxo,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"doxo")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_normal_doxo,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_normal_doxo,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_normal_doxo,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_normal_doxo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_normal_doxo,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"doxo")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_normal_doxo,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"doxo")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_normal_doxo,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"doxo")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_normal_doxo,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"doxo")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
hr_adj_normal_doxo = cox2[c(2,4,6,8),2]
lhr_adj_normal_doxo = cox2[c(2,4,6,8),3]
uhr_adj_normal_doxo = cox2[c(2,4,6,8),4]
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_overweight_doxo,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_overweight_doxo,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_overweight_doxo,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_overweight_doxo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_overweight_doxo,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"doxo")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_overweight_doxo,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"doxo")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_overweight_doxo,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"doxo")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_overweight_doxo,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"doxo")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
hr_adj_over_doxo = cox3[c(2,4,6,8),2]
lhr_adj_over_doxo = cox3[c(2,4,6,8),3]
uhr_adj_over_doxo = cox3[c(2,4,6,8),4]
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_obese1_doxo,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_obese1_doxo,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_obese1_doxo,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_obese1_doxo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_obese1_doxo,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"doxo")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_obese1_doxo,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"doxo")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_obese1_doxo,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"doxo")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_obese1_doxo,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"doxo")
cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
hr_adj_ob1_doxo = cox4[c(2,4,6,8),2]
lhr_adj_ob1_doxo = cox4[c(2,4,6,8),3]
uhr_adj_ob1_doxo = cox4[c(2,4,6,8),4]
cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_obese2_doxo,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_obese2_doxo,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_obese2_doxo,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_obese2_doxo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_htn <- coxtab2(phs.aim2.data[a2_htn_ef_obese2_doxo,],'cvdrf_htn',prevcvd[c(-2,-6,-7)],"doxo")
cox_diab <- coxtab2(phs.aim2.data[a2_diab_ef_obese2_doxo,],'cvdrf_diab',prevcvd[c(-2,-6,-7)],"doxo")
cox_dyslipid <- coxtab2(phs.aim2.data[a2_dyslipid_ef_obese2_doxo,],'cvdrf_dyslipid',prevcvd[c(-2,-6,-7)],"doxo")
cox_rfcomb <- coxtab2(phs.aim2.data[a2_rfcombo_ef_obese2_doxo,],'cvdrfcombo',prevcvd[c(-2,-6,-7)],"doxo")
cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,],cox_htn[1:2,],cox_diab[1:2,], cox_dyslipid[1:2,], 
              cox_rfcomb[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
hr_adj_ob2_doxo = cox5[c(2,4,6,8),2]
lhr_adj_ob2_doxo = cox5[c(2,4,6,8),3]
uhr_adj_ob2_doxo = cox5[c(2,4,6,8),4]
cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))


cox_all <- data.frame(cbind(cox1, '', cox2, '', cox3, '', cox4, '', cox5))
write.csv(cox_all,'T6e_cox_all.csv')



# summary statistics for dose in table 6
s1 = summary(phs.aim2.data.lvefgr50$anthra_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Underweight"])
s2 = summary(phs.aim2.data.lvefgr50$anthra_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Normal"])
s3 = summary(phs.aim2.data.lvefgr50$anthra_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Overweight"])
s4 = summary(phs.aim2.data.lvefgr50$anthra_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Obese I"])
s5 = summary(phs.aim2.data.lvefgr50$anthra_tot_dose[phs.aim2.data.doxo.lvefgr50$bmicat1 == "Obese II+"])
s.anthra = c(s1[4],s1[3],s1[1],s1[6],s2[4],s2[3],s2[1],s2[6],s3[4],s3[3],s3[1],s3[6],s4[4],s4[3],s4[1],s4[6],s5[4],s5[3],s5[1],s5[6])

s1 = summary(phs.aim2.data.lvefgr50$tras_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Underweight"])
s2 = summary(phs.aim2.data.lvefgr50$tras_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Normal"])
s3 = summary(phs.aim2.data.lvefgr50$tras_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Overweight"])
s4 = summary(phs.aim2.data.lvefgr50$tras_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Obese I"])
s5 = summary(phs.aim2.data.lvefgr50$tras_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Obese II+"])
s.tras = c(s1[4],s1[3],s1[1],s1[6],s2[4],s2[3],s2[1],s2[6],s3[4],s3[3],s3[1],s3[6],s4[4],s4[3],s4[1],s4[6],s5[4],s5[3],s5[1],s5[6])

s1 = summary(phs.aim2.data.lvefgr50$taxane_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Underweight"])
s2 = summary(phs.aim2.data.lvefgr50$taxane_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Normal"])
s3 = summary(phs.aim2.data.lvefgr50$taxane_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Overweight"])
s4 = summary(phs.aim2.data.lvefgr50$taxane_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Obese I"])
s5 = summary(phs.aim2.data.lvefgr50$taxane_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Obese II+"])
s.taxane = c(s1[4],s1[3],s1[1],s1[6],s2[4],s2[3],s2[1],s2[6],s3[4],s3[3],s3[1],s3[6],s4[4],s4[3],s4[1],s4[6],s5[4],s5[3],s5[1],s5[6])


s1 = summary(phs.aim2.data.lvefgr50$cyclo_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Underweight"])
s2 = summary(phs.aim2.data.lvefgr50$cyclo_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Normal"])
s3 = summary(phs.aim2.data.lvefgr50$cyclo_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Overweight"])
s4 = summary(phs.aim2.data.lvefgr50$cyclo_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Obese I"])
s5 = summary(phs.aim2.data.lvefgr50$cyclo_tot_dose[phs.aim2.data.lvefgr50$bmicat1 == "Obese II+"])
s.cyclo = c(s1[4],s1[3],s1[1],s1[6],s2[4],s2[3],s2[1],s2[6],s3[4],s3[3],s3[1],s3[6],s4[4],s4[3],s4[1],s4[6],s5[4],s5[3],s5[1],s5[6])


s1 = summary(phs.aim2.data.doxo.lvefgr50$doxorub_tot_dose[phs.aim2.data.doxo.lvefgr50$bmicat1 == "Underweight"])
s2 = summary(phs.aim2.data.doxo.lvefgr50$doxorub_tot_dose[phs.aim2.data.doxo.lvefgr50$bmicat1 == "Normal"])
s3 = summary(phs.aim2.data.doxo.lvefgr50$doxorub_tot_dose[phs.aim2.data.doxo.lvefgr50$bmicat1 == "Overweight"])
s4 = summary(phs.aim2.data.doxo.lvefgr50$doxorub_tot_dose[phs.aim2.data.doxo.lvefgr50$bmicat1 == "Obese I"])
s5 = summary(phs.aim2.data.doxo.lvefgr50$doxorub_tot_dose[phs.aim2.data.doxo.lvefgr50$bmicat1 == "Obese II+"])
s.doxorub = c(s1[4],s1[3],s1[1],s1[6],s2[4],s2[3],s2[1],s2[6],s3[4],s3[3],s3[1],s3[6],s4[4],s4[3],s4[1],s4[6],s5[4],s5[3],s5[1],s5[6])

s.all = rbind(s.anthra,s.tras,s.taxane,s.cyclo,s.doxorub)

write.csv(data.frame(s.all), "s_all.csv")

t1 = table(phs.aim2.data.lvefgr50$anthra[phs.aim2.data.lvefgr50$bmicat1 == "Underweight"])[2:1]
t2 = table(phs.aim2.data.lvefgr50$anthra[phs.aim2.data.lvefgr50$bmicat1 == "Normal"])[2:1]
t3 = table(phs.aim2.data.lvefgr50$anthra[phs.aim2.data.lvefgr50$bmicat1 == "Overweight"])[2:1]
t4 = table(phs.aim2.data.lvefgr50$anthra[phs.aim2.data.lvefgr50$bmicat1 == "Obese I"])[2:1]
t5 = table(phs.aim2.data.lvefgr50$anthra[phs.aim2.data.lvefgr50$bmicat1 == "Obese II+"])[2:1]
t.anthra = c(t1, t2, t3, t4, t5)

t1 = table(phs.aim2.data.lvefgr50$tras[phs.aim2.data.lvefgr50$bmicat1 == "Underweight"])[2:1]
t2 = table(phs.aim2.data.lvefgr50$tras[phs.aim2.data.lvefgr50$bmicat1 == "Normal"])[2:1]
t3 = table(phs.aim2.data.lvefgr50$tras[phs.aim2.data.lvefgr50$bmicat1 == "Overweight"])[2:1]
t4 = table(phs.aim2.data.lvefgr50$tras[phs.aim2.data.lvefgr50$bmicat1 == "Obese I"])[2:1]
t5 = table(phs.aim2.data.lvefgr50$tras[phs.aim2.data.lvefgr50$bmicat1 == "Obese II+"])[2:1]
t.tras = c(t1, t2, t3, t4, t5)

t1 = table(phs.aim2.data.lvefgr50$taxane[phs.aim2.data.lvefgr50$bmicat1 == "Underweight"])[2:1]
t2 = table(phs.aim2.data.lvefgr50$taxane[phs.aim2.data.lvefgr50$bmicat1 == "Normal"])[2:1]
t3 = table(phs.aim2.data.lvefgr50$taxane[phs.aim2.data.lvefgr50$bmicat1 == "Overweight"])[2:1]
t4 = table(phs.aim2.data.lvefgr50$taxane[phs.aim2.data.lvefgr50$bmicat1 == "Obese I"])[2:1]
t5 = table(phs.aim2.data.lvefgr50$taxane[phs.aim2.data.lvefgr50$bmicat1 == "Obese II+"])[2:1]
t.taxane = c(t1, t2, t3, t4, t5)

t1 = table(phs.aim2.data.lvefgr50$cyclo_yn[phs.aim2.data.lvefgr50$bmicat1 == "Underweight"])[2:1]
t2 = table(phs.aim2.data.lvefgr50$cyclo_yn[phs.aim2.data.lvefgr50$bmicat1 == "Normal"])[2:1]
t3 = table(phs.aim2.data.lvefgr50$cyclo_yn[phs.aim2.data.lvefgr50$bmicat1 == "Overweight"])[2:1]
t4 = table(phs.aim2.data.lvefgr50$cyclo_yn[phs.aim2.data.lvefgr50$bmicat1 == "Obese I"])[2:1]
t5 = table(phs.aim2.data.lvefgr50$cyclo_yn[phs.aim2.data.lvefgr50$bmicat1 == "Obese II+"])[2:1]
t.cyclo = c(t1, t2, t3, t4, t5)

t1 = table(phs.aim2.data.doxo.lvefgr50$doxorub_yn[phs.aim2.data.doxo.lvefgr50$bmicat1 == "Underweight"])[2:1]
t2 = table(phs.aim2.data.doxo.lvefgr50$doxorub_yn[phs.aim2.data.doxo.lvefgr50$bmicat1 == "Normal"])[2:1]
t3 = table(phs.aim2.data.doxo.lvefgr50$doxorub_yn[phs.aim2.data.doxo.lvefgr50$bmicat1 == "Overweight"])[2:1]
t4 = table(phs.aim2.data.doxo.lvefgr50$doxorub_yn[phs.aim2.data.doxo.lvefgr50$bmicat1 == "Obese I"])[2:1]
t5 = table(phs.aim2.data.doxo.lvefgr50$doxorub_yn[phs.aim2.data.doxo.lvefgr50$bmicat1 == "Obese II+"])[2:1]
t.doxorub = c(t1, t2, t3, t4, t5)

t.all = rbind(t.anthra,t.tras,t.taxane,t.cyclo,t.doxorub)
t.all = as.data.frame(t.all)

for (i in 1:dim(t.all)[1]){
  for (j in 1:dim(t.all)[2]){
    t.all[i,j] = paste0("n = ",t.all[i,j])
  }
}

write.csv(data.frame(t.all), "t_all.csv")




png("Figure_BMI_Doxorubicin.png", width = 9, height=5.5, res=300, units = 'in')

bmi = c(rep("Normal",4),rep("Overweight",4),rep("Obese I",4),rep("Obese II+",4))
bmi = factor(bmi, levels = c("Normal","Overweight","Obese I","Obese II+"))
fig.doxo.dat = data.frame(bmi = bmi)
fig.doxo.dat$hr = c(hr_adj_normal_doxo,hr_adj_over_doxo,hr_adj_ob1_doxo,hr_adj_ob2_doxo)
fig.doxo.dat$lower = c(lhr_adj_normal_doxo,lhr_adj_over_doxo,lhr_adj_ob1_doxo,lhr_adj_ob2_doxo)
fig.doxo.dat$upper = c(uhr_adj_normal_doxo,uhr_adj_over_doxo,uhr_adj_ob1_doxo,uhr_adj_ob2_doxo)
fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"),4)
fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"))

pd <- position_dodge(width = 0.3) # move them .2 to the left and right
ggplot(fig.doxo.dat, aes(x = bmi,
                         y = hr, 
                         group=cvd, 
                         color=cvd)) +
  geom_point(size = 2, position=pd) +
  #geom_line(size = 1) +
  geom_errorbar(aes(ymin  =lower, 
                    ymax = upper), 
                width = .1, position=pd) + ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank()) + 
  labs(y = "Adjusted Hazard Ratio", 
       x = "BMI Category")

dev.off()




# Table 7
coxtab3 = function(d,x,covar,trt){
  covars <- paste(covar, collapse = '+')
  surv_object <- Surv(time = floor(d[,paste0(x,'_fu')]/30), event = d[,x])
  fit.coxph1 <- eval(parse(text=(paste0("coxph(surv_object ~ ", trt, ", data = d)"))))
  fit.coxph2 <- eval(parse(text=(paste0("coxph(surv_object ~ ", trt, " + ajcc_stage + diab_bl +  
                                        htn_bl + dyslipid_bl + bmicat1 + smok + charlson + 
                                       agegrp + race + edu_cat + income_cat + ",covars,", data = d)"))))
  sum1 <- cbind(summary(fit.coxph1)$coef[,c(2,5)],summary(fit.coxph1)$conf.int[,3:4])
  row.names(sum1) = c("1","2","3")
  sum2 <- cbind(summary(fit.coxph2)$coef[,c(2,5)],summary(fit.coxph2)$conf.int[,3:4])
  sum3 <- rbind(sum1, sum2)[,c(1,3,4,2)]
  sum4 <- data.frame(cbind(var=row.names(sum3), sum3))
  sum4$var <- factor(sum4$var, levels=row.names(sum3))
  sum4
}
phs.aim2.data$chemo_anthra_tras = 0
phs.aim2.data$chemo_anthra_tras[phs.aim2.data$comb1 == 1] = 1
phs.aim2.data$chemo_anthra_tras[phs.aim2.data$comb5 == 1] = 2
phs.aim2.data$chemo_anthra_tras[phs.aim2.data$comb3 == 1 | phs.aim2.data$comb4 == 1] = 3
phs.aim2.data$chemo_anthra_tras = as.factor(phs.aim2.data$chemo_anthra_tras)

cox_isch <- coxtab3(phs.aim2.data[a2_ihd_ef,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"chemo_anthra_tras")
cox_stroke <- coxtab3(phs.aim2.data[a2_stroke_ef,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"chemo_anthra_tras")
cox_hf_cm <- coxtab3(phs.aim2.data[a2_hf_cm_ef,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"chemo_anthra_tras")
cox_cvdcombo <- coxtab3(phs.aim2.data[a2_combo_ef,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"chemo_anthra_tras")
cox1 <- rbind(cox_isch[1:6,],cox_hf_cm[1:6,], cox_stroke[1:6,], 
              cox_cvdcombo[1:6,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))
write.csv(cox1,'T7_cox_all.csv')




########### Model Assesment
coxtab_chemo <- function(d,x,covar,trt){
  #covars <- paste(covar, collapse = '+')
  surv_object <- Surv(time = floor(d[,paste0(x,'_fu')]/30), event = d[,x])
  fit.coxph1 <- eval(parse(text=(paste0("coxph(surv_object ~ ", trt, ", data = d)"))))
  fit.coxph2 <- eval(parse(text=(paste0("coxph(surv_object ~ ", trt, " + ", covar, ", data = d)"))))
  sum1 <- c(summary(fit.coxph1)$coef[c(2,5)],summary(fit.coxph1)$conf.int[3:4])
  sum2 <- cbind(summary(fit.coxph2)$coef[,c(2,5)],summary(fit.coxph2)$conf.int[,3:4])
  sum3 <- rbind(sum1, sum2)[,c(1,3,4,2)]
  sum4 <- data.frame(cbind(var=row.names(sum3), sum3))
  sum4$var <- factor(sum4$var, levels=row.names(sum3))
  sum4[2,]
}
coxtab_chemo(phs.aim2.data[a2_ihd_ef,],'ischemic_heart_disease_grp_inc',"ajcc_stage","anthra")

cov.list = c("ajcc_stage","diab_bl","htn_bl","dyslipid_bl","bmicat1","smok","charlson","agegrp","race","edu_cat","income_cat","prevcvd","horm_yn","rad", "menop")
cox_isch <- do.call(rbind,lapply(cov.list,function(x) coxtab_chemo(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',x,"taxane")))
cox_isch <- data.frame(cbind(cov.list, cox_isch))
cox_isch

selectCox(Surv(ischemic_heart_disease_grp_inc_fu,ischemic_heart_disease_grp_inc) ~ taxane + ajcc_stage + diab_bl +  
            htn_bl + dyslipid_bl + bmicat1 + smok + charlson + agegrp + race + edu_cat + income_cat + prevcvd + horm_yn + rad,
          data=phs.aim2.data[a2_ihd,])

selectCox(Surv(ischemic_heart_disease_grp_inc_fu,ischemic_heart_disease_grp_inc) ~ taxane + diab_bl + htn_bl + dyslipid_bl + charlson,
          data=phs.aim2.data[a2_ihd,])

stepAIC(coxph(Surv(ischemic_heart_disease_grp_inc_fu,ischemic_heart_disease_grp_inc) ~ taxane + ajcc_stage + diab_bl +  
        htn_bl + dyslipid_bl + bmicat1 + smok + charlson + agegrp + race + edu_cat + income_cat + prevcvd + horm_yn + rad,
      data=phs.aim2.data[a2_ihd,]), direction = "both")


coxph(Surv(ischemic_heart_disease_grp_inc_fu,ischemic_heart_disease_grp_inc) ~ taxane + diab_bl +  
        dyslipid_bl + charlson + agegrp + race,
      data=phs.aim2.data[a2_ihd,])

dat2 = phs.aim2.data[,cov.list]
dat2[] <- lapply(dat2, function(x) {
  if(is.factor(x)) as.numeric(x) else x
})
sapply(dat2, class)
require(corrplot)
cor(dat2, use = "complete.obs")
corrplot(cor(dat2,use = "complete.obs"))

# Chemotherapy:  Anthracycline: (yes/no) LVEF > 50
cox_isch <- coxtab(phs.aim2.data[a2_ihd_ef,],'ischemic_heart_disease_grp_inc',prevcvd[-2],"anthra")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_ef,],'stroke_grp_inc',prevcvd[-2],"anthra")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_ef,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[-2],"anthra")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_ef,],'cvdcombo_grp_inc',prevcvd[-2],"anthra")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))


# Chemotherapy:  Trastuzumab: (yes/no) LVEF > 50
cox_isch <- coxtab(phs.aim2.data[a2_ihd_ef,],'ischemic_heart_disease_grp_inc',prevcvd[-2],"tras")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke_ef,],'stroke_grp_inc',prevcvd[-2],"tras")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm_ef,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[-2],"tras")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo_ef,],'cvdcombo_grp_inc',prevcvd[-2],"tras")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))


# Chemotherapy:  Taxane: (yes/no) 	
cox_isch <- coxtab(phs.aim2.data[a2_ihd,],'ischemic_heart_disease_grp_inc',prevcvd[-2],"taxane")
cox_stroke <- coxtab(phs.aim2.data[a2_stroke,],'stroke_grp_inc',prevcvd[-2],"taxane")
cox_hf_cm <- coxtab(phs.aim2.data[a2_hf_cm,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[-2],"taxane")
cox_cvdcombo <- coxtab(phs.aim2.data[a2_combo,],'cvdcombo_grp_inc',prevcvd[-2],"taxane")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
              cox_cvdcombo[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))

cox_all <- data.frame(cbind(cox1, '', cox2,'', cox3))
write_csv(cox_all,'cox_all_chemo.csv')




############
# caculate sample size
# Table 3a1
coxns <- data.frame(t(rbind(c(nrow(phs.aim2.data),
                              sapply(list(a2_ihd,a2_hf_cm,a2_stroke,a2_combo),length)))))
table(phs.aim2.data[a2_ihd,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo,]$cvdcombo_grp_inc)



# Table 3a2
coxns <- data.frame(t(rbind(c(nrow(phs.aim2.data.cvd),
                              sapply(list(a2_ihd_cvd,a2_hf_cm_cvd,a2_stroke_cvd,a2_combo_cvd),length)))))
table(phs.aim2.data[a2_ihd_cvd,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_cvd,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_cvd,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_cvd,]$cvdcombo_grp_inc)



# Table 3b1
coxns <- data.frame(t(rbind(c(nrow(phs.aim2.data_lvefgr50),
                              sapply(list(a2_ihd_ef,a2_hf_cm_ef,a2_stroke_ef,a2_combo_ef),length)))))
table(phs.aim2.data[a2_ihd_ef,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_ef,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_ef,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_ef,]$cvdcombo_grp_inc)


# Table 3b2
coxns <- data.frame(t(rbind(c(nrow(a2_ef_cvd),
                              sapply(list(a2_ihd_ef_cvd,a2_hf_cm_ef_cvd,a2_stroke_ef_cvd,a2_combo_ef_cvd),length)))))
table(phs.aim2.data[a2_ihd_ef_cvd,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_ef_cvd,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_ef_cvd,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_ef_cvd,]$cvdcombo_grp_inc)



# Table 4a
coxns <- data.frame(t(rbind(c(nrow(phs.aim2.data),
                              sapply(list(a2_ihd_post_menop,a2_hf_cm_post_menop,a2_stroke_post_menop,a2_combo_post_menop),length)))))
table(phs.aim2.data[a2_ihd_post_menop,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_post_menop,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_post_menop,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_post_menop,]$cvdcombo_grp_inc)

coxns <- data.frame(t(rbind(c(nrow(phs.aim2.data),
                              sapply(list(a2_ihd_pre_menop,a2_hf_cm_pre_menop,
                                          a2_stroke_pre_menop,a2_combo_pre_menop),length)))))
table(phs.aim2.data[a2_ihd_pre_menop,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_pre_menop,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_pre_menop,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_pre_menop,]$cvdcombo_grp_inc)


# Table 4b
coxns <- data.frame(t(rbind(c(nrow(a2_post_menop_cvd),
                              sapply(list(a2_ihd_post_menop_cvd,a2_hf_cm_post_menop_cvd,a2_stroke_post_menop_cvd,a2_combo_post_menop_cvd),length)))))
table(phs.aim2.data[a2_ihd_post_menop_cvd,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_post_menop_cvd,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_post_menop_cvd,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_post_menop_cvd,]$cvdcombo_grp_inc)

coxns <- data.frame(t(rbind(c(nrow(a2_pre_menop_cvd),
                              sapply(list(a2_ihd_pre_menop_cvd,a2_hf_cm_pre_menop_cvd,
                                          a2_stroke_pre_menop_cvd,a2_combo_pre_menop_cvd),length)))))
table(phs.aim2.data[a2_ihd_pre_menop_cvd,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_pre_menop_cvd,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_pre_menop_cvd,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_pre_menop_cvd,]$cvdcombo_grp_inc)


# Table 5b
coxns <- data.frame(t(rbind(c(nrow(phs.aim2.data.cvd),
                              sapply(list(a2_ihd_cvd,a2_hf_cm_cvd,a2_stroke_cvd,a2_combo_cvd),length)))))
table(phs.aim2.data[a2_ihd_cvd,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_cvd,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_cvd,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_cvd,]$cvdcombo_grp_inc)


coxns <- data.frame(t(rbind(c(nrow(a2_rad_l_cvd),
                              sapply(list(a2_ihd_l_cvd,a2_hf_cm_l_cvd,a2_stroke_l_cvd,a2_combo_l_cvd),length)))))
table(phs.aim2.data[a2_ihd_l_cvd,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_l_cvd,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_l_cvd,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_l_cvd,]$cvdcombo_grp_inc)

# Table 6
coxns <- data.frame(t(rbind(c(nrow(phs.aim2.data),
                              sapply(list(a2_ihd_ef_underweight,a2_hf_cm_ef_underweight,a2_stroke_ef_underweight,a2_combo_ef_underweight),length)))))
table(phs.aim2.data[a2_ihd_ef_underweight,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_ef_underweight,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_ef_underweight,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_ef_underweight,]$cvdcombo_grp_inc)

table(phs.aim2.data.doxo[a2_ihd_ef_underweight_doxo,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data.doxo[a2_hf_cm_ef_underweight_doxo,]$heart_failure_grp_inc)
table(phs.aim2.data.doxo[a2_stroke_ef_underweight_doxo,]$stroke_grp_inc)
table(phs.aim2.data.doxo[a2_combo_ef_underweight_doxo,]$cvdcombo_grp_inc)

coxns <- data.frame(t(rbind(c(nrow(phs.aim2.data),
                              sapply(list(a2_ihd_ef_normal,a2_hf_cm_ef_normal,a2_stroke_ef_normal,a2_combo_ef_normal),length)))))
table(phs.aim2.data[a2_ihd_ef_normal,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_ef_normal,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_ef_normal,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_ef_normal,]$cvdcombo_grp_inc)

table(phs.aim2.data.doxo[a2_ihd_ef_normal_doxo,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data.doxo[a2_hf_cm_ef_normal_doxo,]$heart_failure_grp_inc)
table(phs.aim2.data.doxo[a2_stroke_ef_normal_doxo,]$stroke_grp_inc)
table(phs.aim2.data.doxo[a2_combo_ef_normal_doxo,]$cvdcombo_grp_inc)


coxns <- data.frame(t(rbind(c(nrow(phs.aim2.data),
                              sapply(list(a2_ihd_ef_overweight,a2_hf_cm_ef_overweight,a2_stroke_ef_overweight,a2_combo_ef_overweight),length)))))
table(phs.aim2.data[a2_ihd_ef_overweight,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_ef_overweight,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_ef_overweight,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_ef_overweight,]$cvdcombo_grp_inc)

table(phs.aim2.data.doxo[a2_ihd_ef_overweight_doxo,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data.doxo[a2_hf_cm_ef_overweight_doxo,]$heart_failure_grp_inc)
table(phs.aim2.data.doxo[a2_stroke_ef_overweight_doxo,]$stroke_grp_inc)
table(phs.aim2.data.doxo[a2_combo_ef_overweight_doxo,]$cvdcombo_grp_inc)

coxns <- data.frame(t(rbind(c(nrow(phs.aim2.data),
                              sapply(list(a2_ihd_ef_obese1,a2_hf_cm_ef_obese1,a2_stroke_ef_obese1,a2_combo_ef_obese1),length)))))
table(phs.aim2.data[a2_ihd_ef_obese1,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_ef_obese1,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_ef_obese1,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_ef_obese1,]$cvdcombo_grp_inc)

table(phs.aim2.data.doxo[a2_ihd_ef_obese1_doxo,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data.doxo[a2_hf_cm_ef_obese1_doxo,]$heart_failure_grp_inc)
table(phs.aim2.data.doxo[a2_stroke_ef_obese1_doxo,]$stroke_grp_inc)
table(phs.aim2.data.doxo[a2_combo_ef_obese1_doxo,]$cvdcombo_grp_inc)

coxns <- data.frame(t(rbind(c(nrow(phs.aim2.data),
                              sapply(list(a2_ihd_ef_obese2,a2_hf_cm_ef_obese2,a2_stroke_ef_obese2,a2_combo_ef_obese2),length)))))
table(phs.aim2.data[a2_ihd_ef_obese2,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data[a2_hf_cm_ef_obese2,]$heart_failure_grp_inc)
table(phs.aim2.data[a2_stroke_ef_obese2,]$stroke_grp_inc)
table(phs.aim2.data[a2_combo_ef_obese2,]$cvdcombo_grp_inc)

table(phs.aim2.data.doxo[a2_ihd_ef_obese2_doxo,]$ischemic_heart_disease_grp_inc)
table(phs.aim2.data.doxo[a2_hf_cm_ef_obese2_doxo,]$heart_failure_grp_inc)
table(phs.aim2.data.doxo[a2_stroke_ef_obese2_doxo,]$stroke_grp_inc)
table(phs.aim2.data.doxo[a2_combo_ef_obese2_doxo,]$cvdcombo_grp_inc)

## risk factors
htn_n = c(length(a2_htn_ef_underweight),length(a2_htn_ef_normal),length(a2_htn_ef_overweight),length(a2_htn_ef_obese1),length(a2_htn_ef_obese2))
htn_ev = c(table(phs.aim2.data[a2_htn_ef_underweight,]$cvdrf_htn)[2],
           table(phs.aim2.data[a2_htn_ef_normal,]$cvdrf_htn)[2],
           table(phs.aim2.data[a2_htn_ef_overweight,]$cvdrf_htn)[2],
           table(phs.aim2.data[a2_htn_ef_obese1,]$cvdrf_htn)[2],
           table(phs.aim2.data[a2_htn_ef_obese2,]$cvdrf_htn)[2])
diab_n = c(length(a2_diab_ef_underweight),length(a2_diab_ef_normal),length(a2_diab_ef_overweight),length(a2_diab_ef_obese1),length(a2_diab_ef_obese2))
diab_ev = c(table(phs.aim2.data[a2_diab_ef_underweight,]$cvdrf_diab)[2],
            table(phs.aim2.data[a2_diab_ef_normal,]$cvdrf_diab)[2],
            table(phs.aim2.data[a2_diab_ef_overweight,]$cvdrf_diab)[2],
            table(phs.aim2.data[a2_diab_ef_obese1,]$cvdrf_diab)[2],
            table(phs.aim2.data[a2_diab_ef_obese2,]$cvdrf_diab)[2])
dyslipid_n = c(length(a2_dyslipid_ef_underweight),length(a2_dyslipid_ef_normal),length(a2_dyslipid_ef_overweight),length(a2_dyslipid_ef_obese1),length(a2_dyslipid_ef_obese2))
dyslipid_ev = c(table(phs.aim2.data[a2_dyslipid_ef_underweight,]$cvdrf_dyslipid)[2],
                table(phs.aim2.data[a2_dyslipid_ef_normal,]$cvdrf_dyslipid)[2],
                table(phs.aim2.data[a2_dyslipid_ef_overweight,]$cvdrf_dyslipid)[2],
                table(phs.aim2.data[a2_dyslipid_ef_obese1,]$cvdrf_dyslipid)[2],
                table(phs.aim2.data[a2_dyslipid_ef_obese2,]$cvdrf_dyslipid)[2])
rfcombo_n = c(length(a2_rfcombo_ef_underweight),length(a2_rfcombo_ef_normal),length(a2_rfcombo_ef_overweight),length(a2_rfcombo_ef_obese1),length(a2_rfcombo_ef_obese2))
rfcombo_ev = c(table(phs.aim2.data[a2_rfcombo_ef_underweight,]$cvdrfcombo)[2],
               table(phs.aim2.data[a2_rfcombo_ef_normal,]$cvdrfcombo)[2],
               table(phs.aim2.data[a2_rfcombo_ef_overweight,]$cvdrfcombo)[2],
               table(phs.aim2.data[a2_rfcombo_ef_obese1,]$cvdrfcombo)[2],
               table(phs.aim2.data[a2_rfcombo_ef_obese2,]$cvdrfcombo)[2])


htn_n = c(length(a2_htn_ef_underweight_doxo),length(a2_htn_ef_normal_doxo),length(a2_htn_ef_overweight_doxo),length(a2_htn_ef_obese1_doxo),length(a2_htn_ef_obese2_doxo))
htn_ev = c(table(phs.aim2.data.doxo[a2_htn_ef_underweight_doxo,]$cvdrf_htn)[2],
           table(phs.aim2.data.doxo[a2_htn_ef_normal_doxo,]$cvdrf_htn)[2],
           table(phs.aim2.data.doxo[a2_htn_ef_overweight_doxo,]$cvdrf_htn)[2],
           table(phs.aim2.data.doxo[a2_htn_ef_obese1_doxo,]$cvdrf_htn)[2],
           table(phs.aim2.data.doxo[a2_htn_ef_obese2_doxo,]$cvdrf_htn)[2])
diab_n = c(length(a2_diab_ef_underweight_doxo),length(a2_diab_ef_normal_doxo),length(a2_diab_ef_overweight_doxo),length(a2_diab_ef_obese1_doxo),length(a2_diab_ef_obese2_doxo))
diab_ev = c(table(phs.aim2.data.doxo[a2_diab_ef_underweight_doxo,]$cvdrf_diab)[2],
            table(phs.aim2.data.doxo[a2_diab_ef_normal_doxo,]$cvdrf_diab)[2],
            table(phs.aim2.data.doxo[a2_diab_ef_overweight_doxo,]$cvdrf_diab)[2],
            table(phs.aim2.data.doxo[a2_diab_ef_obese1_doxo,]$cvdrf_diab)[2],
            table(phs.aim2.data.doxo[a2_diab_ef_obese2_doxo,]$cvdrf_diab)[2])
dyslipid_n = c(length(a2_dyslipid_ef_underweight_doxo),length(a2_dyslipid_ef_normal_doxo),length(a2_dyslipid_ef_overweight_doxo),length(a2_dyslipid_ef_obese1_doxo),length(a2_dyslipid_ef_obese2_doxo))
dyslipid_ev = c(table(phs.aim2.data.doxo[a2_dyslipid_ef_underweight_doxo,]$cvdrf_dyslipid)[2],
                table(phs.aim2.data.doxo[a2_dyslipid_ef_normal_doxo,]$cvdrf_dyslipid)[2],
                table(phs.aim2.data.doxo[a2_dyslipid_ef_overweight_doxo,]$cvdrf_dyslipid)[2],
                table(phs.aim2.data.doxo[a2_dyslipid_ef_obese1_doxo,]$cvdrf_dyslipid)[2],
                table(phs.aim2.data.doxo[a2_dyslipid_ef_obese2_doxo,]$cvdrf_dyslipid)[2])
rfcombo_n = c(length(a2_rfcombo_ef_underweight_doxo),length(a2_rfcombo_ef_normal_doxo),length(a2_rfcombo_ef_overweight_doxo),length(a2_rfcombo_ef_obese1_doxo),length(a2_rfcombo_ef_obese2_doxo))
rfcombo_ev = c(table(phs.aim2.data.doxo[a2_rfcombo_ef_underweight_doxo,]$cvdrfcombo)[2],
               table(phs.aim2.data.doxo[a2_rfcombo_ef_normal_doxo,]$cvdrfcombo)[2],
               table(phs.aim2.data.doxo[a2_rfcombo_ef_overweight_doxo,]$cvdrfcombo)[2],
               table(phs.aim2.data.doxo[a2_rfcombo_ef_obese1_doxo,]$cvdrfcombo)[2],
               table(phs.aim2.data.doxo[a2_rfcombo_ef_obese2_doxo,]$cvdrfcombo)[2])

############ Additional Analyses for R01 ####################
coxtab2 <- function(d,x,covar,trt){
  covars <- paste(covar, collapse = '+')
  surv_object <- Surv(time = floor(d[,paste0(x,'_fu')]/30), event = d[,x])
  fit.coxph1 <- eval(parse(text=(paste0("coxph(surv_object ~ ", trt, ", data = d)"))))
  fit.coxph2 <- eval(parse(text=(paste0("coxph(surv_object ~ ", trt, " + ajcc_stage + diab_bl +  
                                        htn_bl + dyslipid_bl + smok + charlson + 
                                       agegrp + race + edu_cat + income_cat + ",covars,", data = d)"))))
  sum1 <- c(summary(fit.coxph1)$coef[c(2,5)],summary(fit.coxph1)$conf.int[3:4])
  sum2 <- cbind(summary(fit.coxph2)$coef[,c(2,5)],summary(fit.coxph2)$conf.int[,3:4])
  sum3 <- rbind(sum1, sum2)[,c(1,3,4,2)]
  sum4 <- data.frame(cbind(var=row.names(sum3), sum3))
  sum4$var <- factor(sum4$var, levels=row.names(sum3))
  sum4
}

prevcvd <- c("prevcvd","chemo_yn","horm_yn","rad_tx_yn", "menop", "lvef_ind", "I(lvef*lvef_ind)")


phs.aim2.data.doxo$doxo = as.factor(phs.aim2.data.doxo$doxorub_yn)

# BMI category
cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_underweight_doxo,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_underweight_doxo,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_underweight_doxo,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_underweight_doxo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
hr_adj_underweight_doxo_bmi = cox1[c(2,4),2]
lhr_adj_underweight_doxo_bmi = cox1[c(2,4),3]
uhr_adj_underweight_doxo_bmi = cox1[c(2,4),4]
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_normal_doxo,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_normal_doxo,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_normal_doxo,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_normal_doxo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
hr_adj_normal_doxo_bmi = cox2[c(2,4),2]
lhr_adj_normal_doxo_bmi = cox2[c(2,4),3]
uhr_adj_normal_doxo_bmi = cox2[c(2,4),4]
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_overweight_doxo,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_overweight_doxo,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_overweight_doxo,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_overweight_doxo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
hr_adj_overweight_doxo_bmi = cox3[c(2,4),2]
lhr_adj_overweight_doxo_bmi = cox3[c(2,4),3]
uhr_adj_overweight_doxo_bmi = cox3[c(2,4),4]
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_obese1_doxo,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_obese1_doxo,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_obese1_doxo,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_obese1_doxo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
hr_adj_obese1_doxo_bmi = cox4[c(2,4),2]
lhr_adj_obese1_doxo_bmi = cox4[c(2,4),3]
uhr_adj_obese1_doxo_bmi = cox4[c(2,4),4]
cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_obese2_doxo,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_obese2_doxo,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_obese2_doxo,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_obese2_doxo,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
hr_adj_obese2_doxo_bmi = cox5[c(2,4),2]
lhr_adj_obese2_doxo_bmi = cox5[c(2,4),3]
uhr_adj_obese2_doxo_bmi = cox5[c(2,4),4]
cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))

bmi = c(rep("Normal",2),rep("Overweight",2),rep("Obese I",2),rep("Obese II+",2))
bmi = factor(bmi, levels = c("Normal","Overweight","Obese I","Obese II+"))
fig.doxo.dat = data.frame(bmi = bmi)
fig.doxo.dat$hr = c(hr_adj_normal_doxo_bmi,hr_adj_overweight_doxo_bmi,hr_adj_obese1_doxo_bmi,hr_adj_obese2_doxo_bmi)
fig.doxo.dat$lower = c(lhr_adj_normal_doxo_bmi,lhr_adj_overweight_doxo_bmi,lhr_adj_obese1_doxo_bmi,lhr_adj_obese2_doxo_bmi)
fig.doxo.dat$upper = c(uhr_adj_normal_doxo_bmi,uhr_adj_overweight_doxo_bmi,uhr_adj_obese1_doxo_bmi,uhr_adj_obese2_doxo_bmi)
fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy"),4)
fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy"))
for (i in 1:dim(fig.doxo.dat)[1]){
  fig.doxo.dat$ci[i] = paste0("(95% CI: ",
                              format(round(fig.doxo.dat$lower[i],2),nsmall=2),"-", format(round(fig.doxo.dat$upper[i],2),nsmall=2),")")
}
fig.doxo.dat$ci[1] = paste0(fig.doxo.dat$ci[1],"*")
fig.doxo.dat$ci[2] = paste0(fig.doxo.dat$ci[2],"*")

p0 = ggplot(fig.doxo.dat, aes(x = bmi,
                              y = hr, 
                              fill=cvd, 
                              color=cvd)) +
  geom_bar(stat="identity", color="white", 
           position=position_dodge()) + geom_hline(yintercept=1, color = "red", size = 1) + scale_fill_manual(values = c("steelblue", "gold")) + 
  #geom_line(size = 1) +
  geom_errorbar(aes(ymin  =lower, 
                    ymax = upper), 
                width = .2, color = "black",
                position=position_dodge(.9), size = 1) + ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
                                                                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                              axis.line.x = element_line(colour = "black"),
                                                                              axis.line.y = element_line(colour = "black"),
                                                                              panel.background = element_blank()) + 
  labs(y = "Adjusted Hazard Ratio", 
       x = "BMI Category") + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))

p0_2 = ggplot(fig.doxo.dat, aes(x = bmi,
                              y = hr, 
                              fill=cvd, 
                              color=cvd)) +
  geom_bar(stat="identity", color="white", 
           position=position_dodge()) + geom_hline(yintercept=1, color = "red", size = 1) + scale_fill_manual(values = c("steelblue", "gold")) + 
  #geom_line(size = 1) +
  geom_text(aes(bmi, hr+0.21, label = format(round(hr,2),nsmall=2), colour = NULL), position=position_dodge(width=0.9), data = fig.doxo.dat) +
  geom_text(aes(bmi, hr+0.11, label = ci, colour = NULL), position=position_dodge(width=0.9), data = fig.doxo.dat, size = 3.5) +
 ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
                                                                                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                        axis.line.x = element_line(colour = "black"),
                                                                                        axis.line.y = element_line(colour = "black"),
                                                                                        panel.background = element_blank()) + 
  labs(y = "Adjusted Hazard Ratio", 
       x = "BMI Category") + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))


png("Figure_Barplots_BMI_2.png", width = 12, height=12, res=300, units = 'in')
p0_2
dev.off()

# HF/CM only
bmi = c(rep("Normal",1),rep("Overweight",1),rep("Obese I",1),rep("Obese II+",1))
bmi = factor(bmi, levels = c("Normal","Overweight","Obese I","Obese II+"))
fig.doxo.dat = data.frame(bmi = bmi)
fig.doxo.dat$hr = c(hr_adj_normal_doxo_bmi[2],hr_adj_overweight_doxo_bmi[2],hr_adj_obese1_doxo_bmi[2],hr_adj_obese2_doxo_bmi[2])
fig.doxo.dat$lower = c(lhr_adj_normal_doxo_bmi[2],lhr_adj_overweight_doxo_bmi[2],lhr_adj_obese1_doxo_bmi[2],lhr_adj_obese2_doxo_bmi[2])
fig.doxo.dat$upper = c(uhr_adj_normal_doxo_bmi[2],uhr_adj_overweight_doxo_bmi[2],uhr_adj_obese1_doxo_bmi[2],uhr_adj_obese2_doxo_bmi[2])
for (i in 1:dim(fig.doxo.dat)[1]){
  fig.doxo.dat$ci[i] = paste0("(95% CI: ",
                              format(round(fig.doxo.dat$lower[i],2),nsmall=2),"-", format(round(fig.doxo.dat$upper[i],2),nsmall=2),")")
}
fig.doxo.dat$ci[1] = paste0(fig.doxo.dat$ci[1],"*")

p0_3 = ggplot(fig.doxo.dat, aes(x = bmi,
                                y = hr)) +
  geom_bar(stat="identity", color="white", 
           position=position_dodge(), fill = "gold") + geom_hline(yintercept=1, color = "red", size = 1) + 
  #geom_line(size = 1) +
  geom_text(aes(bmi, hr+0.30, label = format(round(hr,2),nsmall=2), colour = NULL), position=position_dodge(width=0.9), data = fig.doxo.dat, size = 6.5) +
  geom_text(aes(bmi, hr+0.14, label = ci, colour = NULL), position=position_dodge(width=0.9), data = fig.doxo.dat, size = 6.5) +
  ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 axis.line.x = element_line(colour = "black"),
                                 axis.line.y = element_line(colour = "black"),
                                 panel.background = element_blank(),
                                 axis.text.x = element_text(size = 20, face='bold'),
                                 axis.text.y = element_text(size = 20, face='bold'),
                                 axis.title.x = element_text(size = 22),
                                 axis.title.y = element_text(size = 22),
                                 title = element_text(size=22, face='bold')) + 
  labs(y = "Adjusted Hazard Ratio", 
       x = "") + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))


png("Figure_Barplots_BMI_HF_CM.png", width = 12, height=12, res=300, units = 'in')
p0_3
dev.off()

# BMI
cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_1_doxo_bmi,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_1_doxo_bmi,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_1_doxo_bmi,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_1_doxo_bmi,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
hr_adj_1_doxo_bmi = cox1[c(2,4),2]
lhr_adj_1_doxo_bmi = cox1[c(2,4),3]
uhr_adj_1_doxo_bmi = cox1[c(2,4),4]
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_2_doxo_bmi,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_2_doxo_bmi,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_2_doxo_bmi,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_2_doxo_bmi,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
hr_adj_2_doxo_bmi = cox2[c(2,4),2]
lhr_adj_2_doxo_bmi = cox2[c(2,4),3]
uhr_adj_2_doxo_bmi = cox2[c(2,4),4]
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_3_doxo_bmi,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_3_doxo_bmi,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_3_doxo_bmi,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_3_doxo_bmi,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
hr_adj_3_doxo_bmi = cox3[c(2,4),2]
lhr_adj_3_doxo_bmi = cox3[c(2,4),3]
uhr_adj_3_doxo_bmi = cox3[c(2,4),4]
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_4_doxo_bmi,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_4_doxo_bmi,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_4_doxo_bmi,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_4_doxo_bmi,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
hr_adj_4_doxo_bmi = cox4[c(2,4),2]
lhr_adj_4_doxo_bmi = cox4[c(2,4),3]
uhr_adj_4_doxo_bmi = cox4[c(2,4),4]
cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_5_doxo_bmi,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_5_doxo_bmi,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_5_doxo_bmi,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_5_doxo_bmi,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
hr_adj_5_doxo_bmi = cox5[c(2,4),2]
lhr_adj_5_doxo_bmi = cox5[c(2,4),3]
uhr_adj_5_doxo_bmi = cox5[c(2,4),4]
cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))



bmi = c(rep("Q1",2),rep("Q2",2),rep("Q3",2),rep("Q4",2),rep("Q5",2))
bmi = factor(bmi, levels = c("Q1","Q2","Q3","Q4","Q5"))
fig.doxo.dat = data.frame(bmi = bmi)
fig.doxo.dat$hr = c(hr_adj_1_doxo_bmi,hr_adj_2_doxo_bmi,hr_adj_3_doxo_bmi,hr_adj_4_doxo_bmi,hr_adj_5_doxo_bmi)
fig.doxo.dat$lower = c(lhr_adj_1_doxo_bmi,lhr_adj_2_doxo_bmi,lhr_adj_3_doxo_bmi,lhr_adj_4_doxo_bmi,lhr_adj_5_doxo_bmi)
fig.doxo.dat$upper = c(uhr_adj_1_doxo_bmi,uhr_adj_2_doxo_bmi,uhr_adj_3_doxo_bmi,uhr_adj_4_doxo_bmi,uhr_adj_5_doxo_bmi)
fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy"),5)
fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy"))

p1 = ggplot(fig.doxo.dat, aes(x = bmi,
                         y = hr, 
                         fill=cvd, 
                         color=cvd)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) + geom_hline(yintercept=1, color = "red") + 
  #geom_line(size = 1) +
  geom_errorbar(aes(ymin  =lower, 
                    ymax = upper), 
                width = .2, color = "black",
                position=position_dodge(.9)) + ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
                                                                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                              axis.line.x = element_line(colour = "black"),
                                                                              axis.line.y = element_line(colour = "black"),
                                                                              panel.background = element_blank()) + 
  labs(y = "Adjusted Hazard Ratio", 
       x = "BMI Quintiles") + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9))

# png("Figure_BMI_quintiles_Doxorubicin_2.png", width = 9, height=5.5, res=300, units = 'in')
# 
# bmi = c(rep("0 - 20%",4),rep("20% - 40%",4),rep("40% - 60%",4),rep("60% - 80%",4),rep("80% - 100%",4))
# bmi = factor(bmi, levels = c("0 - 20%","20% - 40%","40% - 60%","60% - 80%","80% - 100%"))
# fig.doxo.dat = data.frame(bmi = bmi)
# fig.doxo.dat$hr = c(hr_adj_1_doxo_bmi,hr_adj_2_doxo_bmi,hr_adj_3_doxo_bmi,hr_adj_4_doxo_bmi,hr_adj_5_doxo_bmi)
# fig.doxo.dat$lower = c(lhr_adj_1_doxo_bmi,lhr_adj_2_doxo_bmi,lhr_adj_3_doxo_bmi,lhr_adj_4_doxo_bmi,lhr_adj_5_doxo_bmi)
# fig.doxo.dat$upper = c(uhr_adj_1_doxo_bmi,uhr_adj_2_doxo_bmi,uhr_adj_3_doxo_bmi,uhr_adj_4_doxo_bmi,uhr_adj_5_doxo_bmi)
# fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"),5)
# fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"))
# 
# ggplot(fig.doxo.dat, aes(x = bmi,
#                          y = hr, 
#                          fill=cvd, 
#                          color=cvd)) +
#   geom_bar(stat="identity", color="black", 
#            position=position_dodge()) + geom_hline(yintercept=1, color = "red") +
#   #geom_line(size = 1) +
#    geom_text(aes(bmi, hr+0.08, label = format(round(hr,2), nsmall=2), colour = NULL), position=position_dodge(width=0.9), data = fig.doxo.dat) + 
#    ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
#                                                                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                               axis.line.x = element_line(colour = "black"),
#                                                                               axis.line.y = element_line(colour = "black"),
#                                                                               panel.background = element_blank()) + 
#   labs(y = "Adjusted Hazard Ratio", 
#        x = "BMI Quintiles") + scale_y_continuous(breaks=c(0,1,2,3,4))
# dev.off()


# PW_BSA_MOSTELLER
cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_1_doxo_pw_bsa_mosteller,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_1_doxo_pw_bsa_mosteller,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_1_doxo_pw_bsa_mosteller,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_1_doxo_pw_bsa_mosteller,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
hr_adj_1_doxo_pw_bsa_mosteller = cox1[c(2,4),2]
lhr_adj_1_doxo_pw_bsa_mosteller = cox1[c(2,4),3]
uhr_adj_1_doxo_pw_bsa_mosteller = cox1[c(2,4),4]
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_2_doxo_pw_bsa_mosteller,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_2_doxo_pw_bsa_mosteller,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_2_doxo_pw_bsa_mosteller,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_2_doxo_pw_bsa_mosteller,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
hr_adj_2_doxo_pw_bsa_mosteller = cox2[c(2,4),2]
lhr_adj_2_doxo_pw_bsa_mosteller = cox2[c(2,4),3]
uhr_adj_2_doxo_pw_bsa_mosteller = cox2[c(2,4),4]
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_3_doxo_pw_bsa_mosteller,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_3_doxo_pw_bsa_mosteller,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_3_doxo_pw_bsa_mosteller,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_3_doxo_pw_bsa_mosteller,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
hr_adj_3_doxo_pw_bsa_mosteller = cox3[c(2,4),2]
lhr_adj_3_doxo_pw_bsa_mosteller = cox3[c(2,4),3]
uhr_adj_3_doxo_pw_bsa_mosteller = cox3[c(2,4),4]
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_4_doxo_pw_bsa_mosteller,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_4_doxo_pw_bsa_mosteller,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_4_doxo_pw_bsa_mosteller,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_4_doxo_pw_bsa_mosteller,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
hr_adj_4_doxo_pw_bsa_mosteller = cox4[c(2,4),2]
lhr_adj_4_doxo_pw_bsa_mosteller = cox4[c(2,4),3]
uhr_adj_4_doxo_pw_bsa_mosteller = cox4[c(2,4),4]
cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_5_doxo_pw_bsa_mosteller,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_5_doxo_pw_bsa_mosteller,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_5_doxo_pw_bsa_mosteller,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_5_doxo_pw_bsa_mosteller,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
hr_adj_5_doxo_pw_bsa_mosteller = cox5[c(2,4),2]
lhr_adj_5_doxo_pw_bsa_mosteller = cox5[c(2,4),3]
uhr_adj_5_doxo_pw_bsa_mosteller = cox5[c(2,4),4]
cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))



bsa = c(rep("Q1",2),rep("Q2",2),rep("Q3",2),rep("Q4",2),rep("Q5",2))
bsa = factor(bsa, levels = c("Q1","Q2","Q3","Q4","Q5"))
fig.doxo.dat = data.frame(bsa = bsa)
fig.doxo.dat$hr = c(hr_adj_1_doxo_pw_bsa_mosteller,hr_adj_2_doxo_pw_bsa_mosteller,hr_adj_3_doxo_pw_bsa_mosteller,hr_adj_4_doxo_pw_bsa_mosteller,hr_adj_5_doxo_pw_bsa_mosteller)
fig.doxo.dat$lower = c(lhr_adj_1_doxo_pw_bsa_mosteller,lhr_adj_2_doxo_pw_bsa_mosteller,lhr_adj_3_doxo_pw_bsa_mosteller,lhr_adj_4_doxo_pw_bsa_mosteller,lhr_adj_5_doxo_pw_bsa_mosteller)
fig.doxo.dat$upper = c(uhr_adj_1_doxo_pw_bsa_mosteller,uhr_adj_2_doxo_pw_bsa_mosteller,uhr_adj_3_doxo_pw_bsa_mosteller,uhr_adj_4_doxo_pw_bsa_mosteller,uhr_adj_5_doxo_pw_bsa_mosteller)
fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy"),5)
fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy"))

p2 = ggplot(fig.doxo.dat, aes(x = bsa,
                         y = hr, 
                         fill=cvd, 
                         color=cvd)) +
  geom_bar(stat="identity", colour = "white",
           position=position_dodge()) + geom_hline(yintercept=1, color = "red", size = 1) + scale_fill_manual(values = c("steelblue", "gold")) + 
  geom_errorbar(aes(ymin  =lower, 
                    ymax = upper), 
                width = .2, color = "black",
                position=position_dodge(.9), size = 1) + ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
                                                                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                              axis.line.x = element_line(colour = "black"),
                                                                              axis.line.y = element_line(colour = "black"),
                                                                              panel.background = element_blank()) + 
  labs(y = "Adjusted Hazard Ratio", 
       x = "BSA MOSTELLER Quintiles") + coord_cartesian(ylim=c(0,10)) + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10))

  

# png("Figure_PW_BSA_MOSTELLER_quintiles_Doxorubicin_2.png", width = 9, height=5.5, res=300, units = 'in')
# 
# bsa = c(rep("0 - 20%",4),rep("20% - 40%",4),rep("40% - 60%",4),rep("60% - 80%",4),rep("80% - 100%",4))
# bsa = factor(bsa, levels = c("0 - 20%","20% - 40%","40% - 60%","60% - 80%","80% - 100%"))
# fig.doxo.dat = data.frame(bsa = bsa)
# fig.doxo.dat$hr = c(hr_adj_1_doxo_pw_bsa_mosteller,hr_adj_2_doxo_pw_bsa_mosteller,hr_adj_3_doxo_pw_bsa_mosteller,hr_adj_4_doxo_pw_bsa_mosteller,hr_adj_5_doxo_pw_bsa_mosteller)
# fig.doxo.dat$lower = c(lhr_adj_1_doxo_pw_bsa_mosteller,lhr_adj_2_doxo_pw_bsa_mosteller,lhr_adj_3_doxo_pw_bsa_mosteller,lhr_adj_4_doxo_pw_bsa_mosteller,lhr_adj_5_doxo_pw_bsa_mosteller)
# fig.doxo.dat$upper = c(uhr_adj_1_doxo_pw_bsa_mosteller,uhr_adj_2_doxo_pw_bsa_mosteller,uhr_adj_3_doxo_pw_bsa_mosteller,uhr_adj_4_doxo_pw_bsa_mosteller,uhr_adj_5_doxo_pw_bsa_mosteller)
# fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"),5)
# fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"))
# 
# ggplot(fig.doxo.dat, aes(x = bsa,
#                          y = hr, 
#                          fill=cvd, 
#                          color=cvd)) +
#   geom_bar(stat="identity", color="black", 
#            position=position_dodge()) + geom_hline(yintercept=1, color = "red") + 
#   #geom_line(size = 1) +
#   geom_text(aes(bsa, hr+0.09, label = format(round(hr,2),nsmall = 2), colour = NULL), position=position_dodge(width=0.9), data = fig.doxo.dat) + 
#   ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
#                                                                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                               axis.line.x = element_line(colour = "black"),
#                                                                               axis.line.y = element_line(colour = "black"),
#                                                                               panel.background = element_blank()) + 
#   labs(y = "Adjusted Hazard Ratio", 
#        x = "BSA MOSTELLER Quintiles") + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11))
# dev.off()


# # PW_BSA_DUBOIS
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_1_doxo_pw_bsa_dubois,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_1_doxo_pw_bsa_dubois,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_1_doxo_pw_bsa_dubois,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_1_doxo_pw_bsa_dubois,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_1_doxo_pw_bsa_dubois = cox1[c(2,4,6,8),2]
# lhr_adj_1_doxo_pw_bsa_dubois = cox1[c(2,4,6,8),3]
# uhr_adj_1_doxo_pw_bsa_dubois = cox1[c(2,4,6,8),4]
# cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))
# 
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_2_doxo_pw_bsa_dubois,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_2_doxo_pw_bsa_dubois,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_2_doxo_pw_bsa_dubois,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_2_doxo_pw_bsa_dubois,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_2_doxo_pw_bsa_dubois = cox2[c(2,4,6,8),2]
# lhr_adj_2_doxo_pw_bsa_dubois = cox2[c(2,4,6,8),3]
# uhr_adj_2_doxo_pw_bsa_dubois = cox2[c(2,4,6,8),4]
# cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))
# 
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_3_doxo_pw_bsa_dubois,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_3_doxo_pw_bsa_dubois,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_3_doxo_pw_bsa_dubois,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_3_doxo_pw_bsa_dubois,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_3_doxo_pw_bsa_dubois = cox3[c(2,4,6,8),2]
# lhr_adj_3_doxo_pw_bsa_dubois = cox3[c(2,4,6,8),3]
# uhr_adj_3_doxo_pw_bsa_dubois = cox3[c(2,4,6,8),4]
# cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))
# 
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_4_doxo_pw_bsa_dubois,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_4_doxo_pw_bsa_dubois,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_4_doxo_pw_bsa_dubois,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_4_doxo_pw_bsa_dubois,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_4_doxo_pw_bsa_dubois = cox4[c(2,4,6,8),2]
# lhr_adj_4_doxo_pw_bsa_dubois = cox4[c(2,4,6,8),3]
# uhr_adj_4_doxo_pw_bsa_dubois = cox4[c(2,4,6,8),4]
# cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))
# 
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_5_doxo_pw_bsa_dubois,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_5_doxo_pw_bsa_dubois,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_5_doxo_pw_bsa_dubois,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_5_doxo_pw_bsa_dubois,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_5_doxo_pw_bsa_dubois = cox5[c(2,4,6,8),2]
# lhr_adj_5_doxo_pw_bsa_dubois = cox5[c(2,4,6,8),3]
# uhr_adj_5_doxo_pw_bsa_dubois = cox5[c(2,4,6,8),4]
# cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))
# 
# 
# png("Figure_PW_BSA_DUBOIS_quintiles_Doxorubicin.png", width = 9, height=5.5, res=300, units = 'in')
# 
# bsa = c(rep("0 - 20%",4),rep("20% - 40%",4),rep("40% - 60%",4),rep("60% - 80%",4),rep("80% - 100%",4))
# bsa = factor(bsa, levels = c("0 - 20%","20% - 40%","40% - 60%","60% - 80%","80% - 100%"))
# fig.doxo.dat = data.frame(bsa = bsa)
# fig.doxo.dat$hr = c(hr_adj_1_doxo_pw_bsa_dubois,hr_adj_2_doxo_pw_bsa_dubois,hr_adj_3_doxo_pw_bsa_dubois,hr_adj_4_doxo_pw_bsa_dubois,hr_adj_5_doxo_pw_bsa_dubois)
# fig.doxo.dat$lower = c(lhr_adj_1_doxo_pw_bsa_dubois,lhr_adj_2_doxo_pw_bsa_dubois,lhr_adj_3_doxo_pw_bsa_dubois,lhr_adj_4_doxo_pw_bsa_dubois,lhr_adj_5_doxo_pw_bsa_dubois)
# fig.doxo.dat$upper = c(uhr_adj_1_doxo_pw_bsa_dubois,uhr_adj_2_doxo_pw_bsa_dubois,uhr_adj_3_doxo_pw_bsa_dubois,uhr_adj_4_doxo_pw_bsa_dubois,uhr_adj_5_doxo_pw_bsa_dubois)
# fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"),5)
# fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"))
# 
# p3 = ggplot(fig.doxo.dat, aes(x = bsa,
#                          y = hr, 
#                          fill=cvd, 
#                          color=cvd)) +
#   geom_bar(stat="identity", color="black", 
#            position=position_dodge()) + geom_hline(yintercept=1, color = "red") + 
#   #geom_line(size = 1) +
#   geom_errorbar(aes(ymin  =lower, 
#                     ymax = upper), 
#                 width = .2, color = "black",
#                 position=position_dodge(.9)) + ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
#                                                                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                               axis.line.x = element_line(colour = "black"),
#                                                                               axis.line.y = element_line(colour = "black"),
#                                                                               panel.background = element_blank()) + 
#   labs(y = "Adjusted Hazard Ratio", 
#        x = "BSA DUBOIS Quintiles") + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11))
# dev.off()
# 
# 
# png("Figure_PW_BSA_DUBOIS_quintiles_Doxorubicin_2.png", width = 9, height=5.5, res=300, units = 'in')
# 
# bsa = c(rep("0 - 20%",4),rep("20% - 40%",4),rep("40% - 60%",4),rep("60% - 80%",4),rep("80% - 100%",4))
# bsa = factor(bsa, levels = c("0 - 20%","20% - 40%","40% - 60%","60% - 80%","80% - 100%"))
# fig.doxo.dat = data.frame(bsa = bsa)
# fig.doxo.dat$hr = c(hr_adj_1_doxo_pw_bsa_dubois,hr_adj_2_doxo_pw_bsa_dubois,hr_adj_3_doxo_pw_bsa_dubois,hr_adj_4_doxo_pw_bsa_dubois,hr_adj_5_doxo_pw_bsa_dubois)
# fig.doxo.dat$lower = c(lhr_adj_1_doxo_pw_bsa_dubois,lhr_adj_2_doxo_pw_bsa_dubois,lhr_adj_3_doxo_pw_bsa_dubois,lhr_adj_4_doxo_pw_bsa_dubois,lhr_adj_5_doxo_pw_bsa_dubois)
# fig.doxo.dat$upper = c(uhr_adj_1_doxo_pw_bsa_dubois,uhr_adj_2_doxo_pw_bsa_dubois,uhr_adj_3_doxo_pw_bsa_dubois,uhr_adj_4_doxo_pw_bsa_dubois,uhr_adj_5_doxo_pw_bsa_dubois)
# fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"),5)
# fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"))
# 
# ggplot(fig.doxo.dat, aes(x = bsa,
#                          y = hr, 
#                          fill=cvd, 
#                          color=cvd)) +
#   geom_bar(stat="identity", color="black", 
#            position=position_dodge()) + geom_hline(yintercept=1, color = "red") + 
#   #geom_line(size = 1) +
#   geom_text(aes(bsa, hr+0.09, label = format(round(hr,2),nsmall = 2), colour = NULL), position=position_dodge(width=0.9), data = fig.doxo.dat) + 
#   ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
#                                                                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                               axis.line.x = element_line(colour = "black"),
#                                                                               axis.line.y = element_line(colour = "black"),
#                                                                               panel.background = element_blank()) + 
#   labs(y = "Adjusted Hazard Ratio", 
#        x = "BSA DUBOIS Quintiles") + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11))
# dev.off()



# Muscle
cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_1_doxo_muscle,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_1_doxo_muscle,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_1_doxo_muscle,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_1_doxo_muscle,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
hr_adj_1_doxo_muscle = cox1[c(2,4),2]
lhr_adj_1_doxo_muscle = cox1[c(2,4),3]
uhr_adj_1_doxo_muscle = cox1[c(2,4),4]
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_2_doxo_muscle,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_2_doxo_muscle,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_2_doxo_muscle,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_2_doxo_muscle,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
hr_adj_2_doxo_muscle = cox2[c(2,4),2]
lhr_adj_2_doxo_muscle = cox2[c(2,4),3]
uhr_adj_2_doxo_muscle = cox2[c(2,4),4]
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_3_doxo_muscle,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_3_doxo_muscle,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_3_doxo_muscle,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_3_doxo_muscle,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
hr_adj_3_doxo_muscle = cox3[c(2,4),2]
lhr_adj_3_doxo_muscle = cox3[c(2,4),3]
uhr_adj_3_doxo_muscle = cox3[c(2,4),4]
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_4_doxo_muscle,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_4_doxo_muscle,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_4_doxo_muscle,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_4_doxo_muscle,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
hr_adj_4_doxo_muscle = cox4[c(2,4),2]
lhr_adj_4_doxo_muscle = cox4[c(2,4),3]
uhr_adj_4_doxo_muscle = cox4[c(2,4),4]
cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_5_doxo_muscle,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_5_doxo_muscle,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_5_doxo_muscle,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_5_doxo_muscle,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
hr_adj_5_doxo_muscle = cox5[c(2,4),2]
lhr_adj_5_doxo_muscle = cox5[c(2,4),3]
uhr_adj_5_doxo_muscle = cox5[c(2,4),4]
cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))



muscle = c(rep("Q1",2),rep("Q2",2),rep("Q3",2),rep("Q4",2),rep("Q5",2))
muscle = factor(muscle, levels = c("Q1","Q2","Q3","Q4","Q5"))
fig.doxo.dat = data.frame(muscle = muscle)
fig.doxo.dat$hr = c(hr_adj_1_doxo_muscle,hr_adj_2_doxo_muscle,hr_adj_3_doxo_muscle,hr_adj_4_doxo_muscle,hr_adj_5_doxo_muscle)
fig.doxo.dat$lower = c(lhr_adj_1_doxo_muscle,lhr_adj_2_doxo_muscle,lhr_adj_3_doxo_muscle,lhr_adj_4_doxo_muscle,lhr_adj_5_doxo_muscle)
fig.doxo.dat$upper = c(uhr_adj_1_doxo_muscle,uhr_adj_2_doxo_muscle,uhr_adj_3_doxo_muscle,uhr_adj_4_doxo_muscle,uhr_adj_5_doxo_muscle)
fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy"),5)
fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy"))
fig.doxo.dat$upper[fig.doxo.dat$upper == "Inf"] = 100

p4 = ggplot(fig.doxo.dat, aes(x = muscle,
                         y = hr, 
                         fill=cvd, 
                         color=cvd)) +
  geom_bar(stat="identity", color="white", 
           position=position_dodge()) + geom_hline(yintercept=1, color = "red", size = 1) + scale_fill_manual(values = c("steelblue", "gold")) + 
  #geom_line(size = 1) +
  geom_errorbar(aes(ymin  =lower, 
                    ymax = upper), 
                width = .2, color = "black",
                position=position_dodge(.9), size = 1) + ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
                                                                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                              axis.line.x = element_line(colour = "black"),
                                                                              axis.line.y = element_line(colour = "black"),
                                                                              panel.background = element_blank()) + 
  labs(y = "Adjusted Hazard Ratio", 
       x = expression(paste("Muscle Area (",cm^2,") Quintiles"))) + 
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)) + coord_cartesian(ylim=c(0,10)) 


# png("Figure_Muscle_Area_quintiles_Doxorubicin_2.png", width = 9, height=5.5, res=300, units = 'in')
# 
# muscle = c(rep("0 - 20%",4),rep("20% - 40%",4),rep("40% - 60%",4),rep("60% - 80%",4),rep("80% - 100%",4))
# muscle = factor(muscle, levels = c("0 - 20%","20% - 40%","40% - 60%","60% - 80%","80% - 100%"))
# fig.doxo.dat = data.frame(muscle = muscle)
# fig.doxo.dat$hr = c(hr_adj_1_doxo_muscle,hr_adj_2_doxo_muscle,hr_adj_3_doxo_muscle,hr_adj_4_doxo_muscle,hr_adj_5_doxo_muscle)
# fig.doxo.dat$lower = c(lhr_adj_1_doxo_muscle,lhr_adj_2_doxo_muscle,lhr_adj_3_doxo_muscle,lhr_adj_4_doxo_muscle,lhr_adj_5_doxo_muscle)
# fig.doxo.dat$upper = c(uhr_adj_1_doxo_muscle,uhr_adj_2_doxo_muscle,uhr_adj_3_doxo_muscle,uhr_adj_4_doxo_muscle,uhr_adj_5_doxo_muscle)
# fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"),5)
# fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"))
# 
# ggplot(fig.doxo.dat, aes(x = muscle,
#                          y = hr, 
#                          fill=cvd, 
#                          color=cvd)) +
#   geom_bar(stat="identity", color="black", 
#            position=position_dodge()) + geom_hline(yintercept=1, color = "red") + 
#   #geom_line(size = 1) +
#   geom_text(aes(muscle, hr+0.07, label = format(round(hr,2),nsmall = 2), colour = NULL), position=position_dodge(width=0.9), data = fig.doxo.dat) + 
#   ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
#                                                                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                               axis.line.x = element_line(colour = "black"),
#                                                                               axis.line.y = element_line(colour = "black"),
#                                                                               panel.background = element_blank()) + 
#   labs(y = "Adjusted Hazard Ratio", 
#        x = expression(paste("Muscle Area (",cm^2,") Quintiles"))) + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13))
# dev.off()



# VF
cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_1_doxo_vf,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_1_doxo_vf,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_1_doxo_vf,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_1_doxo_vf,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
hr_adj_1_doxo_vf = cox1[c(2,4),2]
lhr_adj_1_doxo_vf = cox1[c(2,4),3]
uhr_adj_1_doxo_vf = cox1[c(2,4),4]
cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_2_doxo_vf,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_2_doxo_vf,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_2_doxo_vf,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_2_doxo_vf,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
hr_adj_2_doxo_vf = cox2[c(2,4),2]
lhr_adj_2_doxo_vf = cox2[c(2,4),3]
uhr_adj_2_doxo_vf = cox2[c(2,4),4]
cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_3_doxo_vf,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_3_doxo_vf,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_3_doxo_vf,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_3_doxo_vf,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
hr_adj_3_doxo_vf = cox3[c(2,4),2]
lhr_adj_3_doxo_vf = cox3[c(2,4),3]
uhr_adj_3_doxo_vf = cox3[c(2,4),4]
cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_4_doxo_vf,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_4_doxo_vf,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_4_doxo_vf,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_4_doxo_vf,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
hr_adj_4_doxo_vf = cox4[c(2,4),2]
lhr_adj_4_doxo_vf = cox4[c(2,4),3]
uhr_adj_4_doxo_vf = cox4[c(2,4),4]
cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))

cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_5_doxo_vf,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_5_doxo_vf,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_5_doxo_vf,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_5_doxo_vf,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,])
cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
hr_adj_5_doxo_vf = cox5[c(2,4),2]
lhr_adj_5_doxo_vf = cox5[c(2,4),3]
uhr_adj_5_doxo_vf = cox5[c(2,4),4]
cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))



vf = c(rep("Q1",2),rep("Q2",2),rep("Q3",2),rep("Q4",2),rep("Q5",2))
vf= factor(vf, levels = c("Q1","Q2","Q3","Q4","Q5"))
fig.doxo.dat = data.frame(vf = vf)
fig.doxo.dat$hr = c(hr_adj_1_doxo_vf,hr_adj_2_doxo_vf,hr_adj_3_doxo_vf,hr_adj_4_doxo_vf,hr_adj_5_doxo_vf)
fig.doxo.dat$lower = c(lhr_adj_1_doxo_vf,lhr_adj_2_doxo_vf,lhr_adj_3_doxo_vf,lhr_adj_4_doxo_vf,lhr_adj_5_doxo_vf)
fig.doxo.dat$upper = c(uhr_adj_1_doxo_vf,uhr_adj_2_doxo_vf,uhr_adj_3_doxo_vf,uhr_adj_4_doxo_vf,uhr_adj_5_doxo_vf)
fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy"),5)
fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy"))
fig.doxo.dat$upper[fig.doxo.dat$upper == "Inf"] = 100

p5 = ggplot(fig.doxo.dat, aes(x = vf,
                         y = hr, 
                         fill=cvd, 
                         color=cvd)) +
  geom_bar(stat="identity", color="white", 
           position=position_dodge()) + geom_hline(yintercept=1, color = "red", size = 1) + scale_fill_manual(values = c("steelblue", "gold")) + 
  #geom_line(size = 1) +
  geom_errorbar(aes(ymin  =lower, 
                    ymax = upper), 
                width = .2, color = "black",
                position=position_dodge(.9), size = 1) + ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
                                                                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                              axis.line.x = element_line(colour = "black"),
                                                                              axis.line.y = element_line(colour = "black"),
                                                                              panel.background = element_blank()) + 
  labs(y = "Adjusted Hazard Ratio", 
       x = expression(paste("Visceral Fat (",cm^2, ") Quintiles"))) + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) +
  coord_cartesian(ylim=c(0,10))

# png("Figure_Visceral_Fat_quintiles_Doxorubicin_2.png", width = 9, height=5.5, res=300, units = 'in')
# 
# vf = c(rep("0 - 20%",4),rep("20% - 40%",4),rep("40% - 60%",4),rep("60% - 80%",4),rep("80% - 100%",4))
# vf= factor(vf, levels = c("0 - 20%","20% - 40%","40% - 60%","60% - 80%","80% - 100%"))
# fig.doxo.dat = data.frame(vf = vf)
# fig.doxo.dat$hr = c(hr_adj_1_doxo_vf,hr_adj_2_doxo_vf,hr_adj_3_doxo_vf,hr_adj_4_doxo_vf,hr_adj_5_doxo_vf)
# fig.doxo.dat$lower = c(lhr_adj_1_doxo_vf,lhr_adj_2_doxo_vf,lhr_adj_3_doxo_vf,lhr_adj_4_doxo_vf,lhr_adj_5_doxo_vf)
# fig.doxo.dat$upper = c(uhr_adj_1_doxo_vf,uhr_adj_2_doxo_vf,uhr_adj_3_doxo_vf,uhr_adj_4_doxo_vf,uhr_adj_5_doxo_vf)
# fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"),5)
# fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"))
# 
# ggplot(fig.doxo.dat, aes(x = vf,
#                          y = hr, 
#                          fill=cvd, 
#                          color=cvd)) +
#   geom_bar(stat="identity", color="black", 
#            position=position_dodge()) + geom_hline(yintercept=1, color = "red") + 
#   #geom_line(size = 1) +
#   geom_text(aes(vf, hr+0.10, label = format(round(hr,2),nsmall = 2), colour = NULL), position=position_dodge(width=0.9), data = fig.doxo.dat) + 
#   ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
#                                                                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                               axis.line.x = element_line(colour = "black"),
#                                                                               axis.line.y = element_line(colour = "black"),
#                                                                               panel.background = element_blank()) + 
#   labs(y = "Adjusted Hazard Ratio", 
#        x = expression(paste("Visceral Fat (",cm^2, ") Quintiles"))) + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))
# dev.off()


# SF
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_1_doxo_sf,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_1_doxo_sf,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_1_doxo_sf,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_1_doxo_sf,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_1_doxo_sf = cox1[c(2,4,6,8),2]
# lhr_adj_1_doxo_sf = cox1[c(2,4,6,8),3]
# uhr_adj_1_doxo_sf = cox1[c(2,4,6,8),4]
# cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))
# 
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_2_doxo_sf,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_2_doxo_sf,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_2_doxo_sf,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_2_doxo_sf,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_2_doxo_sf = cox2[c(2,4,6,8),2]
# lhr_adj_2_doxo_sf = cox2[c(2,4,6,8),3]
# uhr_adj_2_doxo_sf = cox2[c(2,4,6,8),4]
# cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))
# 
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_3_doxo_sf,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_3_doxo_sf,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_3_doxo_sf,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_3_doxo_sf,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_3_doxo_sf = cox3[c(2,4,6,8),2]
# lhr_adj_3_doxo_sf = cox3[c(2,4,6,8),3]
# uhr_adj_3_doxo_sf = cox3[c(2,4,6,8),4]
# cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))
# 
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_4_doxo_sf,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_4_doxo_sf,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_4_doxo_sf,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_4_doxo_sf,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_4_doxo_sf = cox4[c(2,4,6,8),2]
# lhr_adj_4_doxo_sf = cox4[c(2,4,6,8),3]
# uhr_adj_4_doxo_sf = cox4[c(2,4,6,8),4]
# cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))
# 
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_5_doxo_sf,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_5_doxo_sf,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_5_doxo_sf,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_5_doxo_sf,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_5_doxo_sf = cox5[c(2,4,6,8),2]
# lhr_adj_5_doxo_sf = cox5[c(2,4,6,8),3]
# uhr_adj_5_doxo_sf = cox5[c(2,4,6,8),4]
# cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))
# 
# 
# png("Figure_Subcutaneous_Fat_quintiles_Doxorubicin.png", width = 9, height=5.5, res=300, units = 'in')
# 
# sf = c(rep("0 - 20%",4),rep("20% - 40%",4),rep("40% - 60%",4),rep("60% - 80%",4),rep("80% - 100%",4))
# sf= factor(sf, levels = c("0 - 20%","20% - 40%","40% - 60%","60% - 80%","80% - 100%"))
# fig.doxo.dat = data.frame(sf = sf)
# fig.doxo.dat$hr = c(hr_adj_1_doxo_sf,hr_adj_2_doxo_sf,hr_adj_3_doxo_sf,hr_adj_4_doxo_sf,hr_adj_5_doxo_sf)
# fig.doxo.dat$lower = c(lhr_adj_1_doxo_sf,lhr_adj_2_doxo_sf,lhr_adj_3_doxo_sf,lhr_adj_4_doxo_sf,lhr_adj_5_doxo_sf)
# fig.doxo.dat$upper = c(uhr_adj_1_doxo_sf,uhr_adj_2_doxo_sf,uhr_adj_3_doxo_sf,uhr_adj_4_doxo_sf,uhr_adj_5_doxo_sf)
# fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"),5)
# fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"))
# 
# p6 = ggplot(fig.doxo.dat, aes(x = sf,
#                          y = hr, 
#                          fill=cvd, 
#                          color=cvd)) +
#   geom_bar(stat="identity", color="black", 
#            position=position_dodge()) + geom_hline(yintercept=1, color = "red") + 
#   #geom_line(size = 1) +
#   geom_errorbar(aes(ymin  =lower, 
#                     ymax = upper), 
#                 width = .2, color = "black",
#                 position=position_dodge(.9)) + ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
#                                                                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                               axis.line.x = element_line(colour = "black"),
#                                                                               axis.line.y = element_line(colour = "black"),
#                                                                               panel.background = element_blank()) +  
#   labs(y = "Adjusted Hazard Ratio", 
#        x = expression(paste("Subcutaneous Fat (",cm^2, ") Quintiles"))) + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23))
# dev.off()
# 
# png("Figure_Subcutaneous_Fat_quintiles_Doxorubicin_2.png", width = 9, height=5.5, res=300, units = 'in')
# 
# sf = c(rep("0 - 20%",4),rep("20% - 40%",4),rep("40% - 60%",4),rep("60% - 80%",4),rep("80% - 100%",4))
# sf= factor(sf, levels = c("0 - 20%","20% - 40%","40% - 60%","60% - 80%","80% - 100%"))
# fig.doxo.dat = data.frame(sf = sf)
# fig.doxo.dat$hr = c(hr_adj_1_doxo_sf,hr_adj_2_doxo_sf,hr_adj_3_doxo_sf,hr_adj_4_doxo_sf,hr_adj_5_doxo_sf)
# fig.doxo.dat$lower = c(lhr_adj_1_doxo_sf,lhr_adj_2_doxo_sf,lhr_adj_3_doxo_sf,lhr_adj_4_doxo_sf,lhr_adj_5_doxo_sf)
# fig.doxo.dat$upper = c(uhr_adj_1_doxo_sf,uhr_adj_2_doxo_sf,uhr_adj_3_doxo_sf,uhr_adj_4_doxo_sf,uhr_adj_5_doxo_sf)
# fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"),5)
# fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"))
# 
# ggplot(fig.doxo.dat, aes(x = sf,
#                          y = hr, 
#                          fill=cvd, 
#                          color=cvd)) +
#   geom_bar(stat="identity", color="black", 
#            position=position_dodge()) + geom_hline(yintercept=1, color = "red") + 
#   #geom_line(size = 1) +
#   geom_text(aes(sf, hr+0.10, label = format(round(hr,2),nsmall = 2), colour = NULL), position=position_dodge(width=0.9), data = fig.doxo.dat) + 
#   ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
#                                                                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                               axis.line.x = element_line(colour = "black"),
#                                                                               axis.line.y = element_line(colour = "black"),
#                                                                               panel.background = element_blank()) +  
#   labs(y = "Adjusted Hazard Ratio", 
#        x = expression(paste("Subcutaneous Fat (",cm^2, ") Quintiles"))) + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23))
# dev.off()
# 
# # Muscle radiodensity, HU
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_1_doxo_mmuscle,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_1_doxo_mmuscle,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_1_doxo_mmuscle,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_1_doxo_mmuscle,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox1 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox1[,-1] <- lapply(cox1[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_1_doxo_mmuscle = cox1[c(2,4,6,8),2]
# lhr_adj_1_doxo_mmuscle = cox1[c(2,4,6,8),3]
# uhr_adj_1_doxo_mmuscle = cox1[c(2,4,6,8),4]
# cox1 <- cbind(round(cox1[,2],2),paste0("(",round(cox1[,3],2),", ", round(cox1[,4],2),")"))
# 
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_2_doxo_mmuscle,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_2_doxo_mmuscle,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_2_doxo_mmuscle,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_2_doxo_mmuscle,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox2 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox2[,-1] <- lapply(cox2[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_2_doxo_mmuscle = cox2[c(2,4,6,8),2]
# lhr_adj_2_doxo_mmuscle = cox2[c(2,4,6,8),3]
# uhr_adj_2_doxo_mmuscle = cox2[c(2,4,6,8),4]
# cox2 <- cbind(round(cox2[,2],2),paste0("(",round(cox2[,3],2),", ", round(cox2[,4],2),")"))
# 
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_3_doxo_mmuscle,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_3_doxo_mmuscle,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_3_doxo_mmuscle,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_3_doxo_mmuscle,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox3 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox3[,-1] <- lapply(cox3[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_3_doxo_mmuscle = cox3[c(2,4,6,8),2]
# lhr_adj_3_doxo_mmuscle = cox3[c(2,4,6,8),3]
# uhr_adj_3_doxo_mmuscle = cox3[c(2,4,6,8),4]
# cox3 <- cbind(round(cox3[,2],2),paste0("(",round(cox3[,3],2),", ", round(cox3[,4],2),")"))
# 
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_4_doxo_mmuscle,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_4_doxo_mmuscle,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_4_doxo_mmuscle,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_4_doxo_mmuscle,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox4 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox4[,-1] <- lapply(cox4[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_4_doxo_mmuscle = cox4[c(2,4,6,8),2]
# lhr_adj_4_doxo_mmuscle = cox4[c(2,4,6,8),3]
# uhr_adj_4_doxo_mmuscle = cox4[c(2,4,6,8),4]
# cox4 <- cbind(round(cox4[,2],2),paste0("(",round(cox4[,3],2),", ", round(cox4[,4],2),")"))
# 
# cox_isch <- coxtab2(phs.aim2.data.doxo[a2_ihd_ef_5_doxo_mmuscle,],'ischemic_heart_disease_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_stroke <- coxtab2(phs.aim2.data.doxo[a2_stroke_ef_5_doxo_mmuscle,],'stroke_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_hf_cm <- coxtab2(phs.aim2.data.doxo[a2_hf_cm_ef_5_doxo_mmuscle,],'heart_failure_cardiomyopathy_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox_cvdcombo <- coxtab2(phs.aim2.data.doxo[a2_combo_ef_5_doxo_mmuscle,],'cvdcombo_grp_inc',prevcvd[c(-2,-6,-7)],"doxo")
# cox5 <- rbind(cox_isch[1:2,],cox_hf_cm[1:2,], cox_stroke[1:2,], 
#               cox_cvdcombo[1:2,])
# cox5[,-1] <- lapply(cox5[,-1], function(x) as.numeric(as.character(x)))
# hr_adj_5_doxo_mmuscle = cox5[c(2,4,6,8),2]
# lhr_adj_5_doxo_mmuscle = cox5[c(2,4,6,8),3]
# uhr_adj_5_doxo_mmuscle = cox5[c(2,4,6,8),4]
# cox5 <- cbind(round(cox5[,2],2),paste0("(",round(cox5[,3],2),", ", round(cox5[,4],2),")"))
# 
# 
# png("Figure_Muscle_Radiodensity_quintiles_Doxorubicin.png", width = 9, height=5.5, res=300, units = 'in')
# 
# mmuscle = c(rep("0 - 20%",4),rep("20% - 40%",4),rep("40% - 60%",4),rep("60% - 80%",4),rep("80% - 100%",4))
# mmuscle = factor(mmuscle , levels = c("0 - 20%","20% - 40%","40% - 60%","60% - 80%","80% - 100%"))
# fig.doxo.dat = data.frame(mmuscle  = mmuscle )
# fig.doxo.dat$hr = c(hr_adj_1_doxo_mmuscle,hr_adj_2_doxo_mmuscle,hr_adj_3_doxo_mmuscle,hr_adj_4_doxo_mmuscle,hr_adj_5_doxo_mmuscle)
# fig.doxo.dat$lower = c(lhr_adj_1_doxo_mmuscle,lhr_adj_2_doxo_mmuscle,lhr_adj_3_doxo_mmuscle,lhr_adj_4_doxo_mmuscle,lhr_adj_5_doxo_mmuscle)
# fig.doxo.dat$upper = c(uhr_adj_1_doxo_mmuscle,uhr_adj_2_doxo_mmuscle,uhr_adj_3_doxo_mmuscle,uhr_adj_4_doxo_mmuscle,uhr_adj_5_doxo_mmuscle)
# fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"),5)
# fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"))
# 
# p7 = ggplot(fig.doxo.dat, aes(x = mmuscle,
#                          y = hr, 
#                          fill=cvd, 
#                          color=cvd)) +
#   geom_bar(stat="identity", color="black", 
#            position=position_dodge()) + geom_hline(yintercept=1, color = "red") + 
#   #geom_line(size = 1) +
#   geom_errorbar(aes(ymin  =lower, 
#                     ymax = upper), 
#                 width = .2, color = "black",
#                 position=position_dodge(.9)) + ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
#                                                                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                               axis.line.x = element_line(colour = "black"),
#                                                                               axis.line.y = element_line(colour = "black"),
#                                                                               panel.background = element_blank()) +  
#   labs(y = "Adjusted Hazard Ratio", 
#        x = "Muscle Radiodensity (HU) Quintiles") + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32))
# dev.off()
# 
# png("Figure_Muscle_Radiodensity_quintiles_Doxorubicin_2.png", width = 9, height=5.5, res=300, units = 'in')
# 
# mmuscle = c(rep("0 - 20%",4),rep("20% - 40%",4),rep("40% - 60%",4),rep("60% - 80%",4),rep("80% - 100%",4))
# mmuscle = factor(mmuscle , levels = c("0 - 20%","20% - 40%","40% - 60%","60% - 80%","80% - 100%"))
# fig.doxo.dat = data.frame(mmuscle  = mmuscle )
# fig.doxo.dat$hr = c(hr_adj_1_doxo_mmuscle,hr_adj_2_doxo_mmuscle,hr_adj_3_doxo_mmuscle,hr_adj_4_doxo_mmuscle,hr_adj_5_doxo_mmuscle)
# fig.doxo.dat$lower = c(lhr_adj_1_doxo_mmuscle,lhr_adj_2_doxo_mmuscle,lhr_adj_3_doxo_mmuscle,lhr_adj_4_doxo_mmuscle,lhr_adj_5_doxo_mmuscle)
# fig.doxo.dat$upper = c(uhr_adj_1_doxo_mmuscle,uhr_adj_2_doxo_mmuscle,uhr_adj_3_doxo_mmuscle,uhr_adj_4_doxo_mmuscle,uhr_adj_5_doxo_mmuscle)
# fig.doxo.dat$cvd = rep(c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"),5)
# fig.doxo.dat$cvd = factor(fig.doxo.dat$cvd, levels = c("Ischemic Heart Disease","Heart Failure and Cardiomyopathy","Stroke","All Primary CVD Outcomes"))
# 
# ggplot(fig.doxo.dat, aes(x = mmuscle,
#                          y = hr, 
#                          fill=cvd, 
#                          color=cvd)) +
#   geom_bar(stat="identity", color="black", 
#            position=position_dodge()) + geom_hline(yintercept=1, color = "red") + 
#   #geom_line(size = 1) +
#   geom_text(aes(mmuscle, hr+0.10, label = format(round(hr,2),nsmall = 2), colour = NULL), position=position_dodge(width=0.9), data = fig.doxo.dat) + 
#   ggtitle("Doxorubicin") + theme(legend.position = "bottom",legend.title=element_blank(),
#                                                                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                               axis.line.x = element_line(colour = "black"),
#                                                                               axis.line.y = element_line(colour = "black"),
#                                                                               panel.background = element_blank()) +  
#   labs(y = "Adjusted Hazard Ratio", 
#        x = "Muscle Radiodensity (HU) Quintiles") + scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32))
# dev.off()


png("Figure_Barplots_BMI.png", width = 18, height=12, res=300, units = 'in')
p0
dev.off()

png("Figure_Barplots_BMI_2.png", width = 12, height=12, res=300, units = 'in')
p0_2
dev.off()

png("Figure_Barplots_BMI_3.png", width = 12, height=12, res=300, units = 'in')
p0_2
dev.off()

png("Figure_Barplots_Others.png", width = 24, height=12, res=300, units = 'in')
ggarrange(p2, p4, p5,
          ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom")
dev.off()


# additional table
tab = function(var, dat, grp){
  var_cat = paste0(var,"_tert")
  tb = table(dat[!is.na(dat[,var]),grp])
  tb.prop = prop.table(tb)*100
  m = aggregate(x = dat[,var],
                by = list(dat[,grp]),FUN = mean,na.rm = T)[,2]
  sd = aggregate(x = dat[,var],
                by = list(dat[,grp]),FUN = sd,na.rm = T)[,2]
  
  mod = summary(aov(dat[,var] ~ dat[,grp]))
  
  sum = cbind(tb[1], paste0(round(tb.prop[1],2),"%"), round(m[1],1), round(sd[1],1), '', 
              tb[2], paste0(round(tb.prop[2],2),"%"), round(m[2],1), round(sd[2],1), '',
              tb[3], paste0(round(tb.prop[3],2),"%"), round(m[3],1), round(sd[3],1), '', mod[[1]][["Pr(>F)"]][1])
  sum
}

tab("sf", phs.aim2.data, "pw_bsa_mosteller_tert")
tab("vf", phs.aim2.data, "pw_bsa_mosteller_tert")
tab("muscle", phs.aim2.data, "pw_bsa_mosteller_tert")
tab("mmuscle", phs.aim2.data, "pw_bsa_mosteller_tert")
tab("bmi_pw", phs.aim2.data, "pw_bsa_mosteller_tert")
tab("vf_to_muscle", phs.aim2.data, "pw_bsa_mosteller_tert")

# Additional figures
# Mosteller BSA and Visceral fat
put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.07,0.975),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02), 
                     bottomcenter = c(0.5525, 0.02), 
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=T, ...)
}

mod <- lm(phs.aim2.data$vf ~ phs.aim2.data$pw_bsa_mosteller)
mod.sum <- summary(mod)
R2 <- mod.sum$r.squared
R2.exp <- expression(paste(" ",R^2 ,"= 0.512"))
quin1 = table(phs.aim2.data$vf_quin)
quin2 = table(phs.aim2.data$pw_bsa_mosteller_quin)
names(quin1) = c("Q1","Q2","Q3","Q4","Q5")
names(quin2) = c("Q1","Q2","Q3","Q4","Q5")


png("Mosteller_BSA_Visceral_fat.png", width = 7, height=6, res=300, units = 'in')
nf = layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE),
       widths = c(2,2), heights = c(1.5,1.5), TRUE)
layout.show(nf)
par(mar = c(4,5,0,2))
plot(phs.aim2.data$pw_bsa_mosteller,phs.aim2.data$vf,xlab = expression(bold(paste("Body surface area (", cm^2,")"))), 
     ylab = expression(bold(paste("Visceral fat (", cm^2,")"))),
     xlim = c(3,6))
abline(mod, lwd = 3)
text(x = 3, y = 400, 
     label = R2.exp,
     pos = 4)
put.fig.letter(label="A", location="topleft", font=2)
par(mar = c(4,5,0,4))
barplot(quin1, ylab = expression(bold("5-FU/VF")), xaxt = "n")
put.fig.letter(label="B", location="topleft", font=2)
par(mar = c(4,5,0,4))
barplot(quin2, ylab = expression(bold("5-FU/BSA")), col = "black", xlab = expression(bold("VF/BSA, Quintiles")), axis.lty	 = 1)
dev.off()


# Mosteller BSA and Muscle area
mod <- lm(phs.aim2.data$muscle ~ phs.aim2.data$pw_bsa_mosteller)
mod.sum <- summary(mod)
R2 <- mod.sum$r.squared
R2.exp <- expression(paste(" ",R^2 ,"= 0.473"))
quin1 = table(phs.aim2.data$muscle_quin)
quin2 = table(phs.aim2.data$pw_bsa_mosteller_quin)
names(quin1) = c("Q1","Q2","Q3","Q4","Q5")
names(quin2) = c("Q1","Q2","Q3","Q4","Q5")


png("Mosteller_BSA_Muscle.png", width = 7, height=6, res=300, units = 'in')
nf = layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE),
            widths = c(2,2), heights = c(1.5,1.5), TRUE)
layout.show(nf)
par(mar = c(4,5,0,2))
plot(phs.aim2.data$pw_bsa_mosteller,phs.aim2.data$muscle,xlab = expression(bold(paste("Body surface area (", cm^2,")"))), 
     ylab = expression(bold(paste("Muscle area (", cm^2,")"))),
     xlim = c(3,6))
abline(mod, lwd = 3)
text(x = 3, y = 200, 
     label = R2.exp,
     pos = 4)
put.fig.letter(label="A", location="topleft", font=2)
par(mar = c(4,5,0,4))
barplot(quin1, ylab = expression(bold("5-FU/Muscle area")), xaxt = "n")
put.fig.letter(label="B", location="topleft", font=2)
par(mar = c(4,5,0,4))
barplot(quin2, ylab = expression(bold("5-FU/BSA")), col = "black", xlab = expression(bold("Muscle area/BSA, Quintiles")), axis.lty	 = 1)
dev.off()


# Mosteller BSA and ratio of visceral fat / muscle area
mod <- lm(phs.aim2.data$vf_to_muscle ~ phs.aim2.data$pw_bsa_mosteller)
mod.sum <- summary(mod)
R2 <- mod.sum$r.squared
R2.exp <- expression(paste(" ",R^2 ,"= 0.372"))
quin1 = table(phs.aim2.data$vf_to_muscle_quin)
quin2 = table(phs.aim2.data$pw_bsa_mosteller_quin)
names(quin1) = c("Q1","Q2","Q3","Q4","Q5")
names(quin2) = c("Q1","Q2","Q3","Q4","Q5")


png("Mosteller_BSA_VF_to_Muscle_Ratio.png", width = 7, height=6, res=300, units = 'in')
nf = layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE),
            widths = c(2,2), heights = c(1.5,1.5), TRUE)
layout.show(nf)
par(mar = c(4,5,0,2))
plot(phs.aim2.data$pw_bsa_mosteller,phs.aim2.data$vf_to_muscle,xlab = expression(bold(paste("Body surface area (", cm^2,")"))), 
     ylab = expression(bold("Visceral fat/Muscle area")),
     xlim = c(3,6))
abline(mod, lwd = 3)
text(x = 3, y = 2.5, 
     label = R2.exp,
     pos = 4)
put.fig.letter(label="A", location="topleft", font=2)
par(mar = c(4,5,0,4))
barplot(quin1, ylab = expression(bold("5-FU/(VF/Muscle)")), xaxt = "n")
put.fig.letter(label="B", location="topleft", font=2)
par(mar = c(4,5,0,4))
barplot(quin2, ylab = expression(bold("5-FU/BSA")), col = "black", xlab = expression(bold("(VF/Muscle)/BSA, Quintiles")), axis.lty	 = 1)
dev.off()

# 
nf = layout(matrix(c(1,2,1,3), 1, 2, byrow = TRUE),
            widths = c(2.5,2.5), heights = c(3.5,3.5), TRUE)
layout.show(nf)
par(mar = c(4,5,0,2))
R2.exp <- expression(paste(" ",R^2 ,"= 0.492"))
mod <- lm(phs.aim2.data$fat_to_muscle ~ phs.aim2.data$pw_bsa_mosteller)
plot(phs.aim2.data$pw_bsa_mosteller,phs.aim2.data$fat_to_muscle,xlab = expression(bold(paste("Body surface area (", cm^2,")"))), 
     ylab = expression(bold("Total fat/muscle ratio")),
     xlim = c(3,6))
abline(mod, lwd = 3)
text(x = 3, y = 6.5, 
     label = R2.exp,
     pos = 4)
par(mar = c(4,5,0,2))
R2.exp <- expression(paste(" ",R^2 ,"= 0.473"))
mod <- lm(phs.aim2.data$muscle ~ phs.aim2.data$pw_bsa_mosteller)
plot(phs.aim2.data$pw_bsa_mosteller,phs.aim2.data$muscle,xlab = expression(bold(paste("Body surface area (", cm^2,")"))), 
     ylab = expression(bold(paste("Muscle area (", cm^2,")"))),
     xlim = c(3,6))
abline(mod, lwd = 3)
text(x = 3, y = 214, 
     label = R2.exp,
     pos = 4)
title('Figure 1: Variablity in visceral fat and muscle area, by body surface area, in \n breast cacner patients at the time of diagnosis (n = 455)', 
      outer = T, line = -3)


############ Model 3: IPW KM Curves and Cox Regression ############################
# Chemotherapy
# Anthracycline
# Primary Outcomes
dat_ihd_anthra$LtAt.data$Y.tplus1 = NULL

shift = function(x, n){
  c(x[-(seq(n))], rep(NA, n))
}

id_uni = unique(dat_ihd_anthra$LtAt.data$ID)
for (i in 1:length(id_uni)){
  dat_ihd_anthra$LtAt.data$Y.tplus1[dat_ihd_anthra$LtAt.data$ID == id_uni[i]] = shift(dat_ihd_anthra$LtAt.data$outcome[dat_ihd_anthra$LtAt.data$ID == id_uni[i]], 1)
}

dat_ihd_anthra$LtAt.data[, ("exposure.set1") := 1L]
dat_ihd_anthra$LtAt.data[, ("exposure.set0") := 0L]

OData = importData(dat_ihd_anthra$LtAt.data, ID = "ID", t = "intnum", covars = covar_chemo, CENS = "censor", TRT = "exposure", OUTCOME = "Y.tplus1")

get_data(OData)[, ("exposure.set1") := 1L]
get_data(OData)[, ("exposure.set0") := 0L]

gform_CENS = "censor ~ I.agegrp +  I.ajcc_stage + I.race + I.bmicat1 + I.smok + I.charlson + I.diab_bl + I.dyslipid_bl + I.htn_bl + I.edu_cat + I.income_cat + I.menop + I.horm_yn + I.rad_tx_yn"
gform_TRT = "exposure ~ I.agegrp +  I.ajcc_stage + I.race + I.bmicat1 + I.smok + I.charlson + I.diab_bl + I.dyslipid_bl + I.htn_bl + I.edu_cat + I.income_cat + I.menop + I.horm_yn + I.rad_tx_yn"
stratify_CENS = list(censor=c("intnum < 12", "intnum == 12"))

OData = fitPropensity(OData,
                      gform_CENS = gform_CENS,
                      gform_TRT = gform_TRT,
                      stratify_CENS = stratify_CENS)


AKME.St.1 = getIPWeights(OData, intervened_TRT = "exposure.set1") %>%
  survNPMSM(OData) %$%
  estimates


head(AKME.St.1[])

IPW.St.1 = getIPWeights(OData, intervened_TRT = "exposure.set1") %>%
  directIPW(OData) %$%
  estimates

head(IPW.St.1[])

wts.DT.1 = getIPWeights(OData = OData, intervened_TRT = "exposure.set1", rule_name = "exposure1")
wts.DT.0 = getIPWeights(OData = OData, intervened_TRT = "exposure.set0", rule_name = "exposure0")
survMSM_res = survMSM(list(wts.DT.1, wts.DT.0), OData, tbreaks = c(0,1,2))

head(survMSM_res[["exposure0"]][["estimates"]])

head(survMSM_res[["exposure1"]][["estimates"]])


###### Estimate the exposure / treatment propensity model with the above defined Super Learner (lrn_sl)
library("sl3")
lrn_xgb <- Lrnr_xgboost$new(nrounds = 5)
lrn_glm <- Lrnr_glm_fast$new()
lrn_glm2 <- Lrnr_glm_fast$new(covariates = c("CVD"))
lrn_glmnet <- Lrnr_glmnet$new(nlambda = 5, family = "binomial")
## Stack the above candidates:
lrn_stack <- Stack$new(lrn_xgb, lrn_glm, lrn_glm2, lrn_glmnet)

lrn_sl <- Lrnr_sl$new(learners = lrn_stack, metalearner = Lrnr_solnp$new())

OData <- fitPropensity(OData,
                       gform_CENS = gform_CENS,
                       gform_TRT = gform_TRT,
                       models_TRT = lrn_sl,
                       stratify_CENS = stratify_CENS)



############ Model 4: Unadjusted Cox Regression with Machine Learning #############






########### Additional analyses for grant 12/21/2020


























