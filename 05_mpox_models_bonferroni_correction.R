#Load in libraries, and install as needed
library(devtools)
#install_github('trendecon/trendecon')
library(trendecon)
library(gtrendsR)
library(tsbox)
library(data.table)
library(lubridate)
library(readr)
library(sf)
library(ggplot2)
library(secr)
library(ggmice)
library(mice)
library(corrplot)
library(MASS)

rm(list=ls())

#set seed
set.seed(1007)
#source data cleaning script
source("~/00_data_import_cleaning.R")

#get a beta estimate from a linear model
set.seed(1007)
get_betas<-function(model_fit, ndraws){
  betas = coef(model_fit)
  vcov = vcov(model_fit)
  betas_simulated = MASS::mvrnorm(ndraws, betas, vcov)
  return(betas_simulated)
}

##### 
#### FULL DATA IMPUTATION, DICHOTOMOUS OUTCOME
#####

#again change the missing "day" data to 999 and then sub out after imputing
full_dat_imp<-full_data[!is.na(yr2022)]

full_dat_imp[,`:=`(date_min=NULL)]
full_dat_imp[is.na(day),day:=999]

imp<-mice(full_dat_imp, method='pmm', m=10)

#analyze at each imputation
df_uni_f<-data.frame()
betas_imputed_uni_f<-data.frame()
ests_imputed_mm1_f<-data.frame()
ests_imputed_mm2_f<-data.frame()
betas_imputed_mm1_f<-data.frame()
betas_imputed_mm2_f<-data.frame()
for(x in 1:10){
  imputed<-complete(imp, x)
  setDT(imputed)
  imputed[day==999, day:=NA]
  imputed[,log_GDP:=log(yr2022)]
  
  imputed[,norm_fem_edu:=(Female_edu_mean_yrs_25_29-min(Female_edu_mean_yrs_25_29, na.rm=T))/(max(Female_edu_mean_yrs_25_29, na.rm=T)-min(Female_edu_mean_yrs_25_29, na.rm=T))]
  imputed[,logit_fem_edu:=log((norm_fem_edu)/(1-norm_fem_edu))]
  imputed[norm_fem_edu %in% c(0,1),logit_fem_edu:=log((norm_fem_edu+(0.5/184))/(1-norm_fem_edu+(0.5/184)))]
  
  imputed[,norm_gender_phone_gap:=(Risk_comms_gender_gap_access_phone_3_6_3a)/100]
  imputed[,logit_gender_phone_gap:=log((norm_gender_phone_gap)/(1-norm_gender_phone_gap))]
  imputed[norm_gender_phone_gap %in% c(0,1),logit_gender_phone_gap:=log((norm_gender_phone_gap+(0.5/184))/(1-norm_gender_phone_gap+(0.5/184)))]
  
  imputed[,norm_gender_internet_gap:=(Risk_comms_gender_gap_access_internet_3_6_4_a)/100]
  imputed[,logit_gender_internet_gap:=log((norm_gender_internet_gap+(0.5/184))/(1-norm_gender_internet_gap+(0.5/184)))]
  imputed[norm_gender_internet_gap %in% c(0,1),logit_gender_internet_gap:=log((norm_gender_internet_gap+(0.5/184))/(1-norm_gender_internet_gap+(0.5/184)))]
  
  imputed[,logit_lib:=log((lib_vdem_owid)/(1-lib_vdem_owid))]
  imputed[lib_vdem_owid %in% c(0,1),logit_lib:=log((lib_vdem_owid+(0.5/184))/(1-lib_vdem_owid+(0.5/184)))]
  
  
  imputed[,factor_pop_inclusion_riskcomm:=factor(Risk_comms_pop_inclusion_3_5_1b, levels=c(0,100),
                                                labels=c('No',"Yes"))]
  imputed[,factor_misinfo_riskcomm:=factor(Risk_comms_leader_share_misinfo_3_5_2b, levels=c(0,100),
                                          labels=c('No',"Yes"))]

  lm1<-glm(mpox_greater~CPI_2022,data=imputed,family = binomial(link = "logit"))
  lm2<-glm(mpox_greater~GAI,data=imputed,family = binomial(link = "logit"))
  lm3<-glm(mpox_greater~HAQI_mean,data=imputed,family = binomial(link = "logit"))
  lm4<-glm(mpox_greater~log_GDP,data=imputed,family = binomial(link = "logit"))
  lm5<-glm(mpox_greater~logit_fem_edu,data=imputed,family = binomial(link = "logit")) 
  lm6<-glm(mpox_greater~Overall,data=imputed,family = binomial(link = "logit")) 
  lm7<-glm(mpox_greater~Risk_comm_3_5,data=imputed,family = binomial(link = "logit"))
  lm8<-glm(mpox_greater~factor_pop_inclusion_riskcomm,data=imputed,family = binomial(link = "logit"))
  lm9<-glm(mpox_greater~factor_misinfo_riskcomm,data=imputed,family = binomial(link = "logit"))
  lm10<-glm(mpox_greater~logit_gender_internet_gap,data=imputed,family = binomial(link = "logit"))
  lm11<-glm(mpox_greater~Risk_coms_mobile_subscribers_3_6_2,data=imputed,family = binomial(link = "logit"))
  lm12<-glm(mpox_greater~logit_gender_phone_gap,data=imputed,family = binomial(link = "logit"))
  lm13<-glm(mpox_greater~Risk_comms_pct_hh_internet_3_6_1a,data=imputed,family = binomial(link = "logit"))
  lm14<-glm(mpox_greater~day,data=imputed,family = binomial(link = "logit"))
  lm15<-glm(mpox_greater~factor(regime_row_owid),data=imputed,family = binomial(link = "logit"))
  lm16<-glm(mpox_greater~electdem_vdem_owid,data=imputed,family = binomial(link = "logit"))
  lm17<-glm(mpox_greater~logit_lib,data=imputed,family = binomial(link = "logit"))
  lm18<-glm(mpox_greater~factor(ever_mpox),data=imputed,family = binomial(link = "logit"))
  
  #test each var if signif:
  df<-data.frame()
  betas_uni<-data.frame()
  for(mod in 1:18){
    model<-paste0('lm',mod)
    rows<-nrow(summary(get(model))$coef)
    if(rows>2){
      rows<-2:rows
    }
    ests<-data.frame(imp=x,
                     var=rownames(summary(get(model))$coef)[rows],
                     est=summary(get(model))$coef[rows,1],
                     lci=summary(get(model))$coef[rows,1]-1.96*summary(get(model))$coef[-1,2],
                     uci=summary(get(model))$coef[rows,1]+1.96*summary(get(model))$coef[-1,2],
                     p=summary(get(model))$coef[rows,4])
    df<-rbind(df,ests) 
    #beta hat estimate for each model
    beta_ests<-data.frame(imp=x,
                          var=rownames(summary(get(model))$coef)[rows],
                          beta=get_betas(get(model),1)[rows])
    betas_uni<-rbind(betas_uni, beta_ests)
  }
  setDT(df)
  setDT(betas_uni)
  
  df_uni_f<-rbind(df_uni_f,df)
  betas_imputed_uni_f<-rbind(betas_imputed_uni_f, betas_uni)
  #remove highly collinear variables
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day', 
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  #cutoff at 0.1 / 18 
  keep<-df[p<0.00556,var]
  
  formula<- paste0("mpox_greater~",paste(keep, collapse = "+"))
  formula<-gsub("Yes","",formula)
  formula<-gsub("[)]1",")", formula)
  formula<-gsub("[)]3",")", formula)
  
  mm1<-glm(formula,data=imputed,family = binomial(link = "logit"))

  mm1_ests<-data.frame(imp=x,
                       var=rownames(summary(mm1)$coef)[-1],
                       est=summary(mm1)$coef[-1,1],
                       lci=summary(mm1)$coef[-1,1]-1.96*summary(mm1)$coef[-1,2],
                       uci=summary(mm1)$coef[-1,1]+1.96*summary(mm1)$coef[-1,2],
                       p=summary(mm1)$coef[-1,4])
  beta_mm1<-data.frame(imp=x,
                       var=rownames(summary(mm1)$coef)[-1],
                       beta=get_betas(mm1,1)[-1])
  ests_imputed_mm1_f<-rbind(ests_imputed_mm1_f, mm1_ests)
  betas_imputed_mm1_f<-rbind(betas_imputed_mm1_f,beta_mm1)

  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022',  'day', 
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]

}

#print out univariable df with the p-values for inclusion into multivariable mode:
df_uni_f[p<0.1/18,]

#print out all variables significant at 0.05 / 18 in the multivariable model 
setDT(ests_imputed_mm1_f)
ests_imputed_mm1_f[p<0.05/18,]

######
### FULL IMPUTATION, CONTINUOUS OUTCOME
#######

#again change the missing "day" data to 999 and then sub out after imputing
full_dat_imp<-full_data[!is.na(yr2022)]

full_dat_imp[,`:=`(date_min=NULL)]
full_dat_imp[is.na(day),day:=999]

imp<-mice(full_dat_imp, method='pmm', m=10)

#analyze at each imputation
df_uni_f<-data.frame()
betas_imputed_uni_f<-data.frame()
ests_imputed_mm1_f<-data.frame()
ests_imputed_mm2_f<-data.frame()
betas_imputed_mm1_f<-data.frame()
betas_imputed_mm2_f<-data.frame()
for(x in 1:10){
  imputed<-complete(imp, x)
  setDT(imputed)
  imputed[day==999, day:=NA]
  imputed[,log_GDP:=log(yr2022)]
  
  imputed[,norm_fem_edu:=(Female_edu_mean_yrs_25_29-min(Female_edu_mean_yrs_25_29, na.rm=T))/(max(Female_edu_mean_yrs_25_29, na.rm=T)-min(Female_edu_mean_yrs_25_29, na.rm=T))]
  imputed[,logit_fem_edu:=log((norm_fem_edu)/(1-norm_fem_edu))]
  imputed[norm_fem_edu %in% c(0,1),logit_fem_edu:=log((norm_fem_edu+(0.5/184))/(1-norm_fem_edu+(0.5/184)))]
  
  imputed[,norm_gender_phone_gap:=(Risk_comms_gender_gap_access_phone_3_6_3a)/100]
  imputed[,logit_gender_phone_gap:=log((norm_gender_phone_gap)/(1-norm_gender_phone_gap))]
  imputed[norm_gender_phone_gap %in% c(0,1),logit_gender_phone_gap:=log((norm_gender_phone_gap+(0.5/184))/(1-norm_gender_phone_gap+(0.5/184)))]
  
  imputed[,norm_gender_internet_gap:=(Risk_comms_gender_gap_access_internet_3_6_4_a)/100]
  imputed[,logit_gender_internet_gap:=log((norm_gender_internet_gap+(0.5/184))/(1-norm_gender_internet_gap+(0.5/184)))]
  imputed[norm_gender_internet_gap %in% c(0,1),logit_gender_internet_gap:=log((norm_gender_internet_gap+(0.5/184))/(1-norm_gender_internet_gap+(0.5/184)))]
  
  imputed[,logit_lib:=log((lib_vdem_owid)/(1-lib_vdem_owid))]
  imputed[lib_vdem_owid %in% c(0,1),logit_lib:=log((lib_vdem_owid+(0.5/184))/(1-lib_vdem_owid+(0.5/184)))]
  
  
  imputed[,factor_pop_inclusion_riskcomm:=factor(Risk_comms_pop_inclusion_3_5_1b, levels=c(0,100),
                                                 labels=c('No',"Yes"))]
  imputed[,factor_misinfo_riskcomm:=factor(Risk_comms_leader_share_misinfo_3_5_2b, levels=c(0,100),
                                           labels=c('No',"Yes"))]

  lm1<-glm(avg_pct_mpox~CPI_2022,data=imputed)
  lm2<-glm(avg_pct_mpox~GAI,data=imputed)
  lm3<-glm(avg_pct_mpox~HAQI_mean,data=imputed)
  lm4<-glm(avg_pct_mpox~log_GDP,data=imputed)
  lm5<-glm(avg_pct_mpox~logit_fem_edu,data=imputed) 
  lm6<-glm(avg_pct_mpox~Overall,data=imputed)
  lm7<-glm(avg_pct_mpox~Risk_comm_3_5,data=imputed)
  lm8<-glm(avg_pct_mpox~factor_pop_inclusion_riskcomm,data=imputed)
  lm9<-glm(avg_pct_mpox~factor_misinfo_riskcomm,data=imputed)
  lm10<-glm(avg_pct_mpox~logit_gender_internet_gap,data=imputed)
  lm11<-glm(avg_pct_mpox~Risk_coms_mobile_subscribers_3_6_2,data=imputed)
  lm12<-glm(avg_pct_mpox~logit_gender_phone_gap,data=imputed)
  lm13<-glm(avg_pct_mpox~Risk_comms_pct_hh_internet_3_6_1a,data=imputed)
  lm14<-glm(avg_pct_mpox~day,data=imputed)
  lm15<-glm(avg_pct_mpox~factor(regime_row_owid),data=imputed)
  lm16<-glm(avg_pct_mpox~electdem_vdem_owid,data=imputed)
  lm17<-glm(avg_pct_mpox~logit_lib,data=imputed)
  lm18<-glm(avg_pct_mpox~factor(ever_mpox),data=imputed)
  
  #test each var if signif:
  df<-data.frame()
  betas_uni<-data.frame()
  for(mod in 1:18){
    model<-paste0('lm',mod)
    rows<-nrow(summary(get(model))$coef)
    if(rows>2){
      rows<-2:rows
    }
    ests<-data.frame(imp=x,
                     var=rownames(summary(get(model))$coef)[rows],
                     est=summary(get(model))$coef[rows,1],
                     lci=summary(get(model))$coef[rows,1]-1.96*summary(get(model))$coef[-1,2],
                     uci=summary(get(model))$coef[rows,1]+1.96*summary(get(model))$coef[-1,2],
                     p=summary(get(model))$coef[rows,4])
    df<-rbind(df,ests) 
    #beta hat estimate for each model
    beta_ests<-data.frame(imp=x,
                          var=rownames(summary(get(model))$coef)[rows],
                          beta=get_betas(get(model),1)[rows])
    betas_uni<-rbind(betas_uni, beta_ests)
  }
  setDT(df)
  setDT(betas_uni)
  
  df_uni_f<-rbind(df_uni_f,df)
  betas_imputed_uni_f<-rbind(betas_imputed_uni_f, betas_uni)
  #remove highly collinear variables
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day', 
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  #cutoff at 0.1 / 18
  keep<-df[p<0.00556,var]

  formula<-paste0("avg_pct_mpox~",paste(keep, collapse = "+"))
  formula<-gsub("Yes","",formula)
  formula<-gsub("[)]1",")", formula)
  formula<-gsub("[)]2",")", formula)
  formula<-gsub("[)]3",")", formula)
  
  mm1<-glm(formula,data=imputed)
  
  mm1_ests<-data.frame(imp=x,
                       var=rownames(summary(mm1)$coef)[-1],
                       est=summary(mm1)$coef[-1,1],
                       lci=summary(mm1)$coef[-1,1]-1.96*summary(mm1)$coef[-1,2],
                       uci=summary(mm1)$coef[-1,1]+1.96*summary(mm1)$coef[-1,2],
                       p=summary(mm1)$coef[-1,4])
  beta_mm1<-data.frame(imp=x,
                       var=rownames(summary(mm1)$coef)[-1],
                       beta=get_betas(mm1,1)[-1])
  ests_imputed_mm1_f<-rbind(ests_imputed_mm1_f, mm1_ests)
  betas_imputed_mm1_f<-rbind(betas_imputed_mm1_f,beta_mm1)

  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day', 
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
}
#print out univariable df with the p-values for inclusion into multivariable model, removing collinear variables:
df_uni<-df_uni_f[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day', 
                           'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
df_uni[p<0.1/18,]

#print out all variables significant at 0.05 / 18 in the multivariable model 
setDT(ests_imputed_mm1_f)
ests_imputed_mm1_f[p<0.05/18,]

      