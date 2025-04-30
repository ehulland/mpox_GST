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
library(pscl)
library(EnvStats)
library(modEvA)
library(glmnet)
rm(list=ls())

#set seed
set.seed(1007)
#source data cleaning script
source("~/mpox_codes/00_data_import_cleaning.R")

##### 
#### FULL DATA IMPUTATION, Dichotomous model, 188 observations
#####

#full imputation
#again change the missing "day" data to 999 and then sub out after imputing
full_dat_imp<-full_data[!is.na(yr2022)]
full_dat_imp[,`:=`(date_min=NULL)]
full_dat_imp[is.na(day),day:=999]

imp<-mice(full_dat_imp, method='pmm', m=10)
set.seed(1007)
get_betas<-function(model_fit, ndraws){
  betas = coef(model_fit)
  vcov = vcov(model_fit)
  betas_simulated = rockchalk::mvrnorm(ndraws, betas, vcov)
  return(betas_simulated)
}
#analyze at each imputation
df_uni_f<-data.frame()
betas_imputed_uni_f<-data.frame()
ests_imputed_mm1_f<-data.frame()
ests_imputed_mm2_f<-data.frame()
betas_imputed_mm1_f<-data.frame()
betas_imputed_mm2_f<-data.frame()

temps<-data.frame()
for(x in 1:10){
  imputed<-complete(imp, x)
  setDT(imputed)
  imputed[day==999, day:=NA]
  imputed[,log_GDP:=log(yr2022)]
  
  imputed[,norm_fem_edu:=(Female_edu_mean_yrs_25_29-min(Female_edu_mean_yrs_25_29, na.rm=T))/(max(Female_edu_mean_yrs_25_29, na.rm=T)-min(Female_edu_mean_yrs_25_29, na.rm=T))]
  imputed[,logit_fem_edu:=log((norm_fem_edu)/(1-norm_fem_edu))]
  imputed[norm_fem_edu %in% c(0,1),logit_fem_edu:=log((norm_fem_edu+(0.5/188))/(1-norm_fem_edu+(0.5/188)))]
  
  imputed[,norm_gender_phone_gap:=(Risk_comms_gender_gap_access_phone_3_6_3a)/100]
  imputed[,logit_gender_phone_gap:=log((norm_gender_phone_gap)/(1-norm_gender_phone_gap))]
  imputed[norm_gender_phone_gap %in% c(0,1),logit_gender_phone_gap:=log((norm_gender_phone_gap+(0.5/188))/(1-norm_gender_phone_gap+(0.5/188)))]
  
  imputed[,norm_gender_internet_gap:=(Risk_comms_gender_gap_access_internet_3_6_4_a)/100]
  imputed[,logit_gender_internet_gap:=log((norm_gender_internet_gap+(0.5/188))/(1-norm_gender_internet_gap+(0.5/188)))]
  imputed[norm_gender_internet_gap %in% c(0,1),logit_gender_internet_gap:=log((norm_gender_internet_gap+(0.5/188))/(1-norm_gender_internet_gap+(0.5/188)))]
  
  imputed[,logit_lib:=log((lib_vdem_owid)/(1-lib_vdem_owid))]
  imputed[lib_vdem_owid %in% c(0,1),logit_lib:=log((lib_vdem_owid+(0.5/188))/(1-lib_vdem_owid+(0.5/188)))]
  
  
  imputed[,factor_pop_inclusion_riskcomm:=factor(Risk_comms_pop_inclusion_3_5_1b, levels=c(0,100),
                                                 labels=c(0,1))]
  imputed[,factor_misinfo_riskcomm:=factor(Risk_comms_leader_share_misinfo_3_5_2b, levels=c(0,100),
                                           labels=c(0,1))]
  
  y<-as.matrix(imputed$mpox_greater)
  X<-as.matrix(imputed[,.(GAI,log_GDP,logit_fem_edu,Overall,Risk_comm_3_5, Risk_comms_pop_inclusion_3_5_1b, Risk_comms_leader_share_misinfo_3_5_2b,
                          logit_gender_internet_gap, Risk_coms_mobile_subscribers_3_6_2,logit_gender_phone_gap,electdem_vdem_owid,ever_mpox)])
  #perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(X, y, alpha = 1)

  
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  best_lambda
  
  mod<-glmnet(X,y,alpha=1,lambda=best_lambda, family = binomial(link = "logit"))
  
  temp<-as.data.frame(as.matrix(coef(mod)))
  temp$imp<-x
  temp$var<-rownames(temp)
  setDT(temp)
  temp<-temp[s0!=0,]
  temp<-temp[var!='(Intercept)',]

  temps<-rbind(temps,temp)
  
}
  
setDT(temps)
temps[,median:=median(s0), by=var]
temps[,q025:=quantile(s0, 0.025), by=var]
temps[,q975:=quantile(s0, 0.975), by=var]
temps2<-temps
temps2[,imp:=NULL]
temps2[,s0:=NULL]

temps_f<-temps2[!duplicated(temps2)]
temps_f[,median_or:=exp(median)]
temps_f[,lor:=exp(q025)]
temps_f[,uor:=exp(q975)]
temps_f[,signif:="Non-significant"]
temps_f[lor>1&uor>1, signif:='Significant - positive']
temps_f[lor<1&uor<1, signif:='Significant - negative']
 
temps_f[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'Risk_comms_leader_share_misinfo_3_5_2b',
                                     'Risk_comms_pop_inclusion_3_5_1b','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'ever_mpox','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBTQ+ Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]
#number of models per variable:

tab<-table(temps$var)
tabs<-as.data.frame(tab)

jpeg(paste0(dir,'/figures/multivar_dichotomous_fullimputation_grouped_lasso2.jpeg'), height=700, width=1000)
ggplot(data=temps_f)+geom_vline(xintercept = 1, col='black',lty=2)+
  geom_point(aes(x=median_or,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lor, xmax=uor, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance', values=c('purple','green'))+
  ylab('')+xlab('Odds Ratio (95% UI)')+
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+coord_cartesian(xlim=c(0,5))+
  theme_classic() +
  theme( panel.spacing=unit(2, "lines")
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20)
  )+
  geom_segment(data=temps_f[uor>5,],aes(x = 5.5, xend = 5.8, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))

dev.off()


tab<-table(temps$var)
tabs<-as.data.frame(tab)
tabs

##### 
#### FULL DATA IMPUTATION, 188 observations
#####

#again change the missing "day" data to 999 and then sub out after imputing
full_dat_imp<-full_data[!is.na(yr2022)]

full_dat_imp[,`:=`(date_min=NULL)]
full_dat_imp[is.na(day),day:=999]

imp<-mice(full_dat_imp, method='pmm', m=10)
set.seed(1007)
get_betas<-function(model_fit, ndraws){
  betas = coef(model_fit)
  vcov = vcov(model_fit)
  betas_simulated = rockchalk::mvrnorm(ndraws, betas, vcov)
  return(betas_simulated)
}
#analyze at each imputation
temps<-data.frame()
for(x in 1:10){
  imputed<-complete(imp, x)
  setDT(imputed)
  imputed[day==999, day:=NA]
  imputed[,log_GDP:=log(yr2022)]
  
  imputed[,norm_fem_edu:=(Female_edu_mean_yrs_25_29-min(Female_edu_mean_yrs_25_29, na.rm=T))/(max(Female_edu_mean_yrs_25_29, na.rm=T)-min(Female_edu_mean_yrs_25_29, na.rm=T))]
  imputed[,logit_fem_edu:=log((norm_fem_edu)/(1-norm_fem_edu))]
  imputed[norm_fem_edu %in% c(0,1),logit_fem_edu:=log((norm_fem_edu+(0.5/188))/(1-norm_fem_edu+(0.5/188)))]
  
  imputed[,norm_gender_phone_gap:=(Risk_comms_gender_gap_access_phone_3_6_3a)/100]
  imputed[,logit_gender_phone_gap:=log((norm_gender_phone_gap)/(1-norm_gender_phone_gap))]
  imputed[norm_gender_phone_gap %in% c(0,1),logit_gender_phone_gap:=log((norm_gender_phone_gap+(0.5/188))/(1-norm_gender_phone_gap+(0.5/188)))]
  
  imputed[,norm_gender_internet_gap:=(Risk_comms_gender_gap_access_internet_3_6_4_a)/100]
  imputed[,logit_gender_internet_gap:=log((norm_gender_internet_gap+(0.5/188))/(1-norm_gender_internet_gap+(0.5/188)))]
  imputed[norm_gender_internet_gap %in% c(0,1),logit_gender_internet_gap:=log((norm_gender_internet_gap+(0.5/188))/(1-norm_gender_internet_gap+(0.5/188)))]
  
  imputed[,logit_lib:=log((lib_vdem_owid)/(1-lib_vdem_owid))]
  imputed[lib_vdem_owid %in% c(0,1),logit_lib:=log((lib_vdem_owid+(0.5/188))/(1-lib_vdem_owid+(0.5/188)))]
  
  
  
  imputed[,factor_pop_inclusion_riskcomm:=factor(Risk_comms_pop_inclusion_3_5_1b, levels=c(0,100),
                                                 labels=c('No',"Yes"))]
  imputed[,factor_misinfo_riskcomm:=factor(Risk_comms_leader_share_misinfo_3_5_2b, levels=c(0,100),
                                           labels=c('No',"Yes"))]
  y<-as.matrix(imputed$avg_pct_mpox)
  X<-as.matrix(imputed[,.(GAI,log_GDP,logit_fem_edu,Overall,Risk_comm_3_5, Risk_comms_pop_inclusion_3_5_1b, Risk_comms_leader_share_misinfo_3_5_2b,
                          logit_gender_internet_gap, Risk_coms_mobile_subscribers_3_6_2,logit_gender_phone_gap,electdem_vdem_owid,ever_mpox)])
  #perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(X, y, alpha = 1)
  
  
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  best_lambda
  
  mod<-glmnet(X,y,alpha=1,lambda=best_lambda)
  temp<-as.data.frame(as.matrix(coef(mod)))
  temp$imp<-x
  temp$var<-rownames(temp)
  setDT(temp)
  temp<-temp[s0!=0,]
  temp<-temp[var!='(Intercept)',]
  temps<-rbind(temps,temp)
  
}

setDT(temps)
temps[,median:=median(s0), by=var]
temps[,q025:=quantile(s0, 0.025), by=var]
temps[,q975:=quantile(s0, 0.975), by=var]
temps2<-temps
temps2[,imp:=NULL]
temps2[,s0:=NULL]

temps_f<-temps2[!duplicated(temps2)]
temps_f[,signif:="Non-significant"]
temps_f[q025>0&q975>0, signif:='Significant - positive']
temps_f[q025<0&q975<0, signif:='Significant - negative']

temps_f[,Variable:=factor(var,
                          levels=c(
                            'logit_fem_edu',
                            'Risk_comms_leader_share_misinfo_3_5_2b',
                            'Risk_comms_pop_inclusion_3_5_1b','factor(regime_row_owid)3',
                            'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                            "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                            'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                            'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                            'ever_mpox','electdem_vdem_owid','day','CPI_2022'),
                          
                          labels=c('Years of education females 25-29 years (logit)',
                                   'Senior leaders used misinformation (Yes vs. No)',
                                   'Risk communications are inclusive (Yes vs. No)', 
                                   'Political regime (Liberal Democracy versus Closed Autocracy)',
                                   'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                   'Percent households with internet',
                                   "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                   'LGBTQ+ Global Acceptance Index', 'Healthcare Access and Quality Index',
                                   'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                   'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                   'Ever had a case of mpox yes vs no','Electoral democracy', 
                                   'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]

jpeg(paste0(dir,'/figures/multivar_continuous_fullimputation_grouped_lasso2.jpeg'), height=700, width=1000)
ggplot(data=temps_f)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=median,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=q025, xmax=q975, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance', values=c('purple','green'))+
  ylab('')+xlab('Beta coefficient (95% UI)')+
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+coord_cartesian(xlim=c(-0.5,0.5))+
  theme_classic() +
  theme( panel.spacing=unit(2, "lines")
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20)
  )
dev.off()

tab<-table(temps$var)
tabs<-as.data.frame(tab)
tabs

save.image(file=paste0(dir,"/univar_fullimpute_LASSO.Rdata"))

