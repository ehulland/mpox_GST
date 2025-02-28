library(devtools)
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
source("~/00_data_import_cleaning_foreign_language.R")

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
full_dat_imp<-full_dat_imp[NAME %in% full_dat$NAME,]

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
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022',  'day',
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  
  keep<-df[p<0.1,var]
  
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

  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day',
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  
  keep<-df[p<0.05,var]
  formula<- paste0("mpox_greater~",paste(keep, collapse = "+"))
  formula<-gsub("Yes","",formula)
  formula<-gsub("[)]1",")", formula)
  mm2<-glm(formula,data=imputed,family = binomial(link = "logit")) #signif
  
  mm2_ests<-data.frame(imp=x,
                       var=rownames(summary(mm2)$coef)[-1],
                       est=summary(mm2)$coef[-1,1],
                       lci=summary(mm2)$coef[-1,1]-1.96*summary(mm2)$coef[-1,2],
                       uci=summary(mm2)$coef[-1,1]+1.96*summary(mm2)$coef[-1,2],
                       p=summary(mm2)$coef[-1,4])
  beta_mm2<-data.frame(imp=x,
                       var=rownames(summary(mm2)$coef)[-1],
                       beta=get_betas(mm2,1)[-1])
  ests_imputed_mm2_f<-rbind(ests_imputed_mm2_f, mm2_ests)
  betas_imputed_mm2_f<-rbind(betas_imputed_mm2_f,beta_mm2)
}
#summarize across imputations
setDT(betas_imputed_mm1_f)
betas_imputed_mm1_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm1_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm1_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm1_f_imputation<-betas_imputed_mm1_f[,c('var','lower_est','upper_est','median_est')]
mm1_f_imputation<-mm1_f_imputation[!duplicated(mm1_f_imputation)]


setDT(betas_imputed_mm2_f)
betas_imputed_mm2_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm2_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm2_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm2_f_imputation<-betas_imputed_mm2_f[,c('var','lower_est','upper_est','median_est')]
mm2_f_imputation<-mm2_f_imputation[!duplicated(mm2_f_imputation)]

setDT(betas_imputed_uni_f)
betas_imputed_uni_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_uni_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_uni_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

uni_f_imputation<-betas_imputed_uni_f[,c('var','lower_est','upper_est','median_est')]
uni_f_imputation<-uni_f_imputation[!duplicated(uni_f_imputation)]

#add color coding & clean up variable names
uni_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
uni_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]

#add color coding & clean up variable names
mm1_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
mm1_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]

#add color coding & clean up variable names
mm2_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
mm2_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]

uni_f_imputation[,included:=ifelse(var %in% c('factor(regime_row_owid)1', 'factor(regime_row_owid)2', 'factor(regime_row_owid)3','HAQI_mean','day',
                                            'logit_lib', 'CPI_2022', 'Risk_comms_pct_hh_internet_3_6_1a'),'Excluded - multicollinearity',
                                 ifelse(var %in% mm1_f_imputation$var, 'Included', 'Excluded - significance'))]

uni_f_imputation<-uni_f_imputation[order(included)]

#Convert to ORs

uni_f_imputation2<-uni_f_imputation
uni_f_imputation2[,`:=`(median_est=exp(median_est), lower_est=exp(lower_est), upper_est=exp(upper_est))]

mm1_f_imputation2<-mm1_f_imputation
mm1_f_imputation2[,`:=`(median_est=exp(median_est), lower_est=exp(lower_est), upper_est=exp(upper_est))]

mm2_f_imputation2<-mm2_f_imputation
mm2_f_imputation2[,`:=`(median_est=exp(median_est), lower_est=exp(lower_est), upper_est=exp(upper_est))]


jpeg(paste0(dir,'/figures/univar_df_dichotomous_foreign_grouped_OR.jpeg'), height=700, width=1000)
ggplot(data=uni_f_imputation2)+geom_vline(xintercept = 1, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance', values=c('Black','Red','Green'))+
  ylab('')+xlab('Odds Ratio (95% UI)')+
  # theme(
  #   panel.border = element_blank(), 
  #   panel.grid.major.y = element_line(color = "light grey", size = 0.3),
  #   panel.grid.major.x = element_line(color = "light grey", size = 0.3),
  #   panel.grid.minor = element_blank(), 
  #   axis.line = element_blank())
  facet_grid(included~., scales='free',switch='y')+
  theme_classic() +coord_cartesian(xlim=c(0,15))+
  theme(panel.spacing=unit(2, "lines"),
        , strip.placement.y = "outside"
        , strip.background = element_blank()
        , strip.text = element_text(face = "bold"),
        text=element_text(size=20))+
  geom_segment(data=uni_f_imputation2[upper_est>15,],aes(x = 15.5, xend = 15.8, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(data=uni_f_imputation2[lower_est<0.01,],aes(x = 0.02, xend = 0.01, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))

dev.off()

jpeg(paste0(dir,'/figures/multivar01_df_dichotomous_foreign_grouped_OR.jpeg'), height=700, width=1000)
ggplot(data=mm1_f_imputation2)+geom_vline(xintercept = 1, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance', values=c('Black','Red','Green'))+
  ylab('')+xlab('Odds Ratio (95% UI)')+ 
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+theme_classic() +coord_cartesian(xlim=c(0,15))+
  theme(panel.spacing=unit(2, "lines"),
        , strip.placement.y = "outside"
        , strip.background = element_blank()
        , strip.text = element_text(face = "bold"),
        text=element_text(size=20))+
  geom_segment(data=mm1_f_imputation2[upper_est>15,],aes(x = 15.5, xend = 15.8, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(data=mm1_f_imputation2[lower_est<0.01,],aes(x = 0.02, xend = 0.01, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))

dev.off()

##### 
#### FULL DATA IMPUTATION, CONTINUOUS OUTCOME
#####
#again change the missing "day" data to 999 and then sub out after imputing
full_dat_imp<-full_data[!is.na(yr2022)]
full_dat_imp<-full_dat_imp[NAME %in% full_dat$NAME,]

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
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022',  "day",
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  
  keep<-df[p<0.1,var]
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
 
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022',"day",
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  
  keep<-df[p<0.05,var]
  formula<- paste0("avg_pct_mpox~",paste(keep, collapse = "+"))
  formula<-gsub("Yes","",formula)
  formula<-gsub("[)]1",")", formula)
  mm2<-glm(formula,data=imputed)
  
  mm2_ests<-data.frame(imp=x,
                       var=rownames(summary(mm2)$coef)[-1],
                       est=summary(mm2)$coef[-1,1],
                       lci=summary(mm2)$coef[-1,1]-1.96*summary(mm2)$coef[-1,2],
                       uci=summary(mm2)$coef[-1,1]+1.96*summary(mm2)$coef[-1,2],
                       p=summary(mm2)$coef[-1,4])
  beta_mm2<-data.frame(imp=x,
                       var=rownames(summary(mm2)$coef)[-1],
                       beta=get_betas(mm2,1)[-1])
  ests_imputed_mm2_f<-rbind(ests_imputed_mm2_f, mm2_ests)
  betas_imputed_mm2_f<-rbind(betas_imputed_mm2_f,beta_mm2)
}
#summarize across imputations
setDT(betas_imputed_mm1_f)
betas_imputed_mm1_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm1_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm1_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm1_f_imputation<-betas_imputed_mm1_f[,c('var','lower_est','upper_est','median_est')]
mm1_f_imputation<-mm1_f_imputation[!duplicated(mm1_f_imputation)]


setDT(betas_imputed_mm2_f)
betas_imputed_mm2_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm2_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm2_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm2_f_imputation<-betas_imputed_mm2_f[,c('var','lower_est','upper_est','median_est')]
mm2_f_imputation<-mm2_f_imputation[!duplicated(mm2_f_imputation)]

setDT(betas_imputed_uni_f)
betas_imputed_uni_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_uni_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_uni_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

uni_f_imputation<-betas_imputed_uni_f[,c('var','lower_est','upper_est','median_est')]
uni_f_imputation<-uni_f_imputation[!duplicated(uni_f_imputation)]

#add color coding & clean up variable names
uni_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
uni_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]
#add color coding & clean up variable names
mm1_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
mm1_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]
#add color coding & clean up variable names
mm2_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
mm2_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]

uni_f_imputation[,included:=ifelse(var %in% c('factor(regime_row_owid)1', 'factor(regime_row_owid)2', 'factor(regime_row_owid)3','HAQI_mean','day',
                                              'logit_lib', 'CPI_2022', 'Risk_comms_pct_hh_internet_3_6_1a'),'Excluded - multicollinearity',
                                   ifelse(var %in% mm1_f_imputation$var, 'Included', 'Excluded - significance'))]

uni_f_imputation<-uni_f_imputation[order(included)]


jpeg(paste0(dir,'/figures/univar_df_continuous_foreign_grouped.jpeg'), height=700, width=1000)
ggplot(data=uni_f_imputation)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance', values=c('Black','Red','Green'))+
  ylab('')+xlab('Point estimate (95% UI)')+
  # theme(
  #   panel.border = element_blank(), 
  #   panel.grid.major.y = element_line(color = "light grey", size = 0.3),
  #   panel.grid.major.x = element_line(color = "light grey", size = 0.3),
  #   panel.grid.minor = element_blank(), 
  #   axis.line = element_blank())
  facet_grid(included~., scales='free',switch='y')+
  theme_classic() +coord_cartesian(xlim=c(-0.5,0.5))+
  theme(panel.spacing=unit(2, "lines"),
        , strip.placement.y = "outside"
        , strip.background = element_blank()
        , strip.text = element_text(face = "bold"),
        text=element_text(size=20))
dev.off()

jpeg(paste0(dir,'/figures/multivar01_df_continuous_foreign_grouped.jpeg'), height=700, width=1000)
ggplot(data=mm1_f_imputation)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance', values=c('Black','Red','Green'))+
  ylab('')+xlab('Point estimate (95% UI)')+ 
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+coord_cartesian(xlim=c(-0.5,0.5))+theme_classic() +
  theme(panel.spacing=unit(2, "lines"),
        , strip.placement.y = "outside"
        , strip.background = element_blank()
        , strip.text = element_text(face = "bold"),
        text=element_text(size=20))
dev.off()

#########################################################################
## Trial run imputing YES and NO for dichotomous analysis for Ecuador and Madagascar!
## The current statistics for Ecuador are that mpox_greater = 1 and avg_pct_mpox is 0.83333. So try a version with Ecuador mpox_greater = 0
## The current statistics for Madagascar are that mpox_greater = 0 and avg_pct_mpox is 0. So try a version with Madagascar mpox_greater = 1.
#############################################################################

# Combos to try: ecuador 1 madagascar 1, ecuador 0 madagascar 1, ecuador 0 madagascar 0 (because ecuador 1 and madagascar 0 is original)

temp1<-temp2<-temp3<-copy(full_data)

######
### TRIAL 1: Both Ecuador and Madagascar have mpox_greater value of 0
 ######

##### 
#### FULL DATA IMPUTATION
#####

#again change the missing "day" data to 999 and then sub out after imputing
temp1[NAME=="Ecuador", mpox_greater:=0]
temp1_imp<-temp1[!is.na(yr2022)]
temp1_imp<-temp1_imp[NAME %in% full_dat$NAME]

temp1_imp[,`:=`(date_min=NULL)]
temp1_imp[is.na(day),day:=999]

imp<-mice(temp1_imp, method='pmm', m=10)
set.seed(1007)
get_betas<-function(model_fit, ndraws){
  betas = coef(model_fit)
  vcov = vcov(model_fit)
  betas_simulated = MASS::mvrnorm(ndraws, betas, vcov)
  return(betas_simulated)
}
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
  
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', "day",
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  
  keep<-df[p<0.1,var]
  
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
  
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', "day",
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  
  keep<-df[p<0.05,var]
  formula<- paste0("mpox_greater~",paste(keep, collapse = "+"))
  formula<-gsub("Yes","",formula)
  formula<-gsub("[)]1",")", formula)
  mm2<-glm(formula,data=imputed,family = binomial(link = "logit"))
  
  mm2_ests<-data.frame(imp=x,
                       var=rownames(summary(mm2)$coef)[-1],
                       est=summary(mm2)$coef[-1,1],
                       lci=summary(mm2)$coef[-1,1]-1.96*summary(mm2)$coef[-1,2],
                       uci=summary(mm2)$coef[-1,1]+1.96*summary(mm2)$coef[-1,2],
                       p=summary(mm2)$coef[-1,4])
  beta_mm2<-data.frame(imp=x,
                       var=rownames(summary(mm2)$coef)[-1],
                       beta=get_betas(mm2,1)[-1])
  ests_imputed_mm2_f<-rbind(ests_imputed_mm2_f, mm2_ests)
  betas_imputed_mm2_f<-rbind(betas_imputed_mm2_f,beta_mm2)
}
setDT(betas_imputed_mm1_f)
betas_imputed_mm1_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm1_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm1_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm1_f_imputation<-betas_imputed_mm1_f[,c('var','lower_est','upper_est','median_est')]
mm1_f_imputation<-mm1_f_imputation[!duplicated(mm1_f_imputation)]


setDT(betas_imputed_mm2_f)
betas_imputed_mm2_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm2_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm2_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm2_f_imputation<-betas_imputed_mm2_f[,c('var','lower_est','upper_est','median_est')]
mm2_f_imputation<-mm2_f_imputation[!duplicated(mm2_f_imputation)]

setDT(betas_imputed_uni_f)
betas_imputed_uni_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_uni_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_uni_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

uni_f_imputation<-betas_imputed_uni_f[,c('var','lower_est','upper_est','median_est')]
uni_f_imputation<-uni_f_imputation[!duplicated(uni_f_imputation)]

#add color coding & clean up variable names
uni_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
uni_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]
#add color coding & clean up variable names
mm1_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
mm1_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]
#add color coding & clean up variable names
mm2_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
mm2_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]

uni_f_imputation[,included:=ifelse(var %in% c('factor(regime_row_owid)1', 'factor(regime_row_owid)2', 'factor(regime_row_owid)3','HAQI_mean','day',
                                              'logit_lib', 'CPI_2022', 'Risk_comms_pct_hh_internet_3_6_1a'),'Excluded - multicollinearity',
                                   ifelse(var %in% mm1_f_imputation$var, 'Included', 'Excluded - significance'))]

uni_f_imputation<-uni_f_imputation[order(included)]


#Convert to ORs

uni_f_imputation2<-uni_f_imputation
uni_f_imputation2[,`:=`(median_est=exp(median_est), lower_est=exp(lower_est), upper_est=exp(upper_est))]

mm1_f_imputation2<-mm1_f_imputation
mm1_f_imputation2[,`:=`(median_est=exp(median_est), lower_est=exp(lower_est), upper_est=exp(upper_est))]

mm2_f_imputation2<-mm2_f_imputation
mm2_f_imputation2[,`:=`(median_est=exp(median_est), lower_est=exp(lower_est), upper_est=exp(upper_est))]



jpeg(paste0(dir,'/figures/univar_df_dichotomous_foreign_both0_grouped_OR.jpeg'), height=700, width=1000)
ggplot(data=uni_f_imputation2)+geom_vline(xintercept = 1, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance', values=c('Black','Green'))+
  ylab('')+xlab('Odds Ratio (95% UI)')+ 
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+coord_cartesian(xlim=c(0,15))+theme_classic() +
  theme(panel.spacing=unit(2, "lines"),
        , strip.placement.y = "outside"
        , strip.background = element_blank()
        , strip.text = element_text(face = "bold"),
        text=element_text(size=20))+
  geom_segment(data=uni_f_imputation2[upper_est>15,],aes(x = 15.5, xend = 15.8, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(data=uni_f_imputation2[lower_est<0.01,],aes(x = 0.02, xend = 0.01, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))

dev.off()

jpeg(paste0(dir,'/figures/multivar01_df_dichotomous_foreign_both0_grouped.jpeg'), height=700, width=1000)
ggplot(data=mm1_f_imputation2)+geom_vline(xintercept = 1, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance', values=c('Black','Red','Green'))+
  ylab('')+xlab('Odds Ratio (95% UI)')+ 
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+coord_cartesian(xlim=c(0,15))+theme_classic() +
  theme(panel.spacing=unit(2, "lines"),
        , strip.placement.y = "outside"
        , strip.background = element_blank()
        , strip.text = element_text(face = "bold"),
        text=element_text(size=20))+
  geom_segment(data=mm1_f_imputation2[upper_est>15,],aes(x = 15.5, xend = 15.8, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(data=mm1_f_imputation2[lower_est<0.01,],aes(x = 0.02, xend = 0.01, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))

dev.off()

######
### TRIAL 2: Both Ecuador and Madagascar have mpox_greater value of 1
######


temp2[NAME=="Madagascar", mpox_greater:=1]
temp2_imp<-temp2[!is.na(yr2022)]
temp2_imp<-temp2_imp[NAME %in% full_dat$NAME,]

temp2_imp[,`:=`(date_min=NULL)]
temp2_imp[is.na(day),day:=999]

imp<-mice(temp2_imp, method='pmm', m=10)
set.seed(1007)
get_betas<-function(model_fit, ndraws){
  betas = coef(model_fit)
  vcov = vcov(model_fit)
  betas_simulated = MASS::mvrnorm(ndraws, betas, vcov)
  return(betas_simulated)
}
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
  
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', "day",
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  
  keep<-df[p<0.1,var]
  
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

  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', "day",
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  
  keep<-df[p<0.05,var]
  formula<- paste0("mpox_greater~",paste(keep, collapse = "+"))
  formula<-gsub("Yes","",formula)
  formula<-gsub("[)]1",")", formula)
  mm2<-glm(formula,data=imputed,family = binomial(link = "logit"))
  
  mm2_ests<-data.frame(imp=x,
                       var=rownames(summary(mm2)$coef)[-1],
                       est=summary(mm2)$coef[-1,1],
                       lci=summary(mm2)$coef[-1,1]-1.96*summary(mm2)$coef[-1,2],
                       uci=summary(mm2)$coef[-1,1]+1.96*summary(mm2)$coef[-1,2],
                       p=summary(mm2)$coef[-1,4])
  beta_mm2<-data.frame(imp=x,
                       var=rownames(summary(mm2)$coef)[-1],
                       beta=get_betas(mm2,1)[-1])
  ests_imputed_mm2_f<-rbind(ests_imputed_mm2_f, mm2_ests)
  betas_imputed_mm2_f<-rbind(betas_imputed_mm2_f,beta_mm2)
}
setDT(betas_imputed_mm1_f)
betas_imputed_mm1_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm1_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm1_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm1_f_imputation<-betas_imputed_mm1_f[,c('var','lower_est','upper_est','median_est')]
mm1_f_imputation<-mm1_f_imputation[!duplicated(mm1_f_imputation)]


setDT(betas_imputed_mm2_f)
betas_imputed_mm2_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm2_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm2_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm2_f_imputation<-betas_imputed_mm2_f[,c('var','lower_est','upper_est','median_est')]
mm2_f_imputation<-mm2_f_imputation[!duplicated(mm2_f_imputation)]

setDT(betas_imputed_uni_f)
betas_imputed_uni_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_uni_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_uni_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

uni_f_imputation<-betas_imputed_uni_f[,c('var','lower_est','upper_est','median_est')]
uni_f_imputation<-uni_f_imputation[!duplicated(uni_f_imputation)]

#add color coding & clean up variable names
uni_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
uni_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]
#add color coding & clean up variable names
mm1_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
mm1_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]
#add color coding & clean up variable names
mm2_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
mm2_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]

uni_f_imputation[,included:=ifelse(var %in% c('factor(regime_row_owid)1', 'factor(regime_row_owid)2', 'factor(regime_row_owid)3','HAQI_mean','day',
                                              'logit_lib', 'CPI_2022', 'Risk_comms_pct_hh_internet_3_6_1a'),'Excluded - multicollinearity',
                                   ifelse(var %in% mm1_f_imputation$var, 'Included', 'Excluded - significance'))]

uni_f_imputation<-uni_f_imputation[order(included)]

#Convert to ORs

uni_f_imputation2<-uni_f_imputation
uni_f_imputation2[,`:=`(median_est=exp(median_est), lower_est=exp(lower_est), upper_est=exp(upper_est))]

mm1_f_imputation2<-mm1_f_imputation
mm1_f_imputation2[,`:=`(median_est=exp(median_est), lower_est=exp(lower_est), upper_est=exp(upper_est))]

mm2_f_imputation2<-mm2_f_imputation
mm2_f_imputation2[,`:=`(median_est=exp(median_est), lower_est=exp(lower_est), upper_est=exp(upper_est))]



jpeg(paste0(dir,'/figures/univar_df_dichotomous_foreign_both1_grouped.jpeg'), height=700, width=1000)
ggplot(data=uni_f_imputation2)+geom_vline(xintercept = 1, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance', values=c('Black','Red','Green'))+
  ylab('')+xlab('Odds Ratio (95% UI)')+ 
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+coord_cartesian(xlim=c(0,15))+theme_classic() +
  theme(panel.spacing=unit(2, "lines"),
        , strip.placement.y = "outside"
        , strip.background = element_blank()
        , strip.text = element_text(face = "bold"),
        text=element_text(size=20))+
  geom_segment(data=uni_f_imputation2[upper_est>15,],aes(x = 15.5, xend = 15.8, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(data=uni_f_imputation2[lower_est<0.01,],aes(x = 0.02, xend = 0.01, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))

dev.off()

jpeg(paste0(dir,'/figures/multivar01_df_dichotomous_foreign_both1_grouped.jpeg'), height=700, width=1000)
ggplot(data=mm1_f_imputation2)+geom_vline(xintercept = 1, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance', values=c('Black','Red','Green'))+
  ylab('')+xlab('Odds Ratio (95% UI)')+ 
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+coord_cartesian(xlim=c(0,15))+theme_classic() +
  theme(panel.spacing=unit(2, "lines"),
        , strip.placement.y = "outside"
        , strip.background = element_blank()
        , strip.text = element_text(face = "bold"),
        text=element_text(size=20))+
  geom_segment(data=mm1_f_imputation2[upper_est>15,],aes(x = 15.5, xend = 15.8, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(data=mm1_f_imputation2[lower_est<0.01,],aes(x = 0.02, xend = 0.01, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))

dev.off()


##### 
## TRIAL 3: Ecuador is 0 for mpox_greater, Madagascar is 1
#####

temp3[NAME=="Madagascar", mpox_greater:=1]
temp3[NAME=="Ecuador", mpox_greater:=0]
temp3_imp<-temp3[!is.na(yr2022)]
temp3_imp<-temp3_imp[NAME %in% full_dat$NAME,]

temp3_imp[,`:=`(date_min=NULL)]
temp3_imp[is.na(day),day:=999]

imp<-mice(temp3_imp, method='pmm', m=10)

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
  
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', "day",
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  
  keep<-df[p<0.1,var]
  
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

  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', "day",
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  
  keep<-df[p<0.05,var]
  formula<- paste0("mpox_greater~",paste(keep, collapse = "+"))
  formula<-gsub("Yes","",formula)
  formula<-gsub("[)]1",")", formula)
  mm2<-glm(formula,data=imputed,family = binomial(link = "logit"))
  
  mm2_ests<-data.frame(imp=x,
                       var=rownames(summary(mm2)$coef)[-1],
                       est=summary(mm2)$coef[-1,1],
                       lci=summary(mm2)$coef[-1,1]-1.96*summary(mm2)$coef[-1,2],
                       uci=summary(mm2)$coef[-1,1]+1.96*summary(mm2)$coef[-1,2],
                       p=summary(mm2)$coef[-1,4])
  beta_mm2<-data.frame(imp=x,
                       var=rownames(summary(mm2)$coef)[-1],
                       beta=get_betas(mm2,1)[-1])
  ests_imputed_mm2_f<-rbind(ests_imputed_mm2_f, mm2_ests)
  betas_imputed_mm2_f<-rbind(betas_imputed_mm2_f,beta_mm2)
}
setDT(betas_imputed_mm1_f)
betas_imputed_mm1_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm1_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm1_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm1_f_imputation<-betas_imputed_mm1_f[,c('var','lower_est','upper_est','median_est')]
mm1_f_imputation<-mm1_f_imputation[!duplicated(mm1_f_imputation)]


setDT(betas_imputed_mm2_f)
betas_imputed_mm2_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm2_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm2_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm2_f_imputation<-betas_imputed_mm2_f[,c('var','lower_est','upper_est','median_est')]
mm2_f_imputation<-mm2_f_imputation[!duplicated(mm2_f_imputation)]

setDT(betas_imputed_uni_f)
betas_imputed_uni_f[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_uni_f[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_uni_f[,median_est:=quantile(beta, c(0.5)),by=c('var')]

uni_f_imputation<-betas_imputed_uni_f[,c('var','lower_est','upper_est','median_est')]
uni_f_imputation<-uni_f_imputation[!duplicated(uni_f_imputation)]

#add color coding & clean up variable names
uni_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
uni_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]
#add color coding & clean up variable names
mm1_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
mm1_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]
#add color coding & clean up variable names
mm2_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]
mm2_f_imputation[,Variable:=factor(var,
                                   levels=c(
                                     'logit_fem_edu',
                                     'factor_misinfo_riskcommYes',
                                     'factor_pop_inclusion_riskcommYes','factor(regime_row_owid)3',
                                     'factor(regime_row_owid)1','factor(regime_row_owid)2','Risk_comms_pct_hh_internet_3_6_1a',
                                     "Risk_coms_mobile_subscribers_3_6_2",'logit_lib',
                                     'GAI','HAQI_mean','Risk_comm_3_5','Overall',
                                     'log_GDP','logit_gender_phone_gap','logit_gender_internet_gap',
                                     'factor(ever_mpox)1','electdem_vdem_owid','day','CPI_2022'),
                                   
                                   labels=c('Years of education females 25-29 years (logit)',
                                            'Senior leaders used misinformation (Yes vs. No)',
                                            'Risk communications are inclusive (Yes vs. No)', 
                                            'Political regime (Liberal Democracy versus Closed Autocracy)',
                                            'Political regime (Electoral versus Closed Autocracy)','Political regime (Electoral Democracy versus Closed Autocracy)',
                                            'Percent households with internet',
                                            "Mobile subscribers per 100 population",'Liberal democracy score (logit)',
                                            'LGBT Global Acceptance Index', 'Healthcare Access and Quality Index',
                                            'GHSI 2021 risk communication score','GHSI 2021 overall score',
                                            'GDP per capita (log)','Female access to mobile phone (logit)','Female access to internet (logit)', 
                                            'Ever had a case of mpox yes vs no','Electoral democracy', 
                                            'Days since first 2022 outbreak case reported','Corruption Perceptions Index 2022'))]

uni_f_imputation[,included:=ifelse(var %in% c('factor(regime_row_owid)1', 'factor(regime_row_owid)2', 'factor(regime_row_owid)3','HAQI_mean','day',
                                              'logit_lib', 'CPI_2022', 'Risk_comms_pct_hh_internet_3_6_1a'),'Excluded - multicollinearity',
                                   ifelse(var %in% mm1_f_imputation$var, 'Included', 'Excluded - significance'))]

uni_f_imputation<-uni_f_imputation[order(included)]

#Convert to ORs

uni_f_imputation2<-uni_f_imputation
uni_f_imputation2[,`:=`(median_est=exp(median_est), lower_est=exp(lower_est), upper_est=exp(upper_est))]

mm1_f_imputation2<-mm1_f_imputation
mm1_f_imputation2[,`:=`(median_est=exp(median_est), lower_est=exp(lower_est), upper_est=exp(upper_est))]

mm2_f_imputation2<-mm2_f_imputation
mm2_f_imputation2[,`:=`(median_est=exp(median_est), lower_est=exp(lower_est), upper_est=exp(upper_est))]



jpeg(paste0(dir,'/figures/univar_df_dichotomous_foreign_e0m1_grouped.jpeg'), height=700, width=1000)
ggplot(data=uni_f_imputation2)+geom_vline(xintercept = 1, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance', values=c('Black','Red','Green'))+
  ylab('')+xlab('Odds Ratio (95% UI)')+ 
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+coord_cartesian(xlim=c(0,15))+theme_classic() +
  theme(panel.spacing=unit(2, "lines"),
        , strip.placement.y = "outside"
        , strip.background = element_blank()
        , strip.text = element_text(face = "bold"),
        text=element_text(size=20))+
  geom_segment(data=uni_f_imputation2[upper_est>15,],aes(x = 15.5, xend = 15.8, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(data=uni_f_imputation2[lower_est<0.01,],aes(x = 0.02, xend = 0.01, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))

dev.off()

jpeg(paste0(dir,'/figures/multivar01_df_dichotomous_foreign_e0m1_grouped.jpeg'), height=700, width=1000)
ggplot(data=mm1_f_imputation2)+geom_vline(xintercept = 1, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance', values=c('Black','Red','Green'))+
  ylab('')+xlab('Odds Ratio (95% UI)')+ 
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+coord_cartesian(xlim=c(0,15))+theme_classic() +
  theme(panel.spacing=unit(2, "lines"),
        , strip.placement.y = "outside"
        , strip.background = element_blank()
        , strip.text = element_text(face = "bold"),
        text=element_text(size=20))+
  geom_segment(data=mm1_f_imputation2[upper_est>15,],aes(x = 15.5, xend = 15.8, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(data=mm1_f_imputation2[lower_est<0.01,],aes(x = 0.02, xend = 0.01, y = Variable, col=signif), arrow = arrow(length = unit(0.3, "cm")))

dev.off()
