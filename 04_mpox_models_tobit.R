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
library(VGAM)

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
#### FULL DATA IMPUTATION, CONTINUOUS OUTCOME
#####

#again change the missing "day" data to 999 and then sub out after imputing
full_dat_imp<-full_data[!is.na(yr2022)]

full_dat_imp[,`:=`( date_min=NULL)]
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
  
  lm1<-vglm(avg_pct_mpox~CPI_2022,data=imputed, family=tobit(Lower=0))
  lm2<-vglm(avg_pct_mpox~GAI,data=imputed, family=tobit(Lower=0))
  lm3<-vglm(avg_pct_mpox~HAQI_mean,data=imputed, family=tobit(Lower=0))
  lm4<-vglm(avg_pct_mpox~log_GDP,data=imputed, family=tobit(Lower=0))
  lm5<-vglm(avg_pct_mpox~logit_fem_edu,data=imputed, family=tobit(Lower=0)) 
  lm6<-vglm(avg_pct_mpox~Overall,data=imputed, family=tobit(Lower=0))
  lm7<-vglm(avg_pct_mpox~Risk_comm_3_5,data=imputed, family=tobit(Lower=0))
  lm8<-vglm(avg_pct_mpox~factor_pop_inclusion_riskcomm,data=imputed, family=tobit(Lower=0))
  lm9<-vglm(avg_pct_mpox~factor_misinfo_riskcomm,data=imputed, family=tobit(Lower=0))
  lm10<-vglm(avg_pct_mpox~logit_gender_internet_gap,data=imputed, family=tobit(Lower=0))
  lm11<-vglm(avg_pct_mpox~Risk_coms_mobile_subscribers_3_6_2,data=imputed, family=tobit(Lower=0))
  lm12<-vglm(avg_pct_mpox~logit_gender_phone_gap,data=imputed, family=tobit(Lower=0))
  lm13<-vglm(avg_pct_mpox~Risk_comms_pct_hh_internet_3_6_1a,data=imputed, family=tobit(Lower=0))
  lm14<-vglm(avg_pct_mpox~day,data=imputed, family=tobit(Lower=0))
  lm15<-vglm(avg_pct_mpox~regime_row_owid,data=imputed, family=tobit(Lower=0))
  lm16<-vglm(avg_pct_mpox~electdem_vdem_owid,data=imputed, family=tobit(Lower=0))
  lm17<-vglm(avg_pct_mpox~logit_lib,data=imputed, family=tobit(Lower=0))
  lm18<-vglm(avg_pct_mpox~factor(ever_mpox),data=imputed, family=tobit(Lower=0))
  
  
  #test each var if signif:
  df<-data.frame()
  betas_uni_f<-data.frame()
  for(mod in 1:18){
    model<-paste0('lm',mod)
    rows<-nrow(coef(summary(get(model))))
    if(rows>3){
      rows<-3:rows
    }
    ests<-data.frame(var=rownames(coef(summary(get(model))))[rows],
                     est=coef(summary(get(model)))[rows,1],
                     lci=coef(summary(get(model)))[rows,1]-1.96*coef(summary(get(model)))[-c(1:2),2],
                     uci=coef(summary(get(model)))[rows,1]+1.96*coef(summary(get(model)))[-c(1,2),2],
                     p=coef(summary(get(model)))[rows,4])
    df<-rbind(df,ests) 
    
    #beta hat estimate for each model
    beta_ests<-data.frame(imp=x,
                          var=rownames(coef(summary(get(model))))[rows],
                          beta=get_betas(get(model),1)[rows])
    beta_ests<-beta_ests[beta_ests$var!='(Intercept):2',]
    betas_uni_f<-rbind(betas_uni_f, beta_ests)
  }
  setDT(df)
  setDT(betas_uni_f)
  
  df_uni_f<-rbind(df_uni_f,df)
  betas_imputed_uni_f<-rbind(betas_imputed_uni_f, betas_uni_f)
  
  # remove due to high collinearity
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day',
                     'logit_lib','regime_row_owid','day')] 
  keep<-df[p<0.1,var]
  formula<- paste0("avg_pct_mpox~",paste(keep, collapse = "+"))
  formula<-gsub('Yes','',formula)
  formula<-gsub(')1',')',formula)
  
  mm1<-vglm(formula,data=imputed, family=tobit(Lower=0))
  mm1_ests_f<-data.frame(imp=x,
                       var=rownames(coef(summary(mm1)))[-c(1,2)],
                       est=coef(summary(mm1))[-c(1,2),1],
                       lci=coef(summary(mm1))[-c(1,2),1]-1.96*coef(summary(mm1))[-c(1,2),2],
                       uci=coef(summary(mm1))[-c(1,2),1]+1.96*coef(summary(mm1))[-c(1,2),2],
                       p=coef(summary(mm1))[-c(1,2),4])
  
  beta_mm1_f<-data.frame(imp=x,
                       var=rownames(coef(summary(mm1)))[-c(1,2)],
                       beta=get_betas(mm1,1)[-c(1,2)])
  ests_imputed_mm1_f<-rbind(ests_imputed_mm1_f, mm1_ests_f)
  betas_imputed_mm1_f<-rbind(betas_imputed_mm1_f,beta_mm1_f)
  #include only p<0.05
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022',  'day',
                     'logit_lib','regime_row_owid')] 
  keep<-df[p<0.05,var]
  formula<- paste0("avg_pct_mpox~",paste(keep, collapse = "+"))
  formula<-gsub('Yes','',formula)
  formula<-gsub(')1',')',formula)
  mm2<-vglm(formula,data=imputed, family=tobit(Lower=0))
  
  mm2_ests_f<-data.frame(imp=x,
                       var=rownames(coef(summary(mm2)))[-c(1,2)],
                       est=coef(summary(mm2))[-c(1,2),1],
                       lci=coef(summary(mm2))[-c(1,2),1]-1.96*coef(summary(mm2))[-c(1,2),2],
                       uci=coef(summary(mm2))[-c(1,2),1]+1.96*coef(summary(mm2))[-c(1,2),2],
                       p=coef(summary(mm2))[-c(1,2),4])
  
  beta_mm2_f<-data.frame(imp=x,
                       var=rownames(coef(summary(mm2)))[-c(1,2)],
                       beta=get_betas(mm2,1)[-c(1,2)])
  ests_imputed_mm2_f<-rbind(ests_imputed_mm2_f, mm2_ests_f)
  betas_imputed_mm2_f<-rbind(betas_imputed_mm2_f,beta_mm2_f)
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
#add color 
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

#make sure regime row isnt included as it is problematic
uni_f_imputation<-uni_f_imputation[var!='regime_row_owid',]
mm1_f_imputation<-mm1_f_imputation[var!='regime_row_owid',]
mm2_f_imputation<-mm2_f_imputation[var!='regime_row_owid',]

uni_f_imputation[,included:=ifelse(var %in% c('factor(regime_row_owid)1', 'factor(regime_row_owid)2', 'factor(regime_row_owid)3','HAQI_mean','day',
                                              'logit_lib', 'CPI_2022', 'Risk_comms_pct_hh_internet_3_6_1a'),'Excluded - multicollinearity',
                                   ifelse(var %in% mm1_f_imputation$var, 'Included', 'Excluded - significance'))]

uni_f_imputation<-uni_f_imputation[order(included)]

jpeg(paste0(dir,'/figures/univar_df_continuous_tobit_grouped.jpeg'), height=700, width=1000)
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
  coord_cartesian(xlim=c(-1,1))+
  facet_grid(included~., scales='free',switch='y')+
  theme_classic() +
  theme( panel.spacing=unit(2, "lines")
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20)
  )
dev.off()

jpeg(paste0(dir,'/figures/multivar01_df_continuous_tobit_grouped.jpeg'), height=700, width=1000)
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
  )+coord_cartesian(xlim=c(-0.5,0.5))+
  theme_classic() +
  theme( panel.spacing=unit(2, "lines")
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20)
  )
dev.off()

