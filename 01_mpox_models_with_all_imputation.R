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

rm(list=ls())

#set seed
set.seed(1007)
#source data cleaning script
source("~/mpox_codes/00_data_import_cleaning.R")
       
#####
## LISTWISE DELETION ANALYSIS, DICHOTOMOUS
######

lm1<-glm(mpox_greater~CPI_2022,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm2<-glm(mpox_greater~GAI,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit")) 
lm3<-glm(mpox_greater~HAQI_mean,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm4<-glm(mpox_greater~log_GDP,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm5<-glm(mpox_greater~logit_fem_edu,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit")) 
lm6<-glm(mpox_greater~Overall,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm7<-glm(mpox_greater~Risk_comm_3_5,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm8<-glm(mpox_greater~factor_pop_inclusion_riskcomm,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm9<-glm(mpox_greater~factor_misinfo_riskcomm,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm10<-glm(mpox_greater~logit_gender_internet_gap,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm11<-glm(mpox_greater~Risk_coms_mobile_subscribers_3_6_2,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm12<-glm(mpox_greater~logit_gender_phone_gap,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm13<-glm(mpox_greater~Risk_comms_pct_hh_internet_3_6_1a,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm14<-glm(mpox_greater~factor(regime_row_owid),data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm15<-glm(mpox_greater~day,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm16<-glm(mpox_greater~electdem_vdem_owid,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm17<-glm(mpox_greater~logit_lib,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))
lm18<-glm(mpox_greater~factor(ever_mpox),data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))

#test each var if significant and pull one beta estimate from multivariable normal model
df<-data.frame()
betas_orig<-data.frame()
for(mod in 1:18){
  model<-paste0('lm',mod)
  rows<-nrow(summary(get(model))$coef)
  if(rows>2){
    rows<-2:rows
  }
  ests<-data.frame(var=rownames(summary(get(model))$coef)[rows],
                   est=summary(get(model))$coef[rows,1],
                   lci=summary(get(model))$coef[rows,1]-1.96*summary(get(model))$coef[-1,2],
                   uci=summary(get(model))$coef[rows,1]+1.96*summary(get(model))$coef[-1,2],
                   p=summary(get(model))$coef[rows,4])
  df<-rbind(df,ests) 
  #beta hat estimate for each model
  beta_ests<-data.frame(var=rep(rownames(summary(get(model))$coef)[rows],10),
                        beta=get_betas(get(model),ndraws=10))
  beta_ests<-beta_ests[,-2]
  names(beta_ests)<-c('var','beta')
  if (length(rows)==3){
  beta_ests2<-melt(beta_ests, id.var='var')
  setDT(beta_ests2)
  beta_ests<-beta_ests2[,variable:=NULL]
  names(beta_ests)<-c('var','beta')
  }
  betas_orig<-rbind(betas_orig, beta_ests)
}

setDT(df)
setDT(betas_orig)
betas_orig[,lci:=quantile(beta, c(0.025)),by=c('var')]
betas_orig[,uci:=quantile(beta, c(0.975)),by=c('var')]
betas_orig[,est:=quantile(beta, c(0.5)),by=c('var')]

df_orig<-betas_orig[,c('var','lci','uci','est')]
df_orig<-df_orig[!duplicated(df_orig)]
#clean up names for a nicer plot
df_orig[,Variable:=factor(var,
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
#color coding for significance
df_orig[,signif:=ifelse(lci<0 & uci<0, 'Significant - negative',
                               ifelse(lci>0 & uci>0, 'Significant - positive',
                                      'Non-significant'))]

#include only variables p<0.1 and p<0.05 in the two models below:
#first remove highly collinear variables from the significance list
df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day',
                   'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]

#keep those significant at 0.1
keep<-df[p<0.1,var]
#Clean up namign for the formula vs results
keep<-gsub("Yes","",keep)
keep<-gsub("[)]1",")", keep)
formula<- paste0("mpox_greater~",paste(keep, collapse = "+"))
mm1<-glm(formula,data=full_dat[NAME %in% listwise_del_ctrys],family = binomial(link = "logit"))

#add in uncertainty 
mm1_uncert<-data.frame(get_betas(mm1,ndraws=10))
mm1_uncert_m<-melt(mm1_uncert[,-1])
names(mm1_uncert_m)<-c('var','beta')
setDT(df)
setDT(mm1_uncert_m)
mm1_uncert_m[,lci:=quantile(beta, c(0.025)),by=c('var')]
mm1_uncert_m[,uci:=quantile(beta, c(0.975)),by=c('var')]
mm1_uncert_m[,est:=quantile(beta, c(0.5)),by=c('var')]

mm1_uncert_m<-mm1_uncert_m[,c('var','lci','uci','est')]
mm1_ests<-mm1_uncert_m[!duplicated(mm1_uncert_m)]

mm1_ests[,Variable:=factor(var,
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

mm1_ests[,signif:=ifelse(lci<0 & uci<0, 'Significant - negative',
                   ifelse(lci>0 & uci>0, 'Significant - positive',
                          'Non-significant'))]
df_orig[,included:=ifelse(var %in% c('factor(regime_row_owid)1', 'factor(regime_row_owid)2', 'factor(regime_row_owid)3','HAQI_mean','day',
                                            'logit_lib', 'CPI_2022', 'Risk_comms_pct_hh_internet_3_6_1a'),'Excluded - multicollinearity',
                                 ifelse(var %in% mm1_ests$var, 'Included', 'Excluded - significance'))]

df_orig<-df_orig[order(included)]
df_orig[included=='Excluded - significance', signif:='Non-significant']


#make a tweak for ever mpox so the line doesn't show up on the figure:
df_orig[var=='factor(ever_mpox)1', `:=`(lci=-20, uci=-20, est=-20)]
jpeg(paste0(dir,'/figures/univar_df_dichotomous_noimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=df_orig)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lci, xmax=uci, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Green'))+
  ylab('')+xlab('Point estimate (95% UI)')+
  coord_cartesian(xlim=c(-5,5))+
  facet_grid(included~., scales='free',switch='y')+
  theme_classic() +
  theme( panel.spacing=unit(2, "lines")
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20)
  )

dev.off()

jpeg(paste0(dir,'/figures/multivar01_df_dichotomous_noimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=mm1_ests)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lci, xmax=uci, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Green'))+
  ylab('')+xlab('Point estimate (95% UI)')+
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+coord_cartesian(xlim=c(-5,5))+
  theme_classic() +
  theme( panel.spacing=unit(2, "lines")
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20)
  )
dev.off()

######
### Imputation for GAI only; 166 countries, dichotomous outcome
######

###PMM imputation for 166 countries (just missing GAI)
# Can use the pre-cleaned dataset here because GAI is the only thing being imputed, and we dont transform or clean GAI
data_final_166<-full_dat[NAME %in% final_166_ctrys]
#drop the date variable from the data so that we can impute (can't impute dates)
data_final_166[,`:=`(date_min=NULL)]
#create a dummy for missing day so we don't impute it as it's not really missing
data_final_166[is.na(day),day:=999]

#use predictive mean matching method and 10 imputations
imp<-mice(data_final_166, method='pmm', m=10)

#plots of raw vs imputed data for GDP vs GAI, as an example (since GDP is fully realized) - very consistent across imputations
ggmice(data_final_166, aes(y=log_GDP, x=GAI))+geom_point()
ggmice(imp, aes(y=log_GDP, x=GAI))+geom_point()
ggmice(imp, aes(y=log_GDP, x=GAI))+geom_point()+facet_grid(~.imp)

#analyze at each imputation then summarize across all imputations
df_uni<-data.frame()
betas_imputed_uni<-data.frame()
ests_imputed_mm1<-data.frame()
ests_imputed_mm2<-data.frame()
betas_imputed_mm1<-data.frame()
betas_imputed_mm2<-data.frame()
for(x in 1:10){
  imputed<-complete(imp, x)
  setDT(imputed)
  #removing the dummy coding for day so that we dont create an erroneous association in our Univar model
  imputed[day==999, day:=NA] 
  
  #logistic regression, univariable models
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

  
  #pull p value, estimates and confidence interval from model for within-imputation multivar model selection
  #pull a draw from a multivariable normal model for summarizing across all imputations
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
  
  df_uni<-rbind(df_uni,df)
  betas_imputed_uni<-rbind(betas_imputed_uni, betas_uni)
  
  #remove highly collinear variables from multivariable model
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day',
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  keep<-df[p<0.1,var]
  formula<- paste0("mpox_greater~",paste(keep, collapse = "+"))
  mm1<-glm(formula,data=imputed,family = binomial(link = "logit"))
  
  mm1_ests<-data.frame(imp=x,
                   var=rownames(summary(mm1)$coef)[-1],
                   est=summary(mm1)$coef[-1,1],
                   lci=summary(mm1)$coef[-1,1]-1.96*summary(mm1)$coef[-1,2],
                   uci=summary(mm1)$coef[-1,1]+1.96*summary(mm1)$coef[-1,2],
                   p=summary(mm1)$coef[-1,4])
  #multivariable model beta hat estimates for each covariate
  beta_mm1<-data.frame(imp=x,
                       var=rownames(summary(mm1)$coef)[-1],
                        beta=get_betas(mm1,1)[-1])
  ests_imputed_mm1<-rbind(ests_imputed_mm1, mm1_ests)
  betas_imputed_mm1<-rbind(betas_imputed_mm1,beta_mm1)

  #again, remove highly collinear variables 
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day',
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  keep<-df[p<0.05,var]
  formula<- paste0("mpox_greater~",paste(keep, collapse = "+"))
  mm2<-glm(formula,data=imputed,family = binomial(link = "logit"))
  #multivariable model beta hat estimates for each covariate
  mm2_ests<-data.frame(imp=x,
                       var=rownames(summary(mm2)$coef)[-1],
                       est=summary(mm2)$coef[-1,1],
                       lci=summary(mm2)$coef[-1,1]-1.96*summary(mm2)$coef[-1,2],
                       uci=summary(mm2)$coef[-1,1]+1.96*summary(mm2)$coef[-1,2],
                       p=summary(mm2)$coef[-1,4])
  beta_mm2<-data.frame(imp=x,
                       var=rownames(summary(mm2)$coef)[-1],
                       beta=get_betas(mm2,1)[-1])
  ests_imputed_mm2<-rbind(ests_imputed_mm2, mm2_ests)
  betas_imputed_mm2<-rbind(betas_imputed_mm2,beta_mm2)
}
setDT(betas_imputed_mm1)
#summarize across all imputations
betas_imputed_mm1[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm1[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm1[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm1_imputation<-betas_imputed_mm1[,c('var','lower_est','upper_est','median_est')]
mm1_imputation<-mm1_imputation[!duplicated(mm1_imputation)]

setDT(betas_imputed_mm2)
betas_imputed_mm2[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm2[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm2[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm2_imputation<-betas_imputed_mm2[,c('var','lower_est','upper_est','median_est')]
mm2_imputation<-mm2_imputation[!duplicated(mm2_imputation)]

setDT(betas_imputed_uni)
betas_imputed_uni[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_uni[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_uni[,median_est:=quantile(beta, c(0.5)),by=c('var')]

uni_imputation<-betas_imputed_uni[,c('var','lower_est','upper_est','median_est')]
uni_imputation<-uni_imputation[!duplicated(uni_imputation)]

#summary stats across imputations to build Table 1
impute_list<-data.frame()
for(x in 1:10){
  impute<-complete(imp,x)
  setDT(impute)
  impute[day==999, day:=NA]
  
  impute_list<-rbind(impute_list,impute)
}

##### 
#### FULL DATA IMPUTATION, Dichotomous model, 184 observations
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
  
  #remove highly collinear variables
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day',
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

  #remove highly collinear variables
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022','day', 
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
setDT(betas_imputed_mm1_f)
#summarize across imputations
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
uni_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                               ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                      'Non-significant'))]
uni_imputation[,Variable:=factor(var,
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
mm1_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                               ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                      'Non-significant'))]
mm1_imputation[,Variable:=factor(var,
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
mm2_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                               ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                      'Non-significant'))]
mm2_imputation[,Variable:=factor(var,
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

uni_imputation[,included:=ifelse(var %in% c('factor(regime_row_owid)1', 'factor(regime_row_owid)2', 'factor(regime_row_owid)3','HAQI_mean','day',
                                              'logit_lib', 'CPI_2022', 'Risk_comms_pct_hh_internet_3_6_1a'),'Excluded - multicollinearity',
                                   ifelse(var %in% mm1_imputation$var, 'Included', 'Excluded - significance'))]

uni_imputation<-uni_imputation[order(included)]
#make a tweak for ever mpox so the line doesn't show up on the figure:
uni_imputation[var=='factor(ever_mpox)1', `:=`(lower_est=-20, upper_est=-20, median_est=-20)]
#GAI imputation only, univariable model
jpeg(paste0(dir,'/figures/univar_df_dichotomous_GAIimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=uni_imputation)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Red','Green'))+
  ylab('')+xlab('Point estimate (95% UI)')+
  # theme(
  #   panel.border = element_blank(), 
  #   panel.grid.major.y = element_line(color = "light grey", size = 0.3),
  #   panel.grid.major.x = element_line(color = "light grey", size = 0.3),
  #   panel.grid.minor = element_blank(), 
  #   axis.line = element_blank())
  coord_cartesian(xlim=c(-5,5))+
  facet_grid(included~., scales='free',switch='y')+
  theme_classic() +
  theme(panel.spacing=unit(2, "lines"),
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20))
dev.off()

jpeg(paste0(dir,'/figures/multivar01_df_dichotomous_GAIimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=mm1_imputation)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Red','Green'))+
  ylab('')+xlab('Point estimate (95% UI)')+
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+coord_cartesian(xlim=c(-5,5))+
  theme_classic() +
  theme(panel.spacing=unit(2, "lines"),
        , strip.placement.y = "outside"
        , strip.background = element_blank()
        , strip.text = element_text(face = "bold"),
        text=element_text(size=20))
dev.off()

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

jpeg(paste0(dir,'/figures/univar_df_dichotomous_fullimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=uni_f_imputation)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Green'))+
  ylab('')+xlab('Point estimate (95% UI)')+
  # theme(
  #   panel.border = element_blank(), 
  #   panel.grid.major.y = element_line(color = "light grey", size = 0.3),
  #   panel.grid.major.x = element_line(color = "light grey", size = 0.3),
  #   panel.grid.minor = element_blank(), 
  #   axis.line = element_blank())
  coord_cartesian(xlim=c(-5,5))+
  facet_grid(included~., scales='free',switch='y')+
  theme_classic() +
  theme( panel.spacing=unit(2, "lines")
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20)
  )

dev.off()

jpeg(paste0(dir,'/figures/multivar01_df_dichotomous_fullimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=mm1_f_imputation)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Red','Green'))+
  ylab('')+xlab('Point estimate (95% UI)')+
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+coord_cartesian(xlim=c(-5,5))+
  theme_classic() +
  theme( panel.spacing=unit(2, "lines")
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20)
  )
dev.off()

jpeg(paste0(dir,'/figures/multivar005_df_dichotomous_fullimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=mm2_f_imputation)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Red','Green'))+
  ylab('')+xlab('Point estimate (95% UI)')+
  theme(
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "light grey", size = 0.3),
    panel.grid.major.x = element_line(color = "light grey", size = 0.3),
    panel.grid.minor = element_blank(), 
    axis.line = element_blank()
  )+coord_cartesian(xlim=c(-5,5))+
  theme_classic() +
  theme( panel.spacing=unit(2, "lines")
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20)
  )
dev.off()

#summary stats across imputations
impute_list2<-data.frame()
for(x in 1:10){
impute<-complete(imp,x)
setDT(impute)
impute[day==999, day:=NA]
impute[,log_GDP:=log(yr2022)]

  impute[,norm_fem_edu:=(Female_edu_mean_yrs_25_29-min(Female_edu_mean_yrs_25_29, na.rm=T))/(max(Female_edu_mean_yrs_25_29, na.rm=T)-min(Female_edu_mean_yrs_25_29, na.rm=T))]
  impute[,logit_fem_edu:=log((norm_fem_edu)/(1-norm_fem_edu))]
  impute[norm_fem_edu %in% c(0,1),logit_fem_edu:=log((norm_fem_edu+(0.5/184))/(1-norm_fem_edu+(0.5/184)))]
  
  impute[,norm_gender_phone_gap:=(Risk_comms_gender_gap_access_phone_3_6_3a)/100]
  impute[,logit_gender_phone_gap:=log((norm_gender_phone_gap)/(1-norm_gender_phone_gap))]
  impute[norm_gender_phone_gap %in% c(0,1),logit_gender_phone_gap:=log((norm_gender_phone_gap+(0.5/184))/(1-norm_gender_phone_gap+(0.5/184)))]
  
  impute[,norm_gender_internet_gap:=(Risk_comms_gender_gap_access_internet_3_6_4_a)/100]
  impute[,logit_gender_internet_gap:=log((norm_gender_internet_gap+(0.5/184))/(1-norm_gender_internet_gap+(0.5/184)))]
  impute[norm_gender_internet_gap %in% c(0,1),logit_gender_internet_gap:=log((norm_gender_internet_gap+(0.5/184))/(1-norm_gender_internet_gap+(0.5/184)))]
  
  impute[,logit_lib:=log((lib_vdem_owid)/(1-lib_vdem_owid))]
  impute[lib_vdem_owid %in% c(0,1),logit_lib:=log((lib_vdem_owid+(0.5/184))/(1-lib_vdem_owid+(0.5/184)))]
  

impute[,factor_pop_inclusion_riskcomm:=factor(Risk_comms_pop_inclusion_3_5_1b, levels=c(0,100),
                                               labels=c('No',"Yes"))]
impute[,factor_misinfo_riskcomm:=factor(Risk_comms_leader_share_misinfo_3_5_2b, levels=c(0,100),
                                         labels=c('No',"Yes"))]

impute_list2<-rbind(impute_list2,impute)
}

######
### CONTINUOUS OUTCOME: % mpox value
#######

######
### LISTWISE DELETION ANALYSIS
#######

lm1<-glm(avg_pct_mpox~CPI_2022,data=full_dat[NAME %in% listwise_del_ctrys])
lm2<-glm(avg_pct_mpox~GAI,data=full_dat[NAME %in% listwise_del_ctrys]) 
lm3<-glm(avg_pct_mpox~HAQI_mean,data=full_dat[NAME %in% listwise_del_ctrys])
lm4<-glm(avg_pct_mpox~log_GDP,data=full_dat[NAME %in% listwise_del_ctrys])
lm5<-glm(avg_pct_mpox~logit_fem_edu,data=full_dat[NAME %in% listwise_del_ctrys]) 
lm6<-glm(avg_pct_mpox~Overall,data=full_dat[NAME %in% listwise_del_ctrys])
lm7<-glm(avg_pct_mpox~Risk_comm_3_5,data=full_dat[NAME %in% listwise_del_ctrys])
lm8<-glm(avg_pct_mpox~factor_pop_inclusion_riskcomm,data=full_dat[NAME %in% listwise_del_ctrys])
lm9<-glm(avg_pct_mpox~factor_misinfo_riskcomm,data=full_dat[NAME %in% listwise_del_ctrys])
lm10<-glm(avg_pct_mpox~logit_gender_internet_gap,data=full_dat[NAME %in% listwise_del_ctrys])
lm11<-glm(avg_pct_mpox~Risk_coms_mobile_subscribers_3_6_2,data=full_dat[NAME %in% listwise_del_ctrys])
lm12<-glm(avg_pct_mpox~logit_gender_phone_gap,data=full_dat[NAME %in% listwise_del_ctrys])
lm13<-glm(avg_pct_mpox~Risk_comms_pct_hh_internet_3_6_1a,data=full_dat[NAME %in% listwise_del_ctrys])
lm14<-glm(avg_pct_mpox~day,data=full_dat[NAME %in% listwise_del_ctrys])
lm15<-glm(avg_pct_mpox~factor(regime_row_owid),data=full_dat[NAME %in% listwise_del_ctrys])
lm16<-glm(avg_pct_mpox~electdem_vdem_owid,data=full_dat[NAME %in% listwise_del_ctrys])
lm17<-glm(avg_pct_mpox~logit_lib,data=full_dat[NAME %in% listwise_del_ctrys])
lm18<-glm(avg_pct_mpox~factor(ever_mpox),data=full_dat[NAME %in% listwise_del_ctrys])

#test each var if signif:
df<-data.frame()
betas_orig<-data.frame()
for(mod in 1:18){
  model<-paste0('lm',mod)
  rows<-nrow(summary(get(model))$coef)
  if(rows>2){
    rows<-2:rows
  }
  ests<-data.frame(var=rownames(summary(get(model))$coef)[rows],
                   est=summary(get(model))$coef[rows,1],
                   lci=summary(get(model))$coef[rows,1]-1.96*summary(get(model))$coef[-1,2],
                   uci=summary(get(model))$coef[rows,1]+1.96*summary(get(model))$coef[-1,2],
                   p=summary(get(model))$coef[rows,4])
  df<-rbind(df,ests) 
  #beta hat estimate for each model
  beta_ests<-data.frame(var=rep(rownames(summary(get(model))$coef)[rows],10),
                        beta=get_betas(get(model),ndraws=10))
  beta_ests<-beta_ests[,-2]
  names(beta_ests)<-c('var','beta')
  if (length(rows)==3){
    beta_ests2<-melt(beta_ests, id.var='var')
    setDT(beta_ests2)
    beta_ests<-beta_ests2[,variable:=NULL]
    names(beta_ests)<-c('var','beta')
  }
  betas_orig<-rbind(betas_orig, beta_ests)
}

setDT(df)
setDT(betas_orig)
betas_orig[,lci:=quantile(beta, c(0.025)),by=c('var')]
betas_orig[,uci:=quantile(beta, c(0.975)),by=c('var')]
betas_orig[,est:=quantile(beta, c(0.5)),by=c('var')]

df_orig<-betas_orig[,c('var','lci','uci','est')]
df_orig<-df_orig[!duplicated(df_orig)]



df_orig[,Variable:=factor(var,
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
df_orig[,signif:=ifelse(lci<0 & uci<0, 'Significant - negative',
                        ifelse(lci>0 & uci>0, 'Significant - positive',
                               'Non-significant'))]

#Drop highly collinear variables
df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022','day',
                   'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]


keep<-df[p<0.1,var]
keep<-gsub("Yes","",keep)
keep<-gsub("[)]1",")", keep)
formula<- paste0("avg_pct_mpox~",paste(keep, collapse = "+"))
mm1<-glm(formula,data=full_dat[NAME %in% listwise_del_ctrys])

#add in uncertainty again
mm1_uncert<-data.frame(get_betas(mm1,ndraws=10))
mm1_uncert_m<-melt(mm1_uncert[,-1])
names(mm1_uncert_m)<-c('var','beta')
setDT(df)
setDT(mm1_uncert_m)
mm1_uncert_m[,lci:=quantile(beta, c(0.025)),by=c('var')]
mm1_uncert_m[,uci:=quantile(beta, c(0.975)),by=c('var')]
mm1_uncert_m[,est:=quantile(beta, c(0.5)),by=c('var')]

mm1_uncert_m<-mm1_uncert_m[,c('var','lci','uci','est')]
mm1_ests<-mm1_uncert_m[!duplicated(mm1_uncert_m)]
mm1_ests[var=='factor.ever_mpox.1',var:='factor(ever_mpox)1']
mm1_ests[,Variable:=factor(var,
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
mm1_ests[,signif:=ifelse(lci<0 & uci<0, 'Significant - negative',
                         ifelse(lci>0 & uci>0, 'Significant - positive',
                                'Non-significant'))]

df_orig[,included:=ifelse(var %in% c('factor(regime_row_owid)1', 'factor(regime_row_owid)2', 'factor(regime_row_owid)3','HAQI_mean','day',
                                     'logit_lib', 'CPI_2022', 'Risk_comms_pct_hh_internet_3_6_1a'),'Excluded - multicollinearity',
                          ifelse(var %in% mm1_ests$var, 'Included', 'Excluded - significance'))]

df_orig<-df_orig[order(included)]
df_orig[included=='Excluded - significance', signif:='Non-significant']

jpeg(paste0(dir,'/figures/univar_df_continuous_noimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=df_orig)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lci, xmax=uci, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Red','Green'))+
  ylab('')+xlab('Point estimate (95% UI)')+
  # theme(
  #   panel.border = element_blank(), 
  #   panel.grid.major.y = element_line(color = "light grey", size = 0.3),
  #   panel.grid.major.x = element_line(color = "light grey", size = 0.3),
  #   panel.grid.minor = element_blank(), 
  #   axis.line = element_blank())
  coord_cartesian(xlim=c(-0.5,0.5))+
  facet_grid(included~., scales='free',switch='y')+
  theme_classic() +
  theme( panel.spacing=unit(2, "lines")
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20)
  )

dev.off()
jpeg(paste0(dir,'/figures/multivar01_df_continuous_noimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=mm1_ests)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lci, xmax=uci, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Red','Green'))+
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

keep<-df[p<0.05,var]
keep<-gsub("Yes","",keep)
keep<-gsub("[)]1",")", keep)
formula<- paste0("avg_pct_mpox~",paste(keep, collapse = "+"))


mm2<-glm(formula,data=full_dat[NAME %in% listwise_del_ctrys])
#add in uncertainty again
mm2_uncert<-data.frame(get_betas(mm2,ndraws=10))
mm2_uncert_m<-melt(mm2_uncert[,-1])
names(mm2_uncert_m)<-c('var','beta')
setDT(df)
setDT(mm2_uncert_m)
mm2_uncert_m[,lci:=quantile(beta, c(0.025)),by=c('var')]
mm2_uncert_m[,uci:=quantile(beta, c(0.975)),by=c('var')]
mm2_uncert_m[,est:=quantile(beta, c(0.5)),by=c('var')]

mm2_uncert_m<-mm2_uncert_m[,c('var','lci','uci','est')]
mm2_ests<-mm2_uncert_m[!duplicated(mm2_uncert_m)]
mm2_ests[var=='factor.ever_mpox.1',var:='factor(ever_mpox)1']

mm2_ests[,Variable:=factor(var,
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
mm2_ests[,signif:=ifelse(lci<0 & uci<0, 'Significant - negative',
                         ifelse(lci>0 & uci>0, 'Significant - positive',
                                'Non-significant'))]


######
##### GAI IMPUTATION,166 countries
#####

data_final_166<-full_dat[NAME %in% final_166_ctrys]
#drop the date variable from the data so that we can impute (can't impute dates)
data_final_166[,`:=`(date_min=NULL)]
#create a dummy for missing day so we don't impute it as it's not really missing
data_final_166[is.na(day),day:=999]

imp<-mice(data_final_166, method='pmm', m=10)

#analyze at each imputation
df_uni<-data.frame()
betas_imputed_uni<-data.frame()
ests_imputed_mm1<-data.frame()
ests_imputed_mm2<-data.frame()
betas_imputed_mm1<-data.frame()
betas_imputed_mm2<-data.frame()
for(x in 1:10){
  imputed<-complete(imp, x)
  setDT(imputed)
  imputed[day==999, day:=NA]
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
  
  df_uni<-rbind(df_uni,df)
  betas_imputed_uni<-rbind(betas_imputed_uni, betas_uni)
  
  # remove highly collinear variables from multivariable model
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day',
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  keep<-df[p<0.1,var]
  formula<- paste0("avg_pct_mpox~",paste(keep, collapse = "+"))
  formula<-gsub('Yes','',formula)
  formula<-gsub(')1',')',formula)
  
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
  ests_imputed_mm1<-rbind(ests_imputed_mm1, mm1_ests)
  betas_imputed_mm1<-rbind(betas_imputed_mm1,beta_mm1)

  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day',
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  keep<-df[p<0.05,var]
  formula<- paste0("avg_pct_mpox~",paste(keep, collapse = "+"))
  formula<-gsub('Yes','',formula)
  formula<-gsub(')1',')',formula)
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
  ests_imputed_mm2<-rbind(ests_imputed_mm2, mm2_ests)
  betas_imputed_mm2<-rbind(betas_imputed_mm2,beta_mm2)
}
setDT(betas_imputed_mm1)
betas_imputed_mm1[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm1[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm1[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm1_imputation<-betas_imputed_mm1[,c('var','lower_est','upper_est','median_est')]
mm1_imputation<-mm1_imputation[!duplicated(mm1_imputation)]


setDT(betas_imputed_mm2)
betas_imputed_mm2[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_mm2[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_mm2[,median_est:=quantile(beta, c(0.5)),by=c('var')]

mm2_imputation<-betas_imputed_mm2[,c('var','lower_est','upper_est','median_est')]
mm2_imputation<-mm2_imputation[!duplicated(mm2_imputation)]

setDT(betas_imputed_uni)
betas_imputed_uni[,lower_est:=quantile(beta, c(0.025)),by=c('var')]
betas_imputed_uni[,upper_est:=quantile(beta, c(0.975)),by=c('var')]
betas_imputed_uni[,median_est:=quantile(beta, c(0.5)),by=c('var')]

uni_imputation<-betas_imputed_uni[,c('var','lower_est','upper_est','median_est')]


##### 
#### FULL DATA IMPUTATION, 184 observations
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
  #remove highly correlated variables
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day',
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
 
  df<-df[!var %in% c('Risk_comms_pct_hh_internet_3_6_1a','HAQI_mean', 'CPI_2022', 'day',
                     'logit_lib','factor(regime_row_owid)1','factor(regime_row_owid)2','factor(regime_row_owid)3')]
  keep<-df[p<0.05,var]
  formula<- paste0("avg_pct_mpox~",paste(keep, collapse = "+"))
  formula<-gsub("Yes","",formula)
  formula<-gsub("[)]1",")", formula)
  mm2<-glm(formula,data=imputed)
  rmm2<-pR2(mm2)
  
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
uni_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                               ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                      'Non-significant'))]
uni_imputation[var=='factor.ever_mpox.1',var:='factor(ever_mpox)1']
mm1_imputation[var=='factor.ever_mpox.1',var:='factor(ever_mpox)1']
mm2_imputation[var=='factor.ever_mpox.1',var:='factor(ever_mpox)1']

uni_imputation[,Variable:=factor(var,
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
mm1_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                               ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                      'Non-significant'))]
mm1_imputation[,Variable:=factor(var,
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
mm2_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                               ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                      'Non-significant'))]
mm2_imputation[,Variable:=factor(var,
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

uni_imputation[,included:=ifelse(var %in% c('factor(regime_row_owid)1', 'factor(regime_row_owid)2', 'factor(regime_row_owid)3','HAQI_mean','day',
                                            'logit_lib', 'CPI_2022', 'Risk_comms_pct_hh_internet_3_6_1a'),'Excluded - multicollinearity',
                                 ifelse(var %in% mm1_imputation$var, 'Included', 'Excluded - significance'))]

uni_imputation<-uni_imputation[order(included)]

jpeg(paste0(dir,'/figures/univar_df_continuous_GAIimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=uni_imputation)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Red','Green'))+
  ylab('')+xlab('Point estimate (95% UI)')+
   coord_cartesian(xlim=c(-0.5,0.5))+
  facet_grid(included~., scales='free',switch='y')+
  theme_classic() +
  theme( panel.spacing=unit(2, "lines")
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20)
  )
dev.off()

jpeg(paste0(dir,'/figures/multivar01_df_continuous_GAIimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=mm1_imputation)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Red','Green'))+
  ylab('')+xlab('Point estimate (95% UI)')+
  coord_cartesian(xlim=c(-0.5,0.5))+
  theme_classic() +
  theme( panel.spacing=unit(2, "lines")
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20)
  )
dev.off()

#color coding & clean up variable names
uni_f_imputation[,signif:=ifelse(lower_est<0 & upper_est<0, 'Significant - negative',
                                 ifelse(lower_est>0 & upper_est>0, 'Significant - positive',
                                        'Non-significant'))]

uni_f_imputation[var=='factor.ever_mpox.1',var:='factor(ever_mpox)1']
mm1_f_imputation[var=='factor.ever_mpox.1',var:='factor(ever_mpox)1']
mm2_f_imputation[var=='factor.ever_mpox.1',var:='factor(ever_mpox)1']

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

jpeg(paste0(dir,'/figures/univar_df_continuous_fullimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=uni_f_imputation)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Red','Green'))+
  ylab('')+xlab('Point estimate (95% UI)')+
   coord_cartesian(xlim=c(-0.5,0.5))+
  facet_grid(included~., scales='free',switch='y')+
  theme_classic() +
  theme( panel.spacing=unit(2, "lines")
         , strip.placement.y = "outside"
         , strip.background = element_blank()
         , strip.text = element_text(face = "bold"),
         text=element_text(size=20)
  )
dev.off()

jpeg(paste0(dir,'/figures/multivar01_df_continuous_fullimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=mm1_f_imputation)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Red','Green'))+
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

jpeg(paste0(dir,'/figures/multivar005_df_continuous_fullimputation_grouped.jpeg'), height=700, width=1000)
ggplot(data=mm2_f_imputation)+geom_vline(xintercept = 0, col='black',lty=2)+
  geom_point(aes(x=median_est,y=Variable, col=signif), cex=6)+
  geom_errorbarh(aes(xmin=lower_est, xmax=upper_est, y=Variable, col=signif), height=0)+
  theme_bw()+scale_color_manual('Significance at p<0.05', values=c('Black','Red','Green'))+
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
