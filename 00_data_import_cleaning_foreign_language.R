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

#specify your repo directory here
dir<-""

#set seed
set.seed(1007)

#get a beta estimate from a linear model
get_betas<-function(model_fit, ndraws){
  betas = coef(model_fit)
  vcov = vcov(model_fit)
  betas_simulated = MASS::mvrnorm(ndraws, betas, vcov)
  return(betas_simulated)
}

#bring in each country individually
fils<-list.files(paste0(dir, "raw data/foreign_official/"), pattern='.csv', full.names = T)
data_foreign<-data.frame(Week=as.Date('11/01/2022'), mpox=NA, monkeypox=NA,monkeypox_foreign=NA, country=NA, language=NA)
for(x in 1:length(fils)){
  ctry<-gsub(paste0(dir, "raw data/foreign_official/"),"",fils[x])
  vect<-strsplit(ctry,"_")[1]
  ctry<-vect[[1]][1]
  language<-vect[[1]][2]
  language<-gsub(".csv","",language)
  dat<-fread(fils[x])
  dat$Week<-as.Date(dat$Date)
  dat[,Date:=NULL]
  dat[,country:=ctry]
  dat[,language:=language]
  names(dat)<-c('monkeypox','mpox','monkeypox_foreign','Week','country','language')
  data_foreign<-rbind(data_foreign, dat)
}
data_foreign<-data_foreign[!is.na(data_foreign$country),]
setDT(data_foreign)
data_foreign[,monkeypox:=NULL]
data_foreign[,monkeypox:=monkeypox_foreign]
data_foreign[mpox=='<1', mpox:=0]
data_foreign[monkeypox=='<1', monkeypox:=0]
data_foreign[,mpox:=as.numeric(mpox)]
data_foreign[,monkeypox:=as.numeric(monkeypox)]


data_foreign[,avg_mpox:=mean(mpox, na.rm=T), by='country']
data_foreign[,avg_monkeypox:=mean(monkeypox, na.rm=T), by='country']

data_foreign[,mpox_greater:=0]
data_foreign[avg_mpox > avg_monkeypox,mpox_greater:=1]

data_foreign[,mpox_non0:=0]
data_foreign[avg_mpox > 0,mpox_non0:=1]

data_foreign[,pct_mpox:=mpox/(mpox+monkeypox)]
#several 0/0, just make sure to do NA.RM for average pct
data_foreign[,avg_pct_mpox:=mean(pct_mpox,na.rm=T), by='country']
data_foreign[is.na(avg_pct_mpox), avg_pct_mpox:=0]

#now bring in the original dataset for everywhere that isn't in the foreign language data set.
filelist<-list.files(paste0(dir, "raw data/cleaned_full_time/"), full.names = TRUE)
data<-data.frame(Week=as.Date('11/01/2022'), mpox=NA, monkeypox=NA, country=NA)
for(x in 1:length(filelist)){
  dat<-fread(filelist[x])
  dat$Week<-as.Date(dat$Week)
  data<-rbind(data, dat)
}
data<-data[!is.na(data$country),]
setDT(data)
data[mpox=='<1', mpox:=0]
data[monkeypox=='<1', monkeypox:=0]
data[,mpox:=as.numeric(mpox)]
data[,monkeypox:=as.numeric(monkeypox)]

#get average over entire time period for each term
data[,avg_mpox:=mean(mpox, na.rm=T), by='country']
data[,avg_monkeypox:=mean(monkeypox, na.rm=T), by='country']

#highlight countries where mpox > monkeypox, if any
data[,mpox_greater:=0]
data[avg_mpox > avg_monkeypox,mpox_greater:=1]

#and just highlight places where mpox at least >0
data[,mpox_non0:=0]
data[avg_mpox > 0,mpox_non0:=1]

#% average mpox of total use
#data[,mpox_noNA:=ifelse(is.na(mpox))]
data[,pct_mpox:=mpox/(mpox+monkeypox)]
#several 0/0, just make sure to do NA.RM for average pct
data[,avg_pct_mpox:=mean(pct_mpox,na.rm=T), by='country']
data[is.na(avg_pct_mpox), avg_pct_mpox:=0]

data_old<-data[!country %in% data_foreign$country]

data_foreign[,monkeypox_foreign:=NULL]
data_foreign[,language:=NULL]
data<-rbind(data_old, data_foreign)

data[,NAME:=country]
data[country=="Congo - Brazzaville", NAME:='Republic of Congo']
data[country=="Congo - Kinshasa", NAME:='Democratic Republic of the Congo']
data[country=="Czechia", NAME:='Czech Republic']
data[country=="MyanmarBurma", NAME:='Myanmar']
data[country=="St. Kitts and Nevis", NAME:='Saint Kitts and Nevis']
data[country=="St. Lucia", NAME:='Saint Lucia']
data[country== "St. Vincent and Grenadines", NAME:='Saint Vincent and the Grenadines']
data[country== "Cote d'Ivoire", NAME:="Ivory Coast"]
data[country== "Côte d’Ivoire", NAME:="Ivory Coast"]

small_dat<-copy(data)
small_dat<-small_dat[,c("NAME",'avg_mpox','avg_monkeypox','avg_pct_mpox')]
small_dat<-small_dat[!duplicated(small_dat)]
#remove Taiwan, Guadeloupe, Martinique, and Reunion from the filled in data
small_dat<-small_dat[!NAME %in% c('Taiwan', 'Guadeloupe','Martinique','Reunion')]



#bring in country level covariates
#corruptions perception index
CPI<-fread(paste0(dir,"raw data/CPI_2022.csv"))
#Standardize names
CPI[,NAME:=Country_name]
CPI[Country_name=='Cabo Verde',NAME:='Cape Verde']
CPI[Country_name=="Cote d'Ivoire",NAME:="Ivory Coast"]
CPI[Country_name=='Congo',NAME:='Republic of Congo']
CPI[Country_name=='Czechia',NAME:='Czech Republic']
CPI[Country_name=='Korea, South',NAME:='South Korea']
CPI[Country_name=='US',NAME:='United States']
CPI[Country_name=='Turkey',NAME:='Turkiye']
CPI[Country_name=='Guinea Bissau',NAME:='Guinea-Bissau']
CPI[,Country_name:=NULL]

#LGBTI acceptance index score
GAI<-fread(paste0(dir,"raw data/LGBTI Acceptance Index score in 2017-2020 - Sheet1.csv"))
#Standardize names
GAI[,NAME:=Country]
GAI[Country=="Cote d’Ivoire",NAME:='Ivory Coast']
GAI[Country=='Macau SAR',NAME:='Macao']
GAI[Country=='Great Britain',NAME:='United Kingdom']
GAI[Country=='Republic of the Congo',NAME:='Republic of Congo']
GAI[Country=='Bosnia Herzegovina',NAME:='Bosnia and Herzegovina']
GAI[Country=='Swaziland',NAME:='Eswatini']
GAI[Country=='Turkey',NAME:='Turkiye']
GAI[Country=='Somaliland',NAME:='Somalia']
GAI[,Country:=NULL]

#IHME Healthcare Access and Quality Index
HAQI<-fread(paste0(dir,"raw data/IHME_GBD_2015_HAQ_INDEX_1990_2015_HAQ_INDEX_AND_VALUES_Y2017M05D18.CSV"))
HAQI<-HAQI[year_id==2015 & indicator_name=='Healthcare Access and Quality',]

#Standardize names
HAQI[,NAME:=location_name]
HAQI[location_name=='Macedonia',NAME:='North Macedonia']
HAQI[location_name=='The Bahamas',NAME:='Bahamas']
HAQI[location_name=='Turkey',NAME:='Turkiye']
HAQI[location_name=='Congo',NAME:='Republic of Congo']
HAQI[location_name=='Swaziland',NAME:='Eswatini']
HAQI[location_name=="Cote d'Ivoire",NAME:="Ivory Coast"]
HAQI[location_name=="The Gambia",NAME:="Gambia"]
HAQI[,location_name:=NULL]
HAQI<-HAQI[,c('NAME','val')]
names(HAQI)<-c('NAME','HAQI_mean')

#GDP per capita
GDP<-fread(paste0(dir,"raw data/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_6224630.csv"))
#use most recent available year of data
GDP[is.na(yr2022), yr2022:=yr2021]
GDP[is.na(yr2022), yr2022:=yr2020] 
GDP[is.na(yr2022), yr2022:=yr2019] 
GDP[is.na(yr2022), yr2022:=yr2018] 
GDP[is.na(yr2022), yr2022:=yr2015] 
GDP[is.na(yr2022), yr2022:=yr2014] #everything else has no data at all after this date
GDP2<-GDP[,c(1,15)]
names(GDP2)<-c('Country','yr2022')
#Standardize names
GDP2[,NAME:=Country]
GDP2[Country=="Yemen, Rep.", NAME:='Yemen']
GDP2[Country=="Viet Nam", NAME:='Vietnam']
GDP2[Country=="Venezuela, RB" , NAME:='Venezuela']
GDP2[Country=="Bahamas, The", NAME:='Bahamas']
GDP2[Country=="Brunei Darussalam", NAME:='Brunei']
GDP2[Country=="Cabo Verde", NAME:="Cape Verde"]
GDP2[Country=="Congo, Rep.", NAME:="Republic of Congo"]
GDP2[Country=="Cote d'Ivoire", NAME:="Ivory Coast"]
GDP2[Country=="Congo, Dem. Rep.", NAME:="Democratic Republic of the Congo"]
GDP2[Country=="Czechia", NAME:="Czech Republic"]
GDP2[Country=="Egypt, Arab Rep.", NAME:="Egypt"]
GDP2[Country=="Gambia, The", NAME:="Gambia"]
GDP2[Country=="Iran, Islamic Rep.", NAME:="Iran" ]
GDP2[Country=="Lao PDR", NAME:= "Laos"]
GDP2[Country=="Hong Kong SAR, China", NAME:="Hong Kong"]
GDP2[Country=="Russian Federation",NAME:="Russia"]
GDP2[Country=="Slovak Republic",NAME:="Slovakia"]
GDP2[Country=="Korea, Rep.",NAME:="South Korea"]
GDP2[Country=="St. Kitts and Nevis", NAME:="Saint Kitts and Nevis"]
GDP2[Country=="St. Vincent and the Grenadines",NAME:="Saint Vincent and the Grenadines"]
GDP2[Country=="Syrian Arab Republic",NAME:="Syria"]
GDP2[Country=="Kyrgyz Republic", NAME:="Kyrgyzstan"]
GDP2[Country=="St. Lucia" , NAME:="Saint Lucia"]
GDP2[Country=="West Bank and Gaza" , NAME:="Palestine"]
GDP2[,Country:=NULL]

#IHME years of education
edu<-fread(paste0(dir,"raw data/IHME_EDUC_DISTRIBUTIONS_1970_2030_Y2020M04D15.csv"))
edu2<-edu[year==2020 & measure=="mean_years",]
edu2<-edu2[,c("location_name","sex","mean")]
edu_c<-dcast(edu2,location_name~sex)
names(edu_c)<-c('Country','Female_edu_mean_yrs_25_29',"Male_edu_mean_yrs_25_29")
#Standardize names
setDT(edu_c)
edu_c[,NAME:=Country]
edu_c[Country=="Cote d'Ivoire" , NAME:="Ivory Coast" ]
edu_c[Country=="Congo", NAME:="Republic of Congo"]
edu_c[Country=="The Bahamas" , NAME:="Bahamas" ]
edu_c[Country=="Swaziland", NAME:="Eswatini" ]
edu_c[Country=="The Gambia" , NAME:="Gambia" ]
edu_c[Country=="Turkey" , NAME:="Turkiye"]
edu_c[Country=="Macedonia"   , NAME:="North Macedonia"]
edu_c[,Country:=NULL]

#GSHI preparedness scores
GHSI<-fread(paste0(dir,"raw data/2021-GHS-Index-April-2022.csv"))
GHSI21<-GHSI[Year==2021,]
GHSI21<-GHSI21[,c(1,3,133,136,140,143,144,147,149)]
names(GHSI21)<-c('Country','Overall','Risk_comm_3_5',
                 'Risk_comms_pop_inclusion_3_5_1b',
                 'Risk_comms_leader_share_misinfo_3_5_2b',
                 'Risk_comms_pct_hh_internet_3_6_1a',
                 'Risk_coms_mobile_subscribers_3_6_2',
                 'Risk_comms_gender_gap_access_phone_3_6_3a',
                 'Risk_comms_gender_gap_access_internet_3_6_4_a')
#Standardize names
GHSI21[,NAME:=Country]
GHSI21[Country=="Antigua & Barbuda",NAME:='Antigua and Barbuda']
GHSI21[Country=="Bosnia and Hercegovina",NAME:="Bosnia and Herzegovina"]
GHSI21[Country=="Cabo Verde",NAME:="Cape Verde"]
GHSI21[Country=="Congo (Brazzaville)",NAME:="Republic of Congo"]
GHSI21[Country=="Congo (Democratic Republic)",NAME:="Democratic Republic of the Congo"]
GHSI21[Country=="Côte d'Ivoire",NAME:="Ivory Coast" ]
GHSI21[Country== "eSwatini",NAME:="Eswatini"]
GHSI21[Country=="Kyrgyz Republic" ,NAME:="Kyrgyzstan"]
GHSI21[Country== "St Kitts & Nevis",NAME:="Saint Kitts and Nevis"]
GHSI21[Country=="St Lucia",NAME:="Saint Lucia"]
GHSI21[Country=="St Vincent & The Grenadines",NAME:="Saint Vincent and the Grenadines"]
GHSI21[Country=="Turkey",NAME:="Turkiye"]
GHSI21[Country=="United States of America",NAME:='United States']
GHSI21[,Country:=NULL]

#democracy indicators
democ<-fread(paste0(dir,"raw data/vdem_row_final.csv"))
#smaller variable subset
democ_all<-democ[year>=2022,c('country_name','regime_row_owid',
                            'electdem_vdem_owid','lib_vdem_owid','region')]
democ_all<-democ_all[!is.na(regime_row_owid)]
#Standardize names
democ_all[,NAME:=country_name]
democ_all[country_name=='Congo',NAME:='Republic of Congo']
democ_all[country_name=="Cote d'Ivoire",NAME:= "Ivory Coast" ]
democ_all[country_name=="Czechia",NAME:='Czech Republic']
democ_all[country_name=="Democratic Republic of Congo",NAME:='Democratic Republic of the Congo']
democ_all[country_name=="Palestine/Gaza",NAME:='Palestine']
democ_all[country_name=="Turkey" ,NAME:='Turkiye']
democ_all[country_name=="Timor" ,NAME:="Timor-Leste"]
democ_all[,country_name:=NULL]

#merge covariates together
covars<-merge(CPI, GAI, by='NAME', all=T)
covars<-merge(covars, HAQI, by='NAME',all=T)
covars<-merge(covars, GDP2, by='NAME',all=T)
covars<-merge(covars, edu_c, by='NAME',all=T)
covars<-merge(covars, GHSI21, by='NAME',all=T)
covars<-merge(covars, democ_all, by='NAME',all=T)


#get static version of data for countries to keep in our covariate dataset for correlations
data_static<-data[,.(NAME,avg_mpox,avg_monkeypox,mpox_greater, mpox_non0, avg_pct_mpox)]
data_static<-data_static[!duplicated(data_static),]

#historic outbreak locations in standardized fashion to match our other covariates
#From occurrence database 
historic_countries<-c('Central African Republic','Nigeria',"Ivory Coast",'South Sudan','Republic of Congo',
                      'Democratic Republic of the Congo','Gabon','Sierra Leone',
                      'Liberia','Cameroon','Ghana', 'United States','United Kingdom','Israel')

#recent outbreak data from OWID:
owid_mpox<-fread(paste0(dir,'raw data/owid-monkeypox-data.csv'))
#keep only those locations in our gst analys, but clean up naming first
owid_mpox[,NAME:=location]
#Standardize names
owid_mpox[location=="Congo", NAME:='Republic of Congo']
owid_mpox[location=="Democratic Republic of Congo", NAME:='Democratic Republic of the Congo']
owid_mpox[location=="Turkey", NAME:='Turkiye']
owid_mpox[location=="Czechia", NAME:='Czech Republic']
owid_mpox[location=="Congo", NAME:='Republic of Congo']

owid_mpox<-owid_mpox[NAME %in% data_static$NAME]
#consider some time-bound covarites - get the first date of an mpox case in a country
owid_mpox[total_cases>0,date_min:=as.Date(min(date, na.rm=T)), by='NAME']
owid_mpox_min<-owid_mpox[,.(NAME, date_min)]
owid_mpox_min<-owid_mpox_min[!duplicated(owid_mpox_min)]
owid_mpox_min<-owid_mpox_min[!is.na(owid_mpox_min$date_min)]
owid_mpox_min<-owid_mpox_min[order(date_min),]
#number of days since the start of the outbreak
owid_mpox_min<-owid_mpox_min[,day:=as.numeric(factor(date_min))]

#Figure out which ones ARENT included in the table above for 2022
dif<-setdiff(historic_countries,owid_mpox_min$NAME)
#create a new entry for these and append to above including some empty rows for 2022 outbreak data
historic_add<-data.frame(NAME=dif,
                         date_min=rep(as.Date('1999-01-01'),length(dif)),
                         day=rep(NA, length(dif)))

owid_mpox_min<-rbind(owid_mpox_min, historic_add)

#subset to only those countries we have GST data for
covars_sub<-covars[NAME %in% data_static$NAME]
#Add a category for not having had mpox
covars_sub<-merge(covars_sub, owid_mpox_min, by='NAME', all.x=T)
#If a country isnt part of the mpox database (old and 2022, we define ever having a case of mpox as NO (0)
covars_sub[,ever_mpox:=0]
covars_sub[NAME %in% owid_mpox_min$NAME, ever_mpox:=1]

#check missingness of all covariates
#remove day because the reason rows are missing is due to there never being an mpox case, not missing (more) at random like these other covariates
covars_temp<-copy(covars_sub)
covars_temp[,`:=`(date_min=NULL, day=NULL,
                  ISO3=NULL, Region=NULL, 
                  region=NULL)]

setnames(covars_temp, c("GAI","HAQI_mean","yr2022","Overall", "regime_row_owid"),
         c("LGBT GAI","HAQI","GDP_2022","GHSI Overall","Political regime"))
#missingness table
plot_pattern(covars_temp[!is.na(GDP_2022),],
             square = TRUE,
             rotate = TRUE
)

# 154 obs have all rows, another 12 missing GAI only remainder all missing more than 1 covariate
listwise_del_ctrys<-covars_sub[!is.na(GAI) & !is.na(regime_row_owid) & !is.na(Overall), NAME]
final_166_ctrys<-covars_sub[!(is.na(regime_row_owid) | is.na(Overall)), NAME]

#merge covariate data with mpox search data
full_data<-merge(data_static, covars_sub, by='NAME',all.x=T)

#Need to drop the countries that have very little data 
#use GDP as the bare-minimum for retention - would drop Guadeloupe, Martinique, Reunion and Taiwan
full_dat<-full_data[!is.na(yr2022)]

#clean up skewed and factor variables:
full_dat[,log_GDP:=log(yr2022)]

full_dat[,norm_fem_edu:=(Female_edu_mean_yrs_25_29-min(Female_edu_mean_yrs_25_29, na.rm=T))/(max(Female_edu_mean_yrs_25_29, na.rm=T)-min(Female_edu_mean_yrs_25_29, na.rm=T))]
full_dat[,logit_fem_edu:=log((norm_fem_edu)/(1-norm_fem_edu))]
full_dat[norm_fem_edu %in% c(0,1),logit_fem_edu:=log((norm_fem_edu+(0.5/184))/(1-norm_fem_edu+(0.5/184)))]

full_dat[,norm_gender_phone_gap:=(Risk_comms_gender_gap_access_phone_3_6_3a)/100]
full_dat[,logit_gender_phone_gap:=log((norm_gender_phone_gap)/(1-norm_gender_phone_gap))]
full_dat[norm_gender_phone_gap %in% c(0,1),logit_gender_phone_gap:=log((norm_gender_phone_gap+(0.5/184))/(1-norm_gender_phone_gap+(0.5/184)))]

full_dat[,norm_gender_internet_gap:=(Risk_comms_gender_gap_access_internet_3_6_4_a)/100]
full_dat[,logit_gender_internet_gap:=log((norm_gender_internet_gap+(0.5/184))/(1-norm_gender_internet_gap+(0.5/184)))]
full_dat[norm_gender_internet_gap %in% c(0,1),logit_gender_internet_gap:=log((norm_gender_internet_gap+(0.5/184))/(1-norm_gender_internet_gap+(0.5/184)))]

full_dat[,logit_lib:=log((lib_vdem_owid)/(1-lib_vdem_owid))]
full_dat[lib_vdem_owid %in% c(0,1),logit_lib:=log((lib_vdem_owid+(0.5/184))/(1-lib_vdem_owid+(0.5/184)))]

full_dat[,factor_pop_inclusion_riskcomm:=factor(Risk_comms_pop_inclusion_3_5_1b, levels=c(0,100),
                                                labels=c('No',"Yes"))]
full_dat[,factor_misinfo_riskcomm:=factor(Risk_comms_leader_share_misinfo_3_5_2b, levels=c(0,100),
                                                labels=c('No',"Yes"))]



rm(list=c('democ_all','edu','edu_c','edu2','GAI','GDP','GDP2','GHSI','GHSI21','HAQI','historic_add','owid_mpox',
          'covars','covars_sub','covars_temp','CPI','democ','owid_mpox_min', 'dif','filelist','x',
          'dat','data','data_static','small_dat', 'historic_countries','vect','data_old','data_foreign','ctry','fils','language'))

