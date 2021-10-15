#author Le GAO
#version 2021-10-15
#library package
library(dplyr)
library(lubridate)
library(readxl)
library(data.table)
library(fasttime)
library(ggplot2)
library(egg)
library(ggsci)
library(exactci)



# 0. generate pesudo index date1 and pesudo index date2----
#read combined DH data
cohort <- readRDS("E:/Cleaned/auto/4.cohort_full.RDS")
#if dod is "" then NA
library(data.table)
cohort[death_date_ymd=="", death_date_ymd := NA]
setDT(cohort)

# 1. Identify cohort (all active patients in HA during study period)
colnames(cohort) <- make.names(colnames(cohort))
library(dplyr)
cohort1 <- cohort %>% filter(!is.na(patient_pssn)) 
cohort_two_dif_vaccine <- cohort1 %>% 
  filter(!is.na(Vaccine.Brand.1st),!is.na(Vaccine.Brand.2nd),Vaccine.Brand.1st!=Vaccine.Brand.2nd)


cohort2 <- cohort1 %>% 
  filter(!patient_pssn%in%cohort_two_dif_vaccine$patient_pssn)
cohort2[Vaccine.Brand.1st=="Sinovac"&Vaccine.Brand.2nd=="Sinovac",vaccine2:="2Sino"]
cohort2[Vaccine.Brand.1st=="Sinovac"&is.na(Vaccine.Brand.2nd),vaccine2:="1Sino"]
cohort2[Vaccine.Brand.1st=="BioNTech/Fosun"&Vaccine.Brand.2nd=="BioNTech/Fosun",vaccine2:="2bion"]
cohort2[Vaccine.Brand.1st=="BioNTech/Fosun"&is.na(Vaccine.Brand.2nd),vaccine2:="1bion"]
cohort2[is.na(Vaccine.Brand.1st)&is.na(Vaccine.Brand.2nd),vaccine2:="none"]
table(cohort2$vaccine2)




# 2. Matching 
dist.match <- merge(cohort2[vaccinated==1, .N, keyby=.(Age, sex)], cohort2[vaccinated==0, .N, keyby=.(Age, sex)], by=c("Age","sex"), all.x=T)[, Ratio:=N.y/N.x]
set.seed(475)
seed <- .Random.seed


#function for generate pseudo id
gp.match <- function(cohort, case.var, match.var, seed) {
  setorderv(cohort, c(case.var, match.var, "PseudoID"))
  cases <- cohort[get(case.var)==1]
  ctrls <- cohort[get(case.var)==0]
  gps <- unique(cases[, match.var, with=F])
  .Random.seed <- seed
  
  res <- list()
  for(i in 1:nrow(gps)) {
    cat("\n", gps[i,Age], gps[i,sex])
    gp1 <- cases[Age==gps[i,Age] & sex==gps[i,sex]]
    gp0 <- ctrls[Age==gps[i,Age] & sex==gps[i,sex]]
    K = floor(nrow(gp0) / nrow(gp1)); cat(" -> 1 :",K)
    gp0_rand <- sample(gp0$PseudoID, nrow(gp0), replace=F)
    
    mat <- matrix(gp0_rand[1:(nrow(gp1)*K)], nrow=nrow(gp1), ncol=K)
    mat <- cbind(mat, c(gp0_rand[-c(1:(nrow(gp1)*K))], rep(NA, nrow(gp1)*(K+1)-nrow(gp0))))
    gpm <- cbind(gp1[, .(PseudoID)], mat)
    res[[i]] <- gpm
  }
  
  return(res)
}

matched <- gp.match(cohort2, "vaccinated", c("Age","sex"), seed)


#function for gen matched cohort
output_matches <- function(md, ct) {
  md_sizes <- unlist(lapply(md, nrow))
  ms <- rbindlist(lapply(1:length(md), function(i) melt(cbind(match=(sum(md_sizes[0:(i-1)]) + (1:nrow(md[[i]]))), md[[i]][, c("PseudoID", grep("^V[0-9]+", names(md[[i]]), value=T)), with=F]), id.vars=c("match"), value.name="PseudoID")))
  ct.m <- merge(ct, ms[, .(PseudoID, match, match_order=as.character(variable))], by="PseudoID")
  return(ct.m)
}

cohort.matched <- output_matches(matched, cohort2)


case_ind <- cohort.matched[match_order=="PseudoID", .(index.date=Date.of.vaccination.1st), keyby=match]
cohort.matched <- merge(cohort.matched, case_ind, by="match")
cohort.matched[, index.date := as.Date(index.date)]


median_vaccine <- cohort.matched %>% 
  filter(!is.na(Vaccine.Brand.1st),!is.na(Vaccine.Brand.2nd)) %>% 
  mutate(time_int=ymd(Date.of.vaccination.2nd)-ymd(Date.of.vaccination.1st)) %>% 
  group_by(Vaccine.Brand.1st) %>% 
  summarise(tt=median(time_int))
#bio 21d sino 28d


cohort.matched1  <- cohort.matched 
cohort.matched1[Vaccine.Brand.1st=="Sinovac"&is.na(Vaccine.Brand.2nd),index.date.2nd:=as.Date(Date.of.vaccination.1st+days(28))]
cohort.matched1[Vaccine.Brand.1st=="Sinovac"&!is.na(Vaccine.Brand.2nd),index.date.2nd:=as.Date(Date.of.vaccination.2nd)]
cohort.matched1[Vaccine.Brand.1st=="BioNTech/Fosun"&is.na(Vaccine.Brand.2nd),index.date.2nd:=as.Date(Date.of.vaccination.1st+days(21))]
cohort.matched1[Vaccine.Brand.1st=="BioNTech/Fosun"&!is.na(Vaccine.Brand.2nd),index.date.2nd:=as.Date(Date.of.vaccination.2nd)]
cohort.matched1[,index.date.2nd_1:=na.omit(index.date.2nd),match]  #na.omit remove na 
col_delete <- c("index.date.2nd","Age.1st","Age.2nd","Sex.1st","Sex.2nd")
cohort.matched1[,(col_delete) :=NULL][]

saveRDS(cohort.matched1,"E:/Cleaned/le_cohort_new/cohort_matched_0714.rds")





# 1. ip primary analysis 1st dose ----

# 1.1 generate data for incidence calculation ----
# read all vaccine data
vaccine <- readRDS("C:/Users/Team2/Documents/gl_dontedit/Le/data_new_0709/cohort_matched_0714.rds")
setkey(vaccine,NULL)
colnames(vaccine) <- make.names(colnames(vaccine))  #3946202
setDT(vaccine)
vaccine1 <- vaccine %>% 
  mutate(death_date_ymd=as.Date(fastPOSIXct(death_date_ymd))) %>% 
  mutate(Date.of.vaccination.2nd=as.Date(fastPOSIXct(Date.of.vaccination.2nd))) %>% 
  mutate(obs_start=if_else(is.na(Date.of.vaccination.1st),as.Date(fastPOSIXct(index.date)),as.Date(fastPOSIXct(Date.of.vaccination.1st)))) %>% 
  filter(is.na(death_date_ymd)|death_date_ymd>=obs_start) %>%    #exclude those death before our obs 
  mutate(obs_end28=pmin(obs_start+days(27),Date.of.vaccination.2nd-days(1),ymd('2021-06-30'),death_date_ymd,na.rm = T)) %>% 
  mutate(obs_end14=pmin(obs_start+days(13),Date.of.vaccination.2nd-days(1),ymd('2021-06-30'),death_date_ymd,na.rm = T)) #3638926


#save vaccine data out
vaccine3 <- vaccine1 %>% 
  select(-PseudoID,-patient_pssn,-pseudo_other_doc_no,-pseudo_hkid)


# read all dx data
all_dx <- readRDS("C:/Users/Team2/Documents/gl_dontedit/Le/data_new_0709/past_history_gopsopipaepx_ranking.RDS")
# table(all_dx$source)
setDT(all_dx)
all_dx[,admdt:=as.Date(fastPOSIXct(date_past_event))]

ip <- all_dx %>% 
  filter(source=="IP",ranking=="P")


# function to find each disease
dis_code <- list(c("Thrombocytopenia","^287.[345]|^279.12|^283.11|^284.1|^446.6|^776.1"),
                 c("Acute aseptic arthritis","^274.0|^696.0|^716.[569]|^712|^711.5"),
                 c("Acute disseminated encephalomyelitis","^323.6|^323.8"),
                 c("Anti-phospholipid antibody syndrome","^795.79|^286.5|^289.8|^D68.61"),
                 c("Guillain-Barré syndrome","^357.0|^357.8|^357.9"),
                 c("Kawasaki disease","^446.1"),
                 c("Narcolepsy","^347|^307.4|^780.5"),
                 c("Psoriatic arthritis all","^696"),
                 c("Reactive arthritis all","^099.3|^711.[13]|^716.[4569]|^719.4|^372.33"),
                 c("Rheumatoid arthritis all","^714"),
                 c("Single organ cutaneous vasculitis","^709.1|^446.2|^287.0"),
                 c("Sjogren syndrome","^710.2|^M35.0"),
                 c("Spondyloarthritis all","^720"),
                 c("Subacute thyroiditis","^245.1"),
                 c("Systemic lupus erythematosus","^710.0"),
                 c("Transverse myelitis","^341.2|^323.[04568]"))


vaccine_short <- vaccine1 %>% 
  select(patient_pssn,PseudoID,obs_start,obs_end28,obs_end14,death_date_ymd,Vaccine.Brand.1st,Vaccine.Brand.2nd)#3938942
vaccine_short[Vaccine.Brand.1st=="Sinovac"&Vaccine.Brand.2nd=="Sinovac",vaccine:="Sino_2"]
vaccine_short[Vaccine.Brand.1st=="Sinovac"&is.na(Vaccine.Brand.2nd),vaccine:="Sino_1"]
vaccine_short[Vaccine.Brand.1st=="BioNTech/Fosun"&Vaccine.Brand.2nd=="BioNTech/Fosun",vaccine:="Bion_2"]
vaccine_short[Vaccine.Brand.1st=="BioNTech/Fosun"&is.na(Vaccine.Brand.2nd),vaccine:="Bion_1"]
vaccine_short[is.na(Vaccine.Brand.1st)&is.na(Vaccine.Brand.2nd),vaccine:="None"]
table(vaccine_short$vaccine)


ip_vaccine <- merge(vaccine_short,ip,by="patient_pssn",all.x=T) %>% 
  filter(is.na(death_date_ymd)|death_date_ymd>=admdt) 
dx_vaccine <- merge(vaccine_short,all_dx,by="patient_pssn",all.x=T)%>% 
  filter(is.na(death_date_ymd)|death_date_ymd>=admdt) 


gl_func <- function(z){
  x <- z[1]
  y <- z[2]
  
  
  #find ip cases of interest disease
  dx_y <- ip_vaccine %>% 
    filter(grepl(y,codes)==T)
  
  dx_y_vaccine28 <- dx_y %>% 
    mutate(case_28=if_else(admdt>=obs_start&admdt<=obs_end28,1,0)) %>% 
    group_by(patient_pssn)%>% 
    arrange(admdt) %>% 
    mutate(case_1st_28=if_else(case_28==1&duplicated(patient_pssn)==F,1,0)) %>% 
    arrange(-case_28,-case_1st_28,admdt) %>% 
    filter(duplicated(patient_pssn)==F) %>%
    ungroup() %>% 
    select(patient_pssn,case_28,case_1st_28,admdt,PseudoID) %>% 
    `colnames<-`(c("patient_pssn","case_28","case_1st_28","admdt28","PseudoID"))
  
  dx_y_vaccine14 <- dx_y %>% 
    mutate(case_14=if_else(admdt>=obs_start&admdt<=obs_end14,1,0)) %>% 
    group_by(patient_pssn)%>% 
    arrange(admdt) %>% 
    mutate(case_1st_14=if_else(case_14==1&duplicated(patient_pssn)==F,1,0)) %>% 
    arrange(-case_14,-case_1st_14,admdt) %>% 
    filter(duplicated(patient_pssn)==F) %>%
    ungroup() %>% 
    select(patient_pssn,case_14,case_1st_14,admdt,PseudoID) %>% 
    `colnames<-`(c("patient_pssn","case_14","case_1st_14","admdt14","PseudoID"))
  
  dx_y_vaccine <- merge(dx_y_vaccine28,dx_y_vaccine14,by=c("patient_pssn","PseudoID"),all=T)
  
  
  #find history in past three years 
  his_y <- dx_vaccine %>% 
    filter(grepl(y,codes)==T)
  
  dx_y_his <- his_y %>% 
    mutate(case_his=if_else(admdt>=obs_start-years(1)&admdt<obs_start,1,0)) %>% 
    group_by(patient_pssn)%>% 
    arrange(-case_his,admdt) %>% 
    filter(duplicated(patient_pssn)==F) %>%
    ungroup() %>% 
    select(patient_pssn,case_his,admdt,PseudoID) %>% 
    `colnames<-` (c("patient_pssn","case_his","admdt_his","PseudoID"))
  
  dx_y_vaccine1 <- merge(dx_y_vaccine,dx_y_his,by=c("patient_pssn","PseudoID"),all=T) %>% 
    mutate(disease=x)
  
  return(dx_y_vaccine1)
}


all_case <- bind_rows(lapply(dis_code,function(x) {print(x);gl_func(x)}))

all_case2 <- all_case
all_case2$id <- openssl::sha256(paste0(all_case2$PseudoID, "----"))

all_case3 <- all_case2 %>% 
  select(-patient_pssn,-PseudoID)

#median day to event without history, first IP happened after first vaccine or first pseudo index date
fun_median <- function(z){
  x <- z[1]
  y <- z[2]
  #find ip cases of interest disease
  dx_y_vaccine_new <- ip_vaccine %>% 
    filter(grepl(y,codes)==T) %>% 
    mutate(case_new=if_else(admdt>=obs_start,1,0)) %>% 
    group_by(patient_pssn)%>% 
    arrange(-case_new,admdt) %>% 
    filter(duplicated(patient_pssn)==F) %>%
    ungroup() %>% 
    select(patient_pssn,case_new,admdt) %>% 
    `colnames<-`(c("patient_pssn","case_new","admdt_new"))
  
  #find history in past three years 
    dx_y_his <- dx_vaccine %>% 
    filter(grepl(y,codes)==T) %>% 
    mutate(case_his=if_else(admdt>=obs_start-years(1)&admdt<obs_start,1,0)) %>% 
    group_by(patient_pssn)%>% 
    arrange(-case_his,admdt) %>% 
    filter(duplicated(patient_pssn)==F) %>%
    ungroup() %>% 
    select(patient_pssn,case_his,admdt) %>% 
    `colnames<-` (c("patient_pssn","case_his","admdt_his"))
  
  dx_y_vaccine1 <- merge(dx_y_vaccine_new,dx_y_his,by="patient_pssn",all = T)
  
  dx_y_vaccine2 <- merge(vaccine_short,dx_y_vaccine1,by="patient_pssn",all.x=T) %>% 
    mutate(case_new=if_else(is.na(case_new),0,case_new)) %>% 
    mutate(case_his=if_else(is.na(case_his),0,case_his)) %>% 
    filter(case_his==0,case_new==1) %>% 
    mutate(fud=admdt_new-obs_start+1) %>% 
    group_by(vaccine) %>% 
    mutate(day_to_case_min=min(fud),day_to_case_max=max(fud)) %>% 
    mutate(day_to_case_25=quantile(fud,0.25),day_to_case_75=quantile(fud,0.75)) %>% 
    mutate(day_to_case_median=median(fud),case_no=sum(case_new)) %>% 
    filter(duplicated(vaccine)==F) %>% 
    select(vaccine,day_to_case_min,day_to_case_max,day_to_case_25,day_to_case_75,day_to_case_median,case_no) %>% 
    mutate(disease=1)
  
  return(dx_y_vaccine2)
  
}
median_case <- bind_rows(lapply(dis_code,function(x) {print(x);fun_median(x)}))



# 1.2 calculate the incidence rate and cumulative incidence---- 
#ref pop no vaccine group
#cut age group
vaccine_new <- vaccine3 %>% 
  mutate(vaccine=if_else(is.na(Vaccine.Brand.1st),"None",Vaccine.Brand.1st)) %>% 
  mutate(age_grp=cut(as.numeric(Age),breaks=c(0,16,20,25,30,35,40,45,50,55,60,65,70,75,80,85,200),right=F,
                     labels = c("0-15","16-19","20-24","25-29","30-34","35-39","40-44","45-49",
                                "50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+"))) %>% 
  filter(age_grp!="0-15")

#calculate the ref pop-- non vaccine as ref group
ref_non <- vaccine_new %>% 
  filter(vaccine=="None") %>% 
  mutate(py_all=as.numeric(sum((obs_end28-obs_start+1)/365.25))) %>% 
  group_by(age_grp) %>% 
  mutate(py_ref=as.numeric(sum((obs_end28-obs_start+1)/365.25))) %>% 
  mutate(py_pro=py_ref/py_all) %>% 
  filter(duplicated(age_grp)==F) %>% 
  select(age_grp,py_pro)

ref_non_sex <- vaccine_new %>% 
  filter(vaccine=="None") %>% 
  group_by(sex) %>% 
  mutate(py_all=as.numeric(sum((obs_end28-obs_start+1)/365.25))) %>%
  ungroup() %>% 
  group_by(age_grp,sex) %>% 
  mutate(py_ref=as.numeric(sum((obs_end28-obs_start+1)/365.25))) %>% 
  mutate(py_pro=py_ref/py_all) %>% 
  filter(duplicated(age_grp)==F) %>% 
  select(age_grp,py_ref,sex,py_pro)


#HK 2021 ref pop
hk_2021 <- read_xlsx("D:/data_no_share_20200506/9.covid_auto_immune/Projected_HKRP_2021.xlsx") %>% 
  mutate(age_grp=cut(as.numeric(Age),breaks=c(0,16,20,25,30,35,40,45,50,55,60,65,70,75,80,85,200),right=F,
                     labels = c("0-15","16-19","20-24","25-29","30-34","35-39","40-44","45-49",
                                "50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")))
hk_2021 <- read_xlsx("C:/Users/Team2/Documents/gl_dontedit/Le/Projected_HKRP_2021.xlsx") %>% 
  mutate(age_grp=cut(as.numeric(Age),breaks=c(0,16,20,25,30,35,40,45,50,55,60,65,70,75,80,85,200),right=F,
                     labels = c("0-15","16-19","20-24","25-29","30-34","35-39","40-44","45-49",
                                "50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")))
ref_pop <-hk_2021 %>% 
  group_by(age_grp) %>% 
  mutate(ref_pop_sum=sum(POP)) %>% 
  filter(duplicated(age_grp)==F) %>%
  ungroup() %>% 
  mutate(ref_all_pro=ref_pop_sum/sum(ref_pop_sum)) %>% 
  select(age_grp,ref_all_pro) %>% 
  filter(age_grp!="0-15")

ref_pop_sex <- hk_2021 %>% 
  group_by(age_grp,SEX) %>% 
  mutate(pop_sum=sum(POP)) %>% 
  filter(duplicated(age_grp)==F) %>%
  ungroup() %>% 
  group_by(SEX) %>%  
  mutate(ref_sex_pro=pop_sum/sum(pop_sum)) %>% 
  select(age_grp,SEX,ref_sex_pro,pop_sum) %>% 
  filter(age_grp!="0-15")

ref_all <- merge(ref_non,ref_pop,by="age_grp")


#read disease data
#find death records
all_case_death <- vaccine3 %>% 
  filter(!is.na(death_date_ymd)) %>% 
  mutate(case_28=if_else(death_date_ymd>=obs_start&death_date_ymd<=obs_end28,1,0)) %>% 
  mutate(case_14=if_else(death_date_ymd>=obs_start&death_date_ymd<=obs_end14,1,0)) %>% 
  mutate(case_1st_28=as.numeric(NA),case_1st_14=as.numeric(NA),case_his=as.numeric(0),disease="Death") %>% 
  mutate(admdt28=if_else(case_28==1,death_date_ymd,as.Date(NA)),admdt14=if_else(case_14==1,death_date_ymd,as.Date(NA)),admdt_his=as.Date(NA)) %>% 
  select(colnames(all_case3))
colnames(all_case_death) <- colnames(all_case3)


all_case4 <- rbind(all_case3,all_case_death) %>% 
  mutate(admdt28=if_else(case_28==1,admdt28,as.Date(NA)))

dis_code1 <- list(c("Thrombocytopenia","^287.[345]|^279.12|^283.11|^284.1|^446.6|^776.1"),
                  c("Acute aseptic arthritis","^274.0|^696.0|^716.[569]|^712|^711.5"),
                  c("Acute disseminated encephalomyelitis","^323.6|^323.8"),
                  c("Anti-phospholipid antibody syndrome","^795.79|^286.5|^289.8|^D68.61"),
                  c("Death","Death date"),
                  c("Guillain-Barré syndrome","^357.0|^357.8|^357.9"),
                  c("Kawasaki disease","^446.1"),
                  c("Narcolepsy","^347|^307.4|^780.5"),
                  c("Psoriatic arthritis all","^696"),
                  c("Reactive arthritis all","^099.3|^711.[13]|^716.[4569]|^719.4|^372.33"),
                  c("Rheumatoid arthritis all","^714"),
                  c("Single organ cutaneous vasculitis","^709.1|^446.2|^287.0"),
                  c("Sjogren syndrome","^710.2|^M35.0"),
                  c("Spondyloarthritis all","^720"),
                  c("Subacute thyroiditis","^245.1"),
                  c("Systemic lupus erythematosus","^710.0"),
                  c("Transverse myelitis","^341.2|^323.0|^323.4|^323.5|^323.6|^323.8"))


vaccine_new1 <- vaccine_new %>% 
  select(id,age_grp,vaccine,sex,obs_end28,obs_end14,obs_start)

fun_inc <- function(z){
  x <- z[1]
  y <- z[2]
  
  x_case <- all_case4 %>% 
    filter(disease==x)
  
  x_case_vaccine <- bind_rows(merge(vaccine_new %>% filter(id %in% x_case$id),x_case,by="id",all.x=T),
                              vaccine_new %>% filter(!id %in% x_case$id)) %>% 
    mutate(case_his=if_else(is.na(case_his),0,case_his)) %>% 
    filter(case_his==0)
  
  
  x_case_denominator <- x_case_vaccine %>% 
    group_by(age_grp,vaccine,sex) %>% 
    mutate(ct_pop_age=length(id),
           py_deno28=as.numeric(sum(pmin(obs_end28,admdt28,na.rm = T)-obs_start+1)/365.25)) %>% 
    filter(duplicated(age_grp)==F) %>% 
    select(age_grp,vaccine,ct_pop_age,sex,py_deno28) %>% 
    arrange(age_grp,vaccine,sex)
  
  x_count <- x_case_vaccine %>% 
    group_by(age_grp,vaccine,sex) %>% 
    mutate(c1=sum(case_28,na.rm = T),c2=sum(case_1st_28,na.rm = T)) %>% 
    filter(duplicated(age_grp)==F) %>% 
    ungroup() %>% 
    select(c1,c2,vaccine,age_grp,sex)
  
  x_ir_byage <- merge(x_count,x_case_denominator,by=c("age_grp","sex","vaccine"),all = T) %>% 
    mutate(c1=ifelse(is.na(c1),0,c1),c2=ifelse(is.na(c2),0,c2)) %>% 
    mutate(disease=x,dis_code=y)
  
  return(x_ir_byage)
  
}

all_ir_cum_age <- bind_rows(lapply(dis_code1,function(x) {print(x);fun_inc(x)}))

#age standardised incidence all sex
all_ir_age_all <- all_ir_cum_age %>% 
  group_by(age_grp,disease,vaccine) %>% 
  mutate(c_all1=sum(c1),c_all2=sum(c2),# count case by age group
         py_deno28_all=sum(py_deno28), # sum the total observed py
         ir_all1=(c_all1/py_deno28_all)*100000,ir_all2=(c_all2/py_deno28_all)*100000,# calculate the crude ir
         pop_all_1=sum(ct_pop_age), # sum the total observed pop
         cuminc_all1=(c_all1/pop_all_1)*100000,cuminc_all2=(c_all2/pop_all_1)*100000) %>% # calculate the crude cum inc 
  filter(duplicated(age_grp)==F) %>% 
  ungroup() %>% 
  select(-sex,-c1,-c2)

std_ir_all <- merge(all_ir_age_all, ref_all, by.x=c("age_grp"),by.y=c("age_grp")) %>% 
  mutate(ir_sum1=ir_all1*py_pro,ir_sum2=ir_all2*py_pro) %>% 
  mutate(cum_inc_sum1=cuminc_all1*ref_all_pro,cum_inc_sum2=cuminc_all2*ref_all_pro) %>% 
  group_by(vaccine,disease) %>% 
  mutate(std_cum_inc1=sum(cum_inc_sum1),std_inc2=sum(cum_inc_sum2),
         ct_case1=sum(c_all1),ct_case2=sum(c_all2),
         std_ir1=sum(ir_sum1),std_ir2=sum(ir_sum2),
         obs_py_28_all=sum(py_deno28_all),
         crud_ir1=ct_case1/obs_py_28_all*100000,crud_ir2=ct_case2/obs_py_28_all*100000,
         obs_pop_all=sum(pop_all_1),
         crud_cum_inc1=ct_case1/obs_pop_all*100000,crud_cum_inc2=ct_case2/obs_pop_all*100000) %>% 
  filter(duplicated(vaccine)==F) %>% 
  select(vaccine,disease,std_cum_inc1,ct_case1,std_ir1,crud_ir1,crud_cum_inc1,obs_pop_all,obs_py_28_all,dis_code) %>% 
  ungroup() %>% 
  arrange(disease,vaccine) %>% 
  mutate(exp_ct1=std_ir1*obs_py_28_all/100000)


#function to calculate the median day to event
fun_median <- function(z){
  x <- z[1]
  y <- z[2]
  
  x_case <- all_case4 %>% 
    mutate(case_his=if_else(is.na(case_his),0,case_his)) %>% 
    filter(disease==x,case_28==1,case_his==0)
  
  x_case_vaccine <- merge(vaccine_new1,x_case,by="id",all.y = T)
  
  
  x_case_median <- x_case_vaccine %>% 
    mutate(day_to_event=admdt28-obs_start+1) %>% 
    select(id,day_to_event,vaccine) %>% 
    mutate(disease=x)
  
  return(x_case_median)
  
}
all_median <- bind_rows(lapply(dis_code1,function(x) {print(x);fun_median(x)}))

median_sum <- all_median %>% 
  group_by(disease,vaccine) %>% 
  mutate(day_to_event_median=median(day_to_event)) %>% 
  mutate(day_to_event_25=quantile(day_to_event,0.25)) %>% 
  mutate(day_to_event_75=quantile(day_to_event,0.75)) %>% 
  mutate(day_to_event_min=min(day_to_event)) %>% 
  mutate(day_to_event_max=max(day_to_event)) %>% 
  filter(duplicated(disease)==F)

#function to calculate risk ratio
rr_fun <- function(z){
  x <- z[1]
  y <- z[2]
  
  c1 <- std_ir_all %>% 
    filter(disease==x,vaccine=="BioNTech/Fosun") %>% 
    mutate(col1=round(exp_ct1),col2=obs_py_28_all) %>% 
    select(col1,col2)
  c2 <- std_ir_all %>% 
    filter(disease==x,vaccine=="Sinovac") %>% 
    mutate(col1=round(exp_ct1),col2=obs_py_28_all) %>% 
    select(col1,col2)
  c3 <- std_ir_all %>% 
    filter(disease==x,vaccine=="None") %>% 
    mutate(col1=round(exp_ct1),col2=obs_py_28_all) %>% 
    select(col1,col2)
  
  # method 1
  out1 <- poisson.test(c(c1$col1,c3$col1),c(c1$col2,c3$col2))
  out.dat1 <- data.frame("IRR_1"=out1$estimate,"LCI_1"=out1$conf.int[1],"UCI_1"=out1$conf.int[2],"P_value_1"=out1$p.value) %>% 
    mutate(vaccine="BioNTech/Fosun",disease=x)
  out2 <- poisson.test(c(c2$col1,c3$col1),c(c2$col2,c3$col2))
  out.dat2 <- data.frame("IRR_1"=out2$estimate,"LCI_1"=out2$conf.int[1],"UCI_1"=out2$conf.int[2],"P_value_1"=out2$p.value) %>% 
    mutate(vaccine="Sinovac",disease=x)
  
  out.data.poisson.test <- rbind(out.dat1,out.dat2)
  
  # method 2
  eout1 <- poisson.exact(c(c1$col1,c3$col1),T=c(c1$col2,c3$col2),
                         r=2,alternative="two.sided",
                         tsmethod = "central",conf.level = 0.95,plot = F,midp = T)
  eout.dat1 <- data.frame("IRR_2"=eout1$estimate,"LCI_2"=eout1$conf.int[1],"UCI_2"=eout1$conf.int[2],"P_value_2"=eout1$p.value) %>% 
    mutate(vaccine="BioNTech/Fosun",disease=x)
  eout2 <- poisson.exact(c(c2$col1,c3$col1),T=c(c2$col2,c3$col2),
                         r=2,alternative="two.sided",
                         tsmethod = "central",conf.level = 0.95,plot = F,midp = T)
  eout.dat2 <- data.frame("IRR_2"=eout2$estimate,"LCI_2"=eout2$conf.int[1],"UCI_2"=eout2$conf.int[2],"P_value_2"=eout2$p.value) %>% 
    mutate(vaccine="Sinovac",disease=x)
  
  out.data.poisson.exact <- rbind(eout.dat1,eout.dat2)
  
  out_all <- merge(out.data.poisson.test,out.data.poisson.exact,by=c("disease","vaccine"))
  
  return(out_all)
}

out_rr <- bind_rows(lapply(dis_code1,rr_fun))



# 2. ip primary analysis 2nd dose----
# same as the 1st dose



# 3. ecology analysis
#read 3946550 HA cohort data 
ha_active <- readRDS("C:/Users/Team2/Documents/gl_dontedit/Le/data_new_0709/4.cohort_full.RDS")
haskey(ha_active)
setkey(ha_active)
colnames(ha_active) <- make.names(colnames(ha_active))

#if dod is "" then NA
ha_active[death_date_ymd=="", death_date_ymd := NA]
setDT(ha_active)

ha_active1 <- ha_active %>% filter(!is.na(patient_pssn))

ha_active2 <- ha_active1
ha_active2$id <- openssl::sha256(paste0(ha_active2$PseudoID, "----"))
ha_active3 <- ha_active2 %>% 
  select(-patient_pssn,-PseudoID,-pseudo_hkid,-pseudo_other_doc_no)

saveRDS(ha_active3,"C:/Users/Team2/Documents/gl_dontedit/Le/data_new_0709/le_ana_0709/ha_active_3946550.rds")


# read all dx data
all_dx <- readRDS("C:/Users/Team2/Documents/gl_dontedit/Le/data_new_0709/past_history_gopsopipaepx_ranking.RDS")
# table(all_dx$source)
setDT(all_dx)
all_dx[,admdt:=as.Date(fastPOSIXct(date_past_event))]

ip <- all_dx %>% 
  filter(source=="IP",ranking=="P")

ip_vaccine <- merge(ha_active1,ip,by="patient_pssn",all.x=T) %>% 
  filter(is.na(death_date_ymd)|death_date_ymd>=admdt) 


dis_code <- list(c("(Idiopathic) Thrombocytopenia","^287.[345]|^279.12|^283.11|^284.1|^446.6|^776.1"),
                 c("Acute aseptic arthritis","^274.0|^696.0|^716.[569]|^712|^711.5"),
                 c("Acute disseminated encephalomyelitis","^323.6|^323.8"),
                 c("Anti-phospholipid antibody syndrome","^795.79|^286.5|^289.8|^D68.61"),
                 c("Guillain-Barr? syndrome","^357.0|^357.8|^357.9"),
                 c("Kawasaki disease","^446.1"),
                 c("Narcolepsy all","^347|^89.1[78]|^307.4|^780.5"),
                 c("Psoriatic arthritis all","^696"),
                 c("Reactive arthritis all","^099.3|^711.[13]|^716.[4569]|^719.4|^372.33"),
                 c("Rheumatoid arthritis all","^714"),
                 c("Single organ cutaneous vasculitis","^709.1|^446.2|^287.0"),
                 c("Sjogren syndrome","^710.2|^M35.0"),
                 c("Spondyloarthritis all","^720"),
                 c("Subacute thyroiditis","^245.1"),
                 c("Systemic lupus erythematosus","^710.0"),
                 c("Transverse myelitis","^341.2|^323.[04568]"))


#HK 2021 ref pop
#D:/data_no_share_20200506/9.covid_auto_immune
hk_2021 <- read_xlsx("C:/Users/Team2/Documents/gl_dontedit/Le/Projected_HKRP_2021.xlsx") %>% 
  mutate(age_grp=cut(as.numeric(Age),breaks=c(0,16,20,25,30,35,40,45,50,55,60,65,70,75,80,85,200),right=F,
                     labels = c("0-15","16-19","20-24","25-29","30-34","35-39","40-44","45-49",
                                "50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")))
ref_pop <-hk_2021 %>% 
  group_by(age_grp) %>% 
  mutate(ref_pop_sum=sum(POP)) %>% 
  filter(duplicated(age_grp)==F) %>%
  ungroup() %>% 
  mutate(ref_all_pro=ref_pop_sum/sum(ref_pop_sum)) %>% 
  select(age_grp,ref_all_pro) %>% 
  filter(age_grp!="0-15")



#historical case

gl_func_his <- function(z,t){
  x <- z[1]
  y <- z[2]
  
  dx_y_case <- ip_vaccine %>% 
    filter(grepl(y,codes)==T) %>%
    mutate(his=if_else(admdt>=ymd(paste0(t,'-03-01'))&admdt<=ymd(paste0(t,'-06-30')),1,0)) %>% 
    arrange(patient_pssn,-his,admdt) %>% 
    filter(duplicated(patient_pssn)==F) %>% 
    mutate(age=t-as.numeric(dob_y)) %>% 
    mutate(age_grp=cut(as.numeric(age),breaks=c(0,16,20,25,30,35,40,45,50,55,60,65,70,75,80,85,200),right=F,
                       labels = c("0-15","16-19","20-24","25-29","30-34","35-39","40-44","45-49",
                                  "50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+"))) %>% 
    filter(age_grp!="0-15") %>% 
    group_by(age_grp) %>% 
    mutate(his_ct=sum(his,na.rm = T)) %>% 
    select(age_grp,his_ct) %>% 
    filter(duplicated(age_grp)==F) %>% 
    ungroup()
  
  dx_y_denominator <- ha_active1 %>% 
    filter(death_date_ymd>=ymd(paste0(t,'-03-01'))|is.na(death_date_ymd)) %>% 
    mutate(age=t-as.numeric(dob_y)) %>% 
    mutate(age_grp=cut(as.numeric(age),breaks=c(0,16,20,25,30,35,40,45,50,55,60,65,70,75,80,85,200),right=F,
                       labels = c("0-15","16-19","20-24","25-29","30-34","35-39","40-44","45-49",
                                  "50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+"))) %>% 
    filter(age_grp!="0-15") %>%
    group_by(age_grp) %>% 
    mutate(deno_ct=length(patient_pssn)) %>% 
    select(age_grp,deno_ct) %>% 
    filter(duplicated(age_grp)==F) %>% 
    ungroup()
  
  dx_y_crude <- merge(dx_y_case,dx_y_denominator,by="age_grp",all = T) %>% 
    mutate(his_ct=ifelse(is.na(his_ct),0,his_ct),
           deno_ct=ifelse(is.na(deno_ct),0,deno_ct),
           c_inc=his_ct/deno_ct*100000)
  
  dx_y_standard <- merge(dx_y_crude,ref_pop,by="age_grp") %>% 
    mutate(c_inc_1=c_inc*ref_all_pro,
           std_cum_inc_his=sum(c_inc_1),
           sum_ct_his=sum(his_ct),
           sum_denominator=sum(deno_ct)) %>% 
    filter(duplicated(sum_denominator)==F) %>% 
    select(std_cum_inc_his,sum_ct_his,sum_denominator) %>% 
    mutate(year=t,disease=x)
  return(dx_y_standard)
}


hist_case2018 <- bind_rows(lapply(dis_code,function(x) {message(x);gl_func_his(x,2018)}))
hist_case2019 <- bind_rows(lapply(dis_code,function(x) {message(x);gl_func_his(x,2019)}))
hist_case2020 <- bind_rows(lapply(dis_code,function(x) {message(x);gl_func_his(x,2020)}))
hist_case2021 <- bind_rows(lapply(dis_code,function(x) {message(x);gl_func_his(x,2021)}))

hist_case <- rbind(hist_case2018,hist_case2019,hist_case2020,hist_case2021)

hist_case_ci <- hist_case %>% 
  mutate(exp_case=round(std_cum_inc_his*sum_denominator/100000),
         lci=NA,
         uci=NA)
for(i in 1: nrow(hist_case_ci)){
  message(i,"/",nrow(hist_case_ci))
  hist_case_ci[i,"lci"] <- binom.test(hist_case_ci$exp_case[i],hist_case_ci$sum_denominator[i])$conf.int[1]*100000
  hist_case_ci[i,"uci"] <- binom.test(hist_case_ci$exp_case[i],hist_case_ci$sum_denominator[i])$conf.int[2]*100000
}



#plot 
setDT(hist_case_ci)

hist_case_ci1 <- hist_case_ci %>% 
  filter(disease %in% c("(Idiopathic) Thrombocytopenia","Acute aseptic arthritis",
                        "Acute disseminated encephalomyelitis","Anti-phospholipid antibody syndrome",
                        "Guillain-Barré syndrome","Kawasaki disease",
                        "Narcolepsy all","Psoriatic arthritis all",
                        "Reactive arthritis all","Rheumatoid arthritis all",
                        "Single organ cutaneous vasculitis",
                        "Sjogren syndrome","Spondyloarthritis all",
                        "Subacute thyroiditis","Systemic lupus erythematosus",
                        "Transverse myelitis")) %>% 
  mutate(`Body system`=as.character(NA))
hist_case_ci1[disease=="(Idiopathic) Thrombocytopenia",disease:="Thrombocytopenia"]
hist_case_ci1[disease=="Narcolepsy all",disease:="Narcolepsy and related disorders"]
hist_case_ci1[disease=="Psoriatic arthritis all",disease:="Psoriatic arthritis"]
hist_case_ci1[disease=="Sjogren syndrome",disease:="Sjogren's syndrome"]
hist_case_ci1[disease=="Spondyloarthritis all",disease:="Spondyloarthritis"]
hist_case_ci1[disease=="Reactive arthritis all",disease:="Reactive arthritis"]
hist_case_ci1[disease=="Rheumatoid arthritis all",disease:="Rheumatoid arthritis"]
colnames(hist_case_ci1)[7] <- "Disease"

plot_code <- list(c("^Kawasaki disease|^Single organ cutaneous vasculitis","Cardiovascular system"),
                  c("^Subacute thyroiditis","Endocrine system"),
                  c("^Anti-phospholipid antibody syndrome|^Thrombocytopenia","Haematological system"),
                  c("^Sjogren's syndrome|^Systemic lupus erythematosus","Multisystem"),
                  c("^Acute aseptic arthritis|^Reactive arthritis|^Rheumatoid arthritis|^Psoriatic arthritis|^Spondyloarthritis","Musculoskeletal system"),
                  c("^Acute disseminated encephalomyelitis|^Guillain-Barré syndrome|^Narcolepsy and related disorders|^Transverse myelitis","Nervous system"))

hist_case_ci1 <- hist_case_ci1 %>% 
  mutate(year=as.factor(year))# if need to change x label name, then need to be factor

plt_fun <- function(z) {
  x <- z[1]
  y <- z[2]
  
  plot1 <- ggplot(data=hist_case_ci1 %>% filter(grepl(x,Disease)==T), 
                  aes(x=year, y=std_cum_inc_his,group=Disease, color=Disease)) +
    geom_line(size = 1.2)+
    geom_errorbar(aes(ymin=lci,ymax=uci),width=.05)+
    geom_point(shape=18,size=5)+
    theme_classic()+
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      axis.text.x = element_text(size = 15,face = "bold"),
      axis.text.y = element_text(size=15,face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size =16,face = "bold"),
      plot.title =element_text(size=24,face = "bold")
    )+
    scale_color_manual(values = pal_nejm()(5))+
    scale_x_discrete(breaks=c("2018","2019","2020","2021"),
                     labels=c("2018 Mar to Jun", "2019 Mar to Jun", "2020 Mar to Jun","2021 Mar to Jun"))+
    scale_y_continuous(limits=c(0,max((hist_case_ci1 %>% filter(grepl(x,Disease)==T) %>% pull(uci))*4/3)))+
    labs(y = "Cunulative incidence /n (per 100,000 persons)")+
    ggtitle(y)+
    theme(legend.position="bottom")+
    guides(color=guide_legend(nrow=2,byrow=TRUE))
  return(plot1)
  
}

plot_end <- lapply(plot_code,plt_fun)

a <- ggarrange(plots =plot_end)
ggsave("D:/OneDrive - connect.hku.hk/Other-share/5.Phd project/20.care_autoimmue/data_code_0722_his_plot/his_1821_ipp.jpeg",a, units="in", width=16, height=16, dpi=600)


