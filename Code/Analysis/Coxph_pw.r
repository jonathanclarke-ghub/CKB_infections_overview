### Libraries ####

library(DT)
library(tidyverse)
library(tidytext)
library(reshape) 
library(dplyr)
library(survival)
library(stringr)
library(tidyr)
library(ggplot2)
library(ckbplotr)
library(plotrix)
library(gridExtra)
library(grid)

### Set paths to data and results ####

path_1<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Raw_data_7/"
path_2<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Project_data_7/"
path_3<-"K:/kadoorie/Staff_Folders/JonathanC/Projects/Serology/Tables_figures/Tables_figures_8/"

### Import participant data ####

data_all<-readRDS(paste0(path_2,"data_serology_2024_7.rds"))

table(data_all$multiplexserology_cancer_b1_subcohort,data_all$cancer_diag,useNA="ifany")

### Pathogen/antigen seroprevalence ####

seropositivity_fun<-function(pathogen_antigen_pos,code,output){  
  
  tab<-list()
  
  data<-data_all %>% filter(multiplexserology_cancer_b1_subcohort==1)
  
  results<-list()
  
  for(ii in pathogen_antigen_pos){
    
    # Calculate summary statistics
    
    n <- as.numeric(length(!is.na(data[[ii]])))
    prop<-mean(data[[ii]],rm.na=TRUE)
    se<-std.error(data[[ii]])
    
    # Calculate margin of error - two-sided 95% CI
    
    margin <- qnorm(0.975)*sqrt(prop*(1-prop)/n)
    
    # Calculate lower and upper bounds of confidence interval
    
    ci_lower <- (prop - margin)*100
    ci_upper <- (prop + margin)*100
    prop     <- prop*100
    se<-se*100
    pathogen<-ii
    
    results[[ii]]<-c(pathogen,n,prop,se,ci_lower,ci_upper)
    
    
  }
  
  tab<-as.data.frame(do.call(rbind,results))
  names(tab)<-c("pathogen","n","prop","stderr","ci_lower","ci_upper")
  
  ### Format table ####
  
  tab$ci_low<-sprintf(as.numeric(tab$ci_lower), fmt = '%#.1f')  
  tab$ci_up<-sprintf(as.numeric(tab$ci_upper), fmt = '%#.1f') 
  tab$prev<-sprintf(as.numeric(tab$prop), fmt = '%#.1f') 
  tab$se<-sprintf(as.numeric(tab$stderr), fmt = '%#.2f') 
  tab$prev_se<-paste0(tab$prev," (",tab$se,")") 
  tab$prev_ci<-paste0(tab$prev," (",tab$ci_low," to ",tab$ci_up,")") 
  tab$pathogen<-str_remove_all(tab$pathogen,code)
  
  saveRDS(tab,paste0(path_2,output,".rds"))
  
}

seropositivity_fun(grep("_pp$",names(data_all),value=TRUE),"_pp","tab_pathogen_pos")
seropositivity_fun(grep("_ap$",names(data_all),value=TRUE),"_ap","tab_antigen_pos")

### Format meta data ####

tab_meta<-readRDS(paste0(path_1,"multiplex_cancer_case_cohort_meta.rds"))
tab_meta$pathogen_label_old<-tab_meta$pathogen_label
tab_meta$pathogen_label[36:47]<-""
tab_meta$pathogen_label[36]<-"*C. trachomatis*"
tab_meta$pathogen_label[37:38]<-"*T. gondii*"
tab_meta$pathogen_label[39]<-"*C. burnetii*"
tab_meta$pathogen_label[40:47]<-"*H. pylori*"
tab_meta$antigen_label_long<-tab_meta$antigen_label
tab_meta$antigen_label_long[1:35]<-paste0(tab_meta$pathogen_label," ",tab_meta$antigen_label)[1:35]
tab_meta_pathogen_labels<-tab_meta[!duplicated(tab_meta$pathogen),c("pathogen","pathogen_label","antigen_label_long")]
tab_meta_pathogen_labels$label<-tab_meta_pathogen_labels$pathogen_label
tab_meta_antigen_labels<-tab_meta[!duplicated(tab_meta$variable_name),c("antigen_label","pathogen_label","antigen_label_long","variable_name")]   
tab_meta_antigen_labels<-tab_meta_antigen_labels %>% dplyr::rename(label=antigen_label,pathogen=variable_name)                                             
tab_meta_labels<-rbind(tab_meta_pathogen_labels,tab_meta_antigen_labels)

### Adjustments ####

adjust_for_1<-paste0("+ age + as.factor(is_female) + as.factor(region_code) + as.factor(highest_education)")
adjust_for_2<-paste0(adjust_for_1, "+ bmi_calc + as.factor(ever_reg_smoker) + as.factor(ever_reg_alcohol) + as.factor(family_cancer) + cluster(csid)")

adjust_for<-adjust_for_2

### Coxmodel - prentice weighted ####

model_fun <- function(exposure, outcome){

  data<-data_all  
  
data$outcome<-data[[paste0(outcome,"_ep")]]
data$outcome_date<-data[[paste0(outcome,"_datedeveloped")]]

data$study_date<-as.numeric(as.Date(data$study_date))  
data$outcome_date<-as.numeric(data$outcome_date)   

### Exclude participants experiencing outcome during first two years of follow-up ####

data$study_time_x2<-data$outcome_date - data$study_date         
data$study_time_x2<-data$study_time_x2-365.25*2
data<-data[data$study_time_x2>0,]               

### Keep subcohort participants and cases with site specific cancer (outcome) of interest ####

data<-data[data$multiplexserology_cancer_b1_subcohort==1  | data$outcome==1,]

###################################################################################################################################################

data$study_time<-data$outcome_date-data$study_date

data$study_date<-0

data$study_date[data$outcome==1 & data$multiplexserology_cancer_b1_subcohort==0]<-data$study_time[data$outcome==1 & data$multiplexserology_cancer_b1_subcohort==0]-0.001

model<-coxph(as.formula(paste0("Surv(study_date,study_time,",outcome,"_ep) ~ ",exposure,"+",adjust_for)),data)

###################################################################################################################################################

        results<- as.data.frame(rbind(estimate   = as.data.frame((summary(model))$coef)[1,1],
                                      stderr     = as.data.frame((summary(model))$coef)[1,3],
                                      p_value    = as.data.frame((summary(model))$coef)[1,5],
                                     pathogen    = exposure,
                                     cancer_code = outcome,
                                     n_cases     = as.data.frame(table(data[[paste0(outcome,"_ep")]]))[2,2]))
                                
}

### Analysis 1a ####

exposures<-grep("_pp$",names(data_all),value=TRUE) 
outcomes<-"ep_cancer_all_combined"

### Run function ####

results<-list()

for(outcome in outcomes){
  
  results[[outcome]]<-list()
  
  for(exposure in exposures){  
    
    results[[outcome]][[exposure]] <- model_fun(exposure, outcome)
    
  }}

tab_1a<-as.data.frame(t(as.data.frame(purrr::list_flatten(results))))           
tab_1a$level<-"pathogen"

### Analysis 1b ####

exposures<-grep("_ap$",names(data_all),value=TRUE) 
outcomes<-"ep_cancer_all_combined"

### Run function ####

results<-list()

for(outcome in outcomes){
  
  results[[outcome]]<-list()
  
  for(exposure in exposures){  
    
    results[[outcome]][[exposure]] <- model_fun(exposure, outcome)
    
  }}

tab_1b<-as.data.frame(t(as.data.frame(purrr::list_flatten(results))))           
tab_1b$level<-"antigen"

### Analysis 2a pathogen ####

exposures<-grep("_pp$",names(data_all),value=TRUE)
exposures<-grep("hbv|hcv",exposures,value=TRUE)
outcomes<-"ep_cancer_C22_combined"      

### Run function ####

results<-list()

for(outcome in outcomes){
  
  results[[outcome]]<-list()
  
  for(exposure in exposures){  
    
    results[[outcome]][[exposure]] <- model_fun(exposure, outcome)
    
  }}

tab_2a<-as.data.frame(t(as.data.frame(purrr::list_flatten(results))))      
tab_2a$level<-"pathogen"

### Analysis 2b antigen ####

exposures<-grep("_ap$",names(data_all),value=TRUE)
exposures<-grep("hbv|hcv",exposures,value=TRUE)
outcomes<-"ep_cancer_C22_combined"      

### Run function ####

results<-list()

for(outcome in outcomes){
  
  results[[outcome]]<-list()
  
  for(exposure in exposures){  
    
    results[[outcome]][[exposure]] <- model_fun(exposure, outcome)
    
  }}

tab_2b<-as.data.frame(t(as.data.frame(purrr::list_flatten(results))))      
tab_2b$level<-"antigen"

### Analysis 2c pathogen ###

exposures<-grep("_pp$",names(data_all),value=TRUE)
exposures<-grep("pylori",exposures,value=TRUE)
outcomes<- c("ep_cancer_C16_combined","ep_cancer_C16.0_combined","ep_cancer_C16x.0_combined","ep_cancer_C16x.0.8_combined") 

### Run function ####

results<-list()

for(outcome in outcomes){
  
  results[[outcome]]<-list()
  
  for(exposure in exposures){  
    
    results[[outcome]][[exposure]] <- model_fun(exposure, outcome)
    
  }}

tab_2c<-as.data.frame(t(as.data.frame(purrr::list_flatten(results))))      
tab_2c$level<-"pathogen"

### Analysis 2d antigen ###

exposures<-grep("_ap$",names(data_all),value=TRUE)
exposures<-grep("pylori",exposures,value=TRUE)
outcomes<- c("ep_cancer_C16_combined","ep_cancer_C16.0_combined","ep_cancer_C16x.0_combined","ep_cancer_C16x.0.8_combined") 

### Run function ####

results<-list()

for(outcome in outcomes){
  
  results[[outcome]]<-list()
  
  for(exposure in exposures){  
    
    results[[outcome]][[exposure]] <- model_fun(exposure, outcome)
    
  }}

tab_2d<-as.data.frame(t(as.data.frame(purrr::list_flatten(results))))      
tab_2d$level<-"antigen"

tab_1<-rbind(tab_1a,tab_1b)
tab_1$analysis<-"cancer_all_cox"
tab_2<-rbind(tab_2a,tab_2b,tab_2c,tab_2d)
tab_2$analysis<-"cancer_selected_cox"
tab<-rbind(tab_1,tab_2)

tab$pathogen<-str_remove_all(tab$pathogen,"_pp$") 
tab$pathogen<-str_remove_all(tab$pathogen,"_ap$") 
tab<-merge(tab,tab_meta_labels,by="pathogen")
tab$estimate<-as.numeric(tab$estimate)
tab$stderr<-as.numeric(tab$stderr)

saveRDS(tab,paste0(path_2,"tab_cancer_all_cox.rds"))
#saveRDS(tab,paste0(path_2,"tab_cancer_all_cox_x2.rds"))




