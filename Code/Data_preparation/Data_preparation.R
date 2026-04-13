### Data ###

library(stringr)
library(tidyverse)
library(ggplot2)
library(DescTools)
library(plotrix)
library(flextable)
library(officer)
library(ggplot2)
library(gridExtra)
library(lattice)
library(patchwork)
library(grid)
library(lubridate)
library(tidyverse)
library(Epi)
library(survival)
library(officer)
library(dplyr)
library(plotrix)
library(broom)
require(nnet)
library(ckbplotr)
library(gridExtra)
library(lattice)
library(patchwork)
library(grid)
library(qvcalc)
library(trend)
library(ggpubr)

# PROGRAM TO IMPORT AND CATEGORISE ADIPOSITY BASELINE DATASET ####

## Set paths ####

path_1<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Raw_data_7/"
path_2<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Project_data_7/"
path_3<-"K:/kadoorie/Staff_Folders/JonathanC/Projects/serology/Tables_figures/Tables_figures_7/"

## Import baseline data ####

data_base<-readRDS(paste0(path_1,"data_baseline_questionnaires.rds"))

# Import data dictionary

column_details<-read.csv(paste0(path_1,"column_details.csv"))
value_definitions<-read.csv(paste0(path_1,"value_definitions.csv"))

### Import multiplex data ####

data_multiplex<-readRDS(paste0(path_1,"data_baseline_multiplex_cancer_case_cohort.rds"))

for(ii in names(data_multiplex)[6:53]){print(summary(data_multiplex[[ii]]))}

### Import ascertainments dataset ####

data_ascertainments<-readRDS(paste0(path_1,"data_baseline_ascertainments.rds"))

# Merge datasets 

data_list <- list(data_base,
                  data_multiplex,
                  data_ascertainments)      

data <- data_list %>% reduce(inner_join, by='csid')

# Select participants with serology assay date

data<-data[!is.na(data$assay_date),]

table(data$multiplexserology_cancer_b1_subcohort,useNA="ifany")

# Add 91 missing plate participants to subcohort 

data$multiplexserology_cancer_b1_subcohort[is.na(data$multiplexserology_cancer_b1_subcohort)]<-1

table(data$multiplexserology_cancer_b1_subcohort,useNA="ifany")

### Define sero-positivity for each antigen using pre-defined cutoffs ####  

### Import multiplex meta data ####

tab_meta<-readRDS(paste0(path_1,"multiplex_cancer_case_cohort_meta.rds"))
tab_meta[tab_meta$variable_name=="hpv16_e6","proposed_cutoff"]<-1000      # Updated cut-off

### Define sero-positive status for each antigen ####

antigen<-tab_meta$variable_name

for(ii in antigen) {
print(summary(data[[ii]]))
data[[paste0(ii,"_ap")]]<-ifelse(data[[ii]]>=as.numeric(tab_meta[tab_meta$variable_name==ii,"proposed_cutoff"]),1,0)
print(table(data[[paste0(ii,"_ap")]],useNA="ifany"))
}

### CREATE PATHOGEN VARIABLES ####

### Pathogen sero-positivity definitions ####

grep("_ap",names(data),value=TRUE)

# Select single antigen for HPV16 and HPV18

table(data$hpv16_l1_ap)             
data$hpv16_pp<-data$hpv16_l1_ap
table(data$hpv16_pp)

table(data$hpv18_l1_ap) 
data$hpv18_pp<-data$hpv18_l1_ap  
table(data$hpv18_pp)

# Modify single antigen - pathogen names

data$hsv1_pp<-data$hsv1_gg_ap   
data$hsv2_pp<-data$hsv2_mggunique_ap 
data$bkv_pp<-data$bk_vp1_ap            
data$jcv_pp<-data$jc_vp1_ap            
data$mcv_pp<-data$mc_vp1_ap    
data$vzv_pp<-data$vzv_gegi_ap
data$hhv7_pp<-data$hhv7_u14_ap   
data$ctrachomatis_pp<-data$ctrachomatis_pgp3_ap
data$cburnetii_pp<-data$cburnetii_com1_ap

# Multiple antigen pathogens (algorithmically defined)

# CMV - 2 or more out of 3

data$cmv_index<-data$cmv_pp150n_ap + data$cmv_pp28_ap + data$cmv_pp52_ap
table(data$cmv_index,useNA="ifany")
data$cmv_pp<-ifelse(data$cmv_index > 1,1,0)
table(data$cmv_pp,useNA="ifany")

# HHV6 - 1 or more out of 2

data$hhv6_index  <- data$hhv6_ie1a_ap + data$hhv6_ie1b_ap 
table(data$hhv6_index,useNA="ifany")
data$hhv6_pp  <- ifelse(data$hhv6_index > 0,1,0)
table(data$hhv6_pp,useNA="ifany")

# HIV - 2 out of 2

data$hiv_index  <- data$hiv_env_ap + data$hiv_gag_ap 
table(data$hiv_index,useNA="ifany")
data$hiv_pp  <- ifelse(data$hiv_index > 1,1,0)
table(data$hiv_pp,useNA="ifany")

# HTLV - 2 out of 2

data$htlv1_index  <- data$htlv1_env_ap + data$htlv1_gag_ap 
table(data$htlv1_index,useNA="ifany")
data$htlv1_pp  <- ifelse(data$htlv1_index > 1,1,0)
table(data$htlv1_pp,useNA="ifany")

# Hpylori - 4 or more out of 8

data$hpylori_index<-data$hpylori_cagan_ap + data$hpylori_catalase_ap + data$hpylori_groel_ap + data$hpylori_hcpc_ap + 
  data$hpylori_hopa_ap + data$hpylori_hp0305_ap + data$hpylori_hp1564_ap + data$hpylori_vacac_ap
table(data$hpylori_index,useNA="ifany")
data$hpylori_pp  <- ifelse(data$hpylori_index > 3,1,0)
table(data$hpylori_pp,useNA="ifany")

# Tgondii - 1 out of 2

data$tgondii_index<-data$tgondii_p22trunc_ap + data$tgondii_sag1d1_ap
table(data$tgondii_index,useNA="ifany")
data$tgondii_pp  <- ifelse(data$tgondii_index > 0,1,0)
table(data$tgondii_pp,useNA="ifany")

# HBV - 2 out of 2

data$hbv_index<-data$hbv_hbc_ap + data$hbv_hbe_ap
table(data$hbv_index,useNA="ifany")
data$hbv_pp  <- ifelse(data$hbv_index > 1,1,0)
table(data$hbv_pp,useNA="ifany")

# HCV - 2 out of 2

data$hcv_index<-data$hcv_ns3_ap + data$hcv_core_ap
table(data$hcv_index,useNA="ifany")
data$hcv_pp  <- ifelse(data$hcv_index > 1,1,0)
table(data$hcv_pp,useNA="ifany")

# EBV - 2 or more out of 4

data$ebv_index<-data$ebv_ebna1_ap + data$ebv_vcap18_ap + data$ebv_ead_ap + data$ebv_zebra_ap
table(data$ebv_index,useNA="ifany")

data$ebv_pp<-ifelse(data$ebv_index > 1,1,0)
table(data$ebv_pp,useNA="ifany")

# Create a pathogen index variable

data$pathogen_index<-data$hbv_pp + 
  data$hcv_pp +
  data$hpv16_pp +
  data$hpv18_pp +
  data$hsv1_pp +
  data$hsv2_pp +
  data$vzv_pp +
  data$ebv_pp +
  data$cmv_pp + 
  data$hhv6_pp + 
  data$hhv7_pp + 
  data$bkv_pp +
  data$jcv_pp +
  data$mcv_pp + 
  data$hiv_pp + 
  data$htlv1_pp + 
  data$hpylori_pp + 
  data$ctrachomatis_pp + 
  data$tgondii_pp + 
  data$cburnetii_pp

table(data$pathogen_index,useNA="ifany")  

### Create a > 10 pathogen variable ####

data$pathogen_11_plus<-ifelse(data$pathogen_index>10,1,0)  

table(data$pathogen_11_plus,useNA="ifany")  

### Create an oncogenic pathogen index ####

data$oc_pathogen_index<-data$hbv_pp + 
  data$hcv_pp +
  data$hpv16_pp +
  data$hpv18_pp +
  data$ebv_pp +
  data$hiv_pp + 
  data$htlv1_pp + 
  data$hpylori_pp

table(data$oc_pathogen_index,useNA="ifany")  

data$oc_pathogen_3_plus<-ifelse(data$oc_pathogen_index>2,1,0)  

table(data$oc_pathogen_3_plus,useNA="ifany")  

## Import baseline data ####

#data<-readRDS(paste0(path_1,"data_baseline_questionnaires.rds"))

### Derive baseline variables ####

# Age

data$age <- as.numeric(floor(data$age_at_study_date_x100 / 100))

data$age_cat5<-as.factor(ifelse(data$age<=39,0,ifelse(data$age<=49,1,ifelse(data$age<=59,2,ifelse(data$age<=69,3,4)))))
table(data$age_cat5,useNA="ifany")

# Fasting time

data$fasting_time<-(data$hours_since_last_ate_x10)/10

# Sex

table(data$is_female)

# Year of birth 

table(data$dob_y)
class(data$dob_y)

data$dob_y[data$dob_y<1931]<-1931
data$dob_y[data$dob_y>1971]<-1971

data$birth_cohort<-as.factor(data$dob_y)

# Birth cohort

table(data$birth_cohort)
data$birth_cohort_cat3<-ifelse(data$dob_y<=1949,0,ifelse(data$dob_y<=1959,1,2))
table(data$birth_cohort_cat3,useNA="ifany")

table(data$birth_cohort)
data$birth_cohort_cat4<-ifelse(data$dob_y<=1939,0,ifelse(data$dob_y<=1949,1,ifelse(data$dob_y<=1959,2,3)))
table(data$birth_cohort_cat4,useNA="ifany")

# Area and region

table(data$region_code)

data$region_name <- as.factor(ifelse(data$region_code==12,'Qingdao',
                                     ifelse(data$region_code==16,'Harbin',
                                            ifelse(data$region_code==26,'Haikou',
                                                   ifelse(data$region_code==36,'Suzhou',
                                                          ifelse(data$region_code==46,'Liuzhou',
                                                                 ifelse(data$region_code==52,'Sichuan',
                                                                        ifelse(data$region_code==58,'Gansu',
                                                                               ifelse(data$region_code==88,'Henan',
                                                                                      ifelse(data$region_code==78,'Zhejiang','Hunan'))))))))))

table(data$region_name)
table(data$region_name,data$region_is_urban)

data$region_urban<-as.factor(ifelse(data$region_name %in% c("Harbin","Qingdao","Liuzhou","Haikou","Suzhou"),"Urban","Rural"))

table(data$region_is_urban,useNA="ifany")

# Education (6 years plus)

table(data$highest_education)
data$education_cat2<-ifelse(data$highest_education %in% c(0,1),0,1)
table(data$education_cat2,useNA="ifany")

# Household income (>=20k vs. <20K)

table(data$household_income)
data$household_income_cat2 <- ifelse(data$household_income<=3,0,1)
table(data$household_income_cat2,useNA="ifany")

# Ever regular smoker

table(data$smoking_category)
data$ever_reg_smoker <- ifelse(data$smoking_category %in% c(1,2),0,1)
table(data$ever_reg_smoker,useNA="ifany")
data$ever_reg_smoker_m<-data$ever_reg_smoker
data$ever_reg_smoker_f<-data$ever_reg_smoker

# Ever regular alcohol drinker

table(data$alcohol_category)
data$ever_reg_alcohol <- as.factor(ifelse(data$alcohol_category %in% c(1,3),0,1))
table(data$ever_reg_alcohol,useNA="ifany")

data$ever_reg_alcohol_m<-data$ever_reg_alcohol
data$ever_reg_alcohol_f<-data$ever_reg_alcohol

# Ever regular tea drinker

table(data$tea_weekly)
data$reg_tea_drinker<- ifelse(data$tea_weekly %in% c(0,1),0,1)
table(data$reg_tea_drinker,useNA="ifany")


# BMI NHS (note NA drops out for case-cohort)

summary(data$bmi_calc)
data$bmi_cat4 <- ifelse(data$bmi_calc<18.5,0,
                        ifelse(data$bmi_calc<25,1,
                               ifelse(data$bmi_calc<30,2,3)))
table(data$bmi_cat4,useNA="ifany")

# BMI (Category 3)

summary(data$bmi_calc)
data$bmi_cat3 <- ifelse(data$bmi_calc<25,0,
                               ifelse(data$bmi_calc<30,1,2))
table(data$bmi_cat3,useNA="ifany")

# BMI (Overweight BMI>=25)

summary(data$bmi_calc,useNA="ifany")
data$bmi_overweight<-ifelse(data$bmi_calc<25,0,1)
table(data$bmi_overweight,useNA="ifany")

# Central obesity (WC >=94cm (male) or >=80 cm (women)

summary(data$waist_mm)
data$waist_cm<-data$waist_mm/10
data$waist_cat2_m<-ifelse(data$waist_cm>=94,1,0)
data$waist_cat2_f<-ifelse(data$waist_cm>=80,1,0)
table(data$waist_cat2_m,useNA="ifany")
table(data$waist_cat2_f,useNA="ifany")

# Physical activity (High vs. Low)

hist(data$met)
data$met_high<-ifelse(data$met>40,1,0)
table(data$met_high,useNA="ifany")

# Self-rated poor health

table(data$self_rated_health,useNA="ifany")
data$srh_poor<-ifelse(data$self_rated_health==3,1,0)
table(data$srh_poor,useNA="ifany")

# Famine experience

table(data$diet_had_shortage,useNA="ifany")

# Regular spicy food

table(data$diet_spice,useNA="ifany")
data$reg_spice<-ifelse(data$diet_spice==4,1,0)
table(data$reg_spice,useNA="ifany")

# Regular fresh fruit

table(data$diet_freq_fresh_fruit,useNA="ifany")
data$reg_fresh_fruit<-ifelse(data$diet_freq_fresh_fruit==0,1,0)
table(data$reg_fresh_fruit,useNA="ifany")

# Regular preserved vegetables

table(data$diet_freq_preserved_veg,useNA="ifany")
data$reg_preserved_veg<-ifelse(data$diet_freq_preserved_veg==0,1,0)
table(data$reg_preserved_veg,useNA="ifany")

# Regular meat/poultry

table(data$diet_freq_meat,useNA="ifany")
table(data$diet_freq_poultry,useNA="ifany")
data$reg_meat_poultry<-ifelse(data$diet_freq_meat==0|data$diet_freq_poultry==0,1,0)
table(data$reg_meat_poultry,useNA="ifany")

# CHD or Stroke/TIA  

table(data$chd_diag,useNA="ifany")   
table(data$stroke_or_tia_diag,useNA="ifany")                
data$chd_stroke_diag<-ifelse(data$chd_diag==1|data$stroke_or_tia_diag==1,1,0)              
table(data$chd_stroke_diag,useNA="ifany")

# Hep B Test

table(data$hep_b,useNA="ifany")         
data$hep_b_diag<-ifelse(data$hep_b %in% c(2,3),NA,ifelse(data$hep_b==1,1,0))
table(data$hep_b_diag,useNA="ifany")  

table(data$rheum_heart_dis_diag,useNA="ifany")              # Rheumatic heart disease
table(data$diabetes_diag,useNA="ifany")                     # Diabetes (SR or SD)
table(data$hypertension_diag,useNA="ifany")                 # Hypertension (SR or SD)
table(data$emph_bronc_diag,useNA="ifany")                   # Emphysema/bronchitis
table(data$peptic_ulcer_diag,useNA="ifany")                 # Peptic ulcer
table(data$tb_diag,useNA="ifany")                           # TB
table(data$cirrhosis_hep_diag,useNA="ifany")                # Cirrhosis/Hepatitis
table(data$asthma_diag,useNA="ifany")                       # Asthma
table(data$rheum_arthritis_diag,useNA="ifany")              # Rheumatoid arthritis

# Family history of cancer

table(data$mother_cancer,useNA="ifany")
table(data$father_cancer,useNA="ifany")
table(data$siblings_cancer,useNA="ifany")
data$family_cancer<-ifelse(data$mother_cancer==1|data$father_cancer==1|data$siblings_cancer==1,1,0)
table(data$family_cancer,useNA="ifany")
data$family_cancer[is.na(data$family_cancer)]<-0


# Blood transfusion

table(data$blood_transfusions,useNA="ifany")
data$blood_transfusion<-ifelse(data$blood_transfusions>0,1,0)
table(data$blood_transfusion,useNA="ifany")

# Blood donation

table(data$paid_blood_donations,useNA="ifany")
data$blood_donation<-ifelse(data$paid_blood_donations>0,1,0)
table(data$blood_donation,useNA="ifany")

# Children (women)

table(data$children,useNA="ifany")
data$children_cat2<-ifelse(data$children>1,1,0)
table(data$children_cat2,useNA="ifany")

# Early menarche (< 12 years)

table(data$first_period_age,useNA="ifany")
data$early_menarche<-ifelse(data$first_period_age<12,1,0)
table(data$early_menarche,useNA="ifany")

# Post menopausal

table(data$menopause_status,useNA="ifany")
data$post_menopausal<-ifelse(data$menopause_status==2,1,0)
table(data$post_menopausal,useNA="ifany")

### Depression symptoms ####

table(data$cidi_a1,useNA="ifany")
data$has_cidi_a<-factor(ifelse(data$cidi_a1 %in% c(0,1),1,
                               ifelse(is.na(data$cidi_a1),0,data$cidi_a1)))
table(data$has_cidi_a,useNA="ifany")

### Anxiety symptoms ####

data$has_cidi_b<-factor(ifelse(data$cidi_b1 %in% c(0,1),1,
                               ifelse(is.na(data$cidi_b1),0,data$cidi_b1)))

table(data$has_cidi_b,useNA="ifany")

### MDE ####

table(data$cidi_mdep,useNA="ifany")
data$mdep<-as.factor(ifelse(data$cidi_mdep %in% c(NA,0),0,1))
table(data$mdep,useNA="ifany")

### GAD ####

table(data$cidi_gad,useNA="ifany")
data$gad<-as.factor(ifelse(data$cidi_gad %in% c(NA,0),0,1))
table(data$gad,useNA="ifany")

### Depression (all) ####

data$depression_any<-as.factor(ifelse(data$mdep==1| 
                                        data$gad==1|
                                        data$has_cidi_a==1|
                                        data$has_cidi_b==1,1,0))

table(data$depression_any,useNA="ifany")


data$depression_any_plus<-as.factor(ifelse(data$sad_depressed==1 |  
                                        data$continuous_anxiety==1|
                                        data$mdep==1| 
                                        data$gad==1|
                                        data$has_cidi_a==1|
                                        data$has_cidi_b==1,1,0))

table(data$depression_any_plus,useNA="ifany")

table(data$has_cidi_a,data$sad_depressed)
table(data$cidi_mdep,data$sad_depressed)

data_endpoints<-readRDS(paste0(path_1,"endpoints.rds"))
data_endpoint_definitions<-readRDS(paste0(path_1,"endpoint_definitions.rds"))

data_endpoints<-data_endpoints[,c("csid",
                                  "ep_CKB0001_du_ep","ep_CKB0001_du_datedeveloped",
                                  "ep_CKB0014_combined_ep","ep_CKB0014_combined_datedeveloped",
                                  "ep_CKB0019_combined_ep","ep_CKB0019_combined_datedeveloped")]


data<-left_join(data,data_endpoints,by="csid")

#saveRDS(data,paste0(path_2,"data_ckb_7.rds"))

# Drop non-cases in cancer arm

data$drop<-ifelse(data$multiplexserology_cancer_b1_subcohort==0 & data$ep_CKB0014_combined_ep==0,0,1)
table(data$drop)
data<-data %>% filter(drop==1)

### Import endpoints data ####

data_long<-read.csv(paste0(path_1,"DAR-2025-00313.all_endpoints.csv"))

# Define outcomes of interest

cancer_all        <-unique(grep("^C",data_long$icd10,value=TRUE))
cancer_C22        <-unique(grep("^C22",data_long$icd10,value=TRUE))
cancer_C16        <-unique(grep("^C16",data_long$icd10,value=TRUE))
cancer_C16.0      <-unique(grep("^C16.0",data_long$icd10,value=TRUE))
cancer_C16x.0     <-unique(grep("^C16.1|^C16.2|^C16.3|^C16.4|^C16.5|^C16.6|^C16.7|^C16.8|^C16.9",data_long$icd10,value=TRUE))
cancer_C16x.0.8   <-unique(grep("^C16.1|^C16.2|^C16.3|^C16.4|^C16.5|^C16.6|^C16.7|^C16.8",data_long$icd10,value=TRUE))

outcomes<-c(cancer_all,cancer_C22,cancer_C16,cancer_C16.0,cancer_C16x.0,cancer_C16x.0.8)

# Create event endpoint

data_long$event<-1

# Remove duplicates

data_long<-data_long[!duplicated(data_long[c("csid","icd10")]), ]

# Keep first date developed

data_long$date_developed<-as.Date(data_long$date_developed)

data_long<-data_long %>%
  group_by(csid,icd10) %>%
  filter(date_developed==min(date_developed)) %>%
  ungroup()

data_long<-data_long %>% filter(icd10 %in% unique(grep(str_c(outcomes,collapse="|"),data_long$icd10,value=TRUE)))

# Create endpoint event variables

data_long_event<-data_long
data_long_event$icd10<-paste0("ep_",data_long_event$icd10,"_combined_ep")  # Rename endpoint events
data_long_event<-data_long_event %>% select(csid,icd10,event)
data_wide_event<-data_long_event %>% pivot_wider(names_from = icd10, values_from = event)

# Create endpoint date variables (152)

data_long_date<-data_long
data_long_date$icd10<-paste0("ep_",data_long_date$icd10,"_combined_datedeveloped")     # Rename endpoint date developed
data_long_date<-data_long_date %>% select(csid,icd10,date_developed)
data_wide_date<-data_long_date %>% pivot_wider(names_from = icd10, values_from = date_developed)

# Join event and date data

data_list<-list(data_wide_event,
                data_wide_date,
                data)

data <- data_list %>% reduce(full_join, by='csid')
data<-data[!is.na(data$multiplexserology_cancer_b1_subcohort),]

rm(data_list)
rm(data_long)
rm(data_long_event)
rm(data_long_date)
rm(data_wide_event)
rm(data_wide_date)
rm(data_base)
rm(data_ascertainments)
rm(data_multiplex)

# Replace NULL with zero

endpoints<-grep("_ep",names(data),value=TRUE)

for(ii in endpoints){

  data[[ii]][is.na(data[[ii]])]<-0
  
}

# Replace NULL date developed with censoring date

end_dates<-grep("_datedeveloped",names(data),value=TRUE)

summary(data$censoring_date)

for (ii in end_dates){
  
  data[[ii]][is.na(data[[ii]])]<-data$censoring_date[is.na(data[[ii]])]
  
}

endpoints<-grep("_combined_ep",names(data),value=TRUE)

data<-data %>%
  mutate(ep_cancer_all_combined_ep = ifelse(rowSums(across(grep(str_c(cancer_all,collapse="|"),endpoints,value=TRUE)), na.rm=TRUE)>0,1,0),
         ep_cancer_C22_combined_ep = ifelse(rowSums(across(grep(str_c(cancer_C22,collapse="|"),endpoints,value=TRUE)), na.rm=TRUE)>0,1,0),
         ep_cancer_C16_combined_ep = ifelse(rowSums(across(grep(str_c(cancer_C16,collapse="|"),endpoints,value=TRUE)), na.rm=TRUE)>0,1,0),
         ep_cancer_C16.0_combined_ep = ifelse(rowSums(across(grep(str_c(cancer_C16.0,collapse="|"),endpoints,value=TRUE)), na.rm=TRUE)>0,1,0),  
         ep_cancer_C16x.0_combined_ep = ifelse(rowSums(across(grep(str_c(cancer_C16x.0,collapse="|"),endpoints,value=TRUE)), na.rm=TRUE)>0,1,0),
         ep_cancer_C16x.0.8_combined_ep = ifelse(rowSums(across(grep(str_c(cancer_C16x.0.8,collapse="|"),endpoints,value=TRUE)), na.rm=TRUE)>0,1,0))

endpoints_date<-grep("_combined_datedeveloped",names(data),value=TRUE)

data<-data %>%
  rowwise() %>%
  mutate(ep_cancer_all_combined_datedeveloped = min(c_across(grep(str_c(cancer_all,collapse="|"),endpoints_date,value=TRUE))),
         ep_cancer_C22_combined_datedeveloped = min(c_across(grep(str_c(cancer_C22,collapse="|"),endpoints_date,value=TRUE))),
         ep_cancer_C16_combined_datedeveloped = min(c_across(grep(str_c(cancer_C16,collapse="|"),endpoints_date,value=TRUE))),
         ep_cancer_C16.0_combined_datedeveloped = min(c_across(grep(str_c(cancer_C16.0,collapse="|"),endpoints_date,value=TRUE))),
         ep_cancer_C16x.0_combined_datedeveloped = min(c_across(grep(str_c(cancer_C16x.0,collapse="|"),endpoints_date,value=TRUE))),
         ep_cancer_C16x.0.8_combined_datedeveloped = min(c_across(grep(str_c(cancer_C16x.0.8,collapse="|"),endpoints_date,value=TRUE))))

saveRDS(data,paste0(path_2,"data_serology_2024_7.rds"))

#######################################################################################################################################

### Table - cases #### (add in extra case definitions)

outcomes<-c("cancer_all","cancer_C22","cancer_C16","cancer_C16.0","cancer_C16x.0","cancer_C16x.0.8","CKB0014","CKB0019")

cases_all<-list()

for(ii in outcomes) {

  n<-as.data.frame(table(data[[paste0("ep_",ii,"_combined_ep")]]))[2,2]
  cases_all[[ii]]<-c(ii,n)
}

tab_cases_all<-as.data.frame(do.call(rbind,cases_all))
names(tab_cases_all)<-c("Disease","cases - all")

# Subcohort cases
  
data_sub<-data %>% filter(multiplexserology_cancer_b1_subcohort==1)

cases_sub<-list()

for(ii in outcomes) {
  
  n<-as.data.frame(table(data_sub[[paste0("ep_",ii,"_combined_ep")]]))[2,2]
  cases_sub[[ii]]<-c(ii,n)
}

tab_cases_sub<-as.data.frame(do.call(rbind,cases_sub))
names(tab_cases_sub)<-c("Disease","cases - subcohort")

# Cancer cases

data_cc<-data %>% filter(multiplexserology_cancer_b1_subcohort==0)

cases_cc<-list()

for(ii in outcomes) {
 
  if(ii %in% c("cancer_C22","cancer_C16","cancer_C16.0","cancer_C16x.0","cancer_C16x.0.8","CKB0019")){
   
  n<-as.data.frame(table(data_cc[[paste0("ep_",ii,"_combined_ep")]]))[2,2]
  cases_cc[[ii]]<-c(ii,n)
  
  }
  
  else if(ii %in% c("cancer_all","CKB0014")){
    n<-as.data.frame(table(data_cc[[paste0("ep_",ii,"_combined_ep")]]))[1,2]
    cases_cc[[ii]]<-c(ii,n)
  }
}

tab_cases_cc<-as.data.frame(do.call(rbind,cases_cc))
names(tab_cases_cc)<-c("Disease","cases - cancer cases")

tab_list<-list(tab_cases_sub,tab_cases_cc,tab_cases_all)

tab_cases <- tab_list %>% reduce(inner_join, by='Disease')

write.csv(tab_cases,paste0(path_3,"tab_cases_v7t2_short_list.csv"))

data_long<-read.csv(paste0(path_1,"DAR-2025-00313.all_endpoints.csv"))
data_long<-data_long[!duplicated(data_long[c("csid","icd10")]), ]

tab<-as.data.frame(table(data_long$icd10))

write.csv(tab_cases,paste0(path_3,"tab_cases_v7t_long_list.csv"))

####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################

### Import endpoints data ####

data_long<-read.csv(paste0(path_1,"DAR-2025-00105.all_endpoints.csv"))

data_long<-data_long[!duplicated(data_long[c("csid","icd10")]), ]
data_long<-data_long %>%
  filter(icd10 %in% grep("^C",data_long$icd10,value=TRUE))

data_long<-data_long[!duplicated(data_long[c("csid")]), ]

data<-readRDS(paste0(path_2,"data_serology_2024_5.rds"))
data<-data[data$multiplexserology_cancer_b1_subcohort==0,]

test<-dplyr::left_join(data,data_long,by="csid")

table(test$icd10,useNA="ifany")


