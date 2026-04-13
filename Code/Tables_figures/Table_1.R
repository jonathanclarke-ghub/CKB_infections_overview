
#  CHARACTERISTICS ####

library(DirectStandardisation)
library(flextable)
library(officer)
library(dplyr)
library(tidyverse)
library(plotrix)
library(ltm)
library(DT)
library(tidyverse)
library(tidytext)
library(reshape) 
library(dplyr)
library(flextable)
library(officer)
library(plotrix)
library(ltm)

### Paths ####

path_1<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Raw_data_7/"
path_2<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Project_data_7/"
path_3<-"K:/kadoorie/Staff_Folders/JonathanC/Projects/Serology/Tables_figures/Tables_figures_10/"


## Import dataset #### 

data<-readRDS(paste0(path_2,"data_serology_2024_7.rds"))
data<-as.data.frame(data)

summary(data$random_glucose_x10)                                    # Move to data prep
data$random_glucose<-data$random_glucose_x10/10                     # Move to data prep

data$study_pop<-1

data_ckb<-readRDS(paste0(path_2,"data_ckb_7.rds"))
data_ckb<-as.data.frame(data_ckb)

summary(data_ckb$random_glucose_x10)                                # Move to data prep
data_ckb$random_glucose<-data_ckb$random_glucose_x10/10             # Move to data prep

data_ckb$study_pop<-1

data$region_code<-as.factor(data$region_code)
data_ckb$region_code<-as.factor(data_ckb$region_code)

table(data_ckb$education_cat2)
table(data$education_cat2[data$multiplexserology_cancer_b1_subcohort==1])

table(data_ckb$household_income_cat2)
table(data$household_income_cat2[data$multiplexserology_cancer_b1_subcohort==1])

###  Baseline table adjusted by age, sex and region  by case-subcohort status #### 

### Adjusted means #### 

results_mn<-as.data.frame(t(adjmeans(dataset = data, outcome_vars = c("bmi_calc",
                                                                      "met",
                                                                      "waist_cm",
                                                                      "sbp_mean",
                                                                      "random_glucose"), 
                                     categorical_vars = "multiplexserology_cancer_b1_subcohort",
                                     adjustment_vars = c("age_cat5","is_female","region_code"), adjustment_var_labels = c("Age","Sex","Area"),
                                     title = "")))

results_mn_age<-as.data.frame(t(adjmeans(dataset = data, outcome_vars = c("age"), 
                                     categorical_vars = "multiplexserology_cancer_b1_subcohort",
                                     adjustment_vars = c("is_female","region_code"), adjustment_var_labels = c("Age","Sex","Area"),
                                     title = "")))


results_mn_all<-as.data.frame(t(adjmeans(dataset = data_ckb, outcome_vars = c("bmi_calc",
                                                                              "met",
                                                                              "waist_cm",
                                                                              "sbp_mean",
                                                                              "random_glucose"), 
                                         categorical_vars = "study_pop",      
                                         adjustment_vars = c("age_cat5","is_female","region_code"), adjustment_var_labels = c("Age","Sex","Area"),
                                         title = "")))

results_mn_all_age<-as.data.frame(t(adjmeans(dataset = data_ckb, outcome_vars = c("age"), 
                                         categorical_vars = "study_pop",      
                                         adjustment_vars = c("is_female","region_code"), adjustment_var_labels = c("Age","Sex","Area"),
                                         title = "")))


### Adjusted proportions ####

results_prop<-as.data.frame(t(adjprop(dataset = data, outcome_vars = c("education_cat2",
                                                                       "household_income_cat2",
                                                                       "chd_stroke_diag",
                                                                       "diabetes_diag",
                                                                       "cirrhosis_hep_diag",
                                                                       "emph_bronc_diag",
                                                                       "tb_diag",
                                                                       "peptic_ulcer_diag",
                                                                       "hep_b_diag",
                                                                       "blood_transfusion",
                                                                       "srh_poor",
                                                                       "family_cancer"), 
                      categorical_vars = "multiplexserology_cancer_b1_subcohort",
                      adjustment_vars = c("age_cat5","is_female","region_code"), adjustment_var_labels = c("Age","Sex","Area"),
                      title = "")))

results_prop_all<-t(adjprop(dataset = data_ckb, outcome_vars = c("education_cat2",
                                                              "household_income_cat2",
                                                              "chd_stroke_diag",
                                                              "diabetes_diag",
                                                              "cirrhosis_hep_diag",
                                                              "emph_bronc_diag",
                                                              "tb_diag",
                                                              "peptic_ulcer_diag",
                                                              "hep_b_diag",
                                                              "blood_transfusion",
                                                              "srh_poor",
                                                              "family_cancer"), 
                   categorical_vars = "study_pop",
                   adjustment_vars = c("age_cat5","is_female","region_code"), adjustment_var_labels = c("Age","Sex","Area"),
                   title = ""))

### Adjusted proportions - female ####

results_prop_sex<-as.data.frame(t(adjprop(dataset = data, outcome_vars = c("is_female"), 
                   categorical_vars = "multiplexserology_cancer_b1_subcohort",
                   adjustment_vars = c("age_cat5","region_code"), adjustment_var_labels = c("Age","Area"),
                   title = "")))

results_prop_sex_all<-t(adjprop(dataset = data_ckb, outcome_vars = c("is_female"), 
                   categorical_vars = "study_pop",
                   adjustment_vars = c("age_cat5","region_code"), adjustment_var_labels = c("Age","Area"),
                   title = ""))

### Adjusted proportions - region ####

results_prop_area<-as.data.frame(t(adjprop(dataset = data, outcome_vars = c("region_is_urban"), 
                      categorical_vars = "multiplexserology_cancer_b1_subcohort",
                      adjustment_vars = c("age_cat5","is_female"), adjustment_var_labels = c("Age","Sex"),
                      title = "")))

results_prop_area_all<-t(adjprop(dataset = data_ckb, outcome_vars = c("region_is_urban"), 
                        categorical_vars = "study_pop",
                        adjustment_vars = c("age_cat5","is_female"), adjustment_var_labels = c("Age","Sex"),
                        title = ""))

### Adjusted proportions - men only ####

data_m<-data[data$is_female==0,]
data_ckb_m<-data_ckb[data_ckb$is_female==0,]

results_prop_m<-as.data.frame(t(adjprop(dataset = data_m, outcome_vars = c("ever_reg_smoker",
                                                                           "ever_reg_alcohol"), 
                   categorical_vars = "multiplexserology_cancer_b1_subcohort",
                   adjustment_vars = c("age_cat5","is_female","region_code"), adjustment_var_labels = c("Age","Sex"),
                   title = "")))


results_prop_m_all<-t(adjprop(dataset = data_ckb_m, outcome_vars = c("ever_reg_smoker",
                                                                 "ever_reg_alcohol"), 
                   categorical_vars = "study_pop",
                   adjustment_vars = c("age_cat5","is_female","region_code"), adjustment_var_labels = c("Age","Sex"),
                   title = ""))

### Adjusted proportions - women only ####

data_f<-data[data$is_female==1,]
data_ckb_f<-data_ckb[data_ckb$is_female==1,]

results_prop_f<-as.data.frame(t(adjprop(dataset = data_f, outcome_vars = c("ever_reg_smoker",
                                                                           "ever_reg_alcohol"), 
                                              categorical_vars = "multiplexserology_cancer_b1_subcohort",
                                              adjustment_vars = c("age_cat5","is_female","region_code"), adjustment_var_labels = c("Age","Sex"),
                                              title = "")))


results_prop_f_all<-t(adjprop(dataset = data_ckb_f, outcome_vars = c("ever_reg_smoker",
                                                                 "ever_reg_alcohol"), 
                                    categorical_vars = "study_pop",
                                    adjustment_vars = c("age_cat5","is_female","region_code"), adjustment_var_labels = c("Age","Sex"),
                                    title = ""))

### Proportions summary ####

tab_prop_levels<-rbind(results_prop,
                       results_prop_area,
                       results_prop_sex,
                       results_prop_m,
                       results_prop_f)

tab_prop_all<-rbind(results_prop_all,
                    results_prop_area_all,
                    results_prop_sex_all,
                    results_prop_m_all,
                    results_prop_f_all)

tab_prop<-cbind(tab_prop_levels,tab_prop_all)
tab_prop<-tab_prop[grep("proportion",row.names(tab_prop),value=TRUE),]

level<-names(tab_prop)

for(ii in level){

tab_prop[[ii]]<-as.numeric(tab_prop[[ii]])
tab_prop[[ii]]<-format(round((tab_prop[[ii]]*100),1),digits=1,nsmall=1,scientific=FALSE)
}

labels<-c(NA,"Schooling > 6 years", 
          NA,"Annual household income < 20,000 (yuan)",
          NA,"CHD or Stroke",
          NA,"Diabetes",
          NA,"Cirrhosis/Hepatitis",
          NA,"Emphysema/Bronchitis",
          NA,"Tuberculosis",
          NA,"Peptic ulcer" ,   
          NA,"HBsAg+, %", 
          NA,"Blood transfusion, %"  , 
          NA,"Poor self-rated health",
          NA,"Family history of cancer, %" ,
          NA,"Urban",
          NA,"Women",
          NA,"Men",
          NA,"Men",
          NA,"Women",
          NA,"Women")

tab_prop<-cbind(labels,tab_prop)
tab_prop<-tab_prop %>%
  filter(!is.na(labels))

### Means ####

results_mn<-rbind(results_mn_age,results_mn)
results_mn_all<-rbind(results_mn_all_age,results_mn_all)
tab_mn<-cbind(results_mn,results_mn_all)

tab_mn<-tab_mn[grep("mean$",row.names(tab_mn),value=TRUE),]

level<-names(tab_mn)

for(ii in level){
  
  tab_mn[[ii]]<-as.numeric(tab_mn[[ii]])
  tab_mn[[ii]]<-format(round(tab_mn[[ii]],1),digits=1,nsmall=1,scientific=FALSE)
}


### Standard deviations (not standardised) ####

var_con<-c("age",
           "bmi_calc",
           "met",
           "waist_cm",
           "sbp_mean",
           "random_glucose")

sd<-list()
sd_all<-list()

for(ii in var_con) {
  
  sd[[ii]]<-format(tapply(data[[ii]], data$multiplexserology_cancer_b1_subcohort, sd, na.rm = TRUE),digits=1,nsmall=1,scientific=FALSE)
  sd_all[[ii]]<-format(tapply(data[[ii]], data$study_pop, sd, na.rm = TRUE),digits=1,nsmall=1,scientific=FALSE)
  
}

sd<-do.call(rbind,sd)
sd_all<-do.call(rbind,sd_all)

tab_sd<-cbind(sd,sd_all)

tab_mn_sd<-cbind(tab_mn,tab_sd)

for(ii in 1:length(level)){
  
  tab_mn_sd[[ii]]<-paste0(tab_mn_sd[,ii]," (", tab_mn_sd[,I(ii + length(level))],")")
  
}

tab_mn_sd<-tab_mn_sd[,1:length(level)]

labels<-c("Age, years, Mean (SD)",
          "BMI, kg/m\u00B2, Mean (SD)",
          "Physical activity, MET-hrs/day, mean (SD)",
          "Waist circumference, cm, mean (SD)",
          "SBP, mmHg, mean (SD)",
          "RPG, mmol/L, mean (SD)")

tab_mn_sd<-cbind(labels,tab_mn_sd)

tab<-rbind(tab_mn_sd,tab_prop)

col_names<-c("varnames","Incident Cancer","Subcohort","All")
names(tab)<-col_names

### Number of participants in each arm ####

row_n<-as.vector(results_prop[4,])
row_n<-as.numeric(unlist(row_n,use.names=FALSE))
row_n<-c(row_n,length(data_ckb$study_pop))
row_n<-format(row_n,big.mark=",",scientific = FALSE)

row_n1<-paste0("N=",row_n[1])
row_n2<-paste0("N=",row_n[2])
row_n3<-paste0("N=",row_n[3])

### Table of results ####

sub1<-c("Socio-demographic, lifestyle, physical measures",NA,NA,NA)
sub2<-c("Ever regular smoker, %",NA,NA,NA)
sub3<-c("Ever regular alcohol drinker, %",NA,NA,NA)
sub4<-c("Prior diseases reported at baseline, %",NA,NA,NA)

sub<-as.data.frame(rbind(sub1,sub2,sub3,sub4))

names(sub)<-names(tab)

tab<-rbind(tab,sub)

n<-c(1:length(tab$varnames))
tab<-cbind(n,tab)

### Create bespoke table #### 

tab<-tab[c(25,1,20,19,7,8,26,21,23,27,22,24,3,2,4,5,6,15,16,
           28,9,10,11,12,13,14,17,18),]

tab<-tab[,c(2:5)]

col_names<-names(tab)

### Create flextable ####

set_flextable_defaults(font.family = "Arial")

header <- data.frame(col_keys = col_names,
                     title   =rep("Table 1: Baselines characteristics of CKB participants in multiplex infection cancer sub-study",4),
                     head_1  = c("Characteristics","Cancer cases","Subcohort","All CKB participants"),
                     head_2  = c("Characteristics",paste0("(",row_n1,")"),paste0("(",row_n2,")"),paste0("(",row_n3,")")))

ft<-flextable(tab, col_keys=c(col_names)) %>%
  set_header_df(mapping = header, key = "col_keys") %>%
  merge_h(part = "header") %>%
  merge_v(part = "header") %>%
  bold(i=c(1:2),j=1:4,part = "header") %>% 
  align(i=1:3,j=2:4,align = "center", part = "header") %>%  
  align(i=1:2,j=1,align = "left", part = "header") %>%   
  hline(i=1,j=NULL,border = fp_border(width = 1, color = "black"), part = "header" ) %>%
  padding(i=c(2:6,7,10,13:19,21:27), j=1, padding.left=15) %>%
  padding(i=c(8:9,11:12), j=1, padding.left=30) %>%
  width(j=c(1,2:4),width=c(4,1.2,1.2,1.2)) %>%
  align(i=NULL,j=c(2:4), align="center", part="body") %>%  
  bold(i=c(1,20,28),j=1, part="body") %>%  
  hline_top(border = fp_border(width =1, color = "black"), part = "body" ) %>%  
  hline_bottom(border = fp_border(width =1, color = "black"), part = "body" ) %>%
  add_footer_lines("Adjusted for age (10-year age groups), sex and region (10 regions) where appropriate.") %>% 
  add_footer_lines("MET: Metabolic Equivalent of Task, BMI: Body Mass Index, CHD: Chronic Heart Disease, SBP: Systolic Blood Pressure, \nRPG: Random Plasma Glucose.") %>%  
  line_spacing(i = NULL, j = NULL, space = 0.6, part = "body") %>%
  line_spacing(i = 1, j = NULL, space = 0.6, part = "footer") %>% 
  fontsize(i = NULL, j = NULL, size = 9, part = "footer") 

# Define page size for output

sect_properties <- prop_section(page_size = page_size(orient = "portrait", width = 13, height = 10),
                                type = "continuous", page_margins = page_mar())

save_as_docx(ft,path="K:/kadoorie/Staff_Folders/JonathanC/Projects/Serology/Tables_figures/Tables_figures_10/Table_1.docx",
             pr_section=sect_properties)


