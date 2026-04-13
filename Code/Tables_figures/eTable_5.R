
library(tidyverse)
library(DescTools)
library(plotrix)
library(flextable)
library(officer)
library(stringr)
library(tidyverse)
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
library(ggplot2)
library(patchwork)
library(grid)
library(qvcalc)
library(trend)
library(ggpubr)

### Hetrogeneity test ####

heterogeneity <- function(beta, se) {
  
  # beta must be a vector of estimates
  # se must be a vector of their standard errors
  
  # degrees of freedom	
  df <- length(beta) - 1
  
  # expected_beta: inverse variance weighted average of betas
  expected_beta <- sum(beta / se^2) / sum(1 / se^2)
  
  heterogeneity_test_statistic <- sum(((beta - expected_beta) / se)^2)
  
  p <- pchisq(heterogeneity_test_statistic, df = df, lower = FALSE)
  
  return(list("test statistic" = heterogeneity_test_statistic, "prob" = p))
  
} 

### Paths ####

path_1<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Raw_data_7/"
path_2<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Project_data_7/"
path_3<-"K:/kadoorie/Staff_Folders/JonathanC/Projects/Serology/Tables_figures/Tables_figures_10/"

### Import simulated data ####

data<-readRDS(paste0(path_2,"data_serology_2024_7.rds"))

### Import multiplex meta data ####

tab_meta<-readRDS(paste0(path_1,"multiplex_cancer_case_cohort_meta.rds"))
tab_meta$pathogen_label_old<-tab_meta$pathogen_label
tab_meta$pathogen_label[36:47]<-""
tab_meta$pathogen_label[36]<-"C.trachomatis"
tab_meta$pathogen_label[37:38]<-"T.gondii"
tab_meta$pathogen_label[39]<-"C.burnetii"
tab_meta$pathogen_label[40:47]<-"H.pylori"

### Table - descriptive, by sex and overall ####

pathogen_pos<-grep("_pp$",names(data),value=TRUE)

group<-c("male_sub","female_sub","male_can","female_can","male","female",
         "rural_sub","urban_sub","rural_can","urban_can","rural","urban",
         "40s_sub","50s_sub","60s_sub","40s_can","50s_can","60s_can","40s","50s","60s",
         "sub_all","can_all","all")

data_all<-data

tab<-list()

fun<-function(group){

  if      (group=="male_sub")   {data<-data_all %>% filter(is_female==0 & multiplexserology_cancer_b1_subcohort==1)}
  else if (group=="female_sub") {data<-data_all %>% filter(is_female==1 & multiplexserology_cancer_b1_subcohort==1)}
  else if (group=="male_can")   {data<-data_all %>% filter(is_female==0 & multiplexserology_cancer_b1_subcohort==0)}
  else if (group=="female_can") {data<-data_all %>% filter(is_female==1 & multiplexserology_cancer_b1_subcohort==0)}
  else if (group=="male")       {data<-data_all %>% filter(is_female==0)}
  else if (group=="female")     {data<-data_all %>% filter(is_female==1)}
  else if (group=="rural_sub")  {data<-data_all %>% filter(region_is_urban==0 & multiplexserology_cancer_b1_subcohort==1)}
  else if (group=="urban_sub")  {data<-data_all %>% filter(region_is_urban==1 & multiplexserology_cancer_b1_subcohort==1)}
  else if (group=="rural_can")  {data<-data_all %>% filter(region_is_urban==0 & multiplexserology_cancer_b1_subcohort==0)}
  else if (group=="urban_can")  {data<-data_all %>% filter(region_is_urban==1 & multiplexserology_cancer_b1_subcohort==0)}
  else if (group=="rural")      {data<-data_all %>% filter(region_is_urban==0)}
  else if (group=="urban")      {data<-data_all %>% filter(region_is_urban==1)}
  else if (group=="40s_sub")    {data<-data_all %>% filter(birth_cohort_cat3==0 & multiplexserology_cancer_b1_subcohort==1)}
  else if (group=="50s_sub")    {data<-data_all %>% filter(birth_cohort_cat3==1 & multiplexserology_cancer_b1_subcohort==1)}  
  else if (group=="60s_sub")    {data<-data_all %>% filter(birth_cohort_cat3==2 & multiplexserology_cancer_b1_subcohort==1)}
  else if (group=="40s_can")    {data<-data_all %>% filter(birth_cohort_cat3==0 & multiplexserology_cancer_b1_subcohort==0)}
  else if (group=="50s_can")    {data<-data_all %>% filter(birth_cohort_cat3==1 & multiplexserology_cancer_b1_subcohort==0)}  
  else if (group=="60s_can")    {data<-data_all %>% filter(birth_cohort_cat3==2 & multiplexserology_cancer_b1_subcohort==0)}
  else if (group=="40s")        {data<-data_all %>% filter(birth_cohort_cat3==0)}
  else if (group=="50s")        {data<-data_all %>% filter(birth_cohort_cat3==1)}
  else if (group=="60s")        {data<-data_all %>% filter(birth_cohort_cat3==2)}
  else if (group=="sub_all")    {data<-data_all %>% filter(multiplexserology_cancer_b1_subcohort==1)}
  else if (group=="can_all")    {data<-data_all %>% filter(multiplexserology_cancer_b1_subcohort==0)}
  else if (group=="all")        {data<-data_all}
  
results<-list()
  
for(ii in pathogen_pos){

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

if      (group=="male_sub")   {names(tab)<-paste0(names(tab),"_m_sub")}
else if (group=="female_sub") {names(tab)<-paste0(names(tab),"_f_sub")}
else if (group=="male_can")   {names(tab)<-paste0(names(tab),"_m_can")}
else if (group=="female_can") {names(tab)<-paste0(names(tab),"_f_can")}
else if (group=="male")       {names(tab)<-paste0(names(tab),"_m")}
else if (group=="female")     {names(tab)<-paste0(names(tab),"_f")}
else if (group=="rural_sub")  {names(tab)<-paste0(names(tab),"_r_sub")}
else if (group=="urban_sub")  {names(tab)<-paste0(names(tab),"_u_sub")}
else if (group=="rural_can")  {names(tab)<-paste0(names(tab),"_r_can")}
else if (group=="urban_can")  {names(tab)<-paste0(names(tab),"_u_can")}
else if (group=="rural")      {names(tab)<-paste0(names(tab),"_r")}
else if (group=="urban")      {names(tab)<-paste0(names(tab),"_u")}
else if (group=="40s_sub")    {names(tab)<-paste0(names(tab),"_40s_sub")}
else if (group=="50s_sub")    {names(tab)<-paste0(names(tab),"_50s_sub")}
else if (group=="60s_sub")    {names(tab)<-paste0(names(tab),"_60s_sub")}
else if (group=="40s_can")    {names(tab)<-paste0(names(tab),"_40s_can")}
else if (group=="50s_can")    {names(tab)<-paste0(names(tab),"_50s_can")}
else if (group=="60s_can")    {names(tab)<-paste0(names(tab),"_60s_can")}
else if (group=="40s")        {names(tab)<-paste0(names(tab),"_40s")}
else if (group=="50s")        {names(tab)<-paste0(names(tab),"_50s")}
else if (group=="60s")        {names(tab)<-paste0(names(tab),"_60s")}
else if (group=="sub_all")    {names(tab)<-paste0(names(tab),"_sub_all")}
else if (group=="can_all")    {names(tab)<-paste0(names(tab),"_can_all")}
else if (group=="all")        {names(tab)<-paste0(names(tab),"_all")}

tab[[group]]<-tab

}

tab_res<-lapply(group,fun)
tab<-do.call(cbind,tab_res)
tab$pathogen<-tab$pathogen_r_sub

tab_meta_selected<-tab_meta[!duplicated(tab_meta$pathogen),c("pathogen","pathogen_label","pathogen_label_long")]

### Heterogeneity test (sex subcohort) ####

pp_list<-list()
for(i in 1:dim(tab)[1]) {
  beta<-as.numeric(c(tab[i,c("prop_m_sub","prop_f_sub")]))
  se<-as.numeric(c(tab[i,c("stderr_m_sub","stderr_f_sub")]))
  het <- heterogeneity(beta, se)
  pp_list[[i]]<-het$p
}
tab$p_value_sex_sub<-unlist(pp_list)
sig_sex_sub<-tab$pathogen[tab$p_value_sex_sub<0.05]
sig_sex_sub<-str_remove_all(sig_sex_sub,"_pp") 
sig_sex_sub<-tab_meta_selected$pathogen_label[tab_meta_selected$pathogen %in% sig_sex_sub]

### Heterogeneity test (sex cancer) ####

pp_list<-list()
for(i in 1:dim(tab)[1]) {
  beta<-as.numeric(c(tab[i,c("prop_m_can","prop_f_can")]))
  se<-as.numeric(c(tab[i,c("stderr_m_can","stderr_f_can")]))
  het <- heterogeneity(beta, se)
  pp_list[[i]]<-het$p
}
tab$p_value_sex_can<-unlist(pp_list)
sig_sex_can<-tab$pathogen[tab$p_value_sex_can<0.05]
sig_sex_can<-str_remove_all(sig_sex_can,"_pp") 
sig_sex_can<-tab_meta_selected$pathogen_label[tab_meta_selected$pathogen %in% sig_sex_can]

### Heterogeneity test (sex all) ####

pp_list<-list()
for(i in 1:dim(tab)[1]) {
  beta<-as.numeric(c(tab[i,c("prop_m","prop_f")]))
  se<-as.numeric(c(tab[i,c("stderr_m","stderr_f")]))
  het <- heterogeneity(beta, se)
  pp_list[[i]]<-het$p
}
tab$p_value_sex_all<-unlist(pp_list)
sig_sex_all<-tab$pathogen[tab$p_value_sex_all<0.05]
sig_sex_all<-str_remove_all(sig_sex_all,"_pp") 
sig_sex_all<-tab_meta_selected$pathogen_label[tab_meta_selected$pathogen %in% sig_sex_all]

### Heterogeneity test (region subcohort) ####

pp_list<-list()
for(i in 1:dim(tab)[1]) {
  beta<-as.numeric(c(tab[i,c("prop_r_sub","prop_u_sub")]))
  se<-as.numeric(c(tab[i,c("stderr_r_sub","stderr_u_sub")]))
  het <- heterogeneity(beta, se)
  pp_list[[i]]<-het$p
}

tab$p_value_region_sub<-unlist(pp_list)
sig_region_sub<-tab$pathogen[tab$p_value_region_sub<0.05]
sig_region_sub<-str_remove_all(sig_region_sub,"_pp") 
sig_region_sub<-tab_meta_selected$pathogen_label[tab_meta_selected$pathogen %in% sig_region_sub]

### Heterogeneity test (region cancer) ####

pp_list<-list()
for(i in 1:dim(tab)[1]) {
  beta<-as.numeric(c(tab[i,c("prop_r_can","prop_u_can")]))
  se<-as.numeric(c(tab[i,c("stderr_r_can","stderr_u_can")]))
  het <- heterogeneity(beta, se)
    pp_list[[i]]<-het$p
}
tab$p_value_region_can<-unlist(pp_list)
sig_region_can<-tab$pathogen[tab$p_value_region_can<0.05]
sig_region_can<-str_remove_all(sig_region_can,"_pp") 
sig_region_can<-tab_meta_selected$pathogen_label[tab_meta_selected$pathogen %in% sig_region_can]

### Heterogeneity test (region all) ####

pp_list<-list()
for(i in 1:dim(tab)[1]) {
  beta<-as.numeric(c(tab[i,c("prop_r","prop_u")]))
  se<-as.numeric(c(tab[i,c("stderr_r","stderr_u")]))
  het <- heterogeneity(beta, se)
  pp_list[[i]]<-het$p
}
tab$p_value_region_all<-unlist(pp_list)
sig_region_all<-tab$pathogen[tab$p_value_region_all<0.05]
sig_region_all<-str_remove_all(sig_region_all,"_pp") 
sig_region_all<-tab_meta_selected$pathogen_label[tab_meta_selected$pathogen %in% sig_region_all]

### Heterogeneity test (birth subcohort) ####

pp_list<-list()
for(i in 1:dim(tab)[1]) {
  beta<-as.numeric(c(tab[i,c("prop_40s_sub","prop_50s_sub","prop_60s_sub")]))
  se<-as.numeric(c(tab[i,c("stderr_40s_sub","stderr_50s_sub","stderr_60s_sub")]))
  het <- heterogeneity(beta, se)
  pp_list[[i]]<-het$p
}

tab$p_value_birth_sub<-unlist(pp_list)
sig_birth_sub<-tab$pathogen[tab$p_value_birth_sub<0.05]
sig_birth_sub<-str_remove_all(sig_birth_sub,"_pp") 
sig_birth_sub<-tab_meta_selected$pathogen_label[tab_meta_selected$pathogen %in% sig_birth_sub]

### Heterogeneity test (birth cancer) ####

pp_list<-list()
for(i in 1:dim(tab)[1]) {
  beta<-as.numeric(c(tab[i,c("prop_40s_can","prop_50s_can","prop_60s_can")]))
  se<-as.numeric(c(tab[i,c("stderr_40s_can","stderr_50s_can","stderr_60s_can")]))
  het <- heterogeneity(beta, se)
  pp_list[[i]]<-het$p
}
tab$p_value_birth_can<-unlist(pp_list)
sig_birth_can<-tab$pathogen[tab$p_value_birth_can<0.05]
sig_birth_can<-str_remove_all(sig_birth_can,"_pp") 
sig_birth_can<-tab_meta_selected$pathogen_label[tab_meta_selected$pathogen %in% sig_birth_can]

### Heterogeneity test (birth all) ####

pp_list<-list()
for(i in 1:dim(tab)[1]) {
  beta<-as.numeric(c(tab[i,c("prop_40s","prop_50s","prop_60s")]))
  se<-as.numeric(c(tab[i,c("stderr_40s","stderr_50s","stderr_60s")]))
  het <- heterogeneity(beta, se)
  pp_list[[i]]<-het$p
}
tab$p_value_birth_all<-unlist(pp_list)
sig_birth_all<-tab$pathogen[tab$p_value_birth_all<0.05]
sig_birth_all<-str_remove_all(sig_birth_all,"_pp") 
sig_birth_all<-tab_meta_selected$pathogen_label[tab_meta_selected$pathogen %in% sig_birth_all]


### Heterogeneity test (all) ####

pp_list<-list()
for(i in 1:dim(tab)[1]) {
  beta<-as.numeric(c(tab[i,c("prop_sub_all","prop_can_all")]))
  se<-as.numeric(c(tab[i,c("stderr_sub_all","stderr_can_all")]))
  het <- heterogeneity(beta, se)
  pp_list[[i]]<-het$p
}
tab$p_value_all<-unlist(pp_list)
sig_all<-tab$pathogen[tab$p_value_all<0.05]
sig_all<-str_remove_all(sig_all,"_pp") 
sig_all<-tab_meta_selected$pathogen_label[tab_meta_selected$pathogen %in% sig_all]

tab$pathogen<-str_remove_all(tab$pathogen,"_pp")
tab_meta<-tab_meta[!duplicated(tab_meta$pathogen),c("pathogen_label","pathogen")]
tab<-merge(tab,tab_meta,by="pathogen")

tab_all<-tab[,c("pathogen_label",
                "prop_m","stderr_m","prop_f","stderr_f","p_value_sex_all",
                "prop_r","stderr_r","prop_u","stderr_u","p_value_region_all",
                "prop_40s","stderr_40s","prop_50s","stderr_50s","prop_60s","stderr_60s","p_value_birth_all",
                "prop_sub_all","stderr_sub_all","prop_can_all","stderr_can_all","p_value_all")]

tab<-tab[,c("pathogen_label",grep("prev_se",names(tab),value=TRUE))]

sub_1<-rep(c("Hepatitis virus",NA),c(1,24))
sub_2<-rep(c("Human Papillomavirus",NA),c(1,24))
sub_3<-rep(c("Human Herpes Virus",NA),c(1,24))
sub_4<-rep(c("Human Polyomavirus",NA),c(1,24))       
sub_5<-rep(c("Human Retrovirus",NA),c(1,24))   
sub_6<-rep(c("Bacteria and Parasite",NA),c(1,24))

tab<-rbind(tab,sub_1,sub_2,sub_3,sub_4,sub_5,sub_6)

tab<-tab[c(23,14,15,20,5,3,8,9,
           21,6:7,
           25,16,10,  
           24,1,17,18,
           22,11:12,
           26,4,13,2,19),]

tab$n<-1:dim(tab)[1]

n_sig_sex_sub<-tab$n[tab$pathogen_label %in% sig_sex_sub]
n_sig_sex_can<-tab$n[tab$pathogen_label %in% sig_sex_can]
n_sig_sex_all<-tab$n[tab$pathogen_label %in% sig_sex_all]
n_sig_region_sub<-tab$n[tab$pathogen_label %in% sig_region_sub]
n_sig_region_can<-tab$n[tab$pathogen_label %in% sig_region_can]
n_sig_region_all<-tab$n[tab$pathogen_label %in% sig_region_all]
n_sig_birth_sub<-tab$n[tab$pathogen_label %in% sig_birth_sub]
n_sig_birth_can<-tab$n[tab$pathogen_label %in% sig_birth_can]
n_sig_birth_all<-tab$n[tab$pathogen_label %in% sig_birth_all]

tab<-tab %>% dplyr::select(-n)

gap_1<-rep(c(""),26)
gap_2<-rep(c(""),26)
gap_3<-rep(c(""),26)

tab<-as.data.frame(cbind(tab[,c(1,4:5)],gap_1,tab[,10:11],gap_2,tab[,17:19],gap_3,tab[,24])) # Cancer

col_names<-names(tab)

n_m_sub<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==1,]$is_female))[1,2]
n_f_sub<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==1,]$is_female))[2,2]
n_m_can<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==0,]$is_female))[1,2]
n_f_can<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==0,]$is_female))[2,2]
n_m<-as.data.frame(table(data$is_female))[1,2]
n_f<-as.data.frame(table(data$is_female))[2,2]
n_r_sub<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==1,]$region_is_urban))[1,2]
n_u_sub<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==1,]$region_is_urban))[2,2]
n_r_can<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==0,]$region_is_urban))[1,2]
n_u_can<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==0,]$region_is_urban))[2,2]
n_r<-as.data.frame(table(data$region_is_urban))[1,2]
n_u<-as.data.frame(table(data$region_is_urban))[2,2]
n_4_sub<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==1,]$birth_cohort_cat3))[1,2]
n_5_sub<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==1,]$birth_cohort_cat3))[2,2]
n_6_sub<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==1,]$birth_cohort_cat3))[3,2]
n_4_can<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==0,]$birth_cohort_cat3))[1,2]
n_5_can<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==0,]$birth_cohort_cat3))[2,2]
n_6_can<-as.data.frame(table(data[data$multiplexserology_cancer_b1_subcohort==0,]$birth_cohort_cat3))[3,2]
n_4<-as.data.frame(table(data$birth_cohort_cat3))[1,2]
n_5<-as.data.frame(table(data$birth_cohort_cat3))[2,2]
n_6<-as.data.frame(table(data$birth_cohort_cat3))[3,2]
n_sub_all<-as.data.frame(table(data$multiplexserology_cancer_b1_subcohort))[2,2] 
n_can_all<-as.data.frame(table(data$multiplexserology_cancer_b1_subcohort))[1,2]


n_male_sub<-format(as.numeric(n_m_sub),big.mark=",",scientific=FALSE)
n_female_sub<-format(as.numeric(n_f_sub),big.mark=",",scientific=FALSE)
n_male_can<-format(as.numeric(n_m_can),big.mark=",",scientific=FALSE)
n_female_can<-format(as.numeric(n_f_can),big.mark=",",scientific=FALSE)
n_rural_sub<-format(as.numeric(n_r_sub),big.mark=",",scientific=FALSE)
n_urban_sub<-format(as.numeric(n_u_sub),big.mark=",",scientific=FALSE)
n_rural_can<-format(as.numeric(n_r_can),big.mark=",",scientific=FALSE)
n_urban_can<-format(as.numeric(n_u_can),big.mark=",",scientific=FALSE)
n_40s_sub<-format(as.numeric(n_4_sub),big.mark=",",scientific=FALSE)
n_50s_sub<-format(as.numeric(n_5_sub),big.mark=",",scientific=FALSE)
n_60s_sub<-format(as.numeric(n_6_sub),big.mark=",",scientific=FALSE)
n_40s_can<-format(as.numeric(n_4_can),big.mark=",",scientific=FALSE)
n_50s_can<-format(as.numeric(n_5_can),big.mark=",",scientific=FALSE)
n_60s_can<-format(as.numeric(n_6_can),big.mark=",",scientific=FALSE)
n_sub_all<-format(as.numeric(n_sub_all),big.mark=",",scientific=FALSE)
n_can_all<-format(as.numeric(n_can_all),big.mark=",",scientific=FALSE)

### Create flextable ####

set_flextable_defaults(font.family = "Arial")

header <- data.frame(col_keys   = col_names,
                     title      = rep("eTable 4: Sero-prevalence proportion (SE) for infectious pathogens for incident cancer, by sex, region and birth cohort",12),
                     col_head_1 = rep(c(NA,"Sex","","Area","","Birth cohort","","All"),c(1,2,1,2,1,3,1,1)),
                     col_head_2 = c("Pathogen",  
                                    paste0("Male \n(n=",n_male_can,")"),paste0("Female \n(n=",n_female_can,")"),NA,
                                    paste0("Rural \n(n=",n_rural_can,")"),paste0("Urban \n(n=",n_urban_can,")"),NA,
                                    paste0("< 1950 \n(n=",n_40s_can,")"),paste0("1950-1959 \n(n=",n_50s_can,")"),paste0("> 1959 \n(n=",n_60s_can,")"),NA,
                                    paste0("All \n(n=",n_can_all,")")))
                                   
ft<-flextable(tab, col_keys=c(col_names)) %>%
  set_header_df(mapping = header, key = "col_keys") %>%
  merge_h(part = "header") %>%
  bold(i=c(1:3),j=NULL,part = "header") %>% 
  align(i=1:3,j=c(2:12),align = "center", part = "header") %>%   
  hline(i=1,j=NULL,border = fp_border(width = 1, color = "black"), part = "header" ) %>%
  hline(i=2,j=c(2:3,5:6,8:10,12),border = fp_border(width = 1, color = "black"), part = "header" ) %>%
  hline(i=3,j=c(2:12),border = fp_border(width = 1, color = "black"), part = "header" ) %>%
  padding(i=c(2:8,10:11,13:14,16:18,20:21,23:26), j=1, padding.left=15) %>%
  width(j=NULL,width=c(rep(c(1.75,1.2,0.1,1.2,0.1,1.2,0.1,1.2),c(1,2,1,2,1,3,1,1)))) %>%
  align(i=NULL,j=c(2:12), align="center", part="body") %>%  
  bold(i=n_sig_sex_sub,j=c(2,3), part="body") %>% 
  bold(i=n_sig_region_sub,j=c(5,6), part="body") %>% 
  bold(i=n_sig_birth_sub,j=c(8:10), part="body") %>% 
  hline_top(border = fp_border(width =1, color = "black"), part = "body" ) %>%  
  hline_bottom(border = fp_border(width =1, color = "black"), part = "body" ) %>%
  italic(i = c(23:26), j = 1, italic = TRUE, part = "body") %>%
  add_footer_lines("Bold values denote statistical significance at the p<0.05 level.") %>% 
  add_footer_lines("BKV, BK polyomavirus; CMV, cytomegalovirus; EBV, Epstein-Barr virus; HBV, hepatitis B virus; HCV, hepatitis C virus;") %>%
  add_footer_lines("HTLV, Human T-cell lymphotropic virus type 1; JCV, JC polyomavirus; MCV, Merkel cell polyomarvirus; VZV, varicella zoster virus.") %>%
  line_spacing(space = 1, part = "header") %>%  
  line_spacing(space = 0.3, part = "body") %>% 
  line_spacing(space = 0.2, part = "footer") %>%
  fontsize(i = NULL, j = NULL, size = 9, part = "header") %>%
  fontsize(i = NULL, j = NULL, size = 9, part = "body") %>%
  fontsize(i = NULL, j = NULL, size = 8, part = "footer") 

# Define page size for output

sect_properties <- prop_section(page_size = page_size(orient = "landscape", width =12, height = 8),
                                type = "continuous", page_margins = page_mar())


save_as_docx(ft,path="K:/kadoorie/Staff_Folders/JonathanC/Projects/serology/Tables_figures/Tables_figures_10/eTable_5.docx",
             pr_section=sect_properties)

###################################################################################################################################################
