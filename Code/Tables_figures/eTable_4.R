
### Libraries ###

library(stringr)
library(tidyverse)
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

## Set paths ####

path_1<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Raw_data_7/"
path_2<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Project_data_7/"
path_3<-"K:/kadoorie/Staff_Folders/JonathanC/Projects/Serology/Tables_figures/Tables_figures_10/"

table(data$hep_b_diag,useNA="always")

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

coinfection_fun<-function(var,cat){

  if(var %in% c("ever_reg_alcohol","ever_reg_smoker")) {
    mean_n<-mean(data$pathogen_index[data[[var]]==cat & data$is_female==0],na.rm=TRUE)
    mean_c<-sprintf(mean_n,fmt = '%#.2f')
    sd<-sprintf(sd(data$pathogen_index[data[[var]]==cat & data$is_female==0]),fmt = '%#.2f')
    se<-std.error(data$pathogen_index[data[[var]]==cat & data$is_female==0],na.rm=TRUE)
    mean_sd<-paste0(mean_c," (",sd,")")
    rbind(var,cat,mean_sd,mean_n,se)
    
  }
  
  else{
    mean_n<-mean(data$pathogen_index[data[[var]]==cat],na.rm=TRUE)
    mean_c<-sprintf(mean_n,fmt = '%#.2f')
    sd<-sprintf(sd(data$pathogen_index[data[[var]]==cat],na.rm=TRUE),fmt = '%#.2f')
    se<-std.error(data$pathogen_index[data[[var]]==cat],na.rm=TRUE)
    mean_sd<-paste0(mean_c," (",sd,")")
    rbind(var,cat,mean_sd,mean_n,se)

  }
}

process_data<-function(){

tab<-as.data.frame(rbind(  
  
c("Age (years)","< 40",                          t(coinfection_fun("age_cat5",0))),
c("Age (years)","40 to 49",                      t(coinfection_fun("age_cat5",1))),
c("Age (years)","50 to 59",                      t(coinfection_fun("age_cat5",2))),
c("Age (years)","60 to 69",                      t(coinfection_fun("age_cat5",3))),
c("Age (years)","70 +",                          t(coinfection_fun("age_cat5",4))),
c("Sex","Women",                                 t(coinfection_fun("is_female",1))),
c("Sex","Men",                                   t(coinfection_fun("is_female",0))),
c("Area","Urban",                                t(coinfection_fun("region_is_urban",1))),
c("Area","Rural",                                t(coinfection_fun("region_is_urban",0))),
c("Schooling","<= 6 years",                     t(coinfection_fun("education_cat2",0))),
c("Schooling","> 6 years",                      t(coinfection_fun("education_cat2",1))),
c("Ever regular smoker (men)","Yes",            t(coinfection_fun("ever_reg_smoker",1))),
c("Ever regular smoker (men)","No",             t(coinfection_fun("ever_reg_smoker",0))),
c("Ever regular alcohol drinker (men)","Yes",   t(coinfection_fun("ever_reg_alcohol",1))),
c("Ever regular alcohol drinker (men)","No",    t(coinfection_fun("ever_reg_alcohol",0))),
c("BMI, kg/m\u00B2","Normal, or below",         t(coinfection_fun("bmi_cat3",0))),
c("BMI, kg/m\u00B2","Overweight",               t(coinfection_fun("bmi_cat3",1))),
c("BMI, kg/m\u00B2","Obese",                    t(coinfection_fun("bmi_cat3",2))),
c("Blood transfusion","Yes",                    t(coinfection_fun("blood_transfusion",1))),
c("Blood transfusion","No",                     t(coinfection_fun("blood_transfusion",0))),
c("HBsAg+","Yes",                               t(coinfection_fun("hep_b_diag",1))),
c("HBsAg+","No",                                t(coinfection_fun("hep_b_diag",0))),
c("CHD or Stroke","Yes",                        t(coinfection_fun("chd_stroke_diag",1))),
c("CHD or Stroke","No",                         t(coinfection_fun("chd_stroke_diag",0))),
c("Diabetes","Yes",                             t(coinfection_fun("diabetes_diag",1))),
c("Diabetes","No",                              t(coinfection_fun("diabetes_diag",0))),
c("Cirrhosis/Hepatitis","Yes",                  t(coinfection_fun("cirrhosis_hep_diag",1))),
c("Cirrhosis/Hepatitis","No",                   t(coinfection_fun("cirrhosis_hep_diag",0))),
c("Emphysema/Bronchitis","Yes",                 t(coinfection_fun("emph_bronc_diag",1))),
c("Emphysema/Bronchitis","No",                  t(coinfection_fun("emph_bronc_diag",0))),
c("Tuberculosis","Yes",                         t(coinfection_fun("tb_diag",1))),
c("Tuberculosis","No",                          t(coinfection_fun("tb_diag",0))),
c("Peptic Ulcer","Yes",                         t(coinfection_fun("peptic_ulcer_diag",1))),
c("Peptic Ulcer","No",                          t(coinfection_fun("peptic_ulcer_diag",0))),
c("Self-rated health","Poor",                   t(coinfection_fun("srh_poor",1))),
c("Self-rated health","Fair, good or excellent",t(coinfection_fun("srh_poor",0))),
c("Family history of cancer","Yes",             t(coinfection_fun("family_cancer",1))),
c("Family history of cancer","No",              t(coinfection_fun("family_cancer",0)))))

names(tab)<-c("varname","varcat","var","varval","mn_sd","mn","se")

return(tab)

}

# Subcohort

data<-readRDS(paste0(path_2,"data_serology_2024_7.rds"))
data<-data %>% filter(multiplexserology_cancer_b1_subcohort==1)
tab<-process_data()


pp_list<-list()
for(ii in tab$varname) {
  
  beta<-as.numeric(tab[tab$varname==ii,"mn"])
  se<-as.numeric(tab[tab$varname==ii,"se"])
  het <- heterogeneity(beta, se)
  pp_list[[ii]]<-het$p
}

pp<-as.data.frame(do.call(rbind,pp_list))
sig_sub<-rownames(pp)[pp$V1 < 0.05]

tab<-tab[,c(1,2,5)]
tab$n<-1:dim(tab)[1]
tab_sub<-tab

# Pre-cancer

data<-readRDS(paste0(path_2,"data_serology_2024_7.rds"))
data<-data %>% filter(multiplexserology_cancer_b1_subcohort==0)
tab<-process_data()

pp_list<-list()
for(ii in tab$varname) {
  
  beta<-as.numeric(tab[tab$varname==ii,"mn"])
  se<-as.numeric(tab[tab$varname==ii,"se"])
  het <- heterogeneity(beta, se)
  pp_list[[ii]]<-het$p
}

pp<-as.data.frame(do.call(rbind,pp_list))
sig_can<-rownames(pp)[pp$V1 < 0.05]

tab<-as.data.frame(tab[,c(5)])
tab$n<-1:dim(tab)[1]
tab_can<-tab

# Merge

tab<-merge(tab_sub,tab_can,by="n")

sig_sub_n<-tab$n[tab$varname %in% sig_sub]
sig_can_n<-tab$n[tab$varname %in% sig_can]

tab<-tab[,2:5]

col_names<-names(tab)

### Create flextable ####

set_flextable_defaults(font.family = "Arial")

header <- data.frame(col_keys = col_names,
                     head_1  = c("Baseline characteristic",NA,"Subcohort","Incident Cancer"))

ft<-flextable(tab, col_keys=c(col_names)) %>%
  set_header_df(mapping = header, key = "col_keys") %>%
  merge_h(part = "header") %>%
  bold(i=c(1),j=1:4,part = "header") %>% 
  align(i=1,j=2:4,align = "center", part = "header") %>%  
  padding(padding.left = 4,padding.right = 4,padding.top = 2.25,padding.bottom = 2.25, part = "all") %>% 
  align(i=1,j=1,align = "left", part = "header") %>%   
  hline_top(border = fp_border(width =1, color = "black"), part = "header" ) %>%  
  hline(i=1,j=NULL,border = fp_border(width = 1, color = "black"), part = "header" ) %>%
  merge_v(part = "body") %>%
  valign(i=NULL,j=1,valign = "top", part = "body") %>% 
  width(j=c(1:4),width=c(3,1.8,1,1)) %>%
  align(i=NULL,j=3:4, align="right", part="body") %>%  
  bold(i=sig_sub_n,j=3,part="body") %>%
  bold(i=sig_can_n,j=4,part="body") %>% 
  add_footer_lines("Bold values denote statistically significant difference at the p<0.05 level.") %>%
  fontsize(i = NULL, j = NULL, size = 10, part = "footer") 

  # Define page size for output

sect_properties <- prop_section(page_size = page_size(orient = "portrait", width = 13, height = 10),
                                type = "continuous", page_margins = page_mar())

save_as_docx(ft,path="K:/kadoorie/Staff_Folders/JonathanC/Projects/Serology/Tables_figures/Tables_figures_10/eTable_4.docx",
             pr_section=sect_properties)


