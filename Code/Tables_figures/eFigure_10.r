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
library(Epi)

### Set paths to data and results ####

path_1<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Raw_data_7/"
path_2<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Project_data_7/"
path_3<-"K:/kadoorie/Staff_Folders/JonathanC/Projects/Serology/Tables_figures/Tables_figures_10/"

### Import participant data ####

data<-readRDS(paste0(path_2,"data_serology_2024_7.rds"))

### Exposure ####

table(data$pathogen_index,useNA="ifany")

table(data$pathogen_index,data$ep_cancer_all_combined_ep)

#data$pathogen_index<-factor(ifelse(data$pathogen_index<6,6,ifelse(data$pathogen_index>13,13,data$pathogen_index)))
data$pathogen_index<-factor(ifelse(data$pathogen_index<5,5,ifelse(data$pathogen_index>15,15,data$pathogen_index)))

table(data$pathogen_index,useNA="ifany")

#data$pathogen_index<-factor(data$pathogen_index,levels=c(9,10,11,12,13,6,7,8))
#data$pathogen_index<-factor(data$pathogen_index,levels=c(10,11,12,13,6,7,8,9))
data$pathogen_index<-factor(data$pathogen_index,levels=c(10,11,12,13,14,15,5,6,7,8,9))

exposure<-"pathogen_index"

### Outcome ####

outcome<-"ep_cancer_all_combined"

### Adjustments ####

adjust_for_1<-paste0("+ age + as.factor(is_female) + as.factor(region_code) + as.factor(highest_education)")
adjust_for_2<-paste0(adjust_for_1, "+ bmi_calc + as.factor(ever_reg_smoker) + as.factor(ever_reg_alcohol) + as.factor(family_cancer)")

adjust_for<-adjust_for_2

### Coxmodel - prentice weighted ####

data$outcome<-data[[paste0(outcome,"_ep")]]
data$outcome_date<-data[[paste0(outcome,"_datedeveloped")]]

data$study_date<-as.numeric(as.Date(data$study_date))  
data$outcome_date<-as.numeric(data$outcome_date)   

### Keep subcohort participants and cases with site specific cancer (outcome) of interest ####

data<-data[data$multiplexserology_cancer_b1_subcohort==1  | data$outcome==1,]

data$study_time<-data$outcome_date-data$study_date
data$study_date<-0
data$study_date[data$outcome==1 & data$multiplexserology_cancer_b1_subcohort==0]<-data$study_time[data$outcome==1 & data$multiplexserology_cancer_b1_subcohort==0]-0.001

model<-coxph(as.formula(paste0("Surv(study_date,study_time,",outcome,"_ep) ~ ",exposure,"+",adjust_for)),data)

results<-data.frame("estimate"          = (float(model)$coef),
                     "stderr"       = sqrt(float(model)$var))


results$index     <- rownames(results)

results$n         <- as.data.frame(table(data$pathogen_index))[,2]
results$n<-format(as.numeric(results$n),big.mark = ",",scientific=FALSE)

results$estimate  <-as.numeric(results$estimate)
results$stderr     <-as.numeric(results$stderr)
results$ci_lower<-exp(results$estimate-1.96*results$stderr)
results$ci_upper<-exp(results$estimate+1.96*results$stderr)

#results$index_cat <-factor(c("9","10","11","12","\u226513","\u22646","7","8"),
#                    levels = c("\u22646","7","8","9","10","11","12","\u226513"))
#results$index_cat <-factor(c("10","11","12","\u226513","\u22646","7","8","9"),
#                          levels = c("\u22646","7","8","9","10","11","12","\u226513"))
results$index_cat <-factor(c("10","11","12","13","14","\u226515","\u22645","6","7","8","9"),
                           levels = c("\u22645","6","7","8","9","10","11","12","13","14","\u226515"))

results$HR<-exp(results$estimate)

write.csv(results,paste0(path_3,"eFig10_data.csv"))

## Shape plot ####

plot<-shape_plot(results,
                       col.x        = "index_cat",
                       col.estimate = "estimate",
                       col.stderr   = "stderr",
                       col.n        = "n", 
                       xlims        = c(0.5,12.5),
                       ylims        = c(0.42,1.26),
                       ybreaks      = c(0.5,0.75,1,1.25),
                       exponentiate = TRUE,
                       scalepoints  = TRUE,
                       ciunder      = TRUE,
                       base_size    = 15,
                       xlab         = "No. of co-infections",
                       ylab         = "HR (95% CI)",
                       width       = unit(12.5, "cm"))


title <- textGrob(substitute(paste(bold("eFigure 11: Adjusted HRs for liver cancer from HBV and HCV and for gastric cancers from "),bolditalic("H. pylori"), bold(" in Chinese adults"))),
                  gp = gpar(cex = 0.8 ,fontface="bold"), 
                  x = unit(0, "npc") + unit(1, "lines"),
                  y = unit(0, "npc") + unit(1, "lines"),
                  just = c(-0.04,-1))

footnote <- textGrob("Adjusted for age, sex, region, education, smoking, alcohol, BMI, family history of cancer",
                     gp = gpar(cex = 0.8), 
                     x = unit(0, "npc") + unit(1, "lines"),
                     y = unit(0, "npc") + unit(1, "lines"),
                     just = c(-0.05,14))



figure<-grid.arrange(plot$plot,
                     footnote,
                     nrow = 2,
                     heights = unit.c(unit(165, "mm"),
                                      unit(1, "mm")))


ggsave(paste0(path_3,"eFigure_10.png"),
       figure,
       bg="white",
       dpi = 450,
       width = 18,                      
       height = 24,                     
       units = "cm")

