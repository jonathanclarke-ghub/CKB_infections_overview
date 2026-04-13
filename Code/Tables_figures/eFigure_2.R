
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

tab_all$pathogen_label<-factor(tab_all$pathogen_label,levels=c("T.gondii","C.burnetii","H.pylori","C.trachomatis",
                                                                        "HPV-18","HPV-16","MCV","JCV","BKV","HIV-1","HTLV-1",
                                                                        "HCV","HBV","HHV-7","HHV-6", "CMV","EBV","VZV","HSV-2","HSV-1"))


tab_m<-tab_all[,c("pathogen_label","prop_m","stderr_m","p_value_sex_all")]
names(tab_m)<-c("pathogen","prop","stderr","p_value")
tab_m$sub_group<-"Male"
tab_m$group<-"Sex"
tab_f<-tab_all[,c("pathogen_label","prop_f","stderr_f","p_value_sex_all")]
names(tab_f)<-c("pathogen","prop","stderr","p_value")
tab_f$sub_group<-"Female"
tab_f$group<-"Sex"
tab_r<-tab_all[,c("pathogen_label","prop_r","stderr_r","p_value_region_all")]
names(tab_r)<-c("pathogen","prop","stderr","p_value")
tab_r$sub_group<-"Rural"
tab_r$group<-"Region"
tab_u<-tab_all[,c("pathogen_label","prop_u","stderr_u","p_value_region_all")]
names(tab_u)<-c("pathogen","prop","stderr","p_value")
tab_u$sub_group<-"Urban"
tab_u$group<-"Region"
tab_40s<-tab_all[,c("pathogen_label","prop_40s","stderr_40s","p_value_birth_all")]
names(tab_40s)<-c("pathogen","prop","stderr","p_value")
tab_40s$sub_group<-"< 1950"
tab_40s$group<-"Birth cohort"
tab_50s<-tab_all[,c("pathogen_label","prop_50s","stderr_50s","p_value_birth_all")]
names(tab_50s)<-c("pathogen","prop","stderr","p_value")
tab_50s$sub_group<-"1950-1959"
tab_50s$group<-"Birth cohort"
tab_60s<-tab_all[,c("pathogen_label","prop_60s","stderr_60s","p_value_birth_all")]
names(tab_60s)<-c("pathogen","prop","stderr","p_value")
tab_60s$sub_group<-"> 1959"
tab_60s$group<-"Birth cohort"
tab_sub<-tab_all[,c("pathogen_label","prop_sub_all","stderr_sub_all","p_value_all")]
names(tab_sub)<-c("pathogen","prop","stderr","p_value")
tab_sub$sub_group<-"Subcohort"
tab_sub$group<-"All"
tab_can<-tab_all[,c("pathogen_label","prop_can_all","stderr_can_all","p_value_all")]
names(tab_can)<-c("pathogen","prop","stderr","p_value")
tab_can$sub_group<-"Incident Cancer"
tab_can$group<-"All"

tab_sum<-rbind(tab_m,tab_f,tab_r,tab_u,tab_40s,tab_50s,tab_60s,tab_sub,tab_can)
tab_sum$group<-factor(tab_sum$group,levels=c("Sex","Region","Birth cohort","All"))
tab_sum$sub_group<-factor(tab_sum$sub_group,levels=c("Male","Female","Rural","Urban","< 1950","1950-1959","> 1959","Incident Cancer","Subcohort"))
tab_sum$prop_1<-sprintf(as.numeric(tab_sum$prop), fmt = '%#.1f') 
tab_sum$adj<-1.7
tab_sum$adj[tab_sum$sub_group %in% c("Subcohort","Female","Urban")]<--0.6
tab_sum$adj[tab_sum$sub_group %in% c("< 1950")]<-2.2
tab_sum$adj[tab_sum$sub_group %in% c("1950-1959")]<-0.4
tab_sum$adj[tab_sum$sub_group %in% c("> 1959")]<--1.5

tab_sum$prop_1[tab_sum$group=="All" & tab_sum$pathogen %in% sig_all]<-paste0(tab_sum$prop_1[tab_sum$group=="All" & tab_sum$pathogen %in% sig_all],"*")

# Create forest plot

plot_all<-ggplot(data=tab_sum[tab_sum$group=="All",], aes(x=pathogen, y=as.numeric(prop),fill=sub_group)) +
  geom_bar(stat = "identity", position = "dodge",colour="black") +
  scale_fill_manual(breaks=c("Subcohort","Incident Cancer"),values= c("grey","white")) +
  scale_y_continuous(expand=c(0,0),limits=c(0,110),breaks=seq(0, 100, 25),labels=seq(0, 100, 25)) +
  xlab("Pathogen") +
  ylab("Seroprevalence, %") +  
  theme_classic() +
  coord_flip() + 
  geom_text(data=tab_sum[tab_sum$sub_group=="Subcohort",],aes(label=prop_1,vjust=adj[1]),hjust=-0.3,size = 6) + 
  geom_text(data=tab_sum[tab_sum$sub_group=="Incident Cancer",],aes(label=prop_1,vjust=adj[1]),hjust=-0.3,size = 6) + 
  theme(axis.text.y  = element_text(face=rep(c("bold.italic","bold"),c(4,16)),size=18),
        axis.text.x  = element_text(face="bold",size=18),
        axis.title.x = element_text(face="bold",vjust=-1,size=18),
        axis.title.y = element_text(vjust=-7.5,face="bold",size=18),
        axis.line.y  = element_blank(),
        legend.title = element_blank(),
        legend.text  = element_text(size=18),
        legend.position = c(0.9,0.05),
        legend.key.size = unit(1.5, units = "cm"),
        plot.title  = element_text(hjust=0,vjust=-7.5,face="bold",size=28))


# Panel of plots

panel_of_plots<-ggarrange(plot_all,ncol=1)
 
title <- textGrob("eFigure 3: Sero-prevalence for infectious pathogens, by study arm and sex, region and birth cohort",
                  gp = gpar(cex = 2 ,fontface="bold"), 
                  x = unit(0, "npc") + unit(1, "lines"),
                  y = unit(0, "npc") + unit(1, "lines"),
                  just = c(0,0))

footnote <- textGrob("",
                     gp = gpar(cex = 1), 
                     x = unit(0, "npc") + unit(1, "lines"),
                     y = unit(0, "npc") + unit(1, "lines"),
                     just = c(0,1.5))

figure<-grid.arrange(panel_of_plots,
                     footnote,
                     nrow = 2,
                     heights = unit.c(unit(500, "mm"),
                                      unit(1, "mm")))

ggsave(paste0(path_3,"eFigure_2.png"),
       figure,
       dpi = 450,
       width = 40,
       height = 58,
       units = "cm")


