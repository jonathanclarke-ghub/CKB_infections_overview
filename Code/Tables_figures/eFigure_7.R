
library(tidyverse)
library(DescTools)
library(plotrix)
library(flextable)
library(officer)
library(stringr)
library(ggplot2)
library(tidyr)

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
tab_meta$pathogen_label[36]<-"C. trachomatis"
tab_meta$pathogen_label[37:38]<-"T. gondii"
tab_meta$pathogen_label[39]<-"C. burnetii"
tab_meta$antigen_label_long<-NA
tab_meta$antigen_label_long[1:38]<-paste0(tab_meta$pathogen_label[1:38]," ",tab_meta$antigen_label[1:38])
tab_meta$antigen_label_long[39:47]<-tab_meta$antigen_label[39:47]

readRDS(paste0(path_1,"multiplex_cancer_case_cohort_meta.rds"))

### Table - descriptive, by sex and overall ####

antigen_pos<-grep("_ap$",names(data),value=TRUE)

data_sub<-data %>% filter(multiplexserology_cancer_b1_subcohort==1) %>% dplyr::select(all_of(antigen_pos))
data_can<-data %>% filter(multiplexserology_cancer_b1_subcohort==0) %>% dplyr::select(all_of(antigen_pos))

names(data_sub)<-str_remove_all(names(data_sub),"_ap")
names(data_can)<-str_remove_all(names(data_can),"_ap")

new_name <- tab_meta$antigen_label_long
old_name <- tab_meta$variable_name

data_sub<-data_sub %>% rename_with(~ new_name, all_of(old_name))
data_can<-data_can %>% rename_with(~ new_name, all_of(old_name))

### Calculate correlation coefficients ####

tab_sub<-as.data.frame(cor(data_sub, method='spearman'))
tab_can<-as.data.frame(cor(data_can, method='spearman'))

row<-list()

for(ii in 2:45){
  
  row[[ii]]<-cbind(tab_can[ii,1:ii],tab_sub[ii,(ii+1):47])
  
}

tab_2_45<-as.data.frame(do.call(rbind,row))

tab_1 <-tab_sub[1,1:47]
tab_46<-cbind(tab_can[46,1:46],tab_sub[46,47])
colnames(tab_46)<-colnames(tab_1)
tab_47<-tab_can[47,1:47]

tab_sub_can<-as.data.frame(rbind(tab_1,tab_2_45,tab_46,tab_47))
tab_sub_can$antigen_sub<-rownames(tab_sub_can)

data_long<-tab_sub_can %>% pivot_longer(cols=colnames(tab_sub_can)[1:47],
                    names_to='antigen_can',
                    values_to="coefficient") 

levels=c("T. gondii sag-1", 
  "T. gondii p22", 
  "CBU-1910 Com1", 
  "HP0887-VacA (N+C)", 
  "HP1564-OMP", 
  "HP0305", 
  "HP-HopA",
  "HP-HcpC",
  "HP0010-GroEL", 
  "HP0875-Catalase", 
  "HP0547-CagA (N+C)", 
  "C. trachomatis pGP3", 
  "HPV-18 L1", 
  "HPV-18 E7", 
  "HPV-18 E6", 
  "HPV-16 L1", 
  "HPV-16 E7", 
  "HPV-16 E6",
  "HPV-16 E2", 
  "HPV-16 E1", 
  "MCV VP1", 
  "JCV VP1", 
  "BKV VP1", 
  "HTLV-1 gag", 
  "HTLV-1 env", 
  "HIV-1 gag",
  "HIV-1 env", 
  "HCV Core", 
  "HCV NS3", 
  "HBV HBe", 
  "HBV HBc", 
  "HHV-7 U14", 
  "HHV-6 IE1B",
  "HHV-6 IE1A", 
  "CMV pp 52",
  "CMV pp 28", 
  "CMV pp150 Nter", 
  "EBV BXLF1",
  "EBV BGLF2", 
  "EBV BFRF1",
  "EBV ZEBRA", 
  "EBV VCA p18", 
  "EBV EA-D", 
  "EBV EBNA-1", 
  "VZV gE/gI", 
  "HSV-2 mgG unique",
  "HSV-1 gg")


data_long$antigen_sub<-factor(data_long$antigen_sub,levels)
data_long$antigen_can<-factor(data_long$antigen_can,levels)

### create Heatmap ####

plot<-ggplot(data_long, aes(antigen_sub, antigen_can)) + 
  geom_tile(aes(fill= coefficient)) + 
  scale_fill_gradient2(low="blue",mid="white", high="red",limits=c(-1,1)) +
  geom_text(aes(label = ifelse(coefficient>-0.01 & coefficient<0.01,"-",round(coefficient, 2))),size=4) +
  labs(
    #title="Correlation between antigen seropositivity by study arms",
       x="Subcohort",
       y="Incident Cancer",
       fill ="Coefficient") +
  theme(plot.title   = element_text(hjust = 0.5,vjust=5,size = 20,face="bold"),
        legend.title = element_text(size=14,vjust = 5),
        legend.text  = element_text(size=14),
        axis.title.x = element_text(size = 18,face="bold"),
        axis.title.y = element_text(size = 18, face="bold"),
        axis.text.x  = element_text(vjust=0.9,hjust=1,angle=45,size=14,face="bold"),
        axis.text.y  = element_text(size=12,face="bold"),
        plot.margin  = unit(c(2,1,1,1), 'cm'))
      
ggsave(paste0(path_3,"eFigure_73.png"),
       plot,
       dpi = 900,
       width = 60,   #50
       height = 62, #52  
       units = "cm")

################################################################################################################################################





