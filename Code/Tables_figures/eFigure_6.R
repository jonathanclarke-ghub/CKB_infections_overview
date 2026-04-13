
library(tidyverse)
library(DescTools)
library(plotrix)
library(flextable)
library(officer)
library(stringr)
library(ggplot2)
library(tidyr)
library(dplyr)

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

data_sub<-data %>% filter(multiplexserology_cancer_b1_subcohort==1) %>% dplyr::select(all_of(pathogen_pos))
data_can<-data %>% filter(multiplexserology_cancer_b1_subcohort==0) %>% dplyr::select(all_of(pathogen_pos))

names(data_sub)<-str_remove_all(names(data_sub),"_pp")
names(data_can)<-str_remove_all(names(data_can),"_pp")

new_name <- tab_meta[!duplicated(tab_meta$pathogen),c("pathogen_label")]
old_name <- tab_meta[!duplicated(tab_meta$pathogen),c("pathogen")]

data_sub<-data_sub %>% rename_with(~ new_name$pathogen_label, all_of(old_name$pathogen))
data_can<-data_can %>% rename_with(~ new_name$pathogen_label, all_of(old_name$pathogen))

### Calculate correlation coefficients ####

tab_sub<-as.data.frame(cor(data_sub, method='spearman'))
tab_can<-as.data.frame(cor(data_can, method='spearman'))

row<-list()

for(ii in 2:18){
  
  row[[ii]]<-cbind(tab_can[ii,1:ii],tab_sub[ii,(ii+1):20])
  
}

tab_2_18<-as.data.frame(do.call(rbind,row))

tab_1 <-tab_sub[1,1:20]
tab_19<-cbind(tab_can[19,1:19],tab_sub[19,20])
colnames(tab_19)<-colnames(tab_1)
tab_20<-tab_can[20,1:20]

tab_sub_can<-as.data.frame(rbind(tab_1,tab_2_18,tab_19,tab_20))
tab_sub_can$pathogen_sub<-rownames(tab_sub_can)

data_long<-tab_sub_can %>% pivot_longer(cols=colnames(tab_sub_can)[1:20],
                    names_to='pathogen_can',
                    values_to="coefficient") 

data_long$pathogen_sub<-factor(data_long$pathogen_sub,levels=c("T.gondii","C.burnetii","H.pylori","C.trachomatis",
                                                              "HPV-18","HPV-16","MCV","JCV","BKV","HIV-1","HTLV-1",
                                                              "HCV","HBV","HHV-7","HHV-6", "CMV","EBV","VZV","HSV-2","HSV-1"))

data_long$pathogen_can<-factor(data_long$pathogen_can,levels=c("T.gondii","C.burnetii","H.pylori","C.trachomatis",
                                                               "HPV-18","HPV-16","MCV","JCV","BKV","HIV-1","HTLV-1",
                                                               "HCV","HBV","HHV-7","HHV-6", "CMV","EBV","VZV","HSV-2","HSV-1"))
                                                               
### create Heatmap ####

plot<-ggplot(data_long, aes(pathogen_sub, pathogen_can, fill= coefficient)) + 
  geom_tile() + 
  scale_fill_gradient2(low="blue",mid="white", high="red",limits=c(-1,1)) +
  geom_text(aes(label = ifelse(coefficient>-0.01 & coefficient<0.01,"<0.01",round(coefficient, 2))),size = 5,face="bold") +
  labs(
    #title="Correlation between pathogen seropositivity by study arms",
       x="Subcohort",
       y="Incident Cancer",
       fill ="Coefficient") +
  theme(plot.title   = element_text(hjust = 0.5,vjust=5,size = 20,face="bold"),
        legend.title = element_text(size=18,vjust = 5),
        legend.text  = element_text(size=16),
        axis.title.x = element_text(size = 20,face="bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        axis.text.x  = element_text(vjust=0.9,hjust=1,angle=45,size=16,face=rep(c("bold.italic","bold"),c(4,16))),
        axis.text.y  = element_text(size=16,face=rep(c("bold.italic","bold"),c(4,16))),
        plot.margin  = unit(c(2,1,1,1), 'cm'))
      
ggsave(paste0(path_3,"eFigure_6.png"),
       plot,
       dpi = 900,
       width = 40,
       height = 42,
       units = "cm")

################################################################################################################################################

