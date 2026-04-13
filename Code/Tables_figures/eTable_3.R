
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
tab_meta$pathogen_label[39]<-"C. burnetti"
tab_meta$pathogen_label[40:47]<-"H. pylori"

tab_meta$antigen_label_long<-paste0(tab_meta$pathogen_label," ",tab_meta$antigen_label)

tab_meta[tab_meta$variable_name=="hpv16_e6","proposed_cutoff"]<-1000      # Updated cut-off

tab_meta_selected<-tab_meta %>% dplyr::select(variable_name,pathogen_label,antigen,proposed_cutoff,algorithm_seropositivity) 
                                     
### Table - descriptive, by sex and overall ####

antigen_pos<-grep("_ap$",names(data),value=TRUE)

results<-list()

for(ii in antigen_pos){

ap_sub<-table(data[[ii]],data$multiplexserology_cancer_b1_subcohort)[2,2]
ap_can<-table(data[[ii]],data$multiplexserology_cancer_b1_subcohort)[2,1]

results[[ii]]<-cbind(ii,ap_sub,ap_can)

}

ap_all<-as.data.frame(do.call(rbind,results))
names(ap_all)[1]<-"variable_name"

ap_all$ap_sub1<-as.numeric(ap_all$ap_sub)/dim(data[data$multiplexserology_cancer_b1_subcohort==1,])[1]*100
ap_all$ap_can1<-as.numeric(ap_all$ap_can)/dim(data[data$multiplexserology_cancer_b1_subcohort==0,])[1]*100
ap_all$ap_sub2<-sprintf(as.numeric(ap_all$ap_sub1), fmt = '%#.1f') 
ap_all$ap_can2<-sprintf(as.numeric(ap_all$ap_can1), fmt = '%#.1f') 
ap_all$variable_name<-str_replace_all(ap_all$variable_name,"_ap","")

tab_meta_selected<-merge(tab_meta_selected,ap_all,by="variable_name")
tab_meta_selected<-tab_meta_selected %>% dplyr::select(pathogen_label,antigen,proposed_cutoff,algorithm_seropositivity,ap_sub2,ap_can2)
tab_meta_selected$pathogen_label[2]<-"C. Burnetti"
tab_meta_selected$pathogen_label[31:38]<-"H. Pylori"
tab_meta_selected$algorithm_seropositivity[7:9]<-"above cutoff"
tab_meta_selected$algorithm_seropositivity[31:38]<-"\u2265 4 out of 8 above cutoff*"
tab_meta_selected$algorithm_seropositivity[c(14:17,21:22,41:42,45:46)]<-"2 out of 2 above cutoff"
tab_meta_selected$algorithm_seropositivity[18:19]<-"\u2265 1 out of 2 above cutoff"
tab_meta_selected$algorithm_seropositivity[3:5]<-"\u2265 2 out of 3 above cutoff"
tab_meta_selected$algorithm_seropositivity[10:13]<-"\u2265 2 out of 4 above cutoff"
tab_meta_selected$algorithm_seropositivity[c(7:9,23:26,28:29)]<-"-"

sub_1<-rep(c("Hepatitis virus",NA),c(1,5))
sub_2<-rep(c("Human Papillomavirus",NA),c(1,5))
sub_3<-rep(c("Human Herpes Virus",NA),c(1,5))
sub_4<-rep(c("Human Polyomavirus",NA),c(1,5))       
sub_5<-rep(c("Human Retrovirus",NA),c(1,5))   
sub_6<-rep(c("Bacteria and Parasite",NA),c(1,5))

tab<-as.data.frame(rbind(sub_1,sub_2,sub_3,sub_4,sub_5,sub_6))
names(tab)<-names(tab_meta_selected)
tab<-rbind(tab_meta_selected,tab)
tab$n<-1:dim(tab)[1]

tab<-tab[c(50,39,40,47,7:13,3:5,18:20,
           48,14:17,
           52,21:22,41:42,  
           51,1,43,44,
           49,23:30,
           53,6,31:38,2,45:46),]

tab<-tab[,c(1:2,5:6)]
col_names<-names(tab)

### Create flextable ####

set_flextable_defaults(font.family = "Arial")

header <- data.frame(col_keys = col_names,
                     head_1  = c("Pathogen","Antigen","Seropositivity","Seropositivity"),
                     head_2  = c("Pathogen","Antigen","Subcohort","Incident Cancer"))

ft<-flextable(tab, col_keys=c(col_names)) %>%
  set_header_df(mapping = header, key = "col_keys") %>%
  merge_h(part = "header") %>%
  merge_v(part = "header") %>%
  merge_v(j=1,part = "body") %>%
  merge_h(i=c(1,18,21,24,27,31)) %>%
  align(i=c(1,18,21,24,27,31),j=NULL,align = "left", part = "body") %>% 
  bold(i=NULL,j=NULL,part = "header") %>% 
  align(i=NULL,j=2:4,align = "center", part = "header") %>%  
  align(i=NULL,j=2:4,align = "center", part = "body") %>% 
  valign(i=NULL,j=3:4,valign = "bottom", part = "header") %>%  
  valign(i=NULL,j=1,valign = "top", part = "body") %>% 
 # valign(i=NULL,j=2:4,valign = "bottom", part = "body") %>% 
  hline_top(border = fp_border(width =1, color = "black"), part = "header" ) %>%  
  hline(i=c(1),j=c(2:3),border = fp_border(width = 1, color = "black"), part = "header" ) %>%
  hline(i=c(2),j=c(2),border = fp_border(width = 1, color = "black"), part = "header" ) %>%
  hline_bottom(border = fp_border(width =1, color = "black"), part = "header") %>%
  padding(padding.left = 4,padding.right = 4,padding.top = 0.8,padding.bottom = 0.8, part = "all") %>% 
  padding(i=c(2:17,19:22,24:27,29:31,33:39,41:53), j=1, padding.left=15) %>%
  width(j=c(1:4),width=c(2,1,1,1)) %>%
  align(i=NULL,j=c(4), align="center", part="body") %>%  
  bold(i=c(1,18,23,28,32,34,36,40,41),j=1, part="body") %>%  
  italic(i=c(42:53),j=1, part="body") %>%  
  fontsize(i = NULL, j = NULL, size = 10, part = "header") %>%
  fontsize(i = NULL, j = NULL, size = 10, part = "body")


# Define page size for output

sect_properties <- prop_section(page_size = page_size(orient = "portrait", width = 8, height = 12),
                                type = "continuous", page_margins = page_mar())

save_as_docx(ft,path="K:/kadoorie/Staff_Folders/JonathanC/Projects/Serology/Tables_figures/Tables_figures_10/eTable_3t.docx",
             pr_section=sect_properties)


################################################################################################################################################
