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
path_3<-"K:/kadoorie/Staff_Folders/JonathanC/Projects/Serology/Tables_figures/Tables_figures_10/"

tab_all<-readRDS(paste0(path_2,"tab_cancer_all_cox.rds"))

tab_pathogen_pos<-readRDS(paste0(path_2,"tab_pathogen_pos.rds"))
tab_antigen_pos<-readRDS(paste0(path_2,"tab_antigen_pos.rds"))

tab_pos<-rbind(tab_antigen_pos,tab_pathogen_pos)

### Pathogens and all cancer #### 

tab_cancer_all<-tab_all %>% filter(cancer_code=="ep_cancer_all_combined", analysis=="cancer_all_cox")
tab_cancer_all<-merge(tab_cancer_all,tab_pos,by="pathogen")

tab_cancer_all<-tab_cancer_all %>% filter(pathogen %in% c("hsv1","hsv1_gg","hsv2","hsv2_mggunique","vzv","vzv_gegi","cmv","cmv_pp150n",
                                       "cmv_pp28","cmv_pp52","hhv6","hhv6_ie1a","hhv6_ie1b","hhv7","hhv7_u14","hbv","hbv_hbc","hbv_hbe","hcv","hcv_core","hcv_ns3","bkv","bk_vp1","jcv","jc_vp1","mcv","mc_vp1",
                                       "hpv16","hpv16_e1","hpv16_e2","hpv16_e6","hpv16_e7","hpv16_l1","hpv18","hpv18_e6","hpv18_e7","hpv18_l1","ctrachomatis","ctrachomatis_pgp3","hpylori","hpylori_cagan",
                                       "hpylori_catalase","hpylori_groel","hpylori_hcpc","hpylori_hopa","hpylori_hp0305","hpylori_hp1564","hpylori_vacac","cburnetii","cburnetii_com1",
                                       "tgondii","tgondii_p22trunc","tgondii_sag1d1"))

tab_cancer_all<-tab_cancer_all %>%
  mutate(estimate = ifelse(pathogen %in% c("hsv1","hsv2","vzv","hhv7","bkv","jcv","mcv","hpv16","hpv18","ctrachomatis","cburnetii"),NA,estimate),
         stderr.x = ifelse(pathogen %in% c("hsv1","hsv2","vzv","hhv7","bkv","jcv","mcv","hpv16","hpv18","ctrachomatis","cburnetii"),NA,stderr.x),
         prev     = ifelse(pathogen %in% c("hsv1","hsv2","vzv","hhv7","bkv","jcv","mcv","hpv16","hpv18","ctrachomatis","cburnetii"),NA,prev))

row_labels<-tab_cancer_all[,c("pathogen","label","pathogen_label","antigen_label_long")]
row_labels<-row_labels %>% arrange(factor(pathogen, 
                            levels = c("hsv1","hsv1_gg","hsv2","hsv2_mggunique","vzv","vzv_gegi","cmv","cmv_pp150n",
                                       "cmv_pp28","cmv_pp52","hhv6","hhv6_ie1a","hhv6_ie1b","hhv7","hhv7_u14","hbv","hbv_hbc","hbv_hbe","hcv","hcv_core","hcv_ns3","bkv","bk_vp1","jcv","jc_vp1","mcv","mc_vp1",
                                       "hpv16","hpv16_e1","hpv16_e2","hpv16_e6","hpv16_e7","hpv16_l1","hpv18","hpv18_e6","hpv18_e7", "hpv18_l1","ctrachomatis","ctrachomatis_pgp3","hpylori","hpylori_cagan",
                                       "hpylori_catalase","hpylori_groel","hpylori_hcpc","hpylori_hopa","hpylori_hp0305","hpylori_hp1564","hpylori_vacac","cburnetii","cburnetii_com1",
                                       "tgondii","tgondii_p22trunc","tgondii_sag1d1")))


row_labels$label[c(2,4,6,8:10,12:13,15,17:18,20:21,23,25,27,29:33,35:37,39,41:48,50,52:53)] <- glue::glue("<span style='color:transparent;'>XX</span>{row_labels$label[c(2,4,6,8:10,12:13,15,17:18,20:21,23,25,27,29:33,35:37,39,41:48,50,52:53)]}")
row_labels$label[c(1,3,5,7,11,16,19,40,51)]<- paste0("<b>", row_labels$label[c(1,3,5,7,11,16,19,40,51)], "<b>")

row_labels<-row_labels[,c(1,2)]


diamond_keys<-c("cmv","hhv6","hbv","hcv","hpylori","tgondii")


plot_cancer_all <- forest_plot(tab_cancer_all, 
                    col.key                       = "pathogen",
                    col.est                       ="estimate",
                    col.stderr                    ="stderr.x",
                    row.labels                    = row_labels, 
                    col.left                      = "prev",
                    col.left.heading              = c("Seropositivity (%) \nin the subcohort"),
                    col.right.heading             = "HR (95% CI)",
                    col.left.hjust                = 0.5,
                    col.right.hjust               = 0.5,
                    xlab                          = "Hazard Ratio",
                    base_size                     = 10,
                    scalepoints                   = TRUE,
                    xlim                          = c(0.49,2.5),
                    xticks                        = c(0.5,1.0,2.0),  
                    diamond                       =c("cmv","hhv6","hbv","hcv","hpylori","tgondii"),
                    nullval                       = 1)
                

title <- textGrob("Figure 3: Adjusted HRs for overall cancer incidence from multiple pathogens in Chinese adults",
                  gp = gpar(cex = 0.8 ,fontface="bold"), 
                  x = unit(0, "npc") + unit(1, "lines"),
                  y = unit(0, "npc") + unit(1, "lines"),
                  just = c(-0.05,0))

footnote <- textGrob("Adjusted for age, sex, region, education, smoking, alcohol, BMI, family history of cancer",
                       gp = gpar(cex = 0.6), 
                       x = unit(0, "npc") + unit(1, "lines"),
                       y = unit(0, "npc") + unit(1, "lines"),
                       just = c(-0.05,3))
                   
figure<-grid.arrange(plot_cancer_all$plot,
                     footnote,
                     nrow = 2,
                     heights = unit.c(unit(250, "mm"),
                               unit(1, "mm")))

ggsave(paste0(path_3,"Figure_2.png"),
       figure,
       bg="white",
       dpi = 450,
       width = 14,                      
       height = 26,                     
       units = "cm")

##############################################################################################################################################

