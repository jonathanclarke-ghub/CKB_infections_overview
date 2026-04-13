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

### Import participant data ####

#tab_all<-readRDS(paste0(path_2,"tab_cancer_all_cch.rds"))
#tab_all<-readRDS(paste0(path_2,"tab_cancer_all_cox.rds"))
#tab_all<-readRDS(paste0(path_2,"tab_cancer_all_cch_x2.rds"))
tab_all<-readRDS(paste0(path_2,"tab_cancer_all_cox_x2.rds"))

tab_pathogen_pos<-readRDS(paste0(path_2,"tab_pathogen_pos.rds"))
tab_antigen_pos<-readRDS(paste0(path_2,"tab_antigen_pos.rds"))

## Antigens/pathogens and site specific cancer 2 ####

#tab_cancer_selected<-tab_all %>% filter(cancer_code %in% c("ep_cancer_C16.0_combined","ep_cancer_C16x.0_combined"), analysis=="cancer_selected_cch")
#tab_cancer_selected<-tab_all %>% filter(cancer_code %in% c("ep_cancer_C16.0_combined","ep_cancer_C16x.0.8_combined"), analysis=="cancer_selected_cch")
#tab_cancer_selected<-tab_all %>% filter(cancer_code %in% c("ep_cancer_C16.0_combined","ep_cancer_C16x.0_combined"), analysis=="cancer_selected_cox")
tab_cancer_selected<-tab_all %>% filter(cancer_code %in% c("ep_cancer_C16.0_combined","ep_cancer_C16x.0.8_combined"), analysis=="cancer_selected_cox")


tab_pathogen_pos_selected<-tab_pathogen_pos[tab_pathogen_pos$pathogen %in% c("hpylori"),]
tab_antigen_pos_selected<-tab_antigen_pos[tab_antigen_pos$pathogen %in% grep("^hpylori",tab_antigen_pos$pathogen,value=TRUE),]
tab_pathogen_antigen_pos_selected<-rbind(tab_pathogen_pos_selected,tab_antigen_pos_selected)
tab_cancer_selected<-merge(tab_cancer_selected,tab_pathogen_antigen_pos_selected,by="pathogen")

tab_cancer_selected$col_key<-tab_cancer_selected$pathogen
#tab_cancer_selected$cancer_code<-factor(tab_cancer_selected$cancer_code, levels=c("ep_cancer_C16.0_combined","ep_cancer_C16x.0_combined"))
tab_cancer_selected$cancer_code<-factor(tab_cancer_selected$cancer_code, levels=c("ep_cancer_C16.0_combined","ep_cancer_C16x.0.8_combined"))
row_labels<-data.frame(col_key=c(tab_cancer_selected$pathogen[c(2,4,6,8,9,12,13,15,17)]),
                       label  =c(tab_cancer_selected$label[c(2,4,6,8,9,12,13,15,17)]))

row_labels$label[c(2:9)] <- glue::glue("<span style='color:transparent;'>XX</span>{row_labels$label[c(2:9)]}")
row_labels$label[c(1)] <- glue::glue("<span style='color:transparent;'>X</span>{row_labels$label[c(1)]}")

cases_C160<-tab_cancer_selected[tab_cancer_selected$cancer_code=="ep_cancer_C16.0_combined","n_cases"][1]
cases_C160<-format(as.numeric(cases_C160),big.mark=",",scientific=FALSE)
#cases_C16x<-tab_cancer_selected[tab_cancer_selected$cancer_code=="ep_cancer_C16x.0_combined","n_cases"][1]
cases_C16x<-tab_cancer_selected[tab_cancer_selected$cancer_code=="ep_cancer_C16x.0.8_combined","n_cases"][1]

cases_C16x<-format(as.numeric(cases_C16x),big.mark=",",scientific=FALSE)


plot_a <- forest_plot(tab_cancer_selected[tab_cancer_selected$cancer_code=="ep_cancer_C16.0_combined",],
                                      col.key                       = "col_key",
                                      col.est                       ="estimate",
                                      col.stderr                    ="stderr.x",
                                      row.labels                    = row_labels, 
                                      col.left                      = "prev",
                                      col.left.heading              = "Seropositivity (%)\nin the subcohort",
                                  #    col.right.hjust               = 1,
                                      col.left.hjust                = 0.5,
                                      col.right.hjust               = 0.5,
                                      xlab                          = "Hazard Ratio",
                                      col.right.heading             = "HR (95% CI)",
                                      right.space                   = unit(0,"mm"), 
                                      scalepoints                   = TRUE,
                                      bold.labels                   = "hpylori",
                                      xlim                          = c(0.75,4.11),
                                      xticks                        = c(1.0,2.0,3.0,4.0),
                                      panel.headings                = paste0("Cardia gastric cancer \n(C16.0), n=",cases_C160),
                                      diamond                       ="hpylori",
                                      nullval                       = 1)


#plot_b <- forest_plot(tab_cancer_selected[tab_cancer_selected$cancer_code=="ep_cancer_C16x.0_combined",],
plot_b <- forest_plot(tab_cancer_selected[tab_cancer_selected$cancer_code=="ep_cancer_C16x.0.8_combined",],
                                      col.key                       = "col_key",
                                      col.est                       ="estimate",
                                      col.stderr                    ="stderr.x",
                                      row.labels                    = row_labels,
                                   #   col.right.hjust               = 1,  
                                      scalepoints                   = TRUE,
                                      right.space                    = unit(20,"mm"),
                                      xlim                          = c(0.73,4.11),
                                      xticks                        = c(1.0,2.0,3.0,4.0),
                                      col.right.hjust               = 0.5,
                                      xlab                          = "Hazard Ratio",
                                      col.right.heading             = "HR (95% CI)",
                                   #   panel.headings                =paste0("Non-cardia gastric cancer \n(C16.1-C16.9), n=",cases_C16x),
                                      panel.headings                =paste0("Non-cardia gastric cancer \n(C16.1-C16.8), n=",cases_C16x),
                                      diamond                       ="hpylori",
                                      nullval                       = 1)

plot_b<-plot_b$plot +
  theme(axis.text.y=element_blank())

panel_1<-gridExtra::grid.arrange(fix_panel(plot_a$plot,width = unit(50,"mm")),
                                 nrow = 1) 

panel_2<-gridExtra::grid.arrange(fix_panel(plot_b,width = unit(50,"mm")),
                                 nrow = 1) 

figure_1<-grid.arrange(panel_1,
                       panel_2,
                       nrow = 1)

title <- textGrob(substitute(paste(bold("eFigure 10: Adjusted HRs for gastric cancer and non-cardia gastric cancers from"),bolditalic("H. pylori"), bold(" in Chinese adults"))),
                  gp = gpar(cex = 0.8 ,fontface="bold"), 
                  x = unit(0, "npc") + unit(1, "lines"),
                  y = unit(0, "npc") + unit(1, "lines"),
                  just = c(-0.04,-1))

footnote <- textGrob("Adjusted for age, sex, region, education, smoking, alcohol, BMI, family history of cancer",
                     gp = gpar(cex = 0.6), 
                     x = unit(0, "npc") + unit(1, "lines"),
                     y = unit(0, "npc") + unit(1, "lines"),
                     just = c(-0.05,2))

figure_2<-grid.arrange(figure_1,
                      footnote,
                      nrow = 2,
                      heights = unit.c(unit(100, "mm"),
                                      unit(1, "mm")))

#ggsave(paste0(path_3,"eFigure_10_cch_new.png"),
 #      figure_2,
  #     bg="white",
   #    dpi = 450,
    #   width = 24,                      
     #  height = 12,                     
      # units = "cm")

#ggsave(paste0(path_3,"eFigure_10_cox_new.png"),
 #      figure_2,
  #     bg="white",
   #    dpi = 450,
    #   width = 24,                      
     #  height = 12,                     
      # units = "cm")

#ggsave(paste0(path_3,"eFigure_10_cch_new_x2_8.png"),
#      figure_2,
#     bg="white",
#    dpi = 450,
#   width = 24,                      
#  height = 12,                     
#  units = "cm")

ggsave(paste0(path_3,"eFigure_11.png"),
       figure_2,
       bg="white",
       dpi = 450,
       width = 24,                      
       height = 12,                     
       units = "cm")

####################################################################################################################################################
