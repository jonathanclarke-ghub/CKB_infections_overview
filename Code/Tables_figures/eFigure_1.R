
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
library(stringr)

path_1<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Raw_data_7/"
path_2<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Project_data_7/"
path_3<-"K:/kadoorie/Staff_Folders/JonathanC/Projects/Serology/Tables_figures/Tables_figures_10/"

### Import raw multiplex data ####

data<-readRDS(paste0(path_2,"data_serology_2024_7.rds"))

data<-as.data.frame(data)

### Import multiplex meta data ####

tab_meta<-readRDS(paste0(path_1,"multiplex_cancer_case_cohort_meta.rds"))
tab_meta[tab_meta$variable_name=="hpv16_e6","proposed_cutoff"]<-1000      # Updated cut-off

tab_meta$pathogen_label_old<-tab_meta$pathogen_label
tab_meta$pathogen_label[36:47]<-""
tab_meta$pathogen_label[36]<-"C.trachomatis"
tab_meta$pathogen_label[37:38]<-"T.gondii"
tab_meta$pathogen_label[39]<-"C.burnetii"
tab_meta$pathogen_label[40:47]<-"H.pylori"
tab_meta<-tab_meta[,-c("antigen")]  # Avoid conflict in names

antigen<-tab_meta$variable_name

###  Barplot to examine sero-positive status/distribution of each antigen #### 

fun<-function(ii,plot_title){

cutoff<-as.numeric(tab_meta[tab_meta$variable_name==antigen[ii],"proposed_cutoff"]) 

plot<-ggplot(data,aes(x=data[,antigen[ii]], fill = ifelse(data[,antigen[ii]]> cutoff,"orange","turquoise"))) +
  geom_histogram(bins=50) + 
  xlab("MFI") +
  ylab("Frequency") +
  ggtitle(plot_title) +
  theme_classic() + 
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(), 
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        plot.title  = element_text(hjust=0.5,face="bold"),
        legend.position = "none")

if(ii %in% 1:35){plot<-plot}
else if(ii %in% 36:47) {plot<-plot + theme(plot.title = element_text(face = "bold.italic"))}

}

plot_1 <-fun(1, "BKV VP1")
plot_2 <-fun(2, "JCV VP1")
plot_3 <-fun(3, "MCV VP1")
plot_4 <-fun(4, "HSV-1 gg")
plot_5 <-fun(5, "HSV-2 mgG unique")
plot_6 <-fun(6, "VZV gE/gl")
plot_7 <-fun(7, "CMV pp150 Nter")
plot_8 <-fun(8, "CMV pp 28")
plot_9 <-fun(9, "CMV pp 52")
plot_10<-fun(10,"EBV EBNA-1")
plot_11<-fun(11,"EBV EA-D")
plot_12<-fun(12,"EBV VCA p18")
plot_13<-fun(13,"EBV ZEBRA")
plot_14<-fun(14,"EBV BFRF1")
plot_15<-fun(15,"EBV BGLF2")
plot_16<-fun(16,"EBV BXLF1")
plot_17<-fun(17,"HHV-6 IE1A")
plot_18<-fun(18,"HHV-6 IE1B")
plot_19<-fun(19,"HHV-7 U14")
plot_20<-fun(20,"HIV env")
plot_21<-fun(21,"HIV gag")
plot_22<-fun(22,"HTLV-1 env")
plot_23<-fun(23,"HTLV-1 gag")
plot_24<-fun(24,"HBV HBc")
plot_25<-fun(25,"HBV HBe")
plot_26<-fun(26,"HCV NS3")
plot_27<-fun(27,"HCV Core")
plot_28<-fun(28,"HPV-16 E1")
plot_29<-fun(29,"HPV-16 E2")
plot_30<-fun(30,"HPV-16 E6")
plot_31<-fun(31,"HPV-16 E7")
plot_32<-fun(32,"HPV-16 L1")
plot_33<-fun(33,"HPV-18 E6")
plot_34<-fun(34,"HPV-18 E7")
plot_35<-fun(35,"HPV-18 L1")
plot_36<-fun(36,"C.trachomatis pGP3")
plot_37<-fun(37,"T.gondii p22")
plot_38<-fun(38,"T.gondii sag-1")
plot_39<-fun(39,"C.burnetii CBU-1910 Com1")
plot_40<-fun(40,"H.pylori HP0547-CagA (N+C)")
plot_41<-fun(41,"H.pylori HP0875-Catalase")
plot_42<-fun(42,"H.pylori HP0010-GroEL")
plot_43<-fun(43,"H.pylori HP-HcpC")
plot_44<-fun(44,"H.pylori HP-HopA")
plot_45<-fun(45,"H.pylori HP0305")
plot_46<-fun(46,"H.pylori HP1564-OMP")
plot_47<-fun(47,"H.pylori VHP0887-VacA (N+C)")

plot<-list(plot_4, plot_5, plot_6, plot_14, plot_15, plot_16, plot_11, plot_10,
           plot_12,plot_13,plot_7,plot_8,plot_9,plot_17,plot_18,plot_19,
           plot_24,plot_25,plot_27,plot_26,plot_20,plot_21,plot_22,plot_23,
           plot_1,plot_2,plot_3,plot_28,plot_29,plot_30,plot_31,plot_32,
           plot_33,plot_34,plot_35,plot_36,plot_40,plot_41,plot_42,plot_43,
           plot_44,plot_45,plot_46,plot_47,plot_39,plot_37,plot_38)

### Panel of plots ####

make_a_plot <-function (i){gridExtra::grid.arrange(fix_panel(plot[[i]],height = unit(50,"mm")),
                                                   nrow = 1)
} 

list_of_plots <- lapply(1:length(tab_meta$variable_name), make_a_plot)

grid_of_plots <- arrangeGrob(grobs = list_of_plots,
                             nrow = 9)

title <- textGrob("eFigure 2: Distribution of median fluorescence intensity (MFI) of antigens measured",
                  gp = gpar(cex = 2 ,fontface="bold"), 
                  x = unit(0, "npc") + unit(1, "lines"),
                  y = unit(0, "npc") + unit(1, "lines"),
                  just = c(0,0))

footnote <- textGrob("For each antigen, orange and blue represent sero-status of positive and negative, respectively. ",
                     gp = gpar(cex = 1.25), 
                     x = unit(0, "npc") + unit(1, "lines"),
                     y = unit(0, "npc") + unit(1, "lines"),
                     just = c(0,-15))

figure<-grid.arrange(grid_of_plots,
                     footnote,
                     nrow = 2,
                     heights = unit.c(unit(700, "mm"),
                                      unit(1, "mm")))

ggsave(paste0(path_3,"eFigure_1t2.png"),
       figure,
       dpi = 450,
       width = 50,
       height = 75,
       units = "cm")


