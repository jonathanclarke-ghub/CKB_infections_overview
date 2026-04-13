
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

## Hetrogeneity test ####

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

path_1<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Raw_data_7/"
path_2<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Project_data_7/"
path_3<-"K:/kadoorie/Staff_Folders/JonathanC/Projects/Serology/Tables_figures/Tables_figures_10/"

### Import simulated data ####

data<-readRDS(paste0(path_2,"data_serology_2024_7.rds"))

data<-as.data.frame(data)

### Import multiplex meta data ####

tab_meta<-readRDS(paste0(path_1,"multiplex_cancer_case_cohort_meta.rds"))
tab_meta$pathogen_label_old<-tab_meta$pathogen_label
tab_meta$pathogen_label[36:47]<-""
tab_meta$pathogen_label[36]<-"C. trachomatis"
tab_meta$pathogen_label[37:38]<-"T. gondii"
tab_meta$pathogen_label[39]<-"C. burnetii"
tab_meta$pathogen_label[40:47]<-"H. Pylori"
tab_meta<-tab_meta[!duplicated(tab_meta$pathogen),c("pathogen_label","pathogen")]

# Vector of pathogen variables

pathogen<-paste0(tab_meta$pathogen,"_pp")
pathogen_label<-tab_meta$pathogen_label

### Barplot - Sero prevalence by region ####

# Vector of regions

region<-unique(data$region_name)

data_all<-data

group<-c("subcohort","cancer")  

results_group<-list()
  
for(i in  group){  
  
  if      (i=="subcohort")   {data<-data_all %>% filter(multiplexserology_cancer_b1_subcohort==1)}
  else if (i=="cancer")      {data<-data_all %>% filter(multiplexserology_cancer_b1_subcohort==0)}
  
### Sero-prevalance by pathogen and region ####

results_pathogen<-list()

for(iii in 1:length(pathogen)){
  
  results<-list()
  
  for(ii in region){
    
    pos<-data[data$region_name==ii,pathogen[iii]]
    
     # Calculate summary statistics
    
    n <- length(pos)
    prop<-mean(pos,rm.na=TRUE)
    se<-std.error(pos)
    
    # Calculate margin of error - two-sided 95% CI
    
    margin <- qnorm(0.975)*sqrt(prop*(1-prop)/n)
    
    # Calculate lower and upper bounds of confidence interval
    
    ci_lower <- (prop - margin)*100
    ci_upper <- (prop + margin)*100
    prop     <- prop*100
    se<-se*100
    area<-ii

    results[[ii]]<-c(n,prop,se,ci_lower,ci_upper,area)
    
  }
  
  results_region<-as.data.frame(do.call(rbind,results))
  results_region$p_var<-pathogen[iii]
  results_region$pathogen_name<-pathogen_label[iii]
  
  results_pathogen[[iii]]<-results_region
  
}

tab<-as.data.frame(do.call(rbind,results_pathogen))
names(tab)<-c("n","prop","stderr","ci_lower","ci_upper","region_name","p_var","pathogen_label")

tab$group<-paste0(i)

results_group[[i]]<-tab

}

tab<- as.data.frame(do.call(rbind,results_group)) 

### Format table ####

tab$ci_low<-sprintf(as.numeric(tab$ci_lower), fmt = '%#.1f')  
tab$ci_up<-sprintf(as.numeric(tab$ci_upper), fmt = '%#.1f') 
tab$prev<-sprintf(as.numeric(tab$prop), fmt = '%#.1f') 
tab$prev_ci<-paste0(tab$prev," (",tab$ci_low," to ",tab$ci_up,")") 

# To order bars by region name

tab$region_name<-factor(tab$region_name,levels=c("Gansu","Henan","Zhejiang","Sichuan","Hunan","Harbin","Qingdao","Liuzhou","Haikou","Suzhou"))

# To colour bars by region

tab$region_urban<-as.factor(ifelse(tab$region_name %in% c("Gansu","Henan","Zhejiang","Sichuan","Hunan"),"Rural","Urban"))
tab$region_urban<-factor(tab$region_urban, levels=c("Rural","Urban"))
tab$prop<-as.numeric(tab$prop)
tab$h_adj<- ifelse(tab$group=="cancer",0.5,0.5)
tab$study_arm<-ifelse(tab$group=="cancer" & tab$region_urban=="Rural","Pre-cancer",
                      ifelse(tab$group=="cancer" & tab$region_urban=="Urban","UC",
                             ifelse(tab$group=="subcohort" & tab$region_urban=="Rural","Subcohort","US")))


### Panel of plots - drop EBV,HIV,HTLV ####

tab<-tab %>% filter(!(p_var %in% c("ebv_pp","hiv_pp","htlv1_pp")))

pathogen<-unique(tab$p_var)

plot_xy<-list()

for(ii in pathogen[c(4,10,13)]) {
  
plot<-ggplot(data=tab[tab$p_var==ii & tab$group=="subcohort",], aes(x=region_name, y=prop,fill=region_urban)) +
  geom_bar(stat = "identity",position = position_dodge(0.9),width=0.9,colour="black",linewidth=0.5) +
  scale_y_continuous(limits=c(0,115),
                     breaks=seq(0, 100, 25),
                     labels=seq(0, 100, 25),expand=c(0,0)) +
  scale_fill_manual(values = c("white","grey")) +
  geom_text(aes(label=prev),vjust=tab[tab$p_var==ii & tab$region_urban=="Rural","h_adj"],hjust=-0.3,size = 4.5) + 
  ylab("Percent") +  
  xlab("Region") +
  ggtitle(paste0(tab[tab$p_var==ii,"pathogen_label"][1])) +
  theme_classic() +
  theme(axis.text.y = element_text(face="bold",size=14),
        axis.line.y =element_blank(),
        axis.text.x = element_text(face="bold",hjust = 0.5,vjust=1,size=14),
        axis.title.y = element_text(face="bold",vjust=0,angle=90,size=16),
        axis.title.x = element_text(face="bold",vjust=0,size=16),
        plot.title  = element_text(hjust=0.5,vjust=1.5,face="bold",size=18),
        legend.position="none") +
  coord_flip()

plot_xy[[ii]]<-plot

}

plot_x<-list()

for(ii in pathogen[c(5:9,11,1:3,12)]) {
  
  plot<-ggplot(data=tab[tab$p_var==ii & tab$group=="subcohort",], aes(x=region_name, y=prop,fill=region_urban)) +
    geom_bar(stat = "identity",position = position_dodge(0.9),width=0.9,colour="black",linewidth=0.5) +
    scale_y_continuous(limits=c(0,115),
                       breaks=seq(0, 100, 25),
                       labels=seq(0, 100, 25),expand=c(0,0)) +
    scale_fill_manual(values = c("white","grey")) +
    geom_text(aes(label=prev),vjust=tab[tab$p_var==ii & tab$region_urban=="Rural","h_adj"],hjust=-0.3,size = 4.5) + 
    ylab("Percent") +
    xlab("Region") +
    ggtitle(paste0(tab[tab$p_var==ii,"pathogen_label"][1])) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y =element_blank(),
          axis.text.x = element_text(face="bold",hjust = 0.5,vjust=1,size=14),
          axis.title.y = element_blank(),
          axis.title.x = element_text(face="bold",vjust=0,size=16),
          plot.title  = element_text(hjust=0.5,vjust=1.5,face="bold",size=18),
          legend.position="none") +
    coord_flip()

    plot_x[[ii]]<-plot
}

plot_xi<-list()

for(ii in pathogen[c(14,17,16)]) {
  
  plot<-ggplot(data=tab[tab$p_var==ii & tab$group=="subcohort",], aes(x=region_name, y=prop,fill=region_urban)) +
    geom_bar(stat = "identity",position = position_dodge(0.9),width=0.9,colour="black",linewidth=0.5) +
    scale_y_continuous(limits=c(0,115),
                       breaks=seq(0, 100, 25),
                       labels=seq(0, 100, 25),expand=c(0,0)) +
    scale_fill_manual(values = c("white","grey")) +
    geom_text(aes(label=prev),vjust=tab[tab$p_var==ii & tab$region_urban=="Rural","h_adj"],hjust=-0.3,size = 4.5) + 
    ylab("Percent") +
    xlab("Region") +
    ggtitle(paste0(tab[tab$p_var==ii,"pathogen_label"][1])) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y =element_blank(),
          axis.text.x = element_text(face="bold",hjust = 0.5,vjust=1,size=14),
          axis.title.y = element_blank(),
          axis.title.x = element_text(face="bold",vjust=0,size=16),
          plot.title  = element_text(hjust=0.5,vjust=1.5,face="bold.italic",size=18),
          legend.position="none") +
    coord_flip()
  
  plot_xi[[ii]]<-plot
}

plot_x_legend<-list()

for(ii in pathogen[c(15)]) {
  
  plot<-ggplot(data=tab[tab$p_var==ii & tab$group=="subcohort",], aes(x=region_name, y=prop,fill=region_urban)) +
    geom_bar(stat = "identity",position = position_dodge(0.9),width=0.9,colour="black",linewidth=0.5) +
    scale_y_continuous(limits=c(0,115),
                       breaks=seq(0, 100, 25),
                       labels=seq(0, 100, 25),expand=c(0,0)) +
    scale_fill_manual(values = c("white","grey")) +
    geom_text(aes(label=prev),vjust=tab[tab$p_var==ii & tab$region_urban=="Rural","h_adj"],hjust=-0.3,size = 4.5) + 
    ylab("Percent") + 
    xlab("Region") + 
    labs(fill="Area") +
    ggtitle(paste0(tab[tab$p_var==ii,"pathogen_label"][1])) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.line.y =element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(face="bold",hjust = 0.5,vjust=1,size=16),
          axis.title.y = element_blank(),
          legend.text = element_text(face="bold",size=14),
          legend.title = element_text(face="bold",size=16),
          axis.title.x = element_text(face="bold",vjust=0,size=16),
          plot.title  = element_text(hjust=0.5,vjust=1.5,face="bold.italic",size=18),
          legend.position=c(0.85,0.85)) +
    coord_flip()
  
  plot_x_legend[[ii]]<-plot
  
}

plot_all1<-c(plot_x[1:5],plot_x[6:10],plot_xi[1:3],plot_x_legend[1])

make_a_plot <-function (i){gridExtra::grid.arrange(fix_panel(plot_all1[[i]],height = unit(90,"mm")),nrow = 1)} 
list_of_plots1 <- lapply(c(1:length(plot_all1)), make_a_plot)

grid_of_plots1 <- arrangeGrob(grobs = list_of_plots1,
                              nrow = 3)

plot_all2<-c(plot_xy[1],plot_xy[2],plot_xy[3])

make_a_plot <-function (i){gridExtra::grid.arrange(fix_panel(plot_all2[[i]],height = unit(90,"mm")),nrow = 1)} 
list_of_plots2 <- lapply(c(1:length(plot_all2)), make_a_plot)


grid_of_plots2 <- arrangeGrob(grobs = list_of_plots2,
                              nrow = 3)

grid_of_plots<-grid.arrange(grid_of_plots2,
                            grid_of_plots1,
                            ncol=2,
                            widths = unit.c(unit(130, "mm"),
                                            unit(500, "mm")))

figure<-grid.arrange(grid_of_plots,
                     nrow = 1,
                     heights = unit.c(unit(375, "mm")))

ggsave(paste0(path_3,"eFigure_3.png"),
       figure,
       dpi = 450,
       width = 64,
       height = 38,
       units = "cm")



########################################################################################################################################


