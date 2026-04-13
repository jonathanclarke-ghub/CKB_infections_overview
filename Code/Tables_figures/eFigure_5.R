
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

### Import simulated data ####

data_all<-readRDS(paste0(path_2,"data_serology_2024_7.rds"))
data_all<-as.data.frame(data_all)

mean(data_all$pathogen_index[data_all$multiplexserology_cancer_b1_subcohort==1])
sd(data_all$pathogen_index[data_all$multiplexserology_cancer_b1_subcohort==1])
mean(data_all$pathogen_index[data_all$multiplexserology_cancer_b1_subcohort==0])
sd(data_all$pathogen_index[data_all$multiplexserology_cancer_b1_subcohort==0])

plot_list<-list()

for(cc in  1:2){

data<-data_all %>% filter(multiplexserology_cancer_b1_subcohort==(cc-1))

### Import multiplex meta data ####

tab_meta<-readRDS(paste0(path_1,"multiplex_cancer_case_cohort_meta.rds"))
tab_meta$pathogen_label_old<-tab_meta$pathogen_label
tab_meta$pathogen_label[36:47]<-""
tab_meta$pathogen_label[36]<-"C.trachomatis"
tab_meta$pathogen_label[37:38]<-"T.gondii"
tab_meta$pathogen_label[39]<-"C.burnetii"
tab_meta$pathogen_label[40:47]<-"H.pylori"
tab_meta<-tab_meta[!duplicated(tab_meta$pathogen),c("pathogen_label","pathogen")]

### Bar plot - all pathogens ####

tab<-as.data.frame(table(data$pathogen_index))
names(tab)<-c("n","freq")

# Insert missing rows (rows with frequency=0)

nn<-as.data.frame(n<-0:20)
names(nn)<-"n"
tab<-merge(nn,tab,by="n",all=TRUE)

# Summary data

tab_mn<-tab
tab_mn$freq[is.na(tab_mn$freq)]<-0
tab_mn$n<-as.numeric(tab_mn$n)
tab_mn$freq<-as.numeric(tab_mn$freq) 
mean<-sum(tab_mn$n*tab_mn$freq)/sum(tab_mn$freq)
mean_1<-round(mean,1)

max_freq<-max(tab$freq[!is.na(tab$freq)])

# Create bar plot

plot_bar<-ggplot(data=tab, aes(x=n, y=freq)) +
  geom_bar(stat="identity",colour="white",fill="grey") +
  scale_x_continuous(breaks = seq(0, 20, by = 1)) +
  ylim(0,max_freq + 10) +
  geom_text(aes(label=freq),vjust=-1,size = 4) +
  geom_segment(x=mean, y=-10, xend=mean, yend=max_freq+10,linetype="dashed",linewidth=1) +
  annotate("text", x = mean_1+0.75, y = 150, label = paste0("Mean \n= ",mean_1) ,size = 4) + 
  xlab("No. of pathogens") +
  ylab("No. of participants") +
  ggtitle("All") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),    
        axis.title.y = element_text(size=12),
        plot.title  = element_text(hjust=0.5,size=18,face="bold"),
        legend.position = "none",
        plot.margin = margin(t=5,r=10,b=0,l=10))

### Scatter plot ####

# Vector of pathogen variables

pathogen<-c("pathogen_11_plus")
pathogen_name <- "Sero-prevalence for > 10 pathogens by year of birth"

# Vector of birth years

cohort<-unique(data$birth_cohort)

### Sero-prevalance by pathogen and birth year ####

tab<-list()

for(iii in 1:length(pathogen)){
  
  results<-list()
  
  for(ii in cohort){
    
    pos<-data[data$birth_cohort==ii,pathogen[iii]]
    
    # Calculate summary statistics
    
    n <- length(pos)
    prop<-mean(pos,rm.na=TRUE)
    se<-std.error(pos)
    
    # Calculate margin of error - two-sided 95% CI
    
    margin <- qnorm(0.975)*sqrt(prop*(1-prop)/n)
    
    # Calculate lower and upper bounds of confidence interval
    
    ci_lower <- (prop - margin)*100
    ci_upper <- (prop + margin)*100
    
    # Convert to percentage
    
    prop<- prop*100
    se<-se*100
    
    birth_cohort<-ii
    
    results[[ii]]<-c(n,prop,se,ci_lower,ci_upper,birth_cohort)
    
  }
  
  res<-as.data.frame(do.call(rbind,results))
  res$birth_cohort<-row.names((tab))
  res$p_var<-pathogen[iii]
  res$pathogen_name<-pathogen_name[iii]
  
  tab[[iii]]<-res
  
}

tab_all<-as.data.frame(do.call(rbind,tab))
names(tab_all)<-c("n","prop","stderr","ci_lower","ci_upper","birth_cohort","p_var","pathogen_name")

plot_scatter<-ggplot(data=tab_all, aes(x=as.numeric(birth_cohort), y=as.numeric(prop))) +
    geom_smooth(method=loess,level=0.95) +
    geom_errorbar(aes(ymin=as.numeric(ci_lower),ymax=as.numeric(ci_upper)),colour="grey",linewidth=1,width=0) +
    geom_point()  + 
    scale_y_continuous(limits = c(0, 105)) + 
    xlab("Birth year") +
    ylab("Percent") +
    ggtitle(paste0(tab_all[,"pathogen_name"][1])) +
    theme_classic() +
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          plot.title  = element_text(hjust=0.5,size=18,face="bold"),
          legend.position = "none",
          plot.margin = margin(t=5,r=10,b=0,l=10))
 
plot_all<-list(plot_bar,plot_scatter)

plot_list[[cc]]<-plot_all

}

plot_list<-purrr::list_flatten(plot_list)

### Panel of plots ####

make_a_plot <-function (i){gridExtra::grid.arrange(fix_panel(plot_list[[i]],height = unit(140,"mm")),nrow = 1)} 
list_of_plots <- lapply(c(1:4), make_a_plot)
                             
grid_of_plots_can <- arrangeGrob(grobs = list_of_plots[1:2],nrow = 1)
grid_of_plots_sub <- arrangeGrob(grobs = list_of_plots[3:4],nrow = 1)

title_sub <- textGrob("Subcohort",
                  gp = gpar(cex = 2 ,fontface="bold"), 
                  x = unit(0, "npc") + unit(1, "lines"),
                  y = unit(0, "npc") + unit(1, "lines"),
                  just = c(-5,2.5))

title_can <- textGrob("Incident Cancer",
                      gp = gpar(cex = 2 ,fontface="bold"), 
                      x = unit(0, "npc") + unit(1, "lines"),
                      y = unit(0, "npc") + unit(1, "lines"),
                      just = c(-3.3,2.5))

footnote <- textGrob("",
                     gp = gpar(cex = 1), 
                     x = unit(0, "npc") + unit(1, "lines"),
                     y = unit(0, "npc") + unit(1, "lines"),
                     just = c(0,0.75))

figure<-grid.arrange(title_sub,
                     grid_of_plots_sub,
                     title_can,
                     grid_of_plots_can,
                     footnote,
                     nrow = 5,
                     heights = unit.c(unit(1, "mm"),
                                      unit(200, "mm"),
                                      unit(1, "mm"),
                                      unit(200, "mm"),
                                      unit(1, "mm")))

ggsave(paste0(path_3,"eFigure_5.png"),
       figure,
       dpi = 450,
       width = 48,
       height = 45,
       units = "cm")

#########################################################################################################################################################


