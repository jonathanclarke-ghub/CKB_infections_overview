
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
tab_meta$pathogen_label[36]<-"C.trachomatis"
tab_meta$pathogen_label[37:38]<-"T.gondii"
tab_meta$pathogen_label[39]<-"C.burnetii"
tab_meta$pathogen_label[40:47]<-"H.pylori"
tab_meta<-tab_meta[!duplicated(tab_meta$pathogen),c("pathogen_label","pathogen")]

# Vector of pathogen variables

pathogen<-paste0(tab_meta$pathogen,"_pp")
pathogen_label<-tab_meta$pathogen_label

data_all<-data

process_data<-function(group,title,name){

    
  if      (group=="subcohort")   {data<-data_all %>% filter(multiplexserology_cancer_b1_subcohort==1)}
  else if (group=="cancer")      {data<-data_all %>% filter(multiplexserology_cancer_b1_subcohort==0)}
  
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
res$pathogen_label<-pathogen_label[iii]

### Annual percentage change in seroprevalence ####

model<-lm(reformulate("age",response=pathogen[iii]),data=data) 

apc_est     <-(-(summary(model)$coefficients[2,1]*100))*10   
apc_stderr  <- summary(model)$coefficients[2,2]*100*10
apc_ci_lower<-sprintf(as.numeric(apc_est - 1.96*apc_stderr), fmt = '%#.2f') 
apc_ci_upper<-sprintf(as.numeric(apc_est + 1.96*apc_stderr), fmt = '%#.2f') 
apc_est     <-sprintf(as.numeric(apc_est), fmt = '%#.2f') 
res$apc<-paste0("Change (%) per 10 years ",apc_est," (",apc_ci_lower,", ",apc_ci_upper,")")
p_value<-summary(model)$coefficients[2,4]            
res$sig<-ifelse(p_value<0.05,1,0)
res$apc<-ifelse(res$sig==1,paste0(res$apc,"*"),res$apc)
tab[[iii]]<-res

}

tab_all<-as.data.frame(do.call(rbind,tab))
names(tab_all)<-c("n","prop","stderr","ci_lower","ci_upper","birth_cohort","p_var","pathogen_label","apc","sig")

### Scatter plot ####

pathogen<-unique(tab_all$p_var)

plot_pathogen<-list()

for(ii in pathogen) {
  
plot<-ggplot(data=tab_all[tab_all$p_var==ii,], aes(x=as.numeric(birth_cohort), y=as.numeric(prop))) +
    geom_smooth(method=loess,level=0.95) +
    geom_errorbar(aes(ymin=as.numeric(ci_lower),ymax=as.numeric(ci_upper)),colour="grey",linewidth=1,width=0) +
    geom_point()  + 
    scale_y_continuous(limits = c(0, 105)) + 
    xlab("Birth year") +
    ylab("Percent") +
    ggtitle(paste0(tab_all[tab_all$p_var==ii,"pathogen_label"][1])) +
#    annotate("text",x=1950,y=105,label=ifelse(tab_all[tab_all$p_var==ii,"sig"][1]==0,NA,tab_all[tab_all$p_var==ii,"apc"][1])) +
    annotate("text",x=1950,y=105,label=tab_all[tab_all$p_var==ii,"apc"][1]) +
    theme_classic() +
    theme(axis.text.x = element_text(face="bold",size=14),
          axis.title.x = element_text(face="bold",size=16),
          axis.text.y = element_text(face="bold",size=14),
          axis.title.y = element_text(face="bold",size=16),
          plot.title  = element_text(hjust=0.5,face="bold",size=18),
          legend.position = "none",
          plot.margin = margin(t=25,r=10,b=10,l=10))

if(ii %in% pathogen[1:16]){plot<-plot}
else if(ii %in% pathogen[17:20]) {plot<-plot + theme(plot.title = element_text(face = "bold.italic"))}




plot_pathogen[[ii]]<-plot

}

### Panel of plots ####

make_a_plot <-function (i){gridExtra::grid.arrange(fix_panel(plot_pathogen[[i]],height = unit(120,"mm")),
                                                   nrow = 1)
} 

list_of_plots <- lapply(c(4,5,6,7,9,10,13,14,1,2,3,15,16,17,20,19,18), make_a_plot)

grid_of_plots <- arrangeGrob(grobs = list_of_plots,
                             nrow = 3)   


title <- textGrob(title,
                  gp = gpar(cex = 2 ,fontface="bold"), 
                  x = unit(0, "npc") + unit(1, "lines"),
                  y = unit(0, "npc") + unit(1, "lines"),
                  just = c(0,0))

footnote <- textGrob("Dots and grey lines represent sero-prevalence estimates and 95% CI for every single birth year. The blue line represents estimates for these sero-prevalences derived from LOESS method, 
and the shaded areas show the 95% CI. Due to few participants that were born before 1932 or after 1971, those who were born between 1927 and 1931, or 1971 and 1976 were combined, separately. 
LOESS, locally weighted scatterplot smoothing",
                       gp = gpar(cex = 1), 
                       x = unit(0, "npc") + unit(1, "lines"),
                       y = unit(0, "npc") + unit(1, "lines"),
                       just = c(0,1.5))

figure<-grid.arrange(grid_of_plots,
                     footnote,
                     nrow = 2,
                     heights = unit.c(unit(460, "mm"),
                                      unit(5, "mm")))

ggsave(paste0(path_3,name,".png"),
       figure,
       dpi = 450,
       width = 76,
       height = 52,
       units = "cm")

}

process_data("subcohort","eFigure 4: Seroprevalance of multiple pathogens by year of birth among subcohort","eFigure_4")
process_data("cancer","eFigure 9: Seroprevalance of multiple pathogens by year of birth among incident cancer cases","eFigure_9")













