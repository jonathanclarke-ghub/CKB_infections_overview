### Libraries ####

library(DT)
library(tidyverse)
library(tidytext)
library(reshape) 
library(dplyr)
library(flextable)
library(officer)
library(plotrix)
library(ltm)

### Paths ####

path_1<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Raw_data_7/"
path_2<-"K:/ckb_data/Staff_Folders/jonathan_clarke/Serology/Project_data_7/"
path_3<-"K:/kadoorie/Staff_Folders/JonathanC/Projects/Serology/Tables_figures/Tables_figures_10/"

### Import participant data ####

data<-readRDS(paste0(path_2,"data_serology_2024_7.rds"))

data<-data %>% filter(multiplexserology_cancer_b1_subcohort==1)

### Import meta data ####

tab_meta<-readRDS(paste0(path_1,"multiplex_cancer_case_cohort_meta.rds"))
tab_meta$pathogen_label_old<-tab_meta$pathogen_label
tab_meta$pathogen_label[36:47]<-""
tab_meta$pathogen_label[36]<-"C. trachomatis"
tab_meta$pathogen_label[37:38]<-"T. gondii"
tab_meta$pathogen_label[39]<-"C. burnetii"
tab_meta$pathogen_label[40:47]<-"H. pylori"

# Import data dictionary

column_details<-read.csv(paste0(path_1,"column_details.csv"))
value_definitions<-read.csv(paste0(path_1,"value_definitions.csv"))

### Exposures ####

#data$age_cat5<-as.factor(data$age_cat5)              # Age
data$is_female<-as.factor(data$is_female)             # Sex 
data$region_code<-as.factor(data$region_code)         # Region

table(data$region_is_urban,useNA="ifany")
table(data$education_cat2,useNA="ifany")               # Education (6 years plus)
table(data$household_income_cat2,useNA="ifany")       # Household income (>=20k vs. <20K)

table(data$ever_reg_smoker,useNA="ifany")                   # Ever regular smoker
data$ever_reg_smoker_m<-data$ever_reg_smoker
data$ever_reg_smoker_f<-data$ever_reg_smoker

table(data$ever_reg_alcohol,useNA="ifany")                  # Ever regular alcohol drinker
data$ever_reg_alcohol_m<-data$ever_reg_alcohol         
data$ever_reg_alcohol_f<-data$ever_reg_alcohol

table(data$reg_tea_drinker,useNA="ifany")                   # Ever regular tea drinker
table(data$bmi_overweight,useNA="ifany")                    # BMI (Overweight BMI>=25)
table(data$waist_cat2_m,useNA="ifany")                      # Central obesity (WC >=94cm (male) or >=80 cm (women)
table(data$waist_cat2_f,useNA="ifany")                     # Central obesity (WC >=94cm (male) or >=80 cm (women)
table(data$met_high,useNA="ifany")                          # Physical activity (High vs. Low)
table(data$srh_poor,useNA="ifany")                          # Self-rated poor health
table(data$diet_had_shortage,useNA="ifany")                 # Famine experience
table(data$reg_spice,useNA="ifany")                         # Regular spicy food
table(data$reg_fresh_fruit,useNA="ifany")                   # Regular fresh fruit
table(data$reg_preserved_veg,useNA="ifany")                 # Regular preserved vegetables
table(data$reg_meat_poultry,useNA="ifany")                  # Regular meat/poultry
table(data$chd_stroke_diag,useNA="ifany")                   # CHD or Stroke/TIA 
table(data$hep_b_diag,useNA="ifany")                        # Hep B Test
table(data$rheum_heart_dis_diag,useNA="ifany")              # Rheumatic heart disease
table(data$diabetes_diag,useNA="ifany")                     # Diabetes (SR or SD)
table(data$hypertension_diag,useNA="ifany")                 # Hypertension (SR or SD)
table(data$emph_bronc_diag,useNA="ifany")                   # Emphysema/bronchitis
table(data$peptic_ulcer_diag,useNA="ifany")                 # Peptic ulcer
table(data$tb_diag,useNA="ifany")                           # TB
table(data$cirrhosis_hep_diag,useNA="ifany")                # Cirrhosis/Hepatitis
table(data$asthma_diag,useNA="ifany")                       # Asthma
table(data$rheum_arthritis_diag,useNA="ifany")              # Rheumatoid arthritis
table(data$family_cancer,useNA="ifany")                     # Family history of cancer
table(data$blood_transfusion,useNA="ifany")                 # Blood transfusion
table(data$blood_donation,useNA="ifany")                    # Blood donation
table(data$children_cat2,useNA="ifany")                     # Children (women)
table(data$early_menarche,useNA="ifany")                    # Early menarche (< 12 years)
table(data$post_menopausal,useNA="ifany")                   # Post menopausal

# Define exposures

data$age_cat5_1<-data$age_cat5
data$age_cat5_2<-data$age_cat5
data$age_cat5_3<-data$age_cat5
data$age_cat5_4<-data$age_cat5
data$birth_cohort_cat4_1<-data$birth_cohort_cat4
data$birth_cohort_cat4_2<-data$birth_cohort_cat4
data$birth_cohort_cat4_3<-data$birth_cohort_cat4

exposures_list<-list()
exposures_list[["Age \n(per year increase)"]]                                     <- "age"
exposures_list[["40 to 50"]]                                <- "age_cat5_1"
exposures_list[["50 to 59"]]                                <- "age_cat5_2"
exposures_list[["60 to 69"]]                                <- "age_cat5_3"
exposures_list[["70 +"]]                                    <- "age_cat5_4"
#exposures_list[["1940 to 1949"]]                            <- "birth_cohort_cat4_1"
#exposures_list[["1950 to 1959"]]                            <- "birth_cohort_cat4_2"
#exposures_list[["1960 +"]]                                  <- "birth_cohort_cat4_3"
exposures_age<-exposures_list

exposures_list<-list()
exposures_list[["Women v Men"]]                                  <- "is_female"
exposures_sex<-exposures_list

exposures_list<-list()
exposures_list[["Urban v Rural"]]                                  <- "region_is_urban"
exposures_region<-exposures_list

exposures_list<-list()
#exposures_list[["Women (Smoker)"]]                         <-"ever_reg_smoker_f"
#exposures_list[["Women (Alcohol drinker)"]]                <-"ever_reg_alcohol_f"
#exposures_list[["Waist circumference < 80 cm "]]           <-"waist_cat2_f"
#exposures_list[["No. of children > 1"]]                    <-"children_cat2"
#exposures_list[["Age at menarche < 12 years"]]             <-"early_menarche"
#exposures_list[["Post menopausal"]]                        <-"post_menopausal"
exposures_female<-exposures_list

exposures_list<-list()
exposures_list[["Regular smoker (men)"]]                           <-"ever_reg_smoker_m"
exposures_list[["Regular alcohol drinker (men)"]]                  <-"ever_reg_alcohol_m"
#exposures_list[["Waist circumference < 94 cm "]]           <-"waist_cat2_m"
exposures_male<-exposures_list

exposures_list<-list()
exposures_list[["Education > 6 years"]]                     <-"education_cat2"
#exposures_list[["Annual household income < \u00A520,000"]] <-"household_income_cat2"
#exposures_list[["Daily tea drinker"]]                       <-"reg_tea_drinker"
exposures_list[["Overweight"]]                               <-"bmi_overweight"
exposures_list[["Poor self-rated health"]]                  <-"srh_poor"
#exposures_list[["Experienced famine"]]                      <-"diet_had_shortage"
#exposures_list[["Spices"]]                                  <-"reg_spice"
exposures_list[["Fresh \nfruit"]]                             <-"reg_fresh_fruit"
#exposures_list[["Preserved vegetables"]]                    <-"reg_preserved_veg"
#exposures_list[["Meat or poultry"]]                         <-"reg_meat_poultry"
#exposures_list[["CHD/Stroke"]]                              <-"chd_stroke_diag"
#exposures_list[["Rheumatic heart disease"]]                 <-"rheum_heart_dis_diag"
exposures_list[["Diabetes"]]                                <-"diabetes_diag"
exposures_list[["Hyper-\ntension"]]                            <-"hypertension_diag"
#exposures_list[["Emphysema/Bronchitis"]]                    <-"emph_bronc_diag"
#exposures_list[["Peptic \nulcer"]]                            <-"peptic_ulcer_diag"
#exposures_list[["TB"]]                                      <-"tb_diag"
exposures_list[["Cirrhosis\n/Hepatitis"]]                     <-"cirrhosis_hep_diag"
#exposures_list[["Asthma"]]                                  <-"asthma_diag"
#exposures_list[["Rheumatic arthritis"]]                     <-"rheum_arthritis_diag"
exposures_list[["Family history of cancer"]]                <-"family_cancer"
#exposures_list[["Blood transfusion"]]                       <-"blood_transfusion"
#exposures_list[["Blood donation"]]                          <-"blood_donation"
exposures_list[["HBsAg+"]]                                   <-"hep_b_diag"
exposures_other<-exposures_list

exposures<-c(exposures_age,exposures_sex,exposures_male,exposures_female,exposures_region,exposures_other)

# Define outcomes

outcomes<-grep("_pp$",names(data),value=TRUE)[c(1:13,16:19)] # Remove high/low prevalence pathogens ("HIV-1","HTLV-1","EBV")

### Logistic model function ####

data_all<-data

model_fun <- function(data_all, exposure, outcome){

  exposure<-unlist(exposure)
  
  # Define population and adjustments
  
  if(exposure %in% exposures_age ){
    adjust<-c("is_female","region_code")
    data<-data_all
  }
  else if(exposure %in% exposures_region ){
    adjust<-c("age","is_female")
    data<-data_all
  }  
  else if(exposure %in% exposures_sex ){
    adjust<-c("age","region_code")
    data<-data_all
  }  
  else if(exposure %in% exposures_female ){
    adjust<-c("age","region_code")
    data<-data_all[data_all$is_female==1,]
  }  
  else if(exposure %in% exposures_male ){
    adjust<-c("age","region_code")
    data<-data_all[data_all$is_female==0,]
  }
  else if(exposure %in% exposures_other){
    adjust<-c("age","is_female","region_code")
    data<-data_all
  }
  
  model<-glm(reformulate(c(exposure,adjust),response=outcome),data=data, family=binomial()) 
  
    n<-as.data.frame(table(data[[outcome]]))[2,2]                                          # Number of cases
    adjust_var<-stringr::str_c(adjust,collapse=" + ")                                      # Create a single adjustment variable

    if(exposure %in% c(exposures_sex,exposures_region,exposures_male,exposures_female,exposures_other,"age","age_cat5_1","birth_cohort_cat4_1")) {
    results<-c(summary(model)$coefficients[2,1:4],n,exposure,outcome,adjust_var)           # Model and summary statistics  
    }
    else if(exposure %in% c("age_cat5_2","birth_cohort_cat4_2")){
    results<-c(summary(model)$coefficients[3,1:4],n,exposure,outcome,adjust_var)           # Model and summary statistics  
    }
    else if(exposure %in% c("age_cat5_3","birth_cohort_cat4_3")){
      results<-c(summary(model)$coefficients[4,1:4],n,exposure,outcome,adjust_var)           # Model and summary statistics  
    }
    else if(exposure %in% c("age_cat5_4")){
      results<-c(summary(model)$coefficients[5,1:4],n,exposure,outcome,adjust_var)           # Model and summary statistics  
    }
}

### Run function ####

  results<-list()
  
  for(exposure in exposures){
    
    results[[exposure]]<-list()
    
    for(outcome in outcomes){  
      {
        
        results[[exposure]][[outcome]] <- model_fun(data_all, exposure, outcome)
        
      }}}
  
tab<-as.data.frame(t(as.data.frame(purrr::list_flatten(results))))
names(tab)<-c("est","stderr","t_stat","p_value","n","exposure","pathogen","adjust")
tab$pathogen<-str_remove_all(tab$pathogen,"_pp")

tab$sig<-ifelse(as.numeric(tab$p_value)<0.05 & as.numeric(tab$est)>=0,1,
                ifelse(as.numeric(tab$p_value)<0.05 & as.numeric(tab$est)<0,2,0))   
tab$sig_pp<-paste0(tab$pathogen,"_",tab$exposure,"_",tab$sig)
sig_plus<-grep("_1$",tab$sig_pp,value=TRUE)
sig_minus<-grep("_2$",tab$sig_pp,value=TRUE)
tab$odds_ratio<-exp(as.numeric(tab$est))
tab$ci_lower<-exp(as.numeric(tab$est) - 1.96*as.numeric(tab$stderr))
tab$ci_upper<-exp(as.numeric(tab$est) + 1.96*as.numeric(tab$stderr))
tab$ci_upper[tab$ci_upper>100]<-100
tab_ci_upper<-as.numeric(tab$ci_upper)
tab$ci_low<-sprintf(tab$ci_lower, fmt = '%#.2f')  
tab$ci_up<-sprintf(tab$ci_upper, fmt = '%#.2f') 
tab$or<-sprintf(tab$odds_ratio, fmt = '%#.2f') 
tab$or_ci<-paste0(tab$or,"\n (",tab$ci_low,", ",tab$ci_up,")") 

tab_long<-tab %>% dplyr::select(pathogen,exposure,or_ci)
tab_wide<-tab_long %>% tidyr::pivot_wider(names_from = exposure, values_from = or_ci)

tab_long_sig<-tab %>% dplyr::select(pathogen,exposure,sig)
tab_wide_sig<-tab_long_sig %>% tidyr::pivot_wider(names_from = exposure, values_from = sig)

# Assign pathogen labels

tab_meta_selected<-tab_meta[!duplicated(tab_meta$pathogen),c("pathogen_label","pathogen")]
tab_meta_selected<-tab_meta_selected %>% filter(pathogen %in% tab_wide$pathogen)
tab_wide<-merge(tab_wide,tab_meta_selected, by="pathogen")
tab_wide_sig<-merge(tab_wide_sig,tab_meta_selected, by="pathogen")


tab_wide<-tab_wide %>% arrange(factor(pathogen, levels = c("hsv1","hsv2","vzv","cmv","hhv6","hhv7",
                                                               "hbv","hcv","bkv","jcv","mcv","hpv16","hpv18",
                                                               "ctrachomatis","hpylori","cburnetii","tgondii")))

tab_wide_sig<-tab_wide_sig %>% arrange(factor(pathogen, levels = c("hsv1","hsv2","vzv","cmv","hhv6","hhv7",
                                                           "hbv","hcv","bkv","jcv","mcv","hpv16","hpv18",
                                                           "ctrachomatis","hpylori","cburnetii","tgondii")))

# Assign exposure labels

new_name<-rownames(as.data.frame(do.call(rbind,exposures)))
old_name<-names(tab_wide)[2:(dim(tab_wide)[2]-1)]
tab_wide<-tab_wide %>% rename_with(~ new_name, all_of(old_name))
tab_wide_sig<-tab_wide_sig %>% rename_with(~ new_name, all_of(old_name))

tab_wide<-tab_wide[,c(20,2,7,10,11,8,9,12,14,19,13,18)]
tab_wide_sig<-tab_wide_sig[,c(20,2,7,10,11,8,9,12,14,19,13,18)]

gap_1<-rep(c(""),17)
gap_2<-rep(c(""),17)

tab_wide<-cbind(tab_wide[,1:5],gap_1,tab_wide[,6:9],gap_2,tab_wide[,10:12])
tab_wide_sig<-cbind(tab_wide_sig[,1:5],gap_1,tab_wide_sig[,6:9],gap_2,tab_wide_sig[,10:12])

# Identify significant associations

tab_wide_sig$n<-1:dim(tab_wide_sig)[1]
exposures<-names(tab_wide_sig)[2:14]

sig_plus<-list()
sig_minus<-list()
for(ii in exposures) {sig_plus[[ii]]<-(tab_wide_sig$n[tab_wide_sig[[ii]]==1])}
for(ii in exposures) {sig_minus[[ii]]<-(tab_wide_sig$n[tab_wide_sig[[ii]]==2])}

### Create flextable ####

col_names<-names(tab_wide)

header <- data.frame(col_keys = col_names,
                             head_1  = rep(c(NA,"Demographic & socio-economic factors",NA,"Lifestyle",NA,"Medical history"),c(1,4,1,4,1,3)),
                             head_2 = c(NA,names(tab_wide)[2:5],NA,names(tab_wide)[7:10],NA,names(tab_wide)[12:14]))


ft<-flextable(tab_wide, col_keys=c(col_names)) %>%
    set_header_df(mapping = header, key = "col_keys") %>%
    merge_h(part = "header") %>%
    width(j=1,width=0.9) %>%
    width(j=c(2:5,7:10,12:14),width=0.9) %>%
    width(j=c(6,11),width=0.1) %>%
  bold(i=as.numeric(unlist(sig_plus[1])),j=2,part="body") %>% 
  bold(i=as.numeric(unlist(sig_minus[1])),j=2,part="body") %>% 
  bold(i=as.numeric(unlist(sig_plus[2])),j=3,part="body") %>% 
  bold(i=as.numeric(unlist(sig_minus[2])),j=3,part="body") %>% 
  bold(i=as.numeric(unlist(sig_plus[3])),j=4,part="body") %>% 
  bold(i=as.numeric(unlist(sig_minus[3])),j=4,part="body") %>% 
  bold(i=as.numeric(unlist(sig_plus[4])),j=5,part="body") %>% 
  bold(i=as.numeric(unlist(sig_minus[4])),j=5,part="body") %>% 
  bold(i=as.numeric(unlist(sig_plus[5])),j=6,part="body") %>% 
  bold(i=as.numeric(unlist(sig_minus[5])),j=6,part="body") %>% 
  bold(i=as.numeric(unlist(sig_plus[6])),j=7,part="body") %>% 
  bold(i=as.numeric(unlist(sig_minus[6])),j=7,part="body") %>% 
  bold(i=as.numeric(unlist(sig_plus[7])),j=8,part="body") %>% 
  bold(i=as.numeric(unlist(sig_minus[7])),j=8,part="body") %>% 
  bold(i=as.numeric(unlist(sig_plus[8])),j=9,part="body") %>% 
  bold(i=as.numeric(unlist(sig_minus[8])),j=9,part="body") %>% 
  bold(i=as.numeric(unlist(sig_plus[9])),j=10,part="body") %>% 
  bold(i=as.numeric(unlist(sig_minus[9])),j=10,part="body") %>% 
  bold(i=as.numeric(unlist(sig_plus[10])),j=11,part="body") %>% 
  bold(i=as.numeric(unlist(sig_minus[10])),j=11,part="body") %>% 
  bold(i=as.numeric(unlist(sig_plus[11])),j=12,part="body") %>% 
  bold(i=as.numeric(unlist(sig_minus[11])),j=12,part="body") %>% 
  bold(i=as.numeric(unlist(sig_plus[12])),j=13,part="body") %>% 
  bold(i=as.numeric(unlist(sig_minus[12])),j=13,part="body") %>% 
  bold(i=as.numeric(unlist(sig_plus[13])),j=14,part="body") %>% 
  bold(i=as.numeric(unlist(sig_minus[13])),j=14,part="body") %>%   
  color(i=as.numeric(unlist(sig_plus[1])),j=2,color="red",part="body") %>% 
  color(i=as.numeric(unlist(sig_minus[1])),j=2,color="blue",part="body") %>% 
  color(i=as.numeric(unlist(sig_plus[2])),j=3,color="red",part="body") %>% 
  color(i=as.numeric(unlist(sig_minus[2])),j=3,color="blue",part="body") %>% 
  color(i=as.numeric(unlist(sig_plus[3])),j=4,color="red",part="body") %>% 
  color(i=as.numeric(unlist(sig_minus[3])),j=4,color="blue",part="body") %>% 
  color(i=as.numeric(unlist(sig_plus[4])),j=5,color="red",part="body") %>% 
  color(i=as.numeric(unlist(sig_minus[4])),j=5,color="blue",part="body") %>% 
  color(i=as.numeric(unlist(sig_plus[5])),j=6,color="red",part="body") %>% 
  color(i=as.numeric(unlist(sig_minus[5])),j=6,color="blue",part="body") %>% 
  color(i=as.numeric(unlist(sig_plus[6])),j=7,color="red",part="body") %>% 
  color(i=as.numeric(unlist(sig_minus[6])),j=7,color="blue",part="body") %>% 
  color(i=as.numeric(unlist(sig_plus[7])),j=8,color="red",part="body") %>% 
  color(i=as.numeric(unlist(sig_minus[7])),j=8,color="blue",part="body") %>% 
  color(i=as.numeric(unlist(sig_plus[8])),j=9,color="red",part="body") %>% 
  color(i=as.numeric(unlist(sig_minus[8])),j=9,color="blue",part="body") %>% 
  color(i=as.numeric(unlist(sig_plus[9])),j=10,color="red",part="body") %>% 
  color(i=as.numeric(unlist(sig_minus[9])),j=10,color="blue",part="body") %>% 
  color(i=as.numeric(unlist(sig_plus[10])),j=11,color="red",part="body") %>% 
  color(i=as.numeric(unlist(sig_minus[10])),j=11,color="blue",part="body") %>% 
  color(i=as.numeric(unlist(sig_plus[11])),j=12,color="red",part="body") %>% 
  color(i=as.numeric(unlist(sig_minus[11])),j=12,color="blue",part="body") %>% 
  color(i=as.numeric(unlist(sig_plus[12])),j=13,color="red",part="body") %>% 
  color(i=as.numeric(unlist(sig_minus[12])),j=13,color="blue",part="body") %>% 
  color(i=as.numeric(unlist(sig_plus[13])),j=14,color="red",part="body") %>% 
  color(i=as.numeric(unlist(sig_minus[13])),j=14,color="blue",part="body") %>% 
   italic(i=c(14:17),j=1, part="body") %>% 
    italic(i=c(14:17),j=1, part="body") %>% 
    bold(i=NULL,j=NULL, part="header") %>% 
    bold(i=NULL,j=1, part="body") %>% 
    padding(padding.left = 4,padding.right = 4,padding.top = 2,padding.bottom = 2, part = "all") %>% 
    hline_top(border = fp_border(width =1.5, color = "black"), part = "header" ) %>%  
    hline(i=1,j=c(2:5,7:10,12:14),border = fp_border(width = 1, color = "black"), part = "header" ) %>%
    hline_bottom(border = fp_border(width =1.5, color = "black"), part = "header" ) %>%
    align(i=NULL,j=NULL,align="center",part="header") %>%
    align(i=NULL,j=c(2:14),align="center",part="body") %>%
    add_footer_lines("Note: ORs: odds ratios; Adjusted for age (10-year age groups), sex and region (10 regions) where appropriate; Compared Yes vs. No groups for each specific demographic & socio-economic, lifestyle factor and medical history;\nBold values denote statistical significance at the p<0.05 level, positive in red and negative in blue") %>%
  #  line_spacing(space = 0.78, part = "header") %>%  
   # line_spacing(space = 0.78, part = "body") %>% 
    #line_spacing(space = 0.1, part = "footer") %>%
    fontsize(i = NULL, j = NULL, size = 8,part="body") %>%
    fontsize(i = NULL, j = NULL, size = 8,part="header") %>% 
    fontsize(i = NULL, j = NULL, size = 7,part="footer")

sect_properties <- prop_section(page_size = page_size(orient = "landscape", width = 12, height = 8),
                                type = "continuous", page_margins = page_mar())

save_as_docx(ft,path="K:/kadoorie/Staff_Folders/JonathanC/Projects/serology/Tables_figures/Tables_figures_10/Table_3.docx",
             pr_section=sect_properties)

 









