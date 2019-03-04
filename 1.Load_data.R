# Damien Beillouin
# Méta-synthèse sur la divsersification des cultures: code d'analyse
# travail initial: Octobre 2018, reprise Février 2019

### SCRIPT POUR CHARGER LES DONNEES #######

# Initialisation -----------------------------------------------------------------------------------------------

# Load packages
  x<-c("magrittr", "tidyverse", "ggplot2","lme4","metafor","broom","ggpubr","reshape2","tidyr","dplyr")
  lapply(x, require, character.only = TRUE)

# set wd
  setwd("~/Documents/shinyapp")

# Load functions
    
    # Fisher'Z -> Cohen s d
    FUN_Z_to_r<-function(Z){
      r=(exp(2*Z)-1)/(exp(2*Z)+1)
      return(r)}
    
    FUN_r_to_d<-function(r){
      d= 2*r/(sqrt(1-r*r))
      return(d)}
    
    FUN_var_from_CI<-function(mean, CI){
      Var= (abs(CI-mean)/1.96)^2
      return(Var)
    }
    
    FUN_var_d<-function(Var_r, r){
      Var_d= 4*Var_r/((1- r*r)^3)
      return(Var_d)
    }
    
    FUN_CI_from_var<-function(Var){
      CI= sqrt(Var)/2*1.96
      return(CI)
    }
    
    # Hedges’ g -> Cohen s d
    #To convert from d to Hedges’ g we use a correction factor, which is
    #called J. Hedges (1981) gives the exact formula for J, but in common practice
    #researchers use an approximation
    FUN_J<-function(df){
      1-3/(4*df-1)
    }
    FUN_Var_from_CI<-function(CI){
      Var= (CI/1.96*2)^2
      return(CI)
    }


# Load data ---------------------------------------------------------------------------------------------------
  
  # Primary studies and description of the MA
    IND_STU  <- read.csv("~/Documents/shinyapp/FINAL/Primary_Studies.csv", sep=";")
    Des_MA   <- read.csv("~/Documents/shinyapp/FINAL/Description_Meta-Analyses.csv", sep=";", comment.char="#")     IND_STU  <-merge(IND_STU, Des_MA[,c("ID","Diversification_types")],by="ID")
      
  ##### Quality of the MA######
    QUAL <- read.csv("~/Documents/shinyapp/FINAL/Quality.csv", sep=";", comment.char="#")  %>%
      select("ID", "X",starts_with("Criterion"))                                           %>%
      setNames(.,c("ID","MA",paste0("Item_",seq(1:20))))                                   %>%
      mutate_at(.vars = vars(contains("Item")), .funs = funs(ifelse(.=="YES", 1, 0)))
    
  ####### Effect size#########
    ES<-read.csv("~/Documents/shinyapp/FINAL/Effect_Size.csv", sep=";")    %>%
      mutate (ES=as.numeric(gsub(",", ".", ES)),
              ES_upper_CI = as.numeric(gsub(",", ".", ES_upper_CI)),
              ES_lower_CI = as.numeric(gsub(",", ".", ES_lower_CI)),
              ES_SE       = as.numeric(gsub(",", ".", ES_SE)))

# data transformation ----------------------------------------------------------------------------------------
  
    ES$ES_minus_SE    <- as.numeric(as.character(ES$ES_minus_SE))
    ES$ES_Plus_SE     <- as.numeric(as.character(ES$ES_Plus_SE))
    ES$ES_SE          <- as.numeric(as.character(ES$ES_SE))
    ES$ES_P_Value     <- as.numeric(as.character(ES$ES_P_Value))
    ES$Nb_Paired_Data <- as.numeric(as.character(ES$Nb_Paired_Data))
    
  
  # On calcule la moyenne des SE (pour les données qui n'ont pas de CI)
    ES$SE2          <- (abs(ES$ES- ES$ES_minus_SE)+ abs(ES$ES_Plus_SE - ES$ES))/2
    ES              <- ES %>% mutate(ES_lower_CI= ifelse(is.na(ES_lower_CI),ES-1.96*SE2,ES_lower_CI))
    ES              <- ES %>% mutate(ES_upper_CI= ifelse(is.na(ES_upper_CI),ES-1.96*SE2,ES_upper_CI))
  
  # Pour les données avec uniquement un SE fourni, on calcule l'IC
    ES            <-ES %>% mutate(ES_lower_CI= ifelse(is.na(ES_lower_CI),ES-1.96*ES_SE,ES_lower_CI))
    ES            <-ES %>% mutate(ES_upper_CI= ifelse(is.na(ES_upper_CI),ES+1.96*ES_SE,ES_upper_CI))
  
  
  ## calculate CI from P value 
    ES <- ES %>% mutate(ES_lower_CI= 
                          ifelse(is.na(ES_lower_CI) & Metric_Output %in% c('log(Ratio)'),
                                 ES -1.96*ES/abs(-0.862 + sqrt(0.743 - 2.404*log(ES_P_Value))),
                                 ES_lower_CI))
    ES <- ES %>% mutate(ES_upper_CI= 
                          ifelse(is.na(ES_upper_CI) & Metric_Output %in% c('log(Ratio)'),
                                 ES +1.96*ES/abs(-0.862 + sqrt(0.743 - 2.404*log(ES_P_Value))),
                                 ES_upper_CI))
  
    # Relative différence -> Ratio
    ES <- ES%>% mutate(ES          = 
                         ifelse(Metric_Output %in% c("Relative difference"), 
                                ES+1,
                                ES),
                       ES_lower_CI = 
                         ifelse(Metric_Output %in% c("Relative difference"),
                                ES_lower_CI+1, 
                                ES_lower_CI),
                       ES_upper_CI =
                         ifelse(Metric_Output %in% c("Relative difference"),
                                ES_upper_CI+1,
                                ES_upper_CI),
                       ES_SE       =
                         ifelse(Metric_Output %in% c("Relative difference"),
                                ES_SE+1,
                                ES_SE)     ) %>%
      mutate(Metric_Output = forcats::fct_recode(Metric_Output, "Ratio" = "Relative difference"))
  
    # Ratio -> log (Ratio)
    ES <- ES%>% mutate(ES          =
                         ifelse(Metric_Output %in% c("Ratio") & !Transformation %in% c("log", "ln"), 
                                log(ES),
                                ES),
                       ES_lower_CI =
                         ifelse(Metric_Output %in% c("Ratio")& !Transformation %in% c("log", "ln"), 
                                log(ES_lower_CI),
                                ES_lower_CI),
                       ES_upper_CI = 
                         ifelse(Metric_Output %in% c("Ratio")& !Transformation %in% c("log", "ln"),
                                log(ES_upper_CI),
                                ES_upper_CI)) %>%
      mutate(Metric_Output = forcats::fct_recode(Metric_Output, "log(Ratio)" = "Ratio"))
  
  
    # Attention à calculer les CI avant la moyenne....
    ES <- ES %>%
      mutate(ES_lower_CI= 
               ifelse(Metric_Output=="Fisher's Z",
                      ES - FUN_CI_from_var(Var= FUN_var_d(Var_r = FUN_var_from_CI(mean = ES,CI=ES_lower_CI),
                      r = ES )),
                      ES_lower_CI))
    
    ES <- ES %>% 
      mutate(ES_lower_CI= 
               ifelse(Metric_Output=="Fisher's Z",
                      ES+FUN_CI_from_var(Var= FUN_var_d(Var_r = FUN_var_from_CI(mean=ES,CI=ES_upper_CI),
                      r     = ES )),
                     ES_lower_CI))
    
    ES <- ES %>%
      mutate(ES= 
               ifelse(Metric_Output=="Fisher's Z" ,
                      FUN_r_to_d(FUN_Z_to_r(ES)),
                      ES))
    #ES<- ES%>% mutate(Metric_Output = fct_recode(Metric_Output, "Cohen's d" = "Fisher's Z"))
    
    ES <- ES %>% 
      mutate(ES= 
               ifelse(Metric_Output %in% c('Hedges g','Hedge s d'),
                      ES/FUN_J(df=2*Nb_Paired_Data-2),
                      ES))
    
    ES <- ES %>%
      mutate(ES_lower_CI=
               ifelse(Metric_Output %in% c('Hedges g','Hedge s d'),
               ES - FUN_CI_from_var(FUN_Var_from_CI(abs(ES-ES_lower_CI))/FUN_J(df=2*Nb_Paired_Data-2)^2),
               ES_lower_CI))
    
    ES <- ES %>%
      mutate(ES_upper_CI=
               ifelse(Metric_Output %in% c('Hedges g','Hedge s d'),
                ES + FUN_CI_from_var(FUN_Var_from_CI(abs(ES_upper_CI-ES))/FUN_J(df=2*Nb_Paired_Data-2)^2),
                ES_upper_CI))
    
    #ES<- ES%>% mutate(Metric_Output = fct_recode(Metric_Output, "Cohen's d" = "Hedge's g"))
    #ES<- ES%>% mutate(Metric_Output = fct_recode(Metric_Output, "Cohen's d" = "Hedge's d"))
    ES$Metric_Output<-factor(ES$Metric_Output)
    
  
  #PART2: mean of the CI 
  ## On moyenne les deux intervalles de confiances et on attribue la moyenne valeur pour les bornes
    ES$CI_moy<- (abs(ES$ES- ES$ES_lower_CI)+ abs(ES$ES_upper_CI - ES$ES))/2
    ES$ES_lower_CI<- ES$ES - ES$CI_moy
    ES$ES_upper_CI<- ES$ES + ES$CI_moy

