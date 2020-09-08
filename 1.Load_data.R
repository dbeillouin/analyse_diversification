# Damien Beillouin
# Méta-synthèse sur la divsersification des cultures: code d'analyse
# travail initial: Octobre 2018, reprise Février 2019, reprise AVRIl 2020

### SCRIPT POUR CHARGER LES DONNEES #######

# Initialisation -----------------------------------------------------------------------------------------------

# Load packages
  x<-c("rlang", "magrittr", "tidyverse", "ggplot2","lme4","metafor","broom",
       "ggpubr","reshape2","tidyr","dplyr", "readxl")
  lapply(x, require, character.only = TRUE)

# set wd
  setwd("~/Documents/analyse_diversification")

# Load functions
    
    # Fisher'Z -> Cohen s d
    FUN_Z_to_r<-function(Z){
      r=(exp(2*Z)-1)/(exp(2*Z)+1)
      return(r)}
    
    FUN_r_to_d<-function(r){
      d= 2*r/(sqrt(1-r*r))
      return(d)}
    
    #http://handbook-5-1.cochrane.org/chapter_7/7_7_3_3_obtaining_standard_deviations_from_standard_errors.htm
    FUN_var_from_CI<-function(mean, CI){
      Var= (abs(CI-mean)/1.96)^2
      return(Var)
    }
    
    FUN_var_r_to_d<-function(Var_r, r){
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

    
    ## from LOR to d
    FUN_LOR_to_d<-function(LOR){
      d= LOR*sqrt(3)/3.14159
      return(d)
    }
    
    FUN_var_LOR_to_d<-function(Var_LOR){
      Var_d= Var_LOR *3/((3.14159)^2)
      return(Var_d)
    }

# Load data ---------------------------------------------------------------------------------------------------
  
  # Primary studies and description of the MA
    IND_STU  <-       read_excel("~/Documents/divers2020/PS_MAI.xlsx", sheet="PS_MAI")
   
    RR<- IND_STU %>% group_by(Species) %>% summarise(n= n_distinct(DOI)) %>% filter(n>5)
      length(unique(RR$Country))
     Des_MA <- read_excel("~/Documents/Meta_synthesis_GLOBAL/Description_MA_2019_2020.xlsx", sheet = "Feuil1") 
      
      
  ##### Quality of the MA######
    
    
    QUAL <- read_excel("~/Documents/Meta_synthesis_GLOBAL/Quality_2019SAUV_V3.xlsx", sheet = "Feuil1") %>%
      select("ID", "Score",starts_with("Criterion"))                                           %>%
      setNames(.,c("ID","Score",paste0("Item_",seq(1:20))))                                   %>%
      mutate_at(.vars = vars(contains("Item")), .funs = funs(ifelse(.=="YES", 1, 0)))
    
  ####### Effect size#########
    
    ES<- read_excel("~/Documents/Meta_synthesis_GLOBAL/Effect_sizes2019_FEV2020V4.xlsx",sheet = "Feuil1")   %>%
              mutate (ES=as.numeric(gsub(",", ".", ES)),
              ES_upper_95 = as.numeric(gsub(",", ".", ES_upper_95)),
              ES_lower_95 = as.numeric(gsub(",", ".", ES_lower_95)),
              ES_SE       = as.numeric(gsub(",", ".", ES_SE)),
              ES_minus_SE = as.numeric(as.character(ES_minus_SE)), 
              ES_Plus_SE  = as.numeric(as.character(ES_Plus_SE)),
              ES_SE       = as.numeric(as.character(ES_SE)), 
              ES_P_Value  = as.numeric(as.character(ES_P_Value)), 
              Nb_Paired_Data  = as.numeric(as.character(Nb_Paired_Data)),
              Reg_Slope        = as.numeric(as.character(Reg_Slope)),
              Reg_Intercept        = as.numeric(as.character(Reg_Intercept)),
              Reg_Slope_SE        = as.numeric(as.character(Reg_Slope_SE)),
              Sub_Cat_Output   = tolower(Sub_Cat_Output))
    
    dim(ES)
    library(plyr)
              ES <- ES  %>% mutate(Sub_Cat_Output = revalue(Sub_Cat_Output,
                                                c("grain yield" = "Production", 
                                                 "ler"  = "Production",
                                   "grass production"  = "Production")))
              detach("package:plyr", unload = TRUE)
    
    TMP<- ES %>% filter(Metric_Output== "Ratio")
    boxplot(TMP$ES)
    
  # On calcule la moyenne des SE (pour les données qui n'ont pas de CI)
    # Pour les données avec uniquement un SE fourni, on calcule l'IC
    ES$SE2          <- (abs(ES$ES- ES$ES_minus_SE)+ abs(ES$ES_Plus_SE - ES$ES))/2
    ES          %<>%
      mutate(ES_lower_95= ifelse(is.na(ES_lower_95),ES-1.96*SE2,ES_lower_95), 
             ES_upper_95= ifelse(is.na(ES_upper_95),ES-1.96*SE2,ES_upper_95), 
             ES_lower_95= ifelse(is.na(ES_lower_95),ES-1.96*ES_SE,ES_lower_95), 
             ES_upper_95= ifelse(is.na(ES_upper_95),ES+1.96*ES_SE,ES_upper_95))
  VERIF <- ES %>% filter(DIVERS=="Rotation")  %>% 
          filter(Sub_Cat_Output== "products quality")
  
    TMP<- ES %>% filter(Metric_Output== "Ratio")
    boxplot(TMP$ES)
    
  ## calculate CI from P value 
    ES <- ES %>% 
      mutate(ES_lower_95= 
               ifelse(is.na(ES_lower_95) & Metric_Output %in% c('log(Ratio)'),
                     ES -1.96*ES/abs(-0.862 + sqrt(0.743 - 2.404*log(ES_P_Value))),
                     ES_lower_95), 
             ES_upper_95= 
                ifelse(is.na(ES_upper_95) & Metric_Output %in% c('log(Ratio)'),
                      ES +1.96*ES/abs(-0.862 + sqrt(0.743 - 2.404*log(ES_P_Value))),
                     ES_upper_95))
  
    TMP<- ES %>% filter(Metric_Output== "Ratio")
    boxplot(TMP$ES)
    
  
  # Relative différence -> Ratio
    ES <- ES%>% mutate(ES          = 
                         ifelse(Metric_Output %in% c("Percent change"), 
                                ES+1,
                                ES),
                       ES_lower_95 = 
                         ifelse(Metric_Output %in% c("Percent change"),
                                ES_lower_95+1, 
                                ES_lower_95),
                       ES_upper_95 =
                         ifelse(Metric_Output %in% c("Percent change"),
                                ES_upper_95+1,
                                ES_upper_95),
                       ES_SE       =
                         ifelse(Metric_Output %in% c("Percent change"),
                                ES_SE+1,
                                ES_SE)     ) %>%
      mutate(Metric_Output = forcats::fct_recode(Metric_Output, "Ratio" = "Percent change"))
  
    TMP<- ES %>% filter(Metric_Output== "Ratio")
    boxplot(TMP$ES)
    TMP<- TMP %>% filter(ES>100)
    
    
    `%notin%` <- Negate(`%in%`)
    TTT<- ES %>% filter(ID=="MA_A2")
    # Ratio -> log (Ratio)
    ES <- ES%>% mutate(ES          =
                         ifelse(Metric_Output %in% c("Ratio") & Log_scale %notin% c("YES"), 
                                log(ES),
                                ES),
                       ES_lower_95 =
                         ifelse(Metric_Output %in% c("Ratio")& Log_scale %notin% c("YES"), 
                                log(ES_lower_95),
                                ES_lower_95),
                       ES_upper_95 = 
                         ifelse(Metric_Output %in% c("Ratio")& Log_scale %notin% c("YES"),
                                log(ES_upper_95),
                                ES_upper_95)) 
  
    ## A faire après la conversion en log ratio!!!
    ## on change l'ordre des ES si nécéssaire (= par exemple pour la bulk densité, des faibles valeurs c'est mlieux)
    
    ES <- ES          %>%
      mutate(ES= ifelse(is.na(Change_order_ratio),ES,-ES))
    ES <- ES          %>%
      mutate(ES_lower_95 = ifelse(is.na(Change_order_ratio),ES_lower_95,-ES_lower_95))
    ES <- ES          %>%
      mutate(ES_upper_95 = ifelse(is.na(Change_order_ratio),ES_upper_95,-ES_upper_95))
    
    
    
    TMP<- ES %>% filter(Metric_Output== "Ratio")
    boxplot(TMP$ES)
    TMP<- TMP %>% filter(ES>4)
  
    #########------------------------------##################@
    # Attention à calculer les CI avant la moyenne....
  
    #1. conversion from Fisher to cohen's
    
    ES  <- ES  %>%
      mutate(ES_lower_95= 
               ifelse(Metric_Output=="Fisher's Z",
                      ES - FUN_CI_from_var(
                      Var= FUN_var_r_to_d(Var_r = 
                                          FUN_var_from_CI(mean = ES,CI=ES_lower_95),
                      r = ES )),
                      ES_lower_95))
    
    ES <- ES %>% 
      mutate(ES_upper_95= 
               ifelse(Metric_Output=="Fisher's Z",
                      ES+FUN_CI_from_var(
                        Var= FUN_var_r_to_d(Var_r = 
                                         FUN_var_from_CI(mean=ES,CI=ES_upper_95),
                      r     = ES )),
                      ES_upper_95))
    
    ES <- ES %>%
      mutate(ES= 
               ifelse(Metric_Output=="Fisher's Z" ,
                      FUN_r_to_d(FUN_Z_to_r(ES)),
                      ES))
    ES<- ES%>% mutate(Metric_Output = fct_recode(Metric_Output, 
                                                 "Cohen's d" = "Fisher's Z"))
    
    
    ## 2. Conversion Hedge to cohen
    
    ES <- ES %>% 
      mutate(ES= 
               ifelse(Metric_Output %in% c('Hedges g','Hedge s d'),
                      ES/FUN_J(df=2*Nb_Paired_Data-2),
                      ES))
    
    ES <- ES %>%
      mutate(ES_lower_95=
               ifelse(Metric_Output %in% c('Hedges g','Hedge s d'),
               ES - FUN_CI_from_var(FUN_Var_from_CI(abs(ES-ES_lower_95))/FUN_J(df=2*Nb_Paired_Data-2)^2),
               ES_lower_95))
    
    ES <- ES %>%
      mutate(ES_upper_95=
               ifelse(Metric_Output %in% c('Hedges g','Hedge s d'),
                ES + FUN_CI_from_var(FUN_Var_from_CI(abs(ES_upper_95-ES))/FUN_J(df=2*Nb_Paired_Data-2)^2),
                ES_upper_95))
    
    ES<- ES%>% mutate(Metric_Output = fct_recode(Metric_Output, "Cohen's d" = "Hedge's g"))
    ES<- ES%>% mutate(Metric_Output = fct_recode(Metric_Output, "Cohen's d" = "Hedge's d"))
    ES$Metric_Output<-factor(ES$Metric_Output)
    
    
    #3. conversion LOR to cohen's 
    
    ES <- ES %>%
      mutate(ES= 
               ifelse(Metric_Output=="Odds Ratio" ,
                      FUN_LOR_to_d(ES),
                      ES))
    
    ES <- ES %>%
      mutate(ES_lower_95=
               ifelse(Metric_Output %in% c('Odds Ratio'),
                      ES - FUN_CI_from_var(
                        FUN_var_LOR_to_d(FUN_Var_from_CI(abs(ES-ES_lower_95)))),
                      ES_lower_95))
    
    ES <- ES %>%
      mutate(ES_upper_95=
               ifelse(Metric_Output %in% c('Odds Ratio'),
                      ES - FUN_CI_from_var(
                        FUN_var_LOR_to_d(FUN_Var_from_CI(abs(ES-ES_upper_95)))),
                      ES_upper_95))
    
    ES<- ES%>% mutate(Metric_Output = fct_recode(Metric_Output, "Cohen's d" = "Odds Ratio"))
    ES$Metric_Output<-factor(ES$Metric_Output)
    
    
  #PART2: mean of the CI 
  ## On moyenne les deux intervalles de confiances et on attribue la moyenne valeur pour les bornes
    # ES$CI_moy<- (abs(ES$ES- ES$ES_lower_95)+ abs(ES$ES_upper_95 - ES$ES))/2
    # ES$ES_lower_95<- ES$ES - ES$CI_moy
    # ES$ES_upper_95<- ES$ES + ES$CI_moy
    ESSAI <- ES %>% filter(ID=="MA2020_V1")
