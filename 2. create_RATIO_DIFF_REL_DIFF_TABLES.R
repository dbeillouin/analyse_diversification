# Damien BEillouin
# Méta-synthèse sur la divsersification des cultures: code d'analyse
# travail initial: Octobre 2018, reprise Février 2019

### FOCUS SUR LES DONNEES DE RATIO #######

# Initialisation -----------------------------------------------------------------------------------------------

  # Load packages and scripts
    x<-c("magrittr", "tidyverse", "ggplot2","lme4","metafor","broom","ggpubr","reshape2","tidyr","dplyr")
    lapply(x, require, character.only = TRUE)
    source('~/Documents/analyse_diversification/analyse_diversification/Load_data.R', echo=TRUE)
    
  # set wd
   setwd("~/Documents/shinyapp/")


# Select main  effect-sizes ------------------------------------------------------------------------
    
    ES_TMP<-ES                                 %>%  
      filter(Type_Of_Analysis=="effect size")  %>%  # We do not analyse regressions, boxplots,...
      filter(KEEP=="YES")                      %>%  # We analyse only "main" effect sizes (not sub-scenarios)
      select(ID,contains("Scenario_"),DIVERS,code_Type_Output,
             Sub_Cat_Output,Sub_Sub_Cat_Output,Metric_Output,
             Nb_Paired_Data,ES_lower_CI, ES, ES_upper_CI,
             contains("Reg"),contains("Treat_"))%>% # select some columns
      mutate(UN= 1)                             %>%
      group_by(ID)                              %>%
      mutate(tid= cumsum(UN))                   %>%
      mutate(vi= FUN_var_from_CI(mean=ES, CI=ES_upper_CI))
    
    ## Check Data
       ggplot(aes(ES,vi),data=ES_TMP)+geom_point()+ theme_bw()+ facet_wrap(~Metric_Output,scales="free")
      # Etrange pour les differences, l'augmentation de variance!!!
    
    
# Split ratio, relative differences and differences ----------------------------------------------

    # RATIO
      RATIO<- ES_TMP %>%  filter(Metric_Output=="log(Ratio)")
      RATIO$Scenario_ID<-as.factor(RATIO$Scenario_ID)
      
    ##DIFFERENCE RELATIVES
     DIFF_REL<-ES %>%  filter(Metric_Output=="cohen's d")
    
    ####DIFFERENCE
     DIFF<-ES %>% filter(Metric_Output=="unit")
    
        # pPour Angus, on transforme les régression en mean différence, en se basant sur la moyenne des valeurs observées
        DIFF  <-DIFF %>% mutate(ES= ifelse(ID=="Angus_2015" & Scenario_ID ==1,
                                           Reg_Intercept+ 3.46* Reg_Slope
                                           ,ES))
        DIFF  <-DIFF %>% mutate(ES_lower_CI= ifelse(ID=="Angus_2015" & Scenario_ID ==1,
                                          (Reg_Intercept+ 3.46* Reg_Slope)-(Reg_SE/2*1.96),
                                          ES_lower_CI))
        DIFF  <-DIFF %>% mutate(ES_upper_CI= ifelse(ID=="Angus_2015" & Scenario_ID ==1,
                                          (Reg_Intercept+ 3.46* Reg_Slope)-(Reg_SE/2*1.96),
                                          ES_upper_CI))
        
        
        DIFF  <-DIFF %>% mutate(ES= ifelse(ID=="Angus_2015" & Scenario_ID ==2,Reg_Intercept+ 3.62* Reg_Slope,ES))
        DIFF  <-DIFF %>% mutate(ES_lower_CI= ifelse(ID=="Angus_2015" & Scenario_ID ==2,(Reg_Intercept+ 3.62* Reg_Slope)-(Reg_SE/2*1.96),ES_lower_CI))
        DIFF  <-DIFF %>% mutate(ES_upper_CI= ifelse(ID=="Angus_2015" & Scenario_ID ==2,(Reg_Intercept+ 3.62* Reg_Slope)-(Reg_SE/2*1.96),ES_upper_CI))
        
        DIFF  <-DIFF %>% mutate(ES= ifelse(ID=="Angus_2015" & Scenario_ID ==3,Reg_Intercept+ 3.30* Reg_Slope,ES))
        DIFF  <-DIFF %>% mutate(ES_lower_CI= ifelse(ID=="Angus_2015" & Scenario_ID ==3,(Reg_Intercept+ 3.30* Reg_Slope)-(Reg_SE/2*1.96),ES_lower_CI))
        DIFF  <-DIFF %>% mutate(ES_upper_CI= ifelse(ID=="Angus_2015" & Scenario_ID ==3,(Reg_Intercept+ 3.30* Reg_Slope)-(Reg_SE/2*1.96),ES_upper_CI))
        
        DIFF  <-DIFF %>% mutate(ES= ifelse(ID=="Angus_2015" & Scenario_ID ==4,Reg_Intercept+ 3.14* Reg_Slope,ES))
        DIFF  <-DIFF %>% mutate(ES_lower_CI= ifelse(ID=="Angus_2015" & Scenario_ID ==4,(Reg_Intercept+ 3.14* Reg_Slope)-(Reg_SE/2*1.96),ES_lower_CI))
        DIFF  <-DIFF %>% mutate(ES_upper_CI= ifelse(ID=="Angus_2015" & Scenario_ID ==4,(Reg_Intercept+ 3.14* Reg_Slope)-(Reg_SE/2*1.96),ES_upper_CI))
        
        DIFF  <-DIFF %>% mutate(ES= ifelse(ID=="Angus_2015" & Scenario_ID ==5,Reg_Intercept+ 3.43* Reg_Slope,ES))
        DIFF  <-DIFF %>% mutate(ES_lower_CI= ifelse(ID=="Angus_2015" & Scenario_ID ==5,(Reg_Intercept+ 3.43* Reg_Slope)-(Reg_SE/2*1.96),ES_lower_CI))
        DIFF  <-DIFF %>% mutate(ES_upper_CI= ifelse(ID=="Angus_2015" & Scenario_ID ==5,(Reg_Intercept+ 3.43* Reg_Slope)-(Reg_SE/2*1.96),ES_upper_CI))
        
        DIFF  <-DIFF %>% mutate(ES= ifelse(ID=="Angus_2015" & Scenario_ID %in% c(6,7,8),Reg_Intercept+ 2.34* Reg_Slope,ES))
        DIFF  <-DIFF %>% mutate(ES_lower_CI= ifelse(ID=="Angus_2015" & Scenario_ID %in% c(6,7,8),(Reg_Intercept+ 2.34* Reg_Slope)-(Reg_SE/2*1.96),ES_lower_CI))
        DIFF  <-DIFF %>% mutate(ES_upper_CI= ifelse(ID=="Angus_2015" & Scenario_ID %in% c(6,7,8),(Reg_Intercept+ 2.34* Reg_Slope)-(Reg_SE/2*1.96),ES_upper_CI))
        
        DIFF  <-DIFF %>% mutate(ES= ifelse(ID=="Angus_2015" & Scenario_ID %in% c(9,10),Reg_Intercept+ 2.34* Reg_Slope,ES))
        DIFF  <-DIFF %>% mutate(ES_lower_CI= ifelse(ID=="Angus_2015" & Scenario_ID %in% c(9,10),(Reg_Intercept+ 1.68* Reg_Slope)-(Reg_SE/2*1.96),ES_lower_CI))
        DIFF  <-DIFF %>% mutate(ES_upper_CI= ifelse(ID=="Angus_2015" & Scenario_ID %in% c(9,10),(Reg_Intercept+ 1.68* Reg_Slope)-(Reg_SE/2*1.96),ES_upper_CI))
        
        # on convertit Himmelsetein de Kg ha-1 à en T.ha-1
        DIFF  <-DIFF %>% mutate(ES= ifelse(ID=="Himmelstein_2017",ES/1000,ES))
        DIFF  <-DIFF %>% mutate(ES_lower_CI= ifelse(ID=="Himmelstein_2017",ES_lower_CI/1000,ES_lower_CI))
        DIFF  <-DIFF %>% mutate(ES_upper_CI= ifelse(ID=="Himmelstein_2017",ES_upper_CI/1000,ES_upper_CI))
        
        


# Final table for analyses ----------------------------------------------------------------------------------------

    RATIO2<-RATIO %>%select(ID, ES, vi, Sub_Cat_Output,DIVERS, Scenario_ID) %>% ungroup() %>%
      mutate(NUM= as.character(1: length(ID)))
    
    RATIO2<-RATIO2[!RATIO2$vi==0,]
    RATIO2<-RATIO2[!is.na(RATIO2$ES),]
    RATIO2$ID<-as.factor(as.character(RATIO2$ID))
