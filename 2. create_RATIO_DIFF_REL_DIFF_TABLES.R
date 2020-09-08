# Damien BEillouin
# Méta-synthèse sur la divsersification des cultures: code d'analyse
# travail initial: Octobre 2018, reprise Février 2019, mise à jour AVRIL 2020

### FOCUS SUR LES DONNEES DE RATIO #######

# Initialisation -----------------------------------------------------------------------------------------------

  # Load packages and scripts
    x<-c("magrittr", "tidyverse", "ggplot2","lme4","metafor","broom","ggpubr","reshape2","tidyr","dplyr")
    lapply(x, require, character.only = TRUE)
    #source('~/Documents/analyse_diversification/analyse_diversification/Load_data.R', echo=TRUE)
    
  # set wd
    setwd("~/Documents/analyse_diversification")

### Change some effect-sizes
    
    ES  <-ES %>% mutate(ES= ifelse(ID=="MA_A4" & Scenario_ID ==1,Reg_Intercept+ 3.62* Reg_Slope,ES))
    
    ES  <-ES %>% mutate(ES_lower_95= ifelse(ID=="MA_A4" & Scenario_ID ==1,
                                                    (Reg_Intercept+ 3.46* Reg_Slope)-(Reg_Slope_SE/2*1.96),
                                                    ES_lower_95))
    ES  <-ES %>% mutate(ES_upper_95= ifelse(ID=="MA_A4" & Scenario_ID ==1,
                                                    (Reg_Intercept+ 3.46* Reg_Slope)-(Reg_Slope_SE/2*1.96),
                                                    ES_upper_95))
    
    
    ES  <-ES %>% mutate(ES= ifelse(ID=="MA_A4" & Scenario_ID ==2,Reg_Intercept+ 3.62* Reg_Slope,ES))
    ES  <-ES %>% mutate(ES_lower_95= ifelse(ID=="MA_A4" & Scenario_ID ==2,(Reg_Intercept+ 3.62* Reg_Slope)-(Reg_Slope_SE/2*1.96),ES_lower_95))
    ES  <-ES %>% mutate(ES_upper_95= ifelse(ID=="MA_A4" & Scenario_ID ==2,(Reg_Intercept+ 3.62* Reg_Slope)-(Reg_Slope_SE/2*1.96),ES_upper_95))
    
    ES  <-ES %>% mutate(ES= ifelse(ID=="MA_A4" & Scenario_ID ==3,Reg_Intercept+ 3.30* Reg_Slope,ES))
    ES  <-ES %>% mutate(ES_lower_95= ifelse(ID=="MA_A4" & Scenario_ID ==3,(Reg_Intercept+ 3.30* Reg_Slope)-(Reg_Slope_SE/2*1.96),ES_lower_95))
    ES  <-ES %>% mutate(ES_upper_95= ifelse(ID=="MA_A4" & Scenario_ID ==3,(Reg_Intercept+ 3.30* Reg_Slope)-(Reg_Slope_SE/2*1.96),ES_upper_95))
    
    ES  <-ES %>% mutate(ES= ifelse(ID=="MA_A4" & Scenario_ID ==4,Reg_Intercept+ 3.14* Reg_Slope,ES))
    ES  <-ES %>% mutate(ES_lower_95= ifelse(ID=="MA_A4" & Scenario_ID ==4,(Reg_Intercept+ 3.14* Reg_Slope)-(Reg_Slope_SE/2*1.96),ES_lower_95))
    ES  <-ES %>% mutate(ES_upper_95= ifelse(ID=="MA_A4" & Scenario_ID ==4,(Reg_Intercept+ 3.14* Reg_Slope)-(Reg_Slope_SE/2*1.96),ES_upper_95))
    
    ES  <-ES %>% mutate(ES= ifelse(ID=="MA_A4" & Scenario_ID ==5,Reg_Intercept+ 3.43* Reg_Slope,ES))
    ES  <-ES %>% mutate(ES_lower_95= ifelse(ID=="MA_A4" & Scenario_ID ==5,(Reg_Intercept+ 3.43* Reg_Slope)-(Reg_Slope_SE/2*1.96),ES_lower_95))
    ES  <-ES %>% mutate(ES_upper_95= ifelse(ID=="MA_A4" & Scenario_ID ==5,(Reg_Intercept+ 3.43* Reg_Slope)-(Reg_Slope_SE/2*1.96),ES_upper_95))
    
    ES  <-ES %>% mutate(ES= ifelse(ID=="MA_A4" & Scenario_ID %in% c(6,7,8),Reg_Intercept+ 2.34* Reg_Slope,ES))
    ES  <-ES %>% mutate(ES_lower_95= ifelse(ID=="MA_A4" & Scenario_ID %in% c(6,7,8),(Reg_Intercept+ 2.34* Reg_Slope)-(Reg_Slope_SE/2*1.96),ES_lower_95))
    ES  <-ES %>% mutate(ES_upper_95= ifelse(ID=="MA_A4" & Scenario_ID %in% c(6,7,8),(Reg_Intercept+ 2.34* Reg_Slope)-(Reg_Slope_SE/2*1.96),ES_upper_95))
    
    ES  <-ES %>% mutate(ES= ifelse(ID=="MA_A4" & Scenario_ID %in% c(9,10),Reg_Intercept+ 2.34* Reg_Slope,ES))
    ES  <-ES %>% mutate(ES_lower_95= ifelse(ID=="MA_A4" & Scenario_ID %in% c(9,10),(Reg_Intercept+ 1.68* Reg_Slope)-(Reg_Slope_SE/2*1.96),ES_lower_95))
    ES  <-ES %>% mutate(ES_upper_95= ifelse(ID=="MA_A4" & Scenario_ID %in% c(9,10),(Reg_Intercept+ 1.68* Reg_Slope)-(Reg_Slope_SE/2*1.96),ES_upper_95))
    
    # on convertit Himmelsetein de Kg ha-1 à en T.ha-1
    ES  <-ES %>% mutate(ES= ifelse(ID=="MA_H2",ES/1000,ES))
    ES  <-ES %>% mutate(ES_lower_95= ifelse(ID=="MA_H2",ES_lower_95/1000,ES_lower_95))
    ES  <-ES %>% mutate(ES_upper_95= ifelse(ID=="MA_H2",ES_upper_95/1000,ES_upper_95))
    
    
    
    
    
    
# Select main  effect-sizes ------------------------------------------------------------------------
    
    ES_TMP<-ES                                 %>%  
      filter(Keep_global_analyse=="YES")       %>%  # We analyse only "main" effect sizes (not sub-scenarios)
      select(ID,contains("Scenario_"),DIVERS,code_Type_Output,
             Sub_Cat_Output,Sub_sub_Cat_Output,Metric_Output,
             Nb_Paired_Data,ES_lower_95, ES, ES_upper_95,
             contains("Reg"),contains("Treat_"))%>% # select some columns
      mutate(UN= 1)                             %>%
      group_by(ID)                              %>%
      mutate(tid= cumsum(UN))                   %>%
      mutate(vi= FUN_var_from_CI(mean=ES, CI=ES_upper_95))
    
    ## Check Data
       ggplot(aes(ES,vi),data=ES_TMP)+geom_point()+ theme_bw()+ facet_wrap(~Metric_Output,scales="free")
      # Etrange pour les ES_TMPerences, l'augmentation de variance!!!
    
# Split ratio, relative ES_TMPerences and ES_TMPerences ----------------------------------------------

    
    ####ES_TMPERENCE
    # ES_TMP<-ES_TMP %>% filter(Metric_Output=="ES_TMPerence")
               
    
# Final table for analyses ----------------------------------------------------------------------------------------

        ES_TMP<-ES_TMP %>%select(ID, ES, vi, Sub_Cat_Output,DIVERS, Scenario_ID,Metric_Output, Nb_Paired_Data,Sub_sub_Cat_Output) %>% ungroup() %>%
      mutate(NUM= as.character(1: length(ID)))
        ES_TMP$ID<-as.factor(as.character(ES_TMP$ID))
        
        SCORE<- QUAL %>% dplyr::select(ID, Score)
        setdiff(SCORE$ID,ES_TMP$ID)
        setdiff(ES_TMP$ID,SCORE$ID)
        ES_TMP<-left_join(ES_TMP, SCORE)
        
        
    #RATIO2<-RATIO2[!RATIO2$vi==0,]
    #RATIO2<-RATIO2[!is.na(RATIO2$ES),]
    
