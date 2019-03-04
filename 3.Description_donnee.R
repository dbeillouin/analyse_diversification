# Damien BEillouin
# Méta-synthèse sur la divsersification des cultures: code d'analyse
# travail initial: Octobre 2018, reprise Février 2019

### SCRIPT POUR CHARGER LES DONNEES ET FAIRE LE POINT#######

# Initialisation --------------------------------------------------------------------------------------------------
# setwd
setwd("~/Documents/shinyapp/")
# Load packages
x<-c("magrittr", "tidyverse", "ggplot2","lme4","metafor","broom")
lapply(x, require, character.only = TRUE)

# load Data
RATIO2 <- read.csv2("~/Documents/shinyapp/ES_analyse_Paper.csv")
QUAL   <- read.csv("Quality_09_11.csv", sep=";", encoding="UTF-8-MAC", comment.char="#")%>%
  select(ID, Score)

RATIO2<-merge(RATIO2, QUAL)

#RATIO2 %>%filter(vi>0.2)


# Point sur les données -------------------------------------------------------------------------------------------

#1. Globalement
    # nombre d'effect sizes
    RATIO2 %>%
      ungroup() %>%
      group_by(Sub_Cat_Output) %>% count() %>% arrange(desc(n))
    
    # Nombre de meta-analyses
    RATIO2 %>%
      ungroup() %>%
      group_by(Sub_Cat_Output) %>% summarise(Unique_Elements = n_distinct(ID))%>%
      arrange(desc(Unique_Elements))
    
    # qualité moyenne des meta-analyses
    RATIO2$Score<-as.numeric(as.character(RATIO2$Score))
    RATIO2 %>%
      ungroup() %>%
      group_by(Sub_Cat_Output) %>% summarise(QUAL = mean(Score,na.rm=TRUE))%>%
      arrange(desc(QUAL))
    
# Par type de diversification
    
    # Nombre d'effect size
    RATIO2 %>%
      ungroup() %>%
      group_by(Sub_Cat_Output, DIVERS) %>% count()
    
    # Nombre de meta-analyses
    RATIO2 %>%
      ungroup() %>%
      group_by(Sub_Cat_Output,DIVERS) %>% summarise(Unique_Elements = n_distinct(ID))%>%
      arrange(desc(Unique_Elements))
    
    # qualité moyenne des meta-analyses
    RATIO2$Score<-as.numeric(as.character(RATIO2$Score))
    RATIO2 %>%
      ungroup() %>%
      group_by(Sub_Cat_Output,DIVERS) %>% summarise(QUAL = mean(Score,na.rm=TRUE))%>%
      arrange(desc(QUAL))
