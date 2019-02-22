# Damien BEillouin
# Méta-synthèse sur la divsersification des cultures: code d'analyse
# travail initial: Octobre 2018, reprise Février 2019

### ANALYSE DES DONNEES TOUT TYE DE DIVERSIFICATION ENSEMBLE ######


# Initialisation --------------------------------------------------------------------------------------------------
# setwd
setwd("~/Documents/shinyapp/")
# Load packages
x<-c("magrittr", "tidyverse", "ggplot2","lme4","metafor","broom")
lapply(x, require, character.only = TRUE)

# load Data
source('~/Documents/analyse_diversification/1.Description_donnee.R', echo=TRUE)
  
RATIO2<-merge(RATIO2, QUAL)

#RATIO2 %>%filter(vi>0.2)


# 1. Main effect (all diversification strategies together) -----------------------------------------------------------
    
 ##  1.1. Linear models
    
    ## Define some functions (complementary to dplyr)
    tidy.rma <- function(x, ...) {
      Y<-summary(x)
      tibble(
        estimate = Y$b[1],
        std.error = Y$se,
        statistic = Y$zval,
        p.value = Y$pval,
        conf.low = Y$ci.lb,
        conf.high = Y$ci.ub
      )
    }
    glance.rma<-function(x, ...) {
      with(
        x,
        tibble(
          AIC= AIC(x),
          BIC= BIC(x),
          tau2 = tau2,
          se.tau2 = se.tau2,
          k = k,
          p = p,
          m = m,
          QE = QE,
          QEp = QEp,
          QM = QM,
          QMp = QMp,
          I2 = I2,
          H2 = H2,
          # R2 = R2,
          int.only = int.only
        )
      )
    }
    
  ## Define the types of meta-anaytical models
    #1. Function for linear model
    model_lm<-function(DAT) {
      lm(ES~1,weights = 1/((vi)), data=DAT)}
    
    #2. Function for mixed model (random= ID)
    model_lmer_ID<-function(DAT) {
      lmer(ES~1+(1|ID),weights = 1/(vi), data=DAT)}
    
    #3. Function for mixed model (random= NUM)
    model_lmer_NUM<-function(DAT) {
      lmer(ES ~ 1 + (1 | NUM), weights = 1/vi, data=DAT, 
                       control=lmerControl(check.nobs.vs.nlev="ignore", 
                                           check.nobs.vs.nRE="ignore"))}
    
    #4. Function for mixed model (random= NUM and ID)
    model_lmer_ID_NUM<-function(DAT) {
      lmer(ES ~ 1 + (1 |ID/ NUM), weights = 1/vi, data=DAT, 
           control=lmerControl(check.nobs.vs.nlev="ignore",
                               check.nobs.vs.nRE       ="ignore"))}
    
    #5. Function for rma (fixed)
    model_rma_lm<-function(DAT) {
      rma(ES, vi, data=DAT, method="FE")}
    
    #6. Function for rma (random= ID)
    model_rma_ID<-function(DAT) {
      rma(ES, vi,slab= ID, data=DAT, method="REML")}
    
    #7. Function for rma (random= NUM)
    model_rma_NUM<-function(DAT) {
      rma(ES, vi,slab= NUM, data=DAT, method="REML")}
    
    #8. Function for rma (random= NUM and ID)
    model_rma_ID_NUM<-function(DAT) {
      rma.mv(ES, vi, data=DAT,random=~1|ID/NUM, method="REML")  }
    
  ######## Run the linear models (lm and rma) ########
    lm<-RATIO2 %>%
      ungroup() %>%
      nest(-Sub_Cat_Output)%>%                  
      mutate(model_lm       = map(data, model_lm),       # with lme4
             tidied_lm      = map(model_lm, tidy),       # tidy data
             glanced_lm     = map(model_lm, glance),     # analyse models
             augmented_lm   = map(model_lm, augment),    # predict
             model_rma_lm   = map(data, model_rma_lm),   # with metafor
             tidied_rma_lm      = map(model_rma_lm, tidy.rma),
             glanced_rma_lm     = map(model_rma_lm, glance.rma))
            # augmented_rma_lm   = map(model_rma_lm, augment))
    
    
    ESSAI <-RATIO2 %>% filter(Sub_Cat_Output=="Soil quality")
    
    ### Sans considérer les poids ######
    rma(ES, vi, mods= ~DIVERS-1,data=ESSAI, method="FE")

    
  #### Analyse teh results ##########
    
    res_lm  <- lm %>% 
             unnest(tidied_lm, .drop = TRUE) %>%
             mutate(method= "lm")            %>%
             mutate(conf.low    = estimate- 1.96*std.error,
                      conf.high = estimate +1.96*std.error)
    res_rma  <- lm %>%
             unnest(tidied_rma_lm, .drop = TRUE)%>%
             mutate(method= "rma")
    
    res     <- bind_rows(res_lm,res_rma) 
    
    ggplot(data=res)+ 
      geom_point(aes(x=method, y= estimate))+
      geom_errorbar(aes(x=method, ymin=conf.low, ymax=conf.high), width=.1) +
      facet_wrap(~Sub_Cat_Output,scales="free")+ theme_bw()+
      geom_hline(yintercept=0, linetype=2)
    ggsave("Estimates_linear_ALL_strategies.pdf")
    
    ### Vérification des résultats en calculan à la main :
    ESSAI<- RATIO2 %>% filter(Sub_Cat_Output=="Soil quality")
    
    ESSAI$w<-1/ESSAI$vi
    sum(ESSAI$ES*ESSAI$w)/ sum(ESSAI$w)
    sqrt(1/ sum(ESSAI$w))
    
    summary(lm$model_lm[[1]])
    summary(lm(ES~1,weights = 1/((vi)), data=ESSAI))
    summary(lm$model_rma_lm[[1]])
    
    # Adjustements
    Adj_lm    <- lm %>%
      unnest(glanced_lm, .drop = TRUE)      %>% 
      mutate(method= "lm")
    
    Adj_rma   <- lm %>%  
      unnest(glanced_rma_lm, .drop = TRUE)  %>%
      mutate(method= "rma")
    
    Adj       <- bind_rows(Adj_lm,Adj_rma)  %>%
       group_by(Sub_Cat_Output)             %>% 
      mutate(MAX= ifelse(AIC==max(AIC),0,1))
    
    ggplot(Adj, aes(method,Sub_Cat_Output)) + 
      geom_tile(aes(fill =factor(MAX) ), colour = "white") +
      geom_text(aes(label = round(AIC, 0)))+
      scale_fill_manual(values=c("white","darkgreen"))+
      theme_bw()+   theme(legend.position="none")
    
    ggsave("AIC_linear_ALL_strategies.pdf")
    
    # Predictions
    Pred_lm  <- lm%>%  
      unnest(augmented_lm, .drop = TRUE)  %>% 
      mutate(method= "lm")
    

  # 1.2  Mixed models (lmer and rma)
    
    # Select Sub_Cat_Output with more than 2 ID
    A <- RATIO2 %>% group_by(Sub_Cat_Output)       %>% 
                   summarise(length = n(),
                             Nb_ID=n_distinct(ID)) %>%
                   filter(Nb_ID>1)
    
    lmer <- RATIO2 %>%
      filter(Sub_Cat_Output %in% A$Sub_Cat_Output) %>%
      ungroup()                                    %>%
      nest(-Sub_Cat_Output)                        %>% 
      mutate(#NUM = 1:nrow(.),
             model_lmer_ID        = map(data, model_lmer_ID),           # First model
             tidied_lmer_ID        = map(model_lmer_ID, tidy, effects = "fixed", conf.int=TRUE),
             glanced_lmer_ID       = map(model_lmer_ID, glance),
             augmented_lmer_ID     = map(model_lmer_ID, augment),
             model_lmer_NUM        = map(data, model_lmer_NUM),        # second model
             tidied_lmer_NUM       = map(model_lmer_NUM, tidy, effects = "fixed", conf.int=TRUE),
             glanced_lmer_NUM      = map(model_lmer_NUM, glance),
             augmented_lmer_NUM    = map(model_lmer_NUM, augment),
             model_lmer_ID_NUM     = map(data, model_lmer_ID_NUM),    # third model
             tidied_lmer_ID_NUM    = map(model_lmer_ID_NUM, tidy, effects = "fixed", conf.int=TRUE),
             glanced_lmer_ID_NUM   = map(model_lmer_ID_NUM, glance),
             augmented_lmer_ID_NUM = map(model_lmer_ID_NUM, augment),
             model_rma_ID          = map(data, model_rma_ID),         # Fouth model
             tidied_rma_ID         = map(model_rma_ID, tidy.rma, effects = "fixed", conf.int=TRUE),
             glanced_rma_ID        = map(model_rma_ID, glance.rma),
             model_rma_NUM         = map(data, model_rma_NUM),        # fifth model
             tidied_rma_NUM        = map(model_rma_NUM, tidy.rma, effects = "fixed", conf.int=TRUE),
             glanced_rma_NUM       = map(model_rma_NUM, glance.rma),
             model_rma_ID_NUM      = map(data, model_rma_ID_NUM),     # Sixth model
             tidied_rma_ID_NUM     = map(model_rma_ID_NUM, tidy.rma, effects = "fixed", conf.int=TRUE))
             #glanced_rma_ID_NUM   = map(model_rma_ID_NUM, glance.rma))
    
    ## Analyse the results
    
    # Estimates
    res_lmer_ID     <- lmer                    %>%
      unnest(tidied_lmer_ID, .drop = TRUE)     %>%
      mutate(method="lmer_ID")
    
    res_lmer_NUM    <- lmer                    %>%
      unnest(tidied_lmer_NUM, .drop = TRUE)    %>%
      mutate(method="lmer_NUM")
    
    res_lmer_ID_NUM <- lmer                    %>% 
      unnest(tidied_lmer_ID_NUM, .drop = TRUE) %>% 
      mutate(method="lmer_ID_NUM")
    
    res_rma_ID      <- lmer                    %>%  
      unnest(tidied_rma_ID, .drop = TRUE)      %>%
      mutate(method="rma_ID")
    
    res_rma_NUM     <- lmer                    %>%  
      unnest(tidied_rma_NUM, .drop = TRUE)     %>% 
      mutate(method="rma_NUM")
    
    res_rma_ID_NUM  <- lmer                    %>%  
      unnest(tidied_rma_ID_NUM, .drop = TRUE)  %>%
      mutate(method="rma_ID_NUM")
    
    res     <-bind_rows(res_lmer_ID,res_lmer_NUM,
                        res_lmer_ID_NUM,res_rma_ID,
                        res_rma_NUM,res_rma_ID_NUM) 
    
    ggplot(data=res)+ 
                      geom_point(aes(x=method, y= estimate))+
                      geom_errorbar(aes(x=method, ymin=conf.low, ymax=conf.high), width=.1) +
                      facet_wrap(~Sub_Cat_Output)+ theme_bw()+
                      geom_hline(yintercept=0, linetype=2)+
                       theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    ggsave("Estimates_mixed_ALL_strategies.pdf")
    
    # Adjustements
    Adj_lmer_ID      <-lmer                   %>%
      unnest(glanced_lmer_ID, .drop = TRUE)   %>% 
      mutate(method="lmer_ID")
    
    Adj_lmer_NUM     <-lmer                   %>%
      unnest(glanced_lmer_NUM, .drop = TRUE)  %>%
      mutate(method="lmer_NUM")
    
    Adj_lmer_ID_NUM  <-lmer                      %>%
      unnest(glanced_lmer_ID_NUM, .drop = TRUE)  %>%
      mutate(method="lmer_ID_NUM")
    
    Adj_rma_ID       <-lmer                      %>% 
      unnest(glanced_rma_ID, .drop = TRUE)       %>% 
      mutate(method="rma_ID")
    
    Adj_rma_NUM      <-lmer                      %>%
      unnest(glanced_rma_NUM, .drop = TRUE)      %>% 
      mutate(method="rma_NUM")
    
    #Adj_rma_ID_NUM  <-lmer%>%   unnest(glanced_rma_ID_NUM, .drop = TRUE) %>% mutate(method="rma_ID_NUM")
    Adj     <-bind_rows(Adj_lmer_ID,Adj_lmer_NUM,
                        Adj_lmer_ID_NUM,Adj_rma_ID,Adj_rma_NUM) %>%
                group_by(Sub_Cat_Output)                        %>%
                mutate(MAX= ifelse(AIC==min(AIC),1,0))
    
    ggplot(Adj, aes(method,Sub_Cat_Output)) + 
      geom_tile(aes(fill =factor(MAX) ), colour = "white") +
      geom_text(aes(label = round(AIC, 0)))+
      scale_fill_manual(values=c("white","darkgreen"))+
      theme_bw()+theme(legend.position="none")
    
    ggsave("AIC_mixed_ALL_strategies.pdf")
    
    # Predictions 
    PRED_lmer_ID_NUM <- lmer %>%   unnest(augmented_lmer_ID_NUM, .drop = TRUE)
    PRED_lmer_NUM    <- lmer %>%   unnest(augmented_lmer_NUM, .drop = TRUE)
    PRED_lmer_ID     <- lmer %>%   unnest(augmented_lmer_ID, .drop = TRUE)
    #lmer%>%   unnest(glanced_rma_NUM, .drop = TRUE)
    

