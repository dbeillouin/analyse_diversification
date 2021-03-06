# Damien BEillouin
# Méta-synthèse sur la divsersification des cultures: code d'analyse
# travail initial: Octobre 2018, reprise Février 2019

### ANALYSE DES DONNEES TOUT TYPE DE DIVERSIFICATION ENSEMBLE ######


# Initialisation --------------------------------------------------------------------------------------------------
  # setwd
  setwd("~/Documents/analyse_diversification")

  # Load packages
    x<-c("magrittr", "tidyverse", "ggplot2","lme4","metafor","broom","ggpubr")
    lapply(x, require, character.only = TRUE)
  
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
      lm(ES~1+Score,weights = 1/((vi)), data=DAT)}
    
    #2. Function for mixed model (random= ID)
    model_lmer_ID<-function(DAT) {
      lmer(ES~1+Score+(1|ID),weights = 1/(vi), data=DAT)}
    
    #3. Function for mixed model (random= NUM)
    model_lmer_NUM<-function(DAT) {
      lmer(ES ~ 1 + Score+(1 | NUM), weights = 1/vi, data=DAT, 
                       control=lmerControl(check.nobs.vs.nlev="ignore", 
                                           check.nobs.vs.nRE="ignore"))}
    
    #4. Function for mixed model (random= NUM and ID)
    model_lmer_ID_NUM<-function(DAT) {
      lmer(ES ~ 1 + Score+(1 |ID/ NUM), weights = 1/vi, data=DAT, 
           control=lmerControl(check.nobs.vs.nlev="ignore",
                               check.nobs.vs.nRE       ="ignore"))}
    
    #5. Function for rma (fixed)
    model_rma_lm<-function(DAT) {
      rma(ES, vi, mods= Score, data=DAT, method="FE")}
    
    #6. Function for rma (random= ID)
    model_rma_ID<-function(DAT) {
      rma(ES, vi, mods= Score,random=~1|ID, data=DAT, method="REML")}
    
    #7. Function for rma (random= NUM)
    model_rma_NUM<-function(DAT) {
      rma(ES, vi,mods= Score, random=~1|NUM, data=DAT, method="REML")}
    
    #8. Function for rma (random= NUM and ID)
    model_rma_ID_NUM<-function(DAT) {
      rma.mv(ES, vi,mods= Score, data=DAT,random=~1|ID/NUM, method="REML")  }
    
  ######## Run the linear models (lm and rma) ########
    
    RATIO<- ES_TMP %>% filter(Metric_Output=="Ratio")
    RATIO<- RATIO %>% filter(!is.na(vi))
    
    COHEN<- ES_TMP %>% filter(Metric_Output=="Cohen's d")
    COHEN<- COHEN %>% filter(!is.na(vi))
    lm(ES~1,weights = 1/((vi)), data=COHEN)

    lm<-RATIO %>%
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
    
    
    res<- rma(ES, vi, mods= ~DIVERS-1,data=RATIO, method="FE")
    res<-    rma(ES, vi,slab= ID, data=RATIO, method="REML")
    funnel(res)
    
    
    ESSAI <-RATIO %>% filter(Sub_Cat_Output=="Water Quality")
    
    ### Sans considérer les poids ######
    rma(ES, vi, mods= ~DIVERS-1,data=ESSAI, method="FE")

    
  #### Analyse the results ##########
    
    res_lm  <- lm %>% 
             unnest(tidied_lm, .drop = TRUE) %>%
             mutate(method= "lm")            %>%
             mutate(conf.low    = estimate- 1.96*std.error,
                      conf.high = estimate +1.96*std.error) %>% 
           filter(term== '(Intercept)')
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
    A <- RATIO %>% group_by(Sub_Cat_Output, DIVERS)       %>% 
                   summarise(length = n(),
                             Nb_ID=n_distinct(ID))         %>%
      mutate(Diff = length-Nb_ID)   %>% 
                   filter(Diff>1)%>% 
                   filter(Nb_ID>2)  %>% 
      mutate(CODE= paste(Sub_Cat_Output, DIVERS))
    
    
    lmer(ES ~ 1 + Score+(1|NUM), weights = 1/vi, data=lmer[-1,], 
         control=lmerControl(check.nobs.vs.nlev="ignore", 
                             check.nobs.vs.nRE="ignore"))
    
    
    lmer <- RATIO %>%
      mutate(CODE= paste(Sub_Cat_Output, DIVERS))  %>% 
      filter(CODE %in% A$CODE)                     %>%
       ungroup()                                    %>%
      group_by(Sub_Cat_Output, DIVERS)             %>% 
      nest()                                       %>% 
      mutate(#NUM = 1:nrow(.),
             model_lmer_ID        = map(data, model_lmer_ID),           # First model
             tidied_lmer_ID        = map(model_lmer_ID, tidy, effects = "fixed", conf.int=TRUE),
             glanced_lmer_ID       = map(model_lmer_ID, glance),
             augmented_lmer_ID     = map(model_lmer_ID, augment),
             model_lmer_NUM        = map(data, model_lmer_NUM),        # second model
             tidied_lmer_NUM       = map(model_lmer_NUM, tidy, effects = "fixed", conf.int=TRUE),
             glanced_lmer_NUM      = map(model_lmer_NUM, glance),
             augmented_lmer_NUM    = map(model_lmer_NUM, augment),
             model_lmer_ID_NUM     = map(data, model_lmer_ID_NUM),# third model
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
      mutate(method="lmer_ID")                 %>% 
      filter(term == "(Intercept)")
    
    res_lmer_NUM    <- lmer                    %>%
      unnest(tidied_lmer_NUM, .drop = TRUE)    %>%
      mutate(method="lmer_NUM")                %>% 
      filter(term == "(Intercept)")
    
    res_lmer_ID_NUM <- lmer                    %>% 
      unnest(tidied_lmer_ID_NUM, .drop = TRUE) %>% 
      mutate(method="lmer_ID_NUM")             %>% 
      filter(term == "(Intercept)")
    
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
  
    
    ggplot(res, aes(x=interaction(method,Sub_Cat_Output) , y=estimate, group=DIVERS, colour=DIVERS)) +
      geom_point(aes(), na.rm=TRUE, position="dodge") +
      #geom_errorbar(aes(x=method, ymax=conf.high, ymin=conf.low), na.rm=TRUE, position="dodge")+ 
                        theme_bw()+
                      geom_hline(yintercept=0, linetype=2)+
                       theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                       ylim(c(-5,5))
    
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
    

    #######
    ### GRAPHIQUES
    
    # type forest plot classique
    DATA<- res %>% filter(method=="lmer_ID_NUM")
    DATA$CODE<- paste(DATA$Sub_Cat_Output, DATA$DIVERS)
    NB_ES<-RATIO %>% group_by(Sub_Cat_Output, DIVERS)  %>%  count()
    NB_MA<-RATIO %>% group_by(Sub_Cat_Output, DIVERS)  %>% 
      summarise(NB_MA = n_distinct(ID)) 
    NB_Data<-RATIO %>% group_by(Sub_Cat_Output, DIVERS)  %>% 
      summarise(NB_paired = sum(Nb_Paired_Data, na.rm=TRUE)) 
    
    COUNT<- merge(merge(NB_ES,NB_MA), NB_Data)
    COUNT$Label<- paste0(COUNT$Sub_Cat_Output, "%", "(",COUNT$NB_MA,"-", COUNT$n,"-", COUNT$NB_paired, ")")
    COUNT$CODE<- paste(COUNT$Sub_Cat_Output, COUNT$DIVERS)
    COUNT %<>% filter(CODE %in% DATA$CODE) 
    
    DATA<-merge(DATA, COUNT)
    addline_format <- function(x,...){
      gsub('%','\n',x)
    }
  
    ggplot(data=DATA, aes(color=DIVERS, group=DIVERS))+ 
      geom_point(aes(x=reorder(Sub_Cat_Output,estimate), y= estimate),position=position_dodge(width=0.6))+
      geom_errorbar(aes(x=Sub_Cat_Output, ymin=conf.low, ymax=conf.high), width=.1,position=position_dodge(width=0.6)) +
      geom_hline(yintercept=0, linetype=2)+
      scale_x_discrete(breaks=(DATA$Sub_Cat_Output), 
                       labels=addline_format(DATA$Label))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            panel.grid.major = element_blank())+
      labs(title="ALL diversification strategies together", x= "Studied Output", y= "Effect size")+
      labs(x="",y= "Ratio diversified systems/ less diversified")+ 
      scale_colour_viridis_d(option = "plasma")

                             group=Sub_Cat_Output), 
              position=position_dodge(width=1),
              size=3,angle=90)+
    
    # type flower Power 
    library(tidyverse)

    ESSAI<- RATIO %>% select(ES,Sub_Cat_Output)
    ESSAI2<-DATA %>%  select(estimate,Sub_Cat_Output,conf.low,conf.high)
    names(ESSAI2)[1]<-"ES"
    ESSAI3<-RATIO %>% group_by(Sub_Cat_Output)  %>%  count()
    
    ggplot(ESSAI, aes(x = Sub_Cat_Output, y = ES)) +
      geom_hline(yintercept=c(-1,0,1),col=c("gray","black","gray")) + 
      geom_text(data=ESSAI3,aes(Sub_Cat_Output,1.3,label=n, group=Sub_Cat_Output), position=position_dodge(width=1), size=4)+
      geom_violin(aes(fill = Sub_Cat_Output), alpha=0.2) +
      geom_point(data = ESSAI2, color = "red")+
      geom_errorbar(data = ESSAI2,aes(ymin=conf.low, ymax=conf.high), colour="black", width=.1)+
      #coord_polar(theta = "x")+
      theme_pubr()
    
    