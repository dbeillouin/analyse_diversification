# Damien BEillouin
# Méta-synthèse sur la divsersification des cultures: code d'analyse
# travail initial: Octobre 2018, reprise Février 2019

### ANALYSE DES DONNEES TOUT TYE DE DIVERSIFICATION ENSEMBLE ######
## article Evidence Map #####

# Initialisation ------------------------------------------------------------------------------------------------
Effect <-  read.csv("~/Documents/shinyapp/FINAL/Effect_Size.csv", sep=";")

Quality <- read.csv("~/Documents/shinyapp/FINAL/Quality.csv", sep=";", comment.char="#")
Quality %<>%  select(ID, Score)

Effect<- merge(Effect,Quality, by="ID")


  Primary_Studies <- read.csv("~/Documents/shinyapp/FINAL/Primary_Studies.csv", sep=";")
AA <- left_join(Primary_Studies, Effect, by="ID")

AA$DIVERS<- as.character(AA$DIVERS)
s <- strsplit(AA$DIVERS, split = ",| and")
AA<- data.frame(Sub_Cat_Output = rep(AA$Sub_Cat_Output, sapply(s, length)), 
                ID = rep(AA$ID, sapply(s, length)),
                DOI = rep(AA$DOI, sapply(s, length)),
                DIVERS = unlist(s))

AA$DIVERS <- recode_factor(AA$DIVERS, `Buffer` = "Landscape",
                           " crop rotation" = "Rotation",
                           " Rotation"  ="Rotation",
                           "rotation" = "Rotation",
                           "intercropping" = "Intercropping",
                           " intercropping" = "Intercropping",
                            " Intercropping" = "Intercropping",
                           " Associated_plant" = "Associated_plant")

#AA <- AA %>%  filter(DIVERS =="Landscape") 
#AA <- AA %>%  filter( Sub_Cat_Output=="Biodiv. Assoc.") 
COUNT <- AA %>%  group_by(ID, DIVERS,Sub_Cat_Output) %>%  summarise(ALL_studies = n_distinct(DOI)) %>%  group_by(DIVERS,Sub_Cat_Output) %>%  summarise(ALL_studies = sum(ALL_studies)) 
COUNT_unique <- AA %>%  group_by(DIVERS,Sub_Cat_Output) %>%  summarise(Unique_Elements = n_distinct(DOI))


# Par type de diversification

# Nombre d'effect size


Effect$DIVERS<- as.character(Effect$DIVERS)
s <- strsplit(Effect$DIVERS, split = ",| and")
Effect<- data.frame(Sub_Cat_Output = rep(Effect$Sub_Cat_Output, sapply(s, length)), 
                    ID = rep(Effect$ID, sapply(s, length)), 
                    Score = rep(Effect$Score, sapply(s, length)), 
                DIVERS = unlist(s))

Effect$DIVERS <- recode_factor(Effect$DIVERS, `Buffer` = "Landscape",
                           " crop rotation" = "Rotation",
                           " Rotation"  ="Rotation",
                           "rotation" = "Rotation",
                           "intercropping" = "Intercropping",
                           " intercropping" = "Intercropping",
                           " Intercropping" = "Intercropping",
                           " Associated_plant" = "Associated_plant")


ES<- Effect %>% 
  group_by(Sub_Cat_Output, DIVERS) %>% count()

# Nombre de meta-analyses
MA <- Effect %>% 
  group_by(Sub_Cat_Output,DIVERS) %>% summarise(Unique_Elements = n_distinct(ID))%>%
  arrange(desc(Unique_Elements))

names(ES)[3]<- "ES"
names(MA)[3]<- "MA"
AA<-merge(MA, ES, all=TRUE)

# qualité moyenne des meta-analyses
Effect$Score<-as.numeric(as.character(Effect$Score))
QUAL <-Effect %>% 
  group_by(Sub_Cat_Output,DIVERS) %>% summarise(QUAL = mean(Score,na.rm=TRUE))%>%
  arrange(desc(QUAL))


names(COUNT)[3]<-"IND"
COUNT$ID<-NULL

ALL <- merge(merge(merge(AA,QUAL, all=TRUE), COUNT, all=TRUE),COUNT_unique,all=TRUE)
ALL<- ALL %>% filter(!is.na(Sub_Cat_Output))

ALL$cutMA<-cut(ALL$MA, breaks=c(0,2,5,16))
## plot
colors <- colorRampPalette(c("ivory3", "green4"))(30)
colors2<- c("ivory3", "olivedrab3","green4")

new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

ESSAI <- ALL %>% filter(Sub_Cat_Output=="Yield")
P<-ggplot(ALL, aes(Sub_Cat_Output, DIVERS )) +
  geom_tile(aes(fill = MA), color = "black") +
  scale_fill_gradient(low = "ivory3", high = "green4") +
  #scale_fill_manual(values=colors2)+
  geom_point(aes(size=(IND),color=QUAL/5), stroke=1)+
  geom_point(aes(size=(Unique_Elements)),color="black", stroke=1, shape=21)+
  scale_radius(range = c(5, 15))+
 # scale_size(range=c(5,15))+
  scale_color_gradient(low = "red", high = "blue") +
 # scale_fill_gradient(low = "white", high = "steelblue") +
  xlab("Stategies of crop diversificaion ") +
  ylab("outcomes") + theme_pubr(base_size=24)+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Number of Meta-analyses")

P


## Ratio Nb études ind/ nb MA

ALL$Ratio<- ALL$Unique_Elements/ALL$MA
ALL$RATIO_UNI<-ALL$Unique_Elements/ALL$IND*100
