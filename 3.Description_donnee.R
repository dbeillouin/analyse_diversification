# Damien BEillouin
# Méta-synthèse sur la divsersification des cultures: code d'analyse
# travail initial: Octobre 2018, reprise Février 2019

### SCRIPT POUR CHARGER LES DONNEES ET FAIRE LE POINT#######

# Initialisation --------------------------------------------------------------------------------------------------
# setwd
setwd("~/Documents/analyse_diversification")
# Load packages
x<-c("magrittr", "tidyverse", "ggplot2","lme4","metafor","broom")
lapply(x, require, character.only = TRUE)


# load Data

SCORE<- QUAL %>% dplyr::select(ID, Score)
ES_TMP<-merge(ES_TMP, SCORE)
Year<- Des_MA %>% dplyr::select(ID,Publication_Year)
ES_TMP<-merge(ES_TMP, Year)

#RATIO2 %>%filter(vi>0.2)


# Point sur les données -------------------------------------------------------------------------------------------

#1. Globalement
    # nombre d'effect sizes
ES_TMP %>%
      ungroup() %>%
      group_by(Sub_Cat_Output) %>% count() %>% arrange(desc(n))
    
    # Nombre de meta-analyses
ES_TMP %>%
      ungroup() %>%
  filter(!DIVERS=="2 AND MORE") %>% 
      group_by(DIVERS) %>% summarise(Unique_Elements = n_distinct(ID))%>%
      arrange(desc(Unique_Elements)) %>% 
  ggdotchart(x = "DIVERS", y = "Unique_Elements",
             color = "DIVERS",                                # Color by groups
             #palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
             sorting = "descending",                       # Sort value in descending order
             add = "segments",                             # Add segments from y = 0 to dots
             rotate = TRUE,                                # Rotate vertically
             group = "DIVERS",                                # Order by groups
             dot.size = 9,                                 # Large dot size
             #label = mpg),                        # Add mpg values as dot labels
             font.label = list(color = "white", size = 9, 
                               vjust = 0.5),               # Adjust label parameters
             ggtheme = theme_pubr(base_size=24)                        # ggplot2 theme
  ) + theme(legend.position='none')+ labs(y="Number of meta-analyses", x="")+
  scale_colour_manual(values = c("#fe6766", "#993200", "black", "#68cb9b", "#ff9934", "#339866", "#3299fe", "#3266cb",  "red"))


ggsave("NB_MA_par_DIVERS.pdf", width = 20, height = 15, units = "cm")


# Nombre de meta-analyses par année
ES_TMP %>%
  filter(!DIVERS=="2 AND MORE") %>% 
  ungroup() %>%
  group_by(DIVERS, Publication_Year) %>% summarise(Unique_Elements = n_distinct(ID))%>%
  arrange(DIVERS,Publication_Year) %>% 
  ungroup() %>% 
  complete(DIVERS, Publication_Year = min(Publication_Year):2020)%>% 
  replace(., is.na(.), 0) %>% 
  group_by(DIVERS) %>% 
  mutate(cumulative = cumsum(Unique_Elements)) %>% 
  ggplot()+geom_line(aes(x=Publication_Year,
                         y= cumulative, 
                         group=DIVERS,
                         color= DIVERS),
                     size=2)+
  theme(legend.position='none')+ labs(y="Number of meta-analyses", x="")+
  theme_pubr(base_size=24)+ theme(legend.position = c(0.3, 0.8))+
  scale_colour_manual(values = c("#fe6766", "#993200", "black", "#68cb9b", "#ff9934", "#339866", "#3299fe", "#3266cb",  "red"))+
  labs(colour="")+
  scale_x_continuous(name="Years", breaks=seq(1990,2020, 5))

ggsave("NB_MA_par_DIVERS_per_year.pdf", width = 20, height = 20, units = "cm")

# Nombre de meta-analyses par année
TAB1<-ES_TMP %>%
  filter(!DIVERS=="2 AND MORE") %>% 
  ungroup() %>%
  group_by(Publication_Year) %>% summarise(Unique_Elements = n_distinct(ID))%>%
  arrange(Publication_Year) %>% 
  ungroup() %>% 
  complete(Publication_Year = min(Publication_Year):2020)%>% 
  replace(., is.na(.), 0) %>% 
  mutate(cumulative.MA = cumsum(Unique_Elements)) 

TAB2<- IND_STU %>% group_by(Year) %>% summarise(n=n_distinct(DOI))%>% 
  mutate(cumulative.PS = cumsum(n)) 
names(TAB2)[1]<- "Publication_Year"
  
TAB3<-merge(TAB1,TAB2, all=TRUE, by= 'Publication_Year')
TAB3 %>% 
  ggplot()+geom_area(aes(x=Publication_Year,
                         y= n),
                     size=1,
                    stat="identity",
                    fill="#FF9933",
                    color="black")+
  geom_line(aes(x=Publication_Year,
                y= cumulative.MA),
            size=2,
            stat="identity",
            color="darkgreen")+
  theme(legend.position='none')+ labs(y="Number of exp./meta-analyses", x="")+
  theme_pubr(base_size=18)+ theme(legend.position = c(0.3, 0.8))+
 # scale_colour_manual(values = c("#fe6766", "#993200", "black", "#68cb9b", "#ff9934", "#339866", "#3299fe", "#3266cb",  "red"))+
  labs(colour="")+
  scale_x_continuous(name="Years", breaks=seq(1940,2020, 10))+
  scale_y_continuous("Number of exp./meta-analyses", sec.axis = sec_axis(~ . , name = "Number of exp./meta-analyses"))

ggsave("NB_MA_exp_per_year.pdf", width = 20, height = 20, units = "cm")



    # qualité moyenne des meta-analyses
ES_TMP$Score<-as.numeric(as.character(ES_TMP$Score))
ES_TMP %>%
      ungroup() %>%
      group_by(Sub_Cat_Output) %>% 
  ggplot()+ geom_boxplot(aes(x=reorder(DIVERS, Score), y=Score, fill=DIVERS))+
  theme_pubr(base_size=24)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_colour_manual(values = c("#fe6766", "#993200", "black", "#68cb9b", "#ff9934", "#339866", "#3299fe", "#3266cb",  "red"))+
  theme(legend.position='none')+ labs(x="", y="Quality Score (%)")

ggsave("Qual_per_DIVERS.pdf", width = 20, height = 20, units = "cm")


# qualité moyenne des meta-analyses par année
ES_TMP %>%
  ungroup() %>%
  group_by(Publication_Year) %>% 
  ggplot()+ geom_boxplot(aes(x=Publication_Year, y=Score, group= Publication_Year), alpha=0.5)+
  theme_pubr(base_size=24)+
  #scale_color_manual(values = c("#00AFBB"))+
  geom_smooth(aes(x=Publication_Year, y=Score),se=FALSE)+
  theme(axis.text.x = element_text(angle = 00, hjust = 0.5))+
  #scale_colour_manual(values = c("#fe6766", "#993200", "black", "#68cb9b", "#ff9934", "#339866", "#3299fe", "#3266cb",  "red"))+
  theme(legend.position='none')+ labs(y="Quality Score (%)", x="")

ggsave("Qual_per_Year.pdf", width = 15, height = 10, units = "cm")

# Par type de diversification
    
    # Nombre d'effect size
mid<-30

TAB<-ES_TMP %>%
      ungroup() %>%
  filter(!DIVERS=="2 AND MORE") %>% 
      group_by(Sub_Cat_Output, DIVERS) %>% count() %>% 
  filter(!Sub_Cat_Output%in% c("yield and quality", "inputs use efficiency", "root biomass", "yield stability"))

TAB$Sub_Cat_Output<-factor(TAB$Sub_Cat_Output, 
                           levels=c("Production", "pests and diseases regulation",
                                    "products quality", "soil quality",
                                    "biodiv. assoc.", "greenhouses gas emission",
                                    "water quality", "water use",
                                    "profitability"))
levels(TAB$Sub_Cat_Output)[2] <- "Pest regulation"
levels(TAB$Sub_Cat_Output)[6] <- "GHG"

TAB%>%  
  ggplot()+
  geom_tile(aes(y= reorder(DIVERS,n), 
                x= Sub_Cat_Output, fill= n))+
  theme_pubr(base_size=22) +
  geom_text(aes(y= DIVERS, x= Sub_Cat_Output,label = round(n, 1)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_gradient2(midpoint=mid, low="cornsilk", mid="burlywood2",
                       high="darkgoldenrod1", space ="Lab" )+
  labs(x="",y="")+  theme(legend.position='none')
    
ggsave("Matrice_Evidence_map_ES.pdf", width = 20, height = 15, units = "cm")


# Nombre d'effect sizes: FOCUS ON SOIL QUALITY
mid<-30

TAB<-ES_TMP %>%
  ungroup() %>%
  filter(Sub_Cat_Output== "pests and diseases regulation") %>% 
  filter(!DIVERS=="2 AND MORE") %>% 
  group_by(Sub_sub_Cat_Output, DIVERS) %>% count() 
# %>% 
#   filter(!Sub_Cat_Output%in% c("yield and quality", "inputs use efficiency", "root biomass", "yield stability"))

TAB<-TAB[-c(7,16,17),]

TAB%>%  
  ggplot()+
  geom_tile(aes(y= reorder(DIVERS,n), 
                x= Sub_sub_Cat_Output, fill= n))+
  theme_pubr(base_size=22) +
  geom_text(aes(y= DIVERS, x= Sub_sub_Cat_Output,label = round(n, 1)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_gradient2(midpoint=mid, low="cornsilk", mid="burlywood2",
                       high="darkgoldenrod1", space ="Lab" )+
  labs(x="",y="")+  theme(legend.position='none')

ggsave("Matrice_Evidence_mapPEST_ES.pdf", width = 30, height = 15, units = "cm")


    # Nombre de meta-analyses
mid<-8

ES_TMP %>%
      ungroup() %>%
  filter(!DIVERS=="2 AND MORE") %>% 
      group_by(Sub_Cat_Output,DIVERS) %>% summarise(Unique_Elements = n_distinct(ID))%>%
      arrange(desc(Unique_Elements)) %>% 
   ggplot()+
  geom_tile(aes(x= reorder(DIVERS,Unique_Elements), 
                y= reorder(Sub_Cat_Output,Unique_Elements), fill= Unique_Elements))+
  theme_pubr(base_size=16) +
  geom_text(aes(x= DIVERS, y= Sub_Cat_Output,label = round(Unique_Elements, 1)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_gradient2(midpoint=mid, low="cornsilk", mid="burlywood2",
                         high="darkgoldenrod1", space ="Lab" )+
  labs(x="",y="")+  theme(legend.position='none')
  
ggsave("Matrice_Evidence_map_MA.pdf", width = 20, height = 20, units = "cm")


 ## carte


## coccurence et carte du monde
library(tidyverse)
library(tidygraph)
library(ggraph)
library('maps')
pdf("Map_WORLD.pdf", width = 10, height = 10)

map("world", col="grey80", fill=TRUE, bg="white", lwd=0.1)

country_centroids <- read.csv("~/Dropbox/A - Evidence map/Primary studies/plot/country_centroids_az8.csv", sep=";")
country_centroids  %<>% filter(type %in% c("Sovereign country", "Country"))

names(country_centroids)[15]<- "Country"
country_centroids$Country<-recode(country_centroids$Country, "United Kindom" = "UK")
country_centroids$Country<-recode(country_centroids$Country, "United States of America" = "USA")

COUNT<-IND_STU %>%                    
  group_by(Country) %>%     
  summarise(NB = n_distinct(DOI)) 

COUNT2<- merge(COUNT, country_centroids, by= "Country")

points(x=COUNT2$Longitude, y=COUNT2$Latitude,
       pch=21,cex=sqrt(COUNT2$NB/30),  bg = "orange", col = "black")

dev.off()

## map 2
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
names(world)[4]<- "Country"
world$Country<-recode(world$Country, "United Kindom" = "UK")
world$Country<-recode(world$Country, "United States of America" = "USA")


ggplot(data = world) +
  geom_sf(color = "black",fill= "antiquewhite")

MAP2 <- left_join(world, COUNT, by = "Country")
ggplot(data = MAP2) +
  geom_sf(color = "black",aes(fill= NB))+
  scale_fill_viridis_c(trans = "log", alpha = .9,option="inferno")+
  theme_pubr()+
  theme(legend.key.size = unit(0.5, "cm"),
    legend.key.width = unit(2.0,"cm") 
  )

ggsave("MAP_viridis.pdf", width = 20, height = 20, units = "cm")

