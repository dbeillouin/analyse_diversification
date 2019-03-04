# Damien BEillouin
# Méta-synthèse sur la divsersification des cultures: code d'analyse
# travail initial: Octobre 2018, reprise Février 2019

### Evidence map ####
### carte des études/ origine des auteurs ######

# Initialisation --------------------------------------------------------------------------------------------------

# Load packages
x<-c("magrittr", "tidyverse", "ggplot2","lme4","metafor","broom","ggpubr","reshape2","tidyr","dplyr","maps")
lapply(x, require, character.only = TRUE)


IND_STU <- read.csv("~/Documents/shinyapp/FINAL/Primary_Studies.csv", sep=";")
Des_MA  <- read.csv("~/Documents/shinyapp/FINAL/Description_Meta-Analyses.csv", sep=";", comment.char="#")
Type<-Des_MA[, c("ID","Diversification_types")]
IND_STU<-full_join(IND_STU,Type,by="ID") 


### Origine des auteurs
Des_MA$Country_first_author  <-as.character(Des_MA$Country_first_author)
Des_MA$Country_last_author   <-as.character(Des_MA$Country_last_author)
Des_MA$Country_other_authors <-as.character(Des_MA$Country_other_authors)

s <- strsplit(Des_MA$Country_first_author, split = ",")
FIRST<- data.frame(ID = rep(Des_MA$ID, sapply(s, length)), 
                     Reference = rep(Des_MA$Reference, sapply(s, length)),
                   Country = unlist(s))

s <- strsplit(Des_MA$Country_last_author, split = ",")
LAST<- data.frame(ID = rep(Des_MA$ID, sapply(s, length)), 
                   Reference = rep(Des_MA$Reference, sapply(s, length)),
                   Country = unlist(s))

s <- strsplit(Des_MA$Country_other_authors, split = ",")
OTHER<- data.frame(ID = rep(Des_MA$ID, sapply(s, length)), 
                  Reference = rep(Des_MA$Reference, sapply(s, length)),
                  Country = unlist(s))

AUTHORS <- rbind(rbind(FIRST, LAST),OTHER)

# supress white space
AUTHORS$Country <- base::trimws(AUTHORS$Country)
AUTHORS$Country <- tolower(AUTHORS$Country)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
AUTHORS$Country <- firstup(AUTHORS$Country)


COUNT<-AUTHORS %>%
  group_by(Country) %>%
  count()           %>% 
  arrange(desc(n))
COUNT[is.na(COUNT)] <- 0


## PARTIE 2 : nombre d'études individuelles dans chaque pays.
TMP<- IND_STU %>%
  dplyr::mutate(Country= str_trim(as.character(Country))) %>%
  dplyr::count(Country) %>%
  arrange(desc(n))
TMP$Country <- tolower(TMP$Country)
TMP$Country <- firstup(TMP$Country)

## ESSAI

TMP<- IND_STU %>%
  dplyr::mutate(Country= str_trim(as.character(Country))) %>%
  dplyr::group_by(Country,Year) %>%
  summarize(n= n_distinct(DOI)) %>% 
  arrange(desc(n))

TMP %<>% 
  mutate(GROUP = case_when(Country %in% EST_Mid_South_Africa ~ "EST_Mid_South_Africa",
                           Country %in% NORTH_WEST_Africa  ~ "NORTH_WEST_Africa",
                           Country %in% Estern_Asia   ~"Estern_Asia",
                           Country %in% WEstern_south_central_asia ~ "WEstern_south_central_asia",
                           Country %in% South_estern_Asia_POLY ~ "South_estern_Asia_POLY",
                           Country %in% Europe ~"Europe",
                           Country %in% CEN_SOUTH_AME ~"CEN_SOUTH_AME",
                           Country %in% North_America  ~ "North_America",
                           Country %in% AUS_NZ ~ "AUS_NZ"))
TMP <- TMP %>%  filter(!is.na(Year))


AA<-TMP %>% filter(Year>2005)
AA<-AA %>%  group_by(GROUP) %>% summarize(nn=sum(n))

BB<-TMP %>% filter(Year<2005)
BB<- BB %>%  group_by(GROUP) %>% summarize(nn=sum(n))

CC<-cbind(AA,BB)
names(CC)<-c("GOUP","npost","GROUP2","nbefo")
CC$ratio<- CC$npost/(CC$npost+CC$nbefo)


DD<-TMP %>% group_by(GROUP, Year) %>%   summarize(sum= sum(n)) %>%  group_by(GROUP) %>% mutate(cs = cumsum(sum))

ggplot(data=DD)+ geom_line(aes(x=Year, y=cs, group=GROUP, color=GROUP))+
geom_text_repel(
  data = DD %>% group_by(GROUP) %>%  filter(cs== max(cs)),
  aes(x=Year, y=cs,label = GROUP), # don't assign it here, 
  size = 3,
  nudge_x = 2250,
  direction="y",
  hjust=0,
  segment.color = 'grey80',
  colour = "black" # assign it here
)+ theme_pubr()+
  labs(x="Year", y="Number of primary studies") + theme(legend.position="none")

### PARTIE 3 représentation carte

map.world <- map_data("world")

map.world$region <- tolower(map.world$region)
map.world$region <- firstup(map.world$region)

map.world_joined <- left_join(map.world, TAB, by = c('region' = 'Country'))%>%
  mutate(fill_flg = ifelse(is.na(n),F,T)) %>%
  replace(., is.na(.), "") %>%
  mutate(n        = as.numeric(n))

centroids<- read.csv("~/Documents/analyse_diversification/analyse_diversification/country_centroids_az8.csv", sep=";", comment.char="#")
centroids %<>% select(admin, Longitude,Latitude)
DD<- centroids %>%  filter(admin=="France")

names(centroids)[1]<-"region"
INDI<- full_join(centroids,TMP, by = c('region' = 'Country'))

#  map.world_joined$n2<- cut(map.world_joined$n,5)

#         ggplot() +
#          geom_polygon(data = map.world_joined, aes(x = long, y = lat, group = group, fill = n2),color="black") +
#           scale_fill_brewer(palette="Spectral")+
#         #scale_fill_gradient(low = "white", high = "red")+
#          labs(title = paste0("Countries that have carried out meta-analysis"),#as.character(input$variable)),
#               subtitle = "Number of citations in the meta-analyses") +
#          theme(text = element_text( color = "black")
#                ,panel.background = element_rect(fill = "white")
#                ,plot.background = element_rect(fill = "white")
#                ,panel.grid = element_blank()
#                ,plot.title = element_text(size = 30)
#                ,plot.subtitle = element_text(size = 10)
#                ,axis.title = element_blank()
#                ,axis.ticks = element_blank())

###
library(viridis)

#         ggplot(map.world_joined, aes(x = long, y = lat, group = group,
#                               fill = cut_interval(n, 8))) +
#           geom_polygon(color = "grey10", size = 0.2) +
#           coord_equal() +
#           viridis::scale_fill_viridis(discrete = TRUE,direction=-1,option="D",name = "Number of \n primary studies:") +
#           labs(title = "",
#                fill = NULL) +
#           theme_void() +
#           theme(legend.position = "bottom",
#                 panel.background = element_rect(fill = NA, colour = "#cccccc"))+
#           theme(legend.text=element_text(size=20),legend.title=element_text(size=20))
#                                 ###


#         
p<-ggplot(aes(x = long, y = lat), data = map.world_joined) + 
  geom_polygon(aes(group = group, fill=n), colour = "grey65") + 
  geom_point(data=INDI,aes(Longitude, Latitude, size=n),shape=21, color="black", fill="blue")+
  scale_size(range = c(1, 12))+
  scale_y_continuous(limits = c(-60, 85)) + 
  coord_equal() +  theme_bw() +
  scale_fill_gradient(low = "white",high = "aquamarine4")+
  #  viridis::scale_fill_viridis(discrete = TRUE,direction=-1,option="D",name = "Number of \n primary studies:") +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(), 
        panel.border = element_rect(colour = "black"))
#         
p

DD<- INDI %>%  filter(region=="Usa")
