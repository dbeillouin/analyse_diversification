# Damien BEillouin
# Méta-synthèse sur la divsersification des cultures: code d'analyse
# travail initial: Octobre 2018, reprise Février 2019

### ANALYSE DES DONNEES TOUT TYPE DE DIVERSIFICATION ENSEMBLE ######
## article Evidence Map #####

#Plot espèces en focntion des région et outputs en fonction des régions. 


Primary_Studies <- read.csv("~/Documents/shinyapp/FINAL/Primary_Studies.csv", sep=";")

Primary_Studies$Species<- as.character(Primary_Studies$Species)
s <- strsplit(Primary_Studies$Species, split = ";")
SPECIES<- data.frame(ID = rep(Primary_Studies$ID, sapply(s, length)), 
                   Reference = rep(Primary_Studies$Reference, sapply(s, length)),
                   Country = rep(Primary_Studies$Country, sapply(s, length)),
                   Species = unlist(s))


### 

Melanasia<- c("Fiji", "New Caledonia", "Papua New Guinea","Timor-Leste")
Micronesia <- c("Guam", "Kiribati","Palau", "Nauru")
Polynesia<- c("Samoa", "Polynesia", "Wallis")

Eastern_Africa <- c("Burundi", "Comoros", "Djibouti", "Eritrea", "Ethiopia", "Kenya", "Madagascar", "Malawi", "Mauritus", "Mozambique", "Rwanda", "Somalia", "Uganda", "Tanzania", "Zambia", "Zimbabwe" )
Middle_africa <- c("Angola", "Cameroon", "Central African Republic", "Chad", "Republic of Congo", "Guinea", "Gabon", "Sao Tome")
Northern_Africa<- c("Algeria", "Egypt", "Tunisia", "Lybia", "Morocco", "Sudan", "Western sahara")
Southern_africa <-c("Bostwana", "Lesotho", "Namibia", "South Africa", "Swaziland")
Western_Africa <- c( "Benin", "Burkina Faso", "Cape verde", "Ivory Coast", "Gambia", "Ghana", "Guinea", "Liberia", "Mali", "Mauritania", "Niger", "Nigeria", "Senegal", "Sierra Leone", "Togo")
###
#regouper Esterm and Middle and southern:
EST_Mid_South_Africa <- c(Eastern_Africa,Middle_africa,Southern_africa)
NORTH_WEST_Africa <- c(Northern_Africa,Western_Africa)


Estern_Asia <- c("China", "South Korea", "Japan", "Mongolia")
South_Central_Asia <- c("Afghanistan", "Bangladesh", "Bhutan", "India", "Iran", "Kazakhstan", "Kyrgyzstan", "Maldives", "Nepal", "Pakistan", "Sri Lanka", "Tajikistan","Turkmenistan","Uzbekistan" )
South_estern_Asia<- c("Cambodia", "Timor", "Indonesia","Malaysia", "Myanmar","Philippines", "Singapore", "Thailand", "Vietnam")
Western_Asia <- c("Armenia", "Azerbaijan", "Bahrain", "Cyprus", "Syria","Iraq", "Israel", "Jordan", "Kuwait", "Lebanon", "Oman", "Qatar", "Saudi Arabia", "Turkey", "Yemen")
# regrouper Western Asia and SOuth central Asia
WEstern_south_central_asia <- c(South_Central_Asia,Western_Asia)
# regrouper South_estern_Asia et micromnésie,...
South_estern_Asia_POLY<- c(South_estern_Asia,Melanasia,Micronesia,Polynesia)

Eastern_europe<- c("Belarus", "Bulgaria","Romania", "Czech Republic", "Hungary", "Poland", "Moldava", "Russia", "Slovakia", "Ukraine")
Northern_europe<- c("Denmark", "Estonia", "Finland", "Iceland", "Ireland", "Latvia", "Lithunia", "Norway", "Sweden", "UK" )
Southern_europe <- c("Albania", "Andorra", "Croatia", "Serbia", "Greece", "Italy", "Malta", "Portugal", "San Marino", "Slovenia","Spain" )
Western_europe <- c("Austria", "Belgium", "Switzerland","France", "Germany", "Luxembourg", "Monaco", "Netherlands","Swizerland")
# regrouper Europe
Europe<- c(Eastern_europe,Northern_europe,Southern_europe,Western_europe)


Caribbean<-c ( "Bahamas", "Cuba", "Dominican Republic", "Grenada", "Jamaica", "Puerto Rico","Trinidad")
Central_america <- c("Belize", "Costa Rica", "El Salvador", "Guatemala", "Honduras", "Mexico", "Nicaragua", "Panama")
South_America <- c( "Argentina", "Bolivia", "Brazil", "Chile", "Colombia", "Ecuador", "Guyana","Paraguay", "Peru", "Suriname", "Uruguay", "Venezuela")
# regrouper
CEN_SOUTH_AME<- c(Central_america,South_America,Caribbean)

North_America <-c("USA", "Canada", "Bermuda")
AUS_NZ<-c ("Australia", "New Zealand")


## groupes$

SPECIES %<>% 
  mutate(GROUP = case_when(Country %in% EST_Mid_South_Africa ~ "EST_Mid_South_Africa",
                           Country %in% NORTH_WEST_Africa  ~ "NORTH_WEST_Africa",
                           Country %in% Estern_Asia   ~"Estern_Asia",
                           Country %in% WEstern_south_central_asia ~ "WEstern_south_central_asia",
                           Country %in% South_estern_Asia_POLY ~ "South_estern_Asia_POLY",
                           Country %in% Europe ~"Europe",
                           Country %in% CEN_SOUTH_AME ~"CEN_SOUTH_AME",
                           Country %in% North_America  ~ "North_America",
                           Country %in% AUS_NZ ~ "AUS_NZ"))

SPECIES <- SPECIES %>%  filter(!is.na(Species))
SPECIES$Species<- base::trimws(SPECIES$Species)

COUNT <- SPECIES           %>% 
  group_by(GROUP, Species) %>% 
  count()                  %>% 
  group_by(GROUP)        %>% 
  mutate(rank = rank(desc(n), ties.method = "first"))

COUNT <- COUNT %>%  filter(rank<6)

LIST <- c("EST_Mid_South_Africa","NORTH_WEST_Africa","Estern_Asia","WEstern_south_central_asia","South_estern_Asia_POLY","Europe","CEN_SOUTH_AME","North_America", "AUS_NZ")
PLOT<- NULL
for (i in 1: length(LIST)) {
PLOT[[i]] <- ggplot(COUNT %>% filter(GROUP== LIST[i])) +
  geom_bar(aes(x=reorder(Species, desc(n)),y=n), stat="identity")+
  facet_wrap(~GROUP)+
  xlab("")+ ylab("Number od primary studies")
}

library("cowplot", lib.loc="~/Library/R/3.5/library")
plot_grid(PLOT[[1]], PLOT[[2]],PLOT[[3]],PLOT[[4]],PLOT[[5]],PLOT[[6]],PLOT[[7]],PLOT[[8]],PLOT[[9]], labels = c("A", "B","D","E","F","G","H","I","J"))

# même chose avec les outcomes
 #JE pense que ce n'est pas possibl'