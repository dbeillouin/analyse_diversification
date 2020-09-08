## Redondancy

library("readxl");

PS<- read_excel("~/Documents/divers2020/PS_MAI.xlsx",sheet = "PS_MAI")
PS$ID<-tolower(PS$ID)

PS<- PS %>% distinct(ID, DOI) # or distinct(df, var1, var2)
COUNT<-PS %>% group_by(ID) %>% tally()
#PS<- PS %>% filter(ID%in% c("MA_A1", "MA_U1", "MA_B2", "MA_F1", "MA_E1"))

       w <- dcast(melt(PS), ID~DOI)
      # TG<-rowSums(w[,-1])
      # w<-x[,-93]
       x <- as.matrix(w[,-1])
       x[is.na(x)] <- 0
       x <- apply(x, 2,  function(x) as.numeric(x > 0))  #recode as 0/1
       v <- x %*% t(x)                                   #the magic matrix
       #diag(v) <- 0                                      #repalce diagonal
       dimnames(v) <- list(w[, 1], w[,1])                #name the dimensions

       # ESSAI 2

       Nb<- diag(v)
       NAMES<-rownames(v)
       TAB<-data.frame(v)
       colnames(TAB)<-NAMES


       TAB$Nb<-Nb

     #  colnames(TAB)[1]<-"NA"
      # TAB<- TAB %>%
       #  mutate_at(vars(starts_with("MA")), function (x) x/TAB$Nb*100)
       
       #TAB<-TAB %>% select(-Nb)

       TAB$MA<-  NAMES
 
       fwrite(TAB,'Matrice_red.csv')
       
       

 # print(TAB)
       TABM<-as.matrix(TAB)
       diag(TABM)<-0
       TABM<-data.frame(TABM)
       TAB2<- TABM %>% gather(MA)
       TAB2$MA_2<- rep (NAMES, length(TAB2$MA)/135)
     #  print(dim(TAB2))
       
       TAB2$value<-as.numeric(TAB2$value)
      TAB2<- TAB2 %>% filter(!MA=="Nb") %>% filter(value>45)
       
       ggplot(TAB2)+ geom_tile(aes(reorder(MA, value), reorder(MA_2, value), fill=value))+
          scale_fill_gradient2()+ theme_pubr()+
          theme(axis.text.x = element_text(angle = 90, hjust = 1))+
           theme(legend.position = "none") +
          labs(x="",y="")
       
       
       ##
       
       DES<- read_excel("~/Documents/Meta_synthesis_GLOBAL/Description_MA_2019_2020.xlsx")
       DES<- DES %>% select(ID, Publication_Year)
       
       names(DES)[1]<- "MA"
       TAB3<-merge(TAB2, DES, by='MA')
       ggplot(TAB3)+ geom_boxplot(aes(Publication_Year,value, group=Publication_Year ))+ 
          theme_pubr()
       
       
       
       QUA<- read_excel("~/Documents/Meta_synthesis_GLOBAL/Quality_2019SAUV_V3.xlsx")
       QUA<- QUA %>% dplyr::select(ID, Score)
       names(QUA)[1]<- "MA"
       TAB3<-merge(QUA, DES, by='MA')
       ggplot(TAB3)+ geom_boxplot(aes(Publication_Year,Score, group=Publication_Year ))+ 
          theme_pubr()
       
       
          geom_quantile(aes(Publication_Year,value),quantiles = q10,method = "rqss")
       
          geom_smooth(aes(Publication_Year,value), se=FALSE)
       
       MOY<-TAB2 %>% group_by(MA)  %>% 
          summarise_at(vars(value),
                       list(min=min, Q1=~quantile(., probs = 0.25),
                            median=median, Q3=~quantile(., probs = 0.75),
                            max=max))
       
       (mid=muted("blue"), high="red", low="blue")
       
       summary(TAB2$value)
       scale_fill_gradient2(low="orange",mid= "blue", high="red", limits=c(0, 100))
       
       
       names(TAB2)<- c("MA_1","Perc", "MA_2")
       TAB2<-TAB2[!TAB2$MA_2=="Nb",]
       TAB2$Perc<-as.numeric(TAB2$Perc)

      # TAB2<-TAB2[!TAB2$Perc>75,]
       ORDER<-TAB2 %>% group_by(MA_1)%>% summarize(MEAN = mean(Perc, na.rm = TRUE))%>% arrange(MEAN)
      TAB2<-full_join(TAB2,ORDER)

      NBR<-data.frame(TAB$MA,   TAB$Nb)
      names(NBR)<-c("MA_1","Nb")
      print(NBR)
      TAB2<-full_join(TAB2,NBR)

      TAB2$MA_1<-substr(TAB2$MA_1,4,8)
      print("ici")
      print(names(TAB2))
      library("ggrepel")
    ggplot(data=TAB2, aes(x=reorder(MA_1, -Nb), y=Perc))+ geom_point(size=2)+
         labs(title= " ",y="Percent of common individual studies \n
              between the meta-analysis in abscissa and the others",
              x="")+
         theme_pubr()+
         theme(legend.position="none",
               text = element_text(size=20),
               axis.text=element_text(colour="black"),
               axis.text.x = element_text(angle = 90, hjust = 1,size=30),
               axis.title.x = element_text(size=30)) +
              geom_text_repel(data= subset(TAB2, Perc > 15),
                              aes(x=MA_1,y=Perc,label = MA_2),label.size=8,size=8,
                              segment.color = "black", force=0.6) +
       geom_hline(yintercept=c(25,50,75), linetype="dashed", color="gray")+coord_flip()
    
    
    
    
    
    
    # get the path to the shapefile embedded in sf
    # package
    shapepath <- system.file("shape/nc.shp", 
                             package="sf")
    
    # import the shapefile
    nc <- sf::st_read(shapepath)
    # change the projection of the map
    nc <- sf::st_transform(nc, 3857)
    

    # data preparation
    # compute the sudden infant deaths per 1000 births
   # nc$share <- 100000 * nc$SID74 / nc$BIR74
    # quantization breaks of the rate
    # correct the breaks to use the global rate as limit of class 
    global_rate <- sum(nc$SID74) * 100000 / sum(nc$BIR74)
    bks[4] <- global_rate
    
    
    world <- ne_countries(scale = "medium", returnclass = "sf")
    world <- sf::st_transform(world, 3857)
    world<- world  %>% filter(!admin=="Antarctica")
    
    world$admin<-plyr::revalue(world$admin, c("United States of America" = "usa",
                 "United Kingdom" = "uk"))
    
    #levels(world$admin)[levels(world$admin)=="United States of America"] <- "usa"
    #levels(world$admin)[levels(world$admin)=="United Kingdom"] <- "uk"
    
  #  world <- cartogram_cont(world, "pop_est", itermax = 5)
    world = st_transform(world, crs = "+proj=moll")
    tm_shape(world) + tm_polygons("pop_est", style = "jenks") +
       tm_layout(frame = TRUE, legend.position = c("left", "bottom"))
    
    
    
    
    PS<- read_excel("~/Documents/divers2020/PS_MAI.xlsx",sheet = "PS_MAI")
    PS<- PS %>% distinct(DOI, Country) # or distinct(df, var1, var2)
    COUNT<-PS %>% group_by(Country) %>% tally()
    names(COUNT)[1]<- "admin"
    names(COUNT)[2]<- "NB_PS"
    COUNT$admin<-tolower(COUNT$admin)
    world$admin<-tolower(world$admin)
    world<- left_join(world, COUNT,by= "admin")
    world <- sf::st_transform(world, 3857)
    
    world <- cartogram_cont(world, "NB_PS", itermax = 5)
    world = st_transform(world, crs = "+proj=moll")
    tm_shape(world) + tm_polygons("NB_PS", style = "jenks", projection="+proj=robin") +
       tm_layout(frame = FALSE, legend.position = c("left", "bottom"))
    
  
    
    
    world$NB_PS <- replace(   world$NB_PS, is.na(   world$NB_PS), 0)
    summary(world$NB_PS)
    
    
    #####
    world_map <- map_data("world")
    
    
    PS<- read_excel("~/Documents/divers2020/PS_MAI.xlsx",sheet = "PS_MAI")
    PS<- PS %>% distinct(DOI, Country) # or distinct(df, var1, var2)
    
    PS$Country<- as.character(PS$Country)
    s <- strsplit(PS$Country, split = ",")
    PS<- data.frame(DOI = rep(PS$DOI, sapply(s, length)),
                     Country = unlist(s))
    
    PS$Country<-trimws(PS$Country)
    names(PS)[2]   <- "region"
    PS$region<-tolower(PS$region)
    COUNT<-PS %>% group_by(region) %>% tally()
    names(COUNT)[2]<- "NB_PS"
  
    world_map$region<-tolower(world_map$region)
    setdiff(COUNT$region,world_map$region)
    
    world_map<- left_join(world_map, COUNT,by= "region")
    
    world_map$cut_NB<- cut(world_map$NB_PS,breaks= c(0,20,100,200, 400, 1600))
    plot1<-ggplot(world_map, aes(x = long, y = lat, group = group, fill = cut_NB)) +
       geom_polygon(colour = "grey75") +
      # scale_y_continuous("", breaks = (-2:2) * 30) +
      # scale_x_continuous("", breaks = (-4:4) * 45) +
       scale_fill_viridis("", option = "inferno", direction = -1,discrete=TRUE) +
       theme(legend.position = "right",
             axis.text = element_blank()) +
       #coord_map("lambert", lat0=32.5343, lat=142.0095)
     coord_map("orthographic", orientation = c(10, 280, 0))+
        theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "gray80"))+
       labs(x="",y="")+ theme(legend.position = "none") 
       
  plot2<-ggplot(world_map, aes(x = long, y = lat, group = group, fill = cut_NB)) +
       geom_polygon(colour = "grey75") +
       # scale_y_continuous("", breaks = (-2:2) * 30) +
       # scale_x_continuous("", breaks = (-4:4) * 45) +
       scale_fill_viridis("", option = "inferno", direction = -1,discrete=TRUE) +
       theme(legend.position = "right",
             axis.text = element_blank()) +
       #coord_map("lambert", lat0=32.5343, lat=142.0095)
       coord_map("orthographic", orientation = c(10, 100, 0))+
       theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "gray80"))+
     labs(x="",y="")+ theme(legend.position = "none") 
    
  plot3<-
     ggplot(world_map, aes(x = long, y = lat, group = group, fill = cut_NB)) +
     geom_polygon(colour = "grey75") +
     # scale_y_continuous("", breaks = (-2:2) * 30) +
     # scale_x_continuous("", breaks = (-4:4) * 45) +
     scale_fill_viridis("", option = "inferno", direction = -1,,discrete=TRUE) +
     theme(legend.position = "right",
           axis.text = element_blank()) +
     #coord_map("lambert", lat0=32.5343, lat=142.0095)
     coord_map("orthographic", orientation = c(10, 30, 0))+
     theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "gray80"))+
     labs(x="",y="")+ theme(legend.position = "none") 
  
  
    
  plot_grid(plot1, plot2,plot3, labels = c('A', 'B', 'C'), label_size = 12, ncol=3)
    
  
  
  
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
  
  COUNT %<>% 
     mutate(GROUP = case_when(region %in% tolower(EST_Mid_South_Africa) ~ "EST_Mid_South_Africa",
                              region %in% tolower(NORTH_WEST_Africa)  ~ "NORTH_WEST_Africa",
                              region %in% tolower(Estern_Asia)   ~"Estern_Asia",
                              region %in% tolower(WEstern_south_central_asia) ~ "WEstern_south_central_asia",
                              region %in% tolower(South_estern_Asia_POLY) ~ "South_estern_Asia_POLY",
                              region %in% tolower(Europe) ~"Europe",
                              region %in% tolower(CEN_SOUTH_AME) ~"CEN_SOUTH_AME",
                              region %in% tolower(North_America)  ~ "North_America",
                              region %in% tolower(AUS_NZ) ~ "AUS_NZ"))
  
 
  
  library(ggplot2) 
  library(treemapify)
  COUNT2<- COUNT %>% filter(!is.na(GROUP))
  ggplot(COUNT2, aes(area = NB_PS, fill = GROUP,label = substr(region, 1,10), 
                     subgroup=GROUP)) +
     geom_treemap()+ theme(legend.position = "none") +
     geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                       grow = TRUE)+
     geom_treemap_subgroup_border(colour = "blue", size = 1) 
     #geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.5, colour =
      #                             "black", fontface = "italic", min.size = 0)
  
  
  
  
  
  
  PS<- read_excel("~/Documents/divers2020/PS_MAI.xlsx",sheet = "PS_MAI")
  PS<- PS %>% distinct(Country, DOI, divers_MAJ) # or distinct(df, var1, var2)
  
  # PS$Country<- as.character(PS$Country)
  # s <- strsplit(PS$Country, split = ",")
  # PS<- data.frame(DOI = rep(PS$DOI, sapply(s, length)),
  #                 Country = unlist(s))
  
  PS$divers_MAJ<-tolower(PS$divers_MAJ)
  PS$divers_MAJ<-trimws(PS$divers_MAJ)
  names(PS)[2]   <- "region"
  PS$region<-tolower(PS$region)
  COUNT<-PS %>% group_by(divers_MAJ) %>% tally()
  names(COUNT)[2]<- "NB_PS"
  
  COUNT<- COUNT %>% filter(divers_MAJ %in% c("rotation", "associated plant species", "intercropping", "agroforestry", "mixture"))
  ggplot(COUNT, aes(area = NB_PS, fill = divers_MAJ,label = divers_MAJ)) +
     geom_treemap()+ theme(legend.position = "none") +
     geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                       grow = TRUE)
  

  PS<- read_excel("~/Documents/divers2020/PS_MAI.xlsx",sheet = "PS_MAI")
  PS<-PS  %>% filter(divers_MAJ=="rotation")
  PS<- PS %>% distinct(DOI, ID) # or distinct(df, var1, var2)

  ggplot(as.data.frame(UCBAdmissions),
         aes(y = Freq, axis1 = Gender, axis2 = Dept)) +
     geom_alluvium(aes(fill = Admit), width = 1/12) +
     geom_stratum(width = 1/12, fill = "black", color = "grey") +
     geom_label(stat = "stratum", infer.label = TRUE) +
     scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
     scale_fill_brewer(type = "qual", palette = "Set1") +
     ggtitle("UC Berkeley admissions and rejections, by sex and department")
  
  
  +head(UCB_lodes, n = 12)
  
  
  COUNT<-PS %>% group_by(DOI, ID) %>% tally()
  ggplot(as.data.frame(COUNT),
         aes(y = n, axis1 = ID, axis2 = DOI)) +
     geom_alluvium(aes(fill = ID,colour= ID), alpha=0.1, width = 1/40, decreasing = TRUE)+
     geom_stratum(width = 1/40, fill = "black", color = "grey")+
     scale_color_viridis(discrete = TRUE,option = "D")+
     #scale_colour_grey()+
     theme_pubr()+ theme(legend.position='none')
  
  +
     geom_label(stat = "stratum", infer.label = TRUE) 
  +
     scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
     scale_fill_brewer(type = "qual", palette = "Set1") +
     ggtitle("UC Berkeley admissions and rejections, by sex and department")
  
  
  PS<- read_excel("~/Documents/divers2020/PS_MAI.xlsx",sheet = "PS_MAI")
  PS<- PS %>% distinct(Country, DOI, divers_MAJ) # or distinct(df, var1, var2)
  
  # PS$Country<- as.character(PS$Country)
  # s <- strsplit(PS$Country, split = ",")
  # PS<- data.frame(DOI = rep(PS$DOI, sapply(s, length)),
  #                 Country = unlist(s))
  
  PS$divers_MAJ<-tolower(PS$divers_MAJ)
  PS$divers_MAJ<-trimws(PS$divers_MAJ)
  names(PS)[2]   <- "region"
  PS$region<-tolower(PS$region)
 
  PS %<>% 
     mutate(GROUP = case_when(region %in% tolower(EST_Mid_South_Africa) ~ "AFR",
                              region %in% tolower(NORTH_WEST_Africa)  ~ "AFR",
                              region %in% tolower(Estern_Asia)   ~"AS",
                              region %in% tolower(WEstern_south_central_asia) ~ "AS",
                              region %in% tolower(South_estern_Asia_POLY) ~ "AS",
                              region %in% tolower(Europe) ~"Europe",
                              region %in% tolower(CEN_SOUTH_AME) ~"AM",
                              region %in% tolower(North_America)  ~ "AM",
                              region %in% tolower(AUS_NZ) ~ "AS"))
  
  COUNT<-PS %>% group_by(divers_MAJ, GROUP) %>% tally()
  names(COUNT)[3]<- "NB_PS"
  
  COUNT<- COUNT %>% filter(divers_MAJ %in% c("rotation", "associated plant species", "intercropping", "agroforestry", "mixture"))
  AM<- COUNT %>% filter(GROUP =="AFR")
  ggplot(AM, aes(area = NB_PS, fill = divers_MAJ,label = divers_MAJ)) +
     geom_treemap()+ theme(legend.position = "none") +
     geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                       grow = TRUE)
     
  
  
  
  
  head(G20)
  
    nc<-world
    
    world$share <- world$NB_PS / world$pop_est
    
    bks <- getBreaks(v = world$share, method = "quantile", nclass = 5)
    
    # get a color palette
    cols <- carto.pal(pal1 = "green.pal", n1 = 3, 
                      pal2 = "wine.pal", n2 = 2)
    ## Choropleth layer
    
    # set figure margins and background color
    par(mar = c(0,0,1.2,0), bg = "lemonchiffon")
    # display the sudden infant deaths per 1000 births
    choroLayer(x = nc_cont,var = "share",   breaks = bks, col = cols,
               border = "khaki", lwd = 0.5, 
               legend.title.txt = "AAAAA", 
               legend.pos = 'topleft', legend.values.rnd = 0)
# add a title and layout
    layoutLayer(title = "Crop diversification around the world", 
                sources = "", north = TRUE, scale = 50,tabtitle = TRUE,
                theme = "sand.pal", frame = FALSE,  
                author = "*per 100,000 live births. Source: North Carolina SIDS data set")

    nc_cont <- cartogram_cont(world, weight="NB_PS",itermax = 2)    
    
    nc_dorling <- cartogram_dorling(world, "NB_PS", 
                                    k = 5)
    
    plot(nc_dorling)
    
    
    
    
    mtq_pencil <- getPencilLayer(x = plot1)
    # display this MULTILINESTRING layer
    plot(st_geometry(mtq_pencil), col = 1:8)
    plot(st_geometry(nc_cont), col = NA, add=T)
    # and a add the original borders
    
    
    
    # package
    shapepath <- system.file("shape/nc.shp", 
                             package="sf")
    # import the shapefile
    nc <- sf::st_read(shapepath)
    # change the projection of the map
    nc <- sf::st_transform(nc, 3857)
    names(nc)[9]<- "NB_PS"
    names(nc)[10]<- "pop_est"
    
    
    # data preparation
    # compute the sudden infant deaths per 1000 births
    nc$share <- 100000 * nc$pop_est / (nc$NB_PS+1)
    nc$share <- replace(  nc$share, is.na(  nc$share), 3)
    
    
    # quantization breaks of the rate
    bks <- getBreaks(v = nc$share, method = "quantile", nclass = 5)
    # correct the breaks to use the global rate as limit of class 
    # get a color palette
    cols <- carto.pal(pal1 = "green.pal", n1 = 3, 
                      pal2 = "wine.pal", n2 = 2)
    ## Choropleth layer
    
    
    
    # set figure margins and background color
    par(mar = c(0,0,1.2,0), bg = "lemonchiffon")
    # display the sudden infant deaths per 1000 births
    choroLayer(x = nc_cont,var = "share", breaks = bks, col = cols,
               border = "khaki", lwd = 0.5, 
               legend.title.txt = "Sudden infant death syndrome rate*", 
               legend.pos = 'topleft', legend.values.rnd = 0)
    # add a title and layout
    layoutLayer(title = "Sudden Infant Death Syndrome in North Carolina, 1974-1978", 
                sources = "", north = TRUE, scale = 50,tabtitle = TRUE,
                theme = "sand.pal", frame = FALSE,  
                author = "*per 100,000 live births. Source: North Carolina SIDS data set")
    
    
    nc_cont <- cartogram_cont(nc, "pop_est",itermax = 2)
    

    
    
    
    
    library(cartogram)
    library(tmap)
    library(maptools)
    #> Loading required package: sp
    #> Checking rgeos availability: TRUE
    library(cartogram)
    library(tmap)
    library(maptools)
    #> Loading required package: sp
    #> Checking rgeos availability: TRUE
    
    data(wrld_simpl)
    ss<-wrld_simpl@data
    names(COUNT)[1]<- "NAME"
    ss$NAME<-tolower(ss$NAME)
    ss<-left_join(ss, COUNT)
    ss$NB_PS <- replace(   ss$NB_PS, is.na(   ss$NB_PS), 1)
    
    
    afr <- wrld_simpl[wrld_simpl$REGION == 2, ]
    # project the map
    
    
    afr <- spTransform(afr, CRS("+init=epsg:3395"))
    # construct cartogram
    afr_cont <- cartogram_cont(afr, "POP2005", itermax = 5)

    afr$POP2005
    afr$NB_PS<- ss$NB_PS
    # plot it
    tm_shape(afr_cont) + tm_polygons("POP2005", style = "jenks") +
       tm_layout(frame = FALSE, legend.position = c("left", "bottom"))
    
    