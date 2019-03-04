# Damien BEillouin
# Méta-synthèse sur la divsersification des cultures: code d'analyse
# travail initial: Octobre 2018, reprise Février 2019

### Deuxième type de représentation graphique ######
# Comme dans le dataPaper ######

# Initialisation ---------------------------------------------------------------------------------------
    
    setwd("~/Documents/shinyapp/")
    
    source("~/Documents/shinyapp/Load_shiny.R")
    QUAL<-read.csv("Quality_16_10.csv", sep=";", encoding="UTF-8-MAC", comment.char="#")%>%
      dplyr::select(c(MA,ID,Score))
    QUAL <- QUAL[-c(98,99),]
    
    Type<-Des_MA[, c("ID","Diversification_types")]
    QUAL<-full_join(QUAL,Type,by="ID") 
    
    
    TMP<-ES   %>% 
      filter(KEEP=="YES") %>%
      left_join(.,QUAL %>% dplyr::select(ID,Score))%>% 
      # filter(grepl(as.character( input$variable), Diversification_types))%>%
      filter(grepl(as.character("Agroforestry"), DIVERS.1))%>%
      filter(!is.na(ES)) %>%
      dplyr::select(ID,Nb_Ind_Studies,Nb_Paired_Data,ES,ES_lower_CI,ES_upper_CI,Metric_Output,Main_Cat_Output,Score,Sub_Cat_Output,Sub_Sub_Cat_Output)%>%
      mutate(pairs= as.numeric(as.character(Nb_Paired_Data))) %>%
      mutate(ES = ifelse(Metric_Output %in% c("percent difference") & abs(ES)>1, ES/100, ES),
             ES_lower_CI = ifelse(Metric_Output %in% c("percent difference") & abs(ES_lower_CI)>1, ES_lower_CI/100, ES_lower_CI),
             ES_upper_CI = ifelse(Metric_Output %in% c("percent difference") & abs(ES_upper_CI)>1, ES_upper_CI/100, ES_upper_CI)) %>%
      mutate(ES = ifelse(Metric_Output %in% c("percent difference"), 1-ES, ES),
             ES_lower_CI = ifelse(Metric_Output %in% c("percent difference"), 1-ES_lower_CI, ES_lower_CI),
             ES_upper_CI = ifelse(Metric_Output %in% c("percent difference"), 1-ES_upper_CI, ES_upper_CI))
    
    
    Hedge<-TMP %>%
      filter(Metric_Output %in% c("Hedge's d'")) %>%  
      setNames(ifelse(grepl('ES', names(.)),paste("Hedge",gsub(".*S_","",names(TMP)),sep="_"),names(.))) 
    
    Ratio<-TMP %>%
      filter(Metric_Output %in% c("percent difference","Ratio","log(Ratio)")) %>%  
      setNames(ifelse(grepl('ES', names(.)),paste("Ratio",gsub(".*S_","",names(TMP)),sep="_"),names(.)))
    
    
    autres<-TMP %>%
      filter(!Metric_Output %in% c("Hedge","percent difference","Ratio"))
    
    TMP<-full_join(Hedge,Ratio) %>%
      arrange(Main_Cat_Output,Sub_Cat_Output,Hedge_ES,Ratio_ES)
    
    TT<-TMP %>% group_by(Sub_Cat_Output)%>% count()
    TT<-TT[TT$n>0,]
    TMP<-TMP[TMP$Sub_Cat_Output %in% TT$Sub_Cat_Output,]
    

# plot-------------------------------------------------------------------------------------------------
    
    # On défini les paramètres graphiques
    
    MIN <-min(-2, min(as.numeric(as.character(TMP$Ratio_lower_CI)),na.rm=TRUE))
    MAX <-max(2, max(as.numeric(as.character(TMP$Ratio_upper_CI)),na.rm=TRUE))
    RANGE<-MAX-MIN
    
    VECT<-c(TMP$Hedge_ES,TMP$Hedge_lower_CI, TMP$Hedge_upper_CI)
    VECT2<-scales::rescale(VECT, to = c(MIN, MAX))
    
    TMP$Hedge_ES2<-VECT2[1: length(TMP$Hedge_ES)]
    TMP$Hedge_lower_CI2<-VECT2[(length(TMP$Hedge_lower_CI)+1):(2*length(TMP$Hedge_lower_CI))]
    TMP$Hedge_upper_CI2<-VECT2[(2*length(TMP$Hedge_upper_CI)+1):(3*length(TMP$Hedge_upper_CI))]
    
    MIN_AXE_DROITE  <- min(-0.5,min(as.numeric(as.character(TMP$Hedge_lower_CI)),na.rm=TRUE))
    #min(MIN, min(TMP$Hedge_lower_CI2-DECALAGE_POINTS,na.rm=TRUE))
    MAX_AXE_DROITE  <- max(0.5,max(as.numeric(as.character(TMP$Hedge_upper_CI)),na.rm=TRUE))
    RANGE_AXE_DROITE<-  MAX_AXE_DROITE- MIN_AXE_DROITE
    
    DECALAGE_POINTS  <-  ((0-MIN_AXE_DROITE)/RANGE_AXE_DROITE*RANGE+MIN-1)
    DECALAGE_ECHELLE <-  abs(1-MIN)/RANGE*RANGE_AXE_DROITE+MIN_AXE_DROITE
    
    #THEME
    ###### 
    ThemeDB<-     theme_pubr()+
      theme(plot.title = element_text(lineheight=3, face="bold", color="black",size=35),
            legend.text=element_text(size=35),
            legend.title=element_text(size=35),
            plot.margin = unit(c(1,5,1,1), "lines"),
            axis.text.x=element_text(angle=90,colour="grey20",face="bold",size=15),
            axis.text.y=element_text(colour="grey20",face="bold",hjust=1,vjust=0.8,size=20),
            axis.title.x=element_text(colour="grey20",face="bold",size=20),
            axis.title.y=element_text(colour="grey20",face="bold",size=20)) 
    ######  
    
    library(pBrackets) 
    DBFUN<-function(X){
      NEW<-X
      # X2<-cumsum(X)
      NEW[1]<-X[1]/2
      for (i in 2: length(X)){
        NEW[i]<-X[i-1]+(X[i]-X[i-1])/2}
      return(as.numeric(NEW))
    }
    
    addline_format <- function(x,...){
      gsub('\\s','\n',x)
    }
    
    
    #########################
    
    TMP$ORDER<-1:length(TMP$ID)
    
    # TMP<-TMP %>% filter(Type2=="Biomass and grain production")
    
    DBFUN<-function(X){
      NEW<-X
      # X2<-cumsum(X)
      NEW[1]<-X[1]/2
      if (length(X)>2){
        for (i in 2: length(X)){
          NEW[i]<-X[i-1]+(X[i]-X[i-1])/2}} else {}
      return(as.numeric(NEW))
    }
    
    TMP$Nb_Paired_Data[is.na(TMP$Nb_Paired_Data)]<-"-"
    TMP$Nb_Ind_Studies[is.na(TMP$Nb_Ind_Studies)]<-"-"
    p<- ggplot(data = TMP) + 
      geom_point(aes(x = ORDER, y =Ratio_ES, group=interaction(ID,Main_Cat_Output) ,
                     fill     =  Score),
                 size     =  1,
                 shape    =  1,
                 stat     =  "identity",
                 alpha    =  0.6,
                 color    =   "black") +
      geom_text(aes(label=paste0(Nb_Ind_Studies," (",Nb_Paired_Data,")"), y=-1.7, x=ORDER),hjust=+0.1, vjust=0, angle=90)+
      geom_errorbar(aes(x=ORDER, ymin=Ratio_lower_CI, ymax=Ratio_upper_CI), 
                    width    = 0.2, 
                    size     = 0.5, 
                    color    = "black")+
      geom_hline(yintercept = 0,linetype=2,col="gray")+
      ylab("Size effect (log(Ratio))")+ xlab('')+
      ThemeDB
    
    
    TMP$Main_Cat_Output<-factor(TMP$Main_Cat_Output)
    COUNT  <-cumsum(plyr::count(TMP$Main_Cat_Output)[2])
    TMP<-TMP %>% mutate(Sub_Cat_Output = factor(Sub_Cat_Output, Sub_Cat_Output))
    TMP$Sub_Cat_Output<-factor(TMP$Sub_Cat_Output)
    COUNT2 <-cumsum(arrange(plyr::count(TMP$Sub_Cat_Output),unique(TMP$Sub_Cat_Output))[2])
    
    #get the ylim and xlim
    xmin <- min(ggplot_build(p)$layout$panel_ranges[[1]]$x.range) 
    xmax <- max(ggplot_build(p)$layout$panel_ranges[[1]]$x.range)
    ymin <- min(ggplot_build(p)$layout$panel_ranges[[1]]$y.range)
    ymax <- max(ggplot_build(p)$layout$panel_ranges[[1]]$y.range)
    
    
    
    AGROF<-p+ # geom_vline(xintercept=unlist(COUNT[,-dim(COUNT)[1]])+0.5,linetype=2,col="blue")     +
      geom_vline(xintercept=unlist(COUNT2)+0.5,linetype=3,col="black")     +
      # annotate("text", size=4,x = DBFUN(COUNT2$freq), y = 1, # MIN+0.05*RANGE, 
      #         label = addline_format(unique(TMP$Sub_Cat_Output),angle=0))+
      theme(axis.text.y.right = element_text(color = "red"))+
      xlab("")+theme(axis.text.x = element_blank(),
                     axis.title.x = element_text(size=38))+
      theme(legend.position='none')+
      theme(axis.title.y =element_text(size=20),
            axis.title.x =element_text(size=20))
    
    
    
    plot_grid(ROT,AGROF,COVER,align=TRUE,ncol=1)