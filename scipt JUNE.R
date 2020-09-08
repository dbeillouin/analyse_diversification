
  
  #Load global matric of redundancy betwwen MA
MAT<-fread("Matrice_red.csv")


### Analysis by  DIVERS AND OUTPUT:

# plot généraux
AA<- RATIO %>%
  ungroup()       %>%
  mutate(DIVERS = tolower(DIVERS)) %>% 
  filter(!DIVERS %in% c("others","2 and more")) %>% 
  filter(Sub_Cat_Output %in% c("Production","biodiv. assoc.")) %>% 
  group_by(Sub_Cat_Output) %>% 
  mutate(positionInCategory = 1:n()) %>% 
  arrange(Sub_Cat_Output, ES) %>%
  mutate(order = row_number())

ggplot(AA)+ geom_point(aes(order, ES, color= DIVERS))+
  geom_errorbar(aes(order,ES,
                    ymin=ES-sqrt(vi)/2*1.96,
                    ymax=ES+sqrt(vi)/2*1.96), width=.2,
                position=position_dodge(.9))+ theme_pubr()+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray", size=0.4)+
  facet_wrap(~Sub_Cat_Output)


# Structure the data frame by nesting
DATA <- RATIO %>%
  ungroup()       %>%
  mutate(DIVERS = tolower(DIVERS)) %>% 
  filter(!DIVERS %in% c("others","2 and more")) %>% 
  group_by(Sub_Cat_Output, DIVERS)             %>% 
  nest() 
# create new columns
DATA$ES<-0; DATA$lower<-0; DATA$higher<-0
DATA$lowerpred<-0; DATA$higherpred<-0
DATA$ESlm<-0; DATA$lowerlm<-0; DATA$higherlm<-0
DATA$NB_ES<-NA; DATA$NB_MA<-NA; DATA$NB_Data<-NA
DATA$EStf<-0; DATA$lowertf<-0; DATA$highertf<-0
DATA$ROSEN<-0


NB_ES<- data.frame(matrix(ncol=2, nrow=36))
names(NB_ES)<- c("indice","nb")

for (i in 1: 35){
  NB_ES[i,"indice"]<- i
  NB_ES[i,"nb"]<- dim(DATA[[3]][[i]])[1]
}
NB_ESsup1<- NB_ES   %>% filter(nb > 1)

for (i in NB_ESsup1$indice){
  print(i)
TAB<-DATA[[3]][[i]]  # select DIVERS i and Outcome i
TAB$CODE<-1:nrow(TAB)
TAB$Scenario_ID<-1:nrow(TAB)

## Quality folowing XXX et al., 2018
TAB$Qi   <- TAB$Score/max(TAB$Score)
TAB$taui <-(((1-TAB$Qi)/TAB$vi)/(nrow(TAB)-1))

for(j in 1: length(TAB$Score)){
  TAB$Qprimj[j]<- 
  if(TAB$Qi[j]<1){
  (sum(TAB$Qi)*TAB$taui[j])/
    (sum(TAB$taui)*(nrow(TAB)-1))+
    TAB$Qi} else{
 TAB$Qi[j]}
}

for(j in 1: length(TAB$Score)){
    TAB$tauprim[j]<-
      (sum(TAB$taui)*nrow(TAB)*
         (TAB$Qprimj[j]/ (sum(TAB$Qprimj))) -
         TAB$taui[j])
  }
  
TAB$Qv<-TAB$Qi/TAB$vi
TAB$W<-TAB$Qv+TAB$tauprim


### Redundancy  
# create a dataframe with all combination
COMBI <- data.frame(MA1=rep(tolower(TAB$ID), length(TAB$ID)),
                   MA2= rep(tolower(TAB$ID), each=length(TAB$ID)))
rownames(MAT)<-tolower(rownames(MAT))
colnames(MAT)<-tolower(colnames(MAT))

#MAT<-arrange(MAT, ID)
TT <- MAT %>% filter(ma %in% tolower(unique(TAB$ID)))
TT <- TT %>% select( tolower(unique(TAB$ID)))
TT<-subset(TT, select=order(colnames(TT)))

NB<-data.frame(MA=colnames(TT), 
               NB=diag(data.matrix(TT)))

TT     <- data.frame(data.matrix(TT)* 2)
TT$ID  <- colnames(TT)
MAT_spread <- gather(TT, "MA2", "RED", -ID)
names(MAT_spread)[1]<- "MA1"

COMBI$Order<-1: nrow(COMBI)
EE <- merge(COMBI,MAT_spread,by=c("MA1", "MA2"), keep=TRUE)
EE <- merge(EE, NB, by.x = c("MA1"), by.y = c("MA"))
EE <- merge(EE, NB, by.x = c("MA2"), by.y = c("MA"))
EE$SUM <- EE$NB.x+EE$NB.y
EE$REDONDANCY <- EE$RED/EE$SUM
EE<- EE  %>% arrange(Order)
EE<- EE %>% select(-NB.x, -NB.y, -RED, -SUM, -Order)

for (j in 1: length(EE$MA1)){
  if(EE$MA1[j]==EE$MA2[j]){EE$REDONDANCY[j]<-0}
}

LONG <- length(TAB$ID)
d <- split(EE,rep(1:LONG,each=LONG))
TABLE<-data.frame(matrix(ncol=LONG, nrow=LONG))
for (j in 1: LONG){
  TABLE[j,]<- d[[j]]$RED
}

#matrice des variances multipliées entre elles
HH<-data.matrix(TAB$vi)
MM<-sqrt(HH %*% t(HH))

# on multiplie les 2 matrices
NN<-MM*TABLE
#on ajoute les variances pour la diagonale
diag(NN)<- TAB$vi
#diag(NN)<- 1/TAB$W
NN<-data.frame(NN)


# Meta-analytical model:
mod<-rma.mv(ES, vi,data=TAB,random=~1|ID/NUM, method="REML") 

# Publication-bias
plot1<-viz_funnel(x=mod, 
           contours_col = "Greys",
           trim_and_fill = TRUE, trim_and_fill_side = "right", 
           egger = TRUE)
name<- paste("viz_funnel", DATA[[1]][[i]], DATA[[2]][[i]], ".pdf")

ggsave(name, width = 30, height = 15, units = "cm")

# Publication-bias second plot
name<- paste("viz_sunset", DATA[[1]][[i]], DATA[[2]][[i]], ".pdf")
viz_sunset(x=mod, true_effect=mod$b,
           power_contours = "continuous")
ggsave(name, width = 30, height = 15, units = "cm")

# forest plot with all conserved effect-sizes
ggplot(TAB, aes(x=ES, y=reorder(Scenario_ID, ES)))+
  geom_point()+
  geom_rect(aes(xmin=mod$b-sqrt(mod$beta[1])/2*1.96,
                xmax = mod$b+sqrt(mod$beta[1])/2*1.96,
                ymin = 1,
                ymax = length(TAB$ID)),
                fill = 'gray80', alpha = 0.05)+
  geom_point()+
  geom_errorbar(aes(xmin = ES-sqrt(vi)/2*1.96, xmax = ES+sqrt(vi)/2*1.96))+
  geom_vline(aes(xintercept = mod$b))+
  geom_vline(aes(xintercept = 0), linetype=2)+
  theme_pubr()

name<- paste("forest_plot", DATA[[1]][[i]], DATA[[2]][[i]], ".pdf")
ggsave(name, width = 30, height = 15, units = "cm")

summary(mod)

# meta-analytical model 2
mod<-rma.mv(ES, NN,data=TAB,random=~1|ID/NUM, method="REML") 
summary(mod)

modlm<-rma.uni(ES, vi,data=TAB, method="REML") 

if(length(TAB$ID)<3){
modtm<-modlm
} else {
 modtm<-trimfill(modlm)
}

ROSEN<-fsn(yi= TAB$ES, vi=TAB$vi)

# based on mixed effect model
RES             <- summary(mod)
DATA$ES[i]      <- RES$b[1,1]
DATA$lower[i]   <- RES$ci.lb[1]
DATA$higher[i]  <- RES$ci.ub[1]

# based on mixed effect model: prediction intervals
RES             <- predict.rma(mod)
DATA$lowerpred[i]   <- RES$cr.lb[1]
DATA$higherpred[i]  <- RES$cr.ub[1]

mod<-rma.mv(ES, NN,data=TAB,random=~1|ID/NUM, method="REML") 
predict.rma(mod)

# based on linear model
RES             <- summary(modlm)
DATA$ESlm[i]      <- RES$b[1,1]
DATA$lowerlm[i]   <- RES$ci.lb[1]
DATA$higherlm[i]  <- RES$ci.ub[1]

# based on trim and fill
RES              <- summary(modtm)
DATA$EStf[i]      <- RES$b[1,1]
DATA$lowertf[i]   <- RES$ci.lb[1]
DATA$highertf[i]  <- RES$ci.ub[1]

# ROSENfail safe number 
DATA$ROSEN[i] <- ROSEN$fsnum

DATA$NB_ES[i]   <- nrow(TAB)
DATA$NB_MA[i]   <- length(unique(TAB$ID))
DATA$NB_Data[i] <- sum(TAB$Nb_Paired_Data, na.rm=TRUE)
}

# For combination of DIVERS and OUTcomes with only 1, we do not use a model, but directly the estimated effect-sizes
NB_ES1<- NB_ES   %>% filter(nb == 1)

for (i in NB_ES1$indice){
TAB<-DATA[[3]][[i]]
DATA$ES[i]<- TAB$ES
DATA$lower[i]<- TAB$ES-sqrt(TAB$vi)/2*1.96
DATA$higher[i]<- TAB$ES+sqrt(TAB$vi)/2*1.96
DATA$NB_ES[i]<- nrow(TAB)
DATA$NB_MA[i]<- length(unique(TAB$ID))
DATA$NB_Data[i]<- sum(TAB$Nb_Paired_Data, na.rm=TRUE)


}

# PAR output

DATA2 <- RATIO %>%
  mutate(CODE= paste(Sub_Cat_Output))  %>% 
  mutate(DIVERS = tolower(DIVERS)) %>% 
  filter(!DIVERS %in% c("others","2 and more")) %>% 
  ungroup()                            %>%
  group_by(Sub_Cat_Output)             %>% 
  nest() 

DATA2$ES<-0; DATA2$lower<-0; DATA2$higher<-0
DATA2$lowerpred<-0; DATA2$higherpred<-0
DATA2$ESlm<-0; DATA2$lowerlm<-0; DATA2$higherlm<-0
DATA2$NB_ES<-NA; DATA2$NB_MA<-NA; DATA2$NB_Data<-NA
DATA2$EStf<-0; DATA2$lowertf<-0; DATA2$highertf<-0
DATA2$ROSEN<-0


for (i in c(1:13)){
  TAB<-DATA2[[2]][[i]]
  
  mod<-rma.mv(ES, vi, data=TAB,random=~1|ID/NUM, method="REML") 
  RES<-summary(mod)
  DATA2$ES[i]<- RES$b[1,1]
  DATA2$lower[i]<- RES$ci.lb[1]
  DATA2$higher[i]<- RES$ci.ub[1]
  
  
  RES<-predict.rma(mod)
  DATA2$lowerpred[i]<- RES$cr.lb[1]
  DATA2$higherpred[i]<- RES$cr.ub[1]
  
  DATA2$NB_ES[i]<- nrow(TAB)
  DATA2$NB_MA[i]<- length(unique(TAB$ID))
  DATA2$NB_Data[i]<- sum(TAB$Nb_Paired_Data, na.rm=TRUE)
  
  modlm<-rma.uni(ES, vi,data=TAB, method="REML") 
  
  if(length(TAB$ID)<3){
    modtm<-modlm
  } else {
    modtm<-trimfill(modlm)
  }
  
  ROSEN<-fsn(yi= TAB$ES, vi=TAB$vi)
  
  # based on linear model
  RES             <- summary(modlm)
  DATA2$ESlm[i]      <- RES$b[1,1]
  DATA2$lowerlm[i]   <- RES$ci.lb[1]
  DATA2$higherlm[i]  <- RES$ci.ub[1]
  
  # based on trim and fill
  RES              <- summary(modtm)
  DATA2$EStf[i]      <- RES$b[1,1]
  DATA2$lowertf[i]      <- RES$ci.lb[1]
  DATA2$highertf[i]  <- RES$ci.ub[1]
  
  # ROSENfail safe number 
  DATA2$ROSEN[i] <- ROSEN$fsnum
  
  
}

for (i in c(14)){
  TAB<-DATA2[[2]][[i]]
  DATA2$ES[i]<- TAB$ES
  DATA2$lower[i]<- TAB$ES-sqrt(TAB$vi)/2*1.96
  DATA2$higher[i]<- TAB$ES+sqrt(TAB$vi)/2*1.96
  DATA2$NB_ES[i]<- nrow(TAB)
  DATA2$NB_MA[i]<- length(unique(TAB$ID))
  DATA2$NB_Data[i]<- sum(TAB$Nb_Paired_Data, na.rm=TRUE)
}


DATA<- DATA %>% dplyr::select(-data)
DATA2<- DATA2 %>% dplyr::select(-data)
DATA2$DIVERS<-"ALL"
DATA_TOTAL<-merge(DATA,DATA2,all=TRUE)

DATA_TOTAL$DIVERS<-fct_relevel(DATA_TOTAL$DIVERS, "ALL")
#DATA_TOTAL<- DATA_TOTAL %>% filter(!DIVERS =="2 AND MORE")
#DATA_TOTAL$DIVERS<-factor(DATA_TOTAL$DIVERS)

DATA_TOTAL$DIVERS <- relevel(DATA_TOTAL$DIVERS, "ALL")
#DATA_TOTAL<- DATA_TOTAL  %>% filter(!Sub_Cat_Output %in% c("yield stability",
 #                                                          "yield and quality" ,
  #                                                         "root biomass",
   #                                                        "profitability"))

# levels(DATA_TOTAL$DIVERS)[2] <- "AgroF"
# levels(DATA_TOTAL$DIVERS)[3] <- "Assoc. Pl."
# levels(DATA_TOTAL$DIVERS)[4] <- "Intercr."
DATA_TOTAL$DIVERS<- factor(DATA_TOTAL$DIVERS)


d <- data.frame(Sub_Cat_Output= unique(DATA_TOTAL$Sub_Cat_Output),
  image = c("~/Documents/analyse_diversification/images_plot/BIODIV2.jpg",
            "~/Documents/analyse_diversification/images_plot/GHG.jpg",
            "~/Documents/analyse_diversification/images_plot/Nuse efficiency2.jpg",
            "~/Documents/analyse_diversification/images_plot/others.png",
            "~/Documents/analyse_diversification/images_plot/pest2.png",
            "~/Documents/analyse_diversification/images_plot/Production4.png",
            "~/Documents/analyse_diversification/images_plot/Quality.png",
            "~/Documents/analyse_diversification/images_plot/euro.jpg",
            "~/Documents/analyse_diversification/images_plot/root.jpg",
            "~/Documents/analyse_diversification/images_plot/soil_quality.jpg",
            "~/Documents/analyse_diversification/images_plot/Water_Quality.jpg",
            "~/Documents/analyse_diversification/images_plot/Water.png",
            "~/Documents/analyse_diversification/images_plot/Quality.png",
            "~/Documents/analyse_diversification/images_plot/satbility.jpg"))

DATA_TOTAL2<- merge(DATA_TOTAL, d, by= "Sub_Cat_Output")
DATA_TOTAL2$Code<- row.names(DATA_TOTAL2) 
DATA_TOTAL$DIVERS2<- "ALL"
DATA_TOTAL2$DIVERS2<- "ALL"
library(tidyr, warn.conflicts = FALSE)
DATA_TOTAL2 %<>% 
  tidyr::complete(ES, nesting(DIVERS2, Sub_Cat_Output), fill = list(value1 = 0))
detach("package:tidyr", unload = TRUE)


newdata <- DATA2$Sub_Cat_Output[order(-DATA2$ES)]
newdata2<- data.frame(Sub_Cat_Output=c("water quality","water use",
             "pests and diseases regulation","biodiv. assoc." ,
             "inputs use efficiency"  ,  "Production",   "yield and quality" ,
             "root biomass" , "products quality","yield stability" ,
              "soil quality" , "other" ,            "greenhouses gas emission",
               "profitability"), 
             NAME=c("EAU", "EAU", "BIOD", "BIOD", "PROD","PROD","PROD","PROD","PROD","PROD","SOIL_GH","SOIL_GH","SOIL_GH","ECO"))
DATA_TOTAL2<-merge(DATA_TOTAL2, newdata2, by= "Sub_Cat_Output")

DATA_TOTAL2$Sub_Cat_Output <- ordered(DATA_TOTAL2$Sub_Cat_Output,
                                      levels = newdata2$Sub_Cat_Output)


library(ggimage)

DELTA<- 1.5
DATA_TOTAL2$ES<-DATA_TOTAL2$ES+DELTA
DATA_TOTAL2$lower<-DATA_TOTAL2$lower+DELTA
DATA_TOTAL2$higher<-DATA_TOTAL2$higher+DELTA
DATA_TOTAL2$lowerpred<-DATA_TOTAL2$lowerpred+DELTA
DATA_TOTAL2$higherpred<-DATA_TOTAL2$higherpred+DELTA

ggplot(data=DATA_TOTAL2 %>% filter(DIVERS=="ALL"), aes(group=DIVERS))+ 
ylim(c(min(DATA_TOTAL2$higher),max(DATA_TOTAL2$higher)))+
  geom_hline(yintercept=c(0, 0.69, 1.38)+DELTA, linetype=2, color="red")  +
  geom_bar(aes(x=interaction(Sub_Cat_Output, DIVERS2), y= higherpred), fill= "gray90", stat="identity")+
  geom_rect( mapping=aes(xmin=1, xmax=10, ymin=-1+DELTA, ymax=0+DELTA), fill="black", color="black", alpha=0.5) +
  geom_bar(aes(x=interaction(Sub_Cat_Output, DIVERS2), y= higher), fill= "gray80", stat="identity")+
  geom_bar(aes(x=interaction(Sub_Cat_Output, DIVERS2), y= ES+0.005),fill="white", stat="identity")+
  geom_bar(aes(x=interaction(Sub_Cat_Output, DIVERS2), y= ES-0.005),fill="gray80", stat="identity")+
  geom_bar(aes(x=interaction(Sub_Cat_Output, DIVERS2), y= lower), fill="gray90",size=3, stat="identity")+
  geom_bar(aes(x=interaction(Sub_Cat_Output, DIVERS2), y= lowerpred), fill="white",size=3, stat="identity")+
  theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_vline(aes(xintercept =c(1:10)-0.5), linetype=2, color=c("black","gray50","black","gray50","black","gray50","gray50","black","gray50","black"))+
 geom_image(aes(size=0.1, image=image,
                x=interaction(Sub_Cat_Output, DIVERS2),
                y=max(higherpred)+0.2,
               label=Sub_Cat_Output,
               angle=90))+
              scale_size_identity()+
 #annotation_custom(rasterGrob(tiger), xmin = 5, xmax = 6, ymin = 0.2, ymax = 1)+
 coord_polar()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_errorbar(data= DATA_TOTAL2 %>% filter(!DIVERS=="ALL"),
                aes(x=interaction(Sub_Cat_Output, DIVERS2), ymin=lower, ymax=higher,size=DIVERS, color=DIVERS), 
                width=.001,size=1.1, position=position_dodge(width=0.6)) +
  geom_errorbar(data= DATA_TOTAL2 %>% filter(!DIVERS=="ALL"),
                aes(x=interaction(Sub_Cat_Output, DIVERS2), ymin=lowerpred, ymax=higherpred,size=DIVERS, color=DIVERS), 
                width=.001,size=0.6,alpha=0.5, position=position_dodge(width=0.6)) +
  geom_point(data= DATA_TOTAL2 %>% filter(!DIVERS=="ALL"),
             aes(x=interaction(Sub_Cat_Output, DIVERS2), y= ES, fill=DIVERS),shape=21, color="black",
             position=position_dodge(width=0.6),size=2)+
  labs(title="Impacts of diversification strategies", x= "Studied Output", y= "Effect size")+
  labs(x="",y= "Ratio diversified systems/ less diversified")+ 
  #scale_colour_manual(values = c("#339965", "#3299fe", "#fe9833", "#983200", "#65cc98", "#3266cb", "red",  "red"))
#scale_colour_manual(values = c("#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e", "red",  "red"))+
#scale_fill_manual(values = c("#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e", "red",  "red"))+
  scale_colour_manual(values = c("#810f7c", "orange", "#3690c0", "#969696", "red", "red", "red",  "red"))+
  scale_fill_manual(values = c("#810f7c", "orange", "#3690c0", "#969696", "red", "red", "red",  "red"))


DATA_TOTAL2<- DATA_TOTAL2 %>% filter(Sub_Cat_Output %in% c("biodiv. assoc.", "greenhouses gas emission", 
                                                  "inputs use efficiency","pests and diseases regulation",
                                                  "Production"    ,"products quality"  ,"profitability" ,
                                                  "soil quality" ,"water quality"  ,"water use" ))


DATA_TOTAL3<-DATA_TOTAL2 %>% filter(!is.na(DIVERS))
ggplot(data=DATA_TOTAL3, aes(group=DIVERS))+ 
  ylim(c(min(DATA_TOTAL2$higher),max(DATA_TOTAL2$higher)))+
  scale_size_identity()+
  #annotation_custom(rasterGrob(tiger), xmin = 5, xmax = 6, ymin = 0.2, ymax = 1)+
  #coord_polar()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  geom_errorbar(
                aes(x=interaction(Sub_Cat_Output, DIVERS2), ymin=lower, ymax=higher,size=DIVERS, color=DIVERS), 
                width=.001,size=1.1, position=position_dodge(width=0.6)) +
  geom_errorbar(
                aes(x=interaction(Sub_Cat_Output, DIVERS2), ymin=lowerpred, ymax=higherpred,size=DIVERS, color=DIVERS), 
                width=.001,size=0.6,alpha=0.5, position=position_dodge(width=0.6)) +
  geom_point(
             aes(x=interaction(Sub_Cat_Output, DIVERS2), y= ES, fill=DIVERS),shape=21, color="black",
             position=position_dodge(width=0.6),size=2)+
  labs(title="Impacts of diversification strategies", x= "Studied Output", y= "Effect size")+
  labs(x="",y= "Ratio diversified systems/ less diversified")+ 
  geom_text(aes(x= interaction(Sub_Cat_Output, DIVERS2), color=DIVERS ,y= 0.3,
                label=paste0(NB_MA,"-",NB_Data)),
            angle=90,position=position_dodge(width=0.6))+
  scale_fill_manual(values = c("black", "orange", "red", "blue", "purple", "gray70", "red",  "red"))+
  scale_color_manual(values = c("black", "orange", "red", "blue", "purple", "gray70", "red",  "red"))+
  theme_pubr()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  geom_hline(yintercept=c(0)+DELTA, linetype=2, color="black")+
  geom_hline(yintercept=c(0.69, 1.38)+DELTA, linetype=2, color="gray50")


  




TT$higherpred[TT$higherpred == 0] <- NA
TT$lowerpred[TT$lowerpred == 0] <- NA

ggplot(data=TT %>% filter(DIVERS=="ALL"), aes(group=DIVERS))+ 
  ylim(c(min(DATA_TOTAL2$higher),max(DATA_TOTAL2$higher)))+
  geom_bar(aes(x=interaction(Sub_Cat_Output, DIVERS2), y= higherpred), fill= "gray90", stat="identity")+
  geom_bar(aes(x=interaction(Sub_Cat_Output, DIVERS2), y= higher), fill= "gray80", stat="identity")+
  geom_bar(aes(x=interaction(Sub_Cat_Output, DIVERS2), y= ES+0.005),fill="white", stat="identity")+
  geom_bar(aes(x=interaction(Sub_Cat_Output, DIVERS2), y= ES-0.005),fill="gray80", stat="identity")+
  geom_bar(aes(x=interaction(Sub_Cat_Output, DIVERS2), y= lower), fill="gray90",size=3, stat="identity")+
  geom_bar(aes(x=interaction(Sub_Cat_Output, DIVERS2), y= lowerpred), fill="white",size=3, stat="identity")+
  geom_hline(yintercept=c(0)+DELTA, linetype=2, color="red", alpha=0.3)+
  geom_hline(yintercept=c(0, 0.69, 1.38)+DELTA, linetype=2, color="blue", alpha=0.3)+
  theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_vline(aes(xintercept =c(1:10)-0.5), linetype=2, color="gray90")+
  geom_image(aes(size=0.1, image=image,
                 x=interaction(Sub_Cat_Output, DIVERS2),
                 y=4,
                 label=Sub_Cat_Output,
                 angle=90))+
  scale_size_identity()+
  scale_y_continuous(  breaks = c(0+DELTA,0.70+DELTA,1.38+DELTA),
                       labels = c(1,2,4))+
  coord_polar() +
  theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())+
  geom_text(aes(x= interaction(Sub_Cat_Output, DIVERS2),y= 3.5,label=paste(NB_Data)))


### plot général
### sensitivity analysis
# DELTA<- -1
# DATA_TOTAL2$ES<-DATA_TOTAL2$ES+DELTA
# DATA_TOTAL2$lower<-DATA_TOTAL2$lower+DELTA
# DATA_TOTAL2$higher<-DATA_TOTAL2$higher+DELTA


TAB1<-DATA_TOTAL2 %>% dplyr::select( "DIVERS2" ,"Sub_Cat_Output", "DIVERS" ,
                                    "lower" , "ES", "higher" , 
                                    "NB_ES" , "NB_MA" , "NB_Data", "image")
TAB1$method<- "RE"
TAB2<-DATA_TOTAL2 %>% dplyr::select( "DIVERS2" ,"Sub_Cat_Output", "DIVERS" ,
                                     "lowerlm" , "ESlm", "higherlm" , 
                                     "NB_ES" , "NB_MA" , "NB_Data", "image")
TAB2$method<- "FE"
names(TAB2)[c(4,5,6)]<- c("lower", "ES", "higher")

TAB3<-DATA_TOTAL2 %>% dplyr::select( "DIVERS2" ,"Sub_Cat_Output", "DIVERS" ,
                                     "lowertf" , "EStf", "highertf" , 
                                     "NB_ES" , "NB_MA" , "NB_Data", "image")
TAB3$method<- "TF"
names(TAB3)[c(4,5,6)]<- c("lower", "ES", "higher")

TAB_AS<-rbind(TAB1, TAB2, TAB3)
for (i in 1: length(TAB_AS$DIVERS2)){
if(is.na(TAB_AS$lower[i])){
  TAB_AS$ES[i]<-"NA"
}
}
TAB_AS$ES<- as.numeric(as.character(TAB_AS$ES))
ggplot(data=TAB_AS %>% filter(!DIVERS=="ALL"), aes(group=interaction(method)))+ 
geom_errorbar(aes(x=interaction(method,DIVERS),
                  ymin=lower, ymax=higher,size=DIVERS, color=method),
              width=.001,size=0.6, position=position_dodge(width=0.6)) +
  geom_point(aes(x=interaction(method,DIVERS),
                 y= ES, fill=DIVERS),shape=21, color="black",
             position=position_dodge(width=0.6),size=2)+
  facet_grid(~Sub_Cat_Output)+
  labs(title="Impacts of diversification strategies", x= "Studied Output", y= "Effect size")+
  labs(x="",y= "Ratio diversified systems/ less diversified")+ 
  #scale_colour_manual(values = c("#339965", "#3299fe", "#fe9833", "#983200", "#65cc98", "#3266cb", "red",  "red"))
  #scale_colour_manual(values = c("#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e", "red",  "red"))+
  #scale_fill_manual(values = c("#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e", "red",  "red"))+
  scale_fill_manual(values = c("#810f7c", "orange", "#3690c0", "#969696", "red", "red", "red",  "red"))+
  scale_colour_manual(values = c("black", "gray90", "gray80", "#969696", "red", "red", "red",  "red"))+
  scale_size_identity()+
  #annotation_custom(rasterGrob(tiger), xmin = 5, xmax = 6, ymin = 0.2, ymax = 1)+
  #  coord_polar()+
  theme_pubr()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_hline(yintercept=c(0, 0.69, 1.38), linetype=2, color="red")
  

ggplot(data=DATA_TOTAL, aes(color=DIVERS, group=DIVERS))+ 
  geom_point(aes(x=interaction(Sub_Cat_Output, DIVERS,method), y= ES),position=position_dodge(width=0.6),size=3)+
  geom_errorbar(aes(x=interaction(Sub_Cat_Output, DIVERS, method), ymin=lower, ymax=higher,size=DIVERS), 
                width=.1,size=0.9, position=position_dodge(width=0.6)) +
  geom_hline(yintercept=0, linetype=2)

+
 # facet_wrap(~Sub_Cat_Output, drop=TRUE,scales="free")+
  # scale_x_discrete(breaks=(DATA$DIVERS), 
  #                  labels=addline_format(DATA$DIVERS))+
  theme_bw(base_size=12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0 ),
        panel.grid.major = element_blank())+
  labs(title="Impacts of diversification strategies", x= "Studied Output", y= "Effect size")+
  labs(x="",y= "Ratio diversified systems/ less diversified")+ 
  scale_colour_manual(values = c("black", "#fe6766", "#993200", "#68cb9b", "#ff9934", "#339866", "#3299fe", "#3266cb",  "red"))+
  scale_size_manual(values = c(1,rep(0.5,7)))+
  geom_text(aes(x=DIVERS,
                group=DIVERS, y=ES, 
                label=paste("\n","\n" ,NB_Data)), 
           angle=90, size=4,position = position_dodge(width = 1), inherit.aes = FALSE)+
   theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=14))+
  coord_polar()


ggsave("Forest_plot_DETAIL.pdf", width = 20, height = 20, units = "cm")



data <- metaDigitise(dir = "~/Documents/Meta_synthesis_GLOBAL/essai")

