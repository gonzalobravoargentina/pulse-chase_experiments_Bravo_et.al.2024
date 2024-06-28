#NATURAl ISOTOPES VALUES---------
#all station together
DATA <- "DATA"
allSTATIONS <- read.csv(file.path(DATA,"ST435-407-177-Pulse chase experiments.csv"))
allSTATIONS <- subset(allSTATIONS,GROUP!="Nematodes")
naturalisotopes <- subset(allSTATIONS,TREATMENT=="C")
naturalisotopes <- naturalisotopes[!is.na(naturalisotopes$δ15N_air),]
naturalisotopes <- subset(naturalisotopes,GROUP!="Echinodermata")
naturalisotopes <- subset(naturalisotopes,SPECIES!="Dacrydium sp.")
naturalisotopes$FAMILY <- paste(naturalisotopes$FAMILY,naturalisotopes$FEED.MODE,sep=" ")
naturalisotopes <- naturalisotopes[order(naturalisotopes$STATION),]
table(naturalisotopes$STATION)
naturalisotopes$delta15Nsediment <- c(rep(7.885000,26),rep(7.093000,34),rep(6.172333,26))
naturalisotopes$T.L <-   (naturalisotopes$δ15N_air - naturalisotopes$delta15Nsediment)/3.8+1
naturalisotopes$C.N <- naturalisotopes$C_µg /naturalisotopes$N_µg
naturalisotopesMACROFAUNA <- subset(naturalisotopes,GROUP!="Bulk sediments")
#naturalisotopesMACROFAUNA <- subset(naturalisotopes,GROUP="Polychaeta")



##------------------------------------------------------------------------------------------------------------------------------
#MOST FAMILIES AND SEDIMENT--------
library(ggplot2)
naturalisotopes.allFAMILY <- subset(naturalisotopes,GROUP=="Polychaeta"|GROUP=="Bulk sediments"|GROUP=="Bivalve"|GROUP=="Nemertea"|GROUP=="Cnidaria")
library(data.table)
naturalisotopesSTATIONS.FAMILY <- setDT(naturalisotopes)[,list(n=length(Formalin.corrected_δ13C_VPDB),x=mean(Formalin.corrected_δ13C_VPDB), ESx=sd(Formalin.corrected_δ13C_VPDB),RangeX=range(Formalin.corrected_δ13C_VPDB) ,y=mean(δ15N_air),ESy=sd(δ15N_air),Rangey=range(δ15N_air),C.N=mean(C.N),C.N.SD=sd(C.N),T.L=mean(T.L),T.L.SD=sd(T.L)), by=c("STATION","FAMILY")]

#ordination by 15N values low to high
position14c <- c("Bulk sediments ","Diastylidae FF/SDF","Ampeliscidae FF","Onuphidae P/S","Ampharetidae SDF","Polynoidae P/S","Alcyoniidae FF","Nuculidae  FF/SDF","Cirratulidae FF/SDF","Spionidae FF/SDF","Lumbrineridae P/S","Pholoidae  P/S","Thyasiridae FF/SDF","Nemertea P/S","Arcidae FF","Pseudotanaidae SDF","Philomedidae SDF","Yoldiidae FF/SDF","Sipuncula SDF","Priapulidae SDF","Opheliidae SSDF","Chaetodermatidae P/S","Maldanidae SSDF","Syllidae P/S","Nephtyidae  P/S")

#position14c <- c("Bulk sediments ","Ampharetidae SDF","","Pholoidae  P/S","Onuphidae P/S","Polynoidae P/S","Cirratulidae FF/SDF","Spionidae FF/SDF","Syllidae P/S","Opheliidae SSDF","Lumbrineridae P/S","Maldanidae SSDF","Nephtyidae  P/S","Nemertea P/S","Yoldiidae FF/SDF","Thyasiridae FF/SDF","Arcidae FF", "Alcyoniidae FF")


#both 13C and 15C by family
pd <- position_dodge(0.7)
ThesisPalette <- c('gray48','gray76',"gray21")
C13 <- ggplot(naturalisotopesSTATIONS.FAMILY, aes(x=FAMILY, y=x,group=STATION))+
  geom_errorbar(aes(ymin=x-ESx,ymax=x+ESx),width=0.1,na.rm = TRUE,position = pd,alpha = 1/3)+
  geom_point(position=pd,aes(shape=factor(STATION),colour=factor(STATION)),size=3.5,na.rm = TRUE)+  coord_flip()+theme_bw(base_size = 12,base_family = "Helvetica")+
  theme(legend.background = element_rect(fill="transparent",linetype="solid", colour ="black"),axis.text.x = element_text(),axis.title.y=element_blank(),axis.text.y = element_blank(),legend.title = element_blank(),legend.text=element_text(size=10),legend.position=c(0.73,0.15),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+scale_colour_manual(values=ThesisPalette,breaks=c("ST177", "ST407", "ST435"),labels=c("Baffin Bay (Stn 177)", "Amundsen Gulf (Stn 407)", "Beaufort Sea (Stn 435)"))+scale_shape_discrete(breaks=c("ST177", "ST407", "ST435"),labels=c("Baffin Bay (Stn 177)", "Amundsen Gulf (Stn 407)", "Beaufort Sea (Stn 435)"))+
  scale_y_continuous(breaks =c(-27,-26,-25,-24,-23,-22,-21,-20,-19,-18))+
  labs(y=expression({delta}^13*C~'\u2030'),factor="")+scale_x_discrete(limits = position14c,position = "top")

N15<-  ggplot(naturalisotopesSTATIONS.FAMILY, aes(x=FAMILY, y=y,group=STATION))+
  geom_errorbar(aes(ymin=y-ESy,ymax=y+ESy),width=0.1,na.rm = TRUE,position = pd,alpha = 1/3)+
  geom_point(position=pd,aes(shape=factor(STATION),colour=factor(STATION)),size=3.5,na.rm = TRUE,show.legend = F)+
  theme_bw(base_size = 12,base_family = "Helvetica")+
  theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_text(),axis.title.y=element_blank(),axis.text.y = element_text(hjust=0.5),legend.title = element_blank(),legend.text=element_text(size=10),legend.position=c(0.92,0.3),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
  coord_flip() +  
  scale_colour_manual(values=ThesisPalette)+
  scale_y_continuous(limits = c(4, 20),breaks = c(6,8,10,12,14,16,18,20))+
  labs(y=expression({delta}^15*N~'\u2030'))+scale_x_discrete(limits = position14c)

C.N<- ggplot(naturalisotopesSTATIONS.FAMILY, aes(x=FAMILY, y=C.N,group=STATION))+
  geom_errorbar(aes(ymin=C.N-C.N.SD,ymax=C.N+C.N.SD),width=0.1,na.rm = TRUE,position = pd,alpha = 1/3)+
  geom_point(position=pd,aes(shape=factor(STATION),colour=factor(STATION)),size=3.5,na.rm = TRUE,show.legend = F)+
  theme_bw(base_size = 12,base_family = "Helvetica")+
  theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_text(),axis.title.y=element_blank(),axis.text.y = element_text(hjust=0.5),legend.title = element_blank(),legend.text=element_text(size=10),legend.position=c(0.92,0.3),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
  coord_flip() +  
  scale_colour_manual(values=ThesisPalette)+  
  scale_x_discrete(limits = position14c)+ geom_hline(yintercept =c(7.03,6.96,7.56), linetype="dashed",colour=c('gray48','gray76',"gray21"))

library(cowplot)
C.N.Ratios <- plot_grid(C13,N15,C.N,ncol = 3 ,nrow=1, align = "v",labels=c("A","B","C"),hjust=-11.3,vjust = 1.9)
ggsave(filename = "NATURAL.ISOTOPES.png",plot = C.N.Ratios,width=15,height = 8)

library(cowplot)
plot_grid(C13,N15,ncol = 2 ,nrow=1, align = "v")

##------------------------------------------------------------------------------------------------------------------------------
#BOX PLOT FEED MODE by 15N---------
pd <- position_dodge(0.5)
positionsfeedmode <- c("FF/SDF","SDF","SSDF","P/S")
ThesisPalette <- c('gray48','gray76',"gray21")
naturalisotopesMACROFAUNA$STATION.FEEDMODE <- paste(naturalisotopesMACROFAUNA$FEED.MODE,naturalisotopesMACROFAUNA$STATION)

N15_feedmode <- ggplot(naturalisotopesMACROFAUNA, aes(x=FEED.MODE, y=δ15N_air,group=STATION.FEEDMODE))+
  geom_boxplot(aes(fill=STATION,colour=STATION),size=0.5,na.rm = TRUE)+theme_bw(base_size = 12,base_family = "Helvetica")+theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_text(),axis.title.y=element_text(),axis.text.y = element_text(hjust=0.5),legend.title = element_blank(),legend.text=element_text(size=10), panel.grid.minor = element_blank())+ labs(y= expression("15N")) +scale_x_discrete(limits = positionsfeedmode)+  scale_y_continuous(limits = c(6,20), breaks = c(8,10,12,14,16,18))+ labs(y= expression(paste("δ"^"15","N ‰")),x="Feeding groups")+scale_fill_manual(values=ThesisPalette,breaks=c("ST177", "ST407", "ST435"),labels=c("Baffin Bay (Stn 177)", "Amundsen Gulf (Stn 407)", "Beaufort Sea (Stn 435)"))+scale_colour_manual(values=c("black","black","black"),breaks=c("ST177", "ST407", "ST435"),labels=c("Baffin Bay (Stn 177)", "Amundsen Gulf (Stn 407)", "Beaufort Sea (Stn 435)"))

+geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1), dotsize=0.3)

ggsave(filename = "N15_feedmode.png",plot = N15_feedmode,width=10,height = 5)

library(data.table)
naturalisotopesMACROFAUNA2 <- setDT(naturalisotopesMACROFAUNA)[,list(n=length(Formalin.corrected_δ13C_VPDB),x=mean(Formalin.corrected_δ13C_VPDB), ESx=sd(Formalin.corrected_δ13C_VPDB),RangeX=range(Formalin.corrected_δ13C_VPDB),y=mean(δ15N_air),ESy=sd(δ15N_air),Rangey=range(δ15N_air)),by=c("STATION","FEED.MODE")]

N15_feedmode <- ggplot(naturalisotopesMACROFAUNA2, aes(x=FEED.MODE, y=y,group=STATION))+
  geom_errorbar(aes(ymin=y-ESy,ymax=y+ESy),width=0.1,na.rm = TRUE,position = pd,alpha = 1/3)+
  geom_point(position=pd,aes(shape=factor(STATION),colour=factor(STATION)),size=3.5,na.rm = TRUE)+theme_bw(base_size = 12,base_family = "Helvetica")+
  theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_text(),axis.title.y=element_text(),axis.text.y = element_text(hjust=0.5),legend.title = element_blank(),legend.text=element_text(size=10), panel.grid.minor = element_blank())+ labs(y= expression(paste("δ"^"15","N")))+
  scale_colour_manual(values=ThesisPalette)+scale_x_discrete(limits = positionsfeedmode)+  scale_y_continuous(limits = c(9,19), breaks = c(8,10,12,14,16,18))




positionsfeedmode <- c("FF/SDF","SDF","SSDF","P/S")
boxplot(δ15N_air~FEED.MODE*STATION, data=naturalisotopesMACROFAUNA,position = pd, 
        col=(c('gray48','gray76',"gray21")))
##------------------------------------------------------------------------------------------------------------------------------



### delta organism - delta sediment 
naturalisotopesMACROFAUNA[1:22,"delta"] <- naturalisotopesMACROFAUNA[1:22,24] - -22.2
naturalisotopesMACROFAUNA[23:53,"delta"] <- naturalisotopesMACROFAUNA[23:53,24]- -23.8
naturalisotopesMACROFAUNA[53:76,"delta"] <- naturalisotopesMACROFAUNA[53:76,24] - -25.20000

## T.L------
library(data.table)
naturalisotopesSTATIONS.FAMILY <- setDT(naturalisotopesMACROFAUNA)[,list(n=length(Formalin.corrected_δ13C_VPDB),x=mean(Formalin.corrected_δ13C_VPDB), ESx=sd(Formalin.corrected_δ13C_VPDB),RangeX=range(Formalin.corrected_δ13C_VPDB) ,y=mean(δ15N_air),ESy=sd(δ15N_air),Rangey=range(δ15N_air),C.N=mean(C.N),C.N.SD=sd(C.N),T.L=mean(T.L),T.L.SD=sd(T.L),TLrange=range(T.L)), by=c("STATION")]
pd <- position_dodge(0.5)

ggplot(naturalisotopesSTATIONS.FAMILY, aes(x=FAMILY, y=T.L,group=STATION))+
  geom_errorbar(aes(ymin=T.L-T.L.SD,ymax=T.L+T.L.SD),width=0.1,na.rm = TRUE,position = pd,alpha = 1/3)+
  geom_point(position=pd,aes(shape=factor(STATION),colour=factor(STATION)),size=3.5,na.rm = TRUE,show.legend = F)+theme_bw(base_size = 12,base_family = "Helvetica")+
  theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_text(),axis.title.y=element_blank(),axis.text.y = element_text(hjust=0.5),legend.title = element_blank(),legend.text=element_text(size=10),legend.position=c(0.92,0.3),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
  coord_flip() +  
  scale_colour_manual(values=ThesisPalette)+  
  scale_y_continuous(limits = c(0,5))+scale_x_discrete(limits = position14c)

TLbystation <- summaryBy(T.L ~  STATION, data =naturalisotopesSTATIONS.FAMILY, FUN = function(x) { c(m = mean(x), SD = sd(x),Range=range(x)) })
names(naturalisotopesSTATIONS.FAMILY)

#T.L by FEED MODE--------
naturalisotopes.allFAMILY <- subset(naturalisotopes,GROUP=="Polychaeta"|GROUP=="Bulk sediments"|GROUP=="Bivalve"|GROUP=="Nemertea")
library(data.table)
naturalisotopesFEEDMODE <- setDT(naturalisotopes.allFAMILY)[,list(n=length(Formalin.corrected_δ13C_VPDB),x=mean(Formalin.corrected_δ13C_VPDB), ESx=sd(Formalin.corrected_δ13C_VPDB),RangeX=range(Formalin.corrected_δ13C_VPDB) ,y=mean(δ15N_air),ESy=sd(δ15N_air),Rangey=range(δ15N_air),C.N=mean(C.N),C.N.SD=sd(C.N),T.L=mean(T.L),T.L.SD=sd(T.L)), by=c("STATION","FEED.MODE")]

ThesisPalette <- c('gray48','gray76',"gray21")
positionsfeedmode <- c("FF/SDF","SDF","SSDF","P/S")
pd <- position_dodge(0.5)
feedmodeTL <- ggplot(naturalisotopesFEEDMODE, aes(x=FEED.MODE, y=T.L,group=STATION))+
  geom_errorbar(aes(ymin=T.L-T.L.SD,ymax=T.L+T.L.SD),width=0.1,na.rm = TRUE,position = pd,alpha = 1/3)+
  geom_point(position=pd,aes(shape=factor(STATION),colour=factor(STATION)),size=3.5,na.rm = TRUE)+theme_bw(base_size = 12,base_family = "Helvetica")+
  theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_text(),axis.title.y=element_text(),axis.text.y = element_text(hjust=0.5),legend.title = element_blank(),legend.text=element_text(size=10),legend.position=c(0.06,0.1), panel.grid.minor = element_blank())+
  scale_colour_manual(values=ThesisPalette)+scale_x_discrete(limits = positionsfeedmode)+  scale_y_continuous(limits = c(-0.5,5), breaks = c(0,1,2,3,4,5))

ggsave(filename = "feedmode.png",plot = feedmodeTL,width=15,height = 5)

naturalisotopesMACROFAUNA$STATION.FEEDMODE <- paste(naturalisotopesMACROFAUNA$FEED.MODE,naturalisotopesMACROFAUNA$STATION)
boxplot.feedmodeTL <- ggplot(naturalisotopesMACROFAUNA, aes(x=FEED.MODE, y=T.L,group=STATION.FEEDMODE))+
  geom_boxplot(aes(fill=STATION,colour=STATION),size=0.5,na.rm = TRUE)+theme_bw(base_size = 12,base_family = "Helvetica")+theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_text(),axis.title.y=element_text(),axis.text.y = element_text(hjust=0.5),legend.title = element_blank(),legend.text=element_text(size=10), panel.grid.minor = element_blank())+ labs(y= expression("15N")) +scale_x_discrete(limits = positionsfeedmode) + labs(y= expression(paste("δ"^"15","N")))+scale_fill_manual(values=ThesisPalette,breaks=c("ST177", "ST407", "ST435"),labels=c("Stn 177", "Stn 407", "Stn 435"))+scale_colour_manual(values=c("black","black","black"),breaks=c("ST177", "ST407", "ST435"),labels=c("Stn 177", "Stn 407", "Stn 435"))



#mean by layer in each station 
library(doBy)
meanbylayerinstation <- summaryBy(Formalin.corrected_δ13C_VPDB + δ15N_air ~  STATION + LAYER , data =naturalisotopesMACROFAUNA ,FUN = mean)
names(naturalisotopesMACROFAUNA)
abundance.biomasCORES <- summaryBy(COUNTS + ABUNDANCE + BIOMASS + Biomass.mg.C.m2 ~ CORE + FAMILY + STATION  + LAYER , data =subset(allSTATIONS,ABUNDANCE>1|BIOMASS>1) ,FUN = sum,na.rm = TRUE)
tapply(naturalisotopesMACROFAUNA$δ15N_air,naturalisotopesMACROFAUNA$LAYER,mean,na.rm=T)

isotopebylayer <- summaryBy(Formalin.corrected_δ13C_VPDB + δ15N_air ~ GROUP + STATION  + LAYER , data =naturalisotopesMACROFAUNA ,FUN = mean,na.rm = TRUE)

#biplot stations allmacrofauna
library(data.table)
naturalisotopesSTATIONS <- setDT(naturalisotopesMACROFAUNA)[,list(n=length(Formalin.corrected_δ13C_VPDB),x=mean(Formalin.corrected_δ13C_VPDB), ESx=sd(Formalin.corrected_δ13C_VPDB),RangeX=range(Formalin.corrected_δ13C_VPDB) ,y=mean(δ15N_air),ESy=sd(δ15N_air),Rangey=range(δ15N_air)), by=STATION]

#biplot of 14C, 15N isotopes allstations-------------
library(cowplot)#It looks similar to theme_classic(), with a fewdifferences. Most importantly, the default label sizes are larger and the background is transparent.
#colours
ThesisPalette <- c('gray48','gray76',"gray21")

library(ggplot2)
ggplot(naturalisotopesSTATIONS, aes(x=x, y=y))+
  geom_errorbarh(aes(xmin=x-ESx,xmax=x+ESx),height=0.1,na.rm = TRUE)+
  geom_errorbar(aes(ymin=y-ESy,ymax=y+ESy),width=0.1,na.rm = TRUE) + 
  theme_cowplot()+
  theme(legend.position=c(0.85,0.9),legend.title=element_blank(),legend.text=element_text(size=15))+
  geom_point(aes(shape = factor(STATION),colour = factor(STATION)),size=5)+
  scale_shape_manual(values=c(16,16,16))+
  labs(x=expression({delta}^13*C~'\u2030'),y=expression({delta}^15*N~'\u2030'))+
  scale_y_continuous(limits = c(10, 18))+
  scale_x_continuous(limits = c(-27, -20))+
  scale_colour_manual(values=ThesisPalette)

#bivalves, polychaetas and Bulk sediments
naturalisotopesSTATIONS.GROUP <- setDT(naturalisotopes)[,list(n=length(Formalin.corrected_δ13C_VPDB),x=mean(Formalin.corrected_δ13C_VPDB), ESx=sd(Formalin.corrected_δ13C_VPDB),RangeX=range(Formalin.corrected_δ13C_VPDB) ,y=mean(δ15N_air),ESy=sd(δ15N_air),Rangey=range(δ15N_air)), by=c("STATION","GROUP")]
naturalisotopesSTATIONS.GROUP <- subset(naturalisotopesSTATIONS.GROUP, GROUP=="Bulk sediments"|GROUP=="Polychaeta"|GROUP=="Bivalve")
position14c <- c("Bulk sediments","Polychaeta","Bivalve")
ThesisPalette <- c('gray48','gray76',"gray21")
ggplot(naturalisotopesSTATIONS.GROUP, aes(x=x, y=y))+
  geom_errorbarh(aes(xmin=x-ESx,xmax=x+ESx),height=0.1,na.rm = TRUE)+
  geom_errorbar(aes(ymin=y-ESy,ymax=y+ESy),width=0.1,na.rm = TRUE) + 
  theme_cowplot()+
  theme(legend.position=c(0.9,0.3),legend.title=element_blank(),legend.text=element_text(size=15))+
  geom_point(aes(shape = factor(GROUP),colour = factor(STATION)),size=5)+
  labs(x=expression({delta}^13*C~'\u2030'),y=expression({delta}^15*N~'\u2030'))+
  scale_y_continuous(limits = c(4, 19),breaks = c(6,8,10,12,14,16,18))+
  scale_x_continuous(limits = c(-25, -19),breaks = c(-25,-24,-23,-22,-21))+
  scale_colour_manual(values=ThesisPalette)















