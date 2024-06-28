#TOTAL UPTAKE in all station----------------
DATA <- "DATA"
data_STATIONS <- read.csv(file.path(DATA,"ST435-407-177-Pulse chase experiments.csv"))
data_STATIONS <- subset(data_STATIONS,GROUP!="Nematodes")
uptake_data_STATIONS <- subset(data_STATIONS, ug.C.mg.1 > 0 & TREATMENT!="C" )
uptake_data_STATIONS <- subset(data_STATIONS, ug.N.mg.1 > 0)
#uptake_data_STATIONS<- subset(uptake_data_STATIONS, SPECIES!="Amage sp.")
#TOTALuptake
library(doBy)
TOTALuptakeBYstations <-  summaryBy(I.total..ug.C.m2.+I.total..ug.N.m2.~ CORE + TREATMENT +STATION , data = uptake_data_STATIONS, FUN = sum,na.rm = TRUE)
TOTALuptakeBYstations_mean <-  summaryBy(I.total..ug.C.m2..sum + I.total..ug.N.m2..sum ~ TREATMENT + STATION , data = TOTALuptakeBYstations, FUN = function(x) { c(m = mean(x), ES = sqrt(var(x)/length(x))) })

TOTALuptakeBYstations_mean.notreatment <-  summaryBy(I.total..ug.C.m2..sum + I.total..ug.N.m2..sum ~ STATION , data = TOTALuptakeBYstations, FUN = function(x) { c(m = mean(x), ES = sqrt(var(x)/length(x))) })

library(ggplot2)
positionST <- c("ST435","ST407","ST177")
C_TOTALuptakeBYstations<- ggplot(TOTALuptakeBYstations_mean, aes( x=STATION, y=I.total..ug.C.m2..sum.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=I.total..ug.C.m2..sum.m-I.total..ug.C.m2..sum.ES, ymax=I.total..ug.C.m2..sum.m+I.total..ug.C.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Total uptake (",mu,"g C m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.text.y = element_text(size=10),legend.text=element_text(size=8),legend.position=c(0.85,0.95),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=3100, label= "A",size=4) +scale_x_discrete(limits = positionST)+ annotate("text", x=3, y=2500, label= "*",size=6)


positionST <- c("ST435","ST407","ST177")
N_TOTALuptakeBYstations<- ggplot(TOTALuptakeBYstations_mean, aes( x=STATION, y=I.total..ug.N.m2..sum.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=I.total..ug.N.m2..sum.m-I.total..ug.N.m2..sum.ES, ymax=I.total..ug.N.m2..sum.m+I.total..ug.N.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Total uptake (",mu,"g N m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.text.y = element_text(size=10),legend.text=element_text(size=8),legend.position=c(0.85,0.95),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=1000, label= "B",size=4) +scale_x_discrete(limits = positionST)+ annotate("text", x=3, y=500, label= "*",size=6)


#BIOMASS-SPECIFIC uptake in all stations------
library(doBy)
biomassspecificuptakeBYstations <-  summaryBy(ug.C.mg.1+ ug.N.mg.1~ CORE + TREATMENT +STATION , data = uptake_data_STATIONS, FUN = mean,na.rm = TRUE)
biomassspecificuptakeBYstations_mean <-  summaryBy(ug.C.mg.1.mean+ ug.N.mg.1.mean ~ TREATMENT + STATION , data = biomassspecificuptakeBYstations, FUN = function(x) { c(m = mean(x), ES = sqrt(var(x)/length(x))) })

#t-tests
variable <- subset(biomassspecificuptakeBYstations,STATION=="ST177")
variable$ug.C.mg.1.mean.log <- log(variable$ug.C.mg.1.mean+1)
shapiro.test(variable$ug.C.mg.1.mean)
t.test(ug.C.mg.1.mean~ TREATMENT, data=variable, var.equal = FALSE) #var.equal = FALSE for Welch approximation
wilcox.test(ug.C.mg.1.mean ~ TREATMENT, data=variable)



C_biomassspecificuptakeBYstations<- ggplot(biomassspecificuptakeBYstations_mean, aes( x=STATION, y=ug.C.mg.1.mean.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=ug.C.mg.1.mean.m-ug.C.mg.1.mean.ES, ymax=ug.C.mg.1.mean.m+ug.C.mg.1.mean.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Biomass-specific uptake(",mu,"g C mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(axis.title.x=element_blank(),axis.text.x=element_text(),axis.text.y = element_text(size=10),legend.position=c(0.8,0.9),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton")) +guides(fill=FALSE)+guides(fill=FALSE)+ annotate("text", x=0.6, y=8, label= "C",size=4)+ scale_x_discrete(breaks=c("ST435","ST407","ST177"), labels=c("Beaufort Sea (Stn 177)", "Amundsen Gulf (Stn 407)", "Baffin Bay (Stn 435)"),limits = positionST)+scale_y_continuous(limits = c(0, 8))+ annotate("text", x=3, y=3, label= "**",size=6)

N_biomassspecificuptakeBYstations<- ggplot(biomassspecificuptakeBYstations_mean, aes( x=STATION, y=ug.N.mg.1.mean.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=ug.N.mg.1.mean.m-ug.N.mg.1.mean.ES, ymax=ug.N.mg.1.mean.m+ug.N.mg.1.mean.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Biomass-specific uptake (",mu,"g N mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(axis.title.x=element_blank(),axis.text.y = element_text(size=10),legend.position=c(0.78,0.9),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton")) +guides(fill=FALSE)+ annotate("text", x=0.6, y=8, label= "D",size=4) + scale_x_discrete(breaks=c("ST435","ST407","ST177"), labels=c("Beaufort Sea (Stn 177)", "Amundsen Gulf (Stn 407)", "Baffin Bay (Stn 435)"),limits = positionST)+scale_y_continuous(limits = c(0, 8))+ annotate("text", x=3, y=2, label= "**",size=6)

library(cowplot)
totalANDspecifUPTAKEperSTATION <- plot_grid(C_TOTALuptakeBYstations,N_TOTALuptakeBYstations,C_biomassspecificuptakeBYstations,N_biomassspecificuptakeBYstations,ncol = 2, align = "v")
#ggsave(filename = "totalANDspecifUPTAKEperSTATION.png",plot = totalANDspecifUPTAKEperSTATION,width=12,height = 10)

###########################################################################################################

#t-test TOTAL UPTAKES-------
#ST177-ugCm2
ST177.TOTALuptake <-subset(TOTALuptakeBYstations,STATION=="ST177")
t.test(I.total..ug.C.m2..sum~ TREATMENT, data=ST177.TOTALuptake)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(ST177.TOTALuptake$I.total..ug.C.m2..sum)
#ST177-ugNm2
ST177.TOTALuptake <-subset(TOTALuptakeBYstations,STATION=="ST177")
t.test(I.total..ug.N.m2..sum~ TREATMENT, data=ST177.TOTALuptake)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(ST177.TOTALuptake$I.total..ug.N.m2..sum)

#ST407-ugCm2
ST407.TOTALuptake <-subset(TOTALuptakeBYstations,STATION=="ST407")
shapiro.test(ST407.TOTALuptake$I.total..ug.C.m2..sum)
t.test(log(I.total..ug.C.m2..sum)~ TREATMENT, data=ST407.TOTALuptake)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(log(ST407.TOTALuptake$I.total..ug.C.m2..sum))

#ST407-ugNm2
ST407.TOTALuptake <-subset(TOTALuptakeBYstations,STATION=="ST407")
t.test(log(I.total..ug.N.m2..sum)~ TREATMENT, data=ST407.TOTALuptake)
ST407.TOTALuptake$I.total..ug.N.m2..sum.log <- log(ST407.TOTALuptake$I.total..ug.N.m2..sum)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(ST407.TOTALuptake$I.total..ug.N.m2..sum.log)

#ST435-ugCm2
ST435.TOTALuptake <-subset(TOTALuptakeBYstations,STATION=="ST435")
t.test(I.total..ug.C.m2..sum~ TREATMENT, data=ST435.TOTALuptake)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(ST435.TOTALuptake$I.total..ug.C.m2..sum)
#ST435-ugNm2
ST435.TOTALuptake <-subset(TOTALuptakeBYstations,STATION=="ST435")
t.test(I.total..ug.N.m2..sum~ TREATMENT, data=ST435.TOTALuptake)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(ST435.TOTALuptake$I.total..ug.N.m2..sum)

#t-test BIOMASS_SPECIFIC UPTAKES
#ST177 ug C mg
ST177.BIOMASS.SPECIFICuptake <-subset(biomassspecificuptakeBYstations,STATION=="ST177")
t.test(ug.C.mg.1.mean~ TREATMENT, data=ST177.BIOMASS.SPECIFICuptake)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(ST177.BIOMASS.SPECIFICuptake$ug.C.mg.1.mean)

#ST177 ug N mg
ST177.BIOMASS.SPECIFICuptake <-subset(biomassspecificuptakeBYstations,STATION=="ST177")
t.test(ug.N.mg.1.mean~ TREATMENT, data=ST177.BIOMASS.SPECIFICuptake)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(ST177.BIOMASS.SPECIFICuptake$ug.N.mg.1.mean)

#ST407 ug C mg
ST407.BIOMASS.SPECIFICuptake <-subset(biomassspecificuptakeBYstations,STATION=="ST407")
t.test(log(ug.C.mg.1.mean)~ TREATMENT, data=ST407.BIOMASS.SPECIFICuptake)
ST407.BIOMASS.SPECIFICuptake$ug.C.mg.1.mean.log <- log(ST407.BIOMASS.SPECIFICuptake$ug.C.mg.1.mean)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(ST407.BIOMASS.SPECIFICuptake$ug.C.mg.1.mean.log)

#ST407 ug N mg
ST407.BIOMASS.SPECIFICuptake <-subset(biomassspecificuptakeBYstations,STATION=="ST407")
t.test(log(ug.N.mg.1.mean)~ TREATMENT, data=ST407.BIOMASS.SPECIFICuptake)
ST407.BIOMASS.SPECIFICuptake$ug.N.mg.1.mean.log <- log(ST407.BIOMASS.SPECIFICuptake$ug.N.mg.1.mean)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(ST407.BIOMASS.SPECIFICuptake$ug.N.mg.1.mean.log)
wilcox.test(ug.N.mg.1.mean~ TREATMENT, data=ST407.BIOMASS.SPECIFICuptake)

#ST435 ug C mg
ST435.BIOMASS.SPECIFICuptake <-subset(biomassspecificuptakeBYstations,STATION=="ST435")
t.test(ug.C.mg.1.mean~ TREATMENT, data=ST435.BIOMASS.SPECIFICuptake)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(ST435.BIOMASS.SPECIFICuptake$ug.C.mg.1.mean)

#ST435 ug N mg
ST435.BIOMASS.SPECIFICuptake <-subset(biomassspecificuptakeBYstations,STATION=="ST435")
t.test(ug.N.mg.1.mean~ TREATMENT, data=ST435.BIOMASS.SPECIFICuptake)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(ST435.BIOMASS.SPECIFICuptake$ug.N.mg.1.mean)


###########################################################################################################
# TOTAL AND BIOMASS-SPECIFIC uptake by IMPORTANT family in STATIONS -------------------------------------
library(doBy)
#ST177
BS.T.byfamilyST177 <- subset(uptake_data_STATIONS,STATION=="ST177")
BS.T.byfamilyST177 <- subset(BS.T.byfamilyST177,FAMILY=="Yoldiidae"|FAMILY=="Thyasiridae"|FAMILY=="Spionidae"|FAMILY=="Lumbrineridae"|FAMILY=="Paraonidae "|FAMILY=="Nephtyidae "|FAMILY=="Maldanidae")
BS.byfamilyST177 <-  summaryBy(ug.C.mg.1 +ug.N.mg.1 ~ CORE + FAMILY + TREATMENT, data = BS.T.byfamilyST177, FUN = mean,na.rm = TRUE)
BS.byfamily_meanST177 <-  summaryBy(ug.C.mg.1.mean+ ug.N.mg.1.mean ~ FAMILY+TREATMENT , data = BS.byfamilyST177, FUN = function(x) { c(m = mean(x),n=length(x), ES = sqrt(var(x)/length(x))) })


T.byfamilyST177 <-  summaryBy(I.total..ug.C.m2.+I.total..ug.N.m2.~ CORE  + FAMILY+ TREATMENT , data = BS.T.byfamilyST177, FUN = sum,na.rm = TRUE)
T.byfamilyST177 <-  summaryBy(I.total..ug.C.m2..sum + I.total..ug.N.m2..sum ~  FAMILY+TREATMENT , data = T.byfamilyST177, FUN = function(x) { c(m = mean(x), ES = sqrt(var(x)/length(x))) })


positionST177 <- c("Spionidae","Lumbrineridae","Paraonidae ","Nephtyidae ","Maldanidae","Yoldiidae","Thyasiridae")
library(ggplot2)
T.ST177.C<- ggplot(T.byfamilyST177, aes( x=FAMILY, y=I.total..ug.C.m2..sum.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=I.total..ug.C.m2..sum.m-I.total..ug.C.m2..sum.ES, ymax=I.total..ug.C.m2..sum.m+I.total..ug.C.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(title="Baffin Bay Stn",y= expression(paste("Total uptake (ug C m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.text.y = element_text(size=12),legend.text=element_text(size=12),legend.position=c(0.85,0.91),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=2500, label= "A",size=4) +scale_x_discrete(limits = positionST177)

BS.ST177.C <- ggplot(BS.byfamily_meanST177, aes(y=ug.C.mg.1.mean.m, x=FAMILY,fill=TREATMENT)) + geom_col(color="black",position = position_dodge(1),na.rm = TRUE)+geom_errorbar(aes(ymin=ug.C.mg.1.mean.m-ug.C.mg.1.mean.ES, ymax=ug.C.mg.1.mean.m+ug.C.mg.1.mean.ES), width=.2,position=position_dodge(.9))+labs(y= expression(paste("Biomass-specific uptake (ug C mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(axis.title.x=element_blank(),axis.text.x = element_text(size=15,angle = 45,hjust=1),axis.text.y = element_text(size=12),legend.position=c(0.78,0.9),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton")) + annotate("text", x=0.6, y=6.9, label= "C",size=4) + scale_x_discrete(limits = positionST177)+scale_y_continuous(limits = c(-0.17, 7))+ guides(fill=F)+
  annotate("text", x=0.75, y=-0.15, label= "n= 5",size=3)+ annotate("text", x=1.25, y=-0.15, label= "n= 5",size=3)+annotate("text", x=1.75, y=-0.15, label= "n= 4",size=3)+ annotate("text", x=2.25, y=-0.15, label= "n= 4",size=3)+annotate("text", x=2.75, y=-0.15, label= "n= 3",size=3)+ annotate("text", x=3.25, y=-0.15, label= "n= 1",size=3)+annotate("text", x=3.95, y=-0.15, label= "n= 1",size=3)+ annotate("text", x=5, y=-0.15, label= "n= 2",size=3)+
  annotate("text", x=5.75, y=-0.15, label= "n= 5",size=3)+ annotate("text", x=6.25, y=-0.15, label= "n= 5",size=3)+
  annotate("text", x=6.75, y=-0.15, label= "n= 5",size=3)+ annotate("text", x=7.25, y=-0.15, label= "n= 5",size=3)


T.ST177.N<- ggplot(T.byfamilyST177, aes( x=FAMILY, y=I.total..ug.N.m2..sum.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=I.total..ug.N.m2..sum.m-I.total..ug.N.m2..sum.ES, ymax=I.total..ug.N.m2..sum.m+I.total..ug.N.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Total uptake (ug N m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.text.y = element_text(size=12),legend.text=element_text(size=8),legend.position=c(0.85,0.95),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=500, label= "B",size=4) +scale_x_discrete(limits = positionST177)+guides(fill=F)

BS.ST177.N <- ggplot(BS.byfamily_meanST177, aes(y=ug.N.mg.1.mean.m, x=FAMILY,fill=TREATMENT)) + geom_col(color="black",position = position_dodge(1),na.rm = TRUE)+geom_errorbar(aes(ymin=ug.N.mg.1.mean.m-ug.N.mg.1.mean.ES, ymax=ug.N.mg.1.mean.m+ug.N.mg.1.mean.ES), width=.2,position=position_dodge(.9))+labs(y= expression(paste("Biomass-specific uptake (ug N mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(axis.title.x=element_blank(),axis.text.x = element_text(size=15,angle = 45,hjust=1),axis.text.y = element_text(size=12),legend.position=c(0.78,0.9),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton")) + annotate("text", x=0.6, y=6.9, label= "D",size=4) + scale_x_discrete(limits = positionST177)+scale_y_continuous(limits = c(-0.17, 7))+ guides(fill=F)+
  annotate("text", x=0.75, y=-0.15, label= "n= 5",size=3)+ annotate("text", x=1.25, y=-0.15, label= "n= 5",size=3)+annotate("text", x=1.75, y=-0.15, label= "n= 4",size=3)+ annotate("text", x=2.25, y=-0.15, label= "n= 4",size=3)+annotate("text", x=2.75, y=-0.15, label= "n= 3",size=3)+ annotate("text", x=3.25, y=-0.15, label= "n= 1",size=3)+annotate("text", x=3.95, y=-0.15, label= "n= 1",size=3)+ annotate("text", x=5, y=-0.15, label= "n= 2",size=3)+
  annotate("text", x=5.75, y=-0.15, label= "n= 5",size=3)+ annotate("text", x=6.25, y=-0.15, label= "n= 5",size=3)+
  annotate("text", x=6.75, y=-0.15, label= "n= 5",size=3)+ annotate("text", x=7.25, y=-0.15, label= "n= 5",size=3)+ annotate("text", x=6, y=5, label= "*",size=6)

library(cowplot)
TandBS.ST177 <- plot_grid(T.ST177.C,T.ST177.N,BS.ST177.C,BS.ST177.N,ncol = 2, align = "v")
ggsave(filename = "TandBS.ST177.png",plot = TandBS.ST177,width=20,height = 15)

positionST177 <- c("Spionidae","Lumbrineridae","Paraonidae ","Nephtyidae ","Maldanidae","Yoldiidae","Thyasiridae")

#stats
BS.T.byfamilyST177_forTtest <- subset(BS.T.byfamilyST177, FAMILY=="Yoldiidae")

#ug.C.mg.1
t.test(ug.C.mg.1 ~ TREATMENT, data=BS.T.byfamilyST177_forTtest)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(BS.T.byfamilyST177_forTtest$ug.C.mg.1)
shapiro.test(log(BS.T.byfamilyST177_forTtest$ug.C.mg.1))
wilcox.test(ug.C.mg.1 ~ TREATMENT, data=BS.T.byfamilyST177_forTtest)

#ug.N.mg.1
t.test(ug.N.mg.1 ~ TREATMENT, data=BS.T.byfamilyST177_forTtest)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(BS.T.byfamilyST177_forTtest$ug.N.mg.1)
shapiro.test(log(BS.T.byfamilyST177_forTtest$ug.N.mg.1))
wilcox.test(ug.N.mg.1 ~ TREATMENT, data=BS.T.byfamilyST177_forTtest)

#I.total..ug.C.m2.
t.test(I.total..ug.C.m2. ~ TREATMENT, data=BS.T.byfamilyST177_forTtest)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(BS.T.byfamilyST177_forTtest$I.total..ug.C.m2.)
shapiro.test(log(BS.T.byfamilyST177_forTtest$I.total..ug.C.m2.))
wilcox.test(I.total..ug.C.m2. ~ TREATMENT, data=BS.T.byfamilyST177_forTtest)

#I.total..ug.N.m2.
t.test(I.total..ug.N.m2. ~ TREATMENT, data=BS.T.byfamilyST177_forTtest)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(BS.T.byfamilyST177_forTtest$I.total..ug.N.m2.)
shapiro.test(log(BS.T.byfamilyST177_forTtest$I.total..ug.N.m2.))
wilcox.test(I.total..ug.N.m2. ~ TREATMENT, data=BS.T.byfamilyST177_forTtest)


#ST407
BS.T.byfamilyST407 <- subset(uptake_data_STATIONS,STATION=="ST407")
BS.T.byfamilyST407 <- subset(BS.T.byfamilyST407,FAMILY=="Yoldiidae"|FAMILY=="Ampharetidae"|FAMILY=="Spionidae"|FAMILY=="Maldanidae")
BS.byfamilyST407 <-  summaryBy(ug.C.mg.1 +ug.N.mg.1 ~ CORE + FAMILY + TREATMENT, data = BS.T.byfamilyST407, FUN = mean,na.rm = TRUE)
BS.byfamily_meanST407 <-  summaryBy(ug.C.mg.1.mean+ ug.N.mg.1.mean ~ FAMILY+TREATMENT , data = BS.byfamilyST407, FUN = function(x) { c(m = mean(x),n=length(x), ES = sqrt(var(x)/length(x))) })

T.byfamilyST407 <-  summaryBy(I.total..ug.C.m2.+I.total..ug.N.m2.~ CORE  + FAMILY+ TREATMENT , data = BS.T.byfamilyST407, FUN = sum,na.rm = TRUE)
T.byfamilyST407 <-  summaryBy(I.total..ug.C.m2..sum + I.total..ug.N.m2..sum ~  FAMILY+TREATMENT , data = T.byfamilyST407, FUN = function(x) { c(m = mean(x), ES = sqrt(var(x)/length(x))) })


positionST407 <- c("Spionidae","Maldanidae","Ampharetidae","Yoldiidae")
library(ggplot2)
T.ST407.C<- ggplot(T.byfamilyST407, aes( x=FAMILY, y=I.total..ug.C.m2..sum.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=I.total..ug.C.m2..sum.m-I.total..ug.C.m2..sum.ES, ymax=I.total..ug.C.m2..sum.m+I.total..ug.C.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(title="Amundsen Gulf Stn",y= expression(paste("Total uptake (ug C m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.text.y = element_text(size=12),legend.text=element_text(size=12),legend.position=c(0.85,0.91),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=4200, label= "A",size=4) +scale_x_discrete(limits = positionST407)

BS.ST407.C <- ggplot(BS.byfamily_meanST407, aes(y=ug.C.mg.1.mean.m, x=FAMILY,fill=TREATMENT)) + geom_col(color="black",position = position_dodge(1),na.rm = TRUE)+geom_errorbar(aes(ymin=ug.C.mg.1.mean.m-ug.C.mg.1.mean.ES, ymax=ug.C.mg.1.mean.m+ug.C.mg.1.mean.ES), width=.2,position=position_dodge(.9))+labs(y= expression(paste("Biomass-specific uptake (ug C mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(axis.title.x=element_blank(),axis.text.x = element_text(size=15,angle = 45,hjust=1),axis.text.y = element_text(size=12),legend.position=c(0.78,0.9),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton")) + annotate("text", x=0.6, y=35, label= "C",size=4) + scale_x_discrete(limits = positionST407)+scale_y_continuous(limits = c(-1,35))+ guides(fill=F)+
  annotate("text", x=0.75, y=-0.8, label= "n= 2",size=3)+ annotate("text", x=1.25, y=-0.8, label= "n= 2",size=3)+
  annotate("text", x=1.75, y=-0.8, label= "n= 5",size=3)+ annotate("text", x=2.25, y=-0.8, label= "n= 2",size=3)+
  annotate("text", x=2.75, y=-0.8, label= "n= 3",size=3)+ annotate("text", x=3.25, y=-0.8, label= "n= 3",size=3)+
  annotate("text", x=3.74, y=-0.8, label= "n= 2",size=3)+ annotate("text", x=4.25, y=-0.8, label= "n= 1",size=3)


T.ST407.N<- ggplot(T.byfamilyST407, aes( x=FAMILY, y=I.total..ug.N.m2..sum.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=I.total..ug.N.m2..sum.m-I.total..ug.N.m2..sum.ES, ymax=I.total..ug.N.m2..sum.m+I.total..ug.N.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Total uptake (ug N m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.text.y = element_text(size=12),legend.text=element_text(size=8),legend.position=c(0.85,0.95),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=900, label= "B",size=4) +scale_x_discrete(limits = positionST407)+guides(fill=F)

BS.ST407.N <- ggplot(BS.byfamily_meanST407, aes(y=ug.N.mg.1.mean.m, x=FAMILY,fill=TREATMENT)) + geom_col(color="black",position = position_dodge(1),na.rm = TRUE)+geom_errorbar(aes(ymin=ug.N.mg.1.mean.m-ug.N.mg.1.mean.ES, ymax=ug.N.mg.1.mean.m+ug.N.mg.1.mean.ES), width=.2,position=position_dodge(.9))+labs(y= expression(paste("Biomass-specific uptake (ug N mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(axis.title.x=element_blank(),axis.text.x = element_text(size=15,angle = 45,hjust=1),axis.text.y = element_text(size=12),legend.position=c(0.78,0.9),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton")) + annotate("text", x=0.6, y=35, label= "D",size=4) + scale_x_discrete(limits = positionST407)+scale_y_continuous(limits = c(-1, 35))+ guides(fill=F)+
  annotate("text", x=0.75, y=-0.8, label= "n= 2",size=3)+ annotate("text", x=1.25, y=-0.8, label= "n= 2",size=3)+
  annotate("text", x=1.75, y=-0.8, label= "n= 5",size=3)+ annotate("text", x=2.25, y=-0.8, label= "n= 2",size=3)+
  annotate("text", x=2.75, y=-0.8, label= "n= 3",size=3)+ annotate("text", x=3.25, y=-0.8, label= "n= 3",size=3)+
  annotate("text", x=3.74, y=-0.8, label= "n= 2",size=3)+ annotate("text", x=4.25, y=-0.8, label= "n= 1",size=3) 

library(cowplot)
TandBS.ST407 <- plot_grid(T.ST407.C,T.ST407.N,BS.ST407.C,BS.ST407.N,ncol = 2, align = "v")
ggsave(filename = "TandBS.ST407.png",plot = TandBS.ST407,width=20,height = 15)

#stats
BS.T.byfamilyST407_forTtest <- subset(BS.T.byfamilyST407, FAMILY=="Maldanidae")

#ug.C.mg.1
t.test(ug.C.mg.1 ~ TREATMENT, data=BS.T.byfamilyST407_forTtest)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(BS.T.byfamilyST407_forTtest$ug.C.mg.1)
shapiro.test(log(BS.T.byfamilyST407_forTtest$ug.C.mg.1))
wilcox.test(ug.C.mg.1 ~ TREATMENT, data=BS.T.byfamilyST407_forTtest)

#I.total..ug.C.m2.
t.test(I.total..ug..m2. ~ TREATMENT, data=BS.T.byfamilyST407_forTtest)
#normally distributed "the samples come from a Normal distribution"
shapiro.test(BS.T.byfamilyST407_forTtest$I.total..ug.C.m2.)
shapiro.test(log(BS.T.byfamilyST407_forTtest$I.total..ug.C.m2.))
wilcox.test(I.total..ug.C.m2. ~ TREATMENT, data=BS.T.byfamilyST407_forTtest)




#ST435
BS.T.byfamilyST435 <- subset(uptake_data_STATIONS,STATION=="ST435")
BS.T.byfamilyST435 <- subset(BS.T.byfamilyST435,FAMILY=="Spionidae"|FAMILY=="Mytilidae "|FAMILY=="Yoldiidae"|FAMILY=="Alcyoniidae")
BS.byfamilyST435 <-  summaryBy(ug.C.mg.1 +ug.N.mg.1 ~ CORE + FAMILY + TREATMENT, data = BS.T.byfamilyST435, FUN = mean,na.rm = TRUE)
BS.byfamily_meanST435 <-  summaryBy(ug.C.mg.1.mean+ ug.N.mg.1.mean ~ FAMILY+TREATMENT , data = BS.byfamilyST435, FUN = function(x) { c(m = mean(x),n=length(x), ES = sqrt(var(x)/length(x))) })

T.byfamilyST435 <-  summaryBy(I.total..ug.C.m2.+I.total..ug.N.m2.~ CORE  + FAMILY+ TREATMENT , data = BS.T.byfamilyST435, FUN = sum,na.rm = TRUE)
T.byfamilyST435 <-  summaryBy(I.total..ug.C.m2..sum + I.total..ug.N.m2..sum ~  FAMILY+TREATMENT , data = T.byfamilyST435, FUN = function(x) { c(m = mean(x), ES = sqrt(var(x)/length(x))) })


positionST435 <- c("Spionidae","Mytilidae ","Yoldiidae","Alcyoniidae")
library(ggplot2)
T.ST435.C<- ggplot(T.byfamilyST435, aes( x=FAMILY, y=I.total..ug.C.m2..sum.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=I.total..ug.C.m2..sum.m-I.total..ug.C.m2..sum.ES, ymax=I.total..ug.C.m2..sum.m+I.total..ug.C.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(title="Beaufort Sea Stn",y= expression(paste("Total uptake (ug C m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.text.y = element_text(size=12),legend.text=element_text(size=12),legend.position=c(0.85,0.91),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=1000, label= "A",size=4) +scale_x_discrete(limits = positionST435)

BS.ST435.C <- ggplot(BS.byfamily_meanST435, aes(y=ug.C.mg.1.mean.m, x=FAMILY,fill=TREATMENT)) + geom_col(color="black",position = position_dodge(1),na.rm = TRUE)+geom_errorbar(aes(ymin=ug.C.mg.1.mean.m-ug.C.mg.1.mean.ES, ymax=ug.C.mg.1.mean.m+ug.C.mg.1.mean.ES), width=.2,position=position_dodge(.9))+labs(y= expression(paste("Biomass-specific uptake (ug C mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(axis.title.x=element_blank(),axis.text.x = element_text(size=15,angle = 45,hjust=1),axis.text.y = element_text(size=12),legend.position=c(0.78,0.9),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton")) + annotate("text", x=0.6, y=30, label= "C",size=4) + scale_x_discrete(limits = positionST435)+scale_y_continuous(limits = c(-1,30))+ guides(fill=F)+
  annotate("text", x=0.75, y=-0.8, label= "n= 2",size=3)+ annotate("text", x=1.25, y=-0.8, label= "n= 2",size=3)+
  annotate("text", x=1.75, y=-0.8, label= "n= 1",size=3)+ annotate("text", x=2.25, y=-0.8, label= "n= 1",size=3)+
  annotate("text", x=2.75, y=-0.8, label= "n= 3",size=3)+ annotate("text", x=3.25, y=-0.8, label= "n= 2",size=3)+
  annotate("text", x=4, y=-0.8, label= "n= 1",size=3)


T.ST435.N<- ggplot(T.byfamilyST435, aes( x=FAMILY, y=I.total..ug.N.m2..sum.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=I.total..ug.N.m2..sum.m-I.total..ug.N.m2..sum.ES, ymax=I.total..ug.N.m2..sum.m+I.total..ug.N.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Total uptake (ug N m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.text.y = element_text(size=12),legend.text=element_text(size=8),legend.position=c(0.85,0.95),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=200, label= "B",size=4) +scale_x_discrete(limits = positionST435)+guides(fill=F)

BS.ST435.N <- ggplot(BS.byfamily_meanST435, aes(y=ug.N.mg.1.mean.m, x=FAMILY,fill=TREATMENT)) + geom_col(color="black",position = position_dodge(1),na.rm = TRUE)+geom_errorbar(aes(ymin=ug.N.mg.1.mean.m-ug.N.mg.1.mean.ES, ymax=ug.N.mg.1.mean.m+ug.N.mg.1.mean.ES), width=.2,position=position_dodge(.9))+labs(y= expression(paste("Biomass-specific uptake (ug N mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(axis.title.x=element_blank(),axis.text.x = element_text(size=15,angle = 45,hjust=1),axis.text.y = element_text(size=12),legend.position=c(0.78,0.9),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton")) + annotate("text", x=0.6, y=30, label= "D",size=4) + scale_x_discrete(limits = positionST435)+scale_y_continuous(limits = c(-1, 30))+ guides(fill=F)+
  annotate("text", x=0.75, y=-0.8, label= "n= 2",size=3)+ annotate("text", x=1.25, y=-0.8, label= "n= 2",size=3)+
  annotate("text", x=1.75, y=-0.8, label= "n= 1",size=3)+ annotate("text", x=2.25, y=-0.8, label= "n= 1",size=3)+
  annotate("text", x=2.75, y=-0.8, label= "n= 3",size=3)+ annotate("text", x=3.25, y=-0.8, label= "n= 2",size=3)+
  annotate("text", x=4, y=-0.8, label= "n= 1",size=3) 

library(cowplot)
TandBS.ST435 <- plot_grid(T.ST435.C,T.ST435.N,BS.ST435.C,BS.ST435.N,ncol = 2, align = "v")
ggsave(filename = "TandBS.ST435.png",plot = TandBS.ST435,width=20,height = 15)

###########################################################################################################

#Carbon incorporation %-----------------
#S. hyperborea were 22.3 and 3.0%, with a C:N ratio of 3.95
#T. nordenskioeldii, 21.6 and 4.5%, with a C:N ratio of 3.52

TOTALuptakeBYstations_mean$mg.C.algae <- c(545,600,425,545,600,425)
TOTALuptakeBYstations_mean$mg.N.algae <- c(3*545/22.3,3*600/22.3,3*425/22.3,4.5*545/21.6,4.5*600/21.6,4.5*425/21.6)
TOTALuptakeBYstations_mean$I.total.mg.C.m2 <- TOTALuptakeBYstations_mean$I.total..ug.C.m2..sum.m/1000
TOTALuptakeBYstations_mean$I.total.mg.N.m2 <- TOTALuptakeBYstations_mean$I.total..ug.N.m2..sum.m/1000
TOTALuptakeBYstations_mean$I.total.mg.C.m2.ES <- TOTALuptakeBYstations_mean$I.total..ug.C.m2..sum.ES/1000
TOTALuptakeBYstations_mean$I.total.mg.N.m2.ES <- TOTALuptakeBYstations_mean$I.total..ug.N.m2..sum.ES/1000
##percentage of uptake by station and treatment
library(dplyr)
porcentagefunctionC <- function(I.total.mg.C.m2,mg.C.algae){I.total.mg.C.m2*100 /mg.C.algae}
porcentagefunctionN <- function(I.total.mg.N.m2,mg.N.algae){I.total.mg.N.m2*100 /mg.N.algae}
porcentagefunctionCporcen <- function(I.total.mg.C.m2.ES,mg.C.algae){I.total.mg.C.m2.ES*100 /mg.C.algae}
porcentagefunctionNporcen <- function(I.total.mg.N.m2.ES,mg.N.algae){I.total.mg.N.m2.ES*100 /mg.N.algae}
table1 <- mutate(TOTALuptakeBYstations_mean,percentageC= porcentagefunctionC(I.total.mg.C.m2,mg.C.algae),percentageC.ES= porcentagefunctionCporcen(I.total.mg.C.m2.ES,mg.C.algae),percentageN= porcentagefunctionN(I.total.mg.N.m2,mg.N.algae),percentageN.ES= porcentagefunctionNporcen(I.total.mg.N.m2.ES,mg.N.algae))


write.table(table1, file = "table1.txt", sep = ",", quote = FALSE, row.names = F)

##percentage total uptake by family
library(doBy)
TOTALuptakeBYfamily <-  summaryBy(I.total..ug.C.m2.+I.total..ug.N.m2.~ CORE + TREATMENT+ FAMILY +STATION , data = uptake_data_STATIONS, FUN = sum,na.rm = TRUE)
TOTALuptakeBYfamily_mean <-  summaryBy(I.total..ug.C.m2..sum + I.total..ug.N.m2..sum ~ FAMILY + TREATMENT + STATION , data = TOTALuptakeBYfamily, FUN = function(x) { c(m = mean(x), ES = sqrt(var(x)/length(x))) })

library(fossil)
TOTALuptakeBYfamily.MATRIX <- create.matrix(subset(TOTALuptakeBYfamily_mean,STATION=="ST435"), tax.name = "FAMILY", locality = "TREATMENT", abund=TRUE, abund.col="I.total..ug.C.m2..sum.m")
table.TOTALuptakeBYfamily <- prop.table(TOTALuptakeBYfamily.MATRIX,2)*100
colSums(prop.table(percent.abundance.MATRIX,2)*100)

###########################################################################################################

#total and biomass-specif N and C uptake by GROUP----------
uptake_data_STATIONS$GROUP3 <- ifelse(uptake_data_STATIONS$GROUP=="Polychaeta","Polychaeta",ifelse(uptake_data_STATIONS$GROUP=="Bivalve","Bivalve",ifelse(uptake_data_STATIONS$GROUP=="Crustacea","Others",ifelse(uptake_data_STATIONS$GROUP=="Cnidaria","Others","Others"))))
TOTALuptake.groups <-  summaryBy(I.total..ug.C.m2. + I.total..ug.N.m2.~ CORE + TREATMENT +STATION + GROUP3 , data = uptake_data_STATIONS, FUN = sum,na.rm = TRUE)
TOTALuptake.groups.mean <-  summaryBy(I.total..ug.C.m2..sum + I.total..ug.N.m2..sum ~ TREATMENT + STATION + GROUP3 , data = TOTALuptake.groups, FUN = function(x) { c(m = mean(x),n=length(x), ES = sqrt(var(x)/length(x))) })

TOTALuptake.groups.mean$GROUP2 <- paste(TOTALuptake.groups.mean$GROUP3,TOTALuptake.groups.mean$STATION,sep="")

library(ggplot2)
positionST <- c("PolychaetaST435","BivalveST435","OthersST435","PolychaetaST407","BivalveST407","OthersST407","PolychaetaST177","BivalveST177","OthersST177")
A <- ggplot(TOTALuptake.groups.mean, aes(x=GROUP2, y=I.total..ug.C.m2..sum.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=I.total..ug.C.m2..sum.m-I.total..ug.C.m2..sum.ES, ymax=I.total..ug.C.m2..sum.m+I.total..ug.C.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Total uptake (",mu,"g C m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(plot.title = element_text(size=12),legend.background = element_rect(fill="transparent"),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.text.y = element_text(size=15),legend.text=element_text(size=8),legend.position=c(0.85,0.90),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.7, y=3100, label= "A",size=4) +scale_x_discrete(limits = positionST)+ geom_vline(xintercept =c(6.5,10.5), linetype="dashed")+ geom_vline(xintercept =c(3.5,10.5), linetype="dashed")+ggtitle("Beaufort Sea (Stn 177)              Amundsen Gulf (Stn 407)         Baffing Bay (Stn 435)")+ annotate("text", x=7, y=2500, label= "*",size=5)+ annotate("text", x=5, y=200, label= "*",size=6)+ annotate("text", x=2, y=500, label= "*",size=6)


B <- ggplot(TOTALuptake.groups.mean, aes(x=GROUP2, y=I.total..ug.N.m2..sum.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=I.total..ug.N.m2..sum.m-I.total..ug.N.m2..sum.ES, ymax=I.total..ug.N.m2..sum.m+I.total..ug.N.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Total uptake (",mu,"g N m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(plot.title = element_text(size=12),legend.background = element_rect(fill="transparent"),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.text.y = element_text(size=15),legend.text=element_text(size=8),legend.position=c(0.85,0.90),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.7, y=640, label= "B",size=4) +scale_x_discrete(limits = positionST)+ geom_vline(xintercept =c(6.5,10.5), linetype="dashed")+ geom_vline(xintercept =c(3.5,10.5), linetype="dashed")+ggtitle("Beaufort Sea (Stn 177)            Amundsen Gulf (Stn 407)        Baffing Bay (Stn 435)")+ annotate("text", x=7, y=450, label= "*",size=6)+ annotate("text", x=5, y=45, label= "*",size=6)+ annotate("text", x=2, y=180, label= "*",size=6)


biomassspecificuptake.groups <-  summaryBy(ug.C.mg.1+ ug.N.mg.1~ CORE + TREATMENT +STATION + GROUP3 , data = uptake_data_STATIONS, FUN = mean,na.rm = TRUE)
biomassspecificuptake.groups_mean <-  summaryBy(ug.C.mg.1.mean+ ug.N.mg.1.mean ~ TREATMENT + STATION +GROUP3 , data = biomassspecificuptake.groups, FUN = function(x) { c(n=length(x),m = mean(x), ES = sqrt(var(x)/length(x))) })

biomassspecificuptake.groups_mean$GROUP2 <- paste(biomassspecificuptake.groups_mean$GROUP3,biomassspecificuptake.groups_mean$STATION,sep="")


#uptakes <- subset(abundance.biomas_byCORE.STATION.GROUP.TREATMENT,STATION=="ST435"& GROUP2=="Bivalve")
#uptakes$I.total..mg.C.m2..sum.log <- sqrt(uptakes$I.total..mg.C.m2..sum)
#shapiro.test(uptakes$I.total..mg.C.m2..sum.log)
#t.test(I.total..mg.C.m2..sum.log~ TREATMENT, data=uptakes,var.equal = FALSE)
#wilcox.test(I.total..mg.C.m2..sum ~ TREATMENT, data=uptakes)

C<- ggplot(biomassspecificuptake.groups_mean, aes(x=GROUP2, y=ug.C.mg.1.mean.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=ug.C.mg.1.mean.m-ug.C.mg.1.mean.ES, ymax=ug.C.mg.1.mean.m+ug.C.mg.1.mean.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Biomass-specific uptake (",mu,"g C mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(axis.title.x=element_blank(),axis.text.x=element_text(size=15,hjust=0.5),axis.text.y = element_text(size=10),legend.position=c(0.8,0.9),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+guides(fill=FALSE)+guides(fill=FALSE)+ annotate("text", x=0.7, y=22, label= "C",size=4)+ scale_x_discrete(breaks=c("PolychaetaST435","BivalveST435","OthersST435","PolychaetaST407","BivalveST407","OthersST407","PolychaetaST177","BivalveST177","OthersST177"), labels=c("Pol.", "Biv.","Others","Pol.","Biv.","Others","Pol.", "Biv.","Others"),limits = positionST)+ geom_vline(xintercept =c(6.5,10.5), linetype="dashed")+ geom_vline(xintercept =c(3.5,10.5), linetype="dashed")+scale_y_continuous(limits = c(0, 22))+ annotate("text", x=7, y=5, label= "*",size=6)+ annotate("text", x=8, y=4, label= "**",size=6)

D <- ggplot(biomassspecificuptake.groups_mean, aes(x=GROUP2, y=ug.N.mg.1.mean.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=ug.N.mg.1.mean.m-ug.N.mg.1.mean.ES, ymax=ug.N.mg.1.mean.m+ug.N.mg.1.mean.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Biomass-specific uptake (",mu,"g N mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 15,base_family = "Helvetica") + theme(axis.title.x=element_blank(),axis.text.x=element_text(size=15),axis.text.y = element_text(size=10),legend.position=c(0.8,0.9),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+guides(fill=FALSE)+guides(fill=FALSE)+ annotate("text", x=0.7, y=22, label= "D",size=4)+ scale_x_discrete(breaks=c("PolychaetaST435","BivalveST435","OthersST435","PolychaetaST407","BivalveST407","OthersST407","PolychaetaST177","BivalveST177","OthersST177"), labels=c("Pol.", "Biv.","Others","Pol.","Biv.","Others","Pol.", "Biv.","Others"),limits = positionST)+ geom_vline(xintercept =c(6.5,10.5), linetype="dashed")+ geom_vline(xintercept =c(3.5,10.5), linetype="dashed")+scale_y_continuous(limits = c(0, 22))+ annotate("text", x=7, y=3, label= "*",size=6)+ annotate("text", x=8, y=3, label= "*",size=6)


library(cowplot)
fig2 <- plot_grid(A,B,C,D,ncol = 2,nrow = 2, align = "v")
ggsave(filename = "Uptake.by.Groups.png",plot = fig2,width=15,height = 10)

##TOTAL uptake Stats 

uptakeTOTAL <-subset(TOTALuptake.groups,STATION=="ST177"& GROUP3=="Polychaeta")
#uptakeTOTAL <- summaryBy(I.total..ug.C.m2.+I.total..ug.N.m2.~ CORE + TREATMENT +STATION + GROUP3 , data = uptakeTOTAL, FUN = sum,na.rm = TRUE)
#uptakeTOTAL <-subset(uptakeTOTAL,CORE!="PP5")

#normally distributed "the samples come from a Normal distribution"
shapiro.test(uptakeTOTAL$I.total..ug.C.m2..sum)
shapiro.test(uptakeTOTAL$I.total..ug.N.m2..sum)
##t test 
t.test(I.total..ug.C.m2..sum~ TREATMENT, data=uptakeTOTAL)
t.test(I.total..ug.N.m2..sum~ TREATMENT, data=uptakeTOTAL)

#log + normally distributed "the samples come from a Normal distribution"
shapiro.test(log(abs(uptakeTOTAL$I.total..ug.C.m2..sum)))
shapiro.test(log(abs(uptakeTOTAL$I.total..ug.N.m2..sum)))
##t test with log
t.test(log(abs(I.total..ug.C.m2..sum))~ TREATMENT, data=uptakeTOTAL)
t.test(log(abs(I.total..ug.N.m2..sum))~ TREATMENT, data=uptakeTOTAL)

#no parametric Mann-Whitney-Wilcoxon Test
wilcox.test(I.total..ug.C.m2..sum~ TREATMENT, data=uptakeTOTAL)
wilcox.test(I.total..ug.N.m2..sum~ TREATMENT, data=uptakeTOTAL)

###BIOMASS-specif uptake Stats
biomassspecific <-subset(biomassspecificuptake.groups,STATION=="ST407"& GROUP3=="Others")
#biomassspecific <- summaryBy(I.total..ug.C.m2.+I.total..ug.N.m2.~ CORE + TREATMENT +STATION + GROUP3 , data = biomassspecific, FUN = sum,na.rm = TRUE)
#biomassspecific <-subset(biomassspecific,CORE!="PP5")

#normally distributed "the samples come from a Normal distribution"
shapiro.test(biomassspecific$ug.C.mg.1.mean)
shapiro.test(biomassspecific$ug.N.mg.1.mean)
##t test 
t.test(ug.C.mg.1.mean~ TREATMENT, data=biomassspecific)
t.test(ug.N.mg.1.mean~ TREATMENT, data=biomassspecific)

#log + normally distributed "the samples come from a Normal distribution"
shapiro.test(log(abs(biomassspecific$ug.C.mg.1.mean)))
shapiro.test(log(abs(biomassspecific$ug.N.mg.1.mean)))
##t test with log
t.test(log(abs(ug.C.mg.1.mean))~ TREATMENT, data=biomassspecific)
t.test(log(abs(ug.N.mg.1.mean))~ TREATMENT, data=biomassspecific)

#no parametric Mann-Whitney-Wilcoxon Test
wilcox.test(I.total..ug.C.m2..sum~ TREATMENT, data=biomassspecific)
wilcox.test(I.total..ug.N.m2..sum~ TREATMENT, data=biomassspecific)



###########################################################################################################

#Total and biomass especif uptake by polychaeta family---------
TOTALuptakeBYpolyFamily <-  summaryBy(I.total..ug.C.m2.+I.total..ug.N.m2.~ CORE + TREATMENT +STATION + FAMILY , data = subset(uptake_data_STATIONS,GROUP=="Polychaeta"), FUN = sum,na.rm = TRUE)
TOTALuptakeBYpolyFamily_mean <-  summaryBy(I.total..ug.C.m2..sum+ I.total..ug.N.m2..sum ~ TREATMENT + STATION + FAMILY , data = TOTALuptakeBYpolyFamily, FUN = function(x) { c(m = mean(x), ES = sqrt(var(x)/length(x))) })

biomassspecificuptakeBYpolyFamily <-  summaryBy(ug.C.mg.1+ ug.N.mg.1~ CORE + TREATMENT +STATION +FAMILY , data = subset(uptake_data_STATIONS,GROUP=="Polychaeta"), FUN = mean,na.rm = TRUE)
biomassspecificuptakeBYpolyFamily_mean <-  summaryBy(ug.C.mg.1.mean+ ug.N.mg.1.mean ~ TREATMENT + STATION +FAMILY , data = biomassspecificuptakeBYpolyFamily, FUN = function(x) { c(n=length(x),m = mean(x), ES = sqrt(var(x)/length(x))) })

A1 <- ggplot(subset(TOTALuptakeBYpolyFamily_mean, STATION=="ST177"), aes(x=FAMILY, y=I.total..ug.C.m2..sum.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=I.total..ug.C.m2..sum.m-I.total..ug.C.m2..sum.ES, ymax=I.total..ug.C.m2..sum.m+I.total..ug.C.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Total uptake (ug C m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.text.y = element_text(size=10),legend.text=element_text(size=8),legend.position=c(0.8,0.90),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=2500, label= "A",size=4) 

B1 <- ggplot(subset(TOTALuptakeBYpolyFamily_mean, STATION=="ST177"), aes(x=FAMILY, y=I.total..ug.N.m2..sum.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=I.total..ug.N.m2..sum.m-I.total..ug.N.m2..sum.ES, ymax=I.total..ug.N.m2..sum.m+I.total..ug.N.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Total uptake (ug N m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_blank(),axis.title.x=element_blank(),axis.text.y = element_text(size=10),legend.text=element_text(size=8),legend.position=c(0.8,0.90),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=390, label= "B",size=4)

C1 <- ggplot(subset(biomassspecificuptakeBYpolyFamily_mean, STATION=="ST177"), aes(x=FAMILY, y=ug.C.mg.1.mean.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=ug.C.mg.1.mean.m-ug.C.mg.1.mean.ES, ymax=ug.C.mg.1.mean.m+ug.C.mg.1.mean.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Biomass-specific uptake (ug C mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.text.y = element_text(size=10),legend.position=c(0.8,0.9),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+guides(fill=FALSE)+ annotate("text", x=0.6, y=10, label= "C",size=4)


D1 <- ggplot(subset(biomassspecificuptakeBYpolyFamily_mean, STATION=="ST177"), aes(x=FAMILY, y=ug.N.mg.1.mean.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=ug.N.mg.1.mean.m-ug.N.mg.1.mean.ES, ymax=ug.N.mg.1.mean.m+ug.N.mg.1.mean.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Biomass-specific uptake (ug N mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.text.y = element_text(size=10),legend.position=c(0.8,0.9),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+guides(fill=FALSE)+guides(fill=FALSE)+ annotate("text", x=0.6, y=10, label= "D",size=4)

library(cowplot)
plot_grid(A1,B1,C1,D1,ncol = 2, nrow = 2,align = "v",labels ="ST177",hjust = 0, vjust = 1.4)

###########################################################################################################

#TOTAL and BIOMASS SPECIFIC UPTAKES per FEEDING GROUP by each station-------
library(doBy)
totaluptakeBYFEED <-  summaryBy(I.total..ug.C.m2.+I.total..ug.N.m2.~ CORE + FEED.MODE + TREATMENT+ STATION , data = uptake_data_STATIONS, FUN = sum)
totaluptakeBYFEED_mean <-  summaryBy(I.total..ug.C.m2..sum+I.total..ug.N.m2..sum~ TREATMENT + STATION + FEED.MODE, data = totaluptakeBYFEED, FUN = function(x) { c(m = mean(x),n=length(x), ES = sqrt(var(x)/length(x))) })
totaluptakeBYFEED_mean$GROUP.FEED.SATION <- paste(totaluptakeBYFEED_mean$STATION,totaluptakeBYFEED_mean$FEED.MODE)

positionFeedMode <- c("FF/SDF","SDF","SSDF","P/S")
ST435.totalbyFEED <- ggplot(subset(totaluptakeBYFEED_mean,STATION=="ST435"), aes(y=I.total..ug.C.m2..sum.m, x=FEED.MODE,fill=TREATMENT)) + geom_col(color="black",position = position_dodge(1),na.rm = TRUE)+labs(title="Beaufort Sea Stn",y= expression(paste("Total uptake (ug C m"^"-2",")")),fill=element_blank())+ theme_classic()+theme(axis.title.x= element_blank(),axis.text.y = element_text(size=20),axis.text.x = element_text(size=10,hjust=0.5),axis.title.y = element_text(size=20))+geom_errorbar(aes(ymin=I.total..ug.C.m2..sum.m-I.total..ug.C.m2..sum.ES, ymax=I.total..ug.C.m2..sum.m+I.total..ug.C.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Total uptake (ug C m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_text(),axis.title.x=element_blank(),axis.text.y = element_text(size=10),legend.text=element_text(size=8),legend.position=c(0.78,0.90),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=4500, label= "A",size=4)+scale_x_discrete(limits = positionFeedMode)+
  annotate("text", x=0.75, y=-100, label= "n= 5",size=2.5)+ annotate("text", x=1.25, y=-100, label= "n= 4",size=2.5)+
  annotate("text", x=1.75, y=-100, label= "n= 4",size=2.5)+ annotate("text", x=2.25, y=-100, label= "n= 4",size=2.5)+
  annotate("text", x=2.75, y=-100, label= "n= 3",size=2.5)+ annotate("text", x=3.25, y=-100, label= "n= 3",size=2.5)+
  annotate("text", x=3.74, y=-100, label= "n= 1",size=2.5)+ annotate("text", x=4.25, y=-100, label= "n= 4",size=2.5)

ST407.totalbyFEED <- ggplot(subset(totaluptakeBYFEED_mean,STATION=="ST407"), aes(y=I.total..ug.C.m2..sum.m, x=FEED.MODE,fill=TREATMENT)) + geom_col(color="black",position = position_dodge(1),na.rm = TRUE)+labs(title="Amundsen Gulf Stn",y= expression(paste("Total uptake (ug C m"^"-2",")")),fill=element_blank())+ theme_classic()+theme(axis.title.x= element_blank(),axis.text.y = element_text(size=20),axis.text.x = element_text(size=10,hjust=0.5),axis.title.y = element_text(size=20))+geom_errorbar(aes(ymin=I.total..ug.C.m2..sum.m-I.total..ug.C.m2..sum.ES, ymax=I.total..ug.C.m2..sum.m+I.total..ug.C.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Total uptake (ug C m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_text(),axis.title.x=element_blank(),axis.text.y = element_text(size=10),legend.text=element_text(size=8),legend.position=c(0.85,0.90),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=4500, label= "B",size=4)+guides(fill=FALSE)+scale_x_discrete(limits = positionFeedMode)+
  annotate("text", x=0.75, y=-100, label= "n= 4",size=2.5)+ annotate("text", x=1.25, y=-100, label= "n= 5",size=2.5)+
  annotate("text", x=1.75, y=-100, label= "n= 3",size=2.5)+ annotate("text", x=2.25, y=-100, label= "n= 3",size=2.5)+
  annotate("text", x=2.75, y=-100, label= "n= 5",size=2.5)+ annotate("text", x=3.25, y=-100, label= "n= 3",size=2.5)+
  annotate("text", x=3.74, y=-100, label= "n= 2",size=2.5)+ annotate("text", x=4.25, y=-100, label= "n= 3",size=2.5)

ST177.totalbyFEED <- ggplot(subset(totaluptakeBYFEED_mean,STATION=="ST177"), aes(y=I.total..ug.C.m2..sum.m, x=FEED.MODE,fill=TREATMENT)) + geom_col(color="black",position = position_dodge(1),na.rm = TRUE)+labs(title="Baffin Bay Stn",y= expression(paste("Total uptake (ug C m"^"-2",")")),fill=element_blank())+ theme_classic()+theme(axis.title.x= element_blank(),axis.text.y = element_text(size=20),axis.text.x = element_text(size=10,hjust=0.5),axis.title.y = element_text(size=20))+geom_errorbar(aes(ymin=I.total..ug.C.m2..sum.m-I.total..ug.C.m2..sum.ES, ymax=I.total..ug.C.m2..sum.m+I.total..ug.C.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Total uptake (ug C m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_text(),axis.title.x=element_blank(),axis.text.y = element_text(size=10),legend.text=element_text(size=8),legend.position=c(0.85,0.90),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=4500, label= "C",size=4)+guides(fill=FALSE)+scale_x_discrete(limits = positionFeedMode)+
  annotate("text", x=0.75, y=-100, label= "n= 5",size=2.5)+ annotate("text", x=1.25, y=-100, label= "n= 5",size=2.5)+
  annotate("text", x=1.75, y=-100, label= "n= 3",size=2.5)+ annotate("text", x=2.25, y=-100, label= "n= 2",size=2.5)+
  annotate("text", x=3, y=-100, label= "n= 3",size=2.5)+
  annotate("text", x=3.74, y=-100, label= "n= 4",size=2.5)+ annotate("text", x=4.25, y=-100, label= "n= 4",size=2.5)

##############################################################################################################################

#Biomass specif uptake
biomassspecificuptakeBYFEED <-  summaryBy(ug.C.mg.1+ ug.N.mg.1~ CORE + FEED.MODE + TREATMENT+ STATION , data =subset(uptake_data_STATIONS,ug.C.mg.1>0 & ug.N.mg.1>0) , FUN = mean,na.rm = TRUE)
biomassspecificuptakeBYFEED_mean <-  summaryBy(ug.C.mg.1.mean+ ug.N.mg.1.mean ~ TREATMENT + STATION + FEED.MODE, data = biomassspecificuptakeBYFEED, FUN = function(x) { c(m = mean(x),n=length(x), ES = sqrt(var(x)/length(x))) })

ST435.biomassspecificuptakeBYFEED <- ggplot(subset(biomassspecificuptakeBYFEED_mean,STATION=="ST435"), aes(y=ug.C.mg.1.mean.m, x=FEED.MODE,fill=TREATMENT)) + geom_col(color="black",position = position_dodge(1),na.rm = TRUE)+labs(y= expression(paste("Biomass-specific uptake (ug C mg"^"-1",")")),fill=element_blank())+ theme_classic()+theme(axis.title.x= element_blank(),axis.text.y = element_text(size=20),axis.text.x = element_text(size=10,hjust=0.5),axis.title.y = element_text(size=20))+geom_errorbar(aes(ymin=ug.C.mg.1.mean.m-ug.C.mg.1.mean.ES, ymax=ug.C.mg.1.mean.m+ug.C.mg.1.mean.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Biomass-specific uptake (ug C mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_text(),axis.title.x=element_blank(),axis.text.y = element_text(size=10),legend.text=element_text(size=8),legend.position=c(0.78,0.90),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=18, label= "D",size=4)+scale_x_discrete(limits = positionFeedMode)+
  annotate("text", x=0.75, y=-0.4, label= "n= 4",size=2.5)+ annotate("text", x=1.25, y=-0.4, label= "n= 4",size=2.5)+
  annotate("text", x=1.75, y=-0.4, label= "n= 4",size=2.5)+ annotate("text", x=2.25, y=-0.4, label= "n= 4",size=2.5)+
  annotate("text", x=2.75, y=-0.4, label= "n= 3",size=2.5)+ annotate("text", x=3.25, y=-0.4, label= "n= 3",size=2.5)+
  annotate("text", x=4, y=-0.4, label= "n= 4",size=2.5)

ST407.biomassspecificuptakeBYFEED <- ggplot(subset(biomassspecificuptakeBYFEED_mean,STATION=="ST407"), aes(y=ug.C.mg.1.mean.m, x=FEED.MODE,fill=TREATMENT)) + geom_col(color="black",position = position_dodge(1),na.rm = TRUE)+labs(y= expression(paste("Biomass-specific uptake (ug C mg"^"-1",")")),fill=element_blank())+ theme_classic()+theme(axis.title.x= element_blank(),axis.text.y = element_text(size=20),axis.text.x = element_text(size=10,hjust=0.5),axis.title.y = element_text(size=20))+geom_errorbar(aes(ymin=ug.C.mg.1.mean.m-ug.C.mg.1.mean.ES, ymax=ug.C.mg.1.mean.m+ug.C.mg.1.mean.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Biomass-specific uptake (ug C mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_text(),axis.title.x=element_blank(),axis.text.y = element_text(size=10),legend.text=element_text(size=8),legend.position=c(0.78,0.90),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=18, label= "E",size=4)+scale_x_discrete(limits = positionFeedMode)+guides(fill=FALSE)+scale_x_discrete(limits = positionFeedMode)+
  annotate("text", x=0.75, y=-0.4, label= "n= 4",size=2.5)+ annotate("text", x=1.25, y=-0.4, label= "n= 5",size=2.5)+
  annotate("text", x=1.75, y=-0.4, label= "n= 3",size=2.5)+ annotate("text", x=2.25, y=-0.4, label= "n= 3",size=2.5)+
  annotate("text", x=2.75, y=-0.4, label= "n= 5",size=2.5)+ annotate("text", x=3.25, y=-0.4, label= "n= 3",size=2.5)+
  annotate("text", x=3.74, y=-0.4, label= "n= 2",size=2.5)+ annotate("text", x=4.25, y=-0.4, label= "n= 3",size=2.5)

ST177.biomassspecificuptakeBYFEED <- ggplot(subset(biomassspecificuptakeBYFEED_mean,STATION=="ST177"), aes(y=ug.C.mg.1.mean.m, x=FEED.MODE,fill=TREATMENT)) + geom_col(color="black",position = position_dodge(1),na.rm = TRUE)+labs(y= expression(paste("Biomass-specific uptake (ug C mg"^"-1",")")),fill=element_blank())+ theme_classic()+theme(axis.title.x= element_blank(),axis.text.y = element_text(size=20),axis.text.x = element_text(size=10,hjust=0.5),axis.title.y = element_text(size=20))+geom_errorbar(aes(ymin=ug.C.mg.1.mean.m-ug.C.mg.1.mean.ES, ymax=ug.C.mg.1.mean.m+ug.C.mg.1.mean.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Biomass-specific uptake (ug C mg"^"-1",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(legend.background = element_rect(fill="transparent"),axis.text.x = element_text(),axis.title.x=element_blank(),axis.text.y = element_text(size=10),legend.text=element_text(size=8),legend.position=c(0.78,0.90),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton"))+ annotate("text", x=0.6, y=18, label= "F",size=4)+scale_x_discrete(limits = positionFeedMode)+guides(fill=FALSE)+
  annotate("text", x=0.75, y=-0.4, label= "n= 5",size=2.5)+ annotate("text", x=1.25, y=-0.4, label= "n= 5",size=2.5)+
  annotate("text", x=1.75, y=-0.4, label= "n= 3",size=2.5)+ annotate("text", x=2.25, y=-0.4, label= "n= 2",size=2.5)+
  annotate("text", x=3, y=-0.4, label= "n= 3",size=2.5)+
  annotate("text", x=3.74, y=-0.4, label= "n= 4",size=2.5)+ annotate("text", x=4.25, y=-0.4, label= "n= 4",size=2.5)


plot.TandBSuptake_BYFEED <- plot_grid(ST435.totalbyFEED,ST407.totalbyFEED,ST177.totalbyFEED,ST435.biomassspecificuptakeBYFEED,ST407.biomassspecificuptakeBYFEED,ST177.biomassspecificuptakeBYFEED,ncol = 3, align = "v")
ggsave(filename = "plot.TandBSuptake_BYFEED.png",plot = plot.TandBSuptake_BYFEED,width=12,height = 8)



##Stats
#t-test TOTAL UPTAKES BY FEED MODE-------
TOTALuptakeByFEEDMODE <- subset(totaluptakeBYFEED,STATION=="ST435" & FEED.MODE=="P/S")
#TOTALuptakeByFEEDMODE<- subset(TOTALuptakeByFEEDMODE,CORE!="IA2")

#normally distributed "the samples come from a Normal distribution"
shapiro.test(TOTALuptakeByFEEDMODE$I.total..ug.C.m2..sum)
shapiro.test(TOTALuptakeByFEEDMODE$I.total..ug.N.m2..sum)
#t-test
t.test(I.total..ug.C.m2..sum~ TREATMENT, data=TOTALuptakeByFEEDMODE)
t.test(I.total..ug.N.m2..sum~ TREATMENT, data=TOTALuptakeByFEEDMODE)

#log + normally distributed "the samples come from a Normal distribution"
shapiro.test(log(abs(TOTALuptakeByFEEDMODE$I.total..ug.C.m2..sum)))
shapiro.test(log(abs(TOTALuptakeByFEEDMODE$I.total..ug.N.m2..sum)))
#t-test
t.test(log(abs(I.total..ug.C.m2..sum))~ TREATMENT, data=TOTALuptakeByFEEDMODE)
t.test(log(abs(I.total..ug.N.m2..sum))~ TREATMENT, data=TOTALuptakeByFEEDMODE)

#no parametric Mann-Whitney-Wilcoxon Test
wilcox.test(I.total..ug.C.m2..sum~ TREATMENT, data=TOTALuptakeByFEEDMODE)
wilcox.test(abs(I.total..ug.N.m2..sum)~ TREATMENT, data=TOTALuptakeByFEEDMODE)


##Biomass-specific uptake stats
biomassspecificuptakeByFEEDMODE <- subset(biomassspecificuptakeBYFEED,STATION=="ST435" & FEED.MODE=="P/S" )
biomassspecificuptakeByFEEDMODE
#TOTALuptakeByFEEDMODE<- subset(TOTALuptakeByFEEDMODE,CORE!="IA2")

#normally distributed "the samples come from a Normal distribution"
shapiro.test(biomassspecificuptakeByFEEDMODE$ug.C.mg.1.mean)
shapiro.test(biomassspecificuptakeByFEEDMODE$ug.N.mg.1.mean)
#t-test
t.test(ug.C.mg.1.mean~ TREATMENT, data=biomassspecificuptakeByFEEDMODE)
t.test(ug.N.mg.1.mean~ TREATMENT, data=biomassspecificuptakeByFEEDMODE)

#log + normally distributed "the samples come from a Normal distribution"
shapiro.test(log(abs(biomassspecificuptakeByFEEDMODE$ug.C.mg.1.mean)))
shapiro.test(log(abs(biomassspecificuptakeByFEEDMODE$ug.N.mg.1.mean)))
#t-test
t.test(log(abs(ug.C.mg.1.mean))~ TREATMENT, data=biomassspecificuptakeByFEEDMODE)
t.test(log(abs(ug.N.mg.1.mean))~ TREATMENT, data=biomassspecificuptakeByFEEDMODE)

#no parametric Mann-Whitney-Wilcoxon Test
wilcox.test(ug.C.mg.1.mean~ TREATMENT, data=biomassspecificuptakeByFEEDMODE)
wilcox.test(ug.N.mg.1.mean~ TREATMENT, data=biomassspecificuptakeByFEEDMODE)

###########################################################################################################
#TOTAL UPTAKE by spionidae in all stations--------- 
library(doBy)
TOTALuptakeBYspionidae <-  summaryBy(I.total..ug.C.m2.+I.total..ug.N.m2.~ CORE + TREATMENT +STATION+FAMILY , data = uptake_data_STATIONS, FUN = sum,na.rm = TRUE)
TOTALuptakeBYspionidae <-  summaryBy(I.total..ug.C.m2..sum+ I.total..ug.N.m2..sum ~ FAMILY +TREATMENT + STATION , data =TOTALuptakeBYspionidae, FUN = function(x) { c(m = mean(x),n=length(x), ES = sqrt(var(x)/length(x))) })
TOTALuptakeBYspionidae$FAMILY2 <- paste(TOTALuptakeBYspionidae$FAMILY,TOTALuptakeBYspionidae$STATION)

positionsonlySPi <- c("Spionidae ST435","Spionidae ST407","Spionidae ST177")
ggplot(TOTALuptakeBYspionidae, aes( x=FAMILY2, y=I.total..ug.C.m2..sum.m,fill=TREATMENT))+ geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +geom_errorbar(aes(ymin=I.total..ug.C.m2..sum.m-I.total..ug.C.m2..sum.ES, ymax=I.total..ug.C.m2..sum.m+I.total..ug.C.m2..sum.ES), width=.2,position=position_dodge(.9))+ labs(y= expression(paste("Total uptake (ug C m"^"-2",")")),fill=element_blank())+ theme_bw(base_size = 12,base_family = "Helvetica") + theme(axis.text.x = element_blank(),axis.title.x=element_blank(),axis.text.y = element_text(size=10),legend.position=c(0.8,0.9),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values=c('white','gray48'),labels=c("Ice algae", "Phytoplankton")) +guides(fill=FALSE)+scale_x_discrete(limits = positionsonlySPi)


###########################################################################################################
# Biomass-specific C:N uptake ratio--------
library(doBy)
biomassspecificuptake <-  summaryBy(ug.C.mg.1+ ug.N.mg.1~ CORE + GROUP + TREATMENT +STATION , data = uptake_data_STATIONS, FUN = mean,na.rm = TRUE)
biomassspecificuptake_mean <-  summaryBy(ug.C.mg.1.mean+ ug.N.mg.1.mean ~ GROUP+ TREATMENT + STATION , data =biomassspecificuptake, FUN = function(x) { c(m = mean(x), ES = sqrt(var(x)/length(x))) })

biomassspecificuptake_mean$ratio <- biomassspecificuptake_mean$ug.C.mg.1.mean.m/biomassspecificuptake_mean$ug.N.mg.1.mean.m
biomassspecificuptake_mean$TREATMENT.STATION <- paste(biomassspecificuptake_mean$TREATMENT,biomassspecificuptake_mean$STATION)

biomassspecificuptake_mean <- subset(biomassspecificuptake_mean,GROUP!="Echinodermata")
library(ggplot2)
positionRATIO <- c("IA ST435","PP ST435","IA ST407","PP ST407","IA ST177","PP ST177")
ratios <- ggplot(biomassspecificuptake_mean, aes(x=TREATMENT.STATION, y=ratio))+
  theme_bw(base_size = 12,base_family = "Helvetica") + theme(axis.text.x = element_text(),axis.title.x=element_text(),axis.text.y = element_text(size=10),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),legend.text=element_text())+
  geom_point(aes(shape = factor(GROUP)),size=4)+
  scale_shape_manual(values=c(2,6,5,3,1,7,10))+
  labs(y="Biomass-specific C:N uptake",x="Stn 435                                          Stn 407                                      Stn 177")+scale_y_continuous(limits = c(0, 4))+
  scale_x_discrete(breaks=c("IA ST435","PP ST435","","IA ST407","PP ST407","","IA ST177", "PP ST177"), labels=c("Ice algae","Phytoplankton","","Ice algae", "Phytoplankton","","Ice algae","Phytoplankton"),limits = positionRATIO)+
  geom_point(aes(y=3.95,x=1), shape=4,size=4)+
  geom_point(aes(y=3.52,x=2), shape=4,size=4)+  
  geom_point(aes(y=3.95,x=3), shape=4,size=4)+
  geom_point(aes(y=3.52,x=4), shape=4,size=4)+
  geom_point(aes(y=3.95,x=5), shape=4,size=4)+
  geom_point(aes(y=3.52,x=6), shape=4,size=4)


ggsave(filename = "ratios.png",plot = ratios,width=10,height = 5)

#Biomass-specific C:N uptake ratio for species--------- 

Biomassspecific_CN <-  summaryBy(ug.C.mg.1+ ug.N.mg.1~ CORE + TREATMENT +STATION +FAMILY , data = uptake_data_STATIONS, FUN = mean,na.rm = TRUE)
Biomassspecific_CN_mean <-  summaryBy(ug.C.mg.1.mean+ ug.N.mg.1.mean ~ TREATMENT + STATION + FAMILY, data = Biomassspecific_CN, FUN = function(x) { c(m = mean(x), ES = sqrt(var(x)/length(x))) })
Biomassspecific_CN_mean$ratio <- Biomassspecific_CN_mean$ug.C.mg.1.mean.m/Biomassspecific_CN_mean$ug.N.mg.1.mean.m
Biomassspecific_CN_mean <- subset(Biomassspecific_CN_mean, ratio>0)



