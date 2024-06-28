##########ABUNDANCE AND BIOMASS in each STATION------
DATA <- "DATA"
allSTATIONS <- read.csv(file.path(DATA,"ST435-407-177-Pulse chase experiments.csv"))
allSTATIONS <- subset(allSTATIONS,GROUP!="Nematodes")
# using package doBY we calculated the mean, min, max, standar error (ES) and standard deviation(SD) in each core per station, subset with ABUNDANCE.ind.m.2>1|DRY.BIOMASS.mg.m.2>1 was used for not considering NA values in those columns

library(doBy)
##TOTAL ABUNDANCE AND BIOMASS in Station----------------------------------------------------------------------------------------
abundance.biomas_byCORE.STATION <- summaryBy(COUNTS+ ABUNDANCE+ BIOMASS+ Biomass.mg.C.m2 ~ CORE + STATION ,  data =subset(allSTATIONS,BIOMASS>1) ,FUN = sum,na.rm = TRUE)
abundance.biomas_bySTATION <- summaryBy(COUNTS.sum+ ABUNDANCE.sum+ BIOMASS.sum+ Biomass.mg.C.m2.sum ~STATION,   data =abundance.biomas_byCORE.STATION, FUN = function(x) { c(mean = mean(x), SD=sd(x),SE = sqrt(var(x)/length(x))) })

#PERMANOVA
abundance.biomas_byCORE.STATION$CODE <-  paste(abundance.biomas_byCORE.STATION$CORE,abundance.biomas_byCORE.STATION$STATION)
ab <- "density.abundance"
#select density or abundance
DENSITY.STATION.MATRIX<- matrix(abundance.biomas_byCORE.STATION$ABUNDANCE.sum)
DENSITY.STATION.MATRIX <- provideDimnames(DENSITY.STATION.MATRIX, base = list(abundance.biomas_byCORE.STATION$CODE,ab))
STATIONS <- as.character(abundance.biomas_byCORE.STATION$STATION)
library(vegan)
source("pairwise_adonis.R")
Density.station.ED <- vegdist(DENSITY.STATION.MATRIX,method="euclidean")
adonis(Density.station.ED ~ STATIONS, permutations = 999)
pairwise.adonis(DENSITY.STATION.MATRIX, STATIONS, sim.method = "euclidean")

##TOTAL ABUNDANCE AND BIOMASS in Station and layer 0-5 cm----------------------------------------------------------

abundance.biomas_byCORE.STATION <- summaryBy(COUNTS+ ABUNDANCE+ BIOMASS+ Biomass.mg.C.m2 ~ CORE + STATION ,  data =subset(allSTATIONS,BIOMASS>1 & LAYER=="0-5cm") ,FUN = sum,na.rm = TRUE)
abundance.biomas_bySTATION <- summaryBy(COUNTS.sum+ ABUNDANCE.sum+ BIOMASS.sum+ Biomass.mg.C.m2.sum ~STATION,   data =abundance.biomas_byCORE.STATION, FUN = function(x) { c(mean = mean(x), SD=sd(x),SE = sqrt(var(x)/length(x))) })

#PERMANOVA
abundance.biomas_byCORE.STATION$CODE <-  paste(abundance.biomas_byCORE.STATION$CORE,abundance.biomas_byCORE.STATION$STATION)
ab <- "density.abundance"
#select density or abundance
DENSITY.STATION.MATRIX<- matrix(abundance.biomas_byCORE.STATION$Biomass.mg.C.m2.sum)
DENSITY.STATION.MATRIX <- provideDimnames(DENSITY.STATION.MATRIX, base = list(abundance.biomas_byCORE.STATION$CODE,ab))
STATIONS <- as.character(abundance.biomas_byCORE.STATION$STATION)
library(vegan)
source("pairwise_adonis.R")
Density.station.ED <- vegdist(DENSITY.STATION.MATRIX,method="euclidean")
adonis(Density.station.ED ~ STATIONS, permutations = 999)
pairwise.adonis(DENSITY.STATION.MATRIX, STATIONS, sim.method = "euclidean")

##ABUNDANCE AND BIOMASS in Station and Treatments------------------------------------------------------------------
#calculated abundance and biomass at each station and treatment
abundance.biomas_byCORE.STATION.TREATMENT <- summaryBy(COUNTS+ ABUNDANCE+ BIOMASS+ Biomass.mg.C.m2 + I.total..mg.C.m2.+I.total..mg.N.m2. ~ CORE + STATION +TREATMENT ,  data =subset(allSTATIONS,BIOMASS>1 & TREATMENT!="C") ,FUN = sum,na.rm = TRUE)
abundance.biomas_bySTATION.TREATMENT <- summaryBy(COUNTS.sum+ ABUNDANCE.sum+ BIOMASS.sum+ Biomass.mg.C.m2.sum + I.total..mg.C.m2..sum+I.total..mg.N.m2..sum ~ STATION +TREATMENT ,  data =abundance.biomas_byCORE.STATION.TREATMENT, FUN = function(x) { c(mean = mean(x),SD=sd(x),SE = sqrt(var(x)/length(x))) })

#t-tests
variable <- subset(abundance.biomas_byCORE.STATION.TREATMENT,STATION=="ST407")
variable$I.total..mg.N.m2..sum.log <- log(variable$I.total..mg.N.m2..sum)
shapiro.test(variable$I.total..mg.N.m2..sum.log)
t.test(I.total..mg.N.m2..sum ~ TREATMENT, data=variable, var.equal = FALSE) #var.equal = FALSE for Welch approximation
wilcox.test(I.total..mg.N.m2..sum ~ TREATMENT, data=variable)

##ABUNDANCE AND BIOMASS in Station, Groups and Treatments-------------------------------------------------------------------------------
abundance.biomas_byCORE.STATION.GROUP.TREATMENT <- summaryBy( BIOMASS +Biomass.mg.C.m2 + I.total..mg.C.m2.+I.total..mg.N.m2. ~ CORE + STATION+ GROUP +TREATMENT ,  data =subset(allSTATIONS,BIOMASS>1 & TREATMENT!="C" ) ,FUN = sum,na.rm = TRUE)
abundance.biomas_byCORE.STATION.GROUP.TREATMENT$GROUP2 <-  ifelse(abundance.biomas_byCORE.STATION.GROUP.TREATMENT$GROUP=="Polychaeta","Polychaeta",ifelse(abundance.biomas_byCORE.STATION.GROUP.TREATMENT$GROUP=="Bivalve","Bivalve",ifelse(abundance.biomas_byCORE.STATION.GROUP.TREATMENT$GROUP=="Crustacea","Others",ifelse(abundance.biomas_byCORE.STATION.GROUP.TREATMENT$GROUP=="Cnidaria","Others","Others"))))
abundance.biomas_bySTATION.GROUP.TREATMENT <- summaryBy(BIOMASS.sum + Biomass.mg.C.m2.sum + I.total..mg.C.m2..sum +I.total..mg.N.m2..sum ~ STATION +GROUP2 + TREATMENT ,  data =abundance.biomas_byCORE.STATION.GROUP.TREATMENT, FUN = function(x) { c(mean = mean(x),SE = sqrt(var(x)/length(x))) })

#t-test
Biomass <- subset(abundance.biomas_byCORE.STATION.GROUP.TREATMENT,STATION=="ST407"& GROUP2=="Bivalve")
Biomass$Biomass.mg.C.m2.sum.log <- log(Biomass$Biomass.mg.C.m2.sum)
shapiro.test(Biomass$Biomass.mg.C.m2.sum.log)
t.test(Biomass.mg.C.m2.sum.log ~ TREATMENT, data=Biomass,var.equal = FALSE)
wilcox.test(Biomass.mg.C.m2.sum ~ TREATMENT, data=Biomass)

Biomass$BIOMASS.sum.log <- log(Biomass$BIOMASS.sum)
shapiro.test(Biomass$BIOMASS.sum.log)
t.test(BIOMASS.sum.log ~ TREATMENT, data=Biomass,var.equal = FALSE)
wilcox.test(Biomass.mg.C.m2.sum ~ TREATMENT, data=Biomass)

#t-test
uptakes <- subset(abundance.biomas_byCORE.STATION.GROUP.TREATMENT,STATION=="ST435"& GROUP2=="Bivalve")
uptakes$I.total..mg.C.m2..sum.log <- sqrt(uptakes$I.total..mg.C.m2..sum)
shapiro.test(uptakes$I.total..mg.C.m2..sum.log)
t.test(I.total..mg.C.m2..sum.log~ TREATMENT, data=uptakes,var.equal = FALSE)
wilcox.test(I.total..mg.C.m2..sum ~ TREATMENT, data=uptakes)




##TOTAL ABUNDANCE AND BIOMASS in Station by layer---------------------------------------------------------------------------------
abundance.biomas_byCORE.STATION.LAYER <- summaryBy(COUNTS+ ABUNDANCE+ BIOMASS+ Biomass.mg.C.m2 +Biomass.mg.N.m2 ~ CORE + STATION +LAYER,  data =subset(allSTATIONS,BIOMASS>1) ,FUN = sum,na.rm = TRUE)
abundance.biomas_bySTATION.LAYER <- summaryBy(COUNTS.sum+ ABUNDANCE.sum+ BIOMASS.sum+ Biomass.mg.C.m2.sum + Biomass.mg.N.m2.sum ~ STATION + LAYER,  data =abundance.biomas_byCORE.STATION.LAYER, FUN = function(x) { c(mean = mean(x),SD=sd(x),SE = sqrt(var(x)/length(x))) })

#PERMANOVA
abundance.biomas_byCORE.STATION.LAYER$CODE <-  paste(abundance.biomas_byCORE.STATION.LAYER$CORE,abundance.biomas_byCORE.STATION.LAYER$STATION,abundance.biomas_byCORE.STATION.LAYER$LAYER)
ab <- "density.abundance"
DENSITY.STATION.LAYER.MATRIX<- matrix(abundance.biomas_byCORE.STATION.LAYER$Biomass.mg.C.m2.sum)
DENSITY.STATION.LAYER.MATRIX <- provideDimnames(DENSITY.STATION.LAYER.MATRIX, base = list(abundance.biomas_byCORE.STATION.LAYER$CODE,ab))
STATIONS <- as.character(abundance.biomas_byCORE.STATION.LAYER$STATION)
LAYER <- as.character(abundance.biomas_byCORE.STATION.LAYER$LAYER)
FACTORS <- paste(STATIONS,LAYER)
#PERMANOVA DENSITY.LAYER.SATION-------
library(vegan)
source("pairwise_adonis.R")
Density.station.layer.ED <- vegdist(DENSITY.STATION.LAYER.MATRIX,method="euclidean")
adonis(Density.station.layer.ED ~ STATIONS + LAYER, permutations = 999)
pairwise.adonis(DENSITY.STATION.LAYER.MATRIX, FACTORS, sim.method = "euclidean")


## % of total biomass by depth
library(fossil)
percent.biomass.layer.MATRIX <- create.matrix(abundance.biomas_bySTATION.LAYER, tax.name = "LAYER", locality = "STATION", abund=TRUE, abund.col="ABUNDANCE.sum.mean")
table.porcent.biomass.layer <- prop.table(percent.biomass.layer.MATRIX,2)*100
colSums(prop.table(percent.biomass.layer.MATRIX,2)*100)


##ABUNDANCE AND BIOMASS in Station by FEEDMODE-------------------------------------------------------------------------------
abundance.biomas_byCORE.STATION.FEEDMODE <- summaryBy(COUNTS+ ABUNDANCE+ BIOMASS+ Biomass.mg.C.m2 ~ CORE + STATION +FEED.MODE,  data =subset(allSTATIONS,BIOMASS>1) ,FUN = sum,na.rm = TRUE)
abundance.biomas_bySTATION.FEEDMODE <- summaryBy(COUNTS.sum + ABUNDANCE.sum + BIOMASS.sum + Biomass.mg.C.m2.sum ~ FEED.MODE + STATION , data =abundance.biomas_byCORE.STATION.FEEDMODE ,FUN = function(x) { c(mean = mean(x),SD=sd(x)) })



#######################################################################################################################
#abundance and biomass by LAYER in each station PLOT-------
library(ggplot2)
abundance <- ggplot(abundance.biomas_bySTATION.LAYER, aes(y=ABUNDANCE.sum.mean, x=LAYER, fill=STATION)) + 
  geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +
  geom_errorbar(aes(ymin=ABUNDANCE.sum.mean, ymax=ABUNDANCE.sum.mean+ABUNDANCE.sum.SE), width=.2,position=position_dodge(.9),na.rm = TRUE)+geom_bar(position="dodge",stat="identity") +coord_flip()+labs(y= expression(paste("Density (ind. m"^"-2",")")),x="Depth (cm)",fill=element_blank())+scale_x_discrete(limits=c("5-10cm","0-5cm"),labels=c("5-10","0-5"))+ theme_bw(base_size = 12,base_family = "Helvetica")+theme(legend.text=element_text(size=15),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.title.x= element_text(size=20),legend.position=c(0.7,0.2),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20),axis.title.y = element_text(size=20))+ scale_fill_manual(values=c('gray48','gray76',"gray21"),breaks=c("ST435", "ST407", "ST177"),labels=c("Beaufort Sea (Stn 435)", "Amundsen Gulf (Stn 407)", "Baffin Bay (Stn 177)"))+ annotate("text", x=2.5, y=17000, label= "A",size=6) + guides(fill = guide_legend(keywidth = 2, keyheight = 2)) 
#+ annotate("text", x=1.68, y=15000, label= "*",size=8)
#+ annotate("text", x=0.98, y=800, label= "a,b",size=6)+ annotate("text", x=0.7, y=650, label= "a,b",size=6)+ annotate("text", x=1.32, y=500, label= "a",size=6)
#names(abundance.biomasSTATIONS_LAYER)

biomass<- ggplot(abundance.biomas_bySTATION.LAYER, aes(y=BIOMASS.sum.mean, x=LAYER, fill=STATION)) + 
  geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +
  geom_errorbar(aes(ymin=BIOMASS.sum.mean, ymax=BIOMASS.sum.mean+BIOMASS.sum.SE), width=.2,position=position_dodge(.9),na.rm = TRUE)+geom_bar(position="dodge",stat="identity") +coord_flip()+ scale_x_discrete(limits=c("5-10cm","0-5cm"),labels=c("5-10","0-5"))+labs(y= expression(paste("Biomass(mg DW.m"^"-2",")"))) + theme_bw(base_size = 12,base_family = "Helvetica")+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_text(size=20),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20),axis.title.y = element_blank()) + scale_fill_manual(values=c('gray48','gray76',"gray21"))+guides(fill=FALSE)+ annotate("text", x=2.5, y=6500, label= "B",size=6)
+ scale_y_continuous(limits = c(0, 2200))

biomass.mgCm2<- ggplot(abundance.biomas_bySTATION.LAYER, aes(y=Biomass.mg.C.m2.sum.mean, x=LAYER, fill=STATION)) + 
  geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +
  geom_errorbar(aes(ymin=Biomass.mg.C.m2.sum.mean, ymax=Biomass.mg.C.m2.sum.mean+Biomass.mg.C.m2.sum.SE), width=.2,position=position_dodge(.9),na.rm = TRUE)+geom_bar(position="dodge",stat="identity") +coord_flip()+ scale_x_discrete(limits=c("5-10cm","0-5cm"),labels=c("5-10","0-5"))+labs(y= expression(paste("Biomass(mg C m"^"-2",")"))) + theme_bw(base_size = 12,base_family = "Helvetica")+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_text(size=20),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20),axis.title.y = element_blank()) + scale_fill_manual(values=c('gray48','gray76',"gray21"))+guides(fill=FALSE)+ annotate("text", x=2.5, y=2000, label= "B",size=6)

biomass.mgNm2<- ggplot(abundance.biomas_bySTATION.LAYER, aes(y=Biomass.mg.N.m2.sum.mean, x=LAYER, fill=STATION)) + 
  geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +
  geom_errorbar(aes(ymin=Biomass.mg.N.m2.sum.mean, ymax=Biomass.mg.N.m2.sum.mean+Biomass.mg.N.m2.sum.SE), width=.2,position=position_dodge(.9),na.rm = TRUE)+geom_bar(position="dodge",stat="identity") +coord_flip()+ scale_x_discrete(limits=c("5-10cm","0-5cm"),labels=c("5-10","0-5"))+labs(y= expression(paste("Biomass(mg N.m"^"-2",")"))) + theme_bw(base_size = 12,base_family = "Helvetica")+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_text(size=20),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20),axis.title.y = element_blank()) + scale_fill_manual(values=c('gray48','gray76',"gray21"))+guides(fill=FALSE)+ annotate("text", x=2.5, y=500, label= "B",size=6)

library(cowplot)
abundance_biomassBYlayer <- plot_grid(abundance, biomass.mgCm2,ncol = 2, align = "hv",rel_widths=c(2,1.8))
ggsave(filename = "abundance_biomassBYlayer.png",plot = abundance_biomassBYlayer,width=15,height = 5)

#######################################################################################################################

# % of total BIOMMAS by FEED MODE----------
abundance.biomasFEED.MODE <- summaryBy(COUNTS + ABUNDANCE + BIOMASS + Biomass.mg.C.m2 ~ CORE + FEED.MODE + STATION   , data =subset(allSTATIONS,ABUNDANCE>1|BIOMASS>1) ,FUN = sum)
abundance.biomasFEED.MODE <- summaryBy(COUNTS.sum + ABUNDANCE.sum + BIOMASS.sum + Biomass.mg.C.m2.sum ~ FEED.MODE + STATION , data =abundance.biomasFEED.MODE ,FUN = function(x) { c(mean = mean(x),SD=sd(x)) })

library(fossil)
percent.biomassFEEDMODE.MATRIX <- create.matrix(abundance.biomasFEED.MODE, tax.name = "FEED.MODE", locality = "STATION", abund=TRUE, abund.col="Biomass.mg.C.m2.sum.mean")
table.percent.biomassFEEDMODE.MATRIX <- prop.table(percent.biomassFEEDMODE.MATRIX,2)*100
percent.biomassFEEDMODE <-  as.data.frame(table.percent.biomassFEEDMODE.MATRIX)
percent.biomassFEEDMODE$FEED.MODE <- c("FF","FF/SDF","P/S","SDF","SSDF")

stn177 <- ggplot(percent.biomassFEEDMODE, aes(y=ST177, x=FEED.MODE)) + 
  geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +
  coord_flip()+ scale_x_discrete(limits=c("SSDF","SDF","P/S","FF/SDF","FF"),labels=c("SSDF","SDF","P/S","FF/SDF","FF"))+labs(y="% of total C biomass") + theme_bw(base_size = 12,base_family = "Helvetica")+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_text(size=20),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20),axis.title.y = element_blank())+ annotate("text", x=5.2, y=51, label= "Baffing Bay (Stn 177)",size=6)+scale_y_continuous(limits = c(0, 79))

stn407 <- ggplot(percent.biomassFEEDMODE, aes(y=ST407, x=FEED.MODE, colour=FEED.MODE)) + 
  geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +
  coord_flip()+ scale_x_discrete(limits=c("SSDF","SDF","P/S","FF/SDF","FF"),labels=c("SSDF","SDF","P/S","FF/SDF","FF"))+labs(y="% of total C biomass") + theme_bw(base_size = 12,base_family = "Helvetica")+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_text(size=20),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20),axis.title.y = element_blank()) + scale_fill_manual(values=c('gray48','gray76',"gray21","white","black"))+ annotate("text", x=5.2, y=27, label= "Amundsen Gulf (Stn 407)",size=6)

stn435 <- ggplot(percent.biomassFEEDMODE, aes(y=ST435, x=FEED.MODE, colour=FEED.MODE)) + 
  geom_bar(stat="identity", color="black",position=position_dodge(),na.rm = TRUE) +
  coord_flip()+ scale_x_discrete(limits=c("SSDF","SDF","P/S","FF/SDF","FF"),labels=c("SSDF","SDF","P/S","FF/SDF","FF"))+labs(y="% of total C biomass") + theme_bw(base_size = 12,base_family = "Helvetica")+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_text(size=20),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20),axis.title.y = element_blank()) + scale_fill_manual(values=c('gray48','gray76',"gray21","white","black"))+ annotate("text", x=5.2, y=47, label= "Beaufort Sea (Stn 435)",size=6)


library(cowplot)
porcentagefeedmode <-plot_grid(stn435, stn407,stn177 ,ncol = 3, align = "hv")
ggsave(filename = "feedmode.png",plot = porcentagefeedmode,width=15,height = 5)


#######################################################################################################################
####BIODIVERSITY INFO, USING COUNT num-------
allSTATIONS <- subset(allSTATIONS, FAMILY!="unidentified")
library(doBy)
RichnessBYstation <- summaryBy(FAMILY~ STATION, data=subset(allSTATIONS, FAMILY!="Bulk sediments"), FUN= function (x) {(s=length (unique(x)))}) 
RichnessBYgroup.station <- summaryBy(FAMILY~ STATION + GROUP, data=subset(allSTATIONS, FAMILY!=""), FUN= function (x) {(s=length (unique(x)))})                 


##Create matrix of abundances---------
allSTATIONS.no.sediment <- subset(allSTATIONS,FAMILY!="Bulk sediments")
countsBY.STATIONS <- summaryBy(COUNTS~ FAMILY + CORE + STATION, data=allSTATIONS.no.sediment, FUN= sum,na.rm = TRUE) 
countsBY.STATIONS$CODE <- paste(countsBY.STATIONS$CORE,countsBY.STATIONS$STATION)

library(fossil)
FAMILY.MATRIX <- create.matrix(countsBY.STATIONS, tax.name = "FAMILY", locality = "CODE", abund=TRUE, abund.col="COUNTS.sum")
FAMILY.MATRIX <- t(FAMILY.MATRIX)

#add STATIONS GRoups
STATION=rep(c("Stn177","Stn407","Stn435"),15)

#MDS--------

#MDS no transformation
library(vegan)
#bray-curtis
family.bc <- vegdist(FAMILY.MATRIX, method = "bray")
family_NMDS=metaMDS(family.bc,k=2,trymax=100)
family_NMDS$stress
par(xpd=T, mar=par()$mar+c(0,0,0.5,0))
ordiplot(family_NMDS,type="n")
orditorp(family_NMDS,display="sites",cex=0.5,air=0.01)
ordihull(family_NMDS,STATION)
#colors to polygon source 
colors=rep(c('gray48',"gray76","gray21"),15)
for(i in unique(STATION)) {
  ordihull(family_NMDS$point[grep(i,STATION),],draw="polygon",
           groups=STATION[STATION==i],col=colors[grep(i,STATION)],label=F)
} 
text(-0.35314697,0.09,labels="Baffin Bay Stn",col="gray48")
text(0.42,-0.3,labels="Amundsen Gulf Stn",col="gray76")
text(0.48,0.41,labels="Beaufort Sea Stn",col="gray21")
text(-0.55,-0.32766912,labels="2D Stress: 0.18")
text(0.7,0.63,labels="no Transformation")
text(0.7,0.58,labels="Bray–Curtis dissimilarity")

#MDS with transformation-------
library(vegan)
family.trans <- sqrt(sqrt(FAMILY.MATRIX))
#bray-curtis
family.bc <- vegdist(family.trans, method = "bray")
family_NMDS=metaMDS(family.bc,k=2,trymax=100)


par(xpd=T, mar=par()$mar+c(0,0,1,0))
ordiplot(family_NMDS,type="n")
orditorp(family_NMDS,display="sites",cex=0.5,air=0.01)
ordihull(family_NMDS,STATION)
#colors to polygon source 
colors=rep(c('gray48',"gray76","gray21"),15)
for(i in unique(STATION)) {
  ordihull(family_NMDS$point[grep(i,STATION),],draw="polygon",
           groups=STATION[STATION==i],col=colors[grep(i,STATION)],label=F)
} 
text(-0.26890376,0.02111721,labels="Stn177")
text(-0.06390031,-0.01156450,labels="Stn407")
text(0.44712277,0.25583130,labels="Stn435")
text(0.68480793,-0.40077394,labels="2D Stress: 0.19")
text(0.4341107,0.4563443,labels="Transform: Fourth-root")
text(0.4341107,0.4145642,labels="Bray–Curtis dissimilarity")

#######################################################################################################################

#anosim
library(vegan)
family.bc.ano <- anosim(family.bc, STATION)
summary(family.bc.ano)
plot(family.bc.ano)

#PERMANOVA benthic community-------
adonis(family.bc ~ STATION, permutations = 9999)
source("pairwise_adonis.R") # Choisir le fichier pairwise_adonis.R
pairwise.adonis(FAMILY.MATRIX, STATION, sim.method = "bray")

#SIMPER----
sim <- simper(FAMILY.MATRIX, STATION, permutations = 999)
summary(sim)


#diversity index-----
H <- diversity(FAMILY.MATRIX)
simp <- diversity(FAMILY.MATRIX, "simpson")
S <- specnumber(FAMILY.MATRIX)
J <- H/log(S)
#alpha <- fisher.alpha(FAMILY.MATRIX)## Fisher alpha
#unbias.simp <- rarefy(FAMILY.MATRIX, 2) - 1

cnom <- c("H","simp","S","J")
#create matrix with biodiversity index
DIVERSITY <- matrix(c(H,simp,S,J),nrow = 45,ncol=4)
DIVERSITY <- provideDimnames(DIVERSITY, base = list(STATION, cnom))

DIVERSITY.mean <- as.data.frame(DIVERSITY)
DIVERSITY.mean$STATION <- rep(c("Stn177","Stn407","Stn435"),15)
DIVERSITY.mean <- summaryBy(H+simp+S+J ~ STATION, data=DIVERSITY.mean, FUN=mean)

#PERMANOVA DIVERSITY-------
source("pairwise_adonis.R")
DIVERSITY.ED <- vegdist(DIVERSITY,method="euclidean")
adonis(DIVERSITY.ED ~ STATION, permutations = 999)
pairwise.adonis(DIVERSITY, STATION, sim.method = "euclidean")


#######################################################################################################################


#PERCENTAGE of contribution per taxa or family-------
library(doBy)
FAMILYbySTATION_CORE <- summaryBy(COUNTS + ABUNDANCE + BIOMASS + Biomass.mg.C.m2~  FAMILY + STATION + CORE , data =subset(allSTATIONS, ABUNDANCE>1),FUN = sum)
FAMILYbySTATION <- summaryBy(COUNTS.sum +ABUNDANCE.sum+ BIOMASS.sum + Biomass.mg.C.m2.sum ~ STATION + FAMILY, data =FAMILYbySTATION_CORE, FUN = function(x) { c(mean = mean(x),SD=sd(x)) })

##Create matrix 
library(fossil)
percent.abundance.MATRIX <- create.matrix(FAMILYbySTATION, tax.name = "FAMILY", locality = "STATION", abund=TRUE, abund.col="ABUNDANCE.sum.mean")
table.porcent.abundance <- prop.table(percent.abundance.MATRIX,2)*100
colSums(prop.table(percent.abundance.MATRIX,2)*100)


#pie plots abundances-------
cols <- c("grey90","grey50","black","grey30","white","grey20","grey40","grey60","grey70","grey99","grey85","grey90","grey100","grey10","grey8","ghostwhite","darkgray","ivory4")

table.porcent.abundance.dataframe <- as.data.frame(table.porcent.abundance)
table.porcent.abundance.dataframe$FAMILY <-  rownames(table.porcent.abundance)
table.porcent.abundance.dataframe$ST177 <- round(table.porcent.abundance.dataframe$ST177,0)
table.porcent.abundance.dataframe$ST407 <- round(table.porcent.abundance.dataframe$ST407,0)
table.porcent.abundance.dataframe$ST435 <- round(table.porcent.abundance.dataframe$ST435,0)
table.porcent.abundance.dataframe$pielabels.ST177<- paste(table.porcent.abundance.dataframe$ST177, "%", sep="")
table.porcent.abundance.dataframe$pielabels.ST407<- paste(table.porcent.abundance.dataframe$ST407, "%", sep="")
table.porcent.abundance.dataframe$pielabels.ST435<- paste(table.porcent.abundance.dataframe$ST435, "%", sep="")

abundanceST177 <- subset(table.porcent.abundance.dataframe,ST177>2)
abundanceST177 <- abundanceST177[order(abundanceST177$ST177,decreasing = TRUE), ]
labsabundances.ST177<- paste("(",abundanceST177$FAMILY,")", " ", abundanceST177$pielabels.ST177, sep="")
pie(abundanceST177$ST177, main=expression(paste("Abundance (ind.m"^"-2",")")), col=cols, labels=labsabundances.ST177, cex=0.6)

abundanceST407 <- subset(table.porcent.abundance.dataframe,ST407>2)
abundanceST407 <- abundanceST407[order(abundanceST407$ST407,decreasing = TRUE), ]
labsabundances.ST407<- paste("(",abundanceST407$FAMILY,")", " ", abundanceST407$pielabels.ST407, sep="")
pie(abundanceST407$ST407, main=expression(paste("Abundance (ind.m"^"-2",")")), col=cols, labels=labsabundances.ST407, cex=0.6)

abundanceST435 <- subset(table.porcent.abundance.dataframe,ST435>2)
abundanceST435 <- abundanceST435[order(abundanceST435$ST435,decreasing = TRUE), ]
labsabundances.ST435<- paste("(",abundanceST435$FAMILY,")", " ", abundanceST435$pielabels.ST435, sep="")
pie(abundanceST435$ST435, main=expression(paste("Abundance (ind.m"^"-2",")")), col=cols, labels=labsabundances.ST435, cex=0.6)



# pie biomass-----

FAMILYbySTATION_CORE <- summaryBy(COUNTS + ABUNDANCE + BIOMASS + Biomass.mg.C.m2~  FAMILY + STATION + CORE , data =subset(allSTATIONS, BIOMASS>1),FUN = sum)
FAMILYbySTATION <- summaryBy(COUNTS.sum +ABUNDANCE.sum+ BIOMASS.sum + Biomass.mg.C.m2.sum ~ STATION + FAMILY, data =FAMILYbySTATION_CORE, FUN = function(x) { c(mean = mean(x),SD=sd(x)) })
library(fossil)
percent.biomass.MATRIX <- create.matrix(FAMILYbySTATION, tax.name = "FAMILY", locality = "STATION", abund=TRUE, abund.col="BIOMASS.sum.mean")
table.porcent.biomass <- prop.table(percent.biomass.MATRIX,2)*100
colSums(prop.table(percent.biomass.MATRIX,2)*100)
table.porcent.biomass.dataframe <- as.data.frame(table.porcent.biomass)
table.porcent.biomass.dataframe$FAMILY <-  rownames(table.porcent.biomass)
table.porcent.biomass.dataframe$ST177 <- round(table.porcent.biomass.dataframe$ST177,0)
table.porcent.biomass.dataframe$ST407 <- round(table.porcent.biomass.dataframe$ST407,0)
table.porcent.biomass.dataframe$ST435 <- round(table.porcent.biomass.dataframe$ST435,0)
table.porcent.biomass.dataframe$pielabels.ST177<- paste(table.porcent.biomass.dataframe$ST177, "%", sep="")
table.porcent.biomass.dataframe$pielabels.ST407<- paste(table.porcent.biomass.dataframe$ST407, "%", sep="")
table.porcent.biomass.dataframe$pielabels.ST435<- paste(table.porcent.biomass.dataframe$ST435, "%", sep="")

cols2 <- c("grey50","grey90","grey30","grey90","grey90","white","grey20","grey60","grey70","grey99","grey85","black","grey100","grey10","grey8","ghostwhite","darkgray","ivory4")

biomassST177 <- subset(table.porcent.biomass.dataframe,ST177>2)
biomassST177 <- biomassST177[order(biomassST177$ST177,decreasing = TRUE), ]
labsbiomass.ST177<- paste("(",biomassST177$FAMILY,")", " ", biomassST177$pielabels.ST177, sep="")
pie(biomassST177$ST177, main=expression(paste("Biomass(mg DW.m"^"-2",")")), col=cols2, labels=labsbiomass.ST177, cex=0.6)

biomassST407 <- subset(table.porcent.biomass.dataframe,ST407>2)
biomassST407 <- biomassST407[order(biomassST407$ST407,decreasing = TRUE), ]
labsbiomass.ST407<- paste("(",biomassST407$FAMILY,")", " ", biomassST407$pielabels.ST407, sep="")
pie(biomassST407$ST407, main=expression(paste("Biomass(mg DW.m"^"-2",")")), col=cols, labels=labsbiomass.ST407, cex=0.6)

biomassST435 <- subset(table.porcent.biomass.dataframe,ST435>2)
biomassST435 <- biomassST435[order(biomassST435$ST435,decreasing = TRUE), ]
labsbiomass.ST435<- paste("(",biomassST435$FAMILY,")", " ", biomassST435$pielabels.ST435, sep="")
pie(biomassST435$ST435, main=expression(paste("Biomass(mg DW.m"^"-2",")")), col=cols, labels=labsbiomass.ST435, cex=0.6)

par(mfrow = c(3, 2),     # 2x2 layout
    oma = c(0, 0, 1, 0), # c(bottom, left, top, right) margenes fuera del grafico
    mar = c(0, 1.3, 1.3, 0), # c(bottom, left, top, right) separacion entre plots
    mgp = c(0, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = T)            # allow content to protrude into outer margin (and beyond)
pie(abundanceST177$ST177, main=expression(paste("Abundance (ind.m"^"-2",")")), col=cols, labels=labsabundances.ST177, cex=0.8)
mtext("Baffin Bay Stn",side=2)
pie(biomassST177$ST177, main=expression(paste("Biomass(mg DW.m"^"-2",")")), col=cols2, labels=labsbiomass.ST177, cex=0.8)
pie(abundanceST407$ST407, col=cols, labels=labsabundances.ST407, cex=0.8)
mtext("Amundsen Gulf Stn",side=2)
pie(biomassST407$ST407, col=cols, labels=labsbiomass.ST407, cex=0.8)
pie(abundanceST435$ST435, col=cols, labels=labsabundances.ST435, cex=0.8)
mtext("Beaufort Sea Stn",side=2)
pie(biomassST435$ST435, col=cols, labels=labsbiomass.ST435, cex=0.8)



