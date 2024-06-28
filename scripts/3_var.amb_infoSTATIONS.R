
#Enviromental variables
#ST177 (baffin bay), ST407 (Amundsen Gulf) and ST435 (Beaufort Sea)
DATA <- "DATA"
ST177.ST407.ST435grainsize <- read.csv(file.path(DATA,"ST435-407-177-GrainSize.csv"))
rownames(ST177.ST407.ST435grainsize) <- ST177.ST407.ST435grainsize$microns
ST177.ST407.ST435grainsize$microns <- NULL
oxygen <- read.csv(file.path(DATA,"ST435-407-177-O2.csv"))
colnames(oxygen)[2] <- "STATION_NUMBER"
colnames(oxygen)[1] <- "STATION"
DIC <- read.csv(file.path(DATA,"ST435-407-177-DIC.csv"))
nutrients <- read.csv(file.path(DATA,"ST435-177-Nutrients.csv"))

#Librarys 
library(G2Sd)
library(plyr)
library(lattice)
library(sciplot)
library(doBy)
library(car) 
library(cowplot)

#calculate how many are labelled------
allSTATIONSnoControl <- subset(allSTATIONS, TREATMENT!="C")
allSTATIONSnoControl$C.ENRICHED <- ifelse(allSTATIONSnoControl$Formalin.corrected_δ13C_VPDB>-15,"YES.C","NO.C")
allSTATIONSnoControl$N.ENRICHED <- ifelse(allSTATIONSnoControl$δ15N_air>20,"YES.N","NO.N")

library(dplyr) 
dataSummary <- allSTATIONSnoControl%>% group_by(STATION,LAYER)%>%count(C.ENRICHED,N.ENRICHED)
library(data.table)
data <- data.table(dataSummary)
data <- subset(data,C.ENRICHED!="NA")

data[, ntotal := sum(n), by=list(STATION,LAYER)]
data[, percen := sum(n), by=list(STATION,LAYER)]
data[, percen := n/percen]
data

#GRAIN SIZE DATA---------
#Median grain size (µm) using G2sd package 
library(G2Sd)
granstat(ST177.ST407.ST435grainsize,aggr=F)$sed
granstat(ST177.ST407.ST435grainsize ,aggr=F)$stat


#OXYGEN DATA ---------
#Funcion for lineal regression between oxygen.corrected and Time AND returns Slope and Intercept
library(plyr)
regressionFuncion <- function(x){
  results <- lm(oxygen.corrected~TIME, data = x)
  return(summary(results)$coefficients[2:1])
}
regressionR2 <- function(x){
  results <- lm(oxygen.corrected~TIME, data = x)
  return(return(summary(results)$r.squared))
}


# ST177 lineal regression oxygen concentration vs time
ST177regressionR2 <- ddply(subset(oxygen,STATION=="ST177"), .(REPLICATE), regressionR2)
ST177regressionResults <- ddply(subset(oxygen,STATION=="ST177"), .(REPLICATE), regressionFuncion)
#plot each core lineal regresion
library(lattice)
xyplot(oxygen.corrected ~ TIME | REPLICATE,data =subset(oxygen,STATION=="ST177"),par.settings = list(strip.background=list(col="lightgrey")),panel = function(x, y,col="white", ...) {panel.xyplot(x, y, col="black",...)
  panel.abline(lm(y~x), col='black',lineheigth=2)
  },grid = TRUE,main = "ST177")

# ST407 lineal regression oxygen concentration vs time
ST407regressionR2 <- ddply(subset(oxygen,STATION=="ST407"), .(REPLICATE), regressionR2)
ST407regressionResults <- ddply(subset(oxygen,STATION=="ST407"), .(REPLICATE), regressionFuncion)
#plot each core lineal regresion
library(lattice)
xyplot(oxygen.corrected ~ TIME | REPLICATE,data =subset(oxygen,STATION=="ST407"),par.settings = list(strip.background=list(col="lightgrey")),panel = function(x, y,col="white", ...) {panel.xyplot(x, y, col="black",...)
  panel.abline(lm(y~x), col='black',lineheigth=2)
},grid = TRUE,main = "ST407")

# ST435 lineal regression oxygen concentration vs time
ST435regressionR2 <- ddply(subset(oxygen,STATION=="ST435"), .(REPLICATE), regressionR2)
ST435regressionResults <- ddply(subset(oxygen,STATION=="ST435"), .(REPLICATE), regressionFuncion)
#plot each core lineal regresion
library(lattice)
xyplot(oxygen.corrected ~ TIME | REPLICATE,data =subset(oxygen,STATION=="ST435"),par.settings = list(strip.background=list(col="lightgrey")),panel = function(x, y,col="white", ...) {panel.xyplot(x, y, col="black",...)
  panel.abline(lm(y~x), col='black',lineheigth=2)
},grid = TRUE,main = "ST435")

#R squared for cores in each station
r.squared <- data.frame(CORES=ST177regressionR2$REPLICATE,ST177=ST177regressionR2$V1)
r.squared$ST407 <- ST407regressionR2$V1
r.squared$ST435 <- ST435regressionR2$V1
print(r.squared)

#SCOC PLOT------
#(calculation in excel USING SCOC(umol O2 m-2 d-1)=SLOPE* VOLH2O(M3)/AREA (M2))
library(sciplot)
bargraph.CI(STATION,SCOC,TREATMENT, data = oxygen ,legend=T,y.leg=-3,x.leg=10,ylab="SCOC (mmol m–2 d–1)",xlab="TREATMENT",ylim = c(-6,1))
abline(h=0)
text(x=1.5, y=-3, "a", pos=3, cex=1.2)
text(x=2.5, y=-4.6, "b", pos=3, cex=1.2)
text(x=3.5, y=-5.6, "b", pos=3, cex=1.2)
#text(x=5.4, y=0.5, "c", pos=3, cex=1.2)
#text(x=6.4, y=-1, "c", pos=3, cex=1.2)
#text(x=7.4, y=-1.5, "a,c", pos=3, cex=1.2)
text(x=9.4, y=-2, "a", pos=3, cex=1.2)
text(x=10.4, y=-2, "a", pos=3, cex=1.2)
text(x=11.4, y=-2, "a", pos=3, cex=1.2)
#ANOVA SCOC ----

#ST177
oxygenST177 <- subset(oxygen,STATION=="ST177")
result_aovST177 <- aov(SCOC ~ TREATMENT, data = subset(oxygenST177,SCOC!="NA")) 
summary(result_aovST177)
TukeyHSD(result_aovST177,"TREATMENT")
#normally distributed "the samples come from a Normal distribution"
shapiro.test(result_aovST177$residuals)
#Homogeneity of Variances ho; v1=v2=v3
bartlett.test(SCOC ~ TREATMENT, data =oxygen)
fligner.test(SCOC ~ TREATMENT, data = oxygen)
layout(matrix(c(1:6), 2, 3))
plot(result_aovST177, 1:6)
layout(1)

#ST407 (NO DATA )
oxygenST407 <- subset(oxygen,STATION=="ST407")
#result_aovST407 <- aov(SCOC ~ TREATMENT, data = subset(oxygenST407,SCOC!="NA")) 
#summary(result_aovST407)
#TukeyHSD(result_aovST407,"TREATMENT")
#normally distributed "the samples come from a Normal distribution"
#shapiro.test(result_aovST407$residuals)
#Homogeneity of Variances ho; v1=v2=v3
#bartlett.test(SCOC ~ TREATMENT, data =oxygen)
#fligner.test(SCOC ~ TREATMENT, data = oxygen)
#layout(matrix(c(1:6), 2, 3))
#plot(result_aovST407, 1:6)
#layout(1)

#ST435
oxygenST435 <- subset(oxygen,STATION=="ST435")
result_aovST435 <- aov(SCOC ~ TREATMENT, data = subset(oxygenST435,SCOC!="NA")) 
summary(result_aovST435)
TukeyHSD(result_aovST435,"TREATMENT")
#normally distributed "the samples come from a Normal distribution"
shapiro.test(result_aovST435$residuals)
#Homogeneity of Variances ho; v1=v2=v3
bartlett.test(SCOC ~ TREATMENT, data =oxygen)
fligner.test(SCOC ~ TREATMENT, data = oxygen)
layout(matrix(c(1:6), 2, 3))
plot(result_aovST435, 1:6)
layout(1)

#DIC-----
library(sciplot)
par(mfrow = c(1, 3),     # 2x2 layout
    oma = c(2, 1, 1, 1), # c(bottom, left, top, right) margenes fuera del grafico
    mar = c(4, 3.5, 1, 0), # c(bottom, left, top, right) separacion entre plots
    mgp = c(2.6, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = T)            # allow content to protrude into outer margin (and beyond)
bargraph.CI(TIME,DIC.mg.m2.,TREATMENT, data = subset(DIC,STATION=="ST177"),yaxt="n",ylab="DIC (mg C m–2)",ylim=c(0,100),main="ST177")
axis(2,las=1,tck=-.01,col = "gray",cex.axis=0.8)
segments(0.7, 0, 8, 0, col= 'gray')
box(lty = 'solid', col = 'gray')
bargraph.CI(TIME,DIC.mg.m2.,TREATMENT, data = subset(DIC,STATION=="ST407"),xlab="TIME",yaxt="n",ylim=c(0,100),main="ST407")
segments(0.7, 0, 8, 0, col= 'gray')
box(lty = 'solid', col = 'gray')
bargraph.CI(TIME,DIC.mg.m2.,TREATMENT, data = subset(DIC,STATION=="ST435"),yaxt="n",ylim=c(0,100),main="ST435",legend=T)
segments(0.7, 0, 8, 0, col= 'gray')
box(lty = 'solid', col = 'gray')


#NUTRIENTS (ST435,ST177)-----
nutrients <- read.csv(file.path(DATA,"ST435-177-Nutrients.csv"))
#cores water volume and area 
coresVOL.AREA <- read.csv(file.path(DATA,"ST435-407-177-coresVOLUME.AREA.csv"))

#funcion for lineal regresion. nutrient cc vs Time in each core
library(plyr)
regressionFuncion.nutrients <- function(x){
  results <- lm(Nutrients~TIME, data = x)
  return(summary(results)$coefficients[2:1])
}

regressionR2.nutrients <- function(x){
  results <- lm(Nutrients~TIME, data = x)
  return(return(summary(results)$r.squared))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ST177nutrients------
#ST177 Nitrite
ST177nitrite <- nutrients[nutrients$STATION== "ST177" & nutrients$TYPE== "Nitrite",]
ST177nitrite.regression <- ddply(ST177nitrite, .(REPLICATE), regressionFuncion.nutrients)
ST177nitrite.R2 <- ddply(ST177nitrite, .(REPLICATE), regressionR2.nutrients)

#plot each core lineal regresion
library(lattice)
xyplot(Nutrients ~ TIME | REPLICATE,data =ST177nitrite,par.settings = list(strip.background=list(col="lightgrey")),panel = function(x, y,col="white", ...) {panel.xyplot(x, y, col="black",...)
  panel.abline(lm(y~x), col='black',lineheigth=2)
},grid = TRUE,main = "nitriteST177")
#include slope of regresion in the data frame
coresVOL.AREA[1:15,6] <- ST177nitrite.regression$V1 
names(coresVOL.AREA)[6]<-"slope.nitrite"

#ST177 Nitrate
ST177Nitrate <- nutrients[nutrients$STATION== "ST177" & nutrients$TYPE== "Nitrate",]
ST177Nitrate.regression <- ddply(ST177Nitrate, .(REPLICATE), regressionFuncion.nutrients)
#plot each core lineal regresion
library(lattice)
xyplot(Nutrients ~ TIME | REPLICATE,data =ST177Nitrate,par.settings = list(strip.background=list(col="lightgrey")),panel = function(x, y,col="white", ...) {panel.xyplot(x, y, col="black",...)
  panel.abline(lm(y~x), col='black',lineheigth=2)
},grid = TRUE,main = "NitrateST177")
#include slope of regresion in the data frame
coresVOL.AREA[1:15,7] <- ST177Nitrate.regression$V1 
names(coresVOL.AREA)[7]<-"slope.Nitrate"

#ST177 Phosphate
ST177Phosphate <- nutrients[nutrients$STATION== "ST177" & nutrients$TYPE== "Phosphate",]
ST177Phosphate.regression <- ddply(ST177Phosphate, .(REPLICATE), regressionFuncion.nutrients)
#plot each core lineal regresion
library(lattice)
xyplot(Nutrients ~ TIME | REPLICATE,data =ST177Phosphate,par.settings = list(strip.background=list(col="lightgrey")),panel = function(x, y,col="white", ...) {panel.xyplot(x, y, col="black",...)
  panel.abline(lm(y~x), col='black',lineheigth=2)
},grid = TRUE,main = "PhosphateST177")
#include slope of regresion in the data frame
coresVOL.AREA[1:15,8] <- ST177Phosphate.regression$V1 
names(coresVOL.AREA)[8]<-"slope.Phosphate"

#ST177 Silicate
ST177Silicate <- nutrients[nutrients$STATION== "ST177" & nutrients$TYPE== "Silicate",]
ST177Silicate.regression <- ddply(ST177Silicate, .(REPLICATE), regressionFuncion.nutrients)
#plot each core lineal regresion
library(lattice)
xyplot(Nutrients ~ TIME | REPLICATE,data =ST177Silicate,par.settings = list(strip.background=list(col="lightgrey")),panel = function(x, y,col="white", ...) {panel.xyplot(x, y, col="black",...)
  panel.abline(lm(y~x), col='black',lineheigth=2)
},grid = TRUE,main = "SilicateST177")
#include slope of regresion in the data frame
coresVOL.AREA[1:15,9] <- ST177Silicate.regression$V1 
names(coresVOL.AREA)[9]<-"slope.Silicate"

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ST435nutrients---------
#ST435 Nitrite
ST435nitrite <- nutrients[nutrients$STATION== "ST435" & nutrients$TYPE== "Nitrite",]
ST435nitrite.regression <- ddply(ST435nitrite, .(REPLICATE), regressionFuncion.nutrients)
#plot each core lineal regresion
library(lattice)
xyplot(Nutrients ~ TIME | REPLICATE,data =ST435nitrite,par.settings = list(strip.background=list(col="lightgrey")),panel = function(x, y,col="white", ...) {panel.xyplot(x, y, col="black",...)
  panel.abline(lm(y~x), col='black',lineheigth=2)
},grid = TRUE,main = "nitriteST435")
#include slope of regresion in the data frame
coresVOL.AREA[16:30,6] <- ST435nitrite.regression$V1 

#ST435 Nitrate
ST435Nitrate <- nutrients[nutrients$STATION== "ST435" & nutrients$TYPE== "Nitrate",]
ST435Nitrate.regression <- ddply(ST435Nitrate, .(REPLICATE), regressionFuncion.nutrients)
#plot each core lineal regresion
library(lattice)
xyplot(Nutrients ~ TIME | REPLICATE,data =ST435Nitrate,par.settings = list(strip.background=list(col="lightgrey")),panel = function(x, y,col="white", ...) {panel.xyplot(x, y, col="black",...)
  panel.abline(lm(y~x), col='black',lineheigth=2)
},grid = TRUE,main = "NitrateST435")
#include slope of regresion in the data frame
coresVOL.AREA[16:30,7] <- ST435Nitrate.regression$V1 

#ST435 Phosphate
ST435Phosphate <- nutrients[nutrients$STATION== "ST435" & nutrients$TYPE== "Phosphate",]
ST435Phosphate.regression <- ddply(ST435Phosphate, .(REPLICATE), regressionFuncion.nutrients)
#plot each core lineal regresion
library(lattice)
xyplot(Nutrients ~ TIME | REPLICATE,data =ST435Phosphate,par.settings = list(strip.background=list(col="lightgrey")),panel = function(x, y,col="white", ...) {panel.xyplot(x, y, col="black",...)
  panel.abline(lm(y~x), col='black',lineheigth=2)
},grid = TRUE,main = "PhosphateST435")
#include slope of regresion in the data frame
coresVOL.AREA[16:30,8] <- ST435Phosphate.regression$V1 

#ST435 Silicate
ST435Silicate <- nutrients[nutrients$STATION== "ST435" & nutrients$TYPE== "Silicate",]
ST435Silicate.regression <- ddply(ST435Silicate, .(REPLICATE), regressionFuncion.nutrients)
#plot each core lineal regresion
library(lattice)
xyplot(Nutrients ~ TIME | REPLICATE,data =ST435Silicate,par.settings = list(strip.background=list(col="lightgrey")),panel = function(x, y,col="white", ...) {panel.xyplot(x, y, col="black",...)
  panel.abline(lm(y~x), col='black',lineheigth=2)
},grid = TRUE,main = "SilicateST435")
#include slope of regresion in the data frame
coresVOL.AREA[16:30,9] <- ST435Silicate.regression$V1 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#CALCULATIONS of flux in (umol.m-2.d-1) using formula flux= slope * vol(m3)/ area (m2)
names(coresVOL.AREA)
coresVOL.AREA$nitriteFLUX <- ((coresVOL.AREA$slope.nitrite * coresVOL.AREA$VOL.m3)/ coresVOL.AREA$AREA.m2)/4
coresVOL.AREA$NitrateFLUX <- ((coresVOL.AREA$slope.Nitrate * coresVOL.AREA$VOL.m3)/ coresVOL.AREA$AREA.m2)/4
coresVOL.AREA$PhosphateFLUX <- ((coresVOL.AREA$slope.Phosphate * coresVOL.AREA$VOL.m3)/ coresVOL.AREA$AREA.m2)/4
coresVOL.AREA$SilicateFLUX <- ((coresVOL.AREA$slope.Silicate * coresVOL.AREA$VOL.m3)/ coresVOL.AREA$AREA.m2)/4

#plots NUTRIENTS---------
benthicfluxes <- read.csv(file.path(DATA,"benthicflux2.csv"))
library(sciplot)
par(mfrow = c(2, 2),     # 2x2 layout
    oma = c(1, 1, 1, 0), # c(bottom, left, top, right) margenes fuera del grafico
    mar = c(1, 3.5, 1, 1.5), # c(bottom, left, top, right) separacion entre plots
    mgp = c(2.6, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = T)            # allow content to protrude into outer margin (and beyond)
bargraph.CI(STATION,nitriteFLUX,TREATMENT, data = coresVOL.AREA,ylab="Nitrite (umol m–2 d–1)",xaxt="n",ylim = c(-0.0004,0.0003),yaxt="n")
axis(2,las=1,tck=-.01,col = "gray",cex.axis=0.6)
segments(0.7, 0, 8, 0, col= 'gray')
box(lty = 'solid', col = 'gray')
abline(v=4.5,lty = 2)
bargraph.CI(STATION,NitrateFLUX,TREATMENT, data = coresVOL.AREA,ylab="Nitrate (umol m–2 d–1)",xaxt="n",ylim = c(-0.11,0.04),las=1,,yaxt="n")
axis(2,las=1,tck=-.01,col = "gray",cex.axis=0.6)
segments(0.7, 0, 8, 0, col= 'gray')
box(lty = 'solid', col = 'gray')
abline(v=4.5,lty = 2)
bargraph.CI(STATION,PhosphateFLUX,TREATMENT, data = coresVOL.AREA,ylab="Phosphate (umol m–2 d–1)",ylim = c(-0.009,0.004),las=1,,yaxt="n")
axis(2,las=1,tck=-.01,col = "gray",cex.axis=0.6)
segments(0.7, 0, 8, 0, col= 'gray')
box(lty = 'solid', col = 'gray')
abline(v=4.5,lty = 2)
bargraph.CI(STATION,SilicateFLUX,TREATMENT, data = coresVOL.AREA,ylab="Silicate (umol m–2 d–1)",ylim = c(0,0.7),las=1,,yaxt="n")
axis(2,las=1,tck=-.01,col = "gray",cex.axis=0.6)
segments(0.7, 0, 8, 0, col= 'gray')
box(lty = 'solid', col = 'gray')
abline(v=4.5,lty = 2)












