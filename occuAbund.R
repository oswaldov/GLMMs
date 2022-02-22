##Load packages
library(stats)
library(statmod)
library(dplyr)
library(tidyr)
library(rpart)
library(rpart.plot)
library(quantmod)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(scales)
library(glmmTMB)
library(DHARMa)
library(modEvA)
library(performance)
library(MASS)

##Load data

mosqcountcx <-read.csv("Cxmosquito.csv", header=TRUE)
dim(mosqcountcx)
names(mosqcountcx)
head(mosqcountcx)

##Exploratory analysis
mosqcountcx$year <- as.factor(mosqcountcx$year)
mosqcountcx$month <- as.factor(mosqcountcx$month)

hist(mosqcountcx$count, breaks = 20, xlab = "", ylab = "Frequency", main = "")  ## check data distribution

##plot total count of mosquitoes by year and by location

ggplot(mosqcountcx,aes(x = year, y = count, fill = trap_type)) + 
  geom_bar(stat = "identity", width = 0.6) +
  facet_wrap(~site) +
  xlab("Year of sampling") + ylab("Number of Culex mosquitoes")

##plot total count of mosquitoes by year and by location

ggplot(mosqcountcx,aes(x = month, y = count, fill = trap_type)) + 
  geom_bar(stat = "identity", width = 0.6) +
  facet_wrap(~site) +
  xlab("Month of sampling") + ylab("Number of Culex mosquitoes")


##Split data set into training and testing data subsets

split1 <- sample(c(rep(0, 0.7 * nrow(mosqcountcx)), rep(1, 0.3 * nrow(mosqcountcx))))
split1

table(split1)

traindat <- mosqcountcx[split1==0, ]  ## train dataset 

testdat <- mosqcountcx[split1==1, ] ## test dataset


## CART ANALYSIS

##Check for NAs

mosqcount <- traindat[complete.cases(traindat),]

## Preliminary visualizations

par(mfrow=c(2,2))
plot(mosqcount$site, mosqcount$Cxpa)
plot(mosqcount$season, mosqcount$Cxpa)
plot(mosqcount$rainfall_season, mosqcount$Cxpa)
plot(mosqcount$elevation, mosqcount$Cxpa)

##correlation check: geographic
locat <- c(13,14,15)
summary(mosqcount[,locat])
pairs(mosqcount[,locat])

##correlation check: precip variables
precip <- c(13,16:19)
summary(mosqcount[,precip])
pairs(mosqcount[,precip])

## correlation check: temperature variables 
temp<- c(13,14,20:23)
summary(mosqcount[,temp])
pairs(mosqcount[,temp])


##Select variables of interest

preds<-c("Cxpa")  ## response variable
geo <- c("site","trap_type","year","month","rainfall_season") ## meteorological and geografical variables
loc <- c("elevation","distance") ## location variables
precip <- c("precip","preciplag1","preciplag2") ## precipitation variables
temp <- c("tmean","tmax","tmin")
tempmean <- c("tmean","tmeanlag1","tmeanlag2") ## temperature variables


all<- c(preds, geo, loc, precip, tempmean)


dataCX <- mosqcount[,all]
head(dataCX)


## First fit forcing the tree to have lots of branches so we can
## examine it and figure out where to trim
CP=0.0005  ##complexity parameter, lack of fit 
MS=50  ## the minimum number of observations in a node in order for split to be attempted 
cnt<-rpart.control(minsplit=MS, cp=CP, xval=100)

f.null<-rpart(Cxpa ~  . , data=dataCX, method="class", control=cnt)
plotcp(f.null) ## use this to decide on a cp/size for trimming the tree
printcp(f.null)
graphics.off()


x11()
par(mfrow=c(1,1))
#plot(f.null)

rpart.plot(f.null, type = 4, extra = 101, branch.lty = 3)
printcp(f.null)
plotcp(f.null)


## now trimming based on the above: plot(f)
cnt<-rpart.control(minsplit=MS, cp=0.0025, xval=100)
f<-rpart(Cxpa ~  . , data=dataCX, method="class", control=cnt)
plotcp(f)
graphics.off()


x11()
par(mfrow=c(1,1))

rpart.plot(f, type = 4, extra = 101, branch.lty = 3, space = 1.1)

## here's what rpart gives for its trimmed tree if I don't include a
## controler
par(mfrow=c(1,1))
f.d<-rpart(Cxpa ~  . , method="class", data=dataCX)
printcp(f.d)
plotcp(f.d)

par(mfrow=c(1,1))
rpart.plot(f.d,type=4, extra=101, box.palette = "GnBu",branch.lty=3,shadow.col="gray")


## GLM analysis 

## squared variables

dataCX$tmeansq<- dataCX$tmean^2


##models

m.null<-glm(Cxpa ~ 1, family="binomial", data=dataCX)
m.null
m.full<-glm(Cxpa ~ . , family="binomial", data=dataCX)
m.full


modelcx<- step(m.null, scope=formula(m.full), direction="forward", criterion = "BIC")
summary(modelcx)

## Check the residuals: quantile - quantile plots, draes teh correlation between a given sample and the normal distribution
qr2<- qresiduals(modelcx)
length(qr2)
summary(qr2)


par(mfrow=c(1,1))
qqnorm(qr2, ylim = c(-6,6), xlim=c(-5,5), main = "Normal Q-Q Plot", las=1, bty="o"); qqline(qr2)


## lets look at marginal predictions, based on particular predictors

## elevation

## lets look at marginal predictions, based on particular predictors

summary(dataCX$elevation)

o<-order(dataCX$elevation)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")
## plot based off predictions from the glm for data pre-2004
plot(dataCX$elevation[o], dataCX$Cxpa, xlab="Elevation (m)",
     ylab= ey, cex.lab = 1.1, ylim=c(0,1), xlim=c(0,2000),main="GLM fit - Elevation")
fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, elevation=dataCX$elevation[o])
y.m2<-predict(loess(fit~elevation, data=fit.dat.m, span = 1), data.frame(elevation=seq(0, 2000, by=10)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0,2000, by=10), y.m2$fit, col=2, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 


## lets try it with the tree
f.pred<-predict(f, type="vector")
table(f.pred)


#f.pred[f.pred = 1] <- 0 
f.pred1 <- recode(f.pred, '1' = 0, '2' = 1)
table(f.pred1)

plot(jitter(dataCX$elevation[o]), dataCX$Cxpa, ylim=c(0,1),xlim=c(0,2000),
     xlab="Elevation (m)", ylab= ey, cex.lab = 1.1,
     main="CART fit - Elevation")
fittedf<- as.numeric(f.pred1[o])
fit.dat.f<-data.frame(fit=fittedf, elevation=dataCX$elevation[o])
#fit.dat.f
#plot(fit.dat.f$elevation, fit.dat.f$fit)
y.f2<-predict(loess(fit~elevation, data=fit.dat.f, span=1), data.frame(elevation=seq(0, 2000, by=10)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2$fit

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 2000, by=10), y.f2$fit, col=4, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))

## graph legend
legend(x=1900, y=1.2, legend="A", xpd = NA, bty="n", cex = 1.7)



##Distance


o<-order(dataCX$distance)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataCX$distance[o], dataCX$Cxpa, xlab="Distance - anthropogenic features (m)",
     ylab= ey, cex.lab = 1.1, ylim=c(0,1), xlim=c(0,2500),main="GLM fit - Distance")
fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, distance=dataCX$distance[o])
y.m2<-predict(loess(fit~distance, data=fit.dat.m, span=1), data.frame(distance=seq(0, 2500, by=100)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit

##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0,2500, by=100), y.m2$fit, col=2, lwd=3)
legend(200, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$distance[o]), dataCX$Cxpa,ylim=c(0,1),xlim=c(0,2500),
     xlab="Distance - anthropogenic features (m)", ylab= ey, cex.lab = 1.1, 
     main="CART fit - Distance")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, distance=dataCX$distance[o])
fit.dat.f
y.f2<-predict(loess(fit~distance, data=fit.dat.f, span=1), data.frame(distance=seq(0, 2500, by=100)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2$fit
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 2500, by=100), y.f2$fit, col=4, lwd=3)
legend(200, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))

##graph legend
legend(x=2300, y=1.2, legend="B", xpd = NA, bty="n", cex = 1.7)



##Precipitation

o<-order(dataCX$precip)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")
## plot based off predictions from the glm for data pre-2004
plot(dataCX$precip[o], dataCX$Cxpa, xlab="Precipitation (mm)",
     ylab= ey, cex.lab = 1.1, ylim=c(0,1), xlim=c(0,1200),main="GLM fit - Precipitation")
fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, precip=dataCX$precip[o])
y.m2<-predict(loess(fit~precip, data=fit.dat.m, span=1), data.frame(precip=seq(0, 1200, by=1)), se=TRUE)
y.m2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0,1200, by=1), y.m2$fit, col=2, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$precip[o]), dataCX$Cxpa,ylim=c(0,1),xlim=c(0,1200),
     xlab="Precipitation (mm)", ylab = ey, cex.lab = 1.1,
     main="CART fit - Precipitation")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, precip=dataCX$precip[o])
fit.dat.f
y.f2<-predict(loess(fit~precip, data=fit.dat.f, span=1), data.frame(precip=seq(0, 1200, by=1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
y.f2


##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1200, by=1), y.f2$fit, col=4, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


##graph legend
legend(x=1180, y=1.2, legend="E", xpd = NA, bty="n", cex = 1.7)



##Preciplag1

o<-order(dataCX$preciplag1)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataCX$preciplag1[o],dataCX$Cxpa[o], xlab="Precipitation lag1 (mm)",
     ylab= ey, cex.lab= 1.1, xlim=c(0,1200),ylim=c(0,1),
     main="GLM fit - Precipitation lag1")

fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, preciplag1=dataCX$preciplag1[o])
y.m2<-predict(loess(fit~preciplag1, data=fit.dat.m, span=1), data.frame(preciplag1=seq(0, 1200, by=10)), se=TRUE)
y.m2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1200, by=10), y.m2$fit, col=2, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$preciplag1[o]), dataCX$Cxpa[o],
     xlab="Precipitation lag1 (mm)", ylab= ey, cex.lab=1.1, 
     xlim=c(0,1200), ylim=c(0,1),main="CART fit - Precipitation lag1")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, preciplag1=dataCX$preciplag1[o])
fit.dat.f
y.f2<-predict(loess(fit~preciplag1, data=fit.dat.f, span=1), data.frame(preciplag1=seq(0, 1200, by=10)), se=TRUE)
y.f2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1200, by=10), y.f2$fit, col=4, lwd=3)
##lines(jitter(dat$bio09[o]/10), y.f, col=2, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=1200, y=1.2, legend="F", xpd = NA, bty="n", cex = 1.6)



## Precipitation lag2

o<-order(dataCX$preciplag2)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataCX$preciplag2[o],dataCX$Cxpa[o], xlab="Precipitation lag2 (mm)",
     ylab= ey, cex.lab= 1.1, xlim=c(0,1200),ylim=c(0,1),
     main="GLM fit - Precipitation lag2")

fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, preciplag2=dataCX$preciplag2[o])
y.m2<-predict(loess(fit~preciplag2, data=fit.dat.m, span=1), data.frame(preciplag2=seq(0, 1200, by=10)), se=TRUE)
y.m2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1200, by=10), y.m2$fit, col=2, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$preciplag2[o]), dataCX$Cxpa[o],
     xlab="Precipitation lag2 (mm)", ylab= ey, cex.lab=1.1, 
     xlim=c(0,1200), ylim=c(0,1),main="CART fit - Precipitation lag2")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, preciplag2=dataCX$preciplag2[o])
fit.dat.f
y.f2<-predict(loess(fit~preciplag2, data=fit.dat.f, span=1), data.frame(preciplag2=seq(0, 1200, by=10)), se=TRUE)
y.f2
##lines(dat$bio09[o]/10, y.sym, col=2, lwd=3)
lines(seq(0, 1200, by=10), y.f2$fit, col=4, lwd=3)
##lines(jitter(dat$bio09[o]/10), y.f, col=2, lwd=3)
legend(100, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))



legend(x=1200, y=1.2, legend="G", xpd = NA, bty="n", cex = 1.6)


## TEMPERATURE MEAN

o<-order(dataCX$tmean)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataCX$tmean[o], dataCX$Cxpa[o], xlab="Temperature (°C)",
     ylab= ey, cex.lab=1.1, xlim=c(10,28),ylim=c(0,1),
     main="GLM fit - Temperature")

fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, tmean=dataCX$tmean[o])
y.m2<-predict(loess(fit~tmean, data=fit.dat.m, span=1), data.frame(tmean=seq(10,28, by=0.1)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit
##lines(dat$lGDP, y.sym, col=2, lwd=3)
lines(seq(10, 28, by=0.1), y.m2$fit, col=2, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$tmean[o]), dataCX$Cxpa[o],
     xlab="Temperature (°C)", ylab= ey, cex.lab=1.1, 
     xlim=c(10,28), ylim=c(0,1),main="CART fit - Temperature")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, tmean=dataCX$tmean[o])
fit.dat.f
y.f2<-predict(loess(fit~tmean, data=fit.dat.f, span=1), data.frame(tmean=seq(10,28, by=0.1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
##lines
lines(seq(10, 28, by=0.1), y.f2$fit, col=4, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))


legend(x=27, y=1.2, legend="H", xpd = NA, bty="n", cex = 1.7)


## tmean lag1 

o<-order(dataCX$tmeanlag1)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataCX$tmeanlag1[o], dataCX$Cxpa[o], xlab="Temperature lag1(°C)",
     ylab= ey, cex.lab=1.1, xlim=c(10,28),ylim=c(0,1),
     main="GLM fit - Temperature lag1")

fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, tmeanlag1=dataCX$tmeanlag1[o])
y.m2<-predict(loess(fit~tmeanlag1, data=fit.dat.m, span=1), data.frame(tmeanlag1=seq(10,28, by=0.1)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit
##lines(dat$lGDP, y.sym, col=2, lwd=3)
lines(seq(10, 28, by=0.1), y.m2$fit, col=2, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$tmeanlag1[o]), dataCX$Cxpa[o],
     xlab="Temperature lag1(°C)", ylab= ey, cex.lab=1.1, 
     xlim=c(10,28), ylim=c(0,1),main="CART fit - Temperature lag1")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, tmeanlag1=dataCX$tmeanlag1[o])
fit.dat.f
y.f2<-predict(loess(fit~tmeanlag1, data=fit.dat.f, span=1), data.frame(tmeanlag1=seq(10,28, by=0.1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
##lines
lines(seq(10, 28, by=0.1), y.f2$fit, col=4, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))



legend(x=28, y=1.2, legend="I", xpd = NA, bty="n", cex = 1.7)



## tmean lag2

o<-order(dataCX$tmeanlag2)

op <- par(
  mar=c(4,4.5,2,2),
  oma=c(2,2,2,2),
  xpd=TRUE,
  mfrow=c(1,2), bty="n")

## plot based off predictions from the glm for data pre-2004
plot(dataCX$tmeanlag2[o], dataCX$Cxpa[o], xlab="Temperature lag2 (°C)",
     ylab= ey, cex.lab=1.1, xlim=c(10,28),ylim=c(0,1),
     main="GLM fit - Temperature lag2")

fitted<-as.numeric(modelcx$fitted[o])
fitted
fit.dat.m<-data.frame(fit=fitted, tmeanlag2=dataCX$tmeanlag2[o])
y.m2<-predict(loess(fit~tmeanlag2, data=fit.dat.m, span=1), data.frame(tmeanlag2=seq(10,28, by=0.1)), se=TRUE)
y.m2$fit[y.m2$fit<0] <- 0
y.m2$fit
##lines(dat$lGDP, y.sym, col=2, lwd=3)
lines(seq(10, 28, by=0.1), y.m2$fit, col=2, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=1:2, lty=c(0, 1), lwd=c(0,3), pch=c(21, 0)) 



## lets try it with the tree
f.pred<-predict(f, type="vector")
plot(jitter(dataCX$tmeanlag2[o]), dataCX$Cxpa[o],
     xlab="Temperature lag2 (°C)", ylab= ey, cex.lab=1.1, 
     xlim=c(10,28), ylim=c(0,1),main="CART fit - Temperature lag2")
fittedf<- as.numeric(f.pred1[o])
fittedf
fit.dat.f<-data.frame(fit=fittedf, tmeanlag2=dataCX$tmeanlag2[o])
fit.dat.f
y.f2<-predict(loess(fit~tmeanlag2, data=fit.dat.f, span=1), data.frame(tmeanlag2=seq(10,28, by=0.1)), se=TRUE)
y.f2$fit[y.f2$fit<0] <- 0
##lines
lines(seq(10, 28, by=0.1), y.f2$fit, col=4, lwd=3)
legend(10, 0.95, c("fitted values", "ave fitted"),
       col=c(1,4), lty=c(0, 1), lwd=c(0,3), pch=c(21, 0))



legend(x=28, y=1.2, legend="J", xpd = NA, bty="n", cex = 1.7)



## Abundance analysis 


## Negative binomials with dispersion

nbdmd1 <- glmmTMB(count ~ site + (1|year), dispformula = ~month, mosqcount, family = nbinom2)
nbdmd2 <- glmmTMB(count ~ site + month + (1|year), dispformula = ~month, mosqcount, family = nbinom2)

## Zero poisson inflated models 

zipmd1 <- glmmTMB(count ~ site + (1|year), zi=~site, mosqcount, family = poisson)
zipmd2 <- glmmTMB(count ~ site + month + (1|year), zi=~site + month, mosqcount, family = poisson)



## Zero negative binomial inflated models 
zinbmd1 <- glmmTMB(count ~ site + (1|year), zi=~site, mosqcount, family = nbinom2)
zinbmd2 <- glmmTMB(count ~ site + month + (1|year), zi=~site + month, mosqcount, family = nbinom2)
zinbmd10 <- glmmTMB(count ~ site + month + trap_type + rainfall_season + distance + precip +
                     (1|year), zi=~site + month + trap_type + rainfall_season + distance + precip
                      , mosqcount, family = nbinom2)

##Compare models base on AIC values
AIC(nbdmd1, nbdmd2, zipmd1, zipmd2, zinbmd1, zinbmd2, zinbmd10)


## Summary best model, check for residuals and for zero inflation

simulateResiduals(fittedModel = zinbmd10, plot = T)

simulationOutput <- simulateResiduals(fittedModel = zinbmd10)
testZeroInflation(simulationOutput)

##Plot abundance

##Generate new data
newdata0 = newdata = unique(mosqcount[,c("site","month","trap_type","rainfall_season","distance","precip")])

##Set the random effects to zero

X.cond <- model.matrix(lme4::nobars(formula(zinbmd10)[-2]), newdata0)
beta.cond <- fixef(zinbmd10)$cond
pred.cond <- X.cond %*% beta.cond


ziformula = zinbmd10$modelInfo$allForm$ziformula
X.zi <- model.matrix(lme4::nobars(ziformula),newdata0)
beta.zi <- fixef(zinbmd10)$zi
pred.zi <- X.zi %*% beta.zi


##estimates of the linear predictor transformed to teh scale response

pred.ucount <- exp(pred.cond)*(1-plogis(pred.zi))

## Standard errors and CIs

pred.condpar.psim <- mvrnorm(1000, mu=beta.cond, Sigma = vcov(zinbmd10)$cond) 
pred.condpar.psim
pred.cond.psim <- X.cond %*% t(pred.condpar.psim)
pred.zipar.psim <- mvrnorm(1000, mu=beta.zi, Sigma =vcov(zinbmd10)$zi)
pred.zi.psim <- X.zi %*% t(pred.zipar.psim)
pred.ucount.psim <- exp(pred.cond.psim)*(1-plogis(pred.zi.psim))

ci.ucount <- t(apply(pred.ucount.psim,1,quantile,c(0.025,0.975)))
ci.ucount <- data.frame(ci.ucount)
names(ci.ucount) <- c("ucount.low","ucount.high")
pred.ucount <- data.frame(newdata0, pred.ucount, ci.ucount)

## plot abundance

library(plyr)
##observed data
real.count <- ddply(mosquicount, ~site+year+month+trap_type+rainfall_season+distance+precip+tmean+preciplag1, summarize, m=median(count), mu=mean(count))


ggplot(pred.ucount,aes(x=site, y=pred.ucount, colour=trap_type)) + geom_point(shape=1, size=2) +
  geom_errorbar(aes(ymin=ucount.low, ymax=ucount.high))+
  geom_point(data=real.count, aes(x=site, y=m, colour=trap_type), shape=20, size=2)+
  ylab("Predicted mosquito abundace")+
  xlab("Sites")

