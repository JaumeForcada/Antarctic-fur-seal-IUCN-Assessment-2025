################################################################################
### FIGURE 6, QUANTITATIVE ASSESSMENT OF THE ANTARCTIC FUR SEAL - MAIVIKEN PLOTS
##------------------------------------------------------------------------------
### LOAD PACKAGES
library(mgcv)
library(MASS)
library(nlme)
library(shape)
### ADDITIONAL FUNCTIONS
source("code/ffs.obs.R")   # FUNCTION TO FIT LOGISTIC MODEL TO PUP PRODUCTION
source("code/addTrans.R")  # COLOUR TRANSPARENCY FILTER
source("code/glbpf.R")
###-----------------------------------------------------------------------------
##  1   DATA
##  1.1 SSB DAILY NEWBORN PUP COUNTS AND OBTAIN DAILY CUMULATIVE SUM
newpups <- read.csv("raw_data/newpups.csv", header = F); colnames(newpups) <- 1985:2025
sumpup <- apply(newpups, 2, cumsum)
##  1.2 OBSERVED SSB PUP MORTALITY FROM 2009 TO 2025 -day X20yy-
dpups <- read.csv("raw_data/SSBmortality.csv",header=T)
##  1.3 OBSERVED SSB PM FEMALE COUNTS
sfc <- read.csv("raw_data/SSB_total_female_counts.csv",header=T)
##  1.4 OBSERVED MVK COUNTS -year day date males females pups juv-
mmfpjc <- read.csv("raw_data/maiviken.csv", header = T)
##
##  2.  PROCESSED DATA FOR TOP TWO FIGURE PANELS
##  2.1 SET VECTORS OF COUNT YEARS, AND MATRICES OF OBSERVED FEMALE AND PUP COUNTS AT BOTH LOCATIONS
cYears <- 2009:2025
PEAKS <- as.numeric(rep(NA,length(cYears)))
MVKPups <- SSBPups <- MVKFems <- SSBFems <- matrix(NA,nrow=length(cYears),ncol=2)
##  2.2 LOOP OVER YEARS
for(i in 1:length(cYears)){
  # PEAK PUPPING AT SSB
  PEAKS[i] <- as.vector(round(ffs.obs(data.frame(day=1:70,pups=sumpup[,paste(cYears[i])]),plots=F)$coefficients[2]))
  # OBSERVED PUP COUNTS AT SSB AND MVK AT PEAK AND PEAK+7 days
  SSBPups[i,1] <- sumpup[PEAKS[i],paste(cYears[i])]-dpups[PEAKS[i],paste("X",cYears[i],sep="")]
  SSBPups[i,2] <- sumpup[PEAKS[i]+7,paste(cYears[i])]-dpups[PEAKS[i]+7,paste("X",cYears[i],sep="")]
  SSBFems[i,] <- sfc[c(PEAKS[i],PEAKS[i]+7),paste("X",cYears[i],sep="")]
  dat <- subset(mmfpjc,year==cYears[i])
  MVKFems[i,1] <- dat$females[which.min(abs(dat$day-PEAKS[i]))]   
  MVKFems[i,2] <- dat$females[which.min(abs(dat$day-(PEAKS[i]+7)))]
  MVKPups[i,1] <- dat$pups[which.min(abs(dat$day-PEAKS[i]))]   
  MVKPups[i,2] <- dat$pups[which.min(abs(dat$day-(PEAKS[i]+7)))]
}
rm(dat)
##
## 3.  OBTAIN SCALING FACTORS DUE TO HAUL-OUT DIFFERENCES BETWEEN MVK AND SSB
##     - USES FREQUENTIST POISSON GAMs AS THEY PROVIDE SIMILAR RESULTS AS jagam BUT MUCH FASTER
##     -> PREDICT PEAK PUPPING DAYS
newpups <- read.csv("raw_data/newpups.csv",header=F)
colnames(newpups) <- 1985:2025
sumpup <- apply(newpups,2,cumsum)   # CUMULATIVE SUM
## FUNCTION ffs.obs FITS AN NLME (LOGISTIC) MODEL TO OBTAIN THE PEAK PUPPING DAY
source("code/ffs.obs.R")
## YEARS WITH COMPLETE BIRD ISLAND SURVEYS
syears2 <- 2010:2025
pks2 <- ses2 <- numeric(length(syears2))
for(i in 1:length(syears2)){
  pks2[i] <- as.vector(ffs.obs(data.frame(day = 1:70, pups = sumpup[, paste(syears2[i])]), plots = F)$coefficients[2])
  ses2[i] <- sqrt(ffs.obs(data.frame(day = 1:70, pups = sumpup[, paste(syears2[i])]), plots = F)$varBeta[2, 2])
}  
## PEAK PUPPING DATES ARE:
pksS <- round(pks2)
pksM <- round(pks2)+6   # FOR CORRELATED MAIVIKEN COUNTS
##
##  SSB FEMALE COUNTS & MAIVIKEN COUNTS & BI-SSB ABUNDANCE ESTIMATES
Sfemc <- read.csv("raw_data/SSB_total_female_counts.csv",header=T)
Mfemc <- read.csv("raw_data/maiviken.csv", header=T)
preds <- readRDS("processed_data/AFSFemaleSimulations200125.Rds")
#
pksM <- pksS; pksM[1] <- 47
MVKAf <- matrix(NA, 15000, length(syears2))
set.seed(666)
for(i in 1:length(syears2)){
  # SSB
  Sdat <- data.frame(y=Sfemc[1:70,paste("X",syears2[i],sep="")],x0=1:70)
  bS <- gam(y ~ s(x0),family=poisson(link="log"),data=Sdat,method="REML")
  XpS <- predict(bS,newdata=data.frame(x0=1:70),type="lpmatrix")
  brS <- mvrnorm(n=15000,coef(bS),bS$Vp)
  bSsim <- exp(XpS %*% t(brS))
  zSsim <- apply(bSsim,2,FUN=function(x) (x-min(x))/(max(x)-min(x)))
  # MVK
  Mdat <- data.frame(y=Mfemc$females[Mfemc$year==syears2[i]],x0=Mfemc$day[Mfemc$year==syears2[i]])
  bM <- gam(y ~ s(x0),family=poisson(link="log"),data=Mdat,method="REML")
  XpM <- predict(bM,newdata=data.frame(x0=1:70),type="lpmatrix")
  brM <- mvrnorm(n=15000,coef(bM),bM$Vp)
  bMsim <- exp(XpM %*% t(brM))
  zMsim <- apply(bMsim,2,FUN=function(x) (x-min(x))/(max(x)-min(x)))
  ##
  dNsim <- zMsim - zSsim   # zM > zS 
  dPsim <- zSsim - zMsim   # zM < zS 
  fkNsim <- bMsim*dNsim/zMsim
  fkPsim <- bMsim*dPsim/zMsim
  fKNsim <- 1-dNsim/zMsim
  fKPsim <- 1+dPsim/zMsim
   MVKAf[,i] <- round(bMsim[pksM[i],] * fKNsim[pksM[i],] * 
    ((preds[,which(dimnames(preds)[[2]]==paste("NB[",i+9,"]",sep=""))]+preds[,which(dimnames(preds)[[2]]==paste("NN[",i+9,"]",sep=""))])/ 
       bSsim[pksS[i],]))
  print(c(i,mean(fKNsim[pksS[i],]),mean(fKPsim[pksS[i],])))
}
## FIT A GAM TREND AND A LOW-PASS GAUSSIAN FILTER TO REMOVE HIGH FREQUENCY COMPONENTS
# gMVKAf <- lpMVKAf <- matrix(NA, 15000, length(syears2))
#for(i in 1:15000){
#  dat <- data.frame(Nf=MVKAf[i,], year=syears2)
#  tgam <- gam(Nf ~ s(year,fx=T,k=3),family=poisson(link="log"),data=dat,method="REML")
#  gMVKAf[i,] <- predict(tgam, type="response")
#  lpMVKAf[i,] <- glbpf(dat$Nf,8)
#  print(i)
#}
#saveRDS(lpMVKAf,"processed_data/lpMVKAf.RDS")
## -> UNCOMMENT WHEN RUN; BRING SAVED PROCESS_DATA
lpMVKAf <- readRDS("processed_data/lpMVKAf.RDS")
##
##   4. GENERATE PLOTS
pdf("documentation/Figure6.pdf",
    family="Helvetica",width=8.27,height=11.69)
par(omi = c(3.5,1,1,0.9))       ## OUTER MARGINS IN INCHES - c(bottom, left, top, right)
layout(matrix(1:3,ncol=1,byrow=T),respect=F)  # PLOT LAYOUT
par(mai = c(0,1,0,0.6))
#
## 1ST plot - PUPS
plot(cYears[2:17],MVKPups[2:17,2],type="n",pch=19,col=2,axes=F,xlab="",ylab="Peak pup count",ylim=c(200,550))
## gam fit
dat <- data.frame(y=MVKPups[3:17,2],x0=cYears[3:17])
b1 <- gam(y~s(x0,fx=T,k=4),family=poisson(link="log"),data=dat,method="REML")
# b1 <- gam(y~s(x0),family=poisson(link="log"),data=dat,method="REML")
b1pred <- predict(b1,type="response",newdata=data.frame(x0=seq(2011,2025,by=0.1)),se.fit=T)
polygon(c(seq(2011,2025,by=0.1),rev(seq(2011,2025,by=0.1))),c(b1pred$fit+1.96*b1pred$se.fit,rev(b1pred$fit-1.96*b1pred$se.fit)),col=addTrans("blue",80),border=NA)
lines(seq(2011,2025,by=0.1),b1pred$fit,col="blue",lwd=1.3)
points(cYears[2:17],MVKPups[2:17,2],pch=19,col="blue")
## rlm model
modr <- MASS::rlm(MVKPups[3:17,2] ~ cYears[3:17],method="MM")
lines(cYears[3:17],predict(modr),lty=2,col="blue")
box()
axis(1,at=seq(2010,2025,by=5),labels=F,tck=-0.025,cex.axis=.9)
axis(1,at=seq(2010,2025,by=1),labels=F,tck=-0.0125,cex.axis=.9)
axis(2,at=pretty(MVKPups[2:17,2]),labels=T,tck=-0.025,cex.axis=.9,las=2,col.axis="blue",hadj=0.77) # ,col.ticks="red"
# mtext("Peak pup count",side=2,line=2.5,cex=.95,col="blue")
Arrows(2010,350,2010,250,code=2, arr.type="triangle",arr.adj=1,col="darkgrey",arr.lwd=1) #arr.length=0.1,arr.width=0.3/4,
#
## ADD SSB
par(new=T)
dat <- data.frame(y=SSBPups[3:17,1],x0=cYears[3:17])
b0 <- gam(y~s(x0),family=poisson(link="log"),data=dat,method="REML")
b0pred <- predict(b0,type="response",newdata=data.frame(x0=seq(2011,2025,by=0.1)),se.fit=T)
plot(cYears[2:17],SSBPups[2:17,1],type="n",axes=F,xlab="",ylab="",
     ylim=c(min(SSBPups[3:17,1]),max(b0pred$fit+1.96*b0pred$se.fit)))
polygon(c(seq(2011,2025,by=0.1),rev(seq(2011,2025,by=0.1))),c(b0pred$fit+1.96*b0pred$se.fit,rev(b0pred$fit-1.96*b0pred$se.fit)),col=addTrans("red",80),border=NA)
lines(seq(2011,2025,by=0.1),b0pred$fit,col="red",lwd=1.3)
points(cYears[2:17]+0.1,SSBPups[2:17,1],pch=19,col=addTrans("red",100))
points(cYears[2:17]+0.1,SSBPups[2:17,1],pch=1,col="red")
## rlm model
modr <- MASS::rlm(y ~ x0,data=dat,method="MM")
lines(dat$x0,predict(modr),lty=2,col="red")
axis(4,at=pretty(c(min(SSBPups[3:17,1]),max(b0pred$fit+1.96*b0pred$se.fit))),labels=T,tck=-0.025,cex.axis=.9,las=2,col.axis="red",hadj=0.2)
#mtext("Peak pup count at SSB",side=4,line=2.5,cex=.95,col="red")
## LEGEND
legend("topright",title=" ",title.adj=0.2,legend=c("Bird Island","Maiviken"),
       cex=0.8,bty="n",ncol=1,col=c("red","blue"),pch=c(21,19),pt.bg=c(addTrans("red",100),NA),pt.cex=rep(1.1,2),lty=c(1,1))
par(new=F)
#
## 2nd PLOT - FEMALES
#
par(mai = c(0,1,0,0.6))
##
plot(cYears[2:17],MVKFems[2:17,2],type="n",pch=19,col=2,axes=F,xlab="",ylab="Peak female count",ylim=range(MVKFems[2:17,2]))
## gam fit
dat <- data.frame(y=MVKFems[3:17,2],x0=cYears[3:17])
b1 <- gam(y~s(x0,fx=T,k=5),family=poisson(link="log"),data=dat,method="REML")
#b1 <- gam(y~s(x0),family=poisson(link="log"),data=dat,method="REML")
b1pred <- predict(b1,type="response",newdata=data.frame(x0=seq(2011,2025,by=0.1)),se.fit=T)
polygon(c(seq(2011,2025,by=0.1),rev(seq(2011,2025,by=0.1))),c(b1pred$fit+1.96*b1pred$se.fit,rev(b1pred$fit-1.96*b1pred$se.fit)),col=addTrans("blue",80),border=NA)
lines(seq(2011,2025,by=0.1),b1pred$fit,col="blue",lwd=1.3)
points(cYears[2:17],MVKFems[2:17,2],pch=19,col="blue")
## rlm model
modr <- MASS::rlm(y ~ x0,data=dat,method="MM")
lines(dat$x0,predict(modr),lty=2,col="blue")
box()
axis(1,at=seq(2010,2025,by=5),labels=F,tck=-0.025,cex.axis=.9)
axis(1,at=seq(2010,2025,by=1),labels=F,tck=-0.0125,cex.axis=.9)
axis(2,at=pretty(MVKFems[2:17,2]),labels=T,tck=-0.025,cex.axis=.9,las=2,col.axis="blue",hadj=0.77) # ,col.ticks="red"
#mtext("Peak female count",side=2,line=2.5,cex=.95,col="blue")
Arrows(2010,500,2010,400,code=2, arr.type="triangle",arr.adj=1,col="darkgrey",arr.lwd=1) #arr.length=0.1,arr.width=0.3/4,
#
## ADD SSB
par(new=T)
dat <- data.frame(y=SSBFems[3:17,1],x0=cYears[3:17])
b0 <- gam(y~s(x0),family=poisson(link="log"),data=dat,method="REML")
b0pred <- predict(b0,type="response",newdata=data.frame(x0=seq(2011,2025,by=0.1)),se.fit=T)
plot(cYears[2:17],SSBFems[2:17,1],type="n",axes=F,xlab="",ylab="",
     ylim=c(min(SSBFems[2:17,1]),max(b0pred$fit+1.96*b0pred$se.fit)))
polygon(c(seq(2011,2025,by=0.1),rev(seq(2011,2025,by=0.1))),c(b0pred$fit+1.96*b0pred$se.fit,rev(b0pred$fit-1.96*b0pred$se.fit)),col=addTrans("red",80),border=NA)
lines(seq(2011,2025,by=0.1),b0pred$fit,col="red",lwd=1.3)
points(cYears[2:17]+0.1,SSBFems[2:17,1],pch=19,col=addTrans("red",100))
points(cYears[2:17]+0.1,SSBFems[2:17,1],pch=1,col="red")
## rlm model
modr <- MASS::rlm(y ~ x0,data=dat,method="MM")
lines(dat$x0,predict(modr),lty=2,col="red")
axis(4,at=pretty(c(min(SSBFems[2:17,1]),max(b0pred$fit+1.96*b0pred$se.fit))),labels=T,tck=-0.025,cex.axis=.9,las=2,col.axis="red",hadj=0.2)
#mtext("Female count at SSB",side=4,line=2.5,cex=.95,col="red")
## LEGEND
#legend("topright",title="Mature females",title.adj=0.2,legend=c("Bird Island","Maiviken"),
#       cex=0.8,bty="n",ncol=1,col=c("red","blue"),pch=rep(19,2),pt.cex=rep(1.1,2),lty=c(1,1))
par(new=F)
##
##  THIRD PLOT - MAIVIKEN MATURE FEMALE ABUNDANCE
#
par(mai = c(0,1,0,0.6))
##
plot(syears2,colMeans(MVKAf),type="n",ylim=c(600,1700),ylab="Maiviken mature female abundance",xlab="",axes=F)
polygon(c(syears2,rev(syears2)), c(apply(lpMVKAf,2,quantile,0.975),rev(apply(lpMVKAf,2,quantile,0.025))),col=addTrans("magenta",80),border=NA)
lines(syears2,colMeans(lpMVKAf),col="magenta",lwd=1.1)
lines(syears2,predict(lm(colMeans(MVKAf) ~ syears2)),lty=2,col="black")
for(i in 1:length(syears2))
  lines(rep(syears2[i],2),apply(MVKAf,2,quantile,c(0.025,0.975))[,i],col="black")
points(syears2,colMeans(MVKAf),type ="p",pch=21,col="black",bg="magenta",cex=0.85)
box()
axis(1,at=seq(2010,2025,by=5),labels=T,tck=-0.035,cex.axis=.9)
axis(1,at=seq(2010,2025,by=1),labels=F,tck=-0.0125,cex.axis=.9)
axis(2,at=seq(400,2000,by=200),labels=T,tck=-0.025,las=2,hadj=0.825,cex.axis=.85,col.axis="magenta3",col.ticks="magenta3")
axis(2,at=seq(400,2000,by=100),labels=F,tck=-0.0125,las=2,col.axis="magenta3",col.ticks="magenta3")
##
dev.off()



##



