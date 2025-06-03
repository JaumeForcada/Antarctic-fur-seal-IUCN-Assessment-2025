##
## LOAD PACKAGES
library(mgcv)
library(jagsUI)
library(MASS)
library(geodist)
##----------------------------------
##  OBTAIN PROPORTIONS OF THE ESTIMATED SSB MATURE FEMALE POPULATION E([Nt+1|Nt]) THAT REPRESENTS
##  A TOTAL DAILY PM FEMALE COUNT AS SSB FOR THE DAYS IN DECEMBER 2021 WITH SURVEYS AT COLONIES ON 
##  SOUTH GEORGIA MAIN LAND
##---------------------------------- 
## SSB FEMALE INTEGRATED POPULATION MODEL ESTIMATES FROM 2001 TO 2025
## FIT GAMM COMPONENT OF THE IPM LIKELIHOOD L_c(n_s|β,b,λ,φ,ζ) 
#  -> POSTERIOR SIMULATIONS FOR MATURE FEMALES IN 2021-22 
preds <- readRDS("processed_data/AFSFemaleSimulations200125.Rds")
NBIs22 <- preds[,which(dimnames(preds)[[2]]=="NBI[22]")]
NNs22 <- preds[,which(dimnames(preds)[[2]]=="NN[22]")]
##----------------------------------
## FIT GAMM TO ESTIMATE PREDICTED SSB FEMALES COUNTS IN 2022
tfemc <- read.csv("raw_data/SSB_total_female_counts.csv", header = T)
dat <- data.frame(y=tfemc[1:70,paste("X",2022,sep="")],x0=1:70)
b0 <- gam(y ~ s(x0),family=poisson(link="log"),data=dat,method="REML")
jags.file <- "test.jags"
## USE jagam TO GENERATE A DATA OBJECT SUITABLE TO FIT GAMs WITH APPROPRIATE PRIORS FOR JAGS
jd <- jagam(y ~ s(x0),family=poisson(link="log"),data=dat,file=jags.file,sp.prior="gamma",diagonalize = TRUE)
##----------------------------------
## BUGS MODEL RUN WITH JAGS
sink("SSBFem.GAMM") 
cat("
  model{
    eta <- X %*% b ## linear predictor
    for(i in 1:n){
      mu[i] <-  exp(eta[i] + eps[i])
      eps[i] ~ dnorm(0, tau)
    } ## expected response
    for(i in 1:n){
      y[i] ~ dpois(mu[i])
    } ## response 
    ## random effects priors
    tau <- pow(sd, -2)
    sd ~ dunif(0, 5)
    ## Parametric effect priors CHECK tau=1/40^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,0.00062)
    }
    ## prior for s(x0)... 
    for(i in c(2:9)){
      b[i] ~ dnorm(0, lambda[1])
    }
    for(i in c(10)){
      b[i] ~ dnorm(0, lambda[2])
    }
    ## smoothing parameter priors CHECK...
    for (i in 1:2) {
      lambda[i] ~ dgamma(.05,.005)
      rho[i] <- log(lambda[i])
    }
  }
  ", fill = TRUE)
sink()
##
## DATA
jags.data <- list(y=jd$jags.data$y,n=jd$jags.data$n,X=jd$jags.data$X)
## INITIAL VALUES
inits <- function()
  list(b=jd$jags.ini$b,lambda=jd$jags.ini$lambda)
## PARAMETRES MONITORED
parameters <- c("b","rho","mu","sd","eps")
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
### RUN BUGS MODEL
sta <- Sys.time()
out <- jags(jags.data,inits,parameters,"SSBFem.GAMM",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta  
# Time difference of 14.65347 secs
##----------------------------------
## SELECT POSTERIOR SIMULATIONS FOR INFERENCE
predFMu <- rbind(out$samples[[1]],out$samples[[2]],out$samples[[3]])
## INFERENCE FOR PREDICTED SSB COUNTS ON DAYS j (from Nov. 1st): 48, 52, 57
predFMu48 <- predFMu[,which(dimnames(predFMu)[[2]]=="mu[48]")]
predFMu52 <- predFMu[,which(dimnames(predFMu)[[2]]=="mu[52]")]
predFMu57 <- predFMu[,which(dimnames(predFMu)[[2]]=="mu[57]")]
## ESTIMATES OF N^M,t_S/E(n_j^S)
NMf.nc48 <- (NBIs22+NNs22)/predFMu48
NMf.nc52 <- (NBIs22+NNs22)/predFMu52
NMf.nc57 <- (NBIs22+NNs22)/predFMu57
##----------------------------------
## OBTAIN SCALING FACTORS DUE TO HAUL-OUT DIFFERENCES BETWEEN MVK AND SSB
## USES FREQUENTIST POISSON GAMs AS THEY PROVIDE SIMILAR RESULTS AS jagam BUT MUCH FASTER
## -> BIRD ISLAND - SSB
Sfemc <- read.csv("raw_data/SSB_total_female_counts.csv",header=T)
Sdat <- data.frame(y=Sfemc[1:70,paste("X",2022,sep="")],x0=1:70)
bS <- gam(y ~ s(x0),family=poisson(link="log"),data=Sdat,method="REML")
## -> MAIVIKEN
Mfemc <- read.csv("raw_data/maiviken.csv", header = T)
Mdat <- data.frame(y=Mfemc$females[Mfemc$year==2022],x0=Mfemc$day[Mfemc$year==2022])
bM <- gam(y ~ s(x0),family=poisson(link="log"),data=Mdat,method="REML")
#-----------------------------------
## UNCERTAINTY FROM POSTERIOR SIMULATIONS
## SSB
XpS <- predict(bS,newdata=data.frame(x0=1:70),type="lpmatrix")
set.seed(666)
brS <- mvrnorm(n=15000,coef(bS),bS$Vp)
bSsim <- exp(XpS %*% t(brS))
zSsim <- apply(bSsim,2,FUN=function(x) (x-min(x))/(max(x)-min(x)))
## MVK
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
## ("delta arrows" IN FIGURE 3)
#  Hound Bay: 18/12/2021
#  Stromness Harbour: 22/12/2021
#  Husvik Harbour: 22/12/2021
#  Maiviken: 27/12/2021
## FEMALES
## DAYS: 48, 52, 57 -> +, +, -
## KAPPA ESTIMATES - SCALE UP OR DOWN A COLONY COUNT DEPENDING ON HAUL OUT CURVES AT MVK WRT SSB
round(c(mean=mean(fKPsim[48,]),quantile(fKPsim[48,],c(0.025,0.975))),2)
# mean  2.5% 97.5% 
# 1.17  1.05  1.30 
round(c(mean=mean(fKPsim[52,]),quantile(fKPsim[52,],c(0.025,0.975))),2)
# mean  2.5% 97.5% 
# 1.12  1.00  1.24 
round(c(mean=mean(fKNsim[57,]),quantile(fKNsim[57,],c(0.025,0.975))),2)
# mean  2.5% 97.5% 
# 0.81  0.73  0.89 
###-----------------------------------------------------------------------------
## SAME ANALYSIS FOR MALES
##
## SSB MALE INTEGRATED POPULATION MODEL ESTIMATES FROM 1995 TO 2025
## FIT GAMM COMPONENT OF THE IPM LIKELIHOOD L_c(n_s|β,b,λ,φ,ζ) 
#  -> POSTERIOR SIMULATIONS FOR MATURE MALES IN 2021-22 
predsM <- readRDS("processed_data/AFSMaleSimulations19952025.Rds")
NTs28 <- predsM[,which(dimnames(predsM)[[2]]==paste("NT[",19,"]",sep=""))]
NNs28 <- predsM[,which(dimnames(predsM)[[2]]==paste("NN[",19,"]",sep=""))]
##----------------------------------
## FIT GAMM TO ESTIMATE PREDICTED SSB AM MALE COUNTS IN 2022
tmc <- read.csv("raw_data/SSB_total_male_counts.csv",header=T)
dat <- data.frame(y=tmc[1:61,paste("X",2022,sep="")],x0=1:61)
b0 <- gam(y ~ s(x0),family=poisson(link="log"),data=dat,method="REML")
jags.file <- "test.jags"
jd <- jagam(y ~ s(x0),family=poisson(link="log"),data=dat,file=jags.file,sp.prior="gamma",diagonalize = TRUE)
##----------------------------------
## BUGS MODEL
sink("SSBMale.GAMM") 
cat("
  model{
    eta <- X %*% b ## linear predictor
    for(i in 1:n){
      mu[i] <-  exp(eta[i] + eps[i])
      eps[i] ~ dnorm(0, tau)
    } ## expected response
    for(i in 1:n){
      y[i] ~ dpois(mu[i])
    } ## response
    ## random effects priors
    tau <- pow(sd, -2)
    sd ~ dunif(0, 5)
    ## Parametric effect priors CHECK tau=1/32^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,0.001)
    }
    ## prior for s(x0)... 
    for(i in c(2:9)){
      b[i] ~ dnorm(0, lambda[1])
    }
    for(i in c(10)){
      b[i] ~ dnorm(0, lambda[2])
    }
    ## smoothing parameter priors CHECK...
    for(i in 1:2){
      lambda[i] ~ dgamma(.05,.005)
      rho[i] <- log(lambda[i])
    }
  }
  ", fill = TRUE)
sink()
## DATA
jags.data <- list(y=jd$jags.data$y,n=jd$jags.data$n,X=jd$jags.data$X)
## INITIAL VALUES
inits <- function()
  list(b = jd$jags.ini$b, lambda = jd$jags.ini$lambda)
## PARAMETERS MONITORED
parameters <- c("b", "rho", "mu", "sd", "eps")
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
### Run the BUGS model
sta <- Sys.time()
outM <- jags(jags.data,inits,parameters,"SSBMale.GAMM",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
# Time difference of 13.79493 secs
## SELECT POSTERIOR SIMULATIONS FOR INFERENCE
predMMu <- rbind(outM$samples[[1]],outM$samples[[2]],outM$samples[[3]])
## INFERENCE FOR PREDICTION DAYS (from Nov. 1st): 48, 52, 57 
predMMu48 <- predMMu[,which(dimnames(predMMu)[[2]]=="mu[48]")]
predMMu52 <- predMMu[,which(dimnames(predMMu)[[2]]=="mu[52]")]
predMMu57 <- predMMu[,which(dimnames(predMMu)[[2]]=="mu[57]")]
NMm.nc48 <- (NTs28+NNs28)/predMMu48
NMm.nc52 <- (NTs28+NNs28)/predMMu52
NMm.nc57 <- (NTs28+NNs28)/predMMu48
##-------
## SCALING FACTORS DUE TO HAUL-OUT DIFFERENCES BETWEEN MVK AND SSB
## BIRD ISLAND - SSB
Stmc <- read.csv("raw_data/SSB_total_male_counts.csv",header=T)
Sdat <- data.frame(y=Stmc[1:61,paste("X",2022,sep="")],x0=1:61)
bS <- gam(y ~ s(x0),family=poisson(link="log"),data=Sdat,method="REML")
## MAIVIKEN
Mtmc <- read.csv("raw_data/maiviken.csv", header = T)
Mdat <- data.frame(y=Mtmc$males[Mtmc$year==2022],x0=Mtmc$day[Mtmc$year==2022])
bM <- gam(y ~ s(x0),family=poisson(link="log"),data=Mdat,method="REML")
## UNCERTAINTY FROM POSTERIOR SIMULATIONS
## SSB
XpS <- predict(bS,newdata=data.frame(x0=1:61),type="lpmatrix")
set.seed(666)
brS <- mvrnorm(n=15000,coef(bS),bS$Vp)
bSsim <- exp(XpS %*% t(brS))
zSsim <- apply(bSsim,2,FUN=function(x) (x-min(x))/(max(x)-min(x)))
## MVK
XpM <- predict(bM,newdata=data.frame(x0=1:61),type="lpmatrix")
brM <- mvrnorm(n=15000,coef(bM),bM$Vp)
bMsim <- exp(XpM %*% t(brM))
zMsim <- apply(bMsim,2,FUN=function(x) (x-min(x))/(max(x)-min(x)))
##
dNsim <- zMsim - zSsim   # zM > zS 
dPsim <- zSsim - zMsim   # zM < zS 
mkNsim <- bMsim*dNsim/zMsim
mkPsim <- bMsim*dPsim/zMsim
mKNsim <- 1-dNsim/zMsim
mKPsim <- 1+dPsim/zMsim
##
## ("delta arrows" IN FIG 3)
## DAYS: 48, 52, 57 -> -, -, -
## KAPPA ESTIMATES
round(c(mean=mean(mKPsim[48,]),quantile(mKPsim[48,],c(0.025,0.975))),2)
# mean  2.5% 97.5% 
# 0.83  0.67  1.02 
round(c(mean=mean(mKPsim[52,]),quantile(mKPsim[52,],c(0.025,0.975))),2)
# mean  2.5% 97.5% 
# 0.83  0.64  1.08
round(c(mean=mean(mKNsim[57,]),quantile(mKNsim[57,],c(0.025,0.975))),2)
# mean  2.5% 97.5% 
# 0.89  0.61  1.29
##-----------------------------------------------------------------------------------------
##  OBTAINING ABUNDANCE ESTIMATES FOR COMPARABLE SURVEY AREAS BETWEEN 2007 OR 2009 AND 2022
##  2021-22 RPAS COUNT DATA:
#   -> Hound Bay: 18/12/2021 - West: 36.2651099°W 54.3880689°S East: 36.2555069°W 54.3918537°S
#   -> Maiviken: 27/12/2021 - West: 36.5024586°W 54.2531057°S East: 36.4984819°W 54.2443644°S
#   -> Stromness Harbour: 22/12/2021
#   -> Husvik Harbour: 22/12/2021
#   Site no.	Location	        Count Method	No. of Females	No. of Pups	No. of Males
#       06	 Hound Bay	        Manual	                 931	        664	         196
#       09	 Maiviken	          Manual	               1,003	        725	         202
#       23	 Stromness Harbour	Manual	               2,043	      2,043	       1,255
#       24	 Husvik Harbour	    Manual	               2,855	      2,616	       1,503
###----------------------------------------------------------------------------------------
## SOURCE SAVED VECTORS OF POSTERIOR DISTRIBUTION SIMULATIONS FOR LOCATIONS AT SOUTH GEORGIA
#
## -> BIRD ISLAND - OBTAINED WITH SCRIPTS 'BIMatureFemaleAbundance.r' and 'BIMatureMaleAbundance.r'
BI04f <- readRDS("processed_data/BI04f.Rds")
BI04m <- readRDS("processed_data/BI04m.Rds")
BI04 <- BI04f+BI04m
BI09f <- readRDS("processed_data/BI09f.Rds")
BI09m <- readRDS("processed_data/BI09m.Rds")
BI09 <- BI09f+BI09m
BI24f <- readRDS("processed_data/BI24f.Rds")
BI24m <- readRDS("processed_data/BI24m.Rds")
BI24 <- BI24f+BI24m
## NEXT ESTIMATES OBTAINED FROM DRYAD - https://doi.org/10.5061/dryad.2rbnzs7tc
## -> Hound Bay
HB09f <- readRDS("processed_data/HB09.Rds")
HB09m <- readRDS("processed_data/HB09.Rds")
HB09 <- readRDS("processed_data/HB09.Rds")
## -> Maiviken
#  2007
MV07f <- readRDS("processed_data/MV07f.Rds")
MV07m <- readRDS("processed_data/MV07m.Rds")
MV07 <- readRDS("processed_data/MV07.Rds")
## 2009
MV09f <- readRDS("processed_data/MV09f.Rds")
MV09m <- readRDS("processed_data/MV09m.Rds")
MV09 <- readRDS("processed_data/MV09.Rds")
## -> Stromness Harbour
ST09f <- readRDS("processed_data/ST09f.Rds")
ST09m <- readRDS("processed_data/ST09m.Rds")
ST09 <- readRDS("processed_data/ST09.Rds")
## -> Husvik Harbour
HV09f <- readRDS("processed_data/HV09f.Rds")
HV09m <- readRDS("processed_data/HV09m.Rds")
HV09 <- readRDS("processed_data/HV09.Rds")
##---------------------------
##  GENERATE ABUNDANCE ESTIMATES FROM DECEMBER 2021 RPAS COUNTS
## lon - lat COORDINATES FOR THE SPECIAL STUDY AREAS AT BI AND MAIVIKEN
SSB <- cbind(lon=-38.05056,lat=-54.01162)
MVK <- cbind(lon=-36.49972,lat=-54.25139)
## lon - lat COORDINATES FOR 2022 SURVEY AREAS
HOB <- cbind(lon=-36.264531,lat=-54.391703) # Hound Bay
STH <- cbind(lon=-36.71109,lat=-54.15795)   # Stromness Harbour
HVH <- cbind(lon=-36.70928,lat=-54.18021)   # Husvik Harbour
#----------------------------
#   -> HOUND BAY - December 18, 2021
##  DISTANCES FROM HOUND BAY TO SSB AND MVK
w1 <- as.vector(1/(geodist(HOB,SSB,measure="geodesic")/1000))
w2 <- as.vector(1/(geodist(HOB,MVK,measure="geodesic")/1000))
## ABUNDANCE
## FEMALES
HB22f <- (w1*931*NMf.nc48 + w2*931*fKPsim[48,]*NMf.nc48)/(w1+w2)
## MALES
HB22m <- (w1*196*NMm.nc48 + w2*196*mKPsim[48,]*NMm.nc48)/(w1+w2)
HB22 <- HB22f + HB22m
saveRDS(HB22,"processed_data/HB22.Rds")
##
#  -> Maiviken - December 27, 2021
MV22f <- 1003 * fKNsim[57,] * NMf.nc57
## MALES
MV22m <- 202 * mKNsim[57,] * NMm.nc57
MV22 <- MV22f + MV22m
saveRDS(MV22,"processed_data/MV22.Rds")
##
#  -> Stromness Harbour
## DISTANCES TO SSB AND MVK
w1 <- as.vector(1/(geodist(STH,SSB,measure="geodesic")/1000))
w2 <- as.vector(1/(geodist(STH,MVK,measure="geodesic")/1000))
## ABUNDANCE
## FEMALES
ST22f <- (w1*2043*NMf.nc52 + w2*2043*fKPsim[52,]*NMf.nc52)/(w1+w2)
## MALES
ST22m <- (w1*.205*2043*NMm.nc52 + w2*.205*2043*mKPsim[52,]*NMm.nc52)/(w1+w2)
ST22 <- ST22f+ST22m
saveRDS(ST22,"processed_data/ST22.Rds")
##
#  -> Husvik Harbour
w1 <- as.vector(1/(geodist(HVH,SSB,measure="geodesic")/1000))
w2 <- as.vector(1/(geodist(HVH,MVK,measure="geodesic")/1000))
## ABUNDANCE
## FEMALES
HV22f <- (w1*2855*NMf.nc52 + w2*2855*fKPsim[52,]*NMf.nc52)/(w1+w2)
## MALES
HV22m <- (w1*.205*2855*NMm.nc52 + w2*.205*2855*mKPsim[52,]*NMm.nc52)/(w1+w2)
HV22 <- HV22f + HV22m
saveRDS(HV22,"processed_data/HV22.Rds")
##
##------------------------------------------------------------------------------
## GENERATION TIME - OBTAINED FOR 2008-09 WITH SCRIPTS 'AFSFemale_vital_rates.R'
##  AND 'AFSFemale_vital_rates.R'.
##---------
## WEIGHTED MEAN
#  2009 BI MATURE MALE AND FEMALE DATA - GENERATION LENGTH IS A VECTOR OF ONLY 10000
#   ABUNDANCE VECTORS NEED TO BE DOWNSIZED FROM 15000, ACCORDINGLY
BI09m <- readRDS("processed_data/BI09m.Rds")
BI09m <- BI09m[1:10000]
BI09f <- readRDS("processed_data/BI09f.Rds")
BI09f <- BI09f[1:10000]
## SOURCE SAVED GENERATION TIMES
gentChM <- readRDS("processed_data/gentChM.Rds")
gentChF <- readRDS("processed_data/gentChF.Rds")[1:10000]
gentCHw <- (gentChM*BI09m + gentChF*BI09f)/(BI09m + BI09f) 
round(c(mean=mean(gentCHw),SD=sd(gentCHw)),3)
# mean     SD 
# 8.637 0.258
## Save data
saveRDS(gentCHw,"processed_data/gentCHw.Rds")
##
##--------------------------------------------
##  REDUCTION UNDER CRITERION A1/A2, BASED ON ESTIMATES OF POPULATION SIZE FROM 2 YEARS
#   -> THREE GENERATION DECLINE WITH EXPONENTIAL MODEL (CONSTANT RATE OF DECLINE) FOR 
#      EACH SURVEYED COLONY
##-------------------------
##  FUNCTIONS TO CALCULATE REDUCTION AND ANNUAL CHANGE
##  REQUIRED DATA AND INTERNAL CALCULATIONS:
#   ay is assessment year
#   ap is assessment period
#   y3a is year three generations ago
#   ac is annual change
#   yby3y1 is years between 3 generations ago and year1
#   yby2ay is years between year 2 and assessment year
#   cb3y1 is change between 3 generations ago and year 1
#   cby2ay is change between year 2 and assessment year
#   P3 is population 3 generations ago
#   PC is current population
#   C3G is change in three generations (or reduction)
#
##  -> REQUIRE MAIN ARGUMENT yp, WHICH IS A MATRIX WITH 4 COLUMNS X NUMBER OF SIMULATIONS, TYPICALLY 15000
##     e.g. cbind(rep(2007,15000),MV07,rep(2022,15000),MV22)
##
## ANNUAL CHANGE
annualChange <- function(i,yp)
  as.vector((yp[i,4]/yp[i,2])^(1/(yp[i,3]-yp[i,1])))
## REDUCTION
reduction <- function(i,yp, G=gentCHw,ay=2025,short=T,lowc=F){
  ap <- min(c(max(c(10,G[i]*3)),100))
  y3a <- ay - ap
  ac <- (yp[i,4]/yp[i,2])^(1/(yp[i,3]-yp[i,1]))
  yby3y1 <- yp[i,1] - y3a
  yby2ay <- ay - yp[i,3]
  if(lowc)
    cb3y1 <- 0.98^yby3y1
  else
    cb3y1 <- ac^yby3y1
  cby2ay <- ac^yby2ay
  P3 <-yp[i,2]/cb3y1
  PC <- cby2ay*yp[i,4]
  C3G <- (PC-P3)/P3
  if(short)
    C3G
  else(
    c(yp[i,],yby3y1,yby2ay,ac,100*(1-ac),cb3y1,cby2ay,P3,PC,C3G)
  )
}
####
##  CALCULATE ANNUAL CHANGE AND REDUCTIONS
## MAIVIKEN
MVred <- sapply(1:10000,reduction, yp=cbind(rep(2007,10000),MV07[1:10000],rep(2022,10000),MV22[1:10000]),short=F)
MVm <- rowMeans(MVred)
MVsd <- apply(MVred,1,sd)

##
## SOUTG GEORGIA ESTIMATES FOR 2007-09
SGC09 <- readRDS("raw_data/SGC09.Rds")
SGC07 <- readRDS("raw_data/SGC07.Rds")
SGC07Compare <- readRDS("raw_data/SGC07Compare.Rds")
## Ratio female to pup counts
femToPup <- readRDS("processed_data/femToPup.Rds")
femToPup50 <- apply(femToPup, 1, median)
femToPup95 <- apply(femToPup, 1, quantile, 0.975)
## - SSB Integrated model estimates for 2001 to 2022 - FOR REPRODUCIBILITY OF OLD ESTIMATES
preds <- readRDS("processed_data/AFSFemaleSimulations200122Best.Rds")
## Proportions of mature females and breeding females wrt SSB daily count prediction
NMats <- NBIm <- matrix(NA, 7, 15000)
for(i in 1:7){
  NMats[i, ] <-  preds[, which(dimnames(preds)[[2]] == paste("NMat.nc[", i, "]", sep = ""))]
  NBIm[i, ] <-  preds[, which(dimnames(preds)[[2]] == paste("NBI.nc[", i, "]", sep = ""))]
}
## FOR 2006-07 - IF REQUIRED
predF07 <- readRDS("processed_data/FemaleGAMM07.Rds")
nc07 <-  rowMeans(predF07[, c(which(dimnames(predF07)[[2]]=="mu[42]"):which(dimnames(predF07)[[2]]=="mu[44]"))])
NBI07 <- preds[, which(dimnames(preds)[[2]] == "NBI[7]")]
NN07 <- preds[, which(dimnames(preds)[[2]] == "NN[7]")]
NMat.nc07 <- (NBI07 + NN07)/nc07
NBI.nc07 <- NBI07/nc07
##
## -MALE DATA-
predsM <- readRDS("processed_data/AFSMaleSimulations19952007Best.Rds")
## POSTERIOR SIMULATIONS FOR RATIOS OF TOTAL MATURE MALES AND TERRITORIALS TO PREDICTED SSB COUNTS AT DATE
NATot.ncs <- predsM[, which(dimnames(predsM)[[2]] == "NATot.nc[1]"):which(dimnames(predsM)[[2]] == "NATot.nc[7]")]
NT.ncs <- predsM[, which(dimnames(predsM)[[2]] == "NT.nc[1]"):which(dimnames(predsM)[[2]] == "NT.nc[7]")]
## FOR 2006-07  - IF REQUIRED
predM07 <- readRDS("processed_data/MaleGAMM07.Rds")
ncM07 <-  rowMeans(predM07[, c(which(dimnames(predM07)[[2]] == "mu[42]"):which(dimnames(predM07)[[2]] == "mu[44]"))])
NT07 <- predsM[, which(dimnames(predsM)[[2]] == "NT[13]")]
NN07 <- predsM[, which(dimnames(predsM)[[2]] == "NN[13]")]
NATot.nc07 <- (NT07 +  NN07) / ncM07
NT.nc07 <- NT07 / ncM07






##-------------------------------------------------------------------------------------
##  SOUTH GEORGIA PROJECTION
SGC09 <- readRDS("raw_data/SGC09.Rds")
SGC07 <- readRDS("raw_data/SGC07.Rds")
## SOUTH GEORGIA STITCH LON, LAT
x0 <- rbind(SGC09[,c("lon", "lat")],SGC07[,c("lon", "lat")])
## VECTORS OF DISTANCES TO DIFFERENT SURVEY COLONIES - 1 MVK; 2 HB; 3 ST; 4 HV; 5 BI 
x0 <- rbind(SGC09[,c("lon", "lat")],SGC07[,c("lon", "lat")])
x1 <- cbind(lon=-36.49972,lat=-54.25139)
GD1 <- geodist(x0,x1,measure="geodesic")/1000
x2 <- cbind(lon=-36.264531,lat=-54.391703)
GD2 <- geodist(x0,x2,measure="geodesic")/1000
x3 <- cbind(lon=-36.71109,lat=-54.15795)
GD3 <- geodist(x0,x3,measure="geodesic")/1000
x4 <- cbind(lon=-36.70928,lat=-54.18021)
GD4 <- geodist(x0,x4,measure="geodesic")/1000
x5 <- cbind(lon=-1*(38+03/60+02.02/3600),lat=-1*(54+00/60+41.84/3600))
GD5 <- geodist(x0,x5,measure="geodesic")/1000
GD <- cbind(GD1,GD2,GD3,GD4,GD5)
##
## OBTAIN ESTIMATES OF ANNUAL CHANGE FOR SURVEYED COLONIES
MVAC <- sapply(1:10000,annualChange, yp=cbind(rep(2007,10000),MV07[1:10000],rep(2022,10000),MV22[1:10000]))
HBAC <- sapply(1:10000,annualChange,yp=cbind(rep(2009,10000),HB09[1:10000],rep(2022,10000),HB22[1:10000]))
STAC <- sapply(1:10000,annualChange,yp=cbind(rep(2009,10000),ST09[1:10000],rep(2022,10000),ST22[1:10000]))
HVAC <- sapply(1:10000,annualChange,yp=cbind(rep(2009,10000),HV09[1:10000],rep(2022,10000),HV22[1:10000]))
BIAC <- sapply(1:10000,annualChange,yp=cbind(rep(2004,10000),BI04[1:10000],rep(2024,10000),BI24[1:10000]))
##
## GENERATE ANNUAL CHANGE VECTORS FOR SOUTH GEORGIA LOCATIONS WEIGHTED BY DISTANCE TO COLSEST SURVEYED COLONY
sum1 <- (cbind(MVAC)%*%as.vector(GD1)+cbind(HBAC)%*%as.vector(GD2)+cbind(STAC)%*%as.vector(GD3)+
           +cbind(HVAC)%*%as.vector(GD4)+cbind(BIAC)%*%as.vector(GD5))
sum2 <- t(t(sum1)/rowSums(GD))          ## WEIGHTED AVERAGE
EXP <- sum2^(2022-2009)
##
##  FEMALES - GENERATE ABUNDANCE ESTIMATES BY STITCH
source("code/getFemales.R")
source("code/getFemales07.R")
SG09F <- cbind(sapply(1:dim(SGC09)[1], getFemales),sapply(1:dim(SGC07)[1], getFemales07))
SGNTsimf <- SG09F[1:10000,] * EXP
SGNTf <- rowSums(SGNTsimf)
## SG MATURE FEMALE POPULATION IN 2022
round(c(mean(SGNTf),quantile(SGNTf,c(0.025,0.975))))
#          2.5%  97.5% 
# 834319 709490 979666
## TOTAL FEMALE POPULATION IN 2022
preds <- readRDS("processed_data/AFSFemaleSimulations200125.Rds")
NBIs22 <- preds[, which(dimnames(preds)[[2]] == "NBI[22]")]
NNs22 <- preds[, which(dimnames(preds)[[2]] == "NN[22]")]
Npres22 <- preds[, which(dimnames(preds)[[2]] == "Npre[22]")]
preF22 <-  1 + Npres22 / (NBIs22 + NNs22)   ### scaling factor to account for prebreeders; 
round(c(mean(SGNTf*preF22[1:10000]),quantile(SGNTf*preF22[1:10000],c(0.025,0.975))))
#            2.5%   97.5% 
# 1112346  931082 1326136 
##
## MALES - GENERATE ABUNDANCE ESTIMATES BY STITCH
source("code/getMales.R")
source("code/getMales07.R")
SG09M <- cbind(sapply(1:dim(SGC09)[1], getMales),sapply(1:dim(SGC07)[1], getMales07))
SGNTsimm <- SG09M[1:10000,] * EXP
SGNTm <- rowSums(SGNTsimm)
## SG MATURE MALE POPULATION IN 2022
round(c(mean(SGNTm),quantile(SGNTm,c(0.025,0.975))))
#          2.5%  97.5% 
# 190478 152610 237723
## TOTAL MALE POPULATION IN 2022
predsM <- readRDS("processed_data/AFSMaleSimulations19952025.Rds")
NTs28 <- predsM[,which(dimnames(predsM)[[2]]==paste("NT[",29,"]",sep=""))]
NNs28 <- predsM[,which(dimnames(predsM)[[2]]==paste("NN[",29,"]",sep=""))]
Npres28 <- predsM[,which(dimnames(predsM)[[2]]==paste("NPall[",29,"]",sep=""))]
preM28 <-  1 + Npres28 /(NTs28+NNs28)
round(c(mean(SGNTm*preM28[1:10000]),quantile(SGNTm*preM28[1:10000],c(0.025,0.975))))
#           2.5%   97.5% 
# 871342  585158 1186420 
##------
## TOTAL POPULATION IN 2022
round(c(mean(SGNTf*preF22[1:10000]+SGNTm*preM28[1:10000]),quantile(SGNTf*preF22[1:10000]+SGNTm*preM28[1:10000],c(.025,.975))))
#            2.5%   97.5% 
# 1983687 1606408 2411876
##------
## TOTAL MATURE POPULATION IN 2022
SGNT <- SGNTf+SGNTm
round(c(mean(SGNTf+SGNTm),quantile(SGNTf+SGNTm,c(.025,.975))))
#            2.5%   97.5% 
# 1024797  873500 1200048
##-----------------------------------------------------------------------------------------
##  SOUTH GEORGIA REDUCTION
SGred <- sapply(1:10000,reduction,yp=cbind(rep(2009,10000),(rowSums(SG09F)+rowSums(SG09M))[1:10000],rep(2022,10000),SGNT[1:10000]),short=T,lowc=T)
c(mean(SGred),quantile(SGred,c(0.025,0.975)))
#                      SD        2.5%       97.5% 
# -0.56658464  0.04500048 -0.64751108 -0.47227540
#-------------
##  SOUTH GEORGIA ANNUAL CHANGE
SGAC <- sapply(1:10000,annualChange,yp=cbind(rep(2009,10000),(rowSums(SG09F)+rowSums(SG09M))[1:10000],rep(2022,10000),SGNT[1:10000]))
c(mean(SGAC),SD=sd(SGAC),quantile(SGAC,c(0.025,0.975)))
#                      SD        2.5%       97.5% 
# 0.960738959 0.006136052 0.948717024 0.972723656 
## 
c(mean(100*(1-SGAC)),SD=sd(100*(1-SGAC)),quantile(100*(1-SGAC),c(0.025,0.975)))
#                  SD      2.5%     97.5% 
# 3.9261041 0.6136052 2.7276344 5.1282976 
#####
