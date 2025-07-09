library(mgcv)
library(jagsUI)
library(rjags)
load.module('glm')
load.module('dic')
library(pracma)
# install.packages("remotes")
# remotes::install_github("ProcessMiner/nlcor")
library(nlcor)
library(biwavelet)
library(ggplot2)
##------------------------------------------------------------------------------------------## 
##  SET WORKING DIRECTORY AND SOURCE DDITIONAL CODE
# setwd("C:/Users/jfor/OneDrive - NERC/Workstuff/AntarcticFurSeals/IUCN_assessment/Analysis")
source("addTrans.R") # Colour transparency filter
source("glbpf.R")    # Gaussian band pass filter, from Cazelles et al. 2008
##------------------------------------------------------------------------------------------## 
##  LOAD PREVIOUSLY SAVE WORKSPACE IF REQUIRED
# load("AFSDynamicsSigny.RData")
##------------------------------------------------------------------------------------------## 
##  READ, FORMAT AND SUBSET DATA
### 1. SEAL COUNT DATA
##  SOURCE DATA - AND REFORMAT
SignySeals <- read.csv("raw_data/SignySeals.csv", header = T)
##
##  ANTARCTIC FUR SEALS
AFS <- SignySeals[, 1:2]
AFS <- with(AFS, tapply(afs, list(year), sum))
Atemp <- matrix(NA, length(1977:2025), 1)
Atemp[which(1977:2025 %in% as.numeric(dimnames(AFS)[[1]])), ] <- AFS 
## INTERPOLATION
myr <- which(1977:2025 %in% as.numeric(dimnames(AFS)[[1]]) == F)
Atemp[myr, 1] <- spline(1977:2025, Atemp[, 1], xout = (1977:2025)[myr])$y
Atemp <- ifelse(Atemp < 1, 1, Atemp)
AFS <- data.frame(year = 1977:2025, count = round(Atemp))
AFS <- AFS[-c(7:8), ]
## SAVE DATA
saveRDS(AFS, file = "processed_data/AFS.Rds")
AFS <- readRDS("processed_data/AFS.Rds")
##------------------------------------------------------------------------------------------## 
### TREND ANALYSIS
### 2. GAM ANALYSIS TO GENERATE JAGS TEMPLATE MODEL FILE, RUN 'jagam' AND OBTAIN DATA STRUCTURE FOR JAGS
dat <- data.frame(y=round(c(AFS$count[1:6],NA,NA,AFS$count[7:length(AFS$count)])),x0=1977:2025)
b0 <- gam(y ~ s(x0), family = poisson(link = "log"), data = dat, method = "REML")
## INTERPOLATE/PREDICT MISSING
dat$y[which(is.na(dat$y))]<-round(predict(b0,type="response",newdata=data.frame(x0=1977:2025))[which(is.na(dat$y))])
b0 <- gam(y~s(x0),family=poisson(link="log"),data=dat,method="REML")
jags.file <- "test.jags" 
jd <- jagam(y ~ s(x0), family = poisson(link = "log"), data = dat, file = jags.file, sp.prior = "gamma", diagonalize = TRUE)
##--------------------------------
###  BAYESIAN MODELLING WITH JAGS
###  POISSON GAM
sink("AFS.GAM")
cat("
  model{
    eta <- X %*% b ## linear predictor
    for(i in 1:n){
      mu[i] <-  exp(eta[i])
    } ## expected response
    for(i in 1:n){
      y[i] ~ dpois(mu[i])
    } ## response 
    ## Parametric effect priors CHECK
    for(i in 1:1){
      b[i] ~ dnorm(0,0.00012)
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
jags.data <- list(y = jd$jags.data$y, n = jd$jags.data$n, X = jd$jags.data$X)
## INITIAL VALUES
inits <- function()
  list(b = jd$jags.ini$b, lambda = jd$jags.ini$lambda)
## PARAMETRES MONITORED
parameters <- c("b", "rho", "mu")
#
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
#
### Run the BUGS model
sta <- Sys.time()
set.seed(666)
out.test <- jags(jags.data,inits,parameters,"AFS.GAM",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
# Time difference of 21.01203 secs
##
## SAVE DATA FOR PLOTTING
## SIMULATIONS
preds <- rbind(out.test$samples[[1]],out.test$samples[[2]],out.test$samples[[3]])
preds <- preds[,which(dimnames(preds)[[2]]=="mu[1]"):which(dimnames(preds)[[2]]=="mu[49]")]
##
AFSRes <- data.frame(year = dat[ ,2], obs=c(AFS$count[1:6],NA,NA,AFS$count[7:length(AFS$count)]),
                     Afs = dat[ ,1], AGM = colMeans(preds),AGL = apply(preds, 2, quantile, 0.025), 
                     AGH = apply(preds, 2, quantile, 0.975))
## SAVE DATA
saveRDS(AFSRes, file = "AFStrends.Rds")
##
## POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
sink("AFS.GAMM") 
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
    ## Parametric effect priors CHECK
    for(i in 1:1){
      b[i] ~ dnorm(0,0.00012)
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
jags.data <- list(y = jd$jags.data$y, n = jd$jags.data$n, X = jd$jags.data$X)
## INITIAL VALUES
inits <- function()
  list(b = jd$jags.ini$b, lambda = jd$jags.ini$lambda)
## PARAMETRES MONITORED
parameters <- c("b", "rho", "mu", "sd", "eps") #, "scale", "mu")
#
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
#
### Run the BUGS model
sta <- Sys.time()
set.seed(666)
out.test1 <- jags(jags.data,inits,parameters,"AFS.GAMM",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
# Time difference of 20.28994 secs
## OVERDISPERSION
round(c(out.test1$mean$sd,out.test1$q2.5$sd,out.test1$q97.5$sd),2)
# [1] 0.41 0.32 0.52
## SIMULATIONS
preds1 <- rbind(out.test1$samples[[1]],out.test1$samples[[2]],out.test1$samples[[3]])
preds1 <- preds1[,which(dimnames(preds1)[[2]]=="mu[1]"):which(dimnames(preds1)[[2]]=="mu[49]")]
### SAVE SIMULATIONS
saveRDS(preds1,fil="AFSsim.Rds")
##
AFSRes <- data.frame(AFSRes,AGMM=colMeans(preds1),AGML=apply(preds1,2,quantile,0.025),
                     AGMH=apply(preds1,2,quantile,0.975))
## SAVE DATA
saveRDS(AFSRes,file="AFStrends.Rds")
##
################################################################################
##  3. GOMPERTZ POPULATION MODELLING
##  3.1 NEGATIVE BINOMIAL BUGS MODEL
sink("AFSNB.jags")
cat("
  model{
    ## Alpha - Growth rate and Density-dependence
    alpha0 ~ dnorm(0,0.01)T(0,5)
    alpha1 ~ dnorm(0,tauDD)T(-2,2)
    ## PRIORS AND HYPERPARAMETERS
    tauDD <- pow(sigDD,-2)
    sigDD ~ dunif(0,2)
    nuprec  ~ dgamma(1.0E-1,1.0E-2) # dgamma(1.0E-3,1.0E-3)
    r ~ dunif(0,50)
    #
    ## INITIAL COUNTS
    tau1 <- pow(sdy1,-2)
    y0 ~ dnorm(y1,tau1)
    x[1] <- log(max(y0,1.0E-5))
    #
    ## POPULATION PROCESS
    for(j in 2:(T+1)){
      #Density-dependence and sea ice covariates
      Ex[j] <- alpha0+(1+alpha1)*x[j-1]
      x[j] ~ dnorm(Ex[j],nuprec)
    }
    # NEGATIVE BINOMIAL OBSERVATION ERRORS
    for(j in 1:T){
      log(mu[j]) <- x[j]
      p[j] <- r/(r+mu[j])
      y[j] ~ dnegbin(p[j],r)
    }
    varunshared <- 1/nuprec
  }
  ", fill = TRUE)
sink()
##
## Data
y=AFSRes$Afs
jags.data = list(y=y,y1=y[1],sdy1=y[1]*.05,T=length(y))
##
## Initial values
inits <- function() list(nuprec=0.1)
#
## PARAMETERS MONITORED
params <- c("alpha0","alpha1","nuprec","varunshared","r","p","mu")
#
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 200000
nt <- 20
nb <- 100000
nc <- 4
## Call jags and fit model
#(sta <- Sys.time())
#out1 <- jags(jags.data,inits,params,"AFSNB.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T) 
#Sys.time() - sta
# Time difference of 33.87359 secs
##
(sta <- Sys.time())
out1 <- jags(jags.data,inits,params,"AFSNB.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=F) 
Sys.time() - sta
#Time difference of 3.170583 mins
#
##  WAIC
##
samples <- rjags::jags.samples(out1$model, c("WAIC", "deviance"),type="mean",n.iter=ni,n.burnin=nb,n.thin=nt)
samples$p_waic <- samples$WAIC
samples$waic <- samples$deviance + samples$p_waic
tmp <- sapply(samples, sum)
#     WAIC deviance   p_waic     waic 
(waic_temp1 <- round(c(deviance = tmp[["deviance"]], waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]), 2))
# deviance     waic   p_waic 
#   860.31   886.54    26.23
##
## SIMULATIONS
preds2 <- rbind(out1$samples[[1]],out1$samples[[2]],out1$samples[[3]])
preds2 <- preds2[,which(dimnames(preds2)[[2]]=="mu[1]"):which(dimnames(preds2)[[2]]=="mu[49]")]
##
AFSRes <- data.frame(AFSRes,ANBM=colMeans(preds2),ANBL=apply(preds2,2,quantile,0.05),ANBH=apply(preds2,2,quantile,0.95))
## SAVE DATA
saveRDS(AFSRes,file="AFStrends.Rds")
##
## FIT GAM TO MU PREDICTIONS
# NTSimsF <- NTOT
## RERUN THIS BIT ONCE AND SAVE RESULTS
years <- 1977:2025
lpAFSS <- gamAFSS <- matrix(NA,dim(preds2)[1],dim(preds2)[2])
for(i in 1:dim(preds2)[1]){
  lpAFSS[i, ] <- glbpf(preds2[i,],pf=9.6)
  gamAFSS[i,] <- predict(gam(preds2[i,] ~ s(years)))
  print(i)  
}
saveRDS(lpAFSS,file="lpAFSS.Rds")
saveRDS(gamAFSS,file="gamAFSS.Rds")
## UPDATE RESULTS data.frame
AFSRes <- data.frame(AFSRes,lpM=colMeans(lpAFSS),lpL=apply(lpAFSS,2,quantile,0.025),lpH=apply(lpAFSS,2,quantile,0.975))
AFSRes <- data.frame(AFSRes,gamFM=colMeans(gamAFSS),gamFL=apply(gamAFSS,2,quantile,0.025),gamFH=apply(gamAFSS,2,quantile,0.975))
## SAVE DATA
saveRDS(AFSRes,file="AFStrends.Rds")
#-----------------------------------
##  PLOT
AFSRes <- readRDS("AFStrends.Rds")
plot(AFSRes$year,AFSRes$Afs,type="n",xlim=c(1977,2025),ylim=c(1200,37000),xlab="",ylab="",axes=F) #ylim=c(1200,26250)
polygon(c(AFSRes$year,rev(AFSRes$year)),c(AFSRes$gamH,rev(AFSRes$gamL)),col=addTrans("lightblue",150),border=F)
lines(AFSRes$year,AFSRes$gamM,lwd=1,col="blue")
for(i in 1:dim(AFSRes)[1])
  lines(rep(AFSRes$year[i],2),rbind(AFSRes$ANBL,AFSRes$ANBH)[,i],col="blue")
points(AFSRes$year,AFSRes$ANBM,pch=21,col="blue",bg="lightblue",cex=.75)
box()
axis(1,at=seq(1980,2025,by=5),labels=T,tck=-0.025,cex.axis=1,lty=1,padj=-1.0)
axis(1,at=seq(1977,2025,by=1),labels=F,tck=-0.012,las=2,cex.axis=1)
axis(4,at=pretty(c(1200,26250))[1:6],labels=T,tck=-0.02,las=2,cex.axis=.75,col.axis="blue",hadj=0.5)  # col.ticks="blue", hadj=0.7
axis(4,at=seq(0,27000,by=1000),labels=F,tck=-0.01,las=2,cex.axis=.75,col.ticks="blue")
mtext("a", side = 3, adj = -0.22, padj = -1.75, cex = .7, font = 2)
#
#################################################################################
##  4. CROSS-CORRELATION BI MATURE FEMALE ABUNDANCE ~ SI POPULATION INDEX
##  BI MATURE FEMALE DATA
NTSims <- readRDS("NTSims.Rds")
##  CCF ANALYSIS
AFSccf <- ccf(as.vector(colMeans(NTSims)),as.vector(AFSRes$Afs[3:49]),lag.max=10,plot=F)
##
##----------
## FIT NB MODEL WITH BI AS EXPLANATORY VARIABLE WITH LAG -3
##  4.1.1 NEGATIVE BINOMIAL BUGS MODEL - BASIC MODEL, WITHOUT BI EFFECTS
sink("AFSNB1.jags")
cat("
  model{
    ## Alpha - Growth rate and Density-dependence
    alpha0 ~ dnorm(0,0.01)T(0,5)
    alpha1 ~ dnorm(0,tauDD)T(-2,2)
    ## PRIORS AND HYPERPARAMETERS
    tauDD <- pow(sigDD,-2)
    sigDD ~ dunif(0,2)
    nuprec  ~ dgamma(1.0E-1,1.0E-2) # dgamma(1.0E-3,1.0E-3)
    r ~ dunif(0,50)
    #
    ## INITIAL COUNTS
    tau1 <- pow(sdy1,-2)
    y0 ~ dnorm(y1,tau1)
    x[1] <- log(max(y0,1.0E-5))
    #
    ## POPULATION PROCESS
    for(j in 2:(T+1)){
      #Density-dependence and sea ice covariates
      Ex[j] <- alpha0+(1+alpha1)*x[j-1]
      x[j] ~ dnorm(Ex[j],nuprec)
    }
    # NEGATIVE BINOMIAL OBSERVATION ERRORS
    for(j in 1:T){
      log(mu[j]) <- x[j]
      p[j] <- r/(r+mu[j])
      y[j] ~ dnegbin(p[j],r)
    }
    varunshared <- 1/nuprec
  }
  ", fill = TRUE)
sink()
##
## Data
y <- AFSRes$Afs[6:49]
jags.data = list(y=y,y1=y[1],sdy1=y[1]*.05,T=length(y))
##
## Initial values
inits <- function() list(nuprec=0.1)
#
## PARAMETERS MONITORED
params <- c("alpha0","alpha1","nuprec","varunshared","r","p","mu")
#
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 200000
nt <- 20
nb <- 100000
nc <- 4
## Call jags and fit model
(sta <- Sys.time())
out2 <- jags(jags.data,inits,params,"AFSNB1.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=F) 
Sys.time() - sta
##  WAIC
##
samples <- rjags::jags.samples(out2$model, c("WAIC", "deviance"),type="mean",n.iter=ni,n.burnin=nb,n.thin=nt)
samples$p_waic <- samples$WAIC
samples$waic <- samples$deviance + samples$p_waic
tmp <- sapply(samples, sum)
#     WAIC deviance   p_waic     waic 
(waic_temp2 <- round(c(deviance = tmp[["deviance"]], waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]), 2))
# deviance     waic   p_waic 
#   777.93   801.10    23.17 
##----
## SIMULATIONS
preds3 <- rbind(out2$samples[[1]],out2$samples[[2]],out2$samples[[3]])
preds3 <- preds3[,which(dimnames(preds3)[[2]]=="mu[1]"):which(dimnames(preds3)[[2]]=="mu[44]")]
##
AFSRes <- data.frame(AFSRes,ANB2M=c(rep(NA,5),colMeans(preds3)),ANB2L=c(rep(NA,5),apply(preds3,2,quantile,0.05)),
  ANB2H=c(rep(NA,5),apply(preds3,2,quantile,0.95)))
## SAVE DATA
saveRDS(AFSRes,file="AFStrends.Rds")

##----------------------------
##  4.1.2 NEGATIVE BINOMIAL BUGS MODEL WITH BI AS EXPLANATORY VARIABLE WITH LAG -3
sink("AFSNB2.jags")
cat("
  model{
    ## alpha -> Growth rate and Density-dependence 
    alpha0 ~ dnorm(0,0.01)T(0,5)
    alpha1 ~ dnorm(0,tauDD)T(-2,2)
    alpha2 ~ dnorm(0,0.01)
    ## PRIORS AND HYPERPARAMETERS
    tauDD <- pow(sigDD,-2)
    sigDD ~ dunif(0,2)
    nuprec  ~ dgamma(1.0E-1,1.0E-2) # dgamma(1.0E-3,1.0E-3)
    r ~ dunif(0,50)
    #
    ## BIfemales - sBIf SDs and nBIf Ns 
    for(t in 1:T){
      tau.BIf[t] <- pow(sBIf[t],-2)
      mBIf[t] ~ dnorm(nBIf[t],tau.BIf[t])
      BIf[t] <- log(max(mBIf[t],1.0E-5))-6
    }
    #
    ## INITIAL COUNTS
    tau1 <- pow(sdy1,-2)
    y0 ~ dnorm(y1,tau1)
    x[1] <- log(max(y0,1.0E-5))
    #
    ## POPULATION PROCESS
    for(j in 2:(T+1)){
      #Density-dependence and BI female covariate
      Ex[j] <- alpha0+(1+alpha1)*x[j-1]+alpha2*BIf[j-1]
      x[j] ~ dnorm(Ex[j],nuprec)
    }
    # NEGATIVE BINOMIAL OBSERVATION ERRORS
    for(j in 1:T){
      log(mu[j]) <- x[j]
      p[j] <- r/(r+mu[j])
      y[j] ~ dnegbin(p[j],r)
    }
    varunshared <- 1/nuprec
  }
  ", fill = TRUE)
sink()
##
## Data
y <- AFSRes$Afs[6:49]
jags.data <- list(y=y,y1=y[1],sdy1=y[1]*.05,T=length(y),sBIf=as.vector(apply(NTSims[,1:44],2,sd)),
  nBIf=as.vector(colMeans(NTSims[,1:44])))
##
## Initial values
inits <- function() list(nuprec=0.1)
#
## PARAMETERS MONITORED
params <- c("alpha0","alpha1","alpha2","nuprec","varunshared","r","p","mu","BIf")
#
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 200000
nt <- 20
nb <- 100000
nc <- 4
## Call jags and fit model
(sta <- Sys.time())
out3 <- jags(jags.data,inits,params,"AFSNB2.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=F) 
Sys.time() - sta
# Time difference of 2.540582 mins
#
##  WAIC
##
samples <- rjags::jags.samples(out3$model, c("WAIC", "deviance"),type="mean",n.iter=ni,n.burnin=nb,n.thin=nt)
samples$p_waic <- samples$WAIC
samples$waic <- samples$deviance + samples$p_waic
tmp <- sapply(samples, sum)
#     WAIC deviance   p_waic     waic 
(waic_temp3 <- round(c(deviance = tmp[["deviance"]], waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]), 2))
# deviance     waic   p_waic 
#   803.93   828.44    24.51
## SIMULATIONS
preds4 <- rbind(out3$samples[[1]],out3$samples[[2]],out3$samples[[3]])
preds4 <- preds4[,which(dimnames(preds4)[[2]]=="mu[1]"):which(dimnames(preds4)[[2]]=="mu[44]")]
out3res <- rbind(out3$samples[[1]],out3$samples[[2]],out3$samples[[3]])
saveRDS(out3res,file="AFSNegBinBISI.Rds")

##
AFSRes <- data.frame(AFSRes,ANB3M=c(rep(NA,5),colMeans(preds4)),ANB3L=c(rep(NA,5),apply(preds4,2,quantile,0.05)),
  ANB3H=c(rep(NA,5),apply(preds4,2,quantile,0.95)))
## SAVE DATA
saveRDS(AFSRes,file="AFStrends.Rds")
##-----------------------------------
##  4.1.3 NEGATIVE BINOMIAL BUGS MODEL WITH BI AS EXPLANATORY VARIABLE WITH LAG -6
sink("AFSNB26.jags")
cat("
  model{
    ## alpha -> Growth rate and Density-dependence 
    alpha0 ~ dnorm(0,0.01)T(0,5)
    alpha1 ~ dnorm(0,tauDD) # T(-2,2)
    alpha2 ~ dnorm(0,0.01)
    ## PRIORS AND HYPERPARAMETERS
    tauDD <- pow(sigDD,-2)
    sigDD ~ dunif(0,2)
    nuprec  ~ dgamma(1.0E-1,1.0E-2) # dgamma(1.0E-3,1.0E-3)
    r ~ dunif(0,50)
    #
    ## BIfemales - sBIf SDs and nBIf Ns 
    for(t in 1:T){
      tau.BIf[t] <- pow(sBIf[t],-2)
      mBIf[t] ~ dnorm(nBIf[t],tau.BIf[t])
      BIf[t] <- log(max(mBIf[t],1.0E-5))-6
    }
    #
    ## INITIAL COUNTS
    tau1 <- pow(sdy1,-2)
    y0 ~ dnorm(y1,tau1)
    x[1] <- log(max(y0,1.0E-5))
    #
    ## POPULATION PROCESS
    for(j in 2:(T+1)){
      #Density-dependence and BI female covariate
      Ex[j] <- alpha0+(1+alpha1)*x[j-1]+alpha2*BIf[j-1]
      x[j] ~ dnorm(Ex[j],nuprec)
    }
    # NEGATIVE BINOMIAL OBSERVATION ERRORS
    for(j in 1:T){
      log(mu[j]) <- x[j]
      p[j] <- r/(r+mu[j])
      y[j] ~ dnegbin(p[j],r)
    }
    varunshared <- 1/nuprec
  }
  ", fill = TRUE)
sink()
##
## Data
## cor.test(AFSRes$Afs[9:48],as.vector(colMeans(NTSims[,1:40])))
y <- AFSRes$Afs[9:49]
jags.data <- list(y=y,y1=y[1],sdy1=y[1]*.05,T=length(y),sBIf=as.vector(apply(NTSims[,1:41],2,sd)),
  nBIf=as.vector(colMeans(NTSims[,1:41])))
##
## Initial values
inits <- function() list(nuprec=0.1)
##
## PARAMETERS MONITORED
params <- c("alpha0","alpha1","alpha2","nuprec","varunshared","r","p","mu","BIf")
#
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 200000
nt <- 20
nb <- 100000
nc <- 4
## Call jags and fit model
(sta <- Sys.time())
out4 <- jags(jags.data,inits,params,"AFSNB26.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T) 
Sys.time() - sta
# Time difference of 1.210611 mins
##
## SIMULATIONS
preds5 <- rbind(out4$samples[[1]],out4$samples[[2]],out4$samples[[3]])
preds5 <- preds5[,which(dimnames(preds5)[[2]]=="mu[1]"):which(dimnames(preds5)[[2]]=="mu[41]")]
out4res <- rbind(out4$samples[[1]],out4$samples[[2]],out4$samples[[3]])
saveRDS(out4res,file="AFSNegBinBISI6.Rds")
##
AFSRes <- data.frame(AFSRes,ANB4M=c(rep(NA,8),colMeans(preds5)),ANB4L=c(rep(NA,8),apply(preds5,2,quantile,0.05)),
  ANB4H=c(rep(NA,8),apply(preds5,2,quantile,0.95)))
## SAVE DATA
saveRDS(AFSRes,file="AFStrends.Rds")
##-----------------------------------
##  4.3 POPULATION SYNCHRONY ANALYSIS
## lag 0
cor.test(AFSRes$Afs[3:49],as.vector(colMeans(NTSims[,1:47])))$estimate
# 0.2369424
## lag -1
cor.test(AFSRes$Afs[4:49],as.vector(colMeans(NTSims[,1:46])))$estimate
# 0.3641519
## lag -2
cor.test(AFSRes$Afs[5:49],as.vector(colMeans(NTSims[,1:45])))$estimate
# 0.441601 
## lag -3
cor.test(AFSRes$Afs[6:49],as.vector(colMeans(NTSims[,1:44])))$estimate
# 0.5247247
## lag -4
cor.test(AFSRes$Afs[7:49],as.vector(colMeans(NTSims[,1:43])))$estimate
# 0.4565654
## lag -5
cor.test(AFSRes$Afs[8:49],as.vector(colMeans(NTSims[,1:42])))$estimate
# 0.4752209
## lag -6
cor.test(AFSRes$Afs[9:49],as.vector(colMeans(NTSims[,1:41])))$estimate
# 0.6484662
## lag -7
cor.test(AFSRes$Afs[10:49],as.vector(colMeans(NTSims[,1:40])))$estimate
# 0.7445995
##--------------------
##  4.3.1 NEGATIVE BINOMIAL BUGS MODEL WITH SI,BI SYNCHRONY DYNAMICS
### TOTAL VARIANCE MODEL (DENSITY-DEPENDENCE-ONLY MODELS)
### BUGS MODEL
sink("AFSNB30.jags")
cat("
  model{
    ## Alphas - 1 GROWTH RATE 2DENSITY-DEPENDENCE
    for(i in 1:S){
      alpha0[i] ~ dnorm(mualpha0,tau0)
      alpha1[i] ~ dnorm(0,tauDD)
      #alpha0[i] ~ dnorm(mualpha0,tau0)T(0,5)
      #alpha1[i] ~ dnorm(0,tauDD)T(-2,2)
    }
    ## PRIORS AND HYPERPARAMETERS
    mualpha0 ~ dnorm(0,0.01)
    tau0 <- pow(sd0,-2)
    sd0 ~ dunif(0,2)
    tauDD <- pow(sigDD,-2)
    sigDD ~ dunif(0,2)
    #
    ## Random (synchrony) terms or shared variance
    for(j in 1:T){
      delta[j]~dnorm(0,deltaprec)
    }
    deltaprec~dgamma(1.0E-1,1.0E-1)
    ## asynchrony or unshared term
    for(i in 1:S){
      epsilonprec[i]~dgamma(1.0E-1,1.0E-1)  #dgamma(1.0,1.0E-1)
    }
    ## INITIAL COUNTS
    for(i in 1:S){
      tau1[i] <- pow(sdy1[i],-2)
      y0[i] ~ dnorm(y1[i],tau1[i])
      x[1,i] <- log(max(y0[i],1.0E-5))
    }
    ## LIKELIHOOD - SYSTEM PROCESS
    for(i in 1:S){
      for(j in 2:(T+1)){
        #Density-dependence and synchrony component
        Ex[j,i] <- alpha0[i]+(1+alpha1[i])*x[j-1,i]+delta[j-1]
        #Asynchrony component
        x[j,i] ~ dnorm(Ex[j,i],epsilonprec[i])
      }
    }
    #  OBSERVATION PROCESS (On real scale and Negative binomial observation errors)
    ## GENERATE PRIORS FOR r, GET p, AND SAMPLE ABUNDANCE 
    for(i in 1:S){
      r[i] ~ dunif(0,50)
      for(j in 1:T){
        log(mu[j,i]) <- x[j,i]
        p[j,i] <- r[i]/(r[i]+mu[j,i])
        y[j,i] ~ dnegbin(p[j,i],r[i])
      }
    }
    ## Derived variances and ICC
    varshared <- 1/deltaprec
    for(i in 1:S){
      varunshared[i] <- 1/epsilonprec[i]
      ICC[i] <- varshared/(varshared+varunshared[i])
    }
    ICCall <- varshared/(varshared+sum(varunshared[]))      
  }
  ", fill = TRUE)
sink()
##
## Data
y <- cbind(round(as.vector(colMeans(NTSims[,1:47]))),as.vector(AFSRes$Afs[3:49]))
jags.data=list(y=y,y1=y[1,],sdy1=c(as.vector(apply(NTSims[,1:44],2,sd))[1],y[1,1]*.05),T=dim(y)[1],S=dim(y)[2])
##
##
## Initial values
inits <- function() list(deltaprec=0.1,epsilonprec=rep(0.01,2))
##
## Parameters monitored
params <- c("r","alpha0","mualpha0","alpha1","deltaprec","epsilonprec","varshared","varunshared","ICC","ICCall","p")
##
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 200000
nt <- 20
nb <- 100000
nc <- 4
## Call jags and fit model
(sta <- Sys.time())
out30 <- jags(jags.data,inits,params,"AFSNB30.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
## Time difference of 1.597363 mins
#
## GENERATE RESULTS
outT <- out30 
outTdf <- data.frame(par=c("rBI","rBI","mups","upsbi","upssi","nubi","nusi","sigd","sigebi","sigesi","ICCbi","ICCsi","ICC"),
  mean=c(outT$mean$r,outT$mean$mualpha0,outT$mean$alpha0,outT$mean$alpha1,outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
  sd=c(outT$sd$r,outT$sd$mualpha0,outT$sd$alpha0,outT$sd$alpha1,outT$sd$varshared,outT$sd$varunshared,outT$sd$ICC,outT$sd$ICCall))
ouTtable <- as.data.frame(outTdf$par)
names(ouTtable)[1]<-"par"
ouTtable$lag0 <- paste(round(outTdf$mean,2),paste(rep("(",dim(outTdf)[1]),round(outTdf$sd,2),rep(")",dim(outTdf)[1]),sep=""))
##
predsT <- rbind(outT$samples[[1]],outT$samples[[2]],outT$samples[[3]])
varnm <-c("varshared","varunshared[1]","varunshared[2]","ICC[1]","ICC[2]","ICCall")
quant <- quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[1])],c(0.05,0.95))
for(i in 2:6)
  quant <- rbind(quant,quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[i])],c(0.05,0.95)))
ICCtable <- data.frame(par=c("sigdelta","sigepsiBI","sigepsiSI","ICCBI","ICCSI","ICC"),
   mean=c(outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
   lowCI=quant[,1],highCI=quant[,2],lag=rep("lag0",6))
## and clean up
rm(quant,predsT,outT,outTdf)
##---------
##  4.3.2 NEGATIVE BINOMIAL BUGS MODEL WITH SI,BI SYNCHRONY DYNAMICS - SI lagged -1
### BUGS MODEL
sink("AFSNB31.jags")
cat("
  model{
    ## Alphas - 1 GROWTH RATE 2DENSITY-DEPENDENCE
    for(i in 1:S){
      alpha0[i] ~ dnorm(mualpha0,tau0)
      alpha1[i] ~ dnorm(0,tauDD)
      #alpha0[i] ~ dnorm(mualpha0,tau0)T(0,5)
      #alpha1[i] ~ dnorm(0,tauDD)T(-2,2)
    }
    ## PRIORS AND HYPERPARAMETERS
    mualpha0 ~ dnorm(0,0.01)
    tau0 <- pow(sd0,-2)
    sd0 ~ dunif(0,2)
    tauDD <- pow(sigDD,-2)
    sigDD ~ dunif(0,2)
    #
    ## Random (synchrony) terms or shared variance
    for(j in 1:T){
      delta[j]~dnorm(0,deltaprec)
    }
    deltaprec~dgamma(1.0E-1,1.0E-1)
    ## asynchrony or unshared term
    for(i in 1:S){
      epsilonprec[i]~dgamma(1.0E-1,1.0E-1)  #dgamma(1.0,1.0E-1)
    }
    ## INITIAL COUNTS
    for(i in 1:S){
      tau1[i] <- pow(sdy1[i],-2)
      y0[i] ~ dnorm(y1[i],tau1[i])
      x[1,i] <- log(max(y0[i],1.0E-5))
    }
    ## LIKELIHOOD - SYSTEM PROCESS
    for(i in 1:S){
      for(j in 2:(T+1)){
        #Density-dependence and synchrony component
        Ex[j,i] <- alpha0[i]+(1+alpha1[i])*x[j-1,i]+delta[j-1]
        #Asynchrony component
        x[j,i] ~ dnorm(Ex[j,i],epsilonprec[i])
      }
    }
    #  OBSERVATION PROCESS (On real scale and Negative binomial observation errors)
    ## GENERATE PRIORS FOR r, GET p, AND SAMPLE ABUNDANCE 
    for(i in 1:S){
      r[i] ~ dunif(0,50)
      for(j in 1:T){
        log(mu[j,i]) <- x[j,i]
        p[j,i] <- r[i]/(r[i]+mu[j,i])
        y[j,i] ~ dnegbin(p[j,i],r[i])
      }
    }
    ## Derived variances and ICC
    varshared <- 1/deltaprec
    for(i in 1:S){
      varunshared[i] <- 1/epsilonprec[i]
      ICC[i] <- varshared/(varshared+varunshared[i])
    }
    ICCall <- varshared/(varshared+sum(varunshared[]))      
  }
  ", fill = TRUE)
sink()
##
## Data
y <- cbind(round(as.vector(colMeans(NTSims[,1:46]))),as.vector(AFSRes$Afs[4:49]))
jags.data=list(y=y,y1=y[1,],sdy1=c(sd(NTSims[,1]),y[1,1]*.05),T=dim(y)[1],S=dim(y)[2])
##
##
## Initial values
inits <- function() list(deltaprec=0.1,epsilonprec=rep(0.01,2))
##
## Parameters monitored
params <- c("r","alpha0","mualpha0","alpha1","deltaprec","epsilonprec","varshared","varunshared","ICC","ICCall","p")
##
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 200000
nt <- 20
nb <- 100000
nc <- 4
## Call jags and fit model
(sta <- Sys.time())
out31 <- jags(jags.data,inits,params,"AFSNB31.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
## Time difference of 1.588677 mins
#
## GENERATE RESULTS
outT <- out31
outTdf <- data.frame(par=c("rBI","rBI","mups","upsbi","upssi","nubi","nusi","sigd","sigebi","sigesi","ICCbi","ICCsi","ICC"),
                     mean=c(outT$mean$r,outT$mean$mualpha0,outT$mean$alpha0,outT$mean$alpha1,outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
                     sd=c(outT$sd$r,outT$sd$mualpha0,outT$sd$alpha0,outT$sd$alpha1,outT$sd$varshared,outT$sd$varunshared,outT$sd$ICC,outT$sd$ICCall))
ouTtable$lag1 <- paste(round(outTdf$mean,2),paste(rep("(",dim(outTdf)[1]),round(outTdf$sd,2),rep(")",dim(outTdf)[1]),sep=""))
##
predsT <- rbind(outT$samples[[1]],outT$samples[[2]],outT$samples[[3]])
varnm <-c("varshared","varunshared[1]","varunshared[2]","ICC[1]","ICC[2]","ICCall")
quant <- quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[1])],c(0.05,0.95))
for(i in 2:6)
  quant <- rbind(quant,quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[i])],c(0.05,0.95)))

ICCtable <- rbind(ICCtable,data.frame(par=c("sigdelta","sigepsiBI","sigepsiSI","ICCBI","ICCSI","ICC"),
                       mean=c(outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
                       lowCI=quant[,1],highCI=quant[,2],lag=rep("lag1",6)))
## and clean up
rm(quant,predsT,outT,outTdf)
##--------------------------
##  4.3.3 NEGATIVE BINOMIAL BUGS MODEL WITH SI,BI SYNCHRONY DYNAMICS - SI lagged -2
### BUGS MODEL
sink("AFSNB32.jags")
cat("
  model{
    ## Alphas - 1 GROWTH RATE 2DENSITY-DEPENDENCE
    for(i in 1:S){
      alpha0[i] ~ dnorm(mualpha0,tau0)
      alpha1[i] ~ dnorm(0,tauDD)
      #alpha0[i] ~ dnorm(mualpha0,tau0)T(0,5)
      #alpha1[i] ~ dnorm(0,tauDD)T(-2,2)
    }
    ## PRIORS AND HYPERPARAMETERS
    mualpha0 ~ dnorm(0,0.01)
    tau0 <- pow(sd0,-2)
    sd0 ~ dunif(0,2)
    tauDD <- pow(sigDD,-2)
    sigDD ~ dunif(0,2)
    #
    ## Random (synchrony) terms or shared variance
    for(j in 1:T){
      delta[j]~dnorm(0,deltaprec)
    }
    deltaprec~dgamma(1.0E-1,1.0E-1)
    ## asynchrony or unshared term
    for(i in 1:S){
      epsilonprec[i]~dgamma(1.0E-1,1.0E-1)  #dgamma(1.0,1.0E-1)
    }
    ## INITIAL COUNTS
    for(i in 1:S){
      tau1[i] <- pow(sdy1[i],-2)
      y0[i] ~ dnorm(y1[i],tau1[i])
      x[1,i] <- log(max(y0[i],1.0E-5))
    }
    ## LIKELIHOOD - SYSTEM PROCESS
    for(i in 1:S){
      for(j in 2:(T+1)){
        #Density-dependence and synchrony component
        Ex[j,i] <- alpha0[i]+(1+alpha1[i])*x[j-1,i]+delta[j-1]
        #Asynchrony component
        x[j,i] ~ dnorm(Ex[j,i],epsilonprec[i])
      }
    }
    #  OBSERVATION PROCESS (On real scale and Negative binomial observation errors)
    ## GENERATE PRIORS FOR r, GET p, AND SAMPLE ABUNDANCE 
    for(i in 1:S){
      r[i] ~ dunif(0,50)
      for(j in 1:T){
        log(mu[j,i]) <- x[j,i]
        p[j,i] <- r[i]/(r[i]+mu[j,i])
        y[j,i] ~ dnegbin(p[j,i],r[i])
      }
    }
    ## Derived variances and ICC
    varshared <- 1/deltaprec
    for(i in 1:S){
      varunshared[i] <- 1/epsilonprec[i]
      ICC[i] <- varshared/(varshared+varunshared[i])
    }
    ICCall <- varshared/(varshared+sum(varunshared[]))      
  }
  ", fill = TRUE)
sink()
##
## Data
# cor.test(AFSRes$Afs[5:48],as.vector(colMeans(NTSims[,1:44])))$estimate
y <- cbind(round(as.vector(colMeans(NTSims[,1:45]))),as.vector(AFSRes$Afs[5:49]))
jags.data=list(y=y,y1=y[1,],sdy1=c(sd(NTSims[,1]),y[1,1]*.05),T=dim(y)[1],S=dim(y)[2])
##
##
## Initial values
inits <- function() list(deltaprec=0.1,epsilonprec=rep(0.01,2))
##
## Parameters monitored
params <- c("r","alpha0","mualpha0","alpha1","deltaprec","epsilonprec","varshared","varunshared","ICC","ICCall","p")
##
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 200000
nt <- 20
nb <- 100000
nc <- 4
## Call jags and fit model
(sta <- Sys.time())
out32 <- jags(jags.data,inits,params,"AFSNB32.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
## Time difference of 1.86795 mins
#
## GENERATE RESULTS
outT <- out32
outTdf <- data.frame(par=c("rBI","rBI","mups","upsbi","upssi","nubi","nusi","sigd","sigebi","sigesi","ICCbi","ICCsi","ICC"),
                     mean=c(outT$mean$r,outT$mean$mualpha0,outT$mean$alpha0,outT$mean$alpha1,outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
                     sd=c(outT$sd$r,outT$sd$mualpha0,outT$sd$alpha0,outT$sd$alpha1,outT$sd$varshared,outT$sd$varunshared,outT$sd$ICC,outT$sd$ICCall))
ouTtable$lag2 <- paste(round(outTdf$mean,2),paste(rep("(",dim(outTdf)[1]),round(outTdf$sd,2),rep(")",dim(outTdf)[1]),sep=""))
##
predsT <- rbind(outT$samples[[1]],outT$samples[[2]],outT$samples[[3]])
varnm <-c("varshared","varunshared[1]","varunshared[2]","ICC[1]","ICC[2]","ICCall")
quant <- quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[1])],c(0.05,0.95))
for(i in 2:6)
  quant <- rbind(quant,quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[i])],c(0.05,0.95)))

ICCtable <- rbind(ICCtable,data.frame(par=c("sigdelta","sigepsiBI","sigepsiSI","ICCBI","ICCSI","ICC"),
                                      mean=c(outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
                                      lowCI=quant[,1],highCI=quant[,2],lag=rep("lag2",6)))
## and clean up
rm(quant,predsT,outT,outTdf)
##--------------------------
##  4.3.4 NEGATIVE BINOMIAL BUGS MODEL WITH SI,BI SYNCHRONY DYNAMICS - SI lagged -3
### BUGS MODEL
sink("AFSNB33.jags")
cat("
  model{
    ## Alphas - 1 GROWTH RATE 2DENSITY-DEPENDENCE
    for(i in 1:S){
      alpha0[i] ~ dnorm(mualpha0,tau0)
      alpha1[i] ~ dnorm(0,tauDD)
      #alpha0[i] ~ dnorm(mualpha0,tau0)T(0,5)
      #alpha1[i] ~ dnorm(0,tauDD)T(-2,2)
    }
    ## PRIORS AND HYPERPARAMETERS
    mualpha0 ~ dnorm(0,0.01)
    tau0 <- pow(sd0,-2)
    sd0 ~ dunif(0,2)
    tauDD <- pow(sigDD,-2)
    sigDD ~ dunif(0,2)
    #
    ## Random (synchrony) terms or shared variance
    for(j in 1:T){
      delta[j]~dnorm(0,deltaprec)
    }
    deltaprec~dgamma(1.0E-1,1.0E-1)
    ## asynchrony or unshared term
    for(i in 1:S){
      epsilonprec[i]~dgamma(1.0E-1,1.0E-1)  #dgamma(1.0,1.0E-1)
    }
    ## INITIAL COUNTS
    for(i in 1:S){
      tau1[i] <- pow(sdy1[i],-2)
      y0[i] ~ dnorm(y1[i],tau1[i])
      x[1,i] <- log(max(y0[i],1.0E-5))
    }
    ## LIKELIHOOD - SYSTEM PROCESS
    for(i in 1:S){
      for(j in 2:(T+1)){
        #Density-dependence and synchrony component
        Ex[j,i] <- alpha0[i]+(1+alpha1[i])*x[j-1,i]+delta[j-1]
        #Asynchrony component
        x[j,i] ~ dnorm(Ex[j,i],epsilonprec[i])
      }
    }
    #  OBSERVATION PROCESS (On real scale and Negative binomial observation errors)
    ## GENERATE PRIORS FOR r, GET p, AND SAMPLE ABUNDANCE 
    for(i in 1:S){
      r[i] ~ dunif(0,50)
      for(j in 1:T){
        log(mu[j,i]) <- x[j,i]
        p[j,i] <- r[i]/(r[i]+mu[j,i])
        y[j,i] ~ dnegbin(p[j,i],r[i])
      }
    }
    ## Derived variances and ICC
    varshared <- 1/deltaprec
    for(i in 1:S){
      varunshared[i] <- 1/epsilonprec[i]
      ICC[i] <- varshared/(varshared+varunshared[i])
    }
    ICCall <- varshared/(varshared+sum(varunshared[]))      
  }
  ", fill = TRUE)
sink()
##
## Data
# cor.test(AFSRes$Afs[6:48],as.vector(colMeans(NTSims[,1:43])))$estimate
y <- cbind(round(as.vector(colMeans(NTSims[,1:44]))),as.vector(AFSRes$Afs[6:49]))
jags.data=list(y=y,y1=y[1,],sdy1=c(sd(NTSims[,1]),y[1,1]*.05),T=dim(y)[1],S=dim(y)[2])
##
##
## Initial values
inits <- function() list(deltaprec=0.1,epsilonprec=rep(0.01,2))
##
## Parameters monitored
params <- c("r","alpha0","mualpha0","alpha1","deltaprec","epsilonprec","varshared","varunshared","ICC","ICCall","p")
##
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 200000
nt <- 20
nb <- 100000
nc <- 4
## Call jags and fit model
(sta <- Sys.time())
out33 <- jags(jags.data,inits,params,"AFSNB33.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
## Time difference of 1.700351 mins
#
## GENERATE RESULTS
outT <- out33
outTdf <- data.frame(par=c("rBI","rBI","mups","upsbi","upssi","nubi","nusi","sigd","sigebi","sigesi","ICCbi","ICCsi","ICC"),
                     mean=c(outT$mean$r,outT$mean$mualpha0,outT$mean$alpha0,outT$mean$alpha1,outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
                     sd=c(outT$sd$r,outT$sd$mualpha0,outT$sd$alpha0,outT$sd$alpha1,outT$sd$varshared,outT$sd$varunshared,outT$sd$ICC,outT$sd$ICCall))
ouTtable$lag3 <- paste(round(outTdf$mean,2),paste(rep("(",dim(outTdf)[1]),round(outTdf$sd,2),rep(")",dim(outTdf)[1]),sep=""))
##
predsT <- rbind(outT$samples[[1]],outT$samples[[2]],outT$samples[[3]])
varnm <-c("varshared","varunshared[1]","varunshared[2]","ICC[1]","ICC[2]","ICCall")
quant <- quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[1])],c(0.05,0.95))
for(i in 2:6)
  quant <- rbind(quant,quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[i])],c(0.05,0.95)))

ICCtable <- rbind(ICCtable,data.frame(par=c("sigdelta","sigepsiBI","sigepsiSI","ICCBI","ICCSI","ICC"),
                                      mean=c(outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
                                      lowCI=quant[,1],highCI=quant[,2],lag=rep("lag3",6)))
## and clean up
rm(quant,predsT,outT,outTdf)
##--------------------------
##  4.3.5 NEGATIVE BINOMIAL BUGS MODEL WITH SI,BI SYNCHRONY DYNAMICS - SI lagged -4
### BUGS MODEL
sink("AFSNB34.jags")
cat("
  model{
    ## Alphas - 1 GROWTH RATE 2DENSITY-DEPENDENCE
    for(i in 1:S){
      alpha0[i] ~ dnorm(mualpha0,tau0)
      alpha1[i] ~ dnorm(0,tauDD)
      #alpha0[i] ~ dnorm(mualpha0,tau0)T(0,5)
      #alpha1[i] ~ dnorm(0,tauDD)T(-2,2)
    }
    ## PRIORS AND HYPERPARAMETERS
    mualpha0 ~ dnorm(0,0.01)
    tau0 <- pow(sd0,-2)
    sd0 ~ dunif(0,2)
    tauDD <- pow(sigDD,-2)
    sigDD ~ dunif(0,2)
    #
    ## Random (synchrony) terms or shared variance
    for(j in 1:T){
      delta[j]~dnorm(0,deltaprec)
    }
    deltaprec~dgamma(1.0E-1,1.0E-1)
    ## asynchrony or unshared term
    for(i in 1:S){
      epsilonprec[i]~dgamma(1.0E-1,1.0E-1)  #dgamma(1.0,1.0E-1)
    }
    ## INITIAL COUNTS
    for(i in 1:S){
      tau1[i] <- pow(sdy1[i],-2)
      y0[i] ~ dnorm(y1[i],tau1[i])
      x[1,i] <- log(max(y0[i],1.0E-5))
    }
    ## LIKELIHOOD - SYSTEM PROCESS
    for(i in 1:S){
      for(j in 2:(T+1)){
        #Density-dependence and synchrony component
        Ex[j,i] <- alpha0[i]+(1+alpha1[i])*x[j-1,i]+delta[j-1]
        #Asynchrony component
        x[j,i] ~ dnorm(Ex[j,i],epsilonprec[i])
      }
    }
    #  OBSERVATION PROCESS (On real scale and Negative binomial observation errors)
    ## GENERATE PRIORS FOR r, GET p, AND SAMPLE ABUNDANCE 
    for(i in 1:S){
      r[i] ~ dunif(0,50)
      for(j in 1:T){
        log(mu[j,i]) <- x[j,i]
        p[j,i] <- r[i]/(r[i]+mu[j,i])
        y[j,i] ~ dnegbin(p[j,i],r[i])
      }
    }
    ## Derived variances and ICC
    varshared <- 1/deltaprec
    for(i in 1:S){
      varunshared[i] <- 1/epsilonprec[i]
      ICC[i] <- varshared/(varshared+varunshared[i])
    }
    ICCall <- varshared/(varshared+sum(varunshared[]))      
  }
  ", fill = TRUE)
sink()
##
## Data
# cor.test(AFSRes$Afs[7:48],as.vector(colMeans(NTSims[,1:42])))$estimate
y <- cbind(round(as.vector(colMeans(NTSims[,1:43]))),as.vector(AFSRes$Afs[7:49]))
jags.data=list(y=y,y1=y[1,],sdy1=c(sd(NTSims[,1]),y[1,1]*.05),T=dim(y)[1],S=dim(y)[2])
##
##
## Initial values
inits <- function() list(deltaprec=0.1,epsilonprec=rep(0.01,2))
##
## Parameters monitored
params <- c("r","alpha0","mualpha0","alpha1","deltaprec","epsilonprec","varshared","varunshared","ICC","ICCall","p")
##
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 200000
nt <- 20
nb <- 100000
nc <- 4
## Call jags and fit model
(sta <- Sys.time())
out34 <- jags(jags.data,inits,params,"AFSNB34.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
## Time difference of 1.700351 mins
#
## GENERATE RESULTS
outT <- out34
outTdf <- data.frame(par=c("rBI","rBI","mups","upsbi","upssi","nubi","nusi","sigd","sigebi","sigesi","ICCbi","ICCsi","ICC"),
                     mean=c(outT$mean$r,outT$mean$mualpha0,outT$mean$alpha0,outT$mean$alpha1,outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
                     sd=c(outT$sd$r,outT$sd$mualpha0,outT$sd$alpha0,outT$sd$alpha1,outT$sd$varshared,outT$sd$varunshared,outT$sd$ICC,outT$sd$ICCall))
ouTtable$lag4 <- paste(round(outTdf$mean,2),paste(rep("(",dim(outTdf)[1]),round(outTdf$sd,2),rep(")",dim(outTdf)[1]),sep=""))
##
predsT <- rbind(outT$samples[[1]],outT$samples[[2]],outT$samples[[3]])
varnm <-c("varshared","varunshared[1]","varunshared[2]","ICC[1]","ICC[2]","ICCall")
quant <- quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[1])],c(0.05,0.95))
for(i in 2:6)
  quant <- rbind(quant,quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[i])],c(0.05,0.95)))

ICCtable <- rbind(ICCtable,data.frame(par=c("sigdelta","sigepsiBI","sigepsiSI","ICCBI","ICCSI","ICC"),
  mean=c(outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
  lowCI=quant[,1],highCI=quant[,2],lag=rep("lag4",6)))
## and clean up
rm(quant,predsT,outT,outTdf)
##--------------------------
##  4.3.6 NEGATIVE BINOMIAL BUGS MODEL WITH SI,BI SYNCHRONY DYNAMICS - SI lagged -5
### BUGS MODEL
sink("AFSNB35.jags")
cat("
  model{
    ## Alphas - 1 GROWTH RATE 2DENSITY-DEPENDENCE
    for(i in 1:S){
      alpha0[i] ~ dnorm(mualpha0,tau0)
      alpha1[i] ~ dnorm(0,tauDD)
      #alpha0[i] ~ dnorm(mualpha0,tau0)T(0,5)
      #alpha1[i] ~ dnorm(0,tauDD)T(-2,2)
    }
    ## PRIORS AND HYPERPARAMETERS
    mualpha0 ~ dnorm(0,0.01)
    tau0 <- pow(sd0,-2)
    sd0 ~ dunif(0,2)
    tauDD <- pow(sigDD,-2)
    sigDD ~ dunif(0,2)
    #
    ## Random (synchrony) terms or shared variance
    for(j in 1:T){
      delta[j]~dnorm(0,deltaprec)
    }
    deltaprec~dgamma(1.0E-1,1.0E-1)
    ## asynchrony or unshared term
    for(i in 1:S){
      epsilonprec[i]~dgamma(1.0E-1,1.0E-1)  #dgamma(1.0,1.0E-1)
    }
    ## INITIAL COUNTS
    for(i in 1:S){
      tau1[i] <- pow(sdy1[i],-2)
      y0[i] ~ dnorm(y1[i],tau1[i])
      x[1,i] <- log(max(y0[i],1.0E-5))
    }
    ## LIKELIHOOD - SYSTEM PROCESS
    for(i in 1:S){
      for(j in 2:(T+1)){
        #Density-dependence and synchrony component
        Ex[j,i] <- alpha0[i]+(1+alpha1[i])*x[j-1,i]+delta[j-1]
        #Asynchrony component
        x[j,i] ~ dnorm(Ex[j,i],epsilonprec[i])
      }
    }
    #  OBSERVATION PROCESS (On real scale and Negative binomial observation errors)
    ## GENERATE PRIORS FOR r, GET p, AND SAMPLE ABUNDANCE 
    for(i in 1:S){
      r[i] ~ dunif(0,50)
      for(j in 1:T){
        log(mu[j,i]) <- x[j,i]
        p[j,i] <- r[i]/(r[i]+mu[j,i])
        y[j,i] ~ dnegbin(p[j,i],r[i])
      }
    }
    ## Derived variances and ICC
    varshared <- 1/deltaprec
    for(i in 1:S){
      varunshared[i] <- 1/epsilonprec[i]
      ICC[i] <- varshared/(varshared+varunshared[i])
    }
    ICCall <- varshared/(varshared+sum(varunshared[]))      
  }
  ", fill = TRUE)
sink()
##
## Data
# cor.test(AFSRes$Afs[8:48],as.vector(colMeans(NTSims[,1:41])))$estimate
y <- cbind(round(as.vector(colMeans(NTSims[,1:42]))),as.vector(AFSRes$Afs[8:49]))
jags.data=list(y=y,y1=y[1,],sdy1=c(sd(NTSims[,1]),y[1,1]*.05),T=dim(y)[1],S=dim(y)[2])
##
##
## Initial values
inits <- function() list(deltaprec=0.1,epsilonprec=rep(0.01,2))
##
## Parameters monitored
params <- c("r","alpha0","mualpha0","alpha1","deltaprec","epsilonprec","varshared","varunshared","ICC","ICCall","p")
##
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 200000
nt <- 20
nb <- 100000
nc <- 4
## Call jags and fit model
(sta <- Sys.time())
out35 <- jags(jags.data,inits,params,"AFSNB35.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
## Time difference of 1.700351 mins
#
## GENERATE RESULTS
outT <- out35
outTdf <- data.frame(par=c("rBI","rBI","mups","upsbi","upssi","nubi","nusi","sigd","sigebi","sigesi","ICCbi","ICCsi","ICC"),
                     mean=c(outT$mean$r,outT$mean$mualpha0,outT$mean$alpha0,outT$mean$alpha1,outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
                     sd=c(outT$sd$r,outT$sd$mualpha0,outT$sd$alpha0,outT$sd$alpha1,outT$sd$varshared,outT$sd$varunshared,outT$sd$ICC,outT$sd$ICCall))
ouTtable$lag5 <- paste(round(outTdf$mean,2),paste(rep("(",dim(outTdf)[1]),round(outTdf$sd,2),rep(")",dim(outTdf)[1]),sep=""))
##
predsT <- rbind(outT$samples[[1]],outT$samples[[2]],outT$samples[[3]])
varnm <-c("varshared","varunshared[1]","varunshared[2]","ICC[1]","ICC[2]","ICCall")
quant <- quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[1])],c(0.05,0.95))
for(i in 2:6)
  quant <- rbind(quant,quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[i])],c(0.05,0.95)))

ICCtable <- rbind(ICCtable,data.frame(par=c("sigdelta","sigepsiBI","sigepsiSI","ICCBI","ICCSI","ICC"),
  mean=c(outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
  lowCI=quant[,1],highCI=quant[,2],lag=rep("lag5",6)))
## and clean up
rm(quant,predsT,outT,outTdf)
##--------------------------
##  4.3.7 NEGATIVE BINOMIAL BUGS MODEL WITH SI,BI SYNCHRONY DYNAMICS - SI lagged -6
### BUGS MODEL
sink("AFSNB36.jags")
cat("
  model{
    ## Alphas - 1 GROWTH RATE 2DENSITY-DEPENDENCE
    for(i in 1:S){
      alpha0[i] ~ dnorm(mualpha0,tau0)
      alpha1[i] ~ dnorm(0,tauDD)
      #alpha0[i] ~ dnorm(mualpha0,tau0)T(0,5)
      #alpha1[i] ~ dnorm(0,tauDD)T(-2,2)
    }
    ## PRIORS AND HYPERPARAMETERS
    mualpha0 ~ dnorm(0,0.01)
    tau0 <- pow(sd0,-2)
    sd0 ~ dunif(0,2)
    tauDD <- pow(sigDD,-2)
    sigDD ~ dunif(0,2)
    #
    ## Random (synchrony) terms or shared variance
    for(j in 1:T){
      delta[j]~dnorm(0,deltaprec)
    }
    deltaprec~dgamma(1.0E-1,1.0E-1)
    ## asynchrony or unshared term
    for(i in 1:S){
      epsilonprec[i]~dgamma(1.0E-1,1.0E-1)  #dgamma(1.0,1.0E-1)
    }
    ## INITIAL COUNTS
    for(i in 1:S){
      tau1[i] <- pow(sdy1[i],-2)
      y0[i] ~ dnorm(y1[i],tau1[i])
      x[1,i] <- log(max(y0[i],1.0E-5))
    }
    ## LIKELIHOOD - SYSTEM PROCESS
    for(i in 1:S){
      for(j in 2:(T+1)){
        #Density-dependence and synchrony component
        Ex[j,i] <- alpha0[i]+(1+alpha1[i])*x[j-1,i]+delta[j-1]
        #Asynchrony component
        x[j,i] ~ dnorm(Ex[j,i],epsilonprec[i])
      }
    }
    #  OBSERVATION PROCESS (On real scale and Negative binomial observation errors)
    ## GENERATE PRIORS FOR r, GET p, AND SAMPLE ABUNDANCE 
    for(i in 1:S){
      r[i] ~ dunif(0,50)
      for(j in 1:T){
        log(mu[j,i]) <- x[j,i]
        p[j,i] <- r[i]/(r[i]+mu[j,i])
        y[j,i] ~ dnegbin(p[j,i],r[i])
      }
    }
    ## Derived variances and ICC
    varshared <- 1/deltaprec
    for(i in 1:S){
      varunshared[i] <- 1/epsilonprec[i]
      ICC[i] <- varshared/(varshared+varunshared[i])
    }
    ICCall <- varshared/(varshared+sum(varunshared[]))      
  }
  ", fill = TRUE)
sink()
##
## Data
# cor.test(AFSRes$Afs[9:48],as.vector(colMeans(NTSims[,1:40])))$estimate
y <- cbind(round(as.vector(colMeans(NTSims[,1:41]))),as.vector(AFSRes$Afs[9:49]))
jags.data=list(y=y,y1=y[1,],sdy1=c(sd(NTSims[,1]),y[1,1]*.05),T=dim(y)[1],S=dim(y)[2])
##
##
## Initial values
inits <- function() list(deltaprec=0.1,epsilonprec=rep(0.01,2))
##
## Parameters monitored
params <- c("r","alpha0","mualpha0","alpha1","deltaprec","epsilonprec","varshared","varunshared","ICC","ICCall","p")
##
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 200000
nt <- 20
nb <- 100000
nc <- 4
## Call jags and fit model
(sta <- Sys.time())
out36 <- jags(jags.data,inits,params,"AFSNB36.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
## Time difference of 1.447232 mins
#
## GENERATE RESULTS
outT <- out36
outTdf <- data.frame(par=c("rBI","rBI","mups","upsbi","upssi","nubi","nusi","sigd","sigebi","sigesi","ICCbi","ICCsi","ICC"),
                     mean=c(outT$mean$r,outT$mean$mualpha0,outT$mean$alpha0,outT$mean$alpha1,outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
                     sd=c(outT$sd$r,outT$sd$mualpha0,outT$sd$alpha0,outT$sd$alpha1,outT$sd$varshared,outT$sd$varunshared,outT$sd$ICC,outT$sd$ICCall))
ouTtable$lag6 <- paste(round(outTdf$mean,2),paste(rep("(",dim(outTdf)[1]),round(outTdf$sd,2),rep(")",dim(outTdf)[1]),sep=""))
##
predsT <- rbind(outT$samples[[1]],outT$samples[[2]],outT$samples[[3]])
varnm <-c("varshared","varunshared[1]","varunshared[2]","ICC[1]","ICC[2]","ICCall")
quant <- quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[1])],c(0.05,0.95))
for(i in 2:6)
  quant <- rbind(quant,quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[i])],c(0.05,0.95)))

ICCtable <- rbind(ICCtable,data.frame(par=c("sigdelta","sigepsiBI","sigepsiSI","ICCBI","ICCSI","ICC"),
  mean=c(outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
  lowCI=quant[,1],highCI=quant[,2],lag=rep("lag6",6)))
## and clean up
rm(quant,predsT,outT,outTdf)
##--------------------------
##  4.3.8 NEGATIVE BINOMIAL BUGS MODEL WITH SI,BI SYNCHRONY DYNAMICS - SI lagged -7
### BUGS MODEL
sink("AFSNB37.jags")
cat("
  model{
    ## Alphas - 1 GROWTH RATE 2DENSITY-DEPENDENCE
    for(i in 1:S){
      alpha0[i] ~ dnorm(mualpha0,tau0)
      alpha1[i] ~ dnorm(0,tauDD)
      #alpha0[i] ~ dnorm(mualpha0,tau0)T(0,5)
      #alpha1[i] ~ dnorm(0,tauDD)T(-2,2)
    }
    ## PRIORS AND HYPERPARAMETERS
    mualpha0 ~ dnorm(0,0.01)
    tau0 <- pow(sd0,-2)
    sd0 ~ dunif(0,2)
    tauDD <- pow(sigDD,-2)
    sigDD ~ dunif(0,2)
    #
    ## Random (synchrony) terms or shared variance
    for(j in 1:T){
      delta[j]~dnorm(0,deltaprec)
    }
    deltaprec~dgamma(1.0E-1,1.0E-1)
    ## asynchrony or unshared term
    for(i in 1:S){
      epsilonprec[i]~dgamma(1.0E-1,1.0E-1)  #dgamma(1.0,1.0E-1)
    }
    ## INITIAL COUNTS
    for(i in 1:S){
      tau1[i] <- pow(sdy1[i],-2)
      y0[i] ~ dnorm(y1[i],tau1[i])
      x[1,i] <- log(max(y0[i],1.0E-5))
    }
    ## LIKELIHOOD - SYSTEM PROCESS
    for(i in 1:S){
      for(j in 2:(T+1)){
        #Density-dependence and synchrony component
        Ex[j,i] <- alpha0[i]+(1+alpha1[i])*x[j-1,i]+delta[j-1]
        #Asynchrony component
        x[j,i] ~ dnorm(Ex[j,i],epsilonprec[i])
      }
    }
    #  OBSERVATION PROCESS (On real scale and Negative binomial observation errors)
    ## GENERATE PRIORS FOR r, GET p, AND SAMPLE ABUNDANCE 
    for(i in 1:S){
      r[i] ~ dunif(0,50)
      for(j in 1:T){
        log(mu[j,i]) <- x[j,i]
        p[j,i] <- r[i]/(r[i]+mu[j,i])
        y[j,i] ~ dnegbin(p[j,i],r[i])
      }
    }
    ## Derived variances and ICC
    varshared <- 1/deltaprec
    for(i in 1:S){
      varunshared[i] <- 1/epsilonprec[i]
      ICC[i] <- varshared/(varshared+varunshared[i])
    }
    ICCall <- varshared/(varshared+sum(varunshared[]))      
  }
  ", fill = TRUE)
sink()
##
## Data
# cor.test(AFSRes$Afs[10:48],as.vector(colMeans(NTSims[,1:39])))$estimate
y <- cbind(round(as.vector(colMeans(NTSims[,1:40]))),as.vector(AFSRes$Afs[10:49]))
jags.data=list(y=y,y1=y[1,],sdy1=c(sd(NTSims[,1]),y[1,1]*.05),T=dim(y)[1],S=dim(y)[2])
##
##
## Initial values
inits <- function() list(deltaprec=0.1,epsilonprec=rep(0.01,2))
##
## Parameters monitored
params <- c("r","alpha0","mualpha0","alpha1","deltaprec","epsilonprec","varshared","varunshared","ICC","ICCall","p")
##
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 350000
nt <- 25
nb <- 250000
nc <- 4
## Call jags and fit model
(sta <- Sys.time())
out37 <- jags(jags.data,inits,params,"AFSNB37.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
## Time difference of 1.447232 mins
#
## GENERATE RESULTS
outT <- out37
outTdf <- data.frame(par=c("rBI","rBI","mups","upsbi","upssi","nubi","nusi","sigd","sigebi","sigesi","ICCbi","ICCsi","ICC"),
  mean=c(outT$mean$r,outT$mean$mualpha0,outT$mean$alpha0,outT$mean$alpha1,outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
  sd=c(outT$sd$r,outT$sd$mualpha0,outT$sd$alpha0,outT$sd$alpha1,outT$sd$varshared,outT$sd$varunshared,outT$sd$ICC,outT$sd$ICCall))
ouTtable$lag7 <- paste(round(outTdf$mean,2),paste(rep("(",dim(outTdf)[1]),round(outTdf$sd,2),rep(")",dim(outTdf)[1]),sep=""))
##
predsT <- rbind(outT$samples[[1]],outT$samples[[2]],outT$samples[[3]])
varnm <-c("varshared","varunshared[1]","varunshared[2]","ICC[1]","ICC[2]","ICCall")
quant <- quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[1])],c(0.05,0.95))
for(i in 2:6)
  quant <- rbind(quant,quantile(predsT[,which(dimnames(predsT)[[2]]==varnm[i])],c(0.05,0.95)))

ICCtable <- rbind(ICCtable,data.frame(par=c("sigdelta","sigepsiBI","sigepsiSI","ICCBI","ICCSI","ICC"),
  mean=c(outT$mean$varshared,outT$mean$varunshared,outT$mean$ICC,outT$mean$ICCall),
  lowCI=quant[,1],highCI=quant[,2],lag=rep("lag7",6)))
## and clean up
rm(quant,predsT,outT,outTdf)
#
## SAVE RESULT FILES
saveRDS(ouTtable,"NBMSynchTable.Rds")
saveRDS(ICCtable,"ICCtable.Rds")
##
ICCdata <- data.frame(Lag=rep(0:-7,2),rbind(ICCtable[which(ICCtable$par=="ICCBI"),2:4],ICCtable[which(ICCtable$par=="ICCSI"),2:4]),Population=rep(c("Bird Island","Signy Island"),each=8))
names(ICCdata)[2] <- "ICC"

ggplot(ICCdata,aes(Lag,ICC))+
  geom_pointrange(aes(ymin=lowCI,ymax=highCI,color=Population),position=position_dodge(0.3))+
  scale_color_manual(values = c("wheat4", "lightblue"))+
  theme_bw()+
  labs(y=expression("Synchrony"~"("*ICC[P]*")"),x="lag (years)")+
  theme(panel.border=element_rect(color="black",fill=NA,linewidth=1),
        legend.position="top",legend.box.margin=margin(-10,-10,-10,-10),legend.margin=margin(0,0,0,0))+
  scale_x_continuous(breaks=seq(-7,0,by=1),limits=c(-7.25,0.25))+
  scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1),minor_breaks=seq(0,1,0.1))
  

##----------------------------
AFSRes <- readRDS("AFStrends.Rds")
plot(AFSRes$year,AFSRes$Afs,type="n",xlim=c(1977,2026),ylim=c(1200,37000),xlab="",ylab="",axes=F) #ylim=c(1200,26250)
polygon(c(AFSRes$year,rev(AFSRes$year)),c(AFSRes$ANBH,rev(AFSRes$ANBL)),col=addTrans("lightgray",150),border=F)
lines(AFSRes$year,AFSRes$ANBM,lwd=1,col="grey")
points(AFSRes$year,AFSRes$Afs,pch=21,col="blue",bg="lightblue",cex=.75)
box()
axis(1,at=seq(1980,2025,by=5),labels=T,tck=-0.025,cex.axis=1,lty=1,padj=-1.0)
axis(1,at=seq(1977,2025,by=1),labels=F,tck=-0.012,las=2,cex.axis=1)
axis(4,at=pretty(c(1200,26250))[1:6],labels=T,tck=-0.02,las=2,cex.axis=.75,col.axis="blue",hadj=0.5)  # col.ticks="blue", hadj=0.7
axis(4,at=seq(0,27000,by=1000),labels=F,tck=-0.01,las=2,cex.axis=.75,col.ticks="blue")
par(new=T)
##
lpAFS <- readRDS("lpAFS.Rds")
gamAFS <- readRDS("gamAFS.Rds")
NTSims <- readRDS("NTSims.Rds")
##
plot(1979:2025,colMeans(NTSims),type="n",pch=21,bg=2,xlim=c(1977,2026),ylim=c(100,1200),axes=F,ylab="",xlab="")
polygon(c(1979:2025,2025:1979)+0.3, c(apply(NTSims,2,quantile,0.975),rev(apply(NTSims,2,quantile,0.025))),col=addTrans("wheat4",100),border="wheat4")
lines(1979:2025+0.3,colMeans(NTSims),col="wheat4", lwd = 1.1) 
points(1979:2025+0.3,colMeans(NTSims),type ="p",pch=21,col="black",bg="wheat4",cex=.75)
box()
axis(2,at=seq(400,1200,by=200),labels=T,tck=-0.02,las=2,cex.axis=.75,col.axis="black",col.ticks="wheat4",hadj=0.4)  # col.ticks="blue", hadj=0.7
axis(2,at=seq(300,1200,by=100),labels=F,tck=-0.01,las=2,cex.axis=.75,col.ticks="wheat4")
mtext("Signy Island population counts",side=4,line=2.25,cex=.85,col="blue",adj=0.125)
mtext("Bird Island mature females",side=2,line=1.9,cex=.85,col="black",adj=0.75)
##
par(new=F)
##--------------------------------------------------------------------------------------------
## PLOT BI EFFECTS ON SI
## EXTRACT SIMULATIONS OF COEFFICIENTS AND PREDICTOR COVARIATE

alphas0 <- preds4[,which(dimnames(preds4)[[2]]=="alpha0")]
alphas1 <- preds4[,which(dimnames(preds4)[[2]]=="alpha1")]
alphas2 <- preds4[,which(dimnames(preds4)[[2]]=="alpha2")]
varPro <- preds4[,which(dimnames(preds4)[[2]]=="varunshared")]
BIfs <- preds4[,which(dimnames(preds4)[[2]]=="BIf[1]"):which(dimnames(preds4)[[2]]=="BIf[44]")]
mus <- preds4[,which(dimnames(preds4)[[2]]=="mu[1]"):which(dimnames(preds4)[[2]]=="mu[44]")]
BIs <- BIfs[,order(colMeans(BIfs))]
SIs <- mus[,order(colMeans(BIfs))]
SI <- AFSRes$Afs[6:49][order(colMeans(BIfs))]
si <- matrix(rep(predict(lm(SI ~ colMeans(BIs))),15000),ncol=44,byrow=T)
BI <- log(NTSims[,1:44])-6
##
# exp(colMeans(BIs)+6)
# BIs=log(mBIf)-6
##
plot(colMeans(BIs),exp(mean(alphas0)+(1+mean(alphas1))*mean(log(colMeans(SIs)))+mean(alphas2)*colMeans(BIs)),
     type="n",xlim=c(0.03,1.08),ylim=c(1400,26400),xlab="",ylab="",axes=F)
polygon(c(colMeans(BIs),rev(colMeans(BIs))),
        c(apply(exp(varPro+alphas0+(1+alphas1)*log(si)+alphas2*BIs),2,quantile,0.95),
          rev(apply(exp(varPro+alphas0+(1+alphas1)*log(si)+alphas2*BIs),2,quantile,0.05))),
        col=addTrans("lightblue",150),border=F)
lines(colMeans(BIs),colMeans(exp(varPro+alphas0+(1+alphas1)*log(si)+alphas2*BIs)),col="blue")
for(i in 1:dim(BI)[2]){
  lines(apply(BI,2,quantile,c(0.05,0.95))[,i],rep(AFSRes$ANB2M[6:48][i],2),col="wheat4")
  lines(rep(colMeans(BI)[i],2),c(AFSRes$ANB2L[6:48][i],AFSRes$ANB2H[6:48][i]),col=addTrans("blue",100))
}
points(colMeans(BI),AFSRes$ANB2M[6:48],pch=21,col="black",bg="wheat4",cex=.75)
box()
axis(1,at=log(pretty(exp(colMeans(BIs)+6)))-6,labels=paste(pretty(exp(colMeans(BIs)+6))),tck=-0.025,cex.axis=.90,
     lty=1,padj=-1.0,col.ticks="wheat4")
axis(2,at=pretty(c(1400,26400)),labels=T,tck=-0.025,las=2,cex.axis=1,col.axis="black",col.ticks="blue") #,hadj=0.4)
mtext(expression(paste("Bird Island mature females ",group("(",scriptstyle(italic(hat(N)[t-3])),")"),sep="")),side=1,line=2,cex=.95)
##--------------------------------------------------
##  1st PLOT - SSB AND SIGNY AFS SMOOTHED TRENDS
AFSRes <- readRDS("AFStrends.Rds")
plot(AFSRes$year,AFSRes$Afs,type="n",xlim=c(1977,2024),ylim=c(1200,40000),xlab="",ylab="",axes=F) #ylim=c(1200,26250)
polygon(c(AFSRes$year,rev(AFSRes$year)),c(AFSRes$gamFH,rev(AFSRes$gamFL)),col=addTrans("lightblue",150),border=F)
lines(AFSRes$year,AFSRes$gamFM,lwd=1,col="blue")
for(i in 1:dim(AFSRes)[1])
  lines(rep(AFSRes$year[i],2),rbind(AFSRes$ANBL,AFSRes$ANBH)[,i],col="blue")
points(AFSRes$year,AFSRes$ANBM,pch=21,col="blue",bg="lightblue",cex=.75)
box()
axis(1,at=seq(1980,2020,by=5),labels=T,tck=-0.025,cex.axis=1,lty=1,padj=-1.0)
axis(1,at=seq(1977,2024,by=1),labels=F,tck=-0.012,las=2,cex.axis=1)
axis(4,at=pretty(c(1200,26250))[1:6],labels=T,tck=-0.02,las=2,cex.axis=.75,col.axis="blue",hadj=0.5)  # col.ticks="blue", hadj=0.7
axis(4,at=seq(0,27000,by=1000),labels=F,tck=-0.01,las=2,cex.axis=.75,col.ticks="blue")
mtext("a", side = 3, adj = -0.22, padj = -1.75, cex = .7, font = 2)
##
par(new=T)
##
lpAFS <- readRDS("lpAFS.Rds")
gamAFS <- readRDS("gamAFS.Rds")
NTSims <- readRDS("NTSims.Rds")
##
plot(1979:2024,colMeans(NTSims),type="n",pch=21,bg=2,xlim=c(1977,2024),ylim=c(-100,1200),axes=F,ylab="",xlab="")
polygon(c(1979:2024,2024:1979)+0.3, c(apply(gamAFS,2,quantile,0.975),rev(apply(gamAFS,2,quantile,0.025))),col=addTrans("wheat4",100),border="wheat4")
lines(1979:2024+0.3,colMeans(gamAFS),col="wheat4", lwd = 1.1) 
for(i in 1:46)
  lines(rep((1979:2024)[i],2)+0.3,apply(NTSims,2,quantile,c(0.025,0.975))[,i],col="wheat4")
points(1979:2024+0.3,colMeans(NTSims),type ="p",pch=21,col="black",bg="wheat4",cex=.75)
box()
axis(2,at=seq(400,1200,by=200),labels=T,tck=-0.02,las=2,cex.axis=.75,col.axis="black",col.ticks="wheat4",hadj=0.4)  # col.ticks="blue", hadj=0.7
axis(2,at=seq(300,1200,by=100),labels=F,tck=-0.01,las=2,cex.axis=.75,col.ticks="wheat4")
mtext("Signy Island population counts",side=4,line=2.25,cex=.85,col="blue",adj=0.125)
mtext("Bird Island mature females",side=2,line=1.9,cex=.85,col="black",adj=0.75)
##
par(new=F)

## CROSS CORRELATION AND NON-LINEAR CORRELATION
X <- as.vector(colMeans(NTSims[,1:43]))
Y <- as.vector(AFSRes$Afs[6:48])
cor.test(X,Y)
# [1] 0.5136513
NLCORXY <- nlcor(X,Y,plt=T)
NLCORXY$cor.estimate
# [1] 0.5136513
NLCORXY$adjusted.p.value
# [1] 0
print(NLCORXY$cor.plot)
  
par(new=F)
plot(-10:10,rep(0,21),type="n",axes=F,ylim=c(-0.5,0.5),xlim=c(-10,10),
     xlab = "",ylab="")
axis(1, at = seq(-10,10,by=5),labels=T,tck=-0.025,cex.axis=0.95,padj=-1.0)
axis(1, at = seq(-10,10,by=1),labels=F,tck=-0.012,cex.axis=0.95)
axis(2,at=seq(-0.8,0.5,by=0.1),labels=T,tck=-0.02,las=2,cex.axis=.75,hadj=0.6)
mtext(expression(paste("Cross-correlation, ",italic(rho[list(scriptscriptstyle(hat(N)[BI]),scriptscriptstyle(hat(N)[SI]))]),(tau),sep="")),side=2,line=1.5,cex=.95)
title(xlab="lag (years)",line=1.5,cex.lab=1)
box()
mira <- ccf(as.vector(colMeans(NTSims)),as.vector(detrend(AFSRes$Afs[3:48])),lag.max=10,plot=F)

l.mi <- qnorm(0.05)/sqrt(mira$n.used)
h.mi <- qnorm(0.95)/sqrt(mira$n.used)
abline(h=l.mi,lty=2)
abline(h=h.mi,lty=2)
##
#l.mi <- qnorm(0.025)/sqrt(mira$n.used)
#h.mi <- qnorm(0.975)/sqrt(mira$n.used)
#abline(h=l.mi,lty=2)
#abline(h=h.mi,lty=2)
##
abline(h=0)
for(i in -10:10)
  lines(c(i,i), c(0, as.vector(mira$acf)[i+11]))
##
par(new=F)
##------------------------------------
## WAVELET ANALYSES
## FILTER AND TRANSFORM THE DATA - DETREND AND SCALE COUNTS
NTSims <- readRDS("NTSims.Rds")
X <- as.vector(colMeans(NTSims))
Y <- as.vector(AFSRes$Afs[3:48])
X <- X-glbpf(X,9)
X <- cbind(1979:2024,scale(X)) # scale(fs1))
Y <- Y-glbpf(Y,9)
Y <- cbind(1979:2024,scale(Y)) # scale(fs2))
## CALCULATE COHERENCE WAVELET
wtc.AFS <- biwavelet::wtc(X,Y,sig.level=0.85,nrands=4999,max.scale=9)
## saveRDS(wtc.AFS,"wtc.AFS.Rds")
wtc.AFS <- readRDS("wtc.AFS.Rds")

par(xpd = F)
##
plot(wtc.AFS,plot.cb=T,legend.horiz=T,plot.phase=T,ylim=c(2,9),yaxt="n",xaxt="n",xlab="",ylab="", 
     cex.axis=0.75,lwd.sig=1.3,alpha.coi=0.7,xlim=c(1979,2023), arrow.col=addTrans("black",100))
axis(1,seq(1980,2020,by=5),labels=T,tck=-0.025,col=1,cex.axis=0.95,padj=-1.0)
axis(1,seq(1979,2024,by=1),labels=F,tck=-0.0125,col=1,cex.axis=0.95) #,padj=-1.5)
axis(2,log2(2:9),labels=paste(2:9),las=2,tck=-0.025,cex.axis=0.9,col=1,hadj=-0.0)
mtext("Scale (years)",side=2,line=2,cex=.85)
mtext(expression(paste("Coherence, ",scriptstyle(Rho[italic(list(x,y))]),scriptstyle(group("(",list(f,tau),")")),sep="")),side=3,line=0.2,cex=.9)
mtext("d",side=3,adj=-0.19,padj=-1.75,cex=.7,font=2)
box()

par(xpd = F)
xwt.AFS <- biwavelet::xwt(X,Y,sig.level=0.90,max.scale=9)
plot(xwt.AFS,plot.cb=T,legend.horiz=T,plot.phase=T,ylim=c(2,9),yaxt="n",xaxt="n",xlab="",ylab="", 
     cex.axis=0.75,lwd.sig=1.3,alpha.coi=0.7,xlim=c(1979,2023),
     arrow.cutoff=0.99,arrow.col=addTrans("black",100)) # ="white")
axis(1,seq(1980,2020,by=5),labels=T,tck=-0.025,col=1,cex.axis=0.95,padj=-1.0)
axis(1,seq(1979,2024,by=1),labels=F,tck=-0.0125,col=1,cex.axis=0.95) #,padj=-1.5)
axis(2,log2(2:9),labels=paste(2:9),las=2,tck=-0.025,cex.axis=0.9,col=1,hadj=-0.0)
mtext("Scale (years)",side=2,line=2,cex=.85)
mtext(expression(paste("Cross-wavelet, ",scriptstyle(W[italic(list(x,y))]),scriptstyle(group("(",list(f,tau),")")),sep="")),side=3,line=0.2,cex=.9)
mtext("c",side=3,adj=-0.22,padj=-1.75,cex=.7,font=2)
box()


save.image("C:/Users/jfor/OneDrive - NERC/Workstuff/AntarcticFurSeals/IUCN_assessment/Analysis/AFSDynamicsSigny.RData")
##
##-----------------------------------


