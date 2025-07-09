################################################################################
## FEMALE FUR SEAL ABUNDANCE AT BIRD ISLAND'S SPECIAL STUDY BEACH 2001 TO 2025 
##------------------------------------------------------------------------------
## LOAD REQUIRED PACKAGES
library(jagsUI)
library(mgcv)


##  LOAD WORKSPACE IF REQUIRED
# load("processed_data/SSBFemaleIPM19842025.RData")
##------------------------------------------------------------------------------------------##
### SUPPORT FUNCTIONS 
source("addTrans.R") # plotting colour transparency
source("glbpf.R")    # Gaussian band pass filter, from Cazelles et al. 2008
source("getfiltered.R")
source("estBetaParams.R")
## additional function to detrend a time series
detrend <- function(x) resid(lm(x ~ seq(x)))
##------------------------------------------------------------------------------------------## 
##  READ, FORMAT AND SUBSET DATA
### 1.  TOTAL PUP PRODUCTION COUNTS
##  1.1 PRODUCTIVITY: SSB total pup production from 1983/1984 to 2023/2024 
##         84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99  00  01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
Prod <- c(850,678,854,820,777,757,822,545,718,770,714,591,801,733,536,693,446,696,747,713,769,652,464,768,736,782,271,568,547,584,258,525,403,361,496,282,409,368,247,335,303,368)
prod <- Prod
##  1.2. EXTENDED SERIES - 1979/1980 to 2024-25
Prod2 <- c(672,668,669,661,748,Prod)
##----
##  1.3  Import capture-recapture data
##
## EVENT DEFINITION IN RAW DATA FILES - IMPORTANT AS USED IN FILTERING
#  1	pup	
# 11	pup serial weighing	
# 12	mother serial weighing or capture	
#  2	mother with pup that survives	
#  3	mother with pup that is recovered	
#  4	mother with pup of unknown fate	
#  5	female without pup	
#  6	DEAD	
#  7	birth elsewhere	
#  8	seen elsewhere	
####
SSBfemales <- read.csv("raw_data/SSB_females_1984-2025.csv", header = T)
SSBfemales2 <- read.csv("raw_data/SSB_pups_1984-2025.csv", header = T)
###  DEFINE BREEDING STATES: 1 breeder 2 non-breeder
SSBfemales$event <- ifelse(SSBfemales$event ==7,4,SSBfemales2$event)
SSBfemales <- SSBfemales[-which(SSBfemales$event>5),]    ## REMOVE SERIAL WEIGHINGS, DEATS AND OBSERVATIONS OUTSIDE SSB
SSBfemales$event <- ifelse(SSBfemales$event < 5,1,2)
SSBfemales2$event <- ifelse(SSBfemales2$event == 1,0,SSBfemales2$event)
SSBfemales2$event <- ifelse(SSBfemales2$event == 6,0,SSBfemales2$event)
SSBfemales2$event <- ifelse(SSBfemales2$event == 7,1,SSBfemales2$event)
SSBfemales2$event <- ifelse(SSBfemales2$event > 7,0,SSBfemales2$event)
SSBfemales2 <- SSBfemales2[-which(SSBfemales2$event ==0),]
SSBfemales2$event <- ifelse(SSBfemales2$event < 5,1,2)
###  create capture histories matrix
CH <- matrix(nrow = 0, ncol = length(1984:2025))
#### loop for each lifted female
for(i in sort(unique(SSBfemales$num))){
  id <- which(SSBfemales$num == i)
  statevec <- SSBfemales$event[id]
  # include capture history if the seal has ever bred at SSB
  if(any(statevec == 1)){
    ch <- rep(0, length(1984:2025))
    ch[match(SSBfemales$season[id], 1984:2025)] <- statevec
    CH <- rbind(CH, ch)
  }
}
#
CH2 <- matrix(nrow = 0, ncol = length(1984:2025))
#### loop for each pup
for(i in sort(unique(SSBfemales2$num))){
  id <- which(SSBfemales2$num == i)
  statevec <- SSBfemales2$event[id]
  # include capture history if the seal has ever bred at SSB
  if(any(statevec==1)){
    ch <- rep(0,length(1984:2025))
    ch[match(SSBfemales2$season[id],1984:2025)] <- statevec
    ch[1:(min(which(ch==1))-1)]<-0
    CH2 <- rbind(CH2,ch)
  }
}
# merge capture histories
CH <- rbind(CH, CH2)
dimnames(CH)[[1]] <- NULL
CH <- CH[-which(rowSums(CH)==0),]
rm(ch, CH2, i, id, statevec)
#
#
## Get vector with occasion of first capture for each female
get.first <- function(x) min(which(x != 0))
f <- apply(CH, 1, get.first)
### Recode CH matrix to 1 = seen as Breeder , 2 = seen as non-breeder, 3 = not seen
rCH <- CH
rCH[rCH == 0] <- 3
#   Functions borrowed from Kery and Schaub to:
##  1. create known latent states z
known.state.ms <- function(ms, notseen){
  state <- ms
  state[state == notseen] <- NA
  for (i in 1:dim(ms)[1]){
    m <- min(which(!is.na(state[i, ])))
    state[i, m] <- NA
  }
  return(state)
}
## generate initial values for unknown z
ms.init.z <- function(ch, f){
  for(i in 1:dim(ch)[1]) 
    ch[i, 1:f[i]] <- NA
  states <- max(ch, na.rm = TRUE)
  known.states <- 1:(states - 1)
  v <- which(ch == states)
  ch[-v] <- NA
  ch[v] <- sample(known.states, length(v), replace = TRUE)
  return(ch)
}
#####------------------------------------------------------------------------------
##  ADULT BREEDING MULTISTATE MODEL STRUCTURE
##  Multi-state Model Description
### States:
##  B: breeder      1
##  N: non-breeder  2
##  D: dead       (non-observable)   3
##
### Parameters: 
##  phi -> survival (for B and N)
##  psiBN -> transition B -> N
##  psiNB -> transition N -> B
##  p -> recapture (B and N)
##
### Transitions Matrix
#               breeder           non-breeder         dead
# breeder       phiB*(1-psiBN)    phiB*psiBN          1-phiB 
# non-breeder   phiN*psiNB        phiN*(1-psiNB)      1-phiN
# dead          0                 0                   1
#
### Recapture Matrix
#               seen.as.B seen.as.N not-seen
# breeder       pB        0         1 - pB 
# non-breeder   0         pN        1 - pN
# dead          0         0         1
##
#####------------------------------------------------------------------------------
###   MULTISTATE model with correlated random time effects on phi B and phi N, constant pB and 
##    random time effects on pN, random time effects on psiBN and psiNB
###
#### specify BUGS model in jags
sink("SSB.ms") 
cat("
  model{
    # Priors
    ## Initial population sizes
    NN[1] ~ dnorm(310, 0.0001)T(0,)
    NB[1] ~ dnorm(850, 0.0001)T(0,)
    NI[1] ~ dnorm(0, 0.0001)T(0,)
    #
    ## phi, psi, p and iota
    for(t in 1:(n.occasions-1)){
      ## survival and psiBN
      logit(phiB[t]) <- eta.phipsiBN[t, 1]
      logit(phiN[t]) <- eta.phipsiBN[t, 2]
      logit(psiBN[t]) <- eta.phipsiBN[t, 3]
      #
      ## psi NB independent temporal random effects
      logit(psiNB[t]) <- muNB + epsilonNB[t]    # temporal random efects on breeding or deferrers
      epsilonNB[t] ~ dnorm(0, tauNB)T(-15,15)
      # p parameters
      pB[t] <- mean.pB
      logit(pN[t]) <- mupN + epsilonpN[t]       # temporal random efects on recapture of non-breeders
      epsilonpN[t] ~ dnorm(0, taupN)T(-15,15)
      # iota
      log(iota[t]) <- miota + epsilonIo[t]
      epsilonIo[t] ~ dnorm(0, tau.iota)T(-15,15)
    }
    ## survival and psiBN
    for(t in 1:(n.occasions-1)){
      eta.phipsiBN[t, 1:3] ~ dmnorm(mu.phipsiBN[], Omega[,])
    }
    for (u in 1:3){
      mean.phipsiBN[u] ~ dunif(0, 1)    # Priors for mean state-specific survival
      mu.phipsiBN[u] <- log(mean.phipsiBN[u]/(1-mean.phipsiBN[u]))
    }
    Omega[1:3,1:3] ~ dwish(R[, ], 4)   # Priors for variance-covariance matrix
    Sigma.phipsiBN[1:3, 1:3] <- inverse(Omega[,])
    #
    ## state transitions
    muNB <- log(mean.psiNB / (1 - mean.psiNB))
    mean.psiNB ~ dunif(0, 1)                    # Prior for mean return to breeding
    tauNB <- pow(sigmaNB, -2)
    sigmaNB ~ dunif(0, 10)
    sigma2NB <- pow(sigmaNB, 2)
    ##  recapture
    mean.pB ~ dunif(0, 1)      # Priors for mean state-spec. recapture of breeders
    mupN <- log(mean.pN / (1 - mean.pN))
    mean.pN ~ dunif(0, 1)                     # Prior for mean survival of N
    taupN <- pow(sigmapN, -2)
    sigmapN ~ dunif(0, 10)
    sigma2pN <- pow(sigmapN, 2)
    ## mean immigration
    miota ~ dnorm(0, 0.0001)T(-10,10)
    sigma.iota ~ dunif(0, 10)
    tau.iota <- pow(sigma.iota, -2)
    #
    ## derived parameters
    eiota <- exp(miota)             # mean immigration rate
    for(t in 1:n.occasions){
      Ntot[t] <- NBI[t]+NN[t]       # total population size
    }
    for(t in 1:(n.occasions-1)){
      lambda[t] <- Ntot[t+1]/Ntot[t]  # population growth rate
      loglambda[t] <- log(lambda[t])  # log-lambda
    }
    meanlambda <- exp((1/(n.occasions-1))*sum(loglambda[1:(n.occasions-1)]))  # geometric mean
    #
    ## LIKELIHOODS
    ## Productivity
    ### system process
    for(t in 2:n.occasions){
      meanN[t] <- NBI[t-1]*phiB[t-1]*psiBN[t-1]+NN[t-1]*phiN[t-1]*(1-psiNB[t-1])
      NN[t] ~ dpois(meanN[t])
      meanB[t] <- NBI[t-1]*phiB[t-1]*(1-psiBN[t-1])+NN[t-1]*phiN[t-1]*psiNB[t-1]
      NB[t] ~ dpois(meanB[t])
      meanI[t] <- (NBI[t-1]+NN[t-1])*iota[t-1]
      NI[t] ~ dpois(meanI[t])
    }
    ### observation process
    for(t in 1:n.occasions){
      NBI[t] <- NB[t] + NI[t]
      prod[t] ~ dpois(NBI[t])
    }
    ## Capture-recapture
    ## State-transition and observation matrices
    for(i in 1:nind){  
      # Define probabilities of state S(t+1) given S(t)
      for (t in f[i]:(n.occasions-1)){
        ps[1,i,t,1] <- phiB[t] * (1 - psiBN[t])
        ps[1,i,t,2] <- phiB[t] * psiBN[t]
        ps[1,i,t,3] <- 1 - phiB[t]
        ps[2,i,t,1] <- phiN[t] * psiNB[t]
        ps[2,i,t,2] <- phiN[t] * (1 - psiNB[t])
        ps[2,i,t,3] <- 1 - phiN[t]
        ps[3,i,t,1] <- 0
        ps[3,i,t,2] <- 0
        ps[3,i,t,3] <- 1
        #    
        ## Define probabilities of O(t) given S(t)
        po[1,i,t,1] <- pB[t]
        po[1,i,t,2] <- 0
        po[1,i,t,3] <- 1 - pB[t]
        po[2,i,t,1] <- 0
        po[2,i,t,2] <- pN[t]
        po[2,i,t,3] <- 1 - pN[t]
        po[3,i,t,1] <- 0
        po[3,i,t,2] <- 0
        po[3,i,t,3] <- 1
      } #t
    } #i
    # Likelihood 
    for(i in 1:nind){
      # Define latent state at first capture
      z[i, f[i]] <- y[i, f[i]]
      for(t in (f[i] + 1):n.occasions){
        # State process: draw S(t) given S(t-1)
        z[i,t] ~ dcat(ps[z[i, t - 1], i, t - 1, ])
        # Observation process: draw O(t) given S(t)
        y[i,  t] ~ dcat(po[z[i, t], i, t  - 1,  ])
      } #t
    } #i
  }
  ", fill = TRUE)
sink()
#
##
### Generate list with the data for jags
jags.data <- list(y = rCH, prod = prod, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 3), R = matrix(c(5, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3))
#
inits <- function()
  list(Omega = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3), 
       mean.psiNB = runif(1, 0, 1), 
       mean.pB = runif(1, 0, 1),
       mean.pN = runif(1, 0, 1),
       sigmapN = runif(1, 0, 1),
       z = ms.init.z(rCH, f),
       miota = rnorm(1, 0.1, 0.5),
       sigma.iota = runif(1, 0.1, 10),
       NN = round(runif(dim(rCH)[2], 10, 100), 0),
       NB = round(runif(dim(rCH)[2], 200, 700), 0),
       NI = round(runif(dim(rCH)[2], 10, 100), 0))
#
## Parameters monitored
parameters <- c("phiB","phiN","psiBN","psiNB","iota","eiota","Sigma.phipsiBN","mean.psiNB",
  "mean.pB","pN","mean.pN","sigma2NB","sigma2pN","sigma.iota","Ntot","NB","NN","NI","NBI",
  "lambda","loglambda","meanlambda")
#
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
#
### Run the model
sta <- Sys.time()
out.test <- jags(jags.data,inits,parameters,"SSB.ms",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
#
# Time difference of 15.16234 hours
#
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
## WRITE SIMULATIONS
write.csv(preds, "processed_data/AFSFemaleSimulations19842025.csv")
saveRDS(preds,file = "processed_data/AFSFemaleSimulations19842025.Rds")
save.image("processed_data/SSBFemaleIPM19842025.RData")
#
#########
## BACKPROJECT TIME SERIES BASED ON EARLIER (1979-1983) SSB PUP PRODUCTION ESTIMATES
##  -> Without early mark-recapture data (prior to 1984) vital rates are generated
##  with beta priors of the rates obtained with the previous model.
sink("SSBearly.ms") 
cat("
  model{
    # Priors
    ## Initial population sizes
    NN[1] ~ dnorm(263,0.0001)T(0,)
    NB[1] ~ dnorm(668,0.0001)T(0,)
    #NBI[1] ~ dnorm(668,0.0001)T(0,)
    NI[1] ~ dnorm(0,0.0001)T(0,)
    #
    ## LIKELIHOOD
    ### system process
    for(t in 2:6){
      meanN[t] <- NBI[t-1]*phiB[t-1]*psiBN[t-1]+NN[t-1]*phiN[t-1]*(1-psiNB[t-1])
      NN[t] ~ dpois(meanN[t])
      meanB[t] <- NBI[t-1]*phiB[t-1]*(1-psiBN[t-1])+NN[t-1]*phiN[t-1]*psiNB[t-1]
      NB[t] ~ dpois(meanB[t])
      meanI[t] <- (NBI[t-1]+NN[t-1])*iota[t-1]
      NI[t] ~ dpois(meanI[t])
    }
    ### observation process
    for(t in 1:6){
      NBI[t] <- NB[t]+NI[t]
      Prod[t] ~ dpois(NBI[t])
      # NBI[t] <- Prod[t]
    }
    ##
    ## PRIORS
    phiB[1] ~ dbeta(88.64,3.91)
    phiB[2] ~ dbeta(86.24,6.23)
    phiB[3] ~ dbeta(82.06,23.74)
    phiB[4] ~ dbeta(85.94,23.35)
    phiB[5] ~ dbeta(124.27,25.94)
    phiN[1] ~ dbeta(51.89,3.35)
    phiN[2] ~ dbeta(46.23,4.64)
    phiN[3] ~ dbeta(47.14,6.44)
    phiN[4] ~ dbeta(39.36,8.16)
    phiN[5] ~ dbeta(50.29,6.69)
    psiBN[1] ~ dbeta(18.08,75.08)
    psiBN[2] ~ dbeta(26.69,85.32)
    psiBN[3] ~ dbeta(29.14,81.10)
    psiBN[4] ~ dbeta(25.17,90.58)    
    psiBN[5] ~ dbeta(33.41,123.18)
    psiNB[1] ~ dbeta(26.77,16.41)
    psiNB[2] ~ dbeta(26.32,14.08)
    psiNB[3] ~ dbeta(30.05,19.68)
    psiNB[4] ~ dbeta(30.06,21.98)    
    psiNB[5] ~ dbeta(34.45,22.16)
    iota[1] ~ dbeta(4.51,26.22)
    iota[2] ~ dbeta(4.59,44.84)
    iota[3] ~ dbeta(7.10,180.5883)
    iota[4] ~ dbeta(6.66,38.51)
    iota[5] ~ dbeta(9.70,38.21)
    #
    #for(t in 1:5){
      # phiB[t] ~ dbeta(228.5648,42.81928)
      # phiN[t] ~ dbeta(152.9776,21.5971)
      # psiBN[t] ~ dbeta(77.54869,239.3087)
      #psiNB[t] ~ dbeta(101.4443,68.66707)
      #iota[t] ~ dbeta(26.9753,180.5883)
    #}
    for(t in 1:6){
      Ntot[t] <- NBI[t]+NN[t]       # total population size
    }
    for(t in 1:5){
      lambda[t] <- Ntot[t+1]/Ntot[t]  # population growth rate
      loglambda[t] <- log(lambda[t])  # log-lambda
    }
    #meanlambda <- exp((1/5)*sum(loglambda[]))  # geometric mean
  }
  ", fill = TRUE)
sink()
#
jags.data2 <- list(Prod=Prod2[1:6])
#
inits <- function()
  list(NN=round(runif(6,10,100),0),NB=round(runif(6,200,700),0),NI=round(runif(6,10,100),0))
#
## Parameters monitored
parameters <- c("phiB","phiN","psiBN","psiNB","iota","Ntot","NB","NN","NI","NBI", "lambda", "loglambda", "meanlambda")
#
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
#
### Run the model
sta <- Sys.time()
out.early <- jags(jags.data2,inits,parameters,"SSBearly.ms",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
#    
## SAVE RESULTS
preds2 <- rbind(out.early$samples[[1]],out.early$samples[[2]],out.early$samples[[3]])
## write simulations
write.csv(preds2, "processed_data/AFSFemaleEarlySimulations19842025.csv")
saveRDS(preds2,file="processed_data/AFSFemaleEarlySimulations19842025.Rds")
save.image("processed_data/SSBFemaleIPM19842025.RData")
#
################################################################################
## GAUSSIAN LOW PASS FILTER FOR ABUNDANCE REPRESENTATION
## EXTRACT SIMULATIONS OF NT
NTSims <- preds[, which(dimnames(preds)[[2]] == "Ntot[1]"):which(dimnames(preds)[[2]] == "Ntot[42]")]
## RUN THIS BIT ONCE AND SAVE RESULTS
years <- 1984:2025
lpAFS8425 <- gamAFS8425 <- matrix(NA,dim(NTSims)[1],dim(NTSims)[2])
for(i in 1:dim(NTSims)[1]){
  lpAFS8425[i,] <- glbpf(NTSims[i,],pf=9)
  gamAFS8425[i,] <- predict(gam(NTSims[i,] ~ s(years)))
  print(i)  
}
saveRDS(lpAFS8425,file="lpAFS8425.Rds")
saveRDS(gamAFS8425,file="gamAFS8425.Rds")
saveRDS(NTSims,file="NTSims8425.Rds")
##
NTSims2 <- cbind(preds2[,which(dimnames(preds2)[[2]]=="Ntot[1]"):which(dimnames(preds2)[[2]]=="Ntot[5]")],NTSims)
## RERUN THIS BIT ONCE AND SAVE RESULTS
lpAFS2 <- matrix(NA, dim(NTSims2)[1], dim(NTSims2)[2])
for(i in 1:dim(NTSims2)[1]){
  lpAFS2[i, ] <- glbpf(NTSims2[i,],pf=8.2)
  print(i)  
}
saveRDS(lpAFS2, file = "lpAFS2.Rds")
lpAFS2 <- readRDS("lpAFS2.Rds")
###
## SAME WITH GAM
gyears <- 1984:2025
dat <- data.frame(y=round(colMeans(NTSims)),x=gyears)
b0 <- gam(y~s(x),family=poisson(link="log"),data=dat,method="REML")
##
## a plot
plot(1984:2025,out.test$mean$Ntot, type = "n", pch = 21, bg = 2, ylim = c(200, 1400), ylab = "mature females", xlab = "")
polygon(c(1984:2025, 2025:1984),  c(apply(lpAFS, 2, quantile, 0.975), rev(apply(lpAFS, 2, quantile, 0.025))), col =  addTrans("antiquewhite4", 100), border = "antiquewhite4")
lines(1984:2025, colMeans(lpAFS), col = "antiquewhite4", lwd = 1.1) 
for(i in 1:41)
  lines(rep((1984:2025)[i], 2), c(out.test$q2.5$Ntot[i], out.test$q97.5$Ntot[i]), col = "black")
points(1984:2025, out.test$mean$Ntot, type = "b", pch = 21, col = "black", bg = "magenta")
##
NTSims <- cbind(
  preds2[,which(dimnames(preds2)[[2]]=="Ntot[1]"):which(dimnames(preds2)[[2]]=="Ntot[6]")],
  preds[, which(dimnames(preds)[[2]] == "Ntot[2]"):which(dimnames(preds)[[2]] == "Ntot[42]")]
)
saveRDS(NTSims,file="NTSims7925.Rds")
years <- 1979:2025
lpAFS <- gamAFS <- matrix(NA,dim(NTSims)[1],dim(NTSims)[2])
for(i in 1:dim(NTSims)[1]){
  lpAFS[i,] <- glbpf(NTSims[i,],pf=9)
  gamAFS[i,] <- predict(gam(NTSims[i,] ~ s(years)))
  print(i)  
}
saveRDS(lpAFS,file="lpAFS7925.Rds")
saveRDS(gamAFS,file="gamAFS7925.Rds")
lpAFS <- readRDS("lpAFS7925.Rds")
gamAFS <- readRDS("gamAFS7925.Rds")
NTSims <- readRDS("NTSims7925.Rds")
#
LogLambdas <- cbind(
  preds2[,which(dimnames(preds2)[[2]]=="loglambda[1]"):which(dimnames(preds2)[[2]]=="loglambda[5]")],
  preds[, which(dimnames(preds)[[2]] == "loglambda[1]"):which(dimnames(preds)[[2]] == "loglambda[41]")]
)
dimnames(LogLambdas)[[2]] <- paste("LL",1:46,sep="")
## EXTRACT VECTOR OF LOG-LAMBDAS FOR 1995:2006 (12 years)
## exp(sum(colMeans(LogLambdas)[16:27])/length(colMeans(LogLambdas)[16:27]))
## [1] 0.988878
saveRDS(exp(rowSums(LogLambdas[,16:27])/12),"processed_data/earlyLogLambdas.Rds")
##
## LogLambda since 2000
LL0025 <- exp(rowSums(LogLambdas[,21:46])/length(21:46))
round(c(mean=mean(LL0025),SD=sd(LL0025),quantile(LL0025,c(0.025,0.975))),3)
#  mean    SD  2.5% 97.5% 
# 0.973 0.003 0.968 0.978