################################################################################
## FEMALE FUR SEAL ABUNDANCE AT BIRD ISLAND'S SPECIAL STUDY BEACH 2001 TO 2025 
##------------------------------------------------------------------------------
## LOAD REQUIRED PACKAGES
library(jagsUI)
library(mgcv)
#------------------------
##  READ, FORMAT AND SUBSET SSB DATA
### 1.  TOTAL PUP PRODUCTION COUNTS
##  1.1 PRODUCTIVITY: SSB total pup production from 1983/1984 to 2023/2025
##         84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99  00  01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
Prod <- c(850,678,854,820,777,757,822,545,718,770,714,591,801,733,536,693,446,696,747,713,769,652,464,768,736,782,271,568,547,584,258,525,403,361,496,282,409,368,247,335,303,368)
## PRODUCTIVITY FOR 2001:2025
prod0125 <- Prod[18:42]
# NUMBERS OF KNOWN PUPS SEXED AT BIRTH
##            01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
allPups <- c(167,141,124,142,316,271,428,293,337,221,284,278,330,205,326,328,243,321,220,280,300,247,340, 34,368)
femPups <- c( 92, 79, 51, 71,147,135,217,149,164,111,136,139,151,119,164,165,123,168,108,140,160,141,179, 16,171) # Females
malPups <- c( 75, 61, 73, 71,169,136,211,144,173,110,149,139,179, 86,162,163,120,153,112,140,140,106,159, 18,197) # Males
## PROPORTION OF FEMALES (MEAN AND SE)
pFem <- femPups / allPups
pFemSE <- sqrt(pFem*(1-pFem)/allPups)
## EARLY SURVIVAL -from birth to tagging (approx. first 1.5 months)
surToTag <- 1 - c(0.265804598, 0.18875502, 0.130434783, 0.239271782, 0.197247706, 0.155172414, 0.266927083,
                  0.244565217, 0.369565217, 0.034090909, 0.163732394, 0.098540146, 0.233276158, 0.03875969,
                  0.167619048, 0.406947891, 0.089709763, 0.114919355, 0.06405694, 0.102941176, 0.24744898,
                  0.052631579, 0.053731343, 0.122112211, 0.04619565)
## AVERAGE -within year- CV IS 0.1215081
surToTagVAR <- (surToTag*0.1215081)^2
n0 <- prod0125 * pFem * surToTag
#--------------------------------------------------------------------------------------------
##  IMPORT CAPTURE-RECAPTURE DATA
##  FEMALES TAGGED AS PUPS SINCE 2001
SSBfemales <- read.csv("raw_data/SSB_pups_2001-2025.csv", header = T)
##  FILTER OUT PUP SERIAL WEIGHTS (event 11) AND DEATHS (event 6)
SSBfemales$event <- ifelse(SSBfemales$event == 7, 3, SSBfemales$event)
SSBfemales$event <- ifelse(SSBfemales$event >= 6, 0, SSBfemales$event)
SSBfemales <- SSBfemales[-which(SSBfemales$event == 0), ]
## CHANGE ALL BREEDING EVENTS TO 3 AND ALL NON-BREEDING EVENTS TO 2
SSBfemales$event <- ifelse(SSBfemales$event > 1 & SSBfemales$event < 5, 3, SSBfemales$event)
SSBfemales$event <- ifelse(SSBfemales$event == 5, 2, SSBfemales$event)
##
## Observations (O):
## 1 seen as weanling
## 2 seen as pre-breeder
## 3 seen breeder
## 4 not seen
## -------------------------------------------------
## FORMAT DATA FOR THE INTEGRATED MODEL
lastyear <- 2025
CH1 <- matrix(nrow = 0, ncol = length(2001:2025))
### loop for each pup
##  FOR DATA ENTRY ERROR CONTROL
pupNum <- NULL
for(i in sort(unique(SSBfemales$num))){
  id <- which(SSBfemales$num == i)
  pupNum <- c(pupNum, i) 
  obsyears <- SSBfemales$season[id][1]:lastyear
  obsvector <- rep(4, length(obsyears))
  obsvector[match(SSBfemales$season[id], obsyears)] <- SSBfemales$event[id]
  ## CHECK REPETITION OF EVENTS WITH VALUE 1
  if(length(which(obsvector == 1)) > 1){
    warning(cat("pup Num ", i, " has more than 1 event = 1"))
    break
  }
  ## RECRUITED NONBREEDERS MUST BE BREEDERS (i.e. 3; OR ADULT)
  if(any(obsvector == 3)){
    ids2 <- which(obsvector == 2)
    if(any(ids2 > min(which(obsvector == 3))))
      obsvector[ids2[ids2 > min(which(obsvector == 3))]] <- 3
  }
  ch1 <- rep(0, length(2001:lastyear))
  ch1[match(obsyears, 2001:lastyear)] <- obsvector
  CH1 <- rbind(CH1, ch1)
  print(i)
}
dimnames(CH1)[[1]] <- NULL
## cleanup
rm(ch1, obsyears, obsvector) #, lastyear)
##
## VECTOR OF OCCASION AT FIRST CAPTURE OF EACH PUP
get.first <- function(x) min(which(x != 0))
f1 <- apply(CH1, 1, get.first)
## -> SUBSEQUENT CHECKS
## ENSURE THAT ALL PUP HISTORIES START WITH 1
er1 <- NULL
for(i in 1:dim(CH1)[1]){
  if(CH1[i, f1[i]] != 1){
    er1 <- c(er1, i)
    CH1[i, f1[i]] <- 1
  }
} ## er1 should be NULL
## ENSURE THAT AT AGE 2 (1) OR 3 (2) THEY CANNOT BE BREEDERS
## ENSURE THAT AT AGE 9 ALL PREBREEDERS ARE BREEDERS
age <- array(NA, dim = dim(CH1))
for(i in 1:nrow(CH1)){
  for(t in f1[i]:ncol(CH1)){
    age[i, t] <- min(c(t - f1[i] + 1, 9))
  }
}
for(i in 1:nrow(CH1)){
  for(t in f1[i]:ncol(CH1)){
    if(CH1[i, t] == 3 & age[i, t] == 2) CH1[i, t] <- 2
    if(CH1[i, t] == 3 & age[i, t] == 3) CH1[i, t] <- 2
    if(CH1[i, t] == 2 & any(CH1[i, 1:(t - 1)] == 3)) CH1[i, t] <- 3
    if(CH1[i, t] == 2 & age[i, t] >= 9) CH1[i, t] <- 3
  }
}
rm(age)
## REPLACE 0s BY 4s
## 1 = seen as pup, 2 = seen as prebreeder, 3 = seen as breeders, 4 = not seen
rCH1 <- CH1
rCH1[rCH1 == 0] <- 4
## FURTHER CHECKS  (if required)
## ...
## FUNCTION TO GENERATE INITIAL VALUES FOR TRUE LATENT STATE MATRIX z
## (This is for analysis with JAGS, below)
agefirst.init <- function(ch, f){
  age <- array(NA, dim = dim(ch))
  for(i in 1:nrow(ch)){
    for(t in f[i]:ncol(ch)){
      age[i, t] <- min(c(t - f[i] + 1, 9))
    }
  }
  ini <- array(NA, dim=dim(ch))
  for (i in 1:nrow(ch)){
    for (t in f[i]:ncol(ch)){
      if(ch[i, t] == 1) ini[i, t] <- 1
      for(j in 2:8)
        if(ch[i, t] == 2 & age[i, t] == j) ini[i, t] <- j
      for(j in 4:8)
        if(ch[i, t] == 3 & age[i, t] == j) ini[i, t] <- 9
      if(ch[i, t] == 4 & age[i, t] == 9) ini[i, t] <- 9
    }
  }
  ini[which(is.na(ini))] <- age[which(is.na(ini))]
  for (i in 1:nrow(ch)){
    ini[i, f[i]] <- NA
    for (t in f[i]:(ncol(ch) - 1)){
      if(ini[i, t] == 9 & ini[i, t + 1] == 4) ini[i, t + 1] <- 9
      if(ini[i, t] == 9 & ini[i, t + 1] == 5) ini[i, t + 1] <- 9
      if(ini[i, t] == 9 & ini[i, t + 1] == 6) ini[i, t + 1] <- 9
      if(ini[i, t] == 9 & ini[i, t + 1] == 7) ini[i, t + 1] <- 9
      if(ini[i, t] == 9 & ini[i, t + 1] == 8) ini[i, t + 1] <- 9
    }
  }
  return(ini)
}
z1 = agefirst.init(rCH1, f1)
## ----------------------------------------------
##  IMPORT SECOND CAPTURE-RECAPTURE DATA SET
##  FEMALES TAGGED AS BREEDERS 2001 TO 2025
###   Lifted females and females tagged as pups since 2001
SSBfemales <- read.csv("raw_data/SSB_females_1984-2025.csv", header = T)
SSBfemales2 <- read.csv("raw_data/SSB_pups_1984-2025.csv", header = T)
###  define breeding states
SSBfemales$event <- ifelse(SSBfemales$event < 5, 1, 2)
SSBfemales2$event <- ifelse(SSBfemales2$event == 1, 0, SSBfemales2$event)
SSBfemales2$event <- ifelse(SSBfemales2$event >= 11, 0, SSBfemales2$event)
SSBfemales2$event <- ifelse(SSBfemales2$event == 6, 0, SSBfemales2$event)
SSBfemales2$event <- ifelse(SSBfemales2$event > 0 & SSBfemales2$event < 5, 1, 2)
##
## Observations (O):
## 1 seen as breeder
## 2 seen as non-breeder
## 3 not seen
## -------------------------------------------------
###  CREATE CAPTURE-HISTORIES MATRIX
###  create capture histories matrix
CH2 <- matrix(nrow = 0, ncol = length(1984:2025))
#
#### loop for each lifted female
for(i in sort(unique(SSBfemales$num))){
  id <- which(SSBfemales$num == i)
  statevec <- SSBfemales$event[id]
  # include capture history if the seal has ever bred at SSB
  if(any(statevec == 1)){
    ch2 <- rep(0, length(1984:2025))
    ch2[match(SSBfemales$season[id], 1984:2025)] <- statevec
    CH2 <- rbind(CH2, ch2)
  }
}
#
CH3 <- matrix(nrow = 0, ncol = length(1984:2025))
#### loop for each pup
for(i in sort(unique(SSBfemales2$num))){
  id <- which(SSBfemales2$num == i)
  statevec <- SSBfemales2$event[id]
  # include capture history if the seal has ever bred at SSB
  if(any(statevec == 1)){
    ch3 <- rep(0, length(1984:2025))
    ch3[match(SSBfemales2$season[id], 1984:2025)] <- statevec
    ch3[1:(min(which(ch3 == 1)) - 1)] <- 0
    CH3 <- rbind(CH3, ch3)
  }
}
# merge capture histories
CH2 <- rbind(CH2, CH3)
dimnames(CH2)[[1]] <- NULL
rm(ch2, ch3, CH3, i, id, statevec)
#
## reduce capture histories to period 2001:2025
CH2 <- CH2[, 18:42]
CH2 <- CH2[-which(rowSums(CH2) == 0), ]
#
## Get vector with occasion of first capture for each female
f2 <- apply(CH2, 1, get.first)
### Recode CH matrix to 1 = seen as Breeder , 2 = seen as non-breeder, 3 = not seen
rCH2 <- CH2
rCH2[rCH2 == 0] <- 3
#   Functions borrowed from Kery and Schaub to:
##  1. create known latent states z for data input
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
#
####----------------------------------------------------------------------------
#### JAGS MODELLING (MCMC)
#### specify BUGS model in JAGS
sink("SSB.IPM") 
cat("
  model{
    # Priors
    ## Initial population sizes
    N1[1] ~ dnorm(110,0.0001)T(0,)
    N2[1] ~ dnorm(87,0.0001)T(0,)
    NB3[1] ~ dnorm(2.23,0.0001)T(0,)
    NP3[1] ~ dnorm(71,0.0001)T(0,)
    NB4[1] ~ dnorm(15,0.0001)T(0,)
    NP4[1] ~ dnorm(44,0.0001)T(0,)
    NB5[1] ~ dnorm(12,0.0001)T(0,)
    NP5[1] ~ dnorm(22,0.0001)T(0,)
    NB6[1] ~ dnorm(3,0.0001)T(0,)
    NP6[1] ~ dnorm(13,0.0001)T(0,)
    NB7[1] ~ dnorm(1.3,0.0001)T(0,)
    NP7[1] ~ dnorm(6,0.0001)T(0,)
    NP[1] ~ dnorm(3.5,0.0001)T(0,)
    NN[1] ~ dnorm(254,0.0001)T(0,)
    NB[1] ~ dnorm(696,0.0001)T(0,)
    NI[1] ~ dnorm(0,0.0001)T(0,)
    #
    ## Sex ratio
    for(i in 1:T){
      sr[i] <- srp[i]/nFp[i]
      srp[i] ~ dbin(pF[i],nFp[i])
    }
    # perinatal survival
    #for(i in 1:T){
    #  phi0[i] ~ dbeta(alphaS[i], betaS[i])
    #}
    ## PRIORS AND CONSTRAINTS
    for(t in 1:(T - 1)){
      ## Adult survival and psiBN with correlated temporal random efffects
      logit(phiB[t]) <- eta.phipsiBN[t, 1]
      logit(phiN[t]) <- eta.phipsiBN[t, 2]
      logit(psiBN[t]) <- eta.phipsiBN[t, 3]
      ## FOR 1ST DATA SET ANALYSIS
      phi.B[t] <- mean.phiB
      ## PRE-BREEDER COMPONENT
      ## phi - epsilon. are temporal random efects
      logit(phi.1[t]) <- mu.1 + epsilon.1[t]
      epsilon.1[t] ~ dnorm(0, tau.1)T(-15, 15)
      logit(phi.2[t]) <- mu.2 + epsilon.2[t]
      epsilon.2[t] ~ dnorm(0, tau.2)T(-15, 15)
      logit(phi.3[t]) <- mu.3 + epsilon.3[t]
      epsilon.3[t] ~ dnorm(0, tau.3)T(-15, 15)
      logit(phi.4[t]) <- mu.4 + epsilon.4[t]
      epsilon.4[t] ~ dnorm(0, tau.4)T(-15, 15)
      logit(phi.5[t]) <- mu.5 + epsilon.5[t]
      epsilon.5[t] ~ dnorm(0, tau.5)T(-15, 15)
      logit(phi.6[t]) <- mu.6 + epsilon.6[t]
      epsilon.6[t] ~ dnorm(0, tau.6)T(-15, 15)
      logit(phi.7[t]) <- mu.7 + epsilon.7[t]
      epsilon.7[t] ~ dnorm(0, tau.7)T(-15, 15)
      logit(phi.8[t]) <- mu.8 + epsilon.8[t]
      epsilon.8[t] ~ dnorm(0, tau.8)T(-15, 15)
      #
      ## alpha (RECRUITMENT)
      logit(alpha.3[t]) <- mu.a3 + epsilon.a3[t]
      epsilon.a3[t] ~ dnorm(0, tau.a3)T(-15, 15)
      logit(alpha.4[t]) <- mu.a4 + epsilon.a4[t]
      epsilon.a4[t] ~ dnorm(0, tau.a4)T(-15, 15)
      logit(alpha.5[t]) <- mu.a5 + epsilon.a5[t]
      epsilon.a5[t] ~ dnorm(0, tau.a5)T(-15, 15)
      logit(alpha.6[t]) <- mu.a6 + epsilon.a6[t]
      epsilon.a6[t] ~ dnorm(0, tau.a6)T(-15, 15)
      logit(alpha.7[t]) <- mu.a7 + epsilon.a7[t]
      epsilon.a7[t] ~ dnorm(0, tau.a7)T(-15, 15)
      #
      ## BREEDING PROBABILITY FOR SKIPPERS
      ## psi NB independent temporal random effects
      logit(psiNB[t]) <- muNB + epsilonNB[t]
      epsilonNB[t] ~ dnorm(0, tauNB)T(-15, 15)
      #
      ## RECAPTURE (p)
      p.P1[t] <- mean.pP1
      p.P2[t] <- mean.pP2
      pB[t] <- mean.pB
      ## temporal random efects on recapture of skippers
      logit(pN[t]) <- mupN + epsilonpN[t]
      epsilonpN[t] ~ dnorm(0, taupN)T(-15, 15)
      ## IMMIGRATION (iota)
      log(iota[t]) <- miota + epsilonIo[t]
      epsilonIo[t] ~ dnorm(0, tau.iota)T(-15, 15)
    }
    ## survival and psiBN
    for(t in 1:(T-1)){
      eta.phipsiBN[t, 1:3] ~ dmnorm(mu.phipsiBN[], Omega[,])
    }
    for(u in 1:3){
      ## Priors for mean mature state-specific survival
      mean.phipsiBN[u] ~ dunif(0, 1)                    
      mu.phipsiBN[u] <- log(mean.phipsiBN[u] / (1-mean.phipsiBN[u]))
    }
    ## Priors for variance-covariance matrix
    Omega[1:3, 1:3] ~ dwish(R[, ], 4)
    Sigma.phipsiBN[1:3, 1:3] <- inverse(Omega[, ])
    #
    ## PREBREEDER SURVIVAL, priors, means and precisions
    mu.1 <- log(mean.phi.1 / (1 - mean.phi.1))
    mean.phi.1 ~ dunif(0, 1)
    tau.1 <- pow(sigma.1, -2)
    sigma.1 ~ dunif(0, 10)
    mu.2 <- log(mean.phi.2 / (1 - mean.phi.2))
    mean.phi.2 ~ dunif(0, 1)
    tau.2 <- pow(sigma.2, -2)
    sigma.2 ~ dunif(0, 10)
    mu.3 <- log(mean.phi.3 / (1 - mean.phi.3))
    mean.phi.3 ~ dunif(0, 1)
    tau.3 <- pow(sigma.3, -2)
    sigma.3 ~ dunif(0, 10)
    mu.4 <- log(mean.phi.4 / (1 - mean.phi.4))
    mean.phi.4 ~ dunif(0, 1)
    tau.4 <- pow(sigma.4, -2)
    sigma.4 ~ dunif(0, 10)
    mu.5 <- log(mean.phi.5 / (1 - mean.phi.5))
    mean.phi.5 ~ dunif(0, 1)
    tau.5 <- pow(sigma.5, -2)
    sigma.5 ~ dunif(0, 10)
    mu.6 <- log(mean.phi.6 / (1 - mean.phi.6))
    mean.phi.6 ~ dunif(0, 1)
    tau.6 <- pow(sigma.6, -2)
    sigma.6 ~ dunif(0, 10)
    mu.7 <- log(mean.phi.7 / (1 - mean.phi.7))
    mean.phi.7 ~ dunif(0, 1)                            # Prior for mean return to breeding
    tau.7 <- pow(sigma.7, -2)
    sigma.7 ~ dunif(0, 10)
    mu.8 <- log(mean.phi.8 / (1 - mean.phi.8))
    mean.phi.8 ~ dunif(0, 1)
    tau.8 <- pow(sigma.8, -2)
    sigma.8 ~ dunif(0, 10)
    mean.phiB ~ dunif(0, 1)     # Prior for mean breeder survival
    #
    ## RECRUITMENT, priors, means and precisions
    mu.a3 <- log(mean.alpha.3 / (1 - mean.alpha.3))
    mean.alpha.3 ~ dunif(0, 1)
    tau.a3 <- pow(sigma.a3, -2)
    sigma.a3 ~ dunif(0, 10)    
    mu.a4 <- log(mean.alpha.4 / (1 - mean.alpha.4))
    mean.alpha.4 ~ dunif(0, 1)
    tau.a4 <- pow(sigma.a4, -2)
    sigma.a4 ~ dunif(0, 10)  
    mu.a5 <- log(mean.alpha.5 / (1 - mean.alpha.5))
    mean.alpha.5 ~ dunif(0, 1)
    tau.a5 <- pow(sigma.a5, -2)
    sigma.a5 ~ dunif(0, 10)  
    mu.a6 <- log(mean.alpha.6 / (1 - mean.alpha.6))
    mean.alpha.6 ~ dunif(0, 1)
    tau.a6 <- pow(sigma.a6, -2)
    sigma.a6 ~ dunif(0, 10)  
    mu.a7 <- log(mean.alpha.7 / (1 - mean.alpha.7))
    mean.alpha.7 ~ dunif(0, 1)
    tau.a7 <- pow(sigma.a7, -2)
    sigma.a7 ~ dunif(0, 10)  
    #
    ## Breeding probability for skippers
    muNB <- log(mean.psiNB / (1 - mean.psiNB))
    ## Prior for mean return to breeding
    mean.psiNB ~ dunif(0, 1)
    tauNB <- pow(sigmaNB, -2)
    sigmaNB ~ dunif(0, 10)
    #
    ## RECAPTURE
    ## Prior for mean recapture young pre-breeders
    mean.pP1 ~ dunif(0, 1)
    ## Prior for mean recapture old pre-breeders
    mean.pP2 ~ dunif(0, 1)
    ## Priors for mean recapture of breeders
    mean.pB ~ dunif(0, 1)
    mupN <- log(mean.pN / (1 - mean.pN))
    ## Prior for mean survival of N
    mean.pN ~ dunif(0, 1)
    taupN <- pow(sigmapN, -2)
    sigmapN ~ dunif(0, 10)
    #
    ## IMMIGRATION RATE priors
    miota ~ dnorm(0, 0.0001)T(-10, 10)
    sigma.iota ~ dunif(0, 10)
    tau.iota <- pow(sigma.iota, -2)
    #
    ## DERIVED PARAMETERS
    ## mean immigration rate
    eiota <- exp(miota)
    for(t in 1:T){
      ## total mature seals
      Ntot[t] <- NBI[t]+NN[t]
      ## recruits
      Nrec[t] <- NB3[t]+NB4[t]+NB5[t]+NB6[t]+NB7[t]
      ## pre-recruits
      Npre[t] <- N1[t]+N2[t]+NP3[t]+NP4[t]+NP5[t]+NP6[t]+NP7[t]+NP[t]
    }
    for(t in 1:T){
      ## total female population size
      NAll[t] <- Ntot[t] + Npre[t]
      ## pregnancy rate
      PR[t] <- NBI[t] / Ntot[t]
    }
    for(t in 1:(T-1)){
      lambda[t] <- Ntot[t+1] / Ntot[t]  # population growth rate
      loglambda[t] <- log(lambda[t])  # log-lambda
    }
    meanlambda <- exp((1 / (T - 1)) * sum(loglambda[1:(T - 1)]))  # geometric mean
    #
    ## PRODUCTIVITY AND RECRUITMENT
    ### system process
    for(t in 2:T){
      mean1[t] <- sr[t-1]*phi0[t-1]*phi.1[t-1]*NBI[t-1]       # weaned pups surviving to age 1
      N1[t] ~ dpois(mean1[t])                                 # 1-year-olds (widges)
      mean2[t] <- phi.2[t-1]*N1[t-1]
      N2[t] ~ dpois(mean2[t])                                 # 2-year-olds (widges)
      meanb3[t] <- phi.3[t-1]*alpha.3[t-1]*N2[t-1]
      NB3[t] ~ dpois(meanb3[t])                               # 3-year-old recruits
      meanp3[t] <- phi.3[t-1]*(1-alpha.3[t-1])*N2[t-1]
      NP3[t] ~ dpois(meanp3[t])                               # 3-year-old pre-breeders
      meanb4[t] <- phi.4[t-1]*alpha.4[t-1]*NP3[t-1]
      NB4[t] ~ dpois(meanb4[t])                               # 4-year-old recruits
      meanp4[t] <- phi.4[t-1]*(1-alpha.4[t-1])*NP3[t-1]
      NP4[t] ~ dpois(meanp4[t])                               # 4-year-old pre-breeders
      meanb5[t] <- phi.5[t-1]*alpha.5[t-1]*NP4[t-1]
      NB5[t] ~ dpois(meanb5[t])                               # 5-year-old recruits
      meanp5[t] <- phi.5[t-1]*(1-alpha.5[t-1])*NP4[t-1]
      NP5[t] ~ dpois(meanp5[t])                               # 5-year-old pre-breeders
      meanb6[t] <- phi.6[t-1]*alpha.6[t-1]*NP5[t-1]
      NB6[t] ~ dpois(meanb6[t])                               # 6-year-old recruits
      meanp6[t] <- phi.6[t-1]*(1-alpha.6[t-1])*NP5[t-1]
      NP6[t] ~ dpois(meanp6[t])                               # 6-year-old pre-breeders
      meanb7[t] <- phi.7[t-1]*alpha.7[t-1]*NP6[t-1]
      NB7[t] ~ dpois(meanb7[t])                               # 7-year-old recruits
      meanp7[t] <- phi.7[t-1]*(1-alpha.7[t-1])*NP6[t-1]
      NP7[t] ~ dpois(meanp7[t])                               # 7-year-old pre-breeders
      meanNP[t] <- NP7[t-1]*phi.8[t-1]                        # Pre-breeders never observed breeding
      NP[t] ~ dpois(meanNP[t])
      meanB[t] <- NBI[t-1]*phiB[t-1]*(1-psiBN[t-1])+phi.3[t-1]*alpha.3[t-1]*N2[t-1]+phi.4[t-1]*alpha.4[t-1]*NP3[t-1]+phi.5[t-1]*alpha.5[t-1]*NP4[t-1]+phi.6[t-1]*alpha.6[t-1]*NP5[t-1]+phi.7[t-1]*alpha.7[t-1]*NP6[t-1]+NN[t-1]*phiN[t-1]*psiNB[t-1]
      NB[t] ~ dpois(meanB[t])                                 # Breeders
      meanN[t] <- NBI[t-1]*phiB[t-1]*psiBN[t-1]+NN[t-1]*phiN[t-1]*(1-psiNB[t-1])        # breeders that skip breeding
      NN[t] ~ dpois(meanN[t])                                 # Non-Breeders
      meanI[t] <- (NBI[t-1]+NN[t-1])*iota[t-1]                # IMMIGRANTS
      NI[t] ~ dpois(meanI[t])
    }
    #
    ## OBSERVATION PROCESS
    ## observed productivity relates directly to the number of observed breeders (NB + NI)
    for(t in 1:T){
      NBI[t] <- NB[t]+NI[t]
      Prod[t] ~ dpois(NBI[t])
    }
    ## FIRST CAPTURE-RECAPTURE; FEMALES TAGGED AS NEWBORN PUPS
    ## State-transition and observation matrices
     for(i in 1:nind1){
      ## Define probabilities of state S(t+1) given S(t)
      for(t in f1[i]:(T - 1)){
        ps1[1, i, t, 1] <- 0
        ps1[1, i, t, 2] <- phi.1[t]
        ps1[1, i, t, 3] <- 0
        ps1[1, i, t, 4] <- 0
        ps1[1, i, t, 5] <- 0
        ps1[1, i, t, 6] <- 0
        ps1[1, i, t, 7] <- 0
        ps1[1, i, t, 8] <- 0
        ps1[1, i, t, 9] <- 0        
        ps1[1, i, t,10] <- 1 - phi.1[t] 
        ps1[2, i, t, 1] <- 0
        ps1[2, i, t, 2] <- 0
        ps1[2, i, t, 3] <- phi.2[t]
        ps1[2, i, t, 4] <- 0
        ps1[2, i, t, 5] <- 0
        ps1[2, i, t, 6] <- 0
        ps1[2, i, t, 7] <- 0
        ps1[2, i, t, 8] <- 0
        ps1[2, i, t, 9] <- 0        
        ps1[2, i, t,10] <- 1 - phi.2[t] 
        ps1[3, i, t, 1] <- 0
        ps1[3, i, t, 2] <- 0
        ps1[3, i, t, 3] <- 0
        ps1[3, i, t, 4] <- phi.3[t] * (1 - alpha.3[t])
        ps1[3, i, t, 5] <- 0
        ps1[3, i, t, 6] <- 0
        ps1[3, i, t, 7] <- 0
        ps1[3, i, t, 8] <- 0
        ps1[3, i, t, 9] <- phi.3[t] * alpha.3[t]
        ps1[3, i, t,10] <- 1 - phi.3[t] 
        ps1[4, i, t, 1] <- 0
        ps1[4, i, t, 2] <- 0
        ps1[4, i, t, 3] <- 0
        ps1[4, i, t, 4] <- 0
        ps1[4, i, t, 5] <- phi.4[t] * (1 - alpha.4[t])
        ps1[4, i, t, 6] <- 0
        ps1[4, i, t, 7] <- 0
        ps1[4, i, t, 8] <- 0
        ps1[4, i, t, 9] <- phi.4[t] * alpha.4[t]
        ps1[4, i, t,10] <- 1 - phi.4[t] 
        ps1[5, i, t, 1] <- 0
        ps1[5, i, t, 2] <- 0
        ps1[5, i, t, 3] <- 0
        ps1[5, i, t, 4] <- 0
        ps1[5, i, t, 5] <- 0
        ps1[5, i, t, 6] <- phi.5[t] * (1 - alpha.5[t])
        ps1[5, i, t, 7] <- 0
        ps1[5, i, t, 8] <- 0
        ps1[5, i, t, 9] <- phi.5[t] * alpha.5[t]
        ps1[5, i, t,10] <- 1 - phi.5[t] 
        ps1[6, i, t, 1] <- 0
        ps1[6, i, t, 2] <- 0
        ps1[6, i, t, 3] <- 0
        ps1[6, i, t, 4] <- 0
        ps1[6, i, t, 5] <- 0
        ps1[6, i, t, 6] <- 0
        ps1[6, i, t, 7] <- phi.6[t] * (1 - alpha.6[t])
        ps1[6, i, t, 8] <- 0
        ps1[6, i, t, 9] <- phi.6[t] * alpha.6[t]
        ps1[6, i, t,10] <- 1 - phi.6[t] 
        ps1[7, i, t, 1] <- 0
        ps1[7, i, t, 2] <- 0
        ps1[7, i, t, 3] <- 0
        ps1[7, i, t, 4] <- 0
        ps1[7, i, t, 5] <- 0
        ps1[7, i, t, 6] <- 0
        ps1[7, i, t, 7] <- 0
        ps1[7, i, t, 8] <- phi.7[t] * (1 - alpha.7[t])
        ps1[7, i, t, 9] <- phi.7[t] * alpha.7[t]
        ps1[7, i, t,10] <- 1 - phi.7[t] 
        ps1[8, i, t, 1] <- 0
        ps1[8, i, t, 2] <- 0
        ps1[8, i, t, 3] <- 0
        ps1[8, i, t, 4] <- 0
        ps1[8, i, t, 5] <- 0
        ps1[8, i, t, 6] <- 0
        ps1[8, i, t, 7] <- 0
        ps1[8, i, t, 8] <- 0
        ps1[8, i, t, 9] <- phi.8[t]
        ps1[8, i, t,10] <- 1 - phi.8[t] 
        ps1[9, i, t, 1] <- 0
        ps1[9, i, t, 2] <- 0
        ps1[9, i, t, 3] <- 0
        ps1[9, i, t, 4] <- 0
        ps1[9, i, t, 5] <- 0
        ps1[9, i, t, 6] <- 0
        ps1[9, i, t, 7] <- 0
        ps1[9, i, t, 8] <- 0
        ps1[9, i, t, 9] <- phi.B[t]
        ps1[9, i, t,10] <- 1 - phi.B[t] 
        ps1[10, i, t, 1] <- 0
        ps1[10, i, t, 2] <- 0
        ps1[10, i, t, 3] <- 0
        ps1[10, i, t, 4] <- 0
        ps1[10, i, t, 5] <- 0
        ps1[10, i, t, 6] <- 0
        ps1[10, i, t, 7] <- 0
        ps1[10, i, t, 8] <- 0
        ps1[10, i, t, 9] <- 0
        ps1[10, i, t,10] <- 1
        ## Define probabilities of O(t) given S(t)
        po1[1, i, t, 1] <- 0
        po1[1, i, t, 2] <- 0
        po1[1, i, t, 3] <- 0
        po1[1, i, t, 4] <- 1
        po1[2, i, t, 1] <- 0
        po1[2, i, t, 2] <- p.P1[t]
        po1[2, i, t, 3] <- 0
        po1[2, i, t, 4] <- 1 - p.P1[t]
        po1[3, i, t, 1] <- 0
        po1[3, i, t, 2] <- p.P1[t]
        po1[3, i, t, 3] <- 0
        po1[3, i, t, 4] <- 1 - p.P1[t]
        po1[4, i, t, 1] <- 0
        po1[4, i, t, 2] <- p.P1[t]
        po1[4, i, t, 3] <- 0
        po1[4, i, t, 4] <- 1 - p.P1[t]        
        po1[5, i, t, 1] <- 0
        po1[5, i, t, 2] <- p.P2[t]
        po1[5, i, t, 3] <- 0
        po1[5, i, t, 4] <- 1 - p.P2[t]        
        po1[6, i, t, 1] <- 0
        po1[6, i, t, 2] <- p.P2[t]
        po1[6, i, t, 3] <- 0
        po1[6, i, t, 4] <- 1 - p.P2[t]        
        po1[7, i, t, 1] <- 0
        po1[7, i, t, 2] <- p.P2[t]
        po1[7, i, t, 3] <- 0
        po1[7, i, t, 4] <- 1 - p.P2[t]        
        po1[8, i, t, 1] <- 0
        po1[8, i, t, 2] <- p.P2[t]
        po1[8, i, t, 3] <- 0
        po1[8, i, t, 4] <- 1 - p.P2[t]        
        po1[9, i, t, 1] <- 0
        po1[9, i, t, 2] <- 0
        po1[9, i, t, 3] <- pB[t]
        po1[9, i, t, 4] <- 1 - pB[t]        
        po1[10, i, t, 1] <- 0
        po1[10, i, t, 2] <- 0
        po1[10, i, t, 3] <- 0
        po1[10, i, t, 4] <- 1
      } #t
    } #i
    ##   
    ## LIKELIHOOD
    for(i in 1:nind1){
      # Define latent state at first capture
      z1[i, f1[i]] <- y1[i, f1[i]]
      for(t in (f1[i] + 1):T){
        # State process: draw S(t) given S(t - 1)
        z1[i, t] ~ dcat(ps1[z1[i, t - 1], i, t - 1, ])
        # Observation process: draw O(t) given S(t)
        y1[i, t] ~ dcat(po1[z1[i, t], i, t - 1, ])
      } #t
    } #i
    ## SECOND CAPTURE-RECAPTURE DATA SET; LIFTED FEMALES
    ## State-transition and observation matrices
    for(i in 1:nind2){  
      # Define probabilities of state S(t + 1) given S(t)
      for (t in f2[i]:(T - 1)){
        ps2[1, i, t, 1] <- phiB[t] * (1 - psiBN[t])
        ps2[1, i, t, 2] <- phiB[t] * psiBN[t]
        ps2[1, i, t, 3] <- 1 - phiB[t]
        ps2[2, i, t, 1] <- phiN[t] * psiNB[t]
        ps2[2, i, t, 2] <- phiN[t] * (1 - psiNB[t])
        ps2[2, i, t, 3] <- 1 - phiN[t]
        ps2[3, i ,t, 1] <- 0
        ps2[3, i, t, 2] <- 0
        ps2[3, i, t, 3] <- 1
        #    
        ## Define probabilities of O(t) given S(t)
        po2[1, i, t, 1] <- pB[t]
        po2[1, i, t, 2] <- 0
        po2[1, i, t, 3] <- 1 - pB[t]
        po2[2, i, t, 1] <- 0
        po2[2, i, t, 2] <- pN[t]
        po2[2, i, t, 3] <- 1 - pN[t]
        po2[3, i, t, 1] <- 0
        po2[3, i, t, 2] <- 0
        po2[3, i, t, 3] <- 1
      } #t
    } #i
    # Likelihood 
    for(i in 1:nind2){
      # Define latent state at first capture
      z2[i, f2[i]] <- y2[i, f2[i]]
      for(t in (f2[i] + 1):T){
        # State process: draw S(t) given S(t - 1)
        z2[i, t] ~ dcat(ps2[z2[i, t - 1], i, t - 1, ])
        # Observation process: draw O(t) given S(t)
        y2[i, t] ~ dcat(po2[z2[i, t], i, t - 1, ])
      } #t
    } #i
  }
  ", fill = TRUE)
sink()
##
# 
jags.data <- list(y1=rCH1,y2=rCH2,Prod=Prod[18:42],pF=pFem,nFp=prod0125,f1=f1,f2=f2, 
  T=dim(rCH1)[2],nind1=dim(rCH1)[1],nind2=dim(rCH2)[1],phi0=surToTag,
  z2=known.state.ms(rCH2,3),R=matrix(c(5,0,0,0,1,0,0,0,1),ncol=3))
#
inits <- function()
  list(Omega=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3), 
       mean.psiNB = runif(1, 0, 1), sigmaNB=runif(1, 0, 1),
       mean.pB = runif(1, 0, 1), mean.pN = runif(1, 0, 1), sigmapN = runif(1, 0, 1), mean.pP1 = runif(1, 0, 1), mean.pP2 = runif(1, 0, 1),
       mean.phi.1=runif(1, 0, 1), mean.phi.2=runif(1, 0, 1),mean.phi.3=runif(1, 0, 1),mean.phi.4=runif(1, 0, 1),mean.phi.5=runif(1, 0, 1),mean.phi.6=runif(1, 0, 1),mean.phi.7=runif(1, 0, 1),mean.phi.8=runif(1, 0, 1),mean.phiB = runif(1, 0, 1),
       mean.alpha.3=runif(1, 0, 1),mean.alpha.4=runif(1, 0, 1),mean.alpha.5=runif(1, 0, 1),mean.alpha.6=runif(1, 0, 1),mean.alpha.7=runif(1, 0, 1),
       z1 = z1, z2 = ms.init.z(rCH2, f2),
       miota = rnorm(1, 0.1, 0.5),
       sigma.iota = runif(1, 0.1, 10),
       N1 = round(runif(dim(rCH1)[2], 10, 110), 0),
       N2 = round(runif(dim(rCH1)[2], 10, 87), 0),
       NB3 = round(runif(dim(rCH1)[2], 1, 3), 0),
       NP3 = round(runif(dim(rCH1)[2], 5, 71), 0),
       NB4 = round(runif(dim(rCH1)[2], 2, 15), 0),
       NP4 = round(runif(dim(rCH1)[2], 5, 44), 0),
       NB5 = round(runif(dim(rCH1)[2], 2, 12), 0),
       NP5 = round(runif(dim(rCH1)[2], 5, 22), 0),
       NB6 = round(runif(dim(rCH1)[2], 1, 5), 0),
       NP6 = round(runif(dim(rCH1)[2], 2, 15), 0),
       NB7 = round(runif(dim(rCH1)[2], 1, 3), 0),
       NP7 = round(runif(dim(rCH1)[2], 2, 10), 0),
       NP = round(runif(dim(rCH1)[2], 3, 5), 0),
       NN = round(runif(dim(rCH1)[2], 10, 100), 0),
       NB = round(runif(dim(rCH1)[2], 200, 700), 0),
       NI = round(runif(dim(rCH1)[2], 10, 100), 0))
#
## Parameters monitored
parameters <- c("phi.1","phi.2","phi.3","phi.4","phi.5","phi.6","phi.7","phi.8","phi.B","phiB", "phiN",
                "alpha.3","alpha.4","alpha.5","alpha.6","alpha.7","psiBN", "psiNB", "iota", "eiota",
                "Sigma.phipsiBN", "mean.psiNB", "mean.pP1", "mean.pP2", "mean.pB","pN", "mean.pN",
                "sigmaNB", "sigmapN", "sigma.iota", 
                "Ntot", "NAll", "Nrec", "Npre", "NB", "NN", "NI", "NBI",
                "N1", "N2", "NB3", "NP3", "NB4", "NP4", "NB5", "NP5", "NB6", "NP6", "NB7", "NP7","NP",
                "PR", "lambda", "loglambda", "meanlambda", "NMat.nc", "NBI.nc")
#
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
##
### Run the BUGS model
(sta <- Sys.time())
out.test <- jags(jags.data,inits,parameters,"SSB.IPM",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
## Time difference of 14.78084 days
# 
## OBTAIN WAIC
# samples <- rjags::jags.samples(out.test$model, c("WAIC", "deviance"), type = "mean", n.iter = ni, n.burnin = nb, n.thin = nc)
# samples$p_waic <- samples$WAIC
# samples$waic <- samples$deviance + samples$p_waic
# (tmp <- sapply(samples, sum))
#      WAIC   deviance     p_waic       waic 
## 
# (waic_temp <- round(c(deviance = tmp[["deviance"]], waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]), 2))
# deviance     waic   p_waic 
##
################################################################################
## SAVE WORKSPACE AND WRITE POSTERIOR SIMULATIONS
preds <- rbind(out.test$samples[[1]],out.test$samples[[2]],out.test$samples[[3]])
write.csv(preds,"processed_data/AFSFemaleSimulations200125.csv")
saveRDS(preds,file="processed_data/AFSFemaleSimulations200125.Rds")
## SAVE WORKSPACE
save.image("processed_data/SSBFemaleIPM200125.RData")
# load("processed_data/SSBFemaleIPM200125.RData")
##----------------------------------------------------
##  GENERATE SOME OUTPUT
# MATURE POPULATION IS NBI+NN; PREBREEDERS ARE Npre
# NBI (breeders including immigrants); NN (non-breeders); Npre (prebreeders)
NBIs <- preds[,which(dimnames(preds)[[2]]=="NBI[25]")]
NNs <- preds[,which(dimnames(preds)[[2]]=="NN[25]")]
Npres <- preds[,which(dimnames(preds)[[2]]=="Npre[25]")]
Ntots <- NBIs+NNs+Npres
## SCALING FACTOR TO ACCOUNT FOR PREBREEDERS (preF) 
##  e.g. (NBIs + NNs) * preF scales up the number of adult females to total population size
preF <-  1+Npres/(NBIs+NNs) 
#####
## LOGLAMBDA FOR 2010-2025
loglambdas <- preds[, which(dimnames(preds)[[2]]=="loglambda[9]"):which(dimnames(preds)[[2]]=="loglambda[24]")]
c(mean(exp((1/16)*rowSums(loglambdas))),quantile(exp((1/16)*rowSums(loglambdas)),c(0.025,0.975)))
#                 2.5%     97.5% 
## 0.9506341 0.9431450 0.9580636
##
## LOGLAMBDA FOR 2001-2025
loglambdas <- preds[,which(dimnames(preds)[[2]]=="loglambda[1]"):which(dimnames(preds)[[2]]=="loglambda[24]")]
c(mean(exp((1/24)*rowSums(loglambdas))),quantile(exp((1/24)*rowSums(loglambdas)),c(0.025,0.975)))
#                 2.5%     97.5% 
## 0.9696380 0.9619445 0.9776896 
################################################################################
## POST MODEL FIT
## SOME SMOOTHED SERIES FOR PRESENTATION PURPOSES
NTSims <- preds[,which(dimnames(preds)[[2]]=="NAll[1]"):which(dimnames(preds)[[2]]=="NAll[25]")]
## SOURCE FUNCTION TO OBTAIN A LOW PASS GAUSSIAN FILTER
source("code/glbpf.R")
## RUN THIS BIT ONCE AND SAVE RESULTS
years <- 2001:2025
lpAFS <- gamAFS <- matrix(NA,dim(NTSims)[1],dim(NTSims)[2])
for(i in 1:dim(NTSims)[1]){
  lpAFS[i,] <- glbpf(NTSims[i,],pf=9)
  gamAFS[i,] <- predict(gam(NTSims[i,] ~ s(years)))
  print(i)  
}
saveRDS(NTSims,file="processed_data/NTSimsF.Rds")
saveRDS(lpAFS,file="processed_data/lpAFSF.Rds")
saveRDS(gamAFS,file="processed_data/gamAFSF.Rds")
##
## SOME PLOTS
lpAFS <- readRDS("processed_data/lpAFSF.Rds")
gamAFS <- readRDS("processed_data/gamAFSF.Rds")
NTSims <- readRDS("processed_data/NTSimsF.Rds")
##
## SOURCE COLOUR TRANSPARENCY FUNCTION FOR PLOTTING
source("code/addTrans.R")
## LOW PASS FILTER
plot(2001:2025,colMeans(NTSims),type="n",pch=21,bg=2,ylim=c(250,2000),ylab="SSB females",xlab="")
polygon(c(2001:2025,2025:2001), c(apply(lpAFS,2,quantile,0.975),rev(apply(lpAFS,2,quantile,0.025))),col=addTrans("antiquewhite4",100),border="antiquewhite4")
lines(2001:2025,colMeans(lpAFS),col="antiquewhite4", lwd = 1.1) 
for(i in 1:25)
  lines(rep((2001:2025)[i],2),apply(NTSims,2,quantile,c(0.025,0.975))[,i],col="black")
points(2001:2025,colMeans(NTSims),type ="b",pch=21,col="black",bg="magenta")
## GAM
plot(2001:2025,colMeans(NTSims),type="n",pch=21,bg=2,ylim=c(250,2000),ylab="SSB females",xlab="")
polygon(c(2001:2025,2025:2001), c(apply(gamAFS,2,quantile,0.975),rev(apply(gamAFS,2,quantile,0.025))),col=addTrans("antiquewhite4",100),border="antiquewhite4")
lines(2001:2025,colMeans(gamAFS),col="antiquewhite4", lwd = 1.1) 
for(i in 1:25)
  lines(rep((2001:2025)[i],2),apply(NTSims,2,quantile,c(0.025,0.975))[,i],col="black")
points(2001:2025,colMeans(NTSims),type ="b",pch=21,col="black",bg="magenta")
######