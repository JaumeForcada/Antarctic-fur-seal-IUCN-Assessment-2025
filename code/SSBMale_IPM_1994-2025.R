l################################################################################
## MALE FUR SEAL ABUNDANCE AT BIRD ISLAND'S SPECIAL STUDY BEACH 1995 TO 2025 
##------------------------------------------------------------------------------
## LOAD REQUIRED PACKAGES
library(jagsUI)
library(mgcv)
#------------------------
##  READ, FORMAT AND SUBSET SSB DATA
### 1.  TOTAL PUP PRODUCTION COUNTS
##  1.1 PRODUCTIVITY: SSB total pup production from 1983/1984 to 2023/2024 
##         84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99  00  01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
Prod <- c(850,678,854,820,777,757,822,545,718,770,714,591,801,733,536,693,446,696,747,713,769,652,464,768,736,782,271,568,547,584,258,525,403,361,496,282,409,368,247,335,303,368)
## productivity 1995:2021
prod9522 <- Prod[12:39]
#
##  1.2 PROPORTIONS OF SEALS SEXED AT BIRTH FROM KNOWN MOTHERS (1994/95 to 2021/22)
#             95  96  97  98  99  00  01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22
allPups <- c(198,164,175,119,104, 69,167,141,124,142,316,271,428,293,337,221,284,278,330,205,326,328,243,321,220,280,300,247) 
# females
femPups <- c( 91, 83, 81, 58, 52, 29, 92, 79, 51, 71,147,135,217,149,164,111,136,139,151,119,164,165,123,168,108,140,160,141)
# males
malPups <- c(107, 81, 94, 61, 52, 40, 75, 61, 73, 71,169,136,211,144,173,110,149,139,179, 86,162,163,120,153,112,140,140,106)
## ESTIMATED PROPORTION OF MALES
pMal <- malPups/allPups
pMalSE <- sqrt(pMal*(1-pMal)/allPups)
surToTag <- 1 - c(0.176,0.196,0.1473,0.3915,0.1313,0.856,0.265804598,0.18875502,0.130434783,0.239271782,0.197247706, 
  0.155172414,0.266927083,0.244565217,0.369565217,0.034090909,0.163732394,0.098540146,0.233276158,0.03875969,
  0.167619048,0.406947891,0.089709763,0.114919355,0.06405694,0.102941176,0.24744898,0.052631579)
## AVERAGE -within year- CV IS .15
surToTagVAR <- (surToTag*0.156628419)^2
## GENERATE ALPHA AND BETA (SHAPE/SCALE) PARAMETERS FOR for SurToTag
## SOURCE AD HOC FUNCTION
source("code/estBetaParams.R")
STTbeta <- estBetaParams(surToTag, surToTagVAR)
## SURVIVING MALE PUPS FROM BIRTH TO TAGGING (PRODUCTIVITY * SEX RATIO * SURVIVAL TO TAGGING)
n0 <- prod9522*pMal*surToTag
##
## ESTIMATE SHAPE PARAMETERS FOR BETA DISTRIBUTION FOR SURVIVAL AT AGE
## PREVIOUSLY ESTIMATED MEANS AND VARIANCES
SUR <- cbind(Age = 0:14,estimate = c(0.5000000,0.7196673,0.8132897,0.8405607,0.8303057,0.7894741,0.7192355,0.6270048,0.5327680,
  0.4643934,0.4453887,0.4828134,0.5609565,0.6510524,0.7309625),VAR=c(0.000000e+00,2.123932e-04,1.717949e-04,9.768589e-05,1.140018e-04,
  2.814642e-04,5.913856e-04,9.857080e-04,2.711171e-03,9.966926e-03,2.840808e-02,6.010530e-02,9.488311e-02,1.169055e-01,1.215227e-01))
## ESTIMATE SHAPE PARAMETERS FOR BETA DISTRIBUTION
surBetas <- estBetaParams(mu = SUR[c(2:10,10,10),2],var=SUR[c(2:10,10,10),3])
##
## ESTIMATE SHAPE PARAMETERS FOR BETA DISTRIBUTIONS FOR RECRUITMENT AT AGE OF TERRITORIALS
## PREVIOUSLY ESTIMATED MEANS AND VARIANCES
REC <- cbind(age=7:12,median=C(0.01256698,0.08589635,0.22874003,0.33441262,0.56351099,0.56787121),
  var=c(8.468579e-05,5.257120e-04,1.962870e-03,4.575174e-03,1.275403e-02,3.438484e-02))
## ESTIMATE SHAPE PARAMETERS FOR BETA DISTRIBUTION FOR EARLY SURVIVAL
aBeta <- estBetaParams(REC[, 2], REC[, 3])$alpha
bBeta <- estBetaParams(REC[, 2], REC[, 3])$beta
##
## INITIAL POPULATION VALUES
#  PRELIMINARY JUVENILE SURVIVAL 0.7282531 (SE = 0.01143565)
#   Juvenile survival estimate  
SURV <- 0.7282531
##  initial population sizes
N1m <- round(n0[1] * surToTag[1] * SURV^(1:11))
# m is proportion of territorials at age
# age at first breeding 8-9, all breeding by 10
m <-  c(rep(0, 7), 0.2, 0.4, 1)
##
N1 <- round(n0[1] * surToTag[1] * SURV)
##
###-----------------------------------------------------------------------------
### MULTI-STATE JOLLY-SEBER MODEL STRUCTURE
###-----------------------------------------------------------------------------
##  TEMPORARY EMIGRATION JOLLY-SEBER MODEL TO ESTIMATE REAL ABUNDANCE OF 
##  TERRITORIALS AND NUMBERS OF RECRUITED TERRITORIALS THAT ARE ALIVE AND AWAY
##  Jolly-Seber Multi-state Model Description
### States:
##  P: pre-breeder (non-observable)   1
##  B: breeder                        2
##  N: non-breeder (non-observable)   3
##  D: dead        (non-observable)   4
### Parameters: 
##  gamma -> 'removal' entry probability
##  phi -> survival (for B and N)
##  psiBN -> transition B -> N
##  psiNB -> transition N -> B
##  p -> recapture (B only)
### Transitions Matrix
#               pre-breeder         breeder           non-breeder         dead
# pre-breeder   1 - gamma           gamma             0                   0
# breeder       0                   phiT*(1-psiBN)    phiT*psiBN          1-phiT 
# non-breeder   0                   phiN*psiNB        phiN*(1-psiNB)      1-phiN
# dead          0                   0                 0                   1
### Recapture Matrix
#               seen  not-seen
# pre-breeder   0     1
# breeder       p     1 - p 
# non-breeder   0     1
# dead          0     1
##------------------------------------------------------------------------------
## 13 seasons, 1994/1995 to 2006/07 [t]
## 11+ Age classes [a]
## The 11+ group is the NT estimate from the Jolly-Seber multi-state model
## The Jolly-Seber model has a temporary emigration structure and produces
##  NT[t] - Territorial males
##  NN[t] - Non-territorial males alive
##  B[t] - new recruits (including SSB born pups and immigrants)
## The integrated count model produces 
##  N[t, a] - over 11 age classes (nAges = 11)
##  N[t, a] = NP[t, a] + NM[t, a]
##  where NP are pre-territorials and NM are maturing territorials.
##  Nall[t] and NMall[t] are sum over ages of N and NM respectively.
##  NI[t] are new territorial immigrants
##  with B[t] = NMall[t] + NI[t]   
####
##
#   JAGS MODEL
sink("AFSMjsms.jags")
cat("
  model{
    ### PRIORS
    # RECRUITMENT PROBABILITY (mature seals = m)
    for(i in 1:6){
      m[i] <- 0
    }
    for(i in 1:5){
      mm[i] ~ dbeta(aBeta[i], bBeta[i])                 # priors parameters for beta distributions
      m[i + 6] <- mm[i]
    }
    ## INITIAL POPULATION SIZES
    for(i in 1:Nages){
      n[1, i] ~ dnorm(N1m[i], 0.001)T(0, )
      N[1, i] <- round(n[1, i])
    }
    ## SEX RATIO
    for(i in 1:nyears){
      sr[i] <- srp[i] / np[i]
      srp[i] ~ dbin(pM[i], np[i])
    }
    ## PUP AND JUVENILE SURVIVAL PRIORS
    for(t in 1:nyears){
      phi0[t] ~ dbeta(alphaS[t], betaS[t])
      phi1[t] ~ dbeta(267.107,233.487)
      # phij[t] ~ dbeta(1075.389,790.8746)
      for(age in 1:Nages){
        phij[t, age] ~ dbeta(surAlpha[age],surBeta[age])
      }
    }
    ## ESTIMATED POPULATION SIZE OF PRETERRITORIALS AND NEW RECRUITS TOGETHER
    for(t in 1:nyears){
      Nall[t] <- sum(N[t, ])              #### OLD Nestim[t]
      NMall[t] <- inprod(N[t, ], m[])     #### OLD Nestim2[t]
    }
    # Vectors of population size classes by age 
    for(t in 1:nyears){
      for(age in 1:Nages){
        NM[t, age] <- N[t, age] * m[age]     # NM are maturing individuals
      }
    }
    for(t in 1:nyears){
      for(age in 1:Nages){
        NP[t, age] <- N[t, age] - NM[t, age]  # NP are non-maturing individuals
      }
    }
    for(t in 1:nyears){
      NPall[t] <- sum(NP[t, ])
    } 
    #
    ### Likelihoods
    ##  Population count data (state-space model)
    #   System process
    for(t in 2:nyears){
      mean1[t] <- sr[t] * np[t] * phi0[t] * phi1[t]
      N[t, 1] ~ dpois(mean1[t])
      for(age in 2:Nages){
        meanN[t, age] <- phij[t - 1, age - 1] * N[t - 1, age - 1]
        N[t, age] ~ dpois(meanN[t, age])
      }
    }
    #
    ## JOLLY-SEBER MULTISTATE MODEL
    # Priors
    for(t in 1:(nyears - 1)){
      logit(phiT[t]) <- muT + epsilonT[t]       # temporal random efects on survival of territorials
      epsilonT[t] ~ dnorm(0, tauT)T(-15,15)
      phiN[t] <- mean.phiN                      # constant survival of N
      psiTN[t] <- mean.psiTN                    # Constant mean.psiNT
      psiNT[t] <- mean.psiNT                    # Constant mean.psiNT
      gamma[t] ~ dunif(0, 1)                    # Prior for return
      p[t] <- 1
      # $mean.p is [1] 0.9922023 -> fixed to 1
      # logit(p[t]) <- muP + epsilonP[t]        # temporal random efects on recapture of territorials
      # epsilonP[t] ~ dnorm(0, tauP)T(-15,15)
    }
    #
    ## survival and psiTN/NT
    muT <- log(mean.phiT / (1 - mean.phiT))
    mean.phiT ~ dunif(0, 1)                     # Prior for mean survival of B
    tauT <- pow(sigmaT, -2)
    sigmaT ~ dunif(0, 10)                       # Prior on standard deviation
    sigma2T <- pow(sigmaT, 2)                   # Temporal variance
    #
    mean.phiN ~ dunif(0, 1)  
    #
    mean.psiTN ~ dunif(0, 1)                    # Prior for mean temporary emigration
    #
    mean.psiNT ~ dunif(0, 1)                    # Prior for mean temporary immigration
    #
    #muP <- log(mean.p / (1 - mean.p))         # Logit transformation
    # mean.p ~ dunif(0, 1)                      # Prior for mean recapture
    #tauP <- pow(sigmaP, -2)
    #sigmaP ~ dunif(0, 5)                      # Prior for standard deviation
    #sigma2P <- pow(sigmaP, 2)
    #
    ## transition and observation matrices   
    for (i in 1:M){  
      # probabilities of state S(t+1) given S(t)
      for (t in 1:(nyears-1)){
        ps[1, i, t, 1] <- 1 - gamma[t]
        ps[1, i, t, 2] <- gamma[t]
        ps[1, i, t, 3] <- 0
        ps[1, i, t, 4] <- 0
        ps[2, i, t, 1] <- 0
        ps[2, i, t, 2] <- phiT[t] * (1 - psiTN[t])
        ps[2, i, t, 3] <- phiT[t] * psiTN[t]
        ps[2, i, t, 4] <- 1 - phiT[t]
        ps[3, i, t, 1] <- 0
        ps[3, i, t, 2] <- phiN[t] * psiNT[t]
        ps[3, i, t, 3] <- phiN[t] * (1 - psiNT[t])
        ps[3, i, t, 4] <- 1 - phiN[t]
        ps[4, i, t, 1] <- 0
        ps[4, i, t, 2] <- 0
        ps[4, i, t, 3] <- 0
        ps[4, i, t, 4] <- 1
        #
        ## probabilities of O(t) given S(t)
        po[1, i, t, 1] <- 0
        po[1, i, t, 2] <- 1
        po[2, i, t, 1] <- p[t]
        po[2, i, t, 2] <- 1 - p[t]
        po[3, i, t, 1] <- 0
        po[3, i, t, 2] <- 1
        po[4, i, t, 1] <- 0
        po[4, i, t, 2] <- 1
      }
    }
    ## Likelihood 
    for(i in 1:M){
      # latent state at first occasion
      z[i, 1] <- 1   # all M seals are in state 1 (pre-breeder) at t = 1
      for(t in 2:nyears){
        # State process: draw S(t) given S(t-1)
        z[i, t] ~ dcat(ps[z[i, t - 1], i, t - 1, ])
        # Observation process: draw O(t) given S(t)
        y[i, t] ~ dcat(po[z[i, t], i, t - 1, ])
      }
    }
    # Derived population parameters
    for(t in 1:(nyears - 1)){
      qgamma[t] <- 1 - gamma[t]
    }
    cprob[1] <- gamma[1]
    for(t in 2:(nyears - 1)){
      cprob[t] <- gamma[t] * prod(qgamma[1:(t - 1)])
    }
    psi <- sum(cprob[])                          # Inclusion probability
    for(t in 1:(nyears - 1)){
      b[t] <- cprob[t] / psi                     # Entry probability
    }
    for(i in 1:M){
      for(t in 2:nyears){
        alt[i, t - 1] <- equals(z[i, t], 2)
        aln[i, t - 1] <- equals(z[i, t], 3)
      }
      for(t in 1:(nyears - 1)){
        d[i, t] <- equals(z[i, t] - (alt[i, t] + aln[i, t]), 0)
      }
      alive[i] <- sum(alt[i, ] + aln[i, ])
    }
    for(t in 1:(nyears - 1)){
      NT[t] <- sum(alt[, t])                          # Actual population size of territorials
      NN[t] <- sum(aln[, t])                          # Actual population size of non-territorials
      NATot[t] <- NT[t] + NN[t]                       # Actual mature population size
      B[t] <- sum(d[, t])                             # Number of entries - recruitment and immigration-
    }
    for (i in 1:M){
      w[i] <- 1 - equals(alive[i], 0)
    }
    Nsuper <- sum(w[])                                # Superpopulation size
    ## END OF JOLLY-SEBER
    ##
    ## DERIVED PARAMETERS
    for(t in 1:(nyears - 1)){
      NI[t] <- max(B[t],B[t]-NMall[t])                    # immigrants
      Ntot[t] <- NATot[t] + NPall[t]                  # total population size
    }
    for(t in 1:(nyears - 2)){
      lambda[t] <- Ntot[t+1] / Ntot[t]                # population growth rate
      loglambda[t] <- log(lambda[t])                  # log-lambda
    }
    meanlambda <- exp((1/(nyears-2))*sum(loglambda[1:(nyears-2)]))  # geometric mean
    #
  }
  ",fill = TRUE)
sink()
#
###
## Import capture histories file
## 0 before first capture
## 1 seen breeding
## 2 not seen after seen breeding on previous occasions
### individual breeding data from SSB for the period 1994 to 2020 (1995/95 to 2020/21)
CH <- as.matrix(read.csv("AFSmale_histories.csv", header = T))
CH <- CH[which(rowSums(CH) != 0), ]
# then continue with other years...
nrows <- dim(CH)[1]
ncols <- dim(CH)[2]
for(i in 1:nrows){
  id <- as.vector(which(CH[i, ] == 1)[1])
  if(id > 1)
    CH[i, 1:(id-1)] <- NA
}
CH[CH == 0] <- 2                            # change 0s to 2s; i.e. Not seen = 2, seen = 1
CH[is.na(CH)] <- 0                          # change NAs to 0s; i.e. Not seen before breeding starts
CH.du <- cbind(rep(0, dim(CH)[1]), CH)      # Add dummy occasion
nz <- 500                                   # Augment data
CH.ms <- rbind(CH.du, matrix(0, ncol = dim(CH.du)[2], nrow = nz))
CH.ms[CH.ms == 0] <- 2                      # change 0s to 2s; i.e. Not seen = 2, seen = 1
colnames(CH.ms) <- NULL
#
# DATA
##
jags.data <- list(Nages=11,N1m=N1m,aBeta=aBeta,bBeta=bBeta,pM=pMal,np=prod9522,
  alphaS=STTbeta$alpha,betaS=STTbeta$beta,surAlpha=surBetas$alpha,surBeta=surBetas$beta,
  y=CH.ms,nyears=dim(CH.ms)[2],M=dim(CH.ms)[1])
#
## INITIAL VALUES
###  Initial values need to be specified for the latent state z must correspond to the true state; not the same as the observed state.
###  "1" for all latent states before an individual is observed
###  "2" when the individual was observed alive and nesting 
###  "3" between observations
###  "4" after the last observation.
##
# code to generate such values after data augmentation
js.ms.init <- function(ch, nz){
  # define states
  state <- ch
  state[state == 2] <- NA
  # "2" individual observed alive and nesting 
  state[state == 1] <- 2
  for(i in 1:(dim(state)[1] - nz)){
    n1 <- min(which(state[i, ] == 2))
    n2 <- max(which(state[i, ] == 2))
    n3 <- dim(state)[2]
    if(n2 > n1)
      state[i, (n1 + 1):n2] <- ifelse(is.na(state[i, (n1 + 1):n2]), 3, 2)
    state[i, 1:(n1 - 1)] <- 1
    if(n3 > n2)
      state[i, (n2 + 1):n3] <- 4  
  }
  state[which(is.na(state))] <- 1
  state[, 1] <- NA
  return(state)
}
#
## initial values
inits <- function(){
  list(mean.phiT = runif(1, 0, 1),
       mean.phiN = runif(1, 0, 1),
       mean.psiTN = runif(1, 0, 1), 
       mean.psiNT = runif(1, 0, 1), 
       # mean.p = runif(1, 0, 1), 
       sigmaT = runif(1, 0, 1),
       sigmaP = runif(1, 0, 1),
       z = js.ms.init(CH.ms, nz))
}    
#
## MCMC settings
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 100000 
nt <- 6
nb <- 70000
nc <- 3
#
## Parameters monitored
params <- c("phij","phiT","mean.phiT","mean.phiN","mean.psiTN","mean.psiNT",  # "mean.p",
  "sigma2T","m","b","Nsuper","NT","NN","NATot","B","NM","NP","NMall","NPall","Ntot","NI", 
  "loglambda", "meanlambda")
#
sta <- Sys.time()
out <- jags(jags.data,inits,params,"AFSMjsms.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T) 
Sys.time() - sta
##
#  Time difference of 1.784385 days
#####
## WAIC
#samples <- rjags::jags.samples(out.0$model, c("WAIC", "deviance"), type = "mean", n.iter = ni, n.burnin = nb, n.thin = nc)
#Sys.time() - sta
## ...
#samples$p_waic <- samples$WAIC
#samples$waic <- samples$deviance + samples$p_waic
#(tmp <- sapply(samples, sum))
##     WAIC  deviance    p_waic      waic 
# (waic_temp <- round(c(deviance = tmp[["deviance"]], waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]), 2))
## deviance     waic   p_waic 
##------------------------------------------------------------------------------
## POST-MODEL AFTER FIRST MODEL IS FITTED
## -THIS IS NECESSARY DUE TO LACK OF GENOTYPE ANALYSIS SINCE 2020/2021-
##
## EXTENDED VECTOR OF OBSERVED TERRITORIAL MALES
obster <- c(as.vector(colSums(as.matrix(read.csv("AFSmale_histories.csv", header = T)))),c(c(71,90,87)*.85,49))
obster[1:27]/out$mean$Ntot
obster[28:31]/(obster[1:27]/out$mean$Ntot)[26]
##
##EXTENDED PRODUCTIVITY SERIES
prod9525 <- Prod[12:42]
# number of known pups sexed at birth
#             95  96  97  98  99  00  01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
allPups2 <- c(198,164,175,119,104, 69,167,141,124,142,316,271,428,293,337,221,284,278,330,205,326,328,243,321,220,280,300,247,340, 34,368) 
# females
femPups2 <- c( 91, 83, 81, 58, 52, 29, 92, 79, 51, 71,147,135,217,149,164,111,136,139,151,119,164,165,123,168,108,140,160,141,179, 16,171)
# males
malPups2 <- c(107, 81, 94, 61, 52, 40, 75, 61, 73, 71,169,136,211,144,173,110,149,139,179, 86,162,163,120,153,112,140,140,106,159, 18,197)
## ESTIMATED PROPORTION OF MALES
pMal2 <- malPups2 / allPups2
pMalSE2 <- sqrt(pMal2 * (1 - pMal2) / allPups2)
## EARLY SURVIVAL -from birth to tagging (approx. first 1.5 months)
surToTag2 <- 1 - c(0.176,0.196,0.1473,0.3915,0.1313,0.856,0.265804598,0.18875502,0.130434783,0.239271782,0.197247706, 
  0.155172414,0.266927083,0.244565217,0.369565217,0.034090909,0.163732394,0.098540146,0.233276158,0.03875969,
  0.167619048,0.406947891,0.089709763,0.114919355,0.06405694,0.102941176,0.24744898,0.052631579,0.053731343, 0.122112211,0.04619565)
## AVERAGE -within year- CV IS .20
surToTagVAR2 <- (surToTag2*0.1215081)^2
## generate alpha and beta parameters for SurToTag
STTbeta2 <- estBetaParams(surToTag2, surToTagVAR2)
##
STbeta <- estBetaParams(c(out$mean$phiT,esTSm),c(out$sd$phiT^2,esTSv))
###
#   JAGS MODEL
sink("AFSM2.jags")
cat("
  model{
    ### PRIORS
    # RECRUITMENT PROBABILITY (mature seals = m)
    for(i in 1:6){
      m[i] <- 0
    }
    for(i in 1:5){
      mm[i] ~ dbeta(aBeta[i],bBeta[i])                 # priors parameters for beta distributions
      m[i + 6] <- mm[i]
    }
    ## INITIAL POPULATION SIZES
    for(i in 1:Nages){
      n[1, i] ~ dnorm(N1m[i], 0.001)T(0, )
      N[1, i] <- round(n[1, i])
    }
    ## SEX RATIO
    for(i in 1:nyears){
      sr[i] <- srp[i] / np[i]
      srp[i] ~ dbin(pM[i], np[i])
    }
    ## PUP AND JUVENILE SURVIVAL PRIORS
    for(t in 1:nyears){
      phi0[t] ~ dbeta(alphaS[t], betaS[t])
      phi1[t] ~ dbeta(267.107,233.487)
      # phij[t] ~ dbeta(1075.389,790.8746)
      for(age in 1:Nages){
        phij[t, age] ~ dbeta(surAlpha[age],surBeta[age])
      }
    }
    ## ESTIMATED POPULATION SIZE OF PRETERRITORIALS AND NEW RECRUITS TOGETHER
    for(t in 1:nyears){
      Nall[t] <- sum(N[t,])
      NMall[t] <- inprod(N[t,],m[])
    }
    # Vectors of population size classes by age 
    for(t in 1:nyears){
      for(age in 1:Nages){
        NM[t,age] <- N[t,age]*m[age]     # NM are maturing individuals
      }
    }
    for(t in 1:nyears){
      for(age in 1:Nages){
        NP[t,age] <- N[t,age]-NM[t,age]  # NP are non-maturing individuals
      }
    }
    for(t in 1:nyears){
      NPall[t] <- sum(NP[t,])
      NP8[t] <- NP[t,1]+NP[t,2]+NP[t,3]+NP[t,4]+NP[t,5]+NP[t,6]+NP[t,7]+NP[t,8]
    } 
    #
    ## PRIORS FOR AGM COMPONENT
    for(t in 1:(nyears-1)){
      phiT[t] ~ dbeta(STS[t,1],STS[t,2])T(0.001,0.999)
      phiN[t] ~ dbeta(2.943019,2.325703)T(0.001,0.999)
      psiTN[t] ~ dbeta(0.1452947,0.9112372)T(0.001,0.999)
      psiNT[t] ~ dbeta(0.6212444,0.1983089)T(0.001,0.999)
    }
    #
    ## Initial ADULT population sizes
    NN[1] ~ dnorm(11,0.03)T(0,)
    # NT[1] ~ dnorm(127,1)T(0,)
    #
    ### Likelihoods
    ##  Population count data (state-space model)
    #   System process 1
    for(t in 2:nyears){
      mean1[t] <- sr[t]*np[t]*phi0[t]*phi1[t]
      N[t, 1] ~ dpois(mean1[t])
      for(age in 2:Nages){
        meanN[t,age] <- phij[t-1,age-1]*N[t-1,age-1]
        N[t,age] ~ dpois(meanN[t,age])
      }
    }
    ### system process 2
    for(t in 2:nyears){
      meanNN[t] <- NT[t-1]*phiT[t-1]*psiTN[t-1]+NN[t-1]*phiN[t-1]*(1-psiNT[t-1])
      NN[t] ~ dpois(meanNN[t])
    }
    ## derived parameters
    for(t in 1:nyears){
      Ntot[t] <- NT[t]+NN[t]       # total mature population size
      NTOT[t] <- Ntot[t]+NPall[t]
    }
    for(t in 1:(nyears-1)){
      lambda[t] <- Ntot[t+1]/Ntot[t]  # population growth rate
      loglambda[t] <- log(lambda[t])  # log-lambda
    }
    meanlambda <- exp((1/(nyears-1))*sum(loglambda[1:(nyears-1)]))  # geometric mean
  }
  ",fill = TRUE)
sink()
#
##
# DATA
##
jags.data2 <- list(Nages=11,N1m=N1m,aBeta=aBeta,bBeta=bBeta,pM=pMal2,np=prod9525,nyears=length(prod9525),
  alphaS=STTbeta2$alpha,betaS=STTbeta2$beta,surAlpha=surBetas$alpha,surBeta=surBetas$beta,
  STS=matrix(unlist(STbeta),ncol=2),NT=round(obster))
#
## initial values
inits2 <- function(){
  list(mm=REC[1:5,2])
}
#
## MCMC settings
## ITERATIONS (ni-nb)/nt*nc -> posterior samples.
ni <- 100000 
nt <- 6
nb <- 70000
nc <- 3
##
## Parameters monitored
params2 <- c("phiT","phiN","psiTN","psiNT","loglambda","meanlambda",
  "m","NM","NP","NMall","NPall","NTOT","Ntot","NN","NT","phij","NP8")
#
sta <- Sys.time()
out2 <- jags(jags.data2,inits2,params2,"AFSM2.jags",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T) 
Sys.time() - sta
#  Time difference of 1.295167 mins
##
################################################################################
## SAVE WORKSPACE AND WRITE POSTERIOR SIMULATIONS
save.image("processed_data/SSBMaleIPM19942025.RData")
preds <- rbind(out2$samples[[1]],out2$samples[[2]],out2$samples[[3]])
write.csv(preds, "AFSMaleSimulations19952025.csv")
saveRDS(preds,file="AFSMaleSimulations19952025.Rds")
##---------------------------------------------------
## POST MODEL FIT
## SOME SMOOTHED SERIES FOR PRESENTATION PURPOSES
years <- 2001:2025
lpAFS <- gamAFS <- matrix(NA,dim(NTSims)[1],dim(NTSims)[2])
for(i in 1:dim(NTSims)[1]){
  lpAFS[i,] <- glbpf(NTSims[i,],pf=4.8)
  gamAFS[i,] <- predict(gam(NTSims[i,] ~ s(years)))
  print(i)  
}
saveRDS(lpAFS,file="lpAFSM.Rds")
saveRDS(gamAFS,file="gamAFSM.Rds")
lpAFSM <- readRDS("lpAFSM.Rds")
gamAFSM <- readRDS("gamAFSM.Rds")
NTSims <- readRDS("NTSimsM.Rds")
##
NTSims9525 <- preds[,which(dimnames(preds)[[2]]=="Ntot[1]"):which(dimnames(preds)[[2]]=="Ntot[31]")]
saveRDS(NTSims9525,file="NTSims9525.Rds")
years <- 1995:2025
lpAFSM9525 <- gamAFSM9525 <- matrix(NA,dim(NTSims9525)[1],dim(NTSims9525)[2])
for(i in 1:dim(NTSims9525)[1]){
  lpAFSM9525[i,] <- glbpf(NTSims9525[i,],pf=6)
  gamAFSM9525[i,] <- predict(gam(NTSims9525[i,] ~ s(years)))
  print(i)  
}
saveRDS(lpAFSM9525,file="lpAFSM9525.Rds")
saveRDS(gamAFSM9525,file="gamAFSM9525.Rds")
##
## PLOTS
plot(2001:2025,colMeans(NTSims),type="n",pch=21,bg=2,ylim=c(300,900),ylab="SSB males",xlab="")
polygon(c(2001:2025,2025:2001), c(apply(lpAFSM,2,quantile,0.975),rev(apply(lpAFSM,2,quantile,0.025))),col=addTrans("antiquewhite4",100),border="antiquewhite4")
lines(2001:2025,colMeans(lpAFS),col="antiquewhite4", lwd = 1.1) 
for(i in 1:25)
  lines(rep((2001:2025)[i],2),apply(NTSims,2,quantile,c(0.025,0.975))[,i],col="black")
points(2001:2025,colMeans(NTSims),type ="b",pch=21,col="black",bg="magenta")
##
#
plot(2001:2025,colMeans(NTSims),type="n",pch=21,bg=2,ylim=c(300,900),ylab="SSB mature females",xlab="")
polygon(c(2001:2025,2025:2001), c(apply(gamAFSM,2,quantile,0.975),rev(apply(gamAFSM,2,quantile,0.025))),col=addTrans("antiquewhite4",100),border="antiquewhite4")
lines(2001:2025,colMeans(gamAFSM),col="antiquewhite4", lwd = 1.1) 
for(i in 1:25)
  lines(rep((2001:2025)[i],2),apply(NTSims,2,quantile,c(0.025,0.975))[,i],col="black")
points(2001:2025,colMeans(NTSims),type ="b",pch=21,col="black",bg="magenta")
##
##--------------------------
save.image("processed_data/SSBMaleIPM19942025.RData")
##
