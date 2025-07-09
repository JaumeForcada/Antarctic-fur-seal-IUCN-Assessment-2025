#############################################################################################
##  VITAL RATE, LIFE TABLE AND GENERATION TIME (LENGTH) ESTIMATES FROM FEMALE IPM
#--------------------------------------------------------------------------------------------
## LOAD PACKAGES
library(Rage)
##---------------
## LOAD SAVED WORKSPACE WITH RESULTS
load("processed_data/SSBMaleIPM19942025.RData")
#-----------------------------------
##  READ, FORMAT AND SUBSET SSB DATA
### 1.  TOTAL PUP PRODUCTION COUNTS
##  1.1 PRODUCTIVITY: SSB total pup production from 1983/1984 to 2023/2025
##         84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99  00  01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
Prod <- c(850,678,854,820,777,757,822,545,718,770,714,591,801,733,536,693,446,696,747,713,769,652,464,768,736,782,271,568,547,584,258,525,403,361,496,282,409,368,247,335,303,368)
## productivity 1995:2021
prod9522 <- Prod[12:39]
##
## SEX RATIO ESTIMATION
#  NUMBERS OF KNOWN PUPS SEXED AT BIRTH
#             95  96  97  98  99  00  01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
allPups <- c(198,164,175,119,104, 69,167,141,124,142,316,271,428,293,337,221,284,278,330,205,326,328,243,321,220,280,300,247,340, 34,368) 
femPups <- c( 91, 83, 81, 58, 52, 29, 92, 79, 51, 71,147,135,217,149,164,111,136,139,151,119,164,165,123,168,108,140,160,141,179, 16,171)
malPups <- c(107, 81, 94, 61, 52, 40, 75, 61, 73, 71,169,136,211,144,173,110,149,139,179, 86,162,163,120,153,112,140,140,106,159, 18,197)
## ESTIMATED PROPORTION OF MALES
pMal <- malPups / allPups
pMalSE <- sqrt(pMal * (1 - pMal) / allPups)
## 2.  EARLY SURVIVAL -from birth to tagging (approx. first 1.5 months)
surToTag <- 1 - c(0.176,0.196,0.1473,0.3915,0.1313,0.856,0.265804598,0.18875502,0.130434783,0.239271782,0.197247706,0.155172414,
  0.266927083,0.244565217,0.369565217,0.034090909,0.163732394,0.098540146,0.233276158,0.03875969,0.167619048,0.406947891,0.089709763,
  0.114919355,0.06405694,0.102941176,0.24744898,0.052631579,0.053731343,0.122112211,0.04619565)
## AVERAGE -within year- CV IS 0.1215081
surToTagVAR <- (surToTag*0.1215081)^2
##
## GENERATE ALPHA AND BETA PARAMETERS TO SIMULATE SurToTag
source("code/estBetaParams.R")
STTbeta <- estBetaParams(surToTag, surToTagVAR)
##------------------------------------------------------
## GENERATE SEX RATIOS AND phi0s
getPhi0s <- function(i, j = 10000) rbeta(j,STTbeta$alpha[i],STTbeta$beta[i])
phi0 <- sapply(1:24, getPhi0s)
phi0 <- rowMeans(phi0[,1:11])
getSRs <- function(i, j = 10000) rbinom(j,prod9522[i],pMal[i])/prod9522[i]
sr <- sapply(1:25, getSRs)
sr <- rowMeans(sr[,1:11])
## OTHER WITAL RATES
## LOAD MALE IPM RESULTS
preds <- readRDS(file = "processed_data/AFSMaleSimulations19952025.Rds")
for(i in 1:11)
  assign(paste("phi",i,sep=""),pred_rep[i+1,])
phiN <- preds[,c(which(dimnames(preds)[[2]]=="phiN[1]"):which(dimnames(preds)[[2]]=="phiN[30]"))];dimnames(phiN)[[2]]<-NULL
phiN <- rowMeans(phiN)[1:10000]
phiT <- preds[,c(which(dimnames(preds)[[2]]=="phiT[1]"):which(dimnames(preds)[[2]]=="phiT[30]"))];dimnames(phiT)[[2]]<-NULL
phiT <- rowMeans(phiT[,1:11])[1:10000]
for(i in 1:5)
  assign(paste("a",i+6,sep=""),a[i,])
psiTN <- preds[,c(which(dimnames(preds)[[2]]=="psiTN[1]"):which(dimnames(preds)[[2]]=="psiTN[30]"))];dimnames(psiTN)<-NULL
psiTN <- rowMeans(psiTN[,1:11])[1:10000]
psiNT <- preds[,c(which(dimnames(preds)[[2]]=="psiNT[1]"):which(dimnames(preds)[[2]]=="psiNT[30]"))];dimnames(psiNT)<-NULL
psiNT <- rowMeans(psiNT[,1:11])[1:10000]
##
## FUNCTION TO GENERATE PROJECTION MATRIX A - index i IS FOR SIMULATIONS - SEASONS ARE AVERAGED; 1995 TO 2005
getA <- function(i){
  rbind(c(0,0,0,0,0,0,sr[i]*0.93*phi0[i]*phi7[i]*a7[i],
          sr[i]*0.93*phi0[i]*phi8[i]*a8[i],
          sr[i]*0.93*phi0[i]*phi9[i]*a9[i],
          sr[i]*0.93*phi0[i]*phi10[i]*a10[i],
          sr[i]*0.93*phi0[i]*phi11[i]*a11[i],
          sr[i]*0.93*phi0[i]*phiN[i]*psiNT[i],
          sr[i]*0.93*phi0[i]*phiT[i]*(1-psiTN[i])),
        c(phi1[i],rep(0,12)), c(0,phi2[i],rep(0,11)), c(rep(0,2),phi3[i],rep(0,10)), c(rep(0,3),phi4[i],rep(0,9)), c(rep(0,4),phi5[i],rep(0,8)),
        c(rep(0,5),phi6[i],rep(0,7)),c(rep(0,6),phi7[i]*(1-a7[i]),rep(0,6)),c(rep(0,7),phi8[i]*(1-a8[i]),rep(0,5)),c(rep(0,8),phi9[i]*(1-a9[i]),rep(0,4)),
        c(rep(0,9),phi10[i]*(1-a10[i]),phi11[i]*(1-a11[i]),rep(0,2)),c(rep(0,11),phiN[i]*(1-psiNT[i]),phiT[i]*psiTN[i]),
        c(rep(0,6),phi7[i]*a7[i],phi8[i]*a8[i],phi9[i]*a9[i],phi10[i]*a10[i],phi11[i]*a11[i],phiN[i]*psiNT[i],phiT[i]*(1-psiTN[i])))
}
#
## FUNCTION TO SPLIT MATRIX A INTO COMPONENTS T AND F (BORROWED FROM splitA in library(popbio))
TFfromA <- function (A, r = 1, c = -1){
  tm <- A; fm <- A
  if (is.matrix(r)) {
    if (is.logical(r)) {
      tm[r] <- 0
      fm[!r] <- 0
    }else{
      tm <- A - r
      fm <- r
    }
  }else{
    tm[r, c] <- 0
    fm[-r, ] <- 0
    fm[r, -(c)] <- 0
  }
  list(T = tm, F = fm)
}
## FUNCTION TO ESTIMATE LAMBDA
Lambda <- function(At){
  ev <- eigen(At)
  lmax <- which.max(Re(ev$values))
  Re(ev$values[lmax])
}
#
##
####--------------
## IUCN ASSESSMENT
## GET PRE-DECLINE GENERATION TIME
#  matrix(rowMeans(sapply(rw+k, getA, j = sims[m])), ncol = 9)
TFA <- TFfromA(matrix(rowMeans(sapply(1:8,getA)),ncol=13))
gen_time(TFA$T,TFA$F,method="cohort")
Lambda(matrix(rowMeans(sapply(1:8,getA)),ncol=13))
(LT <- mpm_to_table(TFA$T,TFA$F,xmax=19,remove_final=F))
##
## SET MATRICES FOR SIMULATIONS OF LIFE TABLES AND (COHORT) GENERATION TIME
gentCh <- numeric(10000)
LTarray <- array(NA,dim=c(19,6,10000))
sims <- seq(1,10000,by=10)
for(m in 1:10000){
  TFA <- TFfromA(matrix(rowMeans(sapply(m,getA)),ncol=13))
  gentCh[m] <- gen_time(TFA$T,TFA$F,method="cohort")
  bLT <- mpm_to_table(TFA$T,TFA$F,xmax=18,lx_crit=0.001,remove_final=F)
  LTarray[1:dim(bLT)[1],,m] <- as.matrix(cbind(bLT[,c("x","px","lx","mx","lxmx")],xlxmx=bLT[,c("x")]*bLT[,c("lxmx")]))
  print(m)
}
##
## TO TABULATE RESULTS
## MEAN AND SD OF LIFE TABLE
mLT <- apply(LTarray,c(1,2),mean,na.rm=TRUE)
sLT <- apply(LTarray,c(1,2),sd,na.rm=TRUE)
## MEAN GENERATION TIME
round(c(mean=mean(gentCh),SD=sd(gentCh)),5)
#  mean    SD 
# 9.90445 0.58698
##
#  MEAN MATRIX MODEL
round(matrix(rowMeans(sapply(seq(1,10000,by=10),getA)),ncol=13),3)
#  SD of MEAN MATRIX MODEL
round(matrix(apply(sapply(seq(1,10000,by=10),getA),1,sd),ncol=13),3)
##
## Save data
saveRDS(gentCh,"processed_data/gentChM.Rds")
##
## TO FIND WEIGHTED MEAN
# load 2009 BI mature seal abundance data
BI09m <- readRDS("processed_data/BI09m.Rds")
BI09m <- BI09m[1:10000]
BI09f <- readRDS("processed_data/BI09f.Rds")
BI09f <- BI09f[1:10000]
## load generation times
gentChM <- readRDS("processed_data/gentChM.Rds")
gentChF <- readRDS("processed_data/gentChF.Rds")[1:10000]
gentCHw <- (gentChM*BI09m + gentChF*BI09f)/(BI09m + BI09f) 
round(c(mean=mean(gentCHw),SD=sd(gentCHw)),3)
# mean     SD 
# 8.637 0.258
## Save data
saveRDS(gentCHw,"processed_data/gentCHw.Rds")

