#############################################################################################
##  VITAL RATE, LIFE TABLE AND GENERATION TIME (LENGTH) ESTIMATES FROM FEMALE IPM
#--------------------------------------------------------------------------------------------
## LOAD PACKAGES
library(Rage)
##---------------
## LOAD SAVED WORKSPACE WITH RESULTS
load("processed_data/SSBFemaleIPM200125.RData")
#-----------------------------------
##  READ, FORMAT AND SUBSET SSB DATA
### 1.  TOTAL PUP PRODUCTION COUNTS
##  1.1 PRODUCTIVITY: SSB total pup production from 1983/1984 to 2023/2025
##         84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99  00  01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
prod <- c(850,678,854,820,777,757,822,545,718,770,714,591,801,733,536,693,446,696,747,713,769,652,464,768,736,782,271,568,547,584,258,525,403,361,496,282,409,368,247,335,303,368)
prod0125 <- prod[18:42] # pup production vector
##
## SEX RATIO ESTIMATION
#  NUMBERS OF KNOWN PUPS SEXED AT BIRTH
##            01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
allPups <- c(167,141,124,142,316,271,428,293,337,221,284,278,330,205,326,328,243,321,220,280,300,247,340, 34,368)
femPups <- c( 92, 79, 51, 71,147,135,217,149,164,111,136,139,151,119,164,165,123,168,108,140,160,141,179, 16,171) # Females
malPups <- c( 75, 61, 73, 71,169,136,211,144,173,110,149,139,179, 86,162,163,120,153,112,140,140,106,159, 18,197) # Males
## PROPORTION OF FEMALES (MEAN AND SE)
pFem <- femPups/allPups
pFemSE <- sqrt(pFem*(1-pFem)/allPups)
##
## 2.  EARLY SURVIVAL -from birth to tagging (approx. first 1.5 months)
surToTag <- 1 - c(0.265804598, 0.18875502, 0.130434783, 0.239271782, 0.197247706, 0.155172414, 0.266927083,
                  0.244565217, 0.369565217, 0.034090909, 0.163732394, 0.098540146, 0.233276158, 0.03875969,
                  0.167619048, 0.406947891, 0.089709763, 0.114919355, 0.06405694, 0.102941176, 0.24744898,
                  0.052631579, 0.053731343, 0.122112211, 0.04619565)
## AVERAGE -within year- CV IS 0.1215081
surToTagVAR <- (surToTag*0.1215081)^2
##
## GENERATE ALPHA AND BETA PARAMETERS TO SIMULATE SurToTag
source("code/estBetaParams.R")
STTbeta <- estBetaParams(surToTag, surToTagVAR)
##------------------------------------------------------
## GENERATE SEX RATIOS AND phi0s
getPhi0s <- function(i, j = 15000) rbeta(j,STTbeta$alpha[i],STTbeta$beta[i])
phi0 <- sapply(1:24, getPhi0s)
getSRs <- function(i, j = 15000) rbinom(j,prod0125[i],pFem[i])/prod0125[i]
sr <- sapply(1:25, getSRs)
## OTHER VITAL RATES
## LOAD FEMALE IPM RESULTS
preds <- readRDS(file ="processed_data/AFSFemaleSimulations200125.Rds")
phi1 <- preds[, c(which(dimnames(preds)[[2]] == "phi.1[1]"):which(dimnames(preds)[[2]] == "phi.1[24]"))];dimnames(phi1)[[2]] <- NULL
phi2 <- preds[, c(which(dimnames(preds)[[2]] == "phi.2[1]"):which(dimnames(preds)[[2]] == "phi.2[24]"))];dimnames(phi2)[[2]] <- NULL
phi3 <- preds[, c(which(dimnames(preds)[[2]] == "phi.3[1]"):which(dimnames(preds)[[2]] == "phi.3[24]"))];dimnames(phi3)[[2]] <- NULL
phi4 <- preds[, c(which(dimnames(preds)[[2]] == "phi.4[1]"):which(dimnames(preds)[[2]] == "phi.4[24]"))];dimnames(phi4)[[2]] <- NULL
phi5 <- preds[, c(which(dimnames(preds)[[2]] == "phi.5[1]"):which(dimnames(preds)[[2]] == "phi.5[24]"))];dimnames(phi5)[[2]] <- NULL
phi6 <- preds[, c(which(dimnames(preds)[[2]] == "phi.6[1]"):which(dimnames(preds)[[2]] == "phi.6[24]"))];dimnames(phi6)[[2]] <- NULL
phi7 <- preds[, c(which(dimnames(preds)[[2]] == "phi.7[1]"):which(dimnames(preds)[[2]] == "phi.7[24]"))];dimnames(phi7)[[2]] <- NULL
phi.B <- preds[, c(which(dimnames(preds)[[2]] == "phi.B[1]"):which(dimnames(preds)[[2]] == "phi.B[24]"))];dimnames(phi.B)[[2]] <- NULL
phiB <- preds[, c(which(dimnames(preds)[[2]] == "phiB[1]"):which(dimnames(preds)[[2]] == "phiB[24]"))];dimnames(phiB)[[2]] <- NULL
phiB <- (phiB + phi.B) / 2
phiN <- preds[, c(which(dimnames(preds)[[2]] == "phiN[1]"):which(dimnames(preds)[[2]] == "phiN[24]"))];dimnames(phiN)[[2]] <- NULL
a3 <- preds[, c(which(dimnames(preds)[[2]] == "alpha.3[1]"):which(dimnames(preds)[[2]] == "alpha.3[24]"))];dimnames(a3)<-NULL
a4 <- preds[, c(which(dimnames(preds)[[2]] == "alpha.4[1]"):which(dimnames(preds)[[2]] == "alpha.4[24]"))];dimnames(a4)<-NULL
a5 <- preds[, c(which(dimnames(preds)[[2]] == "alpha.5[1]"):which(dimnames(preds)[[2]] == "alpha.5[24]"))];dimnames(a5)<-NULL
a6 <- preds[, c(which(dimnames(preds)[[2]] == "alpha.6[1]"):which(dimnames(preds)[[2]] == "alpha.6[24]"))];dimnames(a6)<-NULL
a7 <- preds[, c(which(dimnames(preds)[[2]] == "alpha.7[1]"):which(dimnames(preds)[[2]] == "alpha.7[24]"))];dimnames(a7)<-NULL
psiBN <- preds[, c(which(dimnames(preds)[[2]] == "psiBN[1]"):which(dimnames(preds)[[2]] == "psiBN[24]"))];dimnames(psiBN)<-NULL
psiNB <- preds[, c(which(dimnames(preds)[[2]] == "psiNB[1]"):which(dimnames(preds)[[2]] == "psiNB[24]"))];dimnames(psiNB)<-NULL
io <- preds[, c(which(dimnames(preds)[[2]] == "iota[1]"):which(dimnames(preds)[[2]] == "iota[24]"))];dimnames(io)<-NULL
##--------------------------------------------
## FUNCTION TO GENERATE PROJECTION MATRIX A - j (rows) IS FOR SIMULATIONS, i (columns) IS FOR SEASONS
getA <- function(i, j = 1){
  rbind(c(0,0, sr[j,i]*phi0[j,i]*phi3[j,i]*a3[j,i],
               sr[j,i]*phi0[j,i]*phi4[j,i]*a4[j,i],
               sr[j,i]*phi0[j,i]*phi5[j,i]*a5[j,i],
               sr[j,i]*phi0[j,i]*phi6[j,i]*a6[j,i],
               sr[j,i]*phi0[j,i]*phi7[j,i]*a7[j,i],
               sr[j,i]*phi0[j,i]*phiN[j,i]*psiNB[j,i],
               sr[j,i]*phi0[j,i]*phiB[j,i]*(1-psiBN[j,i])),
             c(phi1[j,i],rep(0,8)),c(0,phi2[j,i],rep(0,7)),c(0,0,phi3[j,i]*(1-a3[j,i]),rep(0,6)),
             c(rep(0,3),phi4[j,i]*(1-a4[j,i]),rep(0,5)),c(rep(0,4),phi5[j,i]*(1-a5[j,i]),rep(0,4)),
             c(rep(0,5),phi6[j,i]*(1-a6[j,i]),phi7[j,i]*(1-a7[j,i]),rep(0,2)),
             c(rep(0,7),phiN[j,i]*(1-psiNB[j,i]),phiB[j,i]*psiBN[j,i]),
             c(0,0,phi3[j,i]*a3[j,i],phi4[j,i]*a4[j,i],phi5[j,i]*a5[j,i],phi6[j,i]*a6[j,i],phi7[j,i]*a7[j,i],phiN[j,i]*psiNB[j,i],phiB[j,i]*(1-psiBN[j,i])))
}
#
## FUNCTION TO SPLIT MATRIX A INTO COMPONENTS T AND F (BORROWED FROM FUNCTION splitA in library(popbio))
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
##
## additional FUNCTION TO ESTIMATE LAMBDA
Lambda <- function(At){
  ev <- eigen(At)
  lmax <- which.max(Re(ev$values))
  Re(ev$values[lmax])
}
#
##
####--------------
## FOR IUCN ASSESSMENT
## GET PRE-DECLINE GENERATION TIME - YEARS 1 TO 8
TFA <- TFfromA(matrix(rowMeans(sapply(1:8,getA,j=1)),ncol=9))
gen_time(TFA$T,TFA$F,method="cohort")
Lambda(matrix(rowMeans(sapply(1:8,getA,j=1)),ncol=9))
(LT <- mpm_to_table(TFA$T,TFA$F,xmax=25,remove_final=F))
##
## SET MATRICES FOR SIMULATIONS OF LIFE TABLES AND (COHORT) GENERATION TIME
gentCh <- numeric(15000)
LTarray <- array(NA,dim=c(25,6,15000))
meanMat <- matrix(NA,81,15000)
sims <- seq(1,15000,by=1) # seq(1,15000,by=15)
for(m in 1:15000){
  TFA <- TFfromA(matrix(rowMeans(sapply(1:8,getA,j=sims[m])),ncol=9))
  meanMat[,m] <- rowMeans(sapply(1:8,getA,j=sims[m])) 
  gentCh[m] <- gen_time(TFA$T,TFA$F,method="cohort")
  bLT <- mpm_to_table(TFA$T,TFA$F,xmax=25,lx_crit=0.001,remove_final=T)
  LTarray[1:dim(bLT)[1],,m] <- as.matrix(cbind(bLT[,c("x","px","lx","mx","lxmx")],xlxmx=bLT[,c("x")]*bLT[,c("lxmx")]))
  print(m)
}
##
## TO SUMMARISE LIFE TABLE RESULTS
## MEAN AND SD OF LIFE TABLE
mLT <- apply(LTarray,c(1,2),mean,na.rm=TRUE)
dimnames(mLT)[[2]] <- c("x","px","lx","mx","lxmx","xlxmx")
sLT <- apply(LTarray,c(1,2),sd,na.rm=TRUE)
dimnames(sLT)[[2]] <- c("x","px","lx","mx","lxmx","xlxmx")
## MEAN GENERATION TIME
round(c(mean=mean(gentCh),SD=sd(gentCh)),3)
#  mean    SD 
# 8.213 0.270 
###
#  MEAN MATRIX MODEL
(meanMatMod <- matrix(rowMeans(meanMat),ncol=9))
#  SD of MEAN MATRIX MODEL
(sdMat <- matrix(apply(meanMat,1,sd),ncol=9))
round(sdMat,3)
##
## Save data
saveRDS(gentCh,"processed_data/gentChF.Rds")
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
# 8.636 0.260
## Save data
saveRDS(gentCHw,"processed_data/gentCHw.Rds")
##



#paste("|",round(mLT[i,1]),"  |",round(mLT[i,2],4),"(",round(sLT[i,2],4),")","|",
#round(mLT[i,3],4),"(",round(sLT[i,3],4),")","|",round(mLT[i,4],4),
#"(",round(sLT[i,4],4),")","|",round(mLT[i,5],4),"(",round(sLT[i,5],4),")","|",
#round(mLT[i,6],4),"(",round(sLT[i,6],4),")","|");i=1
