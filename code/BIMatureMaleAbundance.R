####################################################################################
## ESTIMATES ABUNDANCE AND PRODUCTIVITY FOR A NUMBER OF YEARS WITH AVAILABLE COUNTS 
##                                                                                  
## Requires:
##  1. SSB male abundance and demographic estimates for 1984 to 2025.
##    -> Script: 'SSBMale_IPM_1994-2025.R'
##    -> Data:   'AFSMaleSimulations19952025.Rds'
##  2. Estimates of predicted male counts at SSB around mean counting dates for 
##     each surveyed season. 
##    -> Script: 'SSBMaleAMCountsForTrendAnalysis.r'
##    -> Data:   'SSBPredictedAMCounty1y2.Rds' (12 files)
##------------------------------------------------------------------------------
library(mgcv)
library(jagsUI)
##------------------------------------------------------------------------------
## READ DATA OBJECTS
## BI - SSB ESTIMATES 1995-2025
pred9525 <- readRDS("processed_data/AFSMaleSimulations19952025.Rds")
## TERRITORIALS - ASSUMES RECAPTURE p IS 0.9999999
NTs <- c(127,143,114,74,80,140,160,160,112,144,124,114,111,116,118,61,105,87,80,40,53,48,52,61,56,75,56,60,76,74,49)
NTs <- matrix(rep(NTs,each=15000),nrow=15000)
NALL <- NULL
for(i in 1:31)
  NALL <- cbind(NALL,pred9525[,paste("Ntot[",i,"]",sep="")])
## POSTERIOR SIMULATIONS OF ESTIMATED MATURE MALES AND PRODUCTIVITY FOR 2008-09
BIMatureMales0809 <- readRDS("processed_data/BIMatureMales0809.Rds")
BITerritorialMales0809 <- readRDS("processed_data/BITerritorialMales0809.Rds")
#
##------------------------------------------------------------------------------
##  GET AM MALE COUNT DATA FROM SSB
##   - TOTAL SSB MALE COUNTS -> 1992-2025 -
##------------------------------------------------------------------------------
tmc <- read.csv("raw_data/SSB_total_male_counts.csv", header = T)
day <- 1:61; tmc <- tmc[day,]
##
##   - SSB PUP PRODUCTIONS -> 1985-2025 -
##------------------------------------------------------------------------------
newpups <- read.csv("raw_data/newpups.csv",header=F)
colnames(newpups) <- 1985:2025
sumpup <- apply(newpups,2,cumsum)
## FUNCTION ffs.obs FITS AN NLME (LOGISTIC) MODEL TO OBTAIN THE PEAK PUPPING DAY
source("CODE/ffs.obs.R")
##
syears <- c(2004,2009,2017,2024)
pks <- ses <- numeric(length(syears))
for(i in 1:length(syears)){
  pks[i] <- as.vector(ffs.obs(data.frame(day=1:70,pups=sumpup[,paste(syears[i])]),plots=F)$coefficients[2])
  ses[i] <- sqrt(ffs.obs(data.frame(day=1:70,pups=sumpup[,paste(syears[i])]),plots=F)$varBeta[2, 2])
}  
## PEAK PUPPING DATES ARE:
round(pks, 1)
# [1] 37.4 36.3 43.4 38.4
predday <- round(pks)
##-------------------------------------------------------------------------------------------------------------------------#
## SURVEY COUNTS
years <- c(2004,2009,2017,2024)
counts <-c(7751,7470,3395,2703)
##-------------------------------------------------------------------------------------------------------------------------#
##  ANALYSE SURVEYED YEARS
##  2003-04
##  BAYESIAN GAMM RESULTS
ffB <- rowMeans(readRDS("processed_data/SSBPredictedAMCount0304.Rds")[, 3:4])
round(c(mean(ffB), quantile(ffB, c(0.025,0.975))), 1)
#       2.5% 97.5% 
# 56.1  51.1  61.3
## TERRITORIALS
est <- NTs[,which(1995:2025 %in% 2004)]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#         CV  2.5% 97.5% 
# 144     0   144   144 
round(c(mean(est/ffB),quantile(est/ffB,c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 2.573 2.350 2.819
mEst1 <- counts[1]*est/ffB
round(c(mean(mEst1), CV=100*sd(mEst1)/mean(mEst1),quantile(mEst1,c(0.025,0.975))))
#          CV  2.5% 97.5% 
# 19945     5 18216 21847
## ALL MATURE
est <- NALL[,which(1995:2025 %in% 2004)]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
#   164    13   144   218 
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 2.925 2.427 3.965 
mEstN1 <- counts[1]*est/ffB
round(c(mean(mEstN1),CV=100*sd(mEstN1)/mean(mEstN1),quantile(mEstN1,c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 22668    14 18809 30731
##----------------------------------------------------------------------------------#
##  2008-09
## - SSB Integrated model estimates for 1995 to 2025
##  (ESTIMATES PROVIDED BY SCRIPT 'SouthGeorgiaAbundance200809.r'
BIMatureMales0809 <- readRDS("processed_data/BIMatureMales0809.Rds")
BITerritorialMales0809
##
mEst2 <- BITerritorialMales0809
round(c(mean(mEst2), CV = 100 * sd(mEst2)/mean(mEst2),quantile(mEst2,c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 21629    10 18281 26254 
mEstN2 <- BIMatureMales0809
round(c(mean(mEstN2),CV=100*sd(mEstN2)/mean(mEstN2),quantile(mEstN2,c(0.025,0.975))))
#          CV  2.5% 97.5% 
# 23861     9 20045 28939
##----------------------------------------------------------------------------------#
##  2016-17
##  BAYESIAN GAMM RESULTS
ffB <- readRDS("processed_data/SSBPredictedAMCount1617.Rds")[,4]
round(c(mean(ffB), quantile(ffB,c(0.025,0.975))), 1)
#       2.5% 97.5% 
# 24.8  21.7  28.2 
est <- NTs[,which(1995:2025 %in% 2017)]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#         CV  2.5% 97.5% 
#   52     0    52    52 
round(c(mean(est/ffB), quantile(est/ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 2.102 1.844 2.394
mEst3 <- counts[3]*est/ffB
round(c(mean(mEst3),CV=100*sd(mEst3)/mean(mEst3),quantile(mEst3,c(0.025,0.975))))
#         CV  2.5% 97.5% 
# 7137     7  6261  8129 
est <- NALL[,which(1995:2025 %in% 2017)]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 60    16    52    85 
round(c(mean(est/ffB),quantile(est/ffB,c(0.025,0.975))),3)
#        2.5% 97.5% 
# 2.434 1.919 3.525 
mEstN3 <- counts[3]*est/ffB
round(c(mean(mEstN3),CV=100*sd(mEstN3)/mean(mEstN3),quantile(mEstN3,c(0.025,0.975))))
#round(c(mean(mEstN3),CV=100*sd(mEstN3)/mean(mEstN3),quantile(mEstN3,c(0.025,0.975))))
#          CV  2.5% 97.5% 
#  8263    17  6515 11967 
##----------------------------------------------------------------------------------#
##  2023-24
##  BAYESIAN GAMM RESULTS
# 36 37 38 39 40 41 42
ffB <- rowMeans(readRDS("processed_data/SSBPredictedAMCount2324.Rds")[,3:6])
round(c(mean(ffB), quantile(ffB, c(0.025,0.975))), 1)
#       2.5% 97.5% 
# 34.7  31.2  38.5  
est <- NTs[,which(1995:2025 %in% 2024)]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#         CV  2.5% 97.5% 
#   74     0    74    74 
round(c(mean(est/ffB),quantile(est/ffB,c(0.025,0.975))),3)
#        2.5% 97.5% 
# 2.141 1.921 2.374 
mEst4 <- counts[4]*est/ffB
round(c(mean(mEst4), CV = 100 * sd(mEst4) / mean(mEst4), quantile(mEst4, c(0.025, 0.975))))
#         CV  2.5% 97.5% 
# 5786     5  5193  6418
est <- NALL[,which(1995:2025 %in% 2024)]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#       CV  2.5% 97.5% 
# 86    17    74   125 
round(c(mean(est/ffB),quantile(est/ffB,c(0.025,0.975))),3)
#        2.5% 97.5% 
# 2.501 1.991 3.650
mEstN4 <- counts[4] * est / ffB
round(c(mean(mEstN4), CV = 100 * sd(mEstN4) / mean(mEstN4), quantile(mEstN4, c(0.025, 0.975))))
#         CV  2.5% 97.5% 
# 6760    18  5380  9865 
##----------------------------------------------------------------------------------#
####
##  DATA FRAMES FOR ABUNDANCES OF TERRITORIALS AND MATURE MALES AT BIRD ISLAND DURING THE SURVEY YEARS
territorials <- data.frame(year=syears,mean=colMeans(cbind(mEst1,mEst2,mEst3,mEst4)),
    low = apply(cbind(mEst1,mEst2,mEst3,mEst4),2,quantile,0.025),
    high = apply(cbind(mEst1,mEst2,mEst3,mEst4),2,quantile,0.975))
##
abundances <- data.frame(year=syears,mean=colMeans(cbind(mEstN1,mEstN2,mEstN3,mEstN4)),
    low = apply(cbind(mEstN1,mEstN2,mEstN3,mEstN4),2,quantile,0.025),
    high = apply(cbind(mEstN1,mEstN2,mEstN3,mEstN4),2,quantile,0.975))
##
##
## SAVE FOR IUCN ASSESSMENT
saveRDS(mEstN1,"processed_data/BI04m.Rds")
saveRDS(mEstN2,"processed_data/BI09m.Rds")
saveRDS(mEstN4,"processed_data/BI24m.Rds")
##
saveRDS(territorials,"processed_data/territorialsBI.Rds")
saveRDS(abundances, "processed_data/abundancesMBI.Rds")
##
## MEAN AND DECLINE IN TERRITORIALS FOR 1995 - 2024
c(mean(rowMeans(cbind(mEst1,mEst2))),
  quantile(rowMeans(cbind(mEst1,mEst2)),c(0.025, 0.975)))
#              2.5%    97.5% 
# 20786.75 18874.13 23230.11 
#
## 2017
c(mean(mEst3),quantile(mEst3,c(0.025,0.975)))
#              2.5%    97.5% 
# 7136.571 6260.899 8128.741 
#
## 2024
c(mean(mEst4),quantile(mEst4,c(0.025,0.975)))
#              2.5%    97.5% 
# 5786.081 5192.853 6417.858
####################
## DECLINE IN TERRITORIALS 2009 -> 2017
## exponential decay
100*c(mean(1-exp(log(mEst3/rowMeans(cbind(mEst1,mEst2)))/8)),
        quantile(1-exp(log(mEst3/rowMeans(cbind(mEst1,mEst2)))/8),c(0.025,0.975)))
#              2.5%    97.5% 
# 12.51246 10.68584 14.38083
## is the annual percent decline
## TOTAL DECREASE
perDecP <- 100 * (rowMeans(cbind(mEst1,mEst2))-mEst3)/rowMeans(cbind(mEst1,mEst2))
round(c(mean(perDecP), quantile(perDecP, c(0.025, 0.975))), 1)
#       2.5% 97.5% 
# 65.6  59.5  71.1 
##
## DECLINE OF TERRITORIALS 2017 -> 2024
## exponential decay
100*c(mean(1-exp(log(mEst4/mEst3)/7)),quantile(1-exp(log(mEst4/mEst3)/7),c(0.025,0.975)))
#               2.5%     97.5% 
# 2.9346702 0.5888945 5.2917786 
## is the annual percent decline
## TOTAL DECREASE 2017 -> 2024
perDecP <- 100*(mEst3-mEst4)/mEst3
round(c(mean(perDecP), quantile(perDecP, c(0.025, 0.975))), 1)
#       2.5% 97.5% 
# 18.6   4.1  31.7 
####################
## DECLINE IN MATURE 2009 -> 2024
## exponential decay
100*c(mean(1-exp(log(mEstN4/rowMeans(cbind(mEstN1,mEstN2)))/15)),
      quantile(1-exp(log(mEstN4/rowMeans(cbind(mEstN1,mEstN2)))/15),c(0.025,0.975)))
#              2.5%    97.5% 
# 7.967022 5.359928 9.787444
## is the annual percent decline
## TOTAL DECREASE
perDecP <- 100 * (rowMeans(cbind(mEstN1,mEstN2))-mEstN4)/rowMeans(cbind(mEstN1,mEstN2))
round(c(mean(perDecP), quantile(perDecP, c(0.025, 0.975))), 1)
#       2.5% 97.5% 
# 70.8  56.2  78.7   
##
## DECLINE INTERRITORIALS 2017 -> 2024
## exponential decay
100*c(mean(1-exp(log(mEstN3/mEstN4)/7)),quantile(1-exp(log(mEstN3/mEstN4)/7),c(0.025,0.975)))
#                2.5%     97.5% 
# -2.972274 -9.923573  3.835422 
## is the annual percent decline
## TOTAL DECREASE 2017 -> 2024
perDecP <- 100*(mEstN3-mEstN4)/mEst3
round(c(mean(perDecP), quantile(perDecP, c(0.025, 0.975))), 1)
#        2.5% 97.5% 
#  20.6 -32.7  76.0 
################################################################################
