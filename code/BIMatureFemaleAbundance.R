#####################################################################################
## ESTIMATES ABUNDANCE AND PRODUCTIVITY FOR A NUMBER OF YEARS WITH AVAILABLE COUNTS  
## Requires:
##  1. SSB female abundance and demographic estimates for 1984 to 2025.
##    -> Script: 'SSBFemale_IPM_1984-2025.r'
##    -> Data:   'AFSFemaleSimulations19842025.Rds'
##  2. Estimates of predicted female counts at SSB around mean counting dates for 
##     each surveyed season. 
##    -> Script: 'SSBFemalePMCountsForTrendAnalysis.r'
##    -> Data:   'SSBPredictedPMCounty1y2.Rds' (9 files)
##  3. Productivity and mature female estimate at Bird Island for 2008-09 Survey
##    ->  Data:   'BIMatureFemales0809.Rds'
##    ->  Data:   'BIProductivity0809.Rds'
##---------------------------------------------------------------------------------
library(mgcv)
library(jagsUI)
##---------------------------------------------------------------------------------
## READ DATA OBJECTS
## BI - SSB ESTIMATES 1984-2025
pred8425 <- readRDS("processed_data/AFSFemaleSimulations19842025.Rds")
## POSTERIOR SIMULATIONS OF ESTIMATED MATURE FEMALES AND PRODUCTIVITY FOR 2008-09
BIMatureFemales0809 <- readRDS("BIMatureFemales0809.Rds")
BIProductivity0809 <- readRDS("BIProductivity0809.Rds")
##---------------------------------------------------------------------------------
## TO GENERATE POSTERIOR SIMULATIONS OF PREDICTED SSB COUNTS
##---------------------------------------------------------------------------------
## GET AFTERNOON FEMALE COUNT DATA FROM SSB
##  - TOTAL SSB FEMALE COUNTS - 1989-2024 -
##  1. model observed daily counts n_j (day j=1,...,60) with Poisson GAM(M)s
##  2. predict expected SSB count for survey days at and around peak
##  GET DATA
tfemc <- read.csv("SSB_total_female_counts.csv", header = T)
day <- 1:60; tfemc <- tfemc[day, ]  ## Adding sequential day starting on November 1st
##---------------------------------------------------------------------------------
##  PREDICT PEAK PUPPING DAYS
##  (change file path as required)
newpups <- read.csv("newpups.csv", header = F)
colnames(newpups) <- 1985:2024
sumpup <- apply(newpups, 2, cumsum)
## FUNCTION ffs.obs FITS AN NLME (LOGISTIC) MODEL TO OBTAIN THE PEAK PUPPING DAY
source("ffs.obs.R")
## Survey years
syears <- c(1989, 1990, 1991, 1995, 1998, 2004, 2009, 2017, 2024)
pks <- numeric(length(syears))
for(i in 1:length(syears))
  pks[i] <- as.vector(ffs.obs(data.frame(day = 1:70, pups = sumpup[, paste(syears[i])]), plots = F)$coefficients[2])
## PEAK PUPPING DATES ARE:
## 1989 1990 1991 1995 1998 2004 2009 2017
## 37.6 34.3 39.4 43.4 44.0 37.4 36.3 43.4
## FOR 2007
as.vector(ffs.obs(data.frame(day = 1:70, pups = sumpup[, paste(2007)]), plots = F)$coefficients[2])
# 37.90519
## BASED ON MEAN DATES WHEN TOTAL COUNTS WERE MADE, BEST GAM PREDICTION DAYS ARE
predday <- c(38.5, 42.5, 44, 40, 46, 42.5, 36, 44, 38)
counts <- c(28796, 34286, 12854, 18622, 13579, 29733, 35303, 13955, 9620)
##---------------------------------------------------------------------------------
##  ANALYSE SURVEYED YEARS 
##  1988-89
##  BAYESIAN GAMM RESULTS
ffB <- readRDS("SSBPredictedPMCount8889.Rds")[, 3]
ffB <- rowMeans(readRDS("processed_data/SSBPredictedPMCount8889.Rds")[, 3:4])
round(c(mean(ffB), quantile(ffB, c(0.025,0.975))), 1)
#        2.5% 97.5% 
# 432.1 404.6 460.3
est <- pred8425[, paste("NBI[", which(1984:2024 %in% 1989), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#        CV  2.5% 97.5% 
# 747     3   697   798
# 759     3   708   812 
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 1.730 1.574 1.894
# 1.758 1.598 1.932 
mEst1 <- counts[1] * est / ffB
round(c(mean(mEst1), CV = 100 * sd(mEst1) / mean(mEst1), quantile(mEst1, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 49828     5 45313 54539
# 50625     5 46010 55631
est <- pred8425[, paste("Ntot[", which(1984:2024%in% 1989), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#         CV  2.5% 97.5% 
# 1078     5   987  1178
#  985     5   901  1077 
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 2.498 2.237 2.783
# 2.282 2.039 2.545
mEstN1 <- counts[1] * est / ffB
round(c(mean(mEstN1), CV = 100 * sd(mEstN1) / mean(mEstN1), quantile(mEstN1, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 71938     6 64422 80125
# 65723     6 58717 73274
##----------------------------------------------------------------------------------#
##  1989-90
##  BAYESIAN GAMM RESULTS
# ffB <- readRDS("processed_data/SSBPredictedPMCount8990.Rds")[, 2]
ffB <- rowMeans(readRDS("processed_data/SSBPredictedPMCount8990.Rds")[, 2:3])
round(c(mean(ffB), quantile(ffB, c(0.025,0.975))), 1)
#        2.5% 97.5% 
# 494.3 464.6 524.9
# 494.3 464.6 524.9
est <- pred8425[, paste("NBI[", which(1984:2024 %in% 1990), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#        CV  2.5% 97.5% 
# 825     3   774   878
# 818     3   764   875 
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 1.671 1.527 1.823
# 1.657 1.509 1.810
mEst2 <- counts[2] * est / ffB
round(c(mean(mEst2), CV = 100 * sd(mEst2) / mean(mEst2), quantile(mEst2, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 57279     4 52343 62500
# 56816     5 51746 62053
est <- pred8425[, paste("Ntot[", which(1984:2024 %in% 1990), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#         CV  2.5% 97.5% 
# 1078     4   991  1171
# 1031     4   951  1117 
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 2.183 1.964 2.420
# 2.088 1.886 2.305 
mEstN2 <- counts[2] * est / ffB
round(c(mean(mEstN2), CV = 100 * sd(mEstN2) / mean(mEstN2), quantile(mEstN2, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 74852     5 67346 82956
# 71581     5 64674 79036
##----------------------------------------------------------------------------------#
##  1990-91
##  BAYESIAN GAMM RESULTS
ffB <- readRDS("processed_data/SSBPredictedPMCount9091.Rds")[, 4]
round(c(mean(ffB), quantile(ffB, c(0.025,0.975))), 1)
#        2.5% 97.5% 
# 173.7 154.2 192.1
# 173.7 154.2 192.1 
est <- pred8425[, paste("NBI[", which(1984:2024 %in% 1991), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#        CV  2.5% 97.5% 
# 561     4   517   607
# 551     4   507   595 
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 3.239 2.840 3.727
# 3.180 2.781 3.644
mEst3 <- counts[3] * est / ffB
round(c(mean(mEst3), CV = 100 * sd(mEst3) / mean(mEst3), quantile(mEst3, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 41639     7 36507 47901
# 40872     7 35753 46835
est <- pred8425[, paste("Ntot[", which(1984:2024 %in% 1991), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#        CV  2.5% 97.5% 
# 865     5   776   962
# 819     5   738   908 
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
# 4.999 4.302 5.837
# 4.733 4.081 5.510
mEstN3 <- counts[3] * est / ffB
round(c(mean(mEstN3), CV = 100 * sd(mEstN3) / mean(mEstN3), quantile(mEstN3, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 64255     8 55294 75026  
# 60838     8 52460 70820
##----------------------------------------------------------------------------------#
##  1994-95
##  BAYESIAN GAMM RESULTS
ffB <- readRDS("processed_data/SSBPredictedPMCount9495.Rds")[, 4]
round(c(mean(ffB), quantile(ffB, c(0.025,0.975))), 1)
#        2.5% 97.5% 
# 233.4 211.8 254.9
# 233.4 211.8 254.9 
est <- pred8425[, paste("NBI[", which(1984:2024 %in% 1995), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#        CV  2.5% 97.5% 
# 594     4   548   641
# 593     4   548   640
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 2.550 2.262 2.869
# 2.546 2.256 2.874
mEst4 <- counts[4] * est / ffB
round(c(mean(mEst4), CV = 100 * sd(mEst4) / mean(mEst4), quantile(mEst4, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 47485     6 42128 53434
# 47409     6 42014 53511
est <- pred8425[, paste("Ntot[", which(1984:2024 %in% 1995), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#        CV  2.5% 97.5% 
# 868     5   785   953
# 795     5   725   869 
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 3.726 3.258 4.253
# 3.413 2.994 3.879 
mEstN4 <- counts[4] * est / ffB
round(c(mean(mEstN4), CV = 100 * sd(mEstN4) / mean(mEstN4), quantile(mEstN4, c(0.025, 0.975))))
#         CV  2.5% 97.5% 
# 69394     7 60676 79194
# 63551     7 55761 72233
##----------------------------------------------------------------------------------#
##  1997-98
##  BAYESIAN GAMM RESULTS
ffB <- readRDS("processed_data/SSBPredictedPMCount9798.Rds")[, 4]
round(c(mean(ffB), quantile(ffB, c(0.025,0.975))), 1)
#        2.5% 97.5% 
# 178.0 155.3 201.4
# 178.0 155.3 201.4
est <- pred8425[, paste("NBI[", which(1984:2024 %in% 1998), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#        CV  2.5% 97.5% 
# 548     4   506   591
# 555     4   513   600
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 3.092 2.655 3.597
# 3.134 2.695 3.647 
mEst5 <- counts[5] * est / ffB
round(c(mean(mEst5), CV = 100 * sd(mEst5) / mean(mEst5), quantile(mEst5, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 41992     8 36056 48837
# 42563     8 36599 49525
est <- pred8425[, paste("Ntot[", which(1984:2024 %in% 1998), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#         CV  2.5% 97.5% 
#  913     5   823  1009
#  848     5   768   938 
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 5.153 4.363 6.068
# 4.786 4.059 5.637
mEstN5 <- counts[5] * est / ffB
round(c(mean(mEstN5), CV = 100 * sd(mEstN5) / mean(mEstN5), quantile(mEstN5, c(0.025, 0.975))))
#         CV  2.5% 97.5% 
# 69976     8 59252 82394
# 64984     8 55121 76541
##----------------------------------------------------------------------------------#
##  2003-04
##  BAYESIAN GAMM RESULTS
ffB <- rowMeans(readRDS("processed_data/SSBPredictedPMCount0304.Rds")[, 3:4])
round(c(mean(ffB), quantile(ffB, c(0.025,0.975))), 1)
#        2.5% 97.5% 
# 413.1 386.0 440.7
# 413.1 386.0 440.7 
est <- pred8425[, paste("NBI[", which(1984:2024 %in% 2004), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#        CV  2.5% 97.5% 
# 764     3   715   816
# 765     4   714   819
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 1.853 1.686 2.032
# 1.854 1.684 2.039
mEst6 <- counts[6] * est / ffB
round(c(mean(mEst6), CV = 100 * sd(mEst6) / mean(mEst6), quantile(mEst6, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 55085     5 50123 60424
# 55125     5 50057 60638
est <- pred8425[, paste("Ntot[", which(1984:2024 %in% 2004), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#         CV  2.5% 97.5% 
# 1063     4   976  1156
#  953     4   880  1034 
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 2.575 2.307 2.863
# 2.311 2.080 2.568
mEstN6 <- counts[6] * est / ffB
round(c(mean(mEstN6), CV = 100 * sd(mEstN6) / mean(mEstN6), quantile(mEstN6, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 76575     5 68589 85132
# 68707     5 61844 76341
##----------------------------------------------------------------------------------#
##  2008-09
## - SSB Integrated model estimates for 2001 to 2022
preds <- readRDS("processed_data/AFSFemaleSimulations200122Best.Rds")
round(c(mean(preds[, which(dimnames(preds)[[2]] == "mu[36]")]), 
  quantile(preds[, which(dimnames(preds)[[2]] == "mu[36]")], c(0.025, 0.975))))
#       2.5% 97.5% 
#  304   272   337 
#  304   272   337
round(c(mean(preds[, which(dimnames(preds)[[2]] == "NBI.nc[4]")]), 
        quantile(preds[, which(dimnames(preds)[[2]] == "NBI.nc[4]")], c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 2.577 2.276 2.918
# 2.577 2.276 2.918
est <- pred8425[, paste("NBI[", which(1984:2024 %in% 2009), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#        CV  2.5% 97.5% 
# 778     3   728   830
# 774     3   722   827
round(c(mean(preds[, which(dimnames(preds)[[2]] == "NMat.nc[4]")]), 
  quantile(preds[, which(dimnames(preds)[[2]] == "NMat.nc[4]")], c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 3.065 2.687 3.496
# 3.065 2.687 3.496
est <- pred8425[, paste("Ntot[", which(1984:2024 %in% 2009), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#         CV  2.5% 97.5% 
# 1027     4   943  1116
# 1005     4   924  1090
##  (ESTIMATES PROVIDED BY SCRIPT 'SouthGeorgiaAbundance200809.r'
mEst7 <- BIProductivity0809
round(c(mean(mEst7), CV = 100 * sd(mEst7) / mean(mEst7), quantile(mEst7, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 60547     8 51413 70947
mEstN7 <- BIMatureFemales0809
round(c(mean(mEstN7), CV = 100 * sd(mEstN7) / mean(mEstN7), quantile(mEstN7, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 72019     9 60913 84906
##----------------------------------------------------------------------------------#
##  2016-17
##  BAYESIAN GAMM RESULTS
ffB <- readRDS("processed_data/SSBPredictedPMCount1617.Rds")[, 4]
round(c(mean(ffB), quantile(ffB, c(0.025,0.975))), 1)
#        2.5% 97.5% 
# 183.0 160.1 208.1
# 183.0 160.1 208.1 
est <- pred8425[, paste("NBI[", which(1984:2024 %in% 2017), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#        CV  2.5% 97.5% 
# 364     5   331   399
# 365     5   330   402 
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 1.999 1.700 2.342
# 2.004 1.692 2.355
mEst8 <- counts[8] * est / ffB
round(c(mean(mEst8), CV = 100 * sd(mEst8) / mean(mEst8), quantile(mEst8, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 27902     8 23717 32684
# 27973     8 23606 32870
est <- pred8425[, paste("Ntot[", which(1984:2024 %in% 2017), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#        CV  2.5% 97.5% 
# 541     5   486   599
# 504     6   451   560
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 2.968 2.504 3.500
# 2.766 2.325 3.265 
mEstN8 <- counts[8] * est / ffB
round(c(mean(mEstN8), CV = 100 * sd(mEstN8) / mean(mEstN8), quantile(mEstN8, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 41422     9 34946 48844
# 38599     9 32440 45556
##----------------------------------------------------------------------------------#
##  2023-24
##  BAYESIAN GAMM RESULTS
# 36 37 38 39 40 41 42
ffB <- rowMeans(readRDS("processed_data/SSBPredictedPMCount2324.Rds")[,3:6])
round(c(mean(ffB), quantile(ffB, c(0.025,0.975))), 1)
#        2.5% 97.5% 
# 151.9 140.5 163.7 
#
est <- pred8425[, paste("NBI[", which(1984:2024 %in% 2024), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#        CV  2.5% 97.5% 
# 304     5   273   336
round(c(mean(est / ffB), quantile(est / ffB, c(0.025, 0.975))), 3)
#        2.5% 97.5% 
# 2.002 1.757 2.272 
mEst9 <- counts[9]*est/ffB
round(c(mean(mEst9), CV = 100 * sd(mEst9) / mean(mEst9), quantile(mEst9, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 19256     7 16901 21856 
est <- pred8425[, paste("Ntot[", which(1984:2024 %in% 2024), "]", sep = "")]
round(c(mean(est), CV = 100 * sd(est) / mean(est), quantile(est, c(0.025, 0.975))))
#         CV  2.5% 97.5% 
#  391     6   343   443
round(c(mean(est/ffB),quantile(est/ffB,c(0.025,0.975))),3)
#        2.5% 97.5% 
# 2.579 2.215 2.976 
mEstN9 <- counts[9] * est / ffB
round(c(mean(mEstN9), CV = 100 * sd(mEstN9) / mean(mEstN9), quantile(mEstN9, c(0.025, 0.975))))
#          CV  2.5% 97.5% 
# 24809     8 21304 28629 
##----------------------------------------------------------------------------------#
####
##  DATA FRAMES FOR PRODUCTIVITY AND ABUNDANCES AT BIRD ISLAND DURING THE SURVEY YEARS
productivity <- data.frame(year=syears,mean=colMeans(cbind(mEst1,mEst2,mEst3,mEst4,mEst5,mEst6,mEst7,mEst8,mEst9)),
    low = apply(cbind(mEst1,mEst2,mEst3,mEst4,mEst5,mEst6,mEst7,mEst8,mEst9), 2, quantile, 0.025),
    high = apply(cbind(mEst1,mEst2,mEst3,mEst4,mEst5,mEst6,mEst7,mEst8,mEst9), 2, quantile, 0.975))
##
abundances <- data.frame(year=syears,mean=colMeans(cbind(mEstN1,mEstN2,mEstN3,mEstN4,mEstN5,mEstN6,mEstN7,mEstN8,mEstN9)),
    low = apply(cbind(mEstN1,mEstN2,mEstN3,mEstN4,mEstN5,mEstN6,mEstN7,mEstN8,mEstN9),2,quantile,0.025),
    high = apply(cbind(mEstN1,mEstN2,mEstN3,mEstN4,mEstN5,mEstN6,mEstN7,mEstN8,mEstN9),2,quantile,0.975))
##
## SAVE FOR IUCN ASSESSMENT
saveRDS(mEstN6,"process_data/BI04f.Rds")
saveRDS(mEstN7,"process_data/BI09f.Rds")
saveRDS(mEstN9,"process_data/BI24f.Rds")
#
##
saveRDS(productivity, "process_data/productivityBI.Rds")
saveRDS(abundances, "process_data/abundancesBI.Rds")
##
## MEAN AND DECLINE IN PRODUCTIVITY FOR 1989 - 2024
c(mean(rowMeans(cbind(mEst1,mEst2,mEst3,mEst4,mEst5,mEst6,mEst7))),
  quantile(rowMeans(cbind(mEst1,mEst2,mEst3,mEst4,mEst5,mEst6,mEst7)),c(0.025, 0.975)))
#               2.5%    97.5% 
#  50565.20 48244.38 53032.51 
## 2017
c(mean(mEst8),quantile(mEst8,c(0.025,0.975)))
#              2.5%    97.5% 
# 27972.63 23606.11 32870.16 
## 2024
c(mean(mEst9),quantile(mEst9,c(0.025,0.975)))
# 19256.19 16901.05 21856.01 
#
## exponential decay
100 * c(mean(1-exp(log(mEst8/rowMeans(cbind(mEst1,mEst2,mEst3,mEst4,mEst5,mEst6,mEst7)))/8)),
  quantile(1-exp(log(mEst8/rowMeans(cbind(mEst1,mEst2,mEst3,mEst4,mEst5,mEst6,mEst7)))/8),c(0.025,0.975)))
#              2.5%    97.5% 
# 7.164997 5.160145 9.132037 
## is the annual percent decline
## TOTAL DECREASE
perDecP <- 100 * (rowMeans(cbind(mEst1, mEst2, mEst3, mEst4, mEst5, mEst6, mEst7)) - mEst8) / rowMeans(cbind(mEst1, mEst2, mEst3, mEst4, mEst5, mEst6, mEst7))
#              2.5%    97.5% 
# 44.76833 34.84092 53.39636
round(c(mean(perDecP), quantile(perDecP, c(0.025, 0.975))), 1)
#        2.5% 97.5% 
#  44.6  34.5  53.5 
## 
## DECLINE IN PRODUCTIVITY 2009 -> 2024
## exponential decay
100 * c(mean(1-exp(log(mEst9/rowMeans(cbind(mEst1,mEst2,mEst3,mEst4,mEst5,mEst6,mEst7)))/8)),
        quantile(1-exp(log(mEst9/rowMeans(cbind(mEst1,mEst2,mEst3,mEst4,mEst5,mEst6,mEst7)))/8),c(0.025,0.975)))
#                2.5%     97.5% 
# 11.385468  9.834385 12.905300 
## is the annual percent decline
## TOTAL DECREASE
perDecP <- 100 * (rowMeans(cbind(mEst1, mEst2, mEst3, mEst4, mEst5, mEst6, mEst7)) - mEst9) / rowMeans(cbind(mEst1, mEst2, mEst3, mEst4, mEst5, mEst6, mEst7))
round(c(mean(perDecP), quantile(perDecP, c(0.025, 0.975))), 1)
#        2.5% 97.5% 
#  61.9  56.3  66.9 
##
## DECLINE IN PRODUCTIVITY 2017 -> 2024
## exponential decay
100 * c(mean(1-exp(log(mEst9/mEst8)/7)),quantile(1-exp(log(mEst9/mEst8)/7),c(0.025,0.975)))
#              2.5%    97.5% 
# 5.165437 2.325238 7.950666 
## is the annual percent decline
## TOTAL DECREASE 2017 -> 2024
perDecP <- 100 * (mEst8 - mEst9) / mEst8
round(c(mean(perDecP), quantile(perDecP, c(0.025, 0.975))), 1)
#        2.5% 97.5% 
#  30.7  15.2  44.0 
################################################################################


