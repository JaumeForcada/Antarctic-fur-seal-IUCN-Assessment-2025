#####################################################################################
## ASSESSMENT OF CRITERION A (A2) - ANNUAL CHANGE AND POPULATION REDUCTION
#----------------
## 1. SOURCE SAVED ABUNDANCE SIMULATIONS FOR SOUTH GEORGIA IN 2007/9 AND 2022
##   . GENERATED WITH SCRIPT 'SouthGeorgiaPopulationAssessment.R'
## 2. OBTAIN ABUNDANCE SIMULATIONS FOR OTHER SUBPOPULATIONS
## 3. OBTAIN REDUCTION ESTIMATES UNDER CRITERION A1/A2, BASED ON ESTIMATES OF POPULATION 
##    SIZE FROM TWO DIFFERENT YEARS, ASSUMING A THREE GENERATION DECLINE WITH EXPONENTIAL 
##    PROCESS (CONSTANT RATE OF DECLINE) FOR EACH SURVEYED COLONY
##   -> 
##------------------------------------------------------------------------------
##
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
##  -> BOTH REQUIRE MAIN ARGUMENT yp, WHICH IS A MATRIX WITH 4 COLUMNS X NUMBER OF SIMULATIONS
##     e.g. cbind(rep(2007,15000),MV07,rep(2022,15000),MV22)
## 
## ANNUAL CHANGE
annualChange <- function(i,yp)
  as.vector((yp[i,4]/yp[i,2])^(1/(yp[i,3]-yp[i,1])))
##
## REDUCTION
reduction <- function(i,yp, G=gentCHw,ay=2025,short=T,lowc=F,lb=eLogLs[1:10000]){
  ap <- min(c(max(c(10,G[i]*3)),100))
  y3a <- ay - ap
  ac <- (yp[i,4]/yp[i,2])^(1/(yp[i,3]-yp[i,1]))
  yby3y1 <- yp[i,1] - y3a
  yby2ay <- ay - yp[i,3]
  if(lowc){
    if(length(lb)>1){
      cb3y1 <- lb[i]^yby3y1
    }else{
      cb3y1 <- lb^yby3y1
    }
  }else{
    cb3y1 <- ac^yby3y1  
  }
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
##
## ARGUMENTS INCLUDE
# i -> vector index over which to apply the function
# yp -> data matrix including simulations for year1, abundance 1, year 2 and abundance 2
# G -> generation length simulations, normally G=gentCHw
# ay -> assessment year, by default 2025
# short -> logic 
#   . T = vector of simulations of reduction
#   . F = matrix of with simulations of yp, yby3y1, yby2ay, ac,100*(1-ac), cb3y1, cby2ay, P3, PC, C3G 
# lowc -> logic
#   . T = assumes low rate of decline for years between 3 generations ago and year1
#   . F = assumes same rate of decline for years between 3 generations ago and year1
# lb -> is the "lambda before" or population growth rate between Year 3G and year t_1, and defaults to 
#    eLogLs[1:10000], which is an object with simulations of the average Log-lambda for South Georgia 
#    between 1995 and 2006, obtained with script 'SSBMatureFemale_IPM_1984-2025.R'.
## ALL OUTCOMES ASSUME THAT THE RATE OF DECLINE CONTINUES TO BE THE SAME FOR years between year 2 and assessment year
####
##  SOURCE AVERAGE GENERATION TIME
gentCHw <- readRDS("processed_data/gentCHw.Rds")
##  SOURCE ANNUAL CHANGE BEFORE YEAR 1 IF YEAR 1 > YEAR 3G AGO
eLogLs <- readRDS("processed_data/earlyLogLambdas.Rds")
##  SOURCE MULTIPLIER TO CONVERT PUP PRODUCTION TO MATURE POPULATION
PtoM <- readRDS("processed_data/PtoM.Rds")
## FUNCTION TO GENERATE LOG-NORMAL DEVIATES FROM A MEAN AND VARIANCE
##  -> Converts standard normal random values to lognormals with the defined 
##     means ('smean') and variances ('SVAR').
lns <- function(ndevs, smean, svar){
  nmean = log(smean) - 0.5*log(svar/smean^2 + 1)
  nvar = log(svar/smean^2 + 1)
  normals = rnorm(ndevs) * sqrt(nvar) + nmean
  exp(normals)
}
#
## SUMMARY RESULTS FUNCTION
redSum <- function(x)
  c(N99=round(mean(x[11,])),sdN99=round(sd(x[11,])),
    t1=unique(x[1,]),N1=round(mean(x[2,])),sdN1=round(sd(x[2,])),
    tT=unique(x[3,]),NT=round(mean(x[4,])),sdNT=round(sd(x[4,])),
    N25=round(mean(x[12,])),sdN25=round(sd(x[12,])),
    AC=round(mean(x[8,]),3),sdAC=round(sd(x[8,]),3),
    R=round(mean(x[13,]),3),sdR=round(sd(x[13,]),3))

#A2Table <- data.frame(N99=double(),sdN99=double(),t1=integer(),N1=double(),sdN1=double(),tT=integer(),NT=double(),sdNT=double(),N25=double(),
#  sdN25=double(),AC=double(),sdAC=double(),R=double(),sdR=double(),A2=character())

##------------------------------------------------------------------------------
##  CALCULATE ANNUAL CHANGE AND REDUCTIONS FOR ALL SUBPOPULATIONS
# ------------------- 
##  1. SOUTH GEORGIA
##  -> SOURCE DATA (FROM SCRIPT 'SouthGeorgiaPopulationAssessment.R')
##  .  subsample to obtain only 10000 simulations
SG09 <- readRDS("processed_data/SG09.Rds")[1:10000]
SG22 <- readRDS("processed_data/SG22.Rds")
##  -> AC and REDUCTION
ACSG <- sapply(1:10000,annualChange,yp=cbind(rep(2009,10000),SG09,rep(2022,10000),SG22[1:10000]))
redSG <- sapply(1:10000,reduction,yp=cbind(rep(2009,10000),SG09,rep(2022,10000),SG22[1:10000]),short=F,lowc=T)
##
redSum(redSG)
#     N99  sdN99   t1      N1   sdN1    tT      NT  sdNT     N25   sdN25     AC   sdAC       R    sdR
# 2136953 153084 2009 1721209  71702  2022 1024797 83301  910114   89753  3.926  0.614  -0.572  0.051
##
#  PROBABILITY OF REDUCTION BEING P <= -0.5
100*length(which(redSG[13,]<=-0.5))/length(redSG[13,])
# [1] 91.61
round(c(mean(redSG[13,]),quantile(redSG[13,],c(0.025,0.975))),3) # -0.567 [-0.648 -0.472] -> ENDANGERED
#          2.5%  97.5% 
# -0.572 -0.663 -0.464 
##
#--------------------
##  2. SOUTH SANDWICH ISLANDS
##  1 -> ONLY ZAVODOVSKI SURVEYED IN 2022 - OTHER ISLANDS ASSUMED STABLE
##  . 1997/98 -> 1250 + 500  = 1750
##  . 2023/24 -> 1250 + 1265 = 2515
##  -> AC and REDUCTION
ACSSD <- sapply(1:10000,annualChange,yp=cbind(rep(1998,10000),1750*PtoM[1:10000],rep(2024,10000),2515*PtoM[1:10000]))
redSSD <- sapply(1:10000,reduction,yp=cbind(rep(1998,10000),1750*PtoM[1:10000],rep(2024,10000),2515*PtoM[1:10000]),short=F,lowc=F)
round(c(mean(redSSD[13,]),quantile(redSSD[13,],c(0.025,0.975))),3)  # #0.490 [0.442 0.532] -> LEAST CONCERN
#        2.5% 97.5% 
# 0.435 0.408 0.468
redSum(redSSD)
#  N99  sdN99     t1   N1  sdN1    tT    NT sdNT  N25 sdN25      AC  sdAC     R cc sdR 
# 2535     93   1998 2595    85  2024  3730  122 3782   124  -1.405 0.000 0.493  0.027
##
##--------------------
#  SOUTH ORKNEYS ARE ASSUME TO HAVE NO PUPS, BASED ON RECENT SIGNY IS AND LAURIE IS. RECORDS
##--------------------
#  3.  SOUTH SHETLANDS PUP PRODUCTION
#  -> 1995/96 9530 (135)
#  -> 2001/02 10057 (142)
#  -> 2007/08 7602 (103) 
#  -> 2024/25 901 (2.38)
##
eL <- ((lns(10000,10057,20164)*PtoM[1:10000])/(lns(10000,9530,18106)*PtoM[1:10000]))^(1/6)
set.seed(666)   ## RANDOM NUMBER GENERATOR FOR lns (LOG-NORMAL) DEVIATES TO SAMPLE SD OF PUP COUNTS - REPRODUCIBILITY
ACSSI <- sapply(1:10000,annualChange,yp=cbind(rep(2002,10000),lns(10000,10057,20164)*PtoM[1:10000],rep(2025,10000),lns(10000,901,5.6644)*PtoM[1:10000]))
redSSI <- sapply(1:10000,reduction,yp=cbind(rep(2002,10000),lns(10000,10057,20164)*PtoM[1:10000],rep(2025,10000),lns(10000,901,5.6644)*PtoM[1:10000]),short=F,lowc=T,lb=eL)
round(c(mean(redSSI[13,]),quantile(redSSI[13,],c(0.025,0.975))),3)  # -0.915 [-0.919 -0.912] -> CRITICALLY ENDANGERED
#          2.5%  97.5% 
# -0.908 -0.911 -0.904 
redSum(redSSI)
#   N99  sdN99   t1    N1 sdN1   tT   NT sdNT   N25 sdN25    AC   sdAC       R    sdR 
# 14531    556 2002 14913  531 2025 1336   44  1336    44 9.957  0.056  -0.908  0.002 
#  PROBABILITY OF REDUCTION BEING P <= -0.5
100*length(which(redSSI[13,] < -0.8))/length(redSSI[13,])
# 100
#
##--------------------
#  3. BOUVETOYA PUP PRODUCTION
#  -> 1998/1999 PUP PRODUCTION 15448 (SD = 370)
#  -> 2000/2001 PUP PRODUCTION 15215 (SD = 353)
#  -> 2001/2002 PUP PRODUCTION 15523 (SD = 480)
#  -> 2007/2008 PUP PRODUCTION THOUGHT TO BE round(47000 / 4.26) = 11033 OR at an annual change of 5.6% between 2001 and 2006: round(15523*(1-5.6/100)^5)=11637
#  -> 2014/2015 - Lowther Pers. comm.- 8080 (SE: 310; 95%CI: 7495-8711) - 
##
set.seed(666)
ACBV <- sapply(1:10000,annualChange,yp=cbind(rep(1999,10000),lns(10000,15448,136900)*PtoM[1:10000],rep(2015,10000),lns(10000,8080,98077)*PtoM[1:10000]))
redBV <- sapply(1:10000,reduction,yp=cbind(rep(1999,10000),lns(10000,15448,136900)*PtoM[1:10000],rep(2015,10000),lns(10000,8080,98077)*PtoM[1:10000]),short=F,lowc=F)
round(c(mean(redBV[13,]),quantile(redBV[13,],c(0.025,0.975))),3)  # [-0.649 -0.699 -0.596]   -> ENDANGERED
#          2.5%  97.5% 
# -0.649 -0.701 -0.593 
redSum(redBV)
#   N99 sdN99    t1    N1   sdN1   tT    NT sdNT   N25 sdN25    AC    sdAC       R    sdR 
# 22832  1141  1999 22905   926  2015 11990  610  8009   582 3.967   0.273  -0.649  0.028 
#  PROBABILITY OF REDUCTION BEING P <= -0.5
100*length(which(redBV[13,]<=-0.5))/length(redBV[13,])
# [1] 100
##------------------------------------------------------------------------------
## EAST ANTARCTICA
##--------------------
#  4. MARION ISLAND + PRINCE EDWARD ISLAND
#  -> 1999/2000 PUP PRODUCTION 343
#  -> 2000/2001 PUP PRODUCTION 464
#  -> 2003/2004 PUP PRODUCTION 744 (5)
#  -> 2006/2007 PUP PRODUCTION 1105 (9)
#  -> 2009/2010 PUP PRODUCTION 1379 (27)
#  -> 2012/2013 PUP PRODUCTION 1553 (64)
#  -> 2015/2016 PUP PRODUCTION 2052 (107)
#  -> 2018/2019 PUP PRODUCTION 2434 (169)
#  -> 2021/2022 PUP PRODUCTION 2528 (134)
#  -> 2024/2025 PUP PRODUCTION 3172 (156)
##
set.seed(666)
ACMI <- sapply(1:10000,annualChange,yp=cbind(rep(2000,10000),343*PtoM[1:10000],rep(2025,10000),lns(10000,3172,24336)*PtoM[1:10000]))
redMI <- sapply(1:10000,reduction,yp=cbind(rep(2000,10000),343*PtoM[1:10000],rep(2025,10000),lns(10000,3172,24336)*PtoM[1:10000]),short=F,lowc=F)
round(c(mean(redMI[13,]),quantile(redMI[13,],c(0.025,0.975))),3)  # [6.232 5.484 7.033]   -> LEAST CONCERN
#        2.5% 97.5% 
# 9.054  7.529 10.946 
redSum(redMI)
# N99 sdN99   t1   N1 sdN1   tT   NT sdNT  N25 sdN25      AC   sdAC     R      sdR 
# 470    36 2000  509   17 2025 4704  278 4704   278  -9.300  0.215  9.054   0.873 
## ANNUAL CHANGE
round(c(100*(mean(redMI[7,])-1),SD=sd(redMI[7,]),100*(quantile(redMI[7,],c(0.025,0.975))-1)),3)
#          SD  2.5% 97.5% 
# 9.300 0.002 8.883 9.719  
##--------------------
## 4.1 PRINCE EDWARD ISLAND
#  -> 2000/2001 PUP PRODUCTION 404
set.seed(666)
PEI <- 404*PtoM[1:10000]
round(c(mean=mean(PEI),SD=sd(PEI),c(quantile(PEI,c(0.025,0.975)))))
# mean    SD  2.5% 97.5% 
#  599    20   563   639 
## LEAST CONCERN (STABLE?)
##
##--------------------
## 5. CROZET, (Possession Islans)
#  ->  1999/2000 PUP PRODUCTION 234
#  ->  2003/2004 PUP PRODUCTION 295
#  ->  2020/2021 PUP PRODUCTION 504
#  ->  2023/2024 PUP PRODUCTION 530
##
set.seed(666)
ACCZ <- sapply(1:10000,annualChange,yp=cbind(rep(2000,10000),234*PtoM[1:10000],rep(2024,10000),530*PtoM[1:10000]))
redCZ <- sapply(1:10000,reduction,yp=cbind(rep(2000,10000),234*PtoM[1:10000],rep(2024,10000),530*PtoM[1:10000]),short=F,lowc=F)
round(c(mean(redCZ[13,]),quantile(redCZ[13,],c(0.025,0.975))),3)  # [1.418 1.306 1.555]   -> LEAST CONCERN
redSum(redCZ)
# N99 sdN99    t1   N1  sdN1   tT   NT sdNT  N25 sdN25      AC   sdAC     R    sdR 
# 337    15  2000  347    11 2024  786   26  813    27  -3.465  0.000 1.418  0.064
## ANNUAL CHANGE
round(c(100*(mean(redCZ[7,])-1),SD=sd(redCZ[7,]),100*(quantile(redCZ[7,],c(0.025,0.975))-1)),3)
#           SD  2.5% 97.5% 
#  3.465 0.000 3.465 3.465 
##--------------------
## 6. KERGUELEN (Courbet Peninsula + Nuageuses Is)
##  -> ONLY POSSESSION HAS BEEN REGULARLY COUNTED - NUAGEUSES ASSUMED STABLE (C. Guinet, personal communication)
##  . 1999/00 -> 1600 + 3600  = 5200
##  . 2023/24 -> 3600 + 4793 = 8393
##
set.seed(666)
ACKG <- sapply(1:10000,annualChange,yp=cbind(rep(2000,10000),5200*PtoM[1:10000],rep(2024,10000),8393*PtoM[1:10000]))
redKG <- sapply(1:10000,reduction,yp=cbind(rep(2000,10000),5200*PtoM[1:10000],rep(2024,10000),8393*PtoM[1:10000]),short=F,lowc=T,lb=1.014)
round(c(mean(redKG[13,]),quantile(redKG[13,],c(0.025,0.975))),3)  # [0.668 0.636 0.706 ]   -> LEAST CONCERN
redSum(redKG)
#  N99   sdN99   t1    N1  sdN1    tT    NT sdNT   N25 sdN25     AC   sdAC     R    sdR 
# 7616     267 2000  7712   252  2024 12448  407 12698   415 -2.015  0.000 0.668  0.018 

## ANNUAL CHANGE
round(c(100*(mean(redKG[7,])-1),SD=sd(redKG[7,]),100*(quantile(redKG[7,],c(0.025,0.975))-1)),3)
#          SD  2.5% 97.5% 
# 2.015 0.000 2.015 2.015 
#
##--------------------
##  7. HEARD ISLAND + McDONALD
##  -> ONLY HEARD ISLAND WAS EVALUATED
##  . 2000/01 -> 1012
##  . 2003/04 -> 1278
set.seed(666)
redHI <- sapply(1:10000,reduction,yp=cbind(rep(2001,10000),1012*PtoM[1:10000],rep(2004,10000),1278*PtoM[1:10000]),short=F,lowc=F)
redSum(redHI)
#  N99  sdN99   t1   N1 sdN1   tT   NT sdNT
# 1296     90 2001 1501   49 2004 1895   62
## ANNUAL CHANGE
ACHI <- sapply(1:10000,annualChange,yp=cbind(rep(2001,10000),1012*PtoM[1:10000],rep(2004,10000),1278*PtoM[1:10000]))
round(c(100*(mean(ACHI)-1),SD=sd(ACHI),100*(quantile(ACHI,c(0.025,0.975))-1)),3)
#          SD  2.5% 97.5% 
# 8.089 0.000 8.089 8.089 
#
##--------------------------
##  8. Macquarie Island
#   -> 1999/00 PUP PRODUCTION 105
#   -> 2011/12 PUP PRODUCTION 189
##
set.seed(666)
ACMQ <- sapply(1:10000,annualChange,yp=cbind(rep(1999,10000),105*PtoM[1:10000],rep(2012,10000),189*PtoM[1:10000]))
redMQ <- sapply(1:10000,reduction,yp=cbind(rep(1999,10000),105*PtoM[1:10000],rep(2012,10000),189*PtoM[1:10000]),short=F,lowc=F)
round(c(mean(redMQ[13,]),quantile(redMQ[13,],c(0.025,0.975))),3)  # [0.576 0.574 0.578]   -> DATA DEFICIENT
redSum(redMQ)
# N99 sdN99   t1   N1 sdN1   tT   NT sdNT  N25 sdN25       AC     sdAC        R      sdR 
# 156     8 1999  156    5 2012  280    9  505    17   -4.625    0.000    2.229    0.113 
round(c(mean=100*(mean(ACMQ)-1),SD=sd(ACMQ),100*(quantile(ACMQ,c(0.025,0.975))-1)),3)
#   mean    SD  2.5% 97.5% 
#  4.625 0.000 4.625 4.625
##
###
###################################################################################################
## SUBPOPULATION AND GLOBAL ASSESSMENT
#
## P3 is population 3 generations ago
## PC is current population
## C3G is population change in 3 generations, or REDUCTION, for declining populations
#
#---------------------
## SOUTH GEORGIA
##--------------------
P3SG <- redSG[11,]
paste(round(mean(P3SG))," (",round(sd(P3SG)),")",sep="")
# [1] "2136953 (153084)"
PCSG <- redSG[12,]
paste(round(mean(PCSG))," (",round(sd(PCSG)),")",sep="")
# [1] "910114 (89753)"
C3GSG <- (PCSG-P3SG)/P3SG
c(mean=mean(C3GSG),sd=sd(C3GSG),quantile(C3GSG,c(0.05,0.95)))
#        mean          sd          5%         95% 
# -0.57207421  0.05096133 -0.64994013 -0.48351286 
#
## PROBABILITY OF REDUCTION BEING P <= -0.5
100*length(which(redSG[13,] <= -0.5))/length(redSG[13,])
# [1] 91.61
##  ENDANGERED
## ANNUAL CHANGE
acSG <- (PCSG/P3SG)^(1/(2025-1999))
round(c(mean=mean(-100*(1-acSG)),SD=sd(-100*(1-acSG)),quantile(-100*(1-acSG),c(0.025,0.975))),3)
#   mean     SD   2.5%  97.5% 
# -3.237  0.441 -4.094 -2.373
##------------------------
## SOUTH SHETLAND ISLANDS
##------------------------
P3SS <- redSSI[11,]
paste(round(mean(P3SS))," (",round(sd(P3SS)),")",sep="")
# [1] "14529 (551)"
PCSS <- redSSI[12,]
paste(round(mean(PCSS))," (",round(sd(PCSS)),")",sep="")
# [1] "1336 (44)"
C3GSS <- (PCSS-P3SS)/P3SS
c(mean=mean(C3GSS),sd=sd(C3GSS),quantile(C3GSS,c(0.025,0.975)))
#        mean           sd         2.5%        97.5% 
# -0.907990342  0.001722019 -0.911149751 -0.904432441 
#
## PROBABILITY OF REDUCTION BEING P <= -0.5
100*length(which(redSSI[13,] <= -0.5))/length(redSSI[13,])
# 100
## PROBABILITY OF REDUCTION BEING P <= -0.8
100*length(which(redSSI[13,] <= -0.8))/length(redSSI[13,])
# [1] 100
##  CRITICALLY ENDANGERED
acSS <- (PCSS/P3SS)^(1/(2025-1999))
round(c(mean=mean(-100*(1-acSS)),sd=sd(-100*(1-acSS)),quantile(-100*(1-acSS),c(0.025,0.975))),3)
#   mean     sd   2.5%  97.5% 
# -8.769  0.066 -8.890 -8.635
##------------------------
##  BOUVETØYA
##------------------------
P3BV <- redBV[11,]
paste(round(mean(P3BV))," (",round(sd(P3BV)),")",sep="")
# [1] "22832 (1141)"
PCBV <- redBV[12,]
paste(round(mean(PCBV))," (",round(sd(PCBV)),")",sep="")
# [1] "8009 (582)"
C3GBV <- (PCBV-P3BV)/P3BV
c(mean=mean(C3GBV),sd=sd(C3GBV),quantile(C3GBV,c(0.05,0.95)))
#       mean          sd          5%         95% 
# -0.64860328  0.02786736 -0.69333435 -0.60141860 
# -0.650       0.027      -0.693      -0.606
## PROBABILITY OF REDUCTION BEING P <= -0.5
100*length(which(redBV[13,] <= -0.5))/length(redBV[13,])
# [1] 100
##  ENDANGERED
acBV <- (PCBV/P3BV)^(1/(2025-1999))
round(c(mean=mean(-100*(1-acBV)),sd=sd(-100*(1-acBV)),quantile(-100*(1-acBV),c(0.025,0.975))),3)
#   mean     sd   2.5%  97.5% 
# -3.954  0.293 -4.540 -3.400
#
##------------------------
##  EAST ANTARCTICA
##------------------------
P3EA <- redMI[11,]+PEI+redCZ[11,]+redKG[11,]+redHI[11,]+redMQ[11,]
PCEA <- redMI[12,]+PEI+redCZ[11,]+redKG[12,]+redHI[4,]+redMQ[4,]
acEA <- (PCEA/P3EA)^(1/(2025-1999))
round(c(mean=mean(-100*(1-acEA)),sd=sd(-100*(1-acEA)),quantile(-100*(1-acEA),c(0.025,0.975))),3)
#  mean    sd  2.5% 97.5% 
# 2.620 0.088 2.453 2.800 
# N_3G
round(c(mean=mean(P3EA),SD=sd(P3EA),quantile(P3EA,c(0.025,0.975))))
#   mean    SD  2.5% 97.5% 
#  10474   410  9716 11318 
# N_25
round(c(mean=mean(PCEA),SD=sd(PCEA),quantile(PCEA,c(0.025,0.975))))
#  mean    SD  2.5% 97.5% 
# 20514   710 19211 21975
##  LEAST CONCERN
##
########################
## GLOBAL POPULATION
##---------------------
P3G <- redSG[11,]+redSSD[11,]+redSSI[11,]+redBV[11,]+P3EA
paste(round(mean(P3G))," (",round(sd(P3G)),")",sep="")
# [1] "2187424 (153141)"
PCG <- redSG[12,]+redSSD[12,]+redSSI[12,]+redBV[12,]+PCEA
paste(round(mean(PCG))," (",round(sd(PCG)),")",sep="")
# [1] "943756 (89727)"
C3GG <- (PCG-P3G)/P3G
c(mean=mean(C3GG),sd=sd(C3GG),quantile(C3GG,c(0.05,0.95)))
#       mean          sd          5%         95% 
# -0.56658609  0.04997881 -0.64326633 -0.48016124 
## PROBABILITY OF REDUCTION BEING P <= -0.5
100*length(which(C3GG <= -0.5))/length(C3GG)
# [1] 90.31
##  ENDANGERED
# ANNUAL CHANGE
acG <- (PCG/P3G)^(1/(2025-1999))
round(c(mean=mean(-100*(1-acG)),sd=sd(-100*(1-acG)),quantile(-100*(1-acG),c(0.025,0.975))),3)
#  mean     sd   2.5%  97.5% 
# -3.188  0.427 -4.015 -2.347 
################################################################################
## PROPORTIONS OF GLOBAL POPULATION
##----------------------------------
## south Georgia
round(c(mean=mean(PCSG/PCG*100),SD=sd(PCSG/PCG*100)),2)
#  mean    SD 
# 96.40  0.37
## South Sandwich Islands
round(c(mean=mean(redSSD[12,]/PCG*100),SD=sd(redSSD[12,]/PCG*100)),2)
#  mean   SD 
#  0.40 0.04 
## South Shetlands
round(c(mean=mean(PCSS/PCG*100),SD=sd(PCSS/PCG*100)),2)
#  mean    SD 
#  0.14  0.01
### Bouvetoya
round(c(mean=mean(PCBV/PCG*100),SD=sd(PCBV/PCG*100)),2)
# mean   SD 
# 0.86 0.11
## EAST ANTARCTICA
round(c(mean=mean(PCEA/PCG*100),SD=sd(PCEA/PCG*100)),2)
# mean   SD 
# 2.19 0.22
## SOUTH ATLANTIC REGION
round(c(mean=mean((PCSG+redSSD[12,]+PCSS+PCBV)/PCG*100),SD=sd((PCSG+redSSD[12,]+PCSS+PCBV)/PCG*100)),2)
# mean    SD 
# 97.81  0.2
