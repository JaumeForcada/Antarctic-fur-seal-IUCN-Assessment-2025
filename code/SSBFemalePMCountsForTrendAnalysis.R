################################################################################
### BIRD ISLAND PRODUCTIVITY TRENDS                                             
##------------------------------------------------------------------------------
### LOAD PACKAGES
library(mgcv)
library(jagsUI)
##-------------------------------------------------------------------------------------------------------------------------#
##  GET AFTERNOON FEMALE COUNT DATA FROM SSB
##   - TOTAL SSB FEMALE COUNTS - 1989-2024 -
##  1. model observed daily counts n_j (day j=1,...,70) with Poisson GAMMs
##  2. predict expected SSB count for survey days at and around peak
##  GET DATA
tfemc <- read.csv("raw_data/SSB_total_female_counts.csv",header = T)
## Adding sequential day starting on November 1st
day <- 1:70
##-------------------------------------------------------------------------------------------------------------------------#
##  PREDICT PEAK PUPPING DAYS
newpups <- read.csv("raw_data/newpups.csv",header=F)
colnames(newpups) <- 1985:2025
sumpup <- apply(newpups,2,cumsum)   # CUMULATIVE SUM
#
## FUNCTION ffs.obs FITS AN NLME (LOGISTIC) MODEL TO OBTAIN THE PEAK PUPPING DAY
source(ffs.obs.R)
## YEARS WITH COMPLETE BIRD ISLAND SURVEYS
syears <- c(1989,1990,1991,1995,1998,2004,2009,2017,2024)
pks <- ses <- numeric(length(syears))
for(i in 1:length(syears)){
  pks[i] <- as.vector(ffs.obs(data.frame(day = day, pups = sumpup[, paste(syears[i])]), plots = F)$coefficients[2])
  ses[i] <- sqrt(ffs.obs(data.frame(day = day, pups = sumpup[, paste(syears[i])]), plots = F)$varBeta[2, 2])
}  
## PEAK PUPPING DATES ARE:
round(pks, 1)
# [1] 37.6 34.3 39.4 43.4 44.0 37.4 36.3 43.4 38.4
predday <- c(38.5, 42.5, 44, 40, 46, 42.5, 36, 44, 38.4)
##------------------------------------------------------------------------------
##  FIT Poisson GAMM MODELS
##  1. select data set for specific survey year and create data.frame to fit a GAMM
##  2. set file name (e.g. "test.jags") for file with jags code
##  3. set up jags code and data with function jagam from package mgcv
##  4. Use BUGS code in text.jags to generate the jags model and run the model
##  5. Extract required parameters, including predictions ("mu") of required counts
##  6. Save results to *.Rds file
##----
##  1988-89
dat <- data.frame(y = tfemc[, paste("X", 1989, sep = "")], x0 = day)
jags.file <- paste("C:/workstuff/AntarcticFurSeals/Fur_seal_Survey/analysis", "/test.jags", sep = "") 
jd <- jagam(y ~ s(x0), family = poisson(link = "log"), data = dat, file = jags.file, sp.prior = "gamma", diagonalize = TRUE)
## SET MODEL IN JAGS
## THEN A POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
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
    ## Parametric effect priors CHECK tau=1/52^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,0.00037) # 0.00026)
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
parameters <- c("b", "rho", "mu", "sd", "eps")
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
### Run the BUGS model
sta <- Sys.time()
out.test <- jags(jags.data, inits, parameters, "SSBFem.GAMM", n.chains = nc, n.thin =  nt, n.iter = ni, n.burnin = nb)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1] 0.172 0.134 0.219
## SELECT PREDICTIONS; peak day is mu[39]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
preds1 <- preds[, which(dimnames(preds)[[2]] == paste("mu[", min(pdays(predday[1])), "]", sep = "")):
                 which(dimnames(preds)[[2]] == paste("mu[", max(pdays(predday[1])), "]", sep = ""))]
saveRDS(preds1, "SSBPredictedPMCount8889.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]] == "mu[1]"):which(dimnames(preds)[[2]] == "mu[70]")]
saveRDS(preds1, "SSBPredictedPMCount8889ALL.Rds")
##
##----
##  1989-90
dat <- data.frame(y = tfemc[, paste("X", 1990, sep = "")], x0 = day)
jags.file <- paste("C:/workstuff/AntarcticFurSeals/Fur_seal_Survey/analysis", "/test.jags", sep = "") 
jd <- jagam(y ~ s(x0), family = poisson(link = "log"), data = dat, file = jags.file, sp.prior = "gamma", diagonalize = TRUE)
## SET MODEL IN JAGS
## THEN A POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
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
    ## Parametric effect priors CHECK tau=1/52^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,0.00038)
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
parameters <- c("b", "rho", "mu", "sd", "eps")
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
### Run the BUGS model
sta <- Sys.time()
out.test <- jags(jags.data, inits, parameters, "SSBFem.GAMM", n.chains = nc, n.thin =  nt, n.iter = ni, n.burnin = nb)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1] [1] 0.199 0.154 0.256
## SELECT PREDICTIONS; peak day is mu[43]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
preds1 <- preds[, which(dimnames(preds)[[2]] == paste("mu[", min(pdays(predday[2])), "]", sep = "")):
                  which(dimnames(preds)[[2]] == paste("mu[", max(pdays(predday[2])), "]", sep = ""))]
saveRDS(preds1, "SSBPredictedPMCount8990.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]] == "mu[1]"):which(dimnames(preds)[[2]] == "mu[64]")]
saveRDS(preds1, "SSBPredictedPMCount8990ALL.Rds")
##
##----
##  1990-91
dat <- data.frame(y = tfemc[, paste("X", 1991, sep = "")], x0 = day)
jags.file <- paste("C:/workstuff/AntarcticFurSeals/Fur_seal_Survey/analysis", "/test.jags", sep = "") 
jd <- jagam(y ~ s(x0), family = poisson(link = "log"), data = dat, file = jags.file, sp.prior = "gamma", diagonalize = TRUE)
## SET MODEL IN JAGS
## THEN A POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
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
    ## Parametric effect priors CHECK tau=1/52^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,0.00049)
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
parameters <- c("b", "rho", "mu", "sd", "eps")
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
### Run the BUGS model
sta <- Sys.time()
out.test <- jags(jags.data, inits, parameters, "SSBFem.GAMM", n.chains = nc, n.thin =  nt, n.iter = ni, n.burnin = nb)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1] 0.064 0.021 0.108
## SELECT PREDICTIONS; peak day is mu[44]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
preds1 <- preds[, which(dimnames(preds)[[2]] == paste("mu[", min(pdays(predday[3])), "]", sep = "")):
                  which(dimnames(preds)[[2]] == paste("mu[", max(pdays(predday[3])), "]", sep = ""))]
saveRDS(preds1, "SSBPredictedPMCount9091.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]] == "mu[1]"):which(dimnames(preds)[[2]] == "mu[70]")]
saveRDS(preds1, "SSBPredictedPMCount9091ALL.Rds")
##
##----
##  1994-95
dat <- data.frame(y = tfemc[, paste("X", 1995, sep = "")], x0 = day)
jags.file <- paste("C:/workstuff/AntarcticFurSeals/Fur_seal_Survey/analysis", "/test.jags", sep = "") 
jd <- jagam(y ~ s(x0), family = poisson(link = "log"), data = dat, file = jags.file, sp.prior = "gamma", diagonalize = TRUE)
## SET MODEL IN JAGS
## THEN A POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
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
    ## Parametric effect priors CHECK tau=1/52^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,0.00048)
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
parameters <- c("b", "rho", "mu", "sd", "eps")
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
### Run the BUGS model
sta <- Sys.time()
out.test <- jags(jags.data, inits, parameters, "SSBFem.GAMM", n.chains = nc, n.thin =  nt, n.iter = ni, n.burnin = nb)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1] 0.057 0.017 0.097
## SELECT PREDICTIONS; peak day is mu[40]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
preds1 <- preds[, which(dimnames(preds)[[2]] == paste("mu[", min(pdays(predday[4])), "]", sep = "")):
                  which(dimnames(preds)[[2]] == paste("mu[", max(pdays(predday[4])), "]", sep = ""))]
saveRDS(preds1, "SSBPredictedPMCount9495.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]] == "mu[1]"):which(dimnames(preds)[[2]] == "mu[61]")]
saveRDS(preds1, "SSBPredictedPMCount9495ALL.Rds")
##
##----
##  1997-98
dat <- data.frame(y = tfemc[, paste("X", 1998, sep = "")], x0 = day)
jags.file <- paste("C:/workstuff/AntarcticFurSeals/Fur_seal_Survey/analysis", "/test.jags", sep = "") 
jd <- jagam(y ~ s(x0), family = poisson(link = "log"), data = dat, file = jags.file, sp.prior = "gamma", diagonalize = TRUE)
## SET MODEL IN JAGS
## THEN A POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
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
    ## Parametric effect priors CHECK tau=1/52^2 is appropriate!
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
jags.data <- list(y = jd$jags.data$y, n = jd$jags.data$n, X = jd$jags.data$X)
## INITIAL VALUES
inits <- function()
  list(b = jd$jags.ini$b, lambda = jd$jags.ini$lambda)
## PARAMETRES MONITORED
parameters <- c("b", "rho", "mu", "sd", "eps")
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
### Run the BUGS model
sta <- Sys.time()
out.test <- jags(jags.data, inits, parameters, "SSBFem.GAMM", n.chains = nc, n.thin =  nt, n.iter = ni, n.burnin = nb)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1] 0.099 0.047 0.157
## SELECT PREDICTIONS; peak day is mu[46]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
preds1 <- preds[, which(dimnames(preds)[[2]] == paste("mu[", min(pdays(predday[5])), "]", sep = "")):
                  which(dimnames(preds)[[2]] == paste("mu[", max(pdays(predday[5])), "]", sep = ""))]
saveRDS(preds1, "SSBPredictedPMCount9798.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]] == "mu[1]"):which(dimnames(preds)[[2]] == "mu[61]")]
saveRDS(preds1, "SSBPredictedPMCount9798ALL.Rds")
##
##----
##  2003-04
dat <- data.frame(y = tfemc[, paste("X", 2004, sep = "")], x0 = day)
jags.file <- paste("C:/workstuff/AntarcticFurSeals/Fur_seal_Survey/analysis", "/test.jags", sep = "") 
jd <- jagam(y ~ s(x0), family = poisson(link = "log"), data = dat, file = jags.file, sp.prior = "gamma", diagonalize = TRUE)
## SET MODEL IN JAGS
## THEN A POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
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
    ## Parametric effect priors CHECK tau=1/52^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,0.00037)
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
parameters <- c("b", "rho", "mu", "sd", "eps")
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
### Run the BUGS model
sta <- Sys.time()
out.test <- jags(jags.data, inits, parameters, "SSBFem.GAMM", n.chains = nc, n.thin =  nt, n.iter = ni, n.burnin = nb)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1] 0.149 0.113 0.195
## SELECT PREDICTIONS; peak day is mu[43]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
preds1 <- preds[, which(dimnames(preds)[[2]] == paste("mu[", min(pdays(predday[6])), "]", sep = "")):
                  which(dimnames(preds)[[2]] == paste("mu[", max(pdays(predday[6])), "]", sep = ""))]
saveRDS(preds1, "SSBPredictedPMCount0304.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]] == "mu[1]"):which(dimnames(preds)[[2]] == "mu[70]")]
saveRDS(preds1, "SSBPredictedPMCount0304ALL.Rds")
##
##----
##  2016-17
dat <- data.frame(y = tfemc[, paste("X", 2017, sep = "")], x0 = day)
jags.file <- paste("C:/workstuff/AntarcticFurSeals/Fur_seal_Survey/analysis", "/test.jags", sep = "") 
jd <- jagam(y ~ s(x0), family = poisson(link = "log"), data = dat, file = jags.file, sp.prior = "gamma", diagonalize = TRUE)
## SET MODEL IN JAGS
## THEN A POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
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
    ## Parametric effect priors CHECK tau=1/52^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,0.00057)
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
parameters <- c("b", "rho", "mu", "sd", "eps")
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
### Run the BUGS model
sta <- Sys.time()
out.test <- jags(jags.data, inits, parameters, "SSBFem.GAMM", n.chains = nc, n.thin =  nt, n.iter = ni, n.burnin = nb)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1] 0.148 0.100 0.208
## SELECT PREDICTIONS; peak day is mu[44]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
preds1 <- preds[, which(dimnames(preds)[[2]] == paste("mu[", min(pdays(predday[8])), "]", sep = "")):
                  which(dimnames(preds)[[2]] == paste("mu[", max(pdays(predday[8])), "]", sep = ""))]
saveRDS(preds1, "SSBPredictedPMCount1617.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]] == "mu[1]"):which(dimnames(preds)[[2]] == "mu[68]")]
saveRDS(preds1, "SSBPredictedPMCount1617ALL.Rds")
##
##----
##  2023-24
## SPLINE INTERPOLATION TO ACCOUNT FOR H5N1 SEAL CESSATION WORK
AFSftmp <- tfemc[1:70,paste("X",2024,sep="")]
aicd <- NULL
aicd <- 32:34
aicd2 <- 32:35
aicd3 <- 31:35
lstrcd <- max(which(!is.na(AFSftmp)))
tmprcs <- round(spline(x=which(!is.na(AFSftmp)),y=AFSftmp[which(!is.na(AFSftmp))],xout=1:lstrcd)$y)
AFSftmp[aicd3] <- tmprcs[aicd3]
dat <- data.frame(y=AFSftmp,x0=day)
rm(AFSftmp,aicd3,tmprcs,aicd,aicd2,lstrcd)
jags.file <- "test.jags"
jd <- jagam(y ~ s(x0),family=poisson(link="log"),data=dat,file=jags.file,sp.prior="gamma",diagonalize=T)
## SET MODEL IN JAGS
## POISSON GAMM
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
    ## Parametric effect priors CHECK tau=1/43^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,0.00054)
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
parameters <- c("b", "rho", "mu", "sd", "eps")
## MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3
### Run the BUGS model
sta <- Sys.time()
out.test <- jags(jags.data, inits, parameters, "SSBFem.GAMM", n.chains = nc, n.thin =  nt, n.iter = ni, n.burnin = nb,parallel=T)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1]  0.178 0.126 0.240
## SELECT PREDICTIONS; peak day is mu[38]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
preds1 <- preds[, which(dimnames(preds)[[2]] == paste("mu[", min(pdays(predday[9])), "]", sep = "")):
                  which(dimnames(preds)[[2]] == paste("mu[", max(pdays(predday[9])), "]", sep = ""))]
saveRDS(preds1, "SSBPredictedPMCount2324.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]] == "mu[1]"):which(dimnames(preds)[[2]] == "mu[68]")]
saveRDS(preds1, "SSBPredictedPMCount2324ALL.Rds")
##