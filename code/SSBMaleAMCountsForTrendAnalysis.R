################################################################################
### BIRD ISLAND - SSB AM MALE COUNT ANALYSIS
##------------------------------------------
### LOAD PACKAGES                           
library(mgcv)
library(jagsUI)
###--------------
##  GET AFTERNOON FEMALE COUNT DATA FROM SSB
##   - TOTAL SSB FEMALE COUNTS - 1989-2022 -
##  1. model observed daily counts n_j (day j=1,...,70) with Poisson GAMMs
##  2. predict expected SSB count for survey days at and around peak
##  GET DATA
tmc <- read.csv("SSB_total_male_counts.csv", header = T)
day <- 1:61; tmc <- tmc[day,]
##----------------------------
##  PREDICT PEAK PUPPING DAYS
newpups <- read.csv("newpups.csv",header=F)
colnames(newpups) <- 1985:2025
sumpup <- apply(newpups,2,cumsum)
## FUNCTION ffs.obs FITS AN NLME (LOGISTIC) MODEL TO OBTAIN THE PEAK PUPPING DAY
source("ffs.obs.R")
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
##----------------------------------------------------------------------------------#
##  FIT Poisson GAMM MODELS
##  1. select data set for specific survey year and create data.frame to fit a GAMM
##  2. set ile name (e.g. "test.jags") for file with jags code
##  3. set up jags code and data.
##  4. Use BUGS code in text.jags to generate the jags model and run the model
##  5. Extract required parameters, including predictions ("mu") of required counts
##  6. Save results to *.Rds file
##----
##  1994-95
dat <- data.frame(y=tmc[,paste("X",1995,sep="")],x0=day)
jags.file <- "test.jags"
jd <- jagam(y ~ s(x0), family = poisson(link = "log"), data = dat, file = jags.file, sp.prior = "gamma", diagonalize = TRUE)
## SET MODEL IN JAGS
## THEN A POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
sink("SSBMal.GAMM") 
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
    ## Parametric effect priors CHECK tau=1/34^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,00087)
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
out.test <- jags(jags.data,inits,parameters,"SSBMal.GAMM",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,parallel=T)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1] 3.19 2.62 3.88
## SELECT PREDICTIONS; peak day is mu[40]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]],out.test$samples[[2]],out.test$samples[[3]])
preds1 <- preds[,which(dimnames(preds)[[2]]==paste("mu[",min(pdays(predday[1])),"]",sep="")):
                  which(dimnames(preds)[[2]]==paste("mu[",max(pdays(predday[1])),"]",sep=""))]
saveRDS(preds1,"SSBPredictedAMCount9495.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]]=="mu[1]"):which(dimnames(preds)[[2]]=="mu[61]")]
saveRDS(preds1, "SSBPredictedAMCount9495ALL.Rds")
##
##----
##  1997-98
dat <- data.frame(y = tmc[,paste("X",1998,sep="")],x0=day)
jags.file <- "test.jags"
jd <- jagam(y ~ s(x0),family=poisson(link="log"),data=dat,file=jags.file,sp.prior="gamma",diagonalize = TRUE)
## SET MODEL IN JAGS
## THEN A POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
sink("SSBMal.GAMM") 
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
    ## Parametric effect priors CHECK tau=1/37^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,0.00075)
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
out.test <- jags(jags.data, inits, parameters, "SSBMal.GAMM", n.chains = nc, n.thin =  nt, n.iter = ni, n.burnin = nb,parallel=T)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1] 0.099 0.047 0.157
## SELECT PREDICTIONS; peak day is mu[46]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
preds1 <- preds[, which(dimnames(preds)[[2]] == paste("mu[", min(pdays(predday[2])), "]", sep = "")):
                  which(dimnames(preds)[[2]] == paste("mu[", max(pdays(predday[2])), "]", sep = ""))]
saveRDS(preds1, "SSBPredictedAMCount9798.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]] == "mu[1]"):which(dimnames(preds)[[2]] == "mu[61]")]
saveRDS(preds1, "SSBPredictedAMCount9798ALL.Rds")
##
##----
##  2003-04
dat <- data.frame(y = tmc[, paste("X", 2004, sep = "")], x0 = day)
jags.file <- "test.jags"
jd <- jagam(y ~ s(x0), family = poisson(link = "log"), data = dat, file = jags.file, sp.prior = "gamma", diagonalize = TRUE)
## SET MODEL IN JAGS
## THEN A POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
sink("SSBMal.GAMM") 
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
    ## Parametric effect priors CHECK tau=1/35^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,8e-04)
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
out.test <- jags(jags.data, inits, parameters, "SSBMal.GAMM", n.chains = nc, n.thin =  nt, n.iter = ni, n.burnin = nb,parallel=T)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1] 0.149 0.113 0.195
## SELECT PREDICTIONS; peak day is mu[43]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
preds1 <- preds[, which(dimnames(preds)[[2]] == paste("mu[", min(pdays(predday[3])), "]", sep = "")):
                  which(dimnames(preds)[[2]] == paste("mu[", max(pdays(predday[3])), "]", sep = ""))]
saveRDS(preds1, "SSBPredictedAMCount0304.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]] == "mu[1]"):which(dimnames(preds)[[2]] == "mu[61]")]
saveRDS(preds1, "SSBPredictedAMCount0304ALL.Rds")
##
##----
##  2008-09
dat <- data.frame(y = tmc[,paste("X",2009,sep="")],x0=tmc$day)
jags.file <- "test.jags"
jd <- jagam(y ~ s(x0), family = poisson(link = "log"), data = dat, file = jags.file, sp.prior = "gamma", diagonalize = TRUE)
## SET MODEL IN JAGS
## THEN A POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
sink("SSBMal.GAMM") 
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
    ## Parametric effect priors CHECK tau=1/34^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,0.00088)
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
out.test <- jags(jags.data, inits, parameters, "SSBMal.GAMM", n.chains = nc, n.thin =  nt, n.iter = ni, n.burnin = nb,parallel=T)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1] 0.149 0.113 0.195
## SELECT PREDICTIONS; peak day is mu[43]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
preds1 <- preds[, which(dimnames(preds)[[2]] == paste("mu[", min(pdays(predday[4])), "]", sep = "")):
                  which(dimnames(preds)[[2]] == paste("mu[", max(pdays(predday[4])), "]", sep = ""))]
saveRDS(preds1, "SSBPredictedAMCount0809.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]] == "mu[1]"):which(dimnames(preds)[[2]] == "mu[61]")]
saveRDS(preds1, "SSBPredictedAMCount0809ALL.Rds")
##
##----
##  2016-17
dat <- data.frame(y = tmc[, paste("X", 2017, sep = "")], x0 = day)
jags.file <- "test.jags"
jd <- jagam(y ~ s(x0), family = poisson(link = "log"), data = dat, file = jags.file, sp.prior = "gamma", diagonalize = TRUE)
## SET MODEL IN JAGS
## THEN A POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
sink("SSBMal.GAMM") 
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
    ## Parametric effect priors CHECK tau=1/30^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,0.0011)
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
out.test <- jags(jags.data, inits, parameters, "SSBMal.GAMM", n.chains = nc, n.thin =  nt, n.iter = ni, n.burnin = nb,parallel=T)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1] 0.148 0.100 0.208
## SELECT PREDICTIONS; peak day is mu[44]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
preds1 <- preds[, which(dimnames(preds)[[2]] == paste("mu[", min(pdays(predday[5])), "]", sep = "")):
                  which(dimnames(preds)[[2]] == paste("mu[", max(pdays(predday[5])), "]", sep = ""))]
saveRDS(preds1, "SSBPredictedAMCount1617.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]] == "mu[1]"):which(dimnames(preds)[[2]] == "mu[61]")]
saveRDS(preds1, "SSBPredictedAMCount1617ALL.Rds")
##
##----
##  2023-24
dat <- data.frame(y = tmc[, paste("X", 2024, sep = "")], x0 = day)
jags.file <- "test.jags"
jd <- jagam(y ~ s(x0), family = poisson(link = "log"), data = dat, file = jags.file, sp.prior = "gamma", diagonalize = TRUE)
## SET MODEL IN JAGS
## THEN A POISSON GAMM, TO ACCOUNT FOR OVERDISPERSION
sink("SSBMal.GAMM") 
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
    ## Parametric effect priors CHECK tau=1/32^2 is appropriate!
    for(i in 1:1){
      b[i] ~ dnorm(0,0.00099)
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
out.test <- jags(jags.data, inits, parameters, "SSBMal.GAMM", n.chains = nc, n.thin =  nt, n.iter = ni, n.burnin = nb,parallel=T)
Sys.time() - sta
## INFERENCE
## OVERDISPERSION
round(c(out.test$mean$sd, out.test$q2.5$sd, out.test$q97.5$sd), 3)
# [1] 0.148 0.100 0.208
## SELECT PREDICTIONS; peak day is mu[44]
pdays <- function(x) (ceiling(x) - 3):(ceiling(x) + 3)
preds <- rbind(out.test$samples[[1]], out.test$samples[[2]], out.test$samples[[3]])
preds1 <- preds[, which(dimnames(preds)[[2]] == paste("mu[", min(pdays(predday[6])), "]", sep = "")):
                  which(dimnames(preds)[[2]] == paste("mu[", max(pdays(predday[6])), "]", sep = ""))]
saveRDS(preds1, "SSBPredictedAMCount2324.Rds")
preds1 <- preds[, which(dimnames(preds)[[2]] == "mu[1]"):which(dimnames(preds)[[2]] == "mu[61]")]
saveRDS(preds1, "SSBPredictedAMCount2324ALL.Rds")
##