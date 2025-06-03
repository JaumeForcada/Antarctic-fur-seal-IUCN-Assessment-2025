##  Gaussian band pass filter
##  From Cazelles et al. 2008
##  For time series of 30 – 40 data points with significant periodic components smaller than 20–25% of the series length.
##  e.g., for a 35 year span, .20 * 35 = 7
#   timescale cicles of more than 6 years need to be filtered; for the Gaussian filter we select pf >= 7 / 2
##  
glbpf <- function(y, pf){ 
  ts <- pf
  ts2 <- ts^2
  ts3 <- ts * 3     # effects of the filtering are limited to 3 ts
  errbar <- 0.5     # rms error of 0.5
  #
  ndon <- length(y)
  hano <- y
  hanolp <- rep(0, ndon)
  #
  for(iter in 1:20){
    sumerr <- 0
    delh <- hano - hanolp
    for(i in 1:ndon){
      t1 <- as.numeric(i)
      sumwt <- 0
      cor <- 0
      for(j in 1:ndon){
        t2 <- as.numeric(j)
        if(t2 <= t1 + ts3){
          dt <- abs(t1 - t2)
          if(dt <= ts3){
            dt2 <- dt^2
            wt <- exp(- dt2 / (1.44 * ts2))
            sumwt <- sumwt + wt
            cor <- cor + wt * delh[j]
          }
        }else{
          j <- ndon + 1
        }
      }
      h1 <- hanolp[i]
      if(sumwt != 0){ 
        hanolp[i] <- hanolp[i] + cor / sumwt
      }
      h2 <- hanolp[i]
      sumerr <- sumerr + (h2 - h1)^2
    }
    rmsdh <- sqrt(sumerr / ndon)
    if(rmsdh <= errbar) 
      i <- ndon + 1
  }
  yf <- hanolp
  return(yf)
}
#
