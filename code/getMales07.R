## FILTER AND ABUNDANCE ESTIMATION FUNCTION FOR MALES
## STAGE 2 IS ALL MATURE MALES (INCLUDES SKIPPERS); STAGE 1 IS TERRITORIAL MALES
## RATIOS OF TOTAL MATURE MALES AND TERRITORIALS TO PREDICTED SSB COUNTS:
#  -> NT.ncs RATIO FOR TERRITORIALS
#  -> NATot.ncs RATIO FOR MATURE
getMales07 <- function(i, stage = 1){
  if(sum(SGC07[i, c(1, 4)]) > 0){
    males <- sum(SGC07[i, c(2, 5)])
    if(stage == 1){
      res <- males * NT.nc07
    }
    if(stage == 2){
      res <- res <- males * NATot.nc07
    }
  }else{
    res <- rep(0, 15000)
  }  
  res  
}