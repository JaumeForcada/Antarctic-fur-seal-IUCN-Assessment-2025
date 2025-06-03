##----------------------------------------------------------------------------------------#
## - FOR CONSISTENCY WITH PUBLISHED ESTIMATES FOR 2007-09, I USE SAME DATA FILES AND CODE -
## FILTER AND ABUNDANCE ESTIMATION FUNCTION FOR FEMALES
## STAGE 1 IS MATURE FEMALES; STAGE 2 IS BREEDING FEMALES ONLY (OR PRODUCTIVITY)
## unique survey days are unique(SGC09[, "day"]); i.e. 3 4 5 6 8 9
getFemales07 <- function(i, stage = 1, females = F){
  pupsi <- sum(SGC07[i, c(1, 4)])
  if(pupsi == 0){
    fems <- rep(0, 15000)
  }else{
    fems <- rep(sum(SGC07[i, c(3, 6)]), 15000)
  }
  if(females){
    res <- fems
  }else{
    if(stage == 1){
      res <- NMat.nc07 * fems
    }
    if(stage == 2){
      res <- NBI.nc07 * fems
    }    
  } 
  res  
}