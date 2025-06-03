##----------------------------------------------------------------------------------------#
## - FOR CONSISTENCY WITH PUBLISHED ESTIMATES FOR 2007-09, I USE SAME DATA FILES AND CODE -
## FILTER AND ABUNDANCE ESTIMATION FUNCTION FOR FEMALES
## STAGE 1 IS MATURE FEMALES; STAGE 2 IS BREEDING FEMALES ONLY (OR PRODUCTIVITY)
## unique survey days are unique(SGC09[, "day"]); i.e. 3 4 5 6 8 9
getFemales <- function(i, stage = 1, females = F){
  pupsi <- sum(SGC09[i, c(1, 4)])
  if(pupsi > 0){
    ftop <- sum(SGC09[i, c(3, 6)]) / pupsi                  ## observed ratio in area i
    FtoP <- femToPup95[SGC09[i, "day"] - 2]                 ## predicted ratio from SSB
    if(ftop > FtoP){
      fems <- femToPup[SGC09[i, "day"] - 2, ] * pupsi
      if(mean(fems) < pupsi & (SGC09[i, "day"] - 2) > 4)
        fems <- rep(pupsi, 15000) 
    }else{
      fems <- rep(sum(SGC09[i, c(3, 6)]), 15000)
    }
  }else{
    fems <- rep(0, 15000)
  }
  if(females){
    res <- fems
  }else{
    if(stage == 1){
      res <- NMats[SGC09[i, "day"] - 2, ] * fems
    }
    if(stage == 2){
      res <- NBIm[SGC09[i, "day"] - 2, ] * fems
    }    
  } 
  res  
}