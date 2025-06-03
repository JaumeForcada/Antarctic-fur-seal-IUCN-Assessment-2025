ffs.obs <- function(lgframe = pups04, add = F, obs = T, tl = 1, plots = T, labs = NULL , pleg = T, bgs = NA, cols = 1){
  fit.obs <- gnls(pups ~ Asym / (1 + exp((xmid - day) / scal)), 
                  data = lgframe, start =  c(Asym = 550 , xmid = 40, scal = 4), na.action = "na.omit")
  new <- data.frame(day = seq(1, 70, by = 1))
  preds <- predict(fit.obs, newdata = new)
  if(plots){
    if(add){
      lines(new$day, predict(fit.obs, newdata = new), lty = tl, col = cols) 
      points(coef(fit.obs)[2], predict(fit.obs, newdata = data.frame(day = coef(fit.obs)[2])), pch = 19, cex = 0.6, col = cols)
      if(missing(labs)) 
        labs <- paste(lgframe[70, 2])
      points(coef(fit.obs)[2] + 40, lgframe[70, 2], pch = 19, cex = 0.6, col = 2)
      text(coef(fit.obs)[2] + 42, lgframe[70, 2], labels = labs, cex = 0.5)
      if(obs) 
        points(lgframe$day, lgframe$pups, pch = 1, col = cols, bg = bgs, cex = .4)
    }
    else {
      plot(new$day, predict(fit.obs, newdata = new), type = 'l', xlab = "",
           ylab = "number of pups born at SSB", axes = F, xlim = c(0, 100), ylim = c(0, 900))
      if(pleg){
        lines(c(73, 90), rep(900, 2), col = 2)
        lines(c(73, 90), rep(150, 2), col = 2)
        for(i in 73:90){
          lines(rep(i, 2), c(895, 900), col = 2)
          lines(rep(i, 2), c(145, 155), col = 2)
        }
        lines(rep(79, 2), c(890, 900), col = 2)
        lines(rep(79, 2), c(145, 160), col = 2)
        text(73, 130, paste("3"), col = 2, cex = .6)
        text(79, 130, paste("9"), col = 2, cex = .6)
        text(85, 130, paste("15"), col = 2, cex = .6)
        text(90, 130, paste("20"), col = 2, cex = .6)
        text(92, 130, paste("Dec"), pos = 4, col = 2, cex = .6)
        text(73, 105, "peak season", pos = 4, col = 2, cex = .7)
        lines(rep(79, 2), c(150, 900), col = 2, lty = 3)
      }        
      points(coef(fit.obs)[2], predict(fit.obs, newdata = data.frame(day=coef(fit.obs)[2])),pch=20,cex=0.6)
      axis(1, at = c(1,15,31,45,62), labels = paste(c("1-Nov","15-Nov","1-Dec","15-Dec","1-Jan")), tck = -0.02, cex.axis = 0.75)
      axis(1, at = seq(1,70, by = 1), labels = F, tck = -0.012, cex.axis = 0.75)
      axis(2, at = seq(0,900, by=100), labels = T, tck = -0.02, cex.axis = 0.75, las = 2)
      axis(2, at = seq(0,900, by=25),labels = F, tck = -0.012, cex.axis = 0.75, las = 2)
      axis(3, at = c(1,15,31,45,62), labels = paste(c("1-Nov","15-Nov","1-Dec","15-Dec","1-Jan")), tck = -0.02, cex.axis = 0.75)
      axis(3, at = seq(1,70, by = 1), labels = F, tck = -0.012, cex.axis = 0.75)
      axis(4, at = seq(0,900, by=100), labels = T, tck = -0.02, cex.axis = 0.75, las = 2)
      axis(4, at = seq(0,900, by=25),labels = F, tck = -0.012, cex.axis = 0.75, las = 2)
      box()
      if(missing(labs)) 
        labs <- paste(lgframe[70, 2])
      points(coef(fit.obs)[2] + 40, lgframe[70, 2], pch = 20, cex = 0.6, col = 2)
      text(coef(fit.obs)[2] + 44, lgframe[70, 2], labels = labs, cex = 0.5)
      #
      if(obs) 
        points(lgframe$day, lgframe$pups, pch = 20)
    }
  }
  summary(fit.obs)
}
