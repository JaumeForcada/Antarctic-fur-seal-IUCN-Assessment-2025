## RATIO OF PUP PRODUCTION TO TOTAL POPULATION TO SCALE UP PUP COUNTS TO MATURE AND TOTAL POPULATION SIZE
##--------------------------
## PUP PRODUCTION 1984:2005
Prod <- c(850,678,854,820,777,757,822,545,718,770,714,591,801,733,536,693,446,696,747,713,769,652,464,768,736,782,271,568,547,584,258,525,403,361,496,282,409,368,247,335,303,368)
## PRODUCTIVITY FOR 2001:2025
cbind(2001:2025,Prod[18:42])
## IPM RESULTS
preds <- readRDS("processed_data/AFSFemaleSimulations200125.Rds")
predsM <- readRDS("processed_data/AFSMaleSimulations19952025.Rds")
#
## GENERATE IPM ABUNDANCE VECTORS FOR 2001:2025
ProdToMature <- ProdToTotal <- matrix(NA,15000,25)
for(i in 1:25){
  ProdToMature[,i] <- (preds[,which(dimnames(preds)[[2]]==paste("NN[",i,"]",sep=""))]+
                         preds[,which(dimnames(preds)[[2]]==paste("NBI[",i,"]",sep=""))]+
                         predsM[,which(dimnames(predsM)[[2]]==paste("NT[",6+i,"]",sep=""))]+
                         predsM[,which(dimnames(predsM)[[2]]==paste("NN[",6+i,"]",sep=""))]) / Prod[17+i]
  ProdToTotal[,i] <- (preds[,which(dimnames(preds)[[2]]==paste("NN[",i,"]",sep=""))]+
                        preds[,which(dimnames(preds)[[2]]==paste("NBI[",i,"]",sep=""))]+
                        preds[,which(dimnames(preds)[[2]]==paste("Npre[",i,"]",sep=""))]+
                        predsM[,which(dimnames(predsM)[[2]]==paste("NT[",6+i,"]",sep=""))]+
                        predsM[,which(dimnames(predsM)[[2]]==paste("NN[",6+i,"]",sep=""))]+
                        predsM[,which(dimnames(predsM)[[2]]==paste("NPall[",6+i,"]",sep=""))]) / Prod[17+i]
  print(i)
}
## MEAN 2001:2005
PtoM <- rowMeans(ProdToMature[,1:5])
PtoT <- rowMeans(ProdToTotal[,1:5])
##
## SUMMARIES
round(c(mean=mean(PtoM),quantile(PtoM,c(0.025,0.975))),3)
#  mean  2.5% 97.5% 
# 1.482 1.392 1.582 
round(c(mean=mean(PtoT),quantile(PtoT,c(0.025,0.975))),3)
#  mean  2.5% 97.5% 
# 3.026 2.767 3.303
##
## SAVE VECTORS
saveRDS(PtoM,"processed_data/PtoM.Rds")
saveRDS(PtoT,"processed_data/PtoT.Rds")
##----------------------------------------------