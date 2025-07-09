################################################################################
### FIGURE 7, SSB MATURE FEMALE ABUNDANCE ESTIMATES WAVELET ANALYSIS
##------------------------------------------------------------------------------
### LOAD PACKAGES
library(biwavelet)
### ADDITIONAL FUNCTIONS
source("code/addTrans.R")  # COLOUR TRANSPARENCY FILTER
source("code/glbpf.R")
###-----------------------------------------------------------------------------
##  1   DATA
##  1.1 SSB MATURE FEMALE ABUNDANCE ESTIMATES WAVELET ANALYSIS
##   -> DATA GENERATED WITH SCRIPT: 
##      'NTSims7925.Rds'
NTSimsF <- readRDS("processed_data/NTSims7925.Rds")
#
##  1.2 TRANSFORM DATA FOR WAVELET ANALYSIS   
X <- as.vector(colMeans(NTSimsF))
X <- X-glbpf(X,9.4)                 #apply low pass Gaussian filter to eliminate very long-term waves
X <- cbind(1979:2025,scale(X))
#
##  1.3 ## SIGNY ISLAND PROCESSED DATA
SIAFStrends <- readRDS("processed_data/AFStrends.Rds")
Y <- as.vector(SIAFStrends$Afs[3:49])
Y <- Y - glbpf(Y,9.4)                 #apply low pass Gaussian filter to eliminate very long-term waves
Y <- cbind(1979:2025,scale(Y))
##
#   2 CALCULATE WAVELETS AND GENERATE PLOTS
pdf("documentation/Figure7.pdf",family="Helvetica",width=8.27,height=11.69)
par(omi = c(3.5,1,1,0.9))       ## OUTER MARGINS IN INCHES - c(bottom, left, top, right)
layout(matrix(1:3,ncol=1,byrow=T),respect=F)  # PLOT LAYOUT
par(mai = c(0.2,0.4,0.2,0.8))
#
par(xpd = F)
wt.AFS <- biwavelet::wt(X,sig.level=0.95,max.scale=9)
plot(wt.AFS,plot.cb=T,legend.horiz=F,ylim=c(2,8),yaxt="n",xaxt="n",xlab="",ylab="",cex.axis=0.75,lwd.sig=1.3,alpha.coi=0.7,xlim=c(1979,2025))
axis(1,seq(1980,2025,by=5),labels=T,tck=-0.025,col=1,cex.axis=0.95,padj=-1.0)
axis(1,seq(1979,2025,by=1),labels=F,tck=-0.0125,col=1,cex.axis=0.95) #,padj=-1.5)
axis(2,log2(2:9),labels=paste(2:9),las=2,tck=-0.025,cex.axis=0.9,col=1,hadj=-0.0)
mtext("Scale (years)",side=2,line=2,cex=.85)
box()
dev.off()



pdf("documentation/Figure7.pdf",family="Helvetica",width=8.27,height=11.69)
#
par(omi = c(4.25,1,1,0.9))       ## OUTER MARGINS IN INCHES - c(bottom, left, top, right)
layout(matrix(1:3,ncol=1,byrow=T),respect=F)  # PLOT LAYOUT
#
## 1st PLOT - SSB-BI AFS WAVELET
par(mai = c(0.00,0.3,0.0,0.9))
par(xpd = F)
wt.AFS <- biwavelet::wt(X,sig.level=0.95,max.scale=9)
plot(wt.AFS,plot.cb=T,legend.horiz=F,ylim=c(2,8),yaxt="n",xaxt="n",xlab="",ylab="",cex.axis=0.75,lwd.sig=1.3,alpha.coi=0.7,xlim=c(1979,2025))
axis(1,seq(1980,2025,by=5),labels=F,tck=-0.025,col=1,cex.axis=0.95,padj=-1.0)
axis(1,seq(1979,2025,by=1),labels=F,tck=-0.0125,col=1,cex.axis=0.95) #,padj=-1.5)
axis(2,log2(2:9),labels=F,las=2,tck=-0.025,cex.axis=0.9,col=1,hadj=-0.0)
axis(2,log2(2:7),labels=paste(2:7),las=2,tck=-0.025,cex.axis=0.9,col=1,hadj=-0.0)
# mtext("Scale (years)",side=2,line=2,cex=.85)
mtext(expression(paste("Scale (years) - ",scriptstyle(W[italic(BI)]),scriptstyle(group("(",list(f,tau),")")),sep="")),side=2,line=2,cex=.85)
mtext("a",side=3,adj=-0.08,padj=1.0,cex=.95,font=2)
box()
#
## 2nd PLOT - SI AFS WAVELET
par(mai = c(0.00,0.3,0.0,0.9))
par(xpd = F)
## WAVELET
wt.SIAFS <- biwavelet::wt(Y,sig.level=0.95,max.scale=9)
plot(wt.SIAFS,plot.cb=T,legend.horiz=F,ylim=c(2,8),yaxt="n",xaxt="n",xlab="",ylab="",cex.axis=0.75,lwd.sig=1.3,alpha.coi=0.7,xlim=c(1979,2025))
axis(3,seq(1980,2025,by=5),labels=F,tck=-0.025,col=1,cex.axis=0.95,padj=-1.0)
axis(3,seq(1979,2025,by=1),labels=F,tck=-0.0125,col=1,cex.axis=0.95) #,padj=-1.5)
axis(2,log2(2:9),labels=F,las=2,tck=-0.025,cex.axis=0.9,col=1,hadj=-0.0)
axis(2,log2(2:7),labels=paste(2:7),las=2,tck=-0.025,cex.axis=0.9,col=1,hadj=-0.0)
#mtext("Scale (years)",side=2,line=2,cex=.85)
mtext(expression(paste("Scale (years) - ",scriptstyle(W[italic(SI)]),scriptstyle(group("(",list(f,tau),")")),sep="")),side=2,line=2,cex=.85)
mtext("b",side=3,adj=-0.08,padj=1.5,cex=.95,font=2)
box()
##
## 3rd PLOT - BI_SSB x SI CROSS WAVELET
par(mai = c(0.02,0.3,0.0,0.9))
par(xpd = F)
## CROSS-WAVELET
xwt.AFS <- biwavelet::xwt(X,Y,sig.level=0.95)
plot(xwt.AFS,plot.cb=T,legend.horiz=F,plot.phase=T,ylim=c(2,8),yaxt="n",xaxt="n",xlab="",ylab="",cex.axis=0.75,lwd.sig=1.3,alpha.coi=0.7,
     xlim=c(1979,2025),arrow.cutoff=0.99,arrow.col=addTrans("black",100))
axis(1,seq(1980,2025,by=5),labels=T,tck=-0.025,col=1,cex.axis=0.95,padj=-1.0)
axis(1,seq(1979,2025,by=1),labels=F,tck=-0.0125,col=1,cex.axis=0.95) #,padj=-1.5)
axis(3,seq(1980,2025,by=5),labels=F,tck=-0.025,col=1,cex.axis=0.95,padj=-1.0)
axis(3,seq(1979,2025,by=1),labels=F,tck=-0.0125,col=1,cex.axis=0.95) #,padj=-1.5)
axis(2,log2(2:9),labels=paste(2:9),las=2,tck=-0.025,cex.axis=0.9,col=1,hadj=-0.0)
#mtext("Scale (years)",side=2,line=2,cex=.85)
mtext(expression(paste("Scale (years) - ",scriptstyle(W[italic(list(BI,SI))]),scriptstyle(group("(",list(f,tau),")")),sep="")),side=2,line=2,cex=.85)
mtext("c",side=3,adj=-0.08,padj=1.5,cex=.95,font=2)
box()
##
dev.off()








## CALCULATE COHERENCE WAVELET
## set.seed(666)
## wtc.AFS <- biwavelet::wtc(X, Y,sig.level=0.95,nrands=4999,max.scale=9)
## saveRDS(wtc.AFS,"processed_data/wtc.AFS.Rds")
wtc.AFS <- readRDS("processed_data/wtc.AFS.Rds")

## COHERENCE WAVELET
## CROSS-WAVELET
par(xpd = F)
plot(wtc.AFS,plot.cb=T,legend.horiz=F,plot.phase=T,ylim=c(2,8),yaxt="n",xaxt="n",xlab="",ylab="",cex.axis=0.75,lwd.sig=1.3,alpha.coi=0.7,
     xlim=c(1979,2025),arrow.cutoff=0.99,arrow.col=addTrans("black",100))
axis(1,seq(1980,2025,by=5),labels=T,tck=-0.025,col=1,cex.axis=0.95,padj=-1.0)
axis(1,seq(1979,2025,by=1),labels=F,tck=-0.0125,col=1,cex.axis=0.95) #,padj=-1.5)
axis(2,log2(2:9),labels=paste(2:9),las=2,tck=-0.025,cex.axis=0.9,col=1,hadj=-0.0)
mtext("Scale (years)",side=2,line=2,cex=.85)
box()

