#Figure 2 for Staffing levels 
library(readr)
library(sm)
library(vioplot)
library(tidyverse)
library(dplyr)
library(readxl)
X1000runs18bedfiniteall <- read_excel("Lofgren/Iterations/18bedfinite/1000runs18bedfiniteall.xlsx")
View(X1000runs18bedfiniteall)
data <- X1000runs18bedfiniteall
X1000runs18bedratiosinf <- read_excel("Lofgren/Iterations/18bedinf/1000runs18bedratiosinf.xlsx")
View(X1000runs18bedratiosinf)
data1 <-X1000runs18bedratiosinf

#Turning Finite Models into 1000 patient days
data$baseline_norm <- (data$Baseline/6570) * 1000
summary(data$baseline_norm)
data$OneOne_norm <- (data$OneOne/6570)* 1000
data$OneTwo_norm <- (data$OneTwo/6570)* 1000
data$OneSix_norm <- (data$OneSix/6570)* 1000
data$OneNine_norm <- (data$OneNine/6570)* 1000
data$TwoOne_norm <- (data$TwoOne/6570)* 1000
data$TwoTwo_norm <- (data$TwoTwo/6570)* 1000
data$TwoThree_norm <- (data$TwoThree/6570)* 1000
data$TwoSix_norm <- (data$TwoSix/6570)* 1000
data$TwoNine_norm <- (data$TwoNine/6570)* 1000
data$ThreeOne_norm <- (data$ThreeOne/6570)* 1000
data$ThreeTwo_norm <- (data$ThreeTwo/6570)* 1000
data$ThreeThree_norm <- (data$ThreeThree/6570)* 1000
data$ThreeSix_norm <- (data$ThreeSix/6570)* 1000
data$ThreeNine_norm <- (data$ThreeNine/6570)* 1000
#Turning InFinite Models into 1000 patient days
data1$baseline_norm <- (data1$Baseline/6570) * 1000
summary(data1$baseline_norm)
data1$OneOne_norm <- (data1$OneOne/6570)* 1000
data1$OneTwo_norm <- (data1$OneTwo/6570)* 1000
data1$OneSix_norm <- (data1$OneSix/6570)* 1000
data1$OneNine_norm <- (data1$OneNine/6570)* 1000
data1$TwoOne_norm <- (data1$TwoOne/6570)* 1000
data1$TwoTwo_norm <- (data1$TwoTwo/6570)* 1000
data1$TwoThree_norm <- (data1$TwoThree/6570)* 1000
data1$TwoSix_norm <- (data1$TwoSix/6570)* 1000
data1$TwoNine_norm <- (data1$TwoNine/6570)* 1000
data1$ThreeOne_norm <- (data1$ThreeOne/6570)* 1000
data1$ThreeTwo_norm <- (data1$ThreeTwo/6570)* 1000
data1$ThreeThree_norm <- (data1$ThreeThree/6570)* 1000
data1$ThreeSix_norm <- (data1$ThreeSix/6570)* 1000
data1$ThreeNine_norm <- (data1$ThreeNine/6570)* 1000
summary(data1$OneOne_norm)




par(mfrow=c(3,2))
infpy <- vioplot(data1$OneOne_norm,data1$OneTwo_norm, data1$baseline_norm,data1$OneSix_norm,data1$OneNine_norm, names=c("1", "2","3","6", "9"),
                 drawRect=FALSE,ylim=c(0,20), col="lightgray")
title(main="Infinte Models")
title(ylab="MRSA Acquisitions",cex.lab=1.35)
title(xlab="(a) Infinite Nursing Ratios with One Doctor",cex.lab=1.35)
#title(xlab="Nursing Ratios",cex.lab=1.35)
segments(0.75, mean(data1$OneOne_norm),1.25,mean(data1$OneOne_norm),col= "black",lwd=2)
segments(1.75, mean(data1$OneTwo_norm), 2.25, mean(data1$OneTwo_norm), col="black",lwd=2)
segments(2.75,mean(data1$baseline_norm),3.25,mean(data1$baseline_norm),col="black",lwd=2)
segments(3.75,mean(data1$OneSix_norm),4.25,mean(data1$OneSix_norm),col="black",lwd=2)
segments(4.75,mean(data1$OneNine_norm),5.25,mean(data1$OneNine_norm),col="black",lwd=2)
legend("topleft", c("Model Mean"), lwd=2,col=c("black"),lty=1, bty='n', cex=1.2)

vioplot(data$OneOne_norm, data$OneTwo_norm,data$baseline_norm,data$OneSix_norm,data$OneNine_norm, names=c("1","2","3","6", "9"),
        drawRect=FALSE,ylim=c(0,90),col="lightgray")
title(main="Finte Models")
title(ylab="MRSA Acquisitions",cex.lab=1.35)
title(xlab="(d) Finite Nursing Ratios with One Doctor",cex.lab=1.35)
segments(0.75, mean(data$OneOne_norm),1.25,mean(data$OneOne_norm),col= "black",lwd=2)
segments(1.75, mean(data$OneTwo_norm), 2.25, mean(data$OneTwo_norm), col="black",lwd=2)
segments(2.75,mean(data$baseline_norm),3.25,mean(data$baseline_norm),col="black",lwd=2)
segments(3.75,mean(data$OneSix_norm),4.25,mean(data$OneSix_norm),col="black",lwd=2)
segments(4.75,mean(data$OneNine_norm),5.25,mean(data$OneNine_norm),col="black",lwd=2)
legend("topleft", c("Model Mean"), lwd=2,col=c("black"),lty=1, bty='n', cex=1.2)

vioplot(data1$TwoOne_norm,data1$TwoTwo_norm,data1$TwoThree_norm,data1$TwoSix_norm,data1$TwoNine_norm, names=c("1", "2","3","6", "9"),
        drawRect=FALSE,ylim=c(0,20), col='lightgray')
title(ylab="MRSA Acquisitions",cex.lab=1.35)
title(xlab="(b) Infinite Nursing Ratios with Two Doctors",cex.lab=1.35)
segments(0.75, mean(data1$TwoOne_norm),1.25,mean(data1$TwoOne_norm),col= "black",lwd=2)
segments(1.75, mean(data1$TwoTwo_norm), 2.25, mean(data1$TwoTwo_norm), col="black",lwd=2)
segments(2.75,mean(data1$TwoThree_norm),3.25,mean(data1$TwoThree_norm),col="black",lwd=2)
segments(3.75,mean(data1$TwoSix_norm),4.25,mean(data1$TwoSix_norm),col="black",lwd=2)
segments(4.75,mean(data1$TwoNine_norm),5.25,mean(data1$TwoNine_norm),col="black",lwd=2)
legend("topleft", c("Model Mean"), lwd=2,col=c("black"),lty=1, bty='n', cex=1.2)

vioplot(data$TwoOne_norm, data$TwoTwo_norm,data$TwoThree_norm,data$TwoSix_norm,data$TwoNine_norm, names=c("1","2","3","6", "9"),
        drawRect=FALSE,ylim=c(0,90),col="lightgray")
title(ylab="MRSA Acquisitions",cex.lab=1.35)
title(xlab="(e) Finite Nursing Ratios with Two Doctors",cex.lab=1.35)
segments(0.75, mean(data$TwoOne_norm),1.25,mean(data$TwoOne_norm),col= "black",lwd=2)
segments(1.75, mean(data$TwoTwo_norm), 2.25, mean(data$TwoTwo_norm), col="black",lwd=2)
segments(2.75,mean(data$TwoThree_norm),3.25,mean(data$TwoThree_norm),col="black",lwd=2)
segments(3.75,mean(data$TwoSix_norm),4.25,mean(data$TwoSix_norm),col="black",lwd=2)
segments(4.75,mean(data$TwoNine_norm),5.25,mean(data$TwoNine_norm),col="black",lwd=2)
legend("topleft", c("Model Mean"), lwd=2,col=c("black"),lty=1, bty='n', cex=1.2)


vioplot(data1$ThreeOne_norm,data1$ThreeTwo_norm,data1$ThreeThree_norm,data1$ThreeSix_norm,data1$ThreeNine_norm, names=c("1", "2","3","6", "9"),
        drawRect=FALSE,ylim=c(0,20), col='lightgray')
title(main="")
title(ylab="MRSA Acquisitions",cex.lab=1.35)
title(xlab="(c) Infinite Nursing Ratios with Three Doctors",cex.lab=1.35)
segments(0.75, mean(data1$ThreeOne_norm),1.25,mean(data1$ThreeOne_norm),col= "black",lwd=2)
segments(1.75, mean(data1$ThreeTwo_norm), 2.25, mean(data1$ThreeTwo_norm), col="black",lwd=2)
segments(2.75,mean(data1$ThreeThree_norm),3.25,mean(data1$ThreeThree_norm),col="black",lwd=2)
segments(3.75,mean(data1$ThreeSix_norm),4.25,mean(data1$ThreeSix_norm),col="black",lwd=2)
segments(4.75,mean(data1$ThreeNine_norm),5.25,mean(data1$ThreeNine_norm),col="black",lwd=2)
legend("topleft", c("Model Mean"), lwd=2,col=c("black"),lty=1, bty='n', cex=1.2)

vioplot(data$ThreeOne_norm,data$ThreeTwo_norm,data$ThreeThree_norm,data$ThreeSix_norm,data$ThreeNine_norm, names=c("1", "2","3","6", "9"),
        drawRect=FALSE,ylim=c(0,90), col='lightgray')
title(main="")
title(ylab="MRSA Acquisitions",cex.lab=1.35)
title(xlab="(f) Finite Nursing Ratios with Three Doctors",cex.lab=1.35)
segments(0.75, mean(data$ThreeOne_norm),1.25,mean(data$ThreeOne_norm),col= "black",lwd=2)
segments(1.75, mean(data$ThreeTwo_norm), 2.25, mean(data$ThreeTwo_norm), col="black",lwd=2)
segments(2.75,mean(data$ThreeThree_norm),3.25,mean(data$ThreeThree_norm),col="black",lwd=2)
segments(3.75,mean(data$ThreeSix_norm),4.25,mean(data$ThreeSix_norm),col="black",lwd=2)
segments(4.75,mean(data$ThreeNine_norm),5.25,mean(data$ThreeNine_norm),col="black",lwd=2)
legend("topleft", c("Model Mean"), lwd=2,col=c("black"),lty=1, bty='n', cex=1.2)

#For Epidemics
par(mfrow=c(3,1))
vioplot(data$OneOne_norm, data$OneTwo_norm,data$baseline_norm,data$OneSix_norm,data$OneNine_norm, names=c("1","2","3","6", "9"),
        drawRect=FALSE,ylim=c(0,90),col="maroon")
title(main="Finite Models")
title(ylab="MRSA Acquisitions per 1000 Patient-Days",cex.lab=1.15)
title(xlab="Finite Nursing Ratios with One Doctor",cex.lab=1.35)
segments(0.75, mean(data$OneOne_norm),1.25,mean(data$OneOne_norm),col= "black",lwd=2)
segments(1.75, mean(data$OneTwo_norm), 2.25, mean(data$OneTwo_norm), col="black",lwd=2)
segments(2.75,mean(data$baseline_norm),3.25,mean(data$baseline_norm),col="black",lwd=2)
segments(3.75,mean(data$OneSix_norm),4.25,mean(data$OneSix_norm),col="black",lwd=2)
segments(4.75,mean(data$OneNine_norm),5.25,mean(data$OneNine_norm),col="black",lwd=2)
legend("topleft", c("Model Mean"), lwd=2,col=c("black"),lty=1, bty='n', cex=1.2)

vioplot(data$TwoOne_norm, data$TwoTwo_norm,data$TwoThree_norm,data$TwoSix_norm,data$TwoNine_norm, names=c("1","2","3","6", "9"),
        drawRect=FALSE,ylim=c(0,90),col="lightgray")
title(ylab="MRSA Acquisitions per 1000 Patient-Days",cex.lab=1.15)
title(xlab="Finite Nursing Ratios with Two Doctors",cex.lab=1.35)
segments(0.75, mean(data$TwoOne_norm),1.25,mean(data$TwoOne_norm),col= "black",lwd=2)
segments(1.75, mean(data$TwoTwo_norm), 2.25, mean(data$TwoTwo_norm), col="black",lwd=2)
segments(2.75,mean(data$TwoThree_norm),3.25,mean(data$TwoThree_norm),col="black",lwd=2)
segments(3.75,mean(data$TwoSix_norm),4.25,mean(data$TwoSix_norm),col="black",lwd=2)
segments(4.75,mean(data$TwoNine_norm),5.25,mean(data$TwoNine_norm),col="black",lwd=2)
legend("topleft", c("Model Mean"), lwd=2,col=c("black"),lty=1, bty='n', cex=1.2)

vioplot(data$ThreeOne_norm,data$ThreeTwo_norm,data$ThreeThree_norm,data$ThreeSix_norm,data$ThreeNine_norm, names=c("1", "2","3","6", "9"),
        drawRect=FALSE,ylim=c(0,90), col='darkgray')
title(main="")
title(ylab="MRSA Acquisitions per 1000 Patient-Days",cex.lab=1.15)
title(xlab="Finite Nursing Ratios with Three Doctors",cex.lab=1.35)
segments(0.75, mean(data$ThreeOne_norm),1.25,mean(data$ThreeOne_norm),col= "black",lwd=2)
segments(1.75, mean(data$ThreeTwo_norm), 2.25, mean(data$ThreeTwo_norm), col="black",lwd=2)
segments(2.75,mean(data$ThreeThree_norm),3.25,mean(data$ThreeThree_norm),col="black",lwd=2)
segments(3.75,mean(data$ThreeSix_norm),4.25,mean(data$ThreeSix_norm),col="black",lwd=2)
segments(4.75,mean(data$ThreeNine_norm),5.25,mean(data$ThreeNine_norm),col="black",lwd=2)
legend("topleft", c("Model Mean"), lwd=2,col=c("black"),lty=1, bty='n', cex=1.2)