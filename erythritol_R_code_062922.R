# Generic function for Cox models. Designed to extract a specific 
## hazard ratio from models with and without other variables

## Authors: Marco Witkowski, MD, Lin Li, PhD, Marc Ferrell
## 29 June 2022

## Dependencies: rms, survival

library(rms)
library(survival) 

## Examples

# GeneBankHR(table, "metabolite", adjustments = c())

# KMquartiles(table, "metabolite", ylim = c(0,100), yticks = c(0,50,100))

GeneBankHR  <- function(GBtable, metabolite, adjustments = c(),
                       Y = "DMS3", DT = "DTDMS3_YU"){
  
  x <- GBtable[,metabolite]
  
  ## Recode metabolite concentrations as quartiles
  
  q<-quantile(x,c(0.25,0.50,0.75),na.rm=T)
  index1<-which(!is.na(x)&x<=q[1])
  index2<-which(!is.na(x)&x<=q[2]&x>q[1])
  index3<-which(!is.na(x)&x<=q[3]&x>q[2])          
  index4<-which(!is.na(x)&x>q[3])
  
  a<-rep(NA,length(x))
  a[index1]<-0;a[index2]<-1;a[index3]<-0;a[index4]<-0
  OH2<-a
  
  a<-rep(NA,length(x))
  a[index1]<-0;a[index2]<-0;a[index3]<-1;a[index4]<-0
  OH3<-a
  
  a<-rep(NA,length(x))
  a[index1]<-0;a[index2]<-0;a[index3]<-0;a[index4]<-1
  OH4<-a
  
  ## Make Unadjusted model
  
  coxph.un <- robcov(coxph(Surv(GBtable[,DT], GBtable[,Y]) ~ OH2+OH3+OH4))
  coxph.un.s <- summary(coxph.un)
  
  ## Make Adjusted model
  
  if(length(adjustments) > 0){
    my.cov <- list()
    for(colname in adjustments){
      t <- GBtable[, colname]
      if(!is.numeric(t)) t <- as.numeric(factor(t)) - 1
      my.cov[[colname]] <- t
    }
    my.cov <- do.call("cbind", my.cov)
    
    coxph.ad <- robcov(coxph(Surv(GBtable[,DT], GBtable[,Y]) ~ 
                               OH2+OH3+OH4 + my.cov))
    coxph.ad.s <- summary(coxph.ad)
  }
  
  
    
  # Construct output table
    
  HR = data.frame(Model = rep(c("Unadjusted", "Adjusted"), each = 3),
                  Quartile = rep(2:4, 2),
                  hr=numeric(2), lo95=numeric(2), hi95 = numeric(2), 
                    log_HR = numeric(2), log_se = numeric(2), P = numeric(2))
  # rownames(HR) = c("Unadjusted", "Adjusted")
  HR[1:3, "hr"] = exp(coef(coxph.un))[1:3]
  HR[1:3, "lo95"]   = exp(confint(coxph.un))[1:3,1]
  HR[1:3, "hi95"]   = exp(confint(coxph.un))[1:3,2]
  HR[1:3, "log_HR"] = coef(coxph.un)[1:3]
  HR[1:3, "log_se"] = coxph.un.s$coefficients[1:3,3]
  HR[1:3, "P"] =    coef(coxph.un.s)[1:3,5]
  
  if(length(adjustments) > 0){
    HR[4:6, "hr"]   = exp(coef(coxph.ad))[1:3]
    HR[4:6, "lo95"]     = exp(confint(coxph.ad))[1:3,1]
    HR[4:6, "hi95"]     = exp(confint(coxph.ad))[1:3,2]
    HR[4:6, "log_HR"] = coef(coxph.ad)[1:3]
    HR[4:6, "log_se"] = coxph.ad.s$coefficients[1:3,3]
    HR[4:6, "P"] =    coef(coxph.ad.s)[1:3,5]
  }
  
  (HR)
  
} 

KMquartiles <- function(GBtable, metabolite, ylim = c(80,100),
                        yticks = c(80, 90, 100),
                        Y = "DMS3", DT = "DTDMS3_YU"){
  
  x <- GBtable[,metabolite]
  
  ## Recode metabolite concentrations as quartiles
  
  q<-quantile(x,c(0.25,0.50,0.75),na.rm=T)
  index1<-which(!is.na(x)&x<=q[1])
  index2<-which(!is.na(x)&x<=q[2]&x>q[1])
  index3<-which(!is.na(x)&x<=q[3]&x>q[2])          
  index4<-which(!is.na(x)&x>q[3])
  
  a<-rep(NA,length(x))
  a[index1]<-1;a[index2]<-2;a[index3]<-3;a[index4]<-4
  OH<-a
  
  ## Plot KM
  
  km<-survfit(Surv(GBtable[,DT] / 365, GBtable[,Y]) ~ as.factor(OH), data=GBtable,
              type='fleming')
  color<-c("green", "black", "blue", "red")
  plot(km, ylab="", xlab="", ylim=0.01 * ylim, xlim=c(0,3.2), yaxt="n", xaxt="n",
       col=color, lwd=5, lty=1, main=paste(metabolite, Y))
  title(xlab=expression(bold("Years")),
        ylab=expression(bold("3-Year Event-free survival (%)")), cex.lab=1.2)
  a<- c(0, 1, 2, 3)
  b<- yticks
  axis(1, at=a, labels=a,     col.axis="black", las=0, cex.axis=1.5)
  axis(2, at=0.01*b, labels=b, col.axis="black", las=0, cex.axis=1.5)
  
} 
