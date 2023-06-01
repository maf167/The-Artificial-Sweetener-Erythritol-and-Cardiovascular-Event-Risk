#
# Custom code used for the article: "The artificial sweetener erythritol and 
# cardiovascular event risk" (https://doi.org/10.1038/s41591-023-02223-9)
#  
# Note custom code was used to generate Table 1, Figure 1, and Figure 2
#
# Dependencies: tableone, rms, survminer, dplyr
#
# Note f1 is a data frame containing the fields referenced below.
#
# Author: Xinmin S Li, Marc Ferrell


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Table 1: Clinical characteristics of the discovery and validation cohorts
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#==============Table 1=============P20 
library(tableone)
myVars  <- c("AgeAtProc","MALE", "DIABETICS", "HxHtn","CurrentSmoker","DMS3",
             "ALLCAD", "HxCHF","newHxCHF", "HxMI", "CKD01")
catVars <- c("MALE", "DIABETICS", "HxHtn","CurrentSmoker","DMS3","ALLCAD", 
             "HxCHF","newHxCHF", "HxMI", "CKD01")
tab1    <- CreateTableOne(vars = myVars, data=f1,factorVars = catVars )
f2      <- print(tab1,  showAllLevels = TRUE, 
                 formatOptions = list(big.mark = ","))
write.table(f2,"1157 GB P20 all subjects-07-08-22.csv",sep=",",col.names=NA)

#=======================================Table 1 GB 2149
library(tableone)
myVars <- c("AgeAtProc","BMI","BPSystolic", "LDL.Priority", "TG.Priority",  
            "HDL.Priority","Chol..Priority", "TG.Priority" )

tab2   <- CreateTableOne(vars = myVars, data=f1, )
f3     <- print(tab2, nonnormal =  myVars,  showAllLevels = TRUE, 
                formatOptions = list(big.mark = ","))
write.table(f3,"1157 GB P20 all subjects-07-08-22_2.csv",sep=",",col.names=NA)
#=======================================

library(tableone)
myVars  <- c("AgeAtProc","MALE", "DIABETICS", "HxHtn","CurrentSmoker","DMS3",
             "ALLCAD", "HxCHF","newHxCHF", "HxMI", "CKD01")
catVars <- c("MALE", "DIABETICS", "HxHtn","CurrentSmoker","DMS3","ALLCAD", 
             "HxCHF","newHxCHF", "HxMI", "CKD01")

tab1    <- CreateTableOne(vars = myVars, data=f1,factorVars = catVars )
f2      <- print(tab1,  showAllLevels = TRUE,
                 formatOptions = list(big.mark = ","))
write.table(f2,"2149 GB all subjects-07-08-22.csv",sep=",",col.names=NA)

library(tableone)
myVars <- c("AgeAtProc","BMI","BPSystolic", "LDL.Priority", "TG.Priority",  
            "HDL.Priority","Chol..Priority", "TG.Priority" )

tab2   <- CreateTableOne(vars = myVars, data=f1, )
f3     <- print(tab2, nonnormal =  myVars,  showAllLevels = TRUE, 
                formatOptions = list(big.mark = ","))
write.table(f3,"2149 GB all subjects-07-08-22_2.csv",sep=",",col.names=NA)

#===================Table 1===========Berlin 833
library(tableone)
myVars   <- c("Age","Sex", "Diabetes","Hypertension","Smoker", "D_MI_StrokeYN",
              "CAD01","HF_22", "History.of.MI", "CKD01")
catVars  <- c("Sex", "Diabetes","Hypertension","Smoker", "D_MI_StrokeYN", 
              "CAD01","HF_22", "History.of.MI", "CKD01")
library(tableone)
tab1     <- CreateTableOne(vars = myVars, data=f2,factorVars = catVars )
f3       <- print(tab1,  showAllLevels = TRUE,
                  formatOptions = list(big.mark = ","))
write.table(f3,"833 Berlin all subjects-07-08-22.csv",sep=",",col.names=NA)
#========================================
library(tableone)
myVars <- c("Age","LDL..mg.dl.", "HDL..mg.dl.",  "Total.Chol..md.dl.", 
            "Triglyceride..mg.dl." )

tab2   <- CreateTableOne(vars = myVars, data=f2)
f4     <- print(tab2, nonnormal =  myVars,  showAllLevels = TRUE, 
                formatOptions = list(big.mark = ","))
write.table(f4,"833 Berlin all subjects-07-08-22_2.csv",sep=",",col.names=NA)
#==========================================================


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 1: Kaplan-Meier estimates and forest plots indicating the risks of 
#           MACE, according to erythritol quartile level.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#==============P20======================================KM with risk table
library(rms)
library(survminer)
x   <- f1$erythritol

q   <- quantile(x,c(0.25,0.50,0.75),na.rm=T)
df1 <- which(!is.na(x)&x<=q[1])
df2 <- which(!is.na(x)&x<=q[2]&x>q[1])
df3 <- which(!is.na(x)&x<=q[3]&x>q[2])
df4 <- which(!is.na(x)&x>q[3])

a   <- rep(NA,length(x))
a[df1] <- 1;a[df2]<-2;a[df3]<-3;a[df4]<-4
f1$OH  <- a

fit <- robcov(coxph(Surv(DTDMS3_YU/365,DMS3)~as.factor(OH),f1))
summary(fit)

km <- survfit(Surv(DTDMS3_YU/365,DMS3)~as.factor(OH),data=f1)

ggsurvplot(km, data=f1, pval=TRUE, size=1, legend.labs = c("Q1","Q2","Q3","Q4"), 
           palette =c("green","black","blue","red"), ylim=c(0.8,1),
           risk.table = TRUE, risk.table.col = "strata",risk.table.height =0.4 )

#=======P20=========Unadjust============ erythritol and DMS3 quartile analysis
library(rms)
x      <- f1$erythritol

q      <- quantile(x,c(0.25,0.50,0.75),na.rm=T)
df1    <- which(!is.na(x)&x<=q[1])
df2    <- which(!is.na(x)&x<=q[2]&x>q[1])
df3    <- which(!is.na(x)&x<=q[3]&x>q[2])
df4    <- which(!is.na(x)&x>q[3])

a      <- rep(NA,length(x))
a[df1] <- 1;a[df2] <- 2 ; a[df3] <-3 ; a[df4] <- 4
f1$OH<-a

fit<-robcov(coxph(Surv(DTDMS3_YU/365,DMS3)~as.factor(OH),f1))
summary(fit)

#=======Adjusted===============+AgeAtProc+BPSystolic+MALE+DIABETICS
#==============================+CurrentSmoker+HDL.Priority+LDL.Priority
#==============================+TG.Priority

library(rms)
x      <- f1$erythritol

q      <- quantile(x,c(0.25,0.50,0.75),na.rm=T)
df1    <- which(!is.na(x)&x<=q[1])
df2    <- which(!is.na(x)&x<=q[2]&x>q[1])
df3    <- which(!is.na(x)&x<=q[3]&x>q[2])
df4    <- which(!is.na(x)&x>q[3])

a      <- rep(NA,length(x))
a[df1] <-1;a[df2] <- 2;a[df3] <- 3;a[df4] <- 4
f1$OH  <-a

fit<-robcov(coxph(Surv(DTDMS3_YU/365,DMS3)~
                    as.factor(OH) + AgeAtProc + BPSystolic +MALE + DIABETICS +
                    CurrentSmoker + HDL.Priority + LDL.Priority + TG.Priority,
                  f1))
summary(fit)

#================================GB2149  KM with risk table
# KM with risk table
library(rms)
library(survminer)
x   <- f1$Erythritol.uM

q   <- quantile(x,c(0.25,0.50,0.75),na.rm=T)
df1 <- which(!is.na(x)&x<=q[1])
df2 <- which(!is.na(x)&x<=q[2]&x>q[1])
df3 <- which(!is.na(x)&x<=q[3]&x>q[2])
df4 <- which(!is.na(x)&x>q[3])

a   <- rep(NA,length(x))
a[df1] <- 1;a[df2] <- 2;a[df3] <- 3;a[df4] <- 4
f1$OH <- a

fit<-robcov(coxph(Surv(DTDMS3_YU/365,DMS3)~ as.factor(OH), f1))
summary(fit)

km<-survfit(Surv(DTDMS3_YU/365,DMS3)~as.factor(OH), data=f1)

ggsurvplot(km, data=f1, pval=TRUE, size=1, legend.labs = c("Q1","Q2","Q3","Q4"),
           palette =c("green","black","blue","red"), ylim=c(0.8,1),
           risk.table = TRUE, risk.table.col = "strata",risk.table.height =0.4)

#==========Unadjusted model: Erythritol.uM and DMS3 quartile analysis
library(rms) 
x      <- f1$Erythritol.uM 

q      <- quantile(x,c(0.25,0.50,0.75),na.rm=T)
df1    <- which(!is.na(x)&x<=q[1])
df2    <- which(!is.na(x)&x<=q[2]&x>q[1])
df3    <- which(!is.na(x)&x<=q[3]&x>q[2])
df4    <- which(!is.na(x)&x>q[3])

a      <- rep(NA,length(x))
a[df1] <-1;a[df2] <- 2;a[df3] <- 3;a[df4] <- 4
f1$OH  <- a

fit    <- robcov(coxph(Surv(DTDMS3_YU/365,DMS3)~as.factor(OH),f1))
summary(fit)

#==============Adjusted model: AgeAtProc+BPSystolic+MALE+DIABETICS+
#==============                CurrentSmoker+HDL.Priority+LDL.Priority+
#==============                TG.Priority

library(rms)
x      <- f1$Erythritol.uM

q      <- quantile(x,c(0.25,0.50,0.75),na.rm=T)
df1    <- which(!is.na(x)&x<=q[1])
df2    <- which(!is.na(x)&x<=q[2]&x>q[1])
df3    <- which(!is.na(x)&x<=q[3]&x>q[2])
df4    <- which(!is.na(x)&x>q[3])

a      <- rep(NA,length(x))
a[df1] <-1;a[df2] <- 2;a[df3] <- 3;a[df4] <- 4
f1$OH  <- a

fit<-robcov(coxph(Surv(DTDMS3_YU/365,DMS3)~as.factor(OH) + AgeAtProc + 
                    BPSystolic + MALE + DIABETICS + CurrentSmoker + 
                    HDL.Priority + LDL.Priority + TG.Priority, f1))
summary(fit)

#=============Berlin Cohort 833  =================KM with risk table===========
library(rms)
library(survminer)
x         <- f2$Erythritol..然.

q         <- quantile(x,c(0.25,0.50,0.75),na.rm=T)
index1    <- which(!is.na(x)&x<=q[1])
index2    <- which(!is.na(x)&x<=q[2]&x>q[1])
index3    <- which(!is.na(x)&x<=q[3]&x>q[2])
index4    <- which(!is.na(x)&x>q[3])

a         <- rep(NA,length(x))
a[index1] <-1;a[index2] <- 2;a[index3] <- 3;a[index4] <- 4
f2$OH     <- a

fit       <- robcov(coxph(Surv(D_MI_Stroke/365,D_MI_StrokeYN)~as.factor(OH),f2))
summary(fit)
km<-survfit(Surv(D_MI_Stroke/365,D_MI_StrokeYN)~as.factor(OH),data=f2)
ggsurvplot(km, data=f2, pval=TRUE, size=1, legend.labs = c("Q1","Q2","Q3","Q4"), 
           palette =c("green","black","blue","red"), ylim=c(0.5,1),
           risk.table = TRUE,risk.table.col = "strata",risk.table.height =0.4 ) 

#======Unadjusted model============= Erythritol..然. and DMS3 quartile analysis
library(rms)
x         <- f2$Erythritol..然.

q         <- quantile(x,c(0.25,0.50,0.75),na.rm=T)
index1    <- which(!is.na(x)&x<=q[1])
index2    <- which(!is.na(x)&x<=q[2]&x>q[1])
index3    <- which(!is.na(x)&x<=q[3]&x>q[2])
index4    <- which(!is.na(x)&x>q[3])

a         <- rep(NA,length(x))
a[index1] <-1;a[index2] <- 2;a[index3] <- 3;a[index4] <- 4
f2$OH     <- a

fit       <- robcov(coxph(Surv(D_MI_Stroke/365,D_MI_StrokeYN)~as.factor(OH),f2))
summary(fit)

#========Adjusted model=========+Age+Sex+Smoker+Hypertension+Diabetes
#===============================+LDL..mg.dl.+HDL..mg.dl.+Triglyceride..mg.dl.

library(rms)
x         <- f2$Erythritol..然.

q         <- quantile(x,c(0.25,0.50,0.75),na.rm=T)
index1    <- which(!is.na(x)&x<=q[1])
index2    <- which(!is.na(x)&x<=q[2]&x>q[1])
index3    <- which(!is.na(x)&x<=q[3]&x>q[2])
index4    <- which(!is.na(x)&x>q[3])

a         <- rep(NA,length(x))
a[index1] <- 1;a[index2] <- 2;a[index3] <- 3;a[index4] <- 4
f2$OH     <- a

fit       <- robcov(coxph(Surv(D_MI_Stroke/365,D_MI_StrokeYN)~as.factor(OH) +
                            Age + Sex + Smoker + Hypertension + Diabetes + 
                            LDL..mg.dl. + HDL..mg.dl. + Triglyceride..mg.dl.,
                          f2))
summary(f2)
#==============================


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 2:  Long-term risk of MACE among patient subgroups.
# Sample
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#==============================Female vs Male

library(rms) 

f2     <- subset(f1, MALE==0)

x      <- df$Erythritol.uM

q      <- quantile(x,c(0.25,0.50,0.75),na.rm=T)
df1    <- which(!is.na(x)&x<=q[1])
df2    <- which(!is.na(x)&x<=q[2]&x>q[1])
df3    <- which(!is.na(x)&x<=q[3]&x>q[2])
df4    <- which(!is.na(x)&x>q[3])

a      <- rep(NA,length(x))
a[df1] <- 1;a[df2] <- 2;a[df3] <- 3;a[df4] <- 4
f2$OH  <- a

fit    <- robcov(coxph(Surv(DTDMS3_YU/365,DMS3)~as.factor(OH),f2))
summary(fit)

f3     <- subset(f1,MALE==1)

library(rms)
x      <- f3$Erythritol.uM

q      <- quantile(x,c(0.25,0.50,0.75),na.rm=T)
df1    <- which(!is.na(x)&x<=q[1])
df2    <- which(!is.na(x)&x<=q[2]&x>q[1])
df3    <- which(!is.na(x)&x<=q[3]&x>q[2])
df4    <- which(!is.na(x)&x>q[3])

a      <- rep(NA,length(x))
a[df1] <- 1;a[df2] <- 2;a[df3] <- 3;a[df4] <- 4
f3$OH  <- a

fit    <- robcov(coxph(Surv(DTDMS3_YU/365,DMS3)~as.factor(OH),f3))
summary(fit)

