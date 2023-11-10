# Code for running bivariate model under a CRT setting - model with treatment effect on endpoint outcome measure only

# Will need to install the following packages (if don't already have)
library(nlme)
library(mvtnorm)
library(numDeriv)
library(sjlabelled)

# Need to source the EM algorithm for fitting the bivariate model
source("EM_uncorrected_BV_HH.R")

# Example dataset - Taken from Yang et al paper
KDD <- read.csv("/Users/Primary outcome_K-DPP trial.csv")
KDD$arm <- ifelse(KDD$arms0=="Control",0,1)  #Created a numerical treatment indicator

# Can't have one of the outcomes missing: here I impute missing with a simple mean
mean.sbp2<-mean(KDD$sysbp2,na.rm=TRUE)
KDD$sysbp2<-ifelse(is.na(KDD$sysbp2)=="TRUE",mean.sbp2,KDD$sysbp2)

# Baseline and endline models of interest
lme.baseline<-lme(sysbp0~1,random=~1|cluster0,data=KDD,control=lmeControl(returnObject=TRUE))
lme.endline<-lme(sysbp2~arm,random=~1|cluster0,data=KDD,control=lmeControl(returnObject=TRUE))

# Fit the bivariate model
param <- EM.estim(data = KDD, fm1=lme.baseline, fm2=lme.endline, maxiter=500, epsilon=1e-4, verbose=FALSE)