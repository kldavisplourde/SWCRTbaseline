library(lme4)
library(doMC)
library(doRNG)

#library("lmeInfo", lib.loc="my_R_libs/")    #Need to uncomment?

#source("~/project/proj3_baseline/SimulationStudy/gendata_copri_varCluster_HoopGir.R")
source("gendata_copri_varCluster_HoopGir.R")
source("EM_uncorrected_BV.R")

args<-commandArgs(trailingOnly = TRUE)
k<-as.integer(args[1])
if (is.na(k)) k <- 1
paste("Scenario:",k)

ncores<-as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK",1))
if (is.na(ncores)) ncores<-1
registerDoMC(cores=ncores)

# define scenarios
scenarios <- read.table("Sim_Params.txt", header=TRUE, sep="")
scenarios <- subset(scenarios, scenario == k)

scenario <- k
t <- scenarios$t
N <- scenarios$N
cs <- scenarios$m
eff<-c(0,0) #c(0,scenarios$delta) #Data generation code needs trt effect for each outcome, therefore, setting baseline effect to 0
rho02<-matrix(c(scenarios$rho0,scenarios$rho0.01,scenarios$rho0.01,scenarios$rho0),2)
rho01<-matrix(c(scenarios$rho1,scenarios$rho1.01,scenarios$rho1.01,scenarios$rho1),2)
rho2<-matrix(c(1,scenarios$rho2.01,scenarios$rho2.01,1),2)
vars<-c(scenarios$total.var,scenarios$total.var)
bs <- 0
beta <- cumsum(c(0.1,0.1*0.5,0.1*(0.5^2),0.1*(0.5^3),0.1*(0.5^4),0.1*(0.5^5)))[1:t-1]
nsim<-1000

#set.seed(3678+k) #First set of 1000
set.seed(9461+k)   #Second set of 1000
fail_count <- 0
max_fail <- 200

simData.HG <- simData.CS <- simData.AC <- simData.ACps <- simData.ACcps <- simData.BV <- NULL
# Loop Index
i<-0
# While Loop (discarding false data)
while(i<nsim){
  i<-i+1
  itemp<-i
  
  dt<-datagen_cont(n=N, m=cs, K=2, cv=0, rho01=rho01, rho02=rho02, rho2=rho2, vars=vars, eff=eff, time.eff=c(beta,beta))
  data<-dt$short
  data$out2minusout1 <- data$out2-data$out1
  mean.baseline<-mean(data$out1) 
  ps.means <- aggregate(data$out1, list(data$time), FUN=mean)
  colnames(ps.means)<-c("time","ps.mean")
  new<-merge(data,ps.means,by="time")
  cps.means <- aggregate(data$out1, list(data$cluster.period), FUN=mean)
  colnames(cps.means)<-c("cluster.period","cps.mean")
  data<-merge(new,cps.means,by="cluster.period")
  data$out1.centered.ps<-data$out1-data$ps.mean
  data$out1.centered.cps<-data$out1-data$cps.mean

  if(t==3){
    lme1<-lmer(out1~time.1+time.2 +(1|cluster) +(1|cluster.period),data=data)
    lme2<-try(lmer(out2~time.1+time.2+arm +(1|cluster) +(1|cluster.period),data=data))
    
    lme.cs<-try(lmer(out2minusout1~time.1+time.2+arm +(1|cluster) +(1|cluster.period),data=data))
    
    lme.ancova<-try(lmer(out2~time.1+time.2+out1+arm +(1|cluster) +(1|cluster.period),data=data))
    lme.ancova.ps<-try(lmer(out2~time.1+time.2+out1.centered.ps+arm +(1|cluster) +(1|cluster.period),data=data))
    lme.ancova.cps<-try(lmer(out2~time.1+time.2+out1.centered.cps+arm +(1|cluster) +(1|cluster.period),data=data))
  }
  
  if(t==4){
    lme1<-lmer(out1~time.1+time.2+time.3 +(1|cluster) +(1|cluster.period),data=data)
    lme2<-try(lmer(out2~time.1+time.2+time.3+arm +(1|cluster) +(1|cluster.period),data=data))
    
    lme.cs<-try(lmer(out2minusout1~time.1+time.2+time.3+arm +(1|cluster) +(1|cluster.period),data=data))
    
    lme.ancova<-try(lmer(out2~time.1+time.2+time.3+out1+arm +(1|cluster) +(1|cluster.period),data=data))
    lme.ancova.ps<-try(lmer(out2~time.1+time.2+time.3+out1.centered.ps+arm +(1|cluster) +(1|cluster.period),data=data))
    lme.ancova.cps<-try(lmer(out2~time.1+time.2+time.3+out1.centered.cps+arm +(1|cluster) +(1|cluster.period),data=data))
  }
  
  if(t==5){
    lme1<-lmer(out1~time.1+time.2+time.3+time.4 +(1|cluster) +(1|cluster.period),data=data)
    lme2<-try(lmer(out2~time.1+time.2+time.3+time.4+arm +(1|cluster) +(1|cluster.period),data=data))
    
    lme.cs<-try(lmer(out2minusout1~time.1+time.2+time.3+time.4+arm +(1|cluster) +(1|cluster.period),data=data))
    
    lme.ancova<-try(lmer(out2~time.1+time.2+time.3+time.4+out1+arm +(1|cluster) +(1|cluster.period),data=data))
    lme.ancova.ps<-try(lmer(out2~time.1+time.2+time.3+time.4+out1.centered.ps+arm +(1|cluster) +(1|cluster.period),data=data))
    lme.ancova.cps<-try(lmer(out2~time.1+time.2+time.3+time.4+out1.centered.cps+arm +(1|cluster) +(1|cluster.period),data=data))
  }
  

  #Bivariate Model
  param.BV<-try(EM.estim(data,lme1,lme2,cluster="cluster",cluster.period="cluster.period",maxiter=500, epsilon=1e-4, verbose=FALSE))
  if(class(param.BV)=="try-error"|class(lme2)=="try-error"|class(lme.cs)=="try-error"|class(lme.ancova)=="try-error"|class(lme.ancova.ps)=="try-error"|class(lme.ancova.cps)=="try-error"){
    i<- i-1;fail_count <-fail_count+1
  }else{
    if(param.BV$SEcheck=="ERROR"|anyNA(param.BV$theta$zeta)==TRUE|anyNA(param.BV$theta$SigmaE)==TRUE|anyNA(param.BV$theta$SigmaPhi)==TRUE|anyNA(param.BV$theta$SigmaPsi)==TRUE){i<- i-1;fail_count <-fail_count+1}
  }
  if(fail_count > max_fail){break}
  if(i<itemp){next}
  
  
  zeta.HG <- as.numeric(fixef(lme2))
  vc2<-as.data.frame(VarCorr(lme2))
  SigmaPhi.HG <- vc2[vc2$grp=="cluster",4]
  SigmaPsi.HG <- vc2[vc2$grp=="cluster.period",4]
  SigmaE.HG <- vc2[vc2$grp=="Residual",4]
  
  
  zeta.CS <- as.numeric(fixef(lme.cs))
  vc2.CS<-as.data.frame(VarCorr(lme.cs))
  SigmaPhi.CS <- vc2.CS[vc2.CS$grp=="cluster",4]
  SigmaPsi.CS <- vc2.CS[vc2.CS$grp=="cluster.period",4]
  SigmaE.CS <- vc2.CS[vc2.CS$grp=="Residual",4]
  
  
  zeta.ancova <- as.numeric(fixef(lme.ancova))
  vc2.ancova<-as.data.frame(VarCorr(lme.ancova))
  SigmaPhi.ancova <- vc2.ancova[vc2.ancova$grp=="cluster",4]
  SigmaPsi.ancova <- vc2.ancova[vc2.ancova$grp=="cluster.period",4]
  SigmaE.ancova <- vc2.ancova[vc2.ancova$grp=="Residual",4]
  
  
  zeta.ancova.ps <- as.numeric(fixef(lme.ancova.ps))
  vc2.ancova.ps<-as.data.frame(VarCorr(lme.ancova.ps))
  SigmaPhi.ancova.ps <- vc2.ancova.ps[vc2.ancova.ps$grp=="cluster",4]
  SigmaPsi.ancova.ps <- vc2.ancova.ps[vc2.ancova.ps$grp=="cluster.period",4]
  SigmaE.ancova.ps <- vc2.ancova.ps[vc2.ancova.ps$grp=="Residual",4]
  
  
  zeta.ancova.cps <- as.numeric(fixef(lme.ancova.cps))
  vc2.ancova.cps<-as.data.frame(VarCorr(lme.ancova.cps))
  SigmaPhi.ancova.cps <- vc2.ancova.cps[vc2.ancova.cps$grp=="cluster",4]
  SigmaPsi.ancova.cps <- vc2.ancova.cps[vc2.ancova.cps$grp=="cluster.period",4]
  SigmaE.ancova.cps <- vc2.ancova.cps[vc2.ancova.cps$grp=="Residual",4]
  
  
  zetaSE.HG <- as.numeric(c(sqrt(diag(vcov(lme2))))) # add SE for random effects? No, Wald standard errors are very poor estimates anyway. Not sure I can extract anything other than Wald...
  results.HG.i<-c(zeta.HG,SigmaPhi.HG,SigmaPsi.HG,SigmaE.HG,zetaSE.HG)
  
  zetaSE.CS <- as.numeric(c(sqrt(diag(vcov(lme.cs))))) # add SE for random effects? No, Wald standard errors are very poor estimates anyway. Not sure I can extract anything other than Wald...
  results.CS.i<-c(zeta.CS,SigmaPhi.CS,SigmaPsi.CS,SigmaE.CS,zetaSE.CS)
  
  zetaSE.ancova <- as.numeric(c(sqrt(diag(vcov(lme.ancova))))) # add SE for random effects? No, Wald standard errors are very poor estimates anyway. Not sure I can extract anything other than Wald...
  results.AC.i<-c(zeta.ancova,SigmaPhi.ancova,SigmaPsi.ancova,SigmaE.ancova,zetaSE.ancova)
  
  zetaSE.ancova.ps <- as.numeric(c(sqrt(diag(vcov(lme.ancova.ps))))) # add SE for random effects? No, Wald standard errors are very poor estimates anyway. Not sure I can extract anything other than Wald...
  results.ACps.i<-c(zeta.ancova.ps,SigmaPhi.ancova.ps,SigmaPsi.ancova.ps,SigmaE.ancova.ps,zetaSE.ancova.ps)
  
  zetaSE.ancova.cps <- as.numeric(c(sqrt(diag(vcov(lme.ancova.cps))))) # add SE for random effects? No, Wald standard errors are very poor estimates anyway. Not sure I can extract anything other than Wald...
  results.ACcps.i<-c(zeta.ancova.cps,SigmaPhi.ancova.cps,SigmaPsi.ancova.cps,SigmaE.ancova.cps,zetaSE.ancova.cps)
  
  SEtheta.BV.i<-sqrt(diag(param.BV$Vtheta))
  results.BV.i<- c(param.BV$theta$zeta,param.BV$theta$SigmaPhi[!lower.tri(param.BV$theta$SigmaPhi)],param.BV$theta$SigmaPsi[!lower.tri(param.BV$theta$SigmaPsi)],param.BV$theta$SigmaE[!lower.tri(param.BV$theta$SigmaE)],SEtheta.BV.i,fail_count)
  
  
  # combining results
  simData.HG<-rbind(simData.HG,results.HG.i)
  simData.CS<-rbind(simData.CS,results.CS.i)
  simData.AC<-rbind(simData.AC,results.AC.i)
  simData.ACps<-rbind(simData.ACps,results.ACps.i)
  simData.ACcps<-rbind(simData.ACcps,results.ACcps.i)
  simData.BV<-rbind(simData.BV,results.BV.i)
}

if(t==3){
  colnames(simData.HG)<-colnames(simData.CS)<-c("Intercept.est","Period2.est","Period3.est","Treatment.est",
                          "SigmaPhi22","SigmaPsi22","SigmaE22",
                          "Intercept.se","Period2.se","Period3.se","Treatment.se")
  
  colnames(simData.AC)<-colnames(simData.ACps)<-colnames(simData.ACcps)<-c("Intercept.est","Period2.est","Period3.est","out1.est","Treatment.est",
                          "SigmaPhi22","SigmaPsi22","SigmaE22",
                          "Intercept.se","Period2.se","Period3.se","out1.se","Treatment.se")
  
  colnames(simData.BV)<-c("Intercept.est1","Period2.est1","Period3.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Treatment.est",
                       "SigmaPhi11","SigmaPhi12","SigmaPhi22",
                       "SigmaPsi11","SigmaPsi12","SigmaPsi22","SigmaE11","SigmaE12","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Treatment.se",
                       "SigmaPhi11.se","SigmaPhi12.se","SigmaPhi22.se",
                       "SigmaPsi11.se","SigmaPsi12.se","SigmaPsi22.se","SigmaE11.se","SigmaE12.se","SigmaE22.se","fail.count")
}

if(t==4){
  colnames(simData.HG)<-colnames(simData.CS)<-c("Intercept.est","Period2.est","Period3.est","Period4.est","Treatment.est",
                          "SigmaPhi22","SigmaPsi22","SigmaE22",
                          "Intercept.se","Period2.se","Period3.se","Period4.se","Treatment.se")
  
  colnames(simData.AC)<-colnames(simData.ACps)<-colnames(simData.ACcps)<-c("Intercept.est","Period2.est","Period3.est","Period4.est","out1.est","Treatment.est",
                                                                           "SigmaPhi22","SigmaPsi22","SigmaE22",
                                                                           "Intercept.se","Period2.se","Period3.se","Period4.se","out1.se","Treatment.se")
  
  colnames(simData.BV)<-c("Intercept.est1","Period2.est1","Period3.est1","Period4.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Period4.est2","Treatment.est",
                       "SigmaPhi11","SigmaPhi12","SigmaPhi22",
                       "SigmaPsi11","SigmaPsi12","SigmaPsi22","SigmaE11","SigmaE12","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Period4.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Period4.se2","Treatment.se",
                       "SigmaPhi11.se","SigmaPhi12.se","SigmaPhi22.se",
                       "SigmaPsi11.se","SigmaPsi12.se","SigmaPsi22.se","SigmaE11.se","SigmaE12.se","SigmaE22.se","fail.count")
}

if(t==5){
  colnames(simData.HG)<-colnames(simData.CS)<-c("Intercept.est","Period2.est","Period3.est","Period4.est","Period5.est","Treatment.est",
                          "SigmaPhi22","SigmaPsi22","SigmaE22",
                          "Intercept.se","Period2.se","Period3.se","Period4.se","Period5.se","Treatment.se")
  
  colnames(simData.AC)<-colnames(simData.ACps)<-colnames(simData.ACcps)<-c("Intercept.est","Period2.est","Period3.est","Period4.est","Period5.est","out1.est","Treatment.est",
                                                                           "SigmaPhi22","SigmaPsi22","SigmaE22",
                                                                           "Intercept.se","Period2.se","Period3.se","Period4.se","Period5.se","out1.se","Treatment.se")
  
  colnames(simData.BV)<-c("Intercept.est1","Period2.est1","Period3.est1","Period4.est1","Period5.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Period4.est2","Period5.est2","Treatment.est",
                       "SigmaPhi11","SigmaPhi12","SigmaPhi22",
                       "SigmaPsi11","SigmaPsi12","SigmaPsi22","SigmaE11","SigmaE12","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Period4.se1","Period5.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Period4.se2","Period5.se2","Treatment.se",
                       "SigmaPhi11.se","SigmaPhi12.se","SigmaPhi22.se",
                       "SigmaPsi11.se","SigmaPsi12.se","SigmaPsi22.se","SigmaE11.se","SigmaE12.se","SigmaE22.se","fail.count")
}




simData.HG <- as.data.frame(simData.HG)
simData.CS <- as.data.frame(simData.CS)
simData.AC <- as.data.frame(simData.AC)
simData.ACps <- as.data.frame(simData.ACps)
simData.ACcps <- as.data.frame(simData.ACcps)
simData.BV <- as.data.frame(simData.BV)

write.table(simData.HG, file=paste("results_error/UncorrectedResults_scenario",scenario,"_HG2.txt",sep=""), sep="\t", row.names=F)
write.table(simData.CS, file=paste("results_error/UncorrectedResults_scenario",scenario,"_CS2.txt",sep=""), sep="\t", row.names=F)
write.table(simData.AC, file=paste("results_error/UncorrectedResults_scenario",scenario,"_AC2.txt",sep=""), sep="\t", row.names=F)
write.table(simData.ACps, file=paste("results_error/UncorrectedResults_scenario",scenario,"_ACps2.txt",sep=""), sep="\t", row.names=F)
write.table(simData.ACcps, file=paste("results_error/UncorrectedResults_scenario",scenario,"_ACcps2.txt",sep=""), sep="\t", row.names=F)
write.table(simData.BV, file=paste("results_error/UncorrectedResults_scenario",scenario,"_BV2.txt",sep=""), sep="\t", row.names=F)