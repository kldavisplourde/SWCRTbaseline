library(nlme)
library(doMC)
library(doRNG)
library(lmeInfo)

setwd("/Users/kdavis07/Documents/GitHub/CoPrimarySWCRT")
source("gendata_copri_varCluster_HH.R")
source("EM_uncorrected_BV_HH.R")

args<-commandArgs(trailingOnly = TRUE)
k<-as.integer(args[1])
if (is.na(k)) k <- 1
paste("Scenario:",k)

ncores<-as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK",1))
if (is.na(ncores)) ncores<-1
registerDoMC(cores=ncores)

# define scenarios
scenarios <- read.table("/Users/kdavis07/Dropbox/SW-CRT Methods Development/2_CoPrimary/RCode/Simulations/HH/Sim_Params.txt", header=TRUE, sep="")
scenarios <- subset(scenarios, scenario == k)

scenario <- k
t <- scenarios$t
N <- scenarios$N
cs <- scenarios$m
eff<-c(0,scenarios$delta2)
rho01<-matrix(c(scenarios$rho01.11,scenarios$rho01.12,scenarios$rho01.12,scenarios$rho01.22),2)
rho2<-matrix(c(1,scenarios$rho2.12,scenarios$rho2.12,1),2)
vars<-c(1,1) #c(scenarios$var1,scenarios$var2)
beta <- rep(0,t-1)
nsim<-1000

set.seed(8374+k)
fail_count <- 0
max_fail <- 200
  
simData <- NULL
# Loop Index
i<-0
# While Loop (discarding false data)
while(i<nsim){
  i<-i+1
  itemp<-i
  #set.seed(i+173)
  
  dt<-datagen_cont(n=N, m=cs, K=2, cv=0, rho01=rho01, rho2=rho2, vars=vars, eff=eff, time.eff=c(beta,beta))
  data<-dt$short

    lme1<-lme(out1~1,random=~1|cluster,data=data,control=lmeControl(returnObject=TRUE))
    lme2<-lme(out2~arm,random=~1|cluster,data=data,control=lmeControl(returnObject=TRUE))

  param<-try(EM.estim(data,lme1,lme2, maxiter=500, epsilon=1e-4, verbose=FALSE))
  if(class(param)=="try-error"){
    i<- i-1;fail_count <-fail_count+1
  }else{
    if(param$SEcheck=="ERROR"|anyNA(param$theta$zeta)==TRUE|anyNA(param$theta$SigmaE)==TRUE|anyNA(param$theta$SigmaPhi)==TRUE){i<- i-1;fail_count <-fail_count+1}
  }
  if(fail_count > max_fail){break}
  if(i<itemp){next}
  
  SEtheta.i<-sqrt(diag(param$Vtheta))
  results.i<- c(param$theta$zeta,c(param$theta$SigmaPhi[!lower.tri(param$theta$SigmaPhi)]),param$theta$SigmaE[!lower.tri(param$theta$SigmaE)],SEtheta.i)
  
  # combining results
  simData<-rbind(simData,results.i)
}
#stopCluster(makeCluster(ncores))

  colnames(simData)<-c("Intercept.est1",
                       "Intercept.est2","Treatment.est2",
                       "SigmaPhi11","SigmaPhi12","SigmaPhi22","SigmaE11","SigmaE12","SigmaE22",
                       "Intercept.se1",
                       "Intercept.se2","Treatment.se2",
                       "SigmaPhi11.se","SigmaPhi12.se","SigmaPhi22.se","SigmaE11.se","SigmaE12.se","SigmaE22.se")
  
simData <- as.data.frame(simData)

write.table(simData, file=paste("/Users/kdavis07/Library/CloudStorage/Dropbox/SW-CRT Methods Development/3_Baseline/RCode/z_CRT/UncorrectedResults_scenario",scenario,".txt",sep=""), sep="\t", row.names=F)

#write.table(simData, file=paste("results/CorrectedResults_",analysis,"_scenario",scenario,".txt",sep=""), sep="\t", row.names=F)