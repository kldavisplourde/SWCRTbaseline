# Example Application of calPower.R using TSOS Study

source("/Users/kdavis07/Documents/GitHub/SWCRTbaseline/calPower.R")

#Inputs
N <- 24     #Number of clusters
t <- 5      #Number of periods
m <- 8      #Number of observations within a cluster-period

#Correlations
rho2.01 <- 0.7            #Autocorrelation
rho0.0 <- rho0.1 <- 0.02  #Assuming equal baseline- and endline-specific ICCs 
rho0.01 <- 0.01           #Assuming between measure ICC is half of measure-specific ICC         
rho1.0 <- rho1.1 <- 0.01  #Assuming CAC = 0.5
rho1.01 <- 0.005          #Assuming between measure ICC is half of measure-specific ICC
vars<-c(15^2,15^2)        #Total variance for baseline and endline outcomes - from protocol

#Generate ICC matrices
rho0<-matrix(c(rho0.0,rho0.01,rho0.01,rho0.1),2)
rho1<-matrix(c(rho1.0,rho1.01,rho1.01,rho1.1),2)
rho2<-matrix(c(1,rho2.01,rho2.01,1),2)

#Generate power results

#ANCOVA with 85% power
delta <- 0.23*sqrt(vars[1])          #Intervention effect
calPower(delta=delta,vars=vars,rho0=rho0,rho1=rho1,rho2=rho2,N=N,t=t,m=m,alpha=0.05,user.spec.df=N-2)