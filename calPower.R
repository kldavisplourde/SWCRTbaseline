library(mvtnorm)
########################################################################################################################################################
##Power/Sample Size Calculation for incorporating a baseline measurement.##
# INPUT
# delta:   treatment effect
# vars:    A vector of marginal variances for baseline and endpoint outcomes
# rho0: a K by K dimensional matrix for the correlation parameters (rho0^k) and (rho0^kk')
# For rho0: Outcome-specific ICCs (within-period)
#           the diagonal elements correspond to rho0^k's 
#           the off-diagonal elements correspond to (rho0^kk')'s 
#           For example, rho0[1,1] corresponds to rho0^1, which is the within-period ICC for the baseline outcome
#                        rho0[1,2] corresponds to rho0^12, which is the within-period correlation of outcomes between subjects on the baseline and endpoint measures    
# rho1: a K by K dimensional matrix for the correlation parameters (rho1^k) and (rho1^kk')
# For rho1: Outcome-specific ICCs (between-period)
#           the diagonal elements correspond to rho1^k's 
#           the off-diagonal elements correspond to (rho1^kk')'s 
#           For example, rho1[1,1] corresponds to rho1^1, which is the between-period ICC for the baseline outcome
#                        rho1[1,2] corresponds to rho1^12, which is the between-period correlation of outcomes between subjects on the baseline and endpoint measures
# rho2: a K by K dimensional matrix for the correlation parameters (rho2^kk')
# For rho2: Intra-subject ICC
#           the diagonal elements are 1
#           the off-diagonal elements correspond to (rho2^kk')'s
#           For example, rho2[1,2] corresponds to rho2^12, which is the correlation of outcomes within same subject on the baseline and endpoint measures
# N: number of clusters
# t: number of time periods
# m: cluster-period size
# alpha: type I error control, default is 0.05
# df: degrees of freedom for JM, CS, and ANCOVA models (Default is N-3), for Hooper/Girling (follow-up) model this is set to N-2
########################################################################################################################################################
####Function to Calculate Power Given Design Configurations based on the t test#######
####Critical value c is set to t_alpha, (1-alpha)th quantile of the t distribution with df = N-2###
calPower <- function(delta,vars,rho0,rho1,rho2,N,t,m,alpha=0.05,df=0)
{
  K<-2
  
  # Create X matrix
  X<-NULL
  trtSeq<-matrix(0,t-1,t)
  trtSeq[upper.tri(trtSeq)]<-1
  g<-N/(t-1) # number of clusters per step
  X=trtSeq[rep(1:nrow(trtSeq),each=g),]
  
  # constant calculation
  U=sum(X)
  V=sum((rowSums(X))^2)
  W=sum((colSums(X))^2)
  
  
  # ICCs
  rho0.0 <- rho0[1,1]
  rho0.1 <- rho0[2,2]
  rho1.0 <- rho1[1,1] 
  rho1.1 <- rho1[2,2]
  rho0.01 <- rho0[1,2]
  rho1.01 <- rho1[1,2]
  rho2.01 <- rho2[1,2]
  
  # Eigenvalues and taus
  lambda1.0 <- 1-rho0.0
  lambda1.1 <- 1-rho0.1
  lambda2.0 <- 1-N*rho1.0+(N-1)*rho0.0
  lambda2.1 <- 1-N*rho1.1+(N-1)*rho0.1
  lambda3.0 <- 1+(T-1)*N*rho1.0+(N-1)*rho0.0
  lambda3.1 <- 1+(T-1)*N*rho1.1+(N-1)*rho0.1
  tau1 <- rho2.01-rho0.01
  tau2 <- rho2.01-N*rho1.01+(N-1)*rho0.01
  tau3 <- rho2.01+(T-1)*N*rho1.01+(N-1)*rho0.01
  
  # Total variances
  sigma2.y0 <- vars[1]
  sigma2.y1 <- vars[2]
  
  # Variance of treatment effect: vard
    ## Joint model
    vard.JM <- ((N*t/m)*sigma2.y1*(lambda2.0*lambda2.1-tau2^2)*(lambda3.0*lambda3.1-tau3^2))/((N*t*U-t*W+U^2-N*V)*lambda2.0*(lambda3.0*lambda3.1-tau3^2) - (U^2-N*V)*lambda3.0*(lambda2.0*lambda2.1-tau2^2))
  
    ## Change score model
    vard.CS <- ((N*t/m)*(sigma2.y1*lambda3.1 + sigma2.y0*lambda3.0 - 2*sqrt(sigma2.y0)*sqrt(sigma2.y1)*tau3)*(sigma2.y1*lambda2.1 + sigma2.y0*lambda2.0 - 2*sqrt(sigma2.y0)*sqrt(sigma2.y1)*tau2))/
      ((N*t*U-t*W+U^2-N*V)*(sigma2.y1*lambda3.1 + sigma2.y0*lambda3.0 - 2*sqrt(sigma2.y0)*sqrt(sigma2.y1)*tau3) - (U^2-N*V)*(sigma2.y1*lambda2.1 + sigma2.y0*lambda2.0 - 2*sqrt(sigma2.y0)*sqrt(sigma2.y1)*tau2))

    ## Follow-up model: Hooper and Girling
    vard.HG <- ((N*t/m)*sigma2.y1*lambda2.1*lambda3.1)/((N*t*U-t*W+U^2-N*V)*lambda3.1 - (U^2-N*V)*lambda2.1)
  
    ## ANCOVA model: Asymptotic variance
    vard.ANCOVA <- ((N*t/m)*sigma2.y1*(lambda2.1-tau1^2/lambda1.0-2*sigma2.y0*tau1*(tau2-tau1)+sigma2.y0^2*tau1^2*(lambda2.0-lambda1.0))*(lambda3.1-tau1^2/lambda1.0-2*sigma2.y0*tau1*(tau3-tau1)+sigma2.y0^2*tau1^2*(lambda3.0-lambda1.0)))/
      ((U^2+N*t*U-t*W-N*V)*(lambda3.1-2*sigma2.y0*tau1*tau3+sigma2.y0^2*tau1^2*lambda3.0)-(U^2-N*V)*(lambda2.1-2*sigma2.y0*tau1*tau2+sigma2.y0^2*tau1^2*lambda2.0)-(N*U-W)*t*(tau1^2/lambda1.0)*(sigma2.y0*lambda1.0-1)^2)
    
    sigmaks.sq <- c(vard.HG,vard.JM,vard.CS,vard.ANCOVA)
  
  #Power
  criticalValue.t.HG <- qt(p=(1-alpha), df=(N-2))
  criticalValue.t <- qt(p=(1-alpha), df=(N-K-1))
  criticalValue.z <- qnorm(p=(1-alpha))
  pred.power.t.HG <- 1-pt(criticalValue.t, df=(N-2), delta/sqrt(vard.HG))
  pred.power.t.JM <- 1-pt(criticalValue.t, df=(N-K-1), delta/sqrt(vard.JM))
  pred.power.t.CS <- 1-pt(criticalValue.t, df=(N-K-1), delta/sqrt(vard.CS))
  pred.power.t.ANCOVA <- 1-pt(criticalValue.t, df=(N-K-1), delta/sqrt(vard.ANCOVA))
  pred.power.z.HG <- 1-pnorm(criticalValue.z, mean=(delta/sqrt(vard.HG)))
  pred.power.z.JM <- 1-pnorm(criticalValue.z, mean=(delta/sqrt(vard.JM)))
  pred.power.z.CS <- 1-pnorm(criticalValue.z, mean=(delta/sqrt(vard.CS)))
  pred.power.z.ANCOVA <- 1-pnorm(criticalValue.z, mean=(delta/sqrt(vard.ANCOVA)))
  
  pred.power.t <- c(pred.power.t.HG,pred.power.t.JM,pred.power.t.CS,pred.power.t.ANCOVA)
  pred.power.z <- c(pred.power.z.HG,pred.power.z.JM,pred.power.z.CS,pred.power.z.ANCOVA)
  
  param <- list(vard=c(sigmaks.sq),pred.power.t=c(pred.power.t),pred.power.z=c(pred.power.z))
  return(param)
}
