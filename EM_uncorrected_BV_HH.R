library(mvtnorm)
library(numDeriv)

# function to perform EM estimation with K=2 outcomes
EM.estim <- function(data, fm1,fm2, maxiter=500,epsilon=1e-4
                     , verbose=FALSE){
  K <- 2
  zeta <- as.numeric(c(fm1$coefficients$fixed, fm2$coefficients$fixed))
  beta1 = as.numeric(fm1$coefficients$fixed)
  beta2 = as.numeric(fm2$coefficients$fixed)
  nvar<-length(as.numeric(fm1$coefficients$fixed))
  TermsX1 <- fm1$terms
  TermsX2 <- fm2$terms
  mfX1 <- model.frame(TermsX1, data = data)[,-1]
  mfX2 <- model.frame(TermsX2, data = data)[,-1]
  
  # vector of cluster sizes
  m <- as.numeric(table(fm1$groups[[1]]))
  
  s2phi1 <- VarCorr(fm1)[1,1]
  s2phi2 <- VarCorr(fm2)[1,1]
  SigmaPhi <- diag(c(s2phi1, s2phi2)) # Initialized assuming independence
  InvS2Phi <- solve(SigmaPhi)
  
  s2e1 <- VarCorr(fm1)[2,1]
  s2e2 <- VarCorr(fm2)[2,1]
  SigmaE <- diag(c(s2e1, s2e2))
  InvS2E <- solve(SigmaE)
  
  Y <- as.matrix(cbind(model.frame(TermsX1, data = data)[,1],model.frame(TermsX2, data = data)[,1]))
  ID <- fm1$groups[[1]]
  n <- length(unique(ID))
  X1 <- as.matrix(cbind(1, mfX1)) # design matrix for baseline
  X2 <- as.matrix(cbind(1, mfX2)) # design matrix for endline
  
  ESSphi1 <- matrix(0,n,K)
  ESSphi2 <- array(0,c(K,K,n))
  
  
  #maxiter=500
  #epsilon=1e-4
  delta = 2*epsilon
  max_modi = 20
  
  converge = 0
  
  # log likelihood
  loglik = function(theta){
    beta1 = theta[1:nvar]
    beta2 = theta[(nvar+1):(2*nvar+1)]
    sphi11 = theta[(2*nvar+2)]
    sphi12 = theta[(2*nvar+3)]
    sphi22 = theta[(2*nvar+4)]
    se11 = theta[(2*nvar+5)]
    se12 = theta[(2*nvar+6)]
    se22 = theta[(2*nvar+7)]
    SigmaPhi = matrix(c(sphi11,sphi12,sphi12,sphi22),2,2)
    SigmaE = matrix(c(se11,se12,se12,se22),2,2)
    InvS2Phi <- solve(SigmaPhi)
    InvS2E <- solve(SigmaE)
    
    temp <- 0
    for(j in 1:n){
      Yj <- Y[ID == j,,drop=FALSE]
      X1j <- X1[ID == j,,drop=FALSE]
      X2j <- X2[ID == j,,drop=FALSE]
      residj <- Yj - cbind(X1j%*%beta1, X2j%*%beta2)
      obs = c(t(residj))
      tm1 <- (m[j]-1)*log(det(SigmaE))+log(det(SigmaE+m[j]*SigmaPhi))
      InvSS2 <- solve(SigmaE+m[j]*SigmaPhi)-InvS2E
      Invj <- kronecker(diag(nrow=m[j]),InvS2E) + 
        kronecker(matrix(1,m[j],m[j]),InvSS2)/m[j]
      tm2 <- c(t(obs) %*% Invj %*% obs)
      temp <- temp-(tm1+tm2)/2
    }
    temp
  }
  thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaE[!lower.tri(SigmaE)]))
  LLold <- loglik(thetah)
  
  
  niter=1
  while((niter <= maxiter) & (abs(delta) > epsilon)){
    
    # Expectation step
    for(j in 1:n){
      Yj <- Y[ID == j,,drop=FALSE]
      X1j <- X1[ID == j,,drop=FALSE]
      X2j <- X2[ID == j,,drop=FALSE]
      residj <- Yj - cbind(X1j%*%beta1, X2j%*%beta2)
      Vj <- solve(InvS2Phi + m[j]*InvS2E)
      Muj <- as.numeric(Vj %*% InvS2E %*% colSums(residj))
      Nujj <- Vj + tcrossprod(Muj)
      ESSphi1[j,] <- Muj
      ESSphi2[,,j] <- Nujj
    }
    
    # Maximization step - phi
    SigmaPhi <- apply(ESSphi2,1:2, sum)/n
    InvS2Phi <- solve(SigmaPhi)
    
    # Maximization step - zeta
    # Let Beta = (baseline betas, endline betas)^T
    DSigmaED <- 0
    DSigmaEResid <- 0
    rzeta1 <- Y[,1]-ESSphi1[ID,1]
    rzeta2 <- Y[,2]-ESSphi1[ID,2]
    for(j in 1:nrow(data)){
      Xj1 <- X1[j,,drop=FALSE]
      Xj2 <- X2[j,,drop=FALSE]
      X<-as.matrix(rbind(c(Xj1,rep(0,length(Xj2))),c(rep(0,length(Xj1)),Xj2)))
      DSigmaED <- t(X)%*%InvS2E%*%X + DSigmaED
      DSigmaEResid <- t(X)%*%InvS2E%*%as.matrix(rbind(rzeta1[j,drop=FALSE],rzeta2[j,drop=FALSE])) + DSigmaEResid
    }
    
    zeta <- solve(DSigmaED)%*%DSigmaEResid
    zeta <- c(zeta)
    beta1 = zeta[1:nvar]
    beta2 = zeta[(nvar+1):(2*nvar+1)]
    
    # Maximization step - epsilon
    re <- Y - cbind(X1%*%beta1, X2%*%beta2)
    rss <- crossprod(re) + rowSums(sweep(ESSphi2,3,m,FUN="*"),dims=2) -
      crossprod(ESSphi1,rowsum(re,ID)) - crossprod(rowsum(re,ID),ESSphi1)
    SigmaE <- rss/sum(m)
    InvS2E <- solve(SigmaE)
    
    # whether the algorithm converges
    thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaE[!lower.tri(SigmaE)]))
    LLnew <- loglik(thetah)
    delta <- abs(LLnew - LLold)
    LLold <- LLnew
    converge = (abs(delta)<=epsilon)
    niter <- niter + 1
    if(verbose) cat(paste('iter=',niter),'\t',paste('param.error=',epsilon),'\t',paste('loglik=',LLnew),'\n');  
    
    #print(niter)
    #print(zeta)
    #print(SigmaPhi)
    #print(SigmaE)
    #print(LLnew)
  }
  
  Vtheta <- matrix(NA,length(thetah),length(thetah))
  try(Vtheta <- solve(-hessian(loglik,thetah))) #, silent = TRUE
  SEcheck <- "GOOD"
  if(anyNA(Vtheta)==TRUE|min(diag(Vtheta))<0) SEcheck <- "ERROR"
  #Vtheta = try(solve(-hessian(loglik,thetah)))
  #SEtheta = sqrt(diag(Vtheta))
  
  param <- list(theta=list(zeta=zeta,SigmaPhi=SigmaPhi,SigmaE=SigmaE),loglik=LLnew,eps=epsilon,iter=niter,
                Vtheta=Vtheta,SEcheck=SEcheck)#SEtheta=SEtheta)
  return(param) 
}