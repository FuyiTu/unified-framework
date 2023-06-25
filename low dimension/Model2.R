Model2_SRS<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=4 || length(table(S1))!=4 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    S = sample(c(1,2,3,4),n,replace = TRUE, prob = c(0.2,0.3,0.3,0.2))
    A = SRS(n,pi)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2)
  Y0 = alpha0 + betavec0[1]*log(X1+1) + betavec0[2]*X1^2 + betavec0[3]*exp(X2) + betavec0[4]/(X2+3) + rnorm(n,sd = sigma0)
  Y1 = alpha1 + betavec1[1]*exp(X1+2) + betavec1[2]/(X1+1) + betavec1[3]*X2^2 + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model2_WEI<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=4 || length(table(S1))!=4 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    S = sample(c(1,2,3,4),n,replace = TRUE, prob = c(0.2,0.3,0.3,0.2))
    A = WEI(t(S),pi)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2)
  Y0 = alpha0 + betavec0[1]*log(X1+1) + betavec0[2]*X1^2 + betavec0[3]*exp(X2) + betavec0[4]/(X2+3) + rnorm(n,sd = sigma0)
  Y1 = alpha1 + betavec1[1]*exp(X1+2) + betavec1[2]/(X1+1) + betavec1[3]*X2^2 + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model2_SBR<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=4 || length(table(S1))!=4 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    S = sample(c(1,2,3,4),n,replace = TRUE, prob = c(0.2,0.3,0.3,0.2))
    A = SBR(t(S),pi)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2)
  Y0 = alpha0 + betavec0[1]*log(X1+1) + betavec0[2]*X1^2 + betavec0[3]*exp(X2) + betavec0[4]/(X2+3) + rnorm(n,sd = sigma0)
  Y1 = alpha1 + betavec1[1]*exp(X1+2) + betavec1[2]/(X1+1) + betavec1[3]*X2^2 + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model2_BCD<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,lambda){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=4 || length(table(S1))!=4 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    S = sample(c(1,2,3,4),n,replace = TRUE, prob = c(0.2,0.3,0.3,0.2))
    A = BCD(t(S),pi,lambda)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2)
  Y0 = alpha0 + betavec0[1]*log(X1+1) + betavec0[2]*X1^2 + betavec0[3]*exp(X2) + betavec0[4]/(X2+3) + rnorm(n,sd = sigma0)
  Y1 = alpha1 + betavec1[1]*exp(X1+2) + betavec1[2]/(X1+1) + betavec1[3]*X2^2 + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model2_PS<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,weight,lambda){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=4 || length(table(S1))!=4 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    S = sample(c(1,2,3,4),n,replace = TRUE, prob = c(0.2,0.3,0.3,0.2))
    A = PocSimue(t(S),pi,weight,lambda)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2)
  Y0 = alpha0 + betavec0[1]*log(X1+1) + betavec0[2]*X1^2 + betavec0[3]*exp(X2) + betavec0[4]/(X2+3) + rnorm(n,sd = sigma0)
  Y1 = alpha1 + betavec1[1]*exp(X1+2) + betavec1[2]/(X1+1) + betavec1[3]*X2^2 + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model2_HH<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,omega,lambda){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=4 || length(table(S1))!=4 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    S = sample(c(1,2,3,4),n,replace = TRUE, prob = c(0.2,0.3,0.3,0.2))
    A = HHue(t(S),pi,omega,lambda)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2)
  Y0 = alpha0 + betavec0[1]*log(X1+1) + betavec0[2]*X1^2 + betavec0[3]*exp(X2) + betavec0[4]/(X2+3) + rnorm(n,sd = sigma0)
  Y1 = alpha1 + betavec1[1]*exp(X1+2) + betavec1[2]/(X1+1) + betavec1[3]*X2^2 + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}