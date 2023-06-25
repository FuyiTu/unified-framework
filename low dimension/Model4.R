Model4_SRS<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=2 || length(table(S1))!=2 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    S = sample(c(1,-1),n,replace = TRUE)
    A = SRS(n,pi)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2)
  Y0 = alpha0 + Xbeta%*%betavec0[c(1,2)]*S + betavec0[3]*log(X1+1)*(S==1) + rnorm(n,sd = sigma0)
  Y1 = alpha1 + Xbeta%*%betavec1[c(1,2)]*S + betavec1[3]*exp(X2)*(S==-1) + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model4_WEI<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=2 || length(table(S1))!=2 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    S = sample(c(1,-1),n,replace = TRUE)
    A = WEI(t(S),pi)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2)
  Y0 = alpha0 + Xbeta%*%betavec0[c(1,2)]*S + betavec0[3]*log(X1+1)*(S==1) + rnorm(n,sd = sigma0)
  Y1 = alpha1 + Xbeta%*%betavec1[c(1,2)]*S + betavec1[3]*exp(X2)*(S==-1) + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model4_SBR<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=2 || length(table(S1))!=2 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    S = sample(c(1,-1),n,replace = TRUE)
    A = SBR(t(S),pi)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2)
  Y0 = alpha0 + Xbeta%*%betavec0[c(1,2)]*S + betavec0[3]*log(X1+1)*(S==1) + rnorm(n,sd = sigma0)
  Y1 = alpha1 + Xbeta%*%betavec1[c(1,2)]*S + betavec1[3]*exp(X2)*(S==-1) + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model4_BCD<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,lambda){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=2 || length(table(S1))!=2 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    S = sample(c(1,-1),n,replace = TRUE)
    A = BCD(t(S),pi,lambda)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2)
  Y0 = alpha0 + Xbeta%*%betavec0[c(1,2)]*S + betavec0[3]*log(X1+1)*(S==1) + rnorm(n,sd = sigma0)
  Y1 = alpha1 + Xbeta%*%betavec1[c(1,2)]*S + betavec1[3]*exp(X2)*(S==-1) + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model4_PS<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,weight,lambda){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=2 || length(table(S1))!=2 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    S = sample(c(1,-1),n,replace = TRUE)
    A = PocSimue(t(S),pi,weight,lambda)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2)
  Y0 = alpha0 + Xbeta%*%betavec0[c(1,2)]*S + betavec0[3]*log(X1+1)*(S==1) + rnorm(n,sd = sigma0)
  Y1 = alpha1 + Xbeta%*%betavec1[c(1,2)]*S + betavec1[3]*exp(X2)*(S==-1) + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}

Model4_HH<-function(n,alpha0,alpha1,betavec0,betavec1,pi,sigma0,sigma1,omega,lambda){
  S0 = rep(0,n)
  S1 = rep(0,n)
  while(length(table(S0))!=2 || length(table(S1))!=2 || min(table(S0)) < 3 || min(table(S1)) < 3){
    X1 = rbeta(n,3,4)
    X2 = runif(n,-2,2)
    S = sample(c(1,-1),n,replace = TRUE)
    A = HHue(t(S),pi,omega,lambda)
    S0 = S[which(A == 0)]
    S1 = S[which(A == 1)]
  }
  Xbeta = cbind(X1,X2)
  Y0 = alpha0 + Xbeta%*%betavec0[c(1,2)]*S + betavec0[3]*log(X1+1)*(S==1) + rnorm(n,sd = sigma0)
  Y1 = alpha1 + Xbeta%*%betavec1[c(1,2)]*S + betavec1[3]*exp(X2)*(S==-1) + rnorm(n,sd = sigma1)
  return(list(A=A,S=S,Xbeta=data.frame(Xbeta),Y=Y0*(1-A)+Y1*A))
}