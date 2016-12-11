HO.2PL<-function(vP, Nq, Rqk,Aq){

  P <- (1/(1 + exp(-(vP[1] * Aq + vP[2]))))  # nnodes x K

  sum(( (Rqk) * log(P) ) + ( (Nq - Rqk) * log(1-P) ) )#nnodes x K

}
HO.Rasch<-function(vP, Nq, Rqk,Aq){
  P <- (1/(1 + exp(-vP-Aq)))  # nnodes x K
  sum(( (Rqk) * log(P) ) + ( (Nq - Rqk) * log(1-P) ) )#nnodes x K

}
HO.loglik <- function(a, b, theta, nnodes, X,K)
{
  # nnodes x K for z
  z <- outer(a,theta) + b
  P <- 1/(1 + exp(-t(z)))  # nnodes x K
  P[P>0.999] <- 0.999
  P[P<0.001] <- 0.001
  return(log(P)%*%t(X)+log(1-P)%*%(1-t(X)))  # log { \prod_k P(alpha_k|Aq)^alpha_k[(1-P(alpha_k|Aq))^(1-alpha_k)] } nnodes x 2^K
}

HO.SE.2 <- function(par,Xloglik,K,nnodes=19,N){
  a <- par[1:K];b <- par[(K+1):(2*K)]
  Aq <- seq(-3, 3, length.out = nnodes)
  X <- alpha(K)
  Wght <- dnorm(Aq)
  HOlike <- exp(t(HO.loglik(a,b,Aq,nnodes,X,K)))

  #Compute Nq - NxL %*% L x nnodes * N x nnodes -> N x nnodes
  sum(log(rowSums(exp(Xloglik) %*% HOlike * matrix(Wght, N, nnodes, byrow = T))))
}
HO.SE.1 <- function(par,Xloglik,K,nnodes=19,N){
  a <- rep(par[1],K);b <- par[2:(1+K)]
  Aq <- seq(-3, 3, length.out = nnodes)
  X <- alpha(K)
  Wght <- dnorm(Aq)
  HOlike <- exp(t(HO.loglik(a,b,Aq,nnodes,X,K)))

  #Compute Nq - NxL %*% L x nnodes * N x nnodes -> N x nnodes
  sum(log(rowSums(exp(Xloglik) %*% HOlike * matrix(Wght, N, nnodes, byrow = T))))
}
HO.SE.R <- function(par,Xloglik,K,nnodes=19,N){
  a <- rep(1,K);b <- par
  Aq <- seq(-3, 3, length.out = nnodes)
  X <- alpha(K)
  Wght <- dnorm(Aq)
  HOlike <- exp(t(HO.loglik(a,b,Aq,nnodes,X,K)))

  #Compute Nq - NxL %*% L x nnodes * N x nnodes -> N x nnodes
  sum(log(rowSums(exp(Xloglik) %*% HOlike * matrix(Wght, N, nnodes, byrow = T))))
}
