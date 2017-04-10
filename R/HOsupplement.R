############################################
#
#  HO 1PL model
#
############################################
HO.1PL <- function(vP, Nq, Rq,Aq,mu.loga,mu.b,sigma2.loga,sigma2.b,prior=FALSE){ # vP: c(log-slope,all bs)
  # vP here is a vector of length K
  # Rq is a Q x K matrix
  # Nq is a vector of length Q
  a <- exp(vP[1])
  b <- vP[-1]
  P <- 1/(1 + exp(-(outer(a*Aq,b,"+")))) # Q x K
  if(!prior) {
    sum(( (Rq) * log(P) ) + ( (Nq - Rq) * log(1-P) ) )
  }else{
    sum(( (Rq) * log(P) ) + ( (Nq - Rq) * log(1-P) ) ) +
      sum(dnorm(vP[-1],mean = mu.b, sd=sqrt(sigma2.b),log = TRUE)) +
    dnorm(vP[1], mean = mu.loga, sd = sqrt(sigma2.loga),log = TRUE)
  }
}



############################################
#
#  HO Rasch model
#
############################################
HO.Rasch <- function(vP, Nq, Rq,Aq,mu,sigma2,prior=FALSE){ # vP: all bs
   # vP here is a vector of length K
  # Rq is a Q x K matrix
  # Nq is a vector of length Q
  P <- 1/(1 + exp(-(outer(Aq,vP,"+")))) # Q x K
  # P <- 1/(1 + exp(-vP-Aq))
  if(!prior) {
    sum(( (Rq) * log(P) ) + ( (Nq - Rq) * log(1-P) ) )
  }else{
    sum(( (Rq) * log(P) ) + ( (Nq - Rq) * log(1-P) ) ) +
      sum(dnorm(vP,mean = mu, sd=sqrt(sigma2),log = TRUE))
  }
}

HO.Rasch.j <- function(vP, Nq, Rqk,Aq,mu,sigma2){ # vP: c(intercept)
  P <- 1/(1 + exp(-vP-Aq))
  if(is.null(mu)||is.null(sigma2)) {
    sum(( (Rqk) * log(P) ) + ( (Nq - Rqk) * log(1-P) ) )
  }else{
    sum(( (Rqk) * log(P) ) + ( (Nq - Rqk) * log(1-P) ) ) + dnorm(vP,mean = mu, sd=sqrt(sigma2),log = TRUE)
  }
}

gr.Rasch.j <- function(vP,Nq,Rqk,Aq,mu,sigma2,a,prior){ # vP: c(intercept)
  P <- 1/(1 + exp(-(a*Aq+vP)))
  if(!prior) {
    (sum(Rqk-P*Nq))^2
  }else{
    (sum(Rqk-P*Nq)-(vP-mu)/sigma2)^2
  }
}

obj_fn_b <- function(b,Nq,Rqk,Aq,mu,sigma2,a=1,prior){ # b: c(intercept)
  P <- 1/(1 + exp(-(Aq+b)))
  if(!prior) {
    (sum(Rqk-P*Nq))^2
  }else{
    (sum(Rqk-P*Nq)-(b-mu)/sigma2)^2
  }
}

############################################
#
#  HO 2PL model
#
############################################
HO.2PL <- function(vP, Nq, Rq,Aq,mu.loga,mu.b,sigma2.loga,sigma2.b,prior=FALSE){ # vP: c(log[slope]1,log[slope]2,...intercept1,...)

  a <- exp(vP[1:(length(vP)/2)])
  b <- vP[(1 + (length(vP)/2)):length(vP)]
  P <- outer(a,Aq) + b # K x nnodes
  P <- 1/(1 + exp(-P))  # K x nnodes
if(!prior){
  sum(( t(Rq) * log(P) ) + ( t(Nq - Rq) * log(1-P) ) )  # sum across Q
}else{
  sum(( t(Rq) * log(P) ) + ( t(Nq - Rq) * log(1-P) ) )  + # sum across Q
  sum(dnorm(b,mean = mu.b, sd= sqrt(sigma2.b),log = TRUE)) +
  sum(dnorm(vP[1:(length(vP)/2)], mean = mu.loga, sd = sqrt(sigma2.loga),log = TRUE))
}


}

HO.2PL.j <- function(vP, Nq, Rqk,Aq,mu.loga,mu.b,sigma2.loga,sigma2.b,prior=FALSE){ # vP: c(log[slope],intercept)

  P <- 1/(1 + exp(-exp(vP[1]) * Aq - vP[2])) # length Q
if(!prior){
  sum(( (Rqk) * log(P) ) + ( (Nq - Rqk) * log(1-P) ) )  # sum across Q
}else{
  sum(( (Rqk) * log(P) ) + ( (Nq - Rqk) * log(1-P) ) )  +
  dnorm(vP[2],mean = mu.b, sd= sqrt(sigma2.b),log = TRUE) +
  dnorm(vP[1], mean = mu.loga, sd = sqrt(sigma2.loga),log = TRUE)
}


}

gr.2PL.j <- function(vP, Nq, Rqk,Aq,mu,sigma){ # vP: c(log[slope],intercept)

  P <- 1/(1 + exp(-exp(vP[1]) * Aq - vP[2])) # length Q

  c(sum(exp(vP[1])*Aq*(Rqk-P*Nq))-(vP[1]-mu[1])/(sigma[1]^2),
    sum(Rqk-P*Nq)-(vP[2]-mu[2])/sigma[2]^2)

}


obj_fn_1PLa<-function(loga, Nq, Rq,Aq,mu=0,sigma2=0.25,b,prior=FALSE){
  # we estimate alpha = log(a), which is assumed normally distributed N(mu,sigma2)
  # b here is a vector of length K
  # Rq is a Q x K matrix
  # Nq is a vector of length Q
  P <- 1/(1 + exp(-(outer(exp(loga) * Aq,b,"+")))) # Q x K
  if(!prior) {
    sum(exp(loga)*Aq*(Rq-P*Nq))^2
  }else{
    (sum(exp(loga)*Aq*(Rq-P*Nq))-(loga-mu)/sigma2)^2
  }

}
obj_fn_a<-function(loga, Nq, Rqk,Aq,mu=0,sigma2=0.25,b,prior=FALSE){
 #we estimate alpha = log(a), which is assumed normally distributed N(mu,sigma2)
  P <- 1/(1 + exp(-(exp(loga) * Aq + b)))
  if(!prior) {
    (sum(exp(loga)*Aq*(Rqk-P*Nq)))^2
  }else{
    (sum(exp(loga)*Aq*(Rqk-P*Nq))-(loga-mu)/sigma2)^2
  }
  # sum(exp(loga)*Aq*(Rqk-P*Nq)-(loga-mu)/sigma2)^2
  # sum(( (Rqk) * log(P) ) + ( (Nq - Rqk) * log(1-P) ) )

}

obj_fn_b <- function(b,Nq,Rqk,Aq,mu,sigma2,a=1,prior=FALSE){ # b: c(intercept)
  P <- 1/(1 + exp(-(Aq+b)))
  if(!prior) {
    (sum(Rqk-P*Nq))^2
  }else{
    (sum(Rqk-P*Nq)-(b-mu)/sigma2)^2
  }
}

# HO.loglik <- function(a, b, theta, nnodes, X,K)
# {
#   # nnodes x K for z
#   z <- outer(a,theta) + b
#   P <- 1/(1 + exp(-t(z)))  # nnodes x K
#   P[P>0.999] <- 0.999
#   P[P<0.001] <- 0.001
#   return(log(P)%*%t(X)+log(1-P)%*%(1-t(X)))  # log { \prod_k P(alpha_k|Aq)^alpha_k[(1-P(alpha_k|Aq))^(1-alpha_k)] } nnodes x 2^K
# }

HO.loglik <- function(a, b, theta, nnodes, X, K) # log P(alpha_c|theta,a,b) - 2^K x nnodes
{
  # X - all alpha pattern 2^K x K
  P <- outer(a,theta) + b # K x nnodes
  P <- 1/(1 + exp(-P))  # K x nnodes
  P[P>0.999] <- 0.999
  P[P<0.001] <- 0.001
  return(X%*%log(P)+(1-X)%*%log(1-P))
  # return(log(P)%*%t(X)+log(1-P)%*%(1-t(X)))  # log { \prod_k P(alpha_k|Aq)^alpha_k[(1-P(alpha_k|Aq))^(1-alpha_k)] } nnodes x 2^K
}

HO.SE.2 <- function(par,Xloglik,K,nnodes=19,N){
  a <- par[1:K];b <- par[(K+1):(2*K)]
  Aq <- seq(-3, 3, length.out = nnodes)
  X <- alpha(K)
  Wght <- dnorm(Aq)
  HOlike <- exp(HO.loglik(a,b,Aq,nnodes,X,K))

  #Compute Nq - NxL %*% L x nnodes * N x nnodes -> N x nnodes
  sum(log(rowSums(exp(Xloglik) %*% HOlike * matrix(Wght, N, nnodes, byrow = T))))
}
HO.SE.1 <- function(par,Xloglik,K,nnodes=19,N){
  a <- rep(par[1],K);b <- par[2:(1+K)]
  Aq <- seq(-3, 3, length.out = nnodes)
  X <- alpha(K)
  Wght <- dnorm(Aq)
  HOlike <- exp(HO.loglik(a,b,Aq,nnodes,X,K))

  #Compute Nq - NxL %*% L x nnodes * N x nnodes -> N x nnodes
  sum(log(rowSums(exp(Xloglik) %*% HOlike * matrix(Wght, N, nnodes, byrow = T))))
}
HO.SE.R <- function(par,Xloglik,K,nnodes=19,N){
  a <- rep(1,K);b <- par
  Aq <- seq(-3, 3, length.out = nnodes)
  X <- alpha(K)
  Wght <- dnorm(Aq)
  HOlike <- exp(HO.loglik(a,b,Aq,nnodes,X,K))

  #Compute Nq - NxL %*% L x nnodes * N x nnodes -> N x nnodes
  sum(log(rowSums(exp(Xloglik) %*% HOlike * matrix(Wght, N, nnodes, byrow = T))))
}
