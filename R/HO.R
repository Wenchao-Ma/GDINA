

HO.est <- function(lambda, AlphaPattern, HOgr, Rl, higher.order)
{

  Aq <- higher.order$QuadNodes
  WAq <- higher.order$QuadWghts
  if(is.null(lambda)){
    lambda = higher.order$lambda
  }
  K <- log2(nrow(AlphaPattern))

   NR <- list()
   for(g in HOgr){
    NR[[g]] <- expectedNR(AlphaPattern, Rl[,g], Aq[,g], WAq[,g], lambda[[g]][,1], lambda[[g]][,2])
   }


  if(higher.order$Type=="SameTheta"){
  for(g in HOgr){
    for(k in seq_len(K)){
      lambda[[g]][k,2] <- rootFinder(f = Lfj_intercept, interval = higher.order$InterceptRange,
                           aj = lambda[[g]][k,1],theta = Aq[,g],r = NR[[g]]$r[,k], n = NR[[g]]$n,
                           prior = higher.order$Prior, mu = higher.order$InterceptPrior[1], sigma = higher.order$InterceptPrior[2])
    }
    if(higher.order$model=="1PL"){
      lambda[[g]][,1] <- rootFinder(f = Lfj_commonslope, interval = higher.order$SlopeRange,
                          d = lambda[[g]][,2],theta = Aq[,g],r = NR[[g]]$r, n = NR[[g]]$n,
                          prior = higher.order$Prior, mu = higher.order$SlopePrior[1], sigma = higher.order$SlopePrior[2])
    }else if(higher.order$model=="2PL"){

      for(k in seq_len(K)) lambda[[g]][k,1] <- rootFinder(f = Lfj_slope, interval = higher.order$SlopeRange,
                                                d = lambda[[g]][k,2],theta = Aq[,g],r = NR[[g]]$r[,k], n = NR[[g]]$n,
                                                prior = higher.order$Prior, mu = higher.order$SlopePrior[1], sigma = higher.order$SlopePrior[2])

    }else if(higher.order$model=="Rasch"){
      lambda[[g]][,1] <- 1
    }else{
      stop("Higher-order model is not correctly specified.",call. = FALSE)
    }
  }

  }else if(higher.order$Type=="SameLambda"){
    n <- r <- 0
    for(g in HOgr){
      n <- n + NR[[g]]$n
      r <- r + NR[[g]]$r
    }
    d <- a <- vector("numeric",K)
    for(k in seq_len(K)){
      d[k] <- rootFinder(f = Lfj_intercept, interval = higher.order$InterceptRange,
                          aj = lambda[[HOgr[1]]][k,1],theta = Aq[,1],r = r[,k], n = n,
                          prior = higher.order$Prior, mu = higher.order$InterceptPrior[1],
                          sigma = higher.order$InterceptPrior[2])
    }
    if(higher.order$model=="1PL"){
      a <- rootFinder(f = Lfj_commonslope, interval = higher.order$SlopeRange,
                             d = lambda[[HOgr[1]]][,2],theta = Aq[,1],r = r, n = n,
                             prior = higher.order$Prior, mu = higher.order$SlopePrior[1],
                             sigma = higher.order$SlopePrior[2])
    }else if(higher.order$model=="2PL"){
      for(k in seq_len(K)) a[k,] <- rootFinder(f = Lfj_slope, interval = higher.order$SlopeRange,
                                               d = lambda[[HOgr[1]]][k,2],theta = Aq[,1],r = r[,k], n = n,
                                               prior = higher.order$Prior, mu = higher.order$SlopePrior[1],
                                               sigma = higher.order$SlopePrior[2])
    }else if(higher.order$model=="Rasch"){
      a <- 1
    }else{
      stop("Higher-order model is not correctly specified.",call. = FALSE)
    }
    for(g in HOgr){
      lambda[[g]][,1] <- a
      lambda[[g]][,2] <- d
    }

    Pg_alpha_c <- ColNormalize(Rl)
    for(gg in 2:max(HOgr)){
      WAq[,gg] <-
        higher.order$QuadWghts[,gg]  <-
        colSums(PostTheta(AlphaPattern, Aq[,gg], WAq[,gg], lambda[[gg]][,1],lambda[[gg]][,2])*Pg_alpha_c[,gg])
    }

  }
logprior <- matrix(0,nrow(AlphaPattern),length(HOgr))
for(g in HOgr){
  logprior[,g] <- logP_AlphaPattern(AlphaPattern, Aq[,g], WAq[,g], lambda[[g]][,1], lambda[[g]][,2])
}
  return(list(logprior=logprior,lambda=lambda,higher.order=higher.order))
}


rootFinder <- function(f,interval,...){
  f1 <- f(interval[1],...)
  f2 <- f(interval[2],...)
  if(f1*f2>0 | f1==f2){
    ret <- ifelse(abs(f1)>abs(f2),interval[2],interval[1])
  }else{
    ret <- uniroot(f = f, interval = interval,...)$root
  }
  ret
}

Lfj_intercept <- function(dj,aj,theta,r,n, prior = FALSE, mu, sigma){
  P <- Pr_2PL_vec(theta = theta, a = aj, b = dj) #nnodes x 1
  if(prior){
    ret <- sum(c(r)-c(n)*c(P)) + (dj - mu) / (sigma^2)
  }else{
    ret <- sum(c(r)-c(n)*c(P))
  }
  ret
}
Lfj_slope <- function(aj,dj,theta,r,n, prior = FALSE, mu, sigma){
  P <- Pr_2PL_vec(theta = theta, a = aj, b = dj) #nnodes x 1
  if(prior){
    ret <- sum(theta*(c(r)-c(n)*c(P))) - (log(aj) - mu + sigma^2) /(aj * sigma^2)
  }else{
    ret <- sum(theta*(c(r)-c(n)*c(P)))
  }
  ret

}
Lfj_commonslope <- function(a,d,theta,r,n, prior = FALSE, mu, sigma){
  avec <- rep(a,length(d))
  P <- Pr_2PL_vec(theta = theta, a = avec, b = d) #nnodes x K

  if(prior){
    ret <- sum(theta*(r-c(n)*P)) - (log(a) - mu + sigma^2) /(a * sigma^2)
  }else{
    ret <- sum(theta*(r-c(n)*P))
  }
  ret
}

# HO.SE <- function(a,d,Xloglik,alphapattern,theta,weight,K,model="1PL",h = 1e-7){
#   # This will overestimate standard errors
#   inp <- incomplogL(a=a, b=d,logL = Xloglik, AlphaPattern = alphapattern, theta = theta, f_theta = weight)
#
#   hess <- matrix(NA,K,2)
#   if(model=="2PL"){
#
#       for(rr in 1:K){
#         aa <- a
#         dd <- d
#         x <- rep(0,K)
#         x[rr] <- h
#         hess[rr,1] <- (incomplogL(a=aa+x, b=d,logL = Xloglik, AlphaPattern = alphapattern, theta = theta, f_theta = weight)- 2*inp +
#                          incomplogL(a=aa-x, b=d,logL = Xloglik, AlphaPattern = alphapattern, theta = theta, f_theta = weight))/(h^2)
#         hess[rr,2] <- (incomplogL(a=a, b=dd+x,logL = Xloglik, AlphaPattern = alphapattern, theta = theta, f_theta = weight)- 2*inp +
#                          incomplogL(a=a, b=dd-x,logL = Xloglik, AlphaPattern = alphapattern, theta = theta, f_theta = weight))/(h^2)
#       }
#
#   }else if(model=="1PL"){
#     hess[,1] <- (incomplogL(a=aa+h, b=d,logL = Xloglik, AlphaPattern = alphapattern, theta = theta, f_theta = weight)- 2*inp +
#                      incomplogL(a=aa-h, b=d,logL = Xloglik, AlphaPattern = alphapattern, theta = theta, f_theta = weight))/(h^2)
#     for(rr in 1:K){
#       dd <- d
#       x <- rep(0,K)
#       x[rr] <- h
#       hess[rr,2] <- (incomplogL(a=a, b=dd+x,logL = Xloglik, AlphaPattern = alphapattern, theta = theta, f_theta = weight)- 2*inp +
#                        incomplogL(a=a, b=dd-x,logL = Xloglik, AlphaPattern = alphapattern, theta = theta, f_theta = weight))/(h^2)
#     }
#   }else if(model=="Rasch"){
#     for(rr in 1:K){
#       dd <- d
#       x <- rep(0,K)
#       x[rr] <- h
#       hess[rr,2] <- (incomplogL(a=a, b=dd+x,logL = Xloglik, AlphaPattern = alphapattern, theta = theta, f_theta = weight)- 2*inp +
#                        incomplogL(a=a, b=dd-x,logL = Xloglik, AlphaPattern = alphapattern, theta = theta, f_theta = weight))/(h^2)
#     }
#   }
#
#   return(hess)
# }

# HO.SE <- function(v,Xloglik,alphapattern,theta,weight,K,model="1PL"){
#
#     if(model=="2PL"){
#     a <- v[1:K]
#     d <- v[(K+1):(2*K)]
#   }else if(model=="1PL"){
#     a <- rep(v[1],K)
#     d <- v[-1]
#   }else if(model=="Rasch"){
#     a <- rep(1,K)
#     d <- v
#   }
#
#   incomplogL(a=a, b=d,logL = Xloglik, AlphaPattern = alphapattern, theta = theta, f_theta = weight)
# }
