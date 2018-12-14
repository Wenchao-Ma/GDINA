structural.parm.mg <- function(AlphaPattern,no.mg,logprior,att.dist,att.str,saturated,initial.logprior,
                               lower.prior,loglinear=2,N,Ng,K,lambda,higher.order){
  prior <- exp(logprior)
  parm <- NULL
  ho <- NULL
  npar <- 0
  if (all(att.dist=="higher.order")){
    if(!is.null(att.str))
      stop("Only att.dist = saturated or fixed is available when attributes have structures.",call. = FALSE)
    HO.out <- HO.est(lambda=lambda,AlphaPattern = AlphaPattern, HOgr = seq_len(no.mg), Rl = rowProd(prior,Ng),
                     higher.order = higher.order)
    logprior <- HO.out$logprior
    lambda <- HO.out$lambda
    ho <- HO.out$higher.order
    npar <- HO.out$npar
  }else{
    for(g in seq_len(no.mg)){
      if (att.dist[g]=="saturated"){
        prior[,g] <- (prior[,g]*Ng[g]+saturated$prior[,g])/sum(prior[,g]*Ng[g]+saturated$prior[,g])
        prior[which(prior[,g]<lower.prior[g]),g] <- lower.prior[g]
        prior[,g] <- prior[,g]/sum(prior[,g])
        logprior[,g] <- log(prior[,g])
        lambda[[g]] <- prior[,g]
        npar <- npar + length(lambda[[g]]) - 1
      }else if (att.dist[g]=="loglinear"){
        if(!is.null(att.str)) stop("Only att.dist = saturated or fixed is available when attributes have structures.",call. = FALSE)
        Z <- designM(K,0)
        if(K<2)
          stop("loglinear smoothing is not available when K < 2.",call. = FALSE)
        if(loglinear[g]>K)
          stop("Argument 'loglinear' cannot be greater than K.",call. = FALSE)
        Z <- Z[,seq_len(1+sum(sapply(seq_len(loglinear[g]),choose,n=K)))]
        prior[prior[,g]<1e-9,g] <- 1e-9
        lambda[[g]] <- parm <- stats::lm.wfit(x=Z,log(N*prior[,g]),w=N*prior[,g])$coefficients
        logprior[,g] <- c(Z%*%parm)-log(sum(exp(Z%*%parm)))
        logprior[logprior[,g]<log(.Machine$double.eps),g] <- log(.Machine$double.eps)
        logprior[logprior[,g]>log(1-.Machine$double.eps),g] <- log(1-.Machine$double.eps)

        npar <- npar + length(lambda[[g]])
      }else if(att.dist[g]=="independent"){
        if(!is.null(att.str))
          stop("Only att.dist = saturated or fixed is available when attributes have structures.",call. = FALSE)
        lambda[[g]] <- pk <- colSums(rowProd(prior,Ng)[,g]*AlphaPattern)/Ng[g] #length of K
        logprior[,g] <- AlphaPattern%*%log(pk) + (1-AlphaPattern)%*%log(1-pk)
        npar <- npar + length(lambda[[g]])
      }else if(att.dist[g]=="fixed"){
        logprior[,g] <- initial.logprior[,g]
        lambda[[g]] <- NA
      }
    }
  }





  return(list(logprior=logprior,lambda=lambda,higher.order = ho,npar=npar))
}



structural.parm.sg <- function(AlphaPattern,no.mg = 1,logprior,att.dist,att.str,saturated,initial.logprior,
                               lower.prior,loglinear=2,N,K,lambda,higher.order){
  prior <- c(exp(logprior))
  parm <- NULL
  ho <- NULL

  if (att.dist=="saturated"){

    prior <- (prior*N+saturated$prior)/sum(prior*N+saturated$prior)
    prior[which(prior<lower.prior)] <- lower.prior
    lambda <- prior <- prior/sum(prior)
    logprior <- log(prior)
    npar <- length(lambda) - 1
  }else if(att.dist=="fixed"){
    logprior <- initial.logprior
    lambda <- NA
    npar <- 0
  }else if (att.dist=="higher.order")
  {
    if(!is.null(att.str)) stop("Only att.dist = saturated or fixed is available when attributes have structures.",call. = FALSE)
    HO.out <- SG.HO.est(lambda=lambda,AlphaPattern = AlphaPattern, HOgr = 1, Rl = N*prior,
                     higher.order = higher.order)
    logprior <- HO.out$logprior
    lambda <- HO.out$lambda
    ho <- HO.out$higher.order
    npar <- HO.out$npar
  }else if (att.dist=="loglinear"){
    if(!is.null(att.str)) stop("Only att.dist = saturated or fixed is available when attributes have structures.",call. = FALSE)
    Z <- designM(K,0)
    if(K<2) stop("loglinear smoothing is not available when K < 2.",call. = FALSE)
    if(loglinear>K) stop("Argument 'loglinear' cannot be greater than K.",call. = FALSE)
    Z <- Z[,seq_len(1+sum(sapply(seq_len(loglinear),choose,n=K)))]
    prior[prior<1e-9] <- 1e-9
    lambda <- parm <- stats::lm.wfit(x=Z,log(N*prior),w=N*prior)$coefficients
    logprior <- c(Z%*%parm)-log(sum(exp(Z%*%parm)))
    logprior[logprior<log(.Machine$double.eps)] <- log(.Machine$double.eps)
    logprior[logprior>log(1-.Machine$double.eps)] <- log(1-.Machine$double.eps)
    npar <- length(lambda)
  }else if(att.dist=="independent"){
    if(!is.null(att.str)) stop("Only att.dist = saturated or fixed is available when attributes have structures.",call. = FALSE)
    # print(AlphaPattern)
    # print(prior)
    lambda <- pk <- colSums(N*prior*AlphaPattern)/N #length of K
    logprior <- AlphaPattern%*%log(pk) + (1-AlphaPattern)%*%log(1-pk)
    npar <- length(lambda)
  }
  return(list(logprior=logprior,lambda=lambda,higher.order = ho,npar=npar))
}
