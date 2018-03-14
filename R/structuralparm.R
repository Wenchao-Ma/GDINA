structural.parm <- function(AlphaPattern,no.mg,logprior,att.dist,att.str,saturated,initial.logprior,
                            lower.prior,loglinear=2,N,Ng,K,lambda,higher.order){
  prior <- exp(logprior)
  parm <- NULL
  ho <- NULL
  if (all(att.dist=="higher.order"))
  {
    if(no.mg==1) higher.order$Type <- "SameTheta"
    HO.out <- HO.est(lambda=lambda,AlphaPattern = AlphaPattern, HOgr = seq_len(no.mg), Rl = rowProd(prior,Ng),
                     higher.order = higher.order)
    logprior <- HO.out$logprior
    lambda <- HO.out$lambda
    ho <- HO.out$higher.order
  }else{
    for(g in seq_len(no.mg)){
      if (att.dist[g]=="saturated"){
        prior[,g] <- (prior[,g]*Ng[g]+saturated$prior[,g])/sum(prior[,g]*Ng[g]+saturated$prior[,g])
        prior[which(prior[,g]<lower.prior[g]),g] <- lower.prior[g]
        prior[,g] <- prior[,g]/sum(prior[,g])
        logprior[,g] <- log(prior[,g])
        lambda[[g]] <- prior[,g]
      }else if (att.dist[g]=="loglinear"){
        Z <- designM(K,0)
        if(K<2) stop("loglinear smoothing is not available when K < 2.",call. = FALSE)
        if(loglinear[g]>K) stop("Argument 'loglinear' cannot be greater than K.",call. = FALSE)
        Z <- Z[,seq_len(1+sum(sapply(seq_len(loglinear[g]),choose,n=K)))]
        prior[prior[,g]<1e-9,g] <- 1e-9
        lambda[[g]] <- parm <- stats::lm.wfit(x=Z,log(N*prior[,g]),w=N*prior[,g])$coefficients
        logprior[,g] <- c(Z%*%parm)-log(sum(exp(Z%*%parm)))
        logprior[logprior[,g]<log(.Machine$double.eps),g] <- log(.Machine$double.eps)
        logprior[logprior[,g]>log(1-.Machine$double.eps),g] <- log(1-.Machine$double.eps)
      }else if(att.dist[g]=="higher.order"){
        HO.out <- HO.est(lambda=lambda,AlphaPattern = AlphaPattern, HOgr = g, Rl = rowProd(prior,Ng),
                         higher.order = higher.order)
        logprior[,g] <- HO.out$logprior
        lambda[[g]] <- HO.out$lambda[[g]]
      }else if(att.dist[g]=="independent"){
        lambda[[g]] <- pk <- colSums(rowProd(prior,Ng)[,g]*AlphaPattern)/Ng[g] #length of K
        logprior[,g] <- AlphaPattern%*%log(pk) + (1-AlphaPattern)%*%log(1-pk)
      }else if(att.dist[g]=="fixed"){
        logprior[,g] <- initial.logprior[,g]
        lambda <- NULL
      }
    }
  }


  return(list(logprior=logprior,lambda=lambda,higher.order = ho))
}
