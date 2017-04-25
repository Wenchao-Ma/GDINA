
HO.est <- function(Rl, K, N, IRTmodel = "Rasch", theta.range = c(-3,3),type="testwise",
                   nnodes = 19, a = NULL, b = NULL, loga.prior = c(0,0.25),b.prior = c(0,1),
                   method = "MMLE", a.bound = c(0.1, 5), b.bound = c(-4, 4))
{
  Aq <- seq(theta.range[1], theta.range[2], length.out = nnodes)
  Wght <- dnorm(Aq)
  Wght <- Wght/sum(Wght)
  if (is.null(a)) a <- rep(1, K)
  X <- alpha(K)
  if (is.null(b)) b <- (-1)*qnorm(0.975 - 0.95 * colSums(X * matrix(rep(Rl, K), ncol = K))/N)
  if (toupper(method)=="BL"){

    obj_fn_2PL <- function(vpar)
    {
      # a1,a2,..aJ,b1,b2,.. bJ
      LXl <- t(HO.loglik(vpar[1:(length(vpar)/2)], vpar[(1 + (length(vpar)/2)):length(vpar)],
                       Aq,nnodes, X,K))
      -1 * sum(Rl * log(colSums(exp(matrix(rep(log(Wght), 2^K), nrow = nnodes) +
                                      LXl))))
    }
    obj_fn_Rasch <- function(vpar)
    {
      LXl <- t(HO.loglik(rep(1,length(vpar)), vpar, Aq,nnodes, X,K))
      -1 * sum(Rl * log(colSums(exp(matrix(rep(log(Wght), 2^K), nrow = nnodes) +
                                      LXl))))
    }
    obj_fn_1PL <- function(vpar)
    {
      LXl <- t(HO.loglik(rep(vpar[1],length(vpar)-1), vpar[2:length(vpar)], Aq,nnodes, X,K))
      -1 * sum(Rl * log(colSums(exp(matrix(rep(log(Wght), 2^K), nrow = nnodes) +
                                      LXl))))
    }
    if (IRTmodel=="2PL"){
      vpar <- c(a, b)
      HOopt <- optim(vpar, fn = obj_fn_2PL, method = "L-BFGS-B",
                     lower = c(rep(a.bound[1], K), rep(b.bound[1], K)),
                     upper = c(rep(a.bound[2], K), rep(b.bound[2], K)))
      lambda = data.frame(slope=HOopt$par[1:(length(vpar)/2)],
                          intercept=HOopt$par[(1 + (length(vpar)/2)):length(vpar)])
    }else if (IRTmodel=="Rasch"){
      vpar <- b
      HOopt <- optim(vpar, fn = obj_fn_Rasch, method = "L-BFGS-B",
                     lower = c(rep(b.bound[1], K)),
                     upper = c(rep(b.bound[2], K)))
       lambda = data.frame(slope=rep(1,K), intercept=HOopt$par)
    }else if (IRTmodel=="1PL"){
      vpar <- c(a[1],b)
      HOopt <- optim(vpar, fn = obj_fn_1PL, method = "L-BFGS-B",
                     lower = c(a.bound[1],rep(b.bound[1], K)),
                     upper = c(a.bound[2],rep(b.bound[2], K)))
      lambda = data.frame(slope=rep(HOopt$par[1],K), intercept=HOopt$par[2:length(vpar)])
    }
  }else{
    HOlike <- exp(HO.loglik(a,b,Aq,nnodes,X,K)) #2^K x nnodes - P(alpha_c|Aq,a,b)
    Nlq <- c(Rl)*HOlike * matrix(Wght, 2^K, nnodes, byrow = T)/rowSums(HOlike * matrix(Wght, 2^K, nnodes, byrow = T))
    # Nq <- colSums(post/rowSums(post)) # length of nnodes - expected number of examinees with Aq - P(Xi|alpha_c)P(alpha_c|Aq,a,b)*p(Aq)/P(Xi)
    Nq <- colSums(Nlq)
    #Compute Rqk (each column) - number of examinees at each node that have mastered attribute k - nnodes x K
    Rq <- (t(Nlq)%*%X)

    lambda <- matrix(1,K,2)
    if (method=="BMLE") prior <- TRUE else if (method=="MMLE") prior <- FALSE
    if(type=="attwise"){
      if(toupper(IRTmodel)=="1PL"){
        for(k in 1:K)   b[k] <- lambda[k,2] <- optimise(f=obj_fn_b,interval = b.bound,
                                                        Nq = Nq, Rqk = Rq[,k], Aq=Aq,
                                                        mu=b.prior[1], sigma2=b.prior[2],
                                                        a = a[k],prior=prior)$minimum
        lambda[,1] <- exp(optimise(f=obj_fn_1PLa,interval = log(a.bound),Nq = Nq, Rq = Rq,Aq=Aq,
                                   mu=loga.prior[1], sigma2=loga.prior[2], b = b,prior=prior)$minimum)


      }else if(toupper(IRTmodel)=="2PL"){
        # HO.2PL.j <- function(vP, Nq, Rqk,Aq,mu.loga,mu.b,sigma2.loga,sigma2.b,prior=FALSE)
        for(k in 1:K)  {
          opt <-optim(c(log(a[k]),b[k]), HO.2PL.j, gr = NULL,
                      lower = c(log(a.bound[1]), b.bound[1]),
                      upper = c(log(a.bound[2]),b.bound[2]),
                      method = "L-BFGS-B", control=list(fnscale=-1),
                      Nq=Nq, Rqk=Rq[,k],Aq=Aq,mu.loga=loga.prior[1],sigma2.loga=loga.prior[2],
                      mu.b=b.prior[1],sigma2.b=b.prior[2],prior=prior)$par
          lambda[k,] <- c(exp(opt[1]),opt[2])

        }

      }else if (IRTmodel=="Rasch"){
        for(k in 1:K)  lambda[k,2] <- optimise(f=gr.Rasch.j,interval = b.bound,
                                               Nq = Nq, Rqk = Rq[,k],Aq=Aq,
                                               mu=b.prior[1], sigma2=b.prior[2],
                                               a = 1,prior = prior)$minimum
      }
    }else if(type=="testwise"){
      if(toupper(IRTmodel)=="2PL"){
        HOopt<-optim(c(log(a),b), HO.2PL, gr = NULL,
                     lower = c(rep(log(a.bound[1]), K), rep(b.bound[1], K)),
                     upper = c(rep(log(a.bound[2]), K), rep(b.bound[2], K)),
                     method = "L-BFGS-B", control=list(fnscale=-1),
                     Nq=Nq, Rq=Rq,Aq=Aq,mu.loga=loga.prior[1],sigma2.loga=loga.prior[2],
                     mu.b=b.prior[1],sigma2.b=b.prior[2],prior=prior)$par
        lambda[,1] <- exp(HOopt[1:(length(HOopt)/2)])
        lambda[,2] <- HOopt[(1 + (length(HOopt)/2)):length(HOopt)]
      }else if (IRTmodel=="Rasch"){
        lambda[,2] <- optim(b, HO.Rasch, gr = NULL,
                            lower = rep(b.bound[1], K),upper = rep(b.bound[2], K),
                            method = "L-BFGS-B", control=list(fnscale=-1),
                            Nq=Nq, Rq=Rq,Aq=Aq,
                            mu=b.prior[1],sigma2=b.prior[2],prior=prior)$par
      }else if(IRTmodel=="1PL"){
        HOopt <- optim(c(log(a[1]),b), HO.1PL, gr = NULL,
                       lower = c(log(a.bound[1]),rep(b.bound[1], K)),upper = c(log(a.bound[2]),rep(b.bound[2], K)),
                       method = "L-BFGS-B", control=list(fnscale=-1),
                       Nq=Nq, Rq=Rq,Aq=Aq,mu.loga=loga.prior[1],sigma2.loga=loga.prior[2],
                       mu.b=b.prior[1],sigma2.b=b.prior[2],prior=prior)$par
        lambda[,1] <- exp(HOopt[1])
        lambda[,2] <- HOopt[-1]
      }

    }
  }

  upd.Lx <- HO.loglik(lambda[,1],lambda[,2], Aq,nnodes, X,K) #2^K x nnodes - P(alpha_c|Aq,a,b)

  logprior <- log(rowSums(exp(matrix(log(Wght), nrow = 2^K, ncol = nnodes,byrow = TRUE) + upd.Lx)))

  return(list(logprior=logprior,lambda=lambda,HO.loglik=upd.Lx))
}
