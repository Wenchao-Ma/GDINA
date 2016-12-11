
HO.est <- function(Rl, loglik, K, N, IRTmodel = "Rasch", theta.range = c(-3,3),
                   nnodes = 19, a = NULL, b = NULL,
                   method = "MMLE", a.bound = c(0.1, 5), b.bound = c(-4, 4))
{
  Aq <- seq(theta.range[1], theta.range[2], length.out = nnodes)
  Wght <- dnorm(Aq)
  Wght <- Wght/sum(Wght)
  if (is.null(a)) a <- rep(1, K)
  X <- alpha(K)
  if (is.null(b)) b <- (-1)*qnorm(0.975 - 0.95 * colSums(X * matrix(rep(Rl, K), ncol = K))/N)
if(IRTmodel=="1PL") method <- "BL"
  if (toupper(method)=="BL"){

    obj_fn_2PL <- function(vpar)
    {
      # a1,a2,..aJ,b1,b2,.. bJ
      LXl <- HO.loglik(vpar[1:(length(vpar)/2)], vpar[(1 + (length(vpar)/2)):length(vpar)],
                  Aq,nnodes, X,K)
      -1 * sum(Rl * log(colSums(exp(matrix(rep(log(Wght), 2^K), nrow = nnodes) +
                                      LXl))))
    }
    obj_fn_Rasch <- function(vpar)
    {
      LXl <- HO.loglik(rep(1,length(vpar)), vpar, Aq,nnodes, X,K)
      -1 * sum(Rl * log(colSums(exp(matrix(rep(log(Wght), 2^K), nrow = nnodes) +
                                      LXl))))
    }
    obj_fn_1PL <- function(vpar)
    {
      LXl <- HO.loglik(rep(vpar[1],length(vpar)-1), vpar[2:length(vpar)], Aq,nnodes, X,K)
      -1 * sum(Rl * log(colSums(exp(matrix(rep(log(Wght), 2^K), nrow = nnodes) +
                                      LXl))))
    }
    if (IRTmodel=="2PL"){
      vpar <- c(a, b)
      HOopt <- optim(vpar, fn = obj_fn_2PL, method = "L-BFGS-B",
                     lower = c(rep(a.bound[1], K), rep(b.bound[1], K)),
                     upper = c(rep(a.bound[2], K), rep(b.bound[2], K)))
      # upd.Lx <- P2PL(HOopt$par[1:(length(vpar)/2)], HOopt$par[(1 + (length(vpar)/2)):length(vpar)],
      #                Aq)
      lambda = data.frame(slope=HOopt$par[1:(length(vpar)/2)],
                          intercept=HOopt$par[(1 + (length(vpar)/2)):length(vpar)])
    }else if (IRTmodel=="Rasch"){
      vpar <- b
      HOopt <- optim(vpar, fn = obj_fn_Rasch, method = "L-BFGS-B",
                     lower = c(rep(b.bound[1], K)),
                     upper = c(rep(b.bound[2], K)))
      # upd.Lx <- P2PL(rep(1,length(b)), HOopt$par, Aq)
      lambda = data.frame(slope=rep(1,K), intercept=HOopt$par)
    }else if (IRTmodel=="1PL"){
      vpar <- c(a[1],b)
      HOopt <- optim(vpar, fn = obj_fn_1PL, method = "L-BFGS-B",
                     lower = c(a.bound[1],rep(b.bound[1], K)),
                     upper = c(a.bound[2],rep(b.bound[2], K)))
      # upd.Lx <- P2PL(rep(1,length(b)), HOopt$par, Aq)
      lambda = data.frame(slope=rep(HOopt$par[1],K), intercept=HOopt$par[2:length(vpar)])
    }
   }
  else if(toupper(method)=="MMLE"){

    HOlike <- exp(t(HO.loglik(a,b,Aq,nnodes,X,K))) #L x nnodes
    postNq <- exp(loglik) %*% HOlike * matrix(Wght, N, nnodes, byrow = T)
    Nq <- colSums(postNq/rowSums(postNq))

    #Compute Rqk - number of examinees at each node that have mastered attribute k
    Rq <- apply(X,2,function(x){
      pq <- (exp(loglik) * matrix(x, N, 2^K, byrow = T)) %*% HOlike * matrix(Wght, N, nnodes, byrow = T)
      colSums(pq/rowSums(postNq))
    })

    if(toupper(IRTmodel)=="2PL"){
      lambda <- t(apply(rbind(a,b,Rq),2,function(x){
        HOopt<-optim(x[1:2], HO.2PL, gr = NULL,
                     lower = c(a.bound[1],b.bound[1]),upper = c(a.bound[2],b.bound[2]),
                     method = "L-BFGS-B", control=list(fnscale=-1),
                     Nq, Rqk=x[3:(nnodes+2)],Aq=Aq)
        HOopt$par
      }))
    }else if (IRTmodel=="Rasch"){
      a <- 1
      lambda <- cbind(a,apply(rbind(a,b,Rq),2,function(x){
        HOopt<-optim(x[2], HO.Rasch, gr = NULL,
                     lower = b.bound[1],upper = b.bound[2],
                     method = "L-BFGS-B", control=list(fnscale=-1),
                     Nq, Rqk=x[3:(nnodes+2)],Aq=Aq)
        HOopt$par
      }))
    }

  }

  upd.Lx <- HO.loglik(lambda[,1],lambda[,2], Aq,nnodes, X,K)

  logprior <- log(colSums(exp(matrix(rep(log(Wght), 2^K), nrow = nnodes) +
                                upd.Lx)))

  return(list(logprior=logprior,lambda=lambda,HO.loglik=t(upd.Lx)))
}
