
Mstep_DTM <- function(delta,R0,Iestpar,K,item.no,item.no.0,type="cumulative",linkfunc="identity",eq.const = FALSE,Qc=NULL,Tmatrix=NULL){
  dsg <- GDINA::designmatrix(K)
  L <- 2^K
  delta2 <- delta
  opt <- vector("list",length(unique(item.no)))
  for(j in unique(item.no)){
    mIndj <- Iestpar[which(item.no==j),,drop=FALSE]
    mparj <- delta[which(item.no==j),,drop=FALSE]
    ncatj <- sum(item.no==j)
    if(!eq.const){
      Indj <- c(t(mIndj))
      parj <- c(t(mparj))[which(Indj==1)]
      if(linkfunc=="logit"){
        if(type!="cumulative"){#adjacent or sequential
          opt[[j]] <- optim(par = parj, fn = obj,method = "L-BFGS-B",
                            Indj = Indj, ncat = ncatj,lower = -9.21,upper = 9.21,
                            R0j = R0[which(item.no.0==j),,drop=FALSE],
                            dsg = dsg,control = list(trace=0),type=type,linkfunc=linkfunc,
                            Qcj=Qc[which(item.no==j),,drop=FALSE],Tmatrixj=Tmatrix[[j]])
          parj <- Indj
          parj[which(Indj==1)] <- opt[[j]]$par
          delta2[which(item.no==j),] <- matrix(parj,nrow = ncatj,byrow = TRUE)
        }else if(type=="cumulative"){
          opt[[j]] <- solnp(pars = parj, fun = obj,
                            ineqfun = hin,
                            ineqLB = rep(0.0001,L*(ncatj+1)),ineqUB = rep(0.9999,L*(ncatj+1)),
                            # ineqLB = rep(0,L*(ncatj*2+1)),ineqUB = rep(1,L*(ncatj*2+1)),
                            Indj = Indj, ncat = ncatj,
                            R0j = R0[which(item.no.0==j),,drop=FALSE],
                            dsg = dsg,control = list(trace=0),type=type,linkfunc=linkfunc)
          parj <- Indj
          parj[which(Indj==1)] <- opt[[j]]$pars
          delta2[which(item.no==j),] <- matrix(parj,nrow = ncatj,byrow = TRUE)
        }
      }else if(linkfunc=="identity"){
        opt[[j]] <- solnp(pars = parj, fun = obj,
                          ineqfun = hin, ineqLB = rep(0,L*(ncatj*2+1)),ineqUB = rep(1,L*(ncatj*2+1)),
                          Indj = Indj, ncat = ncatj,
                          R0j = R0[which(item.no.0==j),,drop=FALSE],
                          dsg = dsg,control = list(trace=0),type=type,linkfunc=linkfunc)
        parj <- Indj
        parj[which(Indj==1)] <- opt[[j]]$pars
        delta2[which(item.no==j),] <- matrix(parj,nrow = ncatj,byrow = TRUE)
      }
    }else if(eq.const){
      parj <- parmtrans(mparj,mIndj)
      if(linkfunc=="logit"){
        if(type!="cumulative"){#adjacent or sequential
          opt[[j]] <- optim(par = parj, fn = obj_eq,method = "BFGS",
                            mIndj = mIndj, ncat = ncatj,
                            R0j = R0[which(item.no.0==j),,drop=FALSE],
                            dsg = dsg,control = list(trace=0),type=type,linkfunc=linkfunc)
          delta2[which(item.no==j),] <- t(invparmtrans(opt[[j]]$par,mIndj))
        }else if(type=="cumulative"){

          if(ncatj==1){
            opt[[j]] <- optim(par = parj, fn = obj_eq,method = "BFGS",
                              mIndj = mIndj, ncat = ncatj,
                              R0j = R0[which(item.no.0==j),,drop=FALSE],
                              dsg = dsg,control = list(trace=0),type=type,linkfunc=linkfunc)
            delta2[which(item.no==j),] <- t(invparmtrans(opt[[j]]$par,mIndj))
          }else{
            opt[[j]] <- alabama::auglag(par = parj, fn = obj_eq,hin = hin_eq,
                                        mIndj = mIndj, ncat = ncatj,
                                        R0j = R0[which(item.no.0==j),,drop=FALSE],
                                        dsg = dsg,type=type,linkfunc=linkfunc,
                                        control.outer = list(trace=FALSE,kkt2.check=FALSE))
            delta2[which(item.no==j),] <- t(invparmtrans(opt[[j]]$par,mIndj))
            # print(opt[[j]]$par)
          }

        }
      }else if(linkfunc=="identity"){
        opt[[j]] <- solnp(pars = parj, fun = obj,
                          ineqfun = hin, ineqLB = rep(0,L*(ncatj*2+1)),
                          ineqUB = rep(1,L*(ncatj*2+1)),
                          eqfun = heq,eqB = rep(0,sum(rowSums(matrix(Indj,ncol = ncatj))-1)),
                          Indj = Indj, ncat = ncatj,
                          R0j = R0[which(item.no.0==j),,drop=FALSE],
                          dsg = dsg,control = list(trace=0),type=type,linkfunc=linkfunc)
        parj <- Indj
        parj[which(Indj==1)] <- opt[[j]]$pars
        delta2[which(item.no==j),] <- matrix(parj,nrow = ncatj,byrow = TRUE)
      }
    }

  }
  return(list(delta=delta2,opt=opt))
}



obj <- function(vd,Indj,ncat,R0j,dsg,type,linkfunc,Tmatrixj=NULL,Qcj=NULL){
  d <- Indj
  d[which(Indj==1)] <- vd
  pi <- apply(matrix(d,ncol = ncat),2,function(x)dsg%*%x)
  # pi <- cbind(1-rowSums(pi),pi) #including cat 0 in the first col
  pj <- v2pj(t(pi),type=type,linkfunc=linkfunc,Qcj=Qcj,Tmatrixj=Tmatrixj) #Sj0 x L
  # print(pj)
  pj[pj < .Machine$double.eps] <- .Machine$double.eps
  pj[pj > (1 - .Machine$double.eps)] <- 1 - .Machine$double.eps
  # print(pj)
  -1*sum(R0j*log(pj))
}

#objective function under equality constr.
obj_eq <- function(vd,mIndj,ncat,R0j,dsg,type=type,linkfunc=linkfunc){
  v <- invparmtrans(vd,mIndj)
  pi <- apply(v,2,function(x)dsg%*%x)
  pj <- v2pj(t(pi),type=type,linkfunc=linkfunc) #Sj0 x L
  pj[pj < .Machine$double.eps] <- .Machine$double.eps
  pj[pj > (1 - .Machine$double.eps)] <- 1 - .Machine$double.eps
  -1*sum(R0j*log(pj))
}


hin_eq <- function(vd,mIndj,ncat,R0j,dsg,type=type,linkfunc=linkfunc){
  v <- invparmtrans(vd,mIndj)
  pi <- apply(v,2,function(x)dsg%*%x)
  pj <- c(v2pj(t(pi),type=type,linkfunc=linkfunc)) #Sj0 x L
  return(c(pj,1-pj)) # for logit link
}

hin <- function(vd,Indj,ncat,R0j,dsg,type=type,linkfunc=linkfunc){
  d <- Indj
  d[which(Indj==1)] <- vd
  pi <- apply(matrix(d,ncol = ncat),2,function(x)dsg%*%x)
  # pi <- cbind(1-rowSums(pi),pi) #including cat 0 in the first col
  pj <- c(v2pj(t(pi),type=type,linkfunc=linkfunc)) #Sj0 x L
  # print(pj)
  if (linkfunc=="logit") {
    return(pj)
  }else{
    return(c(pi,pj))
  }
  #return(pj)
}
# only needed when there are more than 1 cat
heq <- function(vd,Indj,ncat,R0j,dsg){
  d <- Indj

  d[which(Indj==1)] <- 1:sum(Indj)
  dm <- matrix(d,ncol = ncat)
  out <- NULL
  for(r in 1:nrow(dm)){
    dmr <- dm[r,]
    dmr <- dmr[which(dmr>0)]
    if(length(dmr)>1) out <- c(out,diff(vd[dmr]))
  }
  out
}
