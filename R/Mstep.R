
Mstep <- function(delta,item.parm,J,Kj,itr,Rg,Ng,model,ConstrType,linkfunc,
                  correction,lower.p,upper.p,algorithm,ConstrMatrix,solver, item.prior,
                  DesignMatrices,MstepMessage, est.bin,auglag_args,nloptr_args,solnp_args){

opts <- list()
  for(j in which(est.bin)){ #for each item

    optims <- GNLOptim_call(par = delta[[j]],solver=solver[j],modelj=model[j],correction=correction,
             auglag_args=auglag_args,nloptr_args=nloptr_args,solnp_args = solnp_args,item.prior=item.prior,
             Nj=Ng[j,1:2^Kj[j]],Rj=Rg[j,1:2^Kj[j]],designMj=DesignMatrices[[j]],
             uPj=upper.p[j],lPj=lower.p[j],MstepMessage=MstepMessage,itr=itr,j=j,
             ConstrMatrix=ConstrMatrix[[j]],linkfunc=linkfunc[j],eps=1e-16,ConstrType=ConstrType[j])
    if(is.null(optims)) stop(paste("M-step fails for item",j,"at iteration",itr),call. = FALSE)

     delta[[j]] <- optims$delta
     item.parm[j,1:(2^Kj[j])] <- optims$phat
    opts[[j]] <- optims
  }

  return(list(item.parm=item.parm,delta=delta,optims=opts))
}

GNLOptim <- function(par,fn,gr=NULL,hin=NULL,hin.gr=NULL,optimizer=NULL,
                     modelj,correction,Nj,Rj,designMj,uPj,lPj,itr,j,auglag_args,nloptr_args,solnp_args,
                     ConstrMatrix,linkfunc,eps,ConstrType,MstepMessage,...){
  ret <- list()

  if(optimizer=="ClosedForm"){
    if(linkfunc!=1) stop("Closed-form solutions are only available for the identity link GDINA, DINA and DINO models.",call. = FALSE)
    if(modelj==0){
      if (any(Nj<correction[2])){
        Nj[which(Nj<correction[2])] <- Nj[which(Nj<correction[2])] + correction[2]
        Rj[which(Nj<correction[2])] <- Rj[which(Nj<correction[2])] + correction[1]
      }
      phat <- c(Rj/Nj)
    }else if(modelj==1 | modelj==2){
      rNj <- c(rep(sum(Nj[designMj[,2]==0],na.rm = TRUE),sum(designMj[,2]==0)),
               rep(sum(Nj[designMj[,2]==1],na.rm = TRUE),sum(designMj[,2]==1)))
      rRj <- c(rep(sum(Rj[designMj[,2]==0],na.rm = TRUE),sum(designMj[,2]==0)),
               rep(sum(Rj[designMj[,2]==1],na.rm = TRUE),sum(designMj[,2]==1)))
      if (any(rNj<correction[2])){
        rNj[which(rNj<correction[2])] <- rNj[which(rNj<correction[2])] + correction[2]
        rRj[which(rNj<correction[2])] <- rRj[which(rNj<correction[2])] + correction[1]
      }
      phat <- c(rRj/rNj)
    }else{
      stop("Closed-form solutions are only available for the identity link GDINA, DINA and DINO models.",call. = FALSE)
    }

    phat[phat<lPj] <- lPj
    phat[phat>uPj] <- uPj
    ret$delta <- c(Calc_Dj(phat, designMj, linkfunc = linkfunc, boundary = 1))
    ret$opt <- list(convergence=0,message="converged - closed-form solution.")
    ret$phat <- phat
    ret$convergence <- 0
  }else if(optimizer=="BFGS"){
    if(modelj>=3&&modelj<=5){
      par <- initials_optim(par,modelj=modelj,Nj=Nj,Rj=Rj,designMj=designMj,uPj=uPj,ConstrMatrix=ConstrMatrix,
                                 lPj=lPj,linkfunc=linkfunc,ConstrType=ConstrType)
       # print(par)
    }
    Lj <- nrow(designMj)
    if(ConstrType==1){
      ui <- rbind(designMj,-1*designMj)
      cii <- rep(c(lPj,-1*uPj),each=Lj)
      ci <- list(gdina=cii,dina=cii, dino=cii, acdm=cii,
                 llm=rep(c(qlogis(lPj),-1*qlogis(uPj)),each=Lj),
                 rrum=rep(c(log(lPj),-1*log(uPj)),each=Lj))
    }else if(ConstrType==3){
      if(modelj>0){
        ui=rbind(designMj,-1*designMj,diag(length(par)))
      }else if(modelj==0){
        ui=ConstrMatrix%*%designMj
      }

      cii <- c(rep(c(lPj,-1*uPj),each=Lj),c(lPj,rep(0,length(par)-1)))
      ci <- list(gdina=rep(0,nrow(ConstrMatrix)),dina=cii,dino=cii,acdm=cii,
                      llm=c(rep(c(qlogis(lPj),-1*qlogis(uPj)),each=Lj),c(qlogis(lPj),rep(0,length(par)-1))),
                      rrum=c(rep(c(log(lPj),-1*log(uPj)),each=Lj),c(log(lPj),rep(0,length(par)-1))))
    }

    op <- tryCatch(stats::constrOptim(theta = par,f = fn,grad = gr,method = "BFGS",
                             ui=ui,ci=ci[[modelj+1]],control=list(trace=0),
                             Nj=Nj,Rj=Rj,designMj=designMj,uPj=uPj,ConstrMatrix=ConstrMatrix,
                             lPj=lPj,linkfunc=linkfunc,eps=eps,ConstrType=ConstrType),
                   error=function(e){
                     if(MstepMessage){
                       message(paste("M-step optimization failed for item",j,"at iteration",itr))
                       message("Error message from stats::constrOptim:")
                       message(e)
                     }
                     return("error")
                   })

    if(length(op)==1 && op=="error"){
        ret$convergence <- -1
      }else{
      ret$delta <- op$par
      ret$opt <- op
      ret$phat <- c(Calc_Pj(op$par,designMj,linkfunc,boundary = 0))
      if(any(ret$phat<0)||any(ret$phat>1)){
        warning(paste("Estimates are out of bounds for item",j,"at iteration",itr,"[stats::constrOptim]"),call. = FALSE)
        ret$convergence <- -3
      }else{
        ret$convergence <- 0
      }
    }
# print(ret$convergence)
  }else if(optimizer=="auglag"){
    args.list <- c(list(par=par,fn=fn,gr=gr,hin = hin,hin.jac = hin.gr,
                      Nj=Nj,Rj=Rj,designMj=designMj,uPj=uPj,ConstrMatrix=ConstrMatrix,
                      lPj=lPj,linkfunc=linkfunc,eps=eps,ConstrType=ConstrType),auglag_args)
    op <- tryCatch(do.call(alabama::auglag,args.list),
                   error=function(e){
                     if(MstepMessage){
                       message(paste("M-step optimization failed for item",j,"at iteration",itr))
                       message("Error message from alabama::auglag:")
                       message(e)
                     }
                     return("error")
                   })
    if(length(op)==1&&op=="error"){
        ret$convergence <- -1
      }else{
      ret$delta <- op$par
      ret$opt <- op
      ret$phat <- c(Calc_Pj(op$par,designMj,linkfunc,boundary = 0))
      if(any(ret$phat<0)||any(ret$phat>1)){
        warning(paste("Estimates are out of bounds for item",j,"at iteration",itr,"[alabama::auglag]"),call. = FALSE)
        ret$convergence <- -3
      }else{
        ret$convergence <- 0
      }
    }
  }else if(optimizer=="solnp"){
    if(ConstrType==1){nineq <- 2*nrow(designMj)}else if(ConstrType==3){nineq <- 2*nrow(designMj)+nrow(ConstrMatrix)}
    args.list <- list(pars=par,fun=fn,ineqfun = hin,ineqLB = rep(0,nineq),
                      ineqUB = rep(Inf,nineq),control=solnp_args,
                      Nj=Nj,Rj=Rj,designMj=designMj,uPj=uPj,ConstrMatrix=ConstrMatrix,
                      lPj=lPj,linkfunc=linkfunc,eps=eps,ConstrType=ConstrType)
    # print(args.list)
    op <- tryCatch(do.call(Rsolnp::solnp,args.list),
                   error=function(e){
                     if(MstepMessage){
                       message(paste("M-step optimization failed for item",j,"at iteration",itr,"\n"))
                       message("Error message from Rsolnp::solnp:\n")
                       message(e)
                     }
                     return("error")
                   })

    if(length(op)==1&&op=="error"){
      ret$convergence <- -1
    }else{
      print(op$pars)
      ret$delta <- op$pars
      ret$opt <- op
      ret$phat <- c(Calc_Pj(op$pars,designMj,linkfunc,boundary = 0))
      if(any(ret$phat<0)||any(ret$phat>1)){
        message(paste("Estimates are out of bounds for item",j,"at iteration",itr,"[Rsolnp:solnp]"),call. = FALSE)
        ret$convergence <- -3
      }else{
        ret$convergence <- 0
      }
    }

  }else if(optimizer=="slsqp"){

    op <- tryCatch(nloptr::nloptr(x0=par,eval_f=fn,eval_grad_f=gr,eval_g_ineq = hin,eval_jac_g_ineq = hin.gr,
                                  Nj=Nj,Rj=Rj,designMj=designMj,uPj=uPj,ConstrMatrix=ConstrMatrix,greaterthan0 = FALSE,
                                  lPj=lPj,linkfunc=linkfunc,eps=eps,ConstrType=ConstrType,
                                  opts=nloptr_args,...),
                   error=function(e){
                     if(MstepMessage){
                       message(paste("M-step optimization failed for item",j,"at iteration",itr))
                       message("Error message from nloptr::slsqp:")
                       message(e)
                     }
                     return("error")
                   })
    if(length(op)==1&&op=="error"){
      ret$convergence <- -1
    }else{
      ret$delta <- op$solution
      ret$opt <- op
      ret$phat <- c(Calc_Pj(op$solution,designMj,linkfunc,boundary = 0))
      if(any(ret$phat<0)||any(ret$phat>1)){
        message(paste("Estimates are out of bounds for item",j,"at iteration",itr,"[nloptr:slsqp]"),call. = FALSE)
        ret$convergence <- -3
      }else{
        ret$convergence <- 0
      }
    }
  }else if(optimizer=="nloptr"){

    op <- tryCatch(nloptr::nloptr(x0=par,eval_f=fn,eval_grad_f=gr,eval_g_ineq = hin,eval_jac_g_ineq = hin.gr,
                                  Nj=Nj,Rj=Rj,designMj=designMj,uPj=uPj,ConstrMatrix=ConstrMatrix,greaterthan0 = FALSE,
                                  lPj=lPj,linkfunc=linkfunc,eps=eps,ConstrType=ConstrType,
                                  opts=nloptr_args,...),
                   error=function(e){
                     if(MstepMessage){
                       message(paste("M-step optimization failed for item",j,"at iteration",itr))
                       message("Error message from nloptr::nloptr:")
                       message(e)
                     }
                     return("error")
                   })
    if(length(op)==1&&op=="error"){
      ret$convergence <- -1
    }else{
      ret$delta <- op$solution
      ret$opt <- op
      ret$phat <- c(Calc_Pj(op$solution,designMj,linkfunc,boundary = 0))
      if(any(ret$phat<0)||any(ret$phat>1)){
        message(paste("Estimates are out of bounds for item",j,"at iteration",itr,"[nloptr:nloptr]"),call. = FALSE)
        ret$convergence <- -3
      }else{
        ret$convergence <- 0
      }
    }
  }
  ret
}


initials_optim <- function(par,modelj,Nj,Rj,designMj,uPj,lPj, ConstrMatrix,linkfunc,eps=1e-5,ConstrType){
  # print(par)
  np <- length(par)
  pj <- c(Calc_Pj(par,designMj,linkfunc,boundary = 0))
  if(ConstrType==1){

    if(any(pj<=lPj)||any(pj>=uPj)){
      # print(pj)
      pj[which(pj<=lPj)] <- lPj + eps
      pj[which(pj>=uPj)] <- uPj - eps

      tmpar <- Calc_Dj(par = pj[seq_len(np)],
                       designMj = designMj[seq_len(np),],
                       linkfunc = linkfunc, boundary = 0)
      tmpj <- c(Calc_Pj(tmpar,designMj,linkfunc,boundary = 0))
      if(any(tmpj<=lPj)||any(tmpj>=uPj)){
        if(modelj==4){
          pj <- qlogis(pj)
        }else if(modelj==5){
          pj <- log(pj)
        }
        tmpar <- c(min(pj),rep((max(pj)-min(pj))/(np-1),np-1))
      }
      par <- tmpar
      # print(c(Calc_Pj(par,designMj,linkfunc,boundary = 0)))
    }
  }else if(ConstrType==3){
    if(modelj==3){#ACDM
      newlPj <- lPj
      newuPj <- uPj
      newpj <- pj
    }else if (modelj==4){
      newlPj <- qlogis(lPj)
      newuPj <- qlogis(uPj)
      newpj <- qlogis(pj)
    }else if (modelj==5){
      newlPj <- log(lPj)
      newuPj <- log(uPj)
      newpj <- log(pj)
    }
    if(par[1]<=newlPj || any(par[2:np]<=0)){
      par[which(par<=newlPj)] <- newlPj + eps
      par[which(par[2:np]<=0)] <- eps
    }
    if (sum(par)>=newuPj) {
      par <- par - par * (sum(par) - newuPj + eps) / sum(par)
    }
    if (par[1]<=newlPj || any(par[2:np]<=0) || sum(par)>=newuPj) {
      Lpar <- which.max(c(newlPj + eps,min(newpj,na.rm = TRUE)))
      Upar <- which.min(c(newuPj - eps,max(newpj,na.rm = TRUE)))
      par <- c(Lpar,rep((Upar-Lpar)/(np-1),np-1))
    }
  }


  return(par)
}

GNLOptim_call <- function(par,solver,modelj,correction,auglag_args,nloptr_args,solnp_args,item.prior=item.prior,
                          Nj,Rj,designMj,uPj,lPj,MstepMessage,itr,j,ConstrMatrix,linkfunc,eps,ConstrType){

  if(item.prior$on){ # adding prior to item parameters

    ret <- list()
    if(modelj==0){

      phat <- c((Rj+item.prior$beta[1]-1)/(Nj+sum(item.prior$beta)-2))

      phat[phat<lPj] <- lPj
      phat[phat>uPj] <- uPj
      ret$delta <- c(Calc_Dj(phat, designMj, linkfunc = linkfunc, boundary = 1))
      ret$opt <- list(convergence=0,message="converged - closed-form solution.")
      ret$phat <- phat
      ret$convergence <- 0
      return(ret)
    }else if(modelj==1 | modelj==2){
      rNj <- c(rep(sum(Nj[designMj[,2]==0]),sum(designMj[,2]==0)),
               rep(sum(Nj[designMj[,2]==1]),sum(designMj[,2]==1)))
      rRj <- c(rep(sum(Rj[designMj[,2]==0]),sum(designMj[,2]==0)),
               rep(sum(Rj[designMj[,2]==1]),sum(designMj[,2]==1)))

      phat <- c((rRj+item.prior$beta[1]-1)/(rNj+sum(item.prior$beta)-2))

      phat[phat<lPj] <- lPj
      phat[phat>uPj] <- uPj
      ret$delta <- c(Calc_Dj(phat, designMj, linkfunc = linkfunc, boundary = 1))
      ret$opt <- list(convergence=0,message="converged - closed-form solution.")
      ret$phat <- phat
      ret$convergence <- 0
      return(ret)

    }else if(linkfunc==2){ #logit link function

      op <- tryCatch(nloptr::slsqp(x0=par,fn=Mstep_obj_fn_prior,
                                    Nj=Nj,Rj=Rj,designMj=designMj,uPj=uPj,ConstrMatrix=ConstrMatrix,greaterthan0 = FALSE,
                                    lPj=lPj,linkfunc=linkfunc,eps=eps,ConstrType=ConstrType,m=item.prior$normal[1],sd=item.prior$normal[2]),
                     error=function(e){
                       if(MstepMessage){
                         message(paste("M-step optimization failed for item",j,"at iteration",itr))
                         message("Error message from nloptr::slsqp:")
                         message(e)
                       }
                       return("error")
                     })
      if(length(op)==1&&op=="error"){
        ret$convergence <- -1
      }else{
        ret$delta <- op$par
        ret$opt <- op
        ret$phat <- c(Calc_Pj(op$par,designMj,linkfunc,boundary = 0))
        if(any(ret$phat<0)||any(ret$phat>1)){
          message(paste("Estimates are out of bounds for item",j,"at iteration",itr,"[nloptr:slsqp]"),call. = FALSE)
          ret$convergence <- -3
        }else{
          ret$convergence <- 0
        }
      }
      return(ret)
    }else{
      stop("Priors cannot be imposed.",call. = FALSE)
    }
  }else{ # No priors

    if(solver=="auto"){
      if(ConstrType==1){
        if(modelj<3&modelj>=0){
          solver <- "ClosedForm"
        }else if(modelj>=3){
          solver <- c("BFGS","slsqp","auglag","solnp")
        }else if(modelj==-1){
          solver <- c("slsqp","auglag","solnp")
        }
      }else{
        if(modelj>0){
          solver <- c("BFGS","slsqp","auglag","solnp")
        }else if(modelj<=0){
          solver <- c("slsqp","auglag","solnp")
        }
      }
    }

    conv <- vector("numeric",length(solver))
    for(s in solver){
      # print(s)
      if(s=="slsqp"){nloptr_args$algorithm = "NLOPT_LD_SLSQP"}
      optims <- GNLOptim(par = par,fn=Mstep_obj_fn,gr=Mstep_obj_gr,hin=Mstep_ineq_fn,hin.gr=Mstep_ineq_jac,
                         optimizer=s,modelj=modelj,correction=correction,
                         auglag_args=auglag_args,nloptr_args=nloptr_args,solnp_args = solnp_args,
                         Nj=Nj,Rj=Rj,designMj=designMj,
                         uPj=uPj,lPj=lPj,MstepMessage=MstepMessage,itr=itr,j=j,
                         ConstrMatrix=ConstrMatrix,linkfunc=linkfunc,eps=1e-16,ConstrType=ConstrType)
      conv[s] <- optims$convergence
      if(optims$convergence==0) break

    }
    if(all(conv<0)){
      return(NULL)
    }else{
      return(optims)
    }
  }



}
