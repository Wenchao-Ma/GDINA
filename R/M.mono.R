initials_optim_mono <- function(ini,Mj,upper.pj,lower.pj,Kjj,model){
  vdelta <- ini
  bound.eps <- 0.001
  if(model==3){#ACDM
    if (any(ineq_identity(vdelta=vdelta,Mj=Mj,upper.pj=upper.pj,
                          lower.pj=lower.pj)<=bound.eps)|any(vdelta<0) ) {
      vdelta[1] <- vdelta[1] + bound.eps
      vdelta[which.max(vdelta)] <- vdelta[which.max(vdelta)] - 2*bound.eps
      vdelta[which(vdelta<0)] <- bound.eps
    }
    if (any(ineq_identity(vdelta=vdelta,Mj=Mj,upper.pj=upper.pj,
                          lower.pj=lower.pj)<=bound.eps)|any(vdelta<0) ) {
      vdelta <- c(0.1,rep(0.8/(length(vdelta)-1),(length(vdelta)-1)))
    }
  }else if (model==4){
    vdelta[vdelta<qlogis(lower.pj)] <- qlogis(lower.pj) + bound.eps
    vdelta[vdelta>qlogis(upper.pj)] <- qlogis(upper.pj) - bound.eps
    if(vdelta[1]>0) vdelta[1] <- qlogis(lower.pj) + bound.eps
    vdelta[which(vdelta[2:length(vdelta)]<0)+1] <- bound.eps
  }else if(model==5){#RRUM

    if (any(ineq_rrum(vdelta,upper.pj=upper.pj,
                      lower.pj=lower.pj,Mj=Mj)<=bound.eps)|any(vdelta[2:length(vdelta)]<0)){
      #print("item ",j)
      vdelta <- vector("numeric",Kjj+1)
      vdelta[1] <- (-3.0)
      vdelta[Kjj+1] <- 2.8
      if (Kjj>1){
        vdelta[2:(Kjj+1)] <- 2.8/(Kjj)
      }
    }
  }
  return(vdelta)
}


Mstep_optim_mono_gdina <- function(ini,Mj,upper.pj,lower.pj,Kjj,
                                   model,Nj,Rj,optimizer,j,itr,
                                   warning.immediate,warning.print){

   # ini <- initials_optim_mono(ini=ini,Mj=Mj,bound.eps=bound.eps,Kjj=Kjj,model=model)

  ineq_const_gdina_mono <- function(vdelta){
    partialorder <- partial_order2(Kjj)
    Pj <- Mj%*%vdelta
    c(Pj-lower.pj,upper.pj-Pj,Pj[partialorder[,1]]-Pj[partialorder[,2]])
  }
  if(tolower(optimizer)=="slsqp"){
    optims <- nloptr::slsqp(x0=ini,fn=obj_fn,gr=gr_fn,hin = ineq_const_gdina_mono,
                            Nj=Nj,Rj=Rj,Mj=Mj,upper.pj=upper.pj,
                            lower.pj=lower.pj,Kjj=Kjj,model=model)
    if (warning.print&&optims$convergence<0)  {
      warning(paste("Unsuccessful M-step optimization for item",j,"at iteration",itr,"\n optimizer [slsqp from nloptr] message:",optims$message),
              call. = FALSE,immediate. = warning.immediate)
    }
    # cat("\n Item",j,"-slsqp\n")
  }else if(tolower(optimizer)=="solnp"){
    optims <- Rsolnp::solnp(pars=ini,fun = obj_fn,ineqfun = ineq_const_gdina_mono2,
                            ineqLB = rep(0,length(ineq_const_gdina_mono(ini))),
                            ineqUB = rep(1,length(ineq_const_gdina_mono(ini))),
                            control = list(trace=0),
                            Nj=Nj,Rj=Rj,Mj=Mj,upper.pj=upper.pj,
                            lower.pj=lower.pj,Kjj=Kjj,model=model)
    if (warning.print&&optims$convergence>0)  {
      warning(paste("Unsuccessful M-step optimization for item",j,"at iteration",itr,"\n optimizer [solnp from Rsolnp] message: error code -",optims$convergence),
              call. = FALSE,immediate. = warning.immediate)
    }

    # cat("\n Item",j,"-solnp\n")
  }else if(tolower(optimizer)=="auglag"){

    optims <- alabama::auglag(ini,obj_fn,gr=gr_fn,hin = ineq_const_gdina_mono2,
                              Nj=Nj,Rj=Rj,Mj=Mj,upper.pj=upper.pj,
                              lower.pj=lower.pj,Kjj=Kjj,model=model,
                              control.outer = list(trace=FALSE))
    if (warning.print&&optims$convergence>0)  {
      warning(paste("Unsuccessful M-step optimization for item",j,"at iteration",itr,"\n optimizer [auglag from alabama] message: error code -",optims$convergence),
              call. = FALSE,immediate. = warning.immediate)
    }
    # cat("\n Item",j,"-auglag\n")
  }else{
    stop(paste("Optimization method",optimizer,"is not available. Try slsqp, solnp, or auglag."),call. = FALSE)

  }

  return(optims)
}



# This function is for M-step for all reduced models
Mstep_optim_mono <- function(ini,Mj,upper.pj,
                             lower.pj,Kjj,model,
                             Nj,Rj,optimizer,j,itr,
                             warning.immediate,
                             warning.print,optim.control=list()){

  ini <- initials_optim_mono(ini=ini,Mj=Mj,upper.pj=upper.pj,
                             lower.pj=lower.pj,Kjj=Kjj,model=model)

  ineq_fn_mono_1 <- function(ini){
    if(model<=3){
      bd <- upper.pj
    }else if(model==4){
      bd <- qlogis(upper.pj)
    }else if(model==5){
      bd <- log(upper.pj)
    }
    bd-sum(ini)
  }

  LB <- rbind(c(lower.pj,lower.pj,lower.pj,qlogis(lower.pj),log(lower.pj)),
              matrix(0,length(ini)-1,5))

  ui.bfgs=rbind(Mj,-1*Mj,diag(length(ini)))
  ci.bfgs <- list(dina=c(rep(c(lower.pj,-1*upper.pj),each=2^Kjj),c(lower.pj,rep(0,length(ini)-1))),
                  dino=c(rep(c(lower.pj,-1*upper.pj),each=2^Kjj),c(lower.pj,rep(0,length(ini)-1))),
                  acdm=c(rep(c(lower.pj,-1*upper.pj),each=2^Kjj),c(lower.pj,rep(0,length(ini)-1))),
                  llm=c(rep(c(qlogis(lower.pj),-1*qlogis(upper.pj)),each=2^Kjj),
                        c(qlogis(lower.pj),rep(0,length(ini)-1))),
                  rrum=c(rep(c(log(lower.pj),-1*log(upper.pj)),each=2^Kjj),
                         c(log(lower.pj),rep(0,length(ini)-1))))
  if (optimizer %in% c("Nelder-Mead","BFGS","CG")){

    optims <- stats::constrOptim(ini,obj_fn,grad=gr_fn,method=optimizer,
                                 ui=ui.bfgs,ci=ci.bfgs[[model]],
                                 Nj=Nj,Rj=Rj,Mj=Mj,upper.pj=upper.pj,
                                 lower.pj=lower.pj,Kjj=Kjj,model=model,
                                 control = optim.control)
    if (warning.print&&optims$convergence>0)  {
      warning(paste("Unsuccessful M-step optimization for item",j,"at iteration",itr,"\n optimizer [constrOptim using",optimizer,"] message:",optims$message),
              call. = FALSE,immediate. = warning.immediate)
    }
    #check gradient
    #maxLik::compareDerivatives(obj_fn,grad=gr_fn,t0=ini,eps=1e-3,Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kjj,model=model)

  }else if(tolower(optimizer)=="slsqp"){

    optims <- nloptr::slsqp(x0=ini,fn=obj_fn,gr=gr_fn,hin = ineq_fn_mono_1,lower=c(LB[,model]),
                            Nj=Nj,Rj=Rj,Mj=Mj,upper.pj=upper.pj,
                            lower.pj=lower.pj,Kjj=Kjj,model=model,
                            control = optim.control)
    if (warning.print&&optims$convergence<0)  {
      warning(paste("Unsuccessful M-step optimization for item",j,"at iteration",itr,"\n optimizer [slsqp from nloptr] message:",optims$message),
              call. = FALSE,immediate. = warning.immediate)
    }
  }else if(tolower(optimizer)=="solnp"){
    optims <- Rsolnp::solnp(pars=ini,fun = obj_fn,ineqfun = ineq_fn_mono_2,
                            ineqLB = 0,LB=c(LB[,model]),ineqUB = +Inf,
                            control = list(trace=0,tol=1e-6),
                            Nj=Nj,Rj=Rj,Mj=Mj,upper.pj=upper.pj,
                            lower.pj=lower.pj,Kjj=Kjj,model=model,
                            control = optim.control)
    if (warning.print&&optims$convergence>0)  {
      warning(paste("Unsuccessful M-step optimization for item",j,"at iteration",itr,"\n optimizer [solnp from Rsolnp] message: error code -",optims$convergence),
              call. = FALSE,immediate. = warning.immediate)
    }
  }else if(tolower(optimizer)=="auglag"){
    if(is.null(optim.control$trace)) optim.control$trace <- FALSE
    optims <- alabama::auglag(ini,obj_fn,gr=gr_fn,hin = ineq_fn_mono,
                              Nj=Nj,Rj=Rj,Mj=Mj,upper.pj=upper.pj,
                              lower.pj=lower.pj,Kjj=Kjj,model=model,
                              control.outer = optim.control)
    if (warning.print&&optims$convergence>0)  {
      warning(paste("Unsuccessful M-step optimization for item",j,"at iteration",itr,"\n optimizer [auglag from alabama] message: error code -",optims$convergence),
              call. = FALSE,immediate. = warning.immediate)
    }
  }else{
    stop(paste("Optimization method",optimizer,"is not available. Try BFGS, slsqp, solnp, auglag, Nelder-Mead, or CG."),call. = FALSE)
  }
  #if (tolower(optimizer)!="bfgs") print(optimizer)
  return(optims)
}


Moptim_set_mono <- function(ini,Mj,upper.pj,
                            lower.pj,Kjj,model,
                            Nj,Rj,optim.method,j,itr,
                            warning.immediate,
                            warning.print){

  if(model==0){
    tryCatch({
      tryCatch({
        tryCatch({Mstep_optim_mono_gdina(ini=c(ini),Mj=Mj,upper.pj=upper.pj,
                                         lower.pj=lower.pj,Kjj=Kjj,model=model,
                                         Nj=Nj,Rj=Rj,optimizer="slsqp",
                                         warning.immediate=warning.immediate,
                                         warning.print=warning.print)
        },error=function(e){
          Mstep_optim_mono_gdina(ini=c(ini),Mj=Mj,upper.pj=upper.pj,
                                 lower.pj=lower.pj,Kjj=Kjj,model=model,
                                 Nj=Nj,Rj=Rj,optimizer="solnp",
                                 warning.immediate=warning.immediate,
                                 warning.print=warning.print)
        })
      },error=function(e){
        Mstep_optim_mono_gdina(ini=c(ini),Mj=Mj,upper.pj=upper.pj,
                               lower.pj=lower.pj,Kjj=Kjj,model=model,Nj=Nj,
                               Rj=Rj,optimizer="auglag",
                               warning.immediate=warning.immediate,
                               warning.print=warning.print)
      })},error=function(e){
        warning(
          paste(
            "Unsuccessful M-step optimization for item",
            j,
            "at iteration",
            itr,
            "\n Check your solution carefully!\n"
          ),
          call. = FALSE,
          immediate. = warning.immediate
        )
      return(list(par=ini))
    })
    }else{

   tryCatch({tryCatch({
    tryCatch({
      tryCatch({Mstep_optim_mono(ini=c(ini),Mj=Mj,upper.pj=upper.pj,
                                 lower.pj=lower.pj,Kjj=Kjj,model=model,
                                 Nj=Nj,Rj=Rj,optimizer="BFGS",
                                 warning.immediate=warning.immediate,
                                 warning.print=warning.print)
      },error=function(e){
        Mstep_optim_mono(ini=c(ini),Mj=Mj,upper.pj=upper.pj,
                         lower.pj=lower.pj,Kjj=Kjj,model=model,Nj=Nj,
                         Rj=Rj,optimizer="slsqp",
                         warning.immediate=warning.immediate,
                         warning.print=warning.print)
      })
    },error=function(e){
      Mstep_optim_mono(ini=c(ini),Mj=Mj,upper.pj=upper.pj,
                       lower.pj=lower.pj,Kjj=Kjj,model=model,Nj=Nj,
                       Rj=Rj,optimizer="solnp",
                       warning.immediate=warning.immediate,
                       warning.print=warning.print)
    })
  },error=function(e){
    Mstep_optim_mono(ini=c(ini),Mj=Mj,upper.pj=upper.pj,
                     lower.pj=lower.pj,Kjj=Kjj,model=model,
                     Nj=Nj,Rj=Rj,optimizer="auglag",
                     warning.immediate=warning.immediate,
                     warning.print=warning.print)
  })},error=function(e){
    warning(
      paste(
        "Unsuccessful M-step optimization for item",
        j,
        "at iteration",
        itr,
        "\n Check your solution carefully!\n"
      ),
      call. = FALSE,
      immediate. = warning.immediate
    )
    return(list(par=ini))
  })
}
}
