initials_optim <- function(ini,Mj,upper.pj, lower.pj,Kjj,model){
  vdelta <- ini
  bound.eps <- 0.0001
  if(model==3){#ACDM
    if (any(ineq_identity(vdelta=vdelta,Mj=Mj,upper.pj=upper.pj,
                          lower.pj=lower.pj)<=0) ) {
      vdelta[1] <- vdelta[1] + 2 * bound.eps
      vdelta[which.max(vdelta)] <- vdelta[which.max(vdelta)] - 3 * bound.eps
    }
    if (any(ineq_identity(vdelta=vdelta,Mj=Mj,upper.pj=upper.pj,
                          lower.pj=lower.pj)<=0) ) {
      vdelta <- c(0.1,rep(0.8/(length(vdelta)-1),(length(vdelta)-1)))
    }
  }else if(model==5){#RRUM

    if (any(ineq_rrum(vdelta,upper.pj=upper.pj,
                      lower.pj=lower.pj,Mj=Mj)<=0)){
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




Mstep_optim <- function(ini,Mj,upper.pj,lower.pj,Kjj,
                        model,Nj,Rj,optimizer,j,itr,
                        warning.immediate=FALSE,
                        warning.print = FALSE){


  ui <- rbind(Mj,-1*Mj)
  ci <- list(dina=NULL,
             dino=NULL,
             acdm=rep(c(lower.pj,-1*upper.pj),each=2^Kjj),
             llm=rep(c(qlogis(lower.pj),-1*qlogis(upper.pj)),each=2^Kjj),
             rrum=rep(c(log(lower.pj),-1*log(upper.pj)),each=2^Kjj))
# print(model)
  if(optimizer %in% c("Nelder-Mead","BFGS","CG")){
    ini <- initials_optim(ini=ini,Mj=Mj,upper.pj=upper.pj,
                          lower.pj=lower.pj,Kjj=Kjj,model=model)
  optims <- stats::constrOptim(ini,obj_fn,grad=gr_fn,method=optimizer,
                               ui=ui,ci=ci[[model]],
                               Nj=Nj,Rj=Rj,Mj=Mj,upper.pj=upper.pj,
                               lower.pj=lower.pj,Kjj=Kjj,model=model)
  if (warning.print&&optims$convergence>0)  {
    warning(paste("Unsuccessful M-step optimization for item",j,"at iteration",itr,"\n optimizer [constrOptim using",optimizer,"] message:",optims$message),
            call. = FALSE,immediate. = warning.immediate)
  }
}else if(tolower(optimizer)=="slsqp"){
    ineq_fn_acdm <- function(vdelta){
      if(model<=3){
        Pj <- Mj%*%vdelta
      }else if(model==4){
        Pj <- plogis(c(Mj%*%vdelta))
      }else if(model==5){
        Pj <- Mj%*%vdelta
        Pj <- exp(Pj)
      }
      c(Pj-lower.pj,upper.pj-Pj)
    }
    optims <- nloptr::slsqp(x0=ini,fn=obj_fn,gr=gr_fn,hin = ineq_fn_acdm,
                            Nj=Nj,Rj=Rj,Mj=Mj,upper.pj=upper.pj,
                            lower.pj=lower.pj,Kjj=Kjj,model=model)
    if (warning.print&&optims$convergence<0)  {
      warning(paste("Unsuccessful M-step optimization for item",j,"at iteration",itr,"\n optimizer [slsqp from nloptr] message:",optims$message),
              call. = FALSE,immediate. = warning.immediate)
    }
  }else if(tolower(optimizer)=="solnp"){
    optims <- Rsolnp::solnp(pars=ini,fun = obj_fn,ineqfun = ineq_fn,
                            ineqLB = rep(0,2*2^Kjj),ineqUB = rep(1,2*2^Kjj),
                            control = list(trace=0,tol=1e-6),
                            Nj=Nj,Rj=Rj,Mj=Mj,upper.pj=upper.pj,
                            lower.pj=lower.pj,Kjj=Kjj,model=model)
    if (warning.print&&optims$convergence>0)  {
      warning(paste("Unsuccessful M-step optimization for item",j,"at iteration",itr,"\n optimizer [solnp from Rsolnp] message: error code -",optims$convergence),
              call. = FALSE,immediate. = warning.immediate)
    }
  }else if(tolower(optimizer)=="auglag"){

    optims <- alabama::auglag(ini,obj_fn,gr=gr_fn,hin = ineq_fn,
                              Nj=Nj,Rj=Rj,Mj=Mj,upper.pj=upper.pj,
                              lower.pj=lower.pj,Kjj=Kjj,model=model,
                              control.outer = list(trace=FALSE))
    if (warning.print&&optims$convergence>0)  {
      warning(paste("Unsuccessful M-step optimization for item",j,"at iteration",itr,"\n optimizer [auglag from alabama] message: error code -",optims$convergence),
              call. = FALSE,immediate. = warning.immediate)
    }
  }else{
    stop(paste("Optimization method",optimizer,"is not available. Try BFGS, slsqp, solnp, auglag, Nelder-Mead, or CG."),call. = FALSE)
  }

    # print(optims)
    #check gradient
    #maxLik::compareDerivatives(obj_fn,grad=gr_fn,t0=ini,eps=1e-3,Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kjj,model=model)
  return(optims)
}

Mstep_optim_all <- function(ini,
                            Mj,
                            upper.pj,
                            lower.pj,
                            Kjj,
                            j,
                            itr,
                            model,
                            Nj,
                            Rj,
                            warning.immediate,
                            warning.print) {
  tryCatch({
    tryCatch({
      tryCatch({
        tryCatch({
          Mstep_optim(
            ini = c(ini),
            Mj = Mj,
            upper.pj = upper.pj,
            lower.pj = lower.pj,
            Kjj = Kjj,
            model = model,
            Nj = Nj,
            Rj = Rj,
            optimizer = "BFGS",
            warning.immediate = warning.immediate,
            warning.print = warning.print
          )
        }, error = function(e) {
          Mstep_optim(
            ini = c(ini),
            Mj = Mj,
            upper.pj = upper.pj,
            lower.pj = lower.pj,
            Kjj = Kjj,
            model = model,
            Nj = Nj,
            Rj = Rj,
            optimizer = "slsqp",
            warning.immediate = warning.immediate,
            warning.print = warning.print
          )
        })
      }, error = function(e) {
        Mstep_optim(
          ini = c(ini),
          Mj = Mj,
          upper.pj = upper.pj,
          lower.pj = lower.pj,
          Kjj = Kjj,
          model = model,
          Nj = Nj,
          Rj = Rj,
          optimizer = "solnp",
          warning.immediate = warning.immediate,
          warning.print = warning.print
        )
      })
    }, error = function(e) {
      Mstep_optim(
        ini = c(ini),
        Mj = Mj,
        upper.pj = upper.pj,
        lower.pj = lower.pj,
        Kjj = Kjj,
        model = model,
        Nj = Nj,
        Rj = Rj,
        optimizer = "auglag",
        warning.immediate = warning.immediate,
        warning.print = warning.print
      )
    })
  }, error = function(e) {
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
    return(list(par = ini))
  })

}

