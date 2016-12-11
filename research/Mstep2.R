
Mstep2 <- function(Kj,RN,model,itmpar,delta=NULL,constr=F,correction=0.0005,optimizer=NULL){
  # bound.eps <- 1e-4
  # J <- length(Kj)
  # options(warn = -1)
  # if (length(constr)==1){
  #   constr <- rep(constr,J)
  # }
  # opts <- vector("list",J)
  # #conv.code <- vector(length = J)
  # if (is.null(delta)) delta <- calc_delta(itmpar,model,Kj)
  # ## returned par is a 2^Kj[j] x 1 col vector
  # for(j in 1:J){ #for each item
  #   Mj <- designmatrix(Kj[j],model[j])
  #   Rj <- RN$Rg[j,1:2^Kj[j]]
  #   Nj <- RN$Ng[j,1:2^Kj[j]]
  #   ini <- delta[[j]]
  #   #print(j)
  #   if (constr[j]==FALSE){#no constraints are set
  #
  #     if (model[j]==0){#G-DINA
  #
  #       par <- (Rj+correction)/(Nj+2*correction)
  #       delta[[j]] <- c(solve(t(Mj)%*%Mj)%*%t(Mj)%*%par)
  #       optims <- list(convergence=0,message=NULL)
  #     }else if(model[j]==1 | model[j]==2){#DINA or DINO
  #
  #       Tm <- combn_matrix_din(Kj[j],model[j])
  #       par <- (Tm %*% Rj+correction)/(Tm %*% Nj+2*correction)
  #       delta[[j]] <- c(par[1],par[2^Kj[j]]-par[1])
  #       optims <- list(convergence=0,message=NULL)
  #     }else if (model[j]==3){#A-CDM
  #
  #       ini <- initials_optim(ini=ini,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j],model=model[j])
  #       optims <- constrOptim(ini,obj_fn_identity,grad=gr_identity,
  #                             ui=rbind(Mj,-1*Mj),ci=rep(c(bound.eps,bound.eps-1),each=2^Kj[j]),
  #                             Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #
  #       delta[[j]] <- optims$par
  #       par <- Mj%*%delta[[j]]
  #     }else if (model[j]==4){ #LLM
  #
  #       optims <- optim(ini,fn=obj_fn_llm,gr=gr_llm,method="BFGS",
  #                       Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #       delta[[j]] <- optims$par
  #       par <- Mj%*%delta[[j]]
  #       par <- (exp(par)/(1+exp(par)))
  #     }else if (model[j]==5){ #R-RUM
  #
  #       ini <- initials_optim(ini=ini,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j],model=model[j])
  #
  #       # optims <- constrOptim(ini,obj_fn_rrum,grad=gr_rrum,
  #       #                           ui=rbind(Mj,-1*Mj,c(-1,rep(0,length(ini)-1))),
  #       #                           ci=c(rep(c(log(bound.eps),log(1-bound.eps)),each=2^Kj[j]),0),
  #       #                           Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #       #print(ini)
  #       optims <- constrOptim(ini,obj_fn_rrum,grad=gr_rrum,
  #                             ui=rbind(Mj,-1*Mj,diag(1,length(ini),length(ini))),
  #                             ci=c(rep(c(log(bound.eps),log(1-bound.eps)),each=2^Kj[j]),rep(log(bound.eps),length(ini))),
  #                             Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #
  #       # ineq_rrum2 <- function(ini){
  #       #   Pj <- Mj%*%ini #Pj is lenth of 2^Kj[j]
  #       #   c(exp(Pj)-bound.eps,1-exp(Pj)-bound.eps)
  #       # }
  #       #
  #       # optims <- nloptr::auglag(ini,obj_fn_rrum,gr=NULL,hin = ineq_rrum2,localsolver = "SLSQP",
  #       #                Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #
  #       delta[[j]] <- optims$par
  #       par <- exp(Mj%*%delta[[j]])
  #     }
  #
  #   }else if (constr[j]==TRUE){  #constr==T
  #
  #     # Calculate Mj
  #     if (model[j]==0){#G-DINA
  #       ineq_const_gdina_mono <- function(vdelta){
  #         partialorder <- partial_order(Kj[j])
  #         Pj <- Mj%*%vdelta
  #         c(Pj-bound.eps,1-Pj-bound.eps,Pj[partialorder[,1]]-Pj[partialorder[,2]])
  #       }
  #       optims <- nloptr::slsqp(ini,obj_fn_identity,gr=gr_identity,hin = ineq_const_gdina_mono,
  #                               Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #       delta[[j]] <- optims$par
  #       par <- Mj%*%delta[[j]]
  #     }else if(model[j]<4){
  #
  #       ini <- initials_optim_mono(ini,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j],model=model[j])
  #       optims <- constrOptim(ini,obj_fn_identity,grad=gr_identity,ui=rbind(Mj,-1*Mj,diag(length(ini))),
  #                             ci=c(rep(c(bound.eps,bound.eps-1),each=2^Kj[j]),rep(0,length(ini))),
  #                             Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #       delta[[j]] <- optims$par
  #       par <- Mj%*%delta[[j]]
  #     }else if (model[j]==4){ #LLM
  #       ini <- initials_optim_mono(ini,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j],model=model[j])
  #       #print(ini)
  #       # optims <- optim(ini,fn=obj_fn_llm,gr=NULL,
  #       #                     lower=c(qlogis(bound.eps),rep(0,length(ini)-1)),
  #       #                     upper=rep(qlogis(1-bound.eps),length(ini)),method="L-BFGS-B",
  #       #                     Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #       optims <- solnp(ini,fun=obj_fn_llm,LB=c(qlogis(bound.eps),rep(0,length(ini)-1)),
  #                       UB=rep(qlogis(1-bound.eps),length(ini)),control = list(trace=0),
  #                       Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #
  #       delta[[j]] <- optims$par
  #       par <- Mj%*%delta[[j]]
  #       par <- (exp(par)/(1+exp(par)))
  #     }else if (model[j]==5){
  #       ini <- initials_optim_mono(ini,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j],model=model[j])
  #
  #       diagM <- diag(length(ini))
  #       diagM[1,1] <- (-1)*diagM[1,1]
  #       # print(ini)
  #       # print(round(Nj,5))
  #       # print(round(Rj,5))
  #       optims <- tryCatch({constrOptim(ini,obj_fn_rrum,grad=gr_rrum,method = "BFGS",
  #                                       ui=rbind(Mj,-1*Mj,diagM),
  #                                       ci=c(rep(c(log(bound.eps),log(1-bound.eps)),each=2^Kj[j]),rep(bound.eps,length(ini))),
  #                                       Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #       },warning=function(w){
  #         conditionMessage(w)
  #       },error=function(e){
  #         conditionMessage(e)
  #         ineq_rrum2 <- function(ini){
  #           ui=rbind(Mj,-1*Mj,diagM)
  #           ci=c(rep(c(log(bound.eps),log(1-bound.eps)),each=2^Kj[j]),rep(bound.eps,length(ini)))
  #           c(ui%*%ini)-ci
  #         }
  #         nloptr::slsqp(ini,obj_fn_rrum,gr=gr_rrum,hin = ineq_rrum2,Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #       })
  #       #print(optims)
  #
  #       # optims <- spg(par=ini, fn=obj_fn_rrum, gr=NULL, project="projectLinear",quiet = TRUE,control=list(trace=FALSE),
  #       #     projectArgs=list(A=rbind(Mj,-1*Mj,diagM),
  #       #                      b=c(rep(c(log(bound.eps),log(1-bound.eps)),each=2^Kj[j]),rep(bound.eps,length(ini))), meq=0),
  #       #     Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #       #
  #       # ineq_rrum2 <- function(ini){
  #       #   ui=rbind(Mj,-1*Mj,diagM)
  #       #   ci=c(rep(c(log(bound.eps),log(1-bound.eps)),each=2^Kj[j]),rep(bound.eps,length(ini)))
  #       #   c(ui%*%ini)-ci
  #       # }
  #       # optims <- nloptr::slsqp(ini,obj_fn_rrum,gr=gr_rrum,hin = ineq_rrum2,Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #       # optims <- nloptr::auglag(ini,obj_fn_rrum,gr=NULL,hin = ineq_rrum2,localsolver = "LBFGS",
  #       #                Nj=Nj,Rj=Rj,Mj=Mj,bound.eps=bound.eps,Kjj=Kj[j])
  #
  #       delta[[j]] <- optims$par
  #
  #       par <- exp(Mj%*%delta[[j]])
  #     }
  #   }
  #   par[par < bound.eps] <- bound.eps
  #   par[par > 1 - bound.eps] <- 1 - bound.eps
  #   itmpar[j,1:(2^Kj[j])] <- par
  #   opts[[j]] <- optims
  # }
  #
  # options(warn = 0)
  # return(list(item.param=itmpar,delta=delta,optims=opts))
}
