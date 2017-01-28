
Mstep <- function(Kj,RN,model,itmpar,delta,constr,
                  correction,lower.p,optimizer,
                  upper.p,itr,warning.immediate,designmatrices,optim.control=list()){

  J <- length(Kj)


  if (length(constr)==1){
    constr <- rep(constr,J)
  }
  opts <- vector("list",J)

  if (is.null(delta)) delta <- calc_delta(itmpar,model,Kj)

  for(j in 1:J){ #for each item
    Mj <- designmatrices[[j]]
    Rj <- RN$Rg[j,1:2^Kj[j]]
    Nj <- RN$Ng[j,1:2^Kj[j]]
    ini <- delta[[j]]
# print(Nj)

    if (constr[j]==FALSE){#no constraints are set

      if (model[j]==0){#G-DINA
        if (any(Nj<correction[2])){
          Nj[which(Nj<correction[2])] <- Nj[which(Nj<correction[2])] + correction[2]
          Rj[which(Nj<correction[2])] <- Rj[which(Nj<correction[2])] + correction[1]
        }
        par <- c(Rj/Nj)
        par[par<lower.p[j]] <- lower.p[j]
        par[par>upper.p[j]] <- upper.p[j]
        delta[[j]] <- c(solve(t(Mj)%*%Mj)%*%t(Mj)%*%par)
        optims <- list(convergence=0,message=NULL)
      }else if(model[j]==1 | model[j]==2){#DINA or DINO

        Nj <- c(combn_matrix_din(Kj[j],model[j])%*%Nj)
        Rj <- c(combn_matrix_din(Kj[j],model[j])%*%Rj)
        if (any(Nj<correction[2])){
          Nj[which(Nj<correction[2])] <- Nj[which(Nj<correction[2])] + correction[2]
          Rj[which(Nj<correction[2])] <- Rj[which(Nj<correction[2])] + correction[1]
        }
        par <- c(Rj/Nj)

        par[par<lower.p[j]] <- lower.p[j]
        par[par>upper.p[j]] <- upper.p[j]
        delta[[j]] <- c(par[1],par[2^Kj[j]]-par[1])
        optims <- list(convergence=0,message=NULL)
      }else if (model[j]>=3){
        if(tolower(optimizer)=="all"){
          optims <- Mstep_optim_all(ini=c(ini),Mj=Mj,lower.pj=lower.p[j],
                                    upper.pj=upper.p[j],j=j,itr=itr,
                               Kjj=Kj[j],model=model[j],Nj=Nj,Rj=Rj,
                               warning.immediate = warning.immediate,
                               warning.print=FALSE)
        }else{
          optims <- Mstep_optim(ini=c(ini),Mj=Mj,lower.pj=lower.p[j],
                                upper.pj=upper.p[j],j=j,itr=itr,
                                    Kjj=Kj[j],model=model[j],Nj=Nj,Rj=Rj,
                                optimizer=optimizer,
                                    warning.immediate = warning.immediate,
                                warning.print=TRUE)
        }

        delta[[j]] <- optims$par
        par <- Mj%*%delta[[j]]
      }


    }else if (constr[j]){

      if(tolower(optimizer)=="all"){
          optims <- Moptim_set_mono(ini=c(ini),Mj=Mj,lower.pj=lower.p[j],
                                    upper.pj=upper.p[j],j=j,itr=itr,
                                    Kjj=Kj[j],model=model[j],Nj=Nj,Rj=Rj,
                                    warning.immediate = warning.immediate,
                                    warning.print=FALSE)
        }else{
          if (model[j]==0){
            optims <- Mstep_optim_mono_gdina(ini=c(ini),Mj=Mj,lower.pj=lower.p[j],
                                upper.pj=upper.p[j],j=j,itr=itr,
                                Kjj=Kj[j],model=model[j],Nj=Nj,Rj=Rj,
                                optimizer=optimizer,
                                warning.immediate = warning.immediate,
                                warning.print=TRUE)
          }else{
            optims <- Mstep_optim_mono(ini=c(ini),Mj=Mj,lower.pj=lower.p[j],
                                             upper.pj=upper.p[j],j=j,itr=itr,
                                             Kjj=Kj[j],model=model[j],Nj=Nj,Rj=Rj,
                                             optimizer=optimizer,
                                             warning.immediate = warning.immediate,
                                             warning.print=TRUE,optim.control=optim.control)
          }
        }
      delta[[j]] <- optims$par
      par <- Mj%*%delta[[j]]
    }

    #print(delta[[j]])

    if (model[j]==4) par <- plogis(par) #logit(par)
    if (model[j]==5) par <- exp(par)
    itmpar[j,1:(2^Kj[j])] <- par
    opts[[j]] <- optims
  }

  return(list(item.parm=itmpar,delta=delta,optims=opts))
}
