partial_order <- function(Kjj){
  alp <- alpha(Kjj)
  out <- NULL
  for (i in 1:nrow(alp)){
    out <- rbind(out,cbind(which(apply(alp,1,function(x) return(all(x-alp[i,]>=0)))),i))
  }
  out <- out[out[,1]>out[,2],,drop=F]
  colnames(out) <- c("l","s")
  return(out)
}

partial_order2 <- function(Kjj){
  alp <- alpha(Kjj)
  alp <- cbind(c(1:nrow(alp)),rowSums(alp),alp)
  out <- NULL

  for(k in 1:max(alp[,2])){
    for(i in 1:sum(alp[,2]==k-1)){
      alpk_1 <- alp[alp[,2]==k-1,,drop=FALSE]
      alpk <-  alp[alp[,2]==k,,drop=FALSE]
      out <- rbind(out,cbind(alpk[(apply(alpk,1,function(x){all(x-alpk_1[i,]>=0)})),1],alpk_1[i,1]))
      }
  }
  colnames(out) <- c("l","s")
  return(out)
}


combn_matrix_din <- function(Kjj,model){
  if (model==1){
    Tm <- matrix(1,2^Kjj,2^Kjj)
    Tm[,2^Kjj] <- 0
    Tm[2^Kjj,] <- 0
    Tm[2^Kjj,2^Kjj] <- 1
  }else if (model==2){
    Tm <- matrix(1,2^Kjj,2^Kjj)
    Tm[1,2:2^Kjj] <- 0
    Tm[2:2^Kjj,1] <- 0
  }
  return(Tm)
}

obj_fn <- function(vdelta,Nj,Rj,Mj,upper.pj,lower.pj,Kjj,model){ #objective function
  if(model<=3){
    Pj <- Mj%*%vdelta
    # print(Pj)
  }else if(model==4){

    Pj <- Mj%*%vdelta
    Pj <- exp(Pj)/(1+exp(Pj))
    #print(Pj)
  }else if(model==5){
    Pj <- Mj%*%vdelta
    Pj <- exp(Pj)
  }
  -1*sum(log(Pj)*Rj+log(1-Pj)*(Nj-Rj))
}

ineq_fn <- function(vdelta,Nj,Rj,Mj,upper.pj,lower.pj,Kjj,model){
  if(model<=3){
    Pj <- Mj%*%vdelta
  }else if(model==4){
    Pj <- plogis(c(Mj%*%vdelta))
  }else if(model==5){
    Pj <- Mj%*%vdelta
    Pj <- exp(Pj)
  }
  # print(c(Pj-lower.pj,upper.pj-Pj))
  c(Pj-lower.pj,upper.pj-Pj)
}

gr_fn <- function(vdelta,Nj,Rj,Mj,upper.pj,lower.pj,Kjj,model){ #objective function
  if(model<=3){
    Pj <- c(Mj%*%vdelta)
    gr <- -1*colSums(Mj*(Rj/Pj-(Nj-Rj)/(1-Pj)))
  }else if(model==4){
    Pj <- plogis(c(Mj%*%vdelta))
    gr <- -1*colSums(Mj*(Rj-Nj*Pj))
  }else if(model==5){
    Pj <- c(Mj%*%vdelta)
    Pj <- exp(Pj)
    gr <- -1*colSums(Mj*((Rj/Pj-(Nj-Rj)/(1-Pj))*Pj))
  }
  return(gr)
}


ineq_fn_mono <- function(vdelta,Nj,Rj,Mj,upper.pj,lower.pj,Kjj,model){
  if(model<=3){
    Pj <- Mj%*%vdelta
  }else if(model==4){
    Pj <- plogis(c(Mj%*%vdelta))
  }else if(model==5){
    Pj <- Mj%*%vdelta
    Pj <- exp(Pj)
  }
  c(Pj-lower.pj,upper.pj-Pj,vdelta[-1])
}


obj_fn_identity <- function(vdelta,Nj,Rj,Mj,upper.pj,lower.pj,Kjj){ #objective function
  Pj <- Mj%*%vdelta #Pj is lenth of 2^Kj[j]
  -1*sum(log(Pj)*Rj+log(1-Pj)*(Nj-Rj))
}


ineq_identity <- function(vdelta,Nj,Rj,Mj,upper.pj,lower.pj,Kjj){
  c(Mj%*%vdelta-lower.pj,upper.pj-Mj%*%vdelta)
}

gr_identity <- function(vdelta,Nj,Rj,Mj,upper.pj,lower.pj,Kjj){ #gradient
  Pj <- Mj%*%vdelta
  -1*colSums(Mj*matrix(rep(Rj/Pj-(Nj-Rj)/(1-Pj),ncol(Mj)),ncol = ncol(Mj)))
}

obj_fn_llm <- function(vdelta,Nj,Rj,Mj,upper.pj,lower.pj,Kjj){ #objective function
  #vdelta is a vector Kj+1 elements for LLM
  Pj <- Mj%*%vdelta #Pj is lenth of 2^Kj[j]
  Pj <- exp(Pj)/(1+exp(Pj))
  -1*sum(log(Pj)*Rj+log(1-Pj)*(Nj-Rj))
}

gr_llm <- function(vdelta,Nj,Rj,Mj,upper.pj,lower.pj,Kjj){ #gradient
  Pj <- Mj%*%vdelta #Pj is lenth of 2^Kj[j]
  Pj <- exp(Pj)/(1+exp(Pj))
  -1*colSums(Mj*matrix(rep(Rj-Nj*Pj,Kjj+1),ncol = Kjj+1))
}

obj_fn_rrum <- function(vdelta,Nj,Rj,Mj,upper.pj,lower.pj,Kjj){ #objective function
  #vdelta is a vector Kj+1 elements for LLM
  Pj <- Mj%*%vdelta #Pj is lenth of 2^Kj[j]
  Pj <- exp(Pj)
  Pj[Pj<lower.pj] <- lower.pj
  Pj[Pj>upper.pj] <- upper.pj
  -1*sum(log(Pj)*Rj+log(1-Pj)*(Nj-Rj))

}


gr_rrum <- function(vdelta,Nj,Rj,Mj,upper.pj,lower.pj,Kjj){
  Pj <- Mj%*%vdelta #Pj is lenth of 2^Kj[j]
  Pj <- exp(Pj)
  -1*colSums(Mj*matrix(rep((Rj/Pj-(Nj-Rj)/(1-Pj))*Pj,Kjj+1),ncol = Kjj+1))

}
ineq_rrum <- function(ini,upper.pj,lower.pj,Mj){
  Pj <- Mj%*%ini #Pj is lenth of 2^Kj[j]
  c(exp(Pj)-lower.pj,upper.pj-exp(Pj))
}


ineq_const_gdina_mono2 <- function(vdelta,Nj,Rj,Mj,upper.pj,lower.pj,Kjj,model){
  partialorder <- partial_order2(Kjj)
  Pj <- Mj%*%vdelta
  c(Pj-lower.pj,upper.pj-Pj,Pj[partialorder[,1]]-Pj[partialorder[,2]])
}

ineq_fn_mono_2 <- function(ini,Nj,Rj,Mj,upper.pj,lower.pj,Kjj,model){
  if(model<=3){
    bd <- upper.pj
  }else if(model==4){
    bd <- qlogis(upper.pj)
  }else if(model==5){
    bd <- log(upper.pj)
  }
  bd-sum(ini)
}
