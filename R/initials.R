
initials <- function(Q,nstarts=1,randomseed,latent.var="att"){
# set.seed(randomseed)
if(tolower(latent.var)=="att"){
  if(nstarts==1){
    par <- gs2p(Q=Q,gs=matrix(runif(nrow(Q)*2,0.05,0.25),ncol = 2),
                model = rep(3,nrow(Q)),type = "equal",
                mono.constraint = TRUE)$itemprob.matrix
  }else if(nstarts==3){
    par <- list(rep(NA,nstarts))
    for(i in seq_len(nstarts)){
      par[[i]] <- gs2p(Q=Q,gs=matrix(runif(nrow(Q)*2,0.05,0.25),ncol = 2),
                       model = rep(i,nrow(Q)),type = "equal",
                       mono.constraint = rep(TRUE,nrow(Q)))$itemprob.matrix
    }
  }else{
    par <- list(rep(NA,nstarts))
    for(i in 1:nstarts){
      par[[i]] <- gs2p(Q=Q,gs=matrix(runif(nrow(Q)*2,0.05,0.25),ncol = 2),
                       model = rep(0,nrow(Q)), type = "random",
                       mono.constraint = rep(TRUE,nrow(Q)))$itemprob.matrix
    }

  }
}else if(tolower(latent.var)=="bugs"){
  if(nstarts==1){
    ip <- gs2p(Q=Q,gs=matrix(runif(nrow(Q)*2,0.05,0.25),ncol = 2),
                model = rep(3,nrow(Q)),type = "equal",
                mono.constraint = TRUE)$itemprob.parm
    ip <- lapply(ip,function(x)1-x)
    if(is.list(ip)){
      par <- l2m(ip)
    }else if(is.matrix){
      par <- t(ip)
    }

  }else{
    par <- list(rep(NA,nstarts))
    for(i in 1:nstarts){
      ip <- gs2p(Q=Q,gs=matrix(runif(nrow(Q)*2,0.05,0.25),ncol = 2),
                 model = rep(3,nrow(Q)),type = "random",
                 mono.constraint = TRUE)$itemprob.parm
      ip <- lapply(ip,function(x)1-x)
      if(is.list(ip)){
        par[[i]] <- l2m(ip)
      }else if(is.matrix){
        par[[i]] <- t(ip)
      }

    }

  }
  }
  return (par)
}


