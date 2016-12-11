
initials <- function(Q,nstarts=1,randomseed){

set.seed(randomseed)

if(nstarts==1){
  par <- gs2p(Q=Q,gs=matrix(runif(nrow(Q)*2,0.05,0.25),ncol = 2),
              model = rep(3,nrow(Q)),type = "equal",
              mono.constraint = TRUE)$itemprob.matrix
}else{
  par <- list(rep(NA,nstarts))
  for(i in 1:nstarts){
    par[[i]] <- gs2p(Q=Q,gs=matrix(runif(nrow(Q)*2,0.05,0.25),ncol = 2),
                model = rep(0,nrow(Q)),
                mono.constraint = rep(TRUE,nrow(Q)))$itemprob.matrix
  }

}
  return (par)
}


