
itmpar.gen <- function(Q,model="GDINA",quality="M",random=F){
  itm <- matrix(0,nrow(Q),2^(max(rowSums(Q)))) #initialize
  if (random==F){
    if (quality=="H"){
      g <- rep(0.1,nrow(Q))
    }else if(quality=="M"){
      g <- rep(0.2,nrow(Q))
    }else if(quality=="L"){
      g <- rep(0.3,nrow(Q))
    }else if (quality=="XH"){
      g <- rep(0.05,nrow(Q))
    }
  }else{
    if (quality=="H"){
      g <- runif(nrow(Q),0.05,0.15)
    }else if(quality=="M"){
      g <- runif(nrow(Q),0.15,0.25)
    }else if(quality=="L"){
      g <- runif(nrow(Q),0.25,0.35)
    }
  }

  if (model=="DINA"){
    for (j in 1:nrow(Q)){
      itm[j,2^sum(Q[j,])] <- 1 - g[j]
      itm[j,1:(2^sum(Q[j,])-1)] <- g[j]
    }

  }else if(model=="DINO"){
    for (j in 1:nrow(Q)){
      itm[j,1] <- g[j]
      itm[j,2:(2^sum(Q[j,]))] <- 1 - g[j]
    }
  }else if(model=="ACDM"){
    for (j in 1:nrow(Q)){
      delta <- (1-g[j]-g[j])/sum(Q[j,])
      itm[j,1:(2^sum(Q[j,]))] <- c(cbind(1,alpha(sum(Q[j,])))%*%c(g[j],rep(delta,sum(Q[j,]))))
    }
  }else if(model=="GDINA"){
    for (j in 1:nrow(Q)){
      if(sum(Q[j,])==1){
        itm[j,1] <- g[j]
        itm[j,2] <- 1-g[j]
      }else if(sum(Q[j,])==2){
        #p2 <- runif(1,g[j],1-g[j])
        #p3 <- runif(1,p2,1-g[j])
        itm[j,1:4] <- c(g[j],runif(2,g[j],1-g[j]),1-g[j])
      }else if(sum(Q[j,])==3){
        p1 <- runif(3,g[j],1-g[j])
        p2 <- c(runif(1,max(p1[1:2]),1-g[j]),runif(1,max(p1[c(1,3)]),1-g[j]),
                runif(1,max(p1[2:3]),1-g[j]))
        itm[j,1:8] <- c(g[j],p1,p2,1-g[j])
      }
    }
  }
  return (itm)
}


N <- c(1000,2000,4000)
quality <- c("L","M","H")
R <- 100
eps.v <- c(0.8,0.85,0.9,0.95,0.99)
Q <- sim30GDINA$simQ
for (n in N){
  for (q in quality){
    for (r in 1:R){

      itmpar <- itmpar.gen(Q,quality=q)
      sim <- GDINA.sim(n,Q=Q,item.param=itmpar)
      est <- GDINA(sim$dat,sim$Q,SE=FALSE,person.est = NULL,verbose = FALSE)
      Qval <- vector("list",length(eps.v))
      for (e in 1:length(eps.v)){
        Qval[[e]] <- Qvalidation(est,eps=eps.v[e])
      }

      #higher.order.par=list(theta = NULL, lambda = NULL)
    }
  }
}


est <- GDINA(frac20$dat,frac20$Q,higher.order = T,higher.order.model = "1PL")

sim <- GDINA.sim(1000,frac20$Q,item.param = est$itemprob.matrix)
est0 <- GDINA(sim$dat,sim$Q)
