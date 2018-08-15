Louis <- function(GDINA.obj,model = "GDINA",SE.type = 2){
  Q <- extract(GDINA.obj,"Q")
  Eta <- LC2LG(as.matrix(Q))
  if(model=="DINA"){
    Eta <- 1*(Eta==Eta[,ncol(Eta)])+1
  }else if(model=="DINO"){
    Eta <- 1*(Eta!=1)+1
  }

  mv <- matchMatrix(unique(extract(GDINA.obj,"seq.dat")),extract(GDINA.obj,"seq.dat"))
  FirstUnique <- which(!duplicated(mv))

  dat <- as.matrix(extract(GDINA.obj,"seq.dat")[FirstUnique,])
  post <- GDINA.obj$technicals$logposterior.i[FirstUnique,]
  weight <- table(mv)

  np <- apply(Eta,1,max)

  L <- LouisC(mX=as.matrix(dat),
              np=matrix(np,ncol = 1),
              mlogPost = as.matrix(post),
              itemparmLC = as.matrix(GDINA.obj$pf),
              parloc = as.matrix(Eta),
              weight = matrix(weight,ncol=1),SEtype=SE.type)

  if(all(eigen(L$An,only.values = TRUE)$values>1e-8)){
    nonPD <- FALSE
    CovA=L$invAn
    CovRobust <- L$robust
  }else{
    nonPD <- TRUE
    CovA <- solve(as.matrix(Matrix::nearPD(as.matrix(L$An))$mat))
    CovRobust <- CovA%*%L$term3%*%CovA
  }
  swSE <- sqrt(diag(CovRobust))
  AnSE <- sqrt(diag(CovA))
  covB <- solve(L$term3)
  BnSE <- sqrt(diag(covB))
  ind2 <- cumsum(np)
  ind1 <- cumsum(c(1,np))
  SE.A <- SE.B <- SE.robust <- vector("list",length(np))
  for(j in 1:length(np)){
    SE.A[[j]] <- AnSE[ind1[j]:ind2[j]]
    SE.B[[j]] <- BnSE[ind1[j]:ind2[j]]
    SE.robust[[j]] <- swSE[ind1[j]:ind2[j]]
  }
  return(list(SE.A=SE.A,SE.B=SE.B,SE.robust=SE.robust,nonPD = nonPD,
              Cov.robust=CovRobust,Cov.B=covB,Cov.A=CovA,
              SEall=cbind(swSE,AnSE,BnSE),ind=cbind(ind1[-length(ind1)],ind2)))
}
