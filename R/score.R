# score function for probabilities
score_pj <- function(Xj,                   # a vector of item responses to item j
                    parloc.j,         # parameter locations for item j - H by 2^K matrix
                    catprob.j,        # a list with H elements giving the reduced catprob.parm for each nonzero category
                    logpost){          # log posterior N x 2^K

  Xj[is.na(Xj)] <- -1 # remove NA from the data
  Lj <- sapply(catprob.j,length)
  Kj <- log2(Lj)
  post <- exp(logpost) # posterior N x 2^K
  score.p <- vector("list",nrow(parloc.j)) #I(Xij>=x)/sjx - I(Xij=x-1)/[1-sjx]
  for(r in 1:length(score.p)){
    score.p[[r]] <- aggregateCol(post,parloc.j[r,])*(outer((Xj>=r),catprob.j[[r]],"/")-outer((Xj==r-1),(1-catprob.j[[r]]),"/"))
  }

  return(score.p)
}


# input is the estimation object
scorefunc2 <- function(object, item = NULL, ...){
  dat <- extract(object,"dat")

  Qc <- extract(object,"Qc")
  Q <- extract(object,"Q")

  J <- ncol(dat)
  par.loc <- eta(as.matrix(Q))
  post <- exp(indlogPost(object)) # posterior N x 2^K

  # urP <- itemparm(object,"itemprob",digits = 9) #unconditional reduced prob.
  C <- table(Qc[,1]) # number of categories
  S0 <- rep(1:J,C+1) # loc for each item including cat 0

  LC.Prob <- sequP(as.matrix(par.loc),as.matrix(object$catprob.matrix),C)
  c.LC.Prob <- LC.Prob$cPr # no cat 0
  u.LC.Prob <- LC.Prob$uPr # including cat 0
  scorej <- vector("list",J)
  if(is.null(item)) item <- 1:J
  for(j in item){
    Xj <- dat[,j]
    sj <- c.LC.Prob[which(Qc[,1]==j),,drop=FALSE] # processing function for item j
    # print(sj)
    locj <- par.loc[which(Qc[,1]==j),,drop=FALSE] # loc for item j
    Pj <- u.LC.Prob[which(S0==j),][-1,,drop=FALSE] # drop cat 0 - P(Xij=h|alpha)
    Pj0 <- u.LC.Prob[which(S0==j),][1,] # P(Xij=0|alpha)
    Dj <- list() #I(Xij=h)/P(Xij=h|alpha) - I(Xij=0)/P(Xij=0|alpha)
    for(h in seq_len(C[j])) Dj[[h]] <- outer((Xj==h),Pj[h,],"/")-outer((Xj==0),Pj0,"/")
    # print(head(Dj[[h]]))
    for(k in seq_len(C[j])) {
      scorejk <- 0
      for(h in seq_len(C[j])){
        if(k<=h) {
          tmp <- Dj[[h]]*matrix(Pj[h,]/sj[k,],nrow = nrow(Dj[[h]]),ncol = ncol(Dj[[h]]),byrow = TRUE)
        }else if(k==h+1){
          tmp <- Dj[[h]]*matrix((-1)*(Pj[h,]/(1-sj[k,])),nrow = nrow(Dj[[h]]),ncol = ncol(Dj[[h]]),byrow = TRUE)
        }else{
          tmp <- matrix(0,nrow = nrow(Dj[[h]]),ncol = ncol(Dj[[h]]))
        }
        scorejk <- scorejk + tmp
      }
      tmp <- post*scorejk
      # print(tmp[1:5,1:5])
      for(lj in sort(unique(locj[k,]))) scorej[[j]] <- cbind(scorej[[j]],rowSums(tmp[,which(locj[k,]==lj),drop=FALSE]))

    }

  }
  return(scorej)
}

# only suitable for DINA, DINO and GDINA
score_p <- function(object){
  model <- extract(object,"models_numeric")
  if(any(model>2)) stop("Score function for success probabilities only available for the G-DINA, DINA and DINO models.",call. = FALSE)
  Q <- extract(object,"Q")
  Qc <- extract(object,"Qc")
  NC <- nrow(Q)
  N <- extract(object,"nobs")
  J <- extract(object,"nitem")
  Kj <- rowSums(Q>0)
  if(extract(object,"sequential")){
    dat <- seq_coding(extract(object,"dat"),Qc)
  }else{
    dat <- extract(object,"dat")
  }

  pj <- as.matrix(extract(object,"catprob.matrix"))
  Lj <- 2^rowSums(Q>0)
  for(j in which(model %in% c(1,2))){#DINA or DINO
    pj[j,2] <- pj[j,Lj[j]]
    if (Lj[j]>2) pj[j,3:ncol(pj)] <- -1
  }

  scof <- scorefun(mX=as.matrix(dat),
                   mlogPost=as.matrix(extract(object,"logposterior.i")),
                   itmpar=pj,
                   parloc=eta(as.matrix(Q)),
                   model=model)
  scorep <- scof$score
  ind <- scof$index + 1 # col 1: location; col 2: category
  score <- vector("list",NC)
  for (j in 1:NC){
    score[[j]] <- as.matrix(scorep[,ind[which(ind[,2]==j),1]]) #score function for P(\alpha_c) for cateogry j
  }
  names(score) <- extract(object,"item.names")
  lik <- exp(indlogLik(object))

  if(extract(object,"ngroup")>1){

    for(g in sort(unique(extract(object,"group")))){
      tmp <- ((lik-lik[,1])/colSums(c(extract(object,"posterior.prob")[g,])*t(lik)))[,-1]
      l <- length(score)+1
      print(l)
      score[[l]] <- tmp
      score[[l]][which(extract(object,"group")!=g),] <- 0
      names(score)[l] <- paste0("G",g)
    }
  }else{
    score[[length(score)+1]] <- ((lik-lik[,1])/colSums(c(extract(object,"posterior.prob"))*t(lik)))[,-1]
  }
  return(score)
}


score_d <- function(object){
  Q <- extract(object,"Q")
  NC <- nrow(Q)
  N <- extract(object,"nobs")
  Kj <- rowSums(Q>0)
  pj <- extract(object,"catprob.parm")
  linkfunc <- extract(object,"linkfunc")
  des <- extract(object,"designmatrix")
  if(extract(object,"sequential")){
    Qc <- extract(object,"Qc")
    dat <- seq_coding(extract(object,"dat"),Qc)
  }else{
    dat <- extract(object,"dat")
  }
  model <- extract(object,"models_numeric")
  scof <- scorefun(mX=as.matrix(dat),
                   mlogPost=as.matrix(extract(object,"logposterior.i")),
                   itmpar=as.matrix(extract(object,"catprob.matrix")),
                   parloc=eta(as.matrix(Q)),
                   model=rep(0,NC))

  scorep <- scof$score
  ind <- scof$index + 1
  score <- vector("list",NC)
  for (r in 1:NC){
    scorepj <- as.matrix(scorep[,ind[which(ind[,2]==r),1]]) #score function for P(\alpha_c) for cateogry j

    if(linkfunc[r]=="identity"){
      score[[r]] <- scorepj%*%des[[r]]
    }else if(linkfunc[r]=="logit"){
      score[[r]] <- (scorepj*matrix(pj[[r]]*(1-pj[[r]]),nrow = N,ncol = length(pj[[r]]),byrow = TRUE))%*%des[[r]]
    }else if(linkfunc[r]=="log"){
      score[[r]] <- (scorepj*matrix(pj[[r]],nrow = N,ncol = length(pj[[r]]),byrow = TRUE))%*%des[[r]]
    }

  }
  names(score) <- extract(object,"item.names")
  lik <- exp(indlogLik(object))

  if(extract(object,"ngroup")>1){

    for(g in sort(unique(extract(object,"group")))){
      tmp <- ((lik-lik[,1])/colSums(c(extract(object,"posterior.prob")[g,])*t(lik)))[,-1]
      l <- length(score)+1
      # print(l)
      score[[l]] <- tmp
      score[[l]][which(extract(object,"group")!=g),] <- 0
      names(score)[l] <- paste0("G",g)
    }
  }else{
    score[[length(score)+1]] <- ((lik-lik[,1])/colSums(c(extract(object,"posterior.prob"))*t(lik)))[,-1]
  }

  return(score)
}

# variance and SE of delta parameters for all models
OPG_d <- function(object,SE.type){
  Q <- extract(object,"Q")
  Qc <- extract(object,"Qc")
  J <- extract(object,"nitem")
  NC <- nrow(Q)
  NG <- extract(object,"ngroup")
  scorejh <- score_d(object) # a list of score function for delta with elements for each category
  IP.loc <- length(scorejh) - NG

  np <- sapply(scorejh,ncol)[seq_len(IP.loc)] # the last element is for mixing proportions
# print(np)
  if(SE.type == 1){
    scorej <- vector("list",J)
    scorejh <- scorejh[seq_len(IP.loc)]
    for(j in 1:J)  scorej[[j]] <- do.call(cbind,scorejh[which(Qc[,1]==j)])
    vars <- bdiag(lapply(scorej,inverse_crossprod))
  }else if(SE.type == 2){
    scorejh <- scorejh[seq_len(IP.loc)]
    vars <- inverse_crossprod(do.call(cbind,scorejh))
  }else if(SE.type==3){
    vars <- inverse_crossprod(do.call(cbind,scorejh))
    vars <- vars[1:sum(np),1:sum(np)]
  }
  covIndex <- data.frame(item = rep(Qc[,1],np),
                         itemcat = rep(Qc[,2],np),
                         cat = rep(1:NC,np),
                         par = unlist(lapply(np,seq_len)),
                         loc = 1:sum(np),row.names = NULL)
  se.vector <- sqrt(diag(vars))
  se <- vector("list",NC)
  for(i in 1:NC) se[[i]] <- se.vector[covIndex$loc[which(covIndex$cat==i)]]
  return(list(cov=vars,se=se,ind=covIndex))
}



# variance of prob parameters for all models
# delta method is used for additive type models
OPG_p <- function(object,SE.type){
  Q <- extract(object,"Q")
  Qc <- extract(object,"Qc")
  J <- extract(object,"nitem")
  NC <- nrow(Q)
  Kj <- rowSums(Q)
  NG <- extract(object,"ngroup")
  linkfunc <- extract(object,"linkfunc")
  des <- extract(object,"designmatrix")
  m <- extract(object,"models_numeric")



  if(all(m<=2)&all(m>=0)){

    scorejh <- score_p(object) # a list with elements for each category
    np <- sapply(scorejh,ncol)
    IP.loc <- length(scorejh) - NG
    np <- np[-length(np)]

    if(SE.type == 1){
      scorejh <- scorejh[seq_len(IP.loc)]
      scorej <- vector("list",J)
      for(j in 1:J)  scorej[[j]] <- do.call(cbind,scorejh[which(Qc[,1]==j)])
      vars <- bdiag(lapply(scorej,inverse_crossprod))
    }else if(SE.type == 2){
      scorejh <- scorejh[seq_len(IP.loc)]
      vars <- inverse_crossprod(do.call(cbind,scorejh))
    }else if(SE.type==3){
      vars <- inverse_crossprod(do.call(cbind,scorejh))
      vars <- vars[1:sum(np),1:sum(np)]
    }
  }else{
    grad <- vector("list",NC)
    for (nc in 1:NC) {
      if (linkfunc[nc]=="identity") {
        grad[[nc]] <- des[[nc]]
      } else if (linkfunc[nc]=="logit") {
        grad[[nc]] <-diag(extract(object, "catprob.parm")[[nc]] * (1 - extract(object, "catprob.parm")[[nc]])) %*%des[[nc]]
        } else if (linkfunc[nc]=="log") {
        grad[[nc]] <- diag(extract(object, "catprob.parm")[[nc]]) %*% des[[nc]]
        }
    }
    g <- bdiag(grad)
    np <- sapply(grad,nrow)
    vars <- g %*% OPG_d(object,SE.type)$cov %*% t(g)
    }


  covIndex <- data.frame(item = rep(Qc[,1],np),
                         itemcat = rep(Qc[,2],np),
                         cat = rep(1:NC,np),
                         par = unlist(lapply(np,seq_len)),
                         loc = 1:sum(np),row.names = NULL)
  se.vector <- sqrt(diag(vars))
  se <- vector("list",NC)
  if(all(m<=2)&all(m>=0)){
  for(h in 1:NC) {

    if(m[h]==1) {
      se.pjh <- se.vector[covIndex$loc[which(covIndex$cat==h)]]
      se[[h]] <- c(rep(se.pjh[1],(2^Kj[h]-1)),se.pjh[2])
    }else if(m[h]==2){
      se.pjh <- se.vector[covIndex$loc[which(covIndex$cat==h)]]
      se[[h]] <- c(se.pjh[1],rep(se.pjh[2],(2^Kj[h]-1)))
    }else if(m[h]==0){
      se[[h]] <- se.vector[covIndex$loc[which(covIndex$cat==h)]]
    }
  }
    }else{
    for(h in 1:NC)
      se[[h]] <- se.vector[covIndex$loc[which(covIndex$cat==h)]]
  }
  return(list(cov=vars,se=se,ind=covIndex))
}


# SE_seq <- function(object,SE.type = 2, item=NULL,...){
#   dat <- extract(object,"dat")
#   Qc <- extract(object,"Qc")
#   Q <- extract(object,"Q")
#
#   J <- ncol(dat)
#   if (is.null(item)) item <- 1:J
#   if(SE.type>1) item <- 1:J
#   par.loc <- eta.loc(Q)
#   lik <- exp(indlogLik(object)) # posterior N x 2^K
#   scorej <- scorefunc2(object, item = item)
#
#   if(SE.type==1){
#     v <- lapply(scorej,inverse_crossprod)
#     if(length(v)>1) v <- bdiag(v)
#   }else if(SE.type==2){
#     v <- inverse_crossprod(do.call(cbind,scorej))
#   }else if(SE.type==3){
#     scopp <- (lik-lik[,1])/colSums(c(extract(object,"posterior.prob"))*t(post))
#     v <- inverse_crossprod(cbind(do.call(cbind,scorej),scopp[,-1]))
#   }
#
#   index <- data.frame(Cat=rep(1:length(rowSums(Q[Qc[,1]%in%item,,drop=FALSE]) ),2^rowSums(Q[Qc[,1]%in%item,,drop=FALSE])) )
#   index$Column <- seq_len(length(index$Cat))
#   index$Item <- rep(Qc[Qc[,1]%in%item,1],2^rowSums(Q[Qc[,1]%in%item,]))
#   return(list(var=v,index=index))
# }
