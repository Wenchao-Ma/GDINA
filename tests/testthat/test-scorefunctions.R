test_that("checking GDINA score function LOGLIKE(I,J) for POLYTOMOUS data", {

  Qc <- matrix(c(1,1,1,0,
                2,1,0,1,
                2,2,1,1),ncol = 4,byrow = TRUE)
  dat <- matrix(c(1,0,
                  0,1,
                  0,2,
                  1,1,
                  1,2),ncol = 2,byrow = TRUE)
  transition.code <- seq_coding(dat,Qc)
  # > transition.code
  # tmp tmp tmp
  # [1,]   1   0  NA
  # [2,]   0   1   0
  # [3,]   0   1   1
  # [4,]   1   1   0
  # [5,]   1   1   1
  itempar <- list(c(0.1,0.9), #ITEM 1
                  c(0.2,0.7), #ITEM 2 CAT 1
                  c(0.1,0.2,0.3,0.8)) #ITEM 2 CAT 2
  # processing function for each category
  lc <- matrix(c(0.1,0.9,0.1,0.9, #P(X1=1|alpha_c)
                 0.2,0.2,0.7,0.7, #S(X2=1|alpha_c)
                 0.1,0.2,0.3,0.8  #S(X2=2|alpha_c)
                 ),ncol = 4,byrow = TRUE)
  lc2 <- matrix(c(0.9,0.1,0.9,0.1, #P(X1=0|alpha_c)
                  0.1,0.9,0.1,0.9, #P(X1=1|alpha_c)
                  0.8,0.8,0.3,0.3, #P(X2=0|alpha_c)
                 0.18,0.16,.49,.14,  #P(X2=1|alpha_c)
                 0.02,0.04,0.21,0.56 #P(X2=2|alpha_c)
                 ),ncol = 4,byrow = TRUE)

  scorefuncj <- function(Xj, Qc,catprob.matrix,logpost,j, ...){

    Q <- Qc[,-c(1:2)]

    J <- length(unique(Qc[,1]))
    par.loc <- eta(Q)
    post <- exp(logpost) # posterior N x 2^K

    # urP <- itemparm(object,"itemprob",digits = 9) #unconditional reduced prob.
    C <- table(Qc[,1]) # number of categories
    S0 <- rep(1:J,C+1) # loc for each item including cat 0

    LC.Prob <- sequP(as.matrix(par.loc),as.matrix(catprob.matrix),C)
    c.LC.Prob <- LC.Prob$cPr # no cat 0
    u.LC.Prob <- LC.Prob$uPr # including cat 0

    sj <- c.LC.Prob[which(Qc[,1]==j),,drop=FALSE] # processing function for item j
    locj <- par.loc[which(Qc[,1]==j),,drop=FALSE] # loc for item j
    Pj <- u.LC.Prob[which(S0==j),][-1,,drop=FALSE] # drop cat 0 - P(Xij=h|alpha)
    Pj0 <- u.LC.Prob[which(S0==j),][1,] # P(Xij=0|alpha)
    Dj <- list() #I(Xij=h)/P(Xij=h|alpha) - I(Xij=0)/P(Xij=0|alpha)
    for(h in seq_len(C[j])) Dj[[h]] <- outer((Xj==h),Pj[h,],"/")-outer((Xj==0),Pj0,"/")
    scorej <- NULL
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
      for(lj in sort(unique(locj[k,]))) scorej <- cbind(scorej,rowSums(tmp[,which(locj[k,]==lj),drop=FALSE]))

    }


    return(scorej)
  }


  loglik_i <- #N x 2^K matrix
    (dat[,1]==0)%o%log(lc2[1,]) + # item 1
    (dat[,1]==1)%o%log(lc2[2,]) + # item 1
    (dat[,2]==0)%o%log(lc2[3,]) + # item 2
    (dat[,2]==1)%o%log(lc2[4,]) + # item 2
    (dat[,2]==2)%o%log(lc2[5,])  # item 2

  prior <- c(0.1,0.2,0.3,0.4)
  mprior <- matrix(prior,nrow = 5,ncol = 4,byrow = TRUE)
  post <- mprior*exp(loglik_i)
  post <- post/rowSums(post) # standardized posterior
  sco <- scorefuncj(Xj=dat[,2], Qc,catprob.matrix=l2m(itempar),logpost=log(post),j=2)

  Xj=dat[,2]
  scorj1c <- post*((Xj>=1)%o%(1/lc[2,])-(Xj==0)%o%(1/(1-lc[2,])))
scoj1 <- cbind(rowSums(scorj1c[,1:2]),rowSums(scorj1c[,3:4]))

scorj2c <- post*((Xj>=2)%o%(1/lc[3,])-(Xj==1)%o%(1/(1-lc[3,])))
  expect_equal(scoj1,sco[,1:2])

})

test_that("checking GDINA score function LOGLIKE(I,J) for POLYTOMOUS data v2", {

  Qc <- matrix(c(1,1,1,0,
                 2,1,0,1,
                 2,2,1,1),ncol = 4,byrow = TRUE)
  dat <- matrix(c(1,0,
                  0,1,
                  0,2,
                  1,1,
                  1,2),ncol = 2,byrow = TRUE)
  transition.code <- seq_coding(dat,Qc)
  # > transition.code
  # tmp tmp tmp
  # [1,]   1   0  NA
  # [2,]   0   1   0
  # [3,]   0   1   1
  # [4,]   1   1   0
  # [5,]   1   1   1
  parloc <- eta(Qc[,-c(1:2)])
  itempar <- list(c(0.1,0.9), #ITEM 1
                  c(0.2,0.7), #ITEM 2 CAT 1
                  c(0.1,0.2,0.3,0.8)) #ITEM 2 CAT 2
  # processing function for each category
  lc <- matrix(c(0.1,0.9,0.1,0.9, #P(X1=1|alpha_c)
                 0.2,0.2,0.7,0.7, #S(X2=1|alpha_c)
                 0.1,0.2,0.3,0.8  #S(X2=2|alpha_c)
  ),ncol = 4,byrow = TRUE)
  lc2 <- matrix(c(0.9,0.1,0.9,0.1, #P(X1=0|alpha_c)
                  0.1,0.9,0.1,0.9, #P(X1=1|alpha_c)
                  0.8,0.8,0.3,0.3, #P(X2=0|alpha_c)
                  0.18,0.16,.49,.14,  #P(X2=1|alpha_c)
                  0.02,0.04,0.21,0.56 #P(X2=2|alpha_c)
  ),ncol = 4,byrow = TRUE)

  loglik_i <- #N x 2^K matrix
    (dat[,1]==0)%o%log(lc2[1,]) + # item 1
    (dat[,1]==1)%o%log(lc2[2,]) + # item 1
    (dat[,2]==0)%o%log(lc2[3,]) + # item 2
    (dat[,2]==1)%o%log(lc2[4,]) + # item 2
    (dat[,2]==2)%o%log(lc2[5,])  # item 2

  scorefuncj <- function(Xj, Qc,catprob.matrix,logpost,j, ...){

    Q <- Qc[,-c(1:2)]

    J <- length(unique(Qc[,1]))
    par.loc <- eta(Q)
    post <- exp(logpost) # posterior N x 2^K

    # urP <- itemparm(object,"itemprob",digits = 9) #unconditional reduced prob.
    C <- table(Qc[,1]) # number of categories
    S0 <- rep(1:J,C+1) # loc for each item including cat 0

    LC.Prob <- sequP(as.matrix(par.loc),as.matrix(catprob.matrix),C)
    c.LC.Prob <- LC.Prob$cPr # no cat 0
    u.LC.Prob <- LC.Prob$uPr # including cat 0

    sj <- c.LC.Prob[which(Qc[,1]==j),,drop=FALSE] # processing function for item j
    locj <- par.loc[which(Qc[,1]==j),,drop=FALSE] # loc for item j
    Pj <- u.LC.Prob[which(S0==j),][-1,,drop=FALSE] # drop cat 0 - P(Xij=h|alpha)
    Pj0 <- u.LC.Prob[which(S0==j),][1,] # P(Xij=0|alpha)
    Dj <- list() #I(Xij=h)/P(Xij=h|alpha) - I(Xij=0)/P(Xij=0|alpha)
    for(h in seq_len(C[j])) Dj[[h]] <- outer((Xj==h),Pj[h,],"/")-outer((Xj==0),Pj0,"/")
    scorej <- NULL
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
      for(lj in sort(unique(locj[k,]))) scorej <- cbind(scorej,rowSums(tmp[,which(locj[k,]==lj),drop=FALSE]))

    }


    return(scorej)
  }



  scorefuncj2 <- function(Xj,parloc.j,catprob.j,logpost, ...){
    Xj[is.na(Xj)] <- -1 # remove NA from the data
    post <- exp(logpost) # posterior N x 2^K
    scorej <- vector("list",nrow(parloc.j)) #I(Xij>=x)/sjx - I(Xij=x-1)/[1-sjx]
    for(r in 1:length(scorej)){
      scorej[[r]] <- aggregateCol(post,parloc.j[r,])*(outer((Xj>=r),catprob.j[[r]],"/")-outer((Xj==r-1),(1-catprob.j[[r]]),"/"))
    }

    return(scorej)
  }


  prior <- c(0.1,0.2,0.3,0.4)
  mprior <- matrix(prior,nrow = 5,ncol = 4,byrow = TRUE)
  post <- mprior*exp(loglik_i)
  post <- post/rowSums(post) # standardized posterior
  sco <- scorefuncj(Xj=dat[,2], Qc,catprob.matrix=l2m(itempar),logpost=log(post),j=2)
  sco2 <- scorefuncj2(Xj = dat[,2],parloc[2:3,],catprob.j = itempar[2:3],logpost=log(post),j=2)
  Xj=dat[,2]
  scorj1c <- post*((Xj>=1)%o%(1/lc[2,])-(Xj==0)%o%(1/(1-lc[2,])))
  scoj1 <- cbind(rowSums(scorj1c[,1:2]),rowSums(scorj1c[,3:4]))

  scorj2c <- post*((Xj>=2)%o%(1/lc[3,])-(Xj==1)%o%(1/(1-lc[3,])))
  expect_equal(scoj1,sco2[[1]])

})



