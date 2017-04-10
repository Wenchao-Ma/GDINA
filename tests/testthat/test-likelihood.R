test_that("checking likelihood calculation for dichotomous model", {

  Q <- matrix(c(1,0,
                0,1,
                1,1),ncol = 2,byrow = TRUE)
  dat <- matrix(c(1,0,1,
                  0,1,1,
                  0,0,1,
                  1,1,0,
                  1,0,0),ncol = 3,byrow = TRUE)
  itempar <- list(c(0.1,0.9),
                  c(0.2,0.8),
                  c(0.1,0.2,0.3,0.8))
  lc <- matrix(c(0.1,0.9,0.1,0.9,
                 0.2,0.2,0.8,0.8,
                 0.1,0.2,0.3,0.8),ncol = 4,byrow = TRUE)

  loglik_i <- dat%*%log(lc)+(1-dat)%*%log(1-lc)

  est <- GDINA(dat,Q,catprob.parm = itempar,maxitr = 0)

  expect_equivalent(extract(est,"loglikelihood.i"), loglik_i)

})


test_that("checking -2LL calculation", {

  Q <- matrix(c(1,0,
                0,1,
                1,1),ncol = 2,byrow = TRUE)
  dat <- matrix(c(1,0,1,
                  0,1,1,
                  0,0,1,
                  1,1,0,
                  1,0,0),ncol = 3,byrow = TRUE)
  itempar <- list(c(0.1,0.9),
                  c(0.2,0.8),
                  c(0.1,0.2,0.3,0.8))
  lc <- matrix(c(0.1,0.9,0.1,0.9,
                 0.2,0.2,0.8,0.8,
                 0.1,0.2,0.3,0.8),ncol = 4,byrow = TRUE)

  prior <- c(0.1,0.2,0.3,0.4)
  mprior <- matrix(prior,nrow = 5,ncol = 4,byrow = TRUE)

  loglik_i <- dat%*%log(lc)+(1-dat)%*%log(1-lc)

  LL <- sum(log(rowSums(mprior*exp(loglik_i))))
  missInd <- 1L - is.na(dat)
  LL2 <- Lik(mP = lc,
      mX = dat,
      vlogPrior = as.matrix(log(prior)),
      vgroup = as.matrix(rep(1,5)))$LL

  expect_equal(LL,LL2)

})

test_that("checking -2LL calculation - 2", {

  Q <- matrix(c(1,0,
                0,1,
                1,1),ncol = 2,byrow = TRUE)
  dat <- matrix(c(1,0,1,
                  0,1,1,
                  0,0,1,
                  1,1,0,
                  1,0,0),ncol = 3,byrow = TRUE)
  itempar <- list(c(0.1,0.9),
                  c(0.2,0.8),
                  c(0.1,0.2,0.3,0.8))
  lc <- matrix(c(0.1,0.9,0.1,0.9,
                 0.2,0.2,0.8,0.8,
                 0.1,0.2,0.3,0.8),ncol = 4,byrow = TRUE)

  prior <- c(0.1,0.2,0.3,0.4)
  mprior <- matrix(prior,nrow = 5,ncol = 4,byrow = TRUE)

  loglik_i <- dat%*%log(lc)+(1-dat)%*%log(1-lc)

  LL <- sum(log(rowSums(mprior*exp(loglik_i))))

  est <- GDINA(dat,Q,catprob.parm = itempar,maxitr = 0,att.prior = prior)

  expect_equal(LL,as.numeric(logLik(est)))

})

test_that("checking GDI calculation", {

  Q <- matrix(c(1,0,
                0,1,
                1,1),ncol = 2,byrow = TRUE)
  dat <- matrix(c(1,0,1,
                  0,1,1,
                  0,0,1,
                  1,1,0,
                  1,0,0),ncol = 3,byrow = TRUE)
  itempar <- list(c(0.1,0.9),
                  c(0.2,0.8),
                  c(0.1,0.2,0.3,0.8))

  w <- c(0.2,0.2,0.4,0.2)


  est <- GDINA(dat,Q,catprob.parm = itempar,maxitr = 0,att.prior = w,att.dist="fixed")
  w <- c(est$posterior.prob)
  estp <- t(extract(est,"expectedCorrect.LC")/extract(est,"expectedTotal.LC"))
  pbar <- colSums(estp*w)
  #11
  vsig11 <- colSums(estp^2*w)-pbar^2
  #01
  wp <- estp*w
  newp <- rbind((wp[1,]+wp[2,])/(w[1]+w[2]),
  (wp[3,]+wp[4,])/(w[3]+w[4]))
  neww <- c((w[1]+w[2]),(w[3]+w[4]))
  vsig01 <- colSums(newp^2*neww)-pbar^2
  #10
  newp <- rbind((wp[1,]+wp[3,])/(w[1]+w[3]),
                (wp[2,]+wp[4,])/(w[2]+w[4]))
  neww <- c((w[1]+w[3]),(w[2]+w[4]))
  vsig10 <- colSums(newp^2*neww)-pbar^2
  calc <- rbind(vsig10,vsig01,vsig11) #varsigma calculated by hand
  z=Qval(est,digits = 10) # by Qval function

  # print(calc)
  # print(z$varsigma)

  expect_equivalent(calc,as.matrix(z$varsigma))

})


test_that("checking person parameter estimation - MLE", {

  Q <- matrix(c(1,0,
                0,1,
                1,1),ncol = 2,byrow = TRUE)
  dat <- matrix(c(1,0,1,
                  0,1,1,
                  0,0,1,
                  1,1,0,
                  1,0,0),ncol = 3,byrow = TRUE)
  itempar <- list(c(0.3,0.9),
                  c(0.1,0.8),
                  c(0.1,0.2,0.3,0.8))
  lc <- matrix(c(0.3,0.9,0.3,0.9,
                 0.1,0.1,0.8,0.8,
                 0.1,0.2,0.3,0.8),ncol = 4,byrow = TRUE)

  loglik_i <- dat%*%log(lc)+(1-dat)%*%log(1-lc)

  est <- GDINA(dat,Q,catprob.parm = itempar,maxitr = 0)

  expect_equivalent(alpha(2)[apply(loglik_i,1,which.max),], as.matrix(personparm(est,"MLE")[,-3]))

})



test_that("checking LOGLIKE(I,J) for POLYTOMOUS data v1", {

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

  loglik_i <- #N x 2^K matrix
    (dat[,1]==0)%o%log(lc2[1,]) + # item 1
    (dat[,1]==1)%o%log(lc2[2,]) + # item 1
    (dat[,2]==0)%o%log(lc2[3,]) + # item 2
    (dat[,2]==1)%o%log(lc2[4,]) + # item 2
    (dat[,2]==2)%o%log(lc2[5,])  # item 2


  prior <- c(0.1,0.2,0.3,0.4)
  mprior <- matrix(prior,nrow = 5,ncol = 4,byrow = TRUE)
  yes.dat <- no.dat <- transition.code
  yes.dat[is.na(transition.code)] <- 0
  no.dat[is.na(transition.code)] <- 1
  loglik_i2 <- yes.dat%*%log(lc)+(1-no.dat)%*%log(1-lc)


  LL <- sum(log(rowSums(mprior*exp(loglik_i))))
  missInd <- 1L - is.na(transition.code)
  LL2 <- Lik(mP = lc,
             mX = transition.code,
             vlogPrior = as.matrix(log(prior)),
             vgroup = as.matrix(rep(1,5)))

  expect_equal(LL2$loglik,loglik_i2)

})


test_that("checking LOGLIKE(I,J) for POLYTOMOUS data v2", {

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

  loglik_i <- #N x 2^K matrix
    (dat[,1]==0)%o%log(lc2[1,]) + # item 1
    (dat[,1]==1)%o%log(lc2[2,]) + # item 1
    (dat[,2]==0)%o%log(lc2[3,]) + # item 2
    (dat[,2]==1)%o%log(lc2[4,]) + # item 2
    (dat[,2]==2)%o%log(lc2[5,])  # item 2


  prior <- c(0.1,0.2,0.3,0.4)
  mprior <- matrix(prior,nrow = 5,ncol = 4,byrow = TRUE)
  yes.dat <- no.dat <- transition.code
  yes.dat[is.na(transition.code)] <- 0
  no.dat[is.na(transition.code)] <- 1
  loglik_i2 <- yes.dat%*%log(lc)+(1-no.dat)%*%log(1-lc)


  LL <- sum(log(rowSums(mprior*exp(loglik_i))))
  missInd <- 1L - is.na(transition.code)
  LL2 <- Lik(mP = lc,
             mX = transition.code,
             vlogPrior = as.matrix(log(prior)),
             vgroup = as.matrix(rep(1,5)))

  expect_equal(LL2$loglik,loglik_i)

})

test_that("checking likelihood calculation with missing data", {

  Q <- matrix(c(1,0,
                0,1,
                1,1),ncol = 2,byrow = TRUE)
  dat <- data.frame(X1=c(1,0,0,1,1),
                    X2=c(NA,1,0,1,0),
                    X3=c(0,NA,1,0,1))
  itempar <- list(c(0.1,0.9),
                  c(0.2,0.8),
                  c(0.1,0.2,0.3,0.8))
  lc <- matrix(c(0.1,0.9,0.1,0.9,
                 0.2,0.2,0.8,0.8,
                 0.1,0.2,0.3,0.8),ncol = 4,byrow = TRUE)

  prior <- c(0.1,0.2,0.3,0.4)
  mprior <- matrix(prior,nrow = 5,ncol = 4,byrow = TRUE)
  yes.dat <- no.dat <- as.matrix(dat)
  yes.dat[is.na(dat)] <- 0
  no.dat[is.na(dat)] <- 1
  loglik_i <- yes.dat%*%log(lc)+(1-no.dat)%*%log(1-lc)

  LL <- sum(log(rowSums(mprior*exp(loglik_i))))
  missInd <- 1L - is.na(dat)
  LL2 <- Lik(mP = lc,
             mX = as.matrix(dat),
             vlogPrior = as.matrix(log(prior)),
             vgroup = as.matrix(rep(1,5)))$loglik

  expect_equal(loglik_i,LL2)

})
