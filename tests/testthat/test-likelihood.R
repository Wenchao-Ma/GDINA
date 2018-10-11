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

  est <- GDINA(dat,Q,catprob.parm = itempar,control=list(maxitr = 0))

  expect_equivalent(as.matrix(extract(est,"loglikelihood.i")), loglik_i)

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

  est <- GDINA(dat,Q,catprob.parm = itempar,control=list(maxitr = 0), att.prior = prior, att.dist = "fixed")
  expect_equal(-2*LL,est$testfit$Deviance,tolerance = 0.001)

})


test_that("checking likelihood calculation for multiple-group dichotomous model", {

  Q <- matrix(c(1,0,
                0,1,
                1,1),ncol = 2,byrow = TRUE)
  dat <- matrix(c(1,0,1,
                  0,1,1,
                  0,0,1,
                  1,0,1,
                  1,0,0),ncol = 3,byrow = TRUE)
  gr <- c(rep(1,3),rep(2,2))
  itempar <- list(c(0.1,0.9),
                  c(0.2,0.8),
                  c(0.1,0.2,0.3,0.8))
  lc <- matrix(c(0.1,0.9,0.1,0.9,
                 0.2,0.2,0.8,0.8,
                 0.1,0.2,0.3,0.8),ncol = 4,byrow = TRUE)

  loglik_i <- dat%*%log(lc)+(1-dat)%*%log(1-lc)
  prior <- matrix(c(0.1,0.2,0.3,0.4,0.4,0.3,0.2,0.1),ncol = 2)
  mprior <- t(prior[,c(1,1,1,2,2)])

  loglik_i <- dat%*%log(lc)+(1-dat)%*%log(1-lc)
  lik.i <- exp(loglik_i)
  joint <- lik.i*mprior
  logpost.i <- log(joint/rowSums(joint))
  LL <- sum(log(rowSums(mprior*exp(loglik_i))))


  est <- GDINA(dat,Q,catprob.parm = itempar,control=list(maxitr = 0), att.prior = prior, att.dist = "fixed",group = gr)

  expect_equivalent(as.matrix(indlogPost(est)), logpost.i)

})



test_that("checking -2LL for compact data", {

  Q <- matrix(c(1,0,
                0,1,
                1,1),ncol = 2,byrow = TRUE)
  dat <- matrix(c(1,0,1,
                  0,1,1,
                  0,0,1,
                  1,1,0,
                  1,0,0,
                  1,0,1,
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
  mprior <- matrix(prior,nrow = nrow(dat),ncol = 4,byrow = TRUE)

  loglik_i <- dat%*%log(lc)+(1-dat)%*%log(1-lc)

  LL <- sum(log(rowSums(mprior*exp(loglik_i))))

  compdat <-  matrix(c(1,0,1,
                       0,1,1,
                       0,0,1,
                       1,1,0,
                       1,0,0),ncol = 3,byrow = TRUE)

  LikNR(mpar = l2m(itempar),
        mX = as.matrix(compdat[,seq_len(ncol(dat))]),
        vlogPrior = as.matrix(log(prior)),
        vgroup = as.matrix(rep(1,nrow(dat))),
        mloc = eta(Q),
        weights = rep(2,5),
        simplify = 0)

  est <- GDINA(dat,Q,catprob.parm = itempar,control=list(maxitr = 0),att.prior = prior)

  expect_equal(LL,as.numeric(logLik(est)))

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

  est <- GDINA(dat,Q,catprob.parm = itempar,control=list(maxitr = 0))

  expect_equivalent(alpha2(2)[apply(loglik_i,1,which.max),], as.matrix(personparm(est,"MLE")[,-3]))

})



test_that("checking likelihood for POLYTOMOUS data", {

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


est <- GDINA(dat,Qc,sequential = T,catprob.parm = itempar,control=list(maxitr = 0))

expect_equivalent(as.matrix(indlogLik(est)),loglik_i)

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
  est <- GDINA(dat,Q,catprob.parm = itempar,att.prior  = prior,att.dist="fixed",control=list(maxitr = 0))

  expect_equivalent(loglik_i,as.matrix(indlogLik(est)))

})


test_that("checking posterior calculation with missing data", {

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
  lik.i <- exp(loglik_i)
  joint <- lik.i*matrix(prior,nrow(lik.i),length(prior),byrow = T)
  logpost.i <- log(joint/rowSums(joint))

  LL <- sum(log(rowSums(mprior*exp(loglik_i))))
  est <- GDINA(dat,Q,catprob.parm = itempar,att.prior  = prior,att.dist="fixed",control=list(maxitr = 0))

  expect_equivalent(logpost.i,as.matrix(indlogPost(est)))

})

