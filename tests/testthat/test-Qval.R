

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


  est <- GDINA(dat,Q,catprob.parm = itempar,control=list(maxitr = 0),att.prior = w,att.dist="fixed")
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

