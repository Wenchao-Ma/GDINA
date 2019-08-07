

test_that("checking itemfit calculation", {

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

  x <- itemfit(est,randomseed = 123456)
  expect_equal(x$r$rstat[2],.9517)
  expect_equal(x$p$pstat,c(.0556,.1284,.2503))
})
