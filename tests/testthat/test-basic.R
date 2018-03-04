test_that("checking alpha patterns", {

  ret <- matrix(c(0,1,0,0,1,1,0,1,0,0,1,0,1,0,1,1,0,0,0,1,0,1,1,1),ncol = 3)

  expect_equivalent(attributepattern(K = 3), ret)

})

test_that("checking alpha patterns v2", {

  ret <- matrix(c(0,1,0,0,1,1,0,1,0,0,1,0,1,0,1,1,0,0,0,1,0,1,1,1),ncol = 3)

  Q <- ret[1:3,]

  expect_equivalent(attributepattern(Q = Q), ret)

})

test_that("checking polytomous alpha patterns", {

  ret <- matrix(c(0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2),ncol = 3)

  Q <- matrix(c(1,1,0,
                 2,1,1,
                 2,1,2),ncol = 3,byrow = TRUE)

  expect_equivalent(attributepattern(Q = Q), ret)

})
