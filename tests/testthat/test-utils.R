test_that("integer helper utilities classify values correctly", {

  expect_true(GDINA:::is.wholenumber(2 + 1e-9))
  expect_false(GDINA:::is.wholenumber(2.1))

  expect_true(all(GDINA:::is.nonNegativeInteger(c(0, 2, 3))))
  expect_false(any(GDINA:::is.nonNegativeInteger(c(-1, 1.2))))

  expect_true(GDINA:::is.positiveInteger(3))
  expect_false(GDINA:::is.positiveInteger(0))

})


test_that("missingMsg emits clear argument-missing error", {

  expect_error(GDINA:::missingMsg("Q"), '"Q" argument is missing.')

})


test_that("model.transform validates and normalizes model inputs", {

  expect_equal(GDINA:::model.transform("DINA", J = 3), rep(1, 3))
  expect_equal(GDINA:::model.transform(c(0, 1, 2), J = 3), c(0, 1, 2))

  expect_error(
    GDINA:::model.transform(c("GDINA", "DINA"), J = 3),
    "same length as the test"
  )
  expect_error(
    GDINA:::model.transform("NOT_A_MODEL", J = 1),
    "can only be"
  )

})


test_that("model conversion helpers stay consistent", {

  expect_equal(GDINA:::model2numeric(c("GDINA", "DINO", "RRUM")), c(0, 2, 5))
  expect_equal(GDINA:::model2character(c(0, 2, 5)), c("GDINA", "DINO", "RRUM"))

  expect_equal(GDINA:::model2numeric("DINA", J = 2), rep(1, 2))
  expect_equal(GDINA:::model2character(1, J = 2), rep("DINA", 2))

  expect_equal(unname(GDINA:::model2rule(c("DINA", "DINO", "ACDM", "MSDINA"))), c(1, 2, 3, 4))
  expect_equal(unname(GDINA:::model2linkfunc(c("GDINA", "LOGITGDINA", "LOGGDINA", "RRUM"))), c(1, 2, 3, 3))

})


test_that("linkf.numeric handles null, scalar, vector and invalid lengths", {

  models <- c("GDINA", "LOGITGDINA", "LOGGDINA")

  expect_equal(unname(GDINA:::linkf.numeric(NULL, model.vector = models)), c(1, 2, 3))
  expect_equal(GDINA:::linkf.numeric("logit", model.vector = c("GDINA", "DINA")), c(2, 2))
  expect_equal(GDINA:::linkf.numeric(c("identity", "log"), model.vector = c("GDINA", "RRUM")), c(1, 3))

  expect_error(
    GDINA:::linkf.numeric(c("identity", "logit"), model.vector = c("GDINA", "DINA", "DINO")),
    "Length of linkfunc is not correct"
  )

})


test_that("inputcheck.sim validates simulation setup", {

  Q <- matrix(c(1, 0,
                0, 1), ncol = 2, byrow = TRUE)
  gs <- matrix(0.2, nrow = nrow(Q), ncol = 2)

  expect_no_error(
    GDINA:::inputcheck.sim(
      N = 10,
      Q = Q,
      gs.parm = gs,
      sequential = FALSE,
      type = "random"
    )
  )

  expect_error(
    GDINA:::inputcheck.sim(
      N = -1,
      Q = Q,
      gs.parm = gs,
      sequential = FALSE,
      type = "random"
    ),
    "N must be negative integer"
  )

  expect_error(
    GDINA:::inputcheck.sim(
      N = 10,
      Q = Q,
      gs.parm = gs,
      sequential = FALSE,
      type = "invalid"
    ),
    "type must be either random or equal"
  )

  expect_error(
    GDINA:::inputcheck.sim(
      N = 10,
      Q = Q,
      gs.parm = NULL,
      catprob.parm = NULL,
      delta.parm = NULL,
      sequential = FALSE,
      type = "random"
    ),
    "Item parameters must be specified"
  )

  expect_error(
    GDINA:::inputcheck.sim(
      N = 10,
      Q = Q,
      gs.parm = gs,
      catprob.parm = list(c(0.2, 0.8), c(0.3, 0.7)),
      sequential = FALSE,
      type = "random"
    ),
    "Item parameters can only be specified"
  )

})
