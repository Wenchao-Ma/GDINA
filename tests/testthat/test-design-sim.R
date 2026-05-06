test_that("simGDINA requires matrix-like Q", {

  expect_error(simGDINA(N = 10, Q = 1:4, gs.parm = matrix(c(0.1, 0.1), nrow = 1)), "Q-matrix must be")

})


test_that("simGDINA validates required and exclusive item parameter inputs", {

  q <- matrix(c(1, 0,
                0, 1), nrow = 2, byrow = TRUE)
  gs <- matrix(c(0.10, 0.20,
                 0.15, 0.10), nrow = 2, byrow = TRUE)

  expect_error(simGDINA(N = 10, Q = q), "Item parameters must be specified")
  expect_error(
    simGDINA(N = 10, Q = q, catprob.parm = list(c(0.5, 0.5), c(0.5, 0.5)), delta.parm = list(c(0, 0), c(0, 0))),
    "can only be specified"
  )

})


test_that("simGDINA validates type and Q matrix checks", {

  q <- matrix(c(1, 0,
                0, 1), nrow = 2, byrow = TRUE)
  gs <- matrix(c(0.10, 0.20,
                 0.15, 0.10), nrow = 2, byrow = TRUE)

  expect_error(
    simGDINA(N = 10, Q = q, gs.parm = gs, gs.args = list(type = "badtype")),
    "type must be either random or equal"
  )

  q0 <- matrix(0, nrow = 2, ncol = 2)
  expect_error(simGDINA(N = 10, Q = q0, gs.parm = gs), "Some rows of the Q-matrix contain only 0s")

})


test_that("simGDINA validates N and Q value domain", {

  q <- matrix(c(1, 0,
                0, 1), nrow = 2, byrow = TRUE)

  expect_error(simGDINA(N = -1, Q = q, gs.parm = matrix(c(0.1, 0.1, 0.1, 0.1), nrow = 2)), "N must be")
  expect_error(simGDINA(N = 10, Q = matrix(c(1, -1, 0, 1), nrow = 2), gs.parm = matrix(c(0.1, 0.1, 0.1, 0.1), nrow = 2)), "Q matrix can only contain")

})


test_that("simGDINA validates parameter object types", {

  q <- matrix(c(1, 0,
                0, 1), nrow = 2, byrow = TRUE)

  expect_error(simGDINA(N = 10, Q = q, gs.parm = c(0.1, 0.2)), "gs.parm must be")
  expect_error(simGDINA(N = 10, Q = q, catprob.parm = c(0.1, 0.2)), "itemprob.parm must be")
  expect_error(simGDINA(N = 10, Q = q, delta.parm = matrix(c(0, 0), nrow = 1)), "delta.parm must be")

})


test_that("simGDINA validates sequential Q row checks", {

  q_seq0 <- matrix(c(1, 1, 0, 0,
                     2, 1, 0, 0), nrow = 2, byrow = TRUE)
  expect_error(
    simGDINA(N = 10, Q = q_seq0, sequential = TRUE, catprob.parm = list(c(0.2, 0.8), c(0.3, 0.7))),
    "Some rows of the Q-matrix contain only 0s"
  )

})
