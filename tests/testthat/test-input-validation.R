base_inputcheck_args <- function() {
  list(
    dat = matrix(c(0, 1,
                   1, 0,
                   1, 1), ncol = 2, byrow = TRUE),
    Q = matrix(c(1, 0,
                 0, 1), ncol = 2, byrow = TRUE),
    model = c(0, 0),
    sequential = FALSE,
    att.dist = "saturated",
    no.bugs = 0,
    verbose = 0,
    catprob.parm = NULL,
    mono.constraint = c(TRUE, TRUE),
    loglinear = 1,
    att.prior = c(0.25, 0.25, 0.25, 0.25),
    lower.p = c(0, 0),
    upper.p = c(1, 1),
    att.str = NULL,
    nstarts = 1,
    conv.crit = 1e-4,
    maxitr = 1
  )
}


test_that("model.table provides expected mapping metadata", {

  tbl <- GDINA:::model.table()

  expect_equal(ncol(tbl), 5)
  expect_equal(nrow(tbl), 12)
  expect_true(all(c("model.char", "model.num", "linkf.num", "linkf.char", "rule") %in% names(tbl)))

  gdina_row <- tbl[tbl$model.char == "GDINA", ]
  expect_equal(gdina_row$model.num, 0)
  expect_equal(gdina_row$linkf.num, 1)
  expect_equal(gdina_row$rule, 0)

  msdina_row <- tbl[tbl$model.char == "MSDINA", ]
  expect_equal(msdina_row$rule, 4)

})


test_that("model2numeric and model2character handle replication and pass-through", {

  expect_equal(GDINA:::model2numeric(1, J = 3), rep(1, 3))
  expect_equal(GDINA:::model2character("DINA", J = 2), rep("DINA", 2))

  vals <- c(-3, 0, 2, 6)
  expect_equal(GDINA:::model2numeric(vals), vals)
  expect_equal(GDINA:::model2character(vals), c("LOGGDINA", "GDINA", "DINO", "MSDINA"))

})


test_that("model2rule and model2linkfunc helper variants agree", {

  expect_equal(GDINA:::model2rule.j("DINA"), 1)
  expect_equal(GDINA:::model2rule.j(2), 2)

  expect_equal(GDINA:::model2linkfunc.j("LOGITGDINA"), 2)
  expect_equal(GDINA:::model2linkfunc.j(-3), 3)

  expect_equal(unname(GDINA:::model2rule(c("GDINA", "DINA", "DINO", "ACDM"))), c(0, 1, 2, 3))
  expect_equal(unname(GDINA:::model2linkfunc(c("GDINA", "LOGITGDINA", "RRUM"))), c(1, 2, 3))

})


test_that("inputcheck accepts a valid minimal setup", {

  args <- base_inputcheck_args()

  expect_no_error(do.call(GDINA:::inputcheck, args))

})


test_that("inputcheck rejects invalid sequential, dimensions, and probability bounds", {

  args <- base_inputcheck_args()
  args$sequential <- 1
  expect_error(do.call(GDINA:::inputcheck, args), "sequential must be logical")

  args <- base_inputcheck_args()
  args$dat <- c(0, 1, 1)
  expect_error(do.call(GDINA:::inputcheck, args), "Data must be a matrix or data frame")

  args <- base_inputcheck_args()
  args$Q <- matrix(c(1, 0,
                     0, 1,
                     1, 1), ncol = 2, byrow = TRUE)
  expect_error(do.call(GDINA:::inputcheck, args), "does not match")

  args <- base_inputcheck_args()
  args$lower.p <- c(0.7, 0.7)
  args$upper.p <- c(0.6, 0.9)
  expect_error(do.call(GDINA:::inputcheck, args), "lower.p must be less than upper.p")

  args <- base_inputcheck_args()
  args$upper.p <- c(1.2, 0.8)
  expect_error(do.call(GDINA:::inputcheck, args), "upper.p must range from 0 to 1")

})


test_that("inputcheck enforces att.str and mono.constraint constraints", {

  args <- base_inputcheck_args()
  args$mono.constraint <- c(TRUE)
  expect_no_error(do.call(GDINA:::inputcheck, args))

  args <- base_inputcheck_args()
  args$mono.constraint <- c(TRUE, FALSE, TRUE)
  expect_error(do.call(GDINA:::inputcheck, args), "Length of mono.constraint")

  args <- base_inputcheck_args()
  args$att.str <- matrix(c(0, 1,
                           1, 0), ncol = 2, byrow = TRUE)
  args$att.dist <- "independent"
  expect_error(do.call(GDINA:::inputcheck, args), "Independent structure is not allowed")

  args <- base_inputcheck_args()
  args$att.str <- matrix(c(0, 1,
                           1, 0), ncol = 2, byrow = TRUE)
  args$att.dist <- "higher.order"
  expect_error(do.call(GDINA:::inputcheck, args), "Higher-order structure is not allowed")

})


test_that("inputcheck catches sequential and model-specific invalid states", {

  args <- base_inputcheck_args()
  args$dat <- matrix(c(2, 1,
                       1, 0,
                       1, 1), ncol = 2, byrow = TRUE)
  expect_error(do.call(GDINA:::inputcheck, args), "set sequential = TRUE")

  args <- base_inputcheck_args()
  args$model <- c(7, 7)
  args$mono.constraint <- c(TRUE, TRUE)
  expect_error(do.call(GDINA:::inputcheck, args), "Monotonicity is not allowed")

  args <- base_inputcheck_args()
  args$model <- c(8, 8)
  args$no.bugs <- 3
  expect_error(do.call(GDINA:::inputcheck, args), "no.bugs is not correctly specified")

})
