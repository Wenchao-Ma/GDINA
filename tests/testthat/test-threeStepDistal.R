test_that("threeStepDistal fits binary, ordinal, and continuous distal outcomes", {
  skip_on_cran()
  set.seed(2026)

  Q <- matrix(c(
    1, 0,
    1, 0,
    0, 1,
    0, 1,
    1, 1,
    1, 1
  ), byrow = TRUE, ncol = 2)
  gs <- data.frame(guess = rep(0.2, nrow(Q)), slip = rep(0.2, nrow(Q)))
  sim <- simGDINA(300, Q, gs.parm = gs, model = "DINA")
  fit <- GDINA(sim$dat, sim$Q, model = "DINA", verbose = 0)
  alpha <- extract(sim, "attribute")

  y_binary <- rbinom(300, 1, plogis(-0.5 + 1.1 * alpha[, 1] - 0.7 * alpha[, 2]))
  eta_ordinal <- -0.2 + 0.8 * alpha[, 1] - 0.6 * alpha[, 2]
  p_low <- stats::plogis(0.1 - eta_ordinal)
  p_medium <- stats::plogis(1.0 - eta_ordinal) - p_low
  u <- runif(300)
  y_ordinal <- ordered(
    ifelse(u < p_low, "low", ifelse(u < p_low + p_medium, "medium", "high")),
    levels = c("low", "medium", "high")
  )
  y_continuous <- 0.4 + 0.9 * alpha[, 1] - 0.6 * alpha[, 2] + rnorm(300, sd = 0.8)

  out_att <- ThreeStepDistal(fit, y_binary, level = "attribute", attribute = 1:2, method = c("ML", "BCH"))
  out_ord <- ThreeStepDistal(fit, y_ordinal, level = "profile", method = c("ML", "BCH"))
  out_prof <- ThreeStepDistal(fit, y_continuous, level = "profile", method = c("ML", "BCH"))

  expect_s3_class(out_att, "ThreeStepDistal")
  expect_equal(out_att$outcome_type, "binary")
  expect_equal(out_att$results$attribute, 1:2)
  expect_true(all(c("naive", "ML", "BCH", "bch_weights") %in% names(out_att$results)))
  expect_equal(out_att$results$ML$table$term, c("(Intercept)", "Attribute 1", "Attribute 2"))
  expect_equal(out_att$results$ML$table$odds.ratio, exp(out_att$results$ML$table$estimate), tolerance = 1e-10)
  expect_equal(out_att$results$ML$table$or.conf.low, exp(out_att$results$ML$table$conf.low), tolerance = 1e-10)
  expect_equal(out_att$results$ML$table$or.conf.high, exp(out_att$results$ML$table$conf.high), tolerance = 1e-10)

  expect_s3_class(out_ord, "ThreeStepDistal")
  expect_equal(out_ord$outcome_type, "ordinal")
  expect_true(all(c("threshold", "estimate", "conf.low", "conf.high") %in% names(out_ord$results$ML$threshold_table)))
  expect_equal(out_ord$results$ML$table$odds.ratio, exp(out_ord$results$ML$table$estimate), tolerance = 1e-10)

  expect_s3_class(out_prof, "ThreeStepDistal")
  expect_equal(out_prof$outcome_type, "continuous")
  expect_equal(ncol(out_prof$results$ML$posterior), nrow(extract(fit, "attributepattern")))
  expect_true(is.numeric(out_prof$results$ML$residual_sd))
})

test_that("threeStepDistal naive attribute fit matches glm on hard classifications", {
  skip_on_cran()
  set.seed(123)

  Q <- matrix(c(
    1, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 1, 0,
    0, 0, 1,
    0, 0, 1,
    1, 1, 0,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  ), byrow = TRUE, ncol = 3)
  gs <- data.frame(guess = rep(0.3, nrow(Q)), slip = rep(0.24, nrow(Q)))
  sim <- simGDINA(500, Q, gs.parm = gs, model = "GDINA")
  fit <- GDINA(sim$dat, sim$Q, verbose = 0)
  alpha <- extract(sim, "attribute")
  att_map <- extract(fit, "attributepattern")[max.col(extract(fit, "logposterior.i")), , drop = FALSE]

  y_binary <- rbinom(500, 1, plogis(-0.8 + alpha[, 1] - 0.9 * alpha[, 2] + 0.7 * alpha[, 3]))

  out <- ThreeStepDistal(fit, y_binary, level = "attribute", attribute = 1:2)
  truth <- glm(y_binary ~ att_map[, 1] + att_map[, 2], family = binomial())

  expect_equal(unname(out$results$naive$coefficients), unname(coef(truth)), tolerance = 1e-6)
})

test_that("threeStepDistal naive continuous profile fit matches glm on hard classifications", {
  skip_on_cran()
  set.seed(124)

  Q <- matrix(c(
    1, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 1, 0,
    0, 0, 1,
    0, 0, 1,
    1, 1, 0,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  ), byrow = TRUE, ncol = 3)
  gs <- data.frame(guess = rep(0.3, nrow(Q)), slip = rep(0.24, nrow(Q)))
  sim <- simGDINA(500, Q, gs.parm = gs, model = "GDINA")
  fit <- GDINA(sim$dat, sim$Q, verbose = 0)
  alpha <- extract(sim, "attribute")
  pattern <- extract(fit, "attributepattern")
  profile_labels <- apply(pattern, 1, paste, collapse = "")
  profile_map <- factor(
    apply(pattern[max.col(extract(fit, "logposterior.i")), , drop = FALSE], 1, paste, collapse = ""),
    levels = profile_labels
  )
  profile_map <- stats::relevel(profile_map, ref = tail(profile_labels, 1))

  y_continuous <- 0.3 + 0.9 * alpha[, 1] - 1.1 * alpha[, 2] + 0.8 * alpha[, 3] + rnorm(500, sd = 1.1)

  out <- ThreeStepDistal(fit, y_continuous, level = "profile", method = c("ML", "BCH"))
  truth <- glm(y_continuous ~ profile_map, family = gaussian())

  expect_equal(unname(out$results$naive$coefficients), unname(coef(truth)), tolerance = 1e-8)
})

test_that("threeStepDistal naive ordinal profile fit matches polr on hard classifications", {
  skip_on_cran()
  set.seed(126)

  Q <- matrix(c(
    1, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 1, 0,
    0, 0, 1,
    0, 0, 1,
    1, 1, 0,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  ), byrow = TRUE, ncol = 3)
  gs <- data.frame(guess = rep(0.3, nrow(Q)), slip = rep(0.24, nrow(Q)))
  sim <- simGDINA(500, Q, gs.parm = gs, model = "GDINA")
  fit <- GDINA(sim$dat, sim$Q, verbose = 0)
  alpha <- extract(sim, "attribute")
  pattern <- extract(fit, "attributepattern")
  profile_labels <- apply(pattern, 1, paste, collapse = "")
  profile_map <- factor(
    apply(pattern[max.col(extract(fit, "logposterior.i")), , drop = FALSE], 1, paste, collapse = ""),
    levels = profile_labels
  )
  profile_map <- stats::relevel(profile_map, ref = tail(profile_labels, 1))

  eta <- -0.4 + 1.1 * alpha[, 1] - 0.8 * alpha[, 2] + 0.5 * alpha[, 3]
  p_low <- stats::plogis(0.2 - eta)
  p_medium <- stats::plogis(1.0 - eta) - p_low
  u <- runif(500)
  y_ordinal <- ordered(
    ifelse(u < p_low, "low", ifelse(u < p_low + p_medium, "medium", "high")),
    levels = c("low", "medium", "high")
  )

  out <- ThreeStepDistal(fit, y_ordinal, level = "profile", method = "ML")
  truth <- MASS::polr(y_ordinal ~ profile_map, method = "logistic", Hess = TRUE)

  expect_equal(unname(out$results$naive$coefficients), unname(coef(truth)), tolerance = 1e-4)
  expect_equal(unname(out$results$naive$thresholds), unname(truth$zeta), tolerance = 1e-4)
})

test_that("threeStepDistal naive categorical profile fit matches baseline-category glm comparisons", {
  skip_on_cran()
  set.seed(125)

  Q <- matrix(c(
    1, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 1, 0,
    0, 0, 1,
    0, 0, 1,
    1, 1, 0,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  ), byrow = TRUE, ncol = 3)
  gs <- data.frame(guess = rep(0.3, nrow(Q)), slip = rep(0.24, nrow(Q)))
  sim <- simGDINA(500, Q, gs.parm = gs, model = "GDINA")
  fit <- GDINA(sim$dat, sim$Q, verbose = 0)
  alpha <- extract(sim, "attribute")
  pattern <- extract(fit, "attributepattern")
  profile_labels <- apply(pattern, 1, paste, collapse = "")
  profile_map <- factor(
    apply(pattern[max.col(extract(fit, "logposterior.i")), , drop = FALSE], 1, paste, collapse = ""),
    levels = profile_labels
  )
  profile_map <- stats::relevel(profile_map, ref = tail(profile_labels, 1))

  latent_score <- -0.4 + 1.1 * alpha[, 1] - 0.8 * alpha[, 2] + 0.5 * alpha[, 3] + rnorm(500, sd = 0.8)
  cuts <- quantile(latent_score, probs = c(1 / 3, 2 / 3))
  y_categorical <- cut(
    latent_score,
    breaks = c(-Inf, cuts, Inf),
    labels = c("low", "medium", "high")
  )

  out <- ThreeStepDistal(fit, y_categorical, level = "profile", method = "ML")
  low_vs_high <- y_categorical %in% c("low", "high")
  glm_low <- glm(
    as.integer(y_categorical[low_vs_high] == "low") ~ profile_map[low_vs_high],
    family = binomial()
  )
  medium_vs_high <- y_categorical %in% c("medium", "high")
  glm_medium <- glm(
    as.integer(y_categorical[medium_vs_high] == "medium") ~ profile_map[medium_vs_high],
    family = binomial()
  )

  expect_equal(unname(out$results$naive$coefficients[, "low"]), unname(coef(glm_low)), tolerance = 1e-4)
  expect_equal(unname(out$results$naive$coefficients[, "medium"]), unname(coef(glm_medium)), tolerance = 1e-4)
})

test_that("threeStepDistal fits categorical profile distal outcomes", {
  skip_on_cran()
  set.seed(2027)

  Q <- matrix(c(
    1, 0,
    1, 0,
    0, 1,
    0, 1,
    1, 1,
    1, 1
  ), byrow = TRUE, ncol = 2)
  gs <- data.frame(guess = rep(0.18, nrow(Q)), slip = rep(0.18, nrow(Q)))
  sim <- simGDINA(250, Q, gs.parm = gs, model = "DINA")
  fit <- GDINA(sim$dat, sim$Q, model = "DINA", verbose = 0)
  alpha <- extract(sim, "attribute")

  latent_score <- -0.3 + 0.8 * alpha[, 1] - 0.5 * alpha[, 2] + rnorm(250, sd = 0.7)
  cut_points <- quantile(latent_score, probs = c(1 / 3, 2 / 3))
  y_categorical <- cut(latent_score, breaks = c(-Inf, cut_points, Inf), labels = c("low", "medium", "high"))

  out <- ThreeStepDistal(fit, y_categorical, level = "profile", method = "ML")
  out_nominal <- ThreeStepDistal(fit, y_categorical, level = "profile", outcome_type = "nominal", method = "ML")

  expect_s3_class(out, "ThreeStepDistal")
  expect_equal(out$outcome_type, "categorical")
  expect_equal(out_nominal$outcome_type, "categorical")
  expect_true(all(c("outcome_level", "term", "estimate", "odds.ratio", "or.conf.low", "or.conf.high") %in% names(out$results$ML$table)))
  expect_equal(out$results$ML$table$odds.ratio, exp(out$results$ML$table$estimate), tolerance = 1e-10)
  expect_equal(out$results$ML$table$or.conf.low, exp(out$results$ML$table$conf.low), tolerance = 1e-10)
  expect_equal(out$results$ML$table$or.conf.high, exp(out$results$ML$table$conf.high), tolerance = 1e-10)
  expect_equal(sort(unique(out$results$ML$table$outcome_level)), c("low", "medium"))
})

test_that("print.ThreeStepDistal prints coefficient tables", {
  skip_on_cran()
  set.seed(2028)

  Q <- matrix(c(
    1, 0,
    1, 0,
    0, 1,
    0, 1,
    1, 1,
    1, 1
  ), byrow = TRUE, ncol = 2)
  gs <- data.frame(guess = rep(0.18, nrow(Q)), slip = rep(0.18, nrow(Q)))
  sim <- simGDINA(250, Q, gs.parm = gs, model = "DINA")
  fit <- GDINA(sim$dat, sim$Q, model = "DINA", verbose = 0)
  alpha <- extract(sim, "attribute")
  y_binary <- rbinom(250, 1, plogis(-0.3 + 0.8 * alpha[, 1] - 0.5 * alpha[, 2]))

  out <- ThreeStepDistal(fit, y_binary, level = "attribute", attribute = 1:2, method = "ML")
  printed <- paste(capture.output(print(out)), collapse = "\n")

  expect_match(printed, "Coefficient table")
  expect_match(printed, "Attribute 1")
  expect_match(printed, "std.error")
  expect_match(printed, "odds.ratio")
})