test_that("threeStep internal corrections reduce to naive fits under identity misclassification", {
  set.seed(1)

  design <- cbind("(Intercept)" = 1, z = rnorm(300))
  prob <- plogis(design %*% c(-0.4, 0.9))
  observed_binary <- rbinom(300, 1, prob)

  fit_binary_1 <- GDINA:::.three_step_binary_fit_core(design, observed_binary, diag(2))
  fit_binary_2 <- GDINA:::.three_step_binary_fit_core(
    design,
    observed_binary,
    diag(2),
    start = fit_binary_1$coefficients
  )
  expect_lt(max(abs(fit_binary_2$coefficients - fit_binary_1$coefficients)), 1e-6)

  eta <- cbind(design %*% c(0.3, -0.5), design %*% c(-0.2, 0.4), 0)
  prob_class <- exp(eta)
  prob_class <- prob_class / rowSums(prob_class)
  observed_profile <- apply(prob_class, 1, function(p) sample.int(3, 1, prob = p))

  fit_profile_1 <- GDINA:::.three_step_multinomial_fit_core(design, observed_profile, diag(3))
  fit_profile_2 <- GDINA:::.three_step_multinomial_fit_core(
    design,
    observed_profile,
    diag(3),
    start = c(fit_profile_1$coefficients)
  )
  expect_lt(max(abs(c(fit_profile_2$coefficients) - c(fit_profile_1$coefficients))), 1e-5)

  fit_profile_first <- GDINA:::.three_step_multinomial_fit_core(
    design,
    observed_profile,
    diag(3),
    reference = "first"
  )
  expect_equal(fit_profile_first$reference, "Class 1")
  expect_equal(colnames(fit_profile_first$coefficients), c("Class 2", "Class 3"))
  expect_equal(colnames(fit_profile_first$fitted), c("Class 1", "Class 2", "Class 3"))
  expect_equal(as.matrix(fit_profile_first$fitted), as.matrix(fit_profile_1$fitted), tolerance = 1e-5)

  binary_table <- fit_binary_1$table
  expect_equal(binary_table$odds.ratio, exp(binary_table$estimate), tolerance = 1e-10)
  expect_equal(binary_table$or.conf.low, exp(binary_table$conf.low), tolerance = 1e-10)
  expect_equal(binary_table$or.conf.high, exp(binary_table$conf.high), tolerance = 1e-10)

  profile_table <- fit_profile_1$table
  expect_equal(profile_table$odds.ratio, exp(profile_table$estimate), tolerance = 1e-10)
  expect_equal(profile_table$or.conf.low, exp(profile_table$conf.low), tolerance = 1e-10)
  expect_equal(profile_table$or.conf.high, exp(profile_table$conf.high), tolerance = 1e-10)
})

test_that("threeStep corrected attribute regression improves over naive regression in a simple example", {
  skip_on_cran()
  set.seed(314)

  N <- 800
  Z <- data.frame(
    z_cont = rnorm(N),
    z_bin = rbinom(N, 1, 0.45),
    z_cat = factor(sample(c("A", "B", "C"), N, replace = TRUE))
  )

  p1 <- plogis(-0.5 + 1.0 * Z$z_cont + 0.9 * Z$z_bin - 0.7 * (Z$z_cat == "B") +
    0.5 * (Z$z_cat == "C"))
  p2 <- plogis(0.1 - 0.7 * Z$z_cont + 0.6 * Z$z_bin + 0.5 * (Z$z_cat == "B") -
    0.4 * (Z$z_cat == "C"))
  alpha <- cbind(rbinom(N, 1, p1), rbinom(N, 1, p2))

  Q <- matrix(c(
    1, 0,
    1, 0,
    1, 0,
    1, 0,
    0, 1,
    0, 1,
    0, 1,
    0, 1,
    1, 1,
    1, 1
  ), byrow = TRUE, ncol = 2)
  gs <- data.frame(guess = rep(0.2, nrow(Q)), slip = rep(0.2, nrow(Q)))

  sim <- simGDINA(N, Q, gs.parm = gs, model = "DINA", attribute = alpha)
  fit <- GDINA(sim$dat, sim$Q, model = "DINA")

  result <- ThreeStepCov(fit, ~ z_cont + z_bin + z_cat, data = Z, level = "attribute", attribute = 1)
  expect_s3_class(result, "ThreeStepCov")
  truth <- glm(alpha[, 1] ~ z_cont + z_bin + z_cat, data = Z, family = binomial())

  naive_error <- mean(abs(result$results[[1]]$naive$coefficients - coef(truth)))
  corrected_error <- mean(abs(result$results[[1]]$corrected$coefficients - coef(truth)))

  expect_lt(corrected_error, naive_error)
})

test_that("threeStep runs profile-level regression with one observed class per respondent", {
  skip_on_cran()
  set.seed(2718)

  N <- 400
  Z <- data.frame(
    z_cont = rnorm(N),
    z_bin = rbinom(N, 1, 0.4),
    z_cat = factor(sample(c("A", "B", "C"), N, replace = TRUE))
  )

  p1 <- plogis(-0.2 + 0.8 * Z$z_cont + 0.7 * Z$z_bin - 0.4 * (Z$z_cat == "B") +
    0.3 * (Z$z_cat == "C"))
  p2 <- plogis(0.1 - 0.6 * Z$z_cont + 0.5 * Z$z_bin + 0.5 * (Z$z_cat == "B") -
    0.2 * (Z$z_cat == "C"))
  alpha <- cbind(rbinom(N, 1, p1), rbinom(N, 1, p2))

  Q <- matrix(c(
    1, 0,
    1, 0,
    1, 0,
    1, 0,
    0, 1,
    0, 1,
    0, 1,
    0, 1,
    1, 1,
    1, 1
  ), byrow = TRUE, ncol = 2)
  gs <- data.frame(guess = rep(0.2, nrow(Q)), slip = rep(0.2, nrow(Q)))

  sim <- simGDINA(N, Q, gs.parm = gs, model = "DINA", attribute = alpha)
  fit <- GDINA(sim$dat, sim$Q, model = "DINA")

  result <- ThreeStepCov(fit, ~ z_cont + z_bin + z_cat, data = Z, level = "profile")
  result_first <- ThreeStepCov(
    fit,
    ~ z_cont + z_bin + z_cat,
    data = Z,
    level = "profile",
    reference = "first"
  )
  class_labels <- apply(extract(fit, "attributepattern"), 1, paste, collapse = "")

  expect_length(result$results$observed, N)
  expect_equal(nrow(result$results$corrected$fitted), N)
  expect_equal(sort(unique(result$results$observed)), seq_len(nrow(extract(fit, "attributepattern"))))
  expect_equal(result$results$corrected$reference, tail(class_labels, 1))
  expect_equal(result_first$results$corrected$reference, class_labels[1])
  expect_equal(colnames(result_first$results$corrected$coefficients), class_labels[-1])
  expect_equal(colnames(result_first$results$corrected$fitted), class_labels)
})

test_that("print.ThreeStepCov prints corrected coefficient tables", {
  skip_on_cran()
  set.seed(1618)

  N <- 250
  Z <- data.frame(z = rnorm(N))
  p1 <- plogis(-0.2 + 0.7 * Z$z)
  p2 <- plogis(0.1 - 0.5 * Z$z)
  alpha <- cbind(rbinom(N, 1, p1), rbinom(N, 1, p2))

  Q <- matrix(c(
    1, 0,
    1, 0,
    0, 1,
    0, 1,
    1, 1,
    1, 1
  ), byrow = TRUE, ncol = 2)
  gs <- data.frame(guess = rep(0.18, nrow(Q)), slip = rep(0.18, nrow(Q)))
  sim <- simGDINA(N, Q, gs.parm = gs, model = "DINA", attribute = alpha)
  fit <- GDINA(sim$dat, sim$Q, model = "DINA", verbose = 0)

  out <- ThreeStepCov(fit, ~ z, data = Z, level = "attribute")
  printed <- paste(capture.output(print(out)), collapse = "\n")

  expect_match(printed, "Corrected coefficients")
  expect_match(printed, "Attribute 1")
  expect_match(printed, "std.error")
  expect_match(printed, "odds.ratio")
})