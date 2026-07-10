
test_that("0b: extract initial.catprob returns initial parameters", {
  Q <- matrix(c(1,0,
                0,1,
                1,1), ncol = 2, byrow = TRUE)

  itempar <- list(c(0.1, 0.9),
                  c(0.2, 0.8),
                  c(0.1, 0.2, 0.3, 0.8))

  set.seed(12345)
  sim <- simGDINA(200, Q, catprob.parm = itempar, model = "GDINA")
  dat <- extract(sim, "dat")

  est <- GDINA(dat, Q, catprob.parm = itempar, model = "GDINA",
               control = list(maxitr = 2))

  ic <- extract(est, "initial.catprob")
  # Should not be NULL (the old code returned NULL due to $technical typo)
  expect_true(!is.null(ic))
  expect_true(is.list(ic))
  expect_equal(length(ic), nrow(Q))
})

test_that("0c: extract itemprob.parm works for sequential models", {
  Qc <- sim20seqGDINA$simQ
  dat <- sim20seqGDINA$simdat

  # Fit a sequential model with zero iterations (just use supplied params)
  est <- GDINA(dat, Qc, sequential = TRUE,
               control = list(maxitr = 1))

  ip <- extract(est, "itemprob.parm")
  # Should not error (the old code referenced out-of-scope Q)
  expect_true(!is.null(ip))
  expect_true(is.list(ip))
  # Number of items (unique item numbers in Qc column 1)
  n_items <- length(unique(Qc[, 1]))
  expect_equal(length(ip), n_items)
})

test_that("0d: summary.dif prints test results without error", {
  set.seed(12345)
  Q <- matrix(c(1,0,
                0,1,
                1,1), ncol = 2, byrow = TRUE)
  gs1 <- data.frame(guess = c(0.2, 0.2, 0.2), slip = c(0.2, 0.2, 0.2))
  gs2 <- data.frame(guess = c(0.2, 0.2, 0.4), slip = c(0.2, 0.2, 0.1))

  sim1 <- simGDINA(300, Q, gs.parm = gs1, model = "DINA")
  sim2 <- simGDINA(300, Q, gs.parm = gs2, model = "DINA")
  dat <- rbind(extract(sim1, "dat"), extract(sim2, "dat"))
  gr <- rep(c("G1", "G2"), each = 300)

  d <- dif(dat, Q, group = gr, method = "wald")
  # summary.dif should not error (old code accessed nonexistent $CDM1/$CDM2)
  s <- expect_no_error(summary(d))
  # Return value should be the test data frame
  expect_true(is.data.frame(s))
  expect_true("p.value" %in% colnames(s))
})

test_that("0d1: dif supports multi-group Wald with observed factor levels", {
  set.seed(54321)
  Q <- matrix(c(1,0,
                0,1,
                1,1), ncol = 2, byrow = TRUE)
  gs1 <- data.frame(guess = c(0.2, 0.2, 0.2), slip = c(0.2, 0.2, 0.2))
  gs2 <- data.frame(guess = c(0.2, 0.25, 0.35), slip = c(0.2, 0.2, 0.15))
  gs3 <- data.frame(guess = c(0.25, 0.2, 0.25), slip = c(0.15, 0.2, 0.2))

  sim1 <- simGDINA(120, Q, gs.parm = gs1, model = "DINA")
  sim2 <- simGDINA(120, Q, gs.parm = gs2, model = "DINA")
  sim3 <- simGDINA(120, Q, gs.parm = gs3, model = "DINA")
  dat <- rbind(extract(sim1, "dat"), extract(sim2, "dat"), extract(sim3, "dat"))
  gr <- factor(rep(c("G1", "G2", "G3"), each = 120),
               levels = c("G1", "G2", "G3", "G4"))

  d <- expect_no_error(dif(dat, Q, group = gr, method = "wald"))

  expect_s3_class(d, "dif")
  expect_equal(nrow(d$test), nrow(Q))
  expect_equal(colnames(d$test), c("Wald stat.", "df", "p.value", "adj.pvalue"))
})

test_that("0d2: dif supports multi-group LR with observed factor levels", {
  set.seed(54321)
  Q <- matrix(c(1,0,
                0,1,
                1,1), ncol = 2, byrow = TRUE)
  gs <- data.frame(guess = c(0.2, 0.2, 0.2), slip = c(0.2, 0.2, 0.2))

  sim1 <- simGDINA(80, Q, gs.parm = gs, model = "DINA")
  sim2 <- simGDINA(80, Q, gs.parm = gs, model = "DINA")
  sim3 <- simGDINA(80, Q, gs.parm = gs, model = "DINA")
  dat <- rbind(extract(sim1, "dat"), extract(sim2, "dat"), extract(sim3, "dat"))
  gr <- factor(rep(c("G1", "G2", "G3"), each = 80),
               levels = c("G1", "G2", "G3", "G4"))

  d <- expect_no_error(dif(dat, Q, group = gr, method = "LR"))

  expect_s3_class(d, "dif")
  expect_equal(nrow(d$test), nrow(Q))
  expect_equal(colnames(d$test), c("neg2LL", "LRstat", "df", "p.value", "adj.pvalue"))
})

test_that("0d3: pairwiseDIF reuses dif object inputs for post hoc tests", {
  set.seed(54321)
  Q <- matrix(c(1,0,
                0,1,
                1,1), ncol = 2, byrow = TRUE)
  gs1 <- data.frame(guess = c(0.2, 0.2, 0.2), slip = c(0.2, 0.2, 0.2))
  gs2 <- data.frame(guess = c(0.2, 0.2, 0.45), slip = c(0.2, 0.2, 0.1))
  gs3 <- data.frame(guess = c(0.25, 0.2, 0.25), slip = c(0.15, 0.2, 0.2))

  sim1 <- simGDINA(120, Q, gs.parm = gs1, model = "DINA")
  sim2 <- simGDINA(120, Q, gs.parm = gs2, model = "DINA")
  sim3 <- simGDINA(120, Q, gs.parm = gs3, model = "DINA")
  dat <- rbind(extract(sim1, "dat"), extract(sim2, "dat"), extract(sim3, "dat"))
  gr <- factor(rep(c("G1", "G2", "G3"), each = 120),
               levels = c("G1", "G2", "G3", "G4"))

  d <- expect_no_error(dif(dat, Q, group = gr, method = "LR"))

  expect_true(is.list(d$input))
  expect_equal(d$input$gr.label, c("G1", "G2", "G3"))

  pw <- expect_no_error(pairwiseDIF(d, alpha.level = 1))

  expect_s3_class(pw, "pairwiseDIF")
  expect_equal(nrow(pw$test), nrow(Q) * choose(3, 2))
  expect_equal(colnames(pw$test), c("item", "group1", "group2", "Wald stat.", "df", "p.value", "adj.pvalue"))
  expect_equal(colnames(pw$posthoc.rearrangement),
               c("refit.item", "original.item", "item.label", "item.type", "group.no", "group.label"))
  expect_equal(nrow(pw$posthoc.rearrangement), nrow(Q) * 3)
  expect_true(all(pw$posthoc.rearrangement$original.item %in% seq_len(nrow(Q))))
  expect_true(all(pw$posthoc.rearrangement$group.no %in% 1:3))
  expect_no_error(print(pw))
})

test_that("0d4: pairwiseDIF supports manual mode and empty object-mode follow-up", {
  set.seed(54321)
  Q <- matrix(c(1,0,
                0,1,
                1,1), ncol = 2, byrow = TRUE)
  gs1 <- data.frame(guess = c(0.2, 0.2, 0.2), slip = c(0.2, 0.2, 0.2))
  gs2 <- data.frame(guess = c(0.2, 0.2, 0.45), slip = c(0.2, 0.2, 0.1))
  gs3 <- data.frame(guess = c(0.25, 0.2, 0.25), slip = c(0.15, 0.2, 0.2))

  sim1 <- simGDINA(120, Q, gs.parm = gs1, model = "DINA")
  sim2 <- simGDINA(120, Q, gs.parm = gs2, model = "DINA")
  sim3 <- simGDINA(120, Q, gs.parm = gs3, model = "DINA")
  dat <- rbind(extract(sim1, "dat"), extract(sim2, "dat"), extract(sim3, "dat"))
  gr <- factor(rep(c("G1", "G2", "G3"), each = 120),
               levels = c("G1", "G2", "G3", "G4"))

  d <- expect_no_error(dif(dat, Q, group = gr, method = "LR"))

  empty_pw <- expect_warning(pairwiseDIF(d, alpha.level = 0),
                             "No items were flagged for DIF for pairwise follow-up")
  expect_s3_class(empty_pw, "pairwiseDIF")
  expect_equal(nrow(empty_pw$test), 0)

  manual_pw <- expect_no_error(pairwiseDIF(dat = dat, Q = Q, group = gr,
                                           model = "DINA", dif.items = 3,
                                           anchor.items = 1:2))

  expect_s3_class(manual_pw, "pairwiseDIF")
  expect_equal(nrow(manual_pw$test), choose(3, 2))
  expect_true(all(manual_pw$test$item == "Item 3"))
  expect_equal(manual_pw$posthoc.rearrangement$original.item, c(1L, 2L, 3L, 3L, 3L))
  expect_equal(manual_pw$posthoc.rearrangement$group.label,
               c("shared", "shared", "G1", "G2", "G3"))
  expect_no_error(plot(manual_pw, item = 3))
})

test_that("0e: modelcomp LM with only reducedMDINO extracts item names correctly", {
  dat <- sim10GDINA$simdat
  Q <- sim10GDINA$simQ

  dino.fit <- GDINA(dat = dat, Q = Q, model = "DINO")

  # LM test with only reducedMDINO supplied (the branch that was buggy)
  # Must restrict models to "DINO" only, since we only provide reducedMDINO
  mc <- expect_no_error(
    modelcomp(method = "LM", models = "DINO",
              LM.args = list(reducedMDINO = dino.fit))
  )

  stats <- mc$stats
  # Should have item names (not error from wrong object)
  expect_true(!is.null(rownames(stats)))
  # Only items with >1 attribute are tested
  multi_attr_items <- which(rowSums(Q) > 1)
  expect_equal(nrow(stats), length(multi_attr_items))
})

test_that("0f: solnp optimizer does not produce spurious print output", {
  Q <- matrix(c(1,0,
                0,1,
                1,1), ncol = 2, byrow = TRUE)

  itempar <- list(c(0.1, 0.9),
                  c(0.2, 0.8),
                  c(0.1, 0.2, 0.3, 0.8))

  set.seed(12345)
  sim <- simGDINA(200, Q, catprob.parm = itempar, model = "GDINA")
  dat <- extract(sim, "dat")

  # Capture all console output during fitting with solnp
  out <- capture.output({
    est <- GDINA(dat, Q, model = "GDINA", solver = "solnp",
                 control = list(maxitr = 3))
  })

  # Filter out expected iteration messages (Iter = ...)
  unexpected <- out[!grepl("^\\s*Iter\\s*=", out) & nchar(trimws(out)) > 0]
  # No spurious numeric output from print(op$pars) should remain
  expect_equal(length(unexpected), 0,
               info = paste("Unexpected output:", paste(unexpected, collapse = "\n")))
})
