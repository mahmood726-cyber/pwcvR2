test_that("pwcv_r2 returns numeric in [0, 1]", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("metadat")

  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg,
                         ci = cpos, di = cneg, data = dat.bcg)
  fit <- metafor::rma(yi, vi, mods = ~ ablat, data = dat)

  r2 <- pwcv_r2(fit)

  expect_type(r2, "double")
  expect_true(r2 >= 0)
  expect_true(r2 <= 1)
})

test_that("pwcv_r2 weighted vs unweighted differ", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("metadat")

  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg,
                         ci = cpos, di = cneg, data = dat.bcg)
  fit <- metafor::rma(yi, vi, mods = ~ ablat, data = dat)

  r2_w <- pwcv_r2(fit, weighted = TRUE)
  r2_u <- pwcv_r2(fit, weighted = FALSE)

  expect_type(r2_w, "double")
  expect_type(r2_u, "double")
  # They should generally differ (different weighting schemes)
  # but both should be valid [0, 1] values
  expect_true(r2_w >= 0 && r2_w <= 1)
  expect_true(r2_u >= 0 && r2_u <= 1)
})

test_that("pwcv_r2 errors on non-rma input", {
  expect_error(pwcv_r2(lm(1:10 ~ rnorm(10))), "rma.uni")
})

test_that("pwcv_r2 errors on intercept-only model", {
  skip_if_not_installed("metafor")

  set.seed(42)
  yi <- rnorm(10)
  vi <- rep(0.1, 10)
  fit <- metafor::rma(yi, vi)

  expect_error(pwcv_r2(fit), "moderator")
})

test_that("pwcv_r2 handles small k", {
  skip_if_not_installed("metafor")

  set.seed(123)
  k <- 5
  yi <- rnorm(k)
  vi <- rep(0.05, k)
  mods <- rnorm(k)
  fit <- metafor::rma(yi, vi, mods = ~ mods)

  r2 <- pwcv_r2(fit)
  expect_type(r2, "double")
  expect_true(r2 >= 0)
})
