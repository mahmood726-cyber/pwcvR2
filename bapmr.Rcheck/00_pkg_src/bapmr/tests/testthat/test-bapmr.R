test_that("bapmr returns correct structure", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("metafor")

  set.seed(42)
  k <- 15
  yi <- rnorm(k, 0, sqrt(0.3))
  vi <- runif(k, 0.01, 0.1)
  X <- cbind(rnorm(k), rnorm(k))
  colnames(X) <- c("mod1", "mod2")

  result <- bapmr(yi, vi, X, B = 20, seed = 123)

  expect_type(result, "list")
  expect_true("coefficients" %in% names(result))
  expect_true("R2_het" %in% names(result))
  expect_true("tau2_res" %in% names(result))
  expect_true("B_successful" %in% names(result))
  expect_length(result$coefficients, 3)  # intercept + 2 mods
  expect_true(result$R2_het >= 0 && result$R2_het <= 100)
  expect_true(result$tau2_res >= 0)
  expect_true(result$B_successful > 0 && result$B_successful <= 20)
})

test_that("bapmr is reproducible with seed", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("metafor")

  k <- 10
  yi <- rnorm(k, 0, 0.5)
  vi <- rep(0.05, k)
  X <- matrix(rnorm(k * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))

  r1 <- bapmr(yi, vi, X, B = 10, seed = 99)
  r2 <- bapmr(yi, vi, X, B = 10, seed = 99)

  expect_equal(r1$coefficients, r2$coefficients)
  expect_equal(r1$R2_het, r2$R2_het)
})

test_that("bapmr handles single moderator", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("metafor")

  set.seed(1)
  k <- 12
  yi <- rnorm(k)
  vi <- runif(k, 0.01, 0.1)
  X <- matrix(rnorm(k), ncol = 1, dimnames = list(NULL, "mod1"))

  result <- bapmr(yi, vi, X, B = 15, seed = 42)

  expect_type(result, "list")
  expect_length(result$coefficients, 2)  # intercept + 1 mod
})

test_that("bapmr errors without glmnet", {
  skip_if(requireNamespace("glmnet", quietly = TRUE), "glmnet is installed")
  expect_error(bapmr(1:5, rep(0.1, 5), matrix(1:5, ncol = 1)), "glmnet")
})

test_that("bapmr alpha parameter works for Ridge", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("metafor")

  set.seed(7)
  k <- 15
  yi <- rnorm(k, 0, 0.4)
  vi <- runif(k, 0.02, 0.08)
  X <- matrix(rnorm(k * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))

  lasso <- bapmr(yi, vi, X, B = 15, alpha = 1, seed = 42)
  ridge <- bapmr(yi, vi, X, B = 15, alpha = 0, seed = 42)

  expect_type(lasso, "list")
  expect_type(ridge, "list")
  # Ridge coefficients should generally be non-zero (no exact sparsity)
})
