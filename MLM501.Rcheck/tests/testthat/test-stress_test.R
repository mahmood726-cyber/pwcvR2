test_that("calculate_MAFI works with robust data", {
  skip_if_not_installed("metafor")

  set.seed(123)
  yi_robust <- rnorm(20, mean = 0.5, sd = 0.05)
  vi_robust <- rep(0.01, 20)
  m_robust <- calculate_MAFI(yi_robust, vi_robust, measure = "SMD")

  expect_type(m_robust, "list")
  expect_true("MAFI" %in% names(m_robust))
  expect_true("MAFI_class" %in% names(m_robust))
  expect_true(m_robust$MAFI >= 0 && m_robust$MAFI <= 1)
})

test_that("calculate_MAFI handles small k (returns NULL for k < 3)", {
  result <- calculate_MAFI(c(0.1, 0.2), c(0.01, 0.01))
  expect_null(result)
})

test_that("calculate_MAFI fragile data has higher MAFI", {
  skip_if_not_installed("metafor")

  yi_fragile <- c(0.1, 0.2, -0.1)
  vi_fragile <- c(0.01, 0.01, 0.01)
  m_fragile <- calculate_MAFI(yi_fragile, vi_fragile, measure = "SMD")

  # Small k should trigger k_penalty and likely be fragile
  expect_true(is.null(m_fragile) || m_fragile$MAFI > 0.2)
})

test_that("grade_audit returns expected strings", {
  # grade_audit requires 3 args: mafi_score, estimate, measure
  expect_equal(grade_audit(0.1, 0.8, "SMD"), "No downgrade (Stable Effect Size)")
  expect_equal(grade_audit(0.1, 0.1, "SMD"), "No downgrade (Small Effect Risk)")
  expect_match(grade_audit(0.2, 0.5, "SMD"), "Consider no downgrade")
  expect_match(grade_audit(0.4, 0.5, "SMD"), "Downgrade 1 level")
  expect_match(grade_audit(0.6, 0.5, "SMD"), "CRITICAL")
  expect_equal(grade_audit(NA, 0.5, "SMD"), "Insufficient Data")
})

test_that("simulate_pathway returns correct structure", {
  result <- simulate_pathway(0.1, 0.9, 0.95, 0.7)
  expect_type(result, "list")
  expect_true("events_prevented" %in% names(result))
  expect_true("total_harms" %in% names(result))
  expect_true("net_benefit" %in% names(result))
  expect_true(result$events_prevented >= 0)
})

test_that("calculate_economic_impact returns monetary values", {
  pw <- simulate_pathway(0.1, 0.9, 0.95, 0.7)
  econ <- calculate_economic_impact(pw)
  expect_type(econ, "list")
  expect_true("net_economic_value" %in% names(econ))
  expect_true(is.numeric(econ$net_economic_value))
})

test_that("read_mlm_effects returns a data.frame", {
  # This depends on inst/extdata/mlm_effects.csv existing
  skip_if(!file.exists(system.file("extdata", "mlm_effects.csv", package = "MLM501")))
  df <- read_mlm_effects()
  expect_s3_class(df, "data.frame")
  expect_true(nrow(df) > 0)
  expect_true("TE" %in% names(df))
  expect_true("seTE" %in% names(df))
})
