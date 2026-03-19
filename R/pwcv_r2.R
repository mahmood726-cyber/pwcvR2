#' Precision-Weighted Cross-Validation for Meta-Regression R-squared
#'
#' Estimates cross-validated R-squared for heterogeneity in meta-regression
#' using leave-one-out cross-validation with optional precision weighting.
#' Standard (apparent) R-squared in meta-regression is optimistic when the
#' number of studies (k) is small relative to the number of moderators (p).
#' This function provides a less biased estimate by evaluating out-of-sample
#' prediction performance.
#'
#' @param model A \code{rma.uni} object from \code{metafor::rma()} fitted
#'   with at least one moderator.
#' @param weighted Logical. If \code{TRUE} (default), uses precision weights
#'   (1 / (vi + tau2)) for the cross-validated MSE calculation, accounting for
#'   heteroscedasticity. If \code{FALSE}, uses unweighted MSE.
#'
#' @return A numeric value: the cross-validated R-squared for heterogeneity
#'   (proportion of between-study variance explained by the moderators),
#'   truncated at 0 if negative.
#'
#' @details
#' For each study i, the model is refit on the remaining k-1 studies and
#' used to predict the held-out effect size. The cross-validated MSE is
#' compared to the MSE of a null (intercept-only) model fitted via the
#' same LOO procedure. Full-model weights are used for both the full and
#' null MSE to ensure a fair comparison of prediction error reduction.
#'
#' @examples
#' \donttest{
#' library(metafor)
#' data(dat.bcg, package = "metadat")
#' dat <- escalc(measure = "RR", ai = tpos, bi = tneg,
#'               ci = cpos, di = cneg, data = dat.bcg)
#' fit <- rma(yi, vi, mods = ~ ablat, data = dat)
#' pwcv_r2(fit)
#' pwcv_r2(fit, weighted = FALSE)
#' }
#'
#' @export
pwcv_r2 <- function(model, weighted = TRUE) {
  if (!inherits(model, "rma.uni")) {
    stop("Model must be an 'rma.uni' object from metafor.")
  }

  yi <- model$yi
  vi <- model$vi
  X  <- model$X
  k  <- length(yi)

  if (ncol(X) < 2) {
    stop("Model must include at least one moderator (intercept-only models have nothing to cross-validate).")
  }

  # Leave-one-out predictions (full model)
  preds <- numeric(k)
  weights <- numeric(k)

  for (i in 1:k) {
    res_minus_i <- metafor::rma(yi[-i], vi[-i], mods = X[-i, -1, drop = FALSE],
                                method = model$method,
                                control = model$control)
    pred_obj <- predict(res_minus_i, newmods = X[i, -1, drop = FALSE])
    preds[i] <- pred_obj$pred
    weights[i] <- 1 / (vi[i] + res_minus_i$tau2)
  }

  # Leave-one-out predictions (null model)
  null_preds <- numeric(k)
  for (i in 1:k) {
    res_null_minus_i <- metafor::rma(yi[-i], vi[-i], method = model$method,
                                     control = model$control)
    pred_null <- predict(res_null_minus_i)
    null_preds[i] <- pred_null$pred
  }

  # Cross-validated MSE (full-model weights used for both to ensure fair comparison)
  if (weighted) {
    mse_model <- sum(weights * (yi - preds)^2) / sum(weights)
    mse_null  <- sum(weights * (yi - null_preds)^2) / sum(weights)
  } else {
    mse_model <- mean((yi - preds)^2)
    mse_null  <- mean((yi - null_preds)^2)
  }

  r2_cv <- 1 - (mse_model / mse_null)
  max(0, r2_cv)
}
