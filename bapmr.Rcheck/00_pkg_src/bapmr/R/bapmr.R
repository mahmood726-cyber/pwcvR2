#' Bootstrap-Aggregated Penalized Meta-Regression (BAP-MR)
#'
#' Implements a bagged version of LASSO meta-regression to overcome the
#' instability of cross-validation in small sample sizes (the regularization
#' paradox). Coefficients are averaged across B bootstrap resamples to produce
#' stable shrinkage estimates.
#'
#' @param y Vector of effect sizes.
#' @param v Vector of sampling variances.
#' @param X Matrix of moderators (k rows, p columns).
#' @param B Number of bootstrap iterations (default = 100).
#' @param alpha Elastic net mixing parameter (1 = LASSO, 0 = Ridge, default = 1).
#' @param seed Optional integer seed for reproducibility. If NULL, no seed is set.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{coefficients} Named vector of bagged coefficients (intercept + moderators).
#'   \item \code{R2_het} Bagged R-squared for heterogeneity (percentage).
#'   \item \code{tau2_res} Mean residual tau-squared across bootstrap iterations.
#'   \item \code{B_successful} Number of bootstrap iterations that converged.
#'   \item \code{raw_boot_coefs} Matrix of bootstrap coefficient estimates.
#' }
#'
#' @examples
#' \donttest{
#' library(metafor)
#' data(dat.bcg, package = "metadat")
#' dat <- escalc(measure = "RR", ai = tpos, bi = tneg,
#'               ci = cpos, di = cneg, data = dat.bcg)
#' result <- bapmr(y = dat$yi, v = dat$vi,
#'                 X = cbind(dat$ablat, dat$year),
#'                 B = 50, seed = 42)
#' print(result$coefficients)
#' cat("BAP-MR R2_het:", result$R2_het, "%\n")
#' }
#'
#' @importFrom stats var coef
#' @export
bapmr <- function(y, v, X, B = 100, alpha = 1, seed = NULL) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for this function.")
  }
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required for this function.")
  }
  if (!is.null(seed)) set.seed(seed)

  k <- length(y)
  p <- ncol(X)
  
  # Storage for bootstrap coefficients and R2
  boot_coefs <- matrix(NA, nrow = B, ncol = p + 1) # Intercept + p moderators
  boot_tau2 <- numeric(B)
  boot_r2 <- numeric(B)
  
  # Fit base model (intercept only) to get baseline tau2
  base_mod <- tryCatch(metafor::rma(yi = y, vi = v, method = "REML"), error = function(e) NULL)
  tau2_base <- if (!is.null(base_mod)) base_mod$tau2 else var(y)
  
  # Calculate weights for glmnet (inverse variance)
  weights <- 1 / v
  
  for (b in 1:B) {
    # Bootstrap sample with replacement
    idx <- sample(1:k, size = k, replace = TRUE)
    
    y_b <- y[idx]
    X_b <- X[idx, , drop = FALSE]
    w_b <- weights[idx]
    v_b <- v[idx]
    
    # Run CV LASSO
    cv_fit <- tryCatch({
      glmnet::cv.glmnet(x = X_b, y = y_b, weights = w_b, alpha = alpha, 
                        nfolds = min(5, k), grouped = FALSE)
    }, error = function(e) NULL)
    
    if (!is.null(cv_fit)) {
      # Extract coefficients using lambda.1se for more conservative shrinkage
      coef_b <- as.vector(coef(cv_fit, s = "lambda.1se"))
      boot_coefs[b, ] <- coef_b
      
      # Calculate predicted values on out-of-bag or bootstrap sample
      preds <- cbind(1, X_b) %*% coef_b
      residuals <- y_b - preds
      
      # Estimate residual tau2 using DerSimonian-Laird type estimator on residuals
      # Simplified tau2 estimator for speed in bootstrap
      Q <- sum(w_b * residuals^2)
      df <- k - sum(coef_b[-1] != 0) - 1
      df <- max(1, df)
      C <- sum(w_b) - sum(w_b^2) / sum(w_b)
      tau2_res <- max(0, (Q - df) / C)
      
      boot_tau2[b] <- tau2_res
      boot_r2[b] <- max(0, 1 - (tau2_res / tau2_base))
    }
  }
  
  # Remove NAs from failed convergence
  valid_b <- !is.na(boot_coefs[, 1])
  boot_coefs <- boot_coefs[valid_b, , drop = FALSE]
  boot_tau2 <- boot_tau2[valid_b]
  boot_r2 <- boot_r2[valid_b]
  
  # Aggregate (Bagging)
  bagged_coefs <- colMeans(boot_coefs)
  names(bagged_coefs) <- c("Intercept", colnames(X))
  
  bagged_r2 <- mean(boot_r2) * 100 # percentage
  
  list(
    coefficients = bagged_coefs,
    R2_het = bagged_r2,
    tau2_res = mean(boot_tau2),
    B_successful = sum(valid_b),
    raw_boot_coefs = boot_coefs
  )
}
