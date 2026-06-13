# ============================================================
# harness.R -- Truth-recovery yardstick for pwcvR2.
#
# pwcv_r2 is a leave-one-out CROSS-VALIDATED R^2 for meta-regression, meant to
# avoid the optimism (overfitting) of the in-sample R^2 that metafor reports
# (model$R2). The honest test: inject a KNOWN true R^2 (proportion of between-
# study heterogeneity explained by the moderators) and check whether
#   (a) under the NULL (noise moderators, true R^2 = 0) the naive R^2 is biased
#       UP while pwcv_r2 stays ~0, and
#   (b) under a real signal both track the true R^2.
#
# Uses the app's own pwcv_r2 (sourced) on metafor rma.uni fits, run unchanged.
# Truth-first: every number is produced from seeded simulation. Run:
#   Rscript truth-recovery/harness.R 300
# ============================================================
suppressMessages(library(metafor))
this <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(this), "..", "R", "pwcv_r2.R"))

args <- commandArgs(TRUE)
NSIM <- if (length(args) >= 1) as.integer(args[1]) else 300

# Build a meta-regression with a KNOWN true R^2.
# theta_i = X_i %*% beta + u_i ; Var(Xbeta) explained, Var(u)=tau2_resid.
gen <- function(k, p, true_r2, tau2_total, seed) {
  set.seed(seed)
  X <- matrix(rnorm(k * p), k, p)
  if (true_r2 == 0) {
    beta <- rep(0, p); tau2_resid <- tau2_total
  } else {
    # one active moderator carries the signal; scale beta to hit true_r2
    beta <- c(1, rep(0, p - 1))
    var_xb <- var(as.vector(X %*% beta))
    # want var_xb_scaled / (var_xb_scaled + tau2_resid) = true_r2, with total = tau2_total
    explained <- true_r2 * tau2_total; tau2_resid <- (1 - true_r2) * tau2_total
    beta <- beta * sqrt(explained / var_xb)
  }
  mu <- 0.2 + as.vector(X %*% beta)
  theta <- mu + rnorm(k, 0, sqrt(tau2_resid))
  vi <- runif(k, 0.02, 0.12)
  yi <- theta + rnorm(k, 0, sqrt(vi))
  list(yi = yi, vi = vi, X = X)
}

run_cell <- function(k, p, true_r2, tau2_total) {
  naive <- c(); cv <- c()
  for (s in 1:NSIM) {
    d <- gen(k, p, true_r2, tau2_total, 4000 + s)
    m <- tryCatch(rma(d$yi, d$vi, mods = d$X, method = "REML"), error = function(e) NULL)
    if (is.null(m)) next
    naive <- c(naive, max(0, m$R2) / 100)
    cvv <- tryCatch(pwcv_r2(m), error = function(e) NA)
    if (is.finite(cvv)) cv <- c(cv, cvv)
  }
  list(naive = round(mean(naive), 3), cv = round(mean(cv), 3),
       naive_sd = round(sd(naive), 3))
}

cat(sprintf("\n# Truth-recovery yardstick -- pwcvR2  nsim=%d\n", NSIM))
cat("R^2 = proportion of between-study heterogeneity explained by moderators\n\n")
cat(sprintf("%4s %4s %9s | %12s | %10s\n", "k", "p", "true R2", "naive R2", "pwcv R2"))
cells <- list(c(10, 2, 0.0), c(10, 5, 0.0), c(20, 2, 0.0), c(20, 5, 0.0),
              c(15, 2, 0.5), c(30, 2, 0.5))
for (cl in cells) {
  r <- run_cell(cl[1], cl[2], cl[3], 0.08)
  cat(sprintf("%4d %4d %9.2f | %12.3f | %10.3f\n", cl[1], cl[2], cl[3], r$naive, r$cv))
}
cat("\n(under true R2=0 the naive R2 is the OVERFITTING bias; pwcv_r2 should be ~0.\n")
cat(" under true R2=0.5 both should be near 0.5.)\n")
