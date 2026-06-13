# Rscript truth-recovery/test-truth-recovery.R   (exit 0 = all pass)
# Light measured invariants for the pwcvR2 yardstick (LOO refits are slow, so
# nsim is small). Seeded. Full grid in harness.R / REPORT.md.
suppressMessages(library(metafor))
this <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(this), "..", "R", "pwcv_r2.R"))

NSIM <- 80
gen <- function(k, p, true_r2, tau2_total, seed) {
  set.seed(seed)
  X <- matrix(rnorm(k * p), k, p)
  if (true_r2 == 0) { beta <- rep(0, p); tau2_resid <- tau2_total } else {
    beta <- c(1, rep(0, p - 1)); var_xb <- var(as.vector(X %*% beta))
    explained <- true_r2 * tau2_total; tau2_resid <- (1 - true_r2) * tau2_total
    beta <- beta * sqrt(explained / var_xb)
  }
  mu <- 0.2 + as.vector(X %*% beta); theta <- mu + rnorm(k, 0, sqrt(tau2_resid))
  vi <- runif(k, 0.02, 0.12); yi <- theta + rnorm(k, 0, sqrt(vi))
  list(yi = yi, vi = vi, X = X)
}
means <- function(k, p, true_r2) {
  nv <- c(); cv <- c()
  for (s in 1:NSIM) {
    d <- gen(k, p, true_r2, 0.08, 4000 + s)
    m <- tryCatch(rma(d$yi, d$vi, mods = d$X, method = "REML"), error = function(e) NULL)
    if (is.null(m)) next
    nv <- c(nv, max(0, m$R2) / 100)
    cc <- tryCatch(pwcv_r2(m), error = function(e) NA); if (is.finite(cc)) cv <- c(cv, cc)
  }
  list(naive = mean(nv), cv = mean(cv))
}

ok <- TRUE
report <- function(name, cond, detail) { cat(sprintf("%-4s %s  %s\n", if (cond) "PASS" else "FAIL", name, detail)); if (!cond) ok <<- FALSE }

n0 <- means(10, 5, 0.0)   # null, many moderators -> naive overfits
report("under the NULL the naive R2 is badly inflated (overfitting)", n0$naive > 0.15, sprintf("(naive %.3f)", n0$naive))
report("pwcv_r2 correctly stays ~0 under the null (no overfitting)", n0$cv < 0.10, sprintf("(pwcv %.3f)", n0$cv))
report("pwcv_r2 is far closer to the truth (0) than naive R2", n0$cv < n0$naive - 0.1, sprintf("(pwcv %.3f vs naive %.3f)", n0$cv, n0$naive))

s5 <- means(30, 2, 0.5)   # real signal
report("under a real signal pwcv_r2 detects it (clearly above its null level)",
       s5$cv > 0.10, sprintf("(pwcv %.3f at true R2=0.5)", s5$cv))
report("HONEST: but pwcv_r2 UNDER-estimates a real R2 (conservative lower bound)",
       s5$cv < 0.40 && s5$cv < s5$naive, sprintf("(pwcv %.3f vs naive %.3f, true 0.50)", s5$cv, s5$naive))

cat(if (ok) "\nAll measured invariants hold.\n" else "\nSOME INVARIANTS FAILED.\n")
quit(status = if (ok) 0 else 1)
