#' Calculate MAFI (Meta-Analysis Fragility Index)
#' @param yi Vector of effect sizes.
#' @param vi Vector of sampling variances.
#' @param measure Effect size measure (e.g., "SMD", "logOR").
#' @param clinical_threshold Minimum clinically important difference. Auto-set if NULL.
#' @param alpha Significance level (default 0.05).
#' @return A list with MAFI score, class, k, I2, estimate, pval, or NULL if k < 3.
#' @export
calculate_MAFI <- function(yi, vi, measure = "SMD", clinical_threshold = NULL, alpha = 0.05) {
  if (!requireNamespace("metafor", quietly = TRUE)) stop("metafor required")
  k <- length(yi); if (k < 3) return(NULL)
  if (is.null(clinical_threshold)) clinical_threshold <- if (measure == "logOR") log(1.25) else 0.2
  fit <- tryCatch({ metafor::rma(yi = yi, vi = vi, method = "REML") }, error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  estimate <- as.numeric(fit$beta[1]); se <- fit$se; pval <- fit$pval
  if (is.na(estimate) || is.na(se)) return(NULL)
  significant <- pval < alpha; effect_pos <- estimate > 0
  loo <- tryCatch({ metafor::leave1out(fit) }, error = function(e) NULL)
  if (is.null(loo)) return(NULL)
  DFI_rate <- sum((loo$estimate > 0) != effect_pos, na.rm = TRUE) / k
  if (significant) { SFI_rate <- sum(loo$pval >= alpha, na.rm = TRUE) / k } else { SFI_rate <- sum(loo$pval < alpha, na.rm = TRUE) / k }
  CFI_rate <- if (!is.na(clinical_threshold)) { sum((abs(loo$estimate) >= clinical_threshold) != (abs(estimate) >= clinical_threshold), na.rm = TRUE) / k } else 0
  ESI <- min(max(if(abs(estimate) > 0.001) abs(loo$estimate - estimate)/abs(estimate) else abs(loo$estimate - estimate), na.rm=T), 1)
  null_val <- 0; orig_exclude <- (fit$ci.lb > null_val) | (fit$ci.ub < null_val); if (is.na(orig_exclude)) orig_exclude <- FALSE
  CISI_rate <- sum(((loo$ci.lb > null_val) | (loo$ci.ub < null_val)) != orig_exclude, na.rm=T) / k
  het_pen <- (fit$I2 / 100) * 0.2; k_pen <- max(0, (1 - k/20) * 0.3)
  MAFI_core <- (0.30 * DFI_rate + 0.25 * SFI_rate + 0.20 * CFI_rate + 0.15 * ESI + 0.10 * CISI_rate)
  MAFI <- min(1, MAFI_core + het_pen + k_pen)
  return(list(MAFI = round(MAFI, 3), MAFI_class = as.character(cut(MAFI, breaks=c(0,0.15,0.30,0.50,1), labels=c("Robust","Low Fragility","Moderate Fragility","High Fragility"), include.lowest=T)), k = k, I2 = round(fit$I2, 1), estimate = round(estimate, 4), pval = round(pval, 4)))
}

#' Optimized Global Robustness Audit
#' @param df Data frame from read_mlm_effects(). If missing, loaded automatically.
#' @return A data frame with MAFI scores for each analysis.
#' @export
calculate_global_mafi <- function(df) {
  if (missing(df)) df <- read_mlm_effects()
  a_list <- split(df, df$analysis_id)
  message(sprintf("Auditing %d analyses...", length(a_list)))
  results <- lapply(seq_along(a_list), function(i) {
    if (i %% 100 == 0) message(sprintf("  Progress: %d/%d (%.1f%%)", i, length(a_list), (i/length(a_list))*100))
    sub <- a_list[[i]]
    if (nrow(sub) < 3) return(NULL)
    m <- calculate_MAFI(sub$TE, sub$seTE^2, measure = sub$measure[1])
    if (is.null(m)) return(NULL)
    data.frame(review_id=sub$review_id[1], analysis_id=names(a_list)[i], measure=sub$measure[1], k=m$k, MAFI=m$MAFI, MAFI_class=m$MAFI_class, I2=m$I2, estimate=m$estimate, pval=m$pval, stringsAsFactors=F)
  })
  do.call(rbind, results)
}

#' GRADE Imprecision Audit (Enhanced for Reviewer 2)
#' @param mafi_score Numeric MAFI score (0-1).
#' @param estimate Point estimate of the pooled effect.
#' @param measure Effect size measure (currently unused, reserved for future MCID scaling).
#' @return Character string with GRADE downgrade recommendation.
#' @export
grade_audit <- function(mafi_score, estimate, measure) {
  if (is.na(mafi_score)) return("Insufficient Data")
  mcid_buffer <- if(abs(estimate) > 0.5) " (Stable Effect Size)" else " (Small Effect Risk)"
  if (mafi_score < 0.15) return(paste0("No downgrade", mcid_buffer))
  if (mafi_score < 0.30) return(paste0("Consider no downgrade", mcid_buffer))
  if (mafi_score < 0.50) return(paste0("Downgrade 1 level", mcid_buffer))
  return("CRITICAL: Downgrade 1-2 levels (Extreme Fragility)")
}
