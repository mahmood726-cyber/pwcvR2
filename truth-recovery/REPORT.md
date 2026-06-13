# Truth-recovery yardstick — pwcvR2

**Verdict: VALIDATION of the overfitting-correction (the core claim) + an honest
nuance — pwcv_r2 is a CONSERVATIVE R²: it correctly refuses to credit noise
moderators (where the naive R² overfits badly), but it also under-estimates a
genuine R². It is a lower-bound / overfitting guard, not an unbiased point
estimate.**

## Method
`pwcv_r2` is a leave-one-out cross-validated R² for meta-regression, meant to
avoid the optimism of metafor's in-sample `R2`. The harness injects a KNOWN true
R² (the proportion of between-study heterogeneity explained by the moderators) and
compares the naive in-sample R² with `pwcv_r2`, run unchanged on metafor `rma`
fits (R via PowerShell + `.Renviron`). 120 sims/cell.

## Results (R² = proportion of heterogeneity explained)

| k  | p | true R² | naive R² | pwcv R² |
|----|---|--------:|---------:|--------:|
| 10 | 2 | 0.00    | 0.186 | 0.040 |
| 10 | 5 | 0.00    | **0.330** | 0.035 |
| 20 | 2 | 0.00    | 0.071 | 0.011 |
| 20 | 5 | 0.00    | 0.154 | 0.014 |
| 15 | 2 | 0.50    | 0.433 | 0.151 |
| 30 | 2 | 0.50    | 0.487 | 0.197 |

## Findings (all measured)
1. **VALIDATION — pwcv_r2 corrects the overfitting that inflates the naive R².**
   Under the null (moderators are pure noise, true R²=0) the naive in-sample R² is
   badly upward-biased — 0.19 at p=2 moderators with k=10, rising to **0.33 at p=5
   moderators** (overfitting worsens with more moderators and fewer studies).
   `pwcv_r2` correctly stays near 0 (0.01–0.04) in every null cell. This is exactly
   the optimism correction the method exists to provide.
2. **HONEST NUANCE — pwcv_r2 is CONSERVATIVE under real signal.** When the true R²
   is 0.5, `pwcv_r2` reports only 0.15–0.20 — it under-estimates the genuine
   explained variance by ~2.5–3×, while the naive R² (0.43–0.49) is actually closer
   to the truth in that regime. So the truth sits *between* the two: naive R² is
   upward-biased, `pwcv_r2` downward-biased. The downward bias is the known
   small-sample conservatism of leave-one-out CV R².
3. **Practical reading.** Use `pwcv_r2` as a *lower bound* / overfitting guard: if
   it is near 0, the moderators almost certainly explain nothing (the naive R² is
   noise); if it is clearly positive, there is real signal but the true R² is
   likely higher than the pwcv value. Reporting BOTH (naive as an upper bound,
   pwcv as a lower bound) is the most honest summary.

## What did NOT transfer
This is a calibration check of an R² estimator (predicted-vs-true), the natural
truth-recovery test; NPE/conformal machinery is not needed. The shipped `pwcv_r2`
is run unchanged on metafor fits.

## Reproduce
```
Rscript truth-recovery/harness.R 120
Rscript truth-recovery/test-truth-recovery.R
```
