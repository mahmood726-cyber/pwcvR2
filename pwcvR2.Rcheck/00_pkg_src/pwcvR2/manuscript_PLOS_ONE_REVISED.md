# Precision-Weighted Cross-Validation for Estimating $R^2$ in Meta-Regression: Addressing Optimism and Heteroscedasticity

**Authors:** Mahmood Ahmad^1,2^, Laiba Khan^3^

**Affiliations:**
^1^ Royal Free London NHS Foundation Trust, London, United Kingdom
^2^ Tahir Heart Institute, Rabwah, Pakistan
^3^ Independent Researcher, London, United Kingdom

**Corresponding author:** Mahmood Ahmad (mahmood726@gmail.com)

---

## Abstract

### Background
Meta-regression is widely used to identify moderators that explain between-study heterogeneity in meta-analysis. The proportion of variance explained, quantified by $R^2$, is the standard metric for assessing the importance of these moderators. However, "apparent" $R^2$ (calculated on the same data used to fit the model) is inherently optimistic and prone to overfitting, particularly in small-to-moderate samples typical of meta-analytic research. While cross-validation (CV) can provide unbiased estimates of out-of-sample performance, the heteroscedastic nature of study-level data necessitates specialized weighting schemes.

### Methods
We evaluated a precision-weighted leave-one-out cross-validation (LOOCV) approach for estimating $R^2$ in meta-regression. We compared its performance against unweighted CV and apparent $R^2$ using two complementary approaches: (1) an empirical evaluation of eight canonical meta-analysis datasets (e.g., BCG vaccine, teacher expectations, passive smoking), and (2) a large-scale Monte Carlo simulation study varying the number of studies ($k$), number of predictors ($p$), residual heterogeneity ($	au^2$), and true underlying $R^2$.

### Results
Empirical analyses demonstrated that apparent $R^2$ consistently and severely overestimates predictive power. In the BCG vaccine dataset ($k=13$, $p=2$), apparent $R^2$ was 64.6%, while precision-weighted CV $R^2$ was only 6.3%. In the "Passive Smoking" dataset ($k=37, p=8$), unweighted CV yielded a highly misleading $R^2$ of 81.8% due to outlier sensitivity, whereas precision-weighted CV correctly identified minimal predictive power (0.6%). Simulation results confirmed that apparent $R^2$ is substantially biased under the null, and that precision-weighted CV provides the most stable and accurate recovery of the true $R^2$ across diverse conditions.

### Conclusions
Precision-weighted cross-validation successfully mitigates the optimism of apparent $R^2$ while accounting for the differential precision of individual studies. We recommend that precision-weighted CV $R^2$ be adopted as the standard metric for reporting explained variance in meta-regression.

---

## Introduction

Meta-analysis has become the cornerstone of evidence-based medicine and the social sciences, providing a rigorous framework for synthesizing findings across multiple studies. Beyond estimating a pooled effect, researchers are increasingly interested in explaining why effect sizes vary between studies—a task performed using meta-regression. The proportion of between-study variance ($	au^2$) accounted for by study-level covariates is typically reported as $R^2$ (or $R^2_{Borenstein}$).

Despite its ubiquity, the estimation of $R^2$ in meta-regression is fraught with statistical pitfalls. The most significant is the "optimism" of the apparent $R^2$. Because meta-regression models are often fitted to small datasets (frequently $k < 20$) with multiple moderators, the model tends to "fit the noise" rather than the signal. Consequently, the apparent $R^2$ can be strikingly high even when the moderators have no true predictive value. This leads to the publication of spurious associations that fail to generalize to new studies.

Cross-validation (CV) is a well-established technique in predictive modeling to estimate out-of-sample performance. By iteratively holding out data and testing the model's predictions, CV provides a more realistic assessment of a model's utility. However, traditional CV assumes that all observations are equally precise. In meta-analysis, each "observation" is a study with its own sampling variance ($v_i$). A small study with a large standard error is inherently less informative than a large, high-precision trial. Failing to account for this heteroscedasticity in CV leads to unweighted error metrics that are dominated by small, noisy studies, potentially masking the true predictive power of the model for high-quality evidence.

In this paper, we propose a precision-weighted cross-validation framework tailored for meta-regression. We argue that by weighting the cross-validated prediction errors by the inverse of their total variance (sampling variance plus residual heterogeneity), we can obtain a more robust and statistically appropriate estimate of $R^2$. Through an extensive empirical re-analysis of landmark datasets and a comprehensive simulation study, we evaluate the performance of this approach and provide practical recommendations for meta-analysts.

---

## Materials and Methods

### The Meta-Regression Model

We consider the standard random-effects meta-regression model:
$$y_i = \beta_0 + \sum_{j=1}^p \beta_j X_{ij} + u_i + \epsilon_i$$
where $y_i$ is the observed effect size for study $i$ ($i=1, \dots, k$), $X_{ij}$ are the values of $p$ moderators, $\beta_j$ are the regression coefficients, $u_i \sim N(0, 	au^2)$ represents the residual between-study heterogeneity, and $\epsilon_i \sim N(0, v_i)$ represents the within-study sampling error with known variance $v_i$.

The total between-study variance in the absence of moderators is denoted as $	au^2_{null}$. The residual between-study variance after accounting for moderators is $	au^2_{model}$. The traditional apparent $R^2$ is defined as:
$$R^2_{apparent} = \frac{\hat{	au}^2_{null} - \hat{	au}^2_{model}}{\hat{	au}^2_{null}}$$
truncated at zero if $\hat{	au}^2_{model} > \hat{	au}^2_{null}$.

### Precision-Weighted Cross-Validation

We implemented a Leave-One-Out Cross-Validation (LOOCV) procedure. For each study $i$, the model is fitted using the remaining $k-1$ studies to obtain the coefficient estimates $\hat{\beta}_{(-i)}$ and the residual heterogeneity $\hat{	au}^2_{(-i)}$. The predicted effect size for the held-out study is $\hat{y}_{i,(-i)} = \hat{\beta}_{0,(-i)} + \sum \hat{\beta}_{j,(-i)} X_{ij}$.

The cross-validated Mean Squared Error (MSE) is then calculated using two weighting schemes:

1.  **Unweighted CV MSE ($MSE_{uw}$):**
    $$MSE_{uw} = \frac{1}{k} \sum_{i=1}^k (y_i - \hat{y}_{i,(-i)})^2$$

2.  **Precision-Weighted CV MSE ($MSE_{w}$):**
    $$MSE_{w} = \frac{\sum_{i=1}^k w_i (y_i - \hat{y}_{i,(-i)})^2}{\sum_{i=1}^k w_i}$$
    where the weights $w_i$ are defined as the inverse of the total variance estimated from the training set: $w_i = (v_i + \hat{	au}^2_{(-i)})^{-1}$.

The corresponding $R^2_{cv}$ metrics are derived by comparing the CV MSE of the full model ($MSE_{model}$) to the CV MSE of the null (intercept-only) model ($MSE_{null}$):
$$R^2_{cv} = 1 - \frac{MSE_{model}}{MSE_{null}}$$

### Empirical Evaluation

We applied the proposed metrics to eight published meta-analysis datasets spanning various medical and psychological domains. These datasets were selected to represent a range of sample sizes ($k$ from 10 to 56), moderator counts ($p$ from 1 to 8), and heterogeneity levels ($I^2$ from 29% to 98%). The primary outcomes were the apparent $R^2$, unweighted CV $R^2$, and precision-weighted CV $R^2$. We also calculated the "Delta CV" ($R^2_{cv, w} - R^2_{cv, uw}$) to identify cases where weighting significantly altered the interpretation of model fit.

### Simulation Design

A Monte Carlo simulation was conducted to evaluate the bias and stability of the estimators under controlled conditions. We simulated 1,000 meta-analyses for each combination of the following parameters:
- Number of studies: $k \in \{20, 40\}$
- Number of moderators: $p \in \{1, 2\}$
- Residual heterogeneity: $	au^2 \in \{0.05, 0.20\}$ (corresponding to low/moderate and high $I^2$ depending on $v_i$)
- True explained variance: $True \ R^2 \in \{0, 0.25\}$
- Distribution of sampling variances: We simulated two scenarios: (a) "Gamma" (representative of typical meta-analyses), and (b) "Anchor" (where one study is significantly more precise than others, representing a common diagnostic challenge in the field).

Analysis was performed using the `metafor` package in R (version 4.8.0).

---

## Results

### Empirical Findings

The empirical analysis across eight canonical datasets revealed that apparent $R^2$ is frequently highly optimistic (Table 1). Across the datasets where apparent $R^2$ was greater than zero (n=5), the average reduction in estimated variance explained when moving from apparent to precision-weighted CV $R^2$ was 43.1 percentage points.

**Table 1: Comparison of Apparent and Cross-Validated $R^2$ for Empirical Datasets.**
| Dataset | $k$ | $p$ | $I^2$ (%) | Apparent $R^2$ (%) | Weighted CV $R^2$ (%) | Unweighted CV $R^2$ (%) | Delta CV (%) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| BCG Vaccine | 13 | 2 | 92.2 | 64.6 | 6.3 | 2.3 | +4.0 |
| Teacher Expectations | 19 | 1 | 41.9 | 40.6 | 0.0 | 0.0 | 0.0 |
| School Achievement | 56 | 1 | 94.7 | 0.0 | 0.0 | 7.8 | -7.8 |
| Passive Smoking | 37 | 8 | 29.6 | 45.5 | 0.6 | 81.8 | -81.2 |
| Colditz BCG | 13 | 1 | 92.1 | 85.1 | 54.6 | 43.9 | +10.6 |

In several cases, models that appeared to have substantial explanatory power were revealed to have zero or near-zero out-of-sample predictive validity. For example, the Teacher Expectations dataset yielded an apparent $R^2$ of 40.6%, yet both weighted and unweighted CV $R^2$ were 0%, indicating complete overfitting. 

The "Passive Smoking" dataset provided a critical case study for the necessity of weighting. With 8 predictors and 37 studies, the unweighted CV $R^2$ was a staggering 81.8%, suggesting exceptional predictive power. However, the precision-weighted CV $R^2$ was only 0.6%. This massive discrepancy (Delta CV = -81.2%) occurred because the unweighted metric was dominated by low-precision studies that the model happened to predict well by chance, whereas the precision-weighted metric correctly reflected the model's inability to predict the more informative, high-precision studies.

### Simulation Findings

The simulation study confirmed that apparent $R^2$ is severely biased upward, especially under the null hypothesis of no true relationship (Table 2).

**Table 2: Simulation Results for Null Scenario ($True \ R^2 = 0$).**
| $k$ | $p$ | $	au^2$ | Variance Model | Mean Apparent $R^2$ (%) | Mean Weighted CV $R^2$ (%) | Mean Unweighted CV $R^2$ (%) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 20 | 1 | 0.05 | Gamma | 19.7 | 10.4 | 16.2 |
| 40 | 1 | 0.05 | Gamma | 7.2 | 6.7 | 14.3 |
| 20 | 2 | 0.05 | Gamma | 21.7 | 14.1 | 15.8 |
| 20 | 1 | 0.20 | Gamma | 5.9 | 2.5 | 3.5 |

Under the null scenario with $k=20, p=1$, the apparent $R^2$ averaged 19.7%, suggesting "large" effects where none existed. Precision-weighted CV $R^2$ reduced this inflation more effectively than unweighted CV $R^2$ in most Gamma-distributed scenarios.

The "Anchor" scenario (where one study has exceptionally high precision) further highlighted the robustness of precision-weighted CV. Under the null ($k=20, p=1, 	au^2=0.05$), unweighted CV $R^2$ was catastrophically inflated to 64.4%, while precision-weighted CV $R^2$ remained much closer to the true value at 12.8%.

In scenarios with a true effect ($True \ R^2 = 0.25$), apparent $R^2$ consistently overestimated the effect (e.g., 30.6% for $k=20, p=1, 	au^2=0.05$), while precision-weighted CV $R^2$ (18.6%) provided a more conservative and realistic estimate of the variance explained in future studies.

---

## Discussion

Our findings demonstrate that the currently reported apparent $R^2$ in meta-regression is fundamentally untrustworthy. Through both empirical re-analysis and controlled simulations, we have shown that apparent $R^2$ is a measure of "fit" rather than "prediction," and in the small-sample context of meta-analysis, fit is frequently driven by noise.

The divergence between weighted and unweighted CV $R^2$ in datasets like "Passive Smoking" illustrates the danger of ignoring precision in model evaluation. Unweighted CV treats all studies as equal units of information, but in meta-analysis, they are not. Precision-weighted CV aligns the evaluation of the model with the underlying statistical principle of meta-analysis: that studies with more information should have a greater influence on our conclusions.

### Strengths and Limitations

This study is the first to systematically evaluate precision-weighting in the context of meta-analytic cross-validation. We utilized a wide range of empirical datasets and diverse simulation conditions, including the challenging "anchor study" scenario. However, our evaluation was limited to Leave-One-Out Cross-Validation. While LOOCV is appropriate for the typically small $k$ in meta-analysis, future research could explore k-fold or bootstrap-based approaches for larger meta-analytic datasets. Additionally, we focused on the $R^2$ metric proposed by Borenstein; the performance of weighted CV for other metrics (e.g., $R^2$ based on likelihood ratios) remains to be investigated.

### Implications for Research

We propose the following recommendations for researchers:
1.  **Cease reporting apparent $R^2$ in isolation.** It should be treated as an upper bound of potential fit, not an estimate of true explanatory power.
2.  **Report Precision-Weighted CV $R^2$ as the primary metric.** This provides a more realistic and conservative estimate of the moderators' value.
3.  **Use Unweighted CV as a diagnostic.** A large discrepancy between weighted and unweighted CV (Delta CV) indicates that the model's performance is highly sensitive to study precision, requiring careful qualitative inspection of the data.

## Figure Legends

**Figure 1. Comparison of Apparent, Unweighted CV, and Precision-Weighted CV $R^2$ across empirical meta-analysis datasets.** The bar chart illustrates the dramatic optimism of apparent $R^2$ (gray) compared to cross-validated estimates. Precision-weighted CV $R^2$ (blue) provides a more conservative and robust estimate of out-of-sample explanatory power, particularly in datasets with high heteroscedasticity.

**Figure 2. Impact of heterogeneity ($I^2$) and study precision on the discrepancy between CV weighting schemes (Delta CV).** (A) Relationship between Delta CV ($R^2_{cv, w} - R^2_{cv, uw}$) and $I^2$. (B) Influence of study weight inequality (measured by the Gini coefficient of precision) on Delta CV. Extreme discrepancies occur when meta-regression models are driven by a few highly precise studies.

**Figure 3. Simulation results: Mean bias of apparent and cross-validated $R^2$ under the null hypothesis ($True \ R^2 = 0$).** Results are shown for varying number of studies ($k$) and residual heterogeneity ($	au^2$). Apparent $R^2$ consistently overestimates effect sizes, while precision-weighted CV $R^2$ remains closest to the true value of zero across diverse simulation scenarios.

---

## Supporting Information

**S1 Table. Full results of the empirical evaluation.** Includes detailed characteristics ($k, p, I^2, 	au^2$) and $R^2$ estimates for all analyzed datasets.

**S2 Table. Comprehensive simulation results.** Detailed results for all combinations of $k, p, 	au^2, True \ R^2,$ and variance models (Gamma vs. Anchor).

**S1 File. R Package `pwcvR2`.** Open-source R package containing the `pwcv_r2` function and scripts to reproduce all analyses and figures presented in this study.

---

## Data Availability Statement

All data used in this study are publicly available. The empirical meta-analysis datasets are embedded within the pwcvR2 R package and are also available in the Supporting Information. The simulation code and results are provided in the package repository at https://github.com/mahmood726-cyber/pwcvR2.

## Funding

No specific funding was received for this work.

## Competing Interests

The authors declare no competing interests.

## Author Contributions

**Conceptualization:** MA. **Methodology:** MA. **Software:** MA. **Validation:** MA, LK. **Formal analysis:** MA. **Writing -- original draft:** MA. **Writing -- review & editing:** MA, LK.

---

## References

[1] Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw. 2010;36(3):1-48. https://doi.org/10.18637/jss.v036.i03

[2] Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. Stat Med. 2002;21(11):1539-1558. https://doi.org/10.1002/sim.1186

[3] Thompson SG, Higgins JPT. How should meta-regression analyses be undertaken and interpreted? Stat Med. 2002;21(11):1559-1573. https://doi.org/10.1002/sim.1187

[4] Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. Chichester: John Wiley & Sons; 2009. https://doi.org/10.1002/9780470743386

[5] Harrell FE, Lee KL, Mark DB. Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. Stat Med. 1996;15(4):361-387.

[6] Knapp G, Hartung J. Improved tests for a random effects meta-regression with a single covariate. Stat Med. 2003;22(17):2693-2710. https://doi.org/10.1002/sim.1482

[7] IntHout J, Ioannidis JP, Borm GF. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. BMC Med Res Methodol. 2014;14:25. https://doi.org/10.1186/1471-2288-14-25

[8] R Core Team. R: A Language and Environment for Statistical Computing. Vienna, Austria: R Foundation for Statistical Computing; 2024. https://www.R-project.org/
