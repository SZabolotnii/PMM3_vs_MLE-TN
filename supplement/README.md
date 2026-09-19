# Supplementary material

**Paper:** S. Zabolotnii, A. Chepynoha, P. Klopotovskyi, *PMM3 vs MLE-TN: Comparison
of Semiparametric and Parametric Estimation of Linear Regression Parameters under
Two-piece Normal Distributed Errors*, submitted to ITEST 2026.

**Repository:** https://github.com/SZabolotnii/PMM3_vs_MLE-TN

**DOI:** _to be assigned by Zenodo_

The files are stored results of the study pipeline. Nothing here is hand-edited:
`scripts/09_build_supplement.R` copies them from `results/tables/`, and
`MD5SUMS.txt` lets you check them against the repository.

## Reproduction

```r
# from the repository root
source("scripts/run_all.R")               # T1, T2, T6, T7  (about 1–2 min, M = 500)
source("scripts/09_build_supplement.R")   # rebuilds this folder
```

- **Seed:** 20260312 (`CFG$seed` in `R/config.R`).
- **Package versions:** see `SESSION_INFO.txt`. PMM2 and PMM3 come from the
  `EstemPMM` R package.
- The real-data files (`iris_versicolor_results.csv`, `real_data_*.csv`) are
  stored outputs. `run_all.R` does not regenerate them.

## Error model

TN(λ) is the equal mixture ½N(−λ, 1) + ½N(λ, 1), with density
f(z | λ) = φ(z)·exp(−λ²/2)·cosh(λz). It is symmetric and platykurtic for
every λ > 0, and bimodal for λ > 1. Regression model: y = β₀ + β₁x + ε with
β₀ = 1 and β₁ = 2. All efficiency figures refer to β₁.

- **g₃** = Var(PMM3) / Var(OLS). Lower is better.
- **ARE** = Var(OLS) / Var(method). Higher is better, and ARE > 1 beats OLS.

## Files

### `T1_theoretical.csv` — TN moments and theoretical PMM3 efficiency

| Column | Meaning |
|---|---|
| `lambda` | TN shape parameter λ |
| `mu2`, `mu4`, `mu6` | Central moments of order 2, 4, 6 |
| `gamma4`, `gamma6` | Excess kurtosis and sixth-order cumulant coefficient |
| `g3` | Theoretical g₃ = 1 − γ₄² / (6 + 9γ₄ + γ₆) |
| `improvement_pct` | Variance reduction relative to OLS, (1 − g₃)·100 |
| `is_bimodal`, `mode` | Whether TN(λ) is bimodal (λ > 1) |

### `T2_monte_carlo.csv` — full Monte Carlo table

24 scenarios (λ ∈ {0.5, 1, 1.5, 2, 3, 5} × n ∈ {46, 100, 200, 500}) × 4 methods
(OLS, PMM2, PMM3, MLE = MLE-TN), M = 500 replications each.

| Column | Meaning |
|---|---|
| `lambda`, `n`, `method` | Scenario and estimator |
| `bias`, `variance`, `mse` | Of β̂₁ over the valid replications |
| `are` | Var(OLS) / Var(method) |
| `g3_empirical` | Var(PMM3) / Var(OLS); PMM3 rows only |
| `conv_rate` | Share of replications in which the estimator converged |
| `n_valid` | Number of replications with a finite estimate |

### `T6_misspecification.csv` — robustness when the errors are not TN

n = 100, M = 500. MLE-TN still assumes TN errors.

| Column | Meaning |
|---|---|
| `distribution` | True error law: uniform, triangular, logistic, Student t₁₀ |
| `are_mle`, `are_pmm3` | ARE of MLE-TN and PMM3 |
| `bias_mle`, `bias_pmm3` | Bias of β̂₁ |
| `na_pct_mle`, `na_pct_pmm3` | Percent of failed fits |

### `T7_bimodal_test.csv` — PMM3 Newton–Raphson on bimodal TN

n = 100, M = 200, λ ∈ {0.8, 1.0, 1.2, 1.5, 2.0, 3.0}.

| Column | Meaning |
|---|---|
| `lambda`, `is_bimodal` | Scenario |
| `bias_pmm3` | Bias of the PMM3 estimate of β₁ |
| `conv_rate` | Share of replications that converged |
| `mean_iter` | Mean number of Newton–Raphson iterations |
| `na_pct` | Percent of failed fits |

### `iris_versicolor_results.csv` — Iris versicolor, Sepal.Length ~ Sepal.Width (n = 50)

| Column | Meaning |
|---|---|
| `method` | OLS, PMM3, MLE-TN |
| `beta1` | Slope estimate |
| `loo_mse` | Leave-one-out mean squared prediction error |
| `g3_emp` | Bootstrap Var(method) / Var(OLS), B = 2000 |

### `real_data_all_candidates.csv` — screening of the other real datasets

Datasets: whiteside (MASS), trees, cars, faithful, iris setosa.

| Column | Meaning |
|---|---|
| `dataset`, `n` | Dataset, model formula and sample size |
| `gamma3`, `gamma4` | Skewness and excess kurtosis of the OLS residuals |
| `g3_theor`, `improve_pct` | Theoretical PMM3 efficiency and variance reduction |
| `loo_rel_pmm3`, `loo_rel_mle` | LOO-MSE relative to OLS |
| `boot_g3_pmm3`, `boot_g3_mle` | Bootstrap Var(method) / Var(OLS), B = 2000 |

### `real_data_bootstrap.csv`, `real_data_loo.csv` — biaxial fatigue data (n = 46)

Bootstrap (B = 2000) and leave-one-out results for OLS, PMM2, PMM3 and MLE-TN
on the biaxial fatigue dataset. The residuals are skewed (γ₃ ≈ −0.47), which
falls outside the symmetric-error setting of the paper.

| Column | Meaning |
|---|---|
| `mean`, `sd`, `q025`, `q975` | Bootstrap distribution of β̂₁ |
| `na_pct` | Percent of failed bootstrap fits |
| `g3_empirical` | Bootstrap Var(method) / Var(OLS) |
| `loo_mse`, `rel_loo` | LOO-MSE and LOO-MSE relative to OLS |
