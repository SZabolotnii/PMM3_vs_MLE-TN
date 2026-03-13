---
name: pmm-research
description: Research assistant for writing academic papers about Polynomial Maximization Method (PMM) - parameter estimation with asymmetric non-Gaussian distributions
---

# PMM Research & Paper Writing Skill

Comprehensive assistance for conducting research and writing academic papers using the **Polynomial Maximization Method (МПП - Метод Максимізації Поліномів)** for parameter estimation in statistical models with asymmetric non-Gaussian error distributions.

## When to Use This Skill

This skill should be triggered when:
- Writing academic papers about PMM methodology
- Implementing PMM estimators for regression, time series, or other statistical models
- Conducting Monte Carlo simulations to evaluate PMM performance
- Comparing PMM with classical methods (OLS, MLE, CSS)
- Preparing reproducibility packages for journal submissions
- Creating statistical analysis with asymmetric error distributions
- Developing R packages for PMM-based estimation
- Working on papers about variance reduction in parameter estimation

## Key Concepts

### Polynomial Maximization Method (PMM)
- **Core Principle**: Uses polynomials of degree S to estimate distribution parameters
- **PMM1 (S=1)**: Equivalent to classical Ordinary Least Squares (OLS)
- **PMM2 (S=2)**: Provides variance reduction for asymmetric distributions
- **Variance Reduction Formula**: `g = 1 - c₃² / (2 + c₄)` where c₃ is skewness, c₄ is kurtosis
- **Original Developer**: Yu. P. Kunchenko (1992)

### Application Domains
1. **Linear Regression** - Parameter estimation with non-Gaussian errors
2. **Time Series (ARIMA)** - AR, MA, ARMA, ARIMA models with asymmetric innovations
3. **Autoregressive Models** - AR(p) with non-normal residuals
4. **Moving Average Models** - MA(q) with skewed errors
5. **Economic/Financial Data** - WTI crude oil prices, stock returns, volatility modeling

### Performance Characteristics
| Distribution     | Skewness | Kurtosis | Theoretical Efficiency | Actual Improvement |
|-----------------|----------|----------|----------------------|-------------------|
| Gamma (α=0.5)   | 2.83     | 12       | 57%                  | ~50%              |
| Exponential     | 2.00     | 6        | 50%                  | ~45%              |
| Gamma (α=2)     | 1.41     | 3        | 40%                  | ~35%              |
| Lognormal       | 1.00     | 1.5      | 29%                  | ~25%              |

## Quick Reference

### PMM2 Linear Regression (R)

```r
library(EstemPMM)

# 1. Create data with asymmetric errors
n <- 100
x <- rnorm(n)
errors <- rgamma(n, shape = 2, scale = 1) - 2
y <- 2 + 1.5 * x + errors
data <- data.frame(x = x, y = y)

# 2. Fit model using PMM2
fit <- lm_pmm2(y ~ x, data = data)

# 3. Review results and compare with OLS
summary(fit)
compare_with_ols(y ~ x, data)

# 4. Check variance reduction
vf <- pmm2_variance_factor(fit@m2, fit@m3, fit@m4)
cat("Variance reduction factor g:", vf$g, "\n")
```

### PMM2 ARIMA Models (R)

```r
library(EstemPMM)

# AR(1) model with PMM2
ar_fit <- ar_pmm2(y, order = 1, method = "pmm2")

# ARMA(1,1) model
arma_fit <- arma_pmm2(y, order = c(1, 1))

# ARIMA(1,1,1) model
arima_fit <- arima_pmm2(y, order = c(1, 1, 1))

# Compare methods
compare_ar_methods(y, p = 1)
```

### Monte Carlo Simulation Template

```r
# Monte Carlo simulation for PMM2 performance evaluation
run_simulation <- function(n = 100, reps = 1000, dist = "gamma") {
  results <- matrix(NA, nrow = reps, ncol = 4)
  colnames(results) <- c("OLS_beta", "PMM2_beta", "OLS_var", "PMM2_var")

  for (i in 1:reps) {
    # Generate data
    x <- rnorm(n)
    if (dist == "gamma") {
      errors <- rgamma(n, shape = 2) - 2
    } else if (dist == "exp") {
      errors <- rexp(n, rate = 1) - 1
    }
    y <- 1 + 0.5 * x + errors

    # Fit models
    ols <- lm(y ~ x)
    pmm2 <- lm_pmm2(y ~ x, data = data.frame(x, y))

    # Store results
    results[i, 1] <- coef(ols)[2]
    results[i, 2] <- coef(pmm2)[2]
    results[i, 3] <- vcov(ols)[2, 2]
    results[i, 4] <- vcov(pmm2)[2, 2]
  }

  return(results)
}

# Run and analyze
results <- run_simulation(n = 200, reps = 2000, dist = "gamma")
variance_ratio <- mean(results[, 4]) / mean(results[, 3])
cat("Empirical variance ratio (PMM2/OLS):", variance_ratio, "\n")
```

### Paper Structure Template

```markdown
# Paper: Applying PMM to [Your Model Type]

## 1. Introduction
- Motivation: Why asymmetric errors matter
- Research gap: What existing methods miss
- Contribution: How PMM improves estimation

## 2. Methodology
### 2.1 Polynomial Maximization Method
- Theoretical foundation (Kunchenko, 1992)
- Variance reduction formula
- Adaptive estimation procedure

### 2.2 Model Specification
- [Your specific model: ARIMA, regression, etc.]
- Assumptions and conditions
- Parameter space

### 2.3 Implementation
- Algorithm description
- Computational complexity
- R/Python implementation

## 3. Monte Carlo Study
### 3.1 Simulation Design
- Sample sizes: n = {50, 100, 200, 500}
- Replications: 1000-5000
- Distributions: Gamma, Exponential, Lognormal

### 3.2 Performance Metrics
- Bias: E(θ̂) - θ
- Variance: Var(θ̂)
- MSE: E[(θ̂ - θ)²]
- Relative efficiency: Var(OLS)/Var(PMM2)

## 4. Empirical Application
- Real dataset description
- Diagnostic tests (stationarity, residual analysis)
- Model comparison (AIC, BIC, RMSE)
- Out-of-sample forecasting

## 5. Results
- Tables: Parameter estimates, standard errors, efficiency gains
- Figures: Convergence plots, residual diagnostics, QQ-plots

## 6. Conclusion
- Summary of findings
- Practical implications
- Future research directions
```

## Reference Files

### Available Documentation
- **references/README.md** - PMM2-ARIMA reproducibility package guide
- **references/releases.md** - Version history and release notes (2 releases)
- **references/file_structure.md** - Complete repository structure (90 items)

### Key Resources from References
1. **EstemPMM R Package** - Production-ready implementation
2. **PMM2-ARIMA Repository** - ARIMA models with PMM2 estimator
3. **Monte Carlo Framework** - 128,000+ simulation template
4. **WTI Case Study** - Real-world application example
5. **LaTeX Templates** - Article formatting and structure

## Working with This Skill

### For New Papers

**Step 1: Define Research Question**
```
- What model are you estimating? (regression, ARIMA, GLM, etc.)
- What type of asymmetry do you expect? (right-skewed, left-skewed)
- What is your comparison baseline? (OLS, MLE, CSS-ML)
```

**Step 2: Implement PMM Estimator**
```r
# Start with EstemPMM package for regression
library(EstemPMM)

# Or adapt PMM2-ARIMA code for time series
source("scripts/pmm2_estimation.R")

# Or implement from scratch using cumulant formulas
# See references/README.md for mathematical details
```

**Step 3: Design Monte Carlo Study**
```r
# Use template above
# Vary: sample size, distribution, parameters
# Replications: 1000+ for publication quality
# Seed: Set for reproducibility (e.g., seed = 12345)
```

**Step 4: Find Real Dataset**
```r
# Economic: GDP, inflation, oil prices (FRED database)
# Financial: Stock returns, volatility (Yahoo Finance)
# Environmental: Temperature, rainfall (NOAA)
# Medical: Survival times, biomarkers
```

**Step 5: Create Reproducibility Package**
```
data/           # Raw datasets
scripts/        # R scripts for all analyses
results/        # Auto-generated outputs
EstemPMM-lib/   # Bundled package version
README.md       # Replication instructions
DESCRIPTION     # R dependencies
```

### For Code Development

**R Package Structure**
```r
R/
  ├── pmm2_main.R         # Core estimation functions
  ├── pmm2_classes.R      # S4 class definitions
  ├── pmm2_utils.R        # Helper functions
  ├── pmm2_inference.R    # Bootstrap/asymptotic inference
  ├── pmm2_diagnostics.R  # Model diagnostics
  └── pmm2_simulation.R   # Monte Carlo utilities

man/                      # Roxygen2 documentation
tests/                    # Unit tests
vignettes/                # Tutorial vignettes
```

**Essential Functions to Implement**
1. `pmm2_fit()` - Core estimator
2. `pmm2_vcov()` - Variance-covariance matrix
3. `pmm2_cumulants()` - Sample cumulant estimation
4. `pmm2_variance_factor()` - Theoretical efficiency g
5. `compare_methods()` - PMM vs OLS/MLE comparison

### For Paper Writing

**Required Sections**
1. **Abstract** (150-250 words)
   - Problem, method, results, conclusion

2. **Introduction** (3-4 pages)
   - Motivation with real examples
   - Literature review (cite Kunchenko, Zabolotnii papers)
   - Research gap and contribution

3. **Methodology** (4-6 pages)
   - PMM theoretical foundation
   - Your specific model adaptation
   - Computational algorithm

4. **Simulations** (3-5 pages)
   - Design justification
   - Results tables (bias, variance, MSE)
   - Convergence and robustness checks

5. **Application** (3-4 pages)
   - Dataset description and diagnostics
   - Estimation results
   - Model comparison
   - Forecasting performance

6. **Conclusion** (1-2 pages)
   - Summary and practical implications
   - Limitations and extensions

**Key Citations**
```bibtex
@article{kunchenko1992,
  author = {Kunchenko, Yu. P. and Lega, Yu. G.},
  title = {Polynomial parameter estimation of close to Gaussian random variables},
  journal = {Radioelectronics and Communications Systems},
  year = {1992}
}

@inproceedings{zabolotnii2018pmm,
  author = {Zabolotnii, S. and Warsza, Z.L. and Tkachenko, O.},
  title = {Polynomial Estimation of Linear Regression Parameters for the Asymmetric PDF of Errors},
  booktitle = {Automation 2018},
  publisher = {Springer},
  year = {2018},
  doi = {10.1007/978-3-319-77179-3_75}
}

@article{zabolotnii2025arima,
  author = {Zabolotnii, Serhii},
  title = {Applying the Polynomial Maximization Method to Estimate ARIMA Models with Asymmetric Non-Gaussian Innovations},
  journal = {Submitted for publication},
  year = {2025},
  doi = {10.5281/zenodo.17529589}
}
```

## Common Patterns

### Pattern 1: Quick Efficiency Check
```r
# Before full Monte Carlo, check theoretical efficiency
c3 <- 1.5  # Your estimated skewness
c4 <- 3.0  # Your estimated kurtosis
g <- 1 - c3^2 / (2 + c4)
cat("Theoretical variance reduction:", (1 - g) * 100, "%\n")
# If g < 0.8, PMM2 should provide >20% variance reduction
```

### Pattern 2: Reproducibility Workflow
```r
# Master script: scripts/run_full_study.R
source("scripts/00_setup.R")              # Load packages
source("scripts/01_data_preparation.R")   # Import/clean data
source("scripts/02_estimation.R")         # Fit all models
source("scripts/03_monte_carlo.R")        # Run simulations
source("scripts/04_create_tables.R")      # Generate LaTeX tables
source("scripts/05_create_figures.R")     # Create publication plots
source("scripts/06_generate_report.R")    # Compile narrative report
```

### Pattern 3: Model Comparison Table
```r
# Create comparison table for paper
results_df <- data.frame(
  Model = c("OLS", "PMM2", "MLE"),
  Estimate = c(ols_est, pmm2_est, mle_est),
  SE = c(ols_se, pmm2_se, mle_se),
  AIC = c(ols_aic, pmm2_aic, mle_aic),
  BIC = c(ols_bic, pmm2_bic, mle_bic)
)

# Export to LaTeX
library(xtable)
print(xtable(results_df, digits = 4),
      file = "results/method_comparison.tex")
```

## Best Practices

### Research Workflow
1. ✅ **Set seed** for all random operations (reproducibility)
2. ✅ **Version control** with Git (track all changes)
3. ✅ **Document dependencies** in DESCRIPTION/requirements.txt
4. ✅ **Bundle package versions** (e.g., EstemPMM v0.1.1)
5. ✅ **Archive on Zenodo** for DOI and long-term preservation
6. ✅ **Share code publicly** (GitHub) before submission

### Statistical Rigor
1. ✅ **Test stationarity** for time series (ADF, KPSS tests)
2. ✅ **Check residuals** (normality, autocorrelation, heteroscedasticity)
3. ✅ **Report all methods** (CSS, CSS-ML, PMM2 for ARIMA)
4. ✅ **Use multiple metrics** (AIC, BIC, RMSE, MAE)
5. ✅ **Out-of-sample validation** (not just in-sample fit)
6. ✅ **Bootstrap inference** when analytical se not available

### Code Quality
1. ✅ **Modular functions** (one function = one task)
2. ✅ **Roxygen documentation** for all exported functions
3. ✅ **Unit tests** with testthat (>80% coverage)
4. ✅ **Clear variable names** (no `x1`, `temp`, `res2`)
5. ✅ **Vectorized operations** (avoid explicit loops when possible)
6. ✅ **Progress bars** for long-running simulations

## Resources

### Official Repositories
- **EstemPMM**: https://github.com/SZabolotnii/EstemPMM (CRAN package)
- **PMM2-ARIMA**: https://github.com/SZabolotnii/PMM2-ARIMA (Reproducibility package)

### Key Papers
1. Kunchenko & Lega (1992) - Original PMM formulation
2. Zabolotnii et al. (2018) - PMM for asymmetric regression errors
3. Zabolotnii et al. (2022) - PMM for autoregressive models
4. Zabolotnii et al. (2023) - PMM for moving average models
5. Zabolotnii (2025) - PMM for ARIMA models (under review)

### Datasets for Examples
- **WTI Crude Oil Prices**: FRED database (DCOILWTICO)
- **Auto MPG**: UCI Machine Learning Repository
- **Financial Returns**: Yahoo Finance API
- **Economic Indicators**: FRED, World Bank, Eurostat

## Notes

- PMM2 is most effective when skewness |c₃| > 1.0
- For symmetric distributions (c₃ ≈ 0), PMM2 ≈ OLS
- Computational cost is similar to OLS (no major overhead)
- Bootstrap recommended for small samples (n < 50)
- Always report both OLS and PMM2 for transparency

## Updating This Skill

To enhance with new research:
1. Add new PMM applications (GLM, survival analysis, etc.)
2. Update with published papers and citations
3. Include new simulation scenarios
4. Add real-world case studies
5. Expand code templates for other languages (Python, Julia)

---

**Generated by Skill Seeker** | Optimized for PMM Research & Academic Writing
