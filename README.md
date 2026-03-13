# PMM3 vs MLE-TN: Comparative Analysis

**Reference:** Salinas H., Bakouch H., Qarmalah N., Martínez-Flórez G. (2023).  
*A Flexible Class of Two-Piece Normal Distribution with a Regression Illustration  
to Biaxial Fatigue Data.* Mathematics, 11(5), 1271.  
https://doi.org/10.3390/math11051271

**Context:** Part of the EstemPMM R-package development.  
PMM3 = Polynomial Maximization Method (S=3) for symmetric non-Gaussian errors.

---

## Project structure

```
PMM3_vs_MLE-TN/
├── R/
│   ├── config.R               ← edit this to change MC parameters
│   ├── 01_tn_distribution.R   ← rtn, dtn, tn_cumulants, build_table_T1
│   ├── 02_mle_tn.R            ← mle_tn() — competitor method
│   ├── 03_monte_carlo.R       ← run_full_mc_study() with caching
│   ├── 04_real_data.R         ← biaxial fatigue, bootstrap, LOO-CV
│   ├── 05_robustness.R        ← misspecification & bimodal tests
│   ├── 06_visualization.R     ← ggplot2 figures
│   └── pmm3_skeleton.R        ← [if needed] stub lm_pmm3() until EstemPMM updated
├── tests/testthat/
│   ├── test-01-tn-distribution.R
│   ├── test-02-mle-tn.R
│   └── test-03-real-data.R
├── data/
│   ├── raw/                   ← raw data files
│   └── processed/             ← processed data
├── results/
│   ├── tables/                ← CSV output (T1, T2, T6, T7)
│   ├── figures/               ← PNG plots
│   └── mc_cache/              ← RDS cache (excluded from git)
├── reports/
│   └── main_report.Rmd        ← full HTML analysis report
├── scripts/
│   ├── install_deps.R         ← one-time package installation
│   └── run_all.R              ← full study runner
└── docs/                      ← additional documentation
```

---

## Quick start

```r
# 1. Install dependencies (once)
source("scripts/install_deps.R")

# 2. Load project
source("R/config.R")
load_project()

# 3. Smoke test: theoretical table
T1 <- build_table_T1()
print(T1[, c("lambda", "gamma4", "g3", "improvement_pct")])

# 4. Run unit tests
library(testthat)
test_dir("tests/testthat")

# 5. Full study (~5 min with M=500, ~60 min with M=5000)
source("scripts/run_all.R")
```

---

## Study plan (6 blocks)

| Block | Description | Output |
|-------|-------------|--------|
| 1 | Theoretical TN cumulants | Table T1 |
| 2 | Monte Carlo: ARE, g₃ | Table T2 |
| 3 | Biaxial fatigue (n=46) | Table T3, bootstrap, LOO-CV |
| 4 | nlme::Fatigue (fallback) | Table T3' |
| 5 | Misspecification robustness | Table T6 |
| 6 | Bimodal TN (λ > 1) | Table T7 |

---

## Key results (preliminary)

- TN is always platykurtic: γ₄(λ) < 0 → PMM3 always gains over OLS
- Maximum theoretical gain: ~50% variance reduction at large λ
- Biaxial fatigue: γ̂₃ ≈ −0.47 (asymmetric!) → PMM2 may outperform PMM3
- PMM3 expected to be more robust than MLE-TN under misspecification

---

## Notes on PMM3 implementation

If `lm_pmm3` is not yet available in `EstemPMM`, the stub in `R/pmm3_skeleton.R`
is loaded automatically. The stub returns OLS estimates with a warning and marks
MC results as NA to avoid misleading comparisons.

Mathematical specification for the full Newton-Raphson implementation is
documented in `R/pmm3_skeleton.R` (formulas 10, 11a,b,c from Zabolotnyi et al., 2019).
