# PMM Research Scripts

Templates and utilities for writing academic papers using the Polynomial Maximization Method.

## Available Scripts

### 1. `template_monte_carlo.R`
**Purpose**: Run Monte Carlo simulations to evaluate PMM2 performance vs OLS/MLE

**Usage**:
```bash
Rscript scripts/template_monte_carlo.R
```

**Features**:
- Parallel processing (multi-core support)
- Multiple sample sizes (50, 100, 200, 500)
- Various error distributions (Gamma, Exponential, Lognormal)
- Automatic LaTeX table generation
- Calculates: bias, variance, MSE, relative efficiency

**Outputs**:
- `results/monte_carlo_raw.csv` - Raw simulation results
- `results/monte_carlo_summary.csv` - Summary statistics
- `results/mc_comparison.tex` - LaTeX table for paper

**Customization**:
```r
# Edit these parameters in the script:
n_sim = 2000              # Number of replications
sample_sizes = c(50, 100, 200, 500)
seed = 12345             # For reproducibility
n_cores = 4              # Parallel processing cores
```

---

### 2. `template_paper_analysis.R`
**Purpose**: Complete analysis pipeline for empirical applications

**Usage**:
```bash
Rscript scripts/template_paper_analysis.R
```

**Features**:
- Data loading and preprocessing
- Multiple ARIMA model specifications
- CSS-ML vs PMM2 comparison
- Diagnostic tests (Ljung-Box, Jarque-Bera, ARCH)
- Publication-quality plots (ggplot2)
- Automated narrative report generation

**Outputs**:
- `results/tables/descriptive_stats.csv`
- `results/tables/model_comparison.csv`
- `results/tables/diagnostics.csv`
- `results/plots/01_time_series.png`
- `results/plots/02_aic_comparison.png`
- `results/plots/03_qq_plot.png`
- `results/plots/04_acf_residuals.png`
- `results/ANALYSIS_REPORT.md`

**Customization**:
Replace the data loading section with your dataset:
```r
# Load your data
data <- read.csv("data/your_dataset.csv")
y <- data$your_variable

# Or use FRED API
library(fredr)
fredr_set_key("your_api_key")
y <- fredr(series_id = "DCOILWTICO")$value
```

---

## Quick Start Workflow

### For a New Paper:

1. **Setup Project Structure**:
```bash
mkdir -p data scripts results/plots results/tables
cp scripts/template_paper_analysis.R scripts/my_analysis.R
```

2. **Add Your Data**:
```bash
# Place your dataset in data/
cp ~/Downloads/my_data.csv data/
```

3. **Edit Analysis Script**:
```r
# Edit scripts/my_analysis.R
# Update data loading section (STEP 1)
y <- read.csv("data/my_data.csv")$value
```

4. **Run Analysis**:
```bash
Rscript scripts/my_analysis.R
```

5. **Review Outputs**:
```bash
open results/ANALYSIS_REPORT.md
open results/plots/
```

6. **Run Monte Carlo** (optional):
```bash
Rscript scripts/template_monte_carlo.R
```

---

## Requirements

### R Packages (Install Once):
```r
install.packages(c(
  "EstemPMM",      # PMM estimators
  "ggplot2",       # Plots
  "gridExtra",     # Multiple plots
  "dplyr",         # Data manipulation
  "xtable",        # LaTeX tables
  "moments",       # Skewness/kurtosis
  "FinTS",         # ARCH test
  "parallel"       # Parallel processing
))

# EstemPMM from GitHub (if not on CRAN)
remotes::install_github("SZabolotnii/EstemPMM")
```

---

## Script Customization Guide

### Adding New Distribution to Monte Carlo:
```r
# In template_monte_carlo.R, add to distributions list:
distributions = list(
  gamma = function(n) rgamma(n, shape = 2, scale = 1) - 2,
  exp = function(n) rexp(n, rate = 1) - 1,
  # ADD NEW:
  t_dist = function(n) rt(n, df = 5) / sqrt(5/3),  # Standardized t(5)
  weibull = function(n) rweibull(n, shape = 1.5, scale = 1) - gamma(1 + 1/1.5)
)
```

### Adding New Model Specification:
```r
# In template_paper_analysis.R, add to orders list:
orders <- list(
  c(1, 0, 0),  # AR(1)
  c(2, 0, 0),  # AR(2)
  c(1, 0, 1),  # ARMA(1,1)
  c(2, 0, 1),  # ARMA(2,1)
  # ADD NEW:
  c(3, 0, 0),  # AR(3)
  c(1, 0, 2)   # ARMA(1,2)
)
```

### Changing Plot Style:
```r
# Add custom theme to all plots:
my_theme <- theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Apply to plots:
p1 <- p1 + my_theme
```

---

## Troubleshooting

**Error: "EstemPMM package not found"**
```r
# Install from GitHub:
install.packages("remotes")
remotes::install_github("SZabolotnii/EstemPMM")
```

**Error: "Cannot allocate vector of size..."**
```r
# Reduce simulation replications:
n_sim = 500  # Instead of 2000
```

**Plots not displaying correctly**:
```r
# Save to file instead of display:
ggsave("plot.png", plot_object, width = 8, height = 6, dpi = 300)
```

**Parallel processing errors**:
```r
# Reduce cores or disable parallel:
n_cores = 1  # Single-threaded
# Or use: mclapply(..., mc.cores = 1)
```

---

## Contributing

To add new templates or improve existing ones:

1. Create new script: `scripts/template_your_feature.R`
2. Document usage in this README
3. Add example outputs
4. Test on multiple datasets

---

## Citation

If you use these scripts in your research, please cite:

```bibtex
@software{pmm_research_templates,
  title = {PMM Research Script Templates},
  author = {Skill Seeker},
  year = {2025},
  url = {https://github.com/yourusername/your-repo}
}
```

And cite the original PMM papers (see main SKILL.md for references).
