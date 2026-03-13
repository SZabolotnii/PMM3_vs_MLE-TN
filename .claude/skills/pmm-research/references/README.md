# PMM2ARIMA Reproducibility Package

[![GitHub](https://img.shields.io/badge/GitHub-PMM2--ARIMA-blue?logo=github)](https://github.com/SZabolotnii/PMM2-ARIMA)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17529589.svg)](https://doi.org/10.5281/zenodo.17529589)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%E2%89%A54.3-blue)](https://www.r-project.org/)

This repository contains all scripts, data, and tools needed to reproduce the empirical results reported in the manuscript **"Applying the Polynomial Maximization Method to Estimate ARIMA Models with Asymmetric Non-Gaussian Innovations"** by Serhii Zabolotnii.

**Repository:** [https://github.com/SZabolotnii/PMM2-ARIMA](https://github.com/SZabolotnii/PMM2-ARIMA)

The package includes both the WTI crude oil case study and the comprehensive Monte Carlo simulation experiments (128,000+ simulations) demonstrating the PMM2 method's performance across various non-Gaussian innovation distributions.

## Repository Layout

### Public Repository Contents

This public repository contains all materials needed to reproduce the empirical results:

- `data/` – immutable inputs delivered with the package (`DCOILWTICO.csv` from FRED).
- `scripts/` – R entry points that orchestrate analysis, plotting, and reporting.
- `results/` – auto-generated artefacts (CSV tables, RDS stores, plots, Markdown report, Monte Carlo summaries).
- `EstemPMM2-lib/` – packaged version of the EstemPMM library (v0.1.1) for reproducibility.
- `DESCRIPTION` – declared R dependencies for quick inspection/tooling.

### Excluded from Public Repository

The following directories are maintained locally for development but excluded from the public repository:

- `LaTeX/` – draft manuscript versions and LaTeX source files (work in progress).
- `Experiment_R/` – experimental R code and alternative implementations.
- `Experiment_Python/` – Python-based experiments and comparative analyses.
- `data-search/` – preliminary data exploration and source documentation.

## Requirements

- R ≥ 4.3 (tested with 4.5.1).
- R packages: `EstemPMM` (version 0.1.1), `ggplot2`, `gridExtra`, `RColorBrewer`, `tseries`, `knitr`, `MASS`.
  - Install CRAN packages: `install.packages(c("ggplot2", "gridExtra", "RColorBrewer", "tseries", "knitr", "MASS"))`
  - Install `EstemPMM` version 0.1.1 from the included archive:
    `install.packages("EstemPMM2-lib/EstemPMM_0.1.1.tar.gz", repos = NULL, type = "source")`
  - Alternatively, install from GitHub: `remotes::install_github("SZabolotnii/EstemPMM")`
- Optional: `pandoc` if you plan to convert the Markdown report to HTML/PDF, and a LaTeX distribution (`TeXLive`, `TinyTeX`, …) to compile `latex/PMM2_ARIMA.tex`.

## Reproducing the WTI Experiment

From the root of this directory, run:

```bash
Rscript scripts/run_full_study.R
```

The command will:

1. Fit all ARIMA specifications with CSS-ML vs PMM2 (`scripts/comprehensive_study.R`).
2. Regenerate the 10 publication plots (`scripts/create_visualizations.R`).
3. Rebuild the narrative Markdown report (`scripts/generate_report.R`).

You can also execute any step individually:

```bash
# Quick console-only sanity check
Rscript scripts/arima_oil_quick_demo.R

# Full modelling pipeline
Rscript scripts/comprehensive_study.R

# Visuals and report (after the previous step)
Rscript scripts/create_visualizations.R
Rscript scripts/generate_report.R
```

## Output Map

After a successful run you should see:

- `results/full_results.csv` – comprehensive model comparison results for all ARIMA specifications.
- `results/method_comparison.csv` – deltas between CSS-ML and PMM2 methods.
- `results/descriptive_stats.csv` – descriptive statistics for the WTI crude oil series.
- `results/plots/*.png` – publication-quality figures (numbered 1–10).
- `results/ANALYTICAL_REPORT.md` – extended narrative analysis of results.
- `results/monte_carlo/*.csv` – simulation summaries including per-model metrics, relative efficiencies, and residual cumulants.
- `results/monte_carlo/article_comparison.csv` – comparison of simulation results with published values.

All CSV files contain raw outputs from R for full transparency and reproducibility.

## Monte Carlo Simulations

The dedicated driver `scripts/run_monte_carlo.R` reconstructs the simulation study (default: 2000 replications per combination, seed = 12345). Expect the full run to take 20–40 minutes depending on hardware.

```bash
# Full Monte Carlo run
Rscript scripts/run_monte_carlo.R

# Customization (fewer replications for quick verification)
Rscript scripts/run_monte_carlo.R --reps=200 --seed=20250101

# Variation of assumptions
Rscript scripts/run_monte_carlo.R --standardize-innov=false --css-method=CSS-ML
```

Options `--standardize-innov` and `--css-method` allow checking sensitivity to innovation scale and type of classical estimate (available values: `CSS`, `CSS-ML`, `ML`). The `--m-est=false` flag disables additional robust estimates, which are computed by default for AR components.

Key artifacts:

- `monte_carlo_metrics.csv` – complete cross-section of metrics (bias, var, MSE, RE, VR) for each model, parameter, and distribution.
- `arima110_summary.csv`, `arima110_re_vs_sample_size.csv` – ARIMA(1,1,0) results and RE dependence on sample size.
- `arima011_summary.csv`, `arima111_summary.csv`, `arima210_summary.csv` – condensed comparisons for other model configurations.
- `arima110_residual_cumulants.csv` – average cumulants of PMM2 residuals for normality assessment.
- `article_comparison.csv` – comparison of simulation results with published values including deviations.
- The `M-EST` column in summary files contains results of Huber M-estimates for AR components, allowing comparison of PMM2/CSS with robust estimates.

## Bootstrap Confidence Intervals (NEW)

After running Monte Carlo simulations, you can add bootstrap confidence intervals to the results:

```bash
# Add bootstrap CIs to Monte Carlo results
Rscript scripts/add_confidence_intervals.R

# Customization
Rscript scripts/add_confidence_intervals.R --bootstrap-reps=2000 --confidence=0.95
```

This creates `monte_carlo_metrics_with_ci.csv` containing standard errors and confidence intervals for bias, variance, MSE, and relative efficiency. Addresses reviewer request for uncertainty quantification in Monte Carlo tables.

## Out-of-Sample Validation (NEW)

Evaluate forecasting performance on held-out data:

```bash
# Fixed train/test split (80/20) + rolling window forecasts
Rscript scripts/wti_out_of_sample.R

# Customization
Rscript scripts/wti_out_of_sample.R --train-fraction=0.8 --window-size=100
```

Outputs:
- `wti_fixed_split_validation.csv` – RMSE/MAE for train/test split
- `wti_rolling_window_validation.csv` – rolling window forecast metrics

Addresses reviewer request for out-of-sample validation.

## Enhanced Diagnostics (NEW)

Generate comprehensive diagnostic plots and statistics with p-values:

```bash
# Full diagnostic suite for WTI case study
Rscript scripts/wti_diagnostics.R

# For different model specification
Rscript scripts/wti_diagnostics.R --order=1,1,1
```

Outputs:
- `wti_diagnostics_statistics.csv` – Ljung-Box, Jarque-Bera, Shapiro-Wilk, ARCH tests with p-values
- `wti_qq_plots.png` – Q-Q plots for residual normality assessment
- `wti_acf_pacf.png` – autocorrelation diagnostics
- `wti_residual_time_series.png` – residual plots over time
- `wti_residual_histograms.png` – residual distributions with normal overlay

Addresses reviewer request for diagnostic plots with p-values.

## Software Versions

- **R:** Version 4.5.1 (2025-06-13)
- **EstemPMM:** Version 0.1.1 (included in `EstemPMM2-lib/EstemPMM_0.1.1.tar.gz`)
- **Platform:** macOS Sequoia 15.6.1 (aarch64-apple-darwin24.4.0)
- **Random seed:** All scripts use `set.seed(12345)` for reproducibility
- **Session info:** Complete R environment details available in `sessionInfo.txt`

## Citation

If you use this code or data in your research, please cite:

```bibtex
@article{zabolotnii2025pmm2arima,
  title={Applying the Polynomial Maximization Method to Estimate ARIMA Models with Asymmetric Non-Gaussian Innovations},
  author={Zabolotnii, Serhii},
  journal={Submitted for publication},
  year={2025},
  url={https://github.com/SZabolotnii/PMM2-ARIMA},
  doi={10.5281/zenodo.17529589}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Serhii Zabolotnii**
Cherkasy State Business College, Cherkasy, Ukraine
Email: zabolotnii.serhii@csbc.edu.ua

## Acknowledgments

This research builds upon the Polynomial Maximization Method framework and adapts it to time series analysis with non-Gaussian innovations. We thank the R community for the excellent statistical computing tools that made this research possible.

## Housekeeping

- Use `DESCRIPTION` with `renv` or `pak` if you wish to create a fully locked environment before sharing.
