# Releases

Version history for this repository (2 releases).

## v1.1.0: Release v1.1.0
**Published:** 2025-11-05



[View on GitHub](https://github.com/SZabolotnii/PMM2-ARIMA/releases/tag/v1.1.0)

---

## PMM: PMM2-ARIMA
**Published:** 2025-11-05
**Pre-release**

This release provides the complete reproducibility package for the manuscript:

Zabolotnii S. (2025). "Applying the Polynomial Maximization Method to Estimate ARIMA Models with Asymmetric Non-Gaussian Innovations" (under review).

The package includes:
- Full R implementation of the PMM2 estimator for ARIMA(p,d,q) models,
- WTI crude oil case study replication,
- Monte Carlo simulation framework (128,000+ simulations),
- Scripts for generating publication-ready figures and tables.

### Repository Structure
- data/ — immutable WTI series (DCOILWTICO.csv)
- scripts/ — R entry points for experiments, simulations, diagnostics, and report creation
- results/ — automatically generated outputs (CSV, RDS, plots, markdown reports)
- EstemPMM2-lib/ — bundled EstemPMM v0.1.1 for guaranteed reproducibility

### Reproducibility Instructions
Required: R ≥ 4.3

1. Install dependencies:
   install.packages(c("ggplot2","gridExtra","RColorBrewer","tseries","knitr","MASS"))
   install.packages("EstemPMM2-lib/EstemPMM_0.1.1.tar.gz", repos = NULL, type = "source")

2. Reproduce the WTI case study:
   Rscript scripts/run_full_study.R

3. Run Monte Carlo simulations (optional):
   Rscript scripts/run_monte_carlo.R

Outputs (figures, tables, summaries) will appear in `results/`.

### Citation
Please cite both:
- The article (once published), and
- This archived release via DOI.


[View on GitHub](https://github.com/SZabolotnii/PMM2-ARIMA/releases/tag/PMM)

---

