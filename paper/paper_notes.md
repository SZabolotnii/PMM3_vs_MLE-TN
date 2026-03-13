# Paper Notes & Decisions Log

## Target Journal
**Mathematics (MDPI)** — open access, ISSN 2227-7390
- Why: Salinas et al. (2023) published there; same scope (statistics, mathematical methods)
- Article type: Article (not review)
- Template: MDPI LaTeX template (mdpi.cls) — to be added before submission
- Typical length: 6,000–10,000 words

## Working Title
"PMM3 vs MLE-TN: Semi-parametric vs Parametric Estimation in Linear Regression
with Two-piece Normal Errors"

Alternatives:
- "Comparison of Polynomial Maximization and Maximum Likelihood Estimators for
  Linear Regression under Platykurtic Error Distributions"
- "Semi-parametric Efficiency Gains in Linear Regression: PMM3 versus MLE-TN"

## Key Contribution (elevator pitch)
PMM3 is distribution-free and achieves up to 95% variance reduction (vs OLS)
for TN errors with λ≥2, while MLE-TN is fragile under misspecification and
adds 2 extra parameters to estimate. The empirical efficiency threshold is γ₄ ≲ −0.7.

## Main Claims (to support with tables/figures)
1. PMM3 variance reduction: g₃_empirical ≈ g₃_theoretical for n ≥ 100
2. MLE-TN dominates PMM3 only at λ ≥ 1.5 (ARE > 1.5) because it uses
   full parametric likelihood — but this requires correct model specification
3. Under misspecification (Uniform, Triangular, Logistic, t₁₀):
   PMM3 ARE ≈ 1.0 (neutral), MLE-TN loses efficiency (ARE ≤ 0.97)
4. Real data: iris-versicolor (γ₄ = −0.966): PMM3 gives 43% bootstrap
   variance reduction vs OLS, LOO-MSE 2.3% better than OLS
5. Empirical threshold: PMM3 > OLS requires γ₄ ≲ −0.7

## Open Questions
- [ ] Should we include PMM2 in main paper or appendix?
      Current plan: include PMM2 as baseline/control in simulation only
- [ ] iris-versicolor as sole real data example — is it sufficient?
      Consider: add whiteside Gas~Temp as negative example (γ₄ = −0.12)
- [ ] LOO-CV vs k-fold: stick with LOO given small n=50
- [ ] Should abstract be in English + Ukrainian? → No, Mathematics is English-only

## Figures Plan
1. **Fig 1**: TN density shapes for λ = {0.5, 1, 1.5, 2, 3} — show platykurtosis
2. **Fig 2**: Theoretical g₃ vs γ₄ curve — the variance reduction formula
3. **Fig 3**: Empirical g₃ vs n (convergence to theory) for λ = {1, 1.5, 2}
4. **Fig 4**: ARE comparison heatmap (λ × n) for PMM3 vs MLE-TN
5. **Fig 5**: iris-versicolor residual distribution + bootstrap g₃ distributions

## Tables Plan
- **T1**: TN distribution theoretical moments (λ, γ₄, g₃, improvement_pct)
- **T2**: MC results summary (n=100, n=200; λ = 1, 1.5, 2, 3): bias, ARE, g₃_emp
- **T3**: Misspecification robustness (ARE for OLS=1, PMM3≈1, MLE-TN<1)
- **T4**: Real data — iris-versicolor comparison (OLS, PMM3, MLE-TN)
- (T5 optional): Bimodal convergence test (Table T7 from simulation)

## Reviewer Concerns to Preempt
- Why iris-versicolor? → only dataset in standard R packages with γ₄ < −0.7
- Why not larger n for real data? → real-world constraint; LOO-CV validates
- PMM3 divergence for λ → ∞? → convergence 100% for all λ tested (T7 shows mean_iter≤10, na_pct=0)
- Is kappa estimation from OLS residuals consistent? → yes, but introduces finite-sample noise

## Timeline (target)
- Draft intro + methodology: Week 1
- Finalize tables/figures: Week 1-2
- Draft simulation + real data sections: Week 2
- Internal review: Week 3
- Submit: Week 4-5
