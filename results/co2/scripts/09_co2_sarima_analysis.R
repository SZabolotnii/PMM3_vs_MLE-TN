#!/usr/bin/env Rscript
# =============================================================================
# results/co2/scripts/09_co2_sarima_analysis.R
# PMM3-SARIMA applied to Mauna Loa CO2 (R built-in `co2`, n=468)
#
# Models tested:
#   A) ARIMA(0,1,1)(0,1,1)_12  — airline model (most common for CO2-like data)
#   B) ARIMA(1,1,0)(1,1,0)_12  — AR version
# Methods per model: CSS-ML, MLE, PMM3-SARIMA (new)
#
# Comparison metrics:
#   1. Residual diagnostics: γ₃, γ₄, g3_theoretical
#   2. Parameter estimates + shifts vs CSS
#   3. AIC / BIC
#   4. Hold-out RMSE (last 24 obs = 2 years, train on first 444)
#
# Output: results/co2/output/
# =============================================================================

# --- Setup -------------------------------------------------------------------
if (!exists("PROJECT_ROOT")) {
  if (file.exists("R/config.R")) {
    PROJECT_ROOT <- getwd()
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg)) {
      # script is at results/co2/scripts/09_... → root is ../../..
      script_path <- normalizePath(sub("--file=", "", file_arg))
      PROJECT_ROOT <- normalizePath(file.path(dirname(script_path), "../../.."),
                                    mustWork = FALSE)
    } else {
      PROJECT_ROOT <- getwd()
    }
  }
}
setwd(PROJECT_ROOT)
cat("Working directory:", getwd(), "\n")

# Load project (for co2 time-series object)
source("R/config.R")
load_project(verbose = FALSE)

# Load PMM3-SARIMA implementation
source("results/co2/R/pmm3_sarima.R")

OUT <- "results/co2/output"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("\n╔══════════════════════════════════════════╗\n")
cat("║  PMM3-SARIMA: Mauna Loa CO2 analysis    ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

# --- Data --------------------------------------------------------------------
y_full <- as.numeric(co2)      # 468 monthly observations
n_full <- length(y_full)

# Train/test split: hold out last 24 months
n_test  <- 24L
n_train <- n_full - n_test
y_train <- y_full[seq_len(n_train)]
y_test  <- y_full[seq(n_train + 1L, n_full)]

cat(sprintf("Series : n=%d (train=%d, test=%d)\n", n_full, n_train, n_test))
cat(sprintf("Range  : %.1f – %.1f ppm\n\n", min(y_full), max(y_full)))

# ---------------------------------------------------------------------------
# Helper: 1-step-ahead rolling forecasts for stats::arima fit
# ---------------------------------------------------------------------------
rolling_rmse <- function(fit_fn, y_train, y_test, ...) {
  n  <- length(y_train)
  m  <- length(y_test)
  sq <- numeric(m)
  y_all <- c(y_train, y_test)
  for (i in seq_len(m)) {
    fit_i  <- tryCatch(fit_fn(y_all[seq_len(n + i - 1L)], ...),
                       error = function(e) NULL)
    if (is.null(fit_i)) { sq[i] <- NA_real_; next }
    pred <- tryCatch(
      predict(fit_i, n.ahead = 1L)$pred[1L],
      error = function(e) NA_real_
    )
    sq[i] <- (y_all[n + i] - as.numeric(pred))^2
  }
  sqrt(mean(sq, na.rm = TRUE))
}

# ---------------------------------------------------------------------------
# Helper: pretty-print cumulants
# ---------------------------------------------------------------------------
print_cum <- function(cum, prefix = "  ") {
  cat(sprintf("%sγ₃ = %+.4f  γ₄ = %+.4f  γ₆ = %+.4f  g3_theor = %s  (Δ=%.1f%%)\n",
              prefix,
              cum$gamma3, cum$gamma4, cum$gamma6,
              if (is.na(cum$g3_theoretical)) "  NA" else sprintf("%.4f", cum$g3_theoretical),
              if (is.na(cum$improvement_pct)) 0 else cum$improvement_pct))
}

# ---------------------------------------------------------------------------
# Model specifications
# ---------------------------------------------------------------------------
models <- list(
  airline = list(
    name     = "ARIMA(0,1,1)(0,1,1)_12",
    order    = c(0L, 1L, 1L),
    seasonal = c(0L, 1L, 1L, 12L)
  ),
  ar_seas = list(
    name     = "ARIMA(1,1,0)(1,1,0)_12",
    order    = c(1L, 1L, 0L),
    seasonal = c(1L, 1L, 0L, 12L)
  )
)

all_results <- list()   # collect summary rows

# ===========================================================================
# Main loop over models
# ===========================================================================
for (mdl in models) {
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("Model:", mdl$name, "\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

  ord <- mdl$order
  ses <- mdl$seasonal

  # ---- 1. CSS-ML fit --------------------------------------------------------
  css_fit <- stats::arima(y_full, order = ord,
                          seasonal = list(order = ses[1:3], period = ses[4]),
                          method = "CSS-ML")
  cat("CSS-ML coefficients:\n")
  print(round(coef(css_fit), 5))
  eps_css <- as.numeric(residuals(css_fit))
  cum_css <- .sarima_pmm3_cumulants(eps_css[is.finite(eps_css)])
  cat("CSS-ML residual cumulants:\n"); print_cum(cum_css)
  cat(sprintf("  AIC=%.2f  BIC=%.2f  σ²=%.5f\n\n", AIC(css_fit), BIC(css_fit), css_fit$sigma2))

  # ---- 2. MLE fit -----------------------------------------------------------
  mle_fit <- tryCatch(
    stats::arima(y_full, order = ord,
                 seasonal = list(order = ses[1:3], period = ses[4]),
                 method = "ML"),
    error = function(e) { message("MLE failed: ", e$message); NULL }
  )
  if (!is.null(mle_fit)) {
    cat("MLE coefficients:\n")
    print(round(coef(mle_fit), 5))
    eps_mle <- as.numeric(residuals(mle_fit))
    cum_mle <- .sarima_pmm3_cumulants(eps_mle[is.finite(eps_mle)])
    cat("MLE residual cumulants:\n"); print_cum(cum_mle)
    cat(sprintf("  AIC=%.2f  BIC=%.2f  σ²=%.5f\n\n", AIC(mle_fit), BIC(mle_fit), mle_fit$sigma2))
  }

  # ---- 3. PMM3-SARIMA fit ---------------------------------------------------
  cat("Fitting PMM3-SARIMA...\n")
  pmm3_fit <- tryCatch(
    pmm3_sarima(y_full, order = ord, seasonal = ses, trace = TRUE),
    error = function(e) { message("PMM3-SARIMA failed: ", e$message); NULL }
  )

  if (!is.null(pmm3_fit)) {
    cat("\nPMM3-SARIMA coefficients:\n")
    print(round(pmm3_fit$coefficients, 5))
    cat("PMM3-SARIMA residual cumulants (final):\n")
    print_cum(pmm3_fit$cumulants_final)
    cat(sprintf("  AIC=%.2f  BIC=%.2f  σ²=%.5f\n",
                pmm3_fit$aic, pmm3_fit$bic, pmm3_fit$sigma2))
    cat(sprintf("  κ=%.4f  converged=%s  iter=%d\n\n",
                pmm3_fit$kappa, pmm3_fit$convergence, pmm3_fit$iterations))

    # Coefficient shift vs CSS
    shift <- pmm3_fit$coefficients - coef(css_fit)
    cat("PMM3 vs CSS coefficient shifts:\n")
    print(round(shift, 6))
    cat("\n")
  }

  # ---- 4. Hold-out RMSE comparison ------------------------------------------
  cat("Computing hold-out RMSE (last 24 months)...\n")

  rmse_css <- rolling_rmse(
    function(y_i, ...) stats::arima(y_i, order = ord,
      seasonal = list(order = ses[1:3], period = ses[4]), method = "CSS-ML"),
    y_train, y_test
  )

  rmse_mle <- if (!is.null(mle_fit))
    rolling_rmse(
      function(y_i, ...) stats::arima(y_i, order = ord,
        seasonal = list(order = ses[1:3], period = ses[4]), method = "ML"),
      y_train, y_test
    ) else NA_real_

  # PMM3 rolling: re-fit PMM3-SARIMA on each expanding window
  rmse_pmm3 <- if (!is.null(pmm3_fit)) {
    rolling_rmse(
      function(y_i, ...) {
        fit <- pmm3_sarima(y_i, order = ord, seasonal = ses)
        # Return as arima object using the refit strategy:
        # Wrap result so predict() works via stats::arima with fixed coefs
        rf <- tryCatch(
          stats::arima(y_i, order = ord,
                       seasonal = list(order = ses[1:3], period = ses[4]),
                       fixed = fit$coefficients, transform.pars = FALSE,
                       method = "CSS"),
          error = function(e) NULL
        )
        if (is.null(rf)) {
          # fallback: plain CSS
          stats::arima(y_i, order = ord,
                       seasonal = list(order = ses[1:3], period = ses[4]),
                       method = "CSS-ML")
        } else rf
      },
      y_train, y_test
    )
  } else NA_real_

  cat(sprintf("  RMSE CSS-ML  = %.5f\n", rmse_css))
  cat(sprintf("  RMSE MLE     = %.5f\n", rmse_mle))
  cat(sprintf("  RMSE PMM3    = %.5f\n\n", rmse_pmm3))

  # ---- 5. Collect results ---------------------------------------------------
  sc <- function(x) as.numeric(x)[1L]   # safe scalar extraction
  all_results[[mdl$name]] <- data.frame(
    model        = rep(mdl$name, 3L),
    method       = c("CSS-ML", "MLE", "PMM3"),
    gamma4_resid = c(sc(cum_css$gamma4),
                     if (!is.null(mle_fit)) sc(cum_mle$gamma4)  else NA_real_,
                     if (!is.null(pmm3_fit)) sc(pmm3_fit$cumulants_final$gamma4) else NA_real_),
    g3_theor     = c(sc(cum_css$g3_theoretical),
                     if (!is.null(mle_fit)) sc(cum_mle$g3_theoretical) else NA_real_,
                     if (!is.null(pmm3_fit)) sc(pmm3_fit$cumulants_final$g3_theoretical) else NA_real_),
    aic          = c(sc(AIC(css_fit)),
                     if (!is.null(mle_fit)) sc(AIC(mle_fit)) else NA_real_,
                     if (!is.null(pmm3_fit)) sc(pmm3_fit$aic) else NA_real_),
    bic          = c(sc(BIC(css_fit)),
                     if (!is.null(mle_fit)) sc(BIC(mle_fit)) else NA_real_,
                     if (!is.null(pmm3_fit)) sc(pmm3_fit$bic) else NA_real_),
    rmse_24      = c(sc(rmse_css), sc(rmse_mle), sc(rmse_pmm3)),
    stringsAsFactors = FALSE
  )
}

# ===========================================================================
# Seasonal regression model (using existing lm_pmm3)
# ===========================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Model: Seasonal regression  y ~ time + sin(2πt) + cos(2πt)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

t_full <- as.numeric(time(co2))
# Seasonal harmonics: one full cycle per year
# MUST use fractional year (t %% 1), NOT 12*t.
# sin(2*pi*12*t) = sin(2*pi * integer) = 0 for all monthly obs → collinear!
dat_seas <- data.frame(
  y    = y_full,
  x    = t_full,
  sin1 = sin(2 * pi * (t_full %% 1)),
  cos1 = cos(2 * pi * (t_full %% 1))
)
formula_seas <- y ~ x + sin1 + cos1

cat("OLS seasonal regression:\n")
ols_seas <- lm(formula_seas, data = dat_seas)
cat(round(coef(ols_seas), 5), "\n")
cum_ols_seas <- .sarima_pmm3_cumulants(residuals(ols_seas))
cat("OLS seasonal residual cumulants:\n"); print_cum(cum_ols_seas)
cat(sprintf("  R² = %.5f\n\n", summary(ols_seas)$r.squared))

cat("PMM3 seasonal regression:\n")
pmm3_seas <- tryCatch(lm_pmm3(formula_seas, dat_seas),
                      error = function(e) { message("PMM3 reg failed: ", e); NULL })
if (!is.null(pmm3_seas)) {
  cat(round(pmm3_seas@coefficients, 5), "\n")
  cum_pmm3_seas <- .sarima_pmm3_cumulants(pmm3_seas@residuals)
  cat("PMM3 seasonal residual cumulants:\n"); print_cum(cum_pmm3_seas)
  cat(sprintf("  Convergence: %s  iter: %d\n\n",
              pmm3_seas@convergence, pmm3_seas@iterations))

  shift_seas <- pmm3_seas@coefficients - coef(ols_seas)
  cat("PMM3 vs OLS coefficient shifts:\n")
  print(round(shift_seas, 6))
  cat("\n")
}

# ===========================================================================
# PMM2 empirical comparison for seasonal regression
# Question: how much does PMM2 actually help when gamma3 = +0.408?
# ===========================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("PMM2 empirical comparison: y ~ x + sin1 + cos1\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Bootstrap (B=2000): variance of x-slope estimator across methods
# bootstrap_comparison(...)[2] extracts coefficient of x (trend slope)
cat("Bootstrap (B=2000, seed=", CFG$seed, ")...\n")
boot_seas <- bootstrap_comparison(dat_seas, formula = formula_seas, B = 2000)

cat("\nBootstrap results (g3_empirical = Var/Var_OLS for x-slope):\n")
var_ols_boot <- boot_seas$sd[boot_seas$method == "OLS"]^2
boot_seas$delta_pct <- (boot_seas$g3_empirical - 1) * 100
print(boot_seas[, c("method", "sd", "g3_empirical", "delta_pct", "na_pct")],
      digits = 5, row.names = FALSE)
cat("\n")

# LOO-MSE: manual loop to avoid PMM2 dispatch bug
cat("LOO-MSE (n=", nrow(dat_seas), ")...\n")
n_seas <- nrow(dat_seas)

loo_ols_seas <- mean(sapply(seq_len(n_seas), function(i) {
  fit <- lm(formula_seas, data = dat_seas[-i, ])
  (dat_seas$y[i] - as.numeric(predict(fit, newdata = dat_seas[i, ])))^2
}))

loo_pmm2_seas <- mean(sapply(seq_len(n_seas), function(i) {
  tryCatch({
    fit <- EstemPMM::lm_pmm2(formula_seas, dat_seas[-i, ])
    pred <- as.numeric(predict(fit, newdata = dat_seas[i, ]))
    (dat_seas$y[i] - pred)^2
  }, error = function(e) NA_real_)
}), na.rm = TRUE)

loo_pmm3_seas <- mean(sapply(seq_len(n_seas), function(i) {
  tryCatch({
    fit <- lm_pmm3(formula_seas, dat_seas[-i, ])
    pred <- as.numeric(predict(fit, newdata = dat_seas[i, ]))
    (dat_seas$y[i] - pred)^2
  }, error = function(e) NA_real_)
}), na.rm = TRUE)

loo_mltn_seas <- mean(sapply(seq_len(n_seas), function(i) {
  tryCatch({
    fit <- mle_tn(formula_seas, dat_seas[-i, ])
    pred <- as.numeric(predict(fit, newdata = dat_seas[i, ]))
    (dat_seas$y[i] - pred)^2
  }, error = function(e) NA_real_)
}), na.rm = TRUE)

loo_seas_df <- data.frame(
  method  = c("OLS", "PMM2", "PMM3", "MLE_TN"),
  loo_mse = c(loo_ols_seas, loo_pmm2_seas, loo_pmm3_seas, loo_mltn_seas),
  stringsAsFactors = FALSE
)
loo_seas_df$rel_ols  <- loo_seas_df$loo_mse / loo_ols_seas
loo_seas_df$delta_pct <- (loo_seas_df$rel_ols - 1) * 100

cat("\nLOO-MSE results:\n")
print(loo_seas_df, digits = 6, row.names = FALSE)
cat("\n")

# Save combined PMM2 seasonal comparison
pmm2_seas_out <- data.frame(
  model    = "seasonal_regression",
  method   = boot_seas$method,
  gamma3   = cum_ols_seas$gamma3,          # same for all (OLS residuals)
  gamma4   = cum_ols_seas$gamma4,
  g3_theor = cum_ols_seas$g3_theoretical,
  g3_emp   = boot_seas$g3_empirical,
  delta_boot_pct = boot_seas$delta_pct,
  loo_mse  = loo_seas_df$loo_mse[match(boot_seas$method, loo_seas_df$method)],
  rel_ols  = loo_seas_df$rel_ols[match(boot_seas$method, loo_seas_df$method)],
  stringsAsFactors = FALSE
)

write.csv(pmm2_seas_out,
          file.path(OUT, "pmm2_seasonal_comparison.csv"),
          row.names = FALSE)
cat("Saved:", file.path(OUT, "pmm2_seasonal_comparison.csv"), "\n\n")

# Human-readable summary
cat("--- PMM2 seasonal regression: summary ---\n")
best_boot_seas <- boot_seas$method[which.min(boot_seas$g3_empirical)]
best_loo_seas  <- loo_seas_df$method[which.min(loo_seas_df$loo_mse)]
cat(sprintf("Bootstrap winner (g3_emp): %s\n", best_boot_seas))
cat(sprintf("LOO-MSE  winner          : %s\n", best_loo_seas))
cat(sprintf("PMM2 bootstrap g3_emp    : %.4f  (delta = %+.1f%%)\n",
            boot_seas$g3_empirical[boot_seas$method == "PMM2"],
            boot_seas$delta_pct[boot_seas$method == "PMM2"]))
cat(sprintf("PMM2 LOO-MSE rel_OLS     : %.5f  (delta = %+.2f%%)\n",
            loo_seas_df$rel_ols[loo_seas_df$method == "PMM2"],
            loo_seas_df$delta_pct[loo_seas_df$method == "PMM2"]))
cat("\n")

# ===========================================================================
# Postulation E: Two-harmonic seasonal regression
#   y ~ x + sin(2πt) + cos(2πt) + sin(4πt) + cos(4πt)
# Hypothesis: second harmonic absorbs asymmetric overtone → γ₃ ≈ 0,
#             remaining residuals may be platykurtic → PMM3 chance
# ===========================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Postulation E: Two-harmonic seasonal regression\n")
cat("  y ~ x + sin(2πt) + cos(2πt) + sin(4πt) + cos(4πt)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

dat_seas2 <- data.frame(
  y    = y_full,
  x    = t_full,
  sin1 = sin(2 * pi * (t_full %% 1)),
  cos1 = cos(2 * pi * (t_full %% 1)),
  sin2 = sin(4 * pi * (t_full %% 1)),
  cos2 = cos(4 * pi * (t_full %% 1))
)
formula_seas2 <- y ~ x + sin1 + cos1 + sin2 + cos2

cat("OLS two-harmonic regression:\n")
ols_seas2 <- lm(formula_seas2, data = dat_seas2)
cat(round(coef(ols_seas2), 5), "\n")
cum_ols_seas2 <- .sarima_pmm3_cumulants(residuals(ols_seas2))
cat("OLS residual cumulants:\n"); print_cum(cum_ols_seas2)
cat(sprintf("  R² = %.5f\n\n", summary(ols_seas2)$r.squared))

# Assess symmetry before running bootstrap
sym_ok2 <- abs(cum_ols_seas2$gamma3) < 0.3
platy_ok2 <- !is.na(cum_ols_seas2$gamma4) && cum_ols_seas2$gamma4 < -0.5

cat(sprintf("  Symmetry (|γ₃|<0.3): %s  |  Platykurtic (γ₄<-0.5): %s\n\n",
            if (sym_ok2) "YES ✓" else sprintf("NO ✗  (γ₃=%+.3f)", cum_ols_seas2$gamma3),
            if (platy_ok2) sprintf("YES ✓  (γ₄=%.3f)", cum_ols_seas2$gamma4) else
                           sprintf("NO ✗  (γ₄=%.3f)", cum_ols_seas2$gamma4)))

# PMM3 fit
cat("PMM3 two-harmonic regression:\n")
pmm3_seas2 <- tryCatch(lm_pmm3(formula_seas2, dat_seas2),
                       error = function(e) { message("PMM3 failed: ", e); NULL })
if (!is.null(pmm3_seas2)) {
  shift2 <- pmm3_seas2@coefficients - coef(ols_seas2)
  cat("PMM3 vs OLS coefficient shifts:\n")
  print(round(shift2, 6))
  cat(sprintf("  Convergence: %s  iter: %d\n\n",
              pmm3_seas2@convergence, pmm3_seas2@iterations))
}

# Bootstrap + LOO
cat("Bootstrap (B=2000, seed=", CFG$seed, ")...\n")
boot_seas2 <- bootstrap_comparison(dat_seas2, formula = formula_seas2, B = 2000)
boot_seas2$delta_pct <- (boot_seas2$g3_empirical - 1) * 100
cat("\nBootstrap results (g3_empirical = Var/Var_OLS for x-slope):\n")
print(boot_seas2[, c("method", "sd", "g3_empirical", "delta_pct", "na_pct")],
      digits = 5, row.names = FALSE)
cat("\n")

cat("LOO-MSE (n=", nrow(dat_seas2), ")...\n")
n_s2 <- nrow(dat_seas2)

loo_ols_s2 <- mean(sapply(seq_len(n_s2), function(i) {
  fit <- lm(formula_seas2, data = dat_seas2[-i, ])
  (dat_seas2$y[i] - as.numeric(predict(fit, newdata = dat_seas2[i, ])))^2
}))
loo_pmm2_s2 <- mean(sapply(seq_len(n_s2), function(i) {
  tryCatch({
    fit <- EstemPMM::lm_pmm2(formula_seas2, dat_seas2[-i, ])
    (dat_seas2$y[i] - as.numeric(predict(fit, newdata = dat_seas2[i, ])))^2
  }, error = function(e) NA_real_)
}), na.rm = TRUE)
loo_pmm3_s2 <- mean(sapply(seq_len(n_s2), function(i) {
  tryCatch({
    fit <- lm_pmm3(formula_seas2, dat_seas2[-i, ])
    (dat_seas2$y[i] - as.numeric(predict(fit, newdata = dat_seas2[i, ])))^2
  }, error = function(e) NA_real_)
}), na.rm = TRUE)
loo_mltn_s2 <- mean(sapply(seq_len(n_s2), function(i) {
  tryCatch({
    fit <- mle_tn(formula_seas2, dat_seas2[-i, ])
    (dat_seas2$y[i] - as.numeric(predict(fit, newdata = dat_seas2[i, ])))^2
  }, error = function(e) NA_real_)
}), na.rm = TRUE)

loo_s2_df <- data.frame(
  method    = c("OLS", "PMM2", "PMM3", "MLE_TN"),
  loo_mse   = c(loo_ols_s2, loo_pmm2_s2, loo_pmm3_s2, loo_mltn_s2),
  stringsAsFactors = FALSE
)
loo_s2_df$rel_ols   <- loo_s2_df$loo_mse / loo_ols_s2
loo_s2_df$delta_pct <- (loo_s2_df$rel_ols - 1) * 100

cat("\nLOO-MSE results:\n")
print(loo_s2_df, digits = 6, row.names = FALSE)
cat("\n")

# Save
seas2_out <- data.frame(
  model    = "seasonal_regression_2harm",
  method   = boot_seas2$method,
  gamma3   = cum_ols_seas2$gamma3,
  gamma4   = cum_ols_seas2$gamma4,
  g3_theor = cum_ols_seas2$g3_theoretical,
  g3_emp   = boot_seas2$g3_empirical,
  delta_boot_pct = boot_seas2$delta_pct,
  loo_mse  = loo_s2_df$loo_mse[match(boot_seas2$method, loo_s2_df$method)],
  rel_ols  = loo_s2_df$rel_ols[match(boot_seas2$method, loo_s2_df$method)],
  stringsAsFactors = FALSE
)
write.csv(seas2_out,
          file.path(OUT, "seasonal_2harm_comparison.csv"),
          row.names = FALSE)
cat("Saved:", file.path(OUT, "seasonal_2harm_comparison.csv"), "\n\n")

cat("--- Two-harmonic summary ---\n")
cat(sprintf("Bootstrap winner (g3_emp): %s\n",
            boot_seas2$method[which.min(boot_seas2$g3_empirical)]))
cat(sprintf("LOO-MSE  winner          : %s\n",
            loo_s2_df$method[which.min(loo_s2_df$loo_mse)]))
for (m in c("OLS", "PMM2", "PMM3", "MLE_TN")) {
  cat(sprintf("  %-8s  g3_emp=%6.4f (%+.1f%%)  LOO_rel=%.5f (%+.2f%%)\n",
              m,
              boot_seas2$g3_empirical[boot_seas2$method == m],
              boot_seas2$delta_pct[boot_seas2$method == m],
              loo_s2_df$rel_ols[loo_s2_df$method == m],
              loo_s2_df$delta_pct[loo_s2_df$method == m]))
}
cat("\n")

# ===========================================================================
# Save summary table
# ===========================================================================
summary_df <- do.call(rbind, all_results)
rownames(summary_df) <- NULL

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("SUMMARY TABLE\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
print(summary_df, digits = 5)

write.csv(summary_df, file.path(OUT, "sarima_comparison.csv"), row.names = FALSE)
cat("\nSaved:", file.path(OUT, "sarima_comparison.csv"), "\n")

# Save seasonal regression comparison
if (!is.null(pmm3_seas)) {
  seas_df <- data.frame(
    model    = "seasonal_regression",
    method   = c("OLS", "PMM3"),
    gamma4   = c(cum_ols_seas$gamma4, cum_pmm3_seas$gamma4),
    g3_theor = c(cum_ols_seas$g3_theoretical, cum_pmm3_seas$g3_theoretical),
    stringsAsFactors = FALSE
  )
  write.csv(seas_df, file.path(OUT, "seasonal_regression_comparison.csv"),
            row.names = FALSE)
  cat("Saved:", file.path(OUT, "seasonal_regression_comparison.csv"), "\n")
}

cat("\n╔══════════════════════════════════════════╗\n")
cat("║  Analysis complete                       ║\n")
cat(sprintf("║  Output: %-31s ║\n", file.path(OUT)))
cat("╚══════════════════════════════════════════╝\n")
