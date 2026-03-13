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
