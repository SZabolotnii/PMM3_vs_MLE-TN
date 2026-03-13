# =============================================================================
# R/04_real_data.R
# Real data analysis — Plan Blocks 3 & 4
#
# Dataset A: Biaxial Fatigue (Rieck & Nedelman, Technometrics 1991, 33:51-60)
#   n=46, 1% Cr-Mo-V steel
#   x = log(work per cycle),  y = log(cycles to failure)
#   Status: manually restored; NOT available in any CRAN package
#   ⚠ OLS residuals: gamma3 ≈ -0.47 (asymmetric!) → PMM2 may outperform PMM3
#
# Dataset B: nlme::Fatigue  (CRAN, 262 obs, repeated-measures, crack growth)
#   Used as accessible fallback for method comparison
# =============================================================================

# ---------------------------------------------------------------------------
# Dataset A: Biaxial Fatigue
# ---------------------------------------------------------------------------

#' Load biaxial fatigue dataset (Rieck & Nedelman, 1991)
#'
#' Manually restored from the original Technometrics paper.
#' 46 observations, 6 stress levels (8 reps each, minus 2 missing).
#'
#' @return data.frame with columns: x (log work/cycle), y (log cycles), dataset
load_biaxial_fatigue <- function() {
  x <- c(
    -1.66, -1.66, -1.66, -1.66, -1.66, -1.66, -1.66, -1.66,
    -1.39, -1.39, -1.39, -1.39, -1.39, -1.39, -1.39, -1.39,
    -1.11, -1.11, -1.11, -1.11, -1.11, -1.11, -1.11, -1.11,
    -0.85, -0.85, -0.85, -0.85, -0.85, -0.85, -0.85, -0.85,
    -0.60, -0.60, -0.60, -0.60, -0.60, -0.60, -0.60, -0.60,
    -0.27, -0.27, -0.27, -0.27, -0.27, -0.27
  )
  y <- c(
    6.00, 5.90, 5.89, 5.59, 5.78, 5.93, 5.74, 5.30,
    5.58, 5.45, 5.26, 5.76, 5.85, 5.56, 5.90, 5.30,
    5.41, 5.32, 5.23, 5.12, 5.02, 5.31, 5.60, 5.38,
    5.23, 5.27, 5.22, 4.76, 5.00, 5.10, 4.97, 4.78,
    5.03, 4.98, 4.87, 4.85, 4.88, 5.08, 5.03, 5.02,
    4.55, 4.70, 4.60, 4.50, 4.75, 4.80
  )
  data.frame(
    x       = x,
    y       = y,
    dataset = "Biaxial Fatigue (Rieck & Nedelman, 1991)",
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Dataset B: nlme::Fatigue (CRAN fallback)
# ---------------------------------------------------------------------------

#' Load nlme::Fatigue aggregated to simple linear regression
#'
#' Original: 262 obs, 21 metal paths, repeated crack-length measurements.
#' Aggregated: mean(relLength) ~ cycles across paths.
load_nlme_fatigue <- function() {
  if (!requireNamespace("nlme", quietly = TRUE))
    stop("Package 'nlme' required. Run: install.packages('nlme')")
  e <- new.env()
  data("Fatigue", package = "nlme", envir = e)
  agg <- aggregate(relLength ~ cycles, e$Fatigue, mean)
  names(agg) <- c("x", "y")
  agg$dataset <- "nlme::Fatigue (mean across paths)"
  agg
}

# ---------------------------------------------------------------------------
# OLS residual diagnostics — determines which method to recommend
# ---------------------------------------------------------------------------

#' Compute OLS residual diagnostics: gamma3, gamma4, gamma6, g3_theoretical
#'
#' @param data    data.frame with columns y (response) and x (predictor)
#' @param formula R formula (default y ~ x)
#' @return Named list with diagnostics and recommended estimation method
resid_diagnostics <- function(data, formula = y ~ x) {
  r   <- residuals(lm(formula, data = data))
  n   <- length(r)
  mu2 <- mean(r^2); mu3 <- mean(r^3)
  mu4 <- mean(r^4); mu6 <- mean(r^6)

  gamma3 <- mu3 / mu2^1.5
  gamma4 <- mu4 / mu2^2 - 3
  gamma6 <- mu6 / mu2^3 - 15 * (mu4 / mu2^2) + 30

  denom    <- 6 + 9 * gamma4 + gamma6
  g3_theor <- if (!is.na(denom) && denom > 0) 1 - gamma4^2 / denom else NA_real_

  method_rec <- if (abs(gamma3) >= CFG$symmetry_threshold)
                  "PMM2  (asymmetric errors: |gamma3| >= threshold)"
                else if (abs(gamma4) < CFG$gaussian_threshold)
                  "OLS   (near-Gaussian errors: |gamma4| < threshold)"
                else
                  "PMM3  (symmetric non-Gaussian errors)"

  list(
    n                  = n,
    gamma3             = gamma3,
    gamma4             = gamma4,
    gamma6             = gamma6,
    g3_theoretical     = g3_theor,
    improvement_pct    = if (!is.na(g3_theor)) (1 - g3_theor) * 100 else NA_real_,
    is_symmetric       = abs(gamma3) < CFG$symmetry_threshold,
    is_gaussian        = abs(gamma4) < CFG$gaussian_threshold,
    recommended_method = method_rec
  )
}

# ---------------------------------------------------------------------------
# Bootstrap comparison of all four methods
# ---------------------------------------------------------------------------

#' Bootstrap comparison: OLS / PMM2 / PMM3 / MLE-TN
#'
#' @param data     data.frame with y, x
#' @param formula  R formula
#' @param B        Bootstrap replications (default 2000 for speed, use 10000 for paper)
#' @param seed     Random seed
#' @return data.frame: method, mean, sd, q025, q975, na_pct, g3_empirical
bootstrap_comparison <- function(data, formula = y ~ x,
                                 B = 2000, seed = CFG$seed) {
  set.seed(seed)
  n <- nrow(data)

  methods <- list(
    OLS    = function(d) coef(lm(formula, d))[2],
    PMM2   = function(d) EstemPMM::lm_pmm2(formula, d)@coefficients[2],
    PMM3   = function(d) lm_pmm3(formula, d)@coefficients[2],
    MLE_TN = function(d) mle_tn(formula, d)$coefficients[2]
  )

  results <- lapply(names(methods), function(m) {
    message("  Bootstrap | ", sprintf("%-8s", m), " B=", B)
    fn <- methods[[m]]
    b1 <- replicate(B, {
      d <- data[sample(n, n, TRUE), ]
      tryCatch(fn(d), error = function(e) NA_real_)
    })
    data.frame(
      method = m,
      mean   = mean(b1, na.rm = TRUE),
      sd     = sd(b1,   na.rm = TRUE),
      q025   = quantile(b1, 0.025, na.rm = TRUE),
      q975   = quantile(b1, 0.975, na.rm = TRUE),
      na_pct = mean(is.na(b1)) * 100,
      stringsAsFactors = FALSE
    )
  })

  df       <- do.call(rbind, results)
  var_ols  <- df$sd[df$method == "OLS"]^2
  df$g3_empirical <- df$sd^2 / var_ols
  df
}

# ---------------------------------------------------------------------------
# Leave-One-Out cross-validation
# Fair criterion for comparing MLE-TN (parametric) and PMM3 (moment-based)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Dataset C: Mauna Loa CO2 (built-in R dataset)
# ---------------------------------------------------------------------------

#' Load Mauna Loa CO2 as regression dataset: CO2 concentration ~ decimal year
#'
#' Built-in R dataset `co2`: 468 monthly observations, 1959–1997.
#' Response y = atmospheric CO2 (ppm); predictor x = decimal year.
#' OLS of y ~ x captures the linear trend; residuals retain the seasonal
#' oscillation (≈ sinusoid, ±3 ppm), which produces a U-shaped (bimodal-edge)
#' distribution → gamma4 << 0, gamma3 ≈ 0.  Ideal PMM3 scenario.
#'
#' @return data.frame with columns y, x, dataset
load_co2 <- function() {
  data.frame(
    y       = as.numeric(co2),
    x       = as.numeric(time(co2)),
    dataset = "Mauna Loa CO2 (Keeling & Whorf, SIO)",
    stringsAsFactors = FALSE
  )
}

#' LOO cross-validation MSE for all four methods
#'
#' @param data     data.frame with y, x
#' @param formula  R formula
#' @return Named vector: LOO-MSE per method
loo_mse <- function(data, formula = y ~ x) {
  n <- nrow(data)
  methods <- list(
    OLS    = function(d) lm(formula, d),
    PMM2   = function(d) EstemPMM::lm_pmm2(formula, d),
    PMM3   = function(d) lm_pmm3(formula, d),
    MLE_TN = function(d) mle_tn(formula, d)
  )

  message("LOO-CV (n=", n, ")...")
  sapply(names(methods), function(m) {
    fn <- methods[[m]]
    mean(sapply(seq_len(n), function(i) {
      tryCatch({
        fit <- fn(data[-i, ])
        (data$y[i] - predict(fit, newdata = data[i, ]))^2
      }, error = function(e) NA_real_)
    }), na.rm = TRUE)
  })
}
