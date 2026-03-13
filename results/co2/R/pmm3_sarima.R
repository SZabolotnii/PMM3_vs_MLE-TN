# =============================================================================
# results/co2/R/pmm3_sarima.R
# PMM3 estimation for ARIMA/SARIMA models  (4th-cumulant based)
#
# Mathematical basis:
#   Model (after differencing): z_t = X_t' β + ε_t
#     X_t — pseudo-regressors (AR: lagged z; MA: lagged CSS residuals)
#   PMM3 score  (S=3 polynomial, same as pmm3_skeleton.R):
#     Z_p = Σ_t ε_t (κ - ε_t²) x_{tp} = 0,   p = 1,..,k
#   Jacobian:
#     J_{pq} = Σ_t (3ε_t² - κ) x_{tp} x_{tq}   →   J = X' diag(3ε²-κ) X
#   Update:  β ← β - J⁻¹ Z
#   κ computed once from CSS residuals (EstemPMM-style fixed moments)
#
# Conditions for PMM3 efficiency gain:
#   γ₃ ≈ 0 (symmetric innovations)  AND  γ₄ < −0.7  (platykurtic)
# =============================================================================

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Compute PMM3 kappa from residuals
#' kappa = (m6 - 3*m4*m2) / (m4 - 3*m2^2);  NA → near-Gaussian fallback
.sarima_pmm3_kappa <- function(eps) {
  e <- eps[is.finite(eps)]
  m2 <- mean(e^2); m4 <- mean(e^4); m6 <- mean(e^6)
  denom <- m4 - 3 * m2^2
  if (abs(denom) < .Machine$double.eps * 1e6) return(NA_real_)
  (m6 - 3 * m4 * m2) / denom
}

#' Compute standardised cumulants γ₃, γ₄, γ₆ and g3_theoretical
.sarima_pmm3_cumulants <- function(eps) {
  e  <- eps[is.finite(eps)]
  m2 <- mean(e^2); m3 <- mean(e^3); m4 <- mean(e^4); m6 <- mean(e^6)
  g3 <- m3 / m2^1.5
  g4 <- m4 / m2^2 - 3
  g6 <- m6 / m2^3 - 15 * (m4 / m2^2) + 30
  denom  <- 6 + 9 * g4 + g6
  g3_th  <- if (is.finite(denom) && denom > 0) 1 - g4^2 / denom else NA_real_
  list(gamma3 = g3, gamma4 = g4, gamma6 = g6,
       g3_theoretical = g3_th,
       improvement_pct = if (!is.na(g3_th)) (1 - g3_th) * 100 else NA_real_)
}

#' Build SARIMA design matrix from differenced series and CSS residuals
#'
#' For AR(p)/SAR(P): columns = lagged differenced series z
#' For MA(q)/SMA(Q): columns = lagged CSS residuals (fixed — EstemPMM style)
#' Seasonal lags at multiples of s.
.sarima_design_matrix <- function(z, eps_css, p, q, P, Q, s) {
  n   <- length(z)
  cols <- list()
  nms  <- character(0)

  for (i in seq_len(p)) {
    cols[[length(cols) + 1L]] <- c(rep(NA_real_, i), z[seq_len(n - i)])
    nms <- c(nms, paste0("ar", i))
  }
  for (i in seq_len(P)) {
    lag <- i * s
    cols[[length(cols) + 1L]] <- c(rep(NA_real_, lag), z[seq_len(n - lag)])
    nms <- c(nms, paste0("sar", i))
  }
  for (i in seq_len(q)) {
    cols[[length(cols) + 1L]] <- c(rep(NA_real_, i), eps_css[seq_len(n - i)])
    nms <- c(nms, paste0("ma", i))
  }
  for (i in seq_len(Q)) {
    lag <- i * s
    cols[[length(cols) + 1L]] <- c(rep(NA_real_, lag), eps_css[seq_len(n - lag)])
    nms <- c(nms, paste0("sma", i))
  }

  if (length(cols) == 0L)
    stop("PMM3-SARIMA: all orders are 0 — nothing to estimate.")

  X <- do.call(cbind, cols)
  colnames(X) <- nms
  X
}

# ---------------------------------------------------------------------------
# Main estimator
# ---------------------------------------------------------------------------

#' Fit ARIMA/SARIMA by PMM3 (Method of Polynomial Moments, S=3)
#'
#' Implements PMM3 for time-series models. Particularly useful when
#' innovations are *symmetric* and *platykurtic* (γ₄ < −0.7).
#'
#' Algorithm:
#'   1. CSS-ML initialisation  →  b_css, ε_css
#'   2. κ = PMM3 kappa from ε_css  (fixed — no re-estimation)
#'   3. Build design matrix X from z (AR lags) and ε_css (MA lags)
#'   4. Newton-Raphson:  b ← b − J⁻¹ Z
#'   5. Refit via stats::arima(fixed = b_pmm3) to recover proper residuals
#'
#' @param y       Numeric time series (raw, not differenced)
#' @param order   c(p, d, q)
#' @param seasonal c(P, D, Q, s); use c(0,0,0,0) for non-seasonal ARIMA
#' @param maxit   Max Newton-Raphson iterations (default 50)
#' @param tol     Convergence tolerance on ||Δβ||₂ (default 1e-7)
#' @param step_max  Max step-size (default 1.0; guards against divergence)
#' @param trace   Print iteration log (default FALSE)
#'
#' @return Named list:
#'   $method, $order, $seasonal, $coefficients, $coef_names,
#'   $sigma2, $residuals, $fitted,
#'   $kappa, $cumulants_css, $cumulants_final,
#'   $aic, $bic, $loglik,
#'   $convergence (logical), $iterations
pmm3_sarima <- function(y,
                        order    = c(0L, 1L, 1L),
                        seasonal = c(0L, 1L, 1L, 12L),
                        maxit    = 50L,
                        tol      = 1e-7,
                        step_max = 1.0,
                        trace    = FALSE) {

  p  <- as.integer(order[1])
  d  <- as.integer(order[2])
  q  <- as.integer(order[3])
  P  <- as.integer(seasonal[1])
  D  <- as.integer(seasonal[2])
  Q  <- as.integer(seasonal[3])
  s  <- as.integer(seasonal[4])
  use_seasonal <- (s > 0L) && (P > 0L || Q > 0L)

  # ---- Step 1: CSS-ML initialisation ----------------------------------------
  css_args <- list(y, order = c(p, d, q), method = "CSS-ML",
                   include.mean = (d == 0L && (!use_seasonal || D == 0L)))
  if (use_seasonal)
    css_args$seasonal <- list(order = c(P, D, Q), period = s)

  css_fit <- tryCatch(
    do.call(stats::arima, css_args),
    error = function(e) {
      css_args$method <- "CSS"
      tryCatch(do.call(stats::arima, css_args),
               error = function(e2) stop("CSS init failed: ", conditionMessage(e2)))
    }
  )

  b_css   <- as.numeric(stats::coef(css_fit))
  b_names <- names(stats::coef(css_fit))
  eps_css <- as.numeric(stats::residuals(css_fit))
  eps_css[!is.finite(eps_css)] <- 0.0          # burn-in NAs → 0

  # ---- Step 2: κ from CSS residuals -----------------------------------------
  # Use only valid (non-burn-in) residuals
  n_burnin  <- max(p + P * s, q + Q * s, d + D * s, 1L)
  eps_valid <- eps_css[seq(n_burnin + 1L, length(eps_css))]

  kappa <- .sarima_pmm3_kappa(eps_valid)
  cum_css <- .sarima_pmm3_cumulants(eps_valid)

  if (trace) {
    cat(sprintf("[PMM3-SARIMA] CSS init | γ₃=%.3f | γ₄=%.3f | κ=%s\n",
                cum_css$gamma3, cum_css$gamma4,
                if (is.na(kappa)) "NA (near-Gaussian)" else sprintf("%.4f", kappa)))
  }

  # Near-Gaussian fallback — PMM3 gives no gain; return CSS
  if (is.na(kappa)) {
    message("[PMM3-SARIMA] Near-Gaussian innovations. Returning CSS estimates.")
    return(.pmm3_sarima_result(css_fit, b_css, b_names, kappa, cum_css,
                               converged = TRUE, iterations = 0L,
                               note = "near-Gaussian: CSS returned"))
  }

  # ---- Step 3: Design matrix ------------------------------------------------
  # Differenced series z (matches the internal ARIMA state)
  z <- as.numeric(y)
  if (d  > 0L) z <- diff(z, differences = d)
  if (use_seasonal && D > 0L) z <- diff(z, lag = s, differences = D)

  nz <- length(z)
  eps_ma <- tail(eps_css, nz)   # CSS residuals trimmed to match z length

  X_full <- .sarima_design_matrix(z, eps_ma, p, q, P, Q, s)
  valid   <- complete.cases(X_full) & is.finite(z)
  X_v     <- X_full[valid, , drop = FALSE]
  y_v     <- z[valid]

  k <- ncol(X_v)
  if (nrow(X_v) < k + 5L)
    stop("[PMM3-SARIMA] Too few valid rows (", nrow(X_v), ") for k=", k, " parameters.")

  # ---- Step 4: Newton-Raphson -----------------------------------------------
  b <- b_css
  converged <- FALSE
  iter <- 0L

  for (i in seq_len(maxit)) {
    iter <- as.integer(i)
    eps  <- as.numeric(y_v - X_v %*% b)

    # Score: Z = X' (ε (κ − ε²))  — force numeric vector (not matrix)
    Z1 <- eps * (kappa - eps^2)
    Z  <- as.numeric(crossprod(X_v, Z1))   # length-k vector

    # Jacobian: J = X' diag(3ε² − κ) X
    w  <- 3.0 * eps^2 - kappa
    J  <- crossprod(X_v * w, X_v)          # k × k matrix

    # Regularise if near-singular
    if (rcond(J) < 1e-12) J <- J + diag(1e-8, k)

    delta <- tryCatch(as.numeric(solve(J, Z)),
                      error = function(e) rep(NA_real_, k))
    if (anyNA(delta)) { warning("[PMM3-SARIMA] Singular Jacobian — stopping."); break }

    # Step-size limiting
    step_norm <- sqrt(sum(delta^2))
    if (step_norm > step_max) delta <- delta * step_max / step_norm

    b <- b - delta
    if (trace)
      cat(sprintf("  iter %2d | ||Δβ|| = %.3e\n", i, step_norm))

    if (step_norm < tol) { converged <- TRUE; break }
  }

  # Divergence guard (same as pmm3_skeleton)
  if (!all(is.finite(b)) || sqrt(sum((b - b_css)^2)) > 10.0) {
    warning("[PMM3-SARIMA] Diverged — reverting to CSS estimates.")
    b <- b_css;  converged <- FALSE
  }

  # ---- Step 5: Refit with fixed PMM3 coefficients ---------------------------
  .pmm3_sarima_result(css_fit, b, b_names, kappa, cum_css,
                      y = y, use_seasonal = use_seasonal,
                      order = c(p, d, q), seas_order = c(P, D, Q), s = s,
                      converged = converged, iterations = iter)
}

# ---------------------------------------------------------------------------
# Internal: build return object; optionally refit with stats::arima(fixed=)
# ---------------------------------------------------------------------------
.pmm3_sarima_result <- function(css_fit, b, b_names, kappa, cum_css,
                                 y = NULL, use_seasonal = FALSE,
                                 order = NULL, seas_order = NULL, s = NULL,
                                 converged, iterations, note = NULL) {
  names(b) <- b_names

  refit <- NULL
  if (!is.null(y)) {
    refit <- tryCatch({
      rf_args <- list(y, order = order, fixed = b, transform.pars = FALSE,
                      method = "CSS",
                      include.mean = (order[2] == 0L && (!use_seasonal || seas_order[2] == 0L)))
      if (use_seasonal)
        rf_args$seasonal <- list(order = seas_order, period = s)
      do.call(stats::arima, rf_args)
    }, error = function(e) NULL)
  }

  if (is.null(refit)) refit <- css_fit   # fallback

  res_final  <- as.numeric(stats::residuals(refit))
  cum_final  <- .sarima_pmm3_cumulants(res_final[is.finite(res_final)])

  list(
    method          = if (is.null(note)) "PMM3-SARIMA" else paste0("PMM3-SARIMA(", note, ")"),
    order           = if (!is.null(order)) order else css_fit$arma[1:3],
    seasonal        = if (!is.null(seas_order)) c(seas_order, s) else NULL,
    coefficients    = b,
    coef_names      = b_names,
    sigma2          = refit$sigma2,
    residuals       = res_final,
    fitted          = as.numeric(y %||% (refit$x)) - res_final,
    kappa           = kappa,
    cumulants_css   = cum_css,
    cumulants_final = cum_final,
    loglik          = as.numeric(stats::logLik(refit)),
    aic             = stats::AIC(refit),
    bic             = stats::BIC(refit),
    convergence     = converged,
    iterations      = iterations
  )
}

# Minimal NULL-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b
