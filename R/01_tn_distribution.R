# =============================================================================
# R/01_tn_distribution.R
# Two-piece Normal (TN) distribution: rtn, dtn, ptn, moments, Table T1
# Reference: Salinas et al. (2023), Mathematics 11(5), 1271
#
# Exported functions:
#   rtn()             — генератор TN(xi, eta, lambda)
#   dtn()             — щільність
#   ptn()             — функція розподілу
#   tn_raw_moment_2k()— E(Z^{2k}) числовим інтегруванням
#   tn_cumulants()    — gamma4, gamma6, g3 для заданого lambda
#   build_table_T1()  — Таблиця T1 для сітки lambda
# =============================================================================

#' Generate random TN(xi, eta, lambda) variates
#'
#' Algorithm: fold-and-sign (Section 4.2, Salinas et al. 2023)
#'   U ~ N(lambda, 1)  →  Y = |U| ~ FN(lambda, 1)
#'   S ~ Bernoulli(1/2) →  Z = S*Y ~ TN(lambda)
#'   X = xi + eta * Z
#'
#' @param n      Number of observations
#' @param xi     Location parameter
#' @param eta    Scale parameter (> 0)
#' @param lambda Shape parameter (>= 0); lambda=0 gives normal
rtn <- function(n, xi = 0, eta = 1, lambda = 1) {
  stopifnot(n > 0, eta > 0, lambda >= 0)
  U <- rnorm(n, mean = lambda, sd = 1)
  Y <- abs(U)
  S <- sample(c(-1L, 1L), n, replace = TRUE)
  xi + eta * (S * Y)
}

#' TN density f(x; xi, eta, lambda)
#'
#' Formula (7) in Salinas et al. (2023):
#'   f(x) = phi(z) * exp(-lambda^2/2) * cosh(lambda*z) / eta
#' where z = (x - xi) / eta
#'
#' Uses numerically stable log_cosh(u) = |u| + log(1 + exp(-2|u|)) - log(2)
#' to avoid overflow in cosh() for large |lambda*z|.
dtn <- function(x, xi = 0, eta = 1, lambda = 0, log = FALSE) {
  stopifnot(eta > 0, lambda >= 0)
  z   <- (x - xi) / eta
  lz  <- lambda * z
  # Numerically stable: log(cosh(u)) = |u| + log1p(exp(-2*|u|)) - log(2)
  log_cosh <- abs(lz) + log1p(exp(-2 * abs(lz))) - log(2)
  log_dens <- dnorm(z, log = TRUE) - log(eta) - lambda^2 / 2 + log_cosh
  if (log) log_dens else exp(log_dens)
}

#' TN distribution function P(X <= q)
#'
#' Formula (3) in Salinas et al.:
#'   F(x) = [Phi(z - lambda) + Phi(z + lambda)] / 2
ptn <- function(q, xi = 0, eta = 1, lambda = 0) {
  z <- (q - xi) / eta
  (pnorm(z - lambda) + pnorm(z + lambda)) / 2
}

#' E(Z^{2k}) for Z ~ TN(0, 1, lambda) via numerical integration
#'
#' Note: for k=1,2,3 this gives mu2, mu4, mu6 needed for PMM3 coefficients.
#' Numerical integration is more reliable than the 1F1 hypergeometric series
#' for large lambda.
tn_raw_moment_2k <- function(k, lambda) {
  tryCatch(
    integrate(
      function(z) z^(2 * k) * dtn(z, lambda = lambda),
      lower = -Inf, upper = Inf,
      subdivisions = 1000L, rel.tol = 1e-8
    )$value,
    error = function(e) NA_real_
  )
}

#' Compute TN cumulant coefficients and PMM3 theoretical efficiency
#'
#' Returns list:
#'   lambda, mu2, mu4, mu6        — raw central moments
#'   gamma4                        — excess kurtosis = mu4/mu2^2 - 3
#'   gamma6                        — standardised 6th cumulant
#'   g3                            — Var(PMM3)/Var(OLS) = 1 - gamma4^2/(6+9*gamma4+gamma6)
#'   improvement_pct               — (1-g3)*100 %
#'   is_bimodal                    — lambda > 1
#'
#' Key property of TN: gamma4(lambda) = -3*lambda^4/(1+lambda^2)^2 <= 0
#' → TN is always platykurtic → PMM3 always gains over OLS (g3 < 1) when lambda > 0
tn_cumulants <- function(lambda) {
  mu2 <- tn_raw_moment_2k(1, lambda)
  mu4 <- tn_raw_moment_2k(2, lambda)
  mu6 <- tn_raw_moment_2k(3, lambda)

  gamma4 <- mu4 / mu2^2 - 3
  gamma6 <- mu6 / mu2^3 - 15 * (mu4 / mu2^2) + 30

  # Theoretical PMM3 variance reduction (formula 12 from PMM3 theory)
  denom <- 6 + 9 * gamma4 + gamma6
  g3    <- if (!is.na(denom) && denom > 0) 1 - gamma4^2 / denom else NA_real_

  list(
    lambda          = lambda,
    mu2             = mu2,
    mu4             = mu4,
    mu6             = mu6,
    gamma4          = gamma4,
    gamma6          = gamma6,
    g3              = g3,
    improvement_pct = if (!is.na(g3)) (1 - g3) * 100 else NA_real_,
    is_bimodal      = lambda > 1   # unimodal for lambda <= 1, bimodal for lambda > 1
  )
}

#' Build Table T1: theoretical TN cumulants over a grid of lambda values
#'
#' This is the reference table for all Monte Carlo comparisons.
#' g3_theoretical serves as the benchmark for empirical g3 from simulations.
build_table_T1 <- function(lambda_grid = CFG$lambda_grid) {
  message("Building Table T1 (", length(lambda_grid), " lambda values)...")
  rows <- lapply(lambda_grid, tn_cumulants)
  T1   <- do.call(rbind, lapply(rows, as.data.frame))
  rownames(T1) <- NULL
  T1$mode <- ifelse(T1$is_bimodal, "bimodal", "unimodal")
  message("Table T1 complete.")
  T1
}
