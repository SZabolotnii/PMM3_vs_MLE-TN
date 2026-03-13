# =============================================================================
# R/pmm3_skeleton.R
# PMM3 — Method of Polynomial Moments (S=3) for symmetric error distributions
#
# Mathematical basis (Zabolotnyi et al., 2019; Kunchenko, 2002):
#   Model: y = X*beta + eps,  eps ~ symmetric distribution with E(eps)=0
#
# Score equations (estimating equations for S=3):
#   sum_v { eps_v * (kappa - eps_v^2) } * x_v^p  = 0,  p = 0,..,P-1
#   where eps_v = y_v - X[v,] %*% beta  (residuals)
#         kappa = (m6 - 3*m4*m2) / (m4 - 3*m2^2)  (moment ratio, fixed at OLS)
#
# Newton-Raphson Jacobian:
#   d(Z1_v)/d(beta) = -(3*eps_v^2 - kappa) * X[v,]
#   JZs[p,q] = sum_v (3*eps_v^2 - kappa) * x_v^p * x_v^q  = t(X*w) %*% X
#   where w_v = 3*eps_v^2 - kappa  (observation weights)
#
# Theoretical efficiency (Zabolotnyi 2019, eq. 14):
#   g3 = Var(PMM3) / Var(OLS) = 1 - gamma4^2 / (6 + 9*gamma4 + gamma6)
#   where gamma4 = m4/m2^2 - 3,  gamma6 = m6/m2^3 - 15*(m4/m2^2) + 30
#
# Reference implementation: EPF/solvMRs3ex_ref.R
# =============================================================================

STUB <- FALSE   # real Newton-Raphson implemented — MC results will use PMM3

# ---------------------------------------------------------------------------
# Define PMM3fit S4 class (only if not already loaded from EstemPMM)
# ---------------------------------------------------------------------------
if (!isClass("PMM3fit")) {
  setClass("PMM3fit", representation(
    coefficients  = "numeric",
    fitted        = "numeric",
    residuals     = "numeric",
    convergence   = "logical",
    iterations    = "integer",
    call          = "call",
    terms         = "ANY",
    method        = "character",
    m2            = "numeric",
    m4            = "numeric",
    m6            = "numeric",
    gamma4        = "numeric",
    gamma6        = "numeric",
    g_coefficient = "numeric"
  ))
}

# ---------------------------------------------------------------------------
# Internal utility: sample moments of (centred) residuals
# ---------------------------------------------------------------------------
.compute_moments_pmm3 <- function(residuals) {
  r <- residuals - mean(residuals)
  list(mu2 = mean(r^2), mu4 = mean(r^4), mu6 = mean(r^6))
}

# ---------------------------------------------------------------------------
# Internal utility: theoretical g3 efficiency coefficient from raw moments
# ---------------------------------------------------------------------------
.pmm3_g_coefficient <- function(mu2, mu4, mu6) {
  gamma4 <- mu4 / mu2^2 - 3
  gamma6 <- mu6 / mu2^3 - 15 * (mu4 / mu2^2) + 30
  denom  <- 6 + 9 * gamma4 + gamma6
  if (is.na(denom) || denom <= 0) return(1.0)
  max(0, 1 - gamma4^2 / denom)
}

# ---------------------------------------------------------------------------
# Internal utility: warn if residuals appear asymmetric
# ---------------------------------------------------------------------------
.check_symmetry_pmm3 <- function(residuals) {
  mu2    <- mean(residuals^2)
  gamma3 <- mean(residuals^3) / mu2^1.5
  if (abs(gamma3) > CFG$symmetry_threshold)
    warning(sprintf(
      "[PMM3] gamma3 = %.3f (|gamma3| > %.1f). Errors appear asymmetric. Consider PMM2.",
      gamma3, CFG$symmetry_threshold
    ))
  invisible(gamma3)
}

# ---------------------------------------------------------------------------
# Core Newton-Raphson solver (vectorised; based on EPF/solvMRs3ex_ref.R)
#
# Score:    Z1_v  = eps_v * (kappa - eps_v^2)
#           Z_p   = sum_v Z1_v * X[v,p]      (length P vector)
# Jacobian: w_v   = 3*eps_v^2 - kappa
#           JZs   = t(X * w) %*% X            (P x P matrix, = X'diag(w)X)
# Update:   B    -= solve(JZs, Z)
#
# @param B_ols    OLS starting values (length P)
# @param X        design matrix (n x P)
# @param Y        response vector (length n)
# @param m2,m4,m6 2nd/4th/6th moments of OLS residuals (scalars)
# @param adaptive  logical: if TRUE, re-estimate kappa from current residuals
#                  at each NR step (unknown distribution scenario)
# @param tol      convergence tolerance on ||Q||_2
# @param max_iter maximum N-R iterations
# @return list(B, converged, iter, kappa, near_gaussian)
# ---------------------------------------------------------------------------
.pmm3_nr_solver <- function(B_ols, X, Y, m2, m4, m6,
                             adaptive  = FALSE,
                             tol       = 1e-6,
                             max_iter  = 100) {

  P <- length(B_ols)
  B <- B_ols

  # Helper: compute kappa from a moment triplet; returns NA if near-Gaussian
  .kappa_from_moments <- function(m2_, m4_, m6_) {
    denom <- m4_ - 3 * m2_^2
    if (abs(denom) < .Machine$double.eps * 1e6) return(NA_real_)
    (m6_ - 3 * m4_ * m2_) / denom
  }

  # Initial kappa (from OLS moments — used always for fixed mode, and as
  # starting value for adaptive mode)
  kappa <- .kappa_from_moments(m2, m4, m6)
  if (is.na(kappa)) {
    return(list(B = B_ols, converged = TRUE, iter = 0L,
                kappa = NA_real_, near_gaussian = TRUE))
  }

  # In fixed mode, pre-compute Y-dependent polynomial coefficients once
  if (!adaptive) {
    B_c <- -3 * Y
    C_c <- 3 * Y^2 - kappa
    D_c <- kappa * Y - Y^3
  }

  converged <- FALSE
  iter      <- 0L

  for (i in seq_len(max_iter)) {
    iter <- as.integer(i)
    Yx   <- as.vector(X %*% B)
    eps  <- Y - Yx

    # Adaptive mode: re-estimate kappa from current residuals
    if (adaptive) {
      eps_c <- eps - mean(eps)            # centre before moment estimation
      m2_c  <- mean(eps_c^2)
      m4_c  <- mean(eps_c^4)
      m6_c  <- mean(eps_c^6)
      kappa_new <- .kappa_from_moments(m2_c, m4_c, m6_c)
      if (is.na(kappa_new)) { converged <- TRUE; break }  # converged to Gaussian
      kappa <- kappa_new
      B_c <- -3 * Y
      C_c <- 3 * Y^2 - kappa
      D_c <- kappa * Y - Y^3
    }

    # Score Z_p = sum_v eps_v*(kappa - eps_v^2) * x_vp
    Z1 <- Yx^3 + B_c * Yx^2 + C_c * Yx + D_c
    Z  <- colSums(Z1 * X)

    # Jacobian = X' diag(3*eps^2 - kappa) X
    w   <- 3 * eps^2 - kappa
    JZs <- crossprod(X * w, X)

    Q <- tryCatch(solve(JZs, Z), error = function(e) rep(NA_real_, P))
    if (anyNA(Q)) break

    B <- B - Q

    if (sqrt(sum(Q^2)) < tol) {
      converged <- TRUE
      break
    }
  }

  # Divergence guard
  if (sqrt(sum((B - B_ols)^2)) > 10 || anyNA(B)) {
    B         <- B_ols
    converged <- FALSE
  }

  list(B = B, converged = converged, iter = iter, kappa = kappa,
       near_gaussian = FALSE)
}

# ---------------------------------------------------------------------------
#' Fit linear regression by PMM3 (Method of Polynomial Moments, S=3)
#'
#' Solves the PMM3 estimating equations for symmetric non-Gaussian errors
#' using Newton-Raphson iteration with OLS starting values.
#'
#' The theoretical relative efficiency vs OLS is:
#'   g3 = 1 - gamma4^2 / (6 + 9*gamma4 + gamma6)
#' For TN(lambda): g3 decreases from ~1 (lambda=0) toward 0 (lambda→∞).
#'
#' @param formula    R formula, e.g. y ~ x1 + x2
#' @param data       data.frame (or NULL for formula evaluation in parent)
#' @param max_iter   Maximum Newton-Raphson iterations (default 100)
#' @param tol        Convergence tolerance ||Q||_2 (default 1e-6)
#' @param use_ols_start  Always TRUE; OLS provides starting values
#'
#' @return S4 object of class PMM3fit with slots:
#'   @@coefficients — beta estimates (named numeric)
#'   @@fitted       — fitted values
#'   @@residuals    — final residuals
#'   @@convergence  — logical
#'   @@iterations   — integer
#'   @@m2, @@m4, @@m6   — residual moments (from OLS step)
#'   @@gamma4, @@gamma6 — standardised cumulants
#'   @@g_coefficient    — theoretical PMM3/OLS efficiency ratio
# ---------------------------------------------------------------------------
lm_pmm3 <- function(formula, data = NULL,
                    max_iter = 100, tol = 1e-6,
                    use_ols_start = TRUE,
                    adaptive = FALSE) {

  cl  <- match.call()
  mf  <- model.frame(formula, data = data)
  y   <- model.response(mf)
  X   <- model.matrix(attr(mf, "terms"), data = mf)
  n   <- length(y)
  p   <- ncol(X)

  # Step 1: OLS starting values and moment estimation
  ols  <- lm.fit(X, y)
  mom  <- .compute_moments_pmm3(ols$residuals)
  .check_symmetry_pmm3(ols$residuals)

  # Step 2: Newton-Raphson
  nr <- .pmm3_nr_solver(
    B_ols    = unname(ols$coefficients),
    X        = X,
    Y        = y,
    m2       = mom$mu2,
    m4       = mom$mu4,
    m6       = mom$mu6,
    adaptive = adaptive,
    tol      = tol,
    max_iter = max_iter
  )

  # Step 3: Assemble result
  beta_hat   <- setNames(nr$B, colnames(X))
  fitted_val <- as.vector(X %*% beta_hat)
  resid_val  <- y - fitted_val

  gamma4 <- mom$mu4 / mom$mu2^2 - 3
  gamma6 <- mom$mu6 / mom$mu2^3 - 15 * (mom$mu4 / mom$mu2^2) + 30
  g_coef <- .pmm3_g_coefficient(mom$mu2, mom$mu4, mom$mu6)

  new("PMM3fit",
      coefficients  = unname(beta_hat),
      fitted        = fitted_val,
      residuals     = resid_val,
      convergence   = nr$converged,
      iterations    = nr$iter,
      call          = cl,
      terms         = attr(mf, "terms"),
      method        = if (isTRUE(nr$near_gaussian)) "PMM3→OLS (near-Gaussian)"
                      else if (adaptive) "PMM3-adaptive" else "PMM3",
      m2            = mom$mu2,
      m4            = mom$mu4,
      m6            = mom$mu6,
      gamma4        = gamma4,
      gamma6        = gamma6,
      g_coefficient = g_coef)
}

# ---------------------------------------------------------------------------
# S4 methods
# ---------------------------------------------------------------------------

setMethod("show", "PMM3fit", function(object) {
  cat("PMM3 Regression  [", object@method, "]\n", sep = "")
  coef_named <- setNames(object@coefficients,
                         attr(object@terms, "term.labels") |>
                           (\(tl) c("(Intercept)", tl))() |>
                           head(length(object@coefficients)))
  cat("Coefficients:\n"); print(coef_named)
  cat(sprintf(
    "gamma4=%.4f  gamma6=%.4f  g3=%.4f  iter=%d  converged=%s\n",
    object@gamma4, object@gamma6, object@g_coefficient,
    object@iterations, object@convergence))
  invisible(object)
})

setGeneric("coef",    function(object, ...) standardGeneric("coef"))
setGeneric("predict", function(object, ...) standardGeneric("predict"))
setGeneric("residuals", function(object, ...) standardGeneric("residuals"))

setMethod("coef", "PMM3fit", function(object, ...) {
  setNames(object@coefficients,
           tryCatch(
             c("(Intercept)", attr(object@terms, "term.labels")),
             error = function(e) seq_along(object@coefficients)
           ) |> head(length(object@coefficients)))
})

setMethod("predict", "PMM3fit", function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object@fitted)
  X_new <- model.matrix(delete.response(object@terms), newdata)
  as.vector(X_new %*% object@coefficients)
})

setMethod("residuals", "PMM3fit", function(object, ...) object@residuals)

message("[pmm3_skeleton.R] Real PMM3 Newton-Raphson loaded. STUB=FALSE.")
