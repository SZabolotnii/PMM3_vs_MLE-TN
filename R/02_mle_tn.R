# =============================================================================
# R/02_mle_tn.R
# Maximum Likelihood Estimation for TN linear regression
# Competitor method for PMM3 comparison (Salinas et al. 2023, Table 5)
#
# Model: y = X*beta + eps,  eps ~ TN(0, eta, lambda)
# Log-likelihood: formula (9) in Salinas et al.
#
# Exported functions:
#   mle_tn()   — fit TN regression by L-BFGS-B
# =============================================================================

# Internal: negative log-likelihood
# par = c(beta[1..p], log_eta, lambda)
.ll_tn <- function(par, y, X, negative = TRUE) {
  p    <- ncol(X)
  beta <- par[seq_len(p)]
  eta  <- exp(par[p + 1])      # log-parameterised → always positive
  lam  <- par[p + 2]
  if (eta <= 0 || lam < 0) return(if (negative) Inf else -Inf)

  z  <- (y - as.vector(X %*% beta)) / eta
  n  <- length(y)
  lz <- lam * z
  log_cosh_sum <- sum(abs(lz) + log1p(exp(-2 * abs(lz))) - log(2))
  ll <- -n * log(2 * pi) / 2 - n * log(eta) -
        n * lam^2 / 2 - sum(z^2) / 2 + log_cosh_sum
  if (negative) -ll else ll
}

#' Fit TN linear regression by MLE
#'
#' Maximises the TN log-likelihood (formula 9, Salinas et al. 2023)
#' using L-BFGS-B with OLS starting values.
#'
#' @param formula  R formula, e.g. y ~ x1 + x2
#' @param data     data.frame
#' @param lambda_init  Starting value for lambda (default: 1.0)
#' @param maxit    Maximum L-BFGS-B iterations
#'
#' @return List with:
#'   coefficients — beta estimates (named)
#'   eta          — scale estimate
#'   lambda       — shape estimate
#'   loglik       — maximised log-likelihood
#'   se_beta      — SE from Hessian (NA if Hessian singular)
#'   aic, bic     — information criteria (k = p + 2)
#'   converged    — logical
#'   fitted, residuals
mle_tn <- function(formula, data = NULL, lambda_init = 1.0, maxit = 500) {
  cl <- match.call()
  mf <- model.frame(formula, data = data)
  y  <- model.response(mf)
  X  <- model.matrix(attr(mf, "terms"), data = mf)
  n  <- length(y)
  p  <- ncol(X)

  # OLS starting values
  ols  <- lm.fit(X, y)
  par0 <- c(coef(ols), log(sd(ols$residuals)), lambda_init)

  opt <- optim(
    par     = par0,
    fn      = .ll_tn,
    y = y, X = X, negative = TRUE,
    method  = "L-BFGS-B",
    lower   = c(rep(-Inf, p), -6,  0.001),
    upper   = c(rep( Inf, p),  6, 15.0),
    control = list(maxit = maxit, factr = 1e7),
    hessian = TRUE
  )

  # Standard errors via observed Fisher information
  se <- tryCatch(
    sqrt(pmax(diag(solve(opt$hessian)), 0)),
    error = function(e) rep(NA_real_, length(par0))
  )

  beta_hat   <- setNames(opt$par[seq_len(p)], colnames(X))
  fitted_val <- as.vector(X %*% beta_hat)
  terms_obj  <- attr(mf, "terms")

  structure(
    list(
      call         = cl,
      terms        = terms_obj,
      coefficients = beta_hat,
      eta          = unname(exp(opt$par[p + 1])),
      lambda       = unname(opt$par[p + 2]),
      loglik       = -opt$value,
      se_beta      = se[seq_len(p)],
      se_eta       = se[p + 1] * exp(opt$par[p + 1]),
      se_lambda    = se[p + 2],
      converged    = (opt$convergence == 0),
      aic          = 2 * (p + 2) + 2 * opt$value,
      bic          = log(n) * (p + 2) + 2 * opt$value,
      n_params     = p + 2,
      n            = n,
      fitted       = fitted_val,
      residuals    = y - fitted_val
    ),
    class = "mle_tn_fit"
  )
}

# S3 methods for pipeline compatibility with bootstrap / LOO code
coef.mle_tn_fit    <- function(object, ...) object$coefficients
predict.mle_tn_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted)
  X_new <- model.matrix(delete.response(object$terms), newdata)
  as.vector(X_new %*% object$coefficients)
}
print.mle_tn_fit <- function(x, ...) {
  cat("MLE-TN Regression\n")
  cat("Coefficients:\n"); print(x$coefficients)
  cat(sprintf("eta=%.4f  lambda=%.4f  loglik=%.4f  converged=%s\n",
              x$eta, x$lambda, x$loglik, x$converged))
  invisible(x)
}
