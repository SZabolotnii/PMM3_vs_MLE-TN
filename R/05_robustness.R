# =============================================================================
# R/05_robustness.R
# Robustness analysis — Plan Blocks 5 & 6
#
# Block 5: Misspecification experiment
#   DGP ≠ TN (Uniform, Triangular, Logistic, Student-t)
#   Both MLE-TN and PMM3 applied → who degrades more?
#   Hypothesis: PMM3 >= MLE-TN under misspecification (MLE forces wrong shape)
#
# Block 6: Bimodal TN test (lambda > 1)
#   Monitor PMM3 convergence rate and bias as lambda crosses 1
#   Critical: PMM3 designed for unimodal symmetric distributions
# =============================================================================

# ---------------------------------------------------------------------------
# Block 5: Symmetric unit-variance distributions for misspecification test
# ---------------------------------------------------------------------------

# All generators produce E(X)=0, Var(X)=1
.rgen <- list(
  uniform    = function(n) runif(n, -sqrt(3), sqrt(3)),
  triangular = function(n) (runif(n) + runif(n) - 1) * sqrt(6),
  logistic   = function(n) rlogis(n, scale = sqrt(3) / pi),
  student10  = function(n) rt(n, df = 10) / sqrt(10 / 8)
)

#' Misspecification experiment: are PMM3 and MLE-TN robust when DGP ≠ TN?
#'
#' @param dist_name  One of: "uniform", "triangular", "logistic", "student10"
#' @param n          Sample size
#' @param M          Monte Carlo replications
#' @param seed       Random seed
#' @return 1-row data.frame with ARE and bias for MLE-TN and PMM3
run_misspec_block <- function(dist_name, n = 100, M = 1000, seed = CFG$seed) {
  gen <- .rgen[[dist_name]]
  if (is.null(gen)) stop("Unknown distribution: ", dist_name,
                         ". Choose from: ", paste(names(.rgen), collapse=", "))

  set.seed(seed + which(names(.rgen) == dist_name))
  x_fix <- rnorm(n)           # fixed design — same X across replications
  beta  <- CFG$beta_true

  reps <- lapply(seq_len(M), function(i) {
    eps <- gen(n)
    y   <- beta["beta0"] + beta["beta1"] * x_fix + eps
    dat <- data.frame(y = y, x = x_fix)

    b_ols  <- tryCatch(coef(lm(y ~ x, dat))["x"],           error = function(e) NA_real_)
    b_mle  <- tryCatch(mle_tn(y ~ x, dat)$coefficients["x"], error = function(e) NA_real_)
    b_pmm3 <- tryCatch(lm_pmm3(y ~ x, dat)@coefficients[2], error = function(e) NA_real_)

    c(ols = b_ols, mle = b_mle, pmm3 = b_pmm3)
  })

  df <- as.data.frame(do.call(rbind, reps))

  var_ols <- var(df$ols, na.rm = TRUE)
  data.frame(
    distribution = dist_name,
    n            = n,
    are_mle      = var_ols / var(df$mle,  na.rm = TRUE),
    are_pmm3     = var_ols / var(df$pmm3, na.rm = TRUE),
    bias_mle     = mean(df$mle,  na.rm = TRUE) - beta["beta1"],
    bias_pmm3    = mean(df$pmm3, na.rm = TRUE) - beta["beta1"],
    na_pct_mle   = mean(is.na(df$mle))  * 100,
    na_pct_pmm3  = mean(is.na(df$pmm3)) * 100,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Block 6: Bimodal TN convergence test
# ---------------------------------------------------------------------------

#' Test PMM3 behaviour as lambda crosses the unimodal/bimodal boundary (lambda=1)
#'
#' Monitors:
#'   conv_rate   — fraction of PMM3 runs that converge
#'   bias_pmm3   — bias in beta1 estimate
#'   mean_iter   — average Newton-Raphson iterations
#'   na_pct      — fraction of NA results (convergence failure)
#'
#' @param lambda_grid  Grid of lambda values (default includes 0.8, 1.0, 1.2, ...)
#' @param n            Sample size
#' @param M            Monte Carlo replications
run_bimodal_test <- function(lambda_grid = CFG$bimodal_lambda_grid,
                              n = 100, M = 500, seed = CFG$seed) {
  beta_true <- CFG$beta_true["beta1"]

  do.call(rbind, lapply(lambda_grid, function(lam) {
    set.seed(seed + as.integer(lam * 10))
    message(sprintf("  Bimodal test | lambda=%.1f", lam))

    reps <- replicate(M, {
      x   <- rnorm(n)
      y   <- CFG$beta_true["beta0"] + beta_true * x + rtn(n, lambda = lam)
      dat <- data.frame(y = y, x = x)
      tryCatch({
        f <- lm_pmm3(y ~ x, dat)
        c(b1   = unname(f@coefficients[2]),
          conv  = as.integer(f@convergence),
          iter  = as.integer(f@iterations))
      }, error = function(e) c(b1 = NA_real_, conv = 0L, iter = NA_integer_))
    })

    data.frame(
      lambda     = lam,
      is_bimodal = lam > 1,
      bias_pmm3  = mean(reps["b1", ], na.rm = TRUE) - beta_true,
      conv_rate  = mean(reps["conv",], na.rm = TRUE),
      mean_iter  = mean(reps["iter",], na.rm = TRUE),
      na_pct     = mean(is.na(reps["b1",])) * 100,
      stringsAsFactors = FALSE
    )
  }))
}
