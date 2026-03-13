# =============================================================================
# R/03_monte_carlo.R
# Monte Carlo study: PMM3 vs MLE-TN vs OLS vs PMM2
#
# Plan reference: Block 2 → Table T2
#   - g3_empirical vs g3_theoretical
#   - ARE = MSE(OLS) / MSE(method) for all methods
#   - Convergence rates for iterative methods
#
# Key hypothesis (H1): g3_empirical → g3_theoretical as n → ∞
# Control check (H4): ARE(PMM2) ≈ 1.0 for all lambda (gamma3=0 → PMM2 = OLS)
# =============================================================================

# ---------------------------------------------------------------------------
# Single Monte Carlo replication
# Returns a 1-row data.frame with beta1 estimates from all four methods
# ---------------------------------------------------------------------------
.one_mc_rep <- function(lambda, n, beta = CFG$beta_true, eta = CFG$eta_true) {
  x   <- rnorm(n)
  eps <- rtn(n, xi = 0, eta = eta, lambda = lambda)
  y   <- beta["beta0"] + beta["beta1"] * x + eps
  dat <- data.frame(y = y, x = x)

  # Helper: safely extract beta1 coefficient
  safe_b1 <- function(expr) tryCatch(eval(expr), error = function(e) NA_real_)

  # OLS
  b1_ols  <- safe_b1(quote(coef(lm(y ~ x, dat))["x"]))

  # MLE-TN
  b1_mle  <- safe_b1(quote(mle_tn(y ~ x, dat)$coefficients["x"]))
  ok_mle  <- !is.na(b1_mle)

  # PMM3 — EstemPMM::lm_pmm3 or local stub
  pmm3_res <- tryCatch({
    f <- lm_pmm3(y ~ x, data = dat)
    # If stub is active, lm_pmm3 returns OLS values with a warning
    # STUB flag set in pmm3_skeleton.R
    is_stub  <- exists("STUB", envir = globalenv()) && isTRUE(get("STUB"))
    list(b1 = f@coefficients[2], ok = f@convergence, stub = is_stub)
  }, error = function(e) list(b1 = NA_real_, ok = FALSE, stub = FALSE))
  b1_pmm3 <- if (pmm3_res$stub) NA_real_ else pmm3_res$b1
  ok_pmm3 <- if (pmm3_res$stub) FALSE    else pmm3_res$ok

  # PMM2 — control: should be ≈ OLS when gamma3 = 0 (symmetric TN)
  b1_pmm2 <- tryCatch(
    EstemPMM::lm_pmm2(y ~ x, data = dat)@coefficients[2],
    error = function(e) NA_real_
  )

  data.frame(
    lambda  = lambda, n = n,
    ols_b1  = b1_ols,
    mle_b1  = b1_mle,  mle_ok  = ok_mle,
    pmm3_b1 = b1_pmm3, pmm3_ok = ok_pmm3,
    pmm2_b1 = b1_pmm2,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Summary metrics from M replications for one (lambda, n) block
# ---------------------------------------------------------------------------
.compute_metrics <- function(df, beta_true) {
  var_ols  <- var(df$ols_b1, na.rm = TRUE)
  mse_ols  <- mean((df$ols_b1 - beta_true)^2, na.rm = TRUE)

  methods <- list(
    OLS  = list(col = "ols_b1",  ok_col = NULL),
    MLE  = list(col = "mle_b1",  ok_col = "mle_ok"),
    PMM3 = list(col = "pmm3_b1", ok_col = "pmm3_ok"),
    PMM2 = list(col = "pmm2_b1", ok_col = NULL)
  )

  rows <- lapply(names(methods), function(m) {
    col  <- methods[[m]]$col
    vals <- df[[col]]
    vals <- vals[!is.na(vals)]
    if (length(vals) < 10) return(NULL)

    var_m   <- var(vals)
    mse_m   <- mean((vals - beta_true)^2)

    data.frame(
      lambda       = df$lambda[1],
      n            = df$n[1],
      method       = m,
      bias         = mean(vals) - beta_true,
      variance     = var_m,
      mse          = mse_m,
      are          = mse_ols / mse_m,        # > 1 means method beats OLS
      g3_empirical = if (m == "PMM3") var_m / var_ols else NA_real_,
      conv_rate    = if (!is.null(methods[[m]]$ok_col))
                       mean(df[[methods[[m]]$ok_col]], na.rm = TRUE)
                     else 1.0,
      n_valid      = length(vals),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

# ---------------------------------------------------------------------------
# Run MC for one (lambda, n) combination
# ---------------------------------------------------------------------------
run_mc_block <- function(lambda, n, M = CFG$mc_replications, seed = CFG$seed) {
  set.seed(seed + as.integer(lambda * 100) + n)
  message(sprintf("  MC | lambda=%4.1f  n=%4d  M=%d", lambda, n, M))

  reps <- lapply(seq_len(M), function(i) .one_mc_rep(lambda, n))
  df   <- do.call(rbind, reps)
  .compute_metrics(df, beta_true = CFG$beta_true["beta1"])
}

# ---------------------------------------------------------------------------
# Full MC study over lambda x n grid with disk caching
# ---------------------------------------------------------------------------
#' Run complete Monte Carlo study
#'
#' Results cached to results/mc_cache/mc_full_results.rds.
#' Set force_rerun=TRUE to ignore cache and recompute.
#'
#' @param force_rerun  Logical; if TRUE, discard cache and rerun
run_full_mc_study <- function(lambda_grid  = CFG$lambda_grid,
                               sample_sizes = CFG$sample_sizes,
                               M            = CFG$mc_replications,
                               force_rerun  = FALSE) {
  cache_file <- file.path(CFG$paths$mc_cache, "mc_full_results.rds")

  if (!force_rerun && file.exists(cache_file)) {
    message("Loading MC cache: ", cache_file)
    return(readRDS(cache_file))
  }

  n_blocks <- length(lambda_grid) * length(sample_sizes)
  message(sprintf(
    "Starting MC study: %d lambda x %d n = %d blocks, M=%d each",
    length(lambda_grid), length(sample_sizes), n_blocks, M
  ))
  t0 <- proc.time()

  all <- list()
  for (lam in lambda_grid)
    for (n in sample_sizes)
      all[[paste0("l", lam, "_n", n)]] <- run_mc_block(lam, n, M)

  full_df <- do.call(rbind, Filter(Negate(is.null), all))
  rownames(full_df) <- NULL

  dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
  saveRDS(full_df, cache_file)

  elapsed <- round((proc.time() - t0)["elapsed"] / 60, 1)
  message(sprintf("MC study done in %.1f min. Results: %s", elapsed, cache_file))
  full_df
}
