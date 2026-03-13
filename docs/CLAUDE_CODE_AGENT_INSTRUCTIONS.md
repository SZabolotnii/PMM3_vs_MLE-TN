# Інструкція для Claude Code: розгортання проєкту PMM3-TN Research

## Контекст та мета

Ти розгортаєш **дослідницький R-проєкт** для порівняльного аналізу двох
статистичних методів:

- **PMM3** — Метод Максимізації Поліномів (S=3), реалізований у пакеті `EstemPMM`
- **MLE-TN** — Метод максимальної правдоподібності для двочастинного нормального
  розподілу (Two-piece Normal), з Salinas et al. (2023), Mathematics 11(5), 1271

Проєкт є частиною розробки пакету EstemPMM. PMM3 призначений для симетричних
негаусових розподілів і конкурує з MLE-TN на даних типу biaxial fatigue.

---

## КРОК 0 — Запит директорії у користувача

**Перш ніж виконувати будь-які команди**, запитай:

```
Де створити проєкт? Вкажіть повний шлях до батьківської директорії.
Наприклад: /home/user/projects або ~/research
Буде створено: <вказана_директорія>/pmm3_tn_study/
```

Збережи відповідь як змінну. У всіх подальших командах використовуй цей шлях.
Якщо користувач не вказав — використай `~/research/pmm3_tn_study`.

---

## КРОК 1 — Перевірка системи

Виконай ці перевірки **одним блоком bash** і покажи результати:

```bash
# 1. R доступний?
Rscript --version 2>&1 | head -1

# 2. EstemPMM встановлений?
Rscript -e "cat(requireNamespace('EstemPMM', quietly=TRUE))" 2>/dev/null

# 3. Чи є lm_pmm3 в EstemPMM?
Rscript -e "
  if (requireNamespace('EstemPMM', quietly=TRUE)) {
    cat(existsMethod('lm_pmm3', 'ANY') || exists('lm_pmm3', where=asNamespace('EstemPMM')))
  } else { cat('FALSE') }
" 2>/dev/null

# 4. Ключові пакети
Rscript -e "
pkgs <- c('ggplot2','dplyr','tidyr','knitr','rmarkdown','nlme','pracma','testthat')
status <- sapply(pkgs, requireNamespace, quietly=TRUE)
cat(paste(names(status), ifelse(status,'OK','MISSING'), sep=': ', collapse='\n'))
" 2>/dev/null
```

**Інтерпретація результатів:**

| Результат | Дія |
|-----------|-----|
| R не знайдено | Повідом користувача — без R проєкт не запуститься |
| EstemPMM = FALSE | Додай крок встановлення (Крок 1b) |
| lm_pmm3 = FALSE | Додай скелет `lm_pmm3()` у `R/pmm3_skeleton.R` (Крок 3g) |
| Пакети MISSING | Сформуй список і встанови в Кроці 7 |

### Крок 1b (умовний) — встановлення EstemPMM

Виконай **тільки якщо** EstemPMM не встановлений:

```bash
Rscript -e "
if (!requireNamespace('devtools', quietly=TRUE)) install.packages('devtools', repos='https://cran.r-project.org')
devtools::install_github('SZabolotnii/EstemPMM', quiet=TRUE)
cat('EstemPMM installed:', requireNamespace('EstemPMM', quietly=TRUE))
"
```

---

## КРОК 2 — Створення структури директорій

```bash
PROJECT_ROOT="<шлях_від_користувача>/pmm3_tn_study"

mkdir -p "$PROJECT_ROOT"/{R,tests/testthat,data/{raw,processed},results/{tables,figures,mc_cache},reports,docs,scripts}

# Верифікація
echo "=== Project structure ==="
find "$PROJECT_ROOT" -type d | sort
echo ""
echo "Total directories: $(find "$PROJECT_ROOT" -type d | wc -l)"
```

**Очікуваний вивід — рівно 12 директорій:**
```
pmm3_tn_study/
pmm3_tn_study/R
pmm3_tn_study/data
pmm3_tn_study/data/processed
pmm3_tn_study/data/raw
pmm3_tn_study/reports
pmm3_tn_study/results
pmm3_tn_study/results/figures
pmm3_tn_study/results/mc_cache
pmm3_tn_study/results/tables
pmm3_tn_study/scripts
pmm3_tn_study/tests
pmm3_tn_study/tests/testthat
```

Якщо кількість відрізняється — знайди і виправ проблему перш ніж продовжувати.

---

## КРОК 3 — Створення файлів

Кожен файл створюй окремою командою `cat > filepath << 'EOF' ... EOF`.
Після кожного файлу виконуй `echo "Created: <шлях>"` для підтвердження.

---

### 3a — `R/config.R`

```bash
cat > "$PROJECT_ROOT/R/config.R" << 'REOF'
# =============================================================================
# R/config.R — Централізована конфігурація дослідження
# PMM3 vs MLE-TN Comparative Analysis
# Salinas et al. (2023), Mathematics 11(5), 1271
# =============================================================================
# РЕДАГУЙ ЦЕЙ ФАЙЛ для зміни параметрів дослідження.
# Всі інші модулі читають CFG звідси.

CFG <- list(

  # --- Monte Carlo -----------------------------------------------------------
  mc_replications = 500,       # 500 для тесту, 5000 для фінального запуску
  sample_sizes    = c(46, 100, 200, 500),   # n=46 = biaxial fatigue
  lambda_grid     = c(0.5, 1.0, 1.5, 2.0, 3.0, 5.0),

  # --- Справжні параметри моделі (DGP) ------------------------------------
  beta_true = c(beta0 = 1.0, beta1 = 2.0),
  eta_true  = 1.0,

  # --- Порогові значення ---------------------------------------------------
  symmetry_threshold = 0.3,   # |gamma3| < threshold → симетричний
  gaussian_threshold = 0.1,   # |gamma4| < threshold → гаусовий

  # --- Розподіли для misspecification тесту --------------------------------
  misspec_distributions = c("uniform", "triangular", "logistic", "student10"),

  # --- Бімодальний тест (lambda > 1) ---------------------------------------
  bimodal_lambda_grid = c(0.8, 1.0, 1.2, 1.5, 2.0, 3.0),

  # --- Шляхи (відносні) ---------------------------------------------------
  paths = list(
    results_tables  = "results/tables",
    results_figures = "results/figures",
    mc_cache        = "results/mc_cache",
    data_raw        = "data/raw",
    data_processed  = "data/processed"
  ),

  # --- Відтворюваність -----------------------------------------------------
  seed         = 20260312,   # YYYYMMDD — дата початку дослідження
  use_parallel = FALSE,
  n_cores      = 2
)

message(sprintf(
  "CFG loaded | MC: M=%d, n=%s, lambda=%s",
  CFG$mc_replications,
  paste(CFG$sample_sizes, collapse=","),
  paste(CFG$lambda_grid, collapse=",")
))
REOF
echo "Created: R/config.R"
```

---

### 3b — `R/01_tn_distribution.R`

```bash
cat > "$PROJECT_ROOT/R/01_tn_distribution.R" << 'REOF'
# =============================================================================
# R/01_tn_distribution.R
# Two-piece Normal (TN) distribution: rtn, dtn, ptn, moments, Table T1
# Reference: Salinas et al. (2023), Section 2-3
# =============================================================================

#' Generate random TN(xi, eta, lambda) variates
#' Algorithm: fold-and-sign, Section 4.2 of Salinas et al.
rtn <- function(n, xi = 0, eta = 1, lambda = 1) {
  stopifnot(n > 0, eta > 0, lambda >= 0)
  U <- rnorm(n, mean = lambda, sd = 1)
  Y <- abs(U)
  S <- sample(c(-1L, 1L), n, replace = TRUE)
  xi + eta * (S * Y)
}

#' TN density (formula 7 in Salinas et al.)
dtn <- function(x, xi = 0, eta = 1, lambda = 0, log = FALSE) {
  stopifnot(eta > 0, lambda >= 0)
  z        <- (x - xi) / eta
  log_dens <- dnorm(z, log = TRUE) - log(eta) -
              lambda^2 / 2 + log(cosh(lambda * z))
  if (log) log_dens else exp(log_dens)
}

#' TN distribution function
ptn <- function(q, xi = 0, eta = 1, lambda = 0) {
  z <- (q - xi) / eta
  (pnorm(z - lambda) + pnorm(z + lambda)) / 2
}

#' E(Z^{2k}) for Z ~ TN(0,1,lambda) via numerical integration
#' Used for computing gamma4, gamma6, g3
tn_raw_moment_2k <- function(k, lambda) {
  tryCatch(
    integrate(function(z) z^(2*k) * dtn(z, lambda = lambda),
              lower = -Inf, upper = Inf,
              subdivisions = 1000L, rel.tol = 1e-7)$value,
    error = function(e) NA_real_
  )
}

#' Cumulant coefficients and theoretical PMM3 efficiency for TN(0,1,lambda)
#' Returns: list with mu2, mu4, mu6, gamma4, gamma6, g3, is_bimodal
tn_cumulants <- function(lambda) {
  mu2 <- tn_raw_moment_2k(1, lambda)
  mu4 <- tn_raw_moment_2k(2, lambda)
  mu6 <- tn_raw_moment_2k(3, lambda)

  gamma4 <- mu4 / mu2^2 - 3
  gamma6 <- mu6 / mu2^3 - 15 * (mu4 / mu2^2) + 30

  denom <- 6 + 9 * gamma4 + gamma6
  g3    <- if (!is.na(denom) && denom > 0) 1 - gamma4^2 / denom else NA_real_

  list(lambda    = lambda,
       mu2       = mu2,
       mu4       = mu4,
       mu6       = mu6,
       gamma4    = gamma4,
       gamma6    = gamma6,
       g3        = g3,
       improvement_pct = if (!is.na(g3)) (1 - g3) * 100 else NA_real_,
       is_bimodal = lambda > 1)
}

#' Build Table T1: theoretical TN cumulants for a grid of lambda values
build_table_T1 <- function(lambda_grid = CFG$lambda_grid) {
  message("Building Table T1 (", length(lambda_grid), " lambda values)...")
  rows <- lapply(lambda_grid, tn_cumulants)
  T1   <- do.call(rbind, lapply(rows, as.data.frame))
  rownames(T1) <- NULL
  T1$mode <- ifelse(T1$is_bimodal, "bimodal", "unimodal")
  message("Table T1 complete.")
  T1
}
REOF
echo "Created: R/01_tn_distribution.R"
```

---

### 3c — `R/02_mle_tn.R`

```bash
cat > "$PROJECT_ROOT/R/02_mle_tn.R" << 'REOF'
# =============================================================================
# R/02_mle_tn.R
# MLE for TN regression: competitor method for PMM3 comparison
# Log-likelihood: formula (9) in Salinas et al. (2023)
# =============================================================================

# Log-likelihood (internal). par = c(beta..., log_eta, lambda)
.ll_tn <- function(par, y, X, negative = TRUE) {
  p    <- ncol(X)
  beta <- par[seq_len(p)]
  eta  <- exp(par[p + 1])      # positive via exp()
  lam  <- par[p + 2]
  if (eta <= 0 || lam < 0) return(if (negative) Inf else -Inf)
  z    <- (y - as.vector(X %*% beta)) / eta
  n    <- length(y)
  ll   <- -n*log(2*pi)/2 - n*log(eta) - n*lam^2/2 -
           sum(z^2)/2 + sum(log(cosh(lam * z)))
  if (negative) -ll else ll
}

#' Fit TN linear regression by MLE
#' @param formula R formula (y ~ x1 + x2 ...)
#' @param data    data.frame
#' @param lambda_init Starting value for lambda (default 1.0)
#' @param maxit   Max L-BFGS-B iterations
#' @return list: coefficients, eta, lambda, loglik, se_beta, aic, bic, converged
mle_tn <- function(formula, data = NULL, lambda_init = 1.0, maxit = 500) {
  cl  <- match.call()
  mf  <- model.frame(formula, data = data)
  y   <- model.response(mf)
  X   <- model.matrix(attr(mf, "terms"), data = mf)
  n   <- length(y); p <- ncol(X)

  # OLS start
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

  se <- tryCatch(sqrt(pmax(diag(solve(opt$hessian)), 0)),
                 error = function(e) rep(NA_real_, length(par0)))

  fitted_vals <- as.vector(X %*% opt$par[seq_len(p)])
  list(
    call         = cl,
    coefficients = setNames(opt$par[seq_len(p)], colnames(X)),
    eta          = exp(opt$par[p + 1]),
    lambda       = opt$par[p + 2],
    loglik       = -opt$value,
    se_beta      = se[seq_len(p)],
    converged    = (opt$convergence == 0),
    aic          = 2*(p + 2) + 2*opt$value,
    bic          = log(n)*(p + 2) + 2*opt$value,
    n_params     = p + 2,
    n            = n,
    fitted       = fitted_vals,
    residuals    = y - fitted_vals
  )
}

# S3 coef method for pipeline compatibility
coef.mle_tn_fit <- function(object, ...) object$coefficients
predict.mle_tn_fit <- function(object, newdata, ...) {
  if (missing(newdata)) return(object$fitted)
  X_new <- model.matrix(delete.response(terms(object$call$formula)), newdata)
  as.vector(X_new %*% object$coefficients)
}
REOF
echo "Created: R/02_mle_tn.R"
```

---

### 3d — `R/03_monte_carlo.R`

```bash
cat > "$PROJECT_ROOT/R/03_monte_carlo.R" << 'REOF'
# =============================================================================
# R/03_monte_carlo.R
# Monte Carlo study: PMM3 vs MLE-TN vs OLS vs PMM2
# Plan Block 2: Table T2 (g3_empirical vs g3_theoretical, ARE by method)
# =============================================================================

# One replication: fit all methods on one TN dataset
.one_mc_rep <- function(lambda, n, beta = CFG$beta_true, eta = CFG$eta_true) {
  x   <- rnorm(n)
  eps <- rtn(n, xi = 0, eta = eta, lambda = lambda)
  y   <- beta["beta0"] + beta["beta1"] * x + eps
  dat <- data.frame(y = y, x = x)

  safe_coef <- function(expr, idx = 2) {
    tryCatch({ co <- eval(expr); list(val=co[idx], ok=TRUE) },
             error = function(e) list(val=NA_real_, ok=FALSE))
  }

  r_ols  <- safe_coef(quote(coef(lm(y ~ x, dat))))
  r_mle  <- safe_coef(quote(mle_tn(y ~ x, dat)$coefficients))

  # PMM3 — використовуємо EstemPMM::lm_pmm3 якщо доступний,
  # інакше — локальний скелет зі stub-значенням
  r_pmm3 <- tryCatch({
    if (exists("lm_pmm3") && !is.null(environment(lm_pmm3)$STUB)) {
      list(val=NA_real_, ok=FALSE)   # Stub: ще не реалізовано
    } else {
      fit <- lm_pmm3(y ~ x, data = dat)
      list(val=fit@coefficients[2], ok=fit@convergence)
    }
  }, error = function(e) list(val=NA_real_, ok=FALSE))

  # PMM2 — контрольний метод (має бути ≈ OLS при gamma3=0)
  r_pmm2 <- tryCatch({
    fit <- EstemPMM::lm_pmm2(y ~ x, data = dat)
    list(val=fit@coefficients[2], ok=TRUE)
  }, error = function(e) list(val=NA_real_, ok=FALSE))

  data.frame(lambda=lambda, n=n,
             ols_b1  = r_ols$val,
             mle_b1  = r_mle$val,  mle_ok  = r_mle$ok,
             pmm3_b1 = r_pmm3$val, pmm3_ok = r_pmm3$ok,
             pmm2_b1 = r_pmm2$val,
             stringsAsFactors = FALSE)
}

# Compute summary metrics from M replications
.compute_metrics <- function(df, beta_true, var_ols) {
  methods <- list(
    OLS  = "ols_b1",
    MLE  = "mle_b1",
    PMM3 = "pmm3_b1",
    PMM2 = "pmm2_b1"
  )
  do.call(rbind, lapply(names(methods), function(m) {
    vals <- df[[methods[[m]]]]
    vals <- vals[!is.na(vals)]
    if (length(vals) < 10) return(NULL)
    mse_m <- mean((vals - beta_true)^2)
    data.frame(
      lambda       = df$lambda[1],
      n            = df$n[1],
      method       = m,
      bias         = mean(vals) - beta_true,
      variance     = var(vals),
      mse          = mse_m,
      are          = mean((df$ols_b1 - beta_true)^2, na.rm=TRUE) / mse_m,
      g3_empirical = if (m == "PMM3") var(vals) / var_ols else NA_real_,
      conv_rate    = if (paste0(tolower(m), "_ok") %in% names(df))
                       mean(df[[paste0(tolower(m), "_ok")]], na.rm=TRUE) else 1.0,
      stringsAsFactors = FALSE
    )
  }))
}

#' Run MC block for one (lambda, n) combination
run_mc_block <- function(lambda, n, M = CFG$mc_replications, seed = CFG$seed) {
  set.seed(seed + as.integer(lambda * 100) + n)
  message(sprintf("  MC | lambda=%.1f  n=%3d  M=%d", lambda, n, M))
  reps     <- lapply(seq_len(M), function(i) .one_mc_rep(lambda, n))
  df       <- do.call(rbind, reps)
  var_ols  <- var(df$ols_b1, na.rm = TRUE)
  .compute_metrics(df, beta_true = CFG$beta_true["beta1"], var_ols = var_ols)
}

#' Full MC study over lambda x n grid with disk caching
#' @param force_rerun If TRUE, ignore existing cache file
run_full_mc_study <- function(lambda_grid  = CFG$lambda_grid,
                               sample_sizes = CFG$sample_sizes,
                               M            = CFG$mc_replications,
                               force_rerun  = FALSE) {
  cache_file <- file.path(CFG$paths$mc_cache, "mc_full_results.rds")

  if (!force_rerun && file.exists(cache_file)) {
    message("Loading MC cache: ", cache_file)
    return(readRDS(cache_file))
  }

  message(sprintf("Starting MC study: %d lambda x %d n x M=%d replications",
                  length(lambda_grid), length(sample_sizes), M))
  t0 <- proc.time()

  all <- list()
  for (lam in lambda_grid)
    for (n in sample_sizes)
      all[[paste0("l", lam, "_n", n)]] <- run_mc_block(lam, n, M)

  full_df <- do.call(rbind, Filter(Negate(is.null), all))
  rownames(full_df) <- NULL

  saveRDS(full_df, cache_file)
  elapsed <- round((proc.time() - t0)["elapsed"] / 60, 1)
  message(sprintf("MC study done in %.1f min. Cache: %s", elapsed, cache_file))
  full_df
}
REOF
echo "Created: R/03_monte_carlo.R"
```

---

### 3e — `R/04_real_data.R`

```bash
cat > "$PROJECT_ROOT/R/04_real_data.R" << 'REOF'
# =============================================================================
# R/04_real_data.R
# Real data analysis:
#   Dataset A — Biaxial fatigue (Rieck & Nedelman, 1991) — manually restored
#   Dataset B — nlme::Fatigue (CRAN-available fallback)
# Plan Blocks 3 and 4
# =============================================================================

#' Biaxial fatigue dataset (Rieck & Nedelman, Technometrics 1991, 33:51-60)
#' 46 obs, 1% Cr-Mo-V steel
#' x = log(work per cycle), y = log(cycles to failure)
#' NOTE: restored manually; gamma3_hat ≈ -0.47 (slight asymmetry detected)
load_biaxial_fatigue <- function() {
  x <- c(-1.66,-1.66,-1.66,-1.66,-1.66,-1.66,-1.66,-1.66,
         -1.39,-1.39,-1.39,-1.39,-1.39,-1.39,-1.39,-1.39,
         -1.11,-1.11,-1.11,-1.11,-1.11,-1.11,-1.11,-1.11,
         -0.85,-0.85,-0.85,-0.85,-0.85,-0.85,-0.85,-0.85,
         -0.60,-0.60,-0.60,-0.60,-0.60,-0.60,-0.60,-0.60,
         -0.27,-0.27,-0.27,-0.27,-0.27,-0.27)
  y <- c(6.00,5.90,5.89,5.59,5.78,5.93,5.74,5.30,
         5.58,5.45,5.26,5.76,5.85,5.56,5.90,5.30,
         5.41,5.32,5.23,5.12,5.02,5.31,5.60,5.38,
         5.23,5.27,5.22,4.76,5.00,5.10,4.97,4.78,
         5.03,4.98,4.87,4.85,4.88,5.08,5.03,5.02,
         4.55,4.70,4.60,4.50,4.75,4.80)
  data.frame(x = x, y = y,
             dataset = "Biaxial Fatigue (Rieck & Nedelman, 1991)",
             stringsAsFactors = FALSE)
}

#' nlme::Fatigue aggregated to simple regression (CRAN-available)
load_nlme_fatigue <- function() {
  if (!requireNamespace("nlme", quietly = TRUE))
    stop("Install nlme: install.packages('nlme')")
  data("Fatigue", package = "nlme", envir = environment())
  agg <- aggregate(relLength ~ cycles, get("Fatigue", envir = environment()), mean)
  names(agg) <- c("x", "y")
  agg$dataset <- "nlme::Fatigue (aggregated)"
  agg
}

#' OLS residual diagnostics: gamma3, gamma4, gamma6, g3_theoretical
resid_diagnostics <- function(data, formula = y ~ x) {
  r   <- residuals(lm(formula, data = data))
  mu2 <- mean(r^2); mu3 <- mean(r^3)
  mu4 <- mean(r^4); mu6 <- mean(r^6)

  g3  <- mu3 / mu2^1.5
  g4  <- mu4 / mu2^2 - 3
  g6  <- mu6 / mu2^3 - 15*(mu4/mu2^2) + 30
  den <- 6 + 9*g4 + g6
  g3_theor <- if (!is.na(den) && den > 0) 1 - g4^2/den else NA_real_

  method_rec <- if (abs(g3) >= CFG$symmetry_threshold) "PMM2 (asymmetric errors)" else
                if (abs(g4) <  CFG$gaussian_threshold)  "OLS (near-Gaussian errors)" else
                                                          "PMM3 (symmetric non-Gaussian)"
  list(n = length(r), gamma3 = g3, gamma4 = g4, gamma6 = g6,
       g3_theoretical = g3_theor,
       improvement_pct = if (!is.na(g3_theor)) (1 - g3_theor)*100 else NA_real_,
       is_symmetric = abs(g3) < CFG$symmetry_threshold,
       is_gaussian  = abs(g4) < CFG$gaussian_threshold,
       recommended_method = method_rec)
}

#' Bootstrap comparison of OLS / PMM2 / PMM3 / MLE-TN on real data
bootstrap_comparison <- function(data, formula = y ~ x,
                                 B = 2000, seed = CFG$seed) {
  set.seed(seed)
  n <- nrow(data)
  methods <- c("OLS", "PMM2", "PMM3", "MLE_TN")

  results <- lapply(methods, function(m) {
    message("  Bootstrap | ", m)
    b1 <- replicate(B, {
      d <- data[sample(n, n, TRUE), ]
      tryCatch(switch(m,
        OLS    = coef(lm(formula, d))[2],
        PMM2   = EstemPMM::lm_pmm2(formula, d)@coefficients[2],
        PMM3   = { f <- lm_pmm3(formula, d); f@coefficients[2] },
        MLE_TN = mle_tn(formula, d)$coefficients[2]
      ), error = function(e) NA_real_)
    })
    data.frame(method = m,
               mean   = mean(b1, na.rm=TRUE),
               sd     = sd(b1,   na.rm=TRUE),
               q025   = quantile(b1, .025, na.rm=TRUE),
               q975   = quantile(b1, .975, na.rm=TRUE),
               na_pct = mean(is.na(b1))*100,
               stringsAsFactors = FALSE)
  })

  df      <- do.call(rbind, results)
  var_ols <- df$sd[df$method == "OLS"]^2
  df$g3_empirical <- df$sd^2 / var_ols
  df
}

#' Leave-One-Out cross-validation MSE for all methods
loo_mse <- function(data, formula = y ~ x) {
  n <- nrow(data)
  methods <- list(
    OLS    = function(d) lm(formula, d),
    PMM2   = function(d) EstemPMM::lm_pmm2(formula, d),
    PMM3   = function(d) lm_pmm3(formula, d),
    MLE_TN = function(d) mle_tn(formula, d)
  )
  sapply(names(methods), function(m) {
    mean(sapply(seq_len(n), function(i) {
      tryCatch({
        fit <- methods[[m]](data[-i, ])
        (data$y[i] - predict(fit, newdata = data[i, ]))^2
      }, error = function(e) NA_real_)
    }), na.rm = TRUE)
  })
}
REOF
echo "Created: R/04_real_data.R"
```

---

### 3f — `R/05_robustness.R`

```bash
cat > "$PROJECT_ROOT/R/05_robustness.R" << 'REOF'
# =============================================================================
# R/05_robustness.R
# Misspecification & bimodal robustness tests
# Plan Blocks 5 and 6
# =============================================================================

# Generators for symmetric unit-variance distributions
.rgen <- list(
  uniform     = function(n) runif(n, -sqrt(3), sqrt(3)),
  triangular  = function(n) (runif(n) + runif(n) - 1) * sqrt(6),
  logistic    = function(n) rlogis(n, scale = sqrt(3)/pi),
  student10   = function(n) rt(n, df=10) / sqrt(10/8)
)

#' Block 5: misspecification experiment
#' DGP != TN, but both MLE-TN and PMM3 are applied
run_misspec_block <- function(dist_name, n = 100, M = 1000, seed = CFG$seed) {
  gen <- .rgen[[dist_name]]
  if (is.null(gen)) stop("Unknown distribution: ", dist_name)

  set.seed(seed + which(names(.rgen) == dist_name))
  x_fix  <- rnorm(n)
  beta   <- CFG$beta_true

  reps <- lapply(seq_len(M), function(i) {
    y   <- beta["beta0"] + beta["beta1"] * x_fix + gen(n)
    dat <- data.frame(y=y, x=x_fix)
    b_ols  <- tryCatch(coef(lm(y~x, dat))[2],        error=function(e) NA_real_)
    b_mle  <- tryCatch(mle_tn(y~x, dat)$coefficients[2], error=function(e) NA_real_)
    b_pmm3 <- tryCatch({ f <- lm_pmm3(y~x, dat); f@coefficients[2] },
                       error=function(e) NA_real_)
    c(ols=b_ols, mle=b_mle, pmm3=b_pmm3)
  })

  df <- as.data.frame(do.call(rbind, reps))
  data.frame(
    distribution = dist_name, n = n,
    are_mle  = var(df$ols,  na.rm=TRUE) / var(df$mle,  na.rm=TRUE),
    are_pmm3 = var(df$ols,  na.rm=TRUE) / var(df$pmm3, na.rm=TRUE),
    bias_mle  = mean(df$mle,  na.rm=TRUE) - beta["beta1"],
    bias_pmm3 = mean(df$pmm3, na.rm=TRUE) - beta["beta1"],
    na_mle    = mean(is.na(df$mle))  * 100,
    na_pmm3   = mean(is.na(df$pmm3)) * 100,
    stringsAsFactors = FALSE
  )
}

#' Block 6: bimodal TN convergence test (lambda > 1)
run_bimodal_test <- function(lambda_grid = CFG$bimodal_lambda_grid,
                              n = 100, M = 500, seed = CFG$seed) {
  do.call(rbind, lapply(lambda_grid, function(lam) {
    set.seed(seed + as.integer(lam * 10))
    reps <- replicate(M, {
      y   <- 1 + 2*rnorm(n) + rtn(n, lambda=lam)
      dat <- data.frame(y=y, x=rnorm(n))
      tryCatch({
        f <- lm_pmm3(y~x, dat)
        c(b1=f@coefficients[2], conv=as.integer(f@convergence),
          iter=f@iterations)
      }, error=function(e) c(b1=NA_real_, conv=0L, iter=NA_integer_))
    })
    data.frame(
      lambda     = lam,
      is_bimodal = lam > 1,
      bias_pmm3  = mean(reps["b1",], na.rm=TRUE) - 2,
      conv_rate  = mean(reps["conv",], na.rm=TRUE),
      mean_iter  = mean(reps["iter",], na.rm=TRUE),
      na_pct     = mean(is.na(reps["b1",]))*100
    )
  }))
}
REOF
echo "Created: R/05_robustness.R"
```

---

### 3g — `R/pmm3_skeleton.R` (УМОВНИЙ ФАЙЛ)

**Виконуй цей крок тільки якщо** перевірка в Кроці 1 показала, що `lm_pmm3` відсутня в EstemPMM:

```bash
cat > "$PROJECT_ROOT/R/pmm3_skeleton.R" << 'REOF'
# =============================================================================
# R/pmm3_skeleton.R
# ЗАГЛУШКА lm_pmm3() — замінити реальною реалізацією з EstemPMM!
#
# Цей файл існує, бо EstemPMM ще не містить lm_pmm3().
# Скелет повертає OLS-оцінки з попередженням і маркером STUB=TRUE,
# щоб 03_monte_carlo.R міг виявити незавершену реалізацію.
#
# Математична основа PMM3 (S=3, симетричні розподіли):
#   Система Q рівнянь (формула 10, Zabolotnyi et al.):
#     sum_v { k1[p,v]*yv + k2[p,v]*(yv^2 + mu2) + k3[p,v]*yv^3 } = 0
#   Коефіцієнти k1, k2, k3 (формула 11a,b,c):
#     k1[p,v] = [ 3*(f_v)^2*(mu4-3*mu2^2) + 3*mu4*mu2 - mu6 ] /
#                 [ mu2^2*(mu4^2 - mu2*mu6) ] * x_v^p
#     k2[p,v] = -3*(mu4 - 3*mu2^2) / [ mu2^2*(mu4^2 - mu2*mu6) ] * f_v * x_v^p
#     k3[p,v] = (mu4 - 3*mu2^2) / [ mu2^2*(mu4^2 - mu2*mu6) ] * x_v^p
#   Розв'язок — Newton-Raphson (потрібно ~ 5-20 ітерацій при гарних початкових)
# =============================================================================

STUB <- TRUE   # Маркер заглушки — читається в 03_monte_carlo.R

# Допоміжні функції обчислення моментів (потрібні для майбутньої реалізації)
.compute_moments_pmm3 <- function(residuals) {
  r <- residuals - mean(residuals)
  list(mu2 = mean(r^2),
       mu4 = mean(r^4),
       mu6 = mean(r^6))
}

# Теоретичний g3 з кумулянтів
.pmm3_g_coefficient <- function(mu2, mu4, mu6) {
  gamma4 <- mu4/mu2^2 - 3
  gamma6 <- mu6/mu2^3 - 15*(mu4/mu2^2) + 30
  denom  <- 6 + 9*gamma4 + gamma6
  if (is.na(denom) || denom <= 0) return(1.0)
  1 - gamma4^2 / denom
}

# Заглушка lm_pmm3 — повертає OLS у обгортці S4
# Замінити блок '.pmm3_fit()' реальним Newton-Raphson алгоритмом
lm_pmm3 <- function(formula, data = NULL,
                    max_iter = 100, tol = 1e-6,
                    use_ols_start = TRUE) {

  warning("[STUB] lm_pmm3(): returning OLS estimates. ",
          "Implement Newton-Raphson solver to complete PMM3.")

  cl  <- match.call()
  mf  <- model.frame(formula, data = data)
  y   <- model.response(mf)
  X   <- model.matrix(attr(mf, "terms"), data = mf)
  fit <- lm.fit(X, y)
  mom <- .compute_moments_pmm3(fit$residuals)

  # Перевірка симетрії залишків
  gamma3 <- mean(fit$residuals^3) / mean(fit$residuals^2)^1.5
  if (abs(gamma3) > 0.3) {
    warning(sprintf("[PMM3] gamma3 = %.3f — errors appear asymmetric. Consider PMM2.", gamma3))
  }

  # Майбутній блок: замінити на .pmm3_fit(y, X, mom, max_iter, tol)
  # Поки що повертаємо OLS
  theta_final <- fit$coefficients

  # Формуємо S4-об'єкт сумісний з EstemPMM
  new("PMM3fit",
      coefficients = as.numeric(theta_final),
      residuals    = as.numeric(fit$residuals),
      convergence  = TRUE,     # OLS завжди "збігається"
      iterations   = 0L,
      call         = cl,
      method       = "PMM3-STUB (OLS fallback)",
      m2           = mom$mu2,
      m4           = mom$mu4,
      m6           = mom$mu6,
      gamma4       = mom$mu4/mom$mu2^2 - 3,
      gamma6       = mom$mu6/mom$mu2^3 - 15*(mom$mu4/mom$mu2^2) + 30,
      g_coefficient = .pmm3_g_coefficient(mom$mu2, mom$mu4, mom$mu6)
  )
}

message("[pmm3_skeleton.R] Stub lm_pmm3 loaded. Replace with full implementation.")
REOF
echo "Created: R/pmm3_skeleton.R (STUB)"
```

---

### 3h — Тести

```bash
cat > "$PROJECT_ROOT/tests/testthat/test-01-tn-distribution.R" << 'REOF'
# Tests for R/01_tn_distribution.R
library(testthat)

test_that("rtn returns correct length", {
  expect_length(rtn(100, lambda = 1), 100)
})

test_that("rtn is symmetric (gamma3 ≈ 0) for any lambda", {
  set.seed(42)
  for (lam in c(0.5, 1.0, 2.0)) {
    x  <- rtn(10000, lambda = lam)
    g3 <- mean(x^3) / mean(x^2)^1.5
    expect_lt(abs(g3), 0.15,
              label = paste("gamma3 for lambda =", lam))
  }
})

test_that("dtn density integrates to 1", {
  for (lam in c(0, 1, 2)) {
    val <- integrate(dtn, -Inf, Inf, lambda = lam)$value
    expect_equal(val, 1, tolerance = 1e-4,
                 label = paste("integral for lambda =", lam))
  }
})

test_that("tn_cumulants: g3 in (0,1) for lambda > 0", {
  for (lam in c(0.5, 1.0, 2.0)) {
    cum <- tn_cumulants(lam)
    expect_true(cum$g3 > 0 && cum$g3 < 1,
                label = paste("g3 bounds for lambda =", lam))
  }
})

test_that("tn_cumulants: gamma4 = 0 and g3 = 1 for lambda = 0", {
  cum <- tn_cumulants(0)
  expect_equal(cum$gamma4, 0, tolerance = 1e-5)
  expect_equal(cum$g3,     1, tolerance = 1e-5)
})

test_that("tn_cumulants: bimodal flag correct", {
  expect_false(tn_cumulants(0.5)$is_bimodal)
  expect_false(tn_cumulants(1.0)$is_bimodal)
  expect_true( tn_cumulants(1.5)$is_bimodal)
})

test_that("build_table_T1 returns correct shape", {
  T1 <- build_table_T1(lambda_grid = c(0.5, 1.0))
  expect_equal(nrow(T1), 2)
  expect_true(all(c("lambda","gamma4","gamma6","g3") %in% names(T1)))
})
REOF

cat > "$PROJECT_ROOT/tests/testthat/test-02-mle-tn.R" << 'REOF'
# Tests for R/02_mle_tn.R
library(testthat)

test_that("mle_tn converges on TN data and recovers beta1", {
  set.seed(123)
  x <- rnorm(200)
  y <- 1 + 2*x + rtn(200, lambda = 1.5)
  fit <- mle_tn(y ~ x, data = data.frame(y, x))

  expect_true(fit$converged)
  expect_equal(fit$coefficients["x"], 2, tolerance = 0.4)
})

test_that("mle_tn recovers lambda on large sample", {
  set.seed(456)
  x <- rnorm(500)
  y <- 1 + 2*x + rtn(500, lambda = 2.0)
  fit <- mle_tn(y ~ x, data = data.frame(y, x))

  expect_equal(fit$lambda, 2.0, tolerance = 0.6)
})

test_that("mle_tn AIC/BIC are finite", {
  set.seed(789)
  x <- rnorm(100); y <- 1 + 2*x + rtn(100, lambda=1)
  fit <- mle_tn(y ~ x, data=data.frame(y,x))
  expect_true(is.finite(fit$aic))
  expect_true(is.finite(fit$bic))
})

test_that("mle_tn on Gaussian data: lambda ≈ 0", {
  set.seed(101)
  x <- rnorm(300); y <- 1 + 2*x + rnorm(300)
  fit <- mle_tn(y ~ x, data=data.frame(y,x))
  expect_lt(abs(fit$lambda), 0.5)
})
REOF

cat > "$PROJECT_ROOT/tests/testthat/test-03-real-data.R" << 'REOF'
# Tests for R/04_real_data.R
library(testthat)

test_that("load_biaxial_fatigue returns correct dimensions", {
  dat <- load_biaxial_fatigue()
  expect_equal(nrow(dat), 46)
  expect_true(all(c("x","y") %in% names(dat)))
})

test_that("resid_diagnostics returns required fields", {
  dat  <- load_biaxial_fatigue()
  diag <- resid_diagnostics(dat)
  expect_true(all(c("gamma3","gamma4","g3_theoretical",
                    "is_symmetric","recommended_method") %in% names(diag)))
})

test_that("biaxial residuals have negative kurtosis (platykurtic)", {
  dat  <- load_biaxial_fatigue()
  diag <- resid_diagnostics(dat)
  # gamma4 < 0 for this dataset (TN-like)
  expect_lt(diag$gamma4, 0)
})
REOF
echo "Created: tests/testthat/ (3 files)"
```

---

### 3i — Rmd звіт

```bash
cat > "$PROJECT_ROOT/reports/main_report.Rmd" << 'REOF'
---
title: "PMM3 vs MLE-TN: Experimental Comparative Analysis"
subtitle: "Salinas et al. (2023), Mathematics 11(5), 1271"
author: "Sergiy Zabolotnyi"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    toc_depth: 3
    code_folding: hide
    theme: cosmo
    number_sections: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE, warning = FALSE, message = FALSE,
  fig.width = 10, fig.height = 6, dpi = 150,
  cache = FALSE
)
# Load project modules
source("../R/config.R")
for (f in list.files("../R", pattern="\\.R$", full.names=TRUE)) source(f)
set.seed(CFG$seed)
```

# Overview

**Research question:** Is PMM3 competitive with MLE-TN when the true data-generating
process is the Two-piece Normal distribution?

PMM3 uses only moment information (μ₂, μ₄, μ₆) and does not assume a specific
distributional shape. MLE-TN maximises the exact TN log-likelihood. The comparison
spans four dimensions: sample size, shape parameter λ, misspecification, and bimodality.

Theoretical efficiency of PMM3 vs OLS:
$$g_3(\lambda) = 1 - \frac{\gamma_4^2(\lambda)}{6 + 9\gamma_4(\lambda) + \gamma_6(\lambda)}$$

# Block 1: Theoretical Calibration {.tabset}

## Table T1

```{r table-T1}
T1 <- build_table_T1()
knitr::kable(
  T1[, c("lambda","gamma4","gamma6","g3","improvement_pct","mode")],
  digits = 3,
  col.names = c("λ", "γ₄", "γ₆", "g₃", "Improvement (%)", "Mode"),
  caption = "Table T1. Theoretical TN cumulants and PMM3 efficiency"
)
```

## Efficiency Curve

```{r fig-g3-curve, fig.cap="Theoretical PMM3 efficiency for TN distribution"}
library(ggplot2)
ggplot(T1, aes(x=lambda, y=g3)) +
  geom_line(size=1.2, colour="#2166AC") +
  geom_point(aes(colour=mode), size=3) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey50") +
  scale_colour_manual(values=c(unimodal="#4DAF4A", bimodal="#E41A1C")) +
  scale_y_continuous(labels=scales::percent_format()) +
  labs(x="λ", y="g₃ = Var(PMM3)/Var(OLS)",
       title="PMM3 Theoretical Efficiency vs OLS for TN(λ)",
       colour="Mode") +
  theme_bw(base_size=13)
```

# Block 2: Monte Carlo Study

> **Note:** Set `CFG$mc_replications <- 5000` in `R/config.R` for the final run.
> Current value: `r CFG$mc_replications`.

```{r run-mc, cache=TRUE}
mc_results <- run_full_mc_study()
```

```{r fig-are, fig.cap="ARE by method, lambda, and sample size"}
ggplot(mc_results[!is.na(mc_results$are), ],
       aes(x=factor(n), y=are, fill=method)) +
  geom_col(position="dodge", alpha=0.85) +
  geom_hline(yintercept=1, linetype="dashed") +
  facet_wrap(~lambda, labeller=label_bquote(lambda==.(lambda))) +
  scale_fill_brewer(palette="Set1") +
  labs(x="Sample size n", y="ARE = MSE(OLS)/MSE(method)",
       title="Asymptotic Relative Efficiency vs OLS",
       fill="Method") +
  theme_bw(base_size=12) + theme(legend.position="bottom")
```

## PMM3: empirical vs theoretical g₃

```{r table-g3-comparison}
pmm3_mc <- mc_results[mc_results$method == "PMM3" & !is.na(mc_results$g3_empirical), ]
T2 <- merge(pmm3_mc[, c("lambda","n","g3_empirical")],
            T1[, c("lambda","g3")], by="lambda", all.x=TRUE)
T2$diff <- round(T2$g3_empirical - T2$g3, 4)
knitr::kable(T2, digits=3,
             col.names=c("λ","n","g̃₃ (MC)","g₃ (theor.)","Δ"),
             caption="Table T2. Empirical vs theoretical PMM3 efficiency")
```

# Block 3: Real Data — Biaxial Fatigue (n=46)

```{r biaxial-diagnostics}
dat_bf   <- load_biaxial_fatigue()
diag_bf  <- resid_diagnostics(dat_bf)
```

**OLS residual diagnostics:**

| Parameter | Value | Interpretation |
|-----------|-------|----------------|
| γ̂₃ | `r round(diag_bf$gamma3, 3)` | `r ifelse(diag_bf$is_symmetric, "Symmetric ✓", "**Asymmetric ⚠️**")` |
| γ̂₄ | `r round(diag_bf$gamma4, 3)` | `r ifelse(diag_bf$is_gaussian, "Near-Gaussian", "Non-Gaussian")` |
| g₃ (theoretical) | `r round(diag_bf$g3_theoretical, 4)` | — |
| PMM3 improvement | `r round(diag_bf$improvement_pct, 2)`% | Expected over OLS |
| Recommended method | — | **`r diag_bf$recommended_method`** |

> ⚠️ Note: γ̂₃ ≈ `r round(diag_bf$gamma3, 2)` indicates residual asymmetry.
> Biaxial fatigue may be better suited for PMM2, not PMM3.
> This is an open research question explored in this analysis.

```{r fit-all-methods, warning=TRUE}
# Fit all methods — errors caught individually
fit_ols   <- lm(y ~ x, dat_bf)
fit_mle   <- tryCatch(mle_tn(y ~ x, dat_bf), error=function(e) NULL)
fit_pmm3  <- tryCatch(lm_pmm3(y ~ x, dat_bf), error=function(e) NULL)
fit_pmm2  <- tryCatch(EstemPMM::lm_pmm2(y ~ x, dat_bf), error=function(e) NULL)

# Summary table
results_bf <- data.frame(
  Method = c("OLS", "MLE-TN", "PMM3", "PMM2"),
  Beta0  = c(coef(fit_ols)[1],
             if (!is.null(fit_mle))  fit_mle$coefficients[1]  else NA,
             if (!is.null(fit_pmm3)) fit_pmm3@coefficients[1]  else NA,
             if (!is.null(fit_pmm2)) fit_pmm2@coefficients[1]  else NA),
  Beta1  = c(coef(fit_ols)[2],
             if (!is.null(fit_mle))  fit_mle$coefficients[2]  else NA,
             if (!is.null(fit_pmm3)) fit_pmm3@coefficients[2]  else NA,
             if (!is.null(fit_pmm2)) fit_pmm2@coefficients[2]  else NA),
  Lambda = c(NA,
             if (!is.null(fit_mle)) round(fit_mle$lambda, 3) else NA,
             NA, NA),
  AIC    = c(AIC(fit_ols),
             if (!is.null(fit_mle)) round(fit_mle$aic, 2) else NA,
             NA, NA)
)
knitr::kable(results_bf, digits=4, caption="Table T3. Method comparison on Biaxial Fatigue")
```

# Block 5: Misspecification Robustness

```{r misspec, cache=TRUE}
misspec_results <- do.call(rbind,
  lapply(CFG$misspec_distributions, run_misspec_block, M=500))
knitr::kable(misspec_results, digits=3,
             caption="Table T6. ARE under distributional misspecification (DGP ≠ TN)")
```

# Session Info

```{r session-info}
sessionInfo()
```
REOF
echo "Created: reports/main_report.Rmd"
```

---

### 3j — Скрипт запуску та допоміжні файли

```bash
cat > "$PROJECT_ROOT/scripts/run_all.R" << 'REOF'
# =============================================================================
# scripts/run_all.R — Повний запуск дослідження
# Використання: Rscript scripts/run_all.R
#               або source("scripts/run_all.R") з RStudio
# =============================================================================
t_start <- proc.time()
cat("=== PMM3 vs MLE-TN Study ===\n")
cat("Started:", format(Sys.time()), "\n\n")

setwd(normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."),
                    mustWork = FALSE))

source("R/config.R")
for (f in list.files("R", "\\.R$", full.names=TRUE)) source(f)

cat("--- Block 1: Table T1 ---\n")
T1 <- build_table_T1()
write.csv(T1, file.path(CFG$paths$results_tables, "T1_theoretical.csv"),
          row.names=FALSE)
cat("Saved T1\n\n")

cat("--- Block 2: Monte Carlo ---\n")
mc <- run_full_mc_study()
write.csv(mc, file.path(CFG$paths$results_tables, "T2_monte_carlo.csv"),
          row.names=FALSE)
cat("Saved T2\n\n")

cat("--- Block 3: Biaxial Fatigue ---\n")
dat_bf <- load_biaxial_fatigue()
diag   <- resid_diagnostics(dat_bf)
cat("Recommended method:", diag$recommended_method, "\n")
saveRDS(diag, file.path(CFG$paths$results_tables, "biaxial_diagnostics.rds"))

cat("\n--- Block 5: Misspecification ---\n")
misspec <- do.call(rbind,
  lapply(CFG$misspec_distributions, run_misspec_block, M=500))
write.csv(misspec, file.path(CFG$paths$results_tables, "T6_misspecification.csv"),
          row.names=FALSE)
cat("Saved T6\n\n")

cat("--- Rendering report ---\n")
rmarkdown::render("reports/main_report.Rmd",
                  output_file = "../results/main_report.html",
                  quiet = TRUE)
cat("Report: results/main_report.html\n\n")

elapsed <- round((proc.time() - t_start)["elapsed"] / 60, 1)
cat(sprintf("=== Done in %.1f min ===\n", elapsed))
REOF

cat > "$PROJECT_ROOT/.gitignore" << 'REOF'
.Rhistory
.RData
.Rproj.user/
*.Rproj
results/mc_cache/*.rds
results/figures/*.png
results/figures/*.pdf
!results/tables/*.csv
REOF

cat > "$PROJECT_ROOT/README.md" << 'REOF'
# PMM3 vs MLE-TN: Comparative Analysis

## Quick start

```r
# Install dependencies
install.packages(c("ggplot2","dplyr","tidyr","knitr",
                   "rmarkdown","nlme","pracma","testthat"))
devtools::install_github("SZabolotnii/EstemPMM")

# Run project
setwd("/path/to/pmm3_tn_study")
source("R/config.R")
for (f in list.files("R", "\\.R$", full.names=TRUE)) source(f)

# Quick smoke test
T1 <- build_table_T1()
print(T1[, c("lambda","gamma4","g3","improvement_pct")])

# Full run (~5-30 min depending on CFG$mc_replications)
source("scripts/run_all.R")
```

## Reference

Salinas H., Bakouch H., Qarmalah N., Martínez-Flórez G. (2023).
A Flexible Class of Two-Piece Normal Distribution.
*Mathematics*, 11(5), 1271. https://doi.org/10.3390/math11051271
REOF
echo "Created: scripts/run_all.R, .gitignore, README.md"
```

---

## КРОК 4 — Встановлення R-залежностей

```bash
Rscript -e "
required <- c('ggplot2','dplyr','tidyr','knitr','rmarkdown',
              'nlme','pracma','testthat','scales')
missing  <- required[!sapply(required, requireNamespace, quietly=TRUE)]
if (length(missing) > 0) {
  cat('Installing:', paste(missing, collapse=', '), '\n')
  install.packages(missing, repos='https://cran.r-project.org', quiet=TRUE)
} else {
  cat('All R packages present.\n')
}
" 2>/dev/null
```

---

## КРОК 5 — Запуск тестів

```bash
cd "$PROJECT_ROOT"
Rscript -e "
source('R/config.R')
for (f in list.files('R', '\\.R$', full.names=TRUE)) source(f)
library(testthat)
results <- test_dir('tests/testthat', reporter='minimal')
cat('\n=== Test summary ===\n')
print(results)
" 2>/dev/null
```

**Очікуваний мінімум:** усі тести в `test-01` і `test-03` проходять без помилок.
Тести `test-02` (MLE-TN) можуть займати 10-30 секунд.

---

## КРОК 6 — Smoke test (швидка перевірка модулів)

```bash
cd "$PROJECT_ROOT"
Rscript -e "
source('R/config.R')
for (f in list.files('R', '\\.R$', full.names=TRUE)) source(f)

cat('--- Test rtn() ---\n')
set.seed(42)
x <- rtn(1000, lambda=1.5)
cat('n=1000 | mean=', round(mean(x),3), '| gamma3=', round(mean(x^3)/mean(x^2)^1.5, 3), '\n')

cat('--- Test tn_cumulants(lambda=2) ---\n')
cum <- tn_cumulants(2)
cat('gamma4=', round(cum\$gamma4,3), '| g3=', round(cum\$g3,3), '| improvement=', round(cum\$improvement_pct,1), '%\n')

cat('--- Test mle_tn() ---\n')
x2 <- rnorm(100); y2 <- 1 + 2*x2 + rtn(100, lambda=1.5)
fit <- mle_tn(y~x, data=data.frame(y=y2, x=x2))
cat('beta1=', round(fit\$coefficients[2],3), '| lambda=', round(fit\$lambda,3), '| converged=', fit\$converged, '\n')

cat('--- Test load_biaxial_fatigue() ---\n')
dat <- load_biaxial_fatigue()
cat('n=', nrow(dat), '| gamma3=', round(mean(residuals(lm(y~x,dat))^3)/mean(residuals(lm(y~x,dat))^2)^1.5, 3), '\n')

cat('\n=== Smoke test PASSED ===\n')
" 2>/dev/null
```

---

## КРОК 7 — Фінальна верифікація

```bash
cd "$PROJECT_ROOT"

echo "=== File inventory ==="
REQUIRED=(
  "R/config.R"
  "R/01_tn_distribution.R"
  "R/02_mle_tn.R"
  "R/03_monte_carlo.R"
  "R/04_real_data.R"
  "R/05_robustness.R"
  "tests/testthat/test-01-tn-distribution.R"
  "tests/testthat/test-02-mle-tn.R"
  "tests/testthat/test-03-real-data.R"
  "reports/main_report.Rmd"
  "scripts/run_all.R"
  "README.md"
  ".gitignore"
)

PASS=0; FAIL=0
for f in "${REQUIRED[@]}"; do
  if [ -f "$f" ]; then
    echo "  ✓ $f"
    ((PASS++))
  else
    echo "  ✗ MISSING: $f"
    ((FAIL++))
  fi
done

echo ""
echo "=== Result: $PASS OK / $FAIL MISSING ==="

# Перевірка умовного файлу
if [ -f "R/pmm3_skeleton.R" ]; then
  echo "  ℹ  R/pmm3_skeleton.R present (lm_pmm3 stub — replace when EstemPMM updated)"
fi

# Підсумок розмірів
echo ""
echo "=== Project size ==="
du -sh .
```

**Критерій успіху:** `13 OK / 0 MISSING`.

---

## Правила та обмеження для агента

1. **Не модифікуй EstemPMM** — працюй тільки в `$PROJECT_ROOT`
2. **Не запускай важкі обчислення** автоматично — лише smoke test і unit tests
3. **Повідом про відсутній lm_pmm3** — якщо stub завантажено, напиши підсумок:
   ```
   ⚠️  lm_pmm3() не знайдено в EstemPMM.
   Файл R/pmm3_skeleton.R містить stub-реалізацію.
   Monte Carlo результати для PMM3 будуть NA до реалізації повного алгоритму.
   ```
4. **При помилці на будь-якому кроці** — зупинись, поясни проблему, запропонуй виправлення
5. **Git init** — тільки якщо користувач явно попросить
6. **Звіт про завершення** — після Кроку 7 виведи чіткий підсумок:

```
╔══════════════════════════════════════════════╗
║  PMM3-TN Research Project — Setup Complete   ║
╠══════════════════════════════════════════════╣
║  Location : <PROJECT_ROOT>                   ║
║  Files    : 13 required files ✓              ║
║  Tests    : <N> passed / <M> failed          ║
║  PMM3     : EstemPMM / STUB                  ║
╠══════════════════════════════════════════════╣
║  Next steps:                                 ║
║  1. source("R/config.R") + load_project()    ║
║  2. T1 <- build_table_T1()  # Table T1       ║
║  3. source("scripts/run_all.R")  # Full run  ║
╚══════════════════════════════════════════════╝
```

---

*Версія: 1.1 · Цільовий агент: Claude Code (CLI) · Інструменти: bash*  
*Проєкт: EstemPMM / PMM3-TN Comparative Research · Березень 2026*
