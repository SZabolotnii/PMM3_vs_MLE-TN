# Інструкція для агента: розгортання дослідницького проєкту PMM3-TN

**Призначення:** Ця інструкція повністю описує, що саме повинен зробити агент для розгортання нового R-проєкту порівняльного аналізу PMM3 vs. MLE-TN.

**Контекст:** Дослідження проводиться в рамках розробки пакету `EstemPMM`. Референсна стаття — Salinas et al. (2023), Mathematics 11(5), 1271. Теоретичний план — `PMM3_vs_MLTN_comparison_plan.md`.

---

## КРОК 0. Отримання цільової директорії від користувача

Перш ніж щось робити, **запитати у користувача**:

```
Будь ласка, вкажіть повний шлях до директорії,
де потрібно створити проєкт.
Наприклад: /home/user/projects/pmm3_tn_study
```

Зберегти відповідь як `PROJECT_ROOT`. Всі подальші шляхи вказуються відносно `PROJECT_ROOT`.

---

## КРОК 1. Створення структури директорій

Виконати наступні команди (bash або еквівалент):

```bash
mkdir -p $PROJECT_ROOT
cd $PROJECT_ROOT

# Основні директорії
mkdir -p R                      # R-код модулів
mkdir -p tests/testthat         # Unit-тести
mkdir -p data/raw               # Вхідні дані
mkdir -p data/processed         # Оброблені дані
mkdir -p results/tables         # Таблиці результатів (CSV/RDS)
mkdir -p results/figures        # Графіки (PNG/PDF)
mkdir -p results/mc_cache       # Кеш Monte Carlo симуляцій
mkdir -p reports                # Звіти (Rmd → HTML/PDF)
mkdir -p docs                   # Технічна документація
mkdir -p scripts                # Скрипти запуску
```

**Перевірка після створення:**
```bash
find $PROJECT_ROOT -type d | sort
# Повинно вивести всі 10 директорій
```

---

## КРОК 2. Ініціалізація R-проєкту

### 2.1. Файл `DESCRIPTION` (метадані проєкту)

Створити файл `$PROJECT_ROOT/DESCRIPTION`:

```
Title: PMM3 vs MLE-TN Comparative Analysis
Type: Research Project
Version: 0.1.0
Date: 2026-03-12
Authors: Sergiy Zabolotnyi <s.zabolotnyi@research.ua>
Description: Experimental comparison of Polynomial Maximization Method (PMM3)
    with Maximum Likelihood Estimation for the Two-piece Normal Distribution
    (MLE-TN), as proposed by Salinas et al. (2023).
Reference: Salinas H., Bakouch H., Qarmalah N., Martinez-Florez G. (2023).
    A Flexible Class of Two-Piece Normal Distribution with a Regression
    Illustration to Biaxial Fatigue Data. Mathematics, 11(5), 1271.
    https://doi.org/10.3390/math11051271
R-Depends: R (>= 4.1.0)
Imports:
    EstemPMM,
    MASS,
    nlme,
    ggplot2,
    dplyr,
    tidyr,
    knitr,
    rmarkdown,
    pracma
Suggests:
    testthat (>= 3.0.0),
    parallel
License: GPL-3
Encoding: UTF-8
```

### 2.2. Файл `.Rprofile` (налаштування середовища)

Створити файл `$PROJECT_ROOT/.Rprofile`:

```r
# .Rprofile — автоматичне налаштування при відкритті проєкту

# Встановлення робочої директорії
if (file.exists("DESCRIPTION")) {
  local({
    proj_root <- normalizePath(".", mustWork = TRUE)
    message("PMM3-TN Research Project: ", proj_root)
  })
}

# Зручні аліаси для шляхів
.paths <- list(
  root    = ".",
  R       = "R",
  data    = "data",
  results = "results",
  reports = "reports",
  cache   = "results/mc_cache"
)

# Функція для завантаження всіх модулів проєкту
load_project <- function(verbose = TRUE) {
  r_files <- list.files(.paths$R, pattern = "\\.R$", full.names = TRUE)
  for (f in r_files) {
    if (verbose) message("Loading: ", basename(f))
    source(f)
  }
  if (verbose) message("Project loaded. ", length(r_files), " modules.")
  invisible(r_files)
}

# Підключення EstemPMM якщо доступний
if (requireNamespace("EstemPMM", quietly = TRUE)) {
  suppressMessages(library(EstemPMM))
} else {
  message("NOTE: EstemPMM not installed. Run: devtools::install_github('SZabolotnii/EstemPMM')")
}

message("Use load_project() to load all R modules.")
```

### 2.3. Файл `config.R` (конфігурація дослідження)

Створити файл `$PROJECT_ROOT/R/config.R`:

```r
# ===========================================================================
# R/config.R — Централізована конфігурація дослідження PMM3 vs MLE-TN
# ===========================================================================

# ---------------------------------------------------------------------------
# 1. ПАРАМЕТРИ MONTE CARLO
# ---------------------------------------------------------------------------

CFG <- list(

  # Кількість реплікацій (5000 для фінального, 500 для тестування)
  mc_replications = 5000,

  # Розміри вибірок
  # n=46: відповідає biaxial fatigue з Salinas et al. (2023)
  sample_sizes = c(46, 100, 200, 500),

  # Параметри форми TN-розподілу
  # λ=0: нормальний; λ=1: межа уні/бімодальності; λ>1: бімодальний
  lambda_grid = c(0.5, 1.0, 1.5, 2.0, 3.0, 5.0),

  # Параметри лінійної регресії (справжні значення)
  beta_true = c(beta0 = 1.0, beta1 = 2.0),

  # Параметр масштабу TN для регресії
  eta_true = 1.0,

  # ---------------------------------------------------------------------------
  # 2. ПАРАМЕТРИ ЕФЕКТИВНОСТІ
  # ---------------------------------------------------------------------------

  # Поріг для визначення симетричності (|gamma3| < threshold)
  symmetry_threshold = 0.3,

  # Поріг для визначення гаусовості (|gamma4| < threshold)
  gaussian_threshold = 0.1,

  # ---------------------------------------------------------------------------
  # 3. РОЗПОДІЛИ ДЛЯ MISSPECIFICATION ТЕСТУ
  # ---------------------------------------------------------------------------
  # (Блок 5 плану — robustness до хибної специфікації)

  misspec_distributions = list(
    uniform     = list(name = "Uniform(-√3, √3)", gamma4 = -1.20, gamma6 = 6.9),
    trapezoidal = list(name = "Trapezoidal(β=0.75)", gamma4 = -1.10, gamma6 = 6.4),
    triangular  = list(name = "Triangular",          gamma4 = -0.60, gamma6 = 1.7),
    logistic    = list(name = "Logistic",             gamma4 =  1.20, gamma6 = NA),
    student10   = list(name = "Student-t(df=10)",     gamma4 =  1.00, gamma6 = NA)
  ),

  # ---------------------------------------------------------------------------
  # 4. ПАРАМЕТРИ БІМОДАЛЬНОГО ТЕСТУ
  # ---------------------------------------------------------------------------
  # (Блок 6 плану — convergence при λ > 1)

  bimodal_lambda_grid = c(0.8, 1.0, 1.2, 1.5, 2.0, 3.0),

  # ---------------------------------------------------------------------------
  # 5. ШЛЯХИ
  # ---------------------------------------------------------------------------

  paths = list(
    results_tables  = "results/tables",
    results_figures = "results/figures",
    mc_cache        = "results/mc_cache",
    data_raw        = "data/raw",
    data_processed  = "data/processed"
  ),

  # ---------------------------------------------------------------------------
  # 6. ВІДТВОРЮВАНІСТЬ
  # ---------------------------------------------------------------------------

  seed = 20260312,       # Дата початку дослідження (YYYYMMDD)
  use_parallel = FALSE,  # TRUE = використовувати parallel::mclapply
  n_cores = 2            # Кількість ядер при паралельних обчисленнях
)

message("CFG loaded: ", CFG$mc_replications, " MC replications, ",
        length(CFG$lambda_grid), " lambda values, ",
        length(CFG$sample_sizes), " sample sizes.")
```

---

## КРОК 3. Створення R-модулів

### 3.1. Модуль `R/01_tn_distribution.R` — визначення TN-розподілу

Створити файл `$PROJECT_ROOT/R/01_tn_distribution.R`:

```r
# ===========================================================================
# R/01_tn_distribution.R
# Два-частинний нормальний розподіл TN(ξ, η, λ)
# Salinas et al. (2023), Mathematics 11(5), 1271
# ===========================================================================

# ---------------------------------------------------------------------------
# Генератор випадкових чисел TN(ξ, η, λ)
# Алгоритм: Section 4.2 статті
# ---------------------------------------------------------------------------
rtn <- function(n, xi = 0, eta = 1, lambda = 1) {
  stopifnot(eta > 0, lambda >= 0)
  U <- rnorm(n, mean = lambda, sd = 1)       # U ~ N(λ, 1)
  Y <- abs(U)                                 # Y ~ FN(λ, 1)
  S <- sample(c(-1L, 1L), n, replace = TRUE) # S ~ Bernoulli(1/2)
  Z <- S * Y                                  # Z ~ TN(λ)
  xi + eta * Z
}

# ---------------------------------------------------------------------------
# Щільність TN (формула 7 зі статті)
# ---------------------------------------------------------------------------
dtn <- function(x, xi = 0, eta = 1, lambda = 0, log = FALSE) {
  stopifnot(eta > 0, lambda >= 0)
  z <- (x - xi) / eta
  log_dens <- dnorm(z, log = TRUE) - log(eta) -
               lambda^2 / 2 + log(cosh(lambda * z))
  if (log) log_dens else exp(log_dens)
}

# ---------------------------------------------------------------------------
# CDF TN (формула 3 зі статті)
# ---------------------------------------------------------------------------
ptn <- function(q, xi = 0, eta = 1, lambda = 0) {
  z <- (q - xi) / eta
  (pnorm(z - lambda) + pnorm(z + lambda)) / 2
}

# ---------------------------------------------------------------------------
# Аналітичні моменти TN(0, 1, λ)
# E(Z^{2k}) через гіпергеометричну функцію Куммера 1F1
# Потребує: pracma::hypergeo або власна реалізація
# ---------------------------------------------------------------------------
tn_raw_moment_2k <- function(k, lambda) {
  # E(Z^{2k}) = 2^{k-1/2}/sqrt(π) * Γ(k+1/2) * 1F1(-k; 1/2; -λ²/2)
  # Обчислення через числове інтегрування (надійніше ніж 1F1 при великих k)
  integrand <- function(z) z^(2*k) * dtn(z, lambda = lambda)
  result <- tryCatch(
    integrate(integrand, -Inf, Inf, subdivisions = 500L)$value,
    error = function(e) NA_real_
  )
  result
}

# ---------------------------------------------------------------------------
# Кумулянтні коефіцієнти TN(0, 1, λ)
# ---------------------------------------------------------------------------
tn_cumulants <- function(lambda) {
  mu2 <- tn_raw_moment_2k(1, lambda)
  mu4 <- tn_raw_moment_2k(2, lambda)
  mu6 <- tn_raw_moment_2k(3, lambda)

  gamma4 <- mu4 / mu2^2 - 3
  gamma6 <- mu6 / mu2^3 - 15 * (mu4 / mu2^2) + 30

  # Теоретичний коефіцієнт зменшення дисперсії PMM3 vs OLS
  denom <- 6 + 9 * gamma4 + gamma6
  g3 <- if (denom > 0) 1 - gamma4^2 / denom else NA_real_

  list(
    lambda = lambda,
    mu2    = mu2,
    mu4    = mu4,
    mu6    = mu6,
    gamma4 = gamma4,
    gamma6 = gamma6,
    g3     = g3,
    is_bimodal = lambda > 1
  )
}

# ---------------------------------------------------------------------------
# Таблиця T1: теоретичні значення для сітки λ
# ---------------------------------------------------------------------------
build_table_T1 <- function(lambda_grid = CFG$lambda_grid) {
  message("Building Table T1: theoretical TN cumulants...")
  rows <- lapply(lambda_grid, tn_cumulants)
  T1 <- do.call(rbind, lapply(rows, as.data.frame))
  T1$improvement_pct <- round((1 - T1$g3) * 100, 2)
  T1$mode <- ifelse(T1$is_bimodal, "bimodal", "unimodal")
  rownames(T1) <- NULL
  message("Table T1 ready: ", nrow(T1), " rows.")
  T1
}
```

### 3.2. Модуль `R/02_mle_tn.R` — MLE для TN-розподілу

Створити файл `$PROJECT_ROOT/R/02_mle_tn.R`:

```r
# ===========================================================================
# R/02_mle_tn.R
# Maximum Likelihood Estimation для TN-регресії
# Конкурент для PMM3 у порівняльному аналізі
# ===========================================================================

# ---------------------------------------------------------------------------
# Log-likelihood TN регресії (формула 9 зі статті)
# par = c(beta_0, ..., beta_p, log_eta, lambda)
# ---------------------------------------------------------------------------
.ll_tn_regression <- function(par, y, X, negative = TRUE) {
  p    <- ncol(X)
  beta <- par[seq_len(p)]
  eta  <- exp(par[p + 1])       # exp() для гарантії позитивності
  lam  <- par[p + 2]

  if (eta <= 0 || lam < 0) return(if (negative) Inf else -Inf)

  mu  <- as.vector(X %*% beta)
  z   <- (y - mu) / eta
  n   <- length(y)

  ll <- -n * log(2 * pi) / 2 - n * log(eta) -
        n * lam^2 / 2 - sum(z^2) / 2 + sum(log(cosh(lam * z)))

  if (negative) -ll else ll
}

# ---------------------------------------------------------------------------
# Основна функція MLE-TN
# ---------------------------------------------------------------------------
mle_tn <- function(formula, data = NULL,
                   lambda_init = 1.0,
                   method = "L-BFGS-B",
                   maxit = 500) {

  cl  <- match.call()
  mf  <- model.frame(formula, data = data)
  y   <- model.response(mf)
  X   <- model.matrix(attr(mf, "terms"), data = mf)
  n   <- length(y)
  p   <- ncol(X)

  # OLS як стартові значення для beta та eta
  ols     <- lm.fit(X, y)
  beta0   <- coef(ols)
  eta0    <- log(sd(residuals(ols)))  # log-scale
  par0    <- c(beta0, eta0, lambda_init)

  # Межі для L-BFGS-B
  lower <- c(rep(-Inf, p), -5,  0.001)
  upper <- c(rep( Inf, p),  5, 15.0)

  opt <- optim(
    par     = par0,
    fn      = .ll_tn_regression,
    y       = y, X = X,
    negative = TRUE,
    method  = method,
    lower   = lower,
    upper   = upper,
    control = list(maxit = maxit, factr = 1e7),
    hessian = TRUE
  )

  # Виявлення збіжності
  converged <- opt$convergence == 0

  # Стандартні похибки через Hessian
  se <- tryCatch({
    H_inv <- solve(opt$hessian)
    sqrt(pmax(diag(H_inv), 0))
  }, error = function(e) rep(NA_real_, length(par0)))

  list(
    call       = cl,
    coefficients = opt$par[seq_len(p)],
    eta        = exp(opt$par[p + 1]),
    lambda     = opt$par[p + 2],
    loglik     = -opt$value,
    se_beta    = se[seq_len(p)],
    se_eta     = se[p + 1] * exp(opt$par[p + 1]),
    se_lambda  = se[p + 2],
    converged  = converged,
    aic        = 2 * (p + 2) + 2 * opt$value,
    bic        = log(n) * (p + 2) + 2 * opt$value,
    n_params   = p + 2,
    n          = n,
    fitted     = as.vector(X %*% opt$par[seq_len(p)]),
    residuals  = y - as.vector(X %*% opt$par[seq_len(p)])
  )
}

# Обгортка для coef()
coef.mle_tn <- function(object, ...) object$coefficients
```

### 3.3. Модуль `R/03_monte_carlo.R` — Monte Carlo порівняння

Створити файл `$PROJECT_ROOT/R/03_monte_carlo.R`:

```r
# ===========================================================================
# R/03_monte_carlo.R
# Monte Carlo порівняння: PMM3 vs MLE-TN vs OLS vs PMM2
# Відповідає Блоку 2 плану дослідження
# ===========================================================================

# ---------------------------------------------------------------------------
# Одна Monte Carlo реплікація
# ---------------------------------------------------------------------------
.one_mc_rep <- function(lambda, n, beta = CFG$beta_true, eta = CFG$eta_true) {

  # Генерація даних з TN-помилками
  x   <- rnorm(n)
  eps <- rtn(n, xi = 0, eta = eta, lambda = lambda)
  y   <- beta["beta0"] + beta["beta1"] * x + eps
  dat <- data.frame(y = y, x = x)

  result <- list(lambda = lambda, n = n)

  # OLS
  tryCatch({
    fit_ols <- lm(y ~ x, data = dat)
    result$ols_b0 <- coef(fit_ols)[1]
    result$ols_b1 <- coef(fit_ols)[2]
  }, error = function(e) {
    result$ols_b0 <<- NA; result$ols_b1 <<- NA
  })

  # MLE-TN
  tryCatch({
    fit_mle <- mle_tn(y ~ x, data = dat)
    result$mle_b0 <- fit_mle$coefficients[1]
    result$mle_b1 <- fit_mle$coefficients[2]
    result$mle_converged <- fit_mle$converged
  }, error = function(e) {
    result$mle_b0 <<- NA; result$mle_b1 <<- NA
    result$mle_converged <<- FALSE
  })

  # PMM3 (через EstemPMM)
  tryCatch({
    fit_pmm3 <- EstemPMM::lm_pmm3(y ~ x, data = dat)
    result$pmm3_b0 <- fit_pmm3@coefficients[1]
    result$pmm3_b1 <- fit_pmm3@coefficients[2]
    result$pmm3_converged <- fit_pmm3@convergence
    result$pmm3_g <- fit_pmm3@g_coefficient
  }, error = function(e) {
    result$pmm3_b0 <<- NA; result$pmm3_b1 <<- NA
    result$pmm3_converged <<- FALSE; result$pmm3_g <<- NA
  })

  # PMM2 (контроль — має бути ≈ OLS при симетричних даних)
  tryCatch({
    fit_pmm2 <- EstemPMM::lm_pmm2(y ~ x, data = dat)
    result$pmm2_b0 <- fit_pmm2@coefficients[1]
    result$pmm2_b1 <- fit_pmm2@coefficients[2]
  }, error = function(e) {
    result$pmm2_b0 <<- NA; result$pmm2_b1 <<- NA
  })

  as.data.frame(result, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Повний Monte Carlo для однієї пари (λ, n)
# ---------------------------------------------------------------------------
run_mc_block <- function(lambda, n, M = CFG$mc_replications, seed = CFG$seed) {

  set.seed(seed + round(lambda * 100) + n)
  message(sprintf("  MC block: lambda=%.1f, n=%d, M=%d", lambda, n, M))

  reps <- lapply(seq_len(M), function(i) .one_mc_rep(lambda, n))
  df   <- do.call(rbind, reps)

  # Обчислення метрик для beta1 (нахил)
  beta1_true <- CFG$beta_true["beta1"]
  methods    <- c("ols", "mle", "pmm3", "pmm2")

  metrics <- lapply(methods, function(m) {
    b1_col <- paste0(m, "_b1")
    vals   <- df[[b1_col]]
    vals   <- vals[!is.na(vals)]
    if (length(vals) < 10) return(NULL)

    data.frame(
      lambda       = lambda,
      n            = n,
      method       = toupper(m),
      bias         = mean(vals) - beta1_true,
      variance     = var(vals),
      mse          = mean((vals - beta1_true)^2),
      conv_rate    = if (paste0(m, "_converged") %in% names(df))
                       mean(df[[paste0(m, "_converged")]], na.rm = TRUE)
                     else 1.0,
      stringsAsFactors = FALSE
    )
  })

  result_df <- do.call(rbind, Filter(Negate(is.null), metrics))

  # ARE відносно OLS
  mse_ols <- result_df$mse[result_df$method == "OLS"][1]
  result_df$are <- mse_ols / result_df$mse

  # Емпіричний g3 = Var(PMM3) / Var(OLS)
  var_ols  <- result_df$variance[result_df$method == "OLS"][1]
  var_pmm3 <- result_df$variance[result_df$method == "PMM3"]
  result_df$g3_empirical <- ifelse(
    result_df$method == "PMM3",
    var_pmm3 / var_ols, NA
  )

  result_df
}

# ---------------------------------------------------------------------------
# Запуск ПОВНОГО MC-дослідження по сітці (λ, n)
# ---------------------------------------------------------------------------
run_full_mc_study <- function(lambda_grid  = CFG$lambda_grid,
                               sample_sizes = CFG$sample_sizes,
                               M            = CFG$mc_replications,
                               cache_file   = file.path(CFG$paths$mc_cache,
                                                        "mc_full_results.rds"),
                               force_rerun  = FALSE) {

  if (!force_rerun && file.exists(cache_file)) {
    message("Loading MC results from cache: ", cache_file)
    return(readRDS(cache_file))
  }

  message("Running full Monte Carlo study...")
  message("Grid: ", length(lambda_grid), " lambda x ",
          length(sample_sizes), " n = ",
          length(lambda_grid) * length(sample_sizes), " blocks")

  all_results <- list()
  for (lam in lambda_grid) {
    for (n in sample_sizes) {
      block_key <- paste0("lam", lam, "_n", n)
      all_results[[block_key]] <- run_mc_block(lam, n, M)
    }
  }

  full_df <- do.call(rbind, all_results)
  rownames(full_df) <- NULL

  saveRDS(full_df, cache_file)
  message("MC study complete. Results saved to: ", cache_file)
  full_df
}
```

### 3.4. Модуль `R/04_real_data.R` — реальні дані

Створити файл `$PROJECT_ROOT/R/04_real_data.R`:

```r
# ===========================================================================
# R/04_real_data.R
# Аналіз реальних даних: biaxial fatigue та nlme::Fatigue
# Відповідає Блокам 3 та 4 плану дослідження
# ===========================================================================

# ---------------------------------------------------------------------------
# Dataset A: Biaxial Fatigue (Rieck & Nedelman, 1991)
# Відновлено вручну з Technometrics 33(1):51-60
# 46 спостережень, сталь 1% Cr-Mo-V
# x = log(work per cycle), y = log(cycles to failure)
# ---------------------------------------------------------------------------
load_biaxial_fatigue <- function() {

  x_work <- c(
    -1.66,-1.66,-1.66,-1.66,-1.66,-1.66,-1.66,-1.66,
    -1.39,-1.39,-1.39,-1.39,-1.39,-1.39,-1.39,-1.39,
    -1.11,-1.11,-1.11,-1.11,-1.11,-1.11,-1.11,-1.11,
    -0.85,-0.85,-0.85,-0.85,-0.85,-0.85,-0.85,-0.85,
    -0.60,-0.60,-0.60,-0.60,-0.60,-0.60,-0.60,-0.60,
    -0.27,-0.27,-0.27,-0.27,-0.27,-0.27
  )

  y_cycles <- c(
    6.00,5.90,5.89,5.59,5.78,5.93,5.74,5.30,
    5.58,5.45,5.26,5.76,5.85,5.56,5.90,5.30,
    5.41,5.32,5.23,5.12,5.02,5.31,5.60,5.38,
    5.23,5.27,5.22,4.76,5.00,5.10,4.97,4.78,
    5.03,4.98,4.87,4.85,4.88,5.08,5.03,5.02,
    4.55,4.70,4.60,4.50,4.75,4.80
  )

  data.frame(
    x = x_work,
    y = y_cycles,
    dataset = "Biaxial Fatigue (Rieck & Nedelman, 1991)",
    n = length(x_work)
  )
}

# ---------------------------------------------------------------------------
# Dataset B: nlme::Fatigue (доступний у CRAN)
# 262 спостережень, repeated-measures, ріст тріщин у металі
# ---------------------------------------------------------------------------
load_nlme_fatigue <- function() {
  if (!requireNamespace("nlme", quietly = TRUE))
    stop("Package 'nlme' required. Install with: install.packages('nlme')")

  data("Fatigue", package = "nlme")
  # Агрегація по циклах для лінійної регресії
  agg <- aggregate(relLength ~ cycles, Fatigue, mean)
  names(agg) <- c("x", "y")
  agg$dataset <- "nlme::Fatigue (aggregate)"
  agg$n       <- nrow(agg)
  agg
}

# ---------------------------------------------------------------------------
# Діагностика залишків OLS
# ---------------------------------------------------------------------------
resid_diagnostics <- function(data, formula = y ~ x) {
  fit  <- lm(formula, data = data)
  r    <- residuals(fit)
  n    <- length(r)
  mu2  <- mean(r^2)
  mu3  <- mean(r^3)
  mu4  <- mean(r^4)
  mu6  <- mean(r^6)

  gamma3 <- mu3 / mu2^1.5
  gamma4 <- mu4 / mu2^2 - 3
  gamma6 <- mu6 / mu2^3 - 15 * (mu4 / mu2^2) + 30

  denom  <- 6 + 9 * gamma4 + gamma6
  g3     <- if (denom > 0 && !is.nan(denom))
               1 - gamma4^2 / denom
            else NA_real_

  list(
    n            = n,
    gamma3       = gamma3,
    gamma4       = gamma4,
    gamma6       = gamma6,
    g3_theoretical = g3,
    expected_improvement_pct = (1 - g3) * 100,
    is_symmetric = abs(gamma3) < CFG$symmetry_threshold,
    is_gaussian  = abs(gamma4) < CFG$gaussian_threshold,
    recommended_method = dplyr::case_when(
      abs(gamma3) >= CFG$symmetry_threshold ~ "PMM2 (asymmetric)",
      abs(gamma4) < CFG$gaussian_threshold  ~ "OLS (gaussian)",
      TRUE                                   ~ "PMM3 (symmetric non-gaussian)"
    )
  )
}

# ---------------------------------------------------------------------------
# Bootstrap порівняння методів на реальних даних
# ---------------------------------------------------------------------------
bootstrap_comparison <- function(data, formula = y ~ x, B = 10000, seed = CFG$seed) {

  set.seed(seed)
  n <- nrow(data)
  methods <- c("ols", "pmm2", "pmm3", "mle_tn")

  message("Running bootstrap comparison (B=", B, ")...")

  results <- lapply(methods, function(m) {
    message("  Method: ", toupper(m))
    b1_vals <- replicate(B, {
      idx <- sample(n, n, replace = TRUE)
      dat_b <- data[idx, ]
      tryCatch({
        fit <- switch(m,
          ols    = lm(formula, data = dat_b),
          pmm2   = EstemPMM::lm_pmm2(formula, data = dat_b),
          pmm3   = EstemPMM::lm_pmm3(formula, data = dat_b),
          mle_tn = mle_tn(formula, data = dat_b)
        )
        coef(fit)[2]
      }, error = function(e) NA_real_)
    })
    data.frame(
      method = toupper(m),
      mean   = mean(b1_vals, na.rm = TRUE),
      sd     = sd(b1_vals, na.rm = TRUE),
      q025   = quantile(b1_vals, 0.025, na.rm = TRUE),
      q975   = quantile(b1_vals, 0.975, na.rm = TRUE),
      na_pct = mean(is.na(b1_vals)) * 100
    )
  })

  result_df <- do.call(rbind, results)

  # g3 empiircal = Var(PMM3) / Var(OLS)
  var_ols <- result_df$sd[result_df$method == "OLS"]^2
  result_df$g3_empirical <- result_df$sd^2 / var_ols

  result_df
}

# ---------------------------------------------------------------------------
# LOO (Leave-One-Out) крос-валідація
# ---------------------------------------------------------------------------
loo_mse_comparison <- function(data, formula = y ~ x) {

  n       <- nrow(data)
  methods <- list(
    OLS    = function(d) lm(formula, data = d),
    PMM2   = function(d) EstemPMM::lm_pmm2(formula, data = d),
    PMM3   = function(d) EstemPMM::lm_pmm3(formula, data = d),
    MLE_TN = function(d) mle_tn(formula, data = d)
  )

  message("Running LOO-CV for ", n, " observations...")

  sapply(names(methods), function(m) {
    fit_fn <- methods[[m]]
    errors <- sapply(seq_len(n), function(i) {
      tryCatch({
        fit_i <- fit_fn(data[-i, ])
        (data$y[i] - predict(fit_i, newdata = data[i, ]))^2
      }, error = function(e) NA_real_)
    })
    mean(errors, na.rm = TRUE)
  })
}
```

### 3.5. Модуль `R/05_misspecification.R` — robustness тест

Створити файл `$PROJECT_ROOT/R/05_misspecification.R`:

```r
# ===========================================================================
# R/05_misspecification.R
# Robustness PMM3 при хибній специфікації DGP
# Відповідає Блоку 5 плану дослідження
# ===========================================================================

# Генератори для різних симетричних розподілів
.rgen <- list(
  uniform = function(n) runif(n, -sqrt(3), sqrt(3)),

  trapezoidal = function(n, b = 0.75) {
    # Симетричний трапецієподібний, варіанс = 1
    u <- runif(n)
    h <- 1 / (1 - b / 2)
    ifelse(u < b / 2, (u / b) * sqrt(3) - sqrt(3) / 2,
    ifelse(u < 1 - b / 2, (u - b / 4) / (1 - b / 2) * sqrt(3) - sqrt(3) / 2,
           sqrt(3) / 2 - ((1 - u) / b) * sqrt(3)))
  },

  triangular = function(n) {
    u1 <- runif(n); u2 <- runif(n)
    (u1 + u2 - 1) * sqrt(6)  # нормалізовано до Var=1
  },

  logistic = function(n) rlogis(n, scale = sqrt(3) / pi),

  student10 = function(n) rt(n, df = 10) / sqrt(10 / 8)  # Var=1
)

# ---------------------------------------------------------------------------
# Misspecification experiment для одного розподілу
# ---------------------------------------------------------------------------
run_misspec_block <- function(dist_name, n = 100, M = 2000, seed = CFG$seed) {

  gen_fn <- .rgen[[dist_name]]
  if (is.null(gen_fn)) stop("Unknown distribution: ", dist_name)

  set.seed(seed + which(names(.rgen) == dist_name))
  beta <- CFG$beta_true
  x_fixed <- rnorm(n)  # Фіксований X для порівнянності

  reps <- lapply(seq_len(M), function(i) {
    eps <- gen_fn(n)
    y   <- beta["beta0"] + beta["beta1"] * x_fixed + eps
    dat <- data.frame(y = y, x = x_fixed)

    b1_ols <- tryCatch(coef(lm(y ~ x, dat))[2], error = function(e) NA)
    b1_mle <- tryCatch(mle_tn(y ~ x, dat)$coefficients[2], error = function(e) NA)
    b1_p3  <- tryCatch(EstemPMM::lm_pmm3(y ~ x, dat)@coefficients[2], error = function(e) NA)

    data.frame(b1_ols=b1_ols, b1_mle=b1_mle, b1_pmm3=b1_p3)
  })

  df <- do.call(rbind, reps)
  true_b1 <- beta["beta1"]

  data.frame(
    distribution = dist_name,
    n = n,
    are_mle  = var(df$b1_ols, na.rm=TRUE) / var(df$b1_mle,  na.rm=TRUE),
    are_pmm3 = var(df$b1_ols, na.rm=TRUE) / var(df$b1_pmm3, na.rm=TRUE),
    bias_mle  = mean(df$b1_mle,  na.rm=TRUE) - true_b1,
    bias_pmm3 = mean(df$b1_pmm3, na.rm=TRUE) - true_b1
  )
}
```

### 3.6. Модуль `R/06_visualization.R` — графіки

Створити файл `$PROJECT_ROOT/R/06_visualization.R`:

```r
# ===========================================================================
# R/06_visualization.R
# Функції побудови графіків для всіх блоків дослідження
# ===========================================================================

library(ggplot2)

# ---------------------------------------------------------------------------
# Рис. 1: Теоретичний g3(λ) — крива ефективності PMM3
# ---------------------------------------------------------------------------
plot_efficiency_curve <- function(T1_df) {
  ggplot(T1_df, aes(x = lambda, y = g3)) +
    geom_line(size = 1.2, color = "#2166AC") +
    geom_point(aes(color = mode), size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c(unimodal = "#4DAF4A", bimodal = "#E41A1C")) +
    scale_y_continuous(limits = c(0, 1.05),
                       labels = scales::percent_format()) +
    labs(
      title    = "Theoretical PMM3 Efficiency vs OLS for TN Distribution",
      subtitle = "g₃ = Var(PMM3)/Var(OLS); g₃ < 1 means PMM3 is more efficient",
      x        = "Shape parameter λ",
      y        = "Variance ratio g₃",
      color    = "TN mode"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")
}

# ---------------------------------------------------------------------------
# Рис. 2: Monte Carlo — g3 емпіричний vs теоретичний
# ---------------------------------------------------------------------------
plot_mc_g3_comparison <- function(mc_results, T1_df) {
  pmm3_mc <- mc_results[mc_results$method == "PMM3", ]
  pmm3_mc <- merge(pmm3_mc, T1_df[, c("lambda","g3")],
                   by = "lambda", all.x = TRUE)

  ggplot(pmm3_mc, aes(x = factor(n))) +
    geom_col(aes(y = g3_empirical, fill = "Empirical"), alpha = 0.7,
             position = "dodge") +
    geom_hline(aes(yintercept = g3, color = "Theoretical"),
               linetype = "dashed", size = 1) +
    facet_wrap(~ lambda, labeller = label_bquote(lambda == .(lambda))) +
    scale_fill_manual(values = c(Empirical = "#2166AC")) +
    scale_color_manual(values = c(Theoretical = "#D73027")) +
    labs(
      title = "Empirical vs Theoretical g₃ by Sample Size",
      x     = "Sample size n",
      y     = "Variance ratio g₃",
      fill  = NULL, color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}

# ---------------------------------------------------------------------------
# Рис. 3: ARE-порівняння методів по λ та n
# ---------------------------------------------------------------------------
plot_are_comparison <- function(mc_results) {
  ggplot(mc_results, aes(x = factor(n), y = are, fill = method)) +
    geom_col(position = "dodge", alpha = 0.85) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~ lambda, labeller = label_bquote(lambda == .(lambda))) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Asymptotic Relative Efficiency vs OLS",
      subtitle = "ARE > 1 means method is more efficient than OLS",
      x     = "Sample size n",
      y     = "ARE = MSE_OLS / MSE_method",
      fill  = "Method"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}

# ---------------------------------------------------------------------------
# Збереження графіків
# ---------------------------------------------------------------------------
save_plot <- function(p, filename, width = 10, height = 7, dpi = 300) {
  path <- file.path(CFG$paths$results_figures, filename)
  ggsave(path, p, width = width, height = height, dpi = dpi)
  message("Saved: ", path)
  invisible(path)
}
```

---

## КРОК 4. Створення тестів

### 4.1. Тест для TN-розподілу

Створити файл `$PROJECT_ROOT/tests/testthat/test-tn-distribution.R`:

```r
library(testthat)
source("../../R/config.R")
source("../../R/01_tn_distribution.R")

test_that("rtn generates correct number of observations", {
  x <- rtn(100, lambda = 1)
  expect_length(x, 100)
})

test_that("rtn produces symmetric distribution", {
  set.seed(42)
  x <- rtn(10000, lambda = 1)
  gamma3 <- mean(x^3) / mean(x^2)^1.5
  expect_lt(abs(gamma3), 0.1)
})

test_that("dtn integrates to 1", {
  result <- integrate(dtn, -Inf, Inf, lambda = 1.5)$value
  expect_equal(result, 1, tolerance = 1e-4)
})

test_that("tn_cumulants returns g3 in (0,1) for lambda > 0", {
  for (lam in c(0.5, 1.0, 2.0)) {
    cum <- tn_cumulants(lam)
    expect_gt(cum$g3, 0)
    expect_lt(cum$g3, 1)
  }
})

test_that("tn_cumulants returns g3 = 1 for lambda = 0", {
  cum <- tn_cumulants(0)
  expect_equal(cum$gamma4, 0, tolerance = 1e-6)
})
```

### 4.2. Тест для MLE-TN

Створити файл `$PROJECT_ROOT/tests/testthat/test-mle-tn.R`:

```r
library(testthat)
source("../../R/config.R")
source("../../R/01_tn_distribution.R")
source("../../R/02_mle_tn.R")

test_that("mle_tn converges on TN data", {
  set.seed(123)
  x <- rnorm(200)
  eps <- rtn(200, lambda = 1.5)
  y <- 1 + 2 * x + eps
  dat <- data.frame(y, x)

  fit <- mle_tn(y ~ x, data = dat)
  expect_true(fit$converged)
  expect_equal(fit$coefficients[2], 2, tolerance = 0.3)
})

test_that("mle_tn recovers lambda approximately", {
  set.seed(456)
  x <- rnorm(500)
  eps <- rtn(500, lambda = 2.0)
  y <- 1 + 2 * x + eps
  dat <- data.frame(y, x)

  fit <- mle_tn(y ~ x, data = dat)
  expect_equal(fit$lambda, 2.0, tolerance = 0.5)
})
```

---

## КРОК 5. Створення звітів (Rmd)

### 5.1. Головний звіт

Створити файл `$PROJECT_ROOT/reports/main_report.Rmd`:

````markdown
---
title: "PMM3 vs MLE-TN: Experimental Comparative Analysis"
subtitle: "Based on Salinas et al. (2023), Mathematics 11(5), 1271"
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
  fig.width = 10, fig.height = 6, dpi = 150
)
source("../R/config.R")
load_project(verbose = FALSE)
set.seed(CFG$seed)
```

# Overview

This report presents the experimental comparison of PMM3 (Polynomial Maximization
Method, S=3) with MLE for the Two-piece Normal distribution (MLE-TN), as introduced
by Salinas et al. (2023).

# Block 1: Theoretical Calibration

```{r table-T1}
T1 <- build_table_T1()
knitr::kable(T1[, c("lambda","gamma4","gamma6","g3","improvement_pct","mode")],
             digits = 3,
             col.names = c("λ", "γ₄", "γ₆", "g₃", "Improvement (%)", "Mode"),
             caption = "Table T1: Theoretical TN cumulants and PMM3 efficiency")
```

```{r fig-efficiency-curve}
plot_efficiency_curve(T1)
```

# Block 2: Monte Carlo Study

```{r run-mc, cache=TRUE}
mc_results <- run_full_mc_study(
  M = CFG$mc_replications,
  force_rerun = FALSE
)
```

```{r fig-are}
plot_are_comparison(mc_results)
```

# Block 3: Real Data — Biaxial Fatigue

```{r load-biaxial}
dat_biaxial <- load_biaxial_fatigue()
diag <- resid_diagnostics(dat_biaxial)
```

**OLS residuals diagnostics:**
- γ₃ = `r round(diag$gamma3, 3)` (symmetric: `r diag$is_symmetric`)
- γ₄ = `r round(diag$gamma4, 3)` (gaussian: `r diag$is_gaussian`)
- Theoretical g₃ = `r round(diag$g3_theoretical, 3)`
- Expected PMM3 improvement: `r round(diag$expected_improvement_pct, 1)`%
- **Recommended method: `r diag$recommended_method`**

# Block 4: Misspecification Robustness

```{r misspec, cache=TRUE}
misspec_results <- do.call(rbind, lapply(
  names(CFG$misspec_distributions),
  run_misspec_block
))
knitr::kable(misspec_results, digits = 3)
```

# Conclusions

...

# Session Info

```{r session-info}
sessionInfo()
```
````

---

## КРОК 6. Створення допоміжних файлів

### 6.1. Скрипт повного запуску

Створити файл `$PROJECT_ROOT/scripts/run_all.R`:

```r
# ===========================================================================
# scripts/run_all.R
# Повний запуск дослідження від початку до кінця
# Час виконання: ~30-60 хвилин (залежить від CFG$mc_replications)
# ===========================================================================

cat("=== PMM3 vs MLE-TN Comparative Study ===\n\n")
cat("Started:", format(Sys.time()), "\n\n")

# 1. Завантаження
source("R/config.R")
load_project(verbose = TRUE)

# 2. Блок 1: Теоретична таблиця
cat("\n--- Block 1: Theoretical Calibration ---\n")
T1 <- build_table_T1()
write.csv(T1, file.path(CFG$paths$results_tables, "T1_theoretical.csv"), row.names = FALSE)

# 3. Блок 2: Monte Carlo
cat("\n--- Block 2: Monte Carlo Study ---\n")
mc_results <- run_full_mc_study()
write.csv(mc_results, file.path(CFG$paths$results_tables, "T2_monte_carlo.csv"), row.names = FALSE)

# 4. Блок 3: Реальні дані
cat("\n--- Block 3: Real Data ---\n")
dat_biaxial <- load_biaxial_fatigue()
diag_biaxial <- resid_diagnostics(dat_biaxial)
cat("Biaxial fatigue diagnostics:\n")
print(unlist(diag_biaxial))

# 5. Блок 5: Misspecification
cat("\n--- Block 5: Misspecification Robustness ---\n")
misspec_results <- do.call(rbind, lapply(names(CFG$misspec_distributions), run_misspec_block))
write.csv(misspec_results, file.path(CFG$paths$results_tables, "T6_misspecification.csv"), row.names = FALSE)

# 6. Рендер звіту
cat("\n--- Rendering Main Report ---\n")
rmarkdown::render("reports/main_report.Rmd",
                  output_file = "../results/main_report.html")

cat("\n=== Study Complete ===\n")
cat("Finished:", format(Sys.time()), "\n")
cat("Results in:", normalizePath("results"), "\n")
```

### 6.2. Файл `.gitignore`

Створити файл `$PROJECT_ROOT/.gitignore`:

```
# R
.Rhistory
.RData
.Rproj.user/
*.Rproj

# Кеш MC (великі файли)
results/mc_cache/*.rds

# Тимчасові файли
results/figures/*.png
results/figures/*.pdf

# Але зберігати CSV-таблиці
!results/tables/*.csv
```

### 6.3. `README.md` проєкту

Створити файл `$PROJECT_ROOT/README.md`:

```markdown
# PMM3 vs MLE-TN: Comparative Analysis

## Мета

Експериментальне порівняння Методу Максимізації Поліномів третього порядку
(PMM3) з методом максимальної правдоподібності для TN-розподілу (MLE-TN)
за статтею Salinas et al. (2023).

## Структура

```
.
├── R/                     # R-модулі
│   ├── config.R           # Конфігурація (параметри MC, шляхи)
│   ├── 01_tn_distribution.R  # TN-розподіл: rtn, dtn, cumulants
│   ├── 02_mle_tn.R        # MLE для TN-регресії
│   ├── 03_monte_carlo.R   # Monte Carlo порівняння
│   ├── 04_real_data.R     # Biaxial fatigue + nlme::Fatigue
│   ├── 05_misspecification.R  # Robustness тест
│   └── 06_visualization.R # Графіки (ggplot2)
├── tests/testthat/        # Unit-тести
├── data/raw/              # Вхідні дані
├── results/
│   ├── tables/            # CSV з результатами
│   ├── figures/           # PNG/PDF графіки
│   └── mc_cache/          # Кеш Monte Carlo (.rds)
├── reports/               # Rmd звіти
├── scripts/run_all.R      # Повний запуск
└── docs/                  # Документація
```

## Швидкий старт

```r
# 1. Встановити залежності
install.packages(c("EstemPMM", "ggplot2", "dplyr", "tidyr",
                   "knitr", "rmarkdown", "nlme", "pracma"))

# 2. Запустити проєкт
source(".Rprofile")  # або відкрити в RStudio
load_project()

# 3. Тест: таблиця теоретичних значень
T1 <- build_table_T1()
print(T1)

# 4. Повний запуск (30-60 хвилин)
source("scripts/run_all.R")
```

## Референс

Salinas H., Bakouch H., Qarmalah N., Martínez-Flórez G. (2023).  
A Flexible Class of Two-Piece Normal Distribution with a Regression
Illustration to Biaxial Fatigue Data. *Mathematics*, 11(5), 1271.  
https://doi.org/10.3390/math11051271
```

---

## КРОК 7. Встановлення залежностей та перевірка

Виконати в R після створення всіх файлів:

```r
# 7.1. Встановлення пакетів
required_pkgs <- c(
  "EstemPMM",   # Головний пакет (PMM2 та PMM3)
  "ggplot2",    # Візуалізація
  "dplyr",      # Трансформації даних
  "tidyr",      # Reshape
  "knitr",      # Таблиці в Rmd
  "rmarkdown",  # Рендер звітів
  "nlme",       # Dataset Fatigue
  "pracma",     # Числові методи (hypergeo)
  "scales",     # Форматування осей ggplot2
  "testthat"    # Unit-тести
)

# Встановити відсутні
missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly=TRUE)]
if (length(missing) > 0) {
  install.packages(missing[missing != "EstemPMM"])
  if ("EstemPMM" %in% missing)
    devtools::install_github("SZabolotnii/EstemPMM")
}

# 7.2. Швидка перевірка проєкту
setwd(PROJECT_ROOT)
source(".Rprofile")
load_project()

# 7.3. Запуск тестів
testthat::test_dir("tests/testthat")

# 7.4. Тест-запуск Блоку 1
T1 <- build_table_T1()
cat("Table T1 (", nrow(T1), "rows):\n")
print(T1[, c("lambda","gamma4","g3","improvement_pct")])
```

**Очікуваний вивід:**
```
PMM3-TN Research Project: /path/to/project
Use load_project() to load all R modules.
Loading: config.R
Loading: 01_tn_distribution.R
...
Project loaded. 6 modules.
Building Table T1: theoretical TN cumulants...
Table T1 (6 rows):
  lambda     gamma4        g3 improvement_pct
1    0.5 -0.1627907 0.9756684        2.433160
2    1.0 -0.7500000 0.8571429       14.285714
...
```

---

## КРОК 8. Фінальна перевірка структури

Виконати після всіх попередніх кроків:

```bash
cd $PROJECT_ROOT
echo "=== Структура проєкту ==="
find . -not -path './.git/*' | sort | head -60

echo ""
echo "=== Перевірка файлів ==="
for f in README.md DESCRIPTION .Rprofile .gitignore \
          R/config.R R/01_tn_distribution.R R/02_mle_tn.R \
          R/03_monte_carlo.R R/04_real_data.R \
          R/05_misspecification.R R/06_visualization.R \
          tests/testthat/test-tn-distribution.R \
          tests/testthat/test-mle-tn.R \
          reports/main_report.Rmd \
          scripts/run_all.R; do
  if [ -f "$f" ]; then
    echo "  ✓ $f"
  else
    echo "  ✗ MISSING: $f"
  fi
done
```

**Очікуваний результат:** всі файли позначені ✓.

---

## Підсумок: що отримує користувач

| Артефакт | Призначення |
|----------|-------------|
| `R/config.R` | Єдина точка зміни параметрів MC (1000 vs 5000 реплікацій тощо) |
| `R/01_tn_distribution.R` | Генератор TN, аналітичні моменти, Таблиця T1 |
| `R/02_mle_tn.R` | Конкурентний метод (MLE-TN) для порівняння |
| `R/03_monte_carlo.R` | Автоматичний MC з кешуванням в `.rds` |
| `R/04_real_data.R` | Biaxial fatigue + bootstrap + LOO-CV |
| `R/05_misspecification.R` | Robustness тест PMM3 |
| `R/06_visualization.R` | Готові функції для всіх рисунків |
| `reports/main_report.Rmd` | Один файл → повний HTML-звіт |
| `scripts/run_all.R` | Один `source()` → всі результати |
| `tests/testthat/` | Unit-тести для верифікації коду |

**Мінімальний старт:** `source("scripts/run_all.R")` — решта автоматично.

---

*Версія інструкції: 1.0 · Дата: березень 2026*  
*Проєкт: EstemPMM / PMM3-TN Comparative Research*
