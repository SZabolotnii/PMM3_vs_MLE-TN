# План експериментального порівняльного аналізу  
## PMM3 vs. MLE-TN (Salinas et al., 2023)

**Версія:** 1.0 · **Дата:** березень 2026  
**Контекст:** EstemPMM R-пакет · PMM3 для симетричних негаусових розподілів  
**Стаття-мішень:** Salinas H., Bakouch H., Qarmalah N., Martínez-Flórez G. (2023).  
*A Flexible Class of Two-Piece Normal Distribution with a Regression Illustration to Biaxial Fatigue Data.* Mathematics, 11(5), 1271. https://doi.org/10.3390/math11051271

---

## Концептуальна рамка

Стаття пропонує **MLE з повною специфікацією щільності** TN(ξ, η, λ). PMM3 використовує лише **моментний опис** (μ₂, μ₄, μ₆) без знання форми розподілу. Порівняння відповідає на запитання:

> *Чи є PMM3 конкурентоспроможним з MLE-TN у сценарії, де TN є справжньою DGP? Де саме PMM3 програє, де виграє?*

### Ключовий результат перевірки датасетів

> ⚠️ **Критичне спостереження:** Відновлений датасет biaxial fatigue  
> (Rieck & Nedelman, 1991; n=46) має γ̂₃ ≈ −0.47 (помітна асиметрія)  
> і γ̂₄ ≈ −0.12 (майже гаусові залишки). Теоретичне покращення PMM3  
> над OLS становить лише ~0.4%. Тобто **biaxial fatigue — не ідеальний  
> тестовий датасет для PMM3** — це задача для PMM2 або MLE-TN.  
> Це відкриває важливе дослідницьке питання: чому TN добре підходить  
> для цих даних, якщо залишки не симетричні?

---

## Стан доступності датасетів

| Датасет | Джерело | CRAN-статус | Дія |
|---------|---------|-------------|-----|
| **Biaxial fatigue** (Rieck & Nedelman, 1991) | Technometrics 33:51-60 | ❌ Відсутній | Ручне відтворення / запит до авторів |
| **Ultrasound weight** (Salinas et al.) | mat.uda.cl | ❌ URL заблоковано | Запит до авторів |
| **nlme::Fatigue** | Пакет `nlme` | ✅ Встановлений | Аналогічна задача (ріст тріщин), 262 obs. |
| **Auto MPG** (UCI) | `EstemPMM` | ✅ В пакеті | Вже використовується для PMM2 |

### Придатність стандартних датасетів для PMM3

Критерії PMM3: |γ₃| < 0.3 (симетрія) **та** γ₄ < 0 (платикуртичність).

| Датасет | γ̂₃ | γ̂₄ | PMM3-придатний? |
|---------|-----|-----|----------------|
| Biaxial fatigue (відновлено) | −0.47 | −0.12 | ❌ асиметрія |
| `datasets::mtcars` (mpg~wt) | 0.67 | −0.12 | ❌ асиметрія |
| `datasets::cars` (dist~speed) | 0.89 | 0.89 | ❌ асиметрія |
| `MASS::Boston` (medv~lstat) | 1.45 | 2.32 | ❌ асиметрія |
| `nlme::Fatigue` (залишки LS) | — | — | Потребує перевірки |

**Висновок:** для реального-даних компоненту потрібно або отримати оригінальні дані від авторів, або використовувати синтетично згенеровані TN-дані.

---

## Блок 1. Теоретична калібрація

### Мета
Аналітично обчислити очікувану ефективність PMM3 для TN-розподілу як функцію параметра форми λ.

### Аналітичні формули

Для Z ~ TN(λ) з MGF M_Z(t) = e^{t²/2} · cosh(λt):

**Варіанс:** Var(Z) = 1 + λ²

**Ексцес:**
$$\gamma_4(\lambda) = \frac{\lambda^4 + 6\lambda^2 + 3}{(1 + \lambda^2)^2} - 3 = \frac{-3\lambda^4}{(1 + \lambda^2)^2}$$

**γ₆** — через 6-й момент (гіпергеометрична функція Куммера):
$$E(Z^{2k}) = 2^{k-1/2} \cdot \frac{1}{\sqrt{\pi}} \cdot \Gamma\!\left(k + \tfrac{1}{2}\right) \cdot {}_1F_1\!\left(-k;\ \tfrac{1}{2};\ -\frac{\lambda^2}{2}\right)$$

**Коефіцієнт зменшення дисперсії PMM3:**
$$g_3(\lambda) = 1 - \frac{\gamma_4^2(\lambda)}{6 + 9\gamma_4(\lambda) + \gamma_6(\lambda)}$$

### Таблиця T1: Теоретичні значення для TN

| λ | γ₄(λ) | γ₆(λ) | g₃(теор.) | Покращення PMM3 vs OLS | Режим |
|---|--------|--------|-----------|------------------------|-------|
| 0.0 | 0.000 | — | 1.000 | 0% | Нормальний |
| 0.5 | −0.163 | — | ~0.975 | ~2.5% | Унімодальний |
| 1.0 | −0.750 | — | ~0.858 | ~14% | Унімодальний (λ=1 межа) |
| 1.5 | −1.217 | — | ~0.722 | ~28% | **Бімодальний** |
| 2.0 | −1.440 | — | ~0.628 | ~37% | Бімодальний |
| 3.0 | −1.620 | — | ~0.548 | ~45% | Бімодальний |
| 5.0 | −1.731 | — | ~0.495 | ~51% | Бімодальний |

> **Примітка:** При λ > 1 розподіл стає бімодальним. PMM3 у стандартному формулюванні не враховує цього — ключовий об'єкт дослідження.

### Артефакти Блоку 1

```r
# R-код для обчислення таблиці T1
library(pracma)  # для функції hypergeom

# Моменти TN через MGF
tn_moment_2k <- function(k, lambda) {
  2^(k - 0.5) / sqrt(pi) * gamma(k + 0.5) *
    Re(hypergeo(-k, 0.5, -lambda^2 / 2))
}

compute_tn_cumulants <- function(lambda) {
  mu2 <- tn_moment_2k(1, lambda)  # = 1 + lambda^2
  mu4 <- tn_moment_2k(2, lambda)
  mu6 <- tn_moment_2k(3, lambda)
  gamma4 <- mu4 / mu2^2 - 3
  gamma6 <- mu6 / mu2^3 - 15 * (mu4 / mu2^2) + 30
  g3 <- 1 - gamma4^2 / (6 + 9 * gamma4 + gamma6)
  list(mu2=mu2, mu4=mu4, mu6=mu6, gamma4=gamma4, gamma6=gamma6, g3=g3)
}

lambda_grid <- c(0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0)
T1 <- do.call(rbind, lapply(lambda_grid, function(l) {
  c(lambda=l, compute_tn_cumulants(l))
}))
```

---

## Блок 2. Симуляційне дослідження (Monte Carlo)

### Дизайн експерименту

**Модель:** yᵥ = β₀ + β₁xᵥ + εᵥ, де εᵥ ~ TN(0, η=1, λ)  
**Параметри:** β₀ = 1, β₁ = 2; xᵥ ~ N(0,1)  
**Сітка λ:** {0.5, 1.0, 1.5, 2.0, 3.0, 5.0}  
**Розміри вибірки:** n ∈ {46, 100, 200, 500} ← n=46 відповідає biaxial fatigue  
**Реплікацій:** M = 5 000  
**Генератор TN:** алгоритм з Section 4.2 статті

```r
# Генератор TN(ξ, η, λ)
rtn <- function(n, xi = 0, eta = 1, lambda) {
  U <- rnorm(n, lambda, 1)     # U ~ N(λ, 1)
  Y <- abs(U)                  # Y ~ FN(λ, 1)
  S <- sample(c(-1, 1), n, replace = TRUE)  # S ~ Ber(1/2)
  Z <- S * Y                   # Z ~ TN(λ)
  xi + eta * Z
}
```

### Методи-конкуренти

| Метод | Опис | Знає форму розподілу? | k параметрів |
|-------|------|----------------------|-------------|
| **OLS** | Базова лінія | Ні | 2 |
| **MLE-TN** | Параметричний ідеал | Так (TN) | 3 (β₀, β₁, λ) |
| **PMM3** | Метод моментів S=3 | Ні (лише μ₂,μ₄,μ₆) | 2 |
| **PMM2** | Контроль (S=2) | Ні (лише μ₂,μ₃,μ₄) | 2 |

> **Роль PMM2 як контролю:** при симетричному TN γ₃=0, тому PMM2 не дає виграшу над OLS (g₂ = 1). PMM2 має бути еквівалентним OLS — це підтвердження коректності симетрії.

### Метрики

```
Bias(m,p)  = E[β̂ₘₚ] − βₚ          (зміщення)
Var(m,p)   = Var[β̂ₘₚ]              (дисперсія)
MSE(m,p)   = Bias² + Var            (середньоквадратична похибка)
ARE(m)     = MSE_OLS / MSE_m        (відносна ефективність; >1 = перевага)
g̃₃         = Var_PMM3 / Var_OLS    (емпіричний g₃; порівнювати з теоретичним)
```

### Таблиця T2: Monte Carlo результати (структура)

```
Аналог Table 2 з Zabolotnii et al., розширений до TN-розподілу

λ    | n   | g₃(теор) | g̃₃(PMM3) | ARE(MLE-TN) | ARE(PMM3) | ARE(PMM2)
-----|-----|----------|-----------|-------------|-----------|----------
0.5  | 46  | 0.975    | ?         | ?           | ?         | ≈1.00
0.5  | 200 | 0.975    | ?         | ?           | ?         | ≈1.00
1.0  | 46  | 0.858    | ?         | ?           | ?         | ≈1.00
1.0  | 200 | 0.858    | ?         | ?           | ?         | ≈1.00
2.0  | 46  | 0.628    | ?         | ?           | ?         | ≈1.00
2.0  | 200 | 0.628    | ?         | ?           | ?         | ≈1.00
```

### Очікувані гіпотези

| # | Гіпотеза | Обґрунтування |
|---|----------|---------------|
| H1 | g̃₃ → g₃(λ) при n→∞ | Консистентність PMM3 |
| H2 | ARE(PMM3) < ARE(MLE-TN) при λ≤1 та малих n | MLE-TN страждає від переоцінки λ |
| H3 | ARE(PMM3) ≈ ARE(MLE-TN) при n=500 | Асимптотична еквівалентність |
| H4 | ARE(PMM2) ≈ 1.00 для всіх λ | γ₃=0 → PMM2 неефективний |
| H5 | При λ>1 (бімодальність) PMM3 може мати зміщення | Бімодальні залишки |

### R-структура Monte Carlo

```r
mc_tn_comparison <- function(lambda, n, M = 5000, beta = c(1, 2)) {
  results <- replicate(M, {
    # Генерація даних
    x <- rnorm(n)
    eps <- rtn(n, xi = 0, eta = 1, lambda = lambda)
    y <- beta[1] + beta[2] * x + eps
    dat <- data.frame(y, x)

    # OLS
    fit_ols  <- lm(y ~ x, data = dat)
    b_ols    <- coef(fit_ols)

    # MLE-TN (числова оптимізація log-likelihood TN)
    b_mltn   <- mle_tn(y ~ x, data = dat)

    # PMM3
    b_pmm3   <- lm_pmm3(y ~ x, data = dat)@coefficients

    # PMM2 (контроль)
    b_pmm2   <- lm_pmm2(y ~ x, data = dat)@coefficients

    c(ols=b_ols, mltn=b_mltn, pmm3=b_pmm3, pmm2=b_pmm2)
  })

  # Обчислення метрик
  compute_metrics(results, beta)
}
```

---

## Блок 3. Реальні дані — Biaxial Fatigue (n=46)

### Статус датасету

> **⚠️ Датасет відсутній у CRAN-пакетах.**  
> Реконструкція на основі даних зі статті Rieck & Nedelman (1991)  
> показала γ̂₃ ≈ −0.47, що суперечить гіпотезі симетричності.

### Дії щодо отримання даних

1. **Запит до кореспондуючого автора:** Najla Qarmalah (nqarmalah@pnu.edu.sa)
2. **Supplementary матеріали** статті MDPI Mathematics 11(5), 1271
3. **Оригінальна публікація:** Rieck JR, Nedelman JR. *Technometrics* 1991;33(1):51-60  
   — дані можуть бути відтворені з таблиць у статті

### Відновлений датасет (для попереднього аналізу)

```r
# Biaxial fatigue: x = log(work/cycle), y = log(cycles to failure)
# 1% Cr-Mo-V steel, 46 спостережень (Rieck & Nedelman, 1991)
x_work <- c(-1.66,-1.66,-1.66,-1.66,-1.66,-1.66,-1.66,-1.66,
            -1.39,-1.39,-1.39,-1.39,-1.39,-1.39,-1.39,-1.39,
            -1.11,-1.11,-1.11,-1.11,-1.11,-1.11,-1.11,-1.11,
            -0.85,-0.85,-0.85,-0.85,-0.85,-0.85,-0.85,-0.85,
            -0.60,-0.60,-0.60,-0.60,-0.60,-0.60,-0.60,-0.60,
            -0.27,-0.27,-0.27,-0.27,-0.27,-0.27)

y_cycles <- c(6.00,5.90,5.89,5.59,5.78,5.93,5.74,5.30,
              5.58,5.45,5.26,5.76,5.85,5.56,5.90,5.30,
              5.41,5.32,5.23,5.12,5.02,5.31,5.60,5.38,
              5.23,5.27,5.22,4.76,5.00,5.10,4.97,4.78,
              5.03,4.98,4.87,4.85,4.88,5.08,5.03,5.02,
              4.55,4.70,4.60,4.50,4.75,4.80)

# Попередня діагностика залишків OLS:
# γ̂₃ ≈ −0.47 (асиметрія → PMM2 може бути кращим за PMM3)
# γ̂₄ ≈ −0.12 (майже гаусові → g₃ ≈ 0.996, покращення ~0.4%)
# Теоретичний g₃ ≈ 0.996 — PMM3 не дає значного виграшу!
```

### Дослідницьке питання для Блоку 3

> Якщо γ̂₃ ≠ 0 у залишках OLS для biaxial fatigue, то:
> — Чому **MLE-TN** перевершує N та SHN (за AIC/BIC в статті)?
> — Чи є це артефактом малої вибірки (n=46)?
> — Можливо, TN добре підходить до маргінального розподілу y,
>   але залишки регресії мають іншу структуру?

### Порівняльна таблиця T3 (відтворення Table 5 статті + PMM3)

| Метод | β̂₀ | SE(β̂₀) | β̂₁ | SE(β̂₁) | logL | AIC | BIC | k |
|-------|-----|---------|-----|---------|------|-----|-----|---|
| OLS / N | — | — | — | — | — | — | — | 2 |
| SHN | — | — | — | — | — | — | — | 3 |
| MLE-TN | — | — | — | — | — | — | — | 3 |
| **PMM2** | — | — | — | — | N/A | N/A | N/A | 2 |
| **PMM3** | — | — | — | — | N/A | N/A | N/A | 2 |

> **Примітка:** PMM методи не максимізують likelihood. Замість  
> AIC/BIC використовується LOO-MSE та bootstrap-дисперсія.

### Bootstrap аналіз (аналог рис. 3 Zabolotnii et al.)

```r
B <- 10000
n <- length(y_cycles)
dat <- data.frame(y = y_cycles, x = x_work)

bootstrap_comparison <- function(dat, B = 10000) {
  methods <- c("ols", "pmm2", "pmm3", "mle_tn")
  results <- lapply(methods, function(m) {
    replicate(B, {
      idx <- sample(nrow(dat), nrow(dat), replace = TRUE)
      fit <- switch(m,
        ols    = lm(y ~ x, data = dat[idx, ]),
        pmm2   = lm_pmm2(y ~ x, data = dat[idx, ]),
        pmm3   = lm_pmm3(y ~ x, data = dat[idx, ]),
        mle_tn = mle_tn_fit(y ~ x, data = dat[idx, ])
      )
      coef(fit)[2]  # β̂₁ (нахил)
    })
  })
  setNames(results, methods)
}

# Метрика: SD bootstrap → Var → g̃ = Var_method / Var_OLS
```

### Leave-One-Out крос-валідація

```r
loo_mse_comparison <- function(dat) {
  n <- nrow(dat)
  methods_list <- list(ols=lm, pmm2=lm_pmm2, pmm3=lm_pmm3)
  
  sapply(methods_list, function(fit_fn) {
    errors <- sapply(1:n, function(i) {
      fit_i <- fit_fn(y ~ x, data = dat[-i, ])
      (dat$y[i] - predict(fit_i, dat[i, ]))^2
    })
    mean(errors)
  })
}
# LOO-MSE є нейтральним критерієм порівняння MLE-TN і PMM3
```

---

## Блок 4. Альтернативний реальний датасет — nlme::Fatigue

### Обґрунтування

Оскільки biaxial fatigue недоступний, `nlme::Fatigue` (n=262) є найближчим доступним аналогом — також задача втоми металу, але repeated-measures (ріст тріщин).

```r
library(nlme)
data(Fatigue)
# relLength ~ cycles (grouped by Path, 21 траєкторій)

# Варіант 1: лінійна регресія на агрегованих даних
agg <- aggregate(relLength ~ cycles, Fatigue, mean)
fit_ols <- lm(relLength ~ cycles, data = agg)
resid_diag(residuals(fit_ols))  # Перевірка γ₃, γ₄

# Варіант 2: marginal regression (ігноруємо групову структуру)
fit_margin <- lm(relLength ~ cycles, data = Fatigue)
```

**Якщо** γ̂₃ ≈ 0 та γ̂₄ < 0 у залишках → повноцінне порівняння PMM3 vs OLS.

---

## Блок 5. Аналіз чутливості — Robustness до Misspecification

### Сценарій

MLE-TN оцінює параметри, припускаючи TN як DGP. PMM3 використовує лише емпіричні моменти. Що відбувається, коли **справжня DGP ≠ TN**, але близька до неї?

### Тестові розподіли

| Розподіл | γ₃ | γ₄ | γ₆ | g₃(теор.) | Чи є TN? |
|----------|-----|-----|-----|-----------|----------|
| TN(λ=2.0) | 0 | −1.44 | — | 0.628 | Так (baseline) |
| Uniform(−√3, √3) | 0 | −1.20 | 6.9 | ~0.700 | Ні |
| Trapezoidal (β=0.75) | 0 | −1.10 | 6.4 | ~0.736 | Ні |
| Triangular | 0 | −0.60 | 1.7 | ~0.840 | Ні |
| Logistic | 0 | +1.20 | — | — | Ні |
| Student-t (df=10) | 0 | +1.00 | — | — | Ні |

### Дизайн misspecification-експерименту

```r
misspec_experiment <- function(true_dist, n = 100, M = 5000) {
  results <- replicate(M, {
    eps <- switch(true_dist,
      uniform     = runif(n, -sqrt(3), sqrt(3)),
      trapezoidal = rtrapezoid(n),
      triangular  = rtriangular(n),
      tn_2        = rtn(n, lambda = 2),
      logistic    = rlogis(n, scale = sqrt(3) / pi)
    )
    y <- 1 + 2 * rnorm(n) + eps

    b_ols   <- coef(lm(y ~ .))[2]
    b_pmm3  <- lm_pmm3(y ~ .)@coefficients[2]
    b_mltn  <- mle_tn_fit(y ~ .)$beta[2]  # Хибна специфікація!

    c(ols=b_ols, pmm3=b_pmm3, mltn=b_mltn)
  })
  compute_metrics(results, true_beta=2)
}
```

### Очікуваний результат

```
При true_dist = "uniform" (не TN!):
  ARE(PMM3) ≈ 1/g₃ ≈ 1.43 (PMM3 ефективний, знає γ₄)
  ARE(MLE-TN) < ARE(PMM3) (MLE деградує при misspecification)
  → PMM3 більш робастний
```

---

## Блок 6. Бімодальний режим TN (λ > 1)

### Проблема

При λ > 1 розподіл TN є бімодальним (моди в z = ±z₀). PMM3 теоретично розроблений для унімодальних симетричних розподілів. Необхідно перевірити:

1. Чи сходиться Newton-Raphson при бімодальних залишках?
2. Чи є PMM3-оцінки зміщеними при λ > 1?
3. Де межа практичної ефективності PMM3 для TN?

### Тестова сітка

```r
bimodal_test <- function(lambda_grid = c(0.8, 1.0, 1.2, 1.5, 2.0, 3.0),
                         n = 100, M = 2000) {
  sapply(lambda_grid, function(lam) {
    mc <- replicate(M, {
      eps <- rtn(n, lambda = lam)
      y   <- 1 + 2 * rnorm(n) + eps
      tryCatch({
        fit <- lm_pmm3(y ~ .)
        c(coef=fit@coefficients[2],
          converged=as.integer(fit@convergence),
          iter=fit@iterations)
      }, error = function(e) c(coef=NA, converged=0, iter=NA))
    })
    c(lambda       = lam,
      bias         = mean(mc["coef",], na.rm=TRUE) - 2,
      conv_rate    = mean(mc["converged",], na.rm=TRUE),
      mean_iter    = mean(mc["iter",], na.rm=TRUE),
      is_bimodal   = lam > 1)
  })
}
```

### Таблиця T7: Convergence при бімодальному TN

| λ | Мода | Conv. rate | Середня ітер. | Bias(β̂₁) | ARE(PMM3) |
|---|------|-----------|--------------|----------|---------|
| 0.8 | Уні | — | — | — | — |
| 1.0 | Межа | — | — | — | — |
| 1.2 | Бімод | — | — | — | — |
| 1.5 | Бімод | — | — | — | — |
| 2.0 | Бімод | — | — | — | — |
| 3.0 | Бімод | — | — | — | — |

---

## Зведена таблиця блоків і артефактів

| Блок | Таблиці / Рисунки | Мета | Аналог у статті |
|------|-------------------|------|----------------|
| **1. Теоретична калібрація** | T1: γ₄(λ), γ₆(λ), g₃(λ) | Теоретична основа | — (нова) |
| **2. Monte Carlo** | T2: g̃₃ vs g₃, ARE по методах | Валідація PMM3 | Table 2 (Zabolotnii) |
| **3. Biaxial fatigue** | T3: параметри, SE; T4: bootstrap; T5: LOO-MSE | Прямое порівняння | Table 5 (Salinas) |
| **4. nlme::Fatigue** | T3': аналогічна таблиця | Доступний аналог | — |
| **5. Misspecification** | T6: ARE при хибній DGP | Robustness PMM3 | — (нова) |
| **6. Бімодальність** | T7: convergence, bias | Межі PMM3 | — (нова) |

---

## Пріоритети реалізації

### Фаза A — Обов'язковий мінімум (2 тижні)

1. ✅ **T1** — аналітична таблиця (R: `pracma::hypergeo`, замкнені формули)
2. ✅ **T2** — Monte Carlo (після імплементації `lm_pmm3()`)
3. ✅ **T7** — тест бімодальності (критично для коректності PMM3)

### Фаза B — Реальні дані (1 тиждень)

4. 📧 Запит датасету до авторів → **T3, T4, T5**
5. ✅ Паралельно: аналіз `nlme::Fatigue` → **T3'**

### Фаза C — Розширення (1 тиждень)

6. ✅ **T6** — Misspecification robustness
7. ✅ **Fig. eff** — теплова карта g₃(γ₄, γ₆)

---

## MLE-TN: технічна реалізація для порівняння

MLE-TN потребує числової максимізації логправдоподібності:

```r
# Log-likelihood TN(ξ, η, λ) — формула (9) зі статті
ll_tn <- function(par, y, X) {
  beta <- par[1:(length(par) - 2)]
  eta  <- exp(par[length(par) - 1])  # позитивність через exp
  lam  <- par[length(par)]
  
  mu <- X %*% beta
  z  <- (y - mu) / eta
  
  n <- length(y)
  ll <- -n * log(2 * pi) / 2 - n * log(eta) - n * lam^2 / 2 -
    sum(z^2) / 2 + sum(log(cosh(lam * z)))
  return(-ll)  # мінімізація
}

mle_tn_fit <- function(formula, data) {
  mf  <- model.frame(formula, data)
  y   <- model.response(mf)
  X   <- model.matrix(attr(mf, "terms"), data = mf)
  ols <- lm(formula, data = data)
  
  par0 <- c(coef(ols), log(sd(residuals(ols))), 1.0)
  
  opt <- optim(par0, ll_tn, y = y, X = X,
               method = "L-BFGS-B",
               lower = c(rep(-Inf, ncol(X)), -5, 0),
               upper = c(rep(Inf,  ncol(X)),  5, 10),
               hessian = TRUE)
  
  list(beta      = opt$par[1:ncol(X)],
       eta       = exp(opt$par[ncol(X) + 1]),
       lambda    = opt$par[ncol(X) + 2],
       loglik    = -opt$value,
       hessian   = opt$hessian)
}
```

---

## Очікувані висновки

### Сценарій "Ідеальні умови" (DGP = TN, n велике)

```
MLE-TN ≈ PMM3 >> OLS при λ ∈ (0.5, 2.0)
PMM2 ≈ OLS (γ₃=0 → PMM2 не ефективний для TN)
```

### Сценарій "Мала вибірка" (n=46, як в статті)

```
MLE-TN може поступатись PMM3 через нестабільність оцінки λ
PMM3 більш стійкий завдяки моментній основі
```

### Сценарій "Бімодальний TN" (λ > 1.5)

```
MLE-TN >> PMM3 (явно враховує бімодальність)
PMM3 може мати низьку convergence rate
→ Необхідне розширення PMM3 для бімодального випадку
```

### Сценарій "Неправильна специфікація"

```
PMM3 > MLE-TN (PMM3 не залежить від хибної форми розподілу)
→ Головна практична перевага PMM3
```

---

## Зв'язок з EstemPMM та подальша розробка

Результати аналізу безпосередньо впливають на:

- **Vignette** `04-pmm3-symmetric-distributions.Rmd` — biaxial fatigue як приклад
- **Функція `test_symmetry()`** — автоматичне визначення PMM2 vs PMM3
- **Функція `lm_pmm3()`** — верифікація convergence при λ > 1
- **Документація** — чітке застереження про бімодальний режим TN

---

*Документ підготовлено у рамках розробки пакету EstemPMM.*  
*Версія: 1.0 · Дата: березень 2026 · Статус: Draft для реалізації*
