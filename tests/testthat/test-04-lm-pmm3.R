# =============================================================================
# tests/testthat/test-04-lm-pmm3.R
# Unit tests for R/pmm3_skeleton.R — real PMM3 Newton-Raphson implementation
# =============================================================================

library(testthat)

# --- Class and basic structure ------------------------------------------------
test_that("lm_pmm3: returns PMM3fit S4 object", {
  set.seed(1)
  x <- rnorm(100); y <- 1 + 2 * x + rtn(100, lambda = 2)
  fit <- lm_pmm3(y ~ x, data = data.frame(y, x))
  expect_s4_class(fit, "PMM3fit")
})

test_that("lm_pmm3: PMM3fit has required slots", {
  set.seed(2)
  x <- rnorm(100); y <- 1 + 2 * x + rtn(100, lambda = 2)
  fit <- lm_pmm3(y ~ x, data = data.frame(y, x))
  expect_true(isVirtualClass("PMM3fit") || is(fit, "PMM3fit"))
  expect_length(fit@coefficients, 2)
  expect_length(fit@fitted,       100)
  expect_length(fit@residuals,    100)
  expect_true(is.logical(fit@convergence))
  expect_true(is.integer(fit@iterations))
})

test_that("STUB is FALSE (real implementation active)", {
  expect_false(isTRUE(get("STUB", envir = globalenv())))
})

# --- Convergence and coefficient recovery ------------------------------------
test_that("lm_pmm3: converges on TN data (n=200, lambda=1.5)", {
  set.seed(10)
  x <- rnorm(200); y <- 1 + 2 * x + rtn(200, lambda = 1.5)
  fit <- lm_pmm3(y ~ x, data = data.frame(y, x))
  expect_true(fit@convergence)
  expect_gt(fit@iterations, 0L)   # at least one NR step taken
})

test_that("lm_pmm3: recovers beta0 and beta1 on large sample", {
  set.seed(20)
  x <- rnorm(500); y <- 3 + 1.5 * x + rtn(500, lambda = 1.0)
  fit <- lm_pmm3(y ~ x, data = data.frame(y, x))
  expect_equal(fit@coefficients[1], 3,   tolerance = 0.3)
  expect_equal(fit@coefficients[2], 1.5, tolerance = 0.25)
})

test_that("lm_pmm3: beats OLS variance at lambda=3 (ARE > 1)", {
  set.seed(30)
  n <- 300; M <- 200; lambda <- 3
  b1_ols <- b1_pmm3 <- numeric(M)
  for (i in seq_len(M)) {
    x <- rnorm(n); y <- 1 + 2 * x + rtn(n, lambda = lambda)
    dat <- data.frame(y, x)
    b1_ols[i]  <- coef(lm(y ~ x, dat))[2]
    b1_pmm3[i] <- lm_pmm3(y ~ x, data = dat)@coefficients[2]
  }
  are <- var(b1_ols) / var(b1_pmm3)
  expect_gt(are, 1.05,   # PMM3 should be notably more efficient
            label = sprintf("ARE=%.3f at lambda=%g", are, lambda))
})

# --- Near-Gaussian fallback --------------------------------------------------
test_that("lm_pmm3: near-Gaussian data (lambda=0) — coefficients close to OLS", {
  set.seed(40)
  x <- rnorm(300); y <- 1 + 2 * x + rnorm(300)
  fit_ols  <- lm(y ~ x, data = data.frame(y, x))
  fit_pmm3 <- lm_pmm3(y ~ x, data = data.frame(y, x))
  expect_equal(fit_pmm3@coefficients, unname(coef(fit_ols)), tolerance = 0.2)
})

# --- Moments and efficiency coefficient --------------------------------------
test_that("lm_pmm3: gamma4 < 0 for TN errors (platykurtic)", {
  set.seed(50)
  x <- rnorm(300); y <- 1 + 2 * x + rtn(300, lambda = 2.0)
  fit <- lm_pmm3(y ~ x, data = data.frame(y, x))
  expect_lt(fit@gamma4, 0)
})

test_that("lm_pmm3: g_coefficient in (0,1) for TN errors", {
  set.seed(60)
  x <- rnorm(300); y <- 1 + 2 * x + rtn(300, lambda = 2.0)
  fit <- lm_pmm3(y ~ x, data = data.frame(y, x))
  expect_gt(fit@g_coefficient, 0)
  expect_lt(fit@g_coefficient, 1)
})

test_that("lm_pmm3: g_coefficient approaches 1 for Gaussian data", {
  set.seed(70)
  x <- rnorm(500); y <- 1 + 2 * x + rnorm(500)
  fit <- lm_pmm3(y ~ x, data = data.frame(y, x))
  expect_gt(fit@g_coefficient, 0.85)
})

# --- S4 methods ---------------------------------------------------------------
test_that("coef(PMM3fit) returns named numeric vector", {
  set.seed(80)
  x <- rnorm(100); y <- 1 + 2 * x + rtn(100, lambda = 1)
  fit <- lm_pmm3(y ~ x, data = data.frame(y, x))
  co  <- coef(fit)
  expect_true(is.numeric(co))
  expect_length(co, 2)
  expect_named(co)
})

test_that("predict(PMM3fit) without newdata returns fitted values", {
  set.seed(90)
  x <- rnorm(80); y <- 1 + 2 * x + rtn(80, lambda = 1)
  fit <- lm_pmm3(y ~ x, data = data.frame(y, x))
  expect_equal(predict(fit), fit@fitted)
  expect_length(predict(fit), 80)
})

test_that("predict(PMM3fit) with newdata produces out-of-sample predictions", {
  set.seed(100)
  x    <- rnorm(100); y <- 1 + 2 * x + rtn(100, lambda = 1)
  fit  <- lm_pmm3(y ~ x, data = data.frame(y, x))
  nd   <- data.frame(x = c(-1, 0, 1))
  pred <- predict(fit, newdata = nd)
  expect_length(pred, 3)
  expect_true(all(is.finite(pred)))
  # Linear prediction: should be approximately 1 + 2*x
  expect_equal(pred[2], fit@coefficients[1], tolerance = 0.5)
})

test_that("residuals(PMM3fit) sums near zero (unbiased)", {
  set.seed(110)
  x <- rnorm(200); y <- 1 + 2 * x + rtn(200, lambda = 2)
  fit <- lm_pmm3(y ~ x, data = data.frame(y, x))
  expect_equal(mean(residuals(fit)), 0, tolerance = 0.1)
})

# --- Multivariate regression -------------------------------------------------
test_that("lm_pmm3: works with two predictors", {
  set.seed(120)
  x1 <- rnorm(200); x2 <- rnorm(200)
  y  <- 1 + 2 * x1 + (-1) * x2 + rtn(200, lambda = 1.5)
  fit <- lm_pmm3(y ~ x1 + x2, data = data.frame(y, x1, x2))
  expect_s4_class(fit, "PMM3fit")
  expect_length(fit@coefficients, 3)
  expect_equal(fit@coefficients[2],  2,  tolerance = 0.4)
  expect_equal(fit@coefficients[3], -1,  tolerance = 0.4)
})

# --- LOO-CV compatibility (predict with single-row newdata) ------------------
test_that("lm_pmm3: predict works for single-row newdata (LOO-CV pattern)", {
  set.seed(130)
  n  <- 50
  x  <- rnorm(n); y <- 1 + 2 * x + rtn(n, lambda = 1.5)
  dat <- data.frame(y, x)
  fit <- lm_pmm3(y ~ x, data = dat[-1, ])
  pred <- predict(fit, newdata = dat[1, ])
  expect_length(pred, 1)
  expect_true(is.finite(pred))
})
