# =============================================================================
# tests/testthat/test-02-mle-tn.R
# Unit tests for R/02_mle_tn.R
# =============================================================================

source("R/config.R")
source("R/01_tn_distribution.R")
source("R/02_mle_tn.R")

library(testthat)

# --- Basic convergence -------------------------------------------------------
test_that("mle_tn: converges on TN data (n=200, lambda=1.5)", {
  set.seed(123)
  x   <- rnorm(200)
  y   <- 1 + 2 * x + rtn(200, lambda = 1.5)
  fit <- mle_tn(y ~ x, data = data.frame(y, x))

  expect_true(fit$converged)
  expect_s3_class(fit, "mle_tn_fit")
  expect_equal(unname(fit$coefficients["x"]), 2, tolerance = 0.4)
})

test_that("mle_tn: recovers beta0 and beta1 on large sample", {
  set.seed(456)
  x   <- rnorm(500)
  y   <- 3 + 1.5 * x + rtn(500, lambda = 1.0)
  fit <- mle_tn(y ~ x, data = data.frame(y, x))

  expect_equal(unname(fit$coefficients[1]), 3,   tolerance = 0.3)
  expect_equal(unname(fit$coefficients[2]), 1.5, tolerance = 0.25)
})

# --- Lambda recovery --------------------------------------------------------
test_that("mle_tn: recovers lambda=2 on n=500", {
  set.seed(789)
  x   <- rnorm(500)
  y   <- 1 + 2 * x + rtn(500, lambda = 2.0)
  fit <- mle_tn(y ~ x, data = data.frame(y, x))

  expect_equal(fit$lambda, 2.0, tolerance = 0.6)
})

test_that("mle_tn: lambda ≈ 0 on Gaussian data", {
  set.seed(111)
  x   <- rnorm(300)
  y   <- 1 + 2 * x + rnorm(300)
  fit <- mle_tn(y ~ x, data = data.frame(y, x))

  expect_lt(abs(fit$lambda), 0.5)
})

# --- Information criteria ---------------------------------------------------
test_that("mle_tn: AIC/BIC are finite and AIC < BIC for n > 8", {
  set.seed(222)
  x   <- rnorm(100)
  y   <- 1 + 2 * x + rtn(100, lambda = 1.0)
  fit <- mle_tn(y ~ x, data = data.frame(y, x))

  expect_true(is.finite(fit$aic))
  expect_true(is.finite(fit$bic))
  # BIC penalises more for n=100 (log(100) = 4.6 > 2)
  expect_gt(fit$bic, fit$aic)
})

# --- coef / predict methods --------------------------------------------------
test_that("mle_tn: coef() returns named vector", {
  set.seed(333)
  x <- rnorm(100); y <- 1 + 2*x + rtn(100, lambda=1)
  fit <- mle_tn(y ~ x, data=data.frame(y,x))
  co  <- coef(fit)
  expect_named(co)
  expect_length(co, 2)
})

test_that("mle_tn: predict() returns fitted values without newdata", {
  set.seed(444)
  x <- rnorm(80); y <- 1 + 2*x + rtn(80, lambda=1)
  fit <- mle_tn(y ~ x, data=data.frame(y,x))
  expect_length(predict(fit), 80)
})
