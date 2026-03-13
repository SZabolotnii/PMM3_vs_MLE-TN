# =============================================================================
# tests/testthat/test-01-tn-distribution.R
# Unit tests for R/01_tn_distribution.R
# =============================================================================

# Load modules (assumes tests run from PROJECT_ROOT)
source("R/config.R")
source("R/01_tn_distribution.R")

library(testthat)

# --- rtn() ------------------------------------------------------------------
test_that("rtn: returns correct number of observations", {
  expect_length(rtn(100, lambda = 1), 100)
  expect_length(rtn(1,   lambda = 0), 1)
})

test_that("rtn: produces symmetric distribution (gamma3 ≈ 0) for any lambda", {
  set.seed(42)
  for (lam in c(0.5, 1.0, 2.0, 3.0)) {
    x  <- rtn(20000, lambda = lam)
    g3 <- mean(x^3) / mean(x^2)^1.5
    expect_lt(abs(g3), 0.10,
              label = paste("gamma3 for lambda =", lam))
  }
})

test_that("rtn: lambda=0 produces standard normal", {
  set.seed(99)
  x <- rtn(10000, lambda = 0)
  expect_equal(mean(x), 0,  tolerance = 0.05)
  expect_equal(var(x),  1,  tolerance = 0.05)
  expect_equal(mean(x^4) / mean(x^2)^2 - 3, 0, tolerance = 0.15)  # gamma4 ≈ 0
})

test_that("rtn: validates inputs", {
  expect_error(rtn(100, eta = -1))
  expect_error(rtn(0))
})

# --- dtn() ------------------------------------------------------------------
test_that("dtn: density integrates to 1 for lambda in {0, 1, 2, 3}", {
  for (lam in c(0, 1, 2, 3)) {
    val <- integrate(dtn, -Inf, Inf, lambda = lam)$value
    expect_equal(val, 1, tolerance = 1e-4,
                 label = paste("integral for lambda =", lam))
  }
})

test_that("dtn: log argument works correctly", {
  x <- seq(-3, 3, by = 0.5)
  expect_equal(dtn(x, lambda = 1),
               exp(dtn(x, lambda = 1, log = TRUE)),
               tolerance = 1e-10)
})

test_that("dtn: symmetric around xi", {
  for (lam in c(0.5, 1.5)) {
    x <- seq(-3, 3, by = 0.1)
    expect_equal(dtn(x, lambda = lam),
                 dtn(-x, lambda = lam),
                 tolerance = 1e-10,
                 label = paste("symmetry for lambda =", lam))
  }
})

# --- tn_cumulants() ---------------------------------------------------------
test_that("tn_cumulants: g3 in (0,1) for lambda > 0", {
  for (lam in c(0.5, 1.0, 2.0, 5.0)) {
    cum <- tn_cumulants(lam)
    expect_true(cum$g3 > 0,
                label = paste("g3 > 0 for lambda =", lam))
    expect_true(cum$g3 < 1,
                label = paste("g3 < 1 for lambda =", lam))
  }
})

test_that("tn_cumulants: lambda=0 gives gaussian (gamma4=0, g3=1)", {
  cum <- tn_cumulants(0)
  expect_equal(cum$gamma4, 0, tolerance = 1e-5)
  expect_equal(cum$g3,     1, tolerance = 1e-4)
})

test_that("tn_cumulants: gamma4 is always negative for lambda > 0 (platykurtic)", {
  for (lam in c(0.5, 1.0, 2.0, 4.0)) {
    cum <- tn_cumulants(lam)
    expect_lt(cum$gamma4, 0,
              label = paste("gamma4 < 0 for lambda =", lam))
  }
})

test_that("tn_cumulants: bimodal flag correct (boundary at lambda=1)", {
  expect_false(tn_cumulants(0.5)$is_bimodal)
  expect_false(tn_cumulants(1.0)$is_bimodal)   # lambda=1 is unimodal (border)
  expect_true( tn_cumulants(1.5)$is_bimodal)
  expect_true( tn_cumulants(3.0)$is_bimodal)
})

test_that("tn_cumulants: g3 decreases as lambda increases (more non-Gaussian)", {
  g3_vals <- sapply(c(0.5, 1.0, 2.0, 3.0), function(l) tn_cumulants(l)$g3)
  expect_true(all(diff(g3_vals) < 0),
              label = "g3 is monotonically decreasing in lambda")
})

# --- build_table_T1() -------------------------------------------------------
test_that("build_table_T1: returns correct dimensions and columns", {
  T1 <- build_table_T1(lambda_grid = c(0.5, 1.0, 2.0))
  expect_equal(nrow(T1), 3)
  expect_true(all(c("lambda", "gamma4", "gamma6", "g3",
                    "improvement_pct", "mode") %in% names(T1)))
})

test_that("build_table_T1: mode column correct", {
  T1 <- build_table_T1(lambda_grid = c(0.5, 1.5))
  expect_equal(T1$mode[T1$lambda == 0.5], "unimodal")
  expect_equal(T1$mode[T1$lambda == 1.5], "bimodal")
})
