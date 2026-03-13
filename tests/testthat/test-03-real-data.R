# =============================================================================
# tests/testthat/test-03-real-data.R
# Unit tests for R/04_real_data.R
# =============================================================================

source("R/config.R")
source("R/01_tn_distribution.R")
source("R/04_real_data.R")

library(testthat)

# --- load_biaxial_fatigue() --------------------------------------------------
test_that("biaxial dataset: correct dimensions", {
  dat <- load_biaxial_fatigue()
  expect_equal(nrow(dat), 46)
  expect_equal(ncol(dat), 3)   # x, y, dataset
  expect_true(all(c("x", "y", "dataset") %in% names(dat)))
})

test_that("biaxial dataset: x values are 6 distinct stress levels", {
  dat <- load_biaxial_fatigue()
  expect_equal(length(unique(dat$x)), 6)
})

test_that("biaxial dataset: y values in plausible range (log cycles)", {
  dat <- load_biaxial_fatigue()
  expect_true(all(dat$y > 4 & dat$y < 7))
})

# --- resid_diagnostics() -----------------------------------------------------
test_that("resid_diagnostics: returns all required fields", {
  dat  <- load_biaxial_fatigue()
  diag <- resid_diagnostics(dat)
  required_fields <- c("n", "gamma3", "gamma4", "gamma6",
                       "g3_theoretical", "improvement_pct",
                       "is_symmetric", "is_gaussian", "recommended_method")
  expect_true(all(required_fields %in% names(diag)))
})

test_that("resid_diagnostics: biaxial residuals are platykurtic (gamma4 < 0)", {
  dat  <- load_biaxial_fatigue()
  diag <- resid_diagnostics(dat)
  # TN distribution is always platykurtic → gamma4 < 0
  expect_lt(diag$gamma4, 0)
})

test_that("resid_diagnostics: biaxial n matches dataset size", {
  dat  <- load_biaxial_fatigue()
  diag <- resid_diagnostics(dat)
  expect_equal(diag$n, 46)
})

test_that("resid_diagnostics: asymmetric data recommends PMM2", {
  # Biaxial has gamma3 ≈ -0.47 > threshold → should recommend PMM2
  dat  <- load_biaxial_fatigue()
  diag <- resid_diagnostics(dat)
  if (!diag$is_symmetric) {
    expect_match(diag$recommended_method, "PMM2", ignore.case = TRUE)
  }
})

test_that("resid_diagnostics: g3_theoretical in (0,1) for non-Gaussian data", {
  dat  <- load_biaxial_fatigue()
  diag <- resid_diagnostics(dat)
  if (!is.na(diag$g3_theoretical)) {
    expect_gt(diag$g3_theoretical, 0)
    expect_lt(diag$g3_theoretical, 1)
  }
})

# --- load_nlme_fatigue() (conditional) --------------------------------------
test_that("nlme::Fatigue loads and aggregates correctly", {
  skip_if_not_installed("nlme")
  dat <- load_nlme_fatigue()
  expect_true(all(c("x", "y") %in% names(dat)))
  expect_gt(nrow(dat), 5)
})
