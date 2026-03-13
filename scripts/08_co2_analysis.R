#!/usr/bin/env Rscript
# =============================================================================
# scripts/08_co2_analysis.R
# Real data experiment: Mauna Loa CO2 — PMM3 vs MLE-TN vs OLS vs PMM2
#
# Dataset: R built-in `co2`
#   n=468, monthly 1959–1997
#   Model: CO2 (ppm) ~ decimal year  [linear trend]
#   Residuals: seasonal oscillation (~sinusoid) → U-shaped → gamma4 << 0
#   Expected: PMM3 wins (symmetric, strongly platykurtic)
#
# Usage:
#   Rscript scripts/08_co2_analysis.R       (from project root)
#   source("scripts/08_co2_analysis.R")     (from RStudio)
# =============================================================================

# --- Setup -------------------------------------------------------------------
if (!exists("PROJECT_ROOT")) {
  # Prefer: already in project root (R/config.R exists = .Rprofile ran)
  if (file.exists("R/config.R")) {
    PROJECT_ROOT <- getwd()
  } else {
    # Fallback: derive from --file= argument passed by Rscript
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg)) {
      PROJECT_ROOT <- normalizePath(
        file.path(dirname(normalizePath(sub("--file=", "", file_arg))), ".."),
        mustWork = FALSE
      )
    } else {
      PROJECT_ROOT <- getwd()
    }
  }
}
setwd(PROJECT_ROOT)
cat("Working directory:", getwd(), "\n\n")

source("R/config.R")
load_project(verbose = FALSE)

CO2_OUT <- "results/co2"
dir.create(CO2_OUT, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load data ------------------------------------------------------------
cat("=== Mauna Loa CO2: Linear Trend Regression ===\n\n")
dat <- load_co2()
cat(sprintf("n = %d | x range: %.2f – %.2f | y range: %.1f – %.1f ppm\n\n",
            nrow(dat),
            min(dat$x), max(dat$x),
            min(dat$y), max(dat$y)))

# --- 2. OLS fit & residual inspection ----------------------------------------
cat("--- OLS fit ---\n")
ols_fit <- lm(y ~ x, data = dat)
cat(summary(ols_fit)$coefficients |> round(4) |> print(), "\n")
cat(sprintf("R² = %.4f\n\n", summary(ols_fit)$r.squared))

# --- 3. Residual diagnostics -------------------------------------------------
cat("--- Residual diagnostics ---\n")
diag <- resid_diagnostics(dat)
cat(sprintf(
  "  n         = %d\n  gamma3    = %+.4f  (%s)\n  gamma4    = %+.4f\n  gamma6    = %+.4f\n  g3_theor  = %.4f  (PMM3 improves %.1f%% vs OLS)\n  Recommend : %s\n\n",
  diag$n,
  diag$gamma3, if (diag$is_symmetric) "symmetric ✓" else "ASYMMETRIC ✗",
  diag$gamma4,
  diag$gamma6,
  ifelse(is.na(diag$g3_theoretical), NaN, diag$g3_theoretical),
  ifelse(is.na(diag$improvement_pct), NaN, diag$improvement_pct),
  diag$recommended_method
))

# Save diagnostics
saveRDS(diag, file.path(CO2_OUT, "co2_diagnostics.rds"))

# --- 4. Bootstrap comparison (B = 2000) --------------------------------------
cat("--- Bootstrap comparison (B=2000, seed=", CFG$seed, ") ---\n")
boot_res <- bootstrap_comparison(dat, B = 2000)

var_ols <- boot_res$sd[boot_res$method == "OLS"]^2
cat(sprintf("\n%-8s  mean=%7.4f  sd=%7.4f  g3_emp=%6.4f  NA=%.1f%%\n",
            boot_res$method, boot_res$mean, boot_res$sd,
            boot_res$g3_empirical, boot_res$na_pct) |> cat(sep=""))
cat("\n")

write.csv(boot_res,
          file.path(CO2_OUT, "co2_bootstrap.csv"),
          row.names = FALSE)
cat("Saved:", file.path(CO2_OUT, "co2_bootstrap.csv"), "\n\n")

# --- 5. LOO cross-validation -------------------------------------------------
cat("--- LOO cross-validation (n=", nrow(dat), ") ---\n")
# Note: use manual LOO for PMM2 to avoid generic predict dispatch bug
n <- nrow(dat)

loo_ols <- mean(sapply(seq_len(n), function(i) {
  fit <- lm(y ~ x, data = dat[-i, ])
  (dat$y[i] - predict(fit, newdata = dat[i, ]))^2
}))

loo_pmm3 <- mean(sapply(seq_len(n), function(i) {
  tryCatch({
    fit <- lm_pmm3(y ~ x, dat[-i, ])
    (dat$y[i] - predict(fit, newdata = dat[i, ]))^2
  }, error = function(e) NA_real_)
}), na.rm = TRUE)

loo_mltn <- mean(sapply(seq_len(n), function(i) {
  tryCatch({
    fit <- mle_tn(y ~ x, dat[-i, ])
    (dat$y[i] - predict(fit, newdata = dat[i, ]))^2
  }, error = function(e) NA_real_)
}), na.rm = TRUE)

loo_pmm2 <- mean(sapply(seq_len(n), function(i) {
  tryCatch({
    fit <- EstemPMM::lm_pmm2(y ~ x, dat[-i, ])
    pred <- as.numeric(predict(fit, newdata = dat[i, ]))
    (dat$y[i] - pred)^2
  }, error = function(e) NA_real_)
}), na.rm = TRUE)

loo_res <- data.frame(
  method  = c("OLS", "PMM2", "PMM3", "MLE_TN"),
  loo_mse = c(loo_ols, loo_pmm2, loo_pmm3, loo_mltn)
)
loo_res$rel_ols <- loo_res$loo_mse / loo_ols

cat("\nLOO-MSE results:\n")
print(loo_res, digits = 6)

write.csv(loo_res,
          file.path(CO2_OUT, "co2_loo_mse.csv"),
          row.names = FALSE)
cat("\nSaved:", file.path(CO2_OUT, "co2_loo_mse.csv"), "\n")

# --- 6. Summary --------------------------------------------------------------
cat("\n=== SUMMARY ===\n")
cat(sprintf("Dataset  : Mauna Loa CO2 (n=%d, y~x linear trend)\n", nrow(dat)))
cat(sprintf("gamma3   : %+.4f  → %s\n", diag$gamma3,
            if (diag$is_symmetric) "SYMMETRIC" else "asymmetric"))
cat(sprintf("gamma4   : %+.4f  → %s\n", diag$gamma4,
            if (diag$gamma4 < -0.5) "strongly platykurtic" else "mildly platykurtic"))
cat(sprintf("g3_theor : %.4f  → PMM3 expected to reduce slope var by %.1f%%\n",
            ifelse(is.na(diag$g3_theoretical), NaN, diag$g3_theoretical),
            ifelse(is.na(diag$improvement_pct), NaN, diag$improvement_pct)))
best_boot <- boot_res$method[which.min(boot_res$g3_empirical)]
cat(sprintf("Bootstrap winner (g3_emp): %s\n", best_boot))
best_loo  <- loo_res$method[which.min(loo_res$loo_mse)]
cat(sprintf("LOO-MSE winner           : %s\n", best_loo))
