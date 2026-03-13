#!/usr/bin/env Rscript
# =============================================================================
# scripts/run_all.R — Повний запуск дослідження PMM3 vs MLE-TN
#
# Використання:
#   Rscript scripts/run_all.R          (з кореневої директорії проєкту)
#   source("scripts/run_all.R")        (з RStudio)
#
# Час виконання (орієнтовно):
#   CFG$mc_replications = 500   → ~2-5 хв
#   CFG$mc_replications = 5000  → ~30-60 хв
# =============================================================================

t_start <- proc.time()
cat("╔══════════════════════════════════════╗\n")
cat("║  PMM3 vs MLE-TN Comparative Study   ║\n")
cat(sprintf("║  Started: %-26s ║\n", format(Sys.time(), "%Y-%m-%d %H:%M")))
cat("╚══════════════════════════════════════╝\n\n")

# --- Setup -------------------------------------------------------------------
# Set working directory to project root (handles both source() and Rscript)
if (!exists("PROJECT_ROOT")) {
  script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) getwd()
  )
  PROJECT_ROOT <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
}
setwd(PROJECT_ROOT)
cat("Working directory:", getwd(), "\n\n")

source("R/config.R")
load_project(verbose = TRUE)

# --- Block 1: Theoretical Calibration ----------------------------------------
cat("\n--- Block 1: Table T1 (Theoretical) ---\n")
T1 <- build_table_T1()
out_path <- file.path(CFG$paths$results_tables, "T1_theoretical.csv")
write.csv(T1, out_path, row.names = FALSE)
cat("Saved:", out_path, "\n")
print(T1[, c("lambda", "gamma4", "g3", "improvement_pct", "mode")])

# --- Block 2: Monte Carlo ----------------------------------------------------
cat("\n--- Block 2: Monte Carlo Study ---\n")
cat(sprintf("Config: M=%d, lambda=%s, n=%s\n",
            CFG$mc_replications,
            paste(CFG$lambda_grid, collapse=","),
            paste(CFG$sample_sizes, collapse=",")))
mc_results <- run_full_mc_study()
out_path <- file.path(CFG$paths$results_tables, "T2_monte_carlo.csv")
write.csv(mc_results, out_path, row.names = FALSE)
cat("Saved:", out_path, "\n")

# --- Block 3: Real Data -------------------------------------------------------
cat("\n--- Block 3: Biaxial Fatigue Diagnostics ---\n")
dat_bf  <- load_biaxial_fatigue()
diag_bf <- resid_diagnostics(dat_bf)
cat(sprintf("  n=%d | gamma3=%.3f | gamma4=%.3f | g3=%.4f\n",
            diag_bf$n, diag_bf$gamma3, diag_bf$gamma4,
            ifelse(is.na(diag_bf$g3_theoretical), NaN, diag_bf$g3_theoretical)))
cat("  Recommended method:", diag_bf$recommended_method, "\n")
saveRDS(diag_bf, file.path(CFG$paths$results_tables, "biaxial_diagnostics.rds"))

# --- Block 5: Misspecification ------------------------------------------------
cat("\n--- Block 5: Misspecification Robustness (M=500) ---\n")
misspec <- do.call(rbind,
  lapply(CFG$misspec_distributions, function(d) {
    cat("  Distribution:", d, "\n")
    run_misspec_block(d, M = 500)
  })
)
out_path <- file.path(CFG$paths$results_tables, "T6_misspecification.csv")
write.csv(misspec, out_path, row.names = FALSE)
cat("Saved:", out_path, "\n")
print(misspec[, c("distribution","are_mle","are_pmm3","bias_pmm3")])

# --- Block 6: Bimodal Test ---------------------------------------------------
cat("\n--- Block 6: Bimodal Convergence Test (M=200) ---\n")
bimodal <- run_bimodal_test(M = 200)
out_path <- file.path(CFG$paths$results_tables, "T7_bimodal_test.csv")
write.csv(bimodal, out_path, row.names = FALSE)
cat("Saved:", out_path, "\n")
print(bimodal)

# --- Figures -----------------------------------------------------------------
cat("\n--- Saving Figures ---\n")
tryCatch({
  save_plot(plot_efficiency_curve(T1),       "fig1_efficiency_curve.png")
  save_plot(plot_are_comparison(mc_results), "fig2_are_comparison.png")
  save_plot(plot_g3_convergence(mc_results, T1), "fig3_g3_convergence.png")
  save_plot(plot_bimodal_convergence(bimodal),    "fig4_bimodal.png")
  cat("All figures saved.\n")
}, error = function(e) {
  cat("Figure error (ggplot2 installed?): ", conditionMessage(e), "\n")
})

# --- Render Report -----------------------------------------------------------
cat("\n--- Rendering HTML Report ---\n")
tryCatch({
  rmarkdown::render(
    "reports/main_report.Rmd",
    output_file = "../results/main_report.html",
    quiet       = TRUE
  )
  cat("Report: results/main_report.html\n")
}, error = function(e) {
  cat("Render error (rmarkdown installed?): ", conditionMessage(e), "\n")
})

# --- Summary -----------------------------------------------------------------
elapsed <- round((proc.time() - t_start)["elapsed"] / 60, 1)
cat("\n╔══════════════════════════════════════╗\n")
cat(sprintf("║  Study complete in %.1f min          ║\n", elapsed))
cat(sprintf("║  Results: %-27s ║\n", "results/tables/"))
cat("╚══════════════════════════════════════╝\n")
