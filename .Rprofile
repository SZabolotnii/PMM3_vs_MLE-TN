# =============================================================================
# .Rprofile — Автоматичне налаштування при відкритті проєкту в RStudio
# =============================================================================

local({
  cat("\n══════════════════════════════════════════\n")
  cat("  PMM3 vs MLE-TN Research Project\n")
  cat("  Salinas et al. (2023) | EstemPMM\n")
  cat("══════════════════════════════════════════\n\n")

  # Автоматичне встановлення залежностей якщо відсутні
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cat("NOTE: Some packages missing. Run:\n")
    cat("  source('scripts/install_deps.R')\n\n")
  }

  # Завантаження конфігурації
  if (file.exists("R/config.R")) {
    source("R/config.R")
    cat("Use load_project() to load all R modules.\n")
    cat("Use source('scripts/run_all.R') to run the full study.\n\n")
  }
})
